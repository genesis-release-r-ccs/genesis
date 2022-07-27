!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_pme_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_pme_mod

  use sp_energy_pme_opt_1dalltoall_mod
  use sp_energy_pme_noopt_1dalltoall_mod
  use sp_energy_pme_opt_2dalltoall_mod
  use sp_energy_pme_noopt_2dalltoall_mod
  use sp_energy_pme_opt_slab_mod
  use sp_energy_pme_noopt_slab_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  real(dp), public, save :: u_self       ! Ewald self energy
  integer,  public, save :: pme_scheme
  integer,  public, save :: pme_nspline

  ! constants
  integer,  parameter    :: CheckLoopCount = 100

  ! subroutines
  public  :: setup_pme
  public  :: dealloc_pme
  public  :: pme_pre
  public  :: pme_pre_lj
  public  :: pme_recip
  public  :: pme_recip_lj
  private :: select_pme_scheme

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme
  !> @brief        Setup for PME with domain decomposition
  !! @authors      JJ, NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @param[inout] enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme(domain, boundary, enefunc)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(in)    :: enefunc

    logical                  :: calc
    integer                  :: iproc(3), index_x, index_y, index_z


    if (boundary%type /= BoundaryTypePBC) &
      call error_msg('Setup_PME> Error! BoundaryType is not PBC.')

    pme_scheme = enefunc%pme_scheme
    pme_nspline = enefunc%pme_nspline

    ! decide communicator
    !
    iproc(1) = mod(my_city_rank, boundary%num_domain(1))
    iproc(2) = mod(my_city_rank/boundary%num_domain(1),boundary%num_domain(2))
    iproc(3) = my_city_rank/(boundary%num_domain(1)*boundary%num_domain(2))

    index_x  = iproc(2)*nproc_city + iproc(3)
    index_y  = iproc(3)*nproc_city + iproc(1)
    index_z  = iproc(1)*nproc_city + iproc(2)

#ifdef HAVE_MPI_GENESIS
    call mpi_comm_split(mpi_comm_city, index_x,my_city_rank,grid_commx,ierror)
    call mpi_comm_size (grid_commx, nprocx, ierror)
    call mpi_comm_rank (grid_commx, my_x_rank, ierror)
    call mpi_comm_split(mpi_comm_city, index_y,my_city_rank,grid_commy,ierror)
    call mpi_comm_size (grid_commy, nprocy, ierror)
    call mpi_comm_rank (grid_commy, my_y_rank, ierror)
    call mpi_comm_split(mpi_comm_city, index_z,my_city_rank,grid_commz,ierror)
    call mpi_comm_size (grid_commz, nprocz, ierror)
    call mpi_comm_rank (grid_commz, my_z_rank, ierror)
    call mpi_comm_split(mpi_comm_city, iproc(3), my_city_rank, &
                        grid_commxy, ierror)
    call mpi_comm_size(grid_commxy, nprocxy, ierror)
    call mpi_comm_rank(grid_commxy, my_xy_rank, ierror)
#endif

    select case (pme_scheme)

    case (FFT_AutoSelect)
      call select_pme_scheme(domain, boundary, enefunc)

    case (FFT_opt_1dalltoall)
      call setup_pme_opt_1dalltoall(domain, boundary, enefunc, calc)

    case (FFT_noopt_1dalltoall)
      call setup_pme_noopt_1dalltoall(domain, boundary, enefunc)

    case (FFT_opt_2dalltoall)
      call setup_pme_opt_2dalltoall(domain, boundary, enefunc, calc)

    case (FFT_noopt_2dalltoall)
      call setup_pme_noopt_2dalltoall(domain, boundary, enefunc)

    case (FFT_opt_slab)
      call setup_pme_opt_slab(domain, boundary, enefunc, calc)

    case (FFT_noopt_slab)
      call setup_pme_noopt_slab(domain, boundary, enefunc, calc)

    end select

    return

  end subroutine setup_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pme
  !> @brief        deallocate all arrays used in PME
  !! @authors      JJ
  !! @param[in]    enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pme(enefunc)

    type(s_enefunc),         intent(in) :: enefunc


    select case (pme_scheme)

    case (FFT_AutoSelect)
      return

    case (FFT_opt_1dalltoall)
      call dealloc_pme_opt_1dalltoall(enefunc)

    case (FFT_noopt_1dalltoall)
      call dealloc_pme_noopt_1dalltoall(enefunc)

    case (FFT_opt_2dalltoall)
      call dealloc_pme_opt_2dalltoall(enefunc)

    case (FFT_noopt_2dalltoall)
      call dealloc_pme_noopt_2dalltoall(enefunc)

    case (FFT_opt_slab)
      call dealloc_pme_opt_slab

    case (FFT_noopt_slab)
      call dealloc_pme_noopt_slab

    end select

    return

  end subroutine dealloc_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre
  !> @brief        Prepare functions for PME calculation
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary


    select case (pme_scheme)

    case (FFT_AutoSelect)
      return

    case (FFT_opt_1dalltoall)
      call pme_pre_opt_1dalltoall(domain, boundary)
      u_self = u_self_o1

    case (FFT_noopt_1dalltoall)
      call pme_pre_noopt_1dalltoall(domain, boundary)
      u_self = u_self_no1

    case (FFT_opt_2dalltoall)
      call pme_pre_opt_2dalltoall(domain, boundary)
      u_self = u_self_o2

    case (FFT_noopt_2dalltoall)
      call pme_pre_noopt_2dalltoall(domain, boundary)
      u_self = u_self_no2

    case (FFT_opt_slab)
      call pme_pre_opt_slab(domain, boundary)
      u_self = u_self_os

    case (FFT_noopt_slab)
      call pme_pre_noopt_slab(domain, boundary)
      u_self = u_self_nos

    end select

    return

  end subroutine pme_pre

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_lj
  !> @brief        Prepare functions for PME calculation (with LJ)
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_lj(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary


    select case (pme_scheme)

    case (FFT_AutoSelect)
      return

    case (FFT_opt_1dalltoall)
      call pme_pre_opt_1dalltoall_lj(domain, boundary)
      u_self = u_self_o1

    case (FFT_noopt_1dalltoall)
      call pme_pre_noopt_1dalltoall_lj(domain, boundary)
      u_self = u_self_no1

    case (FFT_opt_2dalltoall)
      call pme_pre_opt_2dalltoall_lj(domain, boundary)
      u_self = u_self_o2

    case (FFT_noopt_2dalltoall)
      call pme_pre_noopt_2dalltoall_lj(domain, boundary)
      u_self = u_self_no2

    case (FFT_opt_slab)
      call error_msg('Dispersion_PME with opt_slab is not available now')
!     call pme_pre_opt_slab_lj(domain, boundary)
!     u_self = u_self_os

    case (FFT_noopt_slab)
      call error_msg('Dispersion_PME with noopt_slab is not available now')
!     call pme_pre_noopt_slab_lj(domain, boundary)
!     u_self = u_self_nos

    end select

    return

  end subroutine pme_pre_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)


    select case (pme_scheme)

    case (FFT_AutoSelect)
      return

    case (FFT_opt_1dalltoall)
      if (pme_nspline == 4) then
        call pme_recip_opt_1dalltoall_4(domain, force, virial, eelec)
      else if (pme_nspline == 6) then
        call pme_recip_opt_1dalltoall_6(domain, force, virial, eelec)
      else if (pme_nspline == 8) then
        call pme_recip_opt_1dalltoall_8(domain, force, virial, eelec)
      end if

    case (FFT_noopt_1dalltoall)
      call pme_recip_noopt_1dalltoall(domain, force, virial, eelec)

    case (FFT_opt_2dalltoall)
      if (pme_nspline == 4) then
        call pme_recip_opt_2dalltoall_4(domain, force, virial, eelec)
      else if (pme_nspline == 6) then
        call pme_recip_opt_2dalltoall_6(domain, force, virial, eelec)
      else if (pme_nspline == 8) then
        call pme_recip_opt_2dalltoall_8(domain, force, virial, eelec)
      end if

    case (FFT_noopt_2dalltoall)
      call pme_recip_noopt_2dalltoall(domain, force, virial, eelec)

    case (FFT_opt_slab)
      if (pme_nspline == 4) then
        call pme_recip_opt_slab_4(domain, force, virial, eelec)
      else if (pme_nspline == 6) then
        call pme_recip_opt_slab_6(domain, force, virial, eelec)
      else if (pme_nspline == 8) then
        call pme_recip_opt_slab_8(domain, force, virial, eelec)
      end if

    case (FFT_noopt_slab)
      call pme_recip_noopt_slab(domain, force, virial, eelec)

    end select

    return

  end subroutine pme_recip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_lj
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : energy function
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_lj(domain, enefunc, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)


    select case (pme_scheme)

    case (FFT_AutoSelect)
      return

    case (FFT_opt_1dalltoall)
      if (pme_nspline == 4) then
        call pme_recip_opt_1dalltoall_lj_4 &
                              (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 6) then
        call pme_recip_opt_1dalltoall_lj_6 &
                              (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 8) then
        call pme_recip_opt_1dalltoall_lj_8 &
                              (domain, enefunc, force, virial, eelec, evdw)
      end if

    case (FFT_noopt_1dalltoall)
      call pme_recip_noopt_1dalltoall_lj &
                              (domain, enefunc, force, virial, eelec, evdw)

    case (FFT_opt_2dalltoall)
      if (pme_nspline == 4) then
        call pme_recip_opt_2dalltoall_lj_4 &
                              (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 6) then
        call pme_recip_opt_2dalltoall_lj_6 &
                              (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 8) then
!       call pme_recip_opt_2dalltoall_lj_8 &
!                             (domain, enefunc, force, virial, eelec, evdw)
      end if

    case (FFT_noopt_2dalltoall)
      call pme_recip_noopt_2dalltoall_lj &
                              (domain, enefunc, force, virial, eelec, evdw)

    case (FFT_opt_slab)
      if (pme_nspline == 4) then
!       call pme_recip_opt_slab_lj_4 &
!                             (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 6) then
!       call pme_recip_opt_slab_lj_6 &
!                             (domain, enefunc, force, virial, eelec, evdw)
      else if (pme_nspline == 8) then
!       call pme_recip_opt_slab_lj_8 &
!                             (domain, enefunc, force, virial, eelec, evdw)
      end if

    case (FFT_noopt_slab)
!     call pme_recip_noopt_slab_lj &
!                             (domain, enefunc, force, virial, eelec, evdw)

    end select

    return

  end subroutine pme_recip_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_pme_scheme
  !> @brief        Auto selection of FFT scheme
  !! @authors      NT
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @param[in]    enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_pme_scheme(domain, boundary, enefunc)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(in)    :: enefunc

    ! local variables
    real(wip),   allocatable :: force(:,:,:)
    real(dp),    allocatable :: virial(:,:,:), eelec(:)
    real(wp)                 :: grid_space
    real(dp)                 :: proc_time(7), start_time, minimum_time
    logical                  :: calc_opt, calc_noopt, calc_pme_scheme(7)
    integer                  :: i, ncell, natom, ngridmax
    integer                  :: ngrid(3)
    integer                  :: nlocalx, nlocaly, nlocalz


    if (main_rank) then
      write(MsgOut,'(a)') 'Select_FFT_Scheme> Checking performance of long range interaction operation...'
      write(MsgOut,'(a)') ''
    end if

    ! initialize variables
    !

    proc_time(1:7) = 0.0_dp

    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    allocate(force(3,natom,ncell), &
             virial(3,3,nthread),  &
             eelec(nthread))

    force(1:3,1:natom,1:ncell) = 0.0_wip
    virial(1:3,1:3,1:nthread)  = 0.0_dp
    eelec(1:nthread)           = 0.0_dp

    ngrid(1) = enefunc%pme_ngrid_x
    ngrid(2) = enefunc%pme_ngrid_y
    ngrid(3) = enefunc%pme_ngrid_z

    if (ngrid(1) == 0 .and. ngrid(2) == 0 .and. ngrid(3) ==0) then
      grid_space = enefunc%pme_max_spacing
      ngrid(1) = int(boundary%box_size_x/grid_space)
      ngrid(2) = int(boundary%box_size_y/grid_space)
      ngrid(3) = int(boundary%box_size_z/grid_space)
      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))
    end if

    calc_opt = .true.
    calc_noopt = .true. 
    calc_pme_scheme(1:7) = .true.

    ! Stop the calculation if the pme spline order is not 4, 6, nor 8.
    !
    if (pme_nspline /= 4 .and. pme_nspline /= 6 .and. &
        pme_nspline /= 8)                               &
      call error_msg('PME_nspline should be 4, 6, or 8')

    ! Check if charges in boundary cells are sufficient for charge
    ! grid data in each subdomain (noopt case)
    !
    grid_space = boundary%box_size_x / real(ngrid(1),wp)
!   ngridmax   = int(boundary%cell_size_x / grid_space)
    ngridmax   = ngrid(1) / boundary%num_domain(1)
    grid_space = boundary%box_size_y / real(ngrid(2),wp)
!   ngridmax   = min(ngridmax,int(boundary%cell_size_y/grid_space))
    ngridmax   = min(ngridmax, ngrid(2) / boundary%num_domain(2))
    grid_space = boundary%box_size_z / real(ngrid(3),wp)
!   ngridmax   = min(ngridmax,int(boundary%cell_size_z/grid_space))
    ngridmax   = min(ngridmax, ngrid(3) / boundary%num_domain(3))
    if (ngridmax < pme_nspline) then
      calc_pme_scheme(FFT_noopt_1dalltoall) = .false.
      calc_pme_scheme(FFT_noopt_2dalltoall) = .false.
      calc_pme_scheme(FFT_noopt_slab) = .false.
      calc_noopt = .false.
    end if

    ! Check if grid numbers in each subdomain is sufficent for communication
    ! (opt case)
    !
    ngridmax = ngrid(1) / boundary%num_domain(1) 
    ngridmax = min(ngridmax, ngrid(2)/boundary%num_domain(2))
    ngridmax = min(ngridmax, ngrid(3)/boundary%num_domain(3))
    if (ngridmax < 2*pme_nspline) then
      calc_pme_scheme(FFT_opt_1dalltoall) = .false.
      calc_pme_scheme(FFT_opt_2dalltoall) = .false.
      calc_pme_scheme(FFT_opt_slab) = .false.
      calc_opt = .false.
    end if

    ! Stop the simulation if all 6 cases fail 
    if (.not. calc_opt .and. .not. calc_noopt) &
      call error_msg('Current domain information and grid data ' // &
                     'is not recommended. Please increase the ' // &
                     'PME_ngrid_x/PME_ngrid_y/PME_ngrid_z')


    ! check performance
    !

    ! opt-1dalltoall
    !
    if (calc_opt) then

      call setup_pme_opt_1dalltoall(domain, boundary, enefunc, &
                                    calc_pme_scheme(FFT_opt_1dalltoall))
 
      if (calc_pme_scheme(FFT_opt_1dalltoall)) then

        start_time = get_unix_time()

        do i = 1, CheckLoopCount
          call pme_pre_opt_1dalltoall(domain, boundary)
          if (pme_nspline == 4) then
            call pme_recip_opt_1dalltoall_4(domain, force, virial, eelec)
          else if (pme_nspline == 6) then
            call pme_recip_opt_1dalltoall_6(domain, force, virial, eelec)
          else if (pme_nspline == 8) then
            call pme_recip_opt_1dalltoall_8(domain, force, virial, eelec)
          end if
        end do

        proc_time(FFT_opt_1dalltoall) = (get_unix_time() - start_time)*1.0d-3
        if (main_rank) &
          write(MsgOut,'(a,f12.3,a)') '  OPT_1DALLTOALL   : ', &
                                  proc_time(FFT_opt_1dalltoall), ' (ms)'

      end if

    end if

    ! noopt-1dalltoall
    !
    if (calc_noopt) then

      call setup_pme_noopt_1dalltoall(domain, boundary, enefunc)

      start_time = get_unix_time()
  
      do i = 1, CheckLoopCount
        call pme_pre_noopt_1dalltoall(domain, boundary)
        call pme_recip_noopt_1dalltoall(domain, force, virial, eelec)
      end do

      proc_time(FFT_noopt_1dalltoall) = (get_unix_time() - start_time)*1.0d-3

      if (main_rank) &
        write(MsgOut,'(a,f12.3,a)') '  NOOPT-1DALLTOALL : ', &
                                proc_time(FFT_noopt_1dalltoall), ' (ms)'
  
    end if

    ! opt-2dalltoall
    !
    if (calc_opt) then

      call setup_pme_opt_2dalltoall(domain, boundary, enefunc, &
                                    calc_pme_scheme(FFT_opt_2dalltoall))

      if (calc_pme_scheme(FFT_opt_2dalltoall)) then

        start_time = get_unix_time()
  
        do i = 1, CheckLoopCount
          call pme_pre_opt_2dalltoall(domain, boundary)
          if (pme_nspline == 4) then
            call pme_recip_opt_2dalltoall_4(domain, force, virial, eelec)
          else if (pme_nspline == 6) then
            call pme_recip_opt_2dalltoall_6(domain, force, virial, eelec)
          else if (pme_nspline == 8) then
            call pme_recip_opt_2dalltoall_8(domain, force, virial, eelec)
          end if
        end do
  
        proc_time(FFT_opt_2dalltoall) = (get_unix_time() - start_time)*1.0d-3
    
        if (main_rank) &
          write(MsgOut,'(a,f12.3,a)') '  OPT-2DALLTOALL   : ', &
                                  proc_time(FFT_opt_2dalltoall), ' (ms)'
  
      end if
  
    end if

    ! noopt-2dalltoall
    !
    if (calc_noopt) then

      call setup_pme_noopt_2dalltoall(domain, boundary, enefunc)

      start_time = get_unix_time()
  
      do i = 1, CheckLoopCount
        call pme_pre_noopt_2dalltoall(domain, boundary)
        call pme_recip_noopt_2dalltoall(domain, force, virial, eelec)
      end do
  
      proc_time(FFT_noopt_2dalltoall) = (get_unix_time() - start_time)*1.0d-3
    
      if (main_rank) &
        write(MsgOut,'(a,f12.3,a)') '  NOOPT-2DALLTOALL : ', &
                                proc_time(FFT_noopt_2dalltoall), ' (ms)'
 
    end if

!   ! opt-slab
!   !
!   if (calc_opt) then

!     call setup_pme_opt_slab(domain, boundary, enefunc, &
!                             calc_pme_scheme(FFT_opt_slab))

!     if (calc_pme_scheme(FFT_opt_slab)) then

!       start_time = get_unix_time()

!       do i = 1, CheckLoopCount
!         call pme_pre_opt_slab(domain, boundary)
!         if (pme_nspline == 4) then
!           call pme_recip_opt_slab_4(domain, force, virial, eelec)
!         else if (pme_nspline == 6) then
!           call pme_recip_opt_slab_6(domain, force, virial, eelec)
!         else if (pme_nspline == 8) then
!           call pme_recip_opt_slab_8(domain, force, virial, eelec)
!         end if
!       end do

!       proc_time(FFT_opt_slab) = (get_unix_time() - start_time)*1.0d-3
! 
!       if (main_rank) &
!         write(MsgOut,'(a,f12.3,a)') '  OPT-SLAB         : ', &
!                               proc_time(FFT_opt_slab), ' (ms)'

!     else

!       if (main_rank) &
!         write(MsgOut,'(a)') '  OPT-SLAB is skipped '

!     end if

!   end if

!   ! noopt-slab
!   !
!   if (calc_noopt) then

!     call setup_pme_noopt_slab(domain, boundary, enefunc, &
!                               calc_pme_scheme(FFT_noopt_slab))

!     if (calc_pme_scheme(FFT_noopt_slab)) then

!       start_time = get_unix_time()
! 
!       do i = 1, CheckLoopCount
!         call pme_pre_noopt_slab(domain, boundary)
!         call pme_recip_noopt_slab(domain, force, virial, eelec)
!       end do

!       proc_time(FFT_noopt_slab) = (get_unix_time() - start_time)*1.0d-3
! 
!       if (main_rank) &
!         write(MsgOut,'(a,f12.3,a)') '  NOOPT-SLAB       : ', &
!                                 proc_time(FFT_noopt_slab), ' (ms)'
!
!     else 

!       if (main_rank) &
!         write(MsgOut,'(a)') '  NOOPT-SLAB is skipped '
!     end if

!   end if 

    deallocate(force, virial, eelec)

    ! select scheme
    !

    minimum_time = 1000000000.0_wp


    if (main_rank) then
      do i = FFT_opt_1dalltoall, FFT_noopt_2dalltoall
        if (proc_time(i) < minimum_time .and. calc_pme_scheme(i)) then
          pme_scheme = i
          minimum_time = proc_time(i)
        end if
      end do
    end if
#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(pme_scheme, 1, mpi_integer, 0, mpi_comm_world, ierror)
#endif

    if (main_rank) then
      write(MsgOut,'(a)') ''
      write(MsgOut,'(a,a)') 'Select_FFT_Scheme> selected scheme is ', &
           trim(FFT_Types(pme_scheme))
    end if

    select case (pme_scheme)

    case (FFT_opt_1dalltoall)
      call dealloc_pme_noopt_1dalltoall(enefunc)
      call dealloc_pme_opt_2dalltoall(enefunc)
      call dealloc_pme_noopt_2dalltoall(enefunc)
      call dealloc_pme_opt_slab
      call dealloc_pme_noopt_slab

    case (FFT_noopt_1dalltoall)
      call dealloc_pme_opt_1dalltoall(enefunc)
      call dealloc_pme_opt_2dalltoall(enefunc)
      call dealloc_pme_noopt_2dalltoall(enefunc)
      call dealloc_pme_opt_slab
      call dealloc_pme_noopt_slab

    case (FFT_opt_2dalltoall)
      call dealloc_pme_opt_1dalltoall(enefunc)
      call dealloc_pme_noopt_1dalltoall(enefunc)
      call dealloc_pme_noopt_2dalltoall(enefunc)
      call dealloc_pme_opt_slab
      call dealloc_pme_noopt_slab

    case (FFT_noopt_2dalltoall)
      call dealloc_pme_opt_1dalltoall(enefunc)
      call dealloc_pme_noopt_1dalltoall(enefunc)
      call dealloc_pme_opt_2dalltoall(enefunc)
      call dealloc_pme_opt_slab
      call dealloc_pme_noopt_slab

    case (FFT_opt_slab)
      call dealloc_pme_opt_1dalltoall(enefunc)
      call dealloc_pme_noopt_1dalltoall(enefunc)
      call dealloc_pme_opt_2dalltoall(enefunc)
      call dealloc_pme_noopt_2dalltoall(enefunc)
      call dealloc_pme_noopt_slab

    case (FFT_noopt_slab)
      call dealloc_pme_opt_1dalltoall(enefunc)
      call dealloc_pme_noopt_1dalltoall(enefunc)
      call dealloc_pme_opt_2dalltoall(enefunc)
      call dealloc_pme_noopt_2dalltoall(enefunc)
      call dealloc_pme_opt_slab

    end select

    return

  end subroutine select_pme_scheme

end module sp_energy_pme_mod
