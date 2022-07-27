!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_pme_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Takashi Imai (TI), Jaewoon Jung (JJ), Takaharu Mori (TM), 
!!          Chigusa Kobayashi (CK)
!! @note    TI modified for NPT (2010/11/23)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_pme_mod

  use at_energy_str_mod
  use at_energy_pme_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use at_energy_pme_mod
  use at_energy_mod
  use math_libs_mod
  use molecules_str_mod
  use at_boundary_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters
  integer, parameter      :: NumIndex        = 45
  integer, parameter      :: Index(NumIndex) = (/ &
                              1,   2,   3,   4,   5,   6,   8,   9,  10,  12, &
                             15,  16,  18,  20,  24,  25,  27,  30,  32,  36, &
                             40,  45,  48,  50,  54,  60,  64,  72,  75,  80, &
                             81,  90,  96, 100, 120, 125, 128, 135, 144, 150, &
                            160, 162, 180, 192, 200/)

  ! subroutines
  public  :: define_enefunc_pme
  private :: setup_pme_grid_number
  private :: setup_pme_reciprocal_part
  private :: setup_pme_real_part
  private :: setup_self_energy

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_pme
  !> @brief        Setup for PME
  !! @authors      TM
  !! @param[in]    ene_info : parameters in [ENERGY] section
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                :: ny, nz, nx1, ny1


    if (ene_info%electrostatic == ElectrostaticPME) then
      enefunc%pme_use = .true.
      if (boundary%type /= BoundaryTypePBC) &
        call error_msg('Define_Enefunc_Pme> PME is available under PBC')
    else
      enefunc%pme_use = .false.
      return
    end if

    ! Set common parameters
    !
    enefunc%pme%n_bs            =   ene_info%pme_nspline
    enefunc%pme%alpha           =   ene_info%pme_alpha
    enefunc%pme%alpha2m         = - ene_info%pme_alpha*ene_info%pme_alpha
    enefunc%pme%alpha2sp        = 2.0_wp*ene_info%pme_alpha/sqrt(PI)
    enefunc%pme%pme_max_spacing = ene_info%pme_max_spacing

    ! Setup grid number
    !
    call setup_pme_grid_number(ene_info, boundary, enefunc, ny, nz, nx1, ny1)

    ! Setup parameters for reciprocal part
    !
    if (reciprocal_calc) &
      call setup_pme_reciprocal_part(molecule, ny, nz, nx1, ny1, enefunc)

    ! Setup parameters for real part
    !
    if (real_calc) &
      call setup_pme_real_part(enefunc)

    ! Compute Ewald self energy
    !
    call setup_self_energy(ene_info, molecule, enefunc)

    return

  end subroutine define_enefunc_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_pme
  !> @brief        Setup for PME
  !! @authors      JJ, TM, CK
  !! @param[in]    ene_info : parameters in [ENERGY] section
  !! @param[in]    boundary : boundary conditions information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme_grid_number(ene_info, boundary, enefunc,  &
                                   ny, nz, nx1, ny1)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(out)   :: ny
    integer,                 intent(out)   :: nz
    integer,                 intent(out)   :: nx1
    integer,                 intent(out)   :: ny1

    ! local variables
    real(wp)                  :: grid_space
    integer                   :: i, ngrid(3)
    integer                   :: ngridmax
    integer                   :: remainder, quotient, expo


    ! Setting the number of grids
    !
    ngrid(1) = ene_info%pme_ngrid_x
    ngrid(2) = ene_info%pme_ngrid_y
    ngrid(3) = ene_info%pme_ngrid_z

    if (ngrid(1) == 0 .and. ngrid(2) == 0 .and. ngrid(3) == 0) then
      grid_space = enefunc%pme%pme_max_spacing
      ngrid(1) = int(boundary%box_size_x/grid_space)
      ngrid(2) = int(boundary%box_size_y/grid_space)
      ngrid(3) = int(boundary%box_size_z/grid_space)
    endif

    ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    ! grid partition
    !
#ifdef FFTW
    ny = int(sqrt(real(nproc_city,wp)))
    do while(mod(nproc_city,ny) /= 0)
      ny = ny - 1
    end do
    nz = nproc_city / ny
#else
    if (nproc_city == 1) then
      ny = 1
      nz = 1
    else if (nproc_city == 2) then
      ny = 1
      nz = 2
    else if (nproc_city == 4) then
      ny = 2
      nz = 2
    else if (nproc_city == 8) then
      ny = 2
      nz = 4
    else if (nproc_city == 16) then
      ny = 4
      nz = 4
    else if (nproc_city == 32) then
      ny = 4
      nz = 8
    else if (nproc_city == 64) then
      ny = 8
      nz = 8
    else if (nproc_city == 128) then
      ny = 8
      nz = 16
    else
      call error_msg( &
           'Setup_PME_FFTE> MPI process number for PME should be'//  &
           ' 1, 2, 4, 8, 16, 32, 64, or 128')
    end if
#endif

    nx1 = ny
    ny1 = nz

#ifdef FFTW
    ! Check grid point in x direction
    !
    ngrid(1) = ngrid(1) + mod(ngrid(1),2)


    ! Check grid point in y direction
    !
    remainder = mod(ngrid(2),2*ny*nz)
    if (remainder /= 0) ngrid(2) = ngrid(2) + 2*ny*nz - remainder


    ! Check grid point in z direction
    !
    remainder = mod(ngrid(3),ny*nz)
    if (remainder /= 0) ngrid(3) = ngrid(3) + ny*nz - remainder

#else
    ! Check grid point in x direction
    !
    remainder = mod(ngrid(1),ny)
    if (remainder /= 0) ngrid(1) = ngrid(1) + ny - remainder
    quotient = ngrid(1)/ny

    if (quotient <= Index(NumIndex)) then
      do i = 1, NumIndex
        if (quotient <= Index(i)) exit
      end do
      quotient = Index(i)
      ngrid(1) = ny * quotient
    else
      expo = int(log(real(quotient,wp))/log(real(2,wp)))
      if (2**expo >= quotient) then
        quotient = 2**expo
      else
        quotient = 2**(expo+1)
      end if
      ngrid(1) = ny * quotient
    end if

    ! Check grid point in y direction
    !
    remainder = mod(ngrid(2),nz)
    if (remainder /= 0) ngrid(2) = ngrid(2) + nz - remainder

    quotient = ngrid(2)/nz

    if (quotient <= Index(NumIndex)) then
      do i = 1, NumIndex
        if (quotient <= Index(i)) exit
      end do
      quotient = Index(i)
      ngrid(2) = nz * quotient
    else
      expo = int(log(real(quotient,wp))/log(real(2,wp)))
      if (2**expo >= quotient) then
        quotient = 2**expo
      else
        quotient = 2**(expo+1)
      end if
      ngrid(2) = nz * quotient
    end if

    ! Check grid point in z direction
    !
    remainder = mod(ngrid(3),nz)
    if (remainder /= 0) ngrid(3) = ngrid(3) + nz - remainder

    quotient = ngrid(3)/nz

    if (quotient <= Index(NumIndex)) then
      do i = 1, NumIndex
        if (quotient <= Index(i)) exit
      end do
      quotient = Index(i)
      ngrid(3) = nz * quotient
    else
      expo = int(log(real(quotient,wp))/log(real(2,wp)))
      if (2**expo >= quotient) then
        quotient = 2**expo
      else
        quotient = 2**(expo+1)
      end if
      ngrid(3) = nz * quotient
    end if
#endif

    ! Show warning messages
    !
    if ((ene_info%pme_ngrid_x /= ngrid(1)) .or. &
        (ene_info%pme_ngrid_y /= ngrid(2)) .or. &
        (ene_info%pme_ngrid_z /= ngrid(3))) then
      if (main_rank) then
        write(MsgOut,'(A)') &
          ''
        if (ene_info%pme_ngrid_x == 0 .and. &
            ene_info%pme_ngrid_y == 0 .and. &
            ene_info%pme_ngrid_z == 0) then
          write(MsgOut,'(A)') &
            'Setup_PME> Proper PME grid number was generated automatically'
        else
          write(MsgOut,'(A)') &
            '  WARNING: PME grid number is different from the input'
        end if
        write(MsgOut,'(A20,3I10)') &
          '  pme_ngrid(x,y,z)= ', ngrid(1), ngrid(2), ngrid(3)
        write(MsgOut,'(A)') &
          ''
      end if
    end if

    enefunc%pme%ngrid(1) = ngrid(1)
    enefunc%pme%ngrid(2) = ngrid(2)
    enefunc%pme%ngrid(3) = ngrid(3)

    return

  end subroutine setup_pme_grid_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme_reciprocal_part
  !> @brief        Setup for PME
  !! @authors      TI, JJ, TM
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme_reciprocal_part(molecule, ny, nz, nx1, ny1, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    integer,                 intent(in)    :: ny
    integer,                 intent(in)    :: nz
    integer,                 intent(in)    :: nx1
    integer,                 intent(in)    :: ny1
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                  :: bcos, bsin, fact
    integer                   :: i, j, k, js, n_bs, ngrid(3)
    integer                   :: ngridmax
    integer                   :: nlocalx1, nlocaly1, nlocaly, nlocalz
    integer                   :: index_y, index_z
    integer                   :: istart, iend
    integer,      allocatable :: x_start1(:), x_end1(:)
    integer,      allocatable :: y_start(:),  y_end(:)
    integer,      allocatable :: y_start1(:), y_end1(:)
    integer,      allocatable :: z_start(:),  z_end(:)
    integer,      allocatable :: z_local(:),  y_local(:)
    integer,      allocatable :: x_local1(:), y_local1(:)
    integer,      allocatable :: iproc_y(:),  iproc_z(:)
    integer,      allocatable :: iproc_x1(:), iproc_y1(:)
    real(wp),     allocatable :: bs(:)


    ngrid(1) = enefunc%pme%ngrid(1)
    ngrid(2) = enefunc%pme%ngrid(2)
    ngrid(3) = enefunc%pme%ngrid(3)
    ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    nlocaly  = int((ngrid(2) + ny  - 1)/ny )
    nlocalz  = int((ngrid(3) + nz  - 1)/nz )
    nlocalx1 = int((ngrid(1) + nx1 - 1)/nx1)
    nlocaly1 = int((ngrid(2) + ny1 - 1)/ny1)

    ! allocate local array
    !
    allocate(x_start1(0:nproc_city-1), x_end1  (0:nproc_city-1), &
             y_start (0:nproc_city-1), y_end   (0:nproc_city-1), &
             y_start1(0:nproc_city-1), y_end1  (0:nproc_city-1), &
             z_start (0:nproc_city-1), z_end   (0:nproc_city-1), &
             y_local (0:nproc_city-1), z_local (0:nproc_city-1), &
             x_local1(0:nproc_city-1), y_local1(0:nproc_city-1), &
             iproc_y (0:nproc_city-1), iproc_z (0:nproc_city-1), &
             iproc_x1(0:nproc_city-1), iproc_y1(0:nproc_city-1))

    allocate(enefunc%pme%x_start1(0:nproc_city-1), &
             enefunc%pme%x_end1  (0:nproc_city-1), &
             enefunc%pme%y_start (0:nproc_city-1), &
             enefunc%pme%y_end   (0:nproc_city-1), &
             enefunc%pme%y_local (0:nproc_city-1), &
             enefunc%pme%y_start1(0:nproc_city-1), &
             enefunc%pme%y_end1  (0:nproc_city-1), &
             enefunc%pme%z_start (0:nproc_city-1), &
             enefunc%pme%z_end   (0:nproc_city-1), &
             enefunc%pme%z_local (0:nproc_city-1))

    ! index for each processor
    !
    do i = 0, nproc_city - 1
      iproc_z(i)  = mod(i,nz) + 1
      iproc_y(i)  = (i - iproc_z(i) + 1)/nz
      iproc_y(i)  = mod(iproc_y(i),ny) + 1
      iproc_x1(i) = iproc_y(i)
      iproc_y1(i) = iproc_z(i)

      y_start(i)  = (iproc_y(i)  - 1)*nlocaly  + 1
      z_start(i)  = (iproc_z(i)  - 1)*nlocalz  + 1
      x_start1(i) = (iproc_x1(i) - 1)*nlocalx1 + 1
      y_start1(i) = (iproc_y1(i) - 1)*nlocaly1 + 1
      y_end(i)    = min(ngrid(2), nlocaly *iproc_y (i))
      z_end(i)    = min(ngrid(3), nlocalz *iproc_z (i))
      x_end1(i)   = min(ngrid(1), nlocalx1*iproc_x1(i))
      y_end1(i)   = min(ngrid(2), nlocaly1*iproc_y1(i))
      y_local(i)  = y_end(i)  - y_start (i) + 1
      z_local(i)  = z_end(i)  - z_start (i) + 1
      x_local1(i) = x_end1(i) - x_start1(i) + 1
      y_local1(i) = y_end1(i) - y_start1(i) + 1

      enefunc%pme%x_start1(i) = x_start1(i)
      enefunc%pme%x_end1  (i) = x_end1  (i)
      enefunc%pme%y_start (i) = y_start (i)
      enefunc%pme%y_end   (i) = y_end   (i)
      enefunc%pme%y_local (i) = y_local (i)
      enefunc%pme%y_start1(i) = y_start1(i)
      enefunc%pme%y_end1  (i) = y_end1  (i)
      enefunc%pme%z_start (i) = z_start (i)
      enefunc%pme%z_end   (i) = z_end   (i)
      enefunc%pme%z_local (i) = z_local (i)
    end do

    ! define number of interactions for each processor
    !
    call get_loop_index(molecule%num_atoms, istart, iend)

    enefunc%pme%istart_atom = istart
    enefunc%pme%iend_atom   = iend

    ! new communicator according to the grid index
    !
    index_z = iproc_y(my_city_rank)*ny
    index_y = iproc_z(my_city_rank)*nz
#ifdef HAVE_MPI_GENESIS
    call mpi_comm_split(mpi_comm_city,index_y,my_city_rank,grid_commy,ierror)
    call mpi_comm_size (grid_commy, nprocy, ierror)
    call mpi_comm_rank (grid_commy, my_y_rank, ierror)
    call mpi_comm_split(mpi_comm_city,index_z,my_city_rank,grid_commz,ierror)
    call mpi_comm_size (grid_commz, nprocz, ierror)
    call mpi_comm_rank (grid_commz, my_z_rank, ierror)
#endif

    ! setup bs_fact
    !
    n_bs = enefunc%pme%n_bs
    j = 1
    do i = 1, n_bs - 2
      j = j*(n_bs - i) 
    end do

    enefunc%pme%bs_fact   = 1.0_wp/real(j,wp)
    enefunc%pme%bs_fact3  = enefunc%pme%bs_fact**3
    enefunc%pme%bs_fact3d = enefunc%pme%bs_fact3 * real(n_bs - 1,wp)

    ! Preparing b2=b(h)^2, h shifted
    !
    allocate(bs(n_bs))
    allocate(enefunc%pme%b2(ngridmax,3))

    call b_spline_coef(n_bs, 1.0_wp, bs)

    do i = 1, n_bs - 1
      bs(i) = bs(i)*enefunc%pme%bs_fact
    end do

    do k = 1, 3

      fact = 2.0_wp * PI/real(ngrid(k),wp)
      do j = 1, ngrid(k)/2 + 1
        js = j - 1
        bcos = 0.0_wp
        bsin = 0.0_wp
        do i = 0, n_bs - 2
          bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
          bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
        end do
        enefunc%pme%b2(j, k) = 1.0_wp/(bcos**2 + bsin**2)
      end do

      do j = ngrid(k)/2 + 2, ngrid(k)
        js = j - 1 - ngrid(k)
        bcos = 0.0_wp
        bsin = 0.0_wp
        do i = 0, n_bs - 2
          bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
          bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
        end do
        enefunc%pme%b2(j, k) = 1.0_wp/(bcos**2 + bsin**2)
      end do

    end do

    if (mod(n_bs,2) == 1) then
      do k = 1, 3
        enefunc%pme%b2(ngrid(k)/2 + 1,k) = 0.0_wp
      end do
    end if


    ! Prepareing theta = F^-1[theta](h), h shifted
    !
    allocate(enefunc%pme%gx(ngrid(1)),                                        &
             enefunc%pme%gy(ngrid(2)),                                        &
             enefunc%pme%gz(ngrid(3)),                                        &
             enefunc%pme%vir_fact(ngrid(1), ngrid(2), ngrid(3)),              &
             enefunc%pme%theta   (ngrid(1), ngrid(2), ngrid(3)),              &
             enefunc%pme%qdf     (ngrid(1), ngrid(2), ngrid(3)),              &
             enefunc%pme%f(3, molecule%num_atoms),                            &
             enefunc%pme%v(3, molecule%num_atoms),                            &
             enefunc%pme%ftqdf_localA(ngrid(1),    y_local(0),  z_local(0)),  &
#ifdef FFTW
             enefunc%pme%ftqdf_localB(ngrid(3),  ngrid(1)/2+1,  y_local(0)),  &
#else
             enefunc%pme%ftqdf_localB(x_local1(0), y_local1(0), ngrid(3)  ),  &
#endif
             enefunc%pme%qdf_recv(ngrid(1),y_local(0),z_local(0),nproc_city), &
             enefunc%pme%qdf_send(ngrid(1),y_local(0),z_local(0),nproc_city))

    ! deallocate local array
    !
    deallocate(x_start1, x_end1, y_start, y_end,     &
               y_start1, y_end1, z_start, z_end,     &
               y_local, z_local, x_local1, y_local1, &
               iproc_y, iproc_z, iproc_x1, iproc_y1)
    deallocate(bs)

    call setup_pme(enefunc)

    return

  end subroutine setup_pme_reciprocal_part

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme_real_part
  !> @brief        Setup for PME
  !! @authors      TI, JJ, TM
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme_real_part(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc


    allocate(enefunc%pme%i_list  (  size(enefunc%nonb_excl_list(:))), &
             enefunc%pme%force_wk(3,size(enefunc%nonb_excl_list(:))))

    return

  end subroutine setup_pme_real_part

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_self_energy
  !> @brief        Setup for PME
  !! @authors      TI, TM
  !! @param[in]    ene_info : parameters in [ENERGY] section
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_self_energy(ene_info, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                  :: u_self, el_fact
    integer                   :: i


    el_fact =  ELECOEF/ene_info%dielec_const

    u_self = 0.0_wp
    do i = my_country_rank + 1, molecule%num_atoms, nproc_country
      u_self = u_self + molecule%charge(i)**2
    end do

    enefunc%pme%u_self = - u_self * el_fact * enefunc%pme%alpha/sqrt(PI)

    return

  end subroutine setup_self_energy

end module at_enefunc_pme_mod
