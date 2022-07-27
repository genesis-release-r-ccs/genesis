!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_rpath_mod
!> @brief   String method using replicas
!! @authors Yasuaki Komuro (YK), Takaharu Mori (TM), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_rpath_mod

  use sp_energy_mod
  use sp_energy_str_mod
  use sp_energy_restraints_mod
  use sp_md_leapfrog_mod
  use sp_md_vverlet_mod
  use sp_md_respa_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_dynvars_mod
  use sp_ensemble_str_mod
  use sp_rpath_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_boundary_mod
  use sp_pairlist_str_mod
  use sp_pairlist_mod
  use sp_enefunc_str_mod
  use sp_remd_str_mod
  use sp_enefunc_restraints_mod
  use dihedral_libs_mod
  use timers_mod
  use sp_output_str_mod
  use sp_output_mod
  use sp_domain_str_mod
  use sp_communicate_mod
  use fitting_str_mod
  use fitting_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use fileio_crd_mod
  use fileio_mode_mod
  use fileio_control_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_rpath_info
    integer                           :: dimension    = 0
    integer                           :: rpath_period = 0
    integer                           :: nreplica     = 0
    real(wp)                          :: delta        = 0.1_wp
    real(wp)                          :: smooth       = 0.0_wp
    logical                           :: fix_terminal = .false.
    logical                           :: avoid_shrinkage = .false.
    logical                           :: use_restart  = .true.
    integer,              allocatable :: rest_function(:)
  end type s_rpath_info

  ! subroutines
  public  :: show_ctrl_rpath
  public  :: read_ctrl_rpath
  public  :: define_nreplica
  public  :: setup_rpath
  public  :: run_rpath

  ! private subroutines
  private :: setup_rpath_restart
  private :: assign_condition
  private :: evolve
  private :: reparametrize

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_rpath
  !> @brief        show RPATH section usage
  !! @authors      NT, TM, YK
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_rpath(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode

    if (show_all) then

      select case (run_mode)

      case ('rpath')

        write(MsgOut,'(A)') '[RPATH]'
        write(MsgOut,'(A)') 'nreplica          = 4'
        write(MsgOut,'(A)') 'rpath_period      = 1000'
        write(MsgOut,'(A)') 'delta             = 0.05'
        write(MsgOut,'(A)') 'smooth            = 0.1'
        write(MsgOut,'(A)') 'fix_terminal      = NO'
        write(MsgOut,'(A)') 'avoid_shrinkage   = NO'
        write(MsgOut,'(A)') 'rest_function     = 1 2'
        write(MsgOut,'(A)') 'use_restart       = YES'
        write(MsgOut,'(A)') ''

      end select

    else

      select case (run_mode)

      case ('rpath')

        write(MsgOut,'(A)') '[RPATH]'
        write(MsgOut,'(A)') 'nreplica          = 4'
        write(MsgOut,'(A)') 'rpath_period      = 1000'
        write(MsgOut,'(A)') 'delta             = 0.05'
        write(MsgOut,'(A)') 'smooth            = 0.1'
        write(MsgOut,'(A)') 'fix_terminal      = NO'
        write(MsgOut,'(A)') 'avoid_shrinkage   = NO'
        write(MsgOut,'(A)') 'rest_function     = 1 2'
        write(MsgOut,'(A)') 'use_restart       = YES'
        write(MsgOut,'(A)') ''

      end select

    end if

    return

  end subroutine show_ctrl_rpath
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_rpath
  !> @brief        read RPATH section in the control file
  !! @authors      TM, YK
  !! @param[in]    handle     : unit number of control file
  !! @param[out]   rpath_info : RPATH section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_rpath(handle, rpath_info)

    ! parameters
    character(*),            parameter     :: Section = 'Rpath'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_rpath_info),      intent(inout) :: rpath_info

    ! local variables
    integer                  :: i
    character(30)            :: cdim, partmp, vartmp, fnctmp
    character(MaxLine)       :: rest_function_char


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'nreplica',      &
                               rpath_info%nreplica)
    call read_ctrlfile_integer(handle, Section, 'rpath_period',  &
                               rpath_info%rpath_period)
    call read_ctrlfile_real   (handle, Section, 'delta',         &
                               rpath_info%delta)
    call read_ctrlfile_real   (handle, Section, 'smooth',        &
                               rpath_info%smooth)
    call read_ctrlfile_logical(handle, Section, 'fix_terminal',  &
                               rpath_info%fix_terminal)
    call read_ctrlfile_logical(handle, Section, 'avoid_shrinkage',  &
                            rpath_info%avoid_shrinkage)
    call read_ctrlfile_string (handle, Section, 'rest_function', &
                               rest_function_char)
    call read_ctrlfile_logical(handle, Section, 'use_restart',   &
                               rpath_info%use_restart)

    call end_ctrlfile_section(handle)

    ! error check
    !
    if (rpath_info%nreplica == 0) &
      call error_msg('Read_Ctrl_Rpath> nreplica should be >0 in [RPATH]')

    if (rpath_info%delta < EPS .and. rpath_info%rpath_period > 0) &
    call error_msg('Read_Ctrl_Rpath> delta should be larger than 0 if rpath_period > 0 ')

    ! get dimension from the number of numerics in rest_function_char
    !
    rpath_info%dimension = split_num(rest_function_char)

    ! read rest_function
    !
    allocate(rpath_info%rest_function(rpath_info%dimension))
    rpath_info%rest_function(1:rpath_info%dimension) = 0
    call split(rpath_info%dimension, rpath_info%dimension, &
               rest_function_char, rpath_info%rest_function)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Rpath> Rpath information'
      write(MsgOut,'(A20,I10)') &
        '  dimension       = ', rpath_info%dimension
      write(MsgOut,'(A20,I10)') &
        '  nreplica        = ', rpath_info%nreplica
      write(MsgOut,'(A20,I10)') &
        '  rpath_period    = ', rpath_info%rpath_period
      write(MsgOut,'(A20,F10.5)') &
        '  delta           = ', rpath_info%delta
      if (rpath_info%fix_terminal) then
        write(MsgOut,'(A)') '  fix_terminal    = YES'
      else
        write(MsgOut,'(A)') '  fix_terminal    = NO'
      end if
      if (rpath_info%avoid_shrinkage) then
        write(MsgOut,'(A)') '  avoid_shrinkage = YES'
      end if
      if (rpath_info%use_restart) then
        write(MsgOut,'(A)') '  use_restart     = YES'
      else
        write(MsgOut,'(A)') '  use_restart     = NO'
      end if
      write(MsgOut,'(A20,F10.3)') &
        '  smooth          = ', rpath_info%smooth
      write(MsgOut,'(A20,A)') &
        '  rest_function   = ', trim(rest_function_char)
      write(MsgOut,'(A)') ''

    end if

    return

  end subroutine read_ctrl_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_nreplica
  !> @brief        setup structure of replica
  !! @authors      TM, YK
  !! @param[in]    rpath_info : RPATH section control parameters information
  !! @param[out]   replica    : replica information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_nreplica(rpath_info, rpath)

    ! formal arguments
    type(s_rpath_info),       intent(in)    :: rpath_info
    type(s_rpath),            intent(inout) :: rpath


    rpath%nreplica = rpath_info%nreplica

    return

  end subroutine define_nreplica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath
  !> @brief        setup RPATH
  !! @authors      TM, YK
  !! @param[in]    rpath_info : RPATH section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] domain     : domain information
  !! @param[inout] restraints : restraints information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath(rpath_info,rst, boundary, dynamics, molecule, domain, &
                         restraints, ensemble, enefunc, rpath)

    ! formal arguments
    type(s_rpath_info),      intent(in)    :: rpath_info
    type(s_rst),     target, intent(in)    :: rst
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynamics),target, intent(in)    :: dynamics
    type(s_molecule),target, intent(in)    :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_restraints),      intent(inout) :: restraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_rpath),   target, intent(inout) :: rpath

    ! local variables
    integer                    :: i, j, k, m, dimno, ierr
    integer                    :: max_nreplica
    integer                    :: natmgrp
    character(2000)            :: param
    real(wp),      allocatable :: ddata(:)
    integer,       pointer     :: grouplist(:,:), numatoms(:), atomlist(:,:)
    real(wp),      pointer     :: mass(:)
    integer                    :: iatm, iatm_xyz
    integer                    :: replicaid
    integer                    :: ifunc, igroup, idm, ikind, itmp


    mass          => molecule%mass
    grouplist     => enefunc%restraint_grouplist
    numatoms      => enefunc%restraint_numatoms
    atomlist      => enefunc%restraint_atomlist

    if (dynamics%integrator == IntegratorPMTS ) &
        call error_msg('Setup_Rpath> Rpath+RESPA is not allowed')

    ! setup reparameterization period
    !
    if (rpath_info%rpath_period > 0) then
      if (mod(dynamics%nsteps, rpath_info%rpath_period) /= 0) then
        call error_msg('Setup_Rpath> mod(nsteps,rpath_period) must be zero')
      else
        rpath%equilibration_only = .false.
        rpath%rpath_period = rpath_info%rpath_period
        rpath%ncycle = dynamics%nsteps/rpath%rpath_period
      end if
      if (dynamics%rstout_period > 0) then
        if (mod(dynamics%rstout_period,rpath_info%rpath_period) /= 0) then
          call error_msg('Setup_Rpath> mod(rstout_period, rpath_period)'//&
                         ' must be zero')
        end if
      endif
    else if (rpath_info%rpath_period == 0) then
      rpath%equilibration_only = .true.
      rpath%rpath_period = dynamics%nsteps
      rpath%ncycle = 1
    else
      call error_msg('Setup_Rpath> rpath_period should be non-negative integer')
    end if

    ifunc = rpath_info%rest_function(1)
    ikind = enefunc%restraint_kind(ifunc) 
    do i = 1, rpath_info%dimension
      itmp = rpath_info%rest_function(i)
      if (enefunc%restraint_kind(itmp) /= ikind) then
          call error_msg('Setup_Rpath> [ERROR] multiple types of '//&
          'restraint functions as CV are not allowed')
      end if
      enefunc%restraint_rpath_func(itmp)=i

    end do
    if (ikind /= RestraintsFuncPOSI .and. &
        ikind /= RestraintsFuncDIST .and. &
        ikind /= RestraintsFuncDIHED) then
      call error_msg('Setup_Rpath> [ERROR] This CV is not allowed in spdyn')
    end if

    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then
      enefunc%rpath_pos_func = ifunc
      if (enefunc%fitting_method /= FittingMethodNO) then
        if (.not. enefunc%do_fitting) &
          call error_msg('Setup_Rpath> [ERROR] fitfile should be specified')
        if (enefunc%mass_weight) &
          call error_msg('Setup_Rpath> [ERROR] mass_weight is not allowed')
      end if

      do i = 1, enefunc%num_restraintfuncs
        if (enefunc%restraint_kind(i) == ikind) then
          if (enefunc%restraint_rpath_func(i) == 0) then
            call error_msg('Setup_Rpath> [ERROR] all POSI restraint should'//&
               ' be used as R-path CV')
          end if
        end if
      end do

    end if

    ! setup basic parameters
    !
    rpath%nreplica     = rpath_info%nreplica
    max_nreplica       = rpath_info%nreplica
    rpath%delta        = rpath_info%delta
    rpath%smooth       = rpath_info%smooth
    rpath%fix_terminal = rpath_info%fix_terminal
    rpath%use_restart  = rpath_info%use_restart
    rpath%avoid_shrinkage = rpath_info%avoid_shrinkage


    ! setup force constants and references
    !
    ifunc  = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then

      rpath%dimension = enefunc%restraint_numatoms(igroup) * 3

      call alloc_rpath(rpath, RpathReplicas, rpath%dimension,    &
                    rpath%nreplica, max_nreplica)

      rpath%rest_function(1) = rpath_info%rest_function(1)
      if (rpath%dimension >= 2) rpath%rest_function(2:rpath%dimension) = 0

      call alloc_rpath(rpath, RpathUmbrellas, &
                    rpath%dimension, rpath%nreplica)

      replicaid = my_country_no + 1

      if (rpath%use_restart .and. allocated(rst%rest_reference)) then

        call setup_rpath_restart(rpath_info, rst, enefunc, rpath)

      else

        do i = 1, rpath%dimension
          iatm = enefunc%restraint_atomlist(((i-1)/3)+1,igroup)
          iatm_xyz = mod(i-1, 3) + 1
          do j = 1, rpath%nreplica
            rpath%rest_constants(1:4, i, replicaid) =  &
               enefunc%restraint_const_replica(j, ifunc)
            rpath%rest_reference(1:2, i, replicaid) =  &
               enefunc%restraint_refcoord(iatm_xyz, iatm)
            rpath%rest_reference_prev(i, replicaid) =  &
               rpath%rest_reference(1, i, replicaid)
            rpath%rest_reference_init(i, replicaid) =  &
               rpath%rest_reference(1, i, replicaid)
          end do
        end do

      end if

    else

      rpath%dimension = rpath_info%dimension

      call alloc_rpath(rpath, RpathReplicas, rpath%dimension, &
                       rpath%nreplica, max_nreplica)

      rpath%rest_function(:) = rpath_info%rest_function(:)

      ! split rpath_info%parameters (char) into parameters (real or integer)
      !
      call alloc_rpath(rpath, RpathUmbrellas, &
                       rpath%dimension, rpath%nreplica)

      if (rpath%use_restart .and. allocated(rst%rest_reference)) then

        call setup_rpath_restart(rpath_info, rst, enefunc, rpath)

      else
        do i = 1, rpath%dimension
          ifunc = rpath_info%rest_function(i)
          do j = 1, rpath%nreplica
            rpath%rest_constants(1:4, i, j) =  &
               enefunc%restraint_const_replica(j, ifunc)
            rpath%rest_reference(1:2, i, j) =  &
               enefunc%restraint_ref_replica(j, ifunc)
            rpath%rest_reference_prev(i, j) =  &
               enefunc%restraint_ref_replica(j, ifunc)
            rpath%rest_reference_init(i, j) =  &
               enefunc%restraint_ref_replica(j, ifunc)
          end do
        end do
      end if

    end if

    ! assign restraint parameters for enefunc
    !
    allocate(enefunc%stats_id_atm2cv(1:molecule%num_atoms))

    enefunc%stats_id_atm2cv(:) = 0

    do i = 1, enefunc%restraint_numatoms(igroup)
      j = atomlist(i,igroup)
      !if (j > nsolute) then
      !  call error_msg('Setup_Rpath> [ERROR] too large # of restraint atoms')
      !end if
      enefunc%stats_id_atm2cv(j) = i
    end do

    call assign_condition(rpath, enefunc, domain)

    ! determine the number of atoms involved in collective variable
    !
    ifunc = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    select case(enefunc%restraint_kind(ifunc))

    case (RestraintsFuncPOSI)

      natmgrp = enefunc%restraint_numatoms(igroup)

    case (RestraintsFuncPC:RestraintsFuncPCCOM)

      natmgrp = enefunc%restraint_numatoms(igroup)

    case default

      do i = 1, enefunc%restraint_funcgrp(ifunc)
        if (enefunc%restraint_numatoms(i) > 1) &
          call error_msg('Setup_Rpath> # of atoms should be 1')
      end do
      natmgrp = enefunc%restraint_funcgrp(ifunc)

    end select

    ! setup statistical variables in dynvars
    !
    enefunc%rpath_flag        = .true.
    enefunc%rpath_sum_mf_flag = .false.
    enefunc%stats_count       = 0
    enefunc%stats_natom       = natmgrp
    enefunc%stats_dimension   = rpath%dimension

    allocate(enefunc%stats_delta(1:rpath%dimension), &
             enefunc%stats_grad(1:3, 1:natmgrp, 1:rpath%dimension), &
             enefunc%stats_force(1:rpath%dimension), &
             enefunc%stats_metric(1:rpath%dimension,1:rpath%dimension), &
             enefunc%stats_atom(1:natmgrp, 1:rpath%dimension), &
             enefunc%stats_mass(1:natmgrp, 1:rpath%dimension), &
             enefunc%stats_force_save(1:rpath%dimension),      &
             enefunc%rpath_rest_function(1:rpath%dimension))

    enefunc%stats_delta(:)         = 0.0_wp
    enefunc%stats_grad(:,:,:)      = 0.0_wp
    enefunc%stats_force(:)         = 0.0_dp
    enefunc%stats_metric(:,:)      = 0.0_dp
    enefunc%stats_atom(:,:)        = 0
    enefunc%stats_mass(:,:)        = 0.0_wp
    enefunc%stats_force_save(:)    = 0.0_dp
    enefunc%rpath_rest_function(:) = 0

    ! setup atom index and mass
    !
    ifunc  = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    enefunc%rpath_rest_function(1:rpath%dimension) = &
             rpath%rest_function(1:rpath%dimension)

    select case(enefunc%restraint_kind(ifunc))

    case (RestraintsFuncPOSI)

      do i = 1, natmgrp
        j = atomlist(i,igroup)
        enefunc%stats_atom(i,:) = j
        enefunc%stats_mass(i,:) = mass(j)
      end do

    case (RestraintsFuncPC:RestraintsFuncPCCOM)

      do i = 1, natmgrp
        j = atomlist(i,igroup)
        enefunc%stats_atom(i,:) = j
        enefunc%stats_mass(i,:) = mass(j)
      end do

    case default

      do dimno = 1, rpath%dimension
        idm = rpath_info%rest_function(dimno)
        do i = 1, natmgrp
          j = atomlist(1,grouplist(i,idm))
          enefunc%stats_atom(i,dimno) = j
          enefunc%stats_mass(i,dimno) = mass(j)
        end do
      end do

    end select

    ! write the summary of setup
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Setup_Rpath> Rpath information'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  Restraints'

      if (rpath%use_restart) then
        if ((.not. allocated(rst%rest_reference))) then
          write(MsgOut,'(A)') '  use_restart flag is ignored.'
        else
          write(MsgOut,'(A)') '  References are replaced.'
          write(MsgOut, *)
          do i = 1, rpath%dimension
            !write(MsgOut,'(a,i0,a,a)')'  function', i, ' = ', trim(restraints%function(i))
            write(MsgOut,'(a,i0,a,$)')'  constant', i, ' = '
            do j = 1, rpath%nreplica
              write(MsgOut, '(F10.4,$)') rpath%rest_constants(1, i, j)
            end do
            write(MsgOut, *)
            write(MsgOut,'(a,i0,a,$)')'  reference', i, ' = '
            do j = 1, rpath%nreplica
              write(MsgOut, '(F10.4,$)') rpath%rest_reference(1, i, j)
            end do
            !write(MsgOut,'(a,i0,a,a)')'  select_index', i, ' = ', trim(restraints%select_index(i))
            write(MsgOut,'(A)') ''
            write(MsgOut,'(A)') ''
          end do
        end if
      end if

      ! do i = 1, rpath%dimension
      !   do j = 1, rpath%nreplica
      !     write(MsgOut,'(A,I4,A,I4,A,1F8.3,A,1F8.3)')     &
      !       '    REPLICA =', j, '  DIM = ', i, &
      !       '  Const = ',rpath%rest_constants(1, i, j), &
      !       '  Ref   = ',rpath%rest_reference(1, i, j)
      !   end do
      ! end do

      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_restart
  !> @brief        setup RPATH restart
  !! @authors      KT
  !! @param[in]    rpath_info : RPATH section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_restart(rpath_info, rst, enefunc, rpath)

    ! formal arguments
    type(s_rpath_info),      intent(in)    :: rpath_info
    type(s_rst),     target, intent(in)    :: rst
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_rpath),   target, intent(inout) :: rpath

    ! local variables
    integer                                :: i, j, ifunc
    integer                                :: replicaid


    replicaid = my_country_no + 1

    if (replica_main_rank) then
      call mpi_allgather(rst%rest_reference(1, :, replicaid), &
                         rpath%dimension,                     &
                         mpi_wp_real,                         &
                         rpath%rest_reference(1, :, :),       &
                         rpath%dimension,                     &
                         mpi_wp_real,                         &
                         mpi_comm_airplane, ierror)
    end if

    call mpi_bcast(rpath%rest_reference(1, :, :),  &
                   rpath%dimension*rpath%nreplica, &
                   mpi_wp_real, 0, mpi_comm_country, ierror)

    !write(MsgOut,*) my_country_no + 1

    do i = 1, rpath%dimension

      ifunc = rpath_info%rest_function(1)
      if (enefunc%restraint_kind(ifunc) /= RestraintsFuncPOSI) then
        ifunc = rpath_info%rest_function(i)
      end if

      do j = 1, rpath%nreplica

        rpath%rest_constants(1:4, i, j) = &
                                      enefunc%restraint_const_replica(j, ifunc)
        rpath%rest_reference(2, i, j)   = rpath%rest_reference(1, i, j)

        !write(MsgOut,*) rpath%rest_reference(1, i, j)

        enefunc%restraint_ref_replica(j, ifunc) = rpath%rest_reference(1, i, j)
        rpath%rest_reference_prev(i, j) =  rpath%rest_reference(1, i, j)
        rpath%rest_reference_init(i, j) =  rpath%rest_reference(1, i, j)


      end do

    end do

    !write(MsgOut,*) " "

    return

  end subroutine setup_rpath_restart

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition
  !> @brief        control replica exchange
  !! @authors      TM, YK
  !! @param[in]    rpath       : RPATH information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition(rpath, enefunc, domain)

    ! formal arguments
    type(s_rpath),           intent(in)    :: rpath
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i,id,ix,j,k,ncell_local, ifunc
    integer                  :: replicaid
    integer                  :: iatm, iatm_xyz
    integer(int2),   pointer :: id_g2l(:,:)


    id_g2l      => domain%id_g2l
    ncell_local =  domain%num_cell_local

    replicaid = my_country_no + 1

    ifunc = rpath%rest_function(1)

    select case(enefunc%restraint_kind(ifunc))
    
    case(RestraintsFuncPOSI)
        !enefunc%restraint_const(1, 1) = rpath%rest_constants(1, i, replicaid)
      do i = 1, rpath%dimension
        iatm = enefunc%restraint_atomlist(((i-1)/3)+1, &
                             enefunc%restraint_grouplist(1,ifunc))
        iatm_xyz = mod(i-1, 3) + 1
        enefunc%restraint_refcoord(iatm_xyz, iatm) = &
           rpath%rest_reference(1, i, replicaid)
      end do

      do i = 1, ncell_local

        do ix = 1, enefunc%num_restraint(i)

          j = enefunc%restraint_atom(ix, i)

          k = enefunc%stats_id_atm2cv(j)

          enefunc%restraint_coord(1,ix,i) = &
                                         rpath%rest_reference(1,3*k-2,replicaid)
          enefunc%restraint_coord(2,ix,i) = &
                                         rpath%rest_reference(1,3*k-1,replicaid)
          enefunc%restraint_coord(3,ix,i) = &
                                         rpath%rest_reference(1,3*k  ,replicaid)

        end do

      end do

    case(RestraintsFuncPC:RestraintsFuncPCCOM)

      do i = 1, rpath%dimension
        enefunc%restraint_const(1:4, i) = rpath%rest_constants(1:4, i, replicaid)
        enefunc%restraint_ref  (1:2, i) = rpath%rest_reference(1:2, i, replicaid)
      end do

      do i = 1, rpath%dimension
        enefunc%pc_force(i)  = rpath%rest_constants(1, i, replicaid)
        enefunc%pc_target(i) = rpath%rest_reference(1, i, replicaid)
      end do

    case default

      do i = 1, rpath%dimension
        enefunc%restraint_const(1:4, i) = rpath%rest_constants(1:4, i, replicaid)
        enefunc%restraint_ref  (1:2, i) = rpath%rest_reference(1:2, i, replicaid)
      end do

    end select

    !write(*,*)'##### ', replicaid, enefunc%restraint_const(1, :)
    !write(*,*)'##### ', replicaid, enefunc%restraint_ref(1, :)

    return

  end subroutine assign_condition

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath
  !> @brief        run string method
  !! @authors      TM, YK
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : communicator for domain
  !! @param[inout] rpath       : RPATH information
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath(output, domain, enefunc, dynvars, dynamics, &
                       pairlist, boundary, constraints, ensemble,  &
                       comm, rpath, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm
    type(s_rpath),           intent(inout) :: rpath
    type(s_remd),            intent(inout) :: remd

    ! local variables
    integer                   :: i, j, k
    integer                   :: iloop_start, iloop_end
    integer                   :: replicaid


    ! Open output files
    !
    call open_output(output)
    DynvarsOut = output%logunit

    ! assign conditions to enefunc
    !
    call assign_condition(rpath, enefunc, domain)

    ! run rpath
    !
    do i = 1, rpath%ncycle
      dynamics%istart_step  = (i-1)*rpath%rpath_period + 1
      dynamics%iend_step    =  i   *rpath%rpath_period
      dynamics%initial_time = dynvars%time

      enefunc%rpath_sum_mf_flag = .false.

      ! MD main loop
      !
      if (dynamics%integrator == IntegratorLEAP) then
        call leapfrog_dynamics(output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd)
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd)
      else if (dynamics%integrator == IntegratorVRES) then
        call vverlet_respa_dynamics (output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd)
      end if

      ! calculate mean forces
      !
      call calc_mean_force(enefunc, rpath)

      ! path evolution
      !
      if (.not. rpath%equilibration_only) then
        call evolve(enefunc, rpath)
        call reparametrize(rpath)
        call assign_condition(rpath, enefunc, domain)
      end if

      ! output dynamical variables
      !
      call output_rpath(output, domain, dynamics, dynvars, boundary, rpath, &
                                                                     enefunc)

      ! initialize parameters
      !
      call init_params(enefunc)

    end do
     
    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    evolve
  !> @brief        evolve path
  !! @authors      YK, YM
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine evolve(enefunc, rpath)

    ! formal arguments
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_rpath),    target, intent(inout) :: rpath

    ! local variables
    integer                  :: dimno, dimno_i, dimno_j
    integer                  :: repid
    real(wp)                 :: force_wp

    real(dp),        pointer :: force(:)

    if (.not. replica_main_rank) return

    force  => enefunc%stats_force

    repid = my_country_no + 1

    do dimno_i = 1, rpath%dimension
     
      force_wp = real(force(dimno_i),wp)
      rpath%rest_reference(1:2, dimno_i, repid) = &
        rpath%rest_reference(1:2, dimno_i, repid) + rpath%delta * force_wp
    end do

    return

  end subroutine evolve

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_mean_force
  !> @brief        calculate mean force
  !! @authors      KT
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_mean_force(enefunc, rpath)

    ! formal arguments
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_rpath),    target, intent(inout) :: rpath

    ! local variables
    integer                  :: dimno, dimno_i, dimno_j
    integer                  :: repid, repid_i, repid_j
    real(dp)                 :: dnorm
    real(dp),    allocatable :: tangent_vector(:)
    integer                  :: ifunc

    integer,         pointer :: count
    real(dp),        pointer :: force(:), metric(:,:)
    real(dp),        pointer :: before_gather(:), after_gather(:)


    if (.not. replica_main_rank) return

    count  => enefunc%stats_count
    force  => enefunc%stats_force
    metric => enefunc%stats_metric

    before_gather => rpath%before_gather
    after_gather  => rpath%after_gather

    repid = my_country_no + 1

    do dimno = 1, rpath%dimension
      force(dimno) = force(dimno) / real(count,dp)
    end do

    enefunc%stats_force_save(:) = force(:)

    do dimno_i = 1, rpath%dimension
      do dimno_j = 1, rpath%dimension
        metric(dimno_i,dimno_j) = metric(dimno_i,dimno_j) / real(count,dp)
      end do
    end do

    force = matmul(metric, force)

    if (rpath%fix_terminal) then
      if ((repid == 1) .or. (repid == rpath%nreplica)) then
        force(:) = 0.0_dp
      end if
    else if (rpath%avoid_shrinkage) then
#ifdef HAVE_MPI_GENESIS
      do dimno = 1, rpath%dimension
        before_gather(dimno) = rpath%rest_reference(1, dimno, repid)
      end do
      call mpi_allgather(before_gather, rpath%dimension, MPI_Real8, &
                         after_gather,  rpath%dimension, MPI_Real8, &
                         mpi_comm_airplane, ierror)
      if(repid == 1) then
        repid_i = 1
        repid_j = 2
      else if(repid == rpath%nreplica) then
        repid_i = rpath%nreplica
        repid_j = rpath%nreplica - 1
      end if
      if(((repid == 1) .or. (repid == rpath%nreplica))  &
          .and. (repid_i > 0) .and. (repid_j > 0)) then
        allocate(tangent_vector(rpath%dimension))
        do dimno = 1, rpath%dimension
          tangent_vector(dimno) = after_gather((repid_j-1)*rpath%dimension+dimno) &
                                - after_gather((repid_i-1)*rpath%dimension+dimno)
        end do
        dnorm = 0.0_dp
        do dimno = 1, rpath%dimension
          dnorm = dnorm + tangent_vector(dimno)**2
        end do
        dnorm = sqrt(dnorm)
        do dimno = 1, rpath%dimension
          tangent_vector(dimno) = tangent_vector(dimno)/dnorm
        end do
        do dimno = 1, rpath%dimension
          force(dimno) = force(dimno) - force(dimno)*tangent_vector(dimno)
        end do
        deallocate(tangent_vector)
      end if
#endif
    end if

    return

  end subroutine calc_mean_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reparametrize
  !> @brief        reparametrize path
  !! @authors      YK, YM
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reparametrize(rpath)

    ! formal arguments
    type(s_rpath),    target, intent(inout) :: rpath

    ! local variables
    integer                  :: i, j, k, nrep_int
    integer                  :: dimno, dimno_i, dimno_j
    integer                  :: repid, repid_i, repid_j

    real(wip)                :: smooth_wip
    real(wip)                :: distance_prev, distance_init, dtmp
    real(wip),   allocatable :: path(:,:), path_smooth(:,:), path_reparm(:,:)
    real(wip),   allocatable :: path_leng(:), path_equi(:)
    real(dp),       pointer  :: before_gather(:), after_gather(:)


    before_gather => rpath%before_gather
    after_gather  => rpath%after_gather

    smooth_wip = real(rpath%smooth,wip)

    repid = my_country_no + 1

    allocate(path(rpath%nreplica, rpath%dimension),&
             path_smooth(rpath%nreplica, rpath%dimension),&
             path_reparm(rpath%nreplica, rpath%dimension),&
             path_leng(rpath%nreplica), path_equi(rpath%nreplica))

    do dimno = 1, rpath%dimension
      before_gather(dimno) = real(rpath%rest_reference(1, dimno, repid),dp)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_gather(before_gather, rpath%dimension, MPI_Real8,&
                    after_gather,  rpath%dimension, MPI_Real8,&
                    0, mpi_comm_airplane, ierror)

    if (main_rank) then
      do repid_i = 1, rpath%nreplica
        do dimno = 1, rpath%dimension
          path(repid_i,dimno) = real(after_gather((repid_i-1)*rpath%dimension+dimno),wip)
        end do
      end do

      ! smooth
      path_smooth = path
      do repid_i = 2, rpath%nreplica - 1
        do dimno = 1, rpath%dimension
          path_smooth(repid_i,dimno) = &
            (1.0_wip-smooth_wip)*path(repid_i,dimno) + &
            (smooth_wip*0.5_wip)*(path(repid_i-1,dimno) + path(repid_i+1,dimno))
        end do
      end do

      ! calc path length
      path_leng = 0.0_wip
      do repid_i = 2, rpath%nreplica
        path_leng(repid_i) = path_leng(repid_i-1) + &
          sqrt(sum((path_smooth(repid_i,:)-path_smooth(repid_i-1,:))**2))
        ! write(6,*) 'BEFORE inter_replica', &
        !             sqrt(sum((path_smooth(repid_i,:)-path_smooth(repid_i-1,:))**2))
      end do

      ! equi path length
      path_equi = 0.0_wip
      do repid_i = 2, rpath%nreplica
        nrep_int=rpath%nreplica -1
        path_equi(repid_i) = (repid_i-1)* &
                path_leng(rpath%nreplica)/real(nrep_int,wip)
      end do

      ! reparametrize
      path_reparm = path
      do repid_i = 2, (rpath%nreplica-1)
        do repid_j = 2, rpath%nreplica

          if ((path_leng(repid_j-1) <  path_equi(repid_i)) .and. &
              (path_equi(repid_i))  <= path_leng(repid_j)) then

            path_reparm(repid_i,:) = path_smooth(repid_j-1,:) + &
              (path_equi(repid_i) - path_leng(repid_j-1)) * &
              (path_smooth(repid_j,:) - path_smooth(repid_j-1,:)) / &
              sqrt(sum((path_smooth(repid_j,:) - path_smooth(repid_j-1,:))**2))
            exit
          end if

        end do
      end do
      distance_prev = 0.0_wip
      distance_init = 0.0_wip
      do repid_i = 1, rpath%nreplica
        do dimno = 1, rpath%dimension
          dtmp = path_reparm(repid_i,dimno)
          distance_prev = distance_prev +   &
            (real(rpath%rest_reference_prev(dimno, repid_i),wip) - dtmp)**2
          distance_init = distance_init +   &
            (real(rpath%rest_reference_init(dimno, repid_i),wip) - dtmp)**2

          rpath%rest_reference_prev(dimno, repid_i) = real(dtmp,wp)
        end do
      end do
      distance_prev = sqrt(distance_prev)
      distance_init = sqrt(distance_init)
      rpath%sum_distance = real(path_leng(rpath%nreplica),wp)
      rpath%distance_prev    = real(distance_prev,wp)/ &
                               real((rpath%nreplica), wp)
      rpath%distance_init    = real(distance_init,wp)/ &
                             real((rpath%nreplica), wp)

    end if

    ! broadcast restraint reference
    call mpi_bcast(path_reparm, rpath%dimension*rpath%nreplica,&
                   mpi_wip_real, 0, mpi_comm_world, ierror)
#endif

    do dimno = 1, rpath%dimension
      rpath%rest_reference(1:2, dimno, repid) =  &
                             real(path_reparm(repid, dimno),wp)
    end do

    deallocate(path, path_smooth, path_reparm, path_leng, path_equi)

    return

  end subroutine reparametrize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_params
  !> @brief        initialize parameters
  !! @authors      KT
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_params(enefunc)

    ! formal arguments
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables

    if (.not. replica_main_rank) return

    enefunc%stats_count       = 0
    enefunc%stats_force(:)    = 0.0_dp
    enefunc%stats_metric(:,:) = 0.0_dp

    return

  end subroutine init_params

end module sp_rpath_mod
