!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_minimize_mod
!> @brief   perform energy minimization
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK), Norio Takase (NT)
!!          Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_minimize_mod

  use at_qmmm_mod
  use at_output_mod
  use at_dynvars_mod
  use at_boundary_mod
  use at_pairlist_mod
  use at_energy_mod
  use at_output_str_mod
  use at_minimize_str_mod
  use at_dynvars_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use timers_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use structure_check_mod

  implicit none
  private

  ! structures
  type, public :: s_min_info
    integer            :: method              = MinimizeMethodLBFGS
    integer            :: nsteps              = 100
    integer            :: eneout_period       =  10
    integer            :: crdout_period       =   0
    integer            :: rstout_period       =   0
    integer            :: nbupdate_period     =  10
    character(MaxLine) :: fixatm_select_index = ''
    real(wp)           :: tol_rmsg            = -1.0_wp
    real(wp)           :: tol_maxg            = -1.0_wp
    logical            :: verbose             = .false.

    ! For SD
    real(wp)           :: force_scale_init    = 0.01_wp
    real(wp)           :: force_scale_max     = 0.1_wp
!    real(wp)           :: force_scale_init = 0.00005_wp
!    real(wp)           :: force_scale_max  = 0.0001_wp

    ! For LBFGS
    integer            :: ncorrection         = 10
    logical            :: lbfgs_bnd           = .true.
    logical            :: lbfgs_bnd_qmonly    = .false.
    real(wp)           :: lbfgs_bnd_maxmove   = 0.1_wp

    ! For LBFGS - micro-iteration
    logical            :: macro               = .false.
    integer            :: start_micro         = 2
    integer            :: nsteps_micro        = 100
    real(wp)           :: tol_rmsg_micro      = -1.0_wp
    real(wp)           :: tol_maxg_micro      = -1.0_wp
    character(256)     :: macro_select_index  = ''

    ! For structure check
    logical            :: check_structure      = .true.
    logical            :: fix_ring_error       = .false.
    logical            :: fix_chirality_error  = .false.
    character(MaxLine) :: exclude_ring_grpid   = ''
    character(MaxLine) :: exclude_chiral_grpid = ''

  end type s_min_info

  integer, parameter     :: metric = 10

  ! subroutines
  public  :: show_ctrl_minimize
  public  :: read_ctrl_minimize
  public  :: setup_minimize
  private :: setup_fixatoms
  private :: setup_microatoms
  private :: print_atoms_info
  public  :: run_min
  public  :: steepest_descent
  public  :: minimize_lbfgs
  public  :: macro_iter
  public  :: micro_iter
  public  :: add_fixatm
  private :: delete_element
  private :: calc_rmsg

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_minimize
  !> @brief        show usage of MINIMIZE section
  !! @authors      NT, KY
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_minimize(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('min')

        write(MsgOut,'(A)') '[MINIMIZE]'
        write(MsgOut,'(A)') 'method            = LBFGS    # [SD] or [LBFGS]'
        write(MsgOut,'(A)') 'nsteps            = 100      # number of minimization steps'
        write(MsgOut,'(A)') '# eneout_period   = 10       # energy output period'
        write(MsgOut,'(A)') '# crdout_period   = 0        # coordinates output period'
        write(MsgOut,'(A)') '# rstout_period   = 0        # restart output period'
        write(MsgOut,'(A)') '# nbupdate_period = 10       # nonbond update period'
        write(MsgOut,'(A)') '# fixatm_select_index = 1 2  # fixed atoms in minimization'
        write(MsgOut,'(A)') '# tol_rmsg        = 3.6D-01  # tolerance for RMS gradient'
        write(MsgOut,'(A)') '# tol_maxg        = 5.4D-01  # tolerance for MAX gradient'
        write(MsgOut,'(A)') '# verbose         = NO       # output verbosly'
        write(MsgOut,'(A)') '# for SD'
        write(MsgOut,'(A)') '# force_scale_init = 0.01D+00  # initial scaling for SD'
        write(MsgOut,'(A)') '# force_scale_max  = 0.1D+00   # max scaling for SD'
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# for LBFGS'
        write(MsgOut,'(A)') '# ncorrection       = 10      # parameter in LBFGS'
        write(MsgOut,'(A)') '# lbfgs_bnd         = yes     # add boundary in LBFGS'
        write(MsgOut,'(A)') '# lbfgs_bnd_qmonly  = yes     # add boundary only to QM atoms in LBFGS'
        write(MsgOut,'(A)') '# lbfgs_bnd_maxmove = 0.1D+00 # max move of atoms in [LBFGS]'
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# for micro-iteration'
        write(MsgOut,'(A)') '# macro        = NO    # do macro/micro iteration for QM/MM'
        write(MsgOut,'(A)') '# start_micro  = 2     # step number to start micro iteration'
        write(MsgOut,'(A)') '# nsteps_micro = 100   # number of minimization steps in micro iteration'
        write(MsgOut,'(A)') '# tol_rmsg_micro = 2.7D-01  # tolerance for RMS gradient in micro iteration'
        write(MsgOut,'(A)') '# tol_maxg_micro = 4.1D-01  # tolerance for RMS gradient in micro iteration'
        write(MsgOut,'(A)') '# macro_select_index = 1 2  # atoms optimized in macro-iteration'
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# for structure check'
        write(MsgOut,'(A)') '# check_structure      = YES  # check structure in protein and DNA/RNA'
        write(MsgOut,'(A)') '# fix_ring_error       = NO   # fix ring penetration'
        write(MsgOut,'(A)') '# fix_chirality_error  = NO   # fix chirality error'
        write(MsgOut,'(A)') '# exclude_ring_grpid   =      # exclusion list for ring error fix'
        write(MsgOut,'(A)') '# exclude_chiral_grpid =      # exclusion list for chirality error fix'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('min')

        write(MsgOut,'(A)') '[MINIMIZE]'
        write(MsgOut,'(A)') 'method        = LBFGS     # [SD] of [LBFGS]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of minimization steps'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_minimize
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_minimize
  !> @brief        read MINIMIZE section in the control file
  !! @authors      TM, KY
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   min_info : MINIMIZE section in control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_minimize(handle, min_info)

    ! parameters
    character(*),            parameter     :: Section = 'Minimize'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_min_info),        intent(inout) :: min_info

    ! local
    logical :: found_error


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'method',                   &
                               min_info%method, MinimizeMethodTypes)
    call read_ctrlfile_integer(handle, Section, 'nsteps',                   &
                               min_info%nsteps)
    call read_ctrlfile_integer(handle, Section, 'eneout_period',            &
                               min_info%eneout_period)
    call read_ctrlfile_integer(handle, Section, 'crdout_period',            &
                               min_info%crdout_period)
    call read_ctrlfile_integer(handle, Section, 'rstout_period',            &
                               min_info%rstout_period)
    call read_ctrlfile_integer(handle, Section, 'nbupdate_period',          &
                               min_info%nbupdate_period)
    call read_ctrlfile_string (handle, Section, 'fixatm_select_index',      &
                               min_info%fixatm_select_index)
    call read_ctrlfile_real   (handle, Section, 'tol_rmsg',                 &
                               min_info%tol_rmsg)
    call read_ctrlfile_real   (handle, Section, 'tol_maxg',                 &
                               min_info%tol_maxg)
    call read_ctrlfile_logical(handle, Section, 'verbose',                  &
                               min_info%verbose)

    ! SD
    call read_ctrlfile_real   (handle, Section, 'force_scale_max',          &
                               min_info%force_scale_max)
    call read_ctrlfile_real   (handle, Section, 'force_scale_init',         &
                               min_info%force_scale_init)

    ! LBFGS
    call read_ctrlfile_integer(handle, Section, 'ncorrection',              &
                               min_info%ncorrection)
    call read_ctrlfile_logical(handle, Section, 'lbfgs_bnd',                &
                               min_info%lbfgs_bnd)
    call read_ctrlfile_logical(handle, Section, 'lbfgs_bnd_qmonly',         &
                               min_info%lbfgs_bnd_qmonly)
    call read_ctrlfile_real   (handle, Section, 'lbfgs_bnd_maxmove',        &
                               min_info%lbfgs_bnd_maxmove)

    ! LBFGS - micro-iteration
    call read_ctrlfile_logical(handle, Section, 'macro',                    &
                               min_info%macro)
    call read_ctrlfile_integer(handle, Section, 'start_micro',              &
                               min_info%start_micro)
    call read_ctrlfile_integer(handle, Section, 'nsteps_micro',             &
                               min_info%nsteps_micro)
    call read_ctrlfile_real   (handle, Section, 'tol_rmsg_micro',           &
                               min_info%tol_rmsg_micro)
    call read_ctrlfile_real   (handle, Section, 'tol_maxg_micro',           &
                               min_info%tol_maxg_micro)
    call read_ctrlfile_string (handle, Section, 'macro_select_index',       &
                               min_info%macro_select_index)

    ! For structure check
    call read_ctrlfile_logical(handle, Section, 'check_structure',          &
                               min_info%check_structure)
    call read_ctrlfile_logical(handle, Section, 'fix_ring_error',           &
                               min_info%fix_ring_error)
    call read_ctrlfile_logical(handle, Section, 'fix_chirality_error',      &
                               min_info%fix_chirality_error)
    call read_ctrlfile_string (handle, Section, 'exclude_ring_grpid',       &
                               min_info%exclude_ring_grpid)
    call read_ctrlfile_string (handle, Section, 'exclude_chiral_grpid',     &
                               min_info%exclude_chiral_grpid)

    call end_ctrlfile_section(handle)


    ! set default values for tol_rmsg/maxg
    !
    if (min_info%tol_rmsg < 0.0 .and. min_info%tol_maxg < 0.0) then
      min_info%tol_rmsg = 0.36_wp
      min_info%tol_maxg = min_info%tol_rmsg * 1.5_wp

    else if (min_info%tol_maxg < 0.0) then
      min_info%tol_maxg = min_info%tol_rmsg * 1.5_wp

    else if (min_info%tol_rmsg < 0.0) then
      min_info%tol_rmsg = min_info%tol_maxg * 2.0_wp / 3.0_wp

    end if

    if (min_info%tol_rmsg_micro < 0.0 .and. min_info%tol_maxg_micro < 0.0) then
      min_info%tol_rmsg_micro = min_info%tol_rmsg * 0.75_wp
      min_info%tol_maxg_micro = min_info%tol_maxg * 0.75_wp

    else if (min_info%tol_maxg_micro < 0.0) then
      min_info%tol_maxg_micro = min_info%tol_rmsg_micro * 1.5_wp

    else if (min_info%tol_rmsg_micro < 0.0) then
      min_info%tol_rmsg_micro = min_info%tol_maxg_micro * 2.0_wp / 3.0_wp

    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Parameters of MIN'

      if (min_info%macro .and. min_info%start_micro < 2) then
        write(MsgOut,'(A,I5)')  &
        ' Warning: start_micro =', min_info%start_micro
        write(MsgOut,'(A)')  &
        ' Warning: start_micro must be larger than or equal to 2. It is re-set to 2.'  
        min_info%start_micro = 2
      end if

      write(MsgOut,'(A,A10,A,I10)')       &
            '  method                     = ', MinimizeMethodTypes(min_info%method), &
            '  nsteps                     = ', min_info%nsteps
      write(MsgOut,'(A,I10,A,I10)')       &
            '  eneout_period              = ', min_info%eneout_period, &
            '  crdout_period              = ', min_info%crdout_period
      write(MsgOut,'(A,I10,A,I10)')       &
            '  rstout_period              = ', min_info%rstout_period, &
            '  nbupdate_period            = ', min_info%nbupdate_period
      if(len(trim(min_info%fixatm_select_index)) /= 0) then
        write(MsgOut,'(A,A)')             &
            '  fixatm_select_index        = ', trim(min_info%fixatm_select_index)
      else
        write(MsgOut,'(A)')               &
            '  fixatm_select_index        =       none'
      end if
      write(MsgOut,'(A,E10.2,A,E10.2)')   &
            '  tol_rmsg                   = ', min_info%tol_rmsg, &
            '  tol_maxg                   = ', min_info%tol_maxg
      if (min_info%verbose) then
        write(MsgOut,'(A)')               &
            '  verbose                    =        yes'
      else
        write(MsgOut,'(A)')               &
            '  verbose                    =         no'
      end if

      if (min_info%method == MinimizeMethodSD) then
        write(MsgOut,'(A,F10.3,A,F10.3)') &
            '  force_scale_init           = ', min_info%force_scale_init, &
            '  force_scale_max            = ', min_info%force_scale_max
      end if

      if (min_info%method == MinimizeMethodLBFGS) then
        write(MsgOut,'(A,I10)')           &
            '  ncorrection                = ', min_info%ncorrection
        if (min_info%lbfgs_bnd) then
          if (min_info%lbfgs_bnd_qmonly) then
            write(MsgOut,'(A,A)')         &
            '  lbfgs_bnd                  =        yes', &
            '  lbfgs_bnd_qmonly           =        yes'
          else
            write(MsgOut,'(A,A)')         &
            '  lbfgs_bnd                  =        yes', &
            '  lbfgs_bnd_qmonly           =         no'
          end if
          if (min_info%lbfgs_bnd_maxmove < 0.0D+00) then 
            write(MsgOut,'(A,F10.5,A)')  &
            ' Warning: lbfgs_bnd_maxmove =', min_info%lbfgs_bnd_maxmove, &
            ' is too small. Restoring default value'  
            min_info%lbfgs_bnd_maxmove = 0.1D+00
          end if
          write(MsgOut,'(A,F10.5)')       &
            '  lbfgs_bnd_maxmove          = ', min_info%lbfgs_bnd_maxmove
        else
          write(MsgOut,'(A)')             &
            '  lbfgs_bnd                  =         no'
        end if

      end if

      if (min_info%macro) then
        write(MsgOut, '(A,A,I10)')        &
            '  macro                      =        yes', &
            '  start_micro                = ', min_info%start_micro
        write(MsgOut, '(A,I10)')        &
            '  nsteps_micro               = ', min_info%nsteps_micro
        write(MsgOut,'(A,E10.2,A,E10.2)')         &
            '  tol_rmsg_micro             = ', min_info%tol_rmsg_micro, &
            '  tol_maxg_micro             = ', min_info%tol_maxg_micro
        if (len(trim(min_info%macro_select_index)) /= 0) then
          write(MsgOut,'(A,A10)')         &
            '  macro_select_index         = ', trim(min_info%macro_select_index)
        else
          write(MsgOut,'(A)')             &
            '  macro_select_index         =  qmatm'
        end if

      else
        write(MsgOut, '(A)')              &
            '  macro                      =         no'

      end if

      if (min_info%check_structure) then
        write(MsgOut,'(A,$)') '  check_structure            =        yes'
      else
        write(MsgOut,'(A,$)') '  check_structure            =         no'
      end if
      if (min_info%fix_ring_error) then
        write(MsgOut,'(A)')   '  fix_ring_error             =        yes'
      else
        write(MsgOut,'(A)')   '  fix_ring_error             =         no'
      end if
      if (min_info%fix_chirality_error) then
        write(MsgOut,'(A)')   '  fix_chirality_error        =        yes'
      else
        write(MsgOut,'(A)')   '  fix_chirality_error        =         no'
      end if


      write(MsgOut,'(A)') ' '
    end if


    ! error check
    !
    if (main_rank) then
      found_error = .false.

      if (min_info%eneout_period > 0) then
        if (mod(min_info%nsteps, min_info%eneout_period) /= 0) then
          write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in eneout_period'
          write(MsgOut,'(A)') '  mod(nsteps, eneout_period) is not ZERO'
          found_error = .true.
        end if
      end if

      if (min_info%crdout_period > 0) then
        if (mod(min_info%nsteps, min_info%crdout_period) /= 0) then
          write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in crdout_period'
          write(MsgOut,'(A)') '  mod(nsteps, crdout_period) is not ZERO'
          found_error = .true.
        end if
      end if

      if (min_info%rstout_period > 0) then
        if (mod(min_info%nsteps, min_info%rstout_period) /= 0) then
          write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in rstout_period'
          write(MsgOut,'(A)') '  mod(nsteps, rstout_period) is not ZERO'
          found_error = .true.
        end if
      end if

      if (min_info%nbupdate_period > 0) then
        if (mod(min_info%nsteps, min_info%nbupdate_period) /= 0) then
          write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in nbupdate_period'
          write(MsgOut,'(A)') '  mod(nsteps, nbupdate_period) is not ZERO'
          found_error = .true.
        end if
      end if

      if (found_error) then
        call error_msg
      end if

    end if


    return

  end subroutine read_ctrl_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_minimize
  !> @brief        setup minimize information
  !> @authors      TM, KY
  !! @param[in]    min_info : MINIMIZE section in control parameters
  !! @param[in]    out_info : OUTPUT   section in control parameters
  !! @param[in]    sel_info : SELECTOR section in control parameters
  !! @param[in]    molecule : molecule information
  !! @param[in]    qmmm     : qmmm information in enefunc
  !! @param[in]    boundary : information of boundary condition
  !! @param[out]   minimize : minimize information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_minimize(min_info, out_info, sel_info, molecule, qmmm, &
                            boundary, minimize)

    ! formal arguments
    type(s_min_info),        intent(in)    :: min_info
    type(s_out_info),        intent(in)    :: out_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_boundary),        intent(in)    :: boundary
    type(s_minimize),        intent(inout) :: minimize

    ! local
    integer :: nwork, n3, i, ngroup


    ! initialize variables in minimize
    !
    call init_minimize(minimize)


    ! setup minimize
    !
    minimize%method           = min_info%method
    minimize%nsteps           = min_info%nsteps
    minimize%eneout_period    = min_info%eneout_period
    minimize%crdout_period    = min_info%crdout_period
    minimize%rstout_period    = min_info%rstout_period
    minimize%nbupdate_period  = min_info%nbupdate_period
    minimize%tol_rmsg         = min_info%tol_rmsg
    minimize%tol_maxg         = min_info%tol_maxg
    minimize%verbose          = min_info%verbose
    minimize%check_structure  = min_info%check_structure

    if (out_info%dcdfile /= '' .and. minimize%crdout_period <= 0) &
      minimize%crdout_period = minimize%nsteps
    if (out_info%rstfile /= '' .and. minimize%rstout_period <= 0) &
      minimize%rstout_period = minimize%nsteps

    minimize%eneout           = .false.
    minimize%crdout           = .false.
    minimize%rstout           = .false.
    minimize%eneout_short     = .false.

    minimize%macro            = min_info%macro
    if (.not. qmmm%do_qmmm) then
      if (min_info%macro) then
        if (main_rank) then
          write(MsgOut,'(A)')  &
            "Setup_Minimize> Warning: No QM/MM input is given. Turning off macro."
          write(MsgOut,*)
        end if
        minimize%macro = .false.
      end if
      if (min_info%lbfgs_bnd_qmonly) &
        minimize%lbfgs_bnd_qmonly = .false.
    end if

    ! SD
    if (minimize%method == MinimizeMethodSD) then
      minimize%force_scale_init = min_info%force_scale_init
      minimize%force_scale_max  = min_info%force_scale_max

      if (qmmm%do_qmmm .and. minimize%macro) then
        call error_msg('Setup_Minimize> Error: macro/micro-iteration is not &
                       &available for SD. Use method=LBFGS instead.')
      end if

    end if

    ! LBFGS
    if (minimize%method == MinimizeMethodLBFGS) then

      minimize%ncorrection  = min_info%ncorrection
      if (minimize%ncorrection < 3)  then
        minimize%ncorrection = 3
        if (main_rank) write(MsgOut, '(8x,a,i3/)') &
         "LBFGS-B> # of corrections changed to: ", minimize%ncorrection
      else if (minimize%ncorrection > 20) then
        minimize%ncorrection = 20
        if (main_rank) write(MsgOut, '(8x,a,i3/)') &
         "LBFGS-B> # of corrections changed to: ", minimize%ncorrection
      end if

      ! Is this correct? - KY
      n3    = molecule%num_atoms*3
      nwork = (2*n3 + 11*minimize%ncorrection + 8)*minimize%ncorrection + 5*n3
      call alloc_minimize(minimize, MinimizeLBFGS, n3, nwork)
      ! - KY

      minimize%lbfgs_bnd         = min_info%lbfgs_bnd
      minimize%lbfgs_bnd_qmonly  = min_info%lbfgs_bnd_qmonly
      if(.not. minimize%lbfgs_bnd) minimize%lbfgs_bnd_qmonly  = .false.
          
      minimize%lbfgs_bnd_maxmove = min_info%lbfgs_bnd_maxmove

      ! LBFGS - micro-iteration
      minimize%start_micro      = min_info%start_micro
      minimize%nsteps_micro     = min_info%nsteps_micro
      minimize%tol_rmsg_micro   = min_info%tol_rmsg_micro
      minimize%tol_maxg_micro   = min_info%tol_maxg_micro

      if (minimize%macro) qmmm%qm_get_esp = .true.

    end if

    ngroup = split_num(trim(min_info%fixatm_select_index))
    if (boundary%num_fixatm == 0 .and. ngroup == 0) then
      ! No fixed atoms
      !
      minimize%num_fixatoms = 0
      minimize%num_optatoms = molecule%num_atoms
      allocate(minimize%optatom_id(minimize%num_optatoms))
      do i = 1, molecule%num_atoms
        minimize%optatom_id(i) = i
      end do

    else
      ! Setup fixed atoms
      !
      call setup_fixatoms(min_info, sel_info, molecule, &
                          boundary, minimize)
    end if

    ! Setup macro/micro atoms
    !
    if (minimize%macro) &
      call setup_microatoms(min_info, sel_info, molecule, qmmm, minimize)

    ! print_info
    call print_atoms_info(molecule, minimize)

    return

  end subroutine setup_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fixatoms
  !> @brief        define fixed and opt atoms
  !! @authors      YA, KY
  !! @param[in]    min_info    : minimize input information
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    molecule    : molecular information
  !! @param[in]    boundary    : information of boundary condition
  !! @param[out]   minimize    : minimize parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fixatoms(min_info, sel_info, molecule, boundary, minimize)

    ! formal arguments
    type(s_min_info),        intent(in)    :: min_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_minimize),        intent(inout) :: minimize

    ! local variables

    integer                  :: i, j, nfix, nopt
    integer                  :: ngroup, igroup

    logical,          allocatable :: fixatm(:)
    integer,          allocatable :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)


    allocate(fixatm(molecule%num_atoms))

    ! Fixed atoms at the bounary
    !
    if (boundary%num_fixatm > 0) then
      fixatm = boundary%fixatm
    else
      fixatm = .false.
    end if

    ! Fixed atoms from the input
    !
    ngroup = split_num(trim(min_info%fixatm_select_index))
    if (ngroup /= 0) then
      allocate(group_list(ngroup), selatoms(ngroup))
      call split(ngroup, ngroup, min_info%fixatm_select_index, group_list)

      do i = 1, ngroup
        igroup = group_list(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        do j = 1, size(selatoms(i)%idx)
          fixatm(selatoms(i)%idx) = .true.
        end do
      end do
      deallocate(group_list, selatoms)

    end if

    ! Set up fixed/opt atoms
    minimize%num_fixatoms = 0
    do i = 1, molecule%num_atoms
      if (fixatm(i)) &
        minimize%num_fixatoms = minimize%num_fixatoms + 1
    end do

    minimize%num_optatoms = molecule%num_atoms - minimize%num_fixatoms
    allocate(minimize%fixedatom_id(minimize%num_fixatoms))
    allocate(minimize%optatom_id  (minimize%num_optatoms))

    nfix = 0
    nopt = 0
    do i = 1, molecule%num_atoms
      if (fixatm(i)) then
        nfix = nfix + 1
        minimize%fixedatom_id(nfix) = i
      else
        nopt = nopt + 1
        minimize%optatom_id(nopt) = i
      end if
    end do

    deallocate(fixatm)

    return

  end subroutine setup_fixatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_microatoms
  !> @brief        define macro/micro atoms in opt atoms
  !! @authors      YA, KY
  !! @param[in]    min_info    : minimize input information
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    molecule    : molecular information
  !! @param[in]    qmmm        : QM/MM information
  !! @param[out]   minimize    : minimize parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_microatoms(min_info, sel_info, molecule, qmmm, minimize)

    ! formal arguments
    type(s_min_info),        intent(in)    :: min_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(in)    :: qmmm
    type(s_minimize),        intent(inout) :: minimize

    ! local variables
    integer                       :: i, j
    integer                       :: ntmp, num_macro, num_micro
    integer                       :: igroup, ngroup
    integer                       :: iatom, natom, offset, temp

    integer,          allocatable :: itmp_macro(:), itmp_micro(:)
    integer,          allocatable :: itmp_qmmmbond_list(:)
    integer,          allocatable :: group_list(:)
    type(s_selatoms), allocatable :: selatoms(:)
    integer,          allocatable :: input_id(:)


    ! Set up opt atoms for macro/micro-iterations
    ! 
    ntmp = minimize%num_optatoms
    allocate(itmp_macro(ntmp), itmp_micro(ntmp))
    itmp_micro = minimize%optatom_id
    num_micro  = minimize%num_optatoms
    itmp_macro = 0
    num_macro  = 0

    ! delete QM atoms
    call delete_element(qmmm%qm_natoms, qmmm%qmatom_id, &
                        num_micro, itmp_micro,          &
                        ntmp, itmp_macro(num_macro+1:))
    num_macro = num_macro + ntmp

    ! delete link MM atoms
    if (qmmm%num_qmmmbonds > 0) then
      allocate(itmp_qmmmbond_list(qmmm%num_qmmmbonds))

      itmp_qmmmbond_list = qmmm%qmmmbond_list(2,:)

      ! Sort atom indices in ascending order
      !
      do i = qmmm%num_qmmmbonds, 2, -1
        do j = 1, i - 1
          if (itmp_qmmmbond_list(j) > itmp_qmmmbond_list(j+1)) then
            temp = itmp_qmmmbond_list(j)
            itmp_qmmmbond_list(j)   = itmp_qmmmbond_list(j+1)
            itmp_qmmmbond_list(j+1) = temp
          end if
        end do
      end do

      call delete_element(qmmm%num_qmmmbonds, itmp_qmmmbond_list, &
                          num_micro, itmp_micro,          &
                          ntmp, itmp_macro(num_macro+1:))
      num_macro = num_macro + ntmp

      deallocate(itmp_qmmmbond_list)
    end if

    ! Get input atom ID
    ! 
    ngroup = split_num(trim(min_info%macro_select_index))
    if (ngroup /= 0) then
      allocate(group_list(ngroup), selatoms(ngroup))
      call split(ngroup, ngroup, min_info%macro_select_index, group_list)

      natom = 0
      do i = 1, ngroup
        igroup = group_list(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        natom = natom + size(selatoms(i)%idx)
      end do

      allocate(input_id(natom))

      offset = 0
      do i = 1, ngroup
        iatom = size(selatoms(i)%idx)
        input_id(offset+1:offset+iatom) = selatoms(i)%idx(1:iatom)
        offset = offset + iatom
      end do

      ! Sort input atom indices in ascending order
      !
      do i = natom, 2, -1
        do j = 1, i - 1
          if (input_id(j) > input_id(j+1)) then
            temp = input_id(j)
            input_id(j)   = input_id(j+1)
            input_id(j+1) = temp
          end if
        end do
      end do

      ! delete input atoms
      call delete_element(natom, input_id, num_micro, itmp_micro,   &
                          ntmp, itmp_macro(num_macro+1:))
      num_macro = num_macro + ntmp

      deallocate(group_list, selatoms, input_id)

    end if

    minimize%num_optatoms_micro = num_micro
    allocate(minimize%optatom_micro_id(minimize%num_optatoms_micro))
    do i = 1, minimize%num_optatoms_micro
      minimize%optatom_micro_id(i) = itmp_micro(i)
    end do

    minimize%num_optatoms_macro = num_macro
    allocate(minimize%optatom_macro_id(minimize%num_optatoms_macro))
    do i = 1, minimize%num_optatoms_macro
      minimize%optatom_macro_id(i) = itmp_macro(i)
    end do

    ! Sort atom indices in ascending order
    !
    do i = num_macro, 2, -1
      do j = 1, i - 1
        if (minimize%optatom_macro_id(j) > minimize%optatom_macro_id(j+1)) then
          temp = minimize%optatom_macro_id(j)
          minimize%optatom_macro_id(j)   = minimize%optatom_macro_id(j+1)
          minimize%optatom_macro_id(j+1) = temp
        end if
      end do
    end do

    deallocate(itmp_micro, itmp_macro)

    return

  end subroutine setup_microatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_atoms_info
  !> @brief        print information of atoms
  !! @authors      KY
  !! @param[in]    molecule    : molecular information
  !! @param[out]   minimize    : minimize parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_atoms_info(molecule, minimize)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_minimize),        intent(inout) :: minimize

    ! local variables
    integer :: i, iatom
    integer, parameter :: max_atm_print = 100


    ! Punch out opt atoms info.
    !
    if (main_rank) then
      write(MsgOut,'(a)') "Setup_MinAtoms> Atoms info in minimize"
      if (minimize%num_fixatoms < max_atm_print) then
        do i = 1, minimize%num_fixatoms
          iatom   = minimize%fixedatom_id(i)
          write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            iatom, &
            molecule%segment_name(iatom), &
            molecule%residue_no(iatom),   &
            molecule%residue_name(iatom), &
            molecule%atom_name(iatom),    &
            molecule%atom_cls_name(iatom)
        end do
      end if
      write(MsgOut,'(a,i0)') "  number of fixed atoms     = ", minimize%num_fixatoms
      if (minimize%num_optatoms < max_atm_print) then
        do i = 1, minimize%num_optatoms
          iatom   = minimize%optatom_id(i)
          write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            iatom, &
            molecule%segment_name(iatom), &
            molecule%residue_no(iatom),   &
            molecule%residue_name(iatom), &
            molecule%atom_name(iatom),    &
            molecule%atom_cls_name(iatom)
        end do
      end if
      write(MsgOut,'(a,i0)') "  number of minimized atoms = ", minimize%num_optatoms
      write(MsgOut, '(a)') ' '

      if (minimize%macro) then
        if (minimize%num_optatoms_macro > 0 .and. &
            minimize%num_optatoms_macro < max_atm_print) then
          write(MsgOut,'(a)') "  Atoms optimized during macro-iteration"
          do i = 1, minimize%num_optatoms_macro
            iatom   = minimize%optatom_macro_id(i)
            write(MsgOut,'(i10,x,a4,x,i6,x,a4,x,a6,x,a6)') &
              iatom, &
              molecule%segment_name(iatom), &
              molecule%residue_no(iatom),   &
              molecule%residue_name(iatom), &
              molecule%atom_name(iatom),    &
              molecule%atom_cls_name(iatom)
          end do
        end if
        write(MsgOut,'(a,i0)') "  number of atoms in macro-iteration = ", &
              minimize%num_optatoms_macro
        write(MsgOut,'(a,i0)') "  number of atoms in micro-iteration = ", &
              minimize%num_optatoms_micro
        write(MsgOut, '(a)') ' '
      end if
    end if

    return

  end subroutine print_atoms_info

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_min
  !> @brief        perform energy minimization
  !! @authors      YS
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] minimize    : minimize information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_min(output, molecule, enefunc, dynvars, minimize, &
                     pairlist, boundary)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize),         intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary


    ! Open output files
    !
    call open_output(output)

    ! MIN main loop
    !
    if (minimize%nsteps > 0) then
      select case (minimize%method)

      case (MinimizeMethodSD)

        call steepest_descent (output, molecule, enefunc, dynvars, minimize, &
                               pairlist, boundary)

      case (MinimizeMethodLBFGS)

        if (.not. minimize%macro) then
          call minimize_lbfgs (output, molecule, enefunc, dynvars, minimize, &
                               pairlist, boundary)
        else
          call macro_iter (output, molecule, enefunc, dynvars, minimize, &
                           pairlist, boundary)
        end if

        call dealloc_minimize(minimize, MinimizeLBFGS)

      end select
    end if

    ! close output files
    !
    call close_output(output)

    ! check structure
    !
    call perform_structure_check(dynvars%coord, minimize%check_structure, &
                                 .false., .false., .true.)

    return

  end subroutine run_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    steepest_descent
  !> @brief        steepest descent integrator
  !> @authors      TM, CK, KY
  !! @param[inout] output   : output information
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] minimize : minimize information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine steepest_descent(output, molecule, enefunc, dynvars, minimize, &
                              pairlist, boundary)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)                  :: energy_ref, delta_energy
    real(wp)                  :: delta_r, delta_rini, delta_rmax
    real(wp)                  :: rmsg, maxg
    integer                   :: maxg_id
    integer                   :: i, j, jid, natom, natom_opt, nsteps
    integer                   :: outunit
    character(256)            :: folder,basename
    character                 :: num*5
    logical                   :: conv
    logical                   :: savefile

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),         pointer :: force(:,:)
    real(wp),         pointer :: force_omp(:,:,:)
    integer,          pointer :: optatom_id(:)


    ! output file
    if (output%logout) then
      outunit = output%logunit
    else
      outunit = MsgOut
    end if

    ! use pointers
    !
    coord      => dynvars%coord
    coord_pbc  => dynvars%coord_pbc
    coord_ref  => dynvars%coord_ref
    force      => dynvars%force
    force_omp  => dynvars%force_omp

    nsteps     =  minimize%nsteps
    natom      =  molecule%num_atoms
    natom_opt  =  minimize%num_optatoms
    optatom_id => minimize%optatom_id

    ! initialize
    !
    delta_rini = minimize%force_scale_init
    delta_rmax = minimize%force_scale_max
    delta_r    = delta_rini
    dynvars%step = 0

    ! Compute energy of the initial structure
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter, &
                        dynvars%coord,        &
                        dynvars%trans,        &
                        dynvars%coord_pbc,    &
                        dynvars%energy,       &
                        dynvars%temporary,    &
                        dynvars%force,        &
                        dynvars%force_omp,    &
                        dynvars%virial,       &
                        dynvars%virial_extern)
    dynvars%total_pene = dynvars%energy%total

    call calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id)

    ! print the initial energy
    !
    call output_min(output, molecule, enefunc, minimize, boundary, &
                    delta_r, maxg_id, dynvars)

    ! update structure
    !
    do j = 1, natom_opt
      jid  = optatom_id(j)
      coord(1:3,jid) = coord(1:3,jid) + delta_r*force(1:3,jid)/rmsg
    end do


    ! Main loop of the Steepest descent method
    !
    conv = .false.
    do i = 1, nsteps

      dynvars%step = i

      ! save old energy and coordinates
      !
      do j = 1, natom_opt
        jid  = optatom_id(j)
        coord_ref(1:3,jid) = coord(1:3,jid)
      end do
      energy_ref = dynvars%energy%total

      ! Compute energy and forces
      !
      call compute_energy(molecule, enefunc, pairlist, boundary,  &
                          .true.,               &
                          enefunc%nonb_limiter, &
                          dynvars%coord,        &
                          dynvars%trans,        &
                          dynvars%coord_pbc,    &
                          dynvars%energy,       &
                          dynvars%temporary,    &
                          dynvars%force,        &
                          dynvars%force_omp,    &
                          dynvars%virial,       &
                          dynvars%virial_extern)

      dynvars%total_pene = dynvars%energy%total

      delta_energy = dynvars%energy%total - energy_ref
      if (delta_energy > 0) then
        delta_r = 0.5_wp*delta_r
      else
        delta_r = min(delta_rmax, 1.2_wp*delta_r)
      end if

      ! check RMSG and MAXG convergence
      !
      call calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id)
      if (rmsg <= minimize%tol_rmsg .and. maxg <= minimize%tol_maxg) then
        conv = .true.
        exit

      else if (i == nsteps) then
        conv = .false.
        exit

      end if

      ! output energy and dynamical variables
      !
      call output_min(output, molecule, enefunc, minimize, boundary, &
                      delta_r, maxg_id, dynvars)

      ! update structure
      !
      do j = 1, natom_opt
        jid  = optatom_id(j)
        coord(1:3,jid) = coord(1:3,jid) + delta_r*force(1:3,jid)/rmsg
      end do

      ! update nonbond pairlist
      !
      if (minimize%nbupdate_period > 0) then
        if (mod(i,minimize%nbupdate_period) == 0 .and. real_calc) then

          call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                               coord_pbc, pairlist)


        end if
      end if


    end do

    ! output final energy and dynamical variables
    !
    minimize%eneout = .true.
    minimize%crdout = .true.
    minimize%rstout = .true.
    call output_min(output, molecule, enefunc, minimize, boundary, &
                    delta_r, maxg_id, dynvars)
    if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%save_qminfo) &
      call write_qminfo(enefunc%qmmm, dynvars%energy%qm_ene)

    if (replica_main_rank .and. minimize%eneout_period > 0) then
      write(outunit, '(/,"Final energy = ", F20.10)') dynvars%energy%total
      if (conv) then
        write(outunit, '(/, &
        " >>>>> Minimization converged: RMSG and MAXG sufficiently small.",/)')
      else
        write(outunit, '(/, &
        " >>>>> STOP: Total number of iterations exceeds limit.",/)')
      end if
    end if
    minimize%eneout = .false.
    minimize%crdout = .false.
    minimize%rstout = .false.

    return

  end subroutine steepest_descent

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  minimize_lbfgs
  !> @brief        minimzation with lbfgs
  !> @author       CK, KY, YA
  !! @param[inout] output:     information about output [str]
  !! @param[inout] molecule:   information about molecules [str]
  !! @param[inout] enefunc:    potential energy function [str]
  !! @param[inout] dynvars:    dynamics varibales [str]
  !! @param[inout] minimize:   information about minimize [str]
  !! @param[inout] pairlist:   pairlist [str]
  !! @param[inout] boundary:   boundary condition [str]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine minimize_lbfgs (output, molecule, enefunc, dynvars, minimize, &
                             pairlist, boundary)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)               :: energy1
    real(wp)               :: rmsg, maxg
    real(wp)               :: dsave(29)
    real(wp)               :: factr, pgtol
    real(wp)               :: maxmove
    real(wp)               :: delta_r
    integer                :: i, j
    integer                :: maxg_id
    integer                :: natom
    integer                :: natom_opt, no3
    integer                :: nsteps, ncount
    integer                :: ncorr, nwork
    integer                :: isave(44)
    integer                :: iprint = -1
    integer                :: outunit
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    logical                :: bnd, bnd_qmonly

    real(wp), allocatable  :: lower(:), upper(:)
    real(wp), allocatable  :: vec(:), gradient(:), work_lbfgs(:)
    integer , allocatable  :: list_bound(:), iwork_lbfgs(:)

    integer,       pointer :: optatom_id(:)
    real(wp),      pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),      pointer :: force(:,:)
    real(wp),      pointer :: force_omp(:,:,:)



    ! output file
    !
    if (output%logout) then
      outunit = output%logunit
    else
      outunit = MsgOut
    end if

    ! dummy variable for output_min
    !
    delta_r = 0.0001_wp

    ! use pointers
    !
    coord      => dynvars%coord
    coord_pbc  => dynvars%coord_pbc
    coord_ref  => dynvars%coord_ref
    force      => dynvars%force
    force_omp  => dynvars%force_omp

    nsteps    =  minimize%nsteps
    natom     =  molecule%num_atoms

    ! extract indices of relaxed atoms
    !
    natom_opt = minimize%num_optatoms
    optatom_id => minimize%optatom_id
    no3       = natom_opt*3

    ! variables for LBFGS
    !
    allocate(vec(no3), gradient(no3), lower(no3), upper(no3))
    allocate(list_bound(no3))

    ncorr = minimize%ncorrection
    nwork = (2*no3 + 11*ncorr + 8)*ncorr + 5*no3
    allocate(iwork_lbfgs(no3*3), work_lbfgs(nwork))

    bnd        = minimize%lbfgs_bnd
    bnd_qmonly = minimize%lbfgs_bnd_qmonly
    maxmove    = minimize%lbfgs_bnd_maxmove

    if (minimize%verbose .and. main_rank) then
      iprint = 1
    end if
    factr = 0.0_wp
    pgtol = 0.0_wp

    list_bound = 0
    if (bnd) then
      if (bnd_qmonly) then
        ! set boundary to QM atoms
        do i = 1, natom_opt
          do j = 1, enefunc%qmmm%qm_natoms
             if (optatom_id(i) == enefunc%qmmm%qmatom_id(j)) then
               list_bound(3*i-2:3*i) = 2
               exit
             end if
          end do
        end do

      else
        ! set boundary to all atoms
        list_bound = 2

      end if
    end if

    task = 'START'

    ! set initial coordinates
    do i = 1, natom_opt
      vec(i*3-2:i*3) = coord(1:3,optatom_id(i))
    end do

    ! start LBFGS iteration
    !
    ncount  = -1
    energy1 = 0.0_wp
    do while(task(1:2) .eq. 'FG' .or. task(1:5) .eq. 'NEW_X' .or. &
             task(1:5) .eq. 'START')

      if (bnd) then
        if (bnd_qmonly) then
          ! set boundary to QM atoms
          do i = 1, no3
            if (list_bound(i) == 2) then
              upper(i) = vec(i) + maxmove
              lower(i) = vec(i) - maxmove
            end if
          end do

        else
          ! set boundary to all atoms
          do i = 1, no3
            upper(i) = vec(i) + maxmove
            lower(i) = vec(i) - maxmove
          end do

        end if
      end if

      call setulb(no3, ncorr, vec, lower ,upper, list_bound, energy1, &
                  gradient, factr, pgtol, work_lbfgs, iwork_lbfgs,    &
                  task, iprint, csave, lsave, isave, dsave, 0)
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(vec, no3, mpi_wp_real, 0, mpi_comm_country, ierror)
      call mpi_bcast(work_lbfgs, nwork, mpi_wp_real, 0, mpi_comm_country,ierror)
      call mpi_bcast(dsave, 29, mpi_wp_real, 0, mpi_comm_country, ierror)
#endif

      if (minimize%verbose .and. replica_main_rank) then
        write(outunit,'("task is: ",A)') trim(task)
      end if

      if (task(1:2) .eq. 'FG') then
        ncount       = ncount + 1
        dynvars%step = ncount

        ! update structure
        !
        do i = 1, natom_opt
          coord(1:3,optatom_id(i))     = vec(i*3-2:i*3)
          coord_ref(1:3,optatom_id(i)) = coord(1:3,optatom_id(i))
        end do

        ! update pairlist
        !
        if (minimize%nbupdate_period > 0) then
          if (mod(ncount,minimize%nbupdate_period) == 0 .and. real_calc) then

            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 coord_pbc, pairlist)


          end if
        end if

        ! calc energy and force
        !
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            enefunc%nonb_limiter,                     &
                            coord, dynvars%trans, coord_pbc,          &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

        energy1            = dynvars%energy%total
        dynvars%total_pene = dynvars%energy%total

        ! check RMSG and MAXG convergence
        !
        call calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id)
        if (rmsg <= minimize%tol_rmsg .and. maxg <= minimize%tol_maxg) then
          task = 'STOP: RMSG and MAXG becomes sufficiently small'
          exit

        else if (ncount == nsteps) then
          task = 'STOP: Total number of iterations exceeds limit.'
          exit

        end if

        ! output energy and dynamical variables
        !
        call output_min(output, molecule, enefunc, minimize, boundary, &
                        delta_r, maxg_id, dynvars)

        ! update gradient
        !
        do i = 1, natom_opt
          gradient(i*3-2:i*3) = -force(1:3,optatom_id(i))
        end do

      end if

    end do

    ! output final energy and dynamical variables
    !
    minimize%eneout = .true.
    minimize%crdout = .true.
    minimize%rstout = .true.
    call output_min(output, molecule, enefunc, minimize, boundary, &
                    delta_r, maxg_id, dynvars)
    if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%save_qminfo) &
      call write_qminfo(enefunc%qmmm, dynvars%energy%qm_ene)

    if (replica_main_rank .and. minimize%eneout_period > 0) then
      write(outunit, '(/,"Final energy = ", F20.10)') energy1
      write(outunit, '(/, " >>>>> ", A,/)') trim(task)
    end if
    minimize%eneout = .false.
    minimize%crdout = .false.
    minimize%rstout = .false.

    deallocate(vec, gradient, list_bound, lower, upper)
    deallocate(iwork_lbfgs, work_lbfgs)

    return

  end subroutine minimize_lbfgs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    macro_iter
  !> @brief        macro-iteration for QM atoms with QM/MM energy
  !> @authors      YA, KY
  !! @param[inout] output   : information about output [str]
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] minimize : minimize information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine macro_iter(output, molecule, enefunc, dynvars, minimize, &
                        pairlist, boundary)

    ! formal arguments
    type(s_output),            intent(inout) :: output
    type(s_molecule), target,  intent(inout) :: molecule
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),  target,  intent(inout) :: dynvars
    type(s_minimize), target,  intent(inout) :: minimize
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    ! local variables
    real(wp)               :: energy1
    real(wp)               :: rmsg, maxg, absg
    real(wp)               :: dsave(29)
    real(wp)               :: factr, pgtol
    real(wp)               :: maxmove
    real(wp)               :: energy0corr
    real(wp)               :: delta_r
    integer                :: i, j, k, ii
    integer                :: maxg_id
    integer                :: natom
    integer                :: natom_opt, natom_macro, natom_micro, no3
    integer                :: nsteps, ncount
    integer                :: ncorr, nwork
    integer                :: isave(44)
    integer                :: iprint = -1
    integer                :: outunit
    integer                :: step_save
    integer                :: crdout_period_save, rstout_period_save
    integer                :: eneout_period_save
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    logical                :: bnd, bnd_qmonly
    logical                :: eneout_short_save

    real(wp),  allocatable :: lower(:), upper(:)
    real(wp),  allocatable :: vec(:), gradient(:), work_lbfgs(:)
    integer ,  allocatable :: list_bound(:), iwork_lbfgs(:)
    real(wp),  allocatable :: coord0(:,:), force0corr(:,:)

    integer,       pointer :: optatom_id(:)
    integer,       pointer :: optatom_macro_id(:)
    integer,       pointer :: optatom_micro_id(:)
    real(wp),      pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),      pointer :: force(:,:)
    real(wp),      pointer :: force_omp(:,:,:)


    ! initialize
    !
    minimize%eneout = .false.
    minimize%crdout = .false.
    minimize%rstout = .false.

    ! output file
    !
    if (output%logout) then
      outunit = output%logunit
    else
      outunit = MsgOut
    end if

    ! dummy variable for output_min
    !
    delta_r = 0.0001_wp

    ! use pointers
    !
    coord     => dynvars%coord
    coord_pbc => dynvars%coord_pbc
    coord_ref => dynvars%coord_ref
    force     => dynvars%force
    force_omp => dynvars%force_omp

    nsteps    =  minimize%nsteps
    natom     =  molecule%num_atoms

    ! extract indices of relaxed atoms
    !
    natom_opt         =  minimize%num_optatoms
    optatom_id        => minimize%optatom_id
    natom_macro       =  minimize%num_optatoms_macro
    optatom_macro_id  => minimize%optatom_macro_id
    no3               =  natom_macro*3
    natom_micro       =  minimize%num_optatoms_micro
    optatom_micro_id  => minimize%optatom_micro_id

    ! variables for LBFGS
    !
    allocate(vec(no3), gradient(no3), lower(no3), upper(no3))
    allocate(list_bound(no3))

    ncorr = minimize%ncorrection
    nwork = (2*no3 + 11*ncorr + 8)*ncorr + 5*no3
    allocate(iwork_lbfgs(no3*3), work_lbfgs(nwork))

    bnd        = minimize%lbfgs_bnd
    bnd_qmonly = minimize%lbfgs_bnd_qmonly
    maxmove    = minimize%lbfgs_bnd_maxmove

    if (minimize%verbose .and. replica_main_rank) then
      iprint = 1
    end if
    factr = 0.0_wp
    pgtol = 0.0_wp

    if (bnd) then
      ! set boundary to all atoms in macro region
      list_bound = 2
    else
      list_bound = 0
    end if

    task = 'START'

    ! set initial coordinates
    ii = 1
    do i = 1, natom_macro
      vec(ii:ii+2) = coord(1:3,optatom_macro_id(i))
      ii = ii + 3
    end do

    ! correction for micro_iteration
    allocate(force0corr(3,natom_micro), coord0(3,natom_micro))
    coord0      = 0.0_wp
    energy0corr = 0.0_wp
    force0corr  = 0.0_wp

    ! start macro-iteration
    !
    ncount  = -1
    energy1 = 0.0_wp
    do while(task(1:2) .eq. 'FG' .or. task(1:5) .eq. 'NEW_X' .or. &
             task(1:5) .eq. 'START')

      if (bnd) then
        ! set boundary to all atoms in macro region
        do i = 1, no3
          upper(i) = vec(i) + maxmove
          lower(i) = vec(i) - maxmove
        end do
      end if

      call setulb(no3, ncorr, vec, lower ,upper, list_bound, energy1, &
                  gradient, factr, pgtol, work_lbfgs, iwork_lbfgs,    &
                  task, iprint, csave, lsave, isave, dsave, 0)
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(vec, no3, mpi_wp_real, 0, mpi_comm_country, ierror)
      call mpi_bcast(work_lbfgs, nwork, mpi_wp_real, 0, mpi_comm_country, ierror)
      call mpi_bcast(dsave, 29, mpi_wp_real, 0, mpi_comm_country, ierror)
#endif

      if (minimize%verbose .and. replica_main_rank) then
        write(outunit,'("task is: ",A)') trim(task)
      end if

      if (task(1:2) .eq. 'FG') then
        ncount       = ncount + 1
        dynvars%step = ncount

        ! update structure
        !
        do i = 1, natom_macro
          coord(1:3,optatom_macro_id(i))     = vec(i*3-2:i*3)
          coord_ref(1:3,optatom_macro_id(i)) = coord(1:3,optatom_macro_id(i))
        end do

        ! micro iteration
        !   Skip micro-iteration during the first few iterations
        if (ncount > minimize%start_micro) then

          step_save          = dynvars%step
          eneout_period_save = minimize%eneout_period
          crdout_period_save = minimize%crdout_period
          rstout_period_save = minimize%rstout_period

          ! turn off coordinate output
          minimize%crdout_period = 0
          minimize%rstout_period = 0

          if (minimize%eneout_period > 0) then
            if (mod(dynvars%step,minimize%eneout_period) == 0) then
              ! turn on short energy output
              eneout_short_save  = minimize%eneout_short
            else
              ! turn off energy output
              minimize%eneout_period = 0
            end if
          end if

          ! turn on short energy output
          minimize%eneout_short  = .true.
          call micro_iter(output, molecule, enefunc, dynvars, minimize, &
               pairlist, boundary, coord0, energy0corr, force0corr, ncount)

          dynvars%step           = step_save
          minimize%eneout_period = eneout_period_save
          minimize%crdout_period = crdout_period_save
          minimize%rstout_period = rstout_period_save
          minimize%eneout_short  = eneout_short_save

        end if

        ! update pairlist
        !
        if (minimize%nbupdate_period > 0) then
          if (mod(ncount,minimize%nbupdate_period) == 0 .and. real_calc) then

            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 coord_pbc, pairlist)


          end if
        end if

        ! calc energy and force
        !
        enefunc%qmmm%qm_classical = .false.
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            enefunc%nonb_limiter,                     &
                            coord, dynvars%trans, coord_pbc,          &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

        energy1 = dynvars%energy%total
        dynvars%total_pene = dynvars%energy%total

        ! check RMSG and MAXG convergence
        !
        call calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id)

        if (rmsg <= minimize%tol_rmsg .and. maxg <= minimize%tol_maxg) then
          task = 'STOP: RMSG and MAXG becomes sufficiently small'
          exit

        else if (ncount >= nsteps) then
          task = 'STOP: Total number of iterations exceeds limit.'
          exit

        end if

        ! output energy and dynamical variables
        !
        call output_min(output, molecule, enefunc, minimize, boundary, &
                        delta_r, maxg_id, dynvars)

        ! update gradient
        !
        ii = 1
        do i = 1, natom_macro
          gradient(ii:ii+2) = -force(1:3,optatom_macro_id(i))
          ii = ii+3
        end do

        ! energy and gradient correction terms for micro_iteration
        !
        energy0corr = dynvars%energy%total
        do i = 1, natom_micro
          coord0(1:3,i)     = coord(1:3,optatom_micro_id(i))
          force0corr(1:3,i) = force(1:3,optatom_micro_id(i))
        end do

        enefunc%qmmm%qm_classical = .true.
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            enefunc%nonb_limiter,                     &
                            coord, dynvars%trans, coord_pbc,          &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

        energy0corr = energy0corr - dynvars%energy%total
        do i = 1, natom_micro
          force0corr(1:3,i) = force0corr(1:3,i) - force(1:3,optatom_micro_id(i))
        end do

      end if

    end do
    
    ! output final energy and dynamical variables
    !
    minimize%eneout = .true.
    minimize%crdout = .true.
    minimize%rstout = .true.
    call output_min(output, molecule, enefunc, minimize, boundary, &
                    delta_r, maxg_id, dynvars)
    if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%save_qminfo) &
      call write_qminfo(enefunc%qmmm, dynvars%energy%qm_ene)

    if (replica_main_rank .and. minimize%eneout_period > 0) then
      write(outunit, '(/,"Final energy = ", F20.10)') energy1
      write(outunit, '(/, " >>>>> ", A,/)') trim(task)
    end if

    deallocate(vec, gradient, list_bound, lower, upper)
    deallocate(iwork_lbfgs, work_lbfgs)
    deallocate(force0corr, coord0)

    return

  end subroutine macro_iter

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    micro_iter
  !> @brief        micro-iteration for MM atoms with MM energy
  !> @authors      YA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] minimize    : minimize information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[in]    energy0corr : correction of the energy
  !! @param[in]    force0corr  : correction of the force
  !! @param[in]    iter        : counter of the macro iteration
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine micro_iter(output, molecule, enefunc, dynvars, minimize, &
                        pairlist, boundary, coord0, energy0corr, force0corr, iter)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize), target, intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    real(wp),                 intent(in)    :: coord0(3,*)
    real(wp),                 intent(in)    :: energy0corr
    real(wp),                 intent(in)    :: force0corr(3,*)
    integer, optional,        intent(in)    :: iter

    ! local variables
    real(wp)                 :: rmsg, maxg, absg
    real(wp)                 :: factr, pgtol
    real(wp)                 :: maxmove
    real(wp)                 :: energy1, gdx
    real(wp)                 :: delta_r
    integer                  :: i, j, k
    integer                  :: maxg_id
    integer                  :: natom
    integer                  :: natom_opt, no3
    integer                  :: nsteps, ncount
    integer                  :: ncorr
    integer                  :: iprint = -1
    integer                  :: indx
    integer                  :: outunit
    logical                  :: print_log
    character(len=60)        :: task, csave
    logical                  :: bnd, bnd_qmonly

    real(wp),        pointer :: dsave(:)
    real(wp),        pointer :: lower(:), upper(:)
    real(wp),        pointer :: vec(:), gradient(:), work_lbfgs(:)
    real(wp),        pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),        pointer :: force(:,:)
    real(wp),        pointer :: force_omp(:,:,:)
    integer,         pointer :: optatom_id(:)
    integer,         pointer :: isave(:)
    integer,         pointer :: list_bound(:), iwork_lbfgs(:)
    logical,         pointer :: lsave(:)


    ! output file
    if (output%logout) then
      outunit = output%logunit
    else
      outunit = MsgOut
    end if

    print_log = .false.
    if (replica_main_rank .and. minimize%eneout_period > 0) then
      print_log = .true.
      if (present(iter)) then
        write(outunit, '(A28,i0,A7,/)') &
          "+++++ Start micro iteration ", iter, " +++++ "
      else
        write(outunit, '(A,/)') &
          "+++++ Start micro iteration +++++ "
      end if
    end if


    ! dummy variable for output_min
    !
    delta_r = 0.0001_wp

    ! use pointers
    !
    coord     => dynvars%coord
    coord_pbc => dynvars%coord_pbc
    coord_ref => dynvars%coord_ref
    force     => dynvars%force
    force_omp => dynvars%force_omp

    nsteps  =  minimize%nsteps_micro
    natom   =  molecule%num_atoms

    ! extract indices of relaxed atoms
    !
    natom_opt  =  minimize%num_optatoms_micro
    optatom_id => minimize%optatom_micro_id
    no3        =  natom_opt*3

    ! variables for LBFGS
    csave       =  minimize%csave_micro
    lsave       => minimize%lsave_micro
    isave       => minimize%isave_micro
    dsave       => minimize%dsave_micro
    vec         => minimize%vec_micro
    gradient    => minimize%gradient_micro
    lower       => minimize%lower_micro
    upper       => minimize%upper_micro
    work_lbfgs  => minimize%work_lbfgs_micro
    iwork_lbfgs => minimize%iwork_lbfgs_micro
    list_bound  => minimize%list_bound_micro

    ncorr  = minimize%ncorrection

    bnd        = minimize%lbfgs_bnd
    bnd_qmonly = minimize%lbfgs_bnd_qmonly
    maxmove    = minimize%lbfgs_bnd_maxmove

    if (minimize%verbose .and. replica_main_rank) then
      iprint = 1
    end if
    factr = 0.0_wp
    pgtol = 0.0_wp

    if (bnd .and. .not. bnd_qmonly) then
      ! set boundary to all atoms in micro region
      list_bound = 2
    else
      list_bound = 0
    end if

    task = 'START'

    ! set initial coordinates
    do i = 1, natom_opt
       vec(i*3-2:i*3)     = coord(1:3,optatom_id(i))
    end do

    ! start micro-iteration
    !
    ncount  = 0
    energy1 = 0.0_wp
    do while(task(1:2) .eq. 'FG' .or. task(1:5) .eq. 'NEW_X' .or. &
             task(1:5) .eq. 'START')

      if (bnd .and. .not. bnd_qmonly) then
        ! set boundary to all atoms in micro region
        do i = 1, no3
          upper(i) = vec(i) + maxmove
          lower(i) = vec(i) - maxmove
        end do
      end if

      call setulb(no3, ncorr, vec, lower ,upper, list_bound, energy1, &
                   gradient, factr, pgtol, work_lbfgs, iwork_lbfgs,   &
                   task, iprint, csave, lsave, isave, dsave, 0)
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(vec, no3, mpi_wp_real, 0, mpi_comm_country, ierror)
      call mpi_bcast(work_lbfgs, size(work_lbfgs), mpi_wp_real, 0, mpi_comm_country, ierror)
      call mpi_bcast(dsave, 29, mpi_wp_real, 0, mpi_comm_country, ierror)
#endif

      if (minimize%verbose .and. print_log) then
        write(outunit,'("task is: ",A)') trim(task)
      end if
      
      if (task(1:2) .eq. 'FG') then
        ncount       = ncount + 1
        dynvars%step = ncount

        ! update structure
        !
        do i = 1, natom_opt
          coord(1:3,optatom_id(i))     = vec(i*3-2:i*3) 
          coord_ref(1:3,optatom_id(i)) = coord(1:3,optatom_id(i))
        end do

        ! update pairlist
        if (minimize%nbupdate_period > 0) then
          if (mod(ncount,minimize%nbupdate_period) == 0 .and. real_calc) then

            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 coord_pbc, pairlist)


          end if
        end if

        ! calc energy and force
        !
        enefunc%qmmm%qm_classical = .true.
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            enefunc%nonb_limiter,                     &
                            coord, dynvars%trans, coord_pbc,          &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

        ! correction to energy and force
        !
        gdx = 0.0_wp
        do i = 1, natom_opt
          indx = optatom_id(i)
          gdx = gdx - (force0corr(1,i) * (coord(1,indx) - coord0(1,i)) &
                     + force0corr(2,i) * (coord(2,indx) - coord0(2,i)) &
                     + force0corr(3,i) * (coord(3,indx) - coord0(3,i)))
          force(1:3,indx) = force(1:3,indx) + force0corr(1:3,i)
        end do

        dynvars%energy%total = dynvars%energy%total + energy0corr + gdx
        energy1 = dynvars%energy%total

        ! check RMSG and MAXG convergence
        !
        call calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id, .true.)

        if (rmsg <= minimize%tol_rmsg_micro .and. &
            maxg <= minimize%tol_maxg_micro) then
          task = "STOP: RMSG and MAXG becomes sufficiently small."
          exit

        else if (ncount >= nsteps) then
          task = 'STOP: Total number of iterations exceeds limit.'
          exit

        end if

        ! output energy and dynamical variables
        !
        call output_min(output, molecule, enefunc, minimize, boundary, &
                        delta_r, maxg_id, dynvars)

        ! update gradient
        !
        do i = 1, natom_opt
          gradient(i*3-2:i*3) = -force(1:3,optatom_id(i))
        end do

      end if

    end do
    
    ! output final energy and dynamical variables
    !
    minimize%eneout = .true.
    call output_min(output, molecule, enefunc, minimize, boundary, &
                    delta_r, maxg_id, dynvars)
    if (print_log) then
      write(outunit, '(/,"+++++ ", A,"+++++"/)') trim(task)
    end if
    minimize%eneout = .false.

    return

  end subroutine micro_iter

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    add_fixatm
  !> @brief        add atoms to be fixed in MINIMIZATION
  !> @authors      KY
  !! @param[in]    natoms   : number of atoms
  !! @param[in]    atom_id  : list of atoms
  !! @param[in]    molecule : molecule information
  !! @param[inout] minimize : minimize information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine add_fixatm(natoms, atom_id, molecule, minimize)

    ! formal arguments
    integer,                 intent(in)    :: natoms
    integer,                 intent(in)    :: atom_id(natoms)
    type(s_molecule),        intent(in)    :: molecule
    type(s_minimize),        intent(inout) :: minimize

    ! local variables
    integer                  :: nd, nf0, i, j, k
    integer,     allocatable :: del_id(:), tmp(:)


    nd = minimize%num_optatoms
    allocate(del_id(nd))
    call delete_element(natoms, atom_id, minimize%num_optatoms, &
                        minimize%optatom_id, nd, del_id)

    if (nd == 0) return

    if(main_rank) then
      write(MsgOut,'(a)') "Setup_MinAtoms> Atoms info is changed in minimize"
      write(MsgOut,'(a,i0)') "  number of atoms deleted from opt = ", nd
      write(MsgOut,*)
    end if

    ! reallocate fixedatom_id
    if (minimize%num_fixatoms > 0) then
      allocate(tmp(minimize%num_fixatoms))
      tmp = minimize%fixedatom_id
      nf0 = minimize%num_fixatoms

      deallocate(minimize%fixedatom_id)
      minimize%num_fixatoms = minimize%num_fixatoms + nd
      allocate(minimize%fixedatom_id(minimize%num_fixatoms))

      j = 1
      k = 1
      do i = 1, minimize%num_fixatoms
        if(tmp(j) < del_id(k)) then
          minimize%fixedatom_id(i) = tmp(j)
          j = j + 1
          if (j > nf0) then
            minimize%fixedatom_id(i+1:minimize%num_fixatoms) = &
              del_id(k:nd)
            exit
          end if
        else
          minimize%fixedatom_id(i) = del_id(k)
          k = k + 1
          if (k > nd) then
            minimize%fixedatom_id(i+1:minimize%num_fixatoms) = &
              tmp(j:nf0)
            exit
          end if
        end if
      end do

      deallocate(tmp)

    else
      minimize%num_fixatoms = nd
      allocate(minimize%fixedatom_id(minimize%num_fixatoms))
      minimize%fixedatom_id = del_id

    end if

    if (minimize%macro) then
      call delete_element(natoms, atom_id, minimize%num_optatoms_macro, &
                          minimize%optatom_macro_id, nd, del_id)
      call delete_element(natoms, atom_id, minimize%num_optatoms_micro, &
                          minimize%optatom_micro_id, nd, del_id)
    end if

    call print_atoms_info(molecule, minimize)

    deallocate(del_id)

    return

  end subroutine add_fixatm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    delete_element
  !> @brief        deletes an element of vector b if it is an element of
  !>               vector a. Note that a and b are assumed to be in ascending 
  !>               order, i.e., a(1) < a(2) < ... and b(1) < b(2) < ...
  !> @authors      KY
  !! @param[in]    na     : number of elements of a
  !! @param[in]    a(na)  : vector a
  !! @param[inout] nb     : number of elements of b
  !! @param[inout] b(nb)  : vector b
  !! @param[inout] nd     : number of deleted elements
  !! @param[inout] d(nd)  : vector of deleted elements
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine delete_element(na, a, nb, b, nd, d)

    ! formal arguments
    integer,                 intent(in)    :: na
    integer,                 intent(in)    :: a(na)
    integer,                 intent(inout) :: nb
    integer,                 intent(inout) :: b(*)
    integer,       optional, intent(inout) :: nd
    integer,       optional, intent(inout) :: d(*)

    ! local variables
    integer                  :: i, j, k, nb0
    integer,     allocatable :: tmp(:)
    logical,     allocatable :: del(:)


    allocate(del(nb))
    nd = 0
    do i = 1, nb
      del(i) = .false.
      do j = 1, na
        if (a(j) == b(i)) then
          nd = nd + 1
          del(i) = .true.
          exit
        else if (a(j) > b(i)) then
          exit
        end if
      end do
    end do

    if(nd == 0) then
      deallocate(del)
      return
    end if

    allocate(tmp(nb))
    tmp = b(1:nb)
    nb0 = nb
    nb = nb - nd

    j = 1
    k = 1
    do i = 1, nb0
      if(.not. del(i)) then
        b(j) = tmp(i)
        j = j + 1
      else
        d(k) = tmp(i)
        k = k + 1
      end if
    end do

    deallocate(del,tmp)

    return

  end subroutine delete_element

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_rmsg
  !> @brief        Calculate the RMSG of target atoms for minimization
  !> @authors      KY
  !! @param[inout] minimize : minimize information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] rmsg     : RMSG    
  !! @param[inout] maxg     : MAXG
  !! @param[inout] maxg_id  : The atom ID of MAXG
  !! @param[in]    micro    : true, if micro-iteration
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_rmsg(minimize, dynvars, rmsg, maxg, maxg_id, micro)

    ! formal arguments
    type(s_minimize), target, intent(inout) :: minimize
    type(s_dynvars),  target, intent(inout) :: dynvars
    real(wp),                 intent(inout) :: rmsg, maxg
    integer,                  intent(inout) :: maxg_id
    logical,        optional, intent(in)    :: micro

    ! local variables
    integer  :: i, j, natom_opt
    real(wp) :: absg

    real(wp), pointer :: force(:,:)
    integer,  pointer :: optatom_id(:)


    force      => dynvars%force

    if (present(micro)) then
        if(micro) then
            natom_opt  =  minimize%num_optatoms_micro
            optatom_id => minimize%optatom_micro_id
        end if
    else
      natom_opt  =  minimize%num_optatoms
      optatom_id => minimize%optatom_id
    end if

    rmsg = 0.0_wp
    maxg = 0.0_wp

    do i = 1, natom_opt
      do j = 1, 3
        rmsg = rmsg + force(j,optatom_id(i)) * force(j,optatom_id(i))
        absg = abs(force(j,optatom_id(i)))
        if (absg > maxg) then 
          maxg    = absg
          maxg_id = optatom_id(i)
        end if
      end do
    end do
    rmsg = sqrt(rmsg/real(natom_opt*3,wp))

    dynvars%rms_gradient = rmsg
    dynvars%max_gradient = maxg

    return

  end subroutine calc_rmsg

end module at_minimize_mod
