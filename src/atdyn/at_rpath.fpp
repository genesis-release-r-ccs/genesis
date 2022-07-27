!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_rpath_mod
!> @brief   String method using replicas
!! @authors Yasuaki Komuro (YK), Takaharu Mori (TM), Norio Takase (NT)
!!          Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_rpath_mod

  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_rpath_str_mod
  use at_rpath_mep_mod
  use at_rpath_fep_mod
  use at_restraints_str_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_enefunc_restraints_mod
  use at_minimize_str_mod
  use at_output_str_mod
  use at_output_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_rstmep_mod
  use fileio_control_mod
  use select_mod
  use fitting_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  use at_input_mod
  use molecules_mod
  use at_qmmm_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_rst_mod
  use fileio_psf_mod
  use fileio_gpr_mod
  use fileio_par_mod
  use fileio_top_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_mode_mod
  use fileio_rstmep_mod
  use fileio_eef1_mod
  use fileio_spot_mod  

  implicit none
  private

  ! structures
  type, public :: s_rpath_info
    ! General variables
    integer                           :: rpathmode    = RpathmodeMFEP
    integer                           :: ncycle       = 1000
    integer                           :: nreplica     = 0
    real(wp)                          :: delta        = -1.0_wp
    logical                           :: fix_terminal = .false.
    integer                           :: dimension    = 0

    ! MFEP
    logical                           :: avoid_shrinkage = .false.
    integer                           :: rpath_period = 0
    real(wp)                          :: smooth       = 0.0_wp
    logical                           :: use_restart  = .true.
    integer,              allocatable :: rest_function(:)

    ! MEP
    character(256)                    :: mepatm_select_index = ''
    logical                           :: mep_partial_opt = .true.
    integer                           :: eneout_period = 10
    integer                           :: crdout_period = 10
    integer                           :: rstout_period = 10
    real(wp)                          :: tol_energy = 0.01_wp
    real(wp)                          :: tol_path   = 0.01_wp
    logical                           :: massWeightCoord = .false.

    integer                           :: method = MEPmethod_String

    ! MEP-String 
    !   method = string
    !   delta  = 0.001

    ! MEP-NEB
    real(wp)                          :: k_spring = 10.0_wp
    integer                           :: ncorrection        = 10
    logical                           :: lbfgs_bnd          = .true.
    logical                           :: lbfgs_bnd_qmonly   = .false.
    real(wp)                          :: lbfgs_bnd_maxmove  = 0.1_wp
    logical                           :: climbing_image     = .false.
    logical                           :: verbose            = .false.

    ! to be deprecated
    real(wp)                          :: force_scale_init = 0.01_wp
    real(wp)                          :: force_scale_max  = 0.1_wp
    real(wp)                          :: tol_rmsg = 0.36_wp
    real(wp)                          :: tol_maxg = 0.54_wp
    real(wp)                          :: tol_rmsg_cineb     = -1.0_wp
    character(256)                    :: optimizer = 'GLBFGS'
    ! to be deprecated

    ! FEP
    integer                           :: fep_period = 10
    logical                           :: esp_energy = .true.
    logical                           :: esp_md     = .true.

  end type s_rpath_info

  ! subroutines
  public  :: show_ctrl_rpath
  public  :: read_ctrl_rpath
  public  :: define_nreplica
  public  :: setup_rpath
  public  :: run_rpath

  ! private subroutines
  private :: setup_rpath_mfep
  private :: setup_rpath_restart
  private :: assign_condition
  private :: run_rpath_mfep
  private :: evolve
  private :: reparametrize

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_rpath
  !> @brief        show RPATH section usage
  !! @authors      NT, TM, YK, KY
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
        write(MsgOut,'(A)') 'rpathmode           = MFEP # MFEP or MEP'
        write(MsgOut,'(A)') 'nreplica            = 4    # number of images'
        write(MsgOut,'(A)') 'fix_terminal        = NO'

        write(MsgOut,'(A)') '# For MFEP'
        write(MsgOut,'(A)') 'rpath_period        = 1000'
        write(MsgOut,'(A)') 'delta               = 0.05'
        write(MsgOut,'(A)') 'smooth              = 0.1'
        write(MsgOut,'(A)') 'rest_function       = 1 2'
        write(MsgOut,'(A)') 'use_restart         = YES'
        write(MsgOut,'(A)') 'avoid_shrinkage      = NO'

        write(MsgOut,'(A)') '# For MEP'
        write(MsgOut,'(A)') '#method              = NEB # STRING or NEB'
        write(MsgOut,'(A)') '#ncycle              = 1000'
        write(MsgOut,'(A)') '#mepatm_select_index = 1 2'
        write(MsgOut,'(A)') '#eneout_period       = 10'
        write(MsgOut,'(A)') '#crdout_period       = 10'
        write(MsgOut,'(A)') '#rstout_period       = 10'
        write(MsgOut,'(A)') '#tol_energy          = 0.01'
        write(MsgOut,'(A)') '#tol_path            = 0.01'
        write(MsgOut,'(A)') '#massWeightCoord     = YES'
        write(MsgOut,'(A)') '# For MEP/String'
        write(MsgOut,'(A)') '#delta               = 0.001'
        write(MsgOut,'(A)') '# For MEP/NEB'
        write(MsgOut,'(A)') '#k_spring            = 100.0'

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
  !! @authors      TM, YK, KY
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

    call read_ctrlfile_type   (handle, Section, 'rpathmode',       &
                               rpath_info%rpathmode, RpathmodeTypes)
    call read_ctrlfile_integer(handle, Section, 'ncycle',  &
                               rpath_info%ncycle)
    call read_ctrlfile_integer(handle, Section, 'nreplica',       &
                               rpath_info%nreplica)
    call read_ctrlfile_real   (handle, Section, 'delta',  &
                               rpath_info%delta)
    call read_ctrlfile_logical(handle, Section, 'fix_terminal',  &
                               rpath_info%fix_terminal)

    ! MFEP
    call read_ctrlfile_logical(handle, Section, 'avoid_shrinkage',  &
                               rpath_info%avoid_shrinkage)
    call read_ctrlfile_integer(handle, Section, 'rpath_period', &
                               rpath_info%rpath_period)
    call read_ctrlfile_real   (handle, Section, 'smooth', &
                               rpath_info%smooth)
    call read_ctrlfile_logical(handle, Section, 'use_restart',  &
                               rpath_info%use_restart)
    call read_ctrlfile_string (handle, Section, 'rest_function', &
                               rest_function_char)

    ! MEP
    call read_ctrlfile_string (handle, Section, 'mepatm_select_index',  &
                               rpath_info%mepatm_select_index)
    call read_ctrlfile_logical(handle, Section, 'mep_partial_opt',  &
                               rpath_info%mep_partial_opt)
    call read_ctrlfile_integer(handle, Section, 'eneout_period',  &
                               rpath_info%eneout_period)
    call read_ctrlfile_integer(handle, Section, 'crdout_period',  &
                               rpath_info%crdout_period)
    call read_ctrlfile_integer(handle, Section, 'rstout_period',  &
                               rpath_info%rstout_period)
    call read_ctrlfile_type   (handle, Section, 'method', &
                               rpath_info%method, MEPmethodTypes)

    ! String
    call read_ctrlfile_real   (handle, Section, 'tol_energy',  &
                               rpath_info%tol_energy)
    call read_ctrlfile_real   (handle, Section, 'tol_path',  &
                               rpath_info%tol_path)
    call read_ctrlfile_logical(handle, Section, 'massWeightCoord',  &
                               rpath_info%massWeightCoord)

    ! NEB
    call read_ctrlfile_real   (handle, Section, 'tol_rmsg',  &
                               rpath_info%tol_rmsg)
    call read_ctrlfile_real   (handle, Section, 'tol_maxg',  &
                               rpath_info%tol_maxg)
    call read_ctrlfile_real   (handle, Section, 'k_spring',  &
                               rpath_info%k_spring)
    call read_ctrlfile_logical(handle, Section, 'climbing_image',  &
                               rpath_info%climbing_image)
    call read_ctrlfile_real   (handle, Section, 'tol_rmsg_cineb',  &
                               rpath_info%tol_rmsg_cineb)
    call read_ctrlfile_integer(handle, Section, 'ncorrection',  &
                               rpath_info%ncorrection)
    call read_ctrlfile_logical(handle, Section, 'lbfgs_bnd',  &
                               rpath_info%lbfgs_bnd)
    call read_ctrlfile_logical(handle, Section, 'lbfgs_bnd_qmonly',  &
                               rpath_info%lbfgs_bnd_qmonly)
    call read_ctrlfile_real   (handle, Section, 'lbfgs_bnd_maxmove',  &
                               rpath_info%lbfgs_bnd_maxmove)
    call read_ctrlfile_logical(handle, Section, 'verbose',  &
                               rpath_info%verbose)

    ! to be deprecated
    call read_ctrlfile_string (handle, Section, 'optimizer', &
                               rpath_info%optimizer)
    call read_ctrlfile_real   (handle, Section, 'force_scale_init',  &
                               rpath_info%force_scale_init)
    call read_ctrlfile_real   (handle, Section, 'force_scale_max',  &
                               rpath_info%force_scale_max)
    ! to be deprecated

    ! FEP
    call read_ctrlfile_integer(handle, Section, 'fep_period',  &
                               rpath_info%fep_period)
    call read_ctrlfile_logical(handle, Section, 'esp_energy',  &
                               rpath_info%esp_energy)
    call read_ctrlfile_logical(handle, Section, 'esp_md',  &
                               rpath_info%esp_md)

    call end_ctrlfile_section(handle)

    ! error check
    !
    if ((rpath_info%rpathmode /= RpathmodeMFEP) .and. &
        (rpath_info%rpathmode /= RpathmodeMEP ))      &
      call error_msg('Read_Ctrl_Path> rpathmode should be MFEP or MEP in [RPATH]')
    if (rpath_info%nreplica == 0) &
      call error_msg('Read_Ctrl_Rpath> nreplica should be > 0 in [RPATH]')

    if (rpath_info%rpathmode == RpathmodeMFEP) then
      if (rpath_info%delta < EPS .and. rpath_info%rpath_period > 0) &
      call error_msg('Read_Ctrl_Rpath> delta should be larger than 0 if rpath_period > 0 ')
    
     
      if (main_rank) then
        if (rpath_info%rpath_period == 0) then
          write(MsgOut,'(A)') 'Read_Ctrl_Rpath> delta and smooth are not used'
        end if
      end if

      ! get dimension from the number of numerics in rest_function_char
      !
      rpath_info%dimension = split_num(rest_function_char)

      ! read rest_function
      !
      allocate(rpath_info%rest_function(rpath_info%dimension))
      rpath_info%rest_function(1:rpath_info%dimension) = 0
      call split(rpath_info%dimension, rpath_info%dimension, &
                 rest_function_char, rpath_info%rest_function)

    else
      if (rpath_info%avoid_shrinkage) &
        call error_msg('Read_Ctrl_Rpath> avoid_shrinkage is allowed only in MFEP [RPATH].')


      if (len(trim(rpath_info%mepatm_select_index)) == 0) then
        call error_msg('Read_Ctrl_Rpath> mepatm_select_index is not found in [RPATH].')
      end if

      ! default value for delta
      if (rpath_info%delta < EPS) rpath_info%delta = 0.001_wp
      if (rpath_info%tol_rmsg_cineb < 0.0_wp) rpath_info%tol_rmsg_cineb = rpath_info%tol_rmsg * 3.0_wp

      rpath_info%dimension = 0
      rpath_info%use_restart = .false.

    end if

    if (rpath_info%rpathmode == RpathmodeFEP) then
      if (.not. rpath_info%esp_md .and. rpath_info%esp_energy) then
        if (main_rank) write(MsgOut,'(A)') &
            'Read_Ctrl_Rpath> warning! esp_energy is turned off.'
        rpath_info%esp_energy = .false.
      end if
    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Rpath> Rpath information'
      write(MsgOut,'(A,A10,A,I10)') &
        '  rpathmode         = ', trim(RpathmodeTypes(rpath_info%rpathmode)), &
        '  nreplica          = ', rpath_info%nreplica

      if (rpath_info%fix_terminal) then
        write(MsgOut,'(A)') '  fix_terminal      =        yes'
      else
        write(MsgOut,'(A)') '  fix_terminal      =         no'
      end if

      if (rpath_info%rpathmode == RpathmodeMFEP) then
        if (rpath_info%fix_terminal) then
          write(MsgOut,'(A)') '  avoid_shrinkage =        yes'
        end if
        write(MsgOut,'(A,I10,A,I10)') &
          '  dimension         = ', rpath_info%dimension, &
          '  rpath_period      = ', rpath_info%rpath_period
        write(MsgOut,'(A,F10.5)') &
          '  delta             = ', rpath_info%delta
        if (rpath_info%use_restart) then
          write(MsgOut,'(A)') '  use_restart       =        yes'
        else
          write(MsgOut,'(A)') '  use_restart       =         no'
        end if
        write(MsgOut,'(A,F10.3,A,A10)') &
          '  smooth            = ', rpath_info%smooth,      &
          '  rest_function     = ', trim(rest_function_char)
      end if

      if (rpath_info%rpathmode == RpathmodeMEPMD) then
        write(MsgOut,'(A,A10,A,I10)') &
          '  method            = ', trim(MEPmethodTypes(rpath_info%method)), &
          '  rpath_period      = ', rpath_info%rpath_period
        write(MsgOut,'(A,F10.2,A,F10.2)') &
          '  tol_energy        = ', rpath_info%tol_energy,  &
          '  tol_path          = ', rpath_info%tol_path
        if (rpath_info%massWeightCoord) then
          write(MsgOut,'(A)') &
          '  massweightcoord   =        yes'
        else
          write(MsgOut,'(A)') &
          '  massweightcoord   =         no'
        end if

        if (rpath_info%method == MEPmethod_NEB) then
          ! Todo

        else if (rpath_info%method == MEPmethod_String) then
          write(MsgOut,'(A,F10.2)') &
            '  delta             = ', rpath_info%delta

        end if

        write(MsgOut,'(A,A)') &
          '  mepatm_select_index   = ', trim(rpath_info%mepatm_select_index)

      end if

      if (rpath_info%rpathmode == RpathmodeMEP) then
        write(MsgOut,'(A,A10,A,I10)') &
          '  method            = ', trim(MEPmethodTypes(rpath_info%method)), &
          '  ncycle            = ', rpath_info%ncycle
        write(MsgOut,'(A,F10.2,A,F10.2)') &
          '  tol_energy        = ', rpath_info%tol_energy,  &
          '  tol_path          = ', rpath_info%tol_path
        write(MsgOut,'(A,I10,A,I10)') &
          '  eneout_period     = ', rpath_info%eneout_period, &
          '  crdout_period     = ', rpath_info%crdout_period
        write(MsgOut,'(A,I10)') &
          '  rstout_period     = ', rpath_info%rstout_period
        if (rpath_info%massWeightCoord) then
          write(MsgOut,'(A)') &
          '  massweightcoord   =        yes'
        else
          write(MsgOut,'(A)') &
          '  massweightcoord   =         no'
        end if

        if (rpath_info%method == MEPmethod_NEB) then
          write(MsgOut,'(A,F10.2)')         &
            '  k_spring          = ', rpath_info%k_spring
          write(MsgOut,'(A,I10)') &
            '  ncorrection       = ', rpath_info%ncorrection

          if (rpath_info%lbfgs_bnd) then
            if (rpath_info%lbfgs_bnd_qmonly) then
              write(MsgOut,'(A,A)')         &
              '  lbfgs_bnd         =        yes', &
              '  lbfgs_bnd_qmonly  =        yes'
            else
              write(MsgOut,'(A,A)')         &
              '  lbfgs_bnd         =        yes', &
              '  lbfgs_bnd_qmonly  =         no'
            end if
            if (rpath_info%lbfgs_bnd_maxmove < EPS) then 
              write(MsgOut,'(A,F10.2,A)')  &
              ' Warning: lbfgs_bnd_maxmove =', rpath_info%lbfgs_bnd_maxmove, &
              ' is too small. Restoring default value'  
              rpath_info%lbfgs_bnd_maxmove = 0.1D+00
            end if
            write(MsgOut,'(A,F10.2)')       &
              '  lbfgs_bnd_maxmove = ', rpath_info%lbfgs_bnd_maxmove

          else
            write(MsgOut,'(A)')             &
              '  lbfgs_bnd         =         no'
            rpath_info%lbfgs_bnd_qmonly = .false.

          end if

        else if (rpath_info%method == MEPmethod_String) then
          write(MsgOut,'(A,F10.2)') &
            '  delta             = ', rpath_info%delta

        end if

        write(MsgOut,'(A,A)') &
          '  mepatm_select_index   = ', trim(rpath_info%mepatm_select_index)

      end if

      if (rpath_info%rpathmode == RpathmodeFEP) then
        write(MsgOut,'(A,I10,A,I10)') &
                              '  fep_period        = ', rpath_info%fep_period
        if (rpath_info%esp_md) then
          write(MsgOut,'(A)') '  esp_md            =        yes'
        else
          write(MsgOut,'(A)') '  esp_md            =         no'
        end if
        if (rpath_info%esp_energy) then
          write(MsgOut,'(A)') '  esp_energy        =        yes'
        else
          write(MsgOut,'(A)') '  esp_energy        =         no'
        end if

        write(MsgOut,'(A,A)') &
          '  mepatm_select_index = ', trim(rpath_info%mepatm_select_index)

      end if

      write(MsgOut,*)

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
  !! @authors      KY
  !! @param[in]    rpath_info : RPATH section control parameters information
  !! @param[in]    sel_info   : selector input information
  !! @param[in]    rst        : rst data
  !! @param[in]    rstmep     : rstmep data
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    constraints: information of constraints
  !! @param[in]    molecule   : molecule information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] minimize   : minimize information
  !! @param[inout] dynvars    : dynamic variables
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath(rpath_info, sel_info, rst, rstmep, dynamics, &
                         constraints, molecule, restraints, enefunc,  &
                         minimize, dynvars, rpath)

    ! formal arguments
    type(s_rpath_info),      intent(in)    :: rpath_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_rst),             intent(in)    :: rst
    type(s_rstmep),          intent(in)    :: rstmep
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_constraints),     intent(in)    :: constraints
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_minimize),        intent(inout) :: minimize
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_rpath),           intent(inout) :: rpath


    rpath%rpathmode       = rpath_info%rpathmode

    if (rpath_info%rpathmode == RpathmodeMFEP) then
      call setup_rpath_mfep(rpath_info, rst, dynamics, molecule, &
                      restraints, enefunc, rpath)

    else if (rpath_info%rpathmode == RpathmodeMEP .or.  &
             rpath_info%rpathmode == RpathmodeMEPMD) then

      ! Set rpath variables
      !
      rpath%ncycle          = rpath_info%ncycle
      rpath%nreplica        = rpath_info%nreplica
      rpath%delta           = rpath_info%delta
      rpath%fix_terminal    = rpath_info%fix_terminal
      rpath%rstout_period   = rpath_info%rstout_period
      rpath%eneout_period   = rpath_info%eneout_period
      rpath%crdout_period   = rpath_info%crdout_period
      rpath%mep_partial_opt = rpath_info%mep_partial_opt
      rpath%method          = rpath_info%method

      rpath%massWeightCoord = rpath_info%massWeightCoord
      rpath%tol_energy      = rpath_info%tol_energy
      rpath%tol_path        = rpath_info%tol_path

      rpath%tol_rmsg          = rpath_info%tol_rmsg
      rpath%tol_maxg          = rpath_info%tol_maxg
      rpath%k_spring          = rpath_info%k_spring
      rpath%climbing_image    = rpath_info%climbing_image
      rpath%tol_rmsg_cineb    = rpath_info%tol_rmsg_cineb
      rpath%ncorrection       = rpath_info%ncorrection
      rpath%lbfgs_bnd         = rpath_info%lbfgs_bnd
      rpath%lbfgs_bnd_qmonly  = rpath_info%lbfgs_bnd_qmonly
      rpath%lbfgs_bnd_maxmove = rpath_info%lbfgs_bnd_maxmove
      rpath%verbose           = rpath_info%verbose
      rpath%do_cineb          = .false.

      rpath%optimizer         = rpath_info%optimizer
      rpath%force_scale_init  = rpath_info%force_scale_init
      rpath%force_scale_max   = rpath_info%force_scale_max

      if (nrep_per_proc > 1) then
        rpath%rstout_period = 1
        rpath%crdout_period = 1
      end if

      if (rpath_info%rpathmode == RpathmodeMEP) then
        call setup_rpath_mep(rpath_info%mepatm_select_index, sel_info, rst,   &
                             molecule, enefunc%qmmm, minimize, dynvars, rpath)

      else if (rpath_info%rpathmode == RpathmodeMEPMD) then
        rpath%rpath_period = rpath_info%rpath_period
        call setup_rpath_mepmd(rpath_info%mepatm_select_index, sel_info, rst, &
                             dynamics, constraints, &
                             molecule, enefunc, dynvars, rpath)

      end if

    else if (rpath_info%rpathmode == RpathmodeFEP) then

      ! Set rpath variables
      rpath%nreplica      = rpath_info%nreplica
      rpath%esp_energy    = rpath_info%esp_energy
      rpath%esp_md        = rpath_info%esp_md

      call setup_rpath_fep(rpath_info%fep_period, &
                           rpath_info%mepatm_select_index, &
                           sel_info, rst, rstmep, molecule, &
                           enefunc%qmmm, dynamics, dynvars, rpath)

    end if

    return

  end subroutine setup_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_resetup
  !> @brief        setup MEP atoms for the case "nreplica > nproc"
  !! @authors      SI
  !! @param[in]    molecule   : molecule information
  !! @param[in]    qmmm       : QMMM information
  !! @param[out]   rpath      : rpath parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_resetup(inp_info, out_info, qmmm_info, rpath_info, &
                                 sel_info, dynamics, molecule, &
                                 enefunc, minimize, dynvars, output, &
                                 boundary, top, par, gpr, psf, &
                                 prmtop, grotop, rst, fit, rstmep, &
                                 eef1, spot, rpath)

    ! formal arguments
    type(s_inp_info),        intent(inout) :: inp_info
    type(s_out_info),        intent(inout) :: out_info
    type(s_qmmm_info),       intent(in)    :: qmmm_info
    type(s_rpath_info),      intent(in)    :: rpath_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_minimize),        intent(inout) :: minimize
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_output),          intent(inout) :: output
    type(s_boundary),        intent(inout) :: boundary

    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb)                            :: pdb
    type(s_crd)                            :: crd
    type(s_ambcrd)                         :: ambcrd
    type(s_grocrd)                         :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb)                            :: ref
    type(s_pdb),             intent(inout) :: fit
    type(s_ambcrd)                         :: ambref
    type(s_grocrd)                         :: groref
    type(s_mode)                           :: mode
    type(s_rstmep),          intent(inout) :: rstmep
    type(s_eef1),            intent(inout) :: eef1
    type(s_spot),            intent(inout) :: spot
    type(s_rpath),           intent(inout) :: rpath


    call input_rpath_resetup(inp_info, top, par, gpr, psf, prmtop, &
                             grotop, pdb, crd, ambcrd, grocrd, rst, ref, fit, &
                             ambref, groref, mode, rstmep, eef1, spot, rpath, &
                             out_info%rstfile)
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call setup_dynvars(molecule, rst, dynvars, dynamics)
    if (enefunc%qmmm%do_qmmm) then      
      call setup_qmmm_directory(qmmm_info, enefunc%qmmm)
    end if
    call setup_rpath_mep(rpath_info%mepatm_select_index, &
                         sel_info, rst, molecule, &
                         enefunc%qmmm, minimize, dynvars, rpath, rstmep)
    call setup_output_rpath(out_info, dynamics, rpath, output)

    return

  end subroutine setup_rpath_resetup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath
  !> @brief        run string method
  !! @authors      TM, YK, KY, SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] minimize   : minimize information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath(inp_info, out_info, qmmm_info, rpath_info, sel_info, &
                       output, molecule, enefunc, dynvars, minimize, &
                       dynamics, pairlist, boundary, constraints, ensemble, &
                       rpath)   

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule),         intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_minimize),         intent(inout) :: minimize
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_rpath),            intent(inout) :: rpath

    type(s_inp_info),         intent(inout) :: inp_info  
    type(s_out_info),         intent(inout) :: out_info  
    type(s_qmmm_info),        intent(in)    :: qmmm_info 
    type(s_rpath_info),       intent(in)    :: rpath_info
    type(s_sel_info),         intent(in)    :: sel_info  
    type(s_rst)                             :: rst       
    type(s_rstmep)                          :: rstmep     

    ! local variables   
    integer                                 :: i, j
    logical                                 :: conv = .false.      
    type(s_top)                             :: top
    type(s_par)                             :: par
    type(s_gpr)                             :: gpr
    type(s_psf)                             :: psf
    type(s_prmtop)                          :: prmtop
    type(s_grotop)                          :: grotop
    type(s_pdb)                             :: fit
    type(s_eef1)                            :: eef1
    type(s_spot)                            :: spot


    if (rpath%rpathmode == RpathmodeMFEP) then
      call run_rpath_mfep(output, molecule, enefunc, dynvars, dynamics, &
                    pairlist, boundary, constraints, ensemble, rpath)

    else if (rpath%rpathmode == RpathmodeMEPMD) then
      call run_rpath_mepmd(output, molecule, enefunc, dynvars, dynamics, &
                    pairlist, boundary, constraints, ensemble, rpath, conv)

    else if (rpath%rpathmode == RpathmodeMEP) then
      if (nrep_per_proc == 1) then
        call run_rpath_mep(output, molecule, enefunc, dynvars, minimize, &
             dynamics, pairlist, boundary, rpath, conv)
      else
        ! If nreplica is larger than total MPI proccess        
        !
        i = 0          
        do !! Cycle
          do j = 1, nrep_per_proc
            my_replica_no = nrep_per_proc*(my_country_no) + j

            ! Initialize MEP status again for new replica
            call setup_rpath_resetup(inp_info, out_info, qmmm_info, &
                                     rpath_info, sel_info, dynamics, &
                                     molecule, enefunc, minimize, dynvars, &
                                     output, boundary, top, par, &
                                     gpr, psf, prmtop, grotop, rst, fit, &
                                     rstmep, eef1, spot, rpath)

            ! Enter run_rpath_mep
            call run_rpath_mep(output, molecule, enefunc, dynvars, minimize, &
                               dynamics, pairlist, boundary, rpath, conv, i, j)

            rpath%first_replica = .false.
          end do

          rpath%first_replica  = .true.
          rpath%first_iter     = .false.

          ! If MEP reach limit
          !
          if (rpath%mep_exit) exit
          if (rpath%method == MEPmethod_NEB) then
            if (rpath%neb_cycle == rpath%ncycle + 2) exit
          else
            if (i == rpath%ncycle + 1) exit
          end if  
          i = i + 1
        end do
      end if

    end if

    return

  end subroutine run_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rpath_mfep
  !> @brief        setup RPATH
  !! @authors      TM, YK
  !! @param[in]    rpath_info : RPATH section control parameters information
  !! @param[in]    rst        : rst data
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] rpath      : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rpath_mfep(rpath_info, rst, dynamics, molecule, &
                              restraints, enefunc, rpath)

    ! formal arguments
    type(s_rpath_info),      intent(in)    :: rpath_info
    type(s_rst),     target, intent(in)    :: rst
    type(s_dynamics),target, intent(in)    :: dynamics
    type(s_molecule),target, intent(in)    :: molecule
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_rpath),   target, intent(inout) :: rpath

    ! local variables
    integer                    :: i, j, k, m, dimno, itmp, ikind, ierr
    integer                    :: max_nreplica
    integer                    :: natmgrp
    character(2000)            :: param
    real(wp),      allocatable :: ddata(:)
    integer,       pointer     :: grouplist(:,:), numatoms(:), atomlist(:,:)
    real(wp),      pointer     :: mass(:)
    integer                    :: iatm, iatm_xyz
    integer                    :: replicaid
    integer                    :: ifunc, igroup


    mass          => molecule%mass
    grouplist     => enefunc%restraint_grouplist
    numatoms      => enefunc%restraint_numatoms
    atomlist      => enefunc%restraint_atomlist

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
                         ' must be zero)')
        end if
      end if

    else if (rpath_info%rpath_period == 0) then
      rpath%equilibration_only = .true.
      rpath%rpath_period = dynamics%nsteps
      rpath%ncycle = 1
    else
      call error_msg('Setup_Rpath> rpath_period should be non-negative integer')
    end if

    ! setup basic parameters
    !
    rpath%nreplica     = rpath_info%nreplica
    max_nreplica       = rpath_info%nreplica
    rpath%delta        = rpath_info%delta
    rpath%smooth       = rpath_info%smooth
    rpath%fix_terminal = rpath_info%fix_terminal
    rpath%avoid_shrinkage = rpath_info%avoid_shrinkage
    rpath%use_restart  = rpath_info%use_restart

    ifunc = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    ikind  = enefunc%restraint_kind(ifunc)

    ! setup force constants and references
    !
    do i = 1, rpath_info%dimension
      itmp = rpath_info%rest_function(i)
      if (enefunc%restraint_kind(itmp) .ne. ikind) then
          call error_msg('Setup_Rpath> [ERROR] multiple types of '//&
          'restraint functions as CV are not allowed')
      end if
      enefunc%restraint_rpath_func(itmp)=i

    end do

    if (ikind .ne. RestraintsFuncPOSI .and. &
        ikind .ne. RestraintsFuncDIST .and. &
        ikind .ne. RestraintsFuncDIHED .and. &
        ikind .ne. RestraintsFuncPC) then
      call error_msg('Setup_Rpath> [ERROR] This CV is not allowed in atdyn')
    end if

    ! setup fitting index
    !
    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then
      enefunc%rpath_pos_func = ifunc
      if (enefunc%fitting_method .ne. FittingMethodNO) then
        if (.not. enefunc%do_fitting) &
          call error_msg('Setup_Rpath> [ERROR] fitfile should be specified')
        if (enefunc%mass_weight) &
          call error_msg('Setup_Rpath> [ERROR] mass_weight is not allowed')
      end if

      do i = 1, enefunc%num_restraintfuncs
        if (enefunc%restraint_kind(i) .eq. ikind) then
          if (enefunc%restraint_rpath_func(i) .eq. 0) then
            call error_msg('Setup_Rpath> [ERROR] all POSI restraint should'//&
               ' be used as R-path CV')
          end if
        end if
      end do

    end if

    enefunc%restraint_replica_index(:) = 0
    do i = 1, rpath_info%dimension
      ifunc = rpath_info%rest_function(i)
      enefunc%restraint_replica_index(ifunc) = i
    end do

    if (ikind == RestraintsFuncPOSI) then

      rpath%dimension = numatoms(igroup) * 3

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
          iatm = atomlist(int((i-1)/3)+1,igroup)
          iatm_xyz = mod(i-1, 3) + 1
          do j = 1, rpath%nreplica
            rpath%rest_constants(1:4, i, replicaid) = &
              enefunc%restraint_const_replica(j, ifunc)
            rpath%rest_reference(1:2, i, replicaid) = &
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
    call assign_condition(rpath, enefunc)

    ! determine the number of atoms involved in collective variable
    !
    ifunc = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    select case(enefunc%restraint_kind(ifunc))

    case (RestraintsFuncPOSI)

      natmgrp = numatoms(igroup)

    case (RestraintsFuncPC:RestraintsFuncPCCOM)

      natmgrp = numatoms(igroup)

    case default

      natmgrp = enefunc%max_restraint_numgrps

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
             enefunc%stats_mass(1:natmgrp, 1:rpath%dimension))

    enefunc%stats_delta(:)    = 0.0_wp
    enefunc%stats_grad(:,:,:) = 0.0_wp
    enefunc%stats_force(:)    = 0.0_wp
    enefunc%stats_metric(:,:) = 0.0_wp
    enefunc%stats_atom(:,:)   = 0
    enefunc%stats_mass(:,:)   = 0.0_wp

    ! setup atom index and mass
    !
    ifunc = rpath_info%rest_function(1)
    igroup = grouplist(1, ifunc)
    select case(enefunc%restraint_kind(ifunc))

    case (RestraintsFuncPOSI)

      do i = 1, natmgrp
        j = atomlist(i, 1)
        ! j = atomlist(i,igroup)
        enefunc%stats_atom(i,:) = j
        enefunc%stats_mass(i,:) = mass(j)
      end do

    case (RestraintsFuncPC:RestraintsFuncPCCOM)

      do i = 1, natmgrp
        j = atomlist(i, 1)
        ! j = atomlist(i,igroup)
        enefunc%stats_atom(i,:) = j
        enefunc%stats_mass(i,:) = mass(j)
      end do

    case default

      do dimno = 1, rpath%dimension
        !igroup = grouplist(1, dimno)
        do i = 1, natmgrp
          !j = atomlist(i,igroup)
          j = atomlist(1,grouplist(i,dimno))
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
            write(MsgOut,'(a,i0,a,$)')'  constant', i, ' = '
            do j = 1, rpath%nreplica
              write(MsgOut, '(F10.4,$)') rpath%rest_constants(1, i, j)
            end do
            write(MsgOut, *)
            write(MsgOut,'(a,i0,a,$)')'  reference', i, ' = '
            do j = 1, rpath%nreplica
              write(MsgOut, '(F10.4,$)') rpath%rest_reference(1, i, j)
            end do
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

  end subroutine setup_rpath_mfep

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

#ifdef HAVE_MPI_GENESIS
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
#endif

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
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition(rpath, enefunc)

    ! formal arguments
    type(s_rpath),           intent(in)    :: rpath
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i
    integer                  :: replicaid
    integer                  :: ifunc, igroup
    integer                  :: iatm, iatm_xyz


    replicaid = my_country_no + 1
    ifunc = rpath%rest_function(1)
    igroup = enefunc%restraint_grouplist(1, ifunc)
    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI) then

      do i = 1, rpath%dimension
        iatm = enefunc%restraint_atomlist(int((i-1)/3)+1,igroup)
        iatm_xyz = mod(i-1, 3) + 1
        enefunc%restraint_refcoord(iatm_xyz, iatm) = &
                rpath%rest_reference(1, i, replicaid)
      end do

    else

      do i = 1, rpath%dimension
        ifunc = rpath%rest_function(i)
        enefunc%restraint_const(1:4, ifunc) = &
          rpath%rest_constants(1:4, i, replicaid)
        enefunc%restraint_ref  (1:2, ifunc) = &
          rpath%rest_reference(1:2, i, replicaid)
      end do

    end if

    return

  end subroutine assign_condition

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rpath_mfep
  !> @brief        run string method
  !! @authors      TM, YK
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] rpath       : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rpath_mfep(output, molecule, enefunc, dynvars, dynamics,    &
                            pairlist, boundary, constraints, ensemble, rpath)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_rpath),            intent(inout) :: rpath

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
    call assign_condition(rpath, enefunc)

    ! update boundary and pairlist
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
          enefunc%cg_pairlistdist_126,                     &
          enefunc%cg_pairlistdist_PWMcos,                  &
          enefunc%cg_pairlistdist_DNAbp,                   &
          enefunc%cg_pairlistdist_exv,                     &
          boundary)
    else
      call update_boundary(enefunc%table%table,               &
          enefunc%pairlistdist,                               &
          boundary)
    end if
    if (real_calc) then
      pairlist%allocation = .true.
      call update_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                           dynvars%coord_pbc, pairlist)

    end if

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
        call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      end if

      ! path evolution
      !
      if (.not. rpath%equilibration_only) then
        call evolve(enefunc, rpath)
        call reparametrize(rpath)
        call assign_condition(rpath, enefunc)
      end if

      ! output dynamical variables
      !
      call output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                        boundary, rpath)
    end do
     
    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_rpath_mfep

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
    integer                  :: repid, repid_i, repid_j
    real(dp)                 :: dnorm
    real(dp),    allocatable :: tangent_vector(:)
    integer                  :: ifunc

    integer,         pointer :: count
    real(wp),        pointer :: force(:), metric(:,:)


    if (.not. replica_main_rank) return

    count  => enefunc%stats_count
    force  => enefunc%stats_force
    metric => enefunc%stats_metric

    repid = my_country_no + 1

    do dimno = 1, rpath%dimension
      force(dimno) = force(dimno) / real(count,wp)
    end do

    do dimno_i = 1, rpath%dimension
      do dimno_j = 1, rpath%dimension
        metric(dimno_i,dimno_j) = metric(dimno_i,dimno_j) / real(count,wp)
      end do
    end do

    force = matmul(metric, force)

    if (rpath%fix_terminal) then
      if ((repid == 1) .or. (repid == rpath%nreplica)) then
        force(:) = 0.0_wp
      end if
    else if (rpath%avoid_shrinkage) then
      repid_i = 0
      repid_j = 0
      if (repid == 1) then
        repid_i = 1
        repid_j = 2
      else if (repid == rpath%nreplica) then
        repid_i = rpath%nreplica
        repid_j = rpath%nreplica - 1
      end if
      if (((repid == 1) .or. (repid == rpath%nreplica))  &
          .and. (repid_i > 0) .and. (repid_j > 0)) then
        allocate(tangent_vector(rpath%dimension))
        do dimno = 1, rpath%dimension
          tangent_vector(dimno) = rpath%rest_reference(1,dimno, repid_j) &
                                 -rpath%rest_reference(1,dimno, repid_i)
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
    end if
    
    do dimno_i = 1, rpath%dimension
      rpath%rest_reference(1:2, dimno_i, repid) = &
        rpath%rest_reference(1:2, dimno_i, repid) + rpath%delta * force(dimno_i)
    end do

    count       = 0
    force(:)    = 0.0_wp
    metric(:,:) = 0.0_wp

    return

  end subroutine evolve

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
    integer                  :: i, j, k
    integer                  :: dimno, dimno_i, dimno_j
    integer                  :: repid, repid_i, repid_j
    real(dp)                 :: distance_prev, distance_init, dtmp

    real(dp),    allocatable :: path(:,:), path_smooth(:,:), path_reparm(:,:)
    real(dp),    allocatable :: path_leng(:), path_equi(:)

    real(dp),        pointer :: force(:), metric(:,:)
    real(dp),        pointer :: before_gather(:), after_gather(:)


    force         => rpath%force
    metric        => rpath%metric
    before_gather => rpath%before_gather
    after_gather  => rpath%after_gather

    repid = my_country_no + 1

    allocate(path(rpath%nreplica, rpath%dimension),&
             path_smooth(rpath%nreplica, rpath%dimension),&
             path_reparm(rpath%nreplica, rpath%dimension),&
             path_leng(rpath%nreplica), path_equi(rpath%nreplica))

    do dimno = 1, rpath%dimension
      before_gather(dimno) = rpath%rest_reference(1, dimno, repid)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_gather(before_gather, rpath%dimension, mpi_real8,&
                    after_gather,  rpath%dimension, mpi_real8,&
                    0, mpi_comm_airplane, ierror)

    if (main_rank) then
      do repid_i = 1, rpath%nreplica
        do dimno = 1, rpath%dimension
          path(repid_i,dimno) = after_gather((repid_i-1)*rpath%dimension+dimno)
        end do
      end do

      ! smooth
      path_smooth = path
      do repid_i = 2, rpath%nreplica - 1
        do dimno = 1, rpath%dimension
          path_smooth(repid_i,dimno) = &
               (1.0_dp-rpath%smooth)*path(repid_i,dimno) + &
               (rpath%smooth*0.5_dp)*(path(repid_i-1,dimno) + &
               path(repid_i+1,dimno))
        end do
      end do

      ! calc path length
      path_leng = 0.0_dp
      do repid_i = 2, rpath%nreplica
        path_leng(repid_i) = path_leng(repid_i-1) + &
          sqrt( sum( (path_smooth(repid_i,:)-path_smooth(repid_i-1,:))**2) )
        ! write(6,*) 'BEFORE inter_replica', &
        !             sqrt( sum( (path_smooth(repid_i,:)-path_smooth(repid_i-1,:))**2) )
      end do

      ! equi path length
      path_equi = 0.0_dp
      do repid_i = 2, rpath%nreplica
        path_equi(repid_i) = (repid_i-1)*path_leng(rpath%nreplica)/ &
             real((rpath%nreplica -1),wp)
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
    end if

    distance_prev = 0.0_dp
    distance_init = 0.0_dp
    do repid_i = 1, rpath%nreplica
      do dimno = 1, rpath%dimension
        dtmp = path_reparm(repid_i,dimno)
        distance_prev = distance_prev +   &
          (rpath%rest_reference_prev(dimno, repid_i) - dtmp)**2
        distance_init = distance_init +   &
          (rpath%rest_reference_init(dimno, repid_i) - dtmp)**2

        rpath%rest_reference_prev(dimno, repid_i) = dtmp
      end do
    end do
    distance_prev = sqrt(distance_prev)
    distance_init = sqrt(distance_init)
    rpath%sum_distance = real(path_leng(rpath%nreplica),wp)
    rpath%distance_prev    = real(distance_prev,wp)/ &
                             real((rpath%nreplica), wp)
    rpath%distance_init    = real(distance_init,wp)/ &
                             real((rpath%nreplica), wp)

    ! broadcast restraint reference
    call mpi_bcast(path_reparm, rpath%dimension*rpath%nreplica,&
                   mpi_real8, 0, mpi_comm_world, ierror)
#endif

    do dimno = 1, rpath%dimension
      rpath%rest_reference(1:2, dimno, repid) = path_reparm(repid, dimno)
    end do

    deallocate(path, path_smooth, path_reparm, path_leng, path_equi)

    return

  end subroutine reparametrize

end module at_rpath_mod
