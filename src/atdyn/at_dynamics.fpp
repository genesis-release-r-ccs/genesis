!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_dynamics_mod
!> @brief   perform molecular dynamics simulation
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Norio Takase (NT),
!!          Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
!!          Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_dynamics_mod

  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_md_vverlet_cg_mod
  use at_output_mod
  use at_restraints_mod
  use at_boundary_mod
  use at_output_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_restraints_str_mod
  use at_energy_str_mod
  use at_energy_mod
  use at_qmmm_mod
  use random_mod
  use fileio_rst_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_dyn_info
    integer          :: integrator       =  IntegratorLEAP
    integer          :: nsteps           =    100
    real(wp)         :: timestep         =  0.001_wp
    integer          :: eneout_period    =     10
    integer          :: crdout_period    =      0
    integer          :: velout_period    =      0
    integer          :: rstout_period    =      0
    integer          :: stoptr_period    =     10
    integer          :: nbupdate_period  =     10
    integer          :: iseed            =     -1
    real(wp)         :: initial_time     =    0.0_wp
    integer          :: qmsave_period    =      0

    ! simulated annealing MD
    logical          :: annealing        =  .false.
    integer          :: anneal_period    =      0
    real(wp)         :: dtemperature     =    0.0_wp
    logical          :: verbose          = .false.

    logical          :: target_md        =  .false.
    logical          :: steered_md       =  .false.
    real(wp)         :: initial_value    =  0.0_wp
    real(wp)         :: final_value      =  0.0_wp

    ! box expansion & shrink MD
    logical          :: shrink_box       = .false.
    integer          :: shrink_period    =     0
    real(wp)         :: dbox_x           =   0.0_wp
    real(wp)         :: dbox_y           =   0.0_wp
    real(wp)         :: dbox_z           =   0.0_wp
    ! ESP/MM - MD
    logical          :: esp_mm           = .false.
    integer          :: calc_qm_period   = 0
    logical          :: avg_qm_charge    = .false.

  end type s_dyn_info

  ! subroutines
  public  :: show_ctrl_dynamics
  public  :: read_ctrl_dynamics
  public  :: setup_dynamics
  public  :: setup_dynamics_espmm
  public  :: run_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_dynamics
  !> @brief        show DYNAMICS section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_dynamics(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = LEAP      # [LEAP,VVER, VVER_CG]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of MD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# qmsave_period = 0         # Interval for saving QM files'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# velout_period = 0         # velocities output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# stoptr_period = 10        # remove translational and rotational motions period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# iseed         = -1        # random number seed '
        write(MsgOut,'(A)') '# initial_time  = 0.0       # initial time (ps)'
        write(MsgOut,'(A)') '# annealing     = NO        # simulated annealing'
        write(MsgOut,'(A)') '# anneal_period = 0         # annealing period'
        write(MsgOut,'(A)') '# dtemperature  = 0.0       # temperature change at annealing (K)'
        write(MsgOut,'(A)') '# verbose       = NO        # output verbosly'
        write(MsgOut,'(A)') '# target_md     = no        # targeted MD simulation'
        write(MsgOut,'(A)') '# steered_md    = no        # targeted MD simulation'
        write(MsgOut,'(A)') '# initial_value  = 0.0      # initial TMD target RMSD'
        write(MsgOut,'(A)') '# final_value    = 0.0      # final TMD target RMSD'

        write(MsgOut,'(A)') '# esp_mm         = NO       # ESP/MM MD'
        write(MsgOut,'(A)') '# calc_qm_period = 0        # QM calculation period'
        write(MsgOut,'(A)') '# avg_qm_charge  = YES      # Average QM charges'
        write(MsgOut,'(A)') ' '


      case ('remd', 'rpath')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = LEAP      # [LEAP,VVER]'
        write(MsgOut,'(A)') 'nsteps        = 1000      # number of MD steps in REMD'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# velout_period = 0         # velocities output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# stoptr_period = 10        # remove translational and rotational motions period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# iseed         = -1        # random number seed '
        write(MsgOut,'(A)') '# initial_time  = 0.0       # initial time (ps)'

        write(MsgOut,'(A)') '# esp_mm         = NO       # ESP/MM MD'
        write(MsgOut,'(A)') '# calc_qm_period = 0        # QM calculation period'
        write(MsgOut,'(A)') '# avg_qm_charge  = YES      # Average QM charges'
        write(MsgOut,'(A)') ' '

      case ('bd')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = BDEM      # [BDEM,BD2N,SDMP]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of BD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# velout_period = 0         # velocities output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# stoptr_period = 0         # remove translational and rotational motions period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# iseed         = -1        # random number seed '
        write(MsgOut,'(A)') '# initial_time  = 0.0       # initial time (ps)'
        write(MsgOut,'(A)') '# verbose       = NO        # output verbosly'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'remd', 'rpath')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = LEAP      # [LEAP,VVER]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of MD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') ' '

      case ('bd')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = BDEM      # [BDEM,BD2N,SDMP]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of BD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') ' '

      end select

    end if


    return

  end subroutine show_ctrl_dynamics
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_dynamics
  !> @brief        read DYNAMICS section in the control file
  !! @authors      TM, JJ
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   dyn_info : DYNAMICS section in control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_dynamics(handle, dyn_info)

    ! parameters
    character(*),            parameter     :: Section = 'Dynamics'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_dyn_info),        intent(inout) :: dyn_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'integrator',    &
                               dyn_info%integrator, IntegratorTypes)
    call read_ctrlfile_integer(handle, Section, 'nsteps',        &
                               dyn_info%nsteps)
    call read_ctrlfile_real   (handle, Section, 'timestep',      &
                               dyn_info%timestep)
    call read_ctrlfile_integer(handle, Section, 'eneout_period', &
                               dyn_info%eneout_period)
    call read_ctrlfile_integer(handle, Section, 'crdout_period', &
                               dyn_info%crdout_period)
    call read_ctrlfile_integer(handle, Section, 'velout_period', &
                               dyn_info%velout_period)
    call read_ctrlfile_integer(handle, Section, 'rstout_period', &
                               dyn_info%rstout_period)
    call read_ctrlfile_integer(handle, Section, 'stoptr_period', &
                               dyn_info%stoptr_period)
    call read_ctrlfile_integer(handle, Section, 'nbupdate_period', &
                               dyn_info%nbupdate_period)
    call read_ctrlfile_integer(handle, Section, 'qmsave_period',  &
                               dyn_info%qmsave_period)
    call read_ctrlfile_integer(handle, Section, 'iseed',         &
                               dyn_info%iseed)
    call read_ctrlfile_real   (handle, Section, 'initial_time',  &
                               dyn_info%initial_time)
    call read_ctrlfile_logical(handle, Section, 'annealing',     &
                               dyn_info%annealing)
    call read_ctrlfile_integer(handle, Section, 'anneal_period', &
                               dyn_info%anneal_period)
    call read_ctrlfile_real   (handle, Section, 'dtemperature',  &
                               dyn_info%dtemperature)
    call read_ctrlfile_logical(handle, Section, 'verbose',        &
                               dyn_info%verbose)
    call read_ctrlfile_logical(handle, Section, 'target_md',     &
                               dyn_info%target_md)
    call read_ctrlfile_logical(handle, Section, 'steered_md',    &
                               dyn_info%steered_md)
    call read_ctrlfile_real   (handle, Section, 'initial_value',  &
                               dyn_info%initial_value)
    call read_ctrlfile_real   (handle, Section, 'final_value',    &
                               dyn_info%final_value)
    call read_ctrlfile_logical(handle, Section, 'esp_mm',         &
                               dyn_info%esp_mm)
    call read_ctrlfile_integer(handle, Section, 'calc_qm_period', &
                               dyn_info%calc_qm_period)
    call read_ctrlfile_logical(handle, Section, 'avg_qm_charge',  &
                               dyn_info%avg_qm_charge)

    call read_ctrlfile_logical(handle, Section, 'shrink_box',     &
                               dyn_info%shrink_box)
    call read_ctrlfile_integer(handle, Section, 'shrink_period',  &
                               dyn_info%shrink_period)
    call read_ctrlfile_real   (handle, Section, 'dbox_x',  &
                               dyn_info%dbox_x)
    call read_ctrlfile_real   (handle, Section, 'dbox_y',  &
                               dyn_info%dbox_y)
    call read_ctrlfile_real   (handle, Section, 'dbox_z',  &
                               dyn_info%dbox_z)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Dynamics> Parameters of MD simulation'

      write(MsgOut,'(A20,A10,A20,I10)')                          &
            '  integrator      = ', IntegratorTypes(dyn_info%integrator), &
            '  nsteps          = ', dyn_info%nsteps
      write(MsgOut,'(A20,F10.4,A20,F10.4)')                      &
            '  timestep        = ', dyn_info%timestep,           &
            '  initial_time    = ', dyn_info%initial_time
      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  eneout_period   = ', dyn_info%eneout_period,      &
            '  rstout_period   = ', dyn_info%rstout_period
      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  crdout_period   = ', dyn_info%crdout_period,      &
            '  velout_period   = ', dyn_info%velout_period
      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  nbupdate_period = ', dyn_info%nbupdate_period,    &
            '  stoptr_period   = ', dyn_info%stoptr_period
      write(MsgOut,'(A20,I10)')                                    &
            '  qmsave_period   = ', dyn_info%qmsave_period
      write(MsgOut,'(A20,I10)')                                  &
            '  iseed           = ', dyn_info%iseed

      ! simulated annealing
      !
      if (dyn_info%annealing) then
        write(MsgOut,'(A)') '  annealing       =        yes'
        write(MsgOut,'(A20,I10,A20,F10.3)')                      &
              '  anneal_period   = ', dyn_info%anneal_period,    &
              '  dtemperature    = ', dyn_info%dtemperature
      else
        write(MsgOut,'(A)') '  annealing       =         no'
      end if

      ! verbose
      if (dyn_info%verbose) then
        write(MsgOut,'(A)') '  verbose         =        yes'
      else
        write(MsgOut,'(A)') '  verbose         =         no'
      end if

      if (dyn_info%target_md) then
        write(MsgOut,'(A)') '  target_md       =        yes'
        write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
              '  initial value   = ', dyn_info%initial_value,       &
              '  final value     = ', dyn_info%final_value
      else
        write(MsgOut,'(A)') '  target_md       =         no'
      end if
      if (dyn_info%steered_md) then
        write(MsgOut,'(A)') '  steered_md      =        yes'
        write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
              '  initial value   = ', dyn_info%initial_value,       &
              '  final value     = ', dyn_info%final_value
      else
        write(MsgOut,'(A)') '  steered_md      =         no'
      end if

      if (dyn_info%shrink_box) then
        write(MsgOut,'(A)') '  shrink_box      =        yes'
        write(MsgOut,'(A,I10)') &
              '  shrink_period   = ', dyn_info%shrink_period
        write(MsgOut,'(A,F10.3)') &
              '  dbox_x          = ', dyn_info%dbox_x
        write(MsgOut,'(A,F10.3)') &
              '  dbox_y          = ', dyn_info%dbox_y
        write(MsgOut,'(A,F10.3)') &
              '  dbox_z          = ', dyn_info%dbox_z
      else
        write(MsgOut,'(A)') '  shrink_box      =         no'
      end if

      if (dyn_info%esp_mm) then
        write(MsgOut,'(A)') '  esp_mm          =        yes'
        write(MsgOut,'(A20,I10)')                                  &
              '  calc_qm_period  = ', dyn_info%calc_qm_period
        if (dyn_info%avg_qm_charge) then
          write(MsgOut,'(A)') '  avg_qm_charge   =        yes'
        else
          write(MsgOut,'(A)') '  avg_qm_charge   =         no'
        end if
      else
        write(MsgOut,'(A)') '  esp_mm          =         no'
      end if

      write(MsgOut,'(A)') ' '
    end if


    ! error check
    !
    if (main_rank) then
      if (dyn_info%crdout_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%crdout_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, crdout_period) is not ZERO')
        end if
      end if

      if (dyn_info%velout_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%velout_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, velout_period) is not ZERO')
        end if
      end if

      if (dyn_info%eneout_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%eneout_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, eneout_period) is not ZERO')
        end if
      end if

      if (dyn_info%rstout_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%rstout_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, rstout_period) is not ZERO')
        end if
      end if

      if (dyn_info%stoptr_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%stoptr_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, stoptr_period) is not ZERO')
        end if
      end if

      if (dyn_info%nbupdate_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%nbupdate_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, nbupdate_period) is not ZERO')
        end if
      end if

      if (dyn_info%anneal_period > 0) then
        if (mod(dyn_info%nsteps, dyn_info%anneal_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(nsteps, anneal_period) is not ZERO')
        end if
      end if

      if (dyn_info%annealing .and. dyn_info%anneal_period <= 0) then
        call error_msg('Read_Ctrl_Dynamics> '//                   &
            'annealing = YES, but anneal_period is <= ZERO')
      end if

      if (dyn_info%crdout_period > 0) then
        if (mod(dyn_info%crdout_period,dyn_info%nbupdate_period) /= 0) then
          call error_msg('Read_Ctrl_Dynamics> '//                   &
              'mod(crdout_period,nbupdate_period) is not ZERO')
        end if
      end if
      if (dyn_info%shrink_box) then
        if (dyn_info%shrink_period > 0) then
          if (mod(dyn_info%shrink_period, dyn_info%nbupdate_period) /=0) then
            call error_msg('Read_Ctrl_Dynamics> '//                    &
                           'mod(shrink_period, nbupdate_period) is not ZERO')
          end if
        end if

        if (dyn_info%annealing) then
          call error_msg('Read_Ctrl_Dynamics>'//                    &
                         'annealing and shrink_md are not set at the same time')
        end if

        if (dyn_info%target_md) then
          call error_msg('Read_Ctrl_Dynamics>'//                    &
                         'target_md and shrink_md are not set at the same time')
        end if

        if (dyn_info%steered_md) then
          call error_msg('Read_Ctrl_Dynamics>'//                    &
                         'steered_md and shrink_md are not set at the same time')
        end if
      end if


      if (dyn_info%steered_md .and. dyn_info%target_md) then
        call error_msg('Read_Ctrl_Dynamics> steered_md and target_md are not set at the same time')
      endif

    end if

    return

  end subroutine read_ctrl_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics
  !> @brief        setup dynamics information
  !> @authors      TM
  !! @param[in]    dyn_info   : DYNAMICS  section in control parameters
  !! @param[in]    bound_info : BOUNDARY  section in control parameters
  !! @param[in]    res_info   : RESTRAINT section in control parameters
  !! @param[in]    out_info   : OUTPUT    section in control parameters
  !! @param[in]    rst        : restart information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] dynamics   : dynamics information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics(dyn_info, bound_info, res_info, out_info, &
                            rst, molecule, dynamics)

    ! formal arguments
    type(s_dyn_info),        intent(in)    :: dyn_info
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_res_info),        intent(in)    :: res_info
    type(s_out_info),        intent(in)    :: out_info
    type(s_rst),             intent(in)    :: rst
    type(s_molecule),        intent(inout) :: molecule
    type(s_dynamics),        intent(inout) :: dynamics

    ! local variables
    integer                  :: i
    integer                  :: ierror, iseed


    ! initialize variables in dynamics
    !
    call init_dynamics(dynamics)

    ! setup dynamics
    !
    dynamics%integrator       = dyn_info%integrator
    dynamics%nsteps           = dyn_info%nsteps
    dynamics%timestep         = dyn_info%timestep
    dynamics%crdout_period    = dyn_info%crdout_period
    dynamics%velout_period    = dyn_info%velout_period
    dynamics%eneout_period    = dyn_info%eneout_period
    dynamics%rstout_period    = dyn_info%rstout_period
    dynamics%stoptr_period    = dyn_info%stoptr_period
    dynamics%qmsave_period    = dyn_info%qmsave_period
    dynamics%nbupdate_period  = dyn_info%nbupdate_period
    dynamics%initial_time     = dyn_info%initial_time
    dynamics%annealing        = dyn_info%annealing
    dynamics%anneal_period    = dyn_info%anneal_period
    dynamics%dtemperature     = dyn_info%dtemperature
    dynamics%istart_step      = 1
    dynamics%iend_step        = dynamics%nsteps
    dynamics%verbose          = dyn_info%verbose
    dynamics%target_md        = dyn_info%target_md
    dynamics%steered_md       = dyn_info%steered_md
    dynamics%initial_value    = dyn_info%initial_value
    dynamics%final_value      = dyn_info%final_value
    dynamics%shrink_box       = dyn_info%shrink_box
    dynamics%shrink_period    = dyn_info%shrink_period
    dynamics%dbox_x           = dyn_info%dbox_x
    dynamics%dbox_y           = dyn_info%dbox_y
    dynamics%dbox_z           = dyn_info%dbox_z
    dynamics%esp_mm           = dyn_info%esp_mm
    dynamics%calc_qm_period   = dyn_info%calc_qm_period
    dynamics%avg_qm_charge    = dyn_info%avg_qm_charge
    
    iseed = dyn_info%iseed
    if (rst%rstfile_type == RstfileTypeUndef .or. &
        rst%rstfile_type == RstfileTypeMin ) then
      dynamics%restart        = .false.
      dynamics%multiple_random_seed = .true.
      dynamics%random_restart = .false.
      if (dyn_info%iseed == -1) then
        if (main_rank) &
          call random_seed_initial_time(iseed)
#ifdef HAVE_MPI_GENESIS
        call mpi_bcast(iseed, 1, mpi_integer, 0, mpi_comm_country, ierror)
#endif
      endif
      dynamics%iseed          = iseed + 1000 * my_country_no
    else
      dynamics%restart        = .true.
      dynamics%multiple_random_seed = .false.
      if (dyn_info%iseed == -1) then
        dynamics%iseed        = rst%iseed
      else
        dynamics%iseed        = dyn_info%iseed + 1000 * my_country_no
        dynamics%random_restart  &
                              = .false.
      end if
    end if

    if (bound_info%type == BoundaryTypeNOBC) then
      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .true.

    else if (bound_info%type == BoundaryTypePBC) then
      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .false.

    end if

    do i = 1, res_info%nfunctions
      if (res_info%function(i) == RestraintsFuncPOSI .or. &
          res_info%function(i) == RestraintsFuncEM) then
        dynamics%stop_com_translation = .false.
        dynamics%stop_com_rotation    = .false.
      end if
    end do


    do i = 1, molecule%num_atoms
      if (molecule%mass(i) .lt. EPS) &
        call error_msg('Setup_Dynamics> Mass = 0 is not allowed')
    end do

    if (out_info%dcdfile /= '' .and. dynamics%crdout_period <= 0)    &
      dynamics%crdout_period = dynamics%nsteps
    if (out_info%dcdvelfile /= '' .and. dynamics%velout_period <= 0) &
      dynamics%velout_period = dynamics%nsteps
    if (out_info%rstfile /= '' .and. dynamics%rstout_period <= 0)    &
      dynamics%rstout_period = dynamics%nsteps

    ! update number of degrees of freedom
    !
    if (dynamics%stop_com_translation) then
      call update_num_deg_freedom('After removing translation', &
                                  -3, molecule%num_deg_freedom)
    end if

    if (dynamics%stop_com_rotation) then
      call update_num_deg_freedom('After removing rotation', &
                                  -3, molecule%num_deg_freedom)
    end if

    return

  end subroutine setup_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics_espmm
  !> @brief        setup ESP/MM-MD
  !> @authors      KY
  !! @param[inout] enefunc   : potential energy functions
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    pairlist  : information of nonbonded pair list
  !! @param[in]    dynanmics : dynamics information
  !! @param[inout] dynvars   : information of dynamic variables
  !! @param[inout] qmmm      : qmmm information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics_espmm(enefunc, molecule, pairlist, dynamics, dynvars, qmmm)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_molecule),        intent(inout) :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    type(s_energy)   :: dummy_ene
    real(wp)         :: dummy_force(1,1)


    if (.not. qmmm%do_qmmm) then
        call error_msg('Setup_Dynamics_ESPMM> QMMM must be set '// &
                       'to perform ESP/MM-MD')
    end if

    ! QMMM setings
    qmmm%qm_classical = .true.
    qmmm%qm_get_esp   = .true.
    qmmm%ene_only     = .true.

    ! check for ESP charge
    if (.not. qmmm%is_qm_charge) then

      if(main_rank) &
        write(MsgOut,'(a,/)') &
               "Setup_Dynamics_ESPMM> QM charges are not found. "// &
               "Now, computing QM charges (this may take some time) ..."

      qmmm%qm_count = 0
      qmmm%qm_classical = .false.
      call compute_energy_qmmm(enefunc, molecule, pairlist, dynvars%coord, qmmm, &
                               dummy_ene, dummy_force)
      qmmm%qm_classical = .true.

    else
      if(main_rank) &
        write(MsgOut,'("Setup_Dynamics_ESPMM> QM charge from rstfile",/)') 

    end if
 
    return

  end subroutine setup_dynamics_espmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_md
  !> @brief        perform MD simulation
  !! @authors      YS
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !! @param[inout] constraints : bond constraints
  !! @param[inout] ensemble    : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_md(output, molecule, enefunc, dynvars, dynamics, &
                    pairlist, boundary, constraints, ensemble)

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

    ! GaMD
    !
    if (enefunc%gamd%update_period > 0) then
      enefunc%gamd%gamd_stat = .true.
    else
      enefunc%gamd%gamd_stat = .false.
    end if

    ! Open output files
    !
    call open_output(output)

    ! MD main loop
    !
    select case (dynamics%integrator)

    case (IntegratorLEAP)

      call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble)

    case (IntegratorVVER)

      call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble)

    case (IntegratorVVER_CG)

      call vverlet_dynamics_cg (output, molecule, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble)

    end select

    ! Close output files
    !
    call close_output(output)


    return

  end subroutine run_md

end module at_dynamics_mod
