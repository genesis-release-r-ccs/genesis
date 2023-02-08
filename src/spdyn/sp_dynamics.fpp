!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_dynamics_mod
!> @brief   molecular dynamics simulation
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Norio Takase (NT),
!!          Chigusa Kobayashi(CK), Hiraku Oshima (HO)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_dynamics_mod

  use sp_output_mod
  use sp_md_leapfrog_mod
  use sp_md_vverlet_mod
  use sp_md_respa_mod
  use sp_restraints_mod
  use sp_boundary_mod
  use sp_communicate_mod
  use sp_remd_str_mod
  use sp_output_str_mod
  use sp_parallel_io_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_dynamics_str_mod
  use sp_restraints_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_boundary_str_mod
  use sp_domain_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_control_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_str_mod
  use sp_alchemy_mod
  use sp_fep_energy_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_dyn_info
    integer          :: integrator       =  IntegratorVVER
    integer          :: nsteps           =    100
    real(wp)         :: timestep         =  0.001_wp
    integer          :: eneout_period    =     10
    integer          :: crdout_period    =      0
    integer          :: velout_period    =      0
    integer          :: rstout_period    =      0
    integer          :: stoptr_period    =     10
    integer          :: nbupdate_period  =     10
    integer          :: elec_long_period =      1
    integer          :: iseed            =     -1
    real(wp)         :: initial_time     =    0.0_wp

    ! respa
    logical          :: xi_respa         =  .false.
    logical          :: xo_respa         =  .true.
    integer          :: thermo_period    =      1
    integer          :: baro_period      =      1

    ! simulated annealing MD
    logical          :: annealing        =  .false.
    integer          :: anneal_period    =      0
    real(wp)         :: dtemperature     =    0.0_wp
    logical          :: verbose          =  .false.

    ! targeted MD
    logical          :: target_md        =  .false.
    logical          :: steered_md       =  .false.
    real(wp)         :: initial_rmsd     =  0.0_wp
    real(wp)         :: final_rmsd       =  0.0_wp

    ! hydorgen mass repartitioning
    integer          :: hmr_target       = HmrTargetAll
    logical          :: hydrogen_mr      = .false.
    real(wp)         :: hmr_ratio        = 3.0_wp
    real(wp)         :: hmr_ratio_xh1    = 0.0_wp
    real(wp)         :: hmr_ratio_xh2    = 0.0_wp
    real(wp)         :: hmr_ratio_xh3    = 0.0_wp
    real(wp)         :: hmr_ratio_ring   = 0.0_wp
  end type s_dyn_info

  ! subroutines
  public  :: show_ctrl_dynamics
  public  :: read_ctrl_dynamics
  public  :: setup_dynamics
  public  :: setup_dynamics_pio
  public  :: run_md
  ! FEP
  public  :: run_md_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_dynamics
  !> @brief        show DYNAMICS section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_dynamics(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = VVER      # [VVER,VRES]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of MD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# velout_period = 0         # velocities output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# stoptr_period = 10         # remove translational and rotational motions period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# iseed         = -1        # random number seed '
        write(MsgOut,'(A)') '# initial_time  = 0.0       # initial time (ps)'
        write(MsgOut,'(A)') '# annealing     = NO        # simulated annealing'
        write(MsgOut,'(A)') '# anneal_period = 0         # annealing period'
        write(MsgOut,'(A)') '# dtemperature  = 0.0       # temperature change at annealing (K)'
        write(MsgOut,'(A)') '# target_md     = no        # targeted MD simulation'
        write(MsgOut,'(A)') '# steered_md    = no        # targeted MD simulation'
        write(MsgOut,'(A)') '# hydrogen_mr   = no        # hydrogen mass repartitioning'
        write(MsgOut,'(A)') '# hmr_target    = all       # target of hydrogen mass repartitioning'
        write(MsgOut,'(A)') '# hmr_ratio     = 3.0       # hydrogen mass repartitioning ratio'
        write(MsgOut,'(A)') '# hmr_ratio_xh1 = 3.0       # hydrogen mass repartitioning ratio'
        write(MsgOut,'(A)') '# hmr_ratio_xh2 = 3.0       # hydrogen mass repartitioning ratio'
        write(MsgOut,'(A)') '# hmr_ratio_xh3 = 3.0       # hydrogen mass repartitioning ratio'
        write(MsgOut,'(A)') '# hmr_ratio_ring= 3.0       # hydrogen mass repartitioning ratio'
        write(MsgOut,'(A)') '# initial_rmsd  = 0.0       # initial TMD target RMSD'
        write(MsgOut,'(A)') '# final_rmsd    = 0.0       # final TMD target RMSD'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = VVER      # [VVER,VRES]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of MD steps'
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
  !! @param[out]   dyn_info : DYNAMICS section control parameters information
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
    call read_ctrlfile_integer(handle, Section, 'elec_long_period',&
                               dyn_info%elec_long_period)
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
    call read_ctrlfile_integer(handle, Section, 'thermostat_period', &
                               dyn_info%thermo_period)
    call read_ctrlfile_integer(handle, Section, 'barostat_period', &
                               dyn_info%baro_period)
    call read_ctrlfile_logical(handle, Section, 'verbose',        &
                               dyn_info%verbose)
    call read_ctrlfile_logical(handle, Section, 'target_md',     &
                               dyn_info%target_md)
    call read_ctrlfile_logical(handle, Section, 'steered_md',    &
                               dyn_info%steered_md)
    call read_ctrlfile_real   (handle, Section, 'initial_rmsd',  &
                               dyn_info%initial_rmsd)
    call read_ctrlfile_real   (handle, Section, 'final_rmsd',    &
                               dyn_info%final_rmsd)
    call read_ctrlfile_logical(handle, Section, 'hydrogen_mr',   &
                               dyn_info%hydrogen_mr)
    call read_ctrlfile_type   (handle, Section, 'hmr_target',    &
                               dyn_info%hmr_target, HmrTargetTypes)
    call read_ctrlfile_real   (handle, Section, 'hmr_ratio',     &
                               dyn_info%hmr_ratio)
    call read_ctrlfile_real   (handle, Section, 'hmr_ratio_xh1', &
                               dyn_info%hmr_ratio_xh1)
    call read_ctrlfile_real   (handle, Section, 'hmr_ratio_xh2', &
                               dyn_info%hmr_ratio_xh2)
    call read_ctrlfile_real   (handle, Section, 'hmr_ratio_xh3', &
                               dyn_info%hmr_ratio_xh3)
    call read_ctrlfile_real   (handle, Section, 'hmr_ratio_ring', &
                               dyn_info%hmr_ratio_ring)

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
      write(MsgOut,'(A20,I10)')                                  &
            '  iseed           = ', dyn_info%iseed
      if (dyn_info%hydrogen_mr) then
        write(MsgOut,'(A)') '  hydrogen_mr     =        yes'
        write(MsgOut,'(A20,A10)')                                &
              '  hmr_target      = ', trim(HmrTargetTypes(dyn_info%hmr_target))
        write(MsgOut,'(A20,F10.3)')                              &
          '  hmr_ratio       = ', dyn_info%hmr_ratio
        if (dyn_info%hmr_ratio_xh1 /= dyn_info%hmr_ratio .and.   &
            dyn_info%hmr_ratio_xh1 > EPS)                        &
        write(MsgOut,'(A20,F10.3)')                              &
          '  hmr_ratio_xh1   = ', dyn_info%hmr_ratio_xh1
        if (dyn_info%hmr_ratio_xh2 /= dyn_info%hmr_ratio .and.   &
            dyn_info%hmr_ratio_xh2 > EPS)                        &
        write(MsgOut,'(A20,F10.3)')                              &
          '  hmr_ratio_xh2   = ', dyn_info%hmr_ratio_xh2
        if (dyn_info%hmr_ratio_xh3 /= dyn_info%hmr_ratio .and.   &
            dyn_info%hmr_ratio_xh3 > EPS)                        &
        write(MsgOut,'(A20,F10.3)')                              &
          '  hmr_ratio_xh3   = ', dyn_info%hmr_ratio_xh3
        if (dyn_info%hmr_ratio_ring /= dyn_info%hmr_ratio .and.   &
            dyn_info%hmr_ratio_ring > EPS)                        &
        write(MsgOut,'(A20,F10.3)')                              &
          '  hmr_ratio_ring   = ', dyn_info%hmr_ratio_xh3
      end if

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

      ! respa
      !
      if (dyn_info%elec_long_period > 1) then
        write(MsgOut,'(A21,I8)')                   &
              '  elec_long_period  = ', dyn_info%elec_long_period
        write(MsgOut,'(A21,I8)')                   &
              '  thermostat_period = ', dyn_info%thermo_period
        write(MsgOut,'(A21,I8)')                   &
              '  barostat_period   = ', dyn_info%baro_period
      end if

      ! verbose
      if (dyn_info%verbose) then
        write(MsgOut,'(A)') '  verbose         =        yes'
      else
        write(MsgOut,'(A)') '  verbose         =         no'
      end if

      ! tareted md
      if (dyn_info%target_md) then
        write(MsgOut,'(A)') '  target_md       =        yes'
        write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
              '  initial rmsd    = ', dyn_info%initial_rmsd,       &
              '  final rmsd      = ', dyn_info%final_rmsd
      else if (dyn_info%steered_md) then
        write(MsgOut,'(A)') '  steered_md      =        yes'
        write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
              '  initial rmsd    = ', dyn_info%initial_rmsd,       &
              '  final rmsd      = ', dyn_info%final_rmsd
      else
        write(MsgOut,'(A)') '  target_md       =         no'
        write(MsgOut,'(A)') '  steered_md      =         no'
      end if

      ! adjust barostat_period
      !
      if (mod(dyn_info%baro_period,dyn_info%thermo_period) /= 0) then
        if (main_rank) &
          write(MsgOut,'(A)') '  mod(barostat_period,thermostat_period) is not ZERO'
        do while (.true.)
          dyn_info%baro_period = dyn_info%baro_period + 1
          if (mod(dyn_info%baro_period,dyn_info%thermo_period) == 0) exit
        end do
        if (main_rank) then
          write(MsgOut,'(A)') '  we increase the barostat_period'
          write(MsgOut,'(A20,I10)')  &
            '  barostat_period = ', dyn_info%baro_period
        end if
      end if

      write(MsgOut,'(A)') ' '

    end if

    ! warning with leapfor integrator
    !
    if (dyn_info%integrator == IntegratorLEAP) then
      call error_msg ( &
        'LEAP integrator is not available. '//&
        'Please run VVER integrator instead.')
    end if
      
    ! error check
    !
    if (dyn_info%crdout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%crdout_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in crdout_period'//&
          '  mod(nsteps, crdout_period) is not ZERO')
      end if
    end if

    if (dyn_info%velout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%velout_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in velout_period'//&
          '  mod(nsteps, velout_period) is not ZERO')
      end if
    end if

    if (dyn_info%eneout_period > 0 .and.                        &
        mod(dyn_info%nsteps, dyn_info%eneout_period) /= 0) then
      call error_msg( &
        'Read_Ctrl_Dynamics> Error in eneout_period'//&
        '  mod(nsteps, eneout_period) is not ZERO')
    end if

    if (dyn_info%rstout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%rstout_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in rstout_period'//&
          '  mod(nsteps, rstout_period) is not ZERO')
      end if
    end if

    if (dyn_info%stoptr_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%stoptr_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in stoptr_period'//&
          '  mod(nsteps, stoptr_period) is not ZERO')
      end if
    end if

    if (dyn_info%nbupdate_period > 0 .and.                      &
        mod(dyn_info%nsteps, dyn_info%nbupdate_period) /= 0) then
      call error_msg( &
        'Read_Ctrl_Dynamics> Error in nbupdate_period'//&
        '  mod(nsteps, nbupdate_period) is not ZERO')
    end if

    if (dyn_info%anneal_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%anneal_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in anneal_period'//&
          '  mod(nsteps, anneal_period) is not ZERO')
      end if
    end if

    if (dyn_info%annealing .and. dyn_info%anneal_period <= 0) then
      call error_msg( &
        'Read_Ctrl_Dynamics> Error in anneal_period'//&
        '  annealing = YES, but anneal_period is <= ZERO')
    end if

    if (dyn_info%crdout_period > 0) then
      if (mod(dyn_info%crdout_period,dyn_info%nbupdate_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in crdout_period or nbupdate_period'//&
          '  mod(crdout_period,nbupdate_period) is not ZERO')
      end if
    end if

    if (dyn_info%elec_long_period > 0) then
      if (mod(dyn_info%nbupdate_period,dyn_info%elec_long_period) /= 0) then
        call error_msg( &
          'Read_Ctrl_Dynamics> Error in elec_long_period or nbupdate_period'//&
          '  mod(nbupdate_period,elec_long_period) is not ZERO')
      end if
    end if

    if (dyn_info%elec_long_period > 0) then
      if (mod(dyn_info%thermo_period,dyn_info%elec_long_period) /= 0) then
        call error_msg( &
         'Read_Ctrl_Dynamics> Error in elec_long_period or thermostat_period'//&
         '  mod(thermostat_period,elec_long_period) is not ZERO')
      end if
    end if

    return

  end subroutine read_ctrl_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics
  !> @brief        setup dynamics information
  !> @authors      TM, HO
  !! @param[in]    dyn_info   : DYNAMICS  section control parameters information
  !! @param[in]    bound_info : BOUNDARY  section control parameters information
  !! @param[in]    res_info   : RESTRAINT section control parameters information
  !! @param[in]    alch_info  : ALCHEMY section control parameters information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] dynamics   : dynamics information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics(dyn_info, bound_info, res_info, alch_info, &
                            molecule, dynamics)

    ! formal arguments
    type(s_dyn_info),        intent(in)    :: dyn_info
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_res_info),        intent(in)    :: res_info
    type(s_alch_info),       intent(in)    :: alch_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_dynamics),        intent(inout) :: dynamics

    ! local variables
    integer                  :: i
    integer                  :: ierror, iseed, iseed_init_vel


    ! initialize variables in dynamics
    !
    call init_dynamics(dynamics)

    ! setup variables
    !
    dynamics%integrator        = dyn_info%integrator
    dynamics%nsteps            = dyn_info%nsteps
    dynamics%timestep          = dyn_info%timestep
    dynamics%crdout_period     = dyn_info%crdout_period
    dynamics%velout_period     = dyn_info%velout_period
    dynamics%eneout_period     = dyn_info%eneout_period
    dynamics%rstout_period     = dyn_info%rstout_period
    dynamics%stoptr_period     = dyn_info%stoptr_period
    dynamics%nbupdate_period   = dyn_info%nbupdate_period
    dynamics%elec_long_period  = dyn_info%elec_long_period
    dynamics%initial_time      = dyn_info%initial_time
    dynamics%xi_respa          = dyn_info%xi_respa
    dynamics%xo_respa          = .not.(dynamics%xi_respa)
    dynamics%verbose           = dyn_info%verbose
    dynamics%thermo_period     = dyn_info%thermo_period
    dynamics%baro_period       = dyn_info%baro_period  
    dynamics%istart_step       = 1
    dynamics%iend_step         = dynamics%nsteps
    dynamics%target_md         = dyn_info%target_md
    dynamics%steered_md        = dyn_info%steered_md
    dynamics%initial_rmsd      = dyn_info%initial_rmsd
    dynamics%final_rmsd        = dyn_info%final_rmsd
    dynamics%hydrogen_mr       = dyn_info%hydrogen_mr
    dynamics%hmr_target        = dyn_info%hmr_target
    dynamics%hmr_ratio         = dyn_info%hmr_ratio
    dynamics%hmr_ratio_xh1     = dyn_info%hmr_ratio_xh1
    dynamics%hmr_ratio_xh2     = dyn_info%hmr_ratio_xh2
    dynamics%hmr_ratio_xh3     = dyn_info%hmr_ratio_xh3
    dynamics%hmr_ratio_ring    = dyn_info%hmr_ratio_ring

    ! FEP
    dynamics%fepout_period     = alch_info%fepout_period
    dynamics%equilsteps        = alch_info%equilsteps
    dynamics%nsteps            = dynamics%nsteps + alch_info%equilsteps
    if (mod(alch_info%fepout_period,dynamics%eneout_period) /= 0) then
      call error_msg('Setup_Dynamics> mod(fepout_period,eneout_period)'//&
        ' must be zero)')
    end if

    iseed          = dyn_info%iseed
    iseed_init_vel = iseed
    if (dyn_info%iseed == -1) then
      dynamics%iseed_read=.true.
      if (main_rank) &
        call random_seed_initial_time(iseed)
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(iseed, 1, mpi_integer, 0, mpi_comm_country, ierror)
#endif
      iseed_init_vel = iseed
      iseed = iseed + my_country_rank
    end if
    dynamics%iseed = iseed + 1000 * my_country_no

    dynamics%iseed_init_velocity = iseed_init_vel + 1000 * my_country_no

    if (bound_info%type == BoundaryTypePBC) then

      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .false.

    end if

    do i = 1, res_info%nfunctions
      if (res_info%function(i) == RestraintsFuncPOSI) then
        dynamics%stop_com_translation = .false.
        dynamics%stop_com_rotation    = .false.
      end if
    end do


    ! update number of degrees of freedom
    !
    if (dynamics%stop_com_translation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtract 3 translational degrees of freedom'
        write(MsgOut,'(a)') ' '
      end if

      call update_num_deg_freedom('After removing translation', &
                                  -3, molecule%num_deg_freedom)

    end if

    if (dynamics%stop_com_rotation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtract 3 rotational degrees of freedom'
        write(MsgOut,'(a)') ' '
      end if

      call update_num_deg_freedom('After removing rotation', &
                                  -3, molecule%num_deg_freedom)

    end if

    ! simulated annealing MD
    !
    dynamics%annealing     = dyn_info%annealing
    dynamics%anneal_period = dyn_info%anneal_period
    dynamics%dtemperature  = dyn_info%dtemperature


    ! setup random system
    !
    call random_init(dynamics%iseed)

    return

  end subroutine setup_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics_pio
  !> @brief        setup dynamics information
  !> @authors      NT
  !! @param[in]    dyn_info   : DYNAMICS section control parameters information
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[in]    res_info   : RESTRAINT section control parameters information
  !! @param[inout] domain     : domain information
  !! @param[inout] dynamics   : dynamics information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics_pio(dyn_info, bound_info, res_info, domain, &
                                dynamics)

    ! formal arguments
    type(s_dyn_info),        intent(in)    :: dyn_info
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_res_info),        intent(in)    :: res_info
    type(s_domain),          intent(inout) :: domain
    type(s_dynamics),        intent(inout) :: dynamics

    ! local variables
    integer                  :: i
    integer                  :: iseed, iseed_init_vel


    ! setup variables
    !
    dynamics%integrator        = dyn_info%integrator
    dynamics%nsteps            = dyn_info%nsteps
    dynamics%timestep          = dyn_info%timestep
    dynamics%crdout_period     = dyn_info%crdout_period
    dynamics%velout_period     = dyn_info%velout_period
    dynamics%eneout_period     = dyn_info%eneout_period
    dynamics%rstout_period     = dyn_info%rstout_period
    dynamics%stoptr_period     = dyn_info%stoptr_period
    dynamics%nbupdate_period   = dyn_info%nbupdate_period
    dynamics%elec_long_period  = dyn_info%elec_long_period
    dynamics%initial_time      = dyn_info%initial_time
    dynamics%xi_respa          = dyn_info%xi_respa
    dynamics%xo_respa          = .not.(dynamics%xi_respa)
    dynamics%verbose           = dyn_info%verbose
    dynamics%thermo_period     = dyn_info%thermo_period
    dynamics%baro_period       = dyn_info%baro_period
    dynamics%istart_step       = 1
    dynamics%iend_step         = dynamics%nsteps
    dynamics%target_md         = dyn_info%target_md
    dynamics%steered_md        = dyn_info%steered_md
    dynamics%initial_rmsd      = dyn_info%initial_rmsd
    dynamics%final_rmsd        = dyn_info%final_rmsd
    dynamics%hydrogen_mr       = dyn_info%hydrogen_mr
    dynamics%hmr_target        = dyn_info%hmr_target
    dynamics%hmr_ratio         = dyn_info%hmr_ratio
    dynamics%hmr_ratio_xh1     = dyn_info%hmr_ratio_xh1
    dynamics%hmr_ratio_xh2     = dyn_info%hmr_ratio_xh2
    dynamics%hmr_ratio_xh3     = dyn_info%hmr_ratio_xh3
    dynamics%hmr_ratio_ring    = dyn_info%hmr_ratio_ring

    if (bound_info%type == BoundaryTypePBC) then

      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .false.

    end if

    do i = 1, res_info%nfunctions
      if (res_info%function(i) == RestraintsFuncPOSI) then
        dynamics%stop_com_translation = .false.
        dynamics%stop_com_rotation    = .false.
      end if
    end do

    ! update number of degrees of freedom
    !
    if (dynamics%stop_com_translation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtract 3 translational degrees of freedom'
        write(MsgOut,'(a)') ' '
      end if

      call update_num_deg_freedom('After removing translation', &
                                  -3, domain%num_deg_freedom)

    end if

    if (dynamics%stop_com_rotation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtract 3 rotational degrees of freedom'
        write(MsgOut,'(a)') ' '
      end if

      call update_num_deg_freedom('After removing rotation', &
                                  -3, domain%num_deg_freedom)

    end if

    ! simulated annealing MD
    !
    dynamics%annealing     = dyn_info%annealing
    dynamics%anneal_period = dyn_info%anneal_period
    dynamics%dtemperature  = dyn_info%dtemperature

    iseed          = dyn_info%iseed
    iseed_init_vel = iseed
    if (dyn_info%iseed == -1) then
      if (main_rank) &
        call random_seed_initial_time(iseed)
#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(iseed, 1, mpi_integer, 0, mpi_comm_country, ierror)
#endif
      iseed_init_vel = iseed
      iseed = iseed + my_country_rank
    end if
    dynamics%iseed = iseed + 1000 * my_country_no
    dynamics%iseed_init_velocity = iseed_init_vel + 1000 * my_country_no

    if (.not. dynamics%restart) then
      ! setup random system
      !
      call random_term
      call random_init(dynamics%iseed)
    end if

    ! simulated annealing MD
    !
    dynamics%annealing     = dyn_info%annealing
    dynamics%anneal_period = dyn_info%anneal_period
    dynamics%dtemperature  = dyn_info%dtemperature

    return

  end subroutine setup_dynamics_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_md
  !> @brief        perform MD simulation
  !! @authors      YS, CK
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !! @param[inout] constraints : bond constraints
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_md(output, domain, enefunc, dynvars, dynamics, &
                    pairlist, boundary, constraints, ensemble, comm, remd)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_comm),             intent(inout) :: comm
    type(s_remd),             intent(inout) :: remd


    if (dynamics%target_md) then
      enefunc%restraint_rmsd_target = .true.
    else
      enefunc%restraint_rmsd_target = .false.
    end if

    if (enefunc%gamd%update_period > 0) then
      enefunc%gamd%gamd_stat = .true.
    else
      enefunc%gamd%gamd_stat = .false.
    end if

    call open_output(output)

    select case (dynamics%integrator)

    case (IntegratorLEAP)

      call leapfrog_dynamics(output, domain, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble,  &
                             comm, remd)

    case (IntegratorVVER)

      call vverlet_dynamics (output, domain, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble,  &
                             comm, remd)

    case (IntegratorVRES)

      call vverlet_respa_dynamics (output, domain, enefunc, dynvars, dynamics, &
                             pairlist, boundary, constraints, ensemble, comm,  &
                             remd)

    end select

    call close_output(output)

    return

  end subroutine run_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_md_fep
  !> @brief        perform MD simulation for FEP
  !! @authors      HO
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !! @param[inout] constraints : bond constraints
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !! @param[inout] remd        : remd information
  !! @param[inout] alchemy     : alchemy information for FEP
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_md_fep(output, domain, enefunc, dynvars, dynamics, &
                        pairlist, boundary, constraints, ensemble, comm, remd, &
                        alchemy)
    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_comm),             intent(inout) :: comm
    type(s_remd),             intent(inout) :: remd
    type(s_alchemy),  target, intent(inout) :: alchemy

    ! local variables for FEP
    integer              :: i, num_windows
    integer, allocatable :: window_no(:)
    integer, pointer     :: lambid

    if (dynamics%target_md) then
      enefunc%restraint_rmsd_target = .true.
    else
      enefunc%restraint_rmsd_target = .false.
    end if

    if (enefunc%gamd%update_period > 0) then
      enefunc%gamd%gamd_stat = .true.
    else
      enefunc%gamd%gamd_stat = .false.
    end if

    ! FEP
    lambid => alchemy%lambid

    ! FEP: Set window indices and number of windows
    allocate(window_no(alchemy%num_fep_windows))
    if (enefunc%fep_direction == FEP_Bothsides) then
      num_windows = alchemy%num_fep_windows
      do i = 1, alchemy%num_fep_windows
        window_no(i) = i
      end do
    else if (enefunc%fep_direction == FEP_Forward) then
      num_windows = alchemy%num_fep_windows - 1
      do i = 1, alchemy%num_fep_windows
        window_no(i) = i
      end do
    else if (enefunc%fep_direction == FEP_Reverse) then
      num_windows = alchemy%num_fep_windows - 1
      do i = 1, alchemy%num_fep_windows
        window_no(i) = alchemy%num_fep_windows - i + 1
      end do
    end if

    call open_output(output)

    ! FEP: Loop changing lambda value
    do i = 1, num_windows

      ! Initialize MD time step and output title
      dynamics%initial_time = dynvars%time
      dynamics%istart_step  = 1 + (i-1)*dynamics%nsteps
      dynamics%iend_step    = i*dynamics%nsteps
      dynamics%equilsteps   = alchemy%equilsteps + (i-1)*dynamics%nsteps

      ! If single lambda calculation at state ref_lambid,
      ! other lambda values are skipped.
      if (alchemy%fep_md_type == FEP_Single) then
        if (i /= alchemy%ref_lambid) then
          cycle
        else
          dynamics%istart_step  = 1
          dynamics%iend_step    = dynamics%nsteps
          dynamics%equilsteps   = alchemy%equilsteps
        end if
      end if

      ! Set lambda values
      lambid = window_no(i)
      enefunc%lambljA   = alchemy%lambljA(lambid)
      enefunc%lambljB   = alchemy%lambljB(lambid)
      enefunc%lambelA   = alchemy%lambelA(lambid)
      enefunc%lambelB   = alchemy%lambelB(lambid)
      enefunc%lambbondA = alchemy%lambbondA(lambid)
      enefunc%lambbondB = alchemy%lambbondB(lambid)
      enefunc%lambrest  = alchemy%lambrest(lambid)

      ! REST2-like scaling in FEP
      call assign_lambda(alchemy, domain, enefunc)

      ! Output FEP window index
      if (main_rank) then
        write(MsgOut,'(a,i8)') "FEP window index ", lambid
        if (output%fepout) then
          write(output%fepunit,'(a,i8)') '# FEP window index ', lambid
        end if
      end if

      ! MD main loop
      select case (dynamics%integrator)

      case (IntegratorLEAP)

        call error_msg('LEAP integrator is not available in FEP')
!        call leapfrog_dynamics_fep(output, domain, enefunc, dynvars, dynamics, &
!                               pairlist, boundary, constraints, ensemble,  &
!                               comm, remd, alchemy)

      case (IntegratorVVER)

        call vverlet_dynamics (output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd, alchemy)

      case (IntegratorVRES)

        call vverlet_respa_dynamics(output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble, comm,  &
                               remd, alchemy)

      end select
      
    end do

    call close_output(output)

    return

  end subroutine run_md_fep

end module sp_dynamics_mod
