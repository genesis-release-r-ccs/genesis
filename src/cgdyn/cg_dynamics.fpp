!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_dynamics_mod
!> @brief   Performs molecular dynamics simulations.
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Norio Takase (NT),
!!          Chigusa Kobayashi(CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_dynamics_mod

  use cg_output_mod
  use cg_md_leapfrog_mod
  use cg_md_vverlet_mod
  use cg_restraints_mod
  use cg_boundary_mod
  use cg_communicate_str_mod
  use cg_output_str_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_dynamics_str_mod
  use cg_restraints_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_boundary_str_mod
  use cg_domain_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use fileio_control_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
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
    integer          :: lbupdate_period  =      0
    integer          :: elec_long_period =      1
    integer          :: iseed            =     -1
    real(wp)         :: initial_time     =    0.0_wp

    logical          :: gen_velocity     =  .true.
    ! respa
    integer          :: thermo_period    =      1
    integer          :: baro_period      =      1

    ! simulated annealing MD
    logical          :: annealing        =  .false.
    integer          :: anneal_period    =      0
    real(wp)         :: dtemperature     =    0.0_wp
    logical          :: verbose          =  .false.

    ! targeted MD
    logical          :: target_md        =  .false.
    logical          :: steered_md        =  .false.
    real(wp)         :: initial_rmsd     =  0.0_wp
    real(wp)         :: final_rmsd       =  0.0_wp
  end type s_dyn_info

  ! subroutines
  public  :: show_ctrl_dynamics
  public  :: read_ctrl_dynamics
  public  :: setup_dynamics
  public  :: run_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_dynamics
  !> @brief        Displays DYNAMICS section usage.
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

      case ('md')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = LEAP      # [LEAP,VVER]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of MD steps'
        write(MsgOut,'(A)') 'timestep      = 0.001     # timestep (ps)'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# velout_period = 0         # velocities output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# stoptr_period = 10         # Period for removing translational and &'
        write(MsgOut,'(A)') '#                                     # rotational motions.'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# lbupdate_period = 10      # load balance update period'
        write(MsgOut,'(A)') '# iseed         = -1        # random number seed'
        write(MsgOut,'(A)') '# initial_time  = 0.0       # initial time (ps)'
        write(MsgOut,'(A)') '# annealing     = NO        # simulated annealing'
        write(MsgOut,'(A)') '# anneal_period = 0         # annealing period'
        write(MsgOut,'(A)') '# dtemperature  = 0.0       # temperature change during annealing (K)'
        write(MsgOut,'(A)') '# target_md     = no        # targeted MD simulation'
        write(MsgOut,'(A)') '# steered_md    = no        # steered MD simulation'
        write(MsgOut,'(A)') '# initial_rmsd  = 0.0       # initial TMD target RMSD'
        write(MsgOut,'(A)') '# final_rmsd    = 0.0       # final TMD target RMSD'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[DYNAMICS]'
        write(MsgOut,'(A)') 'integrator    = LEAP      # [LEAP,VVER]'
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
  !> @brief        Reads DYNAMICS section in the control file.
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
    call read_ctrlfile_integer(handle, Section, 'lbupdate_period', &
                               dyn_info%lbupdate_period)
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
    call read_ctrlfile_logical(handle, Section, 'gen_velocity',   &
                               dyn_info%gen_velocity)
    call read_ctrlfile_logical(handle, Section, 'target_md',     &
                               dyn_info%target_md)
    call read_ctrlfile_logical(handle, Section, 'steered_md',    &
                               dyn_info%steered_md)
    call read_ctrlfile_real   (handle, Section, 'initial_rmsd',  &
                               dyn_info%initial_rmsd)
    call read_ctrlfile_real   (handle, Section, 'final_rmsd',    &
                               dyn_info%final_rmsd)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Dynamics> Parameters for MD simulation'

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
            '  lbupdate_period = ', dyn_info%lbupdate_period
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
      if (dyn_info%gen_velocity) write(MsgOut,'(A)') '  gen_velocity    =        yes' 

      ! respa
      !
      if (dyn_info%elec_long_period > 1) then
        write(MsgOut,'(A21,I10,A21,I10,A21,I10)')                   &
              ' elec_long_period  = ', dyn_info%elec_long_period,   &
              ' thermostat_period = ', dyn_info%thermo_period,      &
              ' barostat_period   = ', dyn_info%baro_period
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

      write(MsgOut,'(A)') ' '

    end if


    ! error check
    !
    if (dyn_info%crdout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%crdout_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For crdout_period, '//&
          'nsteps must be a multiple of crdout_period.')
    end if

    if (dyn_info%velout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%velout_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For velout_period, '//&
          'nsteps must be a multiple of velout_period.')
    end if

    if (dyn_info%eneout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%eneout_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For eneout_period, '//&
          'nsteps must be a multiple of eneout_period.')
    end if

    if (dyn_info%rstout_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%rstout_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For rstout_period, '//&
          'nsteps must be a multiple of rstout_period.')
    end if

    if (dyn_info%stoptr_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%stoptr_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For stoptr_period, '//&
          'nsteps must be a multiple of stoptr_period.')
    end if

    if (dyn_info%nbupdate_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%nbupdate_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For nbupdate_period, '//&
          'nsteps must be a multiple of nbupdate_period.')
    end if

    if (dyn_info%lbupdate_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%lbupdate_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For lbupdate_period, '//&
          'nsteps must be a multiple of lbupdate_period.')
    end if

    if (dyn_info%anneal_period > 0) then
      if (mod(dyn_info%nsteps, dyn_info%anneal_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For anneal_period, '//&
          'nsteps must be a multiple of anneal_period.')
    end if

    if (dyn_info%annealing .and. dyn_info%anneal_period <= 0) then
      call error_msg( &
        'Read_Ctrl_Dynamics> Error: Annealing is enabled, but '//&
        'anneal_period is <= 0.')
    end if

    if (dyn_info%annealing .and. dyn_info%integrator /= IntegratorLEAP) then
      call error_msg( &
        'Read_Ctrl_Dynamics> Error: Annealing is only allowed '//&
        'with the leapfrog integrator in cgdyn.')
    end if

    if (dyn_info%crdout_period > 0) then
      if (mod(dyn_info%crdout_period,dyn_info%nbupdate_period) /= 0) &
        call error_msg( &
          'Read_Ctrl_Dynamics> Error: For crdout_period or nbupdate_period, '//&
          'crdout_period must be a multiple of nbupdate_period.')
    end if

    if (dyn_info%elec_long_period > 0) then
      if (mod(dyn_info%nbupdate_period,dyn_info%elec_long_period) /= 0) &
        call error_msg( &
        'Read_Ctrl_Dynamics> Error: For elec_long_period or nbupdate_period,'//&
        ' nbupdate_period must be a multiple of elec_long_period.')
    end if

    return

  end subroutine read_ctrl_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynamics
  !> @brief        Sets up dynamics information.
  !> @authors      TM
  !! @param[in]    dyn_info   : DYNAMICS  section control parameters information
  !! @param[in]    bound_info : BOUNDARY  section control parameters information
  !! @param[in]    res_info   : RESTRAINT section control parameters information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] dynamics   : dynamics information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynamics(dyn_info, bound_info, res_info, molecule, dynamics)

    ! formal arguments
    type(s_dyn_info),        intent(in)    :: dyn_info
    type(s_boundary_info),   intent(in)    :: bound_info
    type(s_res_info),        intent(in)    :: res_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_dynamics),        intent(inout) :: dynamics

    ! local variables
    integer                  :: i, id
    integer                  :: ierror, iseed_init_vel, iseed


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
    dynamics%lbupdate_period   = dyn_info%lbupdate_period
    dynamics%elec_long_period  = dyn_info%elec_long_period
    dynamics%initial_time      = dyn_info%initial_time
    dynamics%verbose           = dyn_info%verbose
    dynamics%thermo_period     = dyn_info%thermo_period
    dynamics%baro_period       = dyn_info%baro_period  
    dynamics%istart_step       = 1
    dynamics%iend_step         = dynamics%nsteps
    dynamics%target_md         = dyn_info%target_md
    dynamics%steered_md        = dyn_info%steered_md
    dynamics%initial_rmsd      = dyn_info%initial_rmsd
    dynamics%final_rmsd        = dyn_info%final_rmsd
    dynamics%gen_velocity      = dyn_info%gen_velocity

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
    !$omp parallel private(id)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    dynamics%iseed_omp(id+1) = iseed + id*nthread + 1000 * my_country_no
    !$omp end parallel

    dynamics%iseed_init_velocity = iseed_init_vel + 1000 * my_country_no

    if (bound_info%type == BoundaryTypePBC) then
      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .false.
    else if (bound_info%type == BoundaryTypeNOBC) then
      dynamics%stop_com_translation = .true.
      dynamics%stop_com_rotation    = .true.
    end if

    do i = 1, res_info%nfunctions
      if (res_info%function(i) == RestraintsFuncPOSI .or. &
          res_info%function(i) == RestraintsFuncEM) then
        dynamics%stop_com_translation = .false.
        dynamics%stop_com_rotation    = .false.
      end if
    end do


    ! update number of degrees of freedom
    !
    if (dynamics%stop_com_translation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtracted 3 translational degrees of freedom.'
        write(MsgOut,'(a)') ' '
      endif

      call update_num_deg_freedom('After removing translation', &
                                  -3, molecule%num_deg_freedom)

    end if

    if (dynamics%stop_com_rotation) then

      if (main_rank)  then
        write(MsgOut,'(a)') &
        'Setup_Dynamics> Subtracted 3 rotational degrees of freedom.'
        write(MsgOut,'(a)') ' '
      endif

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
    call random_init_omp(dynamics%iseed_omp)

    return

  end subroutine setup_dynamics


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_md
  !> @brief        Performs MD simulation.
  !! @authors      YS, CK
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_md(output, grotop, molecule, domain, enefunc, dynvars, &
                    dynamics, pairlist, boundary, ensemble, comm)
    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_grotop),           intent(in   ) :: grotop
    type(s_molecule),         intent(inout) :: molecule
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_comm),             intent(inout) :: comm

    if (dynamics%target_md) then
      enefunc%restraint_rmsd_target = .true.
    else
      enefunc%restraint_rmsd_target = .false.
    end if

    call open_output(output)
    select case (dynamics%integrator)

    case (IntegratorLEAP)

      call leapfrog_dynamics(grotop, molecule, output, domain, enefunc, &
                             dynvars, dynamics, pairlist, boundary,     &
                             ensemble, comm)

    case (IntegratorVVER)

      call vverlet_dynamics (grotop, molecule, output, domain, enefunc, &
                             dynvars, dynamics, pairlist, boundary,     &
                             ensemble, comm)

    end select
    call close_output(output)

    return

  end subroutine run_md

end module cg_dynamics_mod
