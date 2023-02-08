!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_remd_mod
!> @brief   Replica exhcnage molecular dynamics simulation
!! @authors Takaharu Mori (TM), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_remd_mod

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
  use sp_remd_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_constraints_mod
  use sp_boundary_str_mod
  use sp_boundary_mod
  use sp_pairlist_str_mod
  use sp_pairlist_mod
  use sp_enefunc_str_mod
  use sp_enefunc_restraints_mod
  use sp_output_str_mod
  use sp_output_mod
  use sp_domain_str_mod
  use sp_communicate_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_str_mod
  use sp_fep_energy_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
  use sp_energy_pme_mod

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! structures
  type, public :: s_rep_info
    integer                           :: dimension          = 1
    integer                           :: exchange_period    = 100
    integer                           :: iseed              = 3141592
    logical                           :: analysis_grest     = .false.
    integer,              allocatable :: types        (:)
    integer,              allocatable :: nreplicas    (:)
    integer,              allocatable :: nreplicas_x  (:)
    integer,              allocatable :: nreplicas_y  (:)
    integer,              allocatable :: nreplicas_z  (:)
    character(MaxLine),   allocatable :: parameters   (:)
    logical,              allocatable :: cyclic_params(:)
    character(MaxLine),   allocatable :: rest_function(:)
    character(MaxLine),   allocatable :: select_index (:)
    character(MaxLine),   allocatable :: param_type   (:)
    ! auto-adjust
    logical,              allocatable :: param_tuning(:)
    real(wp),             allocatable :: tgt_exc_prob(:)
    real(wp),             allocatable :: mgn_exc_prob(:)
    integer,              allocatable :: trial_freq(:)
    integer,              allocatable :: eq_cycle(:)
    real(wp),             allocatable :: param_grid(:)
    real(wp),             allocatable :: max_param_shift(:)
    character(MaxLine),   allocatable :: fix_terminal(:)
  end type s_rep_info

  ! variables
  integer, public, parameter          :: RemdMaxNdimensions = 999
  integer, public, parameter          :: RemdMaxNumbrellas  = 999

  ! subroutines
  public  :: show_ctrl_remd
  public  :: read_ctrl_remd
  public  :: setup_solute_tempering
  public  :: setup_solute_tempering_pio
  public  :: setup_remd
  public  :: setup_remd_pio
  private :: setup_reus
  private :: setup_remd_solute_tempering
  private :: setup_remd_solute_tempering_pio
  private :: setup_remd_solute_tempering_bonds
  private :: setup_remd_solute_tempering_angles
  private :: setup_remd_solute_tempering_dihedrals
  private :: setup_remd_solute_tempering_rb_dihedrals
  private :: setup_remd_solute_tempering_impropers
  private :: setup_remd_solute_tempering_cmaps
  private :: setup_remd_solute_tempering_contacts
  private :: setup_remd_solute_tempering_lj
  public  :: run_remd
  private :: perform_replica_exchange
  private :: remd_autoadj_displacement
  private :: temperature_remd
  private :: pressure_remd
  private :: surface_tension_remd
  private :: restraint_remd
  private :: solute_tempering_remd
  private :: assign_condition
  private :: assign_condition_solute_tempering
  private :: assign_condition_solute_tempering_internal
  private :: assign_condition_solute_tempering_lj
  private :: get_neighbouring_parmsetid
  private :: get_replica_temperature
  private :: compute_energy_restraints_remd
  ! FEP
  public  :: run_remd_fep
  private :: perform_replica_exchange_fep
  private :: setup_fep_remd
  private :: fep_remd
  private :: assign_condition_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_remd
  !> @brief        show REMD section usage
  !! @authors      NT, TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_remd(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('remd')

        write(MsgOut,'(A)') '[REMD]'
        write(MsgOut,'(A)') 'dimension         = 1'
        write(MsgOut,'(A)') 'exchange_period   = 100'
        write(MsgOut,'(A)') 'iseed             = 3141592'
        write(MsgOut,'(A)') 'analysis_grest    = NO'
        write(MsgOut,'(A)') ''
        write(MsgOut,'(A)') 'type1             = TEMPERATURE'
        write(MsgOut,'(A)') 'nreplica1         = 4'
        write(MsgOut,'(A)') 'parameters1       = 300 310 320 330'
        write(MsgOut,'(A)') 'cyclic_params1    = NO'
        write(MsgOut,'(A)') ''
        write(MsgOut,'(A)') '# type2            = PRESSURE'
        write(MsgOut,'(A)') '# nreplica2        = 4'
        write(MsgOut,'(A)') '# parameters2      = 1 10 100 1000'
        write(MsgOut,'(A)') '# cyclic_params2   = NO'
        write(MsgOut,'(A)') '#'
        write(MsgOut,'(A)') '# type3            = GAMMA'
        write(MsgOut,'(A)') '# nreplica3        = 4'
        write(MsgOut,'(A)') '# parameters3      = 0 6 12 18'
        write(MsgOut,'(A)') '# cyclic_params3   = NO'
        write(MsgOut,'(A)') '#'
        write(MsgOut,'(A)') '# type4            = RESTRAINT'
        write(MsgOut,'(A)') '# nreplica4        = 4'
        write(MsgOut,'(A)') '# cyclic_params4   = NO'
        write(MsgOut,'(A)') '# rest_function4   = 1'
        write(MsgOut,'(A)') '#'
        write(MsgOut,'(A)') '# type5            = REST'
        write(MsgOut,'(A)') '# nreplica5        = 4'
        write(MsgOut,'(A)') '# parameters5      = 300 310 320 330'
        write(MsgOut,'(A)') '# cyclic_params5   = NO'
        write(MsgOut,'(A)') '# select_index5    = 1'
        write(MsgOut,'(A)') '# param_type5      = ALL'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('remd')

        write(MsgOut,'(A)') '[REMD]'
        write(MsgOut,'(A)') 'exchange_period   = 100'
        write(MsgOut,'(A)') 'type1             = TEMPERATURE'
        write(MsgOut,'(A)') 'nreplica1         = 4'
        write(MsgOut,'(A)') 'parameters1       = 300 310 320 330'
        write(MsgOut,'(A)') ' '

      end select


    end if

    return

  end subroutine show_ctrl_remd
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_remd
  !> @brief        read REMD section in the control file
  !! @authors      TM
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   rep_info : REPLICA section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_remd(handle, rep_info)

    ! parameters
    character(*),            parameter     :: Section = 'Remd'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_rep_info),        intent(inout) :: rep_info

    ! local variables
    integer                  :: i, trep
    character(30)            :: cdim, numtmp, partmp, vartmp
    character(30)            :: cnstmp, reftmp, fnctmp, cyctmp, seltmp
    character(30)            :: typtmp, adjtmp


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'dimension',       &
                               rep_info%dimension)
    call read_ctrlfile_integer(handle, Section, 'exchange_period', &
                               rep_info%exchange_period)
    call read_ctrlfile_integer(handle, Section, 'iseed',           &
                               rep_info%iseed)
    call read_ctrlfile_logical(handle, Section, 'analysis_grest',  &
                               rep_info%analysis_grest)

    ! error check
    !
    if (rep_info%dimension > RemdMaxNdimensions .or. rep_info%dimension == 0) &
      call error_msg('Read_Ctrl_Remd> error in dimension in [REMD]')

    ! read allocatable variables
    !
    allocate(rep_info%types        (rep_info%dimension),   &
             rep_info%nreplicas    (rep_info%dimension),   &
             rep_info%nreplicas_x  (rep_info%dimension),   &
             rep_info%nreplicas_y  (rep_info%dimension),   &
             rep_info%nreplicas_z  (rep_info%dimension),   &
             rep_info%parameters   (rep_info%dimension),   &
             rep_info%cyclic_params(rep_info%dimension),   &
             rep_info%rest_function(rep_info%dimension),   &
             rep_info%select_index (rep_info%dimension),   &
             rep_info%param_type   (rep_info%dimension),   &
             rep_info%param_tuning (rep_info%dimension),   &
             rep_info%tgt_exc_prob (rep_info%dimension),   &
             rep_info%mgn_exc_prob (rep_info%dimension),   &
             rep_info%trial_freq   (rep_info%dimension),   &
             rep_info%eq_cycle     (rep_info%dimension),   &
             rep_info%param_grid   (rep_info%dimension),   &
             rep_info%max_param_shift(rep_info%dimension), &
             rep_info%fix_terminal (rep_info%dimension))

    rep_info%types          (1:rep_info%dimension) = RemdTemperature
    rep_info%nreplicas      (1:rep_info%dimension) = 0
    rep_info%nreplicas_x    (1:rep_info%dimension) = 0
    rep_info%nreplicas_y    (1:rep_info%dimension) = 0
    rep_info%nreplicas_z    (1:rep_info%dimension) = 0
    rep_info%parameters     (1:rep_info%dimension) = ''
    rep_info%cyclic_params  (1:rep_info%dimension) = .false.
    rep_info%rest_function  (1:rep_info%dimension) = ''
    rep_info%select_index   (1:rep_info%dimension) = ''
    rep_info%param_type     (1:rep_info%dimension) = 'ALL'
    rep_info%param_tuning   (1:rep_info%dimension) = .false.
    rep_info%tgt_exc_prob   (1:rep_info%dimension) = 0.25_wp
    rep_info%mgn_exc_prob   (1:rep_info%dimension) = 0.05_wp
    rep_info%trial_freq     (1:rep_info%dimension) = 48
    rep_info%eq_cycle       (1:rep_info%dimension) = 2
    rep_info%param_grid     (1:rep_info%dimension) = 0.1_wp
    rep_info%max_param_shift(1:rep_info%dimension) = 20.0_wp
    rep_info%fix_terminal   (1:rep_info%dimension) = 'BOTTOM'

    do i = 1, rep_info%dimension
      write(cdim,'(i0)') i

      vartmp = 'type'          // cdim
      call read_ctrlfile_type   (handle, Section, vartmp,  &
                                 rep_info%types(i), RemdTypes)

      numtmp = 'nreplica'      // cdim
      call read_ctrlfile_integer(handle, Section, numtmp,  &
                                 rep_info%nreplicas(i))

      numtmp = 'nreplicax'      // cdim
      call read_ctrlfile_integer(handle, Section, numtmp,  &
                                 rep_info%nreplicas_x(i))

      numtmp = 'nreplicay'      // cdim
      call read_ctrlfile_integer(handle, Section, numtmp,  &
                                 rep_info%nreplicas_y(i))

      numtmp = 'nreplicaz'      // cdim
      call read_ctrlfile_integer(handle, Section, numtmp,  &
                                 rep_info%nreplicas_z(i))

      partmp = 'parameters'    // cdim
      call read_ctrlfile_string (handle, Section, partmp,  &
                                 rep_info%parameters(i))

      cyctmp = 'cyclic_params' // cdim
      call read_ctrlfile_logical(handle, Section, cyctmp,  &
                                 rep_info%cyclic_params(i))

      fnctmp = 'rest_function' // cdim
      call read_ctrlfile_string (handle, Section, fnctmp,  &
                                 rep_info%rest_function(i))

      seltmp = 'select_index' // cdim
      call read_ctrlfile_string (handle, Section, seltmp,  &
                                 rep_info%select_index(i))

      typtmp = 'param_type'   // cdim
      call read_ctrlfile_string (handle, Section, typtmp,  &
                                 rep_info%param_type(i))

      adjtmp = 'param_tuning' // cdim
      call read_ctrlfile_logical(handle, Section, adjtmp, &
                                 rep_info%param_tuning(i))

      adjtmp = 'tgt_exc_prob' // cdim
      call read_ctrlfile_real   (handle, Section, adjtmp, &
                                 rep_info%tgt_exc_prob(i))

      adjtmp = 'mgn_exc_prob' // cdim
      call read_ctrlfile_real   (handle, Section, adjtmp, &
                                 rep_info%mgn_exc_prob(i))

      adjtmp = 'trial_freq'   // cdim
      call read_ctrlfile_integer(handle, Section, adjtmp, &
                                 rep_info%trial_freq(i))

      adjtmp = 'eq_cycle'     // cdim
      call read_ctrlfile_integer(handle, Section, adjtmp, &
                                 rep_info%eq_cycle(i))

      adjtmp = 'param_grid'   // cdim
      call read_ctrlfile_real   (handle, Section, adjtmp, &
                                 rep_info%param_grid(i))

      adjtmp = 'max_param_shift' // cdim
      call read_ctrlfile_real   (handle, Section, adjtmp, &
                                 rep_info%max_param_shift(i))

      adjtmp = 'fix_terminal' // cdim
      call read_ctrlfile_string (handle, Section, adjtmp, &
                                 rep_info%fix_terminal(i))

    end do


    ! error check
    !
    do i = 1, rep_info%dimension
      if (rep_info%nreplicas(i) == 0) &
        rep_info%nreplicas(i) = rep_info%nreplicas_x(i) &
                               *rep_info%nreplicas_y(i) &
                               *rep_info%nreplicas_z(i)
      if (rep_info%nreplicas_x(i) /= 0 .or. &
          rep_info%nreplicas_y(i) /= 0 .or. &
          rep_info%nreplicas_z(i) /= 0) then
        if (rep_info%nreplicas(i) /= rep_info%nreplicas_x(i)  &
                                    *rep_info%nreplicas_y(i)  &
                                    *rep_info%nreplicas_z(i)) &
          call error_msg('Read_Ctrl_Remd> error in nreplicas in [REMD]')
      end if
    end do

    do i = 1, rep_info%dimension
      if (rep_info%nreplicas(i) == 0) &
        call error_msg('Read_Ctrl_Remd> error in nreplicas in [REMD]')
    end do

    trep = product(rep_info%nreplicas(1:rep_info%dimension))
    if (trep == 0) &
      call error_msg('Read_Ctrl_Remd> error in nreplicas in [REMD]')


    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Remd> Replica information'
      do i = 1, rep_info%dimension

        write(MsgOut,'(A20,I10,A20,I10)') &
          '  dimension       = ', i,      &
          '  nreplica        = ', rep_info%nreplicas(i)

        if (rep_info%types(i) == RemdTemperature) then
          write(MsgOut,'(A)')   '  variable        = TEMPERATURE'
          write(MsgOut,'(A,A)') '  parameters      = ', trim(rep_info%parameters(i))
        else if (rep_info%types(i) == RemdPressure) then
          write(MsgOut,'(A)')   '  variable        = PRESSURE'
          write(MsgOut,'(A,A)') '  parameters      = ', trim(rep_info%parameters(i))
        else if (rep_info%types(i) == RemdGamma) then
          write(MsgOut,'(A)')   '  variable        = GAMMA'
          write(MsgOut,'(A,A)') '  parameters      = ', trim(rep_info%parameters(i))
        else if (rep_info%types(i) == RemdRestraint) then
          write(MsgOut,'(A)')   '  type            = RESTRAINT'
          write(MsgOut,'(A,A)') '  rest_function   = ', trim(rep_info%rest_function(i))
        else if (rep_info%types(i) == RemdSoluteTempering) then
          write(MsgOut,'(A)')   '  variable        = SOLUTE TEMPERATURE'
          write(MsgOut,'(A,A)') '  select_index    = ', trim(rep_info%select_index(i))
          write(MsgOut,'(A,A)') '  param_type      = ', trim(rep_info%param_type(i))
        else if (rep_info%types(i) == RemdAlchemy) then
          write(MsgOut,'(A)')   '  type            = ALCHEMY'
        else if (rep_info%types(i) == RemdAlchemyRest) then
          if ( rep_info%dimension > 1 ) then
            call error_msg('Read_Ctrl_Remd> dimension > 1 is not allowed for FEP/REST')
          end if
          write(MsgOut,'(A)')       '  type            = ALCHEMYREST (FEP/REST)'
          write(MsgOut,'(A20,A30)') '  select_index    = ', rep_info%select_index(i)
          write(MsgOut,'(A20,A30)') '  param_type      = ', rep_info%param_type(i)
        end if

        if (rep_info%cyclic_params(i)) then
          write(MsgOut,'(A)') '  cyclic_params   = yes'
        else
          write(MsgOut,'(A)') '  cyclic_params   = no'
        end if

        write(MsgOut,'(A)') ''
      end do

    end if

    return

  end subroutine read_ctrl_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_tempering
  !> @brief        presetup of replica exchange solute tempering
  !! @authors      MK
  !! @param[in]    rep_info    : REMD section control parameters information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] restraints  : restraints information
  !! @param[in]    cons_info   : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_tempering(rep_info, molecule, restraints, cons_info)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_restraints),      intent(inout) :: restraints
    type(s_cons_info),       intent(in)    :: cons_info

    ! local variables
    integer :: i, j, k, ndata
    logical :: has_rest


    has_rest = .false.
    do i = 1, rep_info%dimension
      if ((rep_info%types(i) == RemdSoluteTempering) .or. &
          (rep_info%types(i) == RemdAlchemyRest)) then
        has_rest = .true.
      end if
    end do
    if (.not. has_rest) return

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Solute_Tempering> preparation of REST'
    end if

    do i = 1, rep_info%dimension
      if ((rep_info%types(i) == RemdSoluteTempering) .or. &
          (rep_info%types(i) == RemdAlchemyRest)) then

        ndata = split_num(rep_info%select_index(i))

        if (ndata /= 1) &
          call error_msg('  Setup_Solute_Tempering> Error in index for REST')

        read(rep_info%select_index(i),*) ndata

        if (restraints%num_groups < ndata) then
          call error_msg('  Setup_Solute_Tempering> group out of range')
        end if

        ! have to change water name to avoid table for solute water mols
        !
        do j = 1, restraints%num_atoms(ndata)

          k = restraints%atomlist(j,ndata)
          if (molecule%residue_name(k) == cons_info%water_model) then
            call error_msg('Error> water molecules must not be in REST solute.')

            molecule%residue_name(k) = DummyWaterName
            if (main_rank) then
              write(MsgOut,'(a,i0)') '  replacing residue name of atom: ', k
            end if
          end if

        end do

      end if
    end do

    if (main_rank) then
      write(MsgOut,*)
    end if

    return

  end subroutine setup_solute_tempering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_tempering_pio
  !> @brief        presetup of replica exchange solute tempering (w/ pio)
  !! @authors      JJ
  !! @param[in]    rep_info    : REMD section control parameters information
  !! @param[inout] restraints  : restraints information
  !! @param[in]    cons_info   : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_tempering_pio(rep_info, restraints, cons_info)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_restraints),      intent(inout) :: restraints
    type(s_cons_info),       intent(in)    :: cons_info

    ! local variables
    integer :: i, j, k, ndata
    logical :: has_rest


    has_rest = .false.
    do i = 1, rep_info%dimension
      if (rep_info%types(i) == RemdSoluteTempering) then
        has_rest = .true.
      end if
    end do
    if (.not. has_rest) return

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Solute_Tempering> preparation of REST'
    end if

!   do i = 1, rep_info%dimension

!     if (rep_info%types(i) == RemdSoluteTempering) then

!       ndata = split_num(rep_info%select_index(i))

!       if (ndata /= 1) &
!         call error_msg('  Setup_Solute_Tempering> Error in index for REST')

!       read(rep_info%select_index(i),*) ndata

!       if (restraints%num_groups < ndata) then
!         call error_msg('  Setup_Solute_Tempering> group out of range')
!       end if

!     end if
!   end do

    if (main_rank) then
      write(MsgOut,*)
    end if

    return

  end subroutine setup_solute_tempering_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd
  !> @brief        setup REMD
  !! @authors      TM
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    dynamics   : dynamics information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] domain     : domain information
  !! @param[inout] restraints : restraints information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd(rep_info, rst, boundary, dynamics, molecule, domain, &
                        restraints, ensemble, enefunc, remd, alchemy)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_rst),             intent(in)    :: rst
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_restraints),      intent(inout) :: restraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),    target, intent(inout) :: remd
    type(s_alchemy),optional,intent(inout) :: alchemy

    ! local variables
    integer                    :: i, j, k, m
    integer                    :: max_nreplicas, cycles
    integer                    :: itmp, jtmp, ktmp
    integer                    :: comp, dsta, dend, found, dtotal
    integer                    :: replicaid, parmsetid
    integer                    :: ncycle, icycle, nlen, ixx
    integer                    :: iref, icomp
    character(MaxLine)         :: param
    character(10)              :: tmp, frmt, cd
    character(18)              :: frmt2
    character(3)               :: ctrep
    logical                    :: keep, double
    real(wp), allocatable      :: ddata(:)
    integer,       allocatable :: ndigit(:), parmid(:)
    integer                    :: funcid
    character(MaxLine)         :: ft
    integer,       allocatable :: sendbuf(:), recvbuf(:,:)
    ! FEP
    integer                    :: lambid
    real(wp)                   :: soltemp

    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Remd> Replica information'
      write(MsgOut,'(A)') ''
    end if


    ! initialize
    !
    remd%dimension       = rep_info%dimension
    remd%total_nreplicas = product(rep_info%nreplicas(1:remd%dimension))
    max_nreplicas        = maxval (rep_info%nreplicas(1:remd%dimension))
    remd%analysis_grest  = rep_info%analysis_grest

    call alloc_remd(remd, RemdReplicas, remd%dimension,      &
                    remd%total_nreplicas, max_nreplicas)
    call alloc_remd(remd, RemdReplicas_rst, remd%dimension,  &
                    remd%total_nreplicas, max_nreplicas)

    do i = 1, remd%dimension
      remd%nreplicas(i)     = rep_info%nreplicas(i)
      remd%types(i)         = rep_info%types(i)
      remd%cyclic_params(i) = rep_info%cyclic_params(i)
      remd%autoadj(i)%param_tuning    = rep_info%param_tuning(i)
      remd%autoadj(i)%tgt_exc_prob    = rep_info%tgt_exc_prob(i)
      remd%autoadj(i)%mgn_exc_prob    = rep_info%mgn_exc_prob(i)
      remd%autoadj(i)%trial_freq      = rep_info%trial_freq(i)
      remd%autoadj(i)%eq_cycle        = rep_info%eq_cycle(i)
      remd%autoadj(i)%param_grid      = rep_info%param_grid(i)
      remd%autoadj(i)%max_param_shift = rep_info%max_param_shift(i)
      ft = rep_info%fix_terminal(i)
      call toupper(ft)
      select case (ft)
        case ("BOTTOM")
          remd%autoadj(i)%fix_terminal = RemdAutoFixBottom
        case ("TOP")
          remd%autoadj(i)%fix_terminal = RemdAutoFixTop
        case default
          call error_msg('Error: fix_terminal must be BOTTOM/TOP.')
      end select
    end do

    ! check ensemble 
    do i = 1, remd%dimension
      if (remd%types(i) == RemdGamma .and.  &
          ensemble%ensemble /= EnsembleNPgT) &
         call error_msg('Setup_Remd> Ensemble should be NPgT in gamma-REMD')
      ! FEP
      if ((remd%types(i) == RemdAlchemy) .or. &
        (remd%types(i) == RemdAlchemyRest)) then
        if (ensemble%tpcontrol /= TpcontrolBussi .and. &
            ensemble%tpcontrol /= TpcontrolNO)         &
            call error_msg('Setup_Remd> tpcontrol should be Bussi or No in FEP')
      end if
    end do

    ! check cyclic_params
    !
    do i = 1, remd%dimension
      if (remd%cyclic_params(i)) then
        if (mod(remd%nreplicas(i),2) /= 0) &
          call error_msg('Setup_Remd> nreplicas must be even,'//&
                         ' if cyclic_params is specified')
      end if
    end do


    ! setup remd_iseed, num_criteria, and num_exchange
    !
    if (rst%rstfile_type == RstfileTypeRemd) then
      remd%iseed = rst%iseed_remd
      do i = 1, remd%total_nreplicas
        remd%repid2parmsetid(i) = rst%repid2parmsetid(i)
        do j = 1, remd%dimension
          remd%num_criteria (i,j,1:2) = rst%num_criteria (i,j,1:2)
          remd%num_exchanges(i,j,1:2) = rst%num_exchanges(i,j,1:2)
        end do
      end do

      ! check the restart information consistency
      !
#ifdef HAVE_MPI_GENESIS
      allocate(sendbuf(remd%total_nreplicas))
      allocate(recvbuf(remd%total_nreplicas, remd%total_nreplicas))

      sendbuf(:) = remd%repid2parmsetid(:)
      call mpi_allgather(sendbuf, remd%total_nreplicas, mpi_integer, &
                         recvbuf, remd%total_nreplicas, mpi_integer, mpi_comm_airplane, ierror)
      do i = 1, remd%total_nreplicas
        iref = recvbuf(i,1)
        do j = 2, remd%total_nreplicas
          icomp = recvbuf(i,j)
          if (iref /= icomp) &
            call error_msg("Setup_Remd> REMD restart information differs between replicas. Perhaps some restart files were not output correctly in the previous run due to a computer accident.")
        end do
      end do

      do k = 1, remd%dimension
        do m = 1, 2
          sendbuf(:) = remd%num_criteria(:,k,m)
          call mpi_allgather(sendbuf, remd%total_nreplicas, mpi_integer, &
                             recvbuf, remd%total_nreplicas, mpi_integer, mpi_comm_airplane, ierror)
          do i = 1, remd%total_nreplicas
            iref = recvbuf(i,1)
            do j = 2, remd%total_nreplicas
              icomp = recvbuf(i,j)
              if (iref /= icomp) &
                call error_msg("Setup_Remd> REMD restart information differs between replicas. Perhaps some restart files were not output correctly in the previous run due to a computer accident.")
            end do
          end do

          sendbuf(:) = remd%num_exchanges(:,k,m)
          call mpi_allgather(sendbuf, remd%total_nreplicas, mpi_integer, &
                             recvbuf, remd%total_nreplicas, mpi_integer, mpi_comm_airplane, ierror)
          do i = 1, remd%total_nreplicas
            iref = recvbuf(i,1)
            do j = 2, remd%total_nreplicas
              icomp = recvbuf(i,j)
              if (iref /= icomp) &
                call error_msg("Setup_Remd> REMD restart information differs between replicas. Perhaps some restart files were not output correctly in the previous run due to a computer accident.")
            end do
          end do

        end do
      end do

      deallocate(sendbuf,recvbuf)
#endif

    else
      remd%iseed = rep_info%iseed
      do i = 1, remd%total_nreplicas
        remd%repid2parmsetid(i) = i
        do j = 1, remd%dimension
          remd%num_criteria (i,j,1:2) = 0
          remd%num_exchanges(i,j,1:2) = 0
        end do
      end do
    end if

    ! setup permutation function
    !
    do i = 1, remd%total_nreplicas
      parmsetid = remd%repid2parmsetid(i)
      remd%parmsetid2repid(parmsetid) = i
    end do


    ! setup replica exchange period
    !
    if (rep_info%exchange_period > 0) then
      if (dynamics%eneout_period > rep_info%exchange_period) then
        call error_msg('Setup_Remd> (eneout_period <= exchange_period)'//&
                       ' is required')
      else
        if (mod(rep_info%exchange_period,dynamics%eneout_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,eneout_period)'//&
                         ' must be zero)')
        end if
        if (mod(rep_info%exchange_period,dynamics%thermo_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,thermo_period)'//&
                         ' must be zero)')
        end if
        if (mod(rep_info%exchange_period,dynamics%baro_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,baro_period)'//&
                         ' must be zero)')
        end if
        ! FEP
        do i = 1, remd%dimension
          if ((remd%types(i) == RemdAlchemy) .or. &
            (remd%types(i) == RemdAlchemyRest)) then
            if (mod(rep_info%exchange_period,dynamics%fepout_period) /= 0) then
              call error_msg('Setup_Remd> mod(exchange_period,fepout_period)'//&
                ' must be zero')
            end if
          end if
        end do
      end if
      cycles = 2*rep_info%exchange_period*remd%dimension
      if (mod(dynamics%nsteps,cycles) /= 0) then
        call error_msg('Setup_Remd> mod(nsteps,(2*exchange_period*dimension))'//&
                       ' must be zero')
      else
        remd%equilibration_only = .false.
        remd%exchange_period = rep_info%exchange_period
        remd%ncycles = dynamics%nsteps/remd%exchange_period
      end if
    else if (rep_info%exchange_period == 0) then
      remd%equilibration_only = .true.
      remd%exchange_period = dynamics%nsteps
      remd%ncycles = 1
    else
      call error_msg('Setup_Remd> error in exchange_period')
    end if

    ! check parameter tuning
    !
    if (any(remd%autoadj(:)%param_tuning)) then
      if (remd%exchange_period <= 0) then
        remd%autoadj(:)%param_tuning = .false.
        if (main_rank) then
          write(MsgOut,'("Warning> param_tuning is disabled &
                          &since exchange is disabled.")')
          write(MsgOut,*)
        end if
      else
        ! clear all the exchange probs
        do i = 1, remd%total_nreplicas
          do j = 1, remd%dimension
            remd%num_criteria (i,j,1:2) = 0
            remd%num_exchanges(i,j,1:2) = 0
          end do
        end do
      end if
    end if


    ! split rep_info%parameters (char) into parameters (real)
    !
    allocate(ddata(max_nreplicas))

    do i = 1, remd%dimension
      if ((remd%types(i) /= RemdRestraint) .and. &
          (remd%types(i) /= RemdAlchemy)) then
        param = rep_info%parameters(i)
        if (split_num(param) /= remd%nreplicas(i)) &
          call error_msg('Setup_Remd> "nreplicas" and number of data in "parameters" are inconsistent')

        call split(remd%nreplicas(i), max_nreplicas, param, ddata)

        do j = 1, remd%nreplicas(i)
          remd%dparameters(i,j) = ddata(j)
        end do
      end if
    end do

    deallocate(ddata)


    ! setup replica-exchange umbrella sampling (REUS)
    !
    call setup_reus(rep_info, rst, restraints, enefunc, remd)


    ! setup replica exchange with solute tempering (REST)
    !
    call setup_remd_solute_tempering(rep_info, molecule, domain, remd, &
                                     restraints, enefunc)

    do i = 1, remd%dimension
      if (main_rank .and. remd%autoadj(i)%param_tuning) then
        write(MsgOut,'(a,i8)') 'Setup_Remd> tuning parameters for dim ', i
        write(MsgOut,'(a,f8.3)') '  tgt_exc_prob    = ', &
                                 remd%autoadj(i)%tgt_exc_prob
        write(MsgOut,'(a,f8.3)') '  mgn_exc_prob    = ', &
                                 remd%autoadj(i)%mgn_exc_prob
        write(MsgOut,'(a,i8)')   '  trial_freq      = ', &
                                 remd%autoadj(i)%trial_freq
        write(MsgOut,'(a,i8)')   '  eq_cycle        = ', &
                                 remd%autoadj(i)%eq_cycle
        write(MsgOut,'(a,f8.3)') '  param_grid      = ', &
                                 remd%autoadj(i)%param_grid
        write(MsgOut,'(a,f8.3)') '  max_param_shift = ', &
                                 remd%autoadj(i)%max_param_shift
        if (remd%autoadj(i)%fix_terminal == RemdAutoFixBottom) then
          write(MsgOut,'(a)')    '  fix_terminal    = BOTTOM'
        else
          write(MsgOut,'(a)')    '  fix_terminal    = TOP'
        end if
        write(MsgOut,*)
      end if
    end do

    if (present(alchemy)) then
      ! FEP: setup lambda-exchange FEP (FEP-REMD) or fep with 
      ! generalized replica exchange with solute tempering (FEP-gREST)
      call setup_fep_remd(rep_info, molecule, domain, remd, &
                          restraints, enefunc, alchemy)
    end if

    ! generate multi-dimensional parameter sets
    !
    allocate(ndigit(remd%dimension),parmid(remd%dimension))

    do i = 1, remd%dimension
      do j = 1, 100
        comp = 10**j
        if (remd%nreplicas(i) < comp) then
          ndigit(i) = j
          exit
        end if
      end do
    end do
    dtotal = sum(ndigit(1:remd%dimension))

    write(cd,'(i10)') dtotal
    frmt = '(i' // trim(adjustl(cd)) // '.' // trim(adjustl(cd)) // ')'

    found     = 0
    parmsetid = 0
    do
      found = found + 1
      write(tmp,frmt) found
      keep = .true.
      dend = 0

      do j = 1, remd%dimension
        dsta = dend + 1
        dend = dsta + ndigit(j) - 1
        read(tmp(dsta:dend),*) itmp
        if (itmp > remd%nreplicas(j) .or. itmp == 0) then
          keep = .false.
          exit
        end if
      end do

      if (keep) then
        parmsetid = parmsetid + 1
        dend = 0
        do j = 1, remd%dimension
          dsta = dend + 1
          dend = dsta + ndigit(j) - 1
          read(tmp(dsta:dend),*) parmid(j)
          remd%parmidsets(parmsetid,j) = parmid(j)
        end do
        if (parmsetid == remd%total_nreplicas) exit
      end if
    end do

    deallocate(ndigit,parmid)


    ! assign itnitial parameters
    !
    replicaid = my_country_no + 1
    parmsetid = remd%repid2parmsetid(replicaid)
    if (present(alchemy)) then
      ! FEP
      call assign_condition_fep(parmsetid, remd, domain, ensemble, enefunc, alchemy)
    else
      call assign_condition(parmsetid, remd, domain, ensemble, enefunc)
    end if


    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') '  ParmsetID'

      write(ctrep,'(I3)') remd%total_nreplicas
      frmt2 = '(I11 , A,' // ctrep // 'I11  )'

      do i = 1, remd%total_nreplicas
        replicaid = i
        parmsetid = remd%repid2parmsetid(replicaid)
        write(MsgOut,frmt2) parmsetid, ' = ', &
          remd%parmidsets(parmsetid,1:remd%dimension)
      end do

      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  Parameters'

      do i = 1, remd%dimension
        if ((remd%types(i) /= RemdRestraint) .and. &
            (remd%types(i) /= RemdAlchemy)   .and. &
            (remd%types(i) /= RemdAlchemyRest)) then
          write(ctrep,'(I3)') remd%nreplicas(i)
          frmt2 = '(A,I4,A,' // ctrep // 'F11.3)'
          write(MsgOut,frmt2) &
            '    Dim', i, ' = ', remd%dparameters(i,1:remd%nreplicas(i))
        else
          write(ctrep,'(I3)') remd%nreplicas(i)
          frmt2 = '(A,I4,A,' // ctrep // 'I11  )'
          write(MsgOut,frmt2) &
            '    Dim', i, ' = ', remd%iparameters(i,1:remd%nreplicas(i))
        end if
      end do

      do i = 1, remd%dimension
        if (remd%types(i) == RemdRestraint) then
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '  Restraints'
          exit
        end if
      end do

      do i = 1, size(remd%umbrid2numfuncs(:))
        write(MsgOut,'(4X,A,I10)') 'umbrella potential index = ', i
          do k = 1, remd%umbrid2numfuncs(i)
            write(MsgOut,'(6X,A,I7,2X,A,F11.3,2X,A,F11.3)')         &
                         'function = ',  remd%umbrid2funclist(i,k), &
                         'constant = ',  remd%rest_constants(i,k),  &
                         'reference = ', remd%rest_reference(i,k)
          end do
        write(MsgOut,'(A)') ''
      end do

      ! FEP: show replica information of fepremd or feprest
      do i = 1, remd%dimension
        if ((remd%types(i) == RemdAlchemy) .or. &
            (remd%types(i) == RemdAlchemyRest)) then
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '  FEP Windows'
          do j = 1, remd%nreplicas(i)

            lambid  = remd%iparameters(i,j)
            soltemp = remd%dparameters(i,j)

            write(MsgOut,'(A,I4)') '    Window index = ', j

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambljA = ', remd%dlambljA(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambljB = ', remd%dlambljB(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambelA = ', remd%dlambelA(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambelB = ', remd%dlambelB(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambbondA = ', remd%dlambbondA(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambbondB = ', remd%dlambbondB(lambid)

            write(MsgOut,'(6X,A,F10.5,5X)') &
              ' lambrest = ', remd%dlambrest(lambid)

            if (remd%types(i) == RemdAlchemyRest) then
              write(MsgOut,'(6X,A,F10.5,5X)') &
                ' solute_temp = ', soltemp
            end if

            write(MsgOut,'(A)') ''
          end do
        end if
      end do

    end if

    return

  end subroutine setup_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_pio
  !> @brief        setup REMD
  !! @authors      TM, JJ
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    dynamics   : dynamics information
  !! @param[inout] domain     : domain information
  !! @param[inout] restraints : restraints information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_pio(rep_info, boundary, dynamics, domain, &
                            restraints, ensemble, enefunc, remd)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_domain),          intent(inout) :: domain
    type(s_restraints),      intent(inout) :: restraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                    :: i, j, k, m
    integer                    :: max_nreplicas, cycles
    integer                    :: itmp, jtmp, ktmp
    integer                    :: comp, dsta, dend, found, dtotal
    integer                    :: replicaid, parmsetid
    integer                    :: ncycle, icycle, nlen, ixx
    character(MaxLine)         :: param
    character(10)              :: tmp, frmt, cd
    character(18)              :: frmt2
    character(3)               :: ctrep
    logical                    :: keep, double
    real(wp), allocatable      :: ddata(:)
    integer,       allocatable :: ndigit(:), parmid(:)
    integer                    :: funcid
    character(MaxLine)         :: ft


    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Remd> Replica information'
      write(MsgOut,'(A)') ''
    end if

    ! initialize
    !
    remd%dimension       = rep_info%dimension
    remd%total_nreplicas = product(rep_info%nreplicas(1:remd%dimension))
    max_nreplicas        = maxval (rep_info%nreplicas(1:remd%dimension))
    remd%analysis_grest  = rep_info%analysis_grest

    call alloc_remd(remd, RemdReplicas, remd%dimension,  &
                    remd%total_nreplicas, max_nreplicas)

    do i = 1, remd%dimension
      remd%nreplicas(i)     = rep_info%nreplicas(i)
      remd%types(i)         = rep_info%types(i)
      remd%cyclic_params(i) = rep_info%cyclic_params(i)
      remd%autoadj(i)%param_tuning    = rep_info%param_tuning(i)
      remd%autoadj(i)%tgt_exc_prob    = rep_info%tgt_exc_prob(i)
      remd%autoadj(i)%mgn_exc_prob    = rep_info%mgn_exc_prob(i)
      remd%autoadj(i)%trial_freq      = rep_info%trial_freq(i)
      remd%autoadj(i)%eq_cycle        = rep_info%eq_cycle(i)
      remd%autoadj(i)%param_grid      = rep_info%param_grid(i)
      remd%autoadj(i)%max_param_shift = rep_info%max_param_shift(i)
      ft = rep_info%fix_terminal(i)
      call toupper(ft)
      select case (ft)
        case ("BOTTOM")
          remd%autoadj(i)%fix_terminal = RemdAutoFixBottom
        case ("TOP")
          remd%autoadj(i)%fix_terminal = RemdAutoFixTop
        case default
          call error_msg('Error: fix_terminal must be BOTTOM/TOP.')
      end select
    end do

    ! check ensemble 
    do i = 1, remd%dimension
      if (remd%types(i) == RemdGamma .and.  &
          ensemble%ensemble /= EnsembleNPgT) &
         call error_msg('Setup_Remd> Ensemble should be NPgT in gamma-REMD')
    end do

    ! check cyclic_params
    !
    do i = 1, remd%dimension
      if (remd%cyclic_params(i)) then
        if (mod(remd%nreplicas(i),2) /= 0) &
          call error_msg('Setup_Remd> nreplicas must be even,'//&
                         ' if cyclic_params is specified')
      end if
    end do

    ! setup remd_iseed, num_criteria, and num_exchange
    !
    if (.not.allocated(remd%num_criteria)) then
      call alloc_remd(remd, RemdReplicas_rst, remd%dimension,  &
                      remd%total_nreplicas)
      remd%iseed = rep_info%iseed
      do i = 1, remd%total_nreplicas
        remd%repid2parmsetid(i) = i
        do j = 1, remd%dimension
          remd%num_criteria (i,j,1:2) = 0
          remd%num_exchanges(i,j,1:2) = 0
        end do
      end do
    end if

    ! setup permutation function
    !
    do i = 1, remd%total_nreplicas
      parmsetid = remd%repid2parmsetid(i)
      remd%parmsetid2repid(parmsetid) = i
    end do


    ! setup replica exchange period
    !
    if (rep_info%exchange_period > 0) then
      if (dynamics%eneout_period > rep_info%exchange_period) then
        call error_msg('Setup_Remd> (eneout_period <= exchange_period)'//&
                       ' is required')
      else
        if (mod(rep_info%exchange_period,dynamics%eneout_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,eneout_period)'//&
                         ' must be zero)')
        end if
        if (mod(rep_info%exchange_period,dynamics%thermo_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,thermo_period)'//&
                         ' must be zero)')
        end if
        if (mod(rep_info%exchange_period,dynamics%baro_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,baro_period)'//&
                         ' must be zero)')
        end if
      end if
      cycles = 2*rep_info%exchange_period*remd%dimension
      if (mod(dynamics%nsteps,cycles) /= 0) then
        call error_msg('Setup_Remd> mod(nsteps,(2*exchange_period*dimension))'//&
                       ' must be zero')
      else
        remd%equilibration_only = .false.
        remd%exchange_period = rep_info%exchange_period
        remd%ncycles = dynamics%nsteps/remd%exchange_period
      end if
    else if (rep_info%exchange_period == 0) then
      remd%equilibration_only = .true.
      remd%exchange_period = dynamics%nsteps
      remd%ncycles = 1
    else
      call error_msg('Setup_Remd> error in exchange_period')
    end if

    ! check parameter tuning
    !
    if (any(remd%autoadj(:)%param_tuning)) then
      if (remd%exchange_period <= 0) then
        remd%autoadj(:)%param_tuning = .false.
        if (main_rank) then
          write(MsgOut,'("Warning> param_tuning is disabled &
                          &since exchange is disabled.")')
          write(MsgOut,*)
        end if
      else
        ! clear all the exchange probs
        do i = 1, remd%total_nreplicas
          do j = 1, remd%dimension
            remd%num_criteria (i,j,1:2) = 0
            remd%num_exchanges(i,j,1:2) = 0
          end do
        end do
      end if
    end if


    ! split rep_info%parameters (char) into parameters (real)
    !
    allocate(ddata(max_nreplicas))

    do i = 1, remd%dimension
      if (remd%types(i) /= RemdRestraint) then
        param = rep_info%parameters(i)
        if (split_num(param) /= remd%nreplicas(i)) &
          call error_msg('Setup_Remd> "nreplicas" and number of data in "parameters" are inconsistent')

        call split(remd%nreplicas(i), max_nreplicas, param, ddata)

        do j = 1, remd%nreplicas(i)
          remd%dparameters(i,j) = ddata(j)
        end do
      end if
    end do

    deallocate(ddata)


    ! setup replica-exchange umbrella sampling (REUS)
    !
!   call setup_reus(rep_info, rst, restraints, enefunc, remd)


    ! setup replica exchange with solute tempering (REST)
    !
    call setup_remd_solute_tempering_pio(rep_info, domain, remd, &
                                         restraints, enefunc)

    do i = 1, remd%dimension
      if (main_rank .and. remd%autoadj(i)%param_tuning) then
        write(MsgOut,'(a,i8)') 'Setup_Remd> tuning parameters for dim ', i
        write(MsgOut,'(a,f8.3)') '  tgt_exc_prob    = ', &
                                 remd%autoadj(i)%tgt_exc_prob
        write(MsgOut,'(a,f8.3)') '  mgn_exc_prob    = ', &
                                 remd%autoadj(i)%mgn_exc_prob
        write(MsgOut,'(a,i8)')   '  trial_freq      = ', &
                                 remd%autoadj(i)%trial_freq
        write(MsgOut,'(a,i8)')   '  eq_cycle        = ', &
                                 remd%autoadj(i)%eq_cycle
        write(MsgOut,'(a,f8.3)') '  param_grid      = ', &
                                 remd%autoadj(i)%param_grid
        write(MsgOut,'(a,f8.3)') '  max_param_shift = ', &
                                 remd%autoadj(i)%max_param_shift
        if (remd%autoadj(i)%fix_terminal == RemdAutoFixBottom) then
          write(MsgOut,'(a)')    '  fix_terminal    = BOTTOM'
        else
          write(MsgOut,'(a)')    '  fix_terminal    = TOP'
        end if
        write(MsgOut,*)
      end if
    end do

    ! generate multi-dimensional parameter sets
    !
    allocate(ndigit(remd%dimension), parmid(remd%dimension))

    do i = 1, remd%dimension
      do j = 1, 100
        comp = 10**j
        if (remd%nreplicas(i) < comp) then
          ndigit(i) = j
          exit
        end if
      end do
    end do
    dtotal = sum(ndigit(1:remd%dimension))

    write(cd,'(i10)') dtotal
    frmt = '(i' // trim(adjustl(cd)) // '.' // trim(adjustl(cd)) // ')'

    found     = 0
    parmsetid = 0
    do
      found = found + 1
      write(tmp,frmt) found
      keep = .true.
      dend = 0

      do j = 1, remd%dimension
        dsta = dend + 1
        dend = dsta + ndigit(j) - 1
        read(tmp(dsta:dend),*) itmp
        if (itmp > remd%nreplicas(j) .or. itmp == 0) then
          keep = .false.
          exit
        end if
      end do

      if (keep) then
        parmsetid = parmsetid + 1
        dend = 0
        do j = 1, remd%dimension
          dsta = dend + 1
          dend = dsta + ndigit(j) - 1
          read(tmp(dsta:dend),*) parmid(j)
          remd%parmidsets(parmsetid,j) = parmid(j)
        end do
        if (parmsetid == remd%total_nreplicas) exit
      end if
    end do

    deallocate(ndigit,parmid)


    ! assign itnitial parameters
    !
    replicaid = my_country_no + 1
    parmsetid = remd%repid2parmsetid(replicaid)
    call assign_condition(parmsetid, remd, domain, ensemble, enefunc)


    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') '  ParmsetID'

      write(ctrep,'(I3)') remd%total_nreplicas
      frmt2 = '(I11 , A,' // ctrep // 'I11  )'

      do i = 1, remd%total_nreplicas
        replicaid = i
        parmsetid = remd%repid2parmsetid(replicaid)
        write(MsgOut,frmt2) parmsetid, ' = ', &
          remd%parmidsets(parmsetid,1:remd%dimension)
      end do

      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  Parameters'

      do i = 1, remd%dimension
        if (remd%types(i) /= RemdRestraint) then
          write(ctrep,'(I3)') remd%nreplicas(i)
          frmt2 = '(A,I4,A,' // ctrep // 'F11.3)'
          write(MsgOut,frmt2) &
            '    Dim', i, ' = ', remd%dparameters(i,1:remd%nreplicas(i))
        else
          write(ctrep,'(I3)') remd%nreplicas(i)
          frmt2 = '(A,I4,A,' // ctrep // 'I11  )'
          write(MsgOut,frmt2) &
            '    Dim', i, ' = ', remd%iparameters(i,1:remd%nreplicas(i))
        end if
      end do

      do i = 1, remd%dimension
        if (remd%types(i) == RemdRestraint) then
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '  Restraints'
          exit
        end if
      end do

      do i = 1, size(remd%umbrid2numfuncs(:))
        write(MsgOut,'(4X,A,I10)') 'umbrella potential index = ', i
!         do k = 1, remd%umbrid2numfuncs(i)
!           write(MsgOut,'(6X,A,I7,2X,A,F11.3,2X,A,F11.3)')         &
!                        'function = ',  remd%umbrid2funclist(i,k), &
!                        'constant = ',  remd%rest_constants(i,k),  &
!                        'reference = ', remd%rest_reference(i,k)
!         end do
        write(MsgOut,'(A)') ''
      end do

    end if

    return

  end subroutine setup_remd_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_reus
  !> @brief        setup REUS
  !! @authors      TM
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_reus(rep_info, rst, restraints, enefunc, remd)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_rst),             intent(in)    :: rst
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                    :: i, j, k, ifound
    integer                    :: ndata, max_ndata
    integer                    :: umbrid, funcid, max_umbrella_id, max_nreplicas
    character(MaxLine)         :: param
    character(10)              :: msg1, msg2
    integer,       allocatable :: idata(:)
    real(wp),      allocatable :: ddata(:)


    ! count numbers for allocation
    !
    max_nreplicas = maxval(rep_info%nreplicas(1:remd%dimension))
    max_ndata = 0
    ifound = 0
    do i = 1, remd%dimension
      if (remd%types(i) == RemdRestraint) then
        param = rep_info%rest_function(i)
        ndata = split_num(param)
        if (ndata > max_ndata) max_ndata = ndata
        ifound = ifound + remd%nreplicas(i)
      end if
    end do
    max_umbrella_id = ifound

    call alloc_remd(remd, RemdUmbrellas, max_umbrella_id, max_ndata, max_nreplicas)


    ! split rep_info%rest_function (char) into parameters (integer)
    !   and make permutation function for REUS
    !
    allocate(idata(max_ndata))

    umbrid = 0
    do i = 1, remd%dimension
      if (remd%types(i) == RemdRestraint) then
        param = rep_info%rest_function(i)
        ndata = split_num(param)

        call split(ndata, max_ndata, param, idata)

        do j = 1, remd%nreplicas(i)
          umbrid = umbrid + 1
          remd%iparameters(i,j) = umbrid
          remd%umbrid2numfuncs(umbrid) = ndata
          remd%umbrid2funclist(umbrid,1:ndata) = idata(1:ndata)
        end do
      end if
    end do

    deallocate(idata)


    ! check consistency between [REMD] and [RESTRAINTS]
    !
    umbrid = 0
    do i = 1, remd%dimension
      if (remd%types(i) == RemdRestraint) then
        do j = 1, remd%nreplicas(i)
          umbrid = umbrid + 1
          do k = 1, remd%umbrid2numfuncs(umbrid)
            funcid = remd%umbrid2funclist(umbrid,k)

            write(msg1,'(i0)') i
            write(msg2,'(i0)') funcid

            if (enefunc%restraint_kind(funcid) == RestraintsFuncPOSI) then
              call error_msg('Setup_Remd> REUS of positional restraints is not allowed')
            end if

            if (enefunc%restraint_kind(funcid) == RestraintsFuncPC .or.   &
                enefunc%restraint_kind(funcid) == RestraintsFuncPCCOM) then
              call error_msg('Setup_Remd> REUS of PC restraints is not allowed')
            end if

            if (split_num(restraints%constant(funcid)) /= remd%nreplicas(i)) then
              call error_msg('Setup_Remd> nreplica' // trim(msg1) // ' in [REMD] and &
                                          constant' // trim(msg2) // ' in [RESTRAINTS] &
                                          have inconsistency')
            end if

            if (enefunc%restraint_kind(funcid) /= RestraintsFuncEM) then
              if (split_num(restraints%reference(funcid)) /= remd%nreplicas(i)) then
                call error_msg('Setup_Remd> nreplica' // trim(msg1) // ' in [REMD] and &
                                            reference'// trim(msg2) // ' in [RESTRAINTS] &
                                            have inconsistency')
              end if
            end if

            if (funcid > restraints%nfunctions .or. funcid <= 0) then
              call error_msg('Setup_Remd> Cannot define index = ' // trim(msg2) // &
                                          ' in rest_function' // trim(msg1) // ' in [REMD]')
            end if

          end do
        end do
      end if
    end do


    ! setup constants and reference in umbrella potentials
    !
    allocate(ddata(max_nreplicas))

    umbrid = 0
    do i = 1, remd%dimension
      if (remd%types(i) == RemdRestraint) then
        do j = 1, remd%nreplicas(i)
          umbrid = umbrid + 1
          do k = 1, remd%umbrid2numfuncs(umbrid)
            funcid = remd%umbrid2funclist(umbrid,k)

            ! force constant
            param = restraints%constant(funcid)
            ndata = split_num(param)
            call split(ndata, max_nreplicas, param, ddata)
            remd%rest_constants(umbrid,k) = ddata(j)

            ! reference value
            param = restraints%reference(funcid)
            ndata = split_num(param)
            call split(ndata, max_nreplicas, param, ddata)
            remd%rest_reference(umbrid,k) = ddata(j)

          end do
        end do
      end if
    end do

    deallocate(ddata)

    ! replace rest_reference values with ones in rst (for REUS after RPATH)
    !
    if ((remd%dimension==1) .and. &
        (remd%types(1) == RemdRestraint) .and. &
        allocated(rst%rest_reference)) then
      do j = 1, remd%nreplicas(1)
        do k = 1, remd%umbrid2numfuncs(j)
          remd%rest_reference(j, k) = rst%rest_reference(1, k, j)
        end do
      end do

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Reus> REUS information'
        write(MsgOut,'(A)') '  References are replaced.'
        write(MsgOut, *)
        do k = 1, remd%umbrid2numfuncs(1)
          write(MsgOut,'(a,i0,a,$)')'  constants', k, ' = '
          do j = 1, remd%nreplicas(1)
            write(MsgOut, '(F10.4,$)') remd%rest_constants(j, k)
          end do
          write(MsgOut, *)
          write(MsgOut,'(a,i0,a,$)')'  reference', k, ' = '
          do j = 1, remd%nreplicas(1)
            write(MsgOut, '(F10.4,$)') remd%rest_reference(j, k)
          end do
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''
        end do
        write(MsgOut,'(A)') ''
      end if
    end if

    return

  end subroutine setup_reus

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering
  !> @brief        main routine for setuping REST
  !! @authors      MK
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] domain     : domain information
  !! @param[inout] remd       : REMD information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering(rep_info, molecule, domain, &
                                         remd, restraints, enefunc)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_remd),            intent(inout) :: remd
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer :: i, j, k, ndata, nrest
    integer :: ierror

    integer, parameter :: max_param_type = 20
    character(MaxLine) :: param_type_str(max_param_type)


    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then

        ! read paramter type
        ndata = split_num(rep_info%param_type(i))
        call split(ndata,max_param_type,rep_info%param_type(i), &
                   param_type_str)

        remd%solute_tempering(i)%sw_bonds        = .false.
        remd%solute_tempering(i)%sw_angles       = .false.
        remd%solute_tempering(i)%sw_ureys        = .false.
        remd%solute_tempering(i)%sw_dihedrals    = .false.
        remd%solute_tempering(i)%sw_rb_dihedrals = .false.
        remd%solute_tempering(i)%sw_impropers    = .false.
        remd%solute_tempering(i)%sw_cmaps        = .false.
        remd%solute_tempering(i)%sw_contacts     = .false.
        remd%solute_tempering(i)%sw_lj           = .false.
        remd%solute_tempering(i)%sw_charge       = .false.

        do j = 1, ndata
          call tolower(param_type_str(j))
          select case (trim(param_type_str(j)))
            case("all")
              remd%solute_tempering(i)%sw_bonds        = .true.
              remd%solute_tempering(i)%sw_angles       = .true.
              remd%solute_tempering(i)%sw_ureys        = .true.
              remd%solute_tempering(i)%sw_dihedrals    = .true.
              remd%solute_tempering(i)%sw_rb_dihedrals = .true.
              remd%solute_tempering(i)%sw_impropers    = .true.
              remd%solute_tempering(i)%sw_cmaps        = .true.
              remd%solute_tempering(i)%sw_contacts     = .true.
              remd%solute_tempering(i)%sw_lj           = .true.
              remd%solute_tempering(i)%sw_charge       = .true.
            case("b", "bond", "bonds")
              remd%solute_tempering(i)%sw_bonds        = .true.
            case("a", "angle", "angles")
              remd%solute_tempering(i)%sw_angles       = .true.
            case("u", "urey", "ureys")
              remd%solute_tempering(i)%sw_ureys        = .true.
            case("d", "dihedral", "dihedrals")
              remd%solute_tempering(i)%sw_dihedrals    = .true.
            case("rd", "rbd", "rb_d", "rb_dihedral", "rb_dihedrals")
              remd%solute_tempering(i)%sw_rb_dihedrals = .true.
            case("i", "improper", "impropers")
              remd%solute_tempering(i)%sw_impropers    = .true.
            case("cm", "cmap", "cmaps")
              remd%solute_tempering(i)%sw_cmaps        = .true.
            case("con", "contact", "contacts")
              remd%solute_tempering(i)%sw_contacts     = .true.
            case("c", "charge", "charges")
              remd%solute_tempering(i)%sw_charge       = .true.
            case("l", "lj", "ljs")
              remd%solute_tempering(i)%sw_lj           = .true.
          end select
        end do

        ! this value is set to true after the first call of assign_condition
        !
        remd%solute_tempering(i)%done_setup = .false.

        remd%solute_tempering(i)%num_solute = 0
        remd%solute_tempering(i)%num_atom_cls = 0

        read(rep_info%select_index(i),*) remd%solute_tempering(i)%mygroup

        nrest = restraints%num_atoms(remd%solute_tempering(i)%mygroup)
        remd%solute_tempering(i)%num_solute = nrest

        allocate(remd%solute_tempering(i)%solute_list(nrest), &
                 remd%solute_tempering(i)%is_solute(enefunc%table%num_all), &
                 remd%solute_tempering(i)%atom_cls_no_org(nrest))
        remd%solute_tempering(i)%is_solute(1:enefunc%table%num_all) = 0

        ! listup solute atoms
        !
        do j = 1, restraints%num_atoms(remd%solute_tempering(i)%mygroup)
          k = restraints%atomlist(j,remd%solute_tempering(i)%mygroup)
          remd%solute_tempering(i)%solute_list(j) = k
          remd%solute_tempering(i)%is_solute(k) = 1
        end do

        if (main_rank) then

          write(MsgOut,'(a)')    'Setup_Remd_Solute_Tempering>'
          write(MsgOut,'(a,i8)') 'Param types for dimension', i
          write(MsgOut,'(a,l8)') '  BONDS         = ', &
                                 remd%solute_tempering(i)%sw_bonds
          write(MsgOut,'(a,l8)') '  ANGLES        = ', &
                                 remd%solute_tempering(i)%sw_angles
          write(MsgOut,'(a,l8)') '  UREYS         = ', &
                                 remd%solute_tempering(i)%sw_ureys
          write(MsgOut,'(a,l8)') '  DIHEDRALS     = ', &
                                 remd%solute_tempering(i)%sw_dihedrals
          write(MsgOut,'(a,l8)') '  RB_DIHEDRALS  = ', &
                                 remd%solute_tempering(i)%sw_rb_dihedrals
          write(MsgOut,'(a,l8)') '  IMPROPERS     = ', &
                                 remd%solute_tempering(i)%sw_impropers
          write(MsgOut,'(a,l8)') '  CMAPS         = ', &
                                 remd%solute_tempering(i)%sw_cmaps
          write(MsgOut,'(a,l8)') '  CONTACTS      = ', &
                                 remd%solute_tempering(i)%sw_contacts
          write(MsgOut,'(a,l8)') '  CHARGE        = ', &
                                 remd%solute_tempering(i)%sw_charge
          write(MsgOut,'(a,l8)') '  LJ            = ', &
                                 remd%solute_tempering(i)%sw_lj
        end if

        if (main_rank) then
          write(MsgOut,'(a,i8)') 'Solute Atoms for dimension', i
          write(MsgOut,'(a,i0)') '  num_solute              = ', nrest

        end if

        ! bond (just count number)
        !
        call setup_remd_solute_tempering_bonds(remd%solute_tempering(i), &
                                               domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_bonds, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_bonds
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_bond         = ', nrest
        end if

        ! angle (just count number)
        !
        call setup_remd_solute_tempering_angles(remd%solute_tempering(i), &
                                                domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_angles, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_angles
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_angles       = ', nrest
        end if
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_ureys, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_ureys
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_ureys        = ', nrest
        end if

        ! dihedral (just count number)
        !
        call setup_remd_solute_tempering_dihedrals(remd%solute_tempering(i), &
                                                   domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_dihedrals
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_dihedrals    = ', nrest
        end if

        ! rb_dihedral (just count number)
        !
        call setup_remd_solute_tempering_rb_dihedrals( &
                        remd%solute_tempering(i), domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_rb_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_rb_dihedrals
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_rb_dihedrals = ', nrest
        end if

        ! improper (just count number)
        !
        call setup_remd_solute_tempering_impropers(remd%solute_tempering(i), &
                                                   domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_impropers, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_impropers
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_impropers    = ', nrest
        end if

        ! cmap (just count number)
        !
        call setup_remd_solute_tempering_cmaps(remd%solute_tempering(i), &
                                               domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_cmaps, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_cmaps
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_cmaps        = ', nrest
        end if

        ! contact (just count number)
        !
        call setup_remd_solute_tempering_contacts(remd%solute_tempering(i), &
                                                  domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_contacts, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_contacts
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_contacts     = ', nrest
        end if

        ! LJ
        !
        call setup_remd_solute_tempering_lj(remd%solute_tempering(i), &
                                            domain, enefunc)

        ! charges are not modified here; to be done in assign_condition

      end if
    end do

    return

  end subroutine setup_remd_solute_tempering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_pio
  !> @brief        main routine for setuping REST
  !! @authors      MK, JJ
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[inout] domain     : domain information
  !! @param[inout] remd       : REMD information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_pio(rep_info, domain, remd, &
                                             restraints, enefunc)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_domain),          intent(inout) :: domain
    type(s_remd),            intent(inout) :: remd
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer :: i, j, k, ndata, nrest
    integer :: ierror

    integer, parameter :: max_param_type = 20
    character(MaxLine) :: param_type_str(max_param_type)


    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then

        ! read paramter type
        ndata = split_num(rep_info%param_type(i))
        call split(ndata,max_param_type,rep_info%param_type(i), &
                   param_type_str)

        remd%solute_tempering(i)%sw_bonds        = .false.
        remd%solute_tempering(i)%sw_angles       = .false.
        remd%solute_tempering(i)%sw_ureys        = .false.
        remd%solute_tempering(i)%sw_dihedrals    = .false.
        remd%solute_tempering(i)%sw_rb_dihedrals = .false.
        remd%solute_tempering(i)%sw_impropers    = .false.
        remd%solute_tempering(i)%sw_cmaps        = .false.
        remd%solute_tempering(i)%sw_contacts     = .false.
        remd%solute_tempering(i)%sw_lj           = .false.
        remd%solute_tempering(i)%sw_charge       = .false.

        do j = 1, ndata
          call tolower(param_type_str(j))
          select case (trim(param_type_str(j)))
            case("all")
              remd%solute_tempering(i)%sw_bonds        = .true.
              remd%solute_tempering(i)%sw_angles       = .true.
              remd%solute_tempering(i)%sw_ureys        = .true.
              remd%solute_tempering(i)%sw_dihedrals    = .true.
              remd%solute_tempering(i)%sw_rb_dihedrals = .true.
              remd%solute_tempering(i)%sw_impropers    = .true.
              remd%solute_tempering(i)%sw_cmaps        = .true.
              remd%solute_tempering(i)%sw_contacts     = .true.
              remd%solute_tempering(i)%sw_lj           = .true.
              remd%solute_tempering(i)%sw_charge       = .true.
            case("b", "bond", "bonds")
              remd%solute_tempering(i)%sw_bonds        = .true.
            case("a", "angle", "angles")
              remd%solute_tempering(i)%sw_angles       = .true.
            case("u", "urey", "ureys")
              remd%solute_tempering(i)%sw_ureys        = .true.
            case("d", "dihedral", "dihedrals")
              remd%solute_tempering(i)%sw_dihedrals    = .true.
            case("rd", "rbd", "rb_d", "rb_dihedral", "rb_dihedrals")
              remd%solute_tempering(i)%sw_rb_dihedrals = .true.
            case("i", "improper", "impropers")
              remd%solute_tempering(i)%sw_impropers    = .true.
            case("cm", "cmap", "cmaps")
              remd%solute_tempering(i)%sw_cmaps        = .true.
            case("con", "contact", "contacts")
              remd%solute_tempering(i)%sw_contacts     = .true.
            case("c", "charge", "charges")
              remd%solute_tempering(i)%sw_charge       = .true.
            case("l", "lj", "ljs")
              remd%solute_tempering(i)%sw_lj           = .true.
          end select
        end do

        ! this value is set to true after the first call of assign_condition
        !
        remd%solute_tempering(i)%done_setup = .false.

        remd%solute_tempering(i)%num_solute = 0
        remd%solute_tempering(i)%num_atom_cls = 0

        read(rep_info%select_index(i),*) remd%solute_tempering(i)%mygroup

        nrest = restraints%num_atoms(remd%solute_tempering(i)%mygroup)
        remd%solute_tempering(i)%num_solute = nrest

        allocate(remd%solute_tempering(i)%solute_list(nrest), &
                 remd%solute_tempering(i)%is_solute(enefunc%table%num_all), &
                 remd%solute_tempering(i)%atom_cls_no_org(nrest))
        remd%solute_tempering(i)%is_solute(1:enefunc%table%num_all) = 0

        ! listup solute atoms
        !
        do j = 1, restraints%num_atoms(remd%solute_tempering(i)%mygroup)
          k = restraints%atomlist(j,remd%solute_tempering(i)%mygroup)
          remd%solute_tempering(i)%solute_list(j) = k
          remd%solute_tempering(i)%is_solute(k) = 1
        end do

        if (main_rank) then

          write(MsgOut,'(a)')    'Setup_Remd_Solute_Tempering>'
          write(MsgOut,'(a,i8)') 'Param types for dimension', i
          write(MsgOut,'(a,l8)') '  BONDS         = ', &
                                 remd%solute_tempering(i)%sw_bonds
          write(MsgOut,'(a,l8)') '  ANGLES        = ', &
                                 remd%solute_tempering(i)%sw_angles
          write(MsgOut,'(a,l8)') '  UREYS         = ', &
                                 remd%solute_tempering(i)%sw_ureys
          write(MsgOut,'(a,l8)') '  DIHEDRALS     = ', &
                                 remd%solute_tempering(i)%sw_dihedrals
          write(MsgOut,'(a,l8)') '  RB_DIHEDRALS  = ', &
                                 remd%solute_tempering(i)%sw_rb_dihedrals
          write(MsgOut,'(a,l8)') '  IMPROPERS     = ', &
                                 remd%solute_tempering(i)%sw_impropers
          write(MsgOut,'(a,l8)') '  CMAPS         = ', &
                                 remd%solute_tempering(i)%sw_cmaps
          write(MsgOut,'(a,l8)') '  CONTACTS      = ', &
                                 remd%solute_tempering(i)%sw_contacts
          write(MsgOut,'(a,l8)') '  CHARGE        = ', &
                                 remd%solute_tempering(i)%sw_charge
          write(MsgOut,'(a,l8)') '  LJ            = ', &
                                 remd%solute_tempering(i)%sw_lj
        end if

        if (main_rank) then
          write(MsgOut,'(a,i8)') 'Solute Atoms for dimension', i
          write(MsgOut,'(a,i0)') '  num_solute              = ', nrest

        end if

        ! bond (just count number)
        !
        call setup_remd_solute_tempering_bonds(remd%solute_tempering(i), &
                                               domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_bonds, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_bonds
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_bond         = ', nrest
        end if

        ! angle (just count number)
        !
        call setup_remd_solute_tempering_angles(remd%solute_tempering(i), &
                                                domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_angles, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_angles
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_angles       = ', nrest
        end if
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_ureys, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_ureys
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_ureys        = ', nrest
        end if

        ! dihedral (just count number)
        !
        call setup_remd_solute_tempering_dihedrals(remd%solute_tempering(i), &
                                                   domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_dihedrals
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_dihedrals    = ', nrest
        end if

        ! rb_dihedral (just count number)
        !
        call setup_remd_solute_tempering_rb_dihedrals( &
                        remd%solute_tempering(i), domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_rb_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_rb_dihedrals
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_rb_dihedrals = ', nrest
        end if

        ! improper (just count number)
        !
        call setup_remd_solute_tempering_impropers(remd%solute_tempering(i), &
                                                   domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_impropers, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_impropers
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_impropers    = ', nrest
        end if

        ! cmap (just count number)
        !
        call setup_remd_solute_tempering_cmaps(remd%solute_tempering(i), &
                                               domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_cmaps, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_cmaps
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_cmaps        = ', nrest
        end if

        ! contact (just count number)
        !
        call setup_remd_solute_tempering_contacts(remd%solute_tempering(i), &
                                                  domain, enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_contacts, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_contacts
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_contacts     = ', nrest
        end if

        ! LJ
        !
        call setup_remd_solute_tempering_lj(remd%solute_tempering(i), &
                                            domain, enefunc)

        ! charges are not modified here; to be done in assign_condition

      end if
    end do

    return

  end subroutine setup_remd_solute_tempering_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_bonds
  !> @brief        REST setup for bonds
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_bonds(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 2
    integer            :: alist(num_unit), i, ix, j, num


    soltemp%num_bonds = 0

    if (.not. soltemp%sw_bonds) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_bond(i)
        alist(1:num_unit) = enefunc%bond_list(1:2,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_bonds = soltemp%num_bonds + 1
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_angles
  !> @brief        REST setup for angles and UBs
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_angles(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 3
    integer            :: alist(num_unit), i, ix, j, num


    soltemp%num_angles = 0
    soltemp%num_ureys = 0

    if (.not. soltemp%sw_angles .and. .not. soltemp%sw_ureys) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_angle(i)
        alist(1:num_unit) = enefunc%angle_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0 .and. soltemp%sw_angles) then
          soltemp%num_angles = soltemp%num_angles + 1
        end if

        if (enefunc%urey_force_const(ix,i) > EPS .and. &
             soltemp%sw_ureys) then
          num = 0
          if (soltemp%is_solute(alist(1)) > 0) num = num + 1
          if (soltemp%is_solute(alist(3)) > 0) num = num + 1

          if (num > 0) then
            soltemp%num_ureys = soltemp%num_ureys + 1
          end if
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_angles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_dihedrals
  !> @brief        REST setup for dihedrals
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_dihedrals(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 4
    integer            :: alist(num_unit), i, ix, j, num

    soltemp%num_dihedrals = 0

    if (.not. soltemp%sw_dihedrals) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_dihedral(i)
        alist(1:num_unit) = enefunc%dihe_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_dihedrals = soltemp%num_dihedrals + 1
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_dihedrals

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_rb_dihedrals
  !> @brief        REST setup for rb_dihedrals
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_rb_dihedrals(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 4
    integer            :: alist(num_unit), i, ix, j, num


    soltemp%num_rb_dihedrals = 0

    if (.not. soltemp%sw_rb_dihedrals) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_rb_dihedral(i)
        alist(1:num_unit) = enefunc%rb_dihe_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_rb_dihedrals = soltemp%num_rb_dihedrals + 1
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_rb_dihedrals

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_impropers
  !> @brief        REST setup for impropers
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_impropers(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 4
    integer            :: alist(num_unit), i, ix, j, num


    soltemp%num_impropers = 0

    if (.not. soltemp%sw_impropers) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_improper(i)
        alist(1:num_unit) = enefunc%impr_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_impropers = soltemp%num_impropers + 1
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_impropers

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_cmaps
  !> @brief        REST setup for cmaps
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_cmaps(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 8
    integer            :: alist(num_unit), i, ix, j, k, num

    integer :: num_new_types, ct, ncmap_type, ierror
    integer,  allocatable :: flags(:,:), gflags(:,:)
    integer,  allocatable :: cmap_resolution_org(:)
    real(wp), allocatable :: cmap_coef_org(:,:,:,:,:)


    soltemp%num_cmaps = 0
    soltemp%num_cmap_type = 0

    if (.not. allocated(enefunc%cmap_coef) .or. &
        .not. soltemp%sw_cmaps) return

    ncmap_type = size(enefunc%cmap_coef(1,1,1,1,:))
    if (ncmap_type == 0) return

    ! need to expand cmap_coef
    allocate(flags(ncmap_type,num_unit), gflags(ncmap_type,num_unit))
    flags(1:ncmap_type,1:num_unit) = 0
    gflags(1:ncmap_type,1:num_unit) = 0

    num_new_types = 0
    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_cmap(i)
        alist(1:num_unit) = enefunc%cmap_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_cmaps = soltemp%num_cmaps + 1
          ct = enefunc%cmap_type(ix,i)
          if (flags(ct,num) == 0) then
            num_new_types = num_new_types + 1
            flags(ct,num) = num_new_types
          end if
        end if

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(flags, gflags, ncmap_type*num_unit, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    num_new_types = 0
    do i = 1, num_unit
      do j = 1, ncmap_type
        if (gflags(j,i) > 0) then
          num_new_types = num_new_types + 1
          flags(j,i) = num_new_types
        end if
      end do
    end do
#endif

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_cmap(i)
        alist(1:num_unit) = enefunc%cmap_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          ct = enefunc%cmap_type(ix,i)
          enefunc%cmap_type(ix,i) = ncmap_type + flags(ct,num)
        end if

      end do
    end do

    if (num_new_types  > 0) then
      allocate(cmap_coef_org(4,4,24,24,ncmap_type), &
               cmap_resolution_org(ncmap_type))
      allocate(soltemp%cmap_weight(num_new_types), &
               soltemp%cmap_type_org(num_new_types))

      cmap_coef_org(1:4,1:4,1:24,1:24,1:ncmap_type) = &
                enefunc%cmap_coef(1:4,1:4,1:24,1:24,1:ncmap_type)
      cmap_resolution_org(1:ncmap_type) = &
                enefunc%cmap_resolution(1:ncmap_type)

      soltemp%cmap_weight(:) = 0.0_wp

      deallocate(enefunc%cmap_coef, &
                 enefunc%cmap_resolution)

      ! expand array for new params
      !
      allocate(enefunc%cmap_coef(4,4,24,24,ncmap_type+num_new_types), &
               enefunc%cmap_resolution(ncmap_type+num_new_types))

      enefunc%cmap_coef(1:4,1:4,1:24,1:24,1:ncmap_type) = &
                cmap_coef_org(1:4,1:4,1:24,1:24,1:ncmap_type)
      enefunc%cmap_resolution(1:ncmap_type) = &
                cmap_resolution_org(1:ncmap_type)

      do k = 1, 8
        do j = 1, ncmap_type
          if (flags(j,k) > 0) then
            soltemp%cmap_weight(flags(j,k)) = real(k,wp) / real(num_unit,wp)
            soltemp%cmap_type_org(flags(j,k)) = j
            enefunc%cmap_resolution(ncmap_type+flags(j,k)) = &
                      enefunc%cmap_resolution(j)
          end if
        end do
      end do

      soltemp%num_cmap_type = num_new_types
      soltemp%istart_cmap_type = ncmap_type + 1

      deallocate(cmap_coef_org)
      deallocate(cmap_resolution_org)

    end if

    deallocate(flags, gflags)

    return

  end subroutine setup_remd_solute_tempering_cmaps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_contacts
  !> @brief        REST setup for contacts
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_contacts(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 2
    integer            :: alist(num_unit), i, ix, j, num


    soltemp%num_contacts = 0

    if (.not. soltemp%sw_contacts) return

    do i = 1, domain%num_cell_local
      do ix = 1, enefunc%num_contact(i)
        alist(1:num_unit) = enefunc%contact_list(1:num_unit,ix,i)

        num = 0
        do j = 1, num_unit
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do

        if (num > 0) then
          soltemp%num_contacts = soltemp%num_contacts + 1
        end if

      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_contacts

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_lj
  !> @brief        REST setup for lj terms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_lj(soltemp, domain, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! temporal allocatable arrays
    real(wp),                allocatable   :: tmp_lj(:,:)
    ! work around for nbfix extension
    logical,                 allocatable   :: check_cls(:)
    integer,                 allocatable   :: lc_atom_cls(:)

    ! local variables
    integer :: i, j, ix, localcount, oldcount, org, new, ncell
    integer :: newcount
#ifdef HAVE_MPI_GENESIS
    integer :: ierror
#endif


    soltemp%num_atom_cls = 0
    if (.not. soltemp%sw_lj) return

    !! check solute atom types
    allocate(check_cls(enefunc%num_atom_cls), &
              lc_atom_cls(enefunc%num_atom_cls))
    check_cls(1:enefunc%num_atom_cls) = .false.
    lc_atom_cls(1:enefunc%num_atom_cls) = -1 ! force to cause error

    do i = 1, domain%num_cell_local + domain%num_cell_boundary
      do ix = 1, domain%num_solute(i)
        if (soltemp%is_solute(domain%id_l2g_solute(ix,i)) > 0) then
          check_cls(domain%atom_cls_no(ix,i)) = .true.
        end if
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(MPI_IN_PLACE, check_cls, enefunc%num_atom_cls, &
                       MPI_LOGICAL, MPI_LOR, mpi_comm_country, ierror)
#endif

    localcount = 0
    do i = 1, enefunc%num_atom_cls
      if (check_cls(i)) then
        localcount = localcount + 1
        lc_atom_cls(i) = localcount
        soltemp%atom_cls_no_org(localcount) = i
      end if
    end do
    soltemp%num_atom_cls = localcount
    soltemp%istart_atom_cls = enefunc%num_atom_cls + 1

    ! loop over local cells; modify atom type
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    do i = 1, ncell
      do ix = 1, domain%num_solute(i)
        if (soltemp%is_solute(domain%id_l2g_solute(ix,i)) > 0) then
          org = lc_atom_cls(domain%atom_cls_no(ix,i))
          domain%atom_cls_no(ix,i) = org + enefunc%num_atom_cls
        end if
      end do
    end do

    oldcount = enefunc%num_atom_cls
    allocate(tmp_lj(oldcount,oldcount))
    enefunc%num_atom_cls = oldcount + localcount
    max_class = enefunc%num_atom_cls

    newcount = oldcount + soltemp%num_atom_cls
    if (allocated(enefunc%nb14_lj6)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nb14_lj6(1:oldcount,1:oldcount)
      deallocate(enefunc%nb14_lj6)
      allocate(enefunc%nb14_lj6(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nb14_lj6)
    end if
    if (allocated(enefunc%nb14_lj12)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nb14_lj12(1:oldcount,1:oldcount)
      deallocate(enefunc%nb14_lj12)
      allocate(enefunc%nb14_lj12(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nb14_lj12)
    end if
    if (allocated(enefunc%nonb_lj6)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nonb_lj6(1:oldcount,1:oldcount)
      deallocate(enefunc%nonb_lj6)
      allocate(enefunc%nonb_lj6(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nonb_lj6)
    end if
    if (allocated(enefunc%nonb_lj12)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nonb_lj12(1:oldcount,1:oldcount)
      deallocate(enefunc%nonb_lj12)
      allocate(enefunc%nonb_lj12(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nonb_lj12)
    end if

    return

  end subroutine setup_remd_solute_tempering_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_lj_each
  !> @brief        REST setup for each LJ term such as LJ6, LJ12
  !! @authors      MK
  !! @param[in]    oldcount   : old lj array count
  !! @@aram[in]    newcount   : new lj array count
  !! @param[in]    soltemp    : REST information
  !! @param[in]    tmp_lj     : work matrix
  !! @param[inout] nbcoeff    : new interaction matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_lj_each(oldcount, newcount, soltemp, &
                                                 tmp_lj, nbcoeff)

    ! formal arguments
    integer,                      intent(in)    :: oldcount
    integer,                      intent(in)    :: newcount
    type(s_soltemp),              intent(in)    :: soltemp
    real(wp),                     intent(in)    :: tmp_lj(oldcount,oldcount)
    real(wp),                     intent(inout) :: nbcoeff(:,:)

    ! local variables
    integer :: i, j, orgi, orgj


    nbcoeff(1:newcount,1:newcount) = 0.0_wp
    nbcoeff(1:oldcount,1:oldcount) = tmp_lj(1:oldcount,1:oldcount)

    do i = 1, soltemp%num_atom_cls
      orgi = soltemp%atom_cls_no_org(i)
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = tmp_lj(j,orgi)
        nbcoeff(i+oldcount,j) = tmp_lj(j,orgi)
      end do
      do j = oldcount + 1, newcount
        orgj = soltemp%atom_cls_no_org(j-oldcount)
        nbcoeff(j,i+oldcount) = tmp_lj(orgj,orgi)
        nbcoeff(i+oldcount,j) = tmp_lj(orgj,orgi)
      end do
    end do

    return

  end subroutine setup_remd_solute_tempering_lj_each

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_remd
  !> @brief        control replica exchange
  !! @authors      TM
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
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_remd(output, domain, enefunc, dynvars, dynamics, &
                      pairlist, boundary, constraints, ensemble,  &
                      comm, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd

    ! local variables
    integer                   :: i, replicaid, parmsetid


    ! Open output files
    !
    call open_output(output)

    ! output restart data
    !
    dynvars%step = 0
    call output_remd(0, output, enefunc, domain, dynamics, dynvars, boundary, &
                     constraints, remd)

    do i = 1, remd%ncycles

      dynamics%istart_step  = (i-1)*remd%exchange_period + 1
      dynamics%iend_step    =  i   *remd%exchange_period
      dynamics%initial_time = dynvars%time

      ! set conditions
      !
      replicaid = my_country_no + 1
      parmsetid = remd%repid2parmsetid(replicaid)
      call assign_condition(parmsetid, remd, domain, ensemble, enefunc)

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
        call vverlet_respa_dynamics(output, domain, enefunc, dynvars,      &
                               dynamics, pairlist, boundary, constraints,  &
                               ensemble, comm, remd)
      end if

      ! perform remd
      !
      if (.not. remd%equilibration_only) then
        call perform_replica_exchange(i, domain, enefunc, dynvars,         &
                                      ensemble, boundary, output, remd,    &
                                      pairlist)
      end if

      ! output restart data
      !
      call output_remd(i, output, enefunc, domain, dynamics, dynvars, &
                       boundary, constraints, remd)


    end do

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    perform_replica_exchange
  !> @brief        perform replica exchange
  !! @authors      TM
  !! @param[in]    icycle      : replica exchange cycle number
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] output      : output information
  !! @param[inout] remd        : REMD information
  !! @param[inout] pairlist    : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine perform_replica_exchange(icycle, domain, enefunc, dynvars, &
                                      ensemble, boundary, output, remd, &
                                      pairlist)

    ! formal arguments
    integer,                  intent(in)    :: icycle
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist

    ! local variables
    integer                   :: exptrn, dimno, itmp, parmsetid
    integer                   :: i, j, jx, k, neighid, id_old
    integer                   :: umbr_0, funcid
    integer                   :: ipara, ipara_ref
    real(wp)                  :: dpara, dpara_ref
    character                 :: ctrep*3, frmt1*12, frmt2*12, aorr*1
    integer,          pointer :: total_nreplicas, iseed
    integer,          pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,          pointer :: repid2parmsetid_ref(:)
    integer,          pointer :: iparameters(:,:)
    real(wp),         pointer :: dparameters(:,:)

    integer, save, allocatable :: counter_skip(:)
    integer, save, allocatable :: counter_tuning(:)

    real(wp) :: temp_i, temp_j, factor


    iseed               => remd%iseed
    total_nreplicas     => remd%total_nreplicas
    repid2parmsetid     => remd%repid2parmsetid
    repid2parmsetid_ref => remd%repid2parmsetid_ref
    parmsetid2repid     => remd%parmsetid2repid
    parmidsets          => remd%parmidsets
    dparameters         => remd%dparameters
    iparameters         => remd%iparameters
    num_criteria        => remd%num_criteria
    num_exchanges       => remd%num_exchanges


    ! decide dimension to be exchanged and exchange pattern
    !   pattern 1: 1<>2  3<>4  5<>6  7<>8 ...
    !   pattern 2: 1  2<>3  4<>5  6<>7  8 ...
    !
    itmp  = mod(icycle,2*remd%dimension)
    dimno = mod(itmp,remd%dimension)+1
    if (itmp+1 <= remd%dimension) then
      exptrn = 1
    else
      exptrn = 2
    end if

    ! parameter tuning
    !
    if (remd%autoadj(dimno)%param_tuning) then
      if (.not. allocated(counter_skip)) then
        allocate(counter_skip(remd%dimension))
        counter_skip(1:remd%dimension) = 0
      end if
      if (.not. allocated(counter_tuning)) then
        allocate(counter_tuning(remd%dimension))
        counter_tuning(1:remd%dimension) = 0
      end if
      counter_skip(dimno) = counter_skip(dimno) + 1
      if (counter_skip(dimno) <= 2 * remd%autoadj(dimno)%eq_cycle) then
        if (main_rank) then
          write(MsgOut,'(A, I12)') &
            "Info> skip exchange at step", dynvars%step
          write(MsgOut,*)
        end if
        return
      end if
    end if


    ! save previous information
    !
    repid2parmsetid_ref(1:total_nreplicas) = repid2parmsetid(1:total_nreplicas)


    ! generate random number
    !
    do i = 1, total_nreplicas-1
      do j = i+1, total_nreplicas
        remd%random_table(i,j) = random_get_legacy(iseed)
        remd%random_table(j,i) = remd%random_table(i,j)
      end do
    end do

    
    ! perform replica exchange
    !
    select case(remd%types(dimno))

    case(RemdTemperature)

      call temperature_remd(dimno, exptrn, domain, dynvars, remd)

    case(RemdPressure)

      call pressure_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdGamma)

      call surface_tension_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdRestraint)

      call restraint_remd(dimno, exptrn, domain, boundary, ensemble, dynvars, &
                          enefunc, remd)

    case(RemdSoluteTempering)

      call solute_tempering_remd(dimno, exptrn, domain, ensemble, dynvars, &
                                 enefunc, remd, pairlist, boundary, output)

    end select

    ! parameter tuning
    !
    if (remd%autoadj(dimno)%param_tuning) then
      counter_tuning(dimno) = counter_tuning(dimno) + 1
      if (counter_tuning(dimno) >= 2 * remd%autoadj(dimno)%trial_freq) then

        ! in case of temperature REMD, we should remember current temp
        if (remd%types(dimno) == RemdTemperature) then
          parmsetid = repid2parmsetid(my_country_no+1)
          temp_i = remd%dparameters(dimno,parmidsets(parmsetid,dimno))
        end if

        ! check current exchange ratio
        call remd_autoadj_displacement(dimno,remd)

        if (remd%types(dimno) == RemdTemperature) then
          temp_j = remd%dparameters(dimno,parmidsets(parmsetid,dimno))
          factor = sqrt(temp_j/temp_i)
          do j = 1, domain%num_cell_local
            do jx = 1, domain%num_atom(j)
              domain%velocity(1:3,jx,j) = domain%velocity(1:3,jx,j) * factor
            end do
          end do

          factor = temp_j/temp_i
          dynvars%thermostat_momentum                 &
            = dynvars%thermostat_momentum * factor
          dynvars%barostat_momentum(1:3)              &
            = dynvars%barostat_momentum(1:3) * factor

        end if

        ! reset counters
        do i = 1, remd%total_nreplicas
          remd%num_criteria (i,dimno,1:2) = 0
          remd%num_exchanges(i,dimno,1:2) = 0
        end do

        if (main_rank) then
          write(MsgOut,'("REMD> New parameter set: ")',advance="NO")
          if (remd%types(dimno) /= RemdRestraint) then
            do i = 1, remd%nreplicas(dimno)
              write(MsgOut,'(2X,F10.5)',advance="NO") remd%dparameters(dimno,i)
            end do
          else
            umbr_0 = remd%iparameters(dimno,parmidsets(1,dimno))
            do k = 1, remd%umbrid2numfuncs(umbr_0)
              funcid = remd%umbrid2funclist(umbr_0,k)
              do i = 1, remd%nreplicas(dimno)
                write(MsgOut,'(2X,F10.5)',advance="NO") remd%rest_reference(i,k)
              end do
            end do
          end if
          write(MsgOut,*)
          write(MsgOut,*)
        end if
        counter_skip(dimno)   = 0
        counter_tuning(dimno) = 0
      end if
    end if


    ! update permutation function
    !
    do i = 1, total_nreplicas
      parmsetid = repid2parmsetid(i)
      parmsetid2repid(parmsetid) = i
    end do


    ! write results
    !
    if (main_rank) then

      write(MsgOut,'(A,I10,3X,A,I4,3X,A,I4)')                           &
        'REMD> Step: ',dynvars%step,'Dimension: ',dimno,                &
        'ExchangePattern: ',exptrn

      write(MsgOut,'(A,6X,A,13X,A,6X,A,7X,A)')     &
        '  Replica', 'ExchangeTrial', 'AcceptanceRatio', 'Before', 'After'
      do i = 1, total_nreplicas
        id_old = repid2parmsetid_ref(i)
        call get_neighbouring_parmsetid(dimno, exptrn, remd, id_old, neighid)

        if (remd%types(dimno) /= RemdRestraint) then

          dpara_ref = dparameters(dimno,parmidsets(repid2parmsetid_ref(i),dimno))
          dpara     = dparameters(dimno,parmidsets(repid2parmsetid(i),    dimno))

          if (abs(dpara_ref - dpara) < EPS) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2F12.3)')         &
               i, id_old, ' > ', neighid, aorr,                           &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ', &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),        &
               dpara_ref, dpara
        else

          ipara_ref = iparameters(dimno,parmidsets(repid2parmsetid_ref(i),dimno))
          ipara     = iparameters(dimno,parmidsets(repid2parmsetid(i),    dimno))

          if (ipara_ref == ipara) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2I12)')           &
               i, id_old, ' > ', neighid, aorr,                           &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ', &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),        &
               ipara_ref, ipara
        end if
      end do

      write(MsgOut,*) ''

      write(ctrep,'(I3)') total_nreplicas
      frmt1 = '(A,' // trim(adjustl(ctrep)) // 'F10.3)'
      frmt2 = '(A,' // trim(adjustl(ctrep)) // 'I10  )'

      if (remd%types(dimno) /= RemdRestraint) then
        write(MsgOut,frmt1) '  Parameter    : ', &
          dparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),dimno))
      else
        write(MsgOut,frmt2) '  Parameter    : ', &
          iparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),dimno))
      end if
      write(MsgOut,frmt2)   '  RepIDtoParmID: ', repid2parmsetid(1:total_nreplicas)
      write(MsgOut,frmt2)   '  ParmIDtoRepID: ', parmsetid2repid(1:total_nreplicas)

      write(MsgOut,*) ''

    end if

    return

  end subroutine perform_replica_exchange

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_autoadj_displacement(dimno, remd)

    integer,                    intent(in)    :: dimno
    type(s_remd),               intent(inout) :: remd

    real(wp), allocatable, save :: disp(:)
    integer  :: i, j, ex, cr, nrep, stride, id
    integer  :: parmsetid, umbrid, numbs, nrep2
    real(wp) :: prob, local_disp, diff, s, s2


    nrep = remd%nreplicas(dimno)
    if (allocated(disp) .and. size(disp) < nrep) then
      deallocate(disp)
    end if

    if (.not. allocated(disp)) then
      allocate(disp(nrep))
    end if
    disp(1:nrep) = 0.0_wp

    ! prep stride
    !
    stride = 1
    if (dimno < remd%dimension) then
      stride = stride * product(remd%nreplicas(dimno+1:remd%dimension))
    end if

    ! only the lowest line will be used for multidimensional case
    !
    do i = 1, nrep - 1
      id = 1 + stride * (i - 1)
      cr = remd%num_criteria (id,dimno,2-mod(i,2))
      ex = remd%num_exchanges(id,dimno,2-mod(i,2))
      prob = real(ex,wp) / real(cr,wp) - remd%autoadj(dimno)%tgt_exc_prob
      if (abs(prob) >= remd%autoadj(dimno)%mgn_exc_prob) then
        prob = sign(1.0_wp,prob) * (abs(prob) - remd%autoadj(dimno)%mgn_exc_prob)
        !! change parameter
        disp(i) = floor(100.0_wp * prob) * remd%autoadj(dimno)%param_grid
        if (abs(disp(i)) > remd%autoadj(dimno)%max_param_shift) then
          disp(i) = sign(1.0_wp,disp(i)) * remd%autoadj(dimno)%max_param_shift
        end if
      end if
    end do

    if (remd%types(dimno) == RemdRestraint) then
      parmsetid = remd%repid2parmsetid(my_country_no+1)
      umbrid    = remd%iparameters(dimno,remd%parmidsets(parmsetid,dimno))
      numbs     = remd%umbrid2numfuncs(umbrid)
      if (remd%autoadj(dimno)%fix_terminal == RemdAutoFixBottom) then
        ! change from the bottom
        do j = 1, numbs
          do i = 1, nrep - 1
            diff = remd%rest_reference(i+1,j) - remd%rest_reference(i,j)
            s = sign(1.0_wp,diff)
            s2 = remd%rest_reference(i+1,j) + s * disp(i) - &
                 remd%rest_reference(i,j) - &
                 s * remd%autoadj(dimno)%param_grid
            s2 = sign(1.0_wp,s2)
            if (s /= s2) then
              ! have to prevent reverting
              disp(i) = (remd%rest_reference(i,j) + &
                         s * remd%autoadj(dimno)%param_grid) - &
                        remd%rest_reference(i+1,j)
            end if
            remd%rest_reference(i+1:nrep,j) = &
                remd%rest_reference(i+1:nrep,j) + s * disp(i)
          end do
        end do
      else if (remd%autoadj(dimno)%fix_terminal == RemdAutoFixTop) then
        ! change from the top
        do j = 1, numbs
          do i = nrep - 1, 1, -1
            diff = remd%rest_reference(i,j) - remd%rest_reference(i+1,j)
            s = sign(1.0_wp,diff)
            s2 = remd%rest_reference(i,j) + s * disp(i) - &
                 remd%rest_reference(i+1,j) - &
                 s * remd%autoadj(dimno)%param_grid
            s2 = sign(1.0_wp,s2)
            if (s /= s2) then
              ! have to prevent reverting
              disp(i) = (remd%rest_reference(i+1,j) - &
                         s * remd%autoadj(dimno)%param_grid) - &
                        remd%rest_reference(i,j)
            end if
            remd%rest_reference(1:i,j) = &
                remd%rest_reference(1:i,j) + s * disp(i)
          end do
        end do
      end if
    else
      ! update "dparameters" according to "disp"
      if (remd%autoadj(dimno)%fix_terminal == RemdAutoFixBottom) then
        ! change from the bottom
        do i = 1, nrep - 1
          diff = remd%dparameters(dimno,i+1) - remd%dparameters(dimno,i)
          s = sign(1.0_wp,diff)
          s2 = remd%dparameters(dimno,i+1) + s * disp(i) - &
               remd%dparameters(dimno,i) - &
               s * remd%autoadj(dimno)%param_grid
          s2 = sign(1.0_wp,s2)
          if (s /= s2) then
            ! have to prevent reverting
            disp(i) = (remd%dparameters(dimno,i) + &
                       s * remd%autoadj(dimno)%param_grid) - &
                      remd%dparameters(dimno,i+1)
          end if
          remd%dparameters(dimno,i+1:nrep) = &
              remd%dparameters(dimno,i+1:nrep) + s * disp(i)
        end do
      else if (remd%autoadj(dimno)%fix_terminal == RemdAutoFixTop) then
        ! change from the top
        do i = nrep - 1, 1, -1
          diff = remd%dparameters(dimno,i) - remd%dparameters(dimno,i+1)
          s = sign(1.0_wp,diff)
          s2 = remd%dparameters(dimno,i) + s * disp(i) - &
               remd%dparameters(dimno,i+1) - &
               s * remd%autoadj(dimno)%param_grid
          s2 = sign(1.0_wp,s2)
          if (s /= s2) then
            ! have to prevent reverting
            disp(i) = (remd%dparameters(dimno,i+1) - &
                       s * remd%autoadj(dimno)%param_grid) - &
                      remd%dparameters(dimno,i)
          end if
          remd%dparameters(dimno,1:i) = &
              remd%dparameters(dimno,1:i) + s * disp(i)
        end do
      end if
    end if

    return

  end subroutine remd_autoadj_displacement

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    temperature_remd
  !> @brief        temperature REMD
  !! @authors      TM
  !! @param[in]    dimno    : dimension number
  !! @param[in]    exptrn   : exchange pattern
  !! @param[inout] domain   : domain information
  !! @param[inout] dynvars  : dynamics variables information
  !! @param[inout] remd     : REMD information
  !! @note         Y.Sugita & Y.Okamoto, Chem.Phys.Lett., 314, 141−151 (1999)
  !!               K.Hukushima & K.Nemoto, J.Phys.Soc.Jpn., 65, 1604-1608 (1996)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074001 (2010)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074003 (2010)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine temperature_remd(dimno, exptrn, domain, dynvars, remd)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_domain),   target, intent(inout) :: domain
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_remd),     target, intent(inout) :: remd

    ! local variables
    integer                   :: i, j, jx
    integer                   :: repid_i, repid_j, neighid
    integer                   :: parmsetid, parmsetid_i, parmsetid_j
    real(wp)                  :: temp_i, temp_j, beta_i, beta_j, factor
    real(wp)                  :: energy_i, energy_j, delta
    real(wp)                  :: before_gather, eneval
    logical                   :: do_trial, do_exchange
    integer,          pointer :: total_nreplicas, iseed
    integer,          pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    real(wp),         pointer :: parameters(:,:)
    real(wp),         pointer :: after_gather(:), random_table(:,:)


    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%dparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    iseed           => remd%iseed
    after_gather    => remd%after_gather
    random_table    => remd%random_table


    ! allgather potential energy
    !
    eneval = dynvars%energy%total
    call mpi_bcast(eneval, 1, mpi_wp_real, 0, mpi_comm_country, ierror)

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%potential_energy(i) = after_gather(i)
    end do

    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if


      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        energy_i = remd%potential_energy(repid_i)
        energy_j = remd%potential_energy(repid_j)
        temp_i   = parameters(dimno,parmidsets(parmsetid_i,dimno))
        temp_j   = parameters(dimno,parmidsets(parmsetid_j,dimno))
        beta_i   = 1.0_wp / (KBOLTZ * temp_i)
        beta_j   = 1.0_wp / (KBOLTZ * temp_j)
        delta    = (beta_j - beta_i)*(energy_i - energy_j)

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if


        if (do_exchange) then

          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1

          ! scale velocities and themostat & barostat momenta
          !
          if (i == my_country_no + 1) then

            factor = sqrt(temp_j/temp_i)
            do j = 1, domain%num_cell_local
              do jx = 1, domain%num_atom(j)
                domain%velocity(1:3,jx,j) = domain%velocity(1:3,jx,j) * factor
              end do
            end do

            factor = temp_j/temp_i
            dynvars%thermostat_momentum                 &
              = dynvars%thermostat_momentum * factor
            dynvars%barostat_momentum(1:3)              &
              = dynvars%barostat_momentum(1:3) * factor

          end if

          repid2parmsetid(i) = parmsetid_j
        end if
 
      end if

    end do

    return

  end subroutine temperature_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pressure_remd
  !> @brief        pressure REMD
  !! @authors      TM
  !! @param[in]    dimno    : dimension number
  !! @param[in]    exptrn   : exchange pattern
  !! @param[in]    ensemble : ensemble information
  !! @param[in]    boundary : boundary information
  !! @param[inout] remd     : REMD information
  !! @note         T.Okabe et al., Chem.Phys.Lett., 335, 435−439 (2001)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074003 (2010)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pressure_remd(dimno, exptrn, ensemble, boundary, remd)

    ! formal arguments
    integer,                 intent(in)    :: dimno
    integer,                 intent(in)    :: exptrn
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),        intent(in)    :: boundary
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k
    integer                  :: repid_i, repid_j, neighid
    integer                  :: parmsetid, parmsetid_i, parmsetid_j
    real(wp)                 :: beta, delta, temperature
    real(wp)                 :: volume_i, volume_j, press_i, press_j
    real(wp)                 :: before_gather
    logical                  :: do_trial, do_exchange
    integer,         pointer :: total_nreplicas
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,         pointer :: parmidsets(:,:)
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    real(wp),        pointer :: parameters(:,:)
    real(wp),        pointer :: after_gather(:), random_table(:,:)


    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%dparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    after_gather    => remd%after_gather
    random_table    => remd%random_table


    ! allgather volume
    !
    before_gather =  boundary%box_size_x &
                   * boundary%box_size_y &
                   * boundary%box_size_z 
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%volume(i) = after_gather(i)
    end do

    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if


      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        call get_replica_temperature(parmsetid_i, remd, ensemble, temperature)

        volume_i = remd%volume(repid_i)
        volume_j = remd%volume(repid_j)
        beta     = 1.0_wp / (KBOLTZ * temperature)
        press_i  = parameters(dimno,parmidsets(parmsetid_i,dimno))*ATMOS_P
        press_j  = parameters(dimno,parmidsets(parmsetid_j,dimno))*ATMOS_P
        delta    = beta*(press_i - press_j)*(volume_j - volume_i)

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (do_exchange) then
          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1

          repid2parmsetid(i) = parmsetid_j
        end if
 
      end if

    end do

    return

  end subroutine pressure_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    surface_tension_remd
  !> @brief        surface-tension REMD
  !! @authors      TM
  !! @param[in]    dimno    : dimension number
  !! @param[in]    exptrn   : exchange pattern
  !! @param[in]    ensemble : ensemble information
  !! @param[in]    boundary : boundary information
  !! @param[inout] remd     : REMD information
  !! @note         T.Mori et al., J.Chem.Theory.Comput., 9, 5629−5640 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine surface_tension_remd(dimno, exptrn, ensemble, boundary, remd)

    ! formal arguments
    integer,                 intent(in)    :: dimno
    integer,                 intent(in)    :: exptrn
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),        intent(in)    :: boundary
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k
    integer                  :: repid_i, repid_j, neighid
    integer                  :: parmsetid, parmsetid_i, parmsetid_j
    real(wp)                 :: beta, factor, delta, temperature
    real(wp)                 :: area_i, area_j, gamma_i, gamma_j
    real(wp)                 :: before_gather
    logical                  :: do_trial, do_exchange
    integer,         pointer :: total_nreplicas
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,         pointer :: parmidsets(:,:)
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    real(wp),        pointer :: parameters(:,:)
    real(wp),        pointer :: after_gather(:), random_table(:,:)


    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%dparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    after_gather    => remd%after_gather
    random_table    => remd%random_table


    ! allgather area
    !
    before_gather =  boundary%box_size_x &
                   * boundary%box_size_y 
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%area(i) = after_gather(i)
    end do

    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if


      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        call get_replica_temperature(parmsetid_i, remd, ensemble, temperature)

        factor  = ATMOS_P*100.0_wp/1.01325_wp
        area_i  = remd%area(repid_i)
        area_j  = remd%area(repid_j)
        gamma_i = parameters(dimno,parmidsets(parmsetid_i,dimno))*factor
        gamma_j = parameters(dimno,parmidsets(parmsetid_j,dimno))*factor
        beta    = 1.0_wp / (KBOLTZ * temperature)
        delta   = -beta*(gamma_i - gamma_j)*(area_j - area_i)

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (do_exchange) then
          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1

          repid2parmsetid(i) = parmsetid_j
        end if
 
      end if

    end do

    return

  end subroutine surface_tension_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    restraint_remd
  !> @brief        Replica exchange umbrella sampling (REUS)
  !! @authors      TM
  !! @param[in]    dimno    : dimension number
  !! @param[in]    exptrn   : exchange pattern
  !! @param[inout] domain   : domain information
  !! @param[inout] boundary : boundary information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamics variables information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] remd     : REMD information
  !! @note         Y.Sugita et al., J.Chem.Phys., 113, 6042-6051 (2000)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine restraint_remd(dimno, exptrn, domain, boundary, ensemble, dynvars,&
                            enefunc, remd)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_domain),   target, intent(inout) :: domain
    type(s_boundary),         intent(inout) :: boundary
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_remd),     target, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k, ix
    integer                  :: natoms
    integer                  :: repid_i, repid_j
    integer                  :: parmsetid, parmsetid_i, parmsetid_j
    integer                  :: neighid
    integer                  :: funcid, umbr_0, umbr_1
    real(wp)                 :: cv
    real(wp)                 :: energy_i_0, energy_j_0
    real(wp)                 :: energy_i_1, energy_j_1
    real(dp)                 :: eneval
    real(wp)                 :: beta, delta, temperature
    real(wp)                 :: before_gather
    logical                  :: do_trial, do_exchange, umbr_exist
    integer,         pointer :: total_nreplicas
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,         pointer :: parmidsets(:,:)
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,         pointer :: parameters(:,:)
    real(wp),        pointer :: after_gather(:), random_table(:,:)
    real(wp),        pointer :: refcoord0(:,:), refcoord1(:,:)


    after_gather    => remd%after_gather
    random_table    => remd%random_table
    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%iparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    refcoord0       => remd%restraint_refcoord0
    refcoord1       => remd%restraint_refcoord1

    ! set target restraint function
    !
    repid_i     = my_country_no + 1
    parmsetid_i = repid2parmsetid(repid_i)
    umbr_0      = parameters(dimno,parmidsets(parmsetid_i,dimno))


    ! calculate restraint potential energy at the current condition
    !
    do k = 1, remd%umbrid2numfuncs(umbr_0)
      funcid = remd%umbrid2funclist(umbr_0,k)
      enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbr_0,k)
      enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbr_0,k)

      if (enefunc%restraint_kind(funcid) == RestraintsFuncRMSD) then
        enefunc%rmsd_force  = enefunc%restraint_const(1,funcid)
        enefunc%rmsd_target = enefunc%restraint_ref  (1,funcid)
      end if
    end do

    call compute_energy_restraints_remd(.true., enefunc, domain, boundary, eneval)
    before_gather = eneval

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%restraint_energy_0(i) = after_gather(i)
    end do


    ! calculate restraint potential energy using the neighboring condition
    !
    call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

    if (neighid /= 0) then
      parmsetid_j = neighid
      repid_j     = parmsetid2repid(parmsetid_j)
      umbr_exist  = .true.
    else
      umbr_exist  = .false.
    end if

    if (umbr_exist) then
      umbr_1 = parameters(dimno,parmidsets(parmsetid_j,dimno))

      do k = 1, remd%umbrid2numfuncs(umbr_1)
        funcid = remd%umbrid2funclist(umbr_1,k)
        enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbr_1,k)
        enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbr_1,k)

        if (enefunc%restraint_kind(funcid) == RestraintsFuncRMSD) then
          enefunc%rmsd_force  = enefunc%restraint_const(1,funcid)
          enefunc%rmsd_target = enefunc%restraint_ref  (1,funcid)
        end if
      end do

      ! compute restraint energy
      !
      call compute_energy_restraints_remd(.false., enefunc, domain, boundary, eneval)

    else
      eneval = 0.0_wp
    end if

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%restraint_energy_1(i) = after_gather(i)
    end do


    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if


      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        call get_replica_temperature(parmsetid_i, remd, ensemble, temperature)

        energy_i_0 = remd%restraint_energy_0(repid_i)
        energy_i_1 = remd%restraint_energy_1(repid_i)
        energy_j_0 = remd%restraint_energy_0(repid_j)
        energy_j_1 = remd%restraint_energy_1(repid_j)
        beta = 1.0_wp / (KBOLTZ * temperature)
        delta = beta*(energy_j_1 - energy_i_0 - energy_j_0 + energy_i_1)

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (do_exchange) then
          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1
          repid2parmsetid(i) = parmsetid_j
        end if

      end if
    end do

    return

  end subroutine restraint_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    solute_tempering_remd
  !> @brief        REST
  !! @authors      MK
  !! @param[in]    dimno      : dimension number
  !! @param[in]    exptrn     : exchange pattern
  !! @param[inout] domain     : domain information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] dynvars    : dynamics variables information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !! @param[inout] pairlist   : pairlist information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solute_tempering_remd(dimno, exptrn, domain, ensemble,  &
                                   dynvars, enefunc, remd, pairlist, &
                                   boundary, output)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_domain),   target, intent(inout) :: domain
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output

    ! local variables
    integer                   :: i, j, jx
    integer                   :: repid_i, repid_j, neighid
    integer                   :: parmsetid, parmsetid_i, parmsetid_j
    real(wp)                  :: beta, factor
    real(wp)                  :: energy_i, energy_j, delta
    real(wp)                  :: before_gather, eneval
    logical                   :: does_exist, do_trial, do_exchange
    integer,          pointer :: total_nreplicas, iseed
    integer,          pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    real(wp),         pointer :: parameters(:,:)
    real(wp),         pointer :: after_gather(:), random_table(:,:)

    ! REST
    integer                   :: nc, nv
    type(s_energy),      save :: energy
    real(wip),    allocatable :: force(:,:,:)
    real(wip),    allocatable :: force_long(:,:,:)
    real(wp),     allocatable :: force_omp(:,:,:,:)
    real(wp),     allocatable :: force_pbc(:,:,:,:)
    real(dp),     allocatable :: virial_cell(:,:)
    real(dp)                  :: virial(3,3)
    real(dp)                  :: virial_long(3,3)
    real(dp)                  :: virial_extern(3,3)

    ! REST allocation
    nc = size(domain%force(3,MaxAtom,:))
    nv = size(domain%virial_cellpair(3,:))
    if (.not. allocated(force)) then
      allocate(force(3,MaxAtom,nc), &
               force_long(3,MaxAtom,nc), &
               force_omp(3,MaxAtom,nc,nthread), &
               virial_cell(3,nv))
    end if

    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      allocate(force_pbc(MaxAtom,3, nc, nthread))

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      allocate(force_pbc(3, MaxAtom*nc,1,nthread))

    else if (domain%nonbond_kernel == NBK_Intel) then

      allocate(force_pbc(MaxAtom*nc,3,1,nthread))

    else if (domain%nonbond_kernel == NBK_GPU) then

      allocate(force_pbc(MaxAtom*nc*4,1,1,nthread))

    end if

    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%dparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    iseed           => remd%iseed
    after_gather    => remd%after_gather
    random_table    => remd%random_table

    ! calculate energy at current coordinate
    call compute_energy(domain, enefunc, pairlist, boundary, domain%coord, &
                        ensemble%use_barostat, .true., .true., .true.,     &
                        enefunc%nonb_limiter,                              &
                        energy, domain%atmcls_pbc,                         &
                        domain%translated, force, force_long,              &
                        force_omp, force_pbc, virial_cell, virial,         &
                        virial_long, virial_extern)

    ! allgather potential energy
    !
    eneval = energy%total
    call mpi_bcast(eneval, 1, mpi_wp_real, 0, mpi_comm_country, ierror)

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%potential_energy(i) = after_gather(i)
    end do

    ! calculate potential energy using the neighbouring condition
    !
    repid_i     = my_country_no + 1
    parmsetid_i = repid2parmsetid(repid_i)

    ! get neighboring parameter set index
    !
    call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

    if (neighid /= 0) then
      parmsetid_j = neighid
      repid_j     = parmsetid2repid(parmsetid_j)
      does_exist  = .true.
    else
      does_exist  = .false.
    end if

    if (does_exist) then
      ! replace enefunc
      !
      call assign_condition_rest(parmsetid_j, remd, domain, ensemble, &
                                 enefunc)

      ! calculate energy
      !
      call compute_energy(domain, enefunc, pairlist, boundary, domain%coord, &
                          ensemble%use_barostat, .true., .true., .true., &
                          .false.,                                       &
                          energy, domain%atmcls_pbc,                     &
                          domain%translated, force, force_long,          &
                          force_omp, force_pbc, virial_cell, virial,     &
                          virial_long, virial_extern)
      eneval = energy%total
      call mpi_bcast(eneval, 1, mpi_wp_real, 0, mpi_comm_country, ierror)
    else
      eneval = 0.0_wp
    end if

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%potential_energy(i) = remd%potential_energy(i) - after_gather(i)
    end do

    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if

      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        beta     = 1.0_wp / (KBOLTZ * ensemble%temperature)
        energy_i = remd%potential_energy(repid_i)
        energy_j = remd%potential_energy(repid_j)
        delta    = - beta * (energy_i + energy_j)

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (replica_main_rank .and. does_exist &
            .and. repid_i == my_country_no + 1) then
          write(output%logunit, &
                '("  REST Exchange> Dim myID tgtID DE(my E) result:")')
          write(output%logunit,'(3(I3,1x),F10.4," (",F10.4,")",L2)') &
            dimno, repid2parmsetid(repid_i), neighid, &
            energy_i + energy_j, energy_i, do_exchange
          write(output%logunit,*)
        end if

        if (do_exchange) then

          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1

          repid2parmsetid(i) = parmsetid_j
        else
          ! have to put back enefunc to the original state
          !
          call assign_condition_rest(parmsetid_i, remd, domain, ensemble, &
                                     enefunc)
        end if

      end if

    end do


    return

  end subroutine solute_tempering_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition
  !> @brief        control replica exchange
  !! @authors      TM
  !! @param[in]    parmsetid   : parameter set id
  !! @param[inout] remd        : REMD information
  !! @param[inout] domain      : domain information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition(parmsetid, remd, domain, ensemble, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_remd),            intent(inout) :: remd
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: umbrid, funcid, atmid


    do i = 1, remd%dimension
      if      (remd%types(i) == RemdTemperature) then
        ensemble%temperature = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdPressure   ) then
        ensemble%pressure    = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdGamma      ) then
        ensemble%gamma       = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdRestraint  ) then
        umbrid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
        do k = 1, remd%umbrid2numfuncs(umbrid)
          funcid = remd%umbrid2funclist(umbrid,k)
          enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbrid,k)
          enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbrid,k)

          if (enefunc%restraint_kind(funcid) == RestraintsFuncRMSD) then
            enefunc%rmsd_force  = enefunc%restraint_const(1,funcid)
            enefunc%rmsd_target = enefunc%restraint_ref  (1,funcid)
          end if

        end do
      end if
    end do

    call assign_condition_rest(parmsetid, remd, domain, ensemble, enefunc)

    return

  end subroutine assign_condition

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_rest
  !> @brief        reassign force constant involving REST
  !! @authors      MK
  !! @param[in]    parmsetid   : parameter set id
  !! @param[inout] remd        : remd information
  !! @param[inout] domain      : domain information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_rest(parmsetid, remd, domain, ensemble, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_remd),            intent(inout) :: remd
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i
    real(wp)                 :: rest_param

    do i = 1, remd%dimension

      if (remd%types(i) == RemdSoluteTempering) then

        if (ensemble%ensemble == EnsembleNVE) then
          call error_msg('Solute_Tempering> REST is not available for NVE!')
        end if

        rest_param = remd%dparameters(i,remd%parmidsets(parmsetid,i))
        call assign_condition_solute_tempering(remd%solute_tempering(i), &
                                               rest_param, domain, &
                                               ensemble, enefunc)
      end if

    end do

    return

  end subroutine assign_condition_rest

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solte_tempering
  !> @brief        reassign force constant involving REST
  !! @authors      MK
  !! @param[inout] soltemp     : REST information
  !! @param[in]    rest_param  : target temperature
  !! @param[inout] domain      : domain information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering(soltemp, rest_param, &
                                               domain, ensemble, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    real(wp),                intent(in)    :: rest_param
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc

    integer  :: alist(8), i, j, ix, num, tgt, n, org, ncell, counter
    real(wp) :: coeff_full, coeff_half, wgt, tmp1, tmp2
    real(wp) :: el_fact, alpha
    real(dp) :: u_self_org


    coeff_full = ensemble%temperature / rest_param
    coeff_half = sqrt(coeff_full)

    if (soltemp%done_setup) then
      tmp1 = coeff_full
      tmp2 = coeff_half

      coeff_full = coeff_full / soltemp%rest_param_full
      coeff_half = coeff_half / soltemp%rest_param_half

      ! remember coeffs for next exchange
      !
      soltemp%rest_param_full = tmp1
      soltemp%rest_param_half = tmp2
    else
      soltemp%rest_param_full = coeff_full
      soltemp%rest_param_half = coeff_half
      soltemp%done_setup = .true.
    end if

    ncell = domain%num_cell_local + domain%num_cell_boundary

    if (soltemp%sw_charge) then
      do i = 1, ncell
        ! charge
        do ix = 1, domain%num_solute(i)
          if (soltemp%is_solute(domain%id_l2g_solute(ix,i)) > 0) then
            domain%charge(ix,i) = coeff_half * domain%charge(ix,i)
          end if
        end do
      end do

      ! need to update PME self energy (not necessary for NPT?)
      u_self = 0.0_dp
      el_fact = ELECOEF / enefunc%dielec_const
      alpha = enefunc%pme_alpha
      do i = 1, domain%num_cell_local
        do ix = 1, domain%num_atom(i)
          u_self = u_self + domain%charge(ix,i)**2
        end do
      end do
      u_self = - u_self * el_fact * alpha / sqrt(PI)
    end if

    do i = 1, domain%num_cell_local

      ! bond
      if (soltemp%sw_bonds) then
        call assign_condition_solute_tempering_internal(soltemp, &
                  2, coeff_full, &
                  enefunc%num_bond(i), enefunc%bond_list(:,:,i), &
                  enefunc%bond_force_const(:,i))
      end if

      ! angle and urey
      do ix = 1, enefunc%num_angle(i)
        alist(1:3) = enefunc%angle_list(1:3,ix,i)
        num = 0
        do j = 1, 3
          if (soltemp%is_solute(alist(j)) > 0) num = num + 1
        end do
        if (num > 0 .and. soltemp%sw_angles) then
          wgt = real(num,wp) / 3.0_wp
          wgt = coeff_full ** wgt
          enefunc%angle_force_const(ix,i) = &
                     wgt * enefunc%angle_force_const(ix,i)
        end if
        if (enefunc%urey_force_const(ix,i) > EPS .and. &
             soltemp%sw_ureys) then
          num = 0
          if (soltemp%is_solute(alist(1)) > 0) num = num + 1
          if (soltemp%is_solute(alist(3)) > 0) num = num + 1
          if (num > 0) then
            wgt = real(num,wp) / 2.0_wp
            wgt = coeff_full ** wgt
            enefunc%urey_force_const(ix,i) = &
                     wgt * enefunc%urey_force_const(ix,i)
          end if
        end if
      end do

      ! dihedral
      if (soltemp%sw_dihedrals) then
        call assign_condition_solute_tempering_internal(soltemp, &
                  4, coeff_full, &
                  enefunc%num_dihedral(i), enefunc%dihe_list(:,:,i), &
                  enefunc%dihe_force_const(:,i))
      end if

      ! rb_dihedral
      if (soltemp%sw_rb_dihedrals) then
        do ix = 1, enefunc%num_rb_dihedral(i)
          alist(1:4) = enefunc%rb_dihe_list(1:4,ix,i)
          num = 0
          do j = 1, 4
            if (soltemp%is_solute(alist(j)) > 0) num = num + 1
          end do
          if (num > 0) then
            wgt = real(num,wp) / real(n,wp)
            wgt = coeff_full ** wgt
            enefunc%rb_dihe_c(1:6,ix,i) = &
                   wgt * enefunc%rb_dihe_c(1:6,ix,i)
          end if
        end do
      end if

      ! impropers
      if (soltemp%sw_impropers) then
        call assign_condition_solute_tempering_internal(soltemp, &
                  4, coeff_full, &
                  enefunc%num_improper(i), enefunc%impr_list(:,:,i), &
                  enefunc%impr_force_const(:,i))
      end if

      ! contact
      if (soltemp%sw_contacts) then
        do ix = 1, enefunc%num_contact(i)
          alist(1:2) = enefunc%contact_list(1:2,ix,i)
          num = 0
          do j = 1, 2
            if (soltemp%is_solute(alist(j)) > 0) num = num + 1
          end do
          if (num > 0) then
            wgt = real(num,wp) / 2.0_wp
            wgt = coeff_full ** wgt
            enefunc%contact_lj6(ix,i)  = wgt * enefunc%contact_lj6(ix,i)
            enefunc%contact_lj12(ix,i) = wgt * enefunc%contact_lj12(ix,i)
          end if
        end do
      end if

    end do

    ! lj
    if (soltemp%sw_lj) then
      n = enefunc%num_atom_cls
      if (allocated(enefunc%nb14_lj6)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj6)
      end if
      if (allocated(enefunc%nb14_lj12)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj12)
      end if
      if (allocated(enefunc%nonb_lj6)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj6)
      end if
      if (allocated(enefunc%nonb_lj12)) then
        call assign_condition_solute_tempering_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj12)
      end if
#ifdef USE_GPU
      if (allocated(enefunc%nonb_lj12) .or. allocated(enefunc%nonb_lj6)) then
        call gpu_upload_lj_coeffs(n, enefunc%nonb_lj12, enefunc%nonb_lj6);
      end if
#endif /* USE_GPU */
    end if

    ! cmap
    if (soltemp%num_cmap_type > 0 .and. soltemp%sw_cmaps) then
      do i = 1, soltemp%num_cmap_type
        tgt = i + soltemp%istart_cmap_type - 1
        org = soltemp%cmap_type_org(i)
        wgt = soltemp%rest_param_full ** soltemp%cmap_weight(i)
        enefunc%cmap_coef(1:4,1:4,1:24,1:24,tgt) = &
          wgt * enefunc%cmap_coef(1:4,1:4,1:24,1:24,org)
      end do
    end if

    return

  end subroutine assign_condition_solute_tempering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solte_tempering_lj
  !> @brief        reassign force constant involving REST for LJ
  !! @authors      MK
  !! @param[in]    n           : array size of interaction matrix
  !! @param[in]    coeff_half  : coefficient for solute-solvent
  !! @param[in]    coeff_full  : coefficient for solute-solute
  !! @param[in]    soltemp     : REST information
  !! @param[inout] nbcoeff     : LJ interaction matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering_lj(n, coeff_half, coeff_full, &
                                                  soltemp, nbcoeff)

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_half
    real(wp),                intent(in)    :: coeff_full
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(inout) :: nbcoeff(n,n)

    ! local variables
    integer :: i, j, oldcount, newcount, org, orgj


    oldcount = soltemp%istart_atom_cls - 1
    newcount = oldcount + soltemp%num_atom_cls

    do i = 1, soltemp%num_atom_cls
      org = soltemp%atom_cls_no_org(i)
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = coeff_half * nbcoeff(j,org)
        nbcoeff(i+oldcount,j) = coeff_half * nbcoeff(j,org)
      end do
      do j = oldcount + 1, newcount
        orgj = soltemp%atom_cls_no_org(j-oldcount)
        nbcoeff(j,i+oldcount) = coeff_full * nbcoeff(orgj,org)
        nbcoeff(i+oldcount,j) = coeff_full * nbcoeff(orgj,org)
      end do
    end do

    return

  end subroutine assign_condition_solute_tempering_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_solte_tempering_internal
  !> @brief        reassign force constant involving REST for internal
  !! @authors      MK
  !! @param[in]    soltemp     : REST information
  !! @param[in]    n           : number of indexes involved
  !! @param[in]    coeff_full  : coefficient
  !! @param[in]    n_internal  : number of terms
  !! @param[in]    aindex      : atom list
  !! @param[inout] fc          : force constant
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering_internal(soltemp, n, &
                                  coeff_full, n_internal, aindex, fc)

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: n_internal
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)

    ! local variables
    integer  :: alist(1:n), ix, j, num
    real(wp) :: wgt


    do ix = 1, n_internal
      alist(1:n) = aindex(1:n,ix)
      num = 0
      do j = 1, n
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if (num > 0) then
        wgt = real(num,wp) / real(n,wp)
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
    end do

    return

  end subroutine assign_condition_solute_tempering_internal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_neighbouring_parmsetid
  !> @brief        get neighbouring parameter set id
  !! @authors      TM
  !! @param[in]    dimno        : dimension number
  !! @param[in]    exptrn       : exchange pattern
  !! @param[in]    remd         : REMD information
  !! @param[in]    parmsetid_i  : parameter set id
  !! @param[out]   neighbour_id : parameter set id neighbouring parmsetid_i
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, &
                                        neighbour_id)

    ! formal arguments
    integer,                 intent(in)    :: dimno
    integer,                 intent(in)    :: exptrn
    type(s_remd),    target, intent(in)    :: remd
    integer,                 intent(in)    :: parmsetid_i
    integer,                 intent(out)   :: neighbour_id

    ! local variables
    integer                  :: j, k
    integer                  :: next, jump, parmid_i
    logical                  :: neighbour1, neighbour2
    logical                  :: cyclic
    integer,         pointer :: dimension, num_parmidsets
    integer,         pointer :: parmidsets(:,:)
    

    cyclic = remd%cyclic_params(dimno)

    dimension      => remd%dimension
    num_parmidsets => remd%total_nreplicas
    parmidsets     => remd%parmidsets

    parmid_i = parmidsets(parmsetid_i,dimno)
    if (exptrn == 1) then
      if (mod(parmid_i,2) == 0) then
        ! 1 <2  3 <4  5 <6  7 <8  ...
        next = -1
        if (cyclic) then
          jump = remd%nreplicas(dimno) - 1
        else
          jump = next
        end if
      else
        ! 1> 2  3> 4  5> 6  7> 8  ...
        next = +1
        if (cyclic) then
          jump = -remd%nreplicas(dimno) + 1
        else
          jump = next
        end if
      end if
    else if (exptrn == 2) then
      if (mod(parmid_i,2) == 0) then
        ! 1  2> 3  4> 5  6> 7  8  ...
        next = +1
        if (cyclic) then
          jump = -remd%nreplicas(dimno) + 1
        else
          jump = next
        end if
      else
        ! 1  2 <3  4 <5  6 <7  8  ...
        next = -1
        if (cyclic) then
          jump = remd%nreplicas(dimno) - 1
        else
          jump = next
        end if
      end if
    end if

    do j = 1, num_parmidsets
      neighbour1 = .true.
      neighbour2 = .false.
      do k = 1, dimension
        if (k /= dimno) then
          if (parmidsets(j,k) /= parmidsets(parmsetid_i,k)) &
            neighbour1 = .false.
        else
          if (parmidsets(j,k) == parmidsets(parmsetid_i,k) + next) then
            neighbour2 = .true.
          else if (parmidsets(j,k) == parmidsets(parmsetid_i,k) + jump) then
            neighbour2 = .true.
          end if
        end if
      end do
      if (neighbour1 .and. neighbour2) then
        neighbour_id = j
        exit
      else
        neighbour_id = 0
      end if
    end do

    return

  end subroutine get_neighbouring_parmsetid

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_replica_temperature
  !> @brief        get temperature of each replica
  !! @authors      TM
  !! @param[in]    parmsetid   : parameter set id
  !! @param[in]    remd        : REMD information
  !! @param[in]    ensemble    : ensemble information
  !! @param[out]   temperature : temperature of each replica
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_replica_temperature(parmsetid, remd, ensemble, temperature)

    ! formal arguments
    integer,          intent(in)  :: parmsetid
    type(s_remd),     intent(in)  :: remd
    type(s_ensemble), intent(in)  :: ensemble
    real(wp),         intent(out) :: temperature

    ! local variables
    integer                       :: i


    temperature = ensemble%temperature

    do i = 1, remd%dimension
      if (remd%types(i) == RemdTemperature) then
        temperature = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      end if
    end do

    return

  end subroutine get_replica_temperature

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_remd
  !> @brief        compute restraint potential energy in REUS routines
  !! @authors      TM
  !! @param[in]    get_coord : flag for whether to get coordinates
  !! @param[inout] enefunc   : enefunc information
  !! @param[inout] domain    : domain information
  !! @param[inout] boundary  : boundary information
  !! @param[out]   eneval    : energy value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_remd(get_coord, enefunc, domain, &
                                            boundary, eneval)

    ! formal arguments
    logical,                  intent(in)    :: get_coord
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_domain),   target, intent(inout) :: domain
    type(s_boundary), target, intent(inout) :: boundary
    real(dp),                 intent(out)   :: eneval

    ! local variables
    integer                  :: i, id
    logical                  :: do_posres
    real(dp)                 :: ebonds, eposi, eposi_mpi, ermsd, eemfit
    real(dp)                 :: rmsd, emcorr
    real(dp)                 :: eposi_omp(nthread)
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(wp),        pointer :: force_omp(:,:,:,:)


    force_omp  => domain%force_omp

    eneval     = 0.0_dp
    ebonds     = 0.0_dp
    ermsd      = 0.0_dp
    eemfit     = 0.0_dp
    rmsd       = 0.0_dp
    emcorr     = 0.0_dp

    eposi_omp(1:nthread) = 0.0_dp

    call compute_energy_restraints(get_coord, .false., domain, boundary,  &
                                   enefunc, domain%coord, force_omp,      &
                                   virial_omp, virial_ext_omp, eposi_omp, &
                                   ermsd, rmsd, ebonds, eemfit, emcorr)

    eneval = ebonds + ermsd + eemfit

    return

  end subroutine compute_energy_restraints_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_remd_fep
  !> @brief        control replica exchange for FEP
  !! @authors      HO
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
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_remd_fep(output, domain, enefunc, dynvars, dynamics, &
                      pairlist, boundary, constraints, ensemble,  &
                      comm, remd, alchemy)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd
    type(s_alchemy),optional,intent(inout) :: alchemy

    ! local variables
    integer                   :: i, replicaid, parmsetid

    ! FEP
    logical                   :: is_feprest

    is_feprest = .false.
    do i = 1, remd%dimension
      if (remd%types(i) == RemdAlchemyRest) then
        is_feprest = .true.
      end if
    end do

    ! Open output files
    !
    call open_output(output)

    ! output restart data
    !
    dynvars%step = 0
    call output_remd(0, output, enefunc, domain, dynamics, dynvars, boundary, &
                     constraints, remd)

    do i = 1, remd%ncycles

      dynamics%istart_step  = (i-1)*remd%exchange_period + 1
      dynamics%iend_step    =  i   *remd%exchange_period
      dynamics%initial_time = dynvars%time

      ! set conditions
      !
      replicaid = my_country_no + 1
      parmsetid = remd%repid2parmsetid(replicaid)
      ! FEP: assign condition
      call assign_condition_fep(parmsetid, remd, domain, ensemble, enefunc, alchemy)

      ! MD main loop
      !
      if (dynamics%integrator == IntegratorLEAP) then
        call error_msg('LEAP integrator is not available in FEP')
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd, alchemy)
      else if (dynamics%integrator == IntegratorVRES) then
        call vverlet_respa_dynamics(output, domain, enefunc, dynvars,      &
                               dynamics, pairlist, boundary, constraints,  &
                               ensemble, comm, remd, alchemy)
      end if

      ! perform remd
      !
      if (.not. remd%equilibration_only) then
        call perform_replica_exchange_fep(i, domain, enefunc, dynvars,         &
                                      ensemble, boundary, output, remd,    &
                                      pairlist, alchemy)
      end if

      ! output restart data
      !
      call output_remd(i, output, enefunc, domain, dynamics, dynvars, &
                       boundary, constraints, remd)


    end do

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_remd_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    perform_replica_exchange_fep
  !> @brief        perform replica exchange for FEP
  !! @authors      HO
  !! @param[in]    icycle      : replica exchange cycle number
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] output      : output information
  !! @param[inout] remd        : REMD information
  !! @param[inout] pairlist    : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine perform_replica_exchange_fep(icycle, domain, enefunc, dynvars, &
                                      ensemble, boundary, output, remd, &
                                      pairlist, alchemy)

    ! formal arguments
    integer,                  intent(in)    :: icycle
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_alchemy),optional, intent(inout) :: alchemy

    ! local variables
    integer                   :: exptrn, dimno, itmp, parmsetid
    integer                   :: i, j, jx, k, neighid, id_old
    integer                   :: umbr_0, funcid
    integer                   :: ipara, ipara_ref
    real(wp)                  :: dpara, dpara_ref
    character                 :: ctrep*3, frmt1*12, frmt2*12, aorr*1
    integer,          pointer :: total_nreplicas, iseed
    integer,          pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,          pointer :: repid2parmsetid_ref(:)
    integer,          pointer :: iparameters(:,:)
    real(wp),         pointer :: dparameters(:,:)

    integer, save, allocatable :: counter_skip(:)
    integer, save, allocatable :: counter_tuning(:)

    real(wp) :: temp_i, temp_j, factor


    iseed               => remd%iseed
    total_nreplicas     => remd%total_nreplicas
    repid2parmsetid     => remd%repid2parmsetid
    repid2parmsetid_ref => remd%repid2parmsetid_ref
    parmsetid2repid     => remd%parmsetid2repid
    parmidsets          => remd%parmidsets
    dparameters         => remd%dparameters
    iparameters         => remd%iparameters
    num_criteria        => remd%num_criteria
    num_exchanges       => remd%num_exchanges


    ! decide dimension to be exchanged and exchange pattern
    !   pattern 1: 1<>2  3<>4  5<>6  7<>8 ...
    !   pattern 2: 1  2<>3  4<>5  6<>7  8 ...
    !
    itmp  = mod(icycle,2*remd%dimension)
    dimno = mod(itmp,remd%dimension)+1
    if (itmp+1 <= remd%dimension) then
      exptrn = 1
    else
      exptrn = 2
    end if

    ! parameter tuning
    !
    if (remd%autoadj(dimno)%param_tuning) then
      if (.not. allocated(counter_skip)) then
        allocate(counter_skip(remd%dimension))
        counter_skip(1:remd%dimension) = 0
      end if
      if (.not. allocated(counter_tuning)) then
        allocate(counter_tuning(remd%dimension))
        counter_tuning(1:remd%dimension) = 0
      end if
      counter_skip(dimno) = counter_skip(dimno) + 1
      if (counter_skip(dimno) <= 2 * remd%autoadj(dimno)%eq_cycle) then
        if (main_rank) then
          write(MsgOut,'(A, I12)') &
            "Info> skip exchange at step", dynvars%step
          write(MsgOut,*)
        end if
        return
      end if
    end if


    ! save previous information
    !
    repid2parmsetid_ref(1:total_nreplicas) = repid2parmsetid(1:total_nreplicas)


    ! generate random number
    !
    do i = 1, total_nreplicas-1
      do j = i+1, total_nreplicas
        remd%random_table(i,j) = random_get_legacy(iseed)
        remd%random_table(j,i) = remd%random_table(i,j)
      end do
    end do

    
    ! perform replica exchange
    !
    select case(remd%types(dimno))

    case(RemdTemperature)

      call temperature_remd(dimno, exptrn, domain, dynvars, remd)

    case(RemdPressure)

      call pressure_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdGamma)

      call surface_tension_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdRestraint)

      call restraint_remd(dimno, exptrn, domain, boundary, ensemble, dynvars, &
                          enefunc, remd)

    case(RemdSoluteTempering)

      call solute_tempering_remd(dimno, exptrn, domain, ensemble, dynvars, &
                                 enefunc, remd, pairlist, boundary, output)

    ! FEP or FEP-gREST
    case(RemdAlchemy, RemdAlchemyRest)

      if (present(alchemy)) then
        call fep_remd(dimno, exptrn, domain, ensemble, dynvars, enefunc, remd, &
                      pairlist, boundary, output, alchemy)
      else
        call error_msg('Alchemy infomartion must be needed')
      end if

    end select

    ! parameter tuning
    !
    if (remd%autoadj(dimno)%param_tuning) then
      counter_tuning(dimno) = counter_tuning(dimno) + 1
      if (counter_tuning(dimno) >= 2 * remd%autoadj(dimno)%trial_freq) then

        ! in case of temperature REMD, we should remember current temp
        if (remd%types(dimno) == RemdTemperature) then
          parmsetid = repid2parmsetid(my_country_no+1)
          temp_i = remd%dparameters(dimno,parmidsets(parmsetid,dimno))
        end if

        ! check current exchange ratio
        call remd_autoadj_displacement(dimno,remd)

        if (remd%types(dimno) == RemdTemperature) then
          temp_j = remd%dparameters(dimno,parmidsets(parmsetid,dimno))
          factor = sqrt(temp_j/temp_i)
          do j = 1, domain%num_cell_local
            do jx = 1, domain%num_atom(j)
              domain%velocity(1:3,jx,j) = domain%velocity(1:3,jx,j) * factor
            end do
          end do

          factor = temp_j/temp_i
          dynvars%thermostat_momentum                 &
            = dynvars%thermostat_momentum * factor
          dynvars%barostat_momentum(1:3)              &
            = dynvars%barostat_momentum(1:3) * factor

        end if

        ! reset counters
        do i = 1, remd%total_nreplicas
          remd%num_criteria (i,dimno,1:2) = 0
          remd%num_exchanges(i,dimno,1:2) = 0
        end do

        if (main_rank) then
          write(MsgOut,'("REMD> New parameter set: ")',advance="NO")
          if (remd%types(dimno) /= RemdRestraint) then
            do i = 1, remd%nreplicas(dimno)
              write(MsgOut,'(2X,F10.5)',advance="NO") remd%dparameters(dimno,i)
            end do
          else
            umbr_0 = remd%iparameters(dimno,parmidsets(1,dimno))
            do k = 1, remd%umbrid2numfuncs(umbr_0)
              funcid = remd%umbrid2funclist(umbr_0,k)
              do i = 1, remd%nreplicas(dimno)
                write(MsgOut,'(2X,F10.5)',advance="NO") remd%rest_reference(i,k)
              end do
            end do
          end if
          write(MsgOut,*)
          write(MsgOut,*)
        end if
        counter_skip(dimno)   = 0
        counter_tuning(dimno) = 0
      end if
    end if


    ! update permutation function
    !
    do i = 1, total_nreplicas
      parmsetid = repid2parmsetid(i)
      parmsetid2repid(parmsetid) = i
    end do


    ! write results
    !
    if (main_rank) then

      write(MsgOut,'(A,I10,3X,A,I4,3X,A,I4)')                           &
        'REMD> Step: ',dynvars%step,'Dimension: ',dimno,                &
        'ExchangePattern: ',exptrn

      write(MsgOut,'(A,6X,A,13X,A,6X,A,7X,A)')     &
        '  Replica', 'ExchangeTrial', 'AcceptanceRatio', 'Before', 'After'
      do i = 1, total_nreplicas
        id_old = repid2parmsetid_ref(i)
        call get_neighbouring_parmsetid(dimno, exptrn, remd, id_old, neighid)

        if ((remd%types(dimno) /= RemdRestraint) .and. &
            (remd%types(dimno) /= RemdAlchemy)) then

          dpara_ref = dparameters(dimno,parmidsets(repid2parmsetid_ref(i),dimno))
          dpara     = dparameters(dimno,parmidsets(repid2parmsetid(i),    dimno))

          if (abs(dpara_ref - dpara) < EPS) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2F12.3)')         &
               i, id_old, ' > ', neighid, aorr,                           &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ', &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),        &
               dpara_ref, dpara
        else

          ipara_ref = iparameters(dimno,parmidsets(repid2parmsetid_ref(i),dimno))
          ipara     = iparameters(dimno,parmidsets(repid2parmsetid(i),    dimno))

          if (ipara_ref == ipara) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2I12)')           &
               i, id_old, ' > ', neighid, aorr,                           &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ', &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),        &
               ipara_ref, ipara
        end if
      end do

      write(MsgOut,*) ''

      write(ctrep,'(I3)') total_nreplicas
      frmt1 = '(A,' // trim(adjustl(ctrep)) // 'F10.3)'
      frmt2 = '(A,' // trim(adjustl(ctrep)) // 'I10  )'

      if ((remd%types(dimno) /= RemdRestraint) .and. &
          (remd%types(dimno) /= RemdAlchemy)) then
        write(MsgOut,frmt1) '  Parameter    : ', &
          dparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),dimno))
      else
        write(MsgOut,frmt2) '  Parameter    : ', &
          iparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),dimno))
      end if
      write(MsgOut,frmt2)   '  RepIDtoParmID: ', repid2parmsetid(1:total_nreplicas)
      write(MsgOut,frmt2)   '  ParmIDtoRepID: ', parmsetid2repid(1:total_nreplicas)

      write(MsgOut,*) ''

    end if

    return

  end subroutine perform_replica_exchange_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fep_remd( rep_info, molecule, domain, &
                             remd, restraints, enefunc, alchemy )

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_domain),          intent(inout) :: domain
    type(s_remd),            intent(inout) :: remd
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    integer :: i, j, k, ndata, nrest
    integer :: ierror
    integer :: lambid, max_nreplicas

    integer, parameter :: max_param_type = 20
    character(MaxLine) :: param_type_str(max_param_type)


    ! allocation
    !
    max_nreplicas = maxval(rep_info%nreplicas(1:remd%dimension))
    call alloc_remd(remd, RemdFEP, remd%dimension, remd%total_nreplicas, max_nreplicas)

    ! make permutation function for FEP
    !
    lambid = 0
    do i = 1, remd%dimension
      if (remd%types(i) == RemdAlchemy .or. &
          remd%types(i) == RemdAlchemyRest) then
        if (alchemy%num_fep_windows /= (remd%nreplicas(i))) &
          call error_msg('Setup_Remd> &
            "nreplicas" and number of FEP windows are inconsistent')
        do j = 1, remd%nreplicas(i)
          lambid = lambid + 1
          remd%iparameters(i,j) = lambid
          remd%dlambljA(lambid) = alchemy%lambljA(j)
          remd%dlambljB(lambid) = alchemy%lambljB(j)
          remd%dlambelA(lambid) = alchemy%lambelA(j)
          remd%dlambelB(lambid) = alchemy%lambelB(j)
          remd%dlambbondA(lambid) = alchemy%lambbondA(j)
          remd%dlambbondB(lambid) = alchemy%lambbondB(j)
          remd%dlambrest(lambid) = alchemy%lambrest(j)
        end do
      end if

      if (remd%types(i) == RemdAlchemyRest) then

        ! read paramter type
        ndata = split_num(rep_info%param_type(i))
        call split(ndata,max_param_type,rep_info%param_type(i), &
                   param_type_str)

        remd%fep_rest(i)%sw_bonds        = .false.
        remd%fep_rest(i)%sw_angles       = .false.
        remd%fep_rest(i)%sw_ureys        = .false.
        remd%fep_rest(i)%sw_dihedrals    = .false.
        remd%fep_rest(i)%sw_impropers    = .false.
        remd%fep_rest(i)%sw_cmaps        = .false.
        remd%fep_rest(i)%sw_lj           = .false.
        remd%fep_rest(i)%sw_charge       = .false.

        do j = 1, ndata
          call tolower(param_type_str(j))
          select case (trim(param_type_str(j)))
            case("all")
              remd%fep_rest(i)%sw_bonds        = .true.
              remd%fep_rest(i)%sw_angles       = .true.
              remd%fep_rest(i)%sw_ureys        = .true.
              remd%fep_rest(i)%sw_dihedrals    = .true.
              remd%fep_rest(i)%sw_impropers    = .true.
              remd%fep_rest(i)%sw_cmaps        = .true.
              remd%fep_rest(i)%sw_lj           = .true.
              remd%fep_rest(i)%sw_charge       = .true.
            case("b", "bond", "bonds")
              remd%fep_rest(i)%sw_bonds        = .true.
            case("a", "angle", "angles")
              remd%fep_rest(i)%sw_angles       = .true.
            case("u", "urey", "ureys")
              remd%fep_rest(i)%sw_ureys        = .true.
            case("d", "dihedral", "dihedrals")
              remd%fep_rest(i)%sw_dihedrals    = .true.
            case("i", "improper", "impropers")
              remd%fep_rest(i)%sw_impropers    = .true.
            case("cm", "cmap", "cmaps")
              remd%fep_rest(i)%sw_cmaps        = .true.
            case("c", "charge", "charges")
              remd%fep_rest(i)%sw_charge       = .true.
            case("l", "lj", "ljs")
              remd%fep_rest(i)%sw_lj           = .true.
          end select
        end do

        ! this value is set to true after the first call of assign_condition
        !
        remd%fep_rest(i)%done_setup = .false.

        remd%fep_rest(i)%num_solute = 0
        remd%fep_rest(i)%num_atom_cls = 0

        read(rep_info%select_index(i),*) remd%fep_rest(i)%mygroup

        nrest = restraints%num_atoms(remd%fep_rest(i)%mygroup)
        remd%fep_rest(i)%num_solute = nrest

        allocate(remd%fep_rest(i)%solute_list(nrest), &
                 remd%fep_rest(i)%is_solute(molecule%num_atoms), &
                 remd%fep_rest(i)%atom_cls_no_org(nrest))
        remd%fep_rest(i)%is_solute(1:molecule%num_atoms) = 0

        ! listup solute atoms
        !
        do j = 1, restraints%num_atoms(remd%fep_rest(i)%mygroup)
          k = restraints%atomlist(j,remd%fep_rest(i)%mygroup)
          remd%fep_rest(i)%solute_list(j) = k
          remd%fep_rest(i)%is_solute(k) = 1
        end do

        if (main_rank) then

          write(MsgOut,'(a)')    'Setup_Remd_Solute_Tempering>'
          write(MsgOut,'(a,i8)') 'Param types for dimension', i
          write(MsgOut,'(a,l8)') '  BONDS         = ', &
                                 remd%fep_rest(i)%sw_bonds
          write(MsgOut,'(a,l8)') '  ANGLES        = ', &
                                 remd%fep_rest(i)%sw_angles
          write(MsgOut,'(a,l8)') '  UREYS         = ', &
                                 remd%fep_rest(i)%sw_ureys
          write(MsgOut,'(a,l8)') '  DIHEDRALS     = ', &
                                 remd%fep_rest(i)%sw_dihedrals
          write(MsgOut,'(a,l8)') '  IMPROPERS     = ', &
                                 remd%fep_rest(i)%sw_impropers
          write(MsgOut,'(a,l8)') '  CMAPS         = ', &
                                 remd%fep_rest(i)%sw_cmaps
          write(MsgOut,'(a,l8)') '  CHARGE        = ', &
                                 remd%fep_rest(i)%sw_charge
          write(MsgOut,'(a,l8)') '  LJ            = ', &
                                 remd%fep_rest(i)%sw_lj
        end if

        if (main_rank) then
          write(MsgOut,'(a,i8)') 'Solute Atoms for dimension', i
          write(MsgOut,'(a,i0)') '  num_solute              = ', nrest
        end if

        ! bond (just count number)
        !
        call setup_remd_solute_tempering_bonds( remd%fep_rest(i), &
                                                domain, enefunc )
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_bonds, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_bonds
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_bond         = ', nrest
        end if

        ! angle (just count number)
        !
        call setup_remd_solute_tempering_angles( remd%fep_rest(i), &
                                                 domain, enefunc )
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_angles, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_angles
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_angles       = ', nrest
        end if
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_ureys, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_ureys
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_ureys        = ', nrest
        end if

        ! dihedral (just count number)
        !
        call setup_remd_solute_tempering_dihedrals( remd%fep_rest(i), &
                                                    domain, enefunc )
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_dihedrals
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_dihedrals    = ', nrest
        end if

        ! improper (just count number)
        !
        call setup_remd_solute_tempering_impropers( remd%fep_rest(i), &
                                                    domain, enefunc )
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_impropers, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_impropers
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_impropers    = ', nrest
        end if

        ! cmap (just count number)
        !
        call setup_remd_solute_tempering_cmaps( remd%fep_rest(i), &
                                                domain, enefunc )
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%fep_rest(i)%num_cmaps, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%fep_rest(i)%num_cmaps
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_cmaps        = ', nrest
        end if

        ! LJ
        !
        call setup_remd_solute_tempering_lj( remd%fep_rest(i), &
                                             domain, enefunc )

        ! charges are not modified here; to be done in assign_condition
        !

      end if
    end do

    return

  end subroutine setup_fep_remd

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fep_remd(dimno, exptrn, domain, ensemble, dynvars, enefunc, remd, &
                      pairlist, boundary, output, alchemy)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_domain),   target, intent(inout) :: domain
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output
    type(s_alchemy),          intent(inout) :: alchemy

    ! local variables
    integer                  :: i, j, k, ix
    integer                  :: natoms
    integer                  :: repid_i, repid_j
    integer                  :: parmsetid, parmsetid_i, parmsetid_j
    integer                  :: neighid
    integer                  :: funcid, umbr_0, umbr_1
    real(wp)                 :: cv
    real(wp)                 :: energy_i_0, energy_j_0
    real(wp)                 :: energy_i_1, energy_j_1
    real(wp)                 :: energy_i, energy_j
    real(wp)                 :: eneval, beta, delta, temperature
    real(wp)                 :: before_gather
    logical                  :: does_exist, do_trial, do_exchange, umbr_exist
    integer,         pointer :: total_nreplicas
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,         pointer :: parmidsets(:,:)
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,         pointer :: parameters(:,:)
    real(wp),        pointer :: after_gather(:), random_table(:,:)
    real(wp),        pointer :: refcoord0(:,:), refcoord1(:,:)

    after_gather    => remd%after_gather
    random_table    => remd%random_table
    total_nreplicas => remd%total_nreplicas
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid
    parmidsets      => remd%parmidsets
    parameters      => remd%iparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    refcoord0       => remd%restraint_refcoord0
    refcoord1       => remd%restraint_refcoord1

    ! allgather energy difference of forward
    !
    eneval = dynvars%energy%deltU_fep(3)
    call mpi_bcast(eneval, 1, mpi_wp_real, 0, mpi_comm_country, ierror)

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%deltU_fwd(i) = after_gather(i)
    end do

    ! allgather energy difference of reverse
    !
    eneval = dynvars%energy%deltU_fep(2)
    call mpi_bcast(eneval, 1, mpi_wp_real, 0, mpi_comm_country, ierror)

    before_gather = eneval
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror)
#endif
    do i = 1, total_nreplicas
      remd%deltU_rev(i) = after_gather(i)
    end do

    ! replica exchange
    !
    do i = 1, total_nreplicas

      repid_i     = i
      parmsetid_i = repid2parmsetid(repid_i)

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(dimno, exptrn, remd, parmsetid_i, neighid)

      if (neighid /= 0) then
        parmsetid_j = neighid
        repid_j     = parmsetid2repid(parmsetid_j)
        do_trial    = .true.
      else
        do_trial    = .false.
      end if

      if (do_trial) then

        num_criteria(parmsetid_i,dimno,exptrn) &
          = num_criteria(parmsetid_i,dimno,exptrn) + 1

        ! compute delta
        !
        call get_replica_temperature(parmsetid_i, remd, ensemble, temperature)

        beta = 1.0_wp / (KBOLTZ * temperature)
        if (parmsetid_i < parmsetid_j) then
          delta = beta*(remd%deltU_fwd(repid_i) + remd%deltU_rev(repid_j))
        else
          delta = beta*(remd%deltU_fwd(repid_j) + remd%deltU_rev(repid_i))
        end if

        ! Metropolis criterion
        !
        if (delta <= 0.0_wp) then
          do_exchange = .true.
        else if (delta > 0.0_wp) then
          if (random_table(repid_i,repid_j) <= exp(-delta)) then
            do_exchange = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (do_exchange) then
          num_exchanges(parmsetid_i,dimno,exptrn) &
            = num_exchanges(parmsetid_i,dimno,exptrn) + 1
          repid2parmsetid(i) = parmsetid_j
        else
          call assign_condition_fep(parmsetid_i, remd, domain, ensemble, &
                                    enefunc, alchemy)
        end if

      end if
    end do

    return

  end subroutine fep_remd

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_fep(parmsetid, remd, domain, ensemble, enefunc, &
                                  alchemy)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_remd),            intent(inout) :: remd
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_alchemy),         intent(in)    :: alchemy

    ! local variables
    integer                  :: i, j, k
    integer                  :: umbrid, funcid, atmid
    integer                  :: lambid
    real(wp)                 :: rest_param

    do i = 1, remd%dimension
      if      (remd%types(i) == RemdTemperature) then
        ensemble%temperature = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdPressure   ) then
        ensemble%pressure    = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdGamma      ) then
        ensemble%gamma       = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdRestraint  ) then
        umbrid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
        do k = 1, remd%umbrid2numfuncs(umbrid)
          funcid = remd%umbrid2funclist(umbrid,k)
          enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbrid,k)
          enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbrid,k)

          if (enefunc%restraint_kind(funcid) == RestraintsFuncRMSD) then
            enefunc%rmsd_force  = enefunc%restraint_const(1,funcid)
            enefunc%rmsd_target = enefunc%restraint_ref  (1,funcid)
          end if

        end do

      else if (remd%types(i) == RemdAlchemy     .or. &
               remd%types(i) == RemdAlchemyRest) then
        ! FEP
        lambid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
        enefunc%lambljA   = remd%dlambljA(lambid)
        enefunc%lambljB   = remd%dlambljB(lambid)
        enefunc%lambelA   = remd%dlambelA(lambid)
        enefunc%lambelB   = remd%dlambelB(lambid)
        enefunc%lambbondA = remd%dlambbondA(lambid)
        enefunc%lambbondB = remd%dlambbondB(lambid)
        enefunc%lambrest  = remd%dlambrest(lambid)

        ! FEP: REST2-like scaling
        call assign_lambda(alchemy, domain, enefunc)

      end if

      if (remd%types(i) == RemdAlchemyRest) then

        if (ensemble%ensemble == EnsembleNVE ) then
          call error_msg('Solute_Tempering> REST is not available for NVE!')
        end if

        rest_param = remd%dparameters(i,remd%parmidsets(parmsetid,i))
        call assign_condition_feprest( remd%fep_rest(i), &
                                       rest_param, domain, &
                                       ensemble, enefunc )

      end if

    end do

    return

  end subroutine assign_condition_fep

end module sp_remd_mod
