!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_remd_mod
!> @brief   Replica exhcnage molecular dynamics simulation
!! @authors Takaharu Mori (TM), Norio Takase (NT), Motoshi Kamiya (MK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_remd_mod

  use at_energy_mod
  use at_energy_str_mod
  use at_energy_restraints_mod
  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_md_vverlet_cg_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_remd_str_mod
  use at_restraints_str_mod
  use at_constraints_str_mod
  use at_constraints_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_enefunc_restraints_mod
  use at_enefunc_gbsa_mod
  use at_output_str_mod
  use at_output_mod
  use molecules_str_mod
  use fileio_rst_mod
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
  type, public :: s_rep_info
    integer                           :: dimension          = 1
    integer                           :: exchange_period    = 100
    integer                           :: iseed              = 3141592
    integer,              allocatable :: types        (:)
    integer,              allocatable :: nreplicas    (:) 
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
  public  :: setup_remd
  public  :: run_remd
  private :: setup_reus
  private :: setup_remd_solute_tempering
  private :: perform_replica_exchange
  private :: remd_autoadj_displacement
  private :: temperature_remd
  private :: pressure_remd
  private :: surface_tension_remd
  private :: restraint_remd
  private :: solute_tempering_remd
  public  :: assign_condition
  private :: assign_condition_rest
  private :: assign_condition_solute_tempering
  private :: assign_condition_solute_tempering_lj
  private :: get_neighbouring_parmsetid
  private :: get_replica_temperature

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

    ! error check
    !
    if (rep_info%dimension > RemdMaxNdimensions .or. rep_info%dimension == 0) &
      call error_msg('Read_Ctrl_Remd> error in dimension in [REMD]')

    ! read allocatable variables
    !
    allocate(rep_info%types        (rep_info%dimension),   &
             rep_info%nreplicas    (rep_info%dimension),   &
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

    call end_ctrlfile_section(handle)


    ! error check
    !
    do i = 1, rep_info%dimension
      if (rep_info%nreplicas(i) == 0) &
        call error_msg('Read_Ctrl_Remd> error in nreplicas in [REMD]')
    end do

    trep = product(rep_info%nreplicas(1:rep_info%dimension))
    if (trep == 0) &
      call error_msg('Read_Ctrl_Remd> error in nreplicas in [REMD]')


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

  subroutine setup_solute_tempering(rep_info, molecule, restraints, &
                                    cons_info)

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
      if (rep_info%types(i) == RemdSoluteTempering) then
        has_rest = .true.
      end if
    end do
    if (.not. has_rest) return

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Solute_Tempering> preparation of REST'
    end if

    do i = 1, rep_info%dimension

      if (rep_info%types(i) == RemdSoluteTempering) then

        ndata = split_num(rep_info%select_index(i))

        if (ndata /= 1) &
          call error_msg('Setup_Solute_Tempering> Error in index for REST')

        read(rep_info%select_index(i),*) ndata

        if (restraints%num_groups < ndata) then
          call error_msg('  Setup_Solute_Tempering> group out of range')
        end if

        ! have to change water name to avoid table for "solute"
        !
        do j = 1, restraints%num_atoms(ndata)

          k = restraints%atomlist(j,ndata)
          if (molecule%residue_name(k) == cons_info%water_model) then
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
  !  Subroutine    setup_remd
  !> @brief        setup REMD
  !! @authors      TM
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] restraints : restraints information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd(rep_info, rst, boundary, dynamics, molecule, &
                        restraints, ensemble, enefunc, remd)

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_rst),             intent(in)    :: rst
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_molecule),        intent(inout) :: molecule
    type(s_restraints),      intent(inout) :: restraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k, m
    integer                  :: max_nreplicas, cycles
    integer                  :: itmp, jtmp, ktmp
    integer                  :: comp, dsta, dend, found, dtotal
    integer                  :: replicaid, parmsetid
    integer                  :: ncycle, icycle, nlen, ixx
    integer                  :: funcid
    integer                  :: iref, icomp
    character(MaxLine)       :: param
    character(10)            :: tmp, frmt, cd
    character(18)            :: frmt2
    character(3)             :: ctrep
    character(MaxLine)       :: ft
    logical                  :: keep, double

    real(wp),      allocatable :: ddata(:)
    integer,       allocatable :: ndigit(:), parmid(:)
    integer,       allocatable :: sendbuf(:), recvbuf(:,:)


    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Remd> Replica information'
      write(MsgOut,'(A)') ''
    end if


    ! initialize
    !
    remd%dimension       = rep_info%dimension
    remd%total_nreplicas = product(rep_info%nreplicas(1:remd%dimension))
    max_nreplicas        = maxval (rep_info%nreplicas(1:remd%dimension))

    remd%rest_mixed_three_atoms = .false.

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
          call error_msg('Setup_Remd> nreplicas must be even,'//           &
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
      if (dynamics%rstout_period > 0) then
        if (mod(dynamics%rstout_period,rep_info%exchange_period) /= 0) then
          call error_msg('Setup_Remd> mod(rstout_period, exchange_period)'//&
                         ' must be zero)')
        end if
      end if
      if (dynamics%eneout_period > rep_info%exchange_period) then
        call error_msg('Setup_Remd> (eneout_period <= exchange_period)'//    &
                       ' is required')
      else
        if (mod(rep_info%exchange_period,dynamics%eneout_period) /= 0) then
          call error_msg('Setup_Remd> mod(exchange_period,eneout_period)'//  &
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
    call setup_reus(rep_info, rst, restraints, enefunc, remd)

    ! setup replica exchange with solute tempering (REST)
    !
    call setup_remd_solute_tempering(rep_info, molecule, remd, restraints, &
                                     enefunc)

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
    call assign_condition(parmsetid, remd, ensemble, molecule, enefunc)


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
        do k = 1, remd%umbrid2numfuncs(i)
          write(MsgOut,'(6X,A,I7,2X,A,F11.3,2X,A,F11.3)')          &
                       'function  = ', remd%umbrid2funclist(i,k),  &
                       'constant  = ', remd%rest_constants(i,k),   &
                       'reference = ', remd%rest_reference(i,k)
          funcid = remd%umbrid2funclist(i,k)
          if (restraints%function(funcid) == RestraintsFuncRG .or. &
              restraints%function(funcid) == RestraintsFuncRGWOMASS) then
            write(MsgOut,'(6X,A,L7,2X,A,F11.3)')                     &
                         'caging    = ', remd%rest_caging(i,k),      &
                         'flat_radi = ', remd%rest_flat_radius(i,k)
          end if
        end do
        write(MsgOut,'(A)') ''
      end do

    end if

    return

  end subroutine setup_remd

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
    integer                  :: i, j, k, ifound
    integer                  :: ndata, max_ndata
    integer                  :: umbrid, funcid, max_umbrella_id, max_nreplicas
    character(MaxLine)       :: param
    character(10)            :: msg1, msg2

    integer,       allocatable :: idata(:)
    real(wp),      allocatable :: ddata(:)
    character(6),  allocatable :: cdata(:)


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
    allocate(cdata(max_nreplicas))

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

            ! caging and flat_radius (Rg restraint)
            if (restraints%function(funcid) == RestraintsFuncRG .or. &
                restraints%function(funcid) == RestraintsFuncRGWOMASS) then
              param = restraints%caging(funcid)
              ndata = split_num(param)
              if (ndata > 0) then
                call split(ndata, max_nreplicas, param, cdata)
                call tolower(cdata(j))
                remd%rest_caging(umbrid,k) = (cdata(j) == 'yes' .or. cdata(j) == 'true')
              end if

              param = restraints%flat_radius(funcid)
              ndata = split_num(param)
              if (ndata > 0) then
                call split(ndata, max_nreplicas, param, ddata)
                remd%rest_flat_radius(umbrid,k) = ddata(j)
              end if
            end if

          end do
        end do
      end if
    end do

    deallocate(ddata)
    deallocate(cdata)

    ! replace rest_reference values with ones in rst (for REUS after RPATH)
    !
    if ((remd%dimension == 1) .and. (remd%types(1) == RemdRestraint) .and. &
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
  !> @brief        main routine for setupping REST
  !! @authors      MK
  !! @param[in]    rep_info   : REMD section control parameters information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] remd       : REMD information
  !! @param[inout] restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering(rep_info, molecule, remd, &
                                         restraints, enefunc )

    ! formal arguments
    type(s_rep_info),        intent(in)    :: rep_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_remd),            intent(inout) :: remd
    type(s_restraints),      intent(inout) :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer :: i, j, k, num, ndata, nrest, nflags
    integer :: ierror

    integer, parameter :: max_param_type = 20
    character(MaxLine) :: param_type_str(max_param_type)


    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then

        ! read parameter type
        ndata = split_num(rep_info%param_type(i))
        call split(ndata,max_param_type,rep_info%param_type(i), &
                   param_type_str)

        remd%solute_tempering(i)%sw_charge    = .false.
        remd%solute_tempering(i)%sw_bonds     = .false.
        remd%solute_tempering(i)%sw_angles    = .false.
        remd%solute_tempering(i)%sw_ureys     = .false.
        remd%solute_tempering(i)%sw_dihedrals = .false.
        remd%solute_tempering(i)%sw_impropers = .false.
        remd%solute_tempering(i)%sw_cmaps     = .false.
        remd%solute_tempering(i)%sw_lj        = .false.
        remd%solute_tempering(i)%sw_contacts  = .false.

        do j = 1, ndata
          call tolower(param_type_str(j))
          select case (trim(param_type_str(j)))
            case("all")
              remd%solute_tempering(i)%sw_charge    = .true.
              remd%solute_tempering(i)%sw_bonds     = .true.
              remd%solute_tempering(i)%sw_angles    = .true.
              remd%solute_tempering(i)%sw_ureys     = .true.
              remd%solute_tempering(i)%sw_dihedrals = .true.
              remd%solute_tempering(i)%sw_impropers = .true.
              remd%solute_tempering(i)%sw_cmaps     = .true.
              remd%solute_tempering(i)%sw_lj        = .true.
              remd%solute_tempering(i)%sw_contacts  = .true.
            case("b", "bond", "bonds")
              remd%solute_tempering(i)%sw_bonds     = .true.
            case("a", "angle", "angles")
              remd%solute_tempering(i)%sw_angles    = .true.
            case("u", "urey", "ureys")
              remd%solute_tempering(i)%sw_ureys     = .true.
            case("d", "dihedral", "dihedrals")
              remd%solute_tempering(i)%sw_dihedrals = .true.
            case("i", "improper", "impropers")
              remd%solute_tempering(i)%sw_impropers = .true.
            case("cm", "cmap", "cmaps")
              remd%solute_tempering(i)%sw_cmaps     = .true.
            case("con", "contact", "contacts")
              remd%solute_tempering(i)%sw_contacts  = .true.
            case("c", "charge", "charges")
              remd%solute_tempering(i)%sw_charge    = .true.
            case("l", "lj", "ljs")
              remd%solute_tempering(i)%sw_lj        = .true.
            case default
              if (main_rank) then
                write(MsgOut,*) 'unknown parameter type:', &
                                trim(param_type_str(j))
              end if
              call error_msg('Setup_solute_tempering> error. plz fix input')
          end select
        end do

        remd%solute_tempering(i)%num_solute    = 0
        remd%solute_tempering(i)%num_bonds     = 0
        remd%solute_tempering(i)%num_angles    = 0
        remd%solute_tempering(i)%num_ureys     = 0
        remd%solute_tempering(i)%num_dihedrals = 0
        remd%solute_tempering(i)%num_impropers = 0
        remd%solute_tempering(i)%num_cmaps     = 0
        remd%solute_tempering(i)%num_cmap_type = 0
        remd%solute_tempering(i)%num_atom_cls  = 0
        remd%solute_tempering(i)%num_contacts  = 0
        remd%solute_tempering(i)%num_mixed_three = 0

        read(rep_info%select_index(i),*) remd%solute_tempering(i)%mygroup

        nrest = restraints%num_atoms(remd%solute_tempering(i)%mygroup)
        remd%solute_tempering(i)%num_solute = nrest
        allocate(remd%solute_tempering(i)%solute_list(nrest), &
                 remd%solute_tempering(i)%is_solute(molecule%num_atoms), &
                 remd%solute_tempering(i)%atom_cls_no_org(nrest))
        remd%solute_tempering(i)%is_solute(1:molecule%num_atoms) = .false.

        ! listup solute atoms
        !
        do j = 1, restraints%num_atoms(remd%solute_tempering(i)%mygroup)
          k = restraints%atomlist(j,remd%solute_tempering(i)%mygroup)
          remd%solute_tempering(i)%solute_list(j) = k
          remd%solute_tempering(i)%is_solute(k) = .true.
        end do

        if (main_rank) then
          write(MsgOut,'(a)')    'Setup_Remd_Solute_Tempering>'
          write(MsgOut,'(a,i8)') 'Param types for dimension', i
          write(MsgOut,'(a,l8)') '  BONDS      = ', &
                                 remd%solute_tempering(i)%sw_bonds
          write(MsgOut,'(a,l8)') '  ANGLES     = ', &
                                 remd%solute_tempering(i)%sw_angles
          write(MsgOut,'(a,l8)') '  UREYS      = ', &
                                 remd%solute_tempering(i)%sw_ureys
          write(MsgOut,'(a,l8)') '  DIHEDRALS  = ', &
                                 remd%solute_tempering(i)%sw_dihedrals
          write(MsgOut,'(a,l8)') '  IMPROPERS  = ', &
                                 remd%solute_tempering(i)%sw_impropers
          write(MsgOut,'(a,l8)') '  CMAPS      = ', &
                                 remd%solute_tempering(i)%sw_cmaps
          write(MsgOut,'(a,l8)') '  CONTACTS   = ', &
                                 remd%solute_tempering(i)%sw_contacts
          write(MsgOut,'(a,l8)') '  CHARGES    = ', &
                                 remd%solute_tempering(i)%sw_charge
          write(MsgOut,'(a,l8)') '  LJ         = ', &
                                 remd%solute_tempering(i)%sw_lj
        end if

        if (main_rank) then
          write(MsgOut,'(a,i8)') 'Solute Atoms for dimension', i
          write(MsgOut,'(a,i0)') '  num_solute          = ', nrest

        end if

        ! bond
        !
        call setup_remd_solute_tempering_bonds(remd%solute_tempering(i), &
                                               enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_bonds, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_bonds
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_bond     = ', nrest
        end if

        ! angle
        !
        nflags = 0
        call setup_remd_solute_tempering_angles(remd%solute_tempering(i), &
                                                enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_angles, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        call mpi_reduce(remd%solute_tempering(i)%num_mixed_three, nflags, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_angles
        nflags = remd%solute_tempering(i)%num_mixed_three
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_angle    = ', nrest
        end if
        if (nflags > 0) remd%rest_mixed_three_atoms = .true.

        ! urey-bradley
        !
        call setup_remd_solute_tempering_ureys(remd%solute_tempering(i), &
                                               enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_ureys, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_ureys
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_urey     = ', nrest
        end if

        ! dihedrals
        !
        nflags = 0
        call setup_remd_solute_tempering_dihedrals(remd%solute_tempering(i), &
                                                   enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_dihedrals, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        call mpi_reduce(remd%solute_tempering(i)%num_mixed_three, nflags, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_dihedrals
        nflags = remd%solute_tempering(i)%num_mixed_three
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_dihedral = ', nrest
        end if
        if (nflags > 0) remd%rest_mixed_three_atoms = .true.

        ! impropers
        !
        nflags = 0
        call setup_remd_solute_tempering_impropers(remd%solute_tempering(i), &
                                                   enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_impropers, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        call mpi_reduce(remd%solute_tempering(i)%num_mixed_three, nflags, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_impropers
        nflags = remd%solute_tempering(i)%num_mixed_three
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_improper = ', nrest
        end if
        if (nflags > 0) remd%rest_mixed_three_atoms = .true.

        ! cmaps
        !
        nflags = 0
        call setup_remd_solute_tempering_cmaps(remd%solute_tempering(i), &
                                               enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_cmaps, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        call mpi_reduce(remd%solute_tempering(i)%num_mixed_three, nflags, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_cmaps
        nflags = remd%solute_tempering(i)%num_mixed_three
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_cmap     = ', nrest
        end if
        if (nflags > 0) remd%rest_mixed_three_atoms = .true.

        ! charge
        !
        call setup_remd_solute_tempering_charges(remd%solute_tempering(i), &
                                                 molecule)

        ! LJ
        !
        call setup_remd_solute_tempering_lj(remd%solute_tempering(i), &
                                            molecule, enefunc)
        nrest = remd%solute_tempering(i)%num_atom_cls

        ! do not gather number for atomcls
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_atomcls  = ', nrest
        end if

        ! native contact
        !
        call setup_remd_solute_tempering_native_contact( &
                          remd%solute_tempering(i), enefunc)
#ifdef HAVE_MPI_GENESIS
        call mpi_reduce(remd%solute_tempering(i)%num_contacts, nrest, &
                        1, mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
#else
        nrest = remd%solute_tempering(i)%num_contacts
#endif
        if (main_rank) then
          write(MsgOut,'(a,i0)') '  num_solute_contact  = ', nrest
          write(MsgOut,*)
        end if

      end if
    end do

    return

  end subroutine setup_remd_solute_tempering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_bonds
  !> @brief        REST setup for bonds
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_bonds(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 2
    integer            :: alist(num_unit), i, j, k, num, ndata
    logical            :: count_only


    if (.not. allocated(enefunc%bond_force_const) .or. &
        .not. soltemp%sw_bonds) return

    count_only = .true.
    do i = 1, 2
      ndata = 0
      do j = enefunc%istart_bond, enefunc%iend_bond
        num = 0
        alist(1:num_unit) = enefunc%bond_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k))) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          if (.not. count_only) then
            soltemp%bond_list(ndata) = j
            soltemp%bond_weight(ndata) = real(num,wp) / real(num_unit,wp)
            soltemp%bond_force_const_org(ndata) = enefunc%bond_force_const(j)
          end if
        end if

      end do

      if (count_only) then
        soltemp%num_bonds = ndata
        allocate(soltemp%bond_list(ndata), &
                 soltemp%bond_weight(ndata), &
                 soltemp%bond_force_const_org(ndata))
      end if

      count_only = .false.

    end do

    return

  end subroutine setup_remd_solute_tempering_bonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_angles
  !> @brief        REST setup for angles
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_angles(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 3
    integer            :: alist(num_unit), i, j, k, num, ndata, nflags
    logical            :: count_only


    if (.not. allocated(enefunc%angl_force_const) .or. &
        .not. soltemp%sw_angles) return

    count_only = .true.
    do i = 1, 2
      ndata = 0
      nflags = 0
      do j = enefunc%istart_angle, enefunc%iend_angle
        num = 0
        alist(1:num_unit) = enefunc%angl_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k)) ) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          if (.not. count_only) then
            soltemp%angl_list(ndata) = j
            soltemp%angl_weight(ndata) = real(num,wp) / real(num_unit,wp)
            soltemp%angl_force_const_org(ndata) = enefunc%angl_force_const(j)
          end if
          if (num /= num_unit) nflags = nflags + 1
        end if

      end do

      if (count_only) then
        soltemp%num_mixed_three = nflags
        soltemp%num_angles = ndata
        allocate(soltemp%angl_list(ndata), &
                 soltemp%angl_weight(ndata), &
                 soltemp%angl_force_const_org(ndata))
      end if
      count_only = .false.
    end do

    return

  end subroutine setup_remd_solute_tempering_angles

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_ureys
  !> @brief        REST setup for Urey-Bradley terms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_ureys(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 2
    integer            :: alist(num_unit), i, j, k, num, ndata
    logical            :: count_only


    if (.not. allocated(enefunc%urey_force_const) .or. &
        .not. soltemp%sw_ureys ) return

    count_only = .true.
    do i = 1, 2
      ndata = 0
      do j = enefunc%istart_urey, enefunc%iend_urey
        num = 0
        alist(1:num_unit) = enefunc%urey_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k)) ) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          if (.not. count_only) then
            soltemp%urey_list(ndata) = j
            soltemp%urey_weight(ndata) = real(num,wp) / real(num_unit,wp)
            soltemp%urey_force_const_org(ndata) = enefunc%urey_force_const(j)
          end if
        end if

      end do

      if (count_only) then
        soltemp%num_ureys = ndata
        allocate(soltemp%urey_list(ndata), &
                 soltemp%urey_weight(ndata), &
                 soltemp%urey_force_const_org(ndata))
      end if
      count_only = .false.
    end do

    return

  end subroutine setup_remd_solute_tempering_ureys

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_dihedrals
  !> @brief        REST setup for dihedral terms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_dihedrals(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 4
    integer            :: alist(num_unit), i, j, k, num, ndata, nflags
    logical            :: count_only


    if (.not. allocated(enefunc%dihe_force_const) .or. &
        .not. soltemp%sw_dihedrals ) return

    count_only = .true.
    do i = 1, 2
      ndata = 0
      nflags = 0
      do j = enefunc%istart_dihedral, enefunc%iend_dihedral
        num = 0
        alist(1:num_unit) = enefunc%dihe_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k)) ) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          if (.not. count_only) then
            soltemp%dihe_list(ndata) = j
            soltemp%dihe_weight(ndata) = real(num,wp) / real(num_unit,wp)
            soltemp%dihe_force_const_org(ndata) = enefunc%dihe_force_const(j)
          end if
          if (num /= num_unit) nflags = nflags + 1
        end if

      end do

      if (count_only) then
        soltemp%num_dihedrals = ndata
        soltemp%num_mixed_three = nflags
        allocate(soltemp%dihe_list(ndata), &
                 soltemp%dihe_weight(ndata), &
                 soltemp%dihe_force_const_org(ndata))
      end if
      count_only = .false.
    end do

    return

  end subroutine setup_remd_solute_tempering_dihedrals

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_impropers
  !> @brief        REST setup for improper terms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_impropers(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 4
    integer            :: alist(num_unit), i, j, k, num, ndata, nflags
    logical            :: count_only


    if (.not. allocated(enefunc%impr_force_const) .or. &
        .not. soltemp%sw_impropers ) return

    count_only = .true.
    do i = 1, 2
      ndata = 0
      nflags = 0
      do j = enefunc%istart_improper, enefunc%iend_improper
        num = 0
        alist(1:num_unit) = enefunc%impr_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k)) ) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          if (.not. count_only) then
            soltemp%impr_list(ndata) = j
            soltemp%impr_weight(ndata) = real(num,wp) / real(num_unit,wp)
            soltemp%impr_force_const_org(ndata) = enefunc%impr_force_const(j)
          end if
          if (num /= num_unit) nflags = nflags + 1
        end if

      end do

      if (count_only) then
        soltemp%num_impropers = ndata
        soltemp%num_mixed_three = nflags
        allocate(soltemp%impr_list(ndata), &
                 soltemp%impr_weight(ndata), &
                 soltemp%impr_force_const_org(ndata))
      end if
      count_only = .false.
    end do

    return

  end subroutine setup_remd_solute_tempering_impropers

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_cmaps
  !> @brief        REST setup for crossterms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_cmaps(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer, parameter :: num_unit = 8
    integer            :: alist(num_unit), i, j, k, num, ndata, nflags
    logical            :: count_only

    integer :: num_new_types, ct, ncmap_type
    integer,  allocatable :: flags(:,:)
    integer,  allocatable :: cmap_resolution_org(:)
    real(wp), allocatable :: cmap_coef_org(:,:,:,:,:)


    if (.not. allocated(enefunc%cmap_coef) .or. &
        .not. soltemp%sw_cmaps ) return

    ncmap_type = size(enefunc%cmap_coef(1,1,1,1,:))
    if (ncmap_type == 0) return

    ! need to expand cmap_coef
    allocate(flags(ncmap_type,8))
    flags(1:ncmap_type,1:8) = 0
    num_new_types = 0

    count_only = .true.
    do i = 1, 2
      ndata = 0
      nflags = 0
      do j = enefunc%istart_cmap, enefunc%iend_cmap

        num = 0
        alist(1:num_unit) = enefunc%cmap_list(1:num_unit,j)
        do k = 1, num_unit
          if (soltemp%is_solute(alist(k)) ) num = num + 1
        end do

        if (num > 0) then
          ndata = ndata + 1
          ct = enefunc%cmap_type(j)
          if (count_only) then
            if (flags(ct,num) == 0) then
              num_new_types = num_new_types + 1
              flags(ct,num) = num_new_types
            end if
          else
            enefunc%cmap_type(j) = ncmap_type + flags(ct,num)
          end if
          if (num /= num_unit) nflags = nflags + 1
        end if

      end do

      if (count_only) then
        soltemp%num_mixed_three = nflags
        if (ndata == 0 ) exit
        soltemp%num_cmaps = ndata
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
      count_only = .false.
    end do

    if (allocated(flags)) deallocate(flags)

    return

  end subroutine setup_remd_solute_tempering_cmaps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_charges
  !> @brief        REST setup for charges
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_charges(soltemp, molecule)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer :: i


    if (soltemp%num_solute <= 0 .or. &
         .not. allocated(molecule%charge) .or. &
         .not. soltemp%sw_charge ) return

    allocate(soltemp%charge_org(soltemp%num_solute))
    do i = 1, soltemp%num_solute
      soltemp%charge_org(i) = molecule%charge(soltemp%solute_list(i))
    end do

    return

  end subroutine setup_remd_solute_tempering_charges

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_lj
  !> @brief        REST setup for lj terms
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_lj(soltemp, molecule, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! temporal allocatable arrays
    real(wp),                allocatable   :: tmp_lj(:,:)
    integer,                 allocatable   :: lc_atom_cls(:)

    ! local variables
    integer :: i, localcount, oldcount, org, new
    integer :: newcount


    soltemp%num_atom_cls = 0
    if (.not. soltemp%sw_lj ) return

    allocate(lc_atom_cls(enefunc%num_atom_cls))
    lc_atom_cls(1:enefunc%num_atom_cls) = 0

    do i = 1, soltemp%num_solute
      lc_atom_cls(molecule%atom_cls_no(soltemp%solute_list(i))) = 1
    end do

    localcount = 0
    do i = 1, enefunc%num_atom_cls
      if (lc_atom_cls(i) > 0) then
        localcount = localcount + 1
        lc_atom_cls(i) = localcount
        soltemp%atom_cls_no_org(localcount) = i
      end if
    end do
    soltemp%num_atom_cls = localcount
    soltemp%istart_atom_cls = enefunc%num_atom_cls + 1

    do i = 1, soltemp%num_solute
      org = molecule%atom_cls_no(soltemp%solute_list(i))
      new = lc_atom_cls(org) + enefunc%num_atom_cls
      molecule%atom_cls_no(soltemp%solute_list(i)) = new
    end do

    oldcount = enefunc%num_atom_cls
    allocate(tmp_lj(oldcount,oldcount))
    enefunc%num_atom_cls = oldcount + localcount

    newcount = oldcount + soltemp%num_atom_cls
    if (allocated(enefunc%nb14_lj6 )) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nb14_lj6(1:oldcount,1:oldcount)
      deallocate(enefunc%nb14_lj6)
      allocate(enefunc%nb14_lj6(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nb14_lj6 )
    end if
    if (allocated(enefunc%nb14_lj10)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nb14_lj10(1:oldcount,1:oldcount)
      deallocate(enefunc%nb14_lj10)
      allocate(enefunc%nb14_lj10(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nb14_lj10 )
    end if
    if (allocated(enefunc%nb14_lj12)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nb14_lj12(1:oldcount,1:oldcount)
      deallocate(enefunc%nb14_lj12)
      allocate(enefunc%nb14_lj12(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nb14_lj12 )
    end if
    if (allocated(enefunc%nonb_lj6)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nonb_lj6(1:oldcount,1:oldcount)
      deallocate(enefunc%nonb_lj6)
      allocate(enefunc%nonb_lj6(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nonb_lj6 )
    end if
    if (allocated(enefunc%nonb_lj10)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nonb_lj10(1:oldcount,1:oldcount)
      deallocate(enefunc%nonb_lj10)
      allocate(enefunc%nonb_lj10(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nonb_lj10 )
    end if
    if (allocated(enefunc%nonb_lj12)) then
      tmp_lj(1:oldcount,1:oldcount) = enefunc%nonb_lj12(1:oldcount,1:oldcount)
      deallocate(enefunc%nonb_lj12)
      allocate(enefunc%nonb_lj12(newcount,newcount))
      call setup_remd_solute_tempering_lj_each( &
                  oldcount, newcount, soltemp, tmp_lj, enefunc%nonb_lj12 )
    end if

    deallocate(tmp_lj, lc_atom_cls)

    return

  end subroutine setup_remd_solute_tempering_lj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_remd_solute_tempering_lj_each
  !> @brief        REST setup for each LJ term such as LJ6, LJ12
  !! @authors      MK
  !! @param[in]    ljsize     : original array size of interaction matrix
  !! @param[inout] soltemp    : REST information
  !! @param[inout] tmp_lj     : work matrix
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
  !  Subroutine    setup_remd_solute_tempering_native_contact
  !> @brief        REST setup for native contact interactions
  !! @authors      MK
  !! @param[inout] soltemp    : REST information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_remd_solute_tempering_native_contact(soltemp, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    type(s_enefunc),         intent(in)    :: enefunc

    ! local variables
    integer, parameter :: num_unit = 2
    integer :: alist(num_unit), i, j, k, num, ndata
    logical :: count_only


    if (enefunc%num_contacts <= 0 .or. &
         .not. soltemp%sw_contacts ) return

    if (allocated(enefunc%contact_list)) then
      count_only = .true.
      do i = 1, 2
        ndata = 0
        do j = enefunc%istart_contact, enefunc%iend_contact
          num = 0
          alist(1:num_unit) = enefunc%contact_list(1:num_unit,j)

          do k = 1, num_unit
            if (soltemp%is_solute(alist(k)) ) num = num + 1
          end do

          if (num > 0) then
            ndata = ndata + 1
            if (.not. count_only) then
              soltemp%contact_list(ndata) = j
              soltemp%contact_weight(ndata) = real(num,wp) / real(num_unit,wp)
              if (allocated(enefunc%contact_lj6) ) &
                soltemp%contact_lj6_org(ndata) = enefunc%contact_lj6(j)
              if (allocated(enefunc%contact_lj10) ) &
                soltemp%contact_lj10_org(ndata) = enefunc%contact_lj10(j)
              if (allocated(enefunc%contact_lj12) ) &
                soltemp%contact_lj12_org(ndata) = enefunc%contact_lj12(j)
            end if
          end if
        end do

        if (count_only) then
          soltemp%num_contacts = ndata
          allocate(soltemp%contact_list(ndata), &
                   soltemp%contact_weight(ndata))
          if (allocated(enefunc%contact_lj6) ) &
            allocate(soltemp%contact_lj6_org(ndata))
          if (allocated(enefunc%contact_lj10) ) &
            allocate(soltemp%contact_lj10_org(ndata))
          if (allocated(enefunc%contact_lj12) ) &
            allocate(soltemp%contact_lj12_org(ndata))
        end if
        count_only = .false.

      end do
    end if

    return

  end subroutine setup_remd_solute_tempering_native_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_remd
  !> @brief        control replica exchange
  !! @authors      TM
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_remd(output, molecule, enefunc, dynvars, dynamics,   &
                      pairlist, boundary, constraints, ensemble, remd)

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
    type(s_remd),             intent(inout) :: remd

    ! local variables
    integer                   :: i, replicaid, parmsetid


    ! Open output files
    !
    call open_output(output)

    ! output restart data
    !
    dynvars%step = 0
    call output_remd(0, output, molecule, dynamics, dynvars, boundary, remd)

    do i = 1, remd%ncycles

      dynamics%istart_step  = (i-1)*remd%exchange_period + 1
      dynamics%iend_step    =  i   *remd%exchange_period
      dynamics%initial_time = dynvars%time

      ! set conditions
      !
      replicaid = my_country_no + 1
      parmsetid = remd%repid2parmsetid(replicaid)
      call assign_condition(parmsetid, remd, ensemble, molecule, enefunc)

      ! MD main loop
      !
      if (dynamics%integrator == IntegratorLEAP) then
        call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      else if (dynamics%integrator == IntegratorVVER_CG) then
        call vverlet_dynamics_cg(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)
      end if

      ! perform remd
      !
      if (.not. remd%equilibration_only) then
        call perform_replica_exchange(i, molecule, enefunc, dynvars,         &
                                      ensemble, boundary, output, remd,      &
                                      pairlist)
      end if

      ! output restart data
      !
      call output_remd(i, output, molecule, dynamics, dynvars, boundary, remd)

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
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine perform_replica_exchange(icycle, molecule, enefunc, dynvars, &
                                      ensemble, boundary, output, remd,   &
                                      pairlist)

    ! formal arguments
    integer,                  intent(in)    :: icycle
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist

    ! local variables
    integer                   :: exptrn, dimno, itmp, parmsetid
    integer                   :: i, j, k, neighid, id_old
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
!        remd%random_table(i,j) = random(iseed)
        remd%random_table(i,j) = random_get_legacy(iseed)
        remd%random_table(j,i) = remd%random_table(i,j)
      end do
    end do

    
    ! perform replica exchange
    !
    select case(remd%types(dimno))

    case(RemdTemperature)

      call temperature_remd(dimno, exptrn, molecule, dynvars, remd)

    case(RemdPressure)

      call pressure_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdGamma)

      call surface_tension_remd(dimno, exptrn, ensemble, boundary, remd)

    case(RemdRestraint)

      call restraint_remd(dimno, exptrn, ensemble, boundary, dynvars, &
                          enefunc, remd)
 
    case(RemdSoluteTempering)

      call solute_tempering_remd(dimno, exptrn, molecule, ensemble, dynvars, &
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
          do j = 1, molecule%num_atoms
            dynvars%velocity(1:3,j) = dynvars%velocity(1:3,j) * factor
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

      write(MsgOut,'(A,I10,3X,A,I4,3X,A,I4)')              &
        'REMD> Step: ',dynvars%step,'Dimension: ',dimno,   &
        'ExchangePattern: ',exptrn

      write(MsgOut,'(A,6X,A,13X,A,6X,A,7X,A)')             &
        '  Replica', 'ExchangeTrial', 'AcceptanceRatio',   &
        'Before', 'After'
      do i = 1, total_nreplicas
        id_old = repid2parmsetid_ref(i)
        call get_neighbouring_parmsetid(dimno, exptrn, remd, id_old, neighid)

        if (remd%types(dimno) /= RemdRestraint) then

          dpara_ref = dparameters(dimno,parmidsets(repid2parmsetid_ref(i),     &
                      dimno))
          dpara     = dparameters(dimno,parmidsets(repid2parmsetid(i),         &
                      dimno))

          if (abs(dpara_ref - dpara) < EPS) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2F12.3)')              &
               i, id_old, ' > ', neighid, aorr,                                &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ',      &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),             &
               dpara_ref, dpara
        else

          ipara_ref = iparameters(dimno,parmidsets(repid2parmsetid_ref(i),     &
                      dimno))
          ipara     = iparameters(dimno,parmidsets(repid2parmsetid(i),         &
                      dimno))

          if (ipara_ref == ipara) then
            aorr = 'R'
          else
            aorr = 'A'
          end if
          if (neighid == 0) aorr = 'N'

          write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2I12)')                &
               i, id_old, ' > ', neighid, aorr,                                &
               num_exchanges(repid2parmsetid_ref(i),dimno,exptrn), ' / ',      &
               num_criteria (repid2parmsetid_ref(i),dimno,exptrn),             &
               ipara_ref, ipara
        end if
      end do

      write(MsgOut,*) ''

      write(ctrep,'(I3)') total_nreplicas
      frmt1 = '(A,' // trim(adjustl(ctrep)) // 'F10.3)'
      frmt2 = '(A,' // trim(adjustl(ctrep)) // 'I10  )'

      if (remd%types(dimno) /= RemdRestraint) then
        write(MsgOut,frmt1) '  Parameter    : ', &
          dparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),     &
          dimno))
      else
        write(MsgOut,frmt2) '  Parameter    : ', &
          iparameters(dimno,parmidsets(repid2parmsetid(1:total_nreplicas),     &
          dimno))
      end if
      write(MsgOut,frmt2)   '  RepIDtoParmID: ',                               &
                            repid2parmsetid(1:total_nreplicas)
      write(MsgOut,frmt2)   '  ParmIDtoRepID: ',                               &
                            parmsetid2repid(1:total_nreplicas)

      write(MsgOut,*) ''

    end if

    return

  end subroutine perform_replica_exchange

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_autoadj_displacement(dimno, remd)

    ! formal arguments
    integer,                    intent(in)    :: dimno
    type(s_remd),               intent(inout) :: remd

    ! local variables
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
  !! @param[inout] molecule : molecule information
  !! @param[inout] remd     : REMD information
  !! @note         Y.Sugita & Y.Okamoto, Chem.Phys.Lett., 314, 141−151 (1999)
  !!               K.Hukushima & K.Nemoto, J.Phys.Soc.Jpn., 65, 1604-1608 (1996)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074001 (2010)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074003 (2010)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine temperature_remd(dimno, exptrn, molecule, dynvars, remd)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_molecule), target, intent(inout) :: molecule
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_remd),     target, intent(inout) :: remd

    ! local variables
    integer                   :: i, j, k
    integer                   :: repid_i, repid_j, neighid
    integer                   :: parmsetid, parmsetid_i, parmsetid_j
    real(wp)                  :: temp_i, temp_j, beta_i, beta_j, factor
    real(wp)                  :: energy_i, energy_j, delta
    real(wp)                  :: before_gather
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
    before_gather = dynvars%energy%total
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
            dynvars%velocity(1:3,1:molecule%num_atoms)  &
              = dynvars%velocity(1:3,1:molecule%num_atoms) * factor

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
    integer,         pointer :: total_nreplicas, iseed
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
    iseed           => remd%iseed
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
    integer,         pointer :: total_nreplicas, iseed
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
    iseed           => remd%iseed
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
  !! @param[in]    ensemble : ensemble information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] dynvars  : dynamics variables information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] remd     : REMD information
  !! @note         Y.Sugita et al., J.Chem.Phys., 113, 6042-6051 (2000)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine restraint_remd(dimno, exptrn, ensemble, boundary, dynvars, &
                            enefunc, remd)

    ! formal arguments
    integer,                 intent(in)    :: dimno
    integer,                 intent(in)    :: exptrn
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),target, intent(in)    :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_remd),    target, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k
    integer                  :: natoms
    integer                  :: repid_i, repid_j
    integer                  :: parmsetid, parmsetid_i, parmsetid_j
    integer                  :: neighid
    integer                  :: funcid, umbr_0, umbr_1
    real(wp)                 :: energy_i_0, energy_j_0
    real(wp)                 :: energy_i_1, energy_j_1
    real(wp)                 :: eneval, beta, delta, temperature
    real(wp)                 :: before_gather
    real(wp)                 :: cv
    logical                  :: do_trial, do_exchange, calc_force, umbr_exist
    integer,         pointer :: total_nreplicas, iseed
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)
    integer,         pointer :: parmidsets(:,:)
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,         pointer :: parameters(:,:)
    real(wp),        pointer :: after_gather(:), random_table(:,:)
    real(wp),        pointer :: force(:,:), virial(:,:), virial_ext(:,:)
    real(wp),        pointer :: force_omp(:,:,:)
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
    iseed           => remd%iseed
    refcoord0       => remd%restraint_refcoord0
    refcoord1       => remd%restraint_refcoord1
    force           => dynvars%force
    force_omp       => dynvars%force_omp
    virial          => dynvars%virial
    virial_ext      => dynvars%virial_extern


    ! set target restraint function
    !
    repid_i     = my_country_no + 1
    parmsetid_i = repid2parmsetid(repid_i)
    umbr_0      = parameters(dimno,parmidsets(parmsetid_i,dimno))

    ! calculate restraint potential energy at the current condition
    !
    before_gather = 0.0_wp
    do k = 1, remd%umbrid2numfuncs(umbr_0)
      funcid = remd%umbrid2funclist(umbr_0,k)
      enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbr_0,k)
      enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbr_0,k)
      calc_force = .false.

      call compute_energy_restraint_each(funcid, calc_force, enefunc,    &
                                         boundary, dynvars%coord, force, &
                                         virial, virial_ext, eneval, cv)

      before_gather = before_gather + eneval
    end do

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
      before_gather = 0.0_wp
      do k = 1, remd%umbrid2numfuncs(umbr_1)
        funcid = remd%umbrid2funclist(umbr_1,k)
        enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbr_1,k)
        enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbr_1,k)
        calc_force = .false.

        call compute_energy_restraint_each(funcid, calc_force, enefunc,    &
                                           boundary, dynvars%coord, force, &
                                           virial, virial_ext, eneval, cv)
        before_gather = before_gather + eneval
      end do
    else
      before_gather = 0.0_wp
    end if

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
  !! @param[inout] molecule   : molecule information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] dynvars    : dynamics variables information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] remd       : REMD information
  !! @param[inout] pairlist   : pairlist information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine solute_tempering_remd(dimno, exptrn, molecule, ensemble, &
                                   dynvars, enefunc, remd, pairlist,  &
                                   boundary, output )

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_molecule), target, intent(inout) :: molecule
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_dynvars),          intent(inout) :: dynvars
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_remd),     target, intent(inout) :: remd
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output

    ! pointers
    real(wp),        pointer :: after_gather(:)
    real(wp),        pointer :: random_table(:,:)
    integer,         pointer :: total_nreplicas
    integer,         pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,         pointer :: parmidsets(:,:)
    real(wp),        pointer :: parameters(:,:)
    integer,         pointer :: num_criteria(:,:,:), num_exchanges(:,:,:)

    ! local variables
    real(wp) :: before_gather
    real(wp) :: energy_i, energy_j, delta, beta
    integer  :: repid_i, parmsetid_i
    integer  :: repid_j, parmsetid_j
    integer  :: neighid
    logical  :: does_exist, do_trial, do_exchange
    integer  :: i, ierror

    ! REST
    integer                     :: nfuncs
    type(s_energy),        save :: energy
    real(wp), allocatable, save :: force(:,:)
    real(wp), allocatable, save :: force_omp(:,:,:)
    real(wp)                    :: virial(3,3)
    real(wp)                    :: virial_ext(3,3)


    ! REST allocation
    nfuncs = enefunc%num_restraintfuncs
    if (.not. allocated(force)) then
      call alloc_energy(energy, EneRestraints, nfuncs)
      allocate(force(3,molecule%num_atoms), &
               force_omp(3,molecule%num_atoms,nthread))
    end if

    total_nreplicas => remd%total_nreplicas
    after_gather    => remd%after_gather
    random_table    => remd%random_table
    parmidsets      => remd%parmidsets
    parameters      => remd%dparameters
    num_criteria    => remd%num_criteria
    num_exchanges   => remd%num_exchanges
    repid2parmsetid => remd%repid2parmsetid
    parmsetid2repid => remd%parmsetid2repid

    ! calculate energy at current coordinate
    !
    call compute_energy(molecule, enefunc, pairlist, boundary,    &
                        .true., enefunc%nonb_limiter,             &
                        dynvars%coord, dynvars%trans,             &
                        dynvars%coord_pbc,                        &
                        energy, dynvars%temporary,                &
                        force, force_omp, virial, virial_ext)

    before_gather = energy%total
#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(before_gather, 1, mpi_wp_real, &
                       after_gather,  1, mpi_wp_real, &
                       mpi_comm_airplane, ierror )
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
      call assign_condition_rest(parmsetid_j, remd, molecule, ensemble, &
                                 enefunc)

      ! calculate energy
      !
      call compute_energy(molecule, enefunc, pairlist, boundary,    &
                          .true., enefunc%nonb_limiter,             &
                          dynvars%coord, dynvars%trans,             &
                          dynvars%coord_pbc,                        &
                          energy, dynvars%temporary,                &
                          force, force_omp, virial, virial_ext)

    else
      energy%total = 0.0_wp
    end if

    before_gather = energy%total
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

      ! get neighbouring parameter set index and replica index
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

        beta     = 1.0_wp / (KBOLTZ * ensemble%temperature)
        energy_i = remd%potential_energy(repid_i)
        energy_j = remd%potential_energy(repid_j)
        delta    = - beta * (energy_i + energy_j)

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
            .and. my_country_no + 1 == i) then
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

          ! nothing to do else; just update replica index
          !
          repid2parmsetid(i) = parmsetid_j

        else
          ! have to put back enefunc to the original state
          !
          call assign_condition_rest(parmsetid_i, remd, molecule, ensemble, &
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
  !! @param[in]    parmsetid   :
  !! @param[i]     remd        : REMD information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition(parmsetid, remd, ensemble, molecule, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_remd),            intent(in)    :: remd
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k
    integer                  :: umbrid, funcid


    do i = 1, remd%dimension
      if      (remd%types(i) == RemdTemperature) then
        ensemble%temperature = remd%dparameters(i,remd%parmidsets(parmsetid,i))

        if (enefunc%eef1_use) then
          call setup_eef1_temperature(ensemble%temperature, enefunc)
        else if (enefunc%gbsa_use) then
          enefunc%gbsa%temperature = ensemble%temperature
        end if
        if (enefunc%cg_ele_calc) then
          enefunc%cg_ele_sol_T = ensemble%temperature
        end if

      else if (remd%types(i) == RemdPressure  ) then
        ensemble%pressure    = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdGamma     ) then
        ensemble%gamma       = remd%dparameters(i,remd%parmidsets(parmsetid,i))
      else if (remd%types(i) == RemdRestraint ) then
        umbrid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
        do k = 1, remd%umbrid2numfuncs(umbrid)
          funcid = remd%umbrid2funclist(umbrid,k)
          enefunc%restraint_const(1:2,funcid) = remd%rest_constants(umbrid,k)
          enefunc%restraint_ref  (1:2,funcid) = remd%rest_reference(umbrid,k)
        end do
      end if
    end do

    call assign_condition_rest(parmsetid, remd, molecule, ensemble, enefunc)

    return

  end subroutine assign_condition

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_condition_rest
  !> @brief        reassign force constant involving REST
  !! @authors      MK
  !! @param[in]    parmsetid   :
  !! @param[in]    remd        : remd information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_rest(parmsetid, remd, molecule, &
                                   ensemble, enefunc)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_remd),            intent(in)    :: remd
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_molecule),        intent(inout) :: molecule
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
        ! scale factorsqrt(ensemble%temperature / solute_temperature)
        call assign_condition_solute_tempering(remd%solute_tempering(i), &
                                               rest_param, molecule, &
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
  !! @param[in]    soltemp     : REST information
  !! @param[in]    rest_param  : target temperature
  !! @param[inout] molecule    : molecule information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_condition_solute_tempering(soltemp, rest_param, &
                                               molecule, ensemble, enefunc)

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(in)    :: rest_param
    type(s_molecule),        intent(inout) :: molecule
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc

    ! local vairables
    integer  :: i, tgt, num_atom_cls, n, org
    real(wp) :: coeff_full, coeff_half, wgt, u_self, el_fact


    coeff_full = ensemble%temperature / rest_param
    coeff_half = sqrt(coeff_full)

    ! charge
    if (soltemp%sw_charge) then
      do i = 1, soltemp%num_solute
        molecule%charge(soltemp%solute_list(i)) = &
          coeff_half * soltemp%charge_org(i)
      end do

      ! update PME self energy
      u_self = 0.0_wp
      el_fact = ELECOEF / enefunc%dielec_const
      do i = my_country_rank + 1, molecule%num_atoms, nproc_country
        u_self = u_self + molecule%charge(i) ** 2
      end do
      enefunc%pme%u_self = - u_self * el_fact * enefunc%pme%alpha/sqrt(PI)
    end if

    ! LJ
    if (soltemp%sw_lj) then
      n = enefunc%num_atom_cls
      if (allocated(enefunc%nb14_lj6)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nb14_lj6 )
      end if
      if (allocated(enefunc%nb14_lj10)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nb14_lj10 )
      end if
      if (allocated(enefunc%nb14_lj12)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nb14_lj12 )
      end if
      if (allocated(enefunc%nonb_lj6)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nonb_lj6 )
      end if
      if (allocated(enefunc%nonb_lj10)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nonb_lj10 )
      end if
      if (allocated(enefunc%nonb_lj12)) then
        call assign_condition_solute_tempering_lj( &
            n, coeff_half, coeff_full, soltemp, enefunc%nonb_lj12 )
      end if
    end if

    ! bonds
    if (soltemp%sw_bonds) then
      do i = 1, soltemp%num_bonds
        tgt = soltemp%bond_list(i)
        wgt = coeff_full ** soltemp%bond_weight(i)
        enefunc%bond_force_const(tgt) = soltemp%bond_force_const_org(i) * wgt
      end do
    end if

    ! angles
    if (soltemp%sw_angles) then
      do i = 1, soltemp%num_angles
        tgt = soltemp%angl_list(i)
        wgt = coeff_full ** soltemp%angl_weight(i)
        enefunc%angl_force_const(tgt) = soltemp%angl_force_const_org(i) * wgt
      end do
    end if

    ! urey
    if (soltemp%sw_ureys) then
      do i = 1, soltemp%num_ureys
        tgt = soltemp%urey_list(i)
        wgt = coeff_full ** soltemp%urey_weight(i)
        enefunc%urey_force_const(tgt) = soltemp%urey_force_const_org(i) * wgt
      end do
    end if

    ! dihedral
    if (soltemp%sw_dihedrals) then
      do i = 1, soltemp%num_dihedrals
        tgt = soltemp%dihe_list(i)
        wgt = coeff_full ** soltemp%dihe_weight(i)
        enefunc%dihe_force_const(tgt) = soltemp%dihe_force_const_org(i) * wgt
      end do
    end if

    ! improper
    if (soltemp%sw_impropers) then
      do i = 1, soltemp%num_impropers
        tgt = soltemp%impr_list(i)
        wgt = coeff_full ** soltemp%impr_weight(i)
        enefunc%impr_force_const(tgt) = soltemp%impr_force_const_org(i) * wgt
      end do
    end if

    ! cmaps
    if (soltemp%sw_cmaps) then
      if (soltemp%num_cmap_type > 0) then
        do i = 1, soltemp%num_cmap_type
          tgt = i + soltemp%istart_cmap_type - 1
          org = soltemp%cmap_type_org(i)
          wgt = coeff_full ** soltemp%cmap_weight(i)
          enefunc%cmap_coef(1:4,1:4,1:24,1:24,tgt) = &
            wgt * enefunc%cmap_coef(1:4,1:4,1:24,1:24,org)
        end do
      end if
    end if

    ! contacts
    if (soltemp%sw_contacts) then
      if (allocated(enefunc%contact_list)) then
        do i = 1, soltemp%num_contacts
          tgt = soltemp%contact_list(i)
          wgt = coeff_full ** soltemp%contact_weight(i)
          if (allocated(soltemp%contact_lj6_org) ) &
            enefunc%contact_lj6(tgt) = soltemp%contact_lj6_org(i) * wgt
          if (allocated(soltemp%contact_lj10_org) ) &
            enefunc%contact_lj10(tgt) = soltemp%contact_lj10_org(i) * wgt
          if (allocated(soltemp%contact_lj12_org) ) &
            enefunc%contact_lj12(tgt) = soltemp%contact_lj12_org(i) * wgt
        end do
      end if
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

end module at_remd_mod
