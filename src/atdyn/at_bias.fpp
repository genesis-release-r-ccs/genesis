!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_bias_mod
!> @brief   Bias exhcnage ABMD simulation
!! @authors Shingo Ito (SI)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_bias_mod

  use at_energy_mod
  use at_energy_str_mod
  use at_energy_restraints_mod
  use at_energy_abmd_mod
  use at_energy_reac_coord_mod
  use at_md_leapfrog_mod
  use at_md_vverlet_mod
  use at_md_vverlet_cg_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_bias_str_mod
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
  type, public :: s_bias_info
    integer                             :: nreplicas          = 1
    integer                             :: nwalkers           = 1
    integer                             :: exchange_period    = 100
    integer                             :: iseed              = 3141592
    integer                             :: type               = BiasExchange
    logical                             :: force_exchange     = .false.
    logical,                allocatable :: cyclic_params(:)
    character(MaxLine),     allocatable :: reac_function(:)
    ! gREST
    !logical                             :: analysis_grest     = .false.
    !character(MaxLine)                  :: param_type         = 'ALL'
    !character(MaxLine)                  :: select_index       = ''
    !character(MaxLineLong)              :: parameters         = ''
    ! auto-adjust
    real(wp),               allocatable :: tgt_exc_prob(:)
    real(wp),               allocatable :: mgn_exc_prob(:)
    integer,                allocatable :: trial_freq(:)
    integer,                allocatable :: eq_cycle(:)
    real(wp),               allocatable :: param_grid(:)
    real(wp),               allocatable :: max_param_shift(:)
  end type s_bias_info

  ! variables
  integer, public, parameter          :: BiasMaxNdimensions = 999

  ! subroutines
  public  :: show_ctrl_bias
  public  :: read_ctrl_bias
  public  :: run_bias
  public  :: setup_bias

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_bias
  !> @brief        show BIAS section usage
  !! @authors      SI
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "bias"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_bias(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('bias')

        write(MsgOut,'(A)') '[BIAS]'
        write(MsgOut,'(A)') 'type              = BIAS'
        write(MsgOut,'(A)') 'exchange_period   = 100'
        write(MsgOut,'(A)') 'nreplica          = 2'
        write(MsgOut,'(A)') 'nwalkers          = 1'
        write(MsgOut,'(A)') 'reac_function1    = 1'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('bias')

        write(MsgOut,'(A)') '[BIAS]'
        write(MsgOut,'(A)') 'type              = BIAS'
        write(MsgOut,'(A)') 'exchange_period   = 100'
        write(MsgOut,'(A)') 'nreplica          = 1'
        write(MsgOut,'(A)') 'nwalkers1         = 1'
        write(MsgOut,'(A)') 'reac_function1    = 1'
        write(MsgOut,'(A)') ' '

      end select


    end if

    return

  end subroutine show_ctrl_bias
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_bias
  !> @brief        read BIAS section in the control file
  !! @authors      SI
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   bias_info : REPLICA section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_bias(handle, bias_info)

    ! parameters
    character(*),            parameter     :: Section = 'Bias'

    ! formal arguments
    integer,                 intent(in)     :: handle
    type(s_bias_info),       intent(inout)  :: bias_info

    ! local variables
    integer                  :: i, trep
    character(30)            :: cdim, numtmp, partmp, vartmp
    character(30)            :: cnstmp, reftmp, fnctmp, cyctmp, seltmp
    character(30)            :: typtmp, adjtmp, rectmp


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'nreplicas',       &
                               bias_info%nreplicas)
    call read_ctrlfile_integer(handle, Section, 'nwalkers',        &
                               bias_info%nwalkers)
    call read_ctrlfile_type   (handle, Section, 'type',            &
                               bias_info%type, BiasTypes)
    !call read_ctrlfile_string (handle, Section, 'parameters',      &
    !                           bias_info%parameters)
    !call read_ctrlfile_string (handle, Section, 'select_index',    &
    !                           bias_info%select_index)
    !call read_ctrlfile_string (handle, Section, 'param_type',      &
    !                           bias_info%param_type)                                                                                                        
    call read_ctrlfile_integer(handle, Section, 'exchange_period', &
                               bias_info%exchange_period)
    call read_ctrlfile_integer(handle, Section, 'iseed',           &
                               bias_info%iseed)
    call read_ctrlfile_logical(handle, Section, 'force_exchange',  &
                               bias_info%force_exchange)      
                               
    ! error check
    !
    if (bias_info%nreplicas > BiasMaxNdimensions .or. bias_info%nreplicas == 0) &
      call error_msg('Read_Ctrl_Bias> Invalid dimension specified in [BIAS] section.')

    ! read allocatable variables
    !
    allocate(bias_info%cyclic_params(bias_info%nreplicas),   &
             bias_info%reac_function(bias_info%nreplicas))

    !bias_info%parameters     (1:bias_info%nreplicas) = ''
    bias_info%cyclic_params  (1:bias_info%nreplicas) = .false.
    bias_info%reac_function  (1:bias_info%nreplicas) = ''
    !bias_info%select_index   (1:bias_info%nreplicas) = ''
    !bias_info%param_type     (1:bias_info%nreplicas) = 'ALL'

    do i = 1, bias_info%nreplicas
      write(cdim,'(i0)') i

      !partmp = 'parameters'    // cdim
      !call read_ctrlfile_string (handle, Section, partmp,  &
      !                           bias_info%parameters(i))

      cyctmp = 'cyclic_params' // cdim
      call read_ctrlfile_logical(handle, Section, cyctmp,  &
                                 bias_info%cyclic_params(i))

      rectmp = 'reac_function' // cdim
      call read_ctrlfile_string (handle, Section, rectmp,  &
                                 bias_info%reac_function(i))

      !seltmp = 'select_index' // cdim
      !call read_ctrlfile_string (handle, Section, seltmp,  &
      !                           bias_info%select_index(i))

      !typtmp = 'param_type'   // cdim
      !call read_ctrlfile_string (handle, Section, typtmp,  &
      !                           bias_info%param_type(i))

    end do

    call end_ctrlfile_section(handle)

    ! error check
    !
    if (bias_info%nwalkers == 0) &
      call error_msg('Read_Ctrl_Bias> Error in nwalkers in [BIAS]')


    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Bias> Bias-Exchange information'

      if (bias_info%type == BiasExchange) then
        write(MsgOut,'(A)')   '  variable        = BIAS'
      end if

      do i = 1, bias_info%nreplicas

        write(MsgOut,'(A20,I10,A20,I10)') &
          '  nreplicas       = ', i,      &
          '  nwalkers        = ', bias_info%nwalkers

        if (bias_info%cyclic_params(i)) then
          write(MsgOut,'(A)') '  cyclic_params   = yes'
        else
          write(MsgOut,'(A)') '  cyclic_params   = no'
        end if

        write(MsgOut,'(A)') ''
      end do

    end if

    return

  end subroutine read_ctrl_bias

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_bias
  !> @brief        setup BE-ABMD
  !! @authors      SI
  !! @param[in]    bias_info  : REMD section control parameters information
  !! @param[in]    rst        : Restart data
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] restraints : restraints information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] bias       : BIAS information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_bias(bias_info, rst, boundary, dynamics, molecule, &
                        restraints, ensemble, enefunc, bias)

    ! formal arguments
    type(s_bias_info),       intent(in)    :: bias_info
    type(s_rst),             intent(in)    :: rst
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_molecule),        intent(inout) :: molecule
    type(s_restraints),      intent(inout) :: restraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_bias),    target, intent(inout) :: bias
    type(s_enefunc_reac_coord)             :: reac_coord

    ! local variables
    integer                    :: i, j, k, m
    integer                    :: max_nwalkers, cycles, nrest
    integer                    :: itmp, jtmp, ktmp
    integer                    :: comp, dsta, dend, found, dtotal
    integer                    :: dockid, replicaid, walkerid
    integer                    :: parmsetid1
    integer                    :: ncycle, icycle, nlen, ixx
    integer                    :: funcid
    integer                    :: ndata, max_ndata
    integer,       allocatable :: idata(:)
    integer                    :: iref, icomp
    character(MaxLine)         :: param
    character(10)              :: tmp, frmt, cd
    character(18)              :: frmt2
    character(3)               :: ctrep
    logical                    :: keep, double
    real(wp), allocatable      :: ddata(:)
    integer,       allocatable :: ndigit(:), parmid(:)
    integer,       allocatable :: sendbuf(:), recvbuf(:,:)

    ! local variables for gREST
    !integer,         parameter :: max_param_type = 20
    !character(MaxLine)         :: param_type_str(max_param_type)
    
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Bias> Bias-Exchange information'
      write(MsgOut,'(A)') ''
    end if


    ! initialize
    !
    enefunc%bias_use     = .true.
    bias%nreplicas       = bias_info%nreplicas
    bias%nwalkers        = bias_info%nwalkers
    bias%type            = bias_info%type
    bias%force_exchange  = bias_info%force_exchange
    bias%ndocks          = bias%nreplicas * bias%nwalkers
    dockid               = my_country_no + 1

    call alloc_bias(bias, ABMDBias, bias%ndocks)

    do i = 1, bias%nreplicas
      bias%cyclic_params(i) = bias_info%cyclic_params(i)
    end do

    ! check cyclic_params
    !
    do i = 1, bias%nreplicas
      if (bias%cyclic_params(i)) then
        if (mod(bias%nwalkers,2) /= 0) &
          call error_msg('Setup_Bias> nwalkers must be even,'// &
                         ' when cyclic_params is specified')
      end if
    end do

    ! setup bias_iseed, num_criteria, and num_exchange
    !
    if (rst%rstfile_type == RstfileTypeRemd) then
      bias%iseed = rst%iseed_remd
      do i = 1, bias%ndocks
        bias%dockid2repsetid(i)       = rst%repid2parmsetid(i)
        bias%dockid2repsetid_ref(i)   = bias%dockid2repsetid(i)
        bias%num_criteria_rep (i,1:2) = rst%num_criteria (i,1,1:2)
        bias%num_exchanges_rep(i,1:2) = rst%num_exchanges(i,1,1:2)
      end do
    else
      bias%iseed = bias_info%iseed
      do i = 1, bias%ndocks
        bias%num_criteria_rep (i,1:2)        = 0
        bias%num_exchanges_rep(i,1:2)        = 0
        bias%dockid2repsetid(i)              = i
        bias%dockid2repsetid_ref(i)          = i  
      end do
    end if
  
    ! setup permutation function
    !
    do i = 1, bias%ndocks
      parmsetid1 = bias%dockid2repsetid(i)
      bias%repsetid2dockid(parmsetid1) = i
    end do

    !k = 1
    !do i = 1, bias%nwalkers
    !  do j = 1, bias%nreplicas
    !    bias%repsetid2paramid(k,1) = j ! CV
    !    bias%repsetid2paramid(k,2) = i ! Temperature
    !    k = k + 1
    !  end do
    !end do

    k = 1
    do i = 1, bias%nwalkers
      do j = 1, bias%nreplicas
        bias%repsetid2paramid(k,1) = j ! CV
        bias%repsetid2paramid(k,2) = i ! Temperature
        k = k + 1
      end do
    end do
        
    ! setup replica exchange period
    !
    if (bias_info%exchange_period > 0) then
      if (dynamics%eneout_period > bias_info%exchange_period) then
        call error_msg('Setup_Bias> eneout_period must be less than or'//    &
                       ' equal to exchange period')
      else
        if (mod(bias_info%exchange_period,dynamics%eneout_period) /= 0) then
          call error_msg('Setup_Bias> exchange_period must be a multiple '//  &
                         'of eneout_period)')
        end if
      end if

      bias%equilibration_only = .false.
      bias%exchange_period = bias_info%exchange_period
      bias%ncycles = dynamics%nsteps/bias%exchange_period

    else if (bias_info%exchange_period == 0) then
      bias%equilibration_only = .true.
      bias%exchange_period = dynamics%nsteps
      bias%ncycles = 1
    else
      call error_msg('Setup_Bias> error in exchange_period')
    end if

    ! count numbers for allocation of reac_function
    !
    max_ndata = 0
    do i = 1, bias%nreplicas
      param = bias_info%reac_function(i)
      ndata = split_num(param)
      if (ndata > max_ndata) max_ndata = ndata
    end do
    
    if (max_ndata > 4) &
     call error_msg('Setup_Bias> The number of reaction functions cannot exceed'//&
                    ' 4')

    ! set up enefunc%reac_coord%nreac_func
    !
    allocate(idata(max_ndata))

    do i = 1, bias%nreplicas
      param = bias_info%reac_function(i)
      ndata = split_num(param)
      j = bias%dockid2repsetid(dockid)
      k = bias%repsetid2paramid(j,1)
      if (i == k) then
        enefunc%reac_coord%nreac_func = ndata
      end if

      call split(ndata, max_ndata, param, idata)

      do j = 1, ndata
        bias%ifunc(i,j) = idata(j)
      end do
    end do

    deallocate(idata)   

    bias%dimension = 1
    
    j = bias%dockid2repsetid(dockid)

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') '  DockID = ReplicaID, WalkerID'

      frmt2 = '(I11, A, I11, I11)'

    end if

    return

  end subroutine setup_bias

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_bias
  !> @brief        control bias exchange
  !! @authors      SI
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : pairlist information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] bias        : BIAS information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_bias(output, molecule, enefunc, dynvars, dynamics,   &
                      pairlist, boundary, constraints, ensemble, bias)

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
    type(s_bias),             intent(inout) :: bias

    ! local variables
    integer                   :: i
    integer                   :: replicaid, walkerid, dockid
    integer                   :: parmsetid1, parmsetid2

    ! Open output files
    !
    call open_output(output)

    ! output restart data
    !
    dynvars%step = 0
    call output_bias(0, output, molecule, dynamics, dynvars, boundary, enefunc, bias)

    do i = 1, bias%ncycles

      dynamics%istart_step  = (i-1)*bias%exchange_period + 1
      dynamics%iend_step    =  i   *bias%exchange_period
      dynamics%initial_time = dynvars%time

      ! set conditions
      !
      dockid     = my_country_no   + 1
      parmsetid1 = enefunc%abmd%dockid2repsetid(dockid)

      ! MD main loop
      !
      if (dynamics%integrator == IntegratorLEAP) then
        call leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble, &
                               bias=bias)
      else if (dynamics%integrator == IntegratorVVER) then
        call vverlet_dynamics (output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble, &
                               bias=bias)
      else if (dynamics%integrator == IntegratorVVER_CG) then
        call vverlet_dynamics_cg(output, molecule, enefunc, dynvars, dynamics, &
                                 pairlist, boundary, constraints, ensemble)
      end if

      if (enefunc%abmd_use .and. enefunc%abmd%stop_when_out_of_range       &
                           .and. enefunc%reac_coord%out_of_range) then
        exit
      endif
      
      ! perform bias-exchange
      !
      if (.not. bias%equilibration_only) then
        call perform_replica_exchange(i, molecule, enefunc, dynvars,         &
                                      ensemble, boundary, output, bias,      &
                                      pairlist)
        !enefunc%abmd%dockid2repsetid = bias%dockid2repsetid
      end if

      ! output restart data
      !
      call output_bias(i, output, molecule, dynamics, dynvars, boundary, enefunc, bias)

      ! output spline coefficients
      !
      call output_spline(dynvars, enefunc, output)

    end do

    ! close output files
    !
    call close_output(output)


    return

  end subroutine run_bias

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    perform_replica_exchange
  !> @brief        perform replica exchange
  !! @authors      SI
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamics variables information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] bias        : BIAS information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine perform_replica_exchange(icycle, molecule, enefunc, dynvars, &
                                      ensemble, boundary, output, bias,   &
                                      pairlist)
  
    ! formal arguments
    integer,                  intent(in)    :: icycle
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_boundary),         intent(inout) :: boundary
    type(s_output),           intent(inout) :: output
    type(s_bias),     target, intent(inout) :: bias
    type(s_pairlist),         intent(inout) :: pairlist
  
    ! local variables
    integer                   :: exptrn, dimno, itmp
    integer                   :: parmsetid1, parmsetid2
    integer                   :: id_old1, id_old2
    integer                   :: i, j, k, neighid
    integer                   :: dpara, dpara_ref
    integer                   :: iseed
    integer,          pointer :: nreplicas, nwalkers, ndocks
    integer,          pointer :: dockid2repsetid(:)
    integer,          pointer :: dockid2repsetid_ref(:)
    integer,          pointer :: num_criteria(:,:), num_exchanges(:,:)
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    integer,          pointer :: repid2parmsetid_ref(:)
    integer,          pointer :: iparameters(:,:)
    real(wp),         pointer :: dparameters(:,:)
    character                 :: ctrep*3, frmt1*12, frmt2*12, aorr*1

    integer, save, allocatable :: counter_skip(:)
    integer, save, allocatable :: counter_tuning(:)
  
    real(wp) :: temp_i, temp_j, factor
  
    nreplicas               => bias%nreplicas
    nwalkers                => bias%nwalkers
    ndocks                  => bias%ndocks
    iseed                   =  3939 ! Temporary
    dockid2repsetid         => enefunc%abmd%dockid2repsetid
    dockid2repsetid_ref     => bias%dockid2repsetid_ref
    num_criteria            => bias%num_criteria_rep
    num_exchanges           => bias%num_exchanges_rep
  
    ! decide dimension to be exchanged and exchange pattern
    !   pattern 1: 1<>2  3<>4  5<>6  7<>8 ...
    !   pattern 2: 1  2<>3  4<>5  6<>7  8 ...
    !
    itmp  = mod(icycle,2)
    dimno = 1

    if (itmp+1 <= 1) then
      exptrn = 1
    else
      exptrn = 2
    end if
  
    ! save previous information
    !
    dockid2repsetid_ref(1:ndocks) = dockid2repsetid(1:ndocks)
  
    ! generate random number
    !
    do i = 1, ndocks-1
      do j = i+1, ndocks
        bias%random_table(i,j) = random_get_legacy(iseed)
        bias%random_table(j,i) = bias%random_table(i,j)
      end do
    end do  
  
    ! perform replica exchange
    !
    select case(bias%type)
  
    case(BiasExchange)
  
      call bias_exchange(dimno, exptrn, molecule, dynvars, ensemble, bias, enefunc)
  
    end select
  
    ! update permutation function
    !
    do i = 1, bias%ndocks ! Dock id
      parmsetid1 = dockid2repsetid(i) ! Replica id
      bias%repsetid2dockid(parmsetid1) = i
    end do
  
    if (main_rank) then

      write(MsgOut,'(A,I10,3X,A,I4,3X,A,I4)')                           &
        'BIAS> Step: ',dynvars%step,'Dimension: ',dimno,                &
        'ExchangePattern: ',exptrn

      write(MsgOut,'(A,6X,A,13X,A,6X,A,7X,A)')     &
        '  Replica', 'ExchangeTrial', 'AcceptanceRatio', 'Before', 'After'
      do i = 1, bias%ndocks
        id_old1 = dockid2repsetid_ref(i)
        neighid = 0        
        call get_neighbouring_parmsetid(exptrn, bias, i, neighid)


        dpara_ref = id_old1
        dpara     = dockid2repsetid(i)

        if (abs(dpara_ref - dpara) < EPS) then
          aorr = 'R'
        else
          aorr = 'A'
        end if
        if (neighid == 0) aorr = 'N'

        write(MsgOut,'(I9,6X,I5,A,I5,3X,A1,3X,I9,A,I9,2I12)')           &
             i, i, ' > ', neighid, aorr,                          &
             num_exchanges(i,exptrn), ' / ',       &
             num_criteria (i,exptrn),              &
             dpara_ref, dpara
 
      end do

      write(MsgOut,*) ''

    end if
  
    return
  
  end subroutine perform_replica_exchange

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bias_exchange
  !> @brief        bias-exchange MD
  !! @authors      SI
  !! @param[in]    dimno    : dimension number
  !! @param[in]    exptrn   : exchange pattern
  !! @param[inout] domain   : domain information
  !! @param[inout] dynvars  : dynamics variables information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] remd     : REMD information
  !! @note         Y.Sugita & Y.Okamoto, Chem.Phys.Lett., 314, 141−151 (1999)
  !!               K.Hukushima & K.Nemoto, J.Phys.Soc.Jpn., 65, 1604-1608 (1996)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074001 (2010)
  !!               Y.Mori & Y.Okamoto, J.Phys.Soc.Jpn., 79, 074003 (2010)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bias_exchange(dimno, exptrn, molecule, dynvars, ensemble, bias, enefunc)

    ! formal arguments
    integer,                  intent(in)    :: dimno
    integer,                  intent(in)    :: exptrn
    type(s_molecule), target, intent(inout) :: molecule
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_enefunc),  target, intent(inout) :: enefunc
    
    type(s_bias),     target, intent(inout) :: bias

    ! local variables
    integer                   :: i, j, k
    integer                   :: dockid_i,  dockid_j, neighid
    integer                   :: parmid_i, parmid_j
    integer                   :: iseed = 3939 ! Temporary for test
    real(wp)                  :: beta, delta, delta1, delta2

    integer                   :: dockid
    integer,          pointer :: nreplicas, ndocks
    integer,          pointer :: dockid2repsetid(:)
    integer,          pointer :: dockid2repsetid_ref(:)
    integer,          pointer :: repsetid2dockid(:)
    integer,          pointer :: num_criteria(:,:), num_exchanges(:,:)
    real(wip),        pointer :: temperature
    real(wp),         pointer :: after_gather(:), random_table(:,:)
    real(wip),        pointer :: coord(:,:)
    real(wp),         pointer :: bonds_coord(:,:)
    real(wp),     allocatable :: xi_gather(:,:)
  
    integer,  pointer         :: nreac_coord, nreac_coord_comb
    integer,  pointer         :: nxi(:), kind(:)
    integer,  pointer         :: bondslist(:,:,:)
    integer,  pointer         :: bondslist_to_atomlist(:,:)
    real(dp), pointer         :: uk(:)
    real(wp), pointer         :: ximin(:), ximax(:), dxi(:), invdxi(:)
    logical,  pointer         :: periodicity(:)

    real(wp)                  :: temp_i, temp_j, beta_i, beta_j, factor
    real(wp)                  :: energy_i, energy_j
    real(wp)                  :: before_gather
    logical                   :: do_trial, do_exchange, do_exchange_bias
    integer,          pointer :: parmidsets(:,:)
    integer,          pointer :: repid2parmsetid(:), parmsetid2repid(:)
    real(wp),         pointer :: dparameters(:)
    real(dp),     allocatable :: u(:), u_gather(:,:)  

    integer,  pointer         :: temp_nreac_coord, temp_nreac_coord_comb
    integer,  pointer         :: temp_nxi(:),   temp_kind(:)
    integer,  pointer         :: temp_bondslist(:,:,:)
    integer,  pointer         :: temp_bondslist_to_atomlist(:,:)
    real(dp), pointer         :: temp_uk(:)
    real(wp), pointer         :: temp_ximin(:), temp_ximax(:)
    real(wp), pointer         :: temp_dxi(:),   temp_invdxi(:)
    logical,  pointer         :: temp_periodicity(:)

    nreplicas                  => bias%nreplicas
    ndocks                     => bias%ndocks
    dockid2repsetid            => enefunc%abmd%dockid2repsetid
    dockid2repsetid_ref        => bias%dockid2repsetid_ref
    repsetid2dockid            => bias%repsetid2dockid
    dparameters                => bias%dparameters
    num_criteria               => bias%num_criteria_rep
    num_exchanges              => bias%num_exchanges_rep
    after_gather               => bias%after_gather
    random_table               => bias%random_table
    bonds_coord                => enefunc%reac_coord%bonds_coord
    coord                      => dynvars%coord

    nreac_coord                => enefunc%abmd%nreac_coord
    nreac_coord_comb           => enefunc%abmd%nreac_coord_comb
    nxi                        => enefunc%abmd%nxi
    kind                       => enefunc%abmd%kind
    ximin                      => enefunc%abmd%ximin
    ximax                      => enefunc%abmd%ximax
    dxi                        => enefunc%abmd%dxi
    invdxi                     => enefunc%abmd%invdxi
    periodicity                => enefunc%abmd%periodicity
    bondslist                  => enefunc%abmd%bondslist
    bondslist_to_atomlist      => enefunc%abmd%bondslist_to_atomlist
    uk                         => enefunc%abmd%uk
    
    temp_nreac_coord           => enefunc%abmd%temp_nreac_coord
    temp_nreac_coord_comb      => enefunc%abmd%temp_nreac_coord_comb
    temp_nxi                   => enefunc%abmd%temp_nxi
    temp_kind                  => enefunc%abmd%temp_kind
    temp_ximin                 => enefunc%abmd%temp_ximin
    temp_ximax                 => enefunc%abmd%temp_ximax
    temp_dxi                   => enefunc%abmd%temp_dxi
    temp_invdxi                => enefunc%abmd%temp_invdxi
    temp_periodicity           => enefunc%abmd%temp_periodicity
    temp_bondslist             => enefunc%abmd%temp_bondslist
    temp_bondslist_to_atomlist => enefunc%abmd%temp_bondslist_to_atomlist
    temp_uk                    => enefunc%abmd%temp_uk

    temperature                => ensemble%temperature

    do_exchange                = .false.
    do_exchange_bias           = .false.

    dockid                     =  my_country_no   + 1

    ! allgather xi when abmd is enabled
    !
    if(enefunc%abmd_use) then
      allocate(xi_gather(enefunc%reac_coord%nreac_coord, ndocks)) ! Xi in all walkers (= replicas, if nwalker = 1)
      allocate(u_gather(ndocks, ndocks))                          ! U in all walkers
      allocate(u(ndocks))                                         ! U_my in all CVs

      ! Compute all xi
      !
      call compute_reac_coord(enefunc, coord)

#ifdef HAVE_MPI_GENESIS
      ! allgather xi
      !
      call mpi_allgather(enefunc%reac_coord%xi,                  &
                         enefunc%reac_coord%nreac_coord,         &
                         mpi_wp_real, xi_gather,                 &
                         enefunc%reac_coord%nreac_coord,         &
                         mpi_wp_real, mpi_comm_airplane, ierror)

      !u(dockid) = enefunc%abmd%u ! Set U_my(R_my)
      do i = 1, ndocks
        ! Compute U_my(R_neigh)
        !
        call compute_bias_potential_abmd(enefunc, &
                                         sumu = u(i), xiex = xi_gather(:,i), &
                                         update_u = .false.)
      enddo

      ! allgater potential energy, U
      !
      call mpi_allgather(u, ndocks, mpi_real8, u_gather, &
                            ndocks, mpi_real8,           &
                            mpi_comm_airplane, ierror)     
#endif
    endif

    ! replica exchange
    !
    do i = 1, ndocks ! Docks

      dockid_i    = i ! Dock id
      parmid_i    = dockid2repsetid_ref(dockid_i) ! Replica = parameter id
      !parmid_i    = repsetid2dockid(repid_i) ! Replica = parameter id

      ! get neighboring parameter set index and replica index
      !
      call get_neighbouring_parmsetid(exptrn, bias, dockid_i, dockid_j)

      if (dockid_j /= 0) then
        parmid_j    = dockid2repsetid_ref(dockid_j) ! Dock id
        do_trial    = .true.
      else
        parmid_j    = 0
        do_trial    = .false.
      end if

      if (do_trial) then

        num_criteria(dockid_i,exptrn) &
          = num_criteria(dockid_i,exptrn) + 1

        ! compute delta
        ! 
        delta  = 0.0_wp 
        beta   = 1.0_wp / (KBOLTZ * temperature)

        delta      = beta * (u_gather(dockid_j, dockid_i) - &
                             u_gather(dockid_i, dockid_i) - &
                             u_gather(dockid_j, dockid_j) + &
                             u_gather(dockid_i, dockid_j))   

        ! Metropolis criterion
        !
        if (exptrn == 1 .and. bias%force_exchange) then !! Debugging
          do_exchange      = .true.
        else if (delta <= 0.0_wp) then
          do_exchange      = .true.
        else if (delta > 0.0_wp) then
          if (random_table(dockid_i,dockid_j) <= exp(-delta)) then
            do_exchange      = .true.
          else
            do_exchange = .false.
          end if
        end if

        if (do_exchange) then

          num_exchanges(dockid_i,exptrn) &
            = num_exchanges(dockid_i,exptrn) + 1

          dockid2repsetid(dockid_i)      = parmid_j
          bias%dockid2repsetid(dockid_i) = parmid_j

          ! Exchange bias potential parameters
          if (enefunc%abmd_use) then
            !call abmd_bias_potential_exchange(repid_i, repid_j, enefunc)
            call abmd_bias_potential_exchange(dockid_i, dockid_j, enefunc, &
                                              bias, do_exchange_bias)
          endif
        endif
 
      end if

    end do

    if (enefunc%abmd_use) then
      if (do_exchange_bias) then
        ! Update data
        nreac_coord           = temp_nreac_coord
        nreac_coord_comb      = temp_nreac_coord_comb                                                                   
        nxi                   = temp_nxi
        kind                  = temp_kind
        ximin                 = temp_ximin
        ximax                 = temp_ximax
        dxi                   = temp_dxi
        invdxi                = temp_invdxi
        periodicity           = temp_periodicity
        bondslist             = temp_bondslist
        bondslist_to_atomlist = temp_bondslist_to_atomlist

        if (allocated(enefunc%abmd%uk)) deallocate(enefunc%abmd%uk)
        allocate(enefunc%abmd%uk(temp_nreac_coord_comb))
        enefunc%abmd%uk = 0.0_wp
        enefunc%abmd%uk = enefunc%abmd%temp_uk

      end if
      
      do_exchange_bias = .false.

      ! Deallocate temp data        
      if(allocated(enefunc%abmd%temp_uk)) deallocate(enefunc%abmd%temp_uk)
      deallocate(xi_gather)
      deallocate(u_gather)
      deallocate(u)
      
    endif

    return

  end subroutine bias_exchange

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

  subroutine get_neighbouring_parmsetid(exptrn, bias, repid_i, repid_j)

    ! formal arguments
    integer,                 intent(in)    :: exptrn
    type(s_bias),    target, intent(in)    :: bias
    integer,                 intent(in)    :: repid_i
    integer,                 intent(out)   :: repid_j

    ! local variables
    integer                  :: j, k
    integer                  :: next, jump 
    logical                  :: cyclic
    
    
    cyclic  = .false.

    if (exptrn == 1) then
      if (mod(repid_i,2) == 0) then
        ! 1 <2  3 <4  5 <6  7 <8  ...
        next = -1
        if (cyclic) then
          !jump = remd%nreplicas(dimno) - 1
        else
          jump = next
        end if

        if(repid_i == 1) then 
          repid_j = 0
          return
        end if

      else
        ! 1> 2  3> 4  5> 6  7> 8  ...
        next = +1
        if (cyclic) then
          !jump = -remd%nreplicas(dimno) + 1
        else
          jump = next
        end if

        if(repid_i == bias%ndocks) then
          repid_j = 0
          return
        end if

      end if
    else if (exptrn == 2) then
      if (mod(repid_i,2) == 0) then
        ! 1  2> 3  4> 5  6> 7  8  ...
        next = +1
        if (cyclic) then
          !jump = -remd%nreplicas(dimno) + 1
        else
          jump = next
        end if

        if(repid_i == bias%ndocks) then
          repid_j = 0
          return
        end if

      else
        ! 1  2 <3  4 <5  6 <7  8  ...
        next = -1
        if (cyclic) then
          !jump = remd%nreplicas(dimno) - 1
        else
          jump = next
        end if

        if(repid_i == 1) then 
          repid_j = 0
          return
        end if

      end if      

    end if

    repid_j = repid_i + next

    return

  end subroutine get_neighbouring_parmsetid

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    abmd_bias_potential_exchange
  !> @brief        exchange abmd bias potential between two replicas
  !! @authors      JK, SI
  !! @param[in]    repi    : the ID for the first replica
  !! @param[in]    repj    : the ID for the second replica
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine abmd_bias_potential_exchange(repi, repj, enefunc, bias, do_exchange_bias)
    integer,                  intent(in)    :: repi, repj
    logical,                  intent(inout) :: do_exchange_bias           
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_bias),     target, intent(inout) :: bias

    integer                                 :: i, j
    integer                                 :: dockid
    integer                                 :: senddock, recvdock, dock2rep
    integer                                 :: size1, size2, size3
    integer                                 :: send_size, ierror
    integer,  pointer                       :: nreac_coord, nreac_coord_comb
    integer,  pointer                       :: nxi(:), kind(:)
    integer,  pointer                       :: bondslist(:,:,:)
    integer,  pointer                       :: bondslist_to_atomlist(:,:)
    integer,  pointer                       :: repsetid2dockid(:)
    integer,  pointer                       :: dockid2repsetid(:)
    real(dp), pointer                       :: uk(:)
    real(wp), pointer                       :: ximin(:), ximax(:), dxi(:), invdxi(:)
    logical,  pointer                       :: periodicity(:)
#ifdef HAVE_MPI_GENESIS
    integer                                 :: istatus(mpi_status_size)
#endif
    integer,  pointer                       :: temp_nreac_coord, temp_nreac_coord_comb
    integer,  pointer                       :: temp_nxi(:),   temp_kind(:)
    integer,  pointer                       :: temp_bondslist(:,:,:)
    integer,  pointer                       :: temp_bondslist_to_atomlist(:,:)
    real(wp), pointer                       :: temp_ximin(:), temp_ximax(:)
    real(wp), pointer                       :: temp_dxi(:),   temp_invdxi(:)
    logical,  pointer                       :: temp_periodicity(:)

    
    real(dp), allocatable                   :: unonzero(:)
    integer,  allocatable                   :: indexnonzero(:)
    integer                                 :: nnonzero
    integer                                 :: irequest0, irequest1, irequest2


    nreac_coord                => enefunc%abmd%nreac_coord
    nreac_coord_comb           => enefunc%abmd%nreac_coord_comb
    nxi                        => enefunc%abmd%nxi
    kind                       => enefunc%abmd%kind
    ximin                      => enefunc%abmd%ximin
    ximax                      => enefunc%abmd%ximax
    dxi                        => enefunc%abmd%dxi
    invdxi                     => enefunc%abmd%invdxi
    periodicity                => enefunc%abmd%periodicity
    bondslist                  => enefunc%abmd%bondslist
    bondslist_to_atomlist      => enefunc%abmd%bondslist_to_atomlist
    uk                         => enefunc%abmd%uk

    temp_nreac_coord           => enefunc%abmd%temp_nreac_coord
    temp_nreac_coord_comb      => enefunc%abmd%temp_nreac_coord_comb
    temp_nxi                   => enefunc%abmd%temp_nxi
    temp_kind                  => enefunc%abmd%temp_kind
    temp_ximin                 => enefunc%abmd%temp_ximin
    temp_ximax                 => enefunc%abmd%temp_ximax
    temp_dxi                   => enefunc%abmd%temp_dxi
    temp_invdxi                => enefunc%abmd%temp_invdxi
    temp_periodicity           => enefunc%abmd%temp_periodicity
    temp_bondslist             => enefunc%abmd%temp_bondslist
    temp_bondslist_to_atomlist => enefunc%abmd%temp_bondslist_to_atomlist

    repsetid2dockid            => bias%repsetid2dockid
    dockid2repsetid            => bias%dockid2repsetid_ref

    dockid = my_country_no + 1

    ! allocate array
    size1 = enefunc%reac_coord%max_atoms
    size2 = enefunc%reac_coord%num_groups
    size3 = maxval(enefunc%reac_coord%num_atoms_bonds_reac_coord)

    senddock = repj
    recvdock = repi

#ifdef HAVE_MPI_GENESIS
    if (dockid == repj) then
      do_exchange_bias = .true.
      ! Exchange repj-1 > repi-1
      ! Synchroniz
      !
      ! 1-1) nreac_coord
      call mpi_send(nreac_coord,      1, mpi_integer, recvdock-1, 0, &
                    mpi_comm_airplane, ierror)

      ! 1-2) nreac_coord_comb
      call mpi_send(nreac_coord_comb,      1, mpi_integer, recvdock-1, 0, &
                    mpi_comm_airplane, ierror)

      ! 2-1) nxi
      call mpi_send(nxi,      4, mpi_integer, recvdock-1, 0, mpi_comm_airplane, ierror)

      ! 2-2) kind
      call mpi_send(kind,      4, mpi_integer, recvdock-1, 0, mpi_comm_airplane, ierror)

      ! 2-3) ximin, ximax, dxi, invdxi
      call mpi_send(ximin,       4, mpi_wp_real, recvdock-1, 0, mpi_comm_airplane, ierror)
      call mpi_send(ximax,       4, mpi_wp_real, recvdock-1, 0, mpi_comm_airplane, ierror)
      call mpi_send(dxi,         4, mpi_wp_real, recvdock-1, 0, mpi_comm_airplane, ierror)
      call mpi_send(invdxi,      4, mpi_wp_real, recvdock-1, 0, mpi_comm_airplane, ierror)

      ! 2-4) periodicity
      call mpi_send(periodicity,      4, mpi_logical, recvdock-1, 0, mpi_comm_airplane, ierror)

      ! 2-5) bondslist, bondslist_to_atomlist
      send_size = size1 * size2 * 4
      call mpi_send(bondslist,      send_size, mpi_integer, recvdock-1, 0, mpi_comm_airplane, ierror)

      send_size = size3 * 4
      call mpi_send(bondslist_to_atomlist,      send_size, mpi_integer, recvdock-1, 0, mpi_comm_airplane, ierror)

      ! 3-1) uk
      call mpi_send(uk, nreac_coord_comb, mpi_double_precision, recvdock-1, 0, mpi_comm_airplane, ierror)
    endif

    if (dockid == repi) then
      do_exchange_bias = .true.

      ! Exchange repi-1 > repj-1
      !
      ! Synchronize
      !
      ! 1-1) nreac_coord
      call mpi_recv(temp_nreac_coord, 1, mpi_integer, senddock-1, 0, &
                    mpi_comm_airplane, istatus, ierror)
      
      ! 1-2) nreac_coord_comb
      call mpi_recv(temp_nreac_coord_comb, 1, mpi_integer, senddock-1, 0, &
                    mpi_comm_airplane, istatus, ierror)

      ! 2-1) nxi
      call mpi_recv(temp_nxi, 4, mpi_integer, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

      ! 2-2) kind
      call mpi_recv(temp_kind, 4, mpi_integer, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

      ! 2-3) ximin, ximax, dxi, invdxi
      call mpi_recv(temp_ximin,  4, mpi_wp_real, senddock-1, 0, mpi_comm_airplane, istatus, ierror)
      call mpi_recv(temp_ximax,  4, mpi_wp_real, senddock-1, 0, mpi_comm_airplane, istatus, ierror)
      call mpi_recv(temp_dxi,    4, mpi_wp_real, senddock-1, 0, mpi_comm_airplane, istatus, ierror)
      call mpi_recv(temp_invdxi, 4, mpi_wp_real, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

      ! 2-4) periodicity
      call mpi_recv(temp_periodicity, 4, mpi_logical, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

      ! 2-5) bondslist, bondslist_to_atomlist
      send_size = size1 * size2 * 4
      call mpi_recv(temp_bondslist, send_size, mpi_integer, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

      send_size = size3 * 4
      call mpi_recv(temp_bondslist_to_atomlist, send_size, mpi_integer, senddock-1, 0, mpi_comm_airplane, istatus, ierror)
      
      ! 3-1) uk
      allocate(enefunc%abmd%temp_uk(temp_nreac_coord_comb))
      call mpi_recv(enefunc%abmd%temp_uk, temp_nreac_coord_comb, &
                    mpi_double_precision, senddock-1, 0, mpi_comm_airplane, istatus, ierror)

    endif
#endif

  end subroutine abmd_bias_potential_exchange

end module at_bias_mod
