!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_remd_str_mod
!> @brief   structure of replica information
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_remd_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_soltemp
    integer                           :: mygroup
    logical,              allocatable :: is_solute(:)

    integer                           :: num_solute
    integer,              allocatable :: solute_list(:)

    logical                           :: sw_charge
    integer                           :: num_mixed_three
    real(wp),             allocatable :: charge_org(:)

    logical                           :: sw_bonds
    integer                           :: num_bonds
    integer,              allocatable :: bond_list(:)
    real(wp),             allocatable :: bond_weight(:)
    real(wp),             allocatable :: bond_force_const_org(:)

    logical                           :: sw_angles
    integer                           :: num_angles
    integer,              allocatable :: angl_list(:)
    real(wp),             allocatable :: angl_weight(:)
    real(wp),             allocatable :: angl_force_const_org(:)

    logical                           :: sw_ureys
    integer                           :: num_ureys
    integer,              allocatable :: urey_list(:)
    real(wp),             allocatable :: urey_weight(:)
    real(wp),             allocatable :: urey_force_const_org(:)

    logical                           :: sw_dihedrals
    integer                           :: num_dihedrals
    integer,              allocatable :: dihe_list(:)
    real(wp),             allocatable :: dihe_weight(:)
    real(wp),             allocatable :: dihe_force_const_org(:)

    logical                           :: sw_impropers
    integer                           :: num_impropers
    integer,              allocatable :: impr_list(:)
    real(wp),             allocatable :: impr_weight(:)
    real(wp),             allocatable :: impr_force_const_org(:)

    logical                           :: sw_cmaps
    integer                           :: num_cmaps
    integer                           :: num_cmap_type
    integer                           :: istart_cmap_type
    integer,              allocatable :: cmap_type_org(:)
    real(wp),             allocatable :: cmap_weight(:)

    logical                           :: sw_lj
    integer                           :: num_atom_cls
    integer                           :: istart_atom_cls
    integer,              allocatable :: atom_cls_no_org(:)

    logical                           :: sw_contacts
    integer                           :: num_contacts
    integer,              allocatable :: contact_list(:)
    real(wp),             allocatable :: contact_weight(:)
    real(wp),             allocatable :: contact_lj6_org(:)
    real(wp),             allocatable :: contact_lj10_org(:)
    real(wp),             allocatable :: contact_lj12_org(:)
  end type s_soltemp

  type, public :: s_adjust_remd
    logical                           :: param_tuning
    real(wp)                          :: tgt_exc_prob
    real(wp)                          :: mgn_exc_prob
    integer                           :: trial_freq
    integer                           :: eq_cycle
    real(wp)                          :: param_grid
    real(wp)                          :: max_param_shift
    integer                           :: fix_terminal
  end type s_adjust_remd

  type, public :: s_remd
    integer                           :: total_nreplicas
    integer                           :: ncycles
    integer                           :: dimension
    integer                           :: exchange_period
    integer                           :: iseed
    logical                           :: equilibration_only
    logical                           :: rest_mixed_three_atoms
    integer,              allocatable :: types(:)
    integer,              allocatable :: nreplicas(:)
    integer,              allocatable :: repid2parmsetid(:)
    integer,              allocatable :: repid2parmsetid_ref(:)
    integer,              allocatable :: parmsetid2repid(:)
    integer,              allocatable :: parmidsets(:,:)
    integer,              allocatable :: num_criteria(:,:,:)
    integer,              allocatable :: num_exchanges(:,:,:)
    integer,              allocatable :: rest_function(:)
    integer,              allocatable :: umbrid2funcid(:)
    integer,              allocatable :: iparameters(:,:)
    integer,              allocatable :: umbrid2numfuncs(:)
    integer,              allocatable :: umbrid2funclist(:,:)
    real(wp),             allocatable :: dparameters(:,:)
    real(wp),             allocatable :: rest_constants(:,:)
    real(wp),             allocatable :: rest_reference(:,:)
    logical,              allocatable :: rest_caging(:,:)
    real(wp),             allocatable :: rest_flat_radius(:,:)
    real(wp),             allocatable :: potential_energy(:)
    real(wp),             allocatable :: restraint_energy_0(:)
    real(wp),             allocatable :: restraint_energy_1(:)
    real(wp),             allocatable :: restraint_refcoord0(:,:)
    real(wp),             allocatable :: restraint_refcoord1(:,:)
    real(wp),             allocatable :: volume(:)
    real(wp),             allocatable :: area(:)
    real(wp),             allocatable :: after_gather (:)
    real(wp),             allocatable :: random_table (:,:)
    logical,              allocatable :: cyclic_params(:)
    type(s_soltemp),      allocatable :: solute_tempering(:)
    type(s_adjust_remd),  allocatable :: autoadj(:)
  end type s_remd

  ! parameters for auto-adjusting
  integer, public, parameter      :: RemdAutoFixBottom = 1
  integer, public, parameter      :: RemdAutoFixTop    = 2

  ! parameters for allocatable variables
  integer, public, parameter      :: RemdReplicas    = 1
  integer, public, parameter      :: RemdUmbrellas   = 2
  integer, public, parameter      :: RemdRefatoms    = 3

  ! parameters for REMD type
  integer, public, parameter      :: RemdTemperature     = 1
  integer, public, parameter      :: RemdPressure        = 2
  integer, public, parameter      :: RemdGamma           = 3
  integer, public, parameter      :: RemdRestraint       = 4
  integer, public, parameter      :: RemdSoluteTempering = 5

  character(5), public, parameter :: DummyWaterName      = 'ZZZZZ'
  character(*), public, parameter :: RemdTypes(5)  = (/'TEMPERATURE',&
                                                       'PRESSURE   ',&
                                                       'GAMMA      ',&
                                                       'RESTRAINT  ',&
                                                       'REST       '/)

  ! subroutines
  public :: alloc_remd
  public :: dealloc_remd
  public :: dealloc_remd_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_remd
  !> @brief        allocate remd information
  !! @authors      TM
  !! @param[inout] remd      : information of remd
  !! @param[in]    variable  : allocatable variable
  !! @param[in]    var_size1 : size of variables (REMD dimension)
  !! @param[in]    var_size2 : size of variables (Total number of replicas)
  !! @param[in]    var_size3 : size of variables (max number of replica)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_remd(remd, variable, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_remd),            intent(inout) :: remd
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,      optional,  intent(in)    :: var_size2
    integer,      optional,  intent(in)    :: var_size3

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RemdReplicas)

      if (allocated(remd%nreplicas)) then
        if (size(remd%types(:)) == var_size1) return
        deallocate(remd%types,               &
                   remd%nreplicas,           &
                   remd%cyclic_params,       &
                   remd%repid2parmsetid,     &
                   remd%repid2parmsetid_ref, &
                   remd%parmsetid2repid,     &
                   remd%parmidsets,          &
                   remd%num_criteria,        &
                   remd%num_exchanges,       &
                   remd%iparameters,         &
                   remd%dparameters,         &
                   remd%rest_function,       &
                   remd%potential_energy,    &
                   remd%restraint_energy_0,  &
                   remd%restraint_energy_1,  &
                   remd%volume,              &
                   remd%area,                &
                   remd%after_gather,        &
                   remd%random_table,        &
                   remd%solute_tempering,    &
                   remd%autoadj,             &
                   stat = dealloc_stat)
      end if

      allocate(remd%types(var_size1),                     &
               remd%nreplicas(var_size1),                 &
               remd%cyclic_params(var_size1),             &
               remd%repid2parmsetid(var_size2),           &
               remd%repid2parmsetid_ref(var_size2),       &
               remd%parmsetid2repid(var_size2),           &
               remd%parmidsets(var_size2,var_size1),      &
               remd%num_criteria(var_size2,var_size1,2),  &
               remd%num_exchanges(var_size2,var_size1,2), & 
               remd%iparameters(var_size1,var_size3),     &
               remd%dparameters(var_size1,var_size3),     &
               remd%rest_function(var_size1),             &
               remd%potential_energy(var_size2),          &
               remd%restraint_energy_0(var_size2),        &
               remd%restraint_energy_1(var_size2),        &
               remd%volume(var_size2),                    &
               remd%area(var_size2),                      &
               remd%after_gather(var_size2),              &
               remd%random_table(var_size2,var_size2),    &
               remd%solute_tempering(var_size1),          &
               remd%autoadj(var_size1),                   &
               stat = alloc_stat)

      remd%types(1:var_size1)                         = 1
      remd%nreplicas(1:var_size1)                     = 0
      remd%cyclic_params(1:var_size1)                 = .false.
      remd%repid2parmsetid(1:var_size2)               = 0
      remd%repid2parmsetid_ref(1:var_size2)           = 0
      remd%parmsetid2repid(1:var_size2)               = 0
      remd%parmidsets(1:var_size2,1:var_size1)        = 0
      remd%num_criteria(1:var_size2,1:var_size1,1:2)  = 0
      remd%num_exchanges(1:var_size2,1:var_size1,1:2) = 0
      remd%iparameters(1:var_size1,1:var_size3)       = 0
      remd%dparameters(1:var_size1,1:var_size3)       = 0.0_wp
      remd%rest_function(1:var_size1)                 = 0

    case(RemdUmbrellas)

      if (allocated(remd%umbrid2numfuncs)) then
        if (size(remd%umbrid2numfuncs(:)) == var_size1) return
        deallocate(remd%umbrid2numfuncs,                  &
                   remd%umbrid2funclist,                  &
                   remd%rest_constants,                   &
                   remd%rest_reference,                   &
                   remd%rest_caging,                      &
                   remd%rest_flat_radius,                 &
                   stat = dealloc_stat)
      end if

      allocate(remd%umbrid2numfuncs(var_size1),            &
               remd%umbrid2funclist(var_size1,var_size2),  &
               remd%rest_constants(var_size1,var_size2),   &
               remd%rest_reference(var_size1,var_size2),   &
               remd%rest_caging(var_size1,var_size2),      &
               remd%rest_flat_radius(var_size1,var_size2), &
               stat = alloc_stat)

      remd%umbrid2numfuncs(1:var_size1)               = 0
      remd%umbrid2funclist(1:var_size1,1:var_size2)   = 0
      remd%rest_constants(1:var_size1,1:var_size2)    = 0_wp
      remd%rest_reference(1:var_size1,1:var_size2)    = 0_wp
      remd%rest_caging(1:var_size1,1:var_size2)       = .false.
      remd%rest_flat_radius(1:var_size1,1:var_size2)  = 0_wp

    case default

      call error_msg('Alloc_Remd> bad variable')

    end select


    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_remd
  !> @brief        deallocate remd information
  !! @authors      TM
  !! @param[inout] remd     : remd information
  !! @param[in]    variable : allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_remd(remd, variable)

    ! formal arguments
    type(s_remd),            intent(inout) :: remd
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (RemdReplicas)

      if (allocated(remd%types)) then
        deallocate(remd%types,               &
                   remd%nreplicas,           &
                   remd%cyclic_params,       &
                   remd%repid2parmsetid,     &
                   remd%repid2parmsetid_ref, &
                   remd%parmsetid2repid,     &
                   remd%parmidsets,          &
                   remd%num_criteria,        &
                   remd%num_exchanges,       &
                   remd%iparameters,         &
                   remd%dparameters,         &
                   remd%rest_function,       &
                   remd%potential_energy,    &
                   remd%restraint_energy_0,  &
                   remd%restraint_energy_1,  &
                   remd%volume,              &
                   remd%area,                &
                   remd%after_gather,        &
                   remd%random_table,        &
                   remd%solute_tempering,    &
                   remd%autoadj,             &
                   stat = dealloc_stat)
      end if

    case (RemdUmbrellas)

      if (allocated(remd%umbrid2funclist)) then
        deallocate(remd%umbrid2funclist,     &
                   remd%umbrid2numfuncs,     &
                   remd%rest_constants,      &
                   remd%rest_reference,      &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Remd> bad variable')

    end select 


    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_remd_all
  !> @brief        deallocate all remd information
  !! @authors      TM
  !! @param[inout] remd : information of remd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_remd_all(remd)

    ! formal arguments
    type(s_remd),         intent(inout) :: remd


    call dealloc_remd(remd, RemdReplicas)
    call dealloc_remd(remd, RemdUmbrellas)

    return

  end subroutine dealloc_remd_all

end module at_remd_str_mod
