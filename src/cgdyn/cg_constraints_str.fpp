!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_constraints_str_mod
!> @brief   Structure of constraints information.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_constraints_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  !
  type, public :: s_constraints
    logical                       :: rigid_bond
    logical                       :: fast_bond
    logical                       :: fast_water
    integer                       :: hydrogen_type
    real(wp)                      :: hydrogen_mass_upper_bound

    ! for SHAKE, RATTLE, and LINCS
    integer                       :: shake_iteration
    real(dp)                      :: shake_tolerance
    integer                       :: num_bonds
    integer                       :: connect
    integer                       :: nh(8)

    !  HBond
    integer,          allocatable :: duplicate(:)
    integer,          allocatable :: H_index(:,:)
    !  BondGroup
    integer,          allocatable :: H_Group(:,:,:)
    !  DomainBond
    integer,          allocatable :: No_HGr(:)
    integer,          allocatable :: HGr_local(:,:)
    integer,          allocatable :: HGr_bond_list(:,:,:,:)
    real(dp),         allocatable :: HGr_bond_dist(:,:,:,:)
    real(dp),         allocatable :: HGr_bond_vector(:,:,:,:,:)
    real(dp),         allocatable :: HGr_shake_force(:,:,:,:)
    !  HGroupMove
    integer,          allocatable :: HGr_move(:,:)
    integer,          allocatable :: HGr_stay(:,:)
    integer,          allocatable :: HGr_move_int(:,:,:,:)
    integer,          allocatable :: HGr_stay_int(:,:,:,:)
    real(dp),         allocatable :: HGr_move_real(:,:,:,:)
    real(dp),         allocatable :: HGr_stay_real(:,:,:,:)

    ! for LINCS
    integer                       :: lincs_iteration
    integer                       :: lincs_order

    ! for SETTLE
    logical                       :: tip4
    character(5)                  :: water_model
    real(dp)                      :: water_rHH
    real(dp)                      :: water_rOH
    real(dp)                      :: water_rOD
    real(dp)                      :: water_massO
    real(dp)                      :: water_massH

  end type s_constraints

  ! parameters for allocatable variables
  integer,      public, parameter :: ConstraintsHBond      = 1
  integer,      public, parameter :: ConstraintsBondGroup  = 2
  integer,      public, parameter :: ConstraintsDomainBond = 3
  integer,      public, parameter :: ConstraintsHGroupMove = 4

  ! parameters
  integer,      public, parameter :: ConstraintAtomName  = 1
  integer,      public, parameter :: ConstraintAtomMass  = 2
  integer,      public, parameter :: ConstraintAtomBoth  = 3

  character(*), public, parameter :: ConstraintAtomType(3)  = (/'NAME', &
                                                                'MASS', &
                                                                'BOTH'/)

  ! parameters
  integer,      public            :: HGroupMax   = 150
  integer,      public            :: HGrpMaxMove = 20

  ! subroutines
  public :: init_constraints
  public :: alloc_constraints
  public :: dealloc_constraints
  public :: dealloc_constraints_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_constraints
  !> @brief        Initializes constraints information.
  !! @authors      JJ
  !! @param[out]   constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_constraints(constraints)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints


    constraints%rigid_bond      = .false.
    constraints%fast_bond       = .false.
    constraints%fast_water      = .false.
    constraints%hydrogen_type   = ConstraintAtomName
    constraints%shake_iteration = 0
    constraints%shake_tolerance = 0.0_dp
    constraints%num_bonds       = 0
    constraints%connect         = 0
    constraints%nh(1:8)         = 0
    constraints%lincs_iteration = 0
    constraints%lincs_order     = 0
    constraints%water_rHH       = 0.0_dp
    constraints%water_rOH       = 0.0_dp
    constraints%water_rOD       = 0.0_dp
    constraints%water_massO     = 0.0_dp
    constraints%water_massH     = 0.0_dp

    return

  end subroutine init_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_constraints
  !> @brief        Allocates constraints information.
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !! @param[in]    variable    : selected variable
  !! @param[in]    var_size    : size of the selected variable
  !! @param[in]    var_size2   : 2nd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_constraints(constraints, variable, var_size, var_size1)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer, optional,       intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: connect


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(ConstraintsHBond)

      if (allocated(constraints%duplicate)) then
        if (size(constraints%duplicate) /= var_size) &
          deallocate(constraints%duplicate, &
                     constraints%H_index,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(constraints%duplicate))    &
        allocate(constraints%duplicate(var_size),    &
                 constraints%H_index  (8, var_size), &
                 stat = alloc_stat)
    
      constraints%duplicate(1:var_size)      = 0
      constraints%H_index  (1:8, 1:var_size) = 0

    case (ConstraintsBondGroup)

      connect = constraints%connect 

      if (allocated(constraints%H_Group)) then
        if (size(constraints%H_Group) /= (connect+1)*var_size*connect) &
          deallocate(constraints%H_Group,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(constraints%H_Group)) &
        allocate(constraints%H_Group  (connect+1, var_size, connect), &
                 stat = alloc_stat)

      constraints%H_Group  (1:connect+1, 1:var_size, 1:connect) = 0

    case (ConstraintsDomainBond)

      if (allocated(constraints%HGr_local)) then
        if (size(constraints%HGr_local) /= var_size*var_size1) &
          deallocate(constraints%No_HGr,          &
                     constraints%HGr_local,       &
                     constraints%HGr_bond_list,   &
                     constraints%HGr_bond_dist,   &
                     constraints%HGr_bond_vector, &
                     constraints%HGr_shake_force, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(constraints%HGr_local)) &
        allocate(constraints%No_HGr                                    &
                      (var_size),                                      &
                 constraints%HGr_local                                 &
                      (var_size1, var_size),                           &
                 constraints%HGr_bond_list                             &
                      (var_size1+1, HGroupMax, var_size1, var_size),   &
                 constraints%HGr_bond_dist                             &
                      (var_size1+1, HGroupMax, var_size1, var_size),   &
                 constraints%HGr_bond_vector                           &
                      (3, var_size1, HGroupMax, var_size1, var_size),  &
                 constraints%HGr_shake_force                           &
                      (var_size1, HGroupMax, var_size1, var_size),     &
                 stat = alloc_stat)

      constraints%No_HGr          &
           (1:var_size)                                             = 0
      constraints%HGr_local       &
           (1:var_size1, 1:var_size)                                = 0
      constraints%HGr_bond_list   &
           (1:var_size1+1, 1:HGroupMax, 1:var_size1, 1:var_size)    = 0
      constraints%HGr_bond_dist   &
           (1:var_size1+1, 1:HGroupMax, 1:var_size1, 1:var_size)    = 0.0_dp
      constraints%HGr_bond_vector &
           (1:3, 1:var_size1, 1:HGroupMax, 1:var_size1, 1:var_size) = 0.0_dp
      constraints%HGr_shake_force &
           (1:var_size1,   1:HGroupMax, 1:var_size1, 1:var_size)    = 0.0_dp

    case (ConstraintsHGroupMove)

      if (allocated(constraints%HGr_move)) then
        if (size(constraints%HGr_move) /= var_size*var_size1) &
          deallocate(constraints%HGr_move,      &
                     constraints%HGr_stay,      &
                     constraints%HGr_move_int,  &
                     constraints%HGr_stay_int,  &
                     constraints%HGr_move_real, &
                     constraints%HGr_stay_real, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(constraints%HGr_move)) &
        allocate(constraints%HGr_move                                    &
                      (var_size, var_size1),                             &
                 constraints%HGr_stay                                    &
                      (var_size, var_size1),                             &
                 constraints%HGr_move_int                                &
                      (2*(var_size+1), HGrpMaxMove, var_size, var_size1),&
                 constraints%HGr_stay_int                                &
                      (2*(var_size+1), HGroupMax, var_size, var_size1),  &
                 constraints%HGr_move_real                               &
                      (9*(var_size+1), HGrpMaxMove, var_size, var_size1),&
                 constraints%HGr_stay_real                               &
                      (9*(var_size+1), HGroupMax, var_size, var_size1),  &
                 stat = alloc_stat)

      constraints%HGr_move      &
           (1:var_size, 1:var_size1)                                   = 0
      constraints%HGr_stay      &
           (1:var_size, 1:var_size1)                                   = 0
      constraints%HGr_move_int  &
           (1:2*(var_size+1), 1:HGrpMaxMove, 1:var_size, 1:var_size1)  = 0
      constraints%HGr_stay_int  &
           (1:2*(var_size+1), 1:HGroupMax, 1:var_size, 1:var_size1)    = 0
      constraints%HGr_move_real &
           (1:9*(var_size+1), 1:HGrpMaxMove, 1:var_size, 1:var_size1)  = 0.0_dp
      constraints%HGr_stay_real &
           (1:9*(var_size+1), 1:HGroupMax, 1:var_size, 1:var_size1)    = 0.0_dp

    case default

      call error_msg('Error: Alloc_Constraints - Invalid variable specified.')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_constraints
  !> @brief        Deallocates constraints information.
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !! @param[in]    variable    : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_constraints(constraints, variable)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(ConstraintsHBond)

      if (allocated(constraints%duplicate)) then
        deallocate(constraints%duplicate, &
                   constraints%H_index,   &
                   stat = dealloc_stat)
      end if

    case (ConstraintsBondGroup)

      if (allocated(constraints%H_Group)) then
        deallocate(constraints%H_Group,   &
                   stat = dealloc_stat)
      end if

    case (ConstraintsDomainBond)

      if (allocated(constraints%HGr_local)) then
        deallocate(constraints%No_HGr,          &
                   constraints%HGr_local,       &
                   constraints%HGr_bond_list,   &
                   constraints%HGr_bond_dist,   &
                   constraints%HGr_bond_vector, &
                   constraints%HGr_shake_force, &
                   stat = dealloc_stat)
      end if

    case (ConstraintsHGroupMove)

      if (allocated(constraints%HGr_move)) then
        deallocate(constraints%HGr_move,      &
                   constraints%HGr_stay,      &
                   constraints%HGr_move_int,  &
                   constraints%HGr_stay_int,  &
                   constraints%HGr_move_real, &
                   constraints%HGr_stay_real, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Error: Dealloc_Constraints - Invalid variable specified.')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_constraints_all
  !> @brief        Deallocates all constraints information.
  !! @authors      JJ
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_constraints_all(constraints)

    ! format arguments
    type(s_constraints),     intent(inout) :: constraints


    call dealloc_constraints(constraints, ConstraintsHBond)
    call dealloc_constraints(constraints, ConstraintsBondGroup)
    call dealloc_constraints(constraints, ConstraintsDomainBond)
    call dealloc_constraints(constraints, ConstraintsHGroupMove)

    return

  end subroutine dealloc_constraints_all

end module cg_constraints_str_mod
