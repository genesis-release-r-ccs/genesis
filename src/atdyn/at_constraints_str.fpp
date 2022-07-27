!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_constraints_str_mod
!> @brief   structure of constraints information
!! @authors Takaharu Mori (TM), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_constraints_str_mod

  use constants_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_constraints
    logical                       :: rigid_bond
    logical                       :: fast_bond
    logical                       :: fast_water
    integer                       :: hydrogen_type

    ! for SHAKE, RATTLE, and LINCS
    integer                       :: shake_iteration
    real(wp)                      :: shake_tolerance
    integer                       :: num_bonds
    integer,          allocatable :: bond_list(:,:)
    real(wp),         allocatable :: bond_dist(:)
    real(dp),         allocatable :: bond_vector(:,:)
    real(wp),         allocatable :: shake_force(:)

    ! for LINCS
    integer                       :: lincs_iteration
    integer                       :: lincs_order
    integer,          allocatable :: num_connected_cons(:)
    integer,          allocatable :: connected_cons_idx(:,:)
    real(wp),         allocatable :: lincs_coef(:,:)
    real(wp),         allocatable :: s_diagonal(:)
    real(wp),         allocatable :: cons_coupl_mat(:,:)
    real(wp),         allocatable :: right_hand_side(:,:)
    real(wp),         allocatable :: solution_array(:)
    real(wp),         allocatable :: bond_length(:)

    ! for SETTLE
    character(5)                  :: water_model
    integer                       :: num_water
    integer,          allocatable :: water_list(:,:)
    real(wp)                      :: water_rHH
    real(wp)                      :: water_rOH
    real(wp)                      :: water_massO
    real(wp)                      :: water_massH
    real(wp),         allocatable :: virial_tmp(:,:,:)

    ! for fixed atoms
    integer                       :: num_fixatm
    logical, allocatable          :: fixatm(:)

  end type s_constraints

  ! parameters for allocatable variables
  integer,      public, parameter :: ConstraintsShake  = 1
  integer,      public, parameter :: ConstraintsLincs  = 2
  integer,      public, parameter :: ConstraintsSettle = 3
  integer,      public, parameter :: ConstraintsFixatm = 4

  ! parameters
  integer,      public, parameter :: ConstraintAtomName  = 1
  integer,      public, parameter :: ConstraintAtomMass  = 2
  integer,      public, parameter :: ConstraintAtomBoth  = 3

  character(*), public, parameter :: ConstraintAtomType(3)  = (/'NAME', &
                                                                'MASS', &
                                                                'BOTH'/)
  ! subroutines
  public :: init_constraints
  public :: alloc_constraints
  public :: dealloc_constraints
  public :: dealloc_constraints_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_constraints
  !> @brief        initialize constraints information
  !! @authors      TM, KY
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
    constraints%shake_tolerance = 0.0_wp
    constraints%num_bonds       = 0
    constraints%lincs_iteration = 0
    constraints%lincs_order     = 0
    constraints%water_model     = ''
    constraints%num_water       = 0
    constraints%water_rHH       = 0.0_wp
    constraints%water_rOH       = 0.0_wp
    constraints%water_massO     = 0.0_wp
    constraints%water_massH     = 0.0_wp
    constraints%num_fixatm      = 0

    return

  end subroutine init_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_constraints
  !> @brief        allocate constraints information
  !! @authors      TM, KY
  !! @param[inout] constraints : constraints information
  !! @param[in]    variable    : selected variable
  !! @param[in]    var_size    : size of the selected variable
  !! @param[in]    var_size2   : 2nd size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_constraints(constraints, variable, var_size1, var_size2)

    ! formal arguments
    type(s_constraints),     intent(inout) :: constraints
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer, optional,       intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (ConstraintsShake)

      if (allocated(constraints%bond_list)) then
        if (size(constraints%bond_list(1,:)) == var_size1) return
        deallocate(constraints%bond_list,   &
                   constraints%bond_dist,   &
                   constraints%bond_vector, &
                   constraints%shake_force, &
                   stat = dealloc_stat)
      end if

      allocate(constraints%bond_list(2,var_size1),   &
               constraints%bond_dist(var_size1),     &
               constraints%bond_vector(3,var_size1), &
               constraints%shake_force(var_size1),   &
               stat = alloc_stat)

      constraints%bond_list  (1:2,1:var_size1) = 0
      constraints%bond_dist      (1:var_size1) = 0.0_wp
      constraints%bond_vector(1:3,1:var_size1) = 0.0_dp
      constraints%shake_force    (1:var_size1) = 0.0_wp

    case (ConstraintsLincs)

      if (present(var_size2)) then

        if (allocated(constraints%num_connected_cons)) then
          if (size(constraints%num_connected_cons(:)) == var_size1) return
          deallocate(constraints%num_connected_cons, &
                     constraints%connected_cons_idx, &
                     constraints%lincs_coef,         &
                     constraints%s_diagonal,         &
                     constraints%cons_coupl_mat,     &
                     constraints%right_hand_side,    &
                     constraints%solution_array,     &
                     constraints%bond_length,        &
                     stat = dealloc_stat)
        end if

        allocate(constraints%num_connected_cons(var_size1),           &
                 constraints%connected_cons_idx(var_size2,var_size1), &
                 constraints%lincs_coef(var_size2,var_size1),         &
                 constraints%s_diagonal(var_size1),                   &
                 constraints%cons_coupl_mat(var_size2,var_size1),     &
                 constraints%right_hand_side(2,var_size1),            &
                 constraints%solution_array(var_size1),               &
                 constraints%bond_length(var_size1),                  &
                 stat = alloc_stat)

        constraints%num_connected_cons(1:var_size1) = 0
        constraints%connected_cons_idx(1:var_size2,1:var_size1) = 0
        constraints%lincs_coef(1:var_size2,1:var_size1) = 0.0_wp
        constraints%s_diagonal(1:var_size1) = 0.0_wp
        constraints%cons_coupl_mat(1:var_size2,1:var_size1) = 0.0_wp
        constraints%right_hand_side(1:2,1:var_size1) = 0.0_wp
        constraints%solution_array(1:var_size1) = 0.0_wp
        constraints%bond_length(1:var_size1) = 0.0_wp

      else

        call error_msg('Alloc_Constraints> allocation error for var_size2')

      end if

    case (ConstraintsSettle)

      if (allocated(constraints%water_list)) then
        if (size(constraints%water_list(1,:)) == var_size1) return
        deallocate(constraints%water_list, &
                   stat = dealloc_stat)
      end if

      allocate(constraints%water_list(3,var_size1), &
               stat = alloc_stat)

      constraints%water_list(1:3,1:var_size1) = 0

    case (ConstraintsFixatm)

      if (allocated(constraints%fixatm)) then
        if (size(constraints%fixatm(:)) == var_size1) return
        deallocate(constraints%fixatm, &
                   stat = dealloc_stat)
      end if

      allocate(constraints%fixatm(var_size1), &
               stat = alloc_stat)

      constraints%fixatm(1:var_size1) = .false.

    case default

      call error_msg('Alloc_Constraints> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc


    return

  end subroutine alloc_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_constraints
  !> @brief        deallocate constraints information
  !! @authors      TM, KY
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

    ! deallocate selected variable
    !
    select case (variable)

    case (ConstraintsShake)

      if (allocated(constraints%bond_list)) then
        deallocate (constraints%bond_list,   &
                    constraints%bond_dist,   &
                    constraints%bond_vector, &
                    constraints%shake_force, &
                    stat = dealloc_stat)
      end if

    case (ConstraintsLincs)

      if (allocated(constraints%num_connected_cons)) then
        deallocate (constraints%num_connected_cons, &
                    constraints%connected_cons_idx, &
                    constraints%lincs_coef,         &
                    constraints%s_diagonal,         &
                    constraints%cons_coupl_mat,     &
                    constraints%right_hand_side,    &
                    constraints%solution_array,     &
                    constraints%bond_length,        &
                    stat = dealloc_stat)
      end if

    case (ConstraintsSettle)

      if (allocated(constraints%water_list)) then
        deallocate (constraints%water_list, &
                    stat = dealloc_stat)
      end if

    case (ConstraintsFixatm)

      if (allocated(constraints%fixatm)) then
        deallocate (constraints%fixatm, &
                    stat = dealloc_stat)
      end if

    case default
      call error_msg('Dealloc_Constraints> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc


    return

  end subroutine dealloc_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_constraints_all
  !> @brief        deallocate all constraints information
  !! @authors      TM
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_constraints_all(constraints)

    ! format arguments
    type(s_constraints),     intent(inout) :: constraints


    call dealloc_constraints(constraints, ConstraintsShake)
    call dealloc_constraints(constraints, ConstraintsLincs)
    call dealloc_constraints(constraints, ConstraintsSettle)
    call dealloc_constraints(constraints, ConstraintsFixatm)

    return

  end subroutine dealloc_constraints_all

end module at_constraints_str_mod
