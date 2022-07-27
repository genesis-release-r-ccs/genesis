!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dr_option_str_mod
!> @brief   structure of option information
!! @authors Chigusa Kobayashi (CK), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module dr_option_str_mod

  use messages_mod
  use string_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    logical                       :: verbose
    integer                       :: num_contact
    integer,allocatable           :: contact_list(:,:)
    integer,allocatable           :: exclude_group_list(:)
    integer,allocatable           :: contact_atoms(:,:)
    integer                       :: num_atoms_group(1:2)
    logical                       :: identical_group
    logical                       :: avoid_bonding
    logical                       :: pbc_correct_setup
    logical                       :: pbc_correct
    logical                       :: ignore_hydrogen
    logical                       :: two_states
    integer                       :: exclude_residues
    real(wp),allocatable          :: contact_dist(:)
    real(wp),allocatable          :: contact_cur_dist(:)
    real(wp)                      :: minimum_difference
    real(wp)                      :: minimum_distance
    real(wp)                      :: maximum_distance
    real(wp)                      :: box_size_cur(1:3)
    real(wp)                      :: box_size_ref(1:3)
    
  end type s_option

  ! parameters for allocatable variables
  integer,        public, parameter :: DA_Contact    = 1
  integer,        public, parameter :: DA_Exclude    = 2
  integer,        public, parameter :: DA_Atoms      = 3

  public :: alloc_option
  public :: dealloc_option
  public :: dealloc_option_all

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_option
  !> @brief        allocate option
  !! @authors      CK
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine alloc_option(option, variable, var_size)

    ! formal argments
    type(s_option),          intent(inout) :: option
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    integer                                :: alloc_stat
    integer                                :: dealloc_stat

    dealloc_stat = 0
    alloc_stat = 0

    select case (variable)

    case (DA_Contact)

      if (allocated(option%contact_dist)) then
        if (size(option%contact_dist) == var_size) return
        deallocate(option%contact_dist, &
                   option%contact_cur_dist, &
                   option%contact_list, &
                   stat = dealloc_stat)
      end if
     
      allocate(option%contact_dist(1:var_size),      &
               option%contact_cur_dist(1:var_size),  &
               option%contact_list(1:2, 1:var_size), &
               stat = alloc_stat)
      option%contact_dist(1:var_size)     = 0.0_wp
      option%contact_cur_dist(1:var_size) = 0.0_wp
      option%contact_list(1:2,1:var_size) = 0

    case (DA_Exclude)
      if (allocated(option%exclude_group_list)) then
        if (size(option%exclude_group_list) == var_size) return
        deallocate(option%exclude_group_list, &
                   stat = dealloc_stat)
      end if
     
      allocate(option%exclude_group_list(1:var_size), &
               stat = alloc_stat)
      option%exclude_group_list(1:var_size) = 0

    case (DA_Atoms)
      if (allocated(option%contact_atoms)) then
        if (size(option%contact_atoms(:,1)) == var_size) return
        deallocate(option%contact_atoms, &
                   stat = dealloc_stat)
      end if
     
      allocate(option%contact_atoms(1:var_size,1:2), &
               stat = alloc_stat)
      option%contact_atoms(1:var_size,1:2) = 0

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      CK
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option,variable)

    ! formal argments
    type(s_option),          intent(inout) :: option
    integer,                 intent(in)    :: variable

    integer                                :: dealloc_stat

    dealloc_stat = 0

    select case (variable)

    case (DA_Contact)

      deallocate(option%contact_dist,     &
                 option%contact_cur_dist, &
                 option%contact_list,     &
                 stat = dealloc_stat)

    case (DA_Exclude)
     
      deallocate(option%exclude_group_list, &
                stat = dealloc_stat)

    case (DA_Atoms)

      deallocate(option%contact_atoms, &
                 stat = dealloc_stat)

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine dealloc_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option_all
  !> @brief        deallocate all
  !! @authors      CK
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option_all(option)
    ! formal argments
    type(s_option),          intent(inout) :: option

    call dealloc_option(option,DA_Atoms)
    call dealloc_option(option,DA_Contact)
    call dealloc_option(option,DA_Exclude)

    return

  end subroutine dealloc_option_all

end module dr_option_str_mod
