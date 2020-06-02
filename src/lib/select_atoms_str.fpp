!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_atoms_str_mod
!> @brief   structure for atom selection
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_atoms_str_mod

  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_selatoms
    integer, allocatable :: idx(:)
  end type s_selatoms

  ! subroutines
  public :: alloc_selatoms
  public :: dealloc_selatoms

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_selatoms
  !> @brief        allocate selatoms
  !! @authors      NT
  !! @param[inout] selatoms  : list of atom index
  !! @param[in]    list_size : size of list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_selatoms(selatoms, list_size)

    ! formal arguments
    type(s_selatoms),        intent(inout) :: selatoms
    integer,                 intent(in)    :: list_size

    ! local variables
    integer                  :: alloc_stat


    call dealloc_selatoms(selatoms)

    allocate(selatoms%idx(list_size), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    selatoms%idx(1:list_size) = 0

    return

  end subroutine alloc_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_selatoms
  !> @brief        deallocate selatoms
  !! @authors      NT
  !! @param[inout] selatoms  : list of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_selatoms(selatoms)

    ! formal arguments
    type(s_selatoms),        intent(inout) :: selatoms

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    if (allocated(selatoms%idx)) &
      deallocate(selatoms%idx, stat = dealloc_stat)

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_selatoms

end module select_atoms_str_mod
