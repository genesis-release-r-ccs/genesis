!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cc_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module cc_option_str_mod

  use constants_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical             :: check_only
    integer             :: trjout_format
    integer             :: trjout_type
    character(MaxLine)  :: trjout_atom_exp
    type(s_selatoms)    :: trjout_atom
    type(s_selatoms)    :: trjout_atom_trj
    logical             :: centering
    character(MaxLine)  :: centering_atom_exp
    type(s_selatoms)    :: centering_atom
    real(wp)            :: center_coord(3)
    logical             :: split_trjpdb
    integer             :: pbcc_mode
    integer             :: num_groups
    integer, allocatable:: trjout_groups(:)
    integer, allocatable:: num_atoms(:)
    integer, allocatable:: atomlist(:,:)
    logical             :: geometric_flag

  end type s_option

  ! subroutines
  public  :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      NT
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option


    call dealloc_selatoms(option%trjout_atom)
    call dealloc_selatoms(option%trjout_atom_trj)
    call dealloc_selatoms(option%centering_atom)

    return

  end subroutine dealloc_option

end module cc_option_str_mod
