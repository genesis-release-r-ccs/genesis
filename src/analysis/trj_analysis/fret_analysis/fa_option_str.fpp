!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fa_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fa_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    ! check only
    logical             :: check_only

    ! parameters for FRET efficiency
    character(MaxLineLong) :: analysis_atom_exp(2)
    type(s_selatoms)       :: analysis_atom(2)
    real(wp)               :: Forster_radius

  end type s_option

  public :: dealloc_option

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      TM
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option


    call dealloc_selatoms(option%analysis_atom(1))
    call dealloc_selatoms(option%analysis_atom(2))

    return

  end subroutine dealloc_option

end module fa_option_str_mod
