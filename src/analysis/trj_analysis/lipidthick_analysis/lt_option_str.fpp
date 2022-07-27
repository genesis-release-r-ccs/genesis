!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   lt_option_str_mod
!> @brief   structure of option information
!! @authors Daisuke Matsuoka (DM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module lt_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  type, public :: s_option
    logical                :: check_only

    type(s_selatoms)       :: membrane_atom
    character(MaxLineLong) :: membrane_atom_exp
  end type s_option

  ! subroutines
  public  :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      DM
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option


    call dealloc_selatoms(option%membrane_atom)

    return

  end subroutine dealloc_option

end module lt_option_str_mod
