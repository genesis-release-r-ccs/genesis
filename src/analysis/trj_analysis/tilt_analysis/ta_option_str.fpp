!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_option_str_mod
!> @brief   structure of option information
!! @authors Daisuke Matsuoka (DM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  type, public :: s_parameter
    logical :: tilt
  end type s_parameter


  type, public :: s_option
    logical                :: check_only

    character(MaxLineLong) :: TM_helix_atom_exp
    type(s_selatoms)       :: TM_helix_atom
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

    call dealloc_selatoms(option%TM_helix_atom)

    return

  end subroutine dealloc_option

end module ta_option_str_mod
