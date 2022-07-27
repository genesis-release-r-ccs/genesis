!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   kc_option_str_mod
!> @brief   structure of option information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kc_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    type(s_selatoms)              :: analysis_atom
    integer                       :: num_clusters 
    integer                       :: max_iteration
    real(wp)                      :: stop_threshold
    integer                       :: num_iterations
    integer                       :: iseed
    type(s_selatoms)              :: trjout_atom
    integer                       :: trjout_format
    integer                       :: trjout_type

  end type s_option

  ! subroutines
  public  :: dealloc_option

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


    call dealloc_selatoms(option%trjout_atom)
    call dealloc_selatoms(option%analysis_atom)

    return

  end subroutine dealloc_option

end module kc_option_str_mod
