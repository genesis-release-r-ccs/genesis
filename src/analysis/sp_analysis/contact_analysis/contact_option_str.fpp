!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   contact_option_str_mod
!> @brief   structure of option information
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module contact_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_contact_option

    integer    :: mode        = 0
    real(wp)   :: ana_range   = 0.0_wp
    integer    :: recenter    = 0

  end type s_contact_option

  ! parameters
  integer,      public, parameter :: ModeMindist = 1
  integer,      public, parameter :: ModeNumber  = 2

  character(*), public, parameter :: ModeTypes(2) = (/'MINDIST', &
                                                      'NUMBER '/)

  ! subroutines
  public :: dealloc_contact_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_contact_option
  !> @brief        deallocate option
  !! @authors      IY
  !! @param[inout] contact_option : Contact option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_contact_option(contact_option)

    ! formal argments
    type(s_contact_option),  intent(inout) :: contact_option

    return

  end subroutine dealloc_contact_option

end module contact_option_str_mod
