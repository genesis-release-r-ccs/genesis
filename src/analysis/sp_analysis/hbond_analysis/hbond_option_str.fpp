!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hbond_option_str_mod
!> @brief   structure of option information
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hbond_option_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_hbond_option

    integer                       :: recenter
    integer                       :: output_type
    integer                       :: analysis_atom
    integer                       :: target_atom
    real(wp)                      :: hb_distance
    real(wp)                      :: dha_angle
    real(wp)                      :: hda_angle
    character(len=6), allocatable :: solvent_list(:)

  end type s_hbond_option

  ! parameters
  integer, public, parameter :: HBOutputModeCountAtom = 1
  integer, public, parameter :: HBOutputModeCountSnap = 2

  character(len=10), public, parameter :: HBOutputType(2) = (/'Count_Atom', &
                                                              'Count_Snap'/)
  ! subroutines
  public :: dealloc_hbond_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_hbond_option
  !> @brief        deallocate hbond_option
  !! @authors      TM
  !! @param[inout] hbond_option : hbond_option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_hbond_option(hbond_option)

    ! formal argments
    type(s_hbond_option), intent(inout) :: hbond_option


    if (allocated(hbond_option%solvent_list))  &
      deallocate(hbond_option%solvent_list)

    return

  end subroutine dealloc_hbond_option

end module hbond_option_str_mod
