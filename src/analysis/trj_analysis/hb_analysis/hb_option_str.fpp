!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hb_option_str_mod
!> @brief   structure of option information
!! @authors Daisuke Matsuoka (DM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hb_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    integer                       :: output_type
    real(wp)                      :: hb_distance
    real(wp)                      :: dha_angle
    real(wp)                      :: hda_angle
    character(MaxLineLong)        :: analysis_atom_exp
    type(s_selatoms)              :: analysis_atom
    character(MaxLineLong)        :: target_atom_exp
    type(s_selatoms)              :: target_atom
    integer                       :: boundary_type
    character(LEN=6), allocatable :: solvent_list(:)

  end type s_option

  ! parameters for result output
  integer, public, parameter :: HBOutputModeCountAtom = 1
  integer, public, parameter :: HBOutputModeCountSnap = 2
  integer, public, parameter :: HBOutputModeLifetime  = 3

  character(len=10), public, parameter :: HBOutputType(3) = (/'COUNT_atom', &
                                                              'COUNT_snap', &
                                                              'LIFETIME  '/)

  ! parameters for boundary type
  integer,      public, parameter :: BoundaryTypeNOBC     = 1
  integer,      public, parameter :: BoundaryTypePBC      = 2

  character(*), public, parameter :: BoundaryTypeTypes(2) = (/'NOBC  ', &
                                                              'PBC   '/)


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


    call dealloc_selatoms(option%analysis_atom)
    call dealloc_selatoms(option%target_atom)

    if (allocated(option%solvent_list))  &
      deallocate(option%solvent_list)

    return

  end subroutine dealloc_option

end module hb_option_str_mod
