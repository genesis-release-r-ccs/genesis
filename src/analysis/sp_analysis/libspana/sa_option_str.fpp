!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_option_str_mod
!> @brief   structure of option information
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical             :: wrap
    integer             :: determine_box

    character(MaxLine), allocatable  :: analysis_atom_exps(:)
    type(s_selatoms),   allocatable  :: analysis_atoms(:)

    ! parameters for rdf, density
    real(wp)            :: buffer

  end type s_option

  ! structures
  type, public :: s_parray
    integer, pointer :: idx(:)
  end type s_parray

  ! parameters
  integer,      public, parameter :: DetermineBoxTrajectory = 1
  integer,      public, parameter :: DetermineBoxManual     = 2
  integer,      public, parameter :: DetermineBoxMax        = 3

  character(*), public, parameter :: DetermineBoxTypes(3) = (/'TRAJECTORY', &
                                                              'MANUAL    ', &
                                                              'MAX       '/)

  ! subroutines
  public :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      IY
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option

    !call dealloc_selatoms(option%analysis_atom)

    return

  end subroutine dealloc_option

end module sa_option_str_mod
