!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_option_str_mod
!> @brief   structure of option information
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_option_str_mod

  use select_atoms_str_mod
  use molecules_str_mod
  use select_molecules_mod

  implicit none
  private

  type :: s_axes
    integer, dimension(:), allocatable :: i
  end type s_axes

  ! structures
  type, public :: s_option

    logical                                         :: check_only
    integer                                         :: delta
    logical                                         :: oversample
    type(s_selmols), dimension(:), allocatable :: analysis_mols
    type(s_one_molecule), dimension(:), allocatable :: all_mols
    type(s_axes), dimension(:), allocatable         :: axes

  end type s_option

  public :: dealloc_option

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      DS
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option

    if (allocated(option%analysis_mols)) deallocate(option%analysis_mols)

    return

  end subroutine dealloc_option

end module ma_option_str_mod
