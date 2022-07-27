!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_option_str_mod
!> @brief   structure of option information
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                :: check_only
                          
    integer                :: calc_mode
    character(MaxLineLong) :: analysis_atom_exp_1
    type(s_selatoms)       :: analysis_atom_1

    character(MaxLineLong) :: analysis_atom_exp_2
    type(s_selatoms)       :: analysis_atom_2
                          
    integer                :: matrix_shape

  end type s_option

  ! parameters

  ! matrix mode
  integer,      public, parameter :: MatrixShapeFull      = 1
  integer,      public, parameter :: MatrixShapeHalf      = 2

  character(*), public, parameter :: MatrixShapeTypes(2)  = &
                                                         (/'FULL','HALF'/)

  ! calculation mode
  integer,      public, parameter :: intramolecule        = 1
  integer,      public, parameter :: intermolecule        = 2

  character(*), public, parameter :: calculation_mode(2)  = &
                                                         (/'INTRA', 'INTER'/)

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


    call dealloc_selatoms(option%analysis_atom_1)
    call dealloc_selatoms(option%analysis_atom_2)

    return

  end subroutine dealloc_option

end module da_option_str_mod
