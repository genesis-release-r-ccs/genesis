!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ca_option_str_mod
!> @brief   structure of option information
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ca_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                :: check_only
    integer                :: analysis_for
    integer                :: pbcc_mode
    integer                :: output_coord
    logical                :: output_atomno
    character(MaxLineLong) :: analysis_atom_exp
    type(s_selatoms)       :: analysis_atom

  end type s_option

  ! parameters
  integer,      public, parameter :: AnalysisForMolecule = 1
  integer,      public, parameter :: AnalysisForAll      = 2

  integer,      public, parameter :: OutputCoordX        = 1
  integer,      public, parameter :: OutputCoordXY       = 2
  integer,      public, parameter :: OutputCoordXYZ      = 3

  character(*), public, parameter :: AnalysisForTypes(2) = &
                                                        (/'MOLECULE',&
                                                          'ALL     '/)

  character(*), public, parameter :: OutputCoordTypes(3) = &
                                                        (/'X  ',&
                                                          'XY ',&
                                                          'XYZ'/)

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

    return

  end subroutine dealloc_option

end module ca_option_str_mod
