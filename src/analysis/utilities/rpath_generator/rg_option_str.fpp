!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    ! RPATH section
    integer                  :: nreplica
    integer                  :: cv_atom = 1
    integer                  :: iseed = 777
    integer                  :: iter_reparam = 1

    ! selection variables
    integer                       :: num_atoms
    type(s_selatoms), allocatable :: selatoms(:)
    character(MaxLine), allocatable :: groups(:)

    ! others
    logical                  :: is_cartesian
    integer                  :: fitting_atom = 1

  end type s_option

end module rg_option_str_mod
