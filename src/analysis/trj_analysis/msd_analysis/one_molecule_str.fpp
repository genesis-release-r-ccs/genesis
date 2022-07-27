!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   one_molecule_str_mod
!> @brief   structure of single molecule
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module one_molecule_str_mod

  use molecules_str_mod

  implicit none

  private

  type, public, extends(s_molecule) :: s_one_molecule

    ! global atom indexes
    integer, dimension(:), allocatable :: atom_idx

    ! global atom index to local tom index mapping
    integer, dimension(:), allocatable :: atom_idx2loc


    ! number of substructures in molecule
    ! usually 1
    integer                            :: num_sub = 1

    ! local substructure index of each atom
    ! all values will usually be one,
    ! unless the molecule is constructed from disconnected substructures
    integer, dimension(:), allocatable :: sub_idx

  end type s_one_molecule

end module one_molecule_str_mod
