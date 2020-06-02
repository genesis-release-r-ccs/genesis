!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_pdb2hm_mod
!> @brief   create huge molecule from PDB data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_pdb2hm_mod

  use pr_huge_molecule_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: hm_create_pdb_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_pdb_coord
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_pdb_coord(pdb_in, is_ref)

    ! formal arguments
    integer,                 intent(in) :: pdb_in
    logical,                 intent(in) :: is_ref

    ! local variables
    real(wip)                :: atom_coord(3)
    integer                  :: iatom
    character(80)            :: line


    if (.not. is_ref) then
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Pdb_Coord> '
    else
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Pdb_Coord (REF)> '
    end if

    iatom = 0

    do while(.true.)

      read(pdb_in, fmt='(a80)', err=900,end=100) line

      if (line(1:4) /= 'ATOM' .and. line(1:4) /= 'HETA') &
        cycle

      read(line, fmt='(30x,3f8.3)', err=900)  &
           atom_coord(1),  &
           atom_coord(2),  &
           atom_coord(3)

      iatom = iatom + 1

      if (.not. is_ref) then
        call hm_set_atom_coord   (1, iatom, atom_coord(1))
        call hm_set_atom_coord   (2, iatom, atom_coord(2))
        call hm_set_atom_coord   (3, iatom, atom_coord(3))
      else
        call hm_set_atom_refcoord(1, iatom, atom_coord(1))
        call hm_set_atom_refcoord(2, iatom, atom_coord(2))
        call hm_set_atom_refcoord(3, iatom, atom_coord(3))
      end if

    end do

100 if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Hm_Create_Pdb_Coord> read ERROR. ')

  end subroutine hm_create_pdb_coord

end module pr_pdb2hm_mod
