!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_pdb2cache_mod
!> @brief   create cached molecule from AMBER data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_pdb2cache_mod

  use pr_cached_molecule_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: cm_create_cache_pdb_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_pdb_coord
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_pdb_coord(pdb_in, is_ref)

    ! formal arguments
    integer,                 intent(in) :: pdb_in
    logical,                 intent(in) :: is_ref

    ! local variables
    integer                  :: iatom, idata, file_out, num
    character(80)            :: line, name

    real(wp),    allocatable :: atom_coord(:,:)


    if (.not. is_ref) then
      write(MsgOut,'(a)') 'Cm_Create_Cache_Pdb_Coord> '
      name = 'atom_coord'
    else
      write(MsgOut,'(a)') 'Cm_Create_Cache_Pdb_Coord (REF)> '
      name = 'atom_refcoord'
    end if

    allocate(atom_coord(3, CMPageSize))

    ! read pdb file
    iatom = 1
    idata = 1

    do while(.true.)

      read(pdb_in, fmt='(a80)', end=200, err=900) line

      if (line(1:4) /= 'ATOM' .and. line(1:4) /= 'HETA') &
        cycle

      read(line, fmt='(30x,3f8.3)', err=900)  &
           atom_coord(1,idata),  &
           atom_coord(2,idata),  &
           atom_coord(3,idata)

      idata = idata + 1

      if (idata > CMPageSize .or. iatom + 1> cm_num_atoms) then

        if (idata - 1 == CMPageSize) then
          num = iatom - CMPageSize + 1
        else
          num = iatom / CMPageSize * CMPageSize + 1
        end if

        ! atom_coord
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, name), &
                              IOFileOutputNew)
        write(file_out) atom_coord(1:3,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

      iatom = iatom + 1

    end do
200 continue

    deallocate(atom_coord)

    write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Cm_Create_Cache_Pdb_Coord> read ERROR. ')

  end subroutine cm_create_cache_pdb_coord

end module pr_pdb2cache_mod
