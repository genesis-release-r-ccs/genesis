!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pd_draw_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pd_draw_mod

  use pd_option_str_mod
  use input_str_mod
  use output_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: draw

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine draw(input, option, output)

    ! formal arguments
    type(s_input),           intent(inout) :: input
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb_in
    integer                  :: i, j, k
    integer                  :: natoms, ilines, nlines, natm2, num_modes
    integer                  :: vec_in, vmd_out, pml_out, pdb_out
    real(wp)                 :: tmp, x, y, z, dx, dy, dz
    real(wp), allocatable    :: vector(:,:)

    if (option%check_only) &
      return


    ! setup coordintes
    !
    call input_pdb(input%pdb_avefile, pdb_in)

    natoms = pdb_in%num_atoms

    allocate(vector(3,natoms))


    ! read VEC file
    !
    call open_file(vec_in, input%vecfile, IOFileInput)

    ilines = 0
    do
      read(vec_in, *, end=999) tmp
      ilines = ilines + 1
    end do

999 nlines = ilines

    ! error check
    !
    natm2 = nlines/9
    if (natm2 /= natoms*natoms) &
      call error_msg('Draw> Number of atoms does not match between PDB and VEC files')

    num_modes = 3*natoms

    rewind(vec_in)

    do i = 1, num_modes
      if (i > option%mode_no) exit
      do j = 1, natoms
        do k = 1, 3
          read(vec_in,*) tmp
          if (i == option%mode_no) then
            vector(k,j) = tmp
          end if
        end do
      end do
    end do

    call close_file(vec_in)

    if (option%arrow_reverse) then
      do i = 1, natoms
        vector(1:3,i) = - vector(1:3,i)
      end do
    end if


    ! output VMD style
    !
    if (output%vmdfile /= '') then

      call open_file(vmd_out, output%vmdfile, IOFileOutputNew)

      write(vmd_out,'(a)') 'draw delete all'
      write(vmd_out,'(a,a)') 'draw color ', option%vector_color_vmd

      do i = 1, natoms

        x  = pdb_in%atom_coord(1,i)
        y  = pdb_in%atom_coord(2,i)
        z  = pdb_in%atom_coord(3,i)

        dx = x + vector(1,i)*option%expand_vector
        dy = y + vector(2,i)*option%expand_vector
        dz = z + vector(3,i)*option%expand_vector

        write(vmd_out,'(a,3f10.3,a,3f10.3,a,f10.3)') &
          'draw cylinder {',x,y,z,'} {',dx,dy,dz,'} radius', option%cylinder_radius

        x  = x + vector(1,i)*option%arrow_length*option%expand_vector
        y  = y + vector(2,i)*option%arrow_length*option%expand_vector
        z  = z + vector(3,i)*option%arrow_length*option%expand_vector
      
        write(vmd_out,'(a,3f10.3,a,3f10.3,a,f10.3)') &
          'draw cone     {',dx,dy,dz,'} {',x,y,z,'} radius', option%cone_radius

      end do

      call close_file(vmd_out)

      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') 'Draw> VMD visualization state file (.vmd) was generated'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  Usage:'
      write(MsgOut,'(A)') '    $ vmd input.pdb -e output.vmd'
      write(MsgOut,'(A)') '  or '
      write(MsgOut,'(A)') '    use "play" command in the VMD command window'
      write(MsgOut,'(A)') '  or '
      write(MsgOut,'(A)') '    File> Load Visualization State> Select and load .vmd file'
      write(MsgOut,'(A)') ''

    end if


    ! output PyMol style
    !
    if (output%pmlfile /= '') then

      call open_file(pml_out, output%pmlfile, IOFileOutputNew)

      write(pml_out,'(a)') 'from pymol.cgo import *'
      write(pml_out,'(a)') 'from pymol import cmd'
      write(pml_out,'(a)') ''
      write(pml_out,'(a)')   '# load INPUT PDB file'
      write(pml_out,'(a,a)') 'load ', trim(input%pdb_avefile)
      write(pml_out,'(a)') ''
      write(pml_out,'(a)') 'obj = [ \'

      do i = 1, natoms
        x  = pdb_in%atom_coord(1,i)
        y  = pdb_in%atom_coord(2,i)
        z  = pdb_in%atom_coord(3,i)

        dx = x + vector(1,i)*option%expand_vector
        dy = y + vector(2,i)*option%expand_vector
        dz = z + vector(3,i)*option%expand_vector

        write(pml_out,'(a,$)') 'CYLINDER, '
        write(pml_out,'(f10.3,a1,$)') x,','
        write(pml_out,'(f10.3,a1,$)') y,','
        write(pml_out,'(f10.3,a1,$)') z,','
        write(pml_out,'(f10.3,a1,$)') dx,','
        write(pml_out,'(f10.3,a1,$)') dy,','
        write(pml_out,'(f10.3,a1,$)') dz,','
        write(pml_out,'(f10.3,a8,$)') option%cylinder_radius,',       '
        write(pml_out,'(a1,a,a1,$)') ' ', trim(option%vector_color_pml),','
        write(pml_out,'(a1,a,a)')   ' ', trim(option%vector_color_pml),', \'

        x  = x + vector(1,i)*option%arrow_length*option%expand_vector
        y  = y + vector(2,i)*option%arrow_length*option%expand_vector
        z  = z + vector(3,i)*option%arrow_length*option%expand_vector

        write(pml_out,'(a,$)') 'CONE,     '
        write(pml_out,'(f10.3,a1,$)') dx,','
        write(pml_out,'(f10.3,a1,$)') dy,','
        write(pml_out,'(f10.3,a1,$)') dz,','
        write(pml_out,'(f10.3,a1,$)') x,','
        write(pml_out,'(f10.3,a1,$)') y,','
        write(pml_out,'(f10.3,a1,$)') z,','
        write(pml_out,'(f10.3,a8,$)') option%cone_radius,', 0.000,'
        write(pml_out,'(a1,a,a1,$)') ' ', trim(option%vector_color_pml),','
        write(pml_out,'(a1,a,a1,$)') ' ', trim(option%vector_color_pml),','
        write(pml_out,'(a)') ' 1.0, 1.0, \'


      end do
      write(pml_out,'(a)') ']'

      write(pml_out,'(a)') 'cmd.load_cgo(obj, "vectors")'

      call close_file(pml_out)

      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') 'Draw> PyMol script file (.pml) was generated'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  Usage:'
      write(MsgOut,'(A)') '    $ pymol output.pml'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') ''

    end if


    deallocate(vector)

    return

  end subroutine draw

end module pd_draw_mod
