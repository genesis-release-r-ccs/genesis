!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   di_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module di_analyze_mod

  use di_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use string_mod
 
  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: output_temporary_pdb
  private :: write_temporary_pdb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: iatom, out_out, tmp_out, idx
    character(MaxLine)       :: command
    character(10)            :: sid

    if (option%check_only) &
      return

    ! check temporary PDB existence
    !
    call open_file (tmp_out, option%temporary_pdb, IOFileOutputNew)
    call close_file(tmp_out)


    ! open output file
    !
    if (output%outfile /= '') then
      call open_file(out_out, output%outfile, IOFileOutputNew)
      write(out_out,'(A)') '# ANALYSIS OF PROTEIN SECONDARY STRUCTURES BY DSSP_INTERFACE'
      write(out_out,'(A)') ''
      call close_file(out_out)
    else
      call error_msg('Analyze> outfile is not specified.')
    end if


    ! analysis loop
    !
    nstru = 0
    num_trjfiles = size(trj_list%md_steps)

    do ifile = 1, num_trjfiles

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        ! read trajectory
        !
        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(ifile)) == 0) then

          nstru = nstru+1
          write(MsgOut,*) '      number of structures = ', nstru
          write(MsgOut,*) ''

          ! write snapshot index
          !
          write (sid,'(i10)') nstru
          command = 'echo "SNAPSHOT' // sid // '" >> ' // trim(output%outfile)
          call system(command)

          ! execute DSSP command
          !
          call output_temporary_pdb(option, molecule, trajectory%coord)

          command = trim(option%dssp_exec)     // '    ' // &
                    trim(option%temporary_pdb) // ' >> ' // &
                    trim(output%outfile)
          call system(command)

          ! insert blank line
          !
          command = 'echo " " >> ' // trim(output%outfile)
          call system(command)

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_output_temporary_pdb 
  !> @brief        a driver subroutine for writing PDB file
  !! @authors      TM
  !! @param[in]    molecule     : structure of molecule
  !! @param[in]    coord        : coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_temporary_pdb(option, molecule, coord)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)

    ! local variables
    integer                  :: file


    ! open PDB file
    !
    call open_file(file, option%temporary_pdb, IOFileOutputReplace)

    ! write coordinate data from MD
    !
    call write_temporary_pdb(file, option, molecule, coord)

    ! close PDB file
    !
    call close_file(file)

    return

  end subroutine output_temporary_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_temporary_pdb
  !> @brief        write temporary PDB file
  !! @authors      CK, TM
  !! @param[in]    file     : unit number of PDB file
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates
  !! @note         copied from at_output.fpp
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_temporary_pdb(file, option, molecule, coord)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)

    ! local variables
    integer                  :: i, j, len, iatom
    character(80)            :: fmt_a, fmt_t
    character(6)             :: crec
    character(4)             :: catm, cres, cstr, cseg
    logical                  :: use_cid
    logical                  :: res_col6, res_col5
    logical                  :: atom_col7

    res_col6  = .false.
    res_col5  = .false.
    atom_col7 = .false.

    if (molecule%num_residues >= 100000) then
      res_col6 = .true.
    else if (molecule%num_residues >= 10000) then
      res_col5 = .true.
    end if

    if (molecule%num_atoms >= 1000000) then
      atom_col7 = .true.
    end if

    if (atom_col7) then
      if (res_col6) then
        use_cid = .false.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,i6,48x,a1)'
      else if (res_col5) then
        use_cid = .false.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,a1,i4,49x,a1)'
      end if
    else
      if (res_col6) then
        use_cid = .false.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,i6,48x,a1)'
      else if (res_col5) then
        use_cid = .false.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,a1,i4,49x,a1)'
      end if
    end if

    do iatom = 1, size(option%analysis_atom%idx)

      i = option%analysis_atom%idx(iatom)

      crec = 'ATOM  '

      read(molecule%atom_name(i), *) cstr
      len = len_trim(cstr)
      if (len < 4) then
        write(catm, fmt='(1x,a3)') cstr
      else
        catm = molecule%atom_name(i)
      end if

      read(molecule%residue_name(i),*) cstr
      len = len_trim(cstr)
      if (len == 2) then
        write(cres, fmt='(1x,a2)') cstr
      else if (len == 1) then
        write(cres, fmt='(2x,a1)') cstr
      else
        cres = cstr
      end if

      cseg = molecule%segment_name(i)

      if (use_cid) then
        write(file, fmt=fmt_a)              &
              crec,                         &
              iatom,                        &
              catm,                         &
              cres,                         &
              molecule%chain_id(i),         &
              molecule%residue_no(i),       &
              (coord(j,i), j=1,3),          &
              molecule%atom_occupancy(i),   &
              molecule%atom_temp_factor(i), &
              cseg
      else
        write(file, fmt=fmt_a)              &
              crec,                         &
              iatom,                        &
              catm,                         &
              cres,                         &
              molecule%residue_no(i),       &
              (coord(j,i), j=1,3),          &
              molecule%atom_occupancy(i),   &
              molecule%atom_temp_factor(i), &
              cseg
      end if

    end do

    write(file, fmt='(a6,69x,a1)') 'END   ', ' '

    return

  end subroutine write_temporary_pdb

end module di_analyze_mod
