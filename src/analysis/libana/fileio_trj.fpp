!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_trj_mod
!> @brief   read/write trajectory file
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK), Donatas Surblys(DS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module fileio_trj_mod

  use trajectory_str_mod
  use molecules_str_mod
  use select_atoms_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: open_trj
  public  :: close_trj
  public  :: read_trj
  public  :: write_trj
  public  :: seek_trj

  public  :: get_num_steps_trj

  private :: open_pdb_trj
  private :: close_pdb_trj
  private :: read_pdb_trj
  private :: write_pdb_trj_sel
  private :: write_pdb_trj

  private :: open_amber_trj
  private :: close_amber_trj
  private :: read_amber_trj
  private :: write_amber_trj_sel
  private :: write_amber_trj
  private :: skip_amber_trj_header
  private :: write_amber_trj_header

  private :: open_dcd_trj
  private :: close_dcd_trj
  private :: read_dcd_trj
  private :: read_dcd_trj_header
  private :: write_dcd_trj_header
  private :: write_dcd_trj_sel
  private :: write_dcd_trj
  private :: seek_dcd_trj

  private :: open_gromacs_trj
  private :: close_gromacs_trj
  private :: read_gromacs_trj
  private :: write_gromacs_trj_sel
  private :: write_gromacs_trj

  private :: to_pbcbox
  private :: to_symmat
  private :: to_deg

  ! parameters
  integer,    private             :: GRXMagic = 1993

  character(4), allocatable, save :: nama(:)
  real(sp),     allocatable, save :: crd4f(:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_trj
  !> @brief        Open trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    filename   : trajectory file name
  !! @param[in]    trj_format : trajectory file format
  !! @param[in]    trj_type   : trajectory file type
  !! @param[in]    in_out     : information for input or output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_trj(file, filename, trj_format, trj_type, in_out)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: trj_format
    integer,                 intent(in)    :: trj_type
    integer,                 intent(in)    :: in_out


    file%name       = filename
    file%trj_format = trj_format
    file%trj_type   = trj_type
    file%in_out     = in_out

    select case (file%trj_format)

    case (TrjFormatPDB)

      call open_pdb_trj(file)

    case (TrjFormatAmber)

      call open_amber_trj(file)

    case (TrjFormatDCD)

      call open_dcd_trj(file)

    case (TrjFormatGromacs)

      call open_gromacs_trj(file)

    end select

    file%step_no = 1

    return

  end subroutine open_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_trj
  !> @brief        Close trajectory file
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file


    select case (file%trj_format)

    case (TrjFormatPDB)

      call close_pdb_trj(file)

    case (TrjFormatAmber)

      call close_amber_trj(file)

    case (TrjFormatDCD)

      call close_dcd_trj(file)
      call update_dcd_trj_header(file)

    case (TrjFormatGromacs)

      call close_gromacs_trj(file)

    end select

    return

  end subroutine close_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_trj
  !> @brief        Read trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_trj(file, trajectory)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(inout) :: trajectory


    select case (file%trj_format)

    case (TrjFormatPDB)

      call read_pdb_trj(file, trajectory)

    case (TrjFormatAmber)

      call skip_amber_trj_header(file)
      call read_amber_trj(file, trajectory)

    case (TrjFormatDCD)

      call read_dcd_trj_header(file)
      call read_dcd_trj(file, trajectory)

    case (TrjFormatGromacs)

       call read_gromacs_trj(file, trajectory)

    end select

    file%step_no = file%step_no + 1

    return

  end subroutine read_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trj
  !> @brief        Write trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selatoms   : atom selection information
  !! @param[in]    molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trj(file, trajectory, selatoms, molecule)

    ! formal arguments
    type(s_trj_file),             intent(inout) :: file
    type(s_trajectory),           intent(in)    :: trajectory
    type(s_selatoms),   optional, intent(in)    :: selatoms
    type(s_molecule),   optional, intent(in)    :: molecule


    select case (file%trj_format)

    case (TrjFormatPDB)

      if (.not. present(molecule)) then
        call error_msg('Write_Trj> Molecular information is necessary '//&
                       'for the output in the PDB format')
      end if

      if (present(selatoms)) then
        call write_pdb_trj_sel(file, trajectory, selatoms, molecule)
      else
        call write_pdb_trj(file, trajectory, molecule)
      end if

    case (TrjFormatAmber)

      call write_amber_trj_header(file)

      if (present(selatoms)) then
        call write_amber_trj_sel(file, trajectory, selatoms)
      else
        call write_amber_trj(file, trajectory)
      end if

    case (TrjFormatDCD)

      if (present(selatoms)) then
        call write_dcd_trj_header(file, size(selatoms%idx))
        call write_dcd_trj_sel(file, trajectory, selatoms)
      else
        call write_dcd_trj_header(file, size(trajectory%coord(1,:)))
        call write_dcd_trj(file, trajectory)
      end if

    case (TrjFormatGromacs)

      if (present(selatoms)) then
        call write_gromacs_trj_sel(file, trajectory, selatoms)
      else
        call write_gromacs_trj(file, trajectory)
      end if

    case (TrjFormatCharmmRst)

      call error_msg('Write_Trj> Charmm restart is not supported.')

    case (TrjFormatNamdRst)

      call error_msg('Write_Trj> NAMD restart is not supported.')

    end select

    file%step_no = file%step_no + 1

    return

  end subroutine write_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    seek_trj
  !> @brief        Seek trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    step_no    : step number 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine seek_trj(file, trajectory, step_no)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(inout) :: trajectory
    integer,                 intent(in)    :: step_no


    select case (file%trj_format)

    case (TrjFormatPDB)

      call error_msg('Seek_Trj> Seek is not supported. [PDB]')

    case (TrjFormatAmber)

      call error_msg('Seek_Trj> Seek is not supported. [AMBER]')

    case (TrjFormatDCD)

      call read_dcd_trj_header(file)
      call seek_dcd_trj(file, trajectory, step_no)

    case (TrjFormatGromacs)

      call error_msg('Seek_Trj> Seek is not supported. [GROMACS]')

    end select

    return

  end subroutine seek_trj

  !*********************************************************************
  ! pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_pdb_trj
  !> @brief        Open pdb trajectory file
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_pdb_trj(file)

    ! formal arugments
    type(s_trj_file),        intent(inout) :: file
    

    call open_file(file%unit_no, file%name, file%in_out)

    return 

  end subroutine open_pdb_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_pdb_trj
  !> @brief        Close pdb trajectory file
  !! @authors      NT
  !! @param[in]    file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_pdb_trj(file)

    ! formal arugments
    type(s_trj_file),        intent(in)    :: file


    call close_file(file%unit_no)

    return

  end subroutine close_pdb_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pdb_trj
  !> @brief        Read pdb trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pdb_trj(file, trajectory)
    
    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(inout) :: trajectory

    ! local variables
    character(80)               :: line
    integer                     :: i, j, natm, unit

    real(wp),           pointer :: crd(:,:)


    unit =  file%unit_no
    crd  => trajectory%coord
    natm =  size(crd(1,:))

    i = 0

    do while (.true.)

      read(unit, '(a80)', end=900, err=900) line

      if ((line(1:4) == 'ATOM') .or. &
          (line(1:4) == 'HETA')) then

        i = i + 1 

        ! atom coordinates
        read(line, '(30x,3f8.3)', err=900) (crd(j,i), j=1,3)
        if (i >= natm) &
          exit

      else if (line(1:4) == 'CRYS') then

        ! crystal
        read(line, '(6x,3f9.3)', err=900) (trajectory%pbc_box(j,j), j=1,3)

      end if

    end do

    return

900 call error_msg('Read_Pdb_Trj> read error')

  end subroutine read_pdb_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pdb_trj_sel
  !> @brief        Write pdb trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selatoms   : atom selection information
  !! @param[in]    molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pdb_trj_sel(file, trajectory, selatoms, molecule)

    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(in)    :: trajectory
    type(s_selatoms),   target, intent(in)    :: selatoms
    type(s_molecule),   target, intent(in)    :: molecule

    ! local varialbles
    integer                     :: i, j, iatm, nsel, unit, molec, molep

    real(wp),           pointer :: crd(:,:)
    integer,            pointer :: numr(:), numm(:), aidx(:)
    character(6),       pointer :: namr(:)
    character(5)                :: nchar_atom
!    character(4),   allocatable :: nama(:)
    character(4),       pointer :: seg(:)
    character(1),       pointer :: chn(:)


    unit =  file%unit_no
    namr => molecule%residue_name
    numr => molecule%residue_no
    numm => molecule%molecule_no
    crd  => trajectory%coord
    seg  => molecule%segment_name
    chn  => molecule%chain_id

!    allocate(nama(molecule%num_atoms))
    if (.not. allocated(nama)) then
      allocate(nama(molecule%num_atoms))
    else
      if (size(nama(:)) /= molecule%num_atoms) then
        deallocate(nama)
        allocate(nama(molecule%num_atoms))
      end if
    end if

    do i = 1, molecule%num_atoms
      if (molecule%atom_name(i)(4:4) == ' ') then
        nama(i) = ' '//trim(molecule%atom_name(i))
      else
        nama(i) = trim(molecule%atom_name(i))
      end if
    end do

    if (file%trj_type == TrjTypeCoorBox) then
      write(unit, '(a6,3(f9.3),a37)') 'CRYST1', &
           (trajectory%pbc_box(j,j), j=1,3), &
           '  90.00  90.00  90.00 P 1           1'
    end if

    write(unit, '(a5,5x,i4)') 'MODEL', file%step_no


    aidx  => selatoms%idx
    nsel  = size(selatoms%idx)

    molep = -1

    do i = 1, nsel
      iatm = aidx(i)

      molec = numm(iatm)
      if (molep /= -1 .and. molep /= molec) then
        write(unit,'(a3)') 'TER'
      end if
      molep = molec
      if (iatm <= 99999) then
        write(nchar_atom,'(i5)') iatm
      else
        write(nchar_atom,'(a5)') '*****'
      end if

      if (numr(iatm) <= 9999) then

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,a1,i4,4x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(iatm), namr(iatm), chn(iatm), numr(iatm), (crd(j,iatm),j=1,3), &
                  0.00, 0.00, seg(iatm)

      else if (numr(iatm) > 999999) then 

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,i6,4x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(iatm), namr(iatm), numr(iatm), (crd(j,iatm),j=1,3), &
                  0.00, 0.00, seg(iatm)

      else

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(iatm), namr(iatm), chn(iatm), numr(iatm), (crd(j,iatm),j=1,3), &
                  0.00, 0.00, seg(iatm)

      end if
    end do

    write(unit,'(a3)') 'TER'
    write(unit,'(a6)') 'ENDMDL'

!    deallocate(nama)

    return

  end subroutine write_pdb_trj_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pdb_trj
  !> @brief        Write pdb trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pdb_trj(file, trajectory, molecule)

    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(in)    :: trajectory
    type(s_molecule),   target, intent(in)    :: molecule

    ! local varialbles
    integer                     :: i, j, natm, unit, molec, molep

    real(wp),           pointer :: crd(:,:)
    integer,            pointer :: numr(:), numm(:)
    character(6),       pointer :: namr(:)
    character(5)                :: nchar_atom
!    character(4),   allocatable :: nama(:)
    character(4),       pointer :: seg(:)
    character(1),       pointer :: chn(:)


    unit =  file%unit_no
    namr => molecule%residue_name
    numr => molecule%residue_no
    numm => molecule%molecule_no
    crd  => trajectory%coord
    seg  => molecule%segment_name
    chn  => molecule%chain_id

!    allocate(nama(molecule%num_atoms))
    if (.not. allocated(nama)) then
      allocate(nama(molecule%num_atoms))
    else
      if (size(nama(:)) /= molecule%num_atoms) then
        deallocate(nama)
        allocate(nama(molecule%num_atoms))
      end if
    end if

    do i = 1, molecule%num_atoms
      if (molecule%atom_name(i)(4:4) == ' ') then
        nama(i) = ' '//trim(molecule%atom_name(i))
      else
        nama(i) = trim(molecule%atom_name(i))
      end if
    end do

    if (file%trj_type == TrjTypeCoorBox) then
      write(unit, '(a6,3(f9.3),a37)') 'CRYST1', &
           (trajectory%pbc_box(j,j), j=1,3),    &
           '  90.00  90.00  90.00 P 1           1'
    end if

    write(unit, '(a5,5x,i4)') 'MODEL', file%step_no


    natm = size(crd(1,:))

    molep = -1

    do i = 1, natm

      molec = numm(i)
      if (molep /= -1 .and. molep /= molec) then
        write(unit,'(a3)') 'TER'
      end if
      molep = molec
      if (i <= 99999) then
        write(nchar_atom,'(i5)') i
      else
        write(nchar_atom,'(a5)') '*****'
      end if

      if (numr(i) <= 9999) then

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,a1,i4,4x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(i), namr(i), chn(i), numr(i), (crd(j,i),j=1,3), &
                  0.00, 0.00, seg(i)

      else if (numr(i) > 999999) then 

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,i6,4x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(i), namr(i), numr(i), (crd(j,i),j=1,3), &
                  0.00, 0.00, seg(i)

      else

        write(unit, '(a4,2x,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,2f6.2,6x,a4)') &
          'ATOM', nchar_atom, nama(i), namr(i), chn(i), numr(i), (crd(j,i),j=1,3), &
                  0.00, 0.00, seg(i)

      end if
    end do
    write(unit,'(a3)') 'TER'
    write(unit,'(a6)') 'ENDMDL'

!    deallocate(nama)

    return

  end subroutine write_pdb_trj

  !*********************************************************************
  ! amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_amber_trj
  !> @brief        Open AMBER trajectory file
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_amber_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file


    call open_file(file%unit_no, file%name, file%in_out)

    file%rw_header = .true.

    return 
    
  end subroutine open_amber_trj
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_amber_trj
  !> @brief        Close AMBER trajectory file
  !! @authors      NT
  !! @param[in]    file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_amber_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(in)    :: file


    call close_file(file%unit_no)

    return

  end subroutine close_amber_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_amber_trj
  !> @brief        Read AMBER trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_amber_trj(file, trajectory)

    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(inout) :: trajectory

    ! local variables
    integer                     :: i, j, natm, unit

    real(wp),           pointer :: crd(:,:)


    unit =  file%unit_no
    crd  => trajectory%coord
    natm =  size(crd(1,:))

    read(unit,*,err=900,end=900) ((crd(j,i), j=1,3), i=1,natm)

    if (file%trj_type == TrjTypeCoorBox) &
      read(unit,*,err=900,end=900) (trajectory%pbc_box(j,j), j=1,3)

    return

900 call error_msg('Read_Amber_Trj> read error ')

  end subroutine read_amber_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_amber_trj_sel
  !> @brief        Write AMBER trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selatoms   : selected atoms information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_amber_trj_sel(file, trajectory, selatoms)

    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(in)    :: trajectory
    type(s_selatoms),   target, intent(in)    :: selatoms

    ! local variables
    integer                     :: i, j, natm, unit

    real(wp),           pointer :: crd(:,:)
    integer,            pointer :: idx(:)


    unit =  file%unit_no
    crd  => trajectory%coord

    natm =  size(selatoms%idx)
    idx  => selatoms%idx

    write(unit,'(10(f8.3))',err=900) ((crd(j, idx(i)), j=1,3), i=1,natm)

    if (file%trj_type == TrjTypeCoorBox) &
      write(unit,'(3(f8.3))',err=900) (trajectory%pbc_box(j,j), j=1,3)

    return

900 call error_msg('Write_Amber_Trj_Sel> write error.')

  end subroutine write_amber_trj_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_amber_trj
  !> @brief        Write AMBER trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_amber_trj(file, trajectory)

    ! formal arguments
    type(s_trj_file),           intent(inout) :: file
    type(s_trajectory), target, intent(in)    :: trajectory

    ! local variables
    integer                     :: i, j, natm, unit

    real(wp),           pointer :: crd(:,:)


    unit =  file%unit_no
    crd  => trajectory%coord
    natm =  size(crd(1,:))

    write(unit,'(10(f8.3))',err=900) ((crd(j, i), j=1,3), i=1,natm)

    if (file%trj_type == TrjTypeCoorBox) &
      write(unit,'(3(f8.3))',err=900) (trajectory%pbc_box(j,j), j=1,3)

    return

900 call error_msg('Write_Amber_Trj> write error.')

  end subroutine write_amber_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    skip_amber_trj_header
  !> @brief        Skip AMBER trajectory file header
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine skip_amber_trj_header(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file

    ! local variables
    character(80)            :: line


    if (.not. file%rw_header) &
      return

    file%rw_header = .false.

    read(file%unit_no,'(a80)',err=900) line

    return 

900 call error_msg('Skip_Amber_Trj_Header> error.')

  end subroutine skip_amber_trj_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_amber_trj_header
  !> @brief        Write AMBER trajectory file header
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_amber_trj_header(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file


    if (.not. file%rw_header) &
      return

    file%rw_header = .false.

    write(file%unit_no,'(a)') ' CREATED BY GENESIS'

    return

  end subroutine write_amber_trj_header

  !*********************************************************************
  ! dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_dcd_trj
  !> @brief        Open DCD trajectory file
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_dcd_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file

    ! local variable 
    integer                  :: unit, endian
    integer(4)               :: i4


    ! check endianess
    !
    endian = IOFileNativeEndian

    if (file%in_out == IOFileInput) then

      unit = get_unit_no()
      open(unit, &
           file    = file%name, &
           status  = 'old',     &
           form    = 'unformatted',   &
           access  = 'stream',        & 
           convert = 'little_endian', &
           err     = 900)

      read(unit) i4
      if (i4 == 84) then
        endian = IOFileLittleEndian
      else
        endian = IOFileBigEndian
      end if

      close(unit)
      call free_unit_no(unit)

    end if

    ! open file
    !
    call open_binary_file(file%unit_no, file%name, file%in_out, endian)

    file%rw_header = .true.
    file%byte_swap = .false.

    return

900 call error_msg('Open_Dcd_Trj> I/O Error.'//trim(file%name))

    return

  end subroutine open_dcd_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_dcd_trj
  !> @brief        Close DCD trajectory file
  !! @authors      NT
  !! @param[in]    file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_dcd_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(in)    :: file


    call close_file(file%unit_no)

    return

  end subroutine close_dcd_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_dcd_trj_header
  !> @brief        Read DCD trajectory file header
  !! @authors      DS, NT
  !! @param[inout] file      : trajectory file handle
  !! @param[out]   num_steps : number of steps in file (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_dcd_trj_header(file, num_steps)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    integer, optional,       intent(out)   :: num_steps

    ! local variables
    integer                  :: i, unit, hdr_size, step_size
    integer                  :: step_size_coorbox, step_size_coor
    integer                  :: num_steps_coorbox, num_steps_coor
    integer                  :: num_steps_header
    integer(8)               :: file_size
    integer(4)               :: icntrl(20), ntitle, natom
    character(4)             :: hdr
    character(80)            :: title(10)

    integer                  :: ftell
    logical                  :: valid_coorbox, valid_coor


    if (.not. file%rw_header) then
      if (present(num_steps)) then
        call error_msg('Read_Dcd_Trj_Header> Can only determine number &
          &of steps before first read.')
      end if
      return
    end if
    file%rw_header = .false.

    unit = file%unit_no


    ! read header,cntrl
    !
    read(unit) hdr, icntrl(1:20)
    num_steps_header = icntrl(1)

    if (hdr /= 'CORD' .and. hdr /= 'VELD') &
      call error_msg('Read_Dcd_Trj_Header> Bad dcd file format.')
      

    ! read title
    !
    read(unit) ntitle,(title(i),i=1,ntitle)


    ! read natom
    !
    read(unit) natom


    ! check pbc box
    !
    ! STANDARD dcd file is not include box information ??
    !
    hdr_size = &
         4 + 4 + 20*4      + 4 + &  ! => read() hdr, icntrl
         4 + 4 + 80*ntitle + 4 + &  ! => read() ntitle, title(:)
         4 + 4             + 4      ! => read() natom

    ! step size without box information
    step_size_coor = (4 + 4 * natom + 4) * 3

    ! step size with box information
    ! box size is described by 6 doubles
    step_size_coorbox = step_size_coor + 4 + 8 * 6 + 4

#ifdef KCOMP
    inquire(unit,flen=file_size)
#else
    inquire(unit,size=file_size)
#endif

    num_steps_coor = (file_size - hdr_size) / step_size_coor
    if (mod(file_size - hdr_size, step_size_coor) == 0) then
      valid_coor = .true.
    else
      valid_coor = .false.
    end if

    num_steps_coorbox = (file_size - hdr_size) / step_size_coorbox
    if (mod(file_size - hdr_size, step_size_coorbox) == 0) then
      valid_coorbox = .true.
    else
      valid_coorbox = .false.
    end if

    ! in case of two positives, check if any step numbers match num_steps_header
    if (valid_coorbox .and. valid_coor .and. num_steps_header /= 0) then

      if (num_steps_coorbox == num_steps_header) then
        valid_coor = .false.
      else if (num_steps_coor == num_steps_header) then
        valid_coorbox = .false.
      end if

    end if

    ! got one positive, can reliably determine TrjType
    if (valid_coorbox .neqv. valid_coor) then

      if (valid_coorbox) then
        file%trj_type = TrjTypeCoorBox
      else if (valid_coor) then
        file%trj_type = TrjTypeCoor
      end if

!TM(181230)
!      if (main_rank) &
!      write(MsgOut, '("Read_Dcd_Trj_Header> Setting dcd file format to ",&
!        &A,".")') trim(TrjTypeTypes(file%trj_type))
!TM(181230)

    ! got two positives, cannot determine TrjType
    else if (valid_coorbox .and. valid_coor) then

      if (main_rank) &
      write(MsgOut, '("Read_Dcd_Trj_Header> Cannot determine dcd file format &
        &due to ambigious size. Using user supplied format: ",A,".")') &
        trim(TrjTypeTypes(file%trj_type))

    ! got no positives, cannot determine TrjType
    else

      if (main_rank) &
      write(MsgOut, '("Read_Dcd_Trj_Header> Cannot determine dcd file format &
        &due to invalid size. Using user supplied format: ",A,".")') &
        trim(TrjTypeTypes(file%trj_type))

    end if

    if (present(num_steps)) then

      select case(file%trj_type)

      case(TrjTypeCoorBox)

        num_steps = num_steps_coorbox

      case(TrjTypeCoor)

        num_steps = num_steps_coor

      case default

        call error_msg("Read_Dcd_Trj_Header> Invalid dcd file format.")

      end select

      if (num_steps_header /= 0 .and. num_steps /= num_steps_header) then
        write(ErrOut, '("Read_Dcd_Trj_Header> Warning: Step number in header ("&
          &,i0,") is different from detected step number (",i0,").")') &
          num_steps_header, num_steps
      end if

    end if

    return

  end subroutine read_dcd_trj_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dcd_trj_header
  !> @brief        Write DCD trajectory file header
  !! @authors      NT
  !! @param[inout] file      : trajectory file handle
  !! @param[in]    num_atoms : number of atom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dcd_trj_header(file, num_atoms)

    ! parameters
    character(4),            parameter     :: hdr = 'CORD'
    integer(4),              parameter     :: ntitle = 4

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    integer,                 intent(in)    :: num_atoms

    ! local variables
    integer                  :: i, unit
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date


    if (.not. file%rw_header) &
      return
    file%rw_header = .false.

    unit = file%unit_no


    ! write header,cntrl
    !
    icntrl(:)  = 0
    icntrl(2)  = file%step_no                 ! first time step     (NPRIV)
    icntrl(3)  = 1                            ! output period       (NSAVC)
    icntrl(4)  = 0                            ! number of time step (NSTEP)
    icntrl(10) = transfer(1.0_sp, icntrl(10)) ! length of time step (DELTA)
    if (file%trj_type == TrjTypeCoorBox) &
    icntrl(11) = 1                            ! flag for with unit-cell
    icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

    write(unit) hdr, icntrl(1:20)


    ! write title
    !
    call fdate (date)
    call getlog(name)

    title(1) = 'REMARKS CREATED BY GENESIS ANALYSIS TOOL'
    title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
    title(3) = ''
    title(4) = ''

    write(unit) ntitle, (title(i),i=1,ntitle)


    ! write natom
    !
    write(unit) num_atoms


    return

  end subroutine write_dcd_trj_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_dcd_trj_header
  !> @brief        Update DCD trajectory file header
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_dcd_trj_header(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file

    ! local variables
    integer                  :: unit


    if (file%in_out /= IOFileInput) then

      unit = get_unit_no()
      open(unit,                   &
           file   = file%name,     &
           status = 'old',         &
           form   = 'unformatted', &
           access = 'direct',      &
           recl   = 4)

      write(unit,rec=3) file%step_no - 1
      write(unit,rec=6) file%step_no - 1

      close(unit)
      call free_unit_no(unit)

    end if

    return

  end subroutine update_dcd_trj_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_dcd_trj
  !> @brief        Read DCD trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_dcd_trj(file, trajectory)
    
    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    real(dp)                 :: symmat(6), pbcbox(6)
    integer                  :: i, natm, unit

!    real(sp), allocatable    :: crd4f(:)


    unit = file%unit_no


    ! read periodic boundary condition
    !
    if (file%trj_type == TrjTypeCoorBox) then

      read(unit) symmat(1:6)

      call to_pbcbox(symmat, pbcbox)

      trajectory%pbc_box(1,1) = real(pbcbox(1),wp)
      trajectory%pbc_box(2,2) = real(pbcbox(2),wp)
      trajectory%pbc_box(3,3) = real(pbcbox(3),wp)

    end if


    ! read coordinates
    !

    natm = size(trajectory%coord(1,:))

!    allocate(crd4f(natm))
    if (.not. allocated(crd4f)) then
      allocate(crd4f(natm))
    else
      if (size(crd4f(:)) /= natm) then
        deallocate(crd4f)
        allocate(crd4f(natm))
      end if
    end if

    read(unit) crd4f(1:natm)
    trajectory%coord(1,1:natm) = real(crd4f(1:natm),wp)
    read(unit) crd4f(1:natm)
    trajectory%coord(2,1:natm) = real(crd4f(1:natm),wp)
    read(unit) crd4f(1:natm)
    trajectory%coord(3,1:natm) = real(crd4f(1:natm),wp)

!    deallocate(crd4f)

    return

  end subroutine read_dcd_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dcd_trj_sel
  !> @brief        Write DCD trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selatoms   : atom selection information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dcd_trj_sel(file, trajectory, selatoms)
    
    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_selatoms),        intent(in)    :: selatoms

    ! local variables
    real(dp)                 :: pbcbox(6), symmat(6)
    integer                  :: natm, unit


    unit = file%unit_no


    ! write periodic boundary condition
    !
    if (file%trj_type == TrjTypeCoorBox) then

      pbcbox(1) = real(trajectory%pbc_box(1,1),dp)
      pbcbox(2) = real(trajectory%pbc_box(2,2),dp)
      pbcbox(3) = real(trajectory%pbc_box(3,3),dp)
      pbcbox(4) = 90.0_dp
      pbcbox(5) = 90.0_dp
      pbcbox(6) = 90.0_dp

      call to_symmat(pbcbox, symmat)

      write(unit) symmat(1:6)

    end if


    ! write atom coordinates
    !
    natm = size(selatoms%idx)

    write(unit) real(trajectory%coord(1,selatoms%idx(1:natm)),sp)
    write(unit) real(trajectory%coord(2,selatoms%idx(1:natm)),sp)
    write(unit) real(trajectory%coord(3,selatoms%idx(1:natm)),sp)

    return

  end subroutine write_dcd_trj_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_dcd_trj
  !> @brief        Write DCD trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_dcd_trj(file, trajectory)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(in)    :: trajectory

    ! local variables
    real(dp)                 :: pbcbox(6), symmat(6)
    integer                  :: natm, unit


    unit = file%unit_no


    ! write periodic boundary condition
    !
    if (file%trj_type == TrjTypeCoorBox) then

      pbcbox(1) = real(trajectory%pbc_box(1,1),dp)
      pbcbox(2) = real(trajectory%pbc_box(2,2),dp)
      pbcbox(3) = real(trajectory%pbc_box(3,3),dp)
      pbcbox(4) = 90.0_dp
      pbcbox(5) = 90.0_dp
      pbcbox(6) = 90.0_dp

      call to_symmat(pbcbox, symmat)

      write(unit) symmat(1:6)

    end if


    ! write atom coordinates
    !
    natm = size(trajectory%coord(1,:))

    write(unit) real(trajectory%coord(1,1:natm),sp)
    write(unit) real(trajectory%coord(2,1:natm),sp)
    write(unit) real(trajectory%coord(3,1:natm),sp)

    return

  end subroutine write_dcd_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    seek_dcd_trj
  !> @brief        Seek DCD trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !! @param[in]    step_no    : step number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine seek_dcd_trj(file, trajectory, step_no)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(inout) :: trajectory
    integer,                 intent(in)    :: step_no

    ! local variables
    real(dp)                 :: symmat(6)
    integer                  :: unit, natm

!    real(4),   allocatable   :: crd4f(:)


    if (file%step_no > step_no) then
      call error_msg('Seek_Dcd_Trj> Seeking to the back is not supported.')
    else if (file%step_no == step_no) then
      return 
    end if

    unit = file%unit_no
    natm = size(trajectory%coord(1,:))

!    allocate(crd4f(natm))
    if (.not. allocated(crd4f)) then
      allocate(crd4f(natm))
    else
      if (size(crd4f(:)) /= natm) then
        deallocate(crd4f)
        allocate(crd4f(natm))
      end if
    end if

    do while (file%step_no /= step_no)

      if (file%trj_type == TrjTypeCoorBox) then
        read(unit) symmat(1:6)
      end if

      read(unit) crd4f(1:natm)
      read(unit) crd4f(1:natm)
      read(unit) crd4f(1:natm)

      file%step_no = file%step_no + 1

    end do
    
    !write(MsgOut,'(A,A,A,I7)') 'Seek_Dcd_Trj> file:', trim(file%name), &
    !                           ' current:', file%step_no

!    deallocate(crd4f)

    return

  end subroutine seek_dcd_trj

  !*********************************************************************
  ! gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_gromacs_trj
  !> @brief        Open Gromacs trajectory file
  !! @authors      NT
  !! @param[inout] file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_gromacs_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file

    ! local variables
    integer                  :: magic


#ifndef RICC
    file%unit_no = get_unit_no()

    if (file%in_out == IOFileInput) then

      open(file%unit_no,              &
           file    = file%name,       &
           status  = 'old',           &
           form    = 'unformatted',   &
           convert = 'little_endian', &
           access  = 'stream',        &
           err     = 900)
      read(file%unit_no) magic

      if (magic /= GRXMagic) then

        close(file%unit_no)
        open(file%unit_no,            &
             file    = file%name,     &
             status  = 'old',         &
             form    = 'unformatted', &
             convert = 'big_endian',  &
             access  = 'stream',      &
             err     = 900)
        read(file%unit_no) magic

        if (magic /= GRXMagic) &
          call error_msg('Open_Gromacs_Trj> unknown file format :'//&
          trim(file%name))

      end if

      rewind(file%unit_no)

    else

      open(file%unit_no,            &
           file    = file%name,     &
           status  = 'new',         &
           form    = 'unformatted', &
           access  = 'stream',      &
           err     = 900)

    end if
#else
#warning
#warning "WARNING: Cannot read/write gromacs dcd file on RICC platform."
#warning
    call error_msg( &
'Open_Gromacs_Trj> ERROR: Cannot read/write gromacs dcd file on RICC platform.')
#endif

    return

900 call error_msg('Open_Gromacs_Trj> open error :'//trim(file%name))

    return

  end subroutine open_gromacs_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_gromacs_trj
  !> @brief        Close Gromacs trajectory file
  !! @authors      NT
  !! @param[in]    file : trajectory file handle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_gromacs_trj(file)

    ! formal arguments
    type(s_trj_file),        intent(in)    :: file


    call close_file(file%unit_no)

    return

  end subroutine close_gromacs_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gromacs_trj
  !> @brief        Read Gromacs trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gromacs_trj(file, trajectory)

    ! parameters
    integer                  :: DIM = 3

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    real(dp)                 :: t_d, lambda_d, v_d
    real(sp)                 :: t_s, lambda_s, v_s
    integer                  :: magic, iverlen1, iverlen2, ir_size, e_size
    integer                  :: box_size, vir_size, pres_size, top_size
    integer                  :: sym_size, x_size, v_size, f_size, natom, step
    integer                  :: nre, real_size, i, j
    character(256)           :: version

    
    ! read header
    !
    read(file%unit_no) magic
    read(file%unit_no) iverlen1
    read(file%unit_no) iverlen2
    read(file%unit_no) version(1:iverlen2)
    read(file%unit_no) ir_size
    read(file%unit_no) e_size
    read(file%unit_no) box_size
    read(file%unit_no) vir_size
    read(file%unit_no) pres_size
    read(file%unit_no) top_size
    read(file%unit_no) sym_size
    read(file%unit_no) x_size
    read(file%unit_no) v_size
    read(file%unit_no) f_size
    read(file%unit_no) natom
    read(file%unit_no) step
    read(file%unit_no) nre

    if (natom /= size(trajectory%coord(1,:))) &
       call error_msg('Read_Gromacs_Trj> difference number of atom count')


    ! check real size
    !
    real_size = 0
    if (box_size /= 0) then
      real_size = box_size / (DIM**2)
    else if (x_size /= 0) then
      real_size = x_size / (natom * DIM)
    else if (v_size /= 0) then
      real_size = v_size / (natom * DIM)
    else if (f_size /= 0) then
      real_size = f_size / (natom * DIM)
    end if


    ! time
    if (real_size == 4) then
      read(file%unit_no) t_s
      read(file%unit_no) lambda_s
    else
      read(file%unit_no) t_d
      read(file%unit_no) lambda_d
    end if

    ! box size
    if (box_size /= 0) then

      if (real_size == 4) then

        do i = 1, 3
          do j = 1, 3
            read(file%unit_no) v_s
            trajectory%pbc_box(j,i) = 10.0_wp * real(v_s,wp)
          end do
        end do

      else

        do i = 1, 3
          do j = 1, 3
            read(file%unit_no) v_d
            trajectory%pbc_box(j,i) = 10.0_wp * real(v_d,wp)
          end do
        end do

      end if

    end if

    ! seek
    if (vir_size /= 0) &
      call fseek(file%unit_no, vir_size, 1)
    if (pres_size /= 0) &
      call fseek(file%unit_no, pres_size, 1)

    ! coordinate
    if (x_size /= 0) then

      if (real_size == 4) then
        
        do i = 1,natom
          do j = 1, 3
            read(file%unit_no) v_s
            trajectory%coord(j,i) = 10.0_wp * real(v_s, wp)
          end do
        end do

      else
        
        do i = 1,natom
          do j = 1, 3
            read(file%unit_no) v_d
            trajectory%coord(j,i) = 10.0_wp * real(v_d, wp)
          end do
        end do

      end if

    end if

    ! seek
    if (v_size /= 0) &
      call fseek(file%unit_no, v_size, 1)
    if (f_size /= 0) &
      call fseek(file%unit_no, f_size, 1)

    return

  end subroutine read_gromacs_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gromacs_trj_sel
  !> @brief        Write Gromacs trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    selatoms   : atom selection information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gromacs_trj_sel(file, trajectory, selatoms)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_selatoms),        intent(in)    :: selatoms

    ! local variables
    real(wp)                 :: t, lambda
    integer                  :: magic, iverlen1, iverlen2, ir_size, e_size
    integer                  :: box_size, vir_size, pres_size, top_size
    integer                  :: sym_size, x_size, v_size, f_size, natom, step
    integer                  :: nre, real_size, i, j
    character(256)           :: version
    

    version = 'GMX_trajectory_file'

    magic = GRXMagic
    write(file%unit_no) magic

    iverlen1 = len_trim(version) + 1
    write(file%unit_no) iverlen1

    iverlen2 = len_trim(version)
    write(file%unit_no) iverlen2

    write(file%unit_no) version(1:iverlen2)

    ir_size = 0
    write(file%unit_no) ir_size

    e_size = 0
    write(file%unit_no) e_size

    if (file%trj_type == TrjTypeCoorBox) then
      box_size = 9 * wp
    else
      box_size = 0
    end if
    write(file%unit_no) box_size

    vir_size = 0
    write(file%unit_no) vir_size

    pres_size = 0
    write(file%unit_no) pres_size

    top_size = 0
    write(file%unit_no) top_size

    sym_size = 0
    write(file%unit_no) sym_size

    x_size = size(selatoms%idx) * wp * 3
    write(file%unit_no) x_size

    v_size = 0
    write(file%unit_no) v_size

    f_size = 0
    write(file%unit_no) f_size

    natom = size(selatoms%idx)
    write(file%unit_no) natom

    step = file%step_no
    write(file%unit_no) step

    nre = 0
    write(file%unit_no) nre

    t = 0.0_wp
    write(file%unit_no) t

    lambda = 0.0_wp
    write(file%unit_no) lambda

    if (box_size /= 0) then
      write(file%unit_no) trajectory%pbc_box(1:3,1:3)/10_wp
    end if

    if (x_size /= 0) then
      write(file%unit_no) trajectory%coord(1:3,selatoms%idx(1:natom))/10.0_wp
    end if

    return

  end subroutine write_gromacs_trj_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gromacs_trj
  !> @brief        Write Gromacs trajectory file
  !! @authors      NT
  !! @param[inout] file       : trajectory file handle
  !! @param[in]    trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gromacs_trj(file, trajectory)

    ! formal arguments
    type(s_trj_file),        intent(inout) :: file
    type(s_trajectory),      intent(in)    :: trajectory

    ! local variables
    type(s_selatoms)         :: selatoms
    integer                  :: i, natom


    natom = size(trajectory%coord(1,:))

    call alloc_selatoms(selatoms, natom)
    do i = 1, natom
      selatoms%idx(i) = i
    end do

    call write_gromacs_trj_sel(file, trajectory, selatoms)

    call dealloc_selatoms(selatoms)

    return

  end subroutine write_gromacs_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    to_pbcbox
  !> @brief        to_pbcbox
  !! @authors      NT
  !! @param[in]    symmat : symmetric shape matrix
  !! @param[inout] pbcbox : periodic boundary condition box
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine to_pbcbox(symmat, pbcbox)

    ! formal arguments 
    real(dp),                intent(in)    :: symmat(:)
    real(dp),                intent(inout) :: pbcbox(:)

    ! local variables
    real(dp)                 :: ab, bc, ca


    pbcbox(1)= sqrt(symmat(1)*symmat(1)+symmat(2)*symmat(2)+symmat(4)*symmat(4))
    pbcbox(2)= sqrt(symmat(2)*symmat(2)+symmat(3)*symmat(3)+symmat(5)*symmat(5))
    pbcbox(3)= sqrt(symmat(4)*symmat(4)+symmat(5)*symmat(5)+symmat(6)*symmat(6))

    ab = symmat(2)*(symmat(1)+symmat(3))+symmat(4)*symmat(5)
    bc = symmat(5)*(symmat(3)+symmat(6))+symmat(2)*symmat(4)
    ca = symmat(4)*(symmat(1)+symmat(6))+symmat(2)*symmat(5)

    pbcbox(4) = to_deg(acos(bc/(pbcbox(2)*pbcbox(3))))
    pbcbox(5) = to_deg(acos(ca/(pbcbox(3)*pbcbox(1))))
    pbcbox(6) = to_deg(acos(ab/(pbcbox(1)*pbcbox(2))))

    return

  end subroutine to_pbcbox

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    to_symmat
  !> @brief        to_symmat
  !! @authors      NT
  !! @param[in]    pbcbox : periodic boundary condition box
  !! @param[inout] symmat : symmetric shape matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine to_symmat(pbcbox, symmat)

    ! formal arguments 
    real(dp),                intent(in)    :: pbcbox(:)
    real(dp),                intent(inout) :: symmat(:)


    if (pbcbox(4) - EPS > 90.0_wp .or. &
        pbcbox(4) + EPS < 90.0_wp .or. &
        pbcbox(5) - EPS > 90.0_wp .or. &
        pbcbox(5) + EPS < 90.0_wp .or. &
        pbcbox(6) - EPS > 90.0_wp .or. &
        pbcbox(6) + EPS < 90.0_wp) &
      write(MsgOut,*) &
      'To_Symmat> WARNING: Non rectangular unit cell is not supported.'

    symmat(1) = pbcbox(1)
    symmat(2) = 0.0_wp   ! angle A and B = 90.0 deg.
    symmat(3) = pbcbox(2)
    symmat(4) = 0.0_wp   ! angle A and C = 90.0 deg.
    symmat(5) = 0.0_wp   ! angle B and C = 90.0 deg.
    symmat(6) = pbcbox(3)

    return

  end subroutine to_symmat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      to_deg
  !> @brief        radian to degree
  !! @authors      NT, CK
  !! @param[in]    radian : radian
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function to_deg(radian)

    ! return value
    real(dp)                  :: to_deg

    ! formal arguments
    real(dp),                 intent(in)    :: radian


    to_deg = radian*RAD/PI

    return

  end function to_deg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_num_steps_trj
  !> @brief        determine number of steps contained in trajectory file
  !! @authors      DS
  !! @param[in]    filename   : trajectory file name
  !! @param[in]    trj_format : trajectory file format
  !! @param[in]    trj_type   : trajectory file type
  !! @return       number of steps
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_num_steps_trj(filename, trj_format, trj_type) result(num_steps)

    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: trj_format
    integer,                 intent(in)    :: trj_type

    integer                                :: num_steps

    type(s_trj_file)                       :: file


    select case(trj_format)

    ! currently, only dcd format is supported
    case (TrjFormatDCD)

      call open_trj(file, filename, trj_format, trj_type, IOFileInput)
      call read_dcd_trj_header(file, num_steps=num_steps)
      call close_trj(file)

    case default

      call error_msg("Get_Num_Steps_Trj> Determining number of steps for " &
        // trim(TrjFormatTypes(trj_format)) // " is not yet supported &
        &and must be provided manually")

    end select

    return

  end function get_num_steps_trj

end module fileio_trj_mod
