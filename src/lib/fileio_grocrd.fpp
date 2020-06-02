!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_grocrd_mod
!> @brief   read coordinate data from GROMACS GRO file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_grocrd_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none

  ! structures
  type s_grocrd

    integer                       :: num_atoms
    real(wp)                      :: pbc_box(3,3) = 0.0_wp

    ! atom
    integer,          allocatable :: residue_no(:)
    character(6),     allocatable :: residue_name(:)
    integer,          allocatable :: atom_no(:)
    character(6),     allocatable :: atom_name(:)
    real(wp),         allocatable :: atom_coord(:,:)
    real(wp),         allocatable :: atom_velocity(:,:)

  end type s_grocrd

  ! parameters for allocatable variables
  integer, parameter :: GrocrdAtom    = 1

  ! subroutines
  public  :: input_grocrd
  public  :: output_grocrd
  public  :: init_grocrd
  public  :: alloc_grocrd
  public  :: dealloc_grocrd
  public  :: dealloc_grocrd_all
  private :: read_grocrd
  private :: write_grocrd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_grocrd
  !> @brief        a driver subroutine for reading GRO CRD file
  !! @authors      NT
  !! @param[in]    grocrd_filename : filename of GRO CRD file
  !! @param[inout] grocrd          : structure of GRO CRD data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_grocrd(grocrd_filename, grocrd)

    ! formal arguments
    character(*),            intent(in)    :: grocrd_filename
    type(s_grocrd),          intent(inout) :: grocrd

    ! local variables
    integer                  :: file


    ! open GROMACS GRO file
    !
    call open_file(file, grocrd_filename, IOFileInput)


    ! read coordinate data from GROMACS GRO file
    !
    call read_grocrd(file, grocrd)


    ! close GROMACS GRO file
    !
    call close_file(file)

    return

  end subroutine input_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_grocrd
  !> @brief        a driver subroutine for writing GRO CRD file
  !! @authors      NT
  !! @param[in]    grocrd_filename : filename of GRO CRD file
  !! @param[in]    grocrd          : structure of GRO CRD data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_grocrd(grocrd_filename, grocrd)

    ! formal arguments
    character(*),            intent(in)    :: grocrd_filename
    type(s_grocrd),          intent(in)    :: grocrd

    ! local variables
    integer                  :: file


    ! open GROMACS GRO file
    !
    call open_file(file, grocrd_filename, IOFileOutputNew)


    ! write coordinate data from GROMACS GRO file
    !
    call write_grocrd(file, grocrd)


    ! close GROMACS GRO file
    !
    call close_file(file)

    return

  end subroutine output_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_grocrd
  !> @brief        initialize GROMACS GRO information
  !! @authors      NT
  !! @param[inout] grocrd          : structure of GRO CRD data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_grocrd(grocrd)

    ! formal arguments
    type(s_grocrd),          intent(inout) :: grocrd


    grocrd%num_atoms        = 0
    grocrd%pbc_box(1:3,1:3) = 0.0_wp

    return

  end subroutine init_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_grocrd
  !> @brief        allocate GROMACS GRO information
  !! @authors      NT
  !! @param[inout] grocrd   : structure of GROMACS GRO information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_grocrd(grocrd, variable, var_size)

    ! formal arguments
    type(s_grocrd),          intent(inout) :: grocrd
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(GrocrdAtom)

      if (allocated(grocrd%residue_no)) then
        if (size(grocrd%residue_no) == var_size) return
        deallocate(grocrd%residue_no,    &
                   grocrd%residue_name,  &
                   grocrd%atom_no,       &
                   grocrd%atom_name,     &
                   grocrd%atom_coord,    &
                   grocrd%atom_velocity, &
                   stat = dealloc_stat   )
      end if

      allocate(grocrd%residue_no(var_size),      &
               grocrd%residue_name(var_size),    &
               grocrd%atom_no(var_size),         &
               grocrd%atom_name(var_size),       &
               grocrd%atom_coord(3,var_size),    &
               grocrd%atom_velocity(3,var_size), &
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Grocrd> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_grocrd
  !> @brief        deallocate GROMACS GRO information
  !! @authors      NT
  !! @param[inout] grocrd   : structure of GROMACS GRO information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_grocrd(grocrd, variable)

    ! formal arguments
    type(s_grocrd),          intent(inout) :: grocrd
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case(GrocrdAtom)

      if (allocated(grocrd%residue_no)) then
        deallocate(grocrd%residue_no,    &
                   grocrd%residue_name,  &
                   grocrd%atom_no,       &
                   grocrd%atom_name,     &
                   grocrd%atom_coord,    &
                   grocrd%atom_velocity, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Grocrd> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_grocrd_all
  !> @brief        deallocate all GROMACS GRO information
  !! @authors      NT
  !! @param[inout] grocrd   : structure of GROMACS GRO information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_grocrd_all(grocrd)

    ! formal arguments
    type(s_grocrd),          intent(inout) :: grocrd


    call dealloc_grocrd(grocrd, GrocrdAtom)

    return

  end subroutine dealloc_grocrd_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_grocrd
  !> @brief        read data from GROMACS GRO file
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[inout] grocrd : structure of GROMACS GRO information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_grocrd(file, grocrd)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grocrd),          intent(inout) :: grocrd

    ! local variables
    integer                  :: i, j, num_atoms
    character(MaxLine)       :: line
    character(5)             :: atom_name
    character(5)             :: resi_name


    read(file,'(a)') line
    read(file,*) num_atoms

    call alloc_grocrd(grocrd, GroCRDAtom, num_atoms)

    do i = 1, num_atoms
      read(file,'(i5,a5,a5,i5,3f8.3,3f8.3)') &
           grocrd%residue_no(i), &
           resi_name, &
           atom_name, &
           grocrd%atom_no(i), &
           (grocrd%atom_coord(j,i),j=1,3), &
           (grocrd%atom_velocity(j,i),j=1,3)

      read(resi_name,*) grocrd%residue_name(i)
      read(atom_name,*) grocrd%atom_name(i)
    end do

    grocrd%num_atoms = num_atoms

    read(file,*) (grocrd%pbc_box(i,i),i=1,3)

    return

  end subroutine read_grocrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_grocrd
  !> @brief        write data from GROMACS GRO file
  !! @authors      NT
  !! @param[in]    file   : file unit number
  !! @param[in]    grocrd : structure of GROMACS GRO information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_grocrd(file, grocrd)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_grocrd),          intent(in)    :: grocrd

    ! local variables
    integer                  :: i, j
    character(5)             :: resi_name, atom_name


    write(file,'(a)') 'Gro file for Gromacs: '//PACKAGE_NAME//': '//          &
                      PACKAGE_VERSION
    write(file,'(i0)') grocrd%num_atoms

    do i = 1, grocrd%num_atoms
      write(resi_name,'(a)')  trim(grocrd%residue_name(i)) ! left align   ..why?
      write(atom_name,'(a5)') trim(grocrd%atom_name(i))    ! right align  ..why?
      write(file,'(i5,a5,a5,i5,3f8.3,3f8.3)') &
           grocrd%residue_no(i), &
           resi_name,            &
           atom_name,            &
           grocrd%atom_no(i),    &
           (grocrd%atom_coord   (j,i),j=1,3), &
           (grocrd%atom_velocity(j,i),j=1,3)
    end do

    write(file,'(3f10.5)') (grocrd%pbc_box(i,i),i=1,3)

    return

  end subroutine write_grocrd

end module fileio_grocrd_mod
