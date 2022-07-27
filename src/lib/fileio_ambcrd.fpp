!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_ambcrd_mod
!> @brief   read coordinate data from AMBER CRD file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_ambcrd_mod

  use fileio_mod
  use messages_mod
  use string_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_ambcrd

    logical                       :: time_rec
    logical                       :: velocity_rec
    logical                       :: box_rec
    logical                       :: box_angl_rec

    character(80)                 :: title
    integer                       :: num_atoms
    real(wp)                      :: time

    real(wp),         allocatable :: atom_coord(:,:)
    real(wp),         allocatable :: atom_velocity(:,:)
    real(wp)                      :: box_size(3)
    real(wp)                      :: box_angl(3)

  end type s_ambcrd

  ! parameters for allocatable variables
  integer,     public, parameter :: PrmcrdAtom = 1

  ! privates
  public  :: init_ambcrd
  public  :: alloc_ambcrd
  public  :: dealloc_ambcrd
  public  :: dealloc_ambcrd_all
  public  :: input_ambcrd
  public  :: output_ambcrd
  private :: read_ambcrd
  private :: write_ambcrd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_ambcrd
  !> @brief        initialize AMBER CRD information
  !! @authors      NT
  !! @param[in]    ambcrd : structure of AMBER CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_ambcrd(ambcrd)

    ! formal arguments
    type(s_ambcrd),          intent(inout) :: ambcrd


    ambcrd%time_rec      = .false.
    ambcrd%velocity_rec  = .false.
    ambcrd%box_rec       = .false.
    ambcrd%box_angl_rec  = .false.

    ambcrd%num_atoms     = 0
    ambcrd%time          = 0.0_wp

    ambcrd%title         = ''

    ambcrd%box_size(1:3) = 0.0_wp
    ambcrd%box_angl(1:3) = 0.0_wp

    return

  end subroutine init_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_ambcrd
  !> @brief        allocate AMBER CRD information
  !! @authors      NT
  !! @param[inout] ambcrd   : structure of AMBER CRD information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_ambcrd(ambcrd, variable, var_size)

    ! formal arguments
    type(s_ambcrd),          intent(inout) :: ambcrd
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

    case(PrmcrdAtom)

      if (allocated(ambcrd%atom_coord)) then
        if (size(ambcrd%atom_coord(1,:)) == var_size) return
        deallocate(ambcrd%atom_coord,    &
                   ambcrd%atom_velocity, &
                   stat = dealloc_stat)
      end if

      allocate(ambcrd%atom_coord   (3,var_size), &
               ambcrd%atom_velocity(3,var_size), &
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Ambcrd> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_ambcrd
  !> @brief        deallocate AMBER CRD information
  !! @authors      NT
  !! @param[inout] ambcrd   : structure of AMBER CRD information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_ambcrd(ambcrd, variable)

    ! formal arguments
    type(s_ambcrd),          intent(inout) :: ambcrd
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variable
    !  
    select case (variable)

    case(PrmcrdAtom)
      if (allocated(ambcrd%atom_coord)) then
        deallocate(ambcrd%atom_coord,    &
                   ambcrd%atom_velocity, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Ambcrd> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_ambcrd_all
  !> @brief        deallocate all AMBER CRD information
  !! @authors      NT
  !! @param[in]    ambcrd : structure of AMBER CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_ambcrd_all(ambcrd)

    ! formal arguments
    type(s_ambcrd),          intent(inout) :: ambcrd


    call dealloc_ambcrd(ambcrd, PrmcrdAtom)

    return

  end subroutine dealloc_ambcrd_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_ambcrd
  !> @brief        a driver subroutine for reading AMBER CRD file
  !! @authors      NT
  !! @param[in]    ambcrd_filename : filename of AMBER CRD file
  !! @param[inout] ambcrd          : structure of AMBER CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_ambcrd(ambcrd_filename, ambcrd)

    ! formal arguments
    character(*),            intent(in)    :: ambcrd_filename
    type(s_ambcrd),          intent(inout) :: ambcrd
   
    ! local variables
    integer                  :: file
    

    ! open AMBER CRD file
    !
    call open_file(file, ambcrd_filename, IOFileInput)


    ! read coordinate data from AMBER CRD file
    !
    call read_ambcrd(file, ambcrd)


    ! close AMBER CRD file
    !
    call close_file(file)

    return

  end subroutine input_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_ambcrd
  !> @brief        a driver subroutine for writing AMBER CRD file
  !! @authors      NT
  !! @param[in]    ambcrd_filename : filename of AMBER CRD file
  !! @param[in]    ambcrd          : structure of AMBER CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_ambcrd(ambcrd_filename, ambcrd)

    ! formal arguments
    character(*),            intent(in)    :: ambcrd_filename
    type(s_ambcrd),          intent(in)    :: ambcrd
   
    ! local variables
    integer                  :: file

    
    ! open AMBER CRD file
    !
    call open_file(file, ambcrd_filename, IOFileOutputNew)


    ! write coordinate data from AMBER CRD file
    !
    call write_ambcrd(file, ambcrd)


    ! close AMBER CRD file
    !
    call close_file(file)

    return

  end subroutine output_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ambcrd
  !> @brief        read data from AMBER CRD file
  !! @authors      NT
  !! @param[in]    file   : unit number of AMBER CRD file
  !! @param[inout] ambcrd : ambcrd data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ambcrd(file, ambcrd)

    ! parameters
    real(wp), parameter      :: CrdEmpty = 9999999.0_wp

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_ambcrd),          intent(inout) :: ambcrd

    ! local variables
    integer                  :: i
    character(MaxLine)       :: line


    ! read ITITL
    !
    read(file,'(a)') line
    ambcrd%title = line


    ! read NATOM, TIME
    !
    ambcrd%num_atoms = 0
    ambcrd%time = 0.0_wp

    read(file,'(a)') line
    read(line,*,end=10,err=900) ambcrd%num_atoms, ambcrd%time

10  if (ambcrd%time /= 0.0_wp) &
      ambcrd%time_rec = .true.

    call alloc_ambcrd(ambcrd, PrmcrdAtom, ambcrd%num_atoms)


    ! read (X(i),Y(i),Z(i), i=1,NATOM)
    !
    read(file,'(6F12.7)') &
         (ambcrd%atom_coord(1:3,i), i = 1, ambcrd%num_atoms)


    ! read (VX(i),VY(i),VZ(i), i=1,NATOM)
    !
    do i = 1, ambcrd%num_atoms
      ambcrd%atom_velocity(1:3,i) = CrdEmpty
    end do

    read(file,'(6F12.7)',end=20,err=900) &
         (ambcrd%atom_velocity(1:3,i), i = 1, ambcrd%num_atoms)

20  if (ambcrd%atom_velocity(1,ambcrd%num_atoms) /= CrdEmpty) &
      ambcrd%velocity_rec = .true.


    ! read (BOX(1),BOX(2),BOX(3)
    !
    ambcrd%box_size(1:3) = 0.0_wp
    ambcrd%box_angl(1:3) = 0.0_wp

    if (ambcrd%velocity_rec) then

      read(file,'(6F12.7)',end=30) ambcrd%box_size(1:3), ambcrd%box_angl(1:3)

    else

      if (ambcrd%atom_velocity(1,1) /= CrdEmpty) &
        ambcrd%box_size(1:3) = ambcrd%atom_velocity(1:3,1)

      if (ambcrd%atom_velocity(1,2) /= CrdEmpty) &
        ambcrd%box_angl(1:3) = ambcrd%atom_velocity(1:3,2)

    end if

30  if (ambcrd%box_size(1) /= 0.0_wp) &
      ambcrd%box_rec = .true.

    if (ambcrd%box_angl(1) /= 0.0_wp) &
      ambcrd%box_angl_rec = .true.


    ! cleanup
    !
    if (.not. ambcrd%velocity_rec) then
      do i = 1, ambcrd%num_atoms
        ambcrd%atom_velocity(1:3,i) = 0.0_wp
      end do
    end if

    return

900 call error_msg('Read_Ambcrd> Format error.')

  end subroutine read_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_ambcrd
  !> @brief        write data to AMBER CRD file
  !! @authors      NT
  !! @param[in]    file   : unit number of AMBER CRD file
  !! @param[in]    ambcrd : ambcrd data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_ambcrd(file, ambcrd)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_ambcrd),          intent(in)    :: ambcrd

    ! local variables
    integer                  :: i


    ! write TITLE
    !

    write(file,'(20A4)') ambcrd%title

    
    ! write NATOM, TIME
    !

    if (ambcrd%time_rec) then
      write(file,'(I5,5E15.7)') ambcrd%num_atoms, ambcrd%time

    else
      write(file,'(I5)') ambcrd%num_atoms

    end if


    ! write (X(i),Y(i),Z(i),i=1,NATOM)
    !

    write(file,'(6F12.7)') (ambcrd%atom_coord(1:3,i),i=1,ambcrd%num_atoms)


    ! write (VX(i),VY(i),VZ(i),i=1,NATOM)
    !

    if (ambcrd%velocity_rec) &
      write(file,'(6F12.7)') (ambcrd%atom_velocity(1:3,i),i=1,ambcrd%num_atoms)


    ! write (BOX(1),BOX(2),BOX(3)
    !

    if (ambcrd%box_rec .and. .not. ambcrd%box_angl_rec) then
      write(file,'(6F12.7)') ambcrd%box_size(1:3)

    else if (ambcrd%box_rec .and. ambcrd%box_angl_rec) then
      write(file,'(6F12.7)') ambcrd%box_size(1:3), ambcrd%box_angl(1:3)

    end if

    return

  end subroutine write_ambcrd

end module fileio_ambcrd_mod
