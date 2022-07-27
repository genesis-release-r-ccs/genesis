!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_atp_mod
!> @brief   read/write GROMACS atom type file
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module fileio_atp_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_atp

    character(4), allocatable  :: name(:)
    real(wp),     allocatable  :: mass(:)

  end type s_atp

  ! parameters for allocatable variables
  integer, parameter :: AtpAtomType = 1

  ! subroutines
  public  :: input_atp
  public  :: alloc_atp
  public  :: dealloc_atp
  public  :: dealloc_atp_all
  public  :: merge_atp
  private :: read_atp

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_atp
  !> @brief        a driver subroutine for reading GROMACS ATP file
  !! @authors      NT
  !! @param[in]    atp_filename : filename of ATP file
  !! @param[inout] atp          : structure of ATP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_atp(atp_filename, atp)

    ! formal arguments
    character(*),            intent(in)    :: atp_filename
    type(s_atp),             intent(inout) :: atp

    ! local variables
    type(s_atp)              :: atp0
    integer                  :: file, i, nfile

    character(MaxFilename), allocatable :: files(:)


    ! parse filename string
    !

    nfile = split_num(atp_filename)
    allocate(files(nfile))

    call split(nfile, nfile, atp_filename, files)


    do i = 1, nfile

      ! open ATP file
      !
      call open_file(file, files(i), IOFileInput)


      ! read coordinate data from ATP file
      !
      call read_atp(file, atp0)


      ! close ATP file
      !
      call close_file(file)


      ! merge ATP data
      !
      call merge_atp(atp, atp0)

      call dealloc_atp_all(atp0)

    end do

    deallocate(files)

    return

  end subroutine input_atp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_atp
  !> @brief        allocate GROMACS ATP data
  !! @authors      NT
  !! @param[inout] atp      : structure of GROMACS ATP data
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_atp(atp, variable, var_size)

    ! formal arguments
    type(s_atp),             intent(inout) :: atp
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

    case(AtpAtomType)

      if (allocated(atp%name)) then
        if (size(atp%name) == var_size) return
        deallocate(atp%name, &
                   atp%mass, &
                   stat = dealloc_stat)
      end if

      allocate(atp%name(var_size),&
               atp%mass(var_size), &
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Atp> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return 

  end subroutine alloc_atp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_atp
  !> @brief        deallocate GROMACS ATP data
  !! @authors      NT
  !! @param[inout] atp      : structure of GROMACS ATP data
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_atp(atp, variable)

    ! formal arguments
    type(s_atp),             intent(inout) :: atp
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(AtpAtomType)

      if (allocated(atp%name)) then
        deallocate(atp%name, &
                   atp%mass, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Atp> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return 

  end subroutine dealloc_atp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_atp_all
  !> @brief        deallocate all GROMACS ATP data
  !! @authors      NT
  !! @param[inout] atp : structure of GROMACS ATP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_atp_all(atp)

    ! formal arguments
    type(s_atp),             intent(inout) :: atp


    call dealloc_atp(atp, AtpAtomType)

    return

  end subroutine dealloc_atp_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_atp
  !> @brief        merge two ATP data
  !! @authors      NT
  !! @param[inout] atp0 : structure of GROMACS ATP data
  !! @param[in]    atp1 : structure of GROMACS ATP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine merge_atp(atp0, atp1)

    ! formal arguments
    type(s_atp),             intent(inout) :: atp0
    type(s_atp),             intent(in)    :: atp1

    ! local variables
    type(s_atp)              :: atp
    integer                  :: i, n, nn


    if (.not. allocated(atp0%name)) then

      n = size(atp1%name)
      call alloc_atp(atp0, AtpAtomType, n)
      atp0%name(1:n) = atp1%name(1:n)
      atp0%mass(1:n) = atp1%mass(1:n)
      return

    end if

    ! copy atp0 -> atp
    !

    call alloc_atp(atp, AtpAtomType, size(atp0%name))

    n = size(atp0%name)
    atp%name(1:n) = atp0%name(1:n)
    atp%mass(1:n) = atp0%mass(1:n)


    ! re-allocate atp0
    !
    call dealloc_atp_all(atp0)

    call alloc_atp(atp0, AtpAtomType, size(atp%name) + size(atp1%name))


    ! copy atp  -> atp0
    !

    n = size(atp%name)
    atp0%name(1:n) = atp%name(1:n)
    atp0%mass(1:n) = atp%mass(1:n)


    ! copy atp1 -> atp0
    !

    nn = size(atp%name)
    n  = size(atp1%name)
    atp0%name(nn+1:nn+n) = atp1%name(1:n)
    atp0%mass(nn+1:nn+n) = atp1%mass(1:n)


    ! deallocate
    !

    call dealloc_atp_all(atp)

    return

  end subroutine merge_atp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_atp
  !> @brief        read data from ATP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] atp  : structure of ATP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_atp(file, atp)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_atp),             intent(inout) :: atp

    ! local variables
    integer                  :: i, ntype
    character(200)           :: line


    ntype = 0
    do while (.true.)
      read(file,'(a)',end=100) line
      ntype = ntype + 1
    end do

100 rewind(file)

    call alloc_atp(atp, AtpAtomType, ntype)

    do i = 1, ntype
      read(file,'(a)') line
      read(line,*) atp%name(i), atp%mass(i)
    end do

    ! write summary of ATP information
    !
    if (main_rank) then

      write(MsgOut,'(a)') 'Read_Atp> Summary of Atp file'
      write(MsgOut,'(a20,i10)') &
           '  num_atomtypes   = ', size(atp%name)
      write(MsgOut,'(a)') ''

    end if

    return 

  end subroutine read_atp

end module fileio_atp_mod
