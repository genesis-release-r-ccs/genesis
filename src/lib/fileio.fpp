!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_mod
!> @brief   utilities for file I/O
!! @authors Yuji Sugita (YS), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_mod

  use messages_mod

  implicit none
  private

  ! constants
  integer, public, parameter :: InvalidUnitNo       = 1
  integer, public, parameter :: FirstUnitNo         = 7
  integer, public, parameter :: MaxUnitCount        = 1000

  integer, public, parameter :: IOFileInput         = 1
  integer, public, parameter :: IOFileOutputNew     = 2
  integer, public, parameter :: IOFileOutputReplace = 3

  integer, public, parameter :: IOFileBigEndian     = 1
  integer, public, parameter :: IOFileLittleEndian  = 2
  integer, public, parameter :: IOFileNativeEndian  = 3


  ! module variables
  logical, private           :: io_initial = .false.
  logical, private           :: io_used   (MaxUnitCount)
  integer, private           :: io_unit_no(MaxUnitCount)


  ! subroutines and functions
  public   :: open_file
  public   :: open_binary_file
  public   :: close_file
  public   :: get_unit_no
  public   :: free_unit_no
  public   :: get_extension

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_file
  !> @brief        open a text file
  !! @authors      YS, NT
  !! @param[out]   unit_no  : unit number of a text file
  !! @param[in]    filename : file name of a text file
  !! @param[in]    in_out   : flag for the file (IOFileInput/IOFileOutput)
  !! @note         get_unit_no is used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_file(unit_no, filename, in_out)

    ! formal arguments
    integer,                 intent(out)   :: unit_no
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out


    unit_no = get_unit_no()
    if (unit_no == InvalidUnitNo) goto 900

    if (in_out == IOFileInput) then
      open(unit_no, file=filename, status='old', form='formatted', err=900)

    else if (in_out == IOFileOutputNew) then
      open(unit_no, file=filename, status='new', form='formatted', err=900)

    else if (in_out == IOFileOutputReplace) then
      open(unit_no, file=filename, status='replace', form='formatted', err=900)

    endif

    return

900 write(ErrOut,'(A)') 'Open_File> ', trim(filename)
    call error_msg('Open_file> Error in opening ' // trim(filename)       // &
                   ' (Input file does not exist, or Output file already ' // &
                   'exists, or other reasons)')

  end subroutine open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_binary_file
  !> @brief        open a binary file
  !! @authors      YS, NT
  !! @param[out]   unit_no   : unit number of a file
  !! @param[in]    filename  : file name of a file
  !! @param[in]    in_out    : flag for the file (IOFileInput/IOFileOutput)
  !! @param[in]    in_endian : endian for the file input
  !! @note         get_unit_no is used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_binary_file(unit_no, filename, in_out, in_endian)

    ! parameters
    character(20),           parameter     :: EndianStr(3) =           &
                                                    (/'big_endian   ', &
                                                      'little_endian', &
                                                      'native       '/)

    ! formal arguments
    integer,                 intent(out)   :: unit_no
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out
    integer,       optional, intent(in)    :: in_endian


    unit_no = get_unit_no()
    if (unit_no == InvalidUnitNo) goto 900

    if (in_out == IOFileInput) then

      if (present(in_endian)) then

        if (in_endian < 0 .or. size(EndianStr) < in_endian) &
          call error_msg('Open_Binary_File> ERROR : bad endian.')

        open(unit_no,                 &
             file    = filename,      &
             status  = 'old',         &
             form    = 'unformatted', &
             err     = 900,           &
             convert = EndianStr(in_endian))
      else

        open(unit_no,                 &
             file    = filename,      &
             status  = 'old',         &
             form    = 'unformatted', &
             err     = 900)

      end if

    else if (in_out == IOFileOutputNew) then

      open(unit_no,                &
           file   = filename,      &
           status = 'new',         &
           form   = 'unformatted', &
           err=900)

    else if (in_out == IOFileOutputReplace) then

      open(unit_no,                &
           file   = filename,      &
           status = 'replace',     &
           form   = 'unformatted', &
           err    = 900)

    endif

    return

900 write(ErrOut,'(A)') 'Open_Binary_File> ', trim(filename)
    call error_msg('Open_Binary_file> Error in opening ' // trim(filename) // &
                   ' (Input file does not exist, or Output file already '  // &
                   'exists, or other reasons)')

  end subroutine open_binary_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_file
  !> @brief        close a file
  !! @authors      YS
  !! @param[in]    unit_no : unit number of a file
  !! @note         free_unit_no is used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_file(unit_no)

    ! formal arguments
    integer,                 intent(in)    :: unit_no

    if((unit_no .eq. MsgOut) .or. (unit_no .eq. ErrOut) &
       .or. (unit_no .eq. 6) .or. (unit_no .eq. 0)) return

    close(unit_no)
    call free_unit_no(unit_no)

    return

  end subroutine close_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_unit_no
  !> @brief        create and get unit number of a file
  !! @authors      YS
  !! @return       unit_no of a file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_unit_no()

    ! return value
    integer :: get_unit_no

    ! local variables
    integer :: i


    ! initialize table for unit no
    !
    if (.not. io_initial) then
      io_initial = .true.
      do i = 1, MaxUnitCount
        io_used(i)    = .false.
        io_unit_no(i) = i + FirstUnitNo - 1
      end do
    end if

    ! get unit no
    !
    do i = 1, MaxUnitCount
      if (.not. io_used(i)) then
        io_used(i)  = .true.
        get_unit_no = io_unit_no(i)
        return
      end if
    end do

    get_unit_no = InvalidUnitNo

  end function get_unit_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    free_unit_no
  !> @brief        free unit number of a file
  !! @authors      YS
  !! @param[in]    unit_no : unit number of an input/output file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine free_unit_no(unit_no)

    ! formal arguments
    integer,                 intent(in)    :: unit_no

    ! local variables
    integer                  :: i


    do i = 1, MaxUnitCount
      if (unit_no == io_unit_no(i)) then
        io_used(i) = .false.
        return
      end if
    end do

    return

  end subroutine free_unit_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_extension
  !> @brief        get the file extension
  !! @authors      NT
  !! @param[in]    filename : file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_extension(filename)

    ! return value
    character(10) :: get_extension

    ! formal arguments
    character(*),            intent(in)    :: filename

    ! local variables
    integer                  :: idx_period, idx_slash


    get_extension = ''

    idx_period = index(filename, '.', back=.true.)
    if (idx_period == 0) &
      return

    idx_slash = index(filename, '/', back=.true.)
    if (idx_slash > idx_period) &
      return

    get_extension = filename(idx_period + 1:)

    return

  end function get_extension

end module fileio_mod
