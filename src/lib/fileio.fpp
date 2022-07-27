!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_mod
!> @brief   utilities for file I/O
!! @authors Yuji Sugita (YS), Norio Takase (NT), Donatas Surblys(DS),
!           Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_mod

  use constants_mod
  use messages_mod
  use string_mod
  use mpi_parallel_mod

  implicit none
  private

  ! constants
  integer, public, parameter :: InvalidUnitNo       = 1
  integer, public, parameter :: FirstUnitNo         = 7
  integer, public, parameter :: MaxUnitCount        = 1000

  integer, public, parameter :: IOFileInput         = 1
  integer, public, parameter :: IOFileOutputNew     = 2
  integer, public, parameter :: IOFileOutputReplace = 3
  integer, public, parameter :: IOFileOutputAppend  = 4

  integer, public, parameter :: IOFileBigEndian     = 1
  integer, public, parameter :: IOFileLittleEndian  = 2
  integer, public, parameter :: IOFileNativeEndian  = 3

  integer, public, parameter :: IOFileSequentialAccess = 1
  integer, public, parameter :: IOFileDirectAccess     = 2
  integer, public, parameter :: IOFileStreamAccess     = 3

  ! module variables
  logical, private           :: io_initial = .false.
  logical, private           :: io_used   (MaxUnitCount)
  integer, private           :: io_unit_no(MaxUnitCount)


  ! file backup options
  character(*), private, parameter :: backup_default_suffix   = "backup"
  integer,      private, parameter :: backup_max_no           = 99
  logical,      save               :: backup_allowed          = .false.

  ! subroutines and functions
  public   :: open_file
  public   :: open_binary_file
  public   :: close_file
  public   :: get_unit_no
  public   :: free_unit_no
  public   :: get_extension
  public   :: setup_backup
  private  :: backup_file
  private  :: on_invalid_unit_number
  private  :: on_input_file_does_not_exist
  private  :: on_output_file_exists
  private  :: on_unknown_io_error


contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_file
  !> @brief        open a text file
  !! @authors      YS, NT, DS, KY
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
    if (unit_no == InvalidUnitNo) goto 901

    if (in_out == IOFileInput) then
      open(unit_no, file=filename, status='old', form='formatted', err=902)

    else if (in_out == IOFileOutputNew) then
      open(unit_no, file=filename, status='new', form='formatted', err=903)

    else if (in_out == IOFileOutputReplace) then
      open(unit_no, file=filename, status='replace', form='formatted', err=904)

    else if (in_out == IOFileOutputAppend) then
      open(unit_no, file=filename, status='unknown', form='formatted', &
                    access='append', err=904)

    endif

    return

901 call on_invalid_unit_number(unit_no, filename, "formatted")
    return

902 call on_input_file_does_not_exist(unit_no, filename, "formatted")
    return

903 call on_output_file_exists(unit_no, filename, "formatted")
    return

904 call on_unknown_io_error(unit_no, filename, "formatted")
    return

  end subroutine open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_binary_file
  !> @brief        open a binary file
  !! @authors      YS, NT, DS, TM
  !! @param[out]   unit_no   : unit number of a file
  !! @param[in]    filename  : file name of a file
  !! @param[in]    in_out    : flag for the file (IOFileInput/IOFileOutput)
  !! @param[in]    in_endian : endian for the file input
  !! @param[in]    in_access : access mode for file io
  !! @note         get_unit_no is used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_binary_file(unit_no, filename, in_out, in_endian, in_access)

    ! parameters
    character(20),           parameter     :: EndianStr(3) =           &
                                                    (/'big_endian   ', &
                                                      'little_endian', &
                                                      'native       '/)
    character(20),           parameter     :: AccessStr(3) =           &
                                                    (/'sequential   ', &
                                                      'direct       ', &
                                                      'stream       '/)

    ! formal arguments
    integer,                 intent(out)   :: unit_no
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out
    integer,       optional, intent(in)    :: in_endian
    integer,       optional, intent(in)    :: in_access


    unit_no = get_unit_no()
    if (unit_no == InvalidUnitNo) goto 901

    if (in_out == IOFileInput) then

      if (present(in_endian)) then

        if (in_endian < 0 .or. size(EndianStr) < in_endian) &
          call error_msg('Open_Binary_File> ERROR : bad endian.')

        if (present(in_access)) then

          if (in_access < 0 .or. size(AccessStr) < in_access) &
            call error_msg('Open_Binary_File> ERROR : bad access mode.')

            open(unit_no,                        &
                 file    = filename,             &
                 status  = 'old',                &
                 form    = 'unformatted',        &
                 err     = 902,                  &
                 access  = AccessStr(in_access), &
                 convert = EndianStr(in_endian))

        else

          open(unit_no,                 &
               file    = filename,      &
               status  = 'old',         &
               form    = 'unformatted', &
               err     = 902,           &
               convert = EndianStr(in_endian))

        end if

      else

        open(unit_no,                 &
             file    = filename,      &
             status  = 'old',         &
             form    = 'unformatted', &
             err     = 902)

      end if

    else if (in_out == IOFileOutputNew) then

      if (present(in_endian)) then

        if (in_endian < 0 .or. size(EndianStr) < in_endian) &
          call error_msg('Open_Binary_File> ERROR : bad endian.')

        if (present(in_access)) then

          if (in_access < 0 .or. size(AccessStr) < in_access) &
            call error_msg('Open_Binary_File> ERROR : bad access mode.')

            open(unit_no,                        &
                 file    = filename,             &
                 status  = 'new',                &
                 form    = 'unformatted',        &
                 err     = 903,                  &
                 access  = AccessStr(in_access), &
                 convert = EndianStr(in_endian))

        else

          open(unit_no,                 &
               file    = filename,      &
               status  = 'new',         &
               form    = 'unformatted', &
               err     = 903,           &
               convert = EndianStr(in_endian))

        end if

      else
        open(unit_no,                &
             file   = filename,      &
             status = 'new',         &
             form   = 'unformatted', &
             err=903)

      end if

    else if (in_out == IOFileOutputReplace) then

      open(unit_no,                &
           file   = filename,      &
           status = 'replace',     &
           form   = 'unformatted', &
           err    = 904)
      
    else if (in_out == IOFileOutputAppend) then

      open(unit_no,                &
           file   = filename,      &
           status = 'unknown',     &
           form   = 'unformatted', &
           position = 'append',      &
           err    = 905)
    endif

    return

901 call on_invalid_unit_number(unit_no, filename, "unformatted")
    return

902 call on_input_file_does_not_exist(unit_no, filename, "unformatted")
    return

903 call on_output_file_exists(unit_no, filename, "unformatted")
    return

904 call on_unknown_io_error(unit_no, filename, "unformatted")
    return

905 call on_unknown_io_error(unit_no, filename, "unformatted")
    return
    
  end subroutine open_binary_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_file
  !> @brief        close a file
  !! @authors      YS, KY
  !! @param[in]    unit_no : unit number of a file
  !! @note         free_unit_no is used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_file(unit_no, stat)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    character(*), optional,  intent(in)    :: stat


    if((unit_no .eq. MsgOut) .or. (unit_no .eq. ErrOut) &
       .or. (unit_no .eq. 6) .or. (unit_no .eq. 0)) return

    if (present(stat)) then
      close(unit_no, status=trim(stat))
    else
      close(unit_no)
    end if
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

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_backup
  !> @brief        define backup yes or no
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_backup(do_backup)

    ! formal arguments
    logical,   intent(in) :: do_backup

    
    if (do_backup) then
      backup_allowed = .true.
    else
      backup_allowed = .false.
    end if

    return

  end subroutine setup_backup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    backup_file
  !> @brief        backup (rename) a file in a predefined pattern
  !! @authors      DS
  !! @param[in]    filename : file name of a file
  !! @param[in]    suffix   : backup file name suffix (optinal)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine backup_file(filename, suffix)

    ! formal arguments
    character(*),           intent(in) :: filename
    character(*), optional, intent(in) :: suffix

    character(MaxLine)      :: backup_suffix
    character(MaxFilename)  :: backup_filename
    character(MaxLine)      :: num_str

    logical                 :: file_exists
    integer                 :: backup_no
    integer                 :: stat, rename


    if (present(suffix)) then
      backup_suffix = suffix
    else
      backup_suffix = backup_default_suffix
    end if
    backup_suffix = "."//trim(backup_suffix)//"."

    do backup_no = 1, backup_max_no

      write(num_str, "(i0)") backup_no

      backup_filename = trim(filename)//trim(backup_suffix)//trim(num_str)

      inquire(file=trim(backup_filename), exist=file_exists)

      if (.not. file_exists) then
        write(MsgOut, '("Backup_file> Backing up ",A," to ",A)') &
          trim(filename), trim(backup_filename)

        stat = rename(trim(filename),trim(backup_filename))

        if (stat /= 0) then
          call error_msg('Backup_file> Error in making a backup file') 
        end if

        return
      end if

    end do

    write(num_str, "(i0)") backup_max_no
    call error_msg('Backup_file> Number of backups ('//trim(num_str) &
      // ') exceeded for ' // trim(filename))

    return

  end subroutine backup_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    on_invalid_unit_number
  !> @brief        invalid unit number error handler
  !! @authors      DS
  !! @param[in]    unit_no  : unit number of a file
  !! @param[in]    filename : file name of a file
  !! @param[in]    form     : file transfer type (formatted/unformatted)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine on_invalid_unit_number(unit_no, filename, form)

    ! formal arguments
    integer,      intent(in) :: unit_no
    character(*), intent(in) :: filename
    character(*), intent(in) :: form


    call error_msg('Open_file> Invalid unit number while opening ' &
      // trim(filename))

    return

  end subroutine on_invalid_unit_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    on_input_file_does_not_exist
  !> @brief        input file does not exist error handler
  !! @authors      DS
  !! @param[in]    unit_no  : unit number of a file
  !! @param[in]    filename : file name of a file
  !! @param[in]    form     : file transfer type (formatted/unformatted)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine on_input_file_does_not_exist(unit_no, filename, form)

    ! formal arguments
    integer,      intent(in) :: unit_no
    character(*), intent(in) :: filename
    character(*), intent(in) :: form


    call error_msg('Open_file> File ' // trim(filename) // ' does not exist')

    return

  end subroutine on_input_file_does_not_exist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    on_output_file_exists
  !> @brief        output file exists error handler
  !! @authors      DS
  !! @param[in]    unit_no  : unit number of a file
  !! @param[in]    filename : file name of a file
  !! @param[in]    form     : file transfer type (formatted/unformatted)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine on_output_file_exists(unit_no, filename, form)

    ! formal arguments
    integer,      intent(in) :: unit_no
    character(*), intent(in) :: filename
    character(*), intent(in) :: form


    if (backup_allowed) then
      call backup_file(filename)
      open(unit_no, file=filename, status="new", form=form)
    else
      call error_msg('Open_file> File ' // trim(filename) // ' already exists')
    end if

    return

  end subroutine on_output_file_exists

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    on_unknown_io_error
  !> @brief        unknown I/O file error handler
  !! @authors      DS
  !! @param[in]    unit_no  : unit number of a file
  !! @param[in]    filename : file name of a file
  !! @param[in]    form     : file transfer type (formatted/unformatted)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine on_unknown_io_error(unit_no, filename, form)

    ! formal arguments
    integer,      intent(in) :: unit_no
    character(*), intent(in) :: filename
    character(*), intent(in) :: form


    write(MsgOut,'("Open_file> rank=",i0)') my_world_rank
    call error_msg('Open_file> Unknown error while opening ' // trim(filename))

  end subroutine on_unknown_io_error

end module fileio_mod
