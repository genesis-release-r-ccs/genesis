!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_data
!> @brief   utilities for reading data from binary file
!! @authors Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_data_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer, public, parameter   :: IOFileDataRead       = 1
  integer, public, parameter   :: IOFileDataWrite      = 2
  integer, public, parameter   :: IOFileDataWriteAscii = 3
  integer, public, parameter   :: IOFileDataAppend     = 4

  integer(8)  :: Zero = 0   ! for gfortran 4.4.7 bug
  integer(8)  :: One = 1    ! for gfortran 4.4.7 bug

  ! structures
  type s_data
    character(40)              :: tag = ''
    character                  :: var_type = ''
    character,    allocatable  :: payload(:)
    type(s_data), pointer      :: next => null()
    integer                    :: read_pos = 1
  end type s_data

  type s_file
    character(MaxFilename)     :: filename = ''
    character(20)              :: mode = ''
    type(s_data), pointer      :: head => null()
    character,    allocatable  :: b(:)
  end type s_file

  integer,        parameter    :: MaxFile = 10
  type(s_file),   target, save :: g_file(MaxFile)

  ! variables
  integer,        parameter    :: SizeHeader     = 8
  integer,        parameter    :: SizeBlockSize  = 8
  integer,        parameter    :: SizeVarType    = 1
  integer,        parameter    :: SizeTag        = 40
  integer,        parameter    :: SizeNextPos    = 8
  integer,        parameter    :: SizeNonPayload = SizeBlockSize + &
                                                   SizeVarType   + &
                                                   SizeTag       + &
                                                   SizeNextPos

  character(8),   parameter    :: Header = 'GD150608'

  ! subroutines
  public  :: open_data
  public  :: close_data
  public  :: write_data_logical_array
  public  :: write_data_logical
  public  :: write_data_byte_array
  public  :: write_data_byte
  public  :: write_data_integer_array
  public  :: write_data_integer
  public  :: write_data_real_wp_array
  public  :: write_data_real_wp
  public  :: write_data_real_wip_array
  public  :: write_data_real_wip
  public  :: write_data_real_dp_array
  public  :: write_data_real_dp
  public  :: write_data_real_sp_array
  public  :: write_data_real_sp
  public  :: read_data_logical_array
  public  :: read_data_logical
  public  :: read_data_byte_array
  public  :: read_data_byte
  public  :: read_data_integer_array
  public  :: read_data_integer
  public  :: read_data_real_wp_array
  public  :: read_data_real_wp
  public  :: read_data_real_wip_array
  public  :: read_data_real_wip
  public  :: read_data_real_dp_array
  public  :: read_data_real_dp
  public  :: read_data_real_sp_array
  public  :: read_data_real_sp
  public  :: get_data_size
  private :: import_data
  private :: import_data_binary
  private :: import_data_ascii
  private :: export_data_binary
  private :: export_data_ascii
  private :: read_data
  private :: write_data
  private :: clear_data
  private :: append_data
  public  :: find_tag
  private :: find_data
  private :: get_data
  private :: get_size
  private :: read_block_size
  private :: read_var_type
  private :: read_tag
  private :: read_payload_rec
  private :: read_next_pos
  private :: read_data_size
  private :: read_data_size_rec
  private :: write_block_size
  private :: write_var_type
  private :: write_tag
  private :: write_payload
  private :: write_next_pos
  private :: transfer_i2b
  private :: transfer_w2b
  private :: transfer_wi2b
  private :: transfer_d2b
  private :: transfer_s2b
  private :: transfer_b2i
  private :: transfer_b2w
  private :: transfer_b2wi
  private :: transfer_b2d
  private :: transfer_b2s
  private :: swap4
  private :: swap8
  private :: fd_open
  private :: fd_close
  private :: fd_flen
  private :: fd_read
  private :: fd_write

  ! File format:
  !
  !  1: Block 1
  !  2: Block 2
  !  3: Block 3
  !       :
  !  N: Block N
  !

  ! Block format:
  !
  !     #byte  contents
  ! -------------------
  !  1:    8   block_size
  !  2:    1   var. type (b:byte, i:integer, w:real(wp), d:real(dp), s:real(sp))
  !  3:   40   tag name
  !  4:    n   payload
  !  5:    8   next common tag block pos. (zero => termiated)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_data
  !> @brief        open data file
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[in]    in_out   : flag for the file
  !! @param[out]   handle   : handle for context
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_data(filename, in_out, handle)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out
    integer,                 intent(out)   :: handle

    ! local variables
    logical                  :: bin


    do handle = 1, MaxFile
      if (g_file(handle)%filename == '') &
        exit
    end do

    if (handle > MaxFile) &
      call error_msg('Open_Data> open file reached the maximum.')

    g_file(handle)%filename = filename

    select case (in_out)

    case (IOFileDataRead)
      call import_data(g_file(handle), bin)
      g_file(handle)%mode = 'read'

    case (IOFileDataWrite)
      g_file(handle)%mode = 'write'

    case (IOFileDataWriteAscii)
      g_file(handle)%mode = 'write_ascii'

    case (IOFileDataAppend)
      call import_data(g_file(handle), bin)
      if (bin) then
        g_file(handle)%mode = 'write'
      else
        g_file(handle)%mode = 'write_ascii'
      end if

    end select

    allocate(g_file(handle)%b(1))

    return

  end subroutine open_data

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_data
  !> @brief        close data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_data(handle)

    ! formal arguments
    integer,                 intent(in)    :: handle


    if (handle <= 0 .or. handle > MaxFile) &
      return

    if (g_file(handle)%filename == '') &
      return

    select case (g_file(handle)%mode)

    case ('write')
      call export_data_binary(g_file(handle))

    case ('write_ascii')
      call export_data_ascii (g_file(handle))

    end select

    call clear_data(g_file(handle))

    deallocate(g_file(handle)%b)

    g_file(handle)%filename = ''
    g_file(handle)%head => null()

    return

  end subroutine close_data

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_logical_array
  !> @brief        write logical array to data file
  !! @authors      KY
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_logical_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    logical,                 intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds
    character, allocatable   :: v(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    allocate(v(ds))

    do i = 1, ds
      if (data(i)) then
        v(i) = 'T'
      else
        v(i) = 'F'
      end if
    end do

    call write_data(handle, 'l', tag, ds, v)
    deallocate(v)

    return

  end subroutine write_data_logical_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_logical
  !> @brief        write logical to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_logical(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    logical,                 intent(in)    :: data

    ! local variables
    character                :: v


    if (data) then
      v = 'T'
    else
      v = 'F'
    end if

    call write_data(handle, 'l', tag, 1, (/v/))

    return

  end subroutine write_data_logical

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_byte_array
  !> @brief        write byte array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_byte_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    character,               intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    call write_data(handle, 'b', tag, ds, data)

    return

  end subroutine write_data_byte_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_byte
  !> @brief        write byte to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_byte(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    character,               intent(in)    :: data


    call write_data(handle, 'b', tag, 1, (/data/))

    return

  end subroutine write_data_byte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_integer_array
  !> @brief        write integer array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_integer_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    integer,                 intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds

    character, pointer :: b(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 4))
    end if

    b => g_file(handle)%b

    call transfer_i2b(data, b, ds)

    call write_data(handle, 'i', tag, ds * 4, b)

    return

  end subroutine write_data_integer_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_integer
  !> @brief        write integer to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_integer(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: data

    ! local variables
    character, pointer :: b(:)


    if (4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(4))
    end if

    b => g_file(handle)%b

    b(1:4) = transfer(data, b(1:4))

    call write_data(handle, 'i', tag, 4, b)

    return

  end subroutine write_data_integer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_wp_array
  !> @brief        write real(wp) array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_wp_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(wp),                intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds

    character, pointer :: b(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call transfer_w2b(data, b, ds)

    call write_data(handle, 'w', tag, ds * 8, b)

    return

  end subroutine write_data_real_wp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_wp
  !> @brief        write real(wp) to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_wp(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(wp),                intent(in)    :: data

    ! local variables
    real(dp)                 :: dd

    character, pointer :: b(:)


    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    dd = real(data,dp)
    b(1:8) = transfer(dd, b(1:8))

    call write_data(handle, 'w', tag, 8, b)

    return

  end subroutine write_data_real_wp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_wip_array
  !> @brief        write real(wip) array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_wip_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(wip),               intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds

    character, pointer :: b(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call transfer_wi2b(data, b, ds)

    call write_data(handle, 'd', tag, ds * 8, b)

    return

  end subroutine write_data_real_wip_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_wip
  !> @brief        write real(wip) to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_wip(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(wip),               intent(in)    :: data

    ! local variables
    real(dp)                 :: dd

    character, pointer :: b(:)


    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    dd = real(data,dp)
    b(1:8) = transfer(dd, b(1:8))

    call write_data(handle, 'd', tag, 8, b)

    return

  end subroutine write_data_real_wip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_dp_array
  !> @brief        write real(dp) array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_dp_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(dp),                intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds

    character, pointer :: b(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call transfer_d2b(data, b, ds)

    call write_data(handle, 'd', tag, ds * 8, b)

    return

  end subroutine write_data_real_dp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_dp
  !> @brief        write real(dp) to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_dp(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(dp),                intent(in)    :: data

    ! local variables
    character, pointer :: b(:)


    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    b(1:8) = transfer(data, b(1:8))

    call write_data(handle, 'd', tag, 8, b)

    return

  end subroutine write_data_real_dp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_sp_array
  !> @brief        write real(sp) array to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    dsize  : size of data
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_sp_array(handle, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(sp),                intent(in)    :: data(*)

    ! local variables
    integer                  :: i, ds

    character, pointer :: b(:)


    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 4))
    end if

    b => g_file(handle)%b

    call transfer_s2b(data, b, ds)

    call write_data(handle, 's', tag, ds * 4, b)

    return

  end subroutine write_data_real_sp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_data_real_sp
  !> @brief        write real(sp) to data file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    tag    : tag name
  !! @param[in]    data   : data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data_real_sp(handle, tag, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(sp),                intent(in)    :: data

    ! local variables
    character, pointer :: b(:)


    if (4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(4))
    end if

    b => g_file(handle)%b

    b(1:4) = transfer(data, b(1:4))

    call write_data(handle, 's', tag, 4, b)

    return

  end subroutine write_data_real_sp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_logical_array
  !> @brief        read logical array from data file
  !! @authors      KY
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_logical_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    logical,                 intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    character, allocatable   :: v(:)
    logical                  :: no_tag, err_stop


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    allocate(v(ds))

    call read_data(handle, 'l', tag, ds, v, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Logical_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
    end if

    do i = 1, ds
      data(i) = (v(i) == 'T')
    end do

    deallocate(v)

    return

  end subroutine read_data_logical_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_logical
  !> @brief        read logical from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_logical(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    logical,                 intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    character(200)           :: msg
    character                :: v(1)
    logical                  :: no_tag, err_stop


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    call read_data(handle, 'l', tag, 1, v, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Logical> indicated tag is not found at file.['// &
            trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    data = (v(1) == 'T')

    return

  end subroutine read_data_logical

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_byte_array
  !> @brief        read byte array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_byte_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    character,               intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    call read_data(handle, 'b', tag, ds, data, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Byte_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
    end if

    return

  end subroutine read_data_byte_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_byte
  !> @brief        read byte from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_byte(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    character,               intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    character(200)           :: msg
    character                :: d(1)
    logical                  :: no_tag, err_stop


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    call read_data(handle, 'b', tag, 1, d, no_tag)
    if (no_tag) then
      msg='Read_Data_Byte> indicated tag is not found at file.['//trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    data = d(1)

    return

  end subroutine read_data_byte

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_integer_array
  !> @brief        read integer array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_integer_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    integer,                 intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 4))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'i', tag, ds * 4, b, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Integer_Array> indicated tag is not found at file.['// &
        trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    call transfer_b2i(b, data, ds)

    return

  end subroutine read_data_integer_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_integer
  !> @brief        read integer from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_integer(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    if (4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(4))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'i', tag, 4, b, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Integer> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    data = transfer(b(1:4), data)

    return

  end subroutine read_data_integer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_wp_array
  !> @brief        read real(wp) array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_wp_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(wp),                intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'w', tag, ds * 8, b, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Real_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    call transfer_b2w(b, data, ds)

    return

  end subroutine read_data_real_wp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_wp
  !> @brief        read real(wp) from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_wp(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(wp),                intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    real(dp)                 :: dd
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'w', tag, 8, b, no_tag)
    if (no_tag) then
      msg='Read_Data_Real> indicated tag is not found at file.['//trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    dd = transfer(b(1:8), dd)
    data = real(dd,wp)

    return

  end subroutine read_data_real_wp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_wip_array
  !> @brief        read real(wip) array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_wip_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(wip),               intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'd', tag, ds * 8, b, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Real_Wip_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    call transfer_b2wi(b, data, ds)

    return

  end subroutine read_data_real_wip_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_wip
  !> @brief        read real(wip) from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_wip(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(wip),               intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    real(dp)                 :: dd
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'd', tag, 8, b, no_tag)
    if (no_tag) then
  msg='Read_Data_Real_Wip> indicated tag is not found at file.['//trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    dd = transfer(b(1:8), dd)
    data = real(dd,wip)

    return

  end subroutine read_data_real_wip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_dp_array
  !> @brief        read real(dp) array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_dp_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(dp),                intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'd', tag, ds * 8, b, no_tag)
    if (no_tag) then
      msg = 'Read_Data_Real_Dp_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    call transfer_b2d(b, data, ds)

    return

  end subroutine read_data_real_dp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_dp
  !> @brief        read real(dp) from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_dp(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(dp),                intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    if (8 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(8))
    end if

    b => g_file(handle)%b

    call read_data(handle, 'd', tag, 8, b, no_tag)

    if (no_tag) then
      msg = 'Read_Data_Real_Dp> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    data = transfer(b(1:8), data)

    return

  end subroutine read_data_real_dp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_sp_array
  !> @brief        read real(sp) array from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[in]    dsize     : size of data
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_sp_array(handle, tag, dsize, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize(:)
    real(sp),                intent(inout) :: data(*)
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    integer                  :: i, ds
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    ds = 1
    do i = 1, size(dsize)
      ds = ds * dsize(i)
    end do

    if (ds * 4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(ds * 4))
    end if

    b => g_file(handle)%b

    call read_data(handle, 's', tag, ds * 4, b, no_tag)

    if (no_tag) then
      msg = 'Read_Data_Real_Sp_Array> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    call transfer_b2s(b, data, ds)

    return

  end subroutine read_data_real_sp_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_data_real_sp
  !> @brief        read real(sp) from data file
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[inout] data      : data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_real_sp(handle, tag, data, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    real(sp),                intent(inout) :: data
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    character(200)           :: msg
    logical                  :: no_tag, err_stop

    character, pointer :: b(:)


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    if (4 > size(g_file(handle)%b)) then
      deallocate(g_file(handle)%b)
      allocate(g_file(handle)%b(4))
    end if

    b => g_file(handle)%b

    call read_data(handle, 's', tag, 4, b, no_tag)

    if (no_tag) then
      msg = 'Read_Data_Real_Sp> indicated tag is not found at file.['// &
           trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      return
    end if

    data = transfer(b(1:4), data)

    return

  end subroutine read_data_real_sp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_data_size
  !> @brief        get the data size for indicated tag
  !! @authors      NT
  !! @param[in]    handle    : handle for file
  !! @param[in]    tag       : tag name
  !! @param[out]   dsize     : size of data
  !! @param[in]    oerr_stop : flag for error-stop or not when tag is not found
  !!                           (optional: default is .true.)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_data_size(handle, tag, dsize, oerr_stop)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag
    integer,                 intent(out)   :: dsize
    logical, optional,       intent(in)    :: oerr_stop

    ! local variables
    type(s_data),    pointer :: dat

    character(200)           :: msg
    logical                  :: err_stop


    err_stop = .true.
    if (present(oerr_stop)) err_stop = oerr_stop

    dat => get_data(g_file(handle), tag)

    if (.not. associated(dat)) then
      msg='Get_Data_Size> indicated tag is not found at file.['//trim(tag)//']'
      if (err_stop) call error_msg(msg)
      if (main_rank) write(MsgOut,'(a,/)') trim(msg)//'  skip..'
      dsize = 0
      return
    end if

    call get_size(dat, dsize)

    return

  end subroutine get_data_size

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine import_data(file, is_binary)

    ! formal arguments
    type(s_file),            intent(inout) :: file
    logical,                 intent(out)   :: is_binary

    ! local variables
    integer                  :: unit_no
    character(8)             :: hdr


    ! check Binary file header
    unit_no = get_unit_no()

    open(unit_no, &
         file = file%filename,    &
         status  = 'old',         &
         form    = 'unformatted', &
         access  = 'direct',      &
         recl    = 8)

    read(unit_no,rec=1) hdr(1:8)

    close(unit_no)

    call free_unit_no(unit_no)

    is_binary = (hdr == Header)

    if (is_binary) then

      call import_data_binary(file)

    else

      call import_data_ascii(file)

    end if

    return

  end subroutine import_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine import_data_binary(file)

    ! formal arguments
    type(s_file),            intent(inout) :: file

    ! local variables
    integer                  :: unit_no
    integer(8)               :: block_size, data_size, pos, flen
    character(SizeTag)       :: tag
    character                :: var_type

    type(s_data),   pointer  :: dat


    call fd_open(file%filename, 1, unit_no)
    if (unit_no == -1) &
      call error_msg('Import_Data_Binary> Open error.')

    call fd_flen(unit_no, flen)

    pos = SizeHeader+1

    do while(.true.)

      if (flen <= pos) &
        exit

      call read_block_size(unit_no, pos, block_size)
      call read_var_type  (unit_no, pos, var_type)
      call read_tag       (unit_no, pos, tag)

      if (.not. find_data(file, tag)) then

        call read_data_size(unit_no, pos, data_size)

        ! allocate
        allocate(dat)
        allocate(dat%payload(data_size))

        ! read payload
        call read_payload_rec(unit_no, pos, dat%payload)

        ! regist
        dat%tag = tag
        dat%var_type = var_type

        call append_data(file, dat)

      end if

      pos = pos + block_size

    end do

    call fd_close(unit_no)

    return

  end subroutine import_data_binary

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine import_data_ascii(file)

    ! formal arguments
    type(s_file),            intent(inout) :: file

    ! local variables
    integer                  :: unit_no, ndata, data_size, i
    character(2000)          :: line
    character(100)           :: tag
    character(1)             :: var_type

    integer,     allocatable :: idata(:)
    real(sp),    allocatable :: sdata(:)
    real(dp),    allocatable :: ddata(:)
    real(wp),    allocatable :: wdata(:)
    type(s_data),    pointer :: dat


    call open_file(unit_no, file%filename, IOFileInput)

    do while(.true.)

      read(unit_no,'(a)',end=100) line
      read(line,*) tag, var_type, ndata

      tag = tag(2:len_trim(tag)-1)

      allocate(dat)
      dat%tag = tag
      dat%var_type = var_type

      select case (var_type)

      case ('b')
        data_size = ndata
        allocate(dat%payload(data_size))
        allocate(idata(ndata))
        read(unit_no,*) (idata(i),i=1,ndata)
        do i = 1, ndata
          dat%payload(i) = char(idata(i))
        end do
        deallocate(idata)

      case ('l')
        data_size = ndata
        allocate(dat%payload(data_size))
        read(unit_no,*) (dat%payload(i),i=1,data_size)

      case ('i')
        data_size = ndata * 4
        allocate(dat%payload(data_size))
        allocate(idata(ndata))
        read(unit_no,*) (idata(i),i=1,ndata)
        call transfer_i2b(idata,dat%payload,ndata)
        deallocate(idata)

      case ('w')
        data_size = ndata * 8
        allocate(dat%payload(data_size))
        allocate(wdata(ndata))
        read(unit_no,*) (wdata(i),i=1,ndata)
        call transfer_w2b(wdata,dat%payload,ndata)
        deallocate(wdata)

      case ('d')
        data_size = ndata * 8
        allocate(dat%payload(data_size))
        allocate(ddata(ndata))
        read(unit_no,*) (ddata(i),i=1,ndata)
        call transfer_d2b(ddata,dat%payload,ndata)
        deallocate(ddata)

      case ('s')
        data_size = ndata * 4
        allocate(dat%payload(data_size))
        allocate(sdata(ndata))
        read(unit_no,*) (sdata(i),i=1,ndata)
        call transfer_s2b(sdata,dat%payload,ndata)
        deallocate(sdata)

      case default
        call error_msg('Import_Data_Ascii> Unknown variable type:'//var_type)

      end select

      call append_data(file, dat)

    end do

100 call close_file(unit_no)

    return

  end subroutine import_data_ascii

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_data_binary(file)

    ! formal arguments
    type(s_file),            intent(inout) :: file

    ! local variables
    integer                  :: unit_no
    integer(8)               :: block_size

    type(s_data),    pointer :: dat


    dat => file%head

    call fd_open(file%filename, 0, unit_no)
    if (unit_no == -1) &
      call error_msg('Export_Data_Binary> Open error.')

    call fd_write(unit_no, Header, int(SizeHeader,8))

    do while(associated(dat))

      block_size = size(dat%payload) + SizeNonPayload

      call write_block_size(unit_no, block_size)
      call write_var_type  (unit_no, dat%var_type)
      call write_tag       (unit_no, dat%tag)
      call write_payload   (unit_no, dat%payload)
      call write_next_pos  (unit_no, Zero)

      dat => dat%next

    end do

    call fd_close(unit_no)

    return

  end subroutine export_data_binary

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_data_ascii(file)

    ! formal arguments
    type(s_file),            intent(inout) :: file

    ! local variables
    integer                  :: unit_no, ndata, i

    integer,     allocatable :: idata(:)
    real(sp),    allocatable :: sdata(:)
    real(dp),    allocatable :: ddata(:)
    real(wp),    allocatable :: wdata(:)
    type(s_data),    pointer :: dat


    dat => file%head

    call open_file(unit_no, file%filename, IOFileOutputReplace)

    do while(associated(dat))

      call get_size(dat, ndata)

      write(unit_no,'("[",a,"]",1x,a,1x,i0)') &
           trim(dat%tag), dat%var_type, ndata

      select case (dat%var_type)

      case ('b')
        write(unit_no,'(10(i3,1x))') (ichar(dat%payload(i)),i=1,ndata)

      case ('l')
        write(unit_no,'(10(a,1x))') (dat%payload(i),i=1,ndata)

      case ('i')
        allocate(idata(ndata))
        call transfer_b2i(dat%payload, idata, ndata)
        write(unit_no,'(6(i12,1x))') (idata(i),i=1,ndata)
        deallocate(idata)

      case ('w')
        allocate(wdata(ndata))
        call transfer_b2w(dat%payload, wdata, ndata)
        write(unit_no,'(3(f26.13,1x))') (wdata(i),i=1,ndata)
        deallocate(wdata)

      case ('d')
        allocate(ddata(ndata))
        call transfer_b2d(dat%payload, ddata, ndata)
        write(unit_no,'(3(f26.13,1x))') (ddata(i),i=1,ndata)
        deallocate(ddata)

      case ('s')
        allocate(sdata(ndata))
        call transfer_b2s(dat%payload, sdata, ndata)
        write(unit_no,'(3(f26.13,1x))') (sdata(i),i=1,ndata)
        deallocate(sdata)

      end select

      dat => dat%next

    end do

    call close_file(unit_no)

    return

  end subroutine export_data_ascii

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data(handle, var_type, tag, dsize, data, no_tag)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character,               intent(in)    :: var_type
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize
    character,               intent(inout) :: data(*)
    logical,                 intent(inout) :: no_tag

    ! local variables
    integer                  :: i, ii, n

    type(s_data),    pointer :: dat


    if (handle <= 0 .or. handle > MaxFile) &
      call error_msg('Read_Data> Invalid handle.')

    if (dsize == 0) then
      no_tag = .false.
      return
    end if

    dat => get_data(g_file(handle), tag)

    if (.not. associated(dat)) then
      no_tag = .true.
      return
    else
      no_tag = .false.
    end if

    n = size(dat%payload)

    if (dat%read_pos + dsize - 1 > n) &
      call error_msg('Read_Data> there is no readable data.')

    ii = 1
    do i = dat%read_pos,dat%read_pos + dsize - 1
      data(ii) = dat%payload(i)
      ii = ii + 1
    end do

    dat%read_pos = dat%read_pos + dsize

    return

  end subroutine read_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_data(handle, var_type, tag, dsize, data)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character,               intent(in)    :: var_type
    character(*),            intent(in)    :: tag
    integer,                 intent(in)    :: dsize
    character,               intent(in)    :: data(*)

    ! local variables
    integer                  :: i, n_old, n_new

    character, allocatable   :: payload(:)
    type(s_data),    pointer :: dat


    if (handle <= 0 .or. handle > MaxFile) &
      return

    dat => get_data(g_file(handle), tag)

    if (.not. associated(dat)) then

      allocate(dat)
      allocate(dat%payload(0))

      dat%tag = tag
      dat%var_type = var_type

      call append_data(g_file(handle), dat)

    else if (dat%var_type /= var_type) then

      call error_msg('Write_Data> different variable type.')

    end if

    n_old = size(dat%payload)

    allocate(payload(n_old))

    do i = 1, n_old
      payload(i) = dat%payload(i)
    end do

    n_new = n_old + dsize

    deallocate(dat%payload)
    allocate(dat%payload(n_new))

    do i = 1, n_old
      dat%payload(i) = payload(i)
    end do

    do i = n_old + 1, n_new
      dat%payload(i) = data(i-n_old)
    end do

    deallocate(payload)

    return

  end subroutine write_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine clear_data(file)

    ! formal arguments
    type(s_file),            intent(inout) :: file

    ! local variables
    type(s_data),    pointer :: d, dd


    d => file%head

    do while(.true.)

      if (.not. associated(d)) &
        exit

      dd => d%next
      deallocate(d%payload)
      deallocate(d)
      d => dd

    end do

    return

  end subroutine clear_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine append_data(file, dat)

    ! formal arguments
    type(s_file),            intent(inout) :: file
    type(s_data),            pointer       :: dat

    ! local variables
    type(s_data),    pointer :: d, dp


    d => file%head
    dp => null()

    do while(associated(d))
      dp => d
      d => d%next
    end do

    if (.not. associated(dp)) then
      file%head => dat
    else
      dp%next => dat
    end if

    return

  end subroutine append_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  function find_tag(handle, tag)

    ! return value
    logical :: find_tag

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: tag

    ! local variables
    type(s_data),            pointer :: d


    if (handle <= 0 .or. handle > MaxFile) &
      call error_msg('Find_Tag> Invalid handle.')

    d => g_file(handle)%head
    do while(associated(d))
      if (d%tag == tag) then
        find_tag = .true.
        return
      end if
      d => d%next
    end do

    find_tag = .false.
    return

  end function find_tag

  !======1=========2=========3=========4=========5=========6=========7=========8

  function find_data(file, tag)

    ! return value
    logical :: find_data

    ! formal arguments
    type(s_file),            intent(in)    :: file
    character(*),            intent(in)    :: tag

    ! local variables
    type(s_data),    pointer :: d


    d => file%head
    do while(associated(d))
      if (d%tag == tag) then
        find_data = .true.
        return
      end if
      d => d%next
    end do

    find_data = .false.
    return

  end function find_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_data(file, tag)

    ! return value
    type(s_data), pointer :: get_data

    ! formal arguments
    type(s_file),            intent(in)    :: file
    character(*),            intent(in)    :: tag

    ! local variables
    type(s_data),    pointer :: d


    d => file%head
    do while(associated(d))
      if (d%tag == tag) then
        get_data => d
        return
      end if
      d => d%next
    end do

    get_data => null()
    return

  end function get_data

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_size(data, dsize)

    ! formal arguments
    type(s_data),            intent(in)    :: data
    integer,                 intent(out)   :: dsize


    dsize = size(data%payload)

    if (dsize > 0) then
      select case (data%var_type)
      case ('b','l')
        dsize = dsize / 1
      case ('i')
        dsize = dsize / 4
      case ('w')
        dsize = dsize / 8
      case ('d')
        dsize = dsize / 8
      case ('s')
        dsize = dsize / 4
      case default
        ! nothing to do.
      end select
    end if

    return

  end subroutine get_size

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_block_size(file, pos, block_size)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    integer(8),              intent(inout) :: block_size

    ! local variables
    integer(8)               :: i
    character                :: b(SizeBlockSize)


    do i = 0, SizeBlockSize - 1
      call fd_read(file, pos+i, b(i+1), One)
    end do

    block_size = transfer(b, block_size)
    return

  end subroutine read_block_size

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_var_type(file, pos, var_type)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    character,               intent(inout) :: var_type

    ! local variables
    integer(8)               :: pos0


    pos0 = pos + SizeBlockSize

    call fd_read(file, pos0, var_type, One)

    return

  end subroutine read_var_type

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_tag(file, pos, tag)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    character(*),            intent(inout) :: tag

    ! local variables
    integer(8)               :: pos0


    pos0 = pos + SizeBlockSize + SizeVarType

    call fd_read(file, pos0, tag(1:SizeTag), int(SizeTag,8))

    return

  end subroutine read_tag

  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine read_payload_rec(file, pos, bytes)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    character,               intent(inout) :: bytes(:)

    ! local variables
    integer(8)               :: block_size, next_pos, payload_size, pos0


    if (pos == 0) &
      return

    call read_block_size(file, pos, block_size)
    call read_next_pos  (file, pos, block_size, next_pos)

    payload_size = block_size - SizeNonPayload

    pos0 = pos + SizeBlockSize + SizeVarType + SizeTag

    call fd_read(file, pos0, bytes(1:payload_size), payload_size)

    call read_payload_rec(file, next_pos, bytes(payload_size+1:))

    return

  end subroutine read_payload_rec

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_next_pos(file, pos, block_size, next_pos)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    integer(8),              intent(in)    :: block_size
    integer(8),              intent(inout) :: next_pos

    ! local variables
    integer(8)               :: i, pos0
    character                :: b(SizeNextPos)


    pos0 = pos + block_size - SizeNextPos

    do i = 0, SizeNextPos - 1
      call fd_read(file, pos0+i, b(i+1), One)
    end do

    next_pos = transfer(b, next_pos)

    return

  end subroutine read_next_pos

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_data_size(file, pos, data_size)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    integer(8),              intent(inout) :: data_size


    data_size = 0

    call read_data_size_rec(file, pos, data_size)

    return

  end subroutine read_data_size

  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine read_data_size_rec(file, pos, data_size)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: pos
    integer(8),              intent(inout) :: data_size

    ! local variables
    integer(8)               :: pos2
    integer(8)               :: block_size


    if (pos == 0) &
      return

    call read_block_size(file, pos, block_size)

    call read_next_pos(file, pos, block_size, pos2)

    data_size = data_size + block_size - SizeNonPayload

    call read_data_size_rec(file, pos2, data_size)

    return

  end subroutine read_data_size_rec

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_block_size(file, block_size)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: block_size

    ! local variables
    integer(8)               :: i
    character                :: b(SizeBlockSize)


    b = transfer(block_size, b)

    do i = 0, SizeBlockSize - 1
      call fd_write(file, b(i+1), One)
    end do

    return

  end subroutine write_block_size

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_var_type(file, var_type)

    ! formal arguments
    integer,                 intent(in)    :: file
    character,               intent(in)    :: var_type


    call fd_write(file, var_type, One)

    return

  end subroutine write_var_type

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_tag(file, tag)

    ! formal arguments
    integer,                 intent(in)    :: file
    character(*),            intent(in)    :: tag

    ! local variables
    character(SizeTag)       :: tagl


    tagl = tag

    call fd_write(file, tagl(1:SizeTag), int(SizeTag,8))

    return

  end subroutine write_tag

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_payload(file, payload)

    ! formal arguments
    integer,                 intent(in)    :: file
    character,               intent(in)    :: payload(:)

    ! local variables
    integer(8)               :: payload_size


    payload_size = size(payload)

    call fd_write(file, payload, payload_size)

    return

  end subroutine write_payload

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_next_pos(file, next_pos)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer(8),              intent(in)    :: next_pos

    ! local variables
    integer(8)               :: i
    character                :: b(SizeNextPos)


    b = transfer(next_pos, b)

    do i = 0, SizeNextPos - 1
      call fd_write(file, b(i+1), One)
    end do

    return

  end subroutine write_next_pos

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_i2b(idata, bdata, ndata)

    ! formal arguments
    integer,                 intent(in)    :: idata(*)
    character,               intent(inout) :: bdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      bdata(ii:ii+3) = transfer(idata(i), bdata(ii:ii+3))
      ii = ii + 4
    end do

    return

  end subroutine transfer_i2b

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_w2b(wdata, bdata, ndata)

    ! formal arguments
    real(wp),                intent(in)    :: wdata(*)
    character,               intent(inout) :: bdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    real(dp)                 :: dd
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      dd = real(wdata(i),dp)
      bdata(ii:ii+7) = transfer(dd, bdata(ii:ii+7))
      ii = ii + 8
    end do

    return

  end subroutine transfer_w2b

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_wi2b(wdata, bdata, ndata)

    ! formal arguments
    real(wip),               intent(in)    :: wdata(*)
    character,               intent(inout) :: bdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    real(dp)                 :: dd
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      dd = real(wdata(i),dp)
      bdata(ii:ii+7) = transfer(dd, bdata(ii:ii+7))
      ii = ii + 8
    end do

    return

  end subroutine transfer_wi2b

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_d2b(ddata, bdata, ndata)

    ! formal arguments
    real(dp),                intent(in)    :: ddata(*)
    character,               intent(inout) :: bdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      bdata(ii:ii+7) = transfer(ddata(i), bdata(ii:ii+7))
      ii = ii + 8
    end do

    return

  end subroutine transfer_d2b

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_s2b(sdata, bdata, ndata)

    ! formal arguments
    real(sp),                intent(in)    :: sdata(*)
    character,               intent(inout) :: bdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      bdata(ii:ii+3) = transfer(sdata(i), bdata(ii:ii+3))
      ii = ii + 4
    end do

    return

  end subroutine transfer_s2b

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_b2i(bdata, idata, ndata)

    ! formal arguments
    character,               intent(inout) :: bdata(*)
    integer,                 intent(inout) :: idata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      idata(i) = transfer(bdata(ii:ii+3), idata(i))
      ii = ii + 4
    end do

    return

  end subroutine transfer_b2i

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_b2w(bdata, wdata, ndata)

    ! formal arguments
    character,               intent(inout) :: bdata(*)
    real(wp),                intent(inout) :: wdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    real(dp)                 :: dd
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      dd = transfer(bdata(ii:ii+7), dd)
      wdata(i) = real(dd,wp)
      ii = ii + 8
    end do

    return

  end subroutine transfer_b2w

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_b2wi(bdata, wdata, ndata)

    ! formal arguments
    character,               intent(inout) :: bdata(*)
    real(wip),               intent(inout) :: wdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    real(dp)                 :: dd
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      dd = transfer(bdata(ii:ii+7), dd)
      wdata(i) = real(dd,wip)
      ii = ii + 8
    end do

    return

  end subroutine transfer_b2wi

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_b2d(bdata, ddata, ndata)

    ! formal arguments
    character,               intent(inout) :: bdata(*)
    real(dp),                intent(inout) :: ddata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      ddata(i) = transfer(bdata(ii:ii+7), ddata(i))
      ii = ii + 8
    end do

    return

  end subroutine transfer_b2d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_b2s(bdata, sdata, ndata)

    ! formal arguments
    character,               intent(inout) :: bdata(*)
    real(sp),                intent(inout) :: sdata(*)
    integer,                 intent(in)    :: ndata

    ! local variables
    integer                  :: i, ii


    ii = 1
    do i = 1, ndata
      sdata(i) = transfer(bdata(ii:ii+3), sdata(i))
      ii = ii + 4
    end do

    return

  end subroutine transfer_b2s

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine swap4(b)

    ! formal arguments
    character,               intent(inout) :: b(4)

    ! local variables
    character                :: tmp


    tmp = b(1)
    b(1) = b(4)
    b(4) = tmp

    tmp = b(2)
    b(2) = b(3)
    b(3) = tmp

    return

  end subroutine swap4

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine swap8(b)

    ! formal arguments
    character,               intent(inout) :: b(8)

    ! local variables
    character                :: tmp


    tmp = b(1)
    b(1) = b(8)
    b(8) = tmp

    tmp = b(2)
    b(2) = b(7)
    b(7) = tmp

    tmp = b(3)
    b(3) = b(6)
    b(6) = tmp

    tmp = b(4)
    b(4) = b(5)
    b(5) = tmp

    return

  end subroutine swap8

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fd_open(filename, read_only, unit_no)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: read_only
    integer,                 intent(inout) :: unit_no


#ifndef KCOMP
    unit_no = get_unit_no()
    if (read_only == 1) then
      open(unit_no, &
           file     = filename,      &
           status   = 'old',         &
           form     = 'unformatted', &
           access   = 'stream')
    else
      open(unit_no, &
           file     = filename,      &
           status   = 'replace',     &
           form     = 'unformatted', &
           access   = 'stream')
    end if
#else
    call fd_open_(filename, read_only, unit_no)
#endif

    return

  end subroutine fd_open

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fd_close(unit_no)

    ! formal arguments
    integer,                 intent(in)    :: unit_no


#ifndef KCOMP
    close(unit_no)
    call free_unit_no(unit_no)
#else
    call fd_close_(unit_no)
#endif

    return

  end subroutine fd_close

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fd_flen(unit_no, flen)

    ! formal arguments
    integer,                 intent(inout) :: unit_no
    integer(8),              intent(inout) :: flen


#ifndef KCOMP
    integer(8)               :: ftell
    call fseek(unit_no, 0, 2)
    flen = ftell(unit_no)
#else
    !!inquire(unit_no,flen=flen)
    call fd_flen_(unit_no, flen)
#endif

    return

  end subroutine fd_flen

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fd_read(unit_no, pos, b, blen)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    integer(8),              intent(in)    :: pos
    character,               intent(inout) :: b(*)
    integer(8),              intent(in)    :: blen


#ifndef KCOMP
    read(unit_no,pos=pos) b(1:blen)
#else
    call fd_read_(unit_no, pos, b, blen)
#endif

    return

  end subroutine fd_read

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fd_write(unit_no, b, blen)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    character,               intent(in)    :: b(*)
    integer(8),              intent(in)    :: blen


#ifndef KCOMP
    write(unit_no) b(1:blen)
#else
    call fd_write_(unit_no, b, blen)
#endif

    return

  end subroutine fd_write

end module fileio_data_mod
