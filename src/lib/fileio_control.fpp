!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_control
!> @brief   utilities for reading control parameters from file
!! @authors Yuji Sugita (YS), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_control_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures

  !
  type s_keyval
    character(MaxLineLong)        :: key
    character(MaxLineLong)        :: val
    type(s_keyval),       pointer :: next         => null()
  end type s_keyval

  type s_val
    character(MaxLine)            :: val
    type(s_val),          pointer :: next         => null()
  end type s_val

  type s_section
    character(MaxLineLong)        :: name
    type(s_keyval),       pointer :: keyval_head  => null()
    type(s_val),          pointer :: val_head     => null()
    type(s_section),      pointer :: next         => null()
  end type s_section

  type s_control
    type(s_section),      pointer :: section_head => null()
  end type s_control

  !
  type s_read_key
    character(MaxLineLong)        :: key
    type(s_read_key),     pointer :: next         => null()
  end type s_read_key

  type s_read_section
    logical                       :: on           = .false.
    character(MaxLineLong)        :: name
    type(s_read_key),     pointer :: key_head     => null()
  end type s_read_section

  ! variables
  integer,              parameter :: MaxControls      = 10
  type(s_control),           save :: g_controls(MaxControls)
  type(s_read_section),      save :: g_read_sections(MaxControls)

  ! subroutines
  !
  public  :: open_ctrlfile
  public  :: close_ctrlfile
  public  :: read_ctrlfile_string
  public  :: read_ctrlfile_integer
  public  :: read_ctrlfile_real
  public  :: read_ctrlfile_logical
  public  :: read_ctrlfile_type
  public  :: read_ctrlfile_section_values
  public  :: find_ctrlfile_section
  public  :: begin_ctrlfile_section
  public  :: end_ctrlfile_section
  private :: parse_ctrlfile
  private :: free_ctrlfile
  private :: print_ctrlfile
  private :: append_section
  private :: append_keyval
  private :: search_keyval
  private :: append_val
  private :: append_read_key
  private :: search_read_key
  private :: free_read_keys

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_ctrlfile
  !> @brief        open control file
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[out]   handle   : unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_ctrlfile(filename, handle)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(out)   :: handle

    ! local variables
    integer                  :: unit_no


    do handle = 1, MaxControls
      if (.not. associated(g_controls(handle)%section_head)) &
        exit
    end do

    if (handle > MaxControls) then
      handle = 0
      return
    end if

    call open_file(unit_no, filename, IOFileInput)

    call parse_ctrlfile(unit_no, g_controls(handle))

    call close_file(unit_no)

    return

  end subroutine open_ctrlfile

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_ctrlfile
  !> @brief        close control file
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_ctrlfile(handle)

    ! formal arguments
    integer,                 intent(in)    :: handle


    if (handle <= 0 .or. handle > MaxControls) &
      return

    if (.not. associated(g_controls(handle)%section_head)) &
      return

    call free_ctrlfile(g_controls(handle))

    return

  end subroutine close_ctrlfile

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_string
  !> @brief        read character value from control file
  !! @authors      YS
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] value   : value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_string(handle, section, key, value)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    character(*),            intent(inout) :: value


    if (handle <= 0 .or. handle > MaxControls) &
      return

    call append_read_key(handle, key)

    if (.not. search_keyval(g_controls(handle), section, key, value)) &
      return

    return

  end subroutine read_ctrlfile_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_integer
  !> @brief        read integer value from control file
  !! @authors      YS
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] value   : value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_integer(handle, section, key, value)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    integer,                 intent(inout) :: value

    ! local variables
    character(MaxLinelong)                 :: cval


    if (handle <= 0 .or. handle > MaxControls) &
      return

    call append_read_key(handle, key)

    if (.not. search_keyval(g_controls(handle), section, key, cval)) &
      return

    read(cval,*,err=90,end=90) value
90  return

  end subroutine read_ctrlfile_integer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_real
  !> @brief        read real value from control file
  !! @authors      YS
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] value   : value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_real(handle, section, key, value)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    real(wp),                intent(inout) :: value

    ! local variables
    character(MaxLinelong)                 :: cval


    if (handle <= 0 .or. handle > MaxControls) &
      return

    call append_read_key(handle, key)

    if (.not. search_keyval(g_controls(handle), section, key, cval)) &
      return

    read(cval,*,err=90,end=90) value
90  return

  end subroutine read_ctrlfile_real

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_logical
  !> @brief        read logical value from control file
  !! @authors      YS
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] value   : value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_logical(handle, section, key, value)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    logical,                 intent(inout) :: value

    ! local variables
    character(MaxLinelong)                 :: cval


    if (handle <= 0 .or. handle > MaxControls) &
      return

    call append_read_key(handle, key)

    if (.not. search_keyval(g_controls(handle), section, key, cval)) &
      return

    call tolower(cval)
    value = (cval == 'yes' .or. cval == 'true')

    return

  end subroutine read_ctrlfile_logical

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_type
  !> @brief        read type(integer) value from control file
  !! @authors      YS
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] value   : value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_type(handle, section, key, value, typeStrs)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    integer,                 intent(inout) :: value
    character(*),            intent(in)    :: typeStrs(:)

    ! local variables
    integer                                :: i
    character(MaxLinelong)                 :: cval, cstr1, cstr2


    if (handle <= 0 .or. handle > MaxControls) &
      return

    call append_read_key(handle, key)

    if (.not. search_keyval(g_controls(handle), section, key, cval)) &
      return

    cstr1 = cval
    call tolower(cstr1)

    do i = 1, size(typeStrs)

      cstr2 = typeStrs(i)
      call tolower(cstr2)

      if (cstr1 == cstr2) then
        value = i
        return
      end if
    end do

    call error_msg('Error in ['//trim(section)//'] '//trim(key)//&
                   ' : Unsupported string: "'//trim(cval)//'"')
    return

  end subroutine read_ctrlfile_type

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrlfile_section_values
  !> @brief        read section values from control file
  !! @authors      NT
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !! @param[inout] values  : value list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrlfile_section_values(handle, section, values)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section
    character(MaxLine),      allocatable   :: values(:)

    ! local variables
    integer                  :: count
    character(MaxLine)       :: lsection

    type(s_section), pointer :: cur_sec
    type(s_val),     pointer :: cur_v


    if (handle <= 0 .or. handle > MaxControls) &
      return

    lsection = section
    call tolower(lsection)

    cur_sec => g_controls(handle)%section_head
    do while(associated(cur_sec))
      if (cur_sec%name == lsection) &
        exit
      cur_sec => cur_sec%next
    end do
    if (.not. associated(cur_sec)) &
      return

    count = 0
    cur_v => cur_sec%val_head
    do while(associated(cur_v))
      count = count + 1
      cur_v => cur_v%next
    end do

    allocate(values(count))

    count = 0
    cur_v => cur_sec%val_head
    do while(associated(cur_v))
      count = count + 1
      values(count) = cur_v%val
      cur_v => cur_v%next
    end do

    return

  end subroutine read_ctrlfile_section_values

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      find_ctrlfile_section
  !> @brief        check whether a control file has indicated section
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[in]    section  : section string
  !! @return       has section or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function find_ctrlfile_section(filename, section)

    ! return values
    logical :: find_ctrlfile_section

    ! formal arguments
    character(*),            intent(in)    :: filename
    character(*),            intent(in)    :: section

    ! local variables
    integer                  :: unit_no, istatus
    integer                  :: idx_sh, idx_bl, idx_br
    integer                  :: ichara, isize, ic
    character(MaxLineLong)   :: line
    character(MaxLineLong)   :: sec_target, sec_control
    character(LineBuffer)    :: chara
    character(1)             :: c


    find_ctrlfile_section = .false.

    sec_target = section
    call tolower(sec_target)

    call open_file(unit_no, filename, IOFileInput)

    do while(.true.)

      line    = ''
      chara   = ''
      istatus = 0
      ichara  = 0
      do while (istatus == 0)
        read(unit_no,'(a)',advance='no',size=isize,                     &
             iostat=istatus,end=100) chara
        do ic = 1, isize
          c = chara(ic:ic)
          ichara = ichara + 1
          if (ichara > MaxLineLong) &
            call error_msg('Error : a line should be less than '//&
                           'MaxLineLong characters')
          line(ichara:ichara) = c
        end do
      end do

      idx_sh = scan(line, '#')
      idx_bl = scan(line, '[')
      idx_br = scan(line, ']')

      if (idx_bl == 0 .or. idx_br == 0 .or. idx_bl + 2 > idx_br) &
        cycle

      if (idx_sh > 0 .and. idx_sh < idx_bl) &
        cycle

      sec_control = adjustl(line(idx_bl+1:idx_br-1))
      call tolower(sec_control)

      if (sec_target == sec_control) then
        find_ctrlfile_section = .true.
        exit
      end if

    end do

100 continue
    call close_file(unit_no)

    return

  end function find_ctrlfile_section

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    begin_ctrlfile_section
  !> @brief        begin of the reading section
  !! @authors      NT
  !! @param[in]    handle  : unit number
  !! @param[in]    section : section string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine begin_ctrlfile_section(handle, section)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: section


    if (handle <= 0 .or. handle > MaxControls) &
      return

    if (g_read_sections(handle)%on) &
      return

    g_read_sections(handle)%on = .true.
    g_read_sections(handle)%name = section

    call tolower(g_read_sections(handle)%name)

    return

  end subroutine begin_ctrlfile_section

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    end_ctrlfile_section
  !> @brief        end of the reading section
  !! @authors      NT
  !! @param[in]    handle  : unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine end_ctrlfile_section(handle)

    ! formal arguments
    integer,                 intent(in)    :: handle

    ! local variables
    type(s_section), pointer :: section
    type(s_keyval),  pointer :: keyval
    character(MaxLineLong)   :: name
    logical                  :: error_stop


    if (handle <= 0 .or. handle > MaxControls) &
      return

    if (.not. g_read_sections(handle)%on) &
      return

    name = g_read_sections(handle)%name

    g_read_sections(handle)%on   = .false.
    g_read_sections(handle)%name = ''


    ! check control file section keys
    !
    error_stop = .false.

    section => g_controls(handle)%section_head

    do while(associated(section))
      if (section%name == name) then

        keyval => section%keyval_head

        do while(associated(keyval))
          if (.not. search_read_key(handle, keyval%key)) then
            if (main_rank) then
              write(MsgOut,'(a)') &
                 'End_Ctrlfile_Section> Unknown parameter: ['// &
                 trim(name)//'] "'//trim(keyval%key)//'":'
            end if
            error_stop = .true.
          end if
          keyval => keyval%next

        end do

      end if
      section => section%next
    end do


    ! clear the read keys
    !
    call free_read_keys(handle)


    if (error_stop .and. main_rank) &
      call error_msg('End_Ctrlfile_Section> Stop.')

    return

  end subroutine end_ctrlfile_section

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_ctrlfile
  !> @brief        parse control file
  !! @authors      NT
  !! @param[in]    unit_no : file unit number
  !! @param[inout] control : control information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_ctrlfile(unit_no, control)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_control),         intent(inout) :: control

    ! local variables
    integer                  :: bb, be, sharp, bslash, equal, i, ichara
    integer                  :: istatus, isize, ic
    character(MaxLineLong)   :: line_pp, val
    character(MaxLineLong)   :: line
    character(MaxLineLong)   :: key, sec
    character(LineBuffer)    :: chara
    character(1)             :: c

    type(s_section), pointer :: cur_sec


    ! create no-name section
    !

    cur_sec => append_section(control, '')

    line_pp = ''

    do while(.true.)
      line    = ''
      chara   = ''
      istatus = 0
      ichara  = 0

      do while (istatus == 0)
        read(unit_no,'(a)',advance='no',size=isize,                     &
             iostat=istatus,end=100) chara
        do ic = 1, isize
          c = chara(ic:ic)
          ichara = ichara + 1
          if (ichara > MaxLineLong)  &
            call error_msg('Error in ['//trim(cur_sec%name)//'] '//&
                           ': a line should be less than MaxLineLong'//&
                           ' characters')
          line(ichara:ichara) = c
        end do
      end do

      ! pre-processing
      !

      ! comment
      sharp = scan(line, '#')
      if (sharp /= 0) &
        line(sharp:) = ''

      if (len_trim(line_pp)+len_trim(line) > MaxLineLong) &
         call error_msg('Error in ['//trim(cur_sec%name)//'] '//&
                   ' : too long line')

      line_pp = trim(line_pp)//line

      ! back-slash
      bslash = scan(line_pp, char(92))
      if (bslash /= 0) then
        line_pp(bslash:) = ''
        cycle
      end if

      ! tab clear
      do i = 1, len(line_pp)
        if (line_pp(i:i) == char(09)) &
          line_pp(i:i) = ' '
      end do

      bb = scan(line_pp, '[')
      be = scan(line_pp, ']', .true.)


      ! read line
      !

      ! section
      if (bb /= 0 .and. be /= 0 .and. bb < be) then

        sec = line_pp(bb+1:be-1)

        call tolower(sec)

        cur_sec => append_section(control, sec)

      else

        equal = scan(line_pp, '=')

        ! key-value statement
        if (equal /= 0) then

          key = adjustl(line_pp(:equal-1))
          val = adjustl(line_pp(equal+1:))

          call tolower(key)

          call append_keyval(cur_sec, key, val)

        ! value statement
        else

          val = adjustl(line_pp)

          call append_val(cur_sec, val)

        end if

      end if

      line_pp = ''

    end do

100 continue
    !!call print_ctrlfile(control)

    return

  end subroutine parse_ctrlfile

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    free_ctrlfile
  !> @brief        free control file
  !! @authors      NT
  !! @param[inout] control : control information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine free_ctrlfile(control)

    ! formal arguments
    type(s_control),         intent(inout) :: control

    ! local variables
    type(s_section), pointer :: cur_sec, sec
    type(s_keyval),  pointer :: cur_kv, kv


    cur_sec => control%section_head
    do while(associated(cur_sec))
      sec => cur_sec

      cur_kv => cur_sec%keyval_head
      do while(associated(cur_kv))
        kv => cur_kv
        cur_kv => cur_kv%next
        deallocate(kv)
      end do

      cur_sec => cur_sec%next
      deallocate(sec)
    end do

    return

  end subroutine free_ctrlfile

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_ctrlfile
  !> @brief        print control file
  !! @authors      NT
  !! @param[inout] control : control information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_ctrlfile(control)

    ! formal arguments
    type(s_control),         intent(inout) :: control

    ! local variables
    type(s_section), pointer :: cur_sec
    type(s_keyval),  pointer :: cur_kv


    cur_sec => control%section_head
    do while(associated(cur_sec))

      write(MsgOut,*) 'SECTION: ', trim(cur_sec%name)

      cur_kv => cur_sec%keyval_head
      do while(associated(cur_kv))
        write(MsgOut,*) '  key:', trim(cur_kv%key), '  val:', trim(cur_kv%val)
        cur_kv => cur_kv%next
      end do

      cur_sec => cur_sec%next
    end do

    return

  end subroutine print_ctrlfile

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      append_section
  !> @brief        allocate and append new section
  !! @authors      NT
  !! @param[inout] control : control information
  !! @param[in]    name    : section name
  !! @return       pointer to section information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function append_section(control, name)

    ! return
    type(s_section), pointer :: append_section

    ! formal arguments
    type(s_control),          intent(inout) :: control
    character(*),             intent(in)    :: name

    ! local variables
    type(s_section), pointer  :: cur_sec, prv_sec


    prv_sec => control%section_head
    cur_sec => control%section_head

    do while(associated(cur_sec))

      if (cur_sec%name == name) then
        append_section => cur_sec
        return
      end if

      prv_sec => cur_sec
      cur_sec => cur_sec%next

    end do

    if (.not. associated(prv_sec)) then
      allocate(control%section_head)
      control%section_head%name = name
      append_section => control%section_head
    else
      allocate(cur_sec)
      cur_sec%name = name
      prv_sec%next => cur_sec
      append_section => cur_sec
    end if

    return

  end function append_section

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    append_keyval
  !> @brief        append new key-value to section
  !! @authors      NT
  !! @param[inout] section : section information
  !! @param[in]    key     : key string
  !! @param[in]    val     : value string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine append_keyval(section, key, val)

    ! formal arguments
    type(s_section),         intent(inout) :: section
    character(*),            intent(in)    :: key
    character(*),            intent(in)    :: val

    ! local variables
    type(s_keyval),  pointer :: cur_kv, prv_kv


    if (key == '') &
      return

    prv_kv => section%keyval_head
    cur_kv => section%keyval_head

    do while(associated(cur_kv))

      if (cur_kv%key == key) then
        if (main_rank)    &
          write(MsgOut,*) &
            '  WARNING : re-defined: '//trim(key)//'  previous was overwote.'
        cur_kv%val = val
        return
      end if

      prv_kv => cur_kv
      cur_kv => cur_kv%next

    end do

    if (.not. associated(prv_kv)) then
      allocate(section%keyval_head)
      section%keyval_head%key = key
      section%keyval_head%val = val
    else
      allocate(cur_kv)
      cur_kv%key = key
      cur_kv%val = val
      prv_kv%next => cur_kv
    end if

    return

  end subroutine append_keyval

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_keyval
  !> @brief        search key-value
  !! @authors      NT
  !! @param[in]    control : control information
  !! @param[in]    section : section string
  !! @param[in]    key     : key string
  !! @param[inout] val     : value string
  !! @return       flag for found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_keyval(control, section, key, val)

    ! return
    logical                  :: search_keyval

    ! formal arguments
    type(s_control),         intent(in)    :: control
    character(*),            intent(in)    :: section
    character(*),            intent(in)    :: key
    character(*),            intent(inout) :: val

    ! local variabls
    character(MaxLineLong)   :: lsection, lkey
    type(s_section), pointer :: cur_sec
    type(s_keyval),  pointer :: cur_kv


    lsection = section
    call tolower(lsection)

    lkey = key
    call tolower(lkey)

    cur_sec => control%section_head
    do while(associated(cur_sec))
      if (cur_sec%name == lsection) then
        cur_kv => cur_sec%keyval_head
        do while(associated(cur_kv))
          if (cur_kv%key == lkey) then
            val = cur_kv%val
            search_keyval = .true.
            return
          end if
          cur_kv => cur_kv%next
        end do
      end if
      cur_sec => cur_sec%next
    end do

    search_keyval = .false.
    return

  end function search_keyval

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    append_val
  !> @brief        append new value to section
  !! @authors      NT
  !! @param[inout] section : section information
  !! @param[in]    val     : value string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine append_val(section, val)

    ! formal arguments
    type(s_section),         intent(inout) :: section
    character(*),            intent(in)    :: val

    ! local variables
    type(s_val),     pointer :: cur_v, prv_v


    if (val == '') &
      return

    prv_v => section%val_head
    cur_v => section%val_head

    do while(associated(cur_v))

      prv_v => cur_v
      cur_v => cur_v%next

    end do

    if (.not. associated(prv_v)) then
      allocate(section%val_head)
      section%val_head%val = val
    else
      allocate(cur_v)
      cur_v%val = val
      prv_v%next => cur_v
    end if

    return

  end subroutine append_val

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    append_read_key
  !> @brief        append the read key name
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    key    : key string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine append_read_key(handle, key)

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: key

    ! local variables
    type(s_read_key), pointer :: cur_rk, prv_rk


    if (.not. g_read_sections(handle)%on) &
      return

    cur_rk => g_read_sections(handle)%key_head
    prv_rk => g_read_sections(handle)%key_head

    do while(associated(cur_rk))
      prv_rk => cur_rk
      cur_rk => cur_rk%next
    end do

    if (.not. associated(prv_rk)) then
      allocate(g_read_sections(handle)%key_head)
      g_read_sections(handle)%key_head%key = key

      call tolower(g_read_sections(handle)%key_head%key)
    else
      allocate(cur_rk)
      cur_rk%key = key
      prv_rk%next => cur_rk

      call tolower(cur_rk%key)
    end if

    return

  end subroutine append_read_key

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      search_read_key
  !> @brief        search the read key name
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !! @param[in]    key    : key string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function search_read_key(handle, key)

    ! return value
    logical                  :: search_read_key

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: key

    ! local variables
    type(s_read_key), pointer :: read_key


    read_key => g_read_sections(handle)%key_head

    do while(associated(read_key))
      if (read_key%key == key) then
        search_read_key = .true.
        return
      end if
      read_key => read_key%next
    end do

    search_read_key = .false.
    return

  end function search_read_key

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    free_read_keys
  !> @brief        free the read keys
  !! @authors      NT
  !! @param[in]    handle : handle for file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine free_read_keys(handle)

    ! formal arguments
    integer,                 intent(in)    :: handle

    ! local variables
    type(s_read_key), pointer :: cur_rk, nxt_rk


    cur_rk => g_read_sections(handle)%key_head

    do while(associated(cur_rk))
      nxt_rk => cur_rk%next
      deallocate(cur_rk)
      cur_rk => nxt_rk
    end do

    g_read_sections(handle)%key_head => null()

    return

  end subroutine free_read_keys

end module fileio_control_mod
