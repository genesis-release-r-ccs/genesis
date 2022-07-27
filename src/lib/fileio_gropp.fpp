!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_gropp_mod
!> @brief   GROMACS topology file preprocessor
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_gropp_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! parameters
  integer, parameter                 :: MaxGroPPFile = 10
  integer, parameter                 :: MaxIncDepth  = 30
  integer, parameter                 :: MaxNstDepth  = 20

  ! structures
  type s_gro_macro
    character(50)                   :: def = ''
    character(MaxLine)              :: val = ''
    integer                         :: def_len = 0
    integer                         :: val_len = 0
    type(s_gro_macro), pointer      :: next => null()
  end type s_gro_macro

  type s_gropp_file
    logical                         :: opened             = .false.
    integer                         :: depth              = 0
    integer                         :: units(MaxIncDepth) = InvalidUnitNo
    character(MaxLine)              :: paths(MaxIncDepth) = '.'
    logical                         :: skips(MaxIncDepth,MaxNstDepth) = .false.
    integer                         :: nstNo(MaxIncDepth) = 1

    integer                         :: restorable_depth   = 0
    integer                         :: rest_reads(MaxIncDepth) = 0
    integer                         :: rest_units(MaxIncDepth) = 0
    logical                         :: rest_flags(MaxIncDepth) = .false.

    type(s_gro_macro), pointer      :: macro_head => null()
    type(s_gro_macro), pointer      :: macro_cur => null()
  end type s_gropp_file

  ! variables
  type(s_gropp_file),  private, save, target :: g_gropp_file(MaxGroPPFile)

  ! subroutines
  public  :: gro_pp_open_file
  public  :: gro_pp_close_file
  public  :: gro_pp_next_line
  public  :: gro_pp_next_line_s
  public  :: gro_pp_record_pos
  public  :: gro_pp_restore_pos
  private :: pp_ifdef
  private :: pp_else
  private :: pp_endif
  private :: dealloc_macro
  private :: predef_macro
  private :: expand_macro
  private :: open_include
  private :: check_skip
  private :: read_line
  private :: is_delimit

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      gro_pp_open_file
  !> @brief        open gromacs file
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[in]    predefs  : pre-defined macros
  !! @param[inout] error    : error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function gro_pp_open_file(filename, predefs, error)

    ! return value
    integer :: gro_pp_open_file

    ! formal arguments
    character(*),            intent(in)    :: filename
    character(*),            intent(in)    :: predefs(:)
    character(*),            intent(inout) :: error

    ! local variables
    integer                  :: file_no, alloc_stat
    type(s_gropp_file), pointer :: file


    ! check file handle
    !
    do file_no = 1, MaxGroPPFile
      if (.not. g_gropp_file(file_no)%opened) then
        exit
      end if
    end do
    if (file_no > MaxGroPPFile) then
      error = 'file open error. [open file > MaxFile]'
      return
    end if

    file => g_gropp_file(file_no)
    file%units(1:MaxIncDepth) = InvalidUnitNo
    file%paths(1:MaxIncDepth) = '.'
    file%skips(1:MaxIncDepth,1:MaxNstDepth) = .false.
    file%nstNo(1:MaxIncDepth) = 1

    file%restorable_depth = 0
    file%rest_reads(1:MaxIncDepth) = 0
    file%rest_units(1:MaxIncDepth) = InvalidUnitNo
    file%rest_flags(1:MaxIncDepth) = .false.


    ! allocate head macro
    !
    allocate(file%macro_head, stat = alloc_stat)
    if (alloc_stat /= 0) then
      error = 'memory allocation error.'
      goto 900
    end if

    file%macro_cur => file%macro_head
    file%macro_cur%def = 'M_A_C_R_O_H_E_A_D'


    ! allocate pre-defined macro
    !
    call predef_macro(predefs, file%macro_cur)


    ! open main file
    !
    call open_include(file, filename)
    if (file%units(1) == InvalidUnitNo) then
        error = 'file open error. ['//trim(filename)//']'
        goto 900
    end if

    file%opened = .true.

    gro_pp_open_file = file_no
    return

900 call dealloc_macro(file%macro_head)
    gro_pp_open_file = 0
    return

  end function gro_pp_open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gro_pp_close_file
  !> @brief        close gromacs file
  !! @authors      NT
  !! @param[in]    file_no : file serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gro_pp_close_file(file_no)

    ! formal arguments
    integer,                 intent(in)    :: file_no

    ! local variables
    integer                  :: i
    type(s_gropp_file), pointer :: file


    if (file_no < 1 .or. file_no > MaxGroPPFile) &
      return

    file => g_gropp_file(file_no)

    if (.not. file%opened) &
      return

    call dealloc_macro(file%macro_head)
    file%macro_head => null()
    file%macro_cur => null()

    do i = 1, MaxIncDepth
      if (file%units(i) > InvalidUnitNo) &
        call close_file(file%units(i))
      if (file%rest_units(i) > InvalidUnitNo) &
        call close_file(file%rest_units(i))
    end do

    file%rest_flags(1:MaxIncDepth) = .false.
    file%rest_units(1:MaxIncDepth) = InvalidUnitNo
    file%rest_reads(1:MaxIncDepth) = 0
    file%restorable_depth = 0

    file%nstNo(1:MaxIncDepth) = 1
    file%skips(1:MaxIncDepth,1:MaxNstDepth) = .false.
    file%paths(1:MaxIncDepth) = ''
    file%units(1:MaxIncDepth) = InvalidUnitNo
    file%depth = 0
    file%opened = .false.

    return

  end subroutine gro_pp_close_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      gro_pp_next_line
  !> @brief        read line from gromacs file
  !! @authors      NT
  !! @param[in]    file_no : file serial number
  !! @param[inout] line    : line strings
  !! @param[inout] error   : error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function gro_pp_next_line(file_no, line, error)

    ! return value
    logical                  :: gro_pp_next_line

    ! formal arguments
    integer,                 intent(in)    :: file_no
    character(*),            intent(inout) :: line
    character(*),            intent(inout) :: error

    ! local variables
    integer                  :: alloc_stat, depth
    character(MaxLine)       :: str, str1, str2, str3, str4
    logical                  :: skip

    type(s_gropp_file), pointer :: file
    type(s_gro_macro),  pointer :: macro, macro2


    gro_pp_next_line = .false.
    error = "end of file."

    if (file_no < 1 .or. file_no > MaxGroPPFile) then
      error = "bad file no."
      return
    end if

    file => g_gropp_file(file_no)
    if (.not. file%opened) then
      error = "file is not opened."
      return
    end if

    if (file%depth < 1) then
      return
    end if

    do while(.true.)

      ! read line without comment
      !
100   if (.not. read_line(file, str)) then

        depth = file%depth
        if (file%nstNo(depth) /= 1) &
          goto 900

        if (file%rest_flags(depth)) then

          if (file%rest_units(depth) /= InvalidUnitNo) &
            call error_msg('Gro_pp_next_line> Unexpected error')
          file%rest_units(depth) = file%units(depth)
          file%rest_flags(depth) = .false.
          file%units(depth) = InvalidUnitNo

        else

          call close_file(file%units(depth))
          file%units(depth) = InvalidUnitNo

        end if

        file%depth = file%depth - 1
        if (file%depth == 0) then
          return
        else
          cycle
        end if

      end if

      str1 = ''
      read(str,*,err=100,end=100) str1

      skip = check_skip(file)

      ! process #ifdef
      !
      if (str1 == '#ifdef') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call pp_ifdef(file, .true., str3)
        cycle


      ! process #ifndef
      !
      else if (str1 == '#ifndef') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call pp_ifdef(file, .false., str3)
        cycle


      ! process #else
      !
      else if (str1 == '#else') then

        call pp_else(file)
        cycle


      ! process #endif
      !
      else if (str1 == '#endif') then

        call pp_endif(file)
        cycle


      ! process #include
      !
      else if (.not. skip .and. str1 == '#include') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call open_include(file, str3)
        if (file%units(file%depth) == InvalidUnitNo) &
          goto 920

        cycle


      ! process #define
      !
      else if (.not. skip .and. str1 == '#define') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3
        str4 = adjustl(str((index(str, trim(str3)) + len_trim(str3)):))

        macro => file%macro_head
        do while(associated(macro))
          if (str3 == macro%def) then
            exit
          end if
          macro => macro%next
        end do

        ! macro expansion in macro
        macro2 => file%macro_head
        do while(associated(macro2))
          call expand_macro(macro2, str4)
          macro2 => macro2%next
        end do

        if (associated(macro)) then

          if (main_rank) &
            write(MsgOut,*) 'Fileio_Gropp> WARNING: macro:'//trim(str3)// &
                 ' re-defined.'

          macro%val = str4
          macro%val_len = len_trim(str4)

        else

          allocate(macro, stat = alloc_stat)
          if (alloc_stat /= 0) &
            goto 910

          macro%def = str3
          macro%val = str4
          macro%def_len = len_trim(str3)
          macro%val_len = len_trim(str4)
          file%macro_cur%next => macro
          file%macro_cur => macro

        end if


      ! unknown pre-processor command
      !
      else if (.not. skip .and. str1(1:1) == '#') then

        if (main_rank) &
          write(MsgOut,*) &
               'Fileio_Gropp> WARNING: unknown pre-process command ['// &
               trim(str1)//']'


      ! process non pre-processor line
      !
      else if (.not. skip) then

        ! macro expansion
        macro => file%macro_head
        do while(associated(macro))
          call expand_macro(macro, str)
          macro => macro%next
        end do

        line = str
        gro_pp_next_line = .true.
        error = ''
        return

      end if

    end do

    return

900 error = 'format error.'
    return
910 error = 'memory allocation error.'
    return
920 error = 'file open error. include file['//trim(str3)//']'
    return

  end function gro_pp_next_line

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gro_pp_next_line_s
  !> @brief        read line from gromacs file
  !! @authors      NT
  !! @param[in]    file_no : file serial number
  !! @param[inout] line    : line strings
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gro_pp_next_line_s(file_no, line)

    ! formal arguments
    integer,                 intent(in)    :: file_no
    character(*),            intent(inout) :: line

    ! local variables
    character(100)           :: error


    if (.not. gro_pp_next_line(file_no, line, error)) &
       call error_msg('Gro_pp_next_line_s> '//error)

    return

  end subroutine gro_pp_next_line_s

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gro_pp_record_pos
  !> @brief        record current position for backspace
  !! @authors      NT
  !! @param[in]    file_no : file serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gro_pp_record_pos(file_no)

    ! formal arguments
    integer,                 intent(in)    :: file_no

    ! local variables
    integer                  :: i
    type(s_gropp_file), pointer :: file


    if (file_no < 1 .or. file_no > MaxGroPPFile) &
      return

    file => g_gropp_file(file_no)
    if (.not. file%opened) then
      return
    end if

    if (file%restorable_depth > 0) &
      call error_msg('Gro_pp_record_pos> Position recording is not possible repeatedly.')

    file%restorable_depth = file%depth
    file%rest_reads(1:MaxIncDepth) = 0
    file%rest_units(1:MaxIncDepth) = InvalidUnitNo
    file%rest_flags(1:file%depth) = .true.

    return

  end subroutine gro_pp_record_pos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gro_pp_restore_pos
  !> @brief        restore recorded position for backspace
  !! @authors      NT
  !! @param[in]    file_no : file serial number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gro_pp_restore_pos(file_no)

    ! formal arguments
    integer,                 intent(in)    :: file_no

    ! local variables
    integer                  :: i, j
    type(s_gropp_file), pointer :: file


    if (file_no < 1 .or. file_no > MaxGroPPFile) &
      return

    file => g_gropp_file(file_no)
    if (.not. file%opened) then
      return
    end if

    if (file%restorable_depth == 0) &
      call error_msg('Gro_pp_restore_pos> There is no position to restore.')

    if (file%depth == 0) &
       file%depth = 1

    do i = 1, file%restorable_depth
      if (file%rest_reads(i) == 0) &
        cycle
      if (file%rest_units(i) > InvalidUnitNo) then
        if (file%units(i) > InvalidUnitNo) then
          call close_file(file%units(i))
        end if
        file%units(i) = file%rest_units(i)
        file%rest_units(i) = InvalidUnitNo
      end if
      do j = 1, file%rest_reads(i)
        backspace(file%units(i))
      end do
    end do

    file%depth = file%restorable_depth

    do i = file%depth + 1, MaxIncDepth
      call close_file(file%units(i))
      file%units(i) = InvalidUnitNo
    end do

    file%restorable_depth = 0
    file%rest_reads(1:MaxIncDepth) = 0
    file%rest_units(1:MaxIncDepth) = 0
    file%rest_flags(1:MaxIncDepth) = .false.

    return

  end subroutine gro_pp_restore_pos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pp_ifdef
  !> @brief        preprocess ifdef
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pp_ifdef(file, positive, def_str)

    ! formal arguments
    type(s_gropp_file),         intent(inout) :: file
    logical,                    intent(in)    :: positive
    character(*),               intent(in)    :: def_str

    ! local variables
    integer                     :: nst_no
    logical                     :: defined, is_true
    type(s_gro_macro), pointer  :: macro


    ! check macro
    !
    defined = .false.

    macro => file%macro_head
    do while(associated(macro))
      if (def_str == macro%def) then
        defined = .true.
        exit
      end if
      macro => macro%next
    end do

    is_true = (      positive .and.       defined .or. &
               .not. positive .and. .not. defined)

    file%nstNo(file%depth) = file%nstNo(file%depth) + 1
    nst_no = file%nstNo(file%depth)
    if (0 < nst_no .and. nst_no <= MaxNstDepth) &
      file%skips(file%depth, nst_no) = .not. is_true
    return

  end subroutine pp_ifdef

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pp_else
  !> @brief        preprocess else
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pp_else(file)

    ! formal arguments
    type(s_gropp_file),         intent(inout) :: file

    ! local variables
    integer                     :: nst_no


    nst_no = file%nstNo(file%depth)
    if (0 < nst_no .and. nst_no <= MaxNstDepth) &
      file%skips(file%depth, nst_no) = .not. file%skips(file%depth, nst_no)

    return

  end subroutine pp_else

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pp_endif
  !> @brief        preprocess endif
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pp_endif(file)

    ! formal arguments
    type(s_gropp_file),         intent(inout) :: file

    ! local variables
    integer                     :: nst_no


    nst_no = file%nstNo(file%depth)
    if (0 < nst_no .and. nst_no <= MaxNstDepth) &
      file%skips(file%depth, nst_no) = .false.

    file%nstNo(file%depth) = file%nstNo(file%depth) - 1

    return

  end subroutine pp_endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_macro
  !> @brief        deallocate macro list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine dealloc_macro(macro_cur)

    ! formal arguments
    type(s_gro_macro), pointer :: macro_cur

    ! local variables
    integer            :: dealloc_stat


    if (.not. associated(macro_cur)) &
      return

    call dealloc_macro(macro_cur%next)

    deallocate(macro_cur, stat = dealloc_stat)

    return

  end subroutine dealloc_macro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    predef_macro
  !> @brief        setup pre-defined macro
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine predef_macro(predefs, macro_cur)

    ! formal arguments
    character(*),               intent(in)    :: predefs(:)
    type(s_gro_macro), pointer, intent(inout) :: macro_cur

    ! local variables
    integer                     :: i, ic, alloc_stat
    character(MaxLine)          :: str1, str2

    type(s_gro_macro), pointer  :: macro


!    do i = 1, size(predefs)
    do i = 2, size(predefs)

      if (predefs(i) == '') &
        cycle

      ic = scan(predefs(i),':')
      if (ic == 0) then
        str1 = predefs(i)
        str2 = ''
      else
        str1 = predefs(i)(:ic-1)
        str2 = predefs(i)(ic+1:)
      end if

      allocate(macro, stat = alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg_alloc

      macro%def = str1
      macro%val = str2
      macro%def_len = len_trim(str1)
      macro%val_len = len_trim(str2)
      macro_cur%next => macro
      macro_cur => macro

    end do

    return

  end subroutine predef_macro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    expand_macro
  !> @brief        expand macro
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine expand_macro(macro, str)

    ! formal arguments
    type(s_gro_macro), pointer, intent(in)    :: macro
    character(*),               intent(inout) :: str

    ! local variables
    integer                     :: idxm, idxb, idxe, len
    character(MaxLine)          :: str2, str3


    if (index(str, macro%def(1:1)) == 0) &
       return

    idxm = index(str, trim(macro%def))
    if (idxm == 0) &
      return

    str2 = ''

100 if (idxm /= 0) then

      len = len_trim(str)
      idxb = idxm - 1
      idxe = idxm + macro%def_len

      str3 = str(1:idxb)
      str  = str(idxe:)

      if (is_delimit(str3(idxb:idxb)) .and. is_delimit(str(1:1))) then
        str2 = trim(str2)//str3(1:idxb)//macro%val(1:macro%val_len)
      else
        str2 = trim(str2)//str3(1:idxb)//macro%def(1:macro%def_len)
      end if

      idxm = index(str, trim(macro%def))
      goto 100

    end if
    str = trim(str2)//str

    return

  end subroutine expand_macro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_include
  !> @brief        open include file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_include(file, filename)

    ! formal arguments
    type(s_gropp_file),      intent(inout) :: file
    character(*),            intent(in)    :: filename

    ! local variables
    integer                  :: i, depth, idx
    character(MaxLine)       :: path, path0


    depth = file%depth
    if (depth + 1 > MaxIncDepth) then
      file%units(depth) = InvalidUnitNo
      return
    end if

    depth = depth + 1

    if (filename(1:1) == '/') then
      call open_file(file%units(depth), filename, IOFileInput)

    else
      path = ""
      do i = 1, depth - 1
        path0 = file%paths(i)
        if (path0(1:1) == '/' .or. i == 1) then
          path = path0
        else
          path = trim(path)//'/'//path0
        end if
      end do

      if (len_trim(path) > 0) then
        path = trim(path)//'/'//filename
      else
        path = filename
      end if
      call open_file(file%units(depth), path, IOFileInput)

    end if

    idx = scan(filename, '/', .true.)
    if (idx == 0) then
      path = '.'
    else
      path = filename(:idx-1)
      if (path == '') &
        path = '/'
    end if

    file%depth = depth
    file%paths(depth) = path

    return

  end subroutine open_include

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check skip
  !> @brief        check skip or not
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_skip(file)

    ! return value
    logical :: check_skip

    ! formal arguments
    type(s_gropp_file),      intent(inout) :: file

    ! local variables
    integer                  :: i, nst_no


    nst_no = file%nstNo(file%depth)

    check_skip = .true.
    if (nst_no <= 0 .or. MaxNstDepth < nst_no) &
         return
    
    do i = 1, nst_no
      if (file%skips(file%depth, i)) return
    end do

    check_skip = .false.
    return

  end function check_skip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_line
  !> @brief        read line
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function read_line(file, str)

    ! return value
    logical :: read_line

    ! formal arguments
    type(s_gropp_file),      intent(inout) :: file
    character(*),            intent(inout) :: str

    ! local variables
    integer                  :: i, idx, unit_no, depth
    character(MaxLine)       :: str1
    logical                  :: bs


    depth = file%depth
    unit_no = file%units(depth)

    str = ''

100 read(unit_no,'(A)',end=900,err=900) str1

    if (file%rest_flags(depth)) &
      file%rest_reads(depth) = file%rest_reads(depth) + 1

    idx = index(str1,';')
    if (idx /= 0) &
      str1(idx:) = ' '

    idx = index(str1,char(9),.true.)  !TAB
    if (idx /= 0) then
      do i = 1, idx
        if (str1(i:i) == char(9)) &
          str1(i:i) = ' '
      end do
    end if

    idx = index(str1,char(92)) !BACKSLASH
    if (idx /= 0) &
      str1(idx:) = ' '
    bs = (idx /= 0)

    if (len_trim(str) == 0) then
      str = str1
    else
      str = trim(str)//' '//str1
    end if

    if (bs) &
      goto 100

    read_line = .true.
    return

900 read_line = .false.
    if (file%rest_flags(depth)) &
      file%rest_reads(depth) = file%rest_reads(depth) + 1
    return

  end function read_line

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_delimit
  !> @brief        check the delimiter
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function is_delimit(c)

    ! return value
    logical :: is_delimit

    ! formal arguments
    character,               intent(in)    :: c


    is_delimit = .not.((c >= 'A' .and. c <= 'Z') .or. &
                       (c >= 'a' .and. c <= 'z') .or. &
                       (c >= '0' .and. c <= '9') .or. &
                       c == '_')

  end function is_delimit

end module fileio_gropp_mod
