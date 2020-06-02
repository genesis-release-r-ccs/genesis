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

  ! structures
  type s_gropp_file
     integer                         :: unit_no            = InvalidUnitNo
     character(MaxLine)              :: filename           = ''
     character(MaxLine)              :: paths(MaxIncDepth) = '.'
     integer                         :: depth              = 1
  end type s_gropp_file

  type s_gro_macro
     character(50)                   :: def = ''
     character(MaxLine)              :: val = ''
     integer                         :: def_len = 0
     integer                         :: val_len = 0
     type(s_gro_macro), pointer      :: next => null()
  end type s_gro_macro

  ! variables
  type(s_gropp_file),  private, save :: g_gropp_file(MaxGroPPFile)

  ! subroutines
  public  :: gro_pp_open_file
  public  :: gro_pp_close_file
  private :: pp_body
  private :: pp_ifdef
  private :: dealloc_macro
  private :: print_macro
  private :: predef_macro
  private :: expand_macro
  private :: open_include
  private :: close_include
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

#if defined(INTEL)
    use ifport, only:getpid

#elif defined(KCOMP) || defined(RICC)
    use service_routines, only:getpid

#else
    integer :: getpid

#endif

    ! return value
    integer :: gro_pp_open_file

    ! formal arguments
    character(*),            intent(in)    :: filename
    character(*),            intent(in)    :: predefs(:)
    character(*),            intent(inout) :: error

    ! local variables
    integer                  :: ifile, file_in, file_out, alloc_stat
    integer                  :: depth
    character(MaxLine)       :: pp_filename

    type(s_gro_macro), pointer :: macro_head
    type(s_gro_macro), pointer :: macro_cur

    
    file_in  = InvalidUnitNo
    file_out = InvalidUnitNo


    ! check file handle
    !
    do ifile = 1, MaxGroPPFile
      if (g_gropp_file(ifile)%unit_no == InvalidUnitNo) then
        g_gropp_file(ifile)%paths(1:MaxIncDepth) = '.'
        g_gropp_file(ifile)%depth                = 1
        exit
      end if
    end do
    if (ifile > MaxGroPPFile) then
      error = 'file open error. [open file > MaxFile]'
      return
    end if


    ! allocate head macro
    !
    allocate(macro_head, stat = alloc_stat)
    if (alloc_stat /= 0) then
      error = 'memory allocation error.'
      goto 900
    end if

    macro_cur => macro_head
    macro_cur%def = 'M_A_C_R_O_H_E_A_D'


    ! allocate pre-defined macro
    !
    call predef_macro(predefs, macro_cur)


    ! open include file (read)
    !
    call open_include(file_in, filename, ifile)
    if (file_in == InvalidUnitNo) then
      error = 'file open error. ['//trim(filename)//']'
      goto 900
    end if


    ! open preprocessed file (write)
    !
    write(pp_filename,'(a,a,i0)') &
         trim(filename),      &
         '__gro_pp_file_pid', &
         getpid()

    call open_file(file_out, pp_filename, IOFileOutputReplace)
    if (file_out == InvalidUnitNo) then
      error = 'file open error. [Temporary-file]'
      goto 900
    end if


    ! preprocessing recursively
    !
    depth = 0
    error = ''

    call pp_body(ifile, file_in, file_out, macro_cur, macro_head, &
                 .true., depth, error)
    if (depth /= 0) then
      error = 'Bad pre-processor command state.'
      goto 900
    end if
    if (error /= '') &
      goto 900

    !DEBUG
    !call print_macro(macro_head)


    ! close file
    !
    call close_file(file_out)
    call close_include(file_in, ifile)


    ! deallocate all macro
    !
    call dealloc_macro(macro_head)


    ! open preprocessed file (read)
    !
    call open_file(gro_pp_open_file, pp_filename, IOFileInput)
    if (gro_pp_open_file == InvalidUnitNo) then
      error = 'file open error. [Preprocessed-file]'
      goto 910
    end if

    g_gropp_file(ifile)%unit_no  = gro_pp_open_file
    g_gropp_file(ifile)%filename = pp_filename

    return

900 call close_file(file_out)
    call close_file(file_in)
    call dealloc_macro(macro_head)
910 call unlink(pp_filename)

    gro_pp_open_file = InvalidUnitNo

    return

  end function gro_pp_open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gro_pp_close_file
  !> @brief        close gromacs file
  !! @authors      NT
  !! @param[in]    unit_no : file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gro_pp_close_file(unit_no)

    ! formal arguments
    integer,                 intent(in)    :: unit_no

    ! local variables
    integer                  :: i


    if (unit_no == InvalidUnitNo) &
      return

    call close_file(unit_no)

    do i = 1, MaxGroPPFile
      if (unit_no == g_gropp_file(i)%unit_no) then
        call unlink(g_gropp_file(i)%filename)
        g_gropp_file(i)%unit_no = InvalidUnitNo
        g_gropp_file(i)%filename = ''
        exit
      end if
    end do

    return

  end subroutine gro_pp_close_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pp_body
  !> @brief        preprocess file body
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine pp_body(ifile, file_in, file_out, &
                               macro_cur, macro_head, read_on, depth, error)

    ! formal arguments
    integer,                    intent(in)    :: ifile
    integer,                    intent(in)    :: file_in
    integer,                    intent(in)    :: file_out
    type(s_gro_macro), pointer                :: macro_cur
    type(s_gro_macro), pointer, intent(in)    :: macro_head
    logical,                    intent(in)    :: read_on
    integer,                    intent(inout) :: depth
    character(*),               intent(inout) :: error

    ! local variables
    integer                     :: file_in_inc, alloc_stat
    character(MaxLine)          :: str, str1, str2, str3, str4

    type(s_gro_macro), pointer  :: macro, macro2


    do while(.true.)

      ! read line without comment
      !
100   if (.not.read_line(file_in, str)) &
        exit

      str1 = ''
      read(str,*,err=100,end=100) str1

      ! process #ifdef
      !
      if (str1 == '#ifdef') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call pp_ifdef(ifile, file_in, file_out, macro_cur, macro_head, &
             str3, read_on, .true., depth, error)
        if (error /= '') &
             return

      ! process #ifndef
      !

      else if (str1 == '#ifndef') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call pp_ifdef(ifile, file_in, file_out, macro_cur, macro_head, &
             str3, read_on, .false., depth, error)
        if (error /= '') &
             return

      ! process #else
      !
      else if (str1 == '#else') then

        depth = depth - 1
        return

      ! process #endif
      !
      else if (str1 == '#endif') then

        depth = depth - 1
        return

      ! process #include
      !
      else if (read_on .and. str1 == '#include') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3

        call open_include(file_in_inc, str3, ifile)
        if (file_in_inc == InvalidUnitNo) &
          goto 920

        call pp_body(ifile, file_in_inc, file_out, macro_cur, macro_head, &
                     read_on, depth, error)

        call close_include(file_in_inc, ifile)

        if (error /= '') &
          return

      ! process #define
      !
      else if (read_on .and. str1 == '#define') then

        str2 = ''
        str3 = ''
        read(str,*,err=900,end=900) str2, str3
        str4 = adjustl(str((index(str, trim(str3)) + len_trim(str3)):))

        macro => macro_head
        do while(associated(macro))
          if (str3 == macro%def) then
            exit
          end if
          macro => macro%next
        end do

        ! macro expansion in macro
        macro2 => macro_head
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
          macro_cur%next => macro
          macro_cur => macro

        end if


      ! unknown pre-processor command
      !
      else if (read_on .and. str1(1:1) == '#') then

        if (main_rank) &
          write(MsgOut,*) &
               'Fileio_Gropp> WARNING: unknown pre-process command ['// &
               trim(str1)//']'


      ! process non pre-processor line
      !
      else if (read_on) then

        ! macro expansion
        macro => macro_head
        do while(associated(macro)) 
          call expand_macro(macro, str)
          macro => macro%next
        end do

        write(file_out,'(A)') trim(str)

      end if

    end do

800 continue

    return

900 error = 'format error.'
    return
910 error = 'memory allocation error.'
    return
920 error = 'file open error. include file['//trim(str3)//']'
    return

  end subroutine pp_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pp_ifdef
  !> @brief        preprocess ifdef
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine pp_ifdef(ifile, file_in, file_out, &
                                 macro_cur, macro_head,   &
                                 def_str, read_on, positive, depth, error)

    ! formal arguments
    integer,                    intent(in)    :: ifile
    integer,                    intent(in)    :: file_in
    integer,                    intent(in)    :: file_out
    type(s_gro_macro), pointer                :: macro_cur
    type(s_gro_macro), pointer, intent(in)    :: macro_head
    character(*),               intent(in)    :: def_str
    logical,                    intent(in)    :: read_on
    logical,                    intent(in)    :: positive
    integer,                    intent(inout) :: depth
    character(*),               intent(inout) :: error

    ! local variables
    logical                     :: defined, is_true
    character(MaxLine)          :: str, str1
    type(s_gro_macro), pointer  :: macro


    ! check macro
    !
    defined = .false.

    macro => macro_head
    do while(associated(macro))
      if (def_str == macro%def) then
        defined = .true.
        exit
      end if
      macro => macro%next
    end do

    is_true = (      positive .and.       defined .or. &
               .not. positive .and. .not. defined)


    ! TRUE block
    !
    depth = depth + 1
    call pp_body(ifile, file_in, file_out, macro_cur, macro_head, &
                 is_true .and. read_on, depth, error)
    if (error /= '') &
      return

    backspace(file_in)

    if (.not. read_line(file_in, str)) &
      goto 900

    str1 = ''
    read(str,*,err=900,end=900) str1

    if (str1 == '#endif') &
      return

    if (str1 /= '#else') &
      goto 900


    ! FALSE block
    !
    depth = depth + 1
    call pp_body(ifile, file_in, file_out, macro_cur, macro_head, &
                 (.not. is_true) .and. read_on, depth, error)

    return

900 error = 'format error.'
    return

  end subroutine pp_ifdef

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
  !  Subroutine    print_macro
  !> @brief        print macro list
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_macro(macro_head)

    ! formal arguments
    type(s_gro_macro), pointer, intent(in)    :: macro_head

    ! local variables
    integer                     :: i

    type(s_gro_macro), pointer  :: macro_cur


    i = 1
    write(MsgOut,*) '----------------------------------------------'
    macro_cur => macro_head
    do while(associated(macro_cur))
      write(MsgOut,'(i4,a,a)') &
            i,' ['//trim(macro_cur%def)//']:"', trim(macro_cur%val)//'"'
      macro_cur => macro_cur%next
      i = i + 1
    end do
    write(MsgOut,*) '----------------------------------------------'

    return

  end subroutine print_macro

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
  !  Function      open_include
  !> @brief        open include file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_include(unit_no, filename, ifile)

    ! formal arguments
    integer,                 intent(out)   :: unit_no
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: ifile

    ! local variables
    integer                  :: i, depth, idx
    character(MaxLine)       :: path, path0


    depth = g_gropp_file(ifile)%depth
    if (depth + 1 > MaxIncDepth) then
      unit_no = InvalidUnitNo
      return
    end if

    if (filename(1:1) == '/') then
      call open_file(unit_no, filename, IOFileInput)

    else
      do i = 1, depth
        path0 = g_gropp_file(ifile)%paths(i)
        if (path0(1:1) == '/' .or. i == 1) then
          path = path0
        else
          path = trim(path)//'/'//path0
        end if
      end do

      path = trim(path)//'/'//filename
      call open_file(unit_no, path, IOFileInput)

    end if

    depth = depth + 1
    
    idx = scan(filename, '/', .true.)
    if (idx == 0) then
      path = '.'
    else
      path = filename(:idx-1)
      if (path == '') &
        path = '/'
    end if

    g_gropp_file(ifile)%depth = depth
    g_gropp_file(ifile)%paths(depth) = path

    return

  end subroutine open_include

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      close_include
  !> @brief        close include file
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_include(unit_no, ifile)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    integer,                 intent(in)    :: ifile


    call close_file(unit_no)

    g_gropp_file(ifile)%depth = g_gropp_file(ifile)%depth - 1

    return

  end subroutine close_include

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_line
  !> @brief        read line
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function read_line(unit_no, str)

    ! return value
    logical :: read_line

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    character(*),            intent(inout) :: str

    ! local variables
    integer                  :: i, idx
    character(MaxLine)       :: str1
    logical                  :: bs


    str = ''
    
100 read(unit_no,'(A)',end=900,err=900) str1
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
