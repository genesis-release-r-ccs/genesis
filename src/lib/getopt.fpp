!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   getopt_mod
!> @brief   option parser
!! @authors Motoshi Kamiya (MK)
!
!  (c) Copyright 2015 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#if defined(__GFORTRAN__) && ( __GNUC__ < 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ <= 7 ) )
#define __NO_DEFERRED_LENGTH_CHARACTER__ 1
#endif

module getopt_mod

  use constants_mod
  use string_mod
  use messages_mod

  implicit none
  private

  character, public, parameter :: getopt_empty_option = ' '
  integer,   public, parameter :: getopt_str_length = MaxLine

  type, public :: genesis_opt_item

    character(getopt_str_length) :: name  = ""
    character(getopt_str_length) :: value

  end type genesis_opt_item

  type, public :: genesis_opt_section

    character(getopt_str_length)        :: name  = ""
    integer                             :: count = 0
    type(genesis_opt_item), allocatable :: items(:)

  end type genesis_opt_section

  !! genesis options
  integer,                   public              :: num_genesis_sections
  type(genesis_opt_section), public, allocatable :: genesis_sections(:)

  !! local interger value for current status
  integer                   :: nargs   = 0
  integer                   :: curpos  = 0
#ifdef __NO_DEFERRED_LENGTH_CHARACTER__
  character(MaxLine)        :: arg
#else
  character(:), allocatable :: arg
#endif
  integer                   :: argpos  = 0
  integer                   :: argvpos = 0

  !! corresponding to *optarg
  character(getopt_str_length), public :: getopt_optarg
  integer,                      public :: getopt_optind ! not used now
  integer,                      public :: getopt_opterr ! not used now
  integer,                      public :: getopt_optopt ! not used now

  public :: getopt_initialize
  public :: getopt
  public :: getopt_getarg
  public :: getopt_genesis

  public :: getopt_atof !! convert string to real(wp)
  public :: getopt_atoi !! convert string to integer
  public :: getopt_atol !! convert string to logical
  public :: getopt_type

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    getopt_initialize
  !> @brief        initialize various variables
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine getopt_initialize()

    ! local variables
    integer :: iargc


    nargs = iargc()
    if ( nargs == 0 ) return
    curpos = 0
    argpos = 0
    argvpos = 0
#ifndef __NO_DEFERRED_LENGTH_CHARACTER__
    if (.not. allocated(arg)) then
      allocate(character(getopt_str_length)::arg)
    end if
#endif

    return
  end subroutine getopt_initialize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt_getarg
  !> @brief        get non-optional argument
  !! @authors      MK
  !  @param[in]    options : list of options (getopt style)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  character(getopt_str_length) function getopt_getarg( options ) result(ret)

    ! formal arguments
    character(*), intent(in) :: options

    ! local variables
    integer   :: iargc
    integer   :: char_code
    logical   :: has_arg
    character :: c


    ret = getopt_empty_option
    if (nargs == 0 .and. iargc() /= 0) then
      call getopt_initialize()
    end if

    if (argvpos == nargs) return

    do
      argvpos = argvpos + 1
      if (argvpos > nargs) return
      call getarg(argvpos, arg)
      if (len_trim(arg) == 0) cycle
      if (arg(1:1) /= '-') then
        !! found non-option related argument
        ret = trim(arg)
        return
      end if
      !! check long option
      if (arg(2:2) == '-') then
        if (is_genesis_arg_without_value(arg(3:len_trim(arg)))) then
          argvpos = argvpos + 1
        end if
        cycle
      end if
      !! found options; check last char
      c = arg(len_trim(arg):len_trim(arg))
      char_code = ichar(c)
      if (.not. valid_char(char_code)) call getopt_error_invalid_opt(c)
      !! if it should have argument
      has_arg = check_opt(options, c)
      if (has_arg) then
        argvpos = argvpos + 1
      end if
    end do

    return

  end function getopt_getarg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt
  !> @brief        get an option from command line
  !! @authors      MK
  !  @param[in]    options : list of options (getopt style)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive character function getopt(options) result(ret)

    ! formal arguments
    character(*), intent(in) :: options

    ! local variables
    integer :: iargc


    ret = getopt_empty_option !! default return value
    getopt_optarg = ""
    if (nargs == 0 .and. iargc() /= 0) then
      call getopt_initialize()
    end if

    !! repeat parsing current position
    if (argpos /= 0) then
      argpos = argpos + 1
      ret = parse_next_arg(options)
      return
    end if

    if (curpos == nargs) return !! already finished

    call get_nextarg()
    !! check whether this is an argument
    if (.not. is_arg()) then
      ret = getopt(options)
      return
    end if

    argpos = 2
    ret = parse_next_arg(options)

    return

  end function getopt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    getopt_genesis
  !> @brief        get options for genesis
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine getopt_genesis()

    !! this subroutine is special for genesis

    ! local variables
    integer :: i
    integer :: n_genesis_args, arglen
    integer :: fid, n
    integer,                      allocatable :: count_params(:)
    character(getopt_str_length), allocatable :: section_names(:)
    character(getopt_str_length)              :: fn
    integer :: iargc


    if ( nargs == 0 .and. iargc() /= 0 ) then
      call getopt_initialize()
    end if

    i = 1
    n_genesis_args = 0
    allocate(section_names(nargs), count_params(nargs))
    section_names(1:nargs) = ""

    !! count number of genesis args
    curpos = 0
    do
      if (curpos >= nargs) exit
      call get_nextarg()
      if (arg(1:2) == '--') then
        if (is_genesis_arg(arg(3:len_trim(arg)))) then
          n_genesis_args = n_genesis_args + 1
          section_names(n_genesis_args) = &
              get_genesis_section_name(arg(3:len_trim(arg)))
        end if
      end if
    end do

    curpos = 0
    if (n_genesis_args == 0) return

    !! it might be longest argument size
    arglen = len(arg)
    if (arglen > getopt_str_length) then
      write(MsgOut,'("Warning(getopt)> there may be very long command line argument.")')
      write(MsgOut,'("                 You should check carefully the output or/and result.")')
    end if

    !! get uniq section types, and also count parameter
    call uniq_section_types(n_genesis_args, section_names, &
                            num_genesis_sections, count_params)
    if (num_genesis_sections == 0) return

    !! allocate memory for arguments
    allocate(genesis_sections(num_genesis_sections))
    do i = 1, num_genesis_sections
      genesis_sections(i)%name = section_names(i)
      allocate(genesis_sections(i)%items(count_params(i)))
    end do

    !! reset position and read once more
    curpos = 0
    n_genesis_args = 0
    do
      if (curpos >= nargs) exit
      call get_nextarg()
      !! please ignore many many stupid wastes <(_ _)>
      if (arg(1:2) == '--') then
        if (is_genesis_arg(arg(3:len_trim(arg)))) then
          !! section name
          fn = get_genesis_section_name(arg(3:len_trim(arg)))
          fid = get_section_id(fn)
          genesis_sections(fid)%count = genesis_sections(fid)%count + 1
          n = genesis_sections(fid)%count
          if (is_genesis_arg_without_value(arg(3:len_trim(arg)))) then
            !! parameter name
            genesis_sections(fid)%items(n)%name = &
                get_genesis_param_name( arg(3:len_trim(arg)))
            !! if this is the last arg, put error
            if (curpos >= nargs) then
              call error_msg('ERROR(getopt)> value for the last parameter is not given')
            end if
            !! value
            call get_nextarg()
            genesis_sections(fid)%items(n)%value = arg
          else
            !! parameter name
            genesis_sections(fid)%items(n)%name = &
                get_genesis_param_name( arg(3:len_trim(arg)) )
            !! value
            genesis_sections(fid)%items(n)%value = &
                get_genesis_param_value( arg(3:len_trim(arg)) )
          end if
          if ( len_trim(genesis_sections(fid)%items(n)%value) == &
               len(genesis_sections(fid)%items(n)%value) ) then
            write(MsgOut,'("Error(getopt)> value for ", a, " might be too long")') &
                genesis_sections(fid)%items(n)%name
            call error_msg('Error(getopt)> failed to parse argument')
          end if
        end if
      end if
    end do

    !! clear when exits
    curpos = 0
    deallocate(section_names, count_params)

    return
  end subroutine getopt_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      parse_next_arg
  !> @brief        get next argument
  !! @authors      MK
  !  @param[in]    options : list of options (getopt style)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  character function parse_next_arg( options ) result(ret)

    ! formal arguments
    character(*), intent(in) :: options

    ! local variables
    integer   :: char_code
    logical   :: has_arg
    character :: c


    ret = getopt_empty_option

    !! normal option found; is it a valid character
    c = arg(argpos:argpos)
    char_code = ichar(c)
    if (.not. valid_char(char_code)) call getopt_error_invalid_opt(c)

    !! check the option is in the requested option
    has_arg = check_opt(options, c)

    ret = c
    if (.not. has_arg) then
      if (len_trim(arg) <= argpos) argpos = 0
      return
    end if

    !! error; require option, but not provided
    if (has_arg .and. len_trim(arg) > argpos) then
      call getopt_error_parse_arg(c)
    end if

    !! get arg and return
    call get_nextarg()
    getopt_optarg = trim(arg)

    return

  end function parse_next_arg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_opt
  !> @brief        check inputted option
  !! @authors      MK
  !  @param[in]    options : list of options (getopt style)
  !  @param[in]    o       : current option name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function check_opt( options, o ) result(ret)

    ! formal arguments
    character(*), intent(in) :: options
    character,    intent(in) :: o

    ! local variables
    integer :: i
    logical :: found


    found = .false.
    ret = .false.
    do i = 1, len_trim(options)
      if (options(i:i) == o) then
        found = .true.
        if (i < len_trim(options)) then
          if (options(i+1:i+1) == ':') then
            ret = .true.
          end if
        end if
        exit !! option found
      end if
    end do

    !! error; not registered option
    if (.not. found) call getopt_error_unknown_opt(o)

    return

  end function check_opt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_nextarg
  !> @brief        get next command line argument
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_nextarg()

    ! local variables
    integer :: s


    curpos = curpos + 1
    argpos = 0
    do
      call getarg(curpos, arg)
      if (len(arg) == len_trim(arg)) then
#ifdef __NO_DEFERRED_LENGTH_CHARACTER__
        call error_msg('ERROR(getopt)> argument is too long')
#else
        !! I'm afraid argument is longer than expected
        s = len(arg)
        deallocate(arg)
        allocate(character(2*s)::arg)
#endif
      else
        exit
      end if
    end do

    return

  end subroutine get_nextarg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_arg
  !> @brief        check whether current argument is an option
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function is_arg() result(ret)

    ret = .false.
    if (arg(1:1) == '-') then
      if (arg(2:2) == '-') then
        !! long option found, ignore it
        if (is_genesis_arg_without_value(arg(2:len_trim(arg)))) then
          curpos = curpos + 1
        end if
        return
      end if
      ret = .true.
    end if

    return

  end function is_arg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_genesis_arg
  !> @brief        check whether current argument is a genesis option
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function is_genesis_arg(myarg) result(ret)

    ! formal arguments
    character(*), intent(in) :: myarg

    ! local variables
    integer :: i


    ret = .false.

    !! genesis arg should have : in the string, which separate
    !! section name and parameter name.
    do i = 1, len_trim(myarg)
      if (myarg(i:i) == ':') then
        ret = .true.
        return
      end if
    end do

    return

  end function is_genesis_arg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_genesis_arg_without_value
  !> @brief        check whether current argument is a genesis option
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function is_genesis_arg_without_value( myarg ) result(ret)

    ! formal arguments
    character(*), intent(in) :: myarg

    ! local variables
    integer :: i, j, colonpos


    ret = .false.

    !! genesis arg should have : in the string, which separate
    !! section name and parameter name.
    colonpos = 0
    do i = 1, len_trim(myarg)
      if (myarg(i:i) == ':') then
        colonpos = i
        !! check repeated colons
        j = colonpos + 1
        do
          if (j >= len_trim(myarg)) exit
          if (myarg(j:j) == ':') then
            colonpos = j
          else
            exit
          end if
          j = j + 1
        end do
        exit
      end if
    end do

    !! what this is? normal long option
    if (colonpos == 0) return

    !! search =
    do j = colonpos + 1, len_trim(myarg)
      !! if there is colon before =, it must not be genesis arg
      if (myarg(j:j) == ':') return
      if (myarg(j:j) == '=') return
    end do

    ret = .TRUE.

    return

  end function is_genesis_arg_without_value

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_genesis_section_name
  !> @brief        get section name
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  character(getopt_str_length) function get_genesis_section_name(a) result(ret)

    ! formal arguments
    character(*), intent(in) :: a

    ! local variables
    integer :: i, ic


    ret = ""
    do i = 1, len_trim(a)
      if (a(i:i) == ':') then
        exit
      end if
      !! tolower
      ic = ichar(a(i:i))
      if (65 <= ic .and. ic <= 90) ic = ic + 32
      ret(i:i) = char(ic)
    end do

    return

  end function get_genesis_section_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_section_id
  !> @brief        get section id
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  integer function get_section_id( fn ) result(ret)

    ! formal arguments
    character(getopt_str_length) :: fn

    ! local variables
    integer :: i


    ret = 0
    do i = 1, num_genesis_sections
      if (fn == genesis_sections(i)%name) then
        ret = i
        return
      end if
    end do

    return
  end function get_section_id

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_genesis_param_name
  !> @brief        get param name
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  character(getopt_str_length) function get_genesis_param_name(a) result(ret)

    ! formal arguments
    character(*), intent(in) :: a

    ! local variables
    integer :: i, ic, n
    logical :: found_colon


    ret = ""
    i = 0
    n = 0
    found_colon = .false.
    do
      i = i + 1
      if (i > len_trim(a)) exit

      if (a(i:i) == ':' .and. .not. found_colon) then
        do
          if (i >= len_trim(a)) exit
          if (a(i+1:i+1) == ':') then
            i = i + 1
            cycle
          end if
          exit
        end do
        found_colon = .true.
        cycle
      end if
      if (.not. found_colon) cycle
      if (a(i:i) == '=') exit

      ic = ichar(a(i:i))
      if (65 <= ic .and. ic <= 90) ic = ic + 32
      n = n + 1
      ret(n:n) = char(ic)
    end do

    return

  end function get_genesis_param_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_genesis_param_value
  !> @brief        get param value
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  character(getopt_str_length) function get_genesis_param_value(a) result(ret)

    ! formal arguments
    character(*), intent(in) :: a

    ! local variables
    integer :: i


    ret = ""
    i = 0
    do
      i = i + 1
      if (i > len_trim(a)) exit

      if (a(i:i) == '=' .and. len_trim(a) > i) then
        ret = a(i+1:len_trim(a))
      end if
    end do

    return

  end function get_genesis_param_value

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    uniq_section_types
  !> @brief        process sections
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine uniq_section_types(n, section_names, count_sections, count_params)

    ! formal arguments
    integer,                      intent(in)    :: n
    character(getopt_str_length), intent(inout) :: section_names(*)
    integer,                      intent(out)   :: count_sections
    integer,                      intent(out)   :: count_params(*)

    ! local variables
    integer :: i, j, filled


    filled = 0
    count_params(1:n) = 0

    do i = 1, n

      if (len_trim(section_names(i)) == 0) cycle
      filled = filled + 1
      if (i /= filled) section_names(filled) = section_names(i)
      count_params(filled) = 1

      !! hello i am stupid!
      do j = i + 1, nargs
        if (section_names(j) == section_names(i)) then
          count_params(filled) = count_params(filled) + 1
          section_names(j) = ""
        end if
      end do
    end do

    count_sections = filled

    return

  end subroutine uniq_section_types

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      valid_char
  !> @brief        check whether current argname is valid
  !! @authors      MK
  !! @params[in]   char_code : character code (integer)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function valid_char( char_code ) result(ret)

    ! formal arguments
    integer, intent(in) :: char_code

    ret = .false.

    if ( (48 <= char_code .and. char_code <= 57) .or. &  !! 0-9
         (65 <= char_code .and. char_code <= 90) .or. &  !! A-Z
         (97 <= char_code .and. char_code <= 122)) then !! a-z
      ret = .true.
      return
    end if

    return

  end function valid_char

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    getopt_error_invalid_opt
  !> @brief        put error and quit
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine getopt_error_invalid_opt(c)

    ! formal arguments
    character, intent(in) :: c


    write(MsgOut,'("Error(getopt)> invalid character for an option found.")')
    write(MsgOut,'("               invalid char is: ", a1 )') c
    call error_msg('Error(getopt)> quit program')

    return

  end subroutine getopt_error_invalid_opt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    getopt_error_unknown_opt
  !> @brief        put error and quit
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine getopt_error_unknown_opt(c)

    ! formal arguments
    character, intent(in) :: c


    write(MsgOut,'("Error(getopt)> unknown option type ", a1, " is found")') c
    call error_msg('Error(getopt)> quit program.')

    return

  end subroutine getopt_error_unknown_opt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    getopt_error_parse_arg
  !> @brief        put error and quit
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine getopt_error_parse_arg(c)

    ! formal arguments
    character, intent(in) :: c


    write(MsgOut,'("Error(getopt)> option ", a1, &
            & " requires argument but not found.")') c
    call error_msg('Error(getopt)> quit program.')

    return

  end subroutine getopt_error_parse_arg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt_atof
  !> @brief        convert given string to real(wp)
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  real(wp) function getopt_atof(c) result(ret)

    ! formal arguments
    character(*), intent(in) :: c


    read(c,*) ret

    return

  end function getopt_atof

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt_atoi
  !> @brief        convert given string to integer
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  integer function getopt_atoi(c) result(ret)

    ! formal arguments
    character(*), intent(in) :: c


    read(c,*) ret

    return

  end function getopt_atoi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt_atol
  !> @brief        convert given string to logical
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  logical function getopt_atol(c) result(ret)

    ! formal arguments
    character(*), intent(in) :: c

    ! local variables
    character(getopt_str_length) :: cvar


    cvar = c;
    call tolower(cvar)
    ret = ( cvar == 'yes' .or. cvar == 'true' )

    return

  end function getopt_atol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      getopt_type
  !> @brief        convert given string to corresponding integer
  !! @authors      MK
  !! @params[in]   c : arg
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  integer function getopt_type(c, typeStrs) result(ret)

    ! formal arguments
    character(*), intent(in) :: c
    character(*), intent(in) :: typeStrs(:)

    ! local variables
    integer :: i
    character(getopt_str_length) :: cstr1, cstr2


    cstr1 = c
    call tolower(cstr1)

    ret = 0
    do i = 1, size(typeStrs)

      cstr2 = typeStrs(i)
      call tolower(cstr2)

      if ( cstr1 == cstr2 ) then
        ret = i
        return
      end if

    end do

    call error_msg('Error> cannot find type')

    return

  end function getopt_type

end module getopt_mod

#undef __WITHOUT_DEFERRED_CHARACTER__
