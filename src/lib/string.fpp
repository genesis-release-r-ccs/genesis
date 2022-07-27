!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   string_mod
!> @brief   utilities for treating characters
!! @authors Yuji Sugita (YS), Chigusa Kobayashi (CK), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module string_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! variables
  integer,      public, parameter :: MaxLine          = 1500
  integer,      public, parameter :: MaxLineLong      = 50000
  integer,      public, parameter :: MaxLineLong_CV   = 100000
  integer,      public, parameter :: LineBuffer       = 128
  integer,      public, parameter :: MaxFilename      = MaxLine
  integer,      public, parameter :: MaxFilenameLong  = MaxLineLong
  integer,      public, parameter :: MaxMultiFilename = MaxLine

  ! subroutines and functions
  !
  public  :: char_line
  public  :: read_int 
  public  :: read_real
  public  :: read_word
  public  :: read_comment
  public  :: read_ndata  
  public  :: split
  public  :: split_num
  public  :: extract
  public  :: tolower
  public  :: toupper
  public  :: parse_tristate_string

  private :: is_blank
  private :: not_blank
  private :: split_string
  private :: split_int
  private :: split_real

  interface split
    module procedure split_string
    module procedure split_int
    module procedure split_real
  end interface split

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    char_line
  !> @brief        remove comments and blank characters in a line
  !! @authors      YS
  !! @param[in]    max_row : max. row in one line
  !! @param[out]   line    : a line containing data
  !! @param[out]   nchar   : number of characters in the line
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine char_line(max_row, line, nchar)

    ! formal arguments
    integer,                 intent(in)    :: max_row
    character(*),            intent(inout) :: line
    integer,                 intent(inout) :: nchar

    ! local variables
    integer                  :: if_dat, nstrow, nenrow, i


    if_dat = -1
    nstrow = -1
    do i = 1, max_row
      if (not_blank(line(i:i)).and.(if_dat == -1)) then 
        nstrow = i
        exit
      end if
    end do

    if_dat = -1
    do i = max_row, 1, -1
      if (not_blank(line(i:i)).and.(if_dat == -1)) then
        nenrow = i
        exit
      end if
    end do

    if (nstrow == -1) then
      nchar = 0
      return 
    end if

    if (line(nstrow:nstrow) == '!') then
      nenrow = 0
      goto 300
    end if

    do i = nstrow+1, nenrow
      if (line(i:i) == '!') then
        nenrow = i-1
        exit
      end if
    end do

300 continue

    nchar = nenrow - nstrow + 1
    if (nstrow < 0) &
      nchar = 0

    if (nchar > 0) then
      if ( is_blank(line(1:nchar))) then
        line(1:nchar) = line(nstrow:nenrow)
     
      else if ( not_blank(line(1:nchar))) then
        line(1:nchar) = line(nstrow:nenrow)
        if (nchar < max_row) then
          nchar = nchar + 1
          line(nchar:nchar) = ' '
        endif
     
      end if
    end if

    return

  end subroutine char_line

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_int
  !> @brief        read an integer in a line
  !! @authors      YS
  !! @param[in]    line : input data in one line 
  !! @param[out]   nsta : starting row in the line
  !! @param[out]   nend : ending row in the line
  !! @param[out]   iii  : data (integer)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_int(line, nsta, nend, iii)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(inout) :: nsta
    integer,                 intent(inout) :: nend
    integer,                 intent(out)   :: iii

    ! local variables
    integer                  :: nend_bk, i
    logical                  :: exist


    exist = .false.
    nend_bk = nend

    do i = nsta, nend
      if (not_blank(line(i:i))) &
        exist = .true.

      if (is_blank(line(i:i)) .and. (exist)) then
        nend = i
        goto 100
      end if
    end do

    if (.not. exist) &
      call error_msg('  read_int> ERROR: no data exist')

100 continue

    read(line(nsta:nend), *) iii
    nsta = nend
    nend = nend_bk

    return

  end subroutine read_int

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_real
  !> @brief        read a real from the line
  !! @authors      YS
  !! @param[in]    line : input data in one line
  !! @param[out]   nsta : starting row in the line
  !! @param[out]   nend : ending row in the line
  !! @param[out]   rrr  : data (real)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_real(line, nsta, nend, rrr)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(inout) :: nsta
    integer,                 intent(inout) :: nend
    real(wp),                intent(out)   :: rrr

    ! local variables
    integer                  :: nend_bk, i
    logical                  :: exist

    
    exist = .false.
    nend_bk = nend

    do i = nsta, nend
      if (not_blank(line(i:i))) &
        exist = .true.

      if (is_blank(line(i:i)) .and. (exist)) then
        nend = i
        goto 100
      end if
    end do

    if (.not. exist) &
      call error_msg('  read_real> ERROR: no data exist')

100 continue

    read(line(nsta:nend), *) rrr
    nsta = nend
    nend = nend_bk

    return

  end subroutine read_real

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_word
  !> @brief        read a word from the line
  !! @authors      YS
  !! @param[in]    line : input data in one line
  !! @param[out]   nsta : starting row in the line
  !! @param[out]   nend : ending row in the line
  !! @param[in]    nc0  : number of characters
  !! @param[out]   ccc  : data (word)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_word(line, nsta, nend, nc0, ccc)

    ! parameter
    integer,       parameter :: NumAlpha = 26

    character(1),  parameter :: csmall(NumAlpha) = &
                                        (/'a','b','c','d','e','f','g', &
                                          'h','i','j','k','l','m','n', &
                                          'o','p','q','r','s','t','u', &
                                          'v','w','x','y','z'/)
    character(1),  parameter :: clarge(NumAlpha) = &
                                        (/'A','B','C','D','E','F','G', &
                                          'H','I','J','K','L','M','N', &
                                          'O','P','Q','R','S','T','U', &
                                          'V','W','X','Y','Z'/)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(inout) :: nsta
    integer,                 intent(inout) :: nend
    integer,                 intent(in)    :: nc0
    character(*),            intent(out)   :: ccc

    ! local variables
    integer                  :: nend_bk, nc, nspace, i, j, k
    logical                  :: exist


    exist   = .false.
    nend_bk = nend
    nc      = 0

    do i = nsta, nend

      if (is_blank(line(i:i)) .and. (exist)) then
        nend = i
        goto 100

      end if

      if (not_blank(line(i:i)) .and. (nc == nc0)) then
        nend = i - 1
        goto 100

      else if (not_blank(line(i:i)) .and. (nc < nc0)) then
        nc = nc + 1
        ccc(nc:nc) = line(i:i)
        exist = .true.

      end if

    end do

    if (.not. exist) &
      call error_msg('  read_word> ERROR: no data exist')

100 continue
    nsta = nend
    nend = nend_bk

    do j = 1, nc
      do k = 1, NumAlpha
        if (ccc(j:j) .eq. csmall(k)) &
          ccc(j:j) = clarge(k)
      end do
    end do

    nspace = nc0 - nc
    if (nspace > 0) then 
      do j = 1, nspace
        ccc(nc+j:nc+j) = ' '
      end do
    end if

    return

  end subroutine read_word

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_comment
  !> @brief        read a comment after "!"
  !! @authors      YS
  !! @param[in]    line : input data in one line
  !! @param[in]    nsta : starting row in the line
  !! @param[in]    nend : ending row in the line
  !! @param[out]   cccc : comment
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_comment(line, nsta, nend, cccc)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nsta
    integer,                 intent(in)    :: nend
    character(*),            intent(out)   :: cccc

    ! local variables
    logical                  :: exist
    integer                  :: nread, i


    exist = .false.
    
    do i = nsta, nend
      if (line(i:i) .eq. '!') then
        nread = i + 1
        exist = .true.
        exit
      end if
    end do

    if (exist) &
      read(line(nread:nend),'(a60)') cccc

    return

  end subroutine read_comment

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ndata
  !> @brief        read the number of data in a line
  !! @authors      YS
  !! @param[in]    line  : input data in one line
  !! @param[in]    nsta  : starting row in the line
  !! @param[in]    nend  : ending row in the line
  !! @param[out]   ndata : number of data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ndata(line, nsta, nend, ndata)

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(in)    :: nsta
    integer,                 intent(in)    :: nend
    integer,                 intent(out)   :: ndata

    ! local variables
    integer                  :: i
    logical                  :: exist


    ndata = 0
    exist = .false.

    do i = nsta, nend
      if (not_blank(line(i:i)) .and. (.not. exist)) then
        exist = .true.
      end if

      if (is_blank(line(i:i)) .and. (exist)) then
        ndata = ndata + 1
        exist = .false.
      end if

    end do
    if (not_blank(line(nend:nend)) .and. exist) then
      ndata = ndata + 1
    endif

    return

  end subroutine read_ndata

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      split_num
  !> @brief        split clone function   
  !! @authors      CK
  !! @return       number of splitted characters
  !! @param[in]    chr : input character
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function split_num(chr)

    ! return value
    integer                  :: split_num

    ! formal arguments
    character(*),            intent(in)    :: chr

    ! local variables
    integer                  :: i, ifg, len


    split_num = 0

    len = len_trim(chr)
    if (len == 0) return
      
    ifg = 1
    split_num = split_num + 1

    do i = 1,len
      if (chr(i:i) .ne. " " .and. chr(i:i) .ne. char(9)) then
        ifg = 0
      else 
        if (ifg /= 1) then
          split_num = split_num+1
        end if
        ifg = 1
      end if
    end do

    if (chr(len:len) .eq. char(9)) split_num = split_num - 1

    return

  end function split_num

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      extract
  !> @brief        extract string from comma separated string list
  !! @authors      NT
  !! @param[in]    list : string list
  !! @param[inout] idx  : index of list
  !! @param[out]   str  : string
  !! @return       flag for string is extract
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function extract(list, idx, str)

    ! return
    logical                  :: extract

    ! formal arguments
    character(*),            intent(in)    :: list
    integer,                 intent(inout) :: idx
    character(*),            intent(inout) :: str

    ! local variables
    integer                  :: i, ib, ie
    character(MaxLineLong)   :: work


    work = list

    ib = 1
    ie = scan(work,',')

    do i = 1, idx
      if (ie == 0) then
        extract = .false.
        return
      end if
      work(ie:ie) = ' '
      ib = ie
      ie = scan(work,',')
    end do

    if (ie == 0) then
      str = adjustl(work(ib:))
    else
      str = adjustl(work(ib:ie-1))
    end if

    idx = idx + 1
    extract = len_trim(str) > 0
    return

  end function extract

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    tolower
  !> @brief        change capital letters into small letters
  !! @authors      YS
  !! @param[inout] str : string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine tolower(str)

    ! formal arguments
    character(*),            intent(inout) :: str

    ! local variables
    integer                  :: i


    do i = 1, len(str)
      if (str(i:i) .ge. 'A' .and. str(i:i) .le. 'Z') then
        str(i:i) = char(ichar(str(i:i)) + 32)
      end if
    end do
       
    return

  end subroutine tolower

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    toupper
  !> @brief        change small letters into capital letters
  !! @authors      YS
  !! @param[inout] str : string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine toupper(str)

    ! formal arguments
    character(*),            intent(inout) :: str

    ! local variables
    integer                  :: i


    do i = 1, len(str)
      if (str(i:i) .ge. 'a' .and. str(i:i) .le. 'z') then
        str(i:i) = char(ichar(str(i:i)) - 32)
      end if
    end do
       
    return

  end subroutine toupper

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    is_blank
  !> @brief        check for character is blank or not
  !! @authors      YS
  !! @param[in]    c : character
  !! @return       flag for blank or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function is_blank(c)

    ! return value
    logical                  :: is_blank

    ! formal arguments
    character,               intent(in)    :: c


    is_blank = (c .eq. ' ' .or. c .eq. char(9))
    return

  end function is_blank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    not_blank
  !> @brief        check for character is blank or not
  !! @authors      YS
  !! @param[in]    c : character
  !! @return       flag for blank or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function not_blank(c)

    ! return value
    logical        :: not_blank

    ! formal arguments
    character,               intent(in) :: c


    not_blank = (c .ne. ' ' .and. c .ne. char(9))
    return

  end function not_blank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      split_string
  !> @brief        split clone function   
  !! @authors      CK, NT
  !! @param[in]    nval       : number of values
  !! @param[in]    nmax       : length of value dimension
  !! @param[in]    char       : input character
  !! @param[out]   value_char : splited values
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine split_string(nval, nmax, chr, value_char)

    ! formal arguments
    integer,                 intent(in)    :: nval
    integer,                 intent(in)    :: nmax
    character(*),            intent(in)    :: chr
    character(*),            intent(out)   :: value_char(1:nmax)

    ! local variables
    integer                  :: i, ib, ie


    ib = 1
    do i = 1, nval

      do while(chr(ib:ib) .eq. ' ')
        ib = ib + 1
        if (ib >= len(chr)) &
          exit
      end do
      if (ib == len(chr)) &
        exit

      ie = ib + 1
      do while(chr(ie:ie) .ne. ' ')
        ie = ie + 1
        if (ie > len(chr)) &
          exit
      end do

      value_char(i) = chr(ib:ie-1)

      ib = ie
      
    end do

    value_char(i:nmax) = ''

    return

999 call error_msg('Split_char> internal-read failed')

  end subroutine split_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      split_int
  !> @brief        split clone function   
  !! @authors      CK
  !! @param[in]    nval      : number of values
  !! @param[in]    nmax      : length of value dimension
  !! @param[in]    char      : input character
  !! @param[out]   value_int : splited values
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine split_int(nval, nmax, chr, value_int)

    ! formal arguments
    integer,                 intent(in)    :: nval
    integer,                 intent(in)    :: nmax
    character(*),            intent(in)    :: chr
    integer,                 intent(out)   :: value_int(1:nmax)

    ! local variables
    integer                  :: i


    if (nval > 0) &
      read(chr,*,err=999) (value_int(i),i=1,nval) 

    do i = nval+1, nmax
      value_int(i) = 0
    end do

    return

999 call error_msg('Split_int> internal-read failed')

  end subroutine split_int

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      split_real
  !> @brief        split clone function   
  !! @authors      CK
  !! @param[in]    nval       : number of values
  !! @param[in]    nmax       : length of value dimension
  !! @param[in]    char       : input character
  !! @param[out]   value_real : splited values (real)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine split_real(nval, nmax, chr, value_real)

    ! formal arguments
    integer,                 intent(in)    :: nval
    integer,                 intent(in)    :: nmax
    character(*),            intent(in)    :: chr
    real(wp),                intent(out)   :: value_real(1:nmax)

    ! local variables
    integer                  :: i


    if (nval > 0) &
      read(chr,*,err=999) (value_real(i),i=1,nval) 

    do i = nval+1,nmax
      value_real(i) = 0.0_wp
    end do

    return

999 call error_msg('Split_real> internal-read failed')

  end subroutine split_real

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      parse_tristate_string
  !> @brief        convert string to tristate value
  !! @authors      DS
  !! @param[in]    string : character string
  !! @return       tristate(integer) value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function parse_tristate_string(string) result(tristate)

    ! return value
    integer                                :: tristate

    ! formal arguments
    character(*),            intent(in)    :: string

    ! local variables
    character(len(string))                 :: string_lower


    if (string == '') then
      tristate = tristate_NOT_SET
      return
    end if

    string_lower = string
    call tolower(string_lower)

    select case(string_lower)
    case ("true", "yes", "1", "y", "t")
      tristate = tristate_TRUE
    case ("false", "no", "0", "n", "f")
      tristate = tristate_FALSE
    case default
      call error_msg('parse_tristate_string> invalid value')
    end select

    return

  end function parse_tristate_string

end module string_mod
