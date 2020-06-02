!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_lexer_mod
!> @brief   lexical analyzer for atom selection
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_lexer_mod

  use string_mod

  implicit none
  private

  ! token type
  !     SLTTString : string token
  !     SLTTOther  : other token
  integer,         public,  parameter :: SLTTString = 1
  integer,         public,  parameter :: SLTTOther  = 2

  ! string token
  integer,         private, parameter :: NumToken = 7
  character(*),    private, parameter :: Token(NumToken) = &
                                        (/'(',')','&','|','!',':','-'/)

  ! space (delimiter)
  integer,         private, parameter :: NumSpace = 2
  character(1),    private, parameter :: Space(NumSpace) = (/' ',char(9)/)

  ! global variables
  integer,                private            :: g_lex_ptr
  integer,                private            :: g_lex_str_len
  character(MaxLineLong), private            :: g_lex_str

  ! subroutines
  public  :: lexer_set
  public  :: lexer_next_token
  public  :: lexer_back

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    lexer_set
  !> @brief        set string for lexical analysis
  !! @authors      NT
  !! @param[in]    str : string for atoms election
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine lexer_set(str)

    ! formal arguments
    character(*),            intent(in)    :: str


    g_lex_ptr = 1
    g_lex_str_len = len_trim(str)
    g_lex_str = str
    
    return

  end subroutine lexer_set

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      lexer_next_token
  !> @brief        set the next token
  !! @authors      NT
  !! @param[inout] token_str  : next token
  !! @param[out]   token_type : type of token (string/others)
  !! @return       flag for token is set or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function lexer_next_token(token_str, token_type)

    ! return value
    logical        :: lexer_next_token

    ! formal arguments
    character(*),            intent(inout) :: token_str
    integer,                 intent(out)   :: token_type

    ! local variables
    integer                  :: i, token_len


    lexer_next_token = .false.
    token_str = ''

    if (g_lex_ptr > g_lex_str_len) &
      return

    do i = g_lex_ptr, g_lex_str_len

      if (.not. check_space(i)) &
        exit

    end do

    g_lex_ptr = i

    if (check_token(g_lex_ptr, token_len)) then

      token_str = g_lex_str(g_lex_ptr:g_lex_ptr + token_len - 1)
      token_type = SLTTString
      g_lex_ptr = g_lex_ptr + token_len

    else

      do i = g_lex_ptr, g_lex_str_len
        if (check_space(i) .or. check_token(i, token_len)) &
          exit
      end do

      token_str = g_lex_str(g_lex_ptr:i - 1)
      token_type = SLTTOther
      g_lex_ptr = i

    end if

    lexer_next_token = .true.
    return

  end function lexer_next_token

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    lexer_back
  !> @brief        set the lexical pointer back
  !! @authors      NT
  !! @param[in]    token_str  : token
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine lexer_back(token_str)

    ! formal arguments
    character(*),            intent(in)    :: token_str


    g_lex_ptr = g_lex_ptr - len_trim(token_str)

    return

  end subroutine lexer_back

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_space
  !> @brief        check space
  !! @authors      NT
  !! @param[in]    iptr  : pointer
  !! @return       flag for check result
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_space(iptr)

    ! return value
    logical        :: check_space

    ! formal arguments
    integer,                 intent(in)    :: iptr

    ! local variable
    integer                  :: i
    character                :: c


    c = g_lex_str(iptr:iptr)
    do i = 1, NumSpace
      if (c == Space(i)) then
        check_space = .true.
        return 
      end if
    end do

    check_space = .false.
    return

  end function check_space

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_token
  !> @brief        check token
  !! @authors      NT
  !! @param[in]    iptr      : pointer
  !! @param[out]   token_len : length of token
  !! @return       flag for check result
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_token(iptr, token_len)

    ! return value
    logical        :: check_token

    ! formal arguments
    integer,                 intent(in)    :: iptr
    integer,                 intent(out)   :: token_len

    ! local variables
    integer                  :: i
    

    do i = 1, NumToken
      token_len = len_trim(Token(i))
      if (g_lex_str(iptr:iptr + token_len - 1) == Token(i)) then
        check_token = .true.
        return
      end if
    end do

    check_token = .false.
    return

  end function check_token

end module select_lexer_mod
