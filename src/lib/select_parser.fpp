!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_parser_mod
!> @brief   expression parser for atom selection
!! @authors Norio Takase (NT), Wataru Nishima (WN)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_parser_mod

  use select_lexer_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structure for parse tree node
  type, public :: s_sp_node

    integer                        :: type     = 0
    integer                        :: priority = 0
    integer                        :: capacity = 0

    type(s_sp_node),       pointer :: parent_node => null()
    type(s_sp_node),       pointer :: child_node1 => null()
    type(s_sp_node),       pointer :: child_node2 => null()

    character(MaxLineLong)         :: param_str = ''
    type(s_selatoms)               :: selatoms

  end type s_sp_node

  ! parameters for node types
  integer,      public,  parameter :: NumNodes      = 19
  integer,      public,  parameter :: NT_Top        = 1
  integer,      public,  parameter :: NT_Bracket    = 2
  integer,      public,  parameter :: NT_AND        = 3
  integer,      public,  parameter :: NT_OR         = 4
  integer,      public,  parameter :: NT_NOT        = 5
  integer,      public,  parameter :: NT_AroundAtom = 6
  integer,      public,  parameter :: NT_AroundRes  = 7
  integer,      public,  parameter :: NT_AroundMol  = 8
  integer,      public,  parameter :: NT_AtomName   = 9
  integer,      public,  parameter :: NT_AtomIdx    = 10
  integer,      public,  parameter :: NT_AtomNo     = 11
  integer,      public,  parameter :: NT_ResName    = 12
  integer,      public,  parameter :: NT_ResNo      = 13
  integer,      public,  parameter :: NT_MolName    = 14
  integer,      public,  parameter :: NT_SegId      = 15
  integer,      public,  parameter :: NT_Hydrogen   = 16
  integer,      public,  parameter :: NT_Backbone   = 17
  integer,      public,  parameter :: NT_Heavy      = 18
  integer,      public,  parameter :: NT_All        = 19

  integer,      public,  parameter :: MaxToken      = 20

  ! subroutines
  public  :: parser_build_tree
  public  :: parser_dealloc_tree
  public  :: parser_print_tree
  public  :: parser_get_tree
  public  :: parser_get_error

  ! node name
  character(*),    private, parameter :: Names(NumNodes) = (/&
                                          'Top       ', &
                                          '()        ', &
                                          '&         ', &
                                          '|         ', &
                                          '!         ', &
                                          'AroundAtom', &
                                          'AroundRes ', &
                                          'AroundMol ', &
                                          'AtomName  ', &
                                          'AtomIdx   ', &
                                          'AtomNo    ', &
                                          'ResName   ', &
                                          'ResNo     ', &
                                          'MolName   ', &
                                          'SegId     ', &
                                          'Hydrogen  ', &
                                          'Backbone  ', &
                                          'Heavy     ', &
                                          'All       '/)

  ! node priority
  integer,         private, parameter :: Priority(NumNodes) = (/&
                                          0, &     ! NT_Top
                                          0, &     ! NT_Bracket
                                          2, &     ! NT_AND
                                          1, &     ! NT_OR
                                          4, &     ! NT_NOT
                                          3, &     ! NT_AroundAtom
                                          3, &     ! NT_AroundRes
                                          3, &     ! NT_AroundMol
                                          0, &     ! NT_AtomName
                                          0, &     ! NT_AtomIdx
                                          0, &     ! NT_AtomNo
                                          0, &     ! NT_ResName
                                          0, &     ! NT_ResNo
                                          0, &     ! NT_MolName
                                          0, &     ! NT_SegId
                                          0, &     ! NT_Hydrogen
                                          0, &     ! NT_Backbone
                                          0, &     ! NT_Heavy
                                          0/)      ! NT_All

  ! node capacity
  integer,         private, parameter :: Capacity(NumNodes) = (/&
                                          1, &     ! NT_Top
                                          1, &     ! NT_Bracket
                                          2, &     ! NT_AND
                                          2, &     ! NT_OR
                                          1, &     ! NT_NOT
                                          1, &     ! NT_AroundAtom
                                          1, &     ! NT_AroundRes
                                          1, &     ! NT_AroundMol
                                          0, &     ! NT_AtomName
                                          0, &     ! NT_AtomIdx
                                          0, &     ! NT_AtomNo
                                          0, &     ! NT_ResName
                                          0, &     ! NT_ResNo
                                          0, &     ! NT_MolName
                                          0, &     ! NT_SegId
                                          0, &     ! NT_Hydrogen
                                          0, &     ! NT_Backbone
                                          0, &     ! NT_Heavy
                                          0/)      ! NT_All

  ! token
  character(*),    private, parameter :: T_BracketL    = '('
  character(*),    private, parameter :: T_BracketR    = ')'
  character(*),    private, parameter :: T_AND         = '&'
  character(*),    private, parameter :: T_AND2        = 'and'
  character(*),    private, parameter :: T_OR          = '|'
  character(*),    private, parameter :: T_OR2         = 'or'
  character(*),    private, parameter :: T_NOT         = '!'
  character(*),    private, parameter :: T_NOT2        = 'not'
  character(*),    private, parameter :: T_AroundAtom  = 'around_atoms'
  character(*),    private, parameter :: T_AroundAtom2 = 'around'
  character(*),    private, parameter :: T_AroundRes   = 'around_residues'
  character(*),    private, parameter :: T_AroundRes2  = 'around_res'
  character(*),    private, parameter :: T_AroundMol   = 'around_molecules'
  character(*),    private, parameter :: T_AroundMol2  = 'around_mol'
  character(*),    private, parameter :: T_AtomName    = 'atomname'
  character(*),    private, parameter :: T_AtomName2   = 'atom_name'
  character(*),    private, parameter :: T_AtomName3   = 'an'
  character(*),    private, parameter :: T_AtomIdx     = 'atomidx'
  character(*),    private, parameter :: T_AtomIdx2    = 'atomindex'
  character(*),    private, parameter :: T_AtomIdx3    = 'ai'
  character(*),    private, parameter :: T_AtomNo      = 'atomno'
  character(*),    private, parameter :: T_AtomNo2     = 'ano'
  character(*),    private, parameter :: T_ResName     = 'resname'
  character(*),    private, parameter :: T_ResName2    = 'residuename'
  character(*),    private, parameter :: T_ResName3    = 'rnam'
  character(*),    private, parameter :: T_ResNo       = 'resno'
  character(*),    private, parameter :: T_ResNo2      = 'residueno'
  character(*),    private, parameter :: T_ResNo3      = 'rno'
  character(*),    private, parameter :: T_MolName     = 'molname'
  character(*),    private, parameter :: T_MolName2    = 'moleculename'
  character(*),    private, parameter :: T_MolName3    = 'mname'
  character(*),    private, parameter :: T_SegId       = 'segid'
  character(*),    private, parameter :: T_SegId2      = 'segmentid'
  character(*),    private, parameter :: T_SegId3      = 'sid'
  character(*),    private, parameter :: T_Hydrogen    = 'hydrogen'
  character(*),    private, parameter :: T_Hydrogen2   = 'hydrogenatom'
  character(*),    private, parameter :: T_Backbone    = 'backbone'
  character(*),    private, parameter :: T_Backbone2   = 'backboneatom'
  character(*),    private, parameter :: T_Heavy       = 'heavy'
  character(*),    private, parameter :: T_Heavy2      = 'heavyatom'
  character(*),    private, parameter :: T_All         = 'all'
  character(*),    private, parameter :: T_All2        = '*'

  ! global variables
  type(s_sp_node), private, target, save :: g_top_node
  type(s_sp_node), private, pointer      :: g_cur_node => null()

  character(100),  private               :: g_error

  ! subroutines
  private :: parse_bracket_l
  private :: parse_bracket_r
  private :: parse_and
  private :: parse_or
  private :: parse_not
  private :: parse_around_atom
  private :: parse_around_res
  private :: parse_around_mol
  private :: parse_atom_name
  private :: parse_atom_idx
  private :: parse_atom_no
  private :: parse_res_name
  private :: parse_res_no
  private :: parse_mol_name
  private :: parse_seg_id
  private :: parse_simple_token
  private :: create_node
  private :: append_node
  private :: dealloc_node
  private :: check_node
  private :: print_node
  private :: is_int
  private :: is_real

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parser_build_tree
  !> @brief        build parse tree
  !! @authors      NT, WN
  !! @param[in]    str : string for atom selection
  !! @param[out]   ret : error flag
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parser_build_tree(str, ret)

    ! formal arguments
    character(*),            intent(in)    :: str
    logical,                 intent(out)   :: ret

    ! local variables
    character(MaxToken)       :: token
    integer                   :: token_type


    g_error = ''

    call lexer_set(str)
    call parser_dealloc_tree

    do while(lexer_next_token(token, token_type))

      call tolower(token)

      if (token == T_BracketL) then
        call parse_bracket_l

      else if (token == T_BracketR) then
        call parse_bracket_r

      else if (token == T_AND .or. token == T_AND2) then
        call parse_and

      else if (token == T_OR .or. token == T_OR2) then
        call parse_or

      else if (token == T_NOT .or. token == T_NOT2) then
        call parse_not

      else if (token == T_AroundAtom .or. token == T_AroundAtom2) then
        call parse_around_atom

      else if (token == T_AroundRes .or. token == T_AroundRes2) then
        call parse_around_res

      else if (token == T_AroundMol .or. token == T_AroundMol2) then
        call parse_around_mol

      else if (token == T_AtomName .or. token == T_AtomName2 .or. &
           token == T_AtomName3) then
        call parse_atom_name

      else if (token == T_AtomIdx .or. token == T_AtomIdx2 .or. &
           token == T_AtomIdx3) then
        call parse_atom_idx

      else if (token == T_AtomNo .or. token == T_AtomNo2) then
        call parse_atom_no

      else if (token == T_ResName .or. token == T_ResName2 .or. &
           token == T_ResName3) then
        call parse_res_name

      else if (token == T_ResNo .or. token == T_ResNo2 .or. &
           token == T_ResNo3) then
        call parse_res_no

      else if (token == T_MolName .or. token == T_MolName2 .or. &
           token == T_MolName3) then
        call parse_mol_name

      else if (token == T_SegId .or. token == T_SegId2 .or. &
           token == T_SegId3) then
        call parse_seg_id

      else if (token == T_Hydrogen .or. token == T_Hydrogen2) then
        call parse_simple_token(NT_Hydrogen)

      else if (token == T_Backbone .or. token == T_backbone2) then
        call parse_simple_token(NT_Backbone)

      else if (token == T_Heavy .or. token == T_heavy2) then
        call parse_simple_token(NT_Heavy)

      else if (token == T_All .or. token == T_All2) then
        call parse_simple_token(NT_All)

      else
         g_error = 'Unknown token: "'// trim(token)//'"'

      end if

      if (g_error /= '') &
        goto 100

    end do

    call check_node(g_top_node%child_node1)


    if (g_error /= '') &
      goto 100

    ret = .true.
    return

100 ret = .false.
    call parser_dealloc_tree
    return

  end subroutine parser_build_tree
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parser_dealloc_tree
  !> @brief        deallocate parse tree
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parser_dealloc_tree

    call dealloc_node(g_top_node%child_node1)
    call dealloc_node(g_top_node%child_node2)

    g_top_node%child_node1 => null()
    g_top_node%child_node2 => null()

    g_top_node%type     = NT_Top
    g_top_node%priority = Priority(NT_Top)
    g_top_node%capacity = Capacity(NT_Top)

    g_cur_node => g_top_node

    return

  end subroutine parser_dealloc_tree

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parser_print_tree
  !> @brief        print parse tree
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parser_print_tree

    write(MsgOut, *) ''

    call print_node(0, g_top_node%child_node1)
    call print_node(0, g_top_node%child_node2)

    write(MsgOut, *) ''

    return

  end subroutine parser_print_tree


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parser_get_tree
  !> @brief        return top node of parse tree
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function parser_get_tree()

    ! return value
    type(s_sp_node), pointer :: parser_get_tree


    parser_get_tree => g_top_node%child_node1
    return

  end function parser_get_tree


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parser_get_error
  !> @brief        get error messages
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function parser_get_error()

    ! return value
    character(100) :: parser_get_error


    parser_get_error = g_error
    return

  end function parser_get_error

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_bracket_l
  !> @brief        parse the character '('
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_bracket_l

    ! local variable
    type(s_sp_node), pointer :: node


    call create_node(NT_Bracket, node)

    ! check preliminary
    if (g_cur_node%capacity == 0) then
      g_error = 'bracket is not a preliminary'
      return
    end if

    ! setup mark
    node%param_str = 'B'

    call append_node(node)

    return

  end subroutine parse_bracket_l

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_bracket_l
  !> @brief        parse the character ')'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_bracket_r

    do while(.true.)

      if (g_cur_node%type == NT_Bracket .and. &
          g_cur_node%param_str == 'B') then

        g_cur_node%param_str = ''
        exit
      end if
      g_cur_node => g_cur_node%parent_node

      ! check unnecessary bracket
      if (.not. associated(g_cur_node)) then
        g_error = 'unnecessary bracket ")"'
        return
      end if

    end do

    ! check bracket inside
    if (g_cur_node%capacity /= 0) then
      g_error = 'bracket is empty'
      return
    end if

    return

  end subroutine parse_bracket_r

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_and
  !> @brief        parse the character '&'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_and

    ! local variable
    type(s_sp_node), pointer :: node


    call create_node(NT_AND, node)
    call append_node(node)

    return

  end subroutine parse_and

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_or
  !> @brief        parse the character '|'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_or

    ! local variable
    type(s_sp_node), pointer :: node


    call create_node(NT_OR, node)
    call append_node(node)

    return

  end subroutine parse_or

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_not
  !> @brief        parse the character '!'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_not

    ! local variable
    type(s_sp_node), pointer :: node


    call create_node(NT_NOT, node)

    ! check preliminary
    if (g_cur_node%capacity == 0) then
      g_error = '['//trim(Names(node%type))//'] not a preliminary'
      return
    end if

    call append_node(node)

    return

  end subroutine parse_not

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_around_atom
  !> @brief        parse the expression 'around_atom'
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_around_atom

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AroundAtom, node)

    ! check postposition
    if (g_cur_node%capacity /= 0) then
      g_error = '['//trim(Names(node%type))//'] not a postposition'
      return
    end if

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_real(token)) then
      g_error = '"'//trim(token)//'": not a real number'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_around_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_around_res
  !> @brief        parse the expression 'around_res'
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_around_res

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AroundRes, node)

    ! check postposition
    if (g_cur_node%capacity /= 0) then
      g_error = '['//trim(Names(node%type))//'] not a postposition'
      return
    end if

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_real(token)) then
      g_error = '"'//trim(token)//'": not a real number'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_around_res

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_around_mol
  !> @brief        parse the expression 'around_mol'
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_around_mol

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AroundMol, node)

    ! check postposition
    if (g_cur_node%capacity /= 0) then
      g_error = '['//trim(Names(node%type))//'] not a postposition'
      return
    end if

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_real(token)) then
      g_error = '"'//trim(token)//'": not a real number'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_around_mol

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_atom_name
  !> @brief        parse the expression 'atom name'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_atom_name

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AtomName, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read name
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no name'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_atom_idx
  !> @brief        parse the expression 'atom_idx'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_atom_idx

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AtomIdx, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read (begin) number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_int(token)) then
      g_error = '"'//trim(token)//'": not a number'
      return
    end if

    node%param_str = token

    ! read "-" (option)
    if (lexer_next_token(token, token_type)) then

      if (token == '-') then

        node%param_str = trim(node%param_str) // token

        ! read end number (option)
        if (lexer_next_token(token, token_type)) then

          if (token_type == SLTTOther) then

            if (is_int(token)) then
              node%param_str = trim(node%param_str) // token
            else
              g_error = '"'//trim(token)//'": not a number'
              return
            end if

          else

            call lexer_back(token)

          end if

        end if

      else

        call lexer_back(token)

      end if

    end if

    call append_node(node)

    return

  end subroutine parse_atom_idx

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_atom_no
  !> @brief        parse the expression 'atom no'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_atom_no

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_AtomNo, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read (begin) number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_int(token)) then
      g_error = '"'//trim(token)//'": not a number'
      return
    end if

    node%param_str = token

    ! read "-" (option)
    if (lexer_next_token(token, token_type)) then

      if (token == '-') then

        node%param_str = trim(node%param_str) // token

        ! read end number (option)
        if (lexer_next_token(token, token_type)) then
             
          if (token_type == SLTTOther) then

            if (is_int(token)) then
              node%param_str = trim(node%param_str) // token
            else
              g_error = '"'//trim(token)//'": not a number'
              return
            end if

          else

            call lexer_back(token)

          end if

        end if

      else

        call lexer_back(token)

      end if

    end if

    call append_node(node)

    return

  end subroutine parse_atom_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_res_name
  !> @brief        parse the expression 'res_name'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_res_name

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_ResName, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read name
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no name'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_res_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_res_no
  !> @brief        parse the expression 'res_no'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_res_no

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_ResNo, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read (begin) number
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no number'
      return
    end if

    if (.not. is_int(token)) then
      g_error = '"'//trim(token)//'": not a number'
      return
    end if

    node%param_str = token

    ! read "-" (option)
    if (lexer_next_token(token, token_type)) then

      if (token == '-') then

        node%param_str = trim(node%param_str) // token

        ! read end number (option)
        if (lexer_next_token(token, token_type)) then
             
          if (token_type == SLTTOther) then

            if (is_int(token)) then
              node%param_str = trim(node%param_str) // token
            else
              g_error = '"'//trim(token)//'": not a number'
              return
            end if

          else

            call lexer_back(token)

          end if

        end if

      else

        call lexer_back(token)

      end if

    end if

    call append_node(node)

    return

  end subroutine parse_res_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_mol_name
  !> @brief        parse the expression 'mol_name'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_mol_name

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_MolName, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read name
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no name'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_mol_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_seg_id
  !> @brief        parse the expression 'seg_id'
  !! @authors      NT, WN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_seg_id

    ! local variable
    integer                  :: token_type
    character(MaxToken)      :: token
    logical                  :: has_token

    type(s_sp_node), pointer :: node


    call create_node(NT_SegId, node)

    ! read ":"
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token /= ':') then
      g_error = '['//trim(Names(node%type))//'] there is no ":"'
      return
    end if

    ! read name
    has_token = lexer_next_token(token, token_type)

    if (.not. has_token .or. token_type == SLTTString) then
      g_error = '['//trim(Names(node%type))//'] there is no id'
      return
    end if

    node%param_str = token

    call append_node(node)

    return

  end subroutine parse_seg_id

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_simple_token
  !> @brief        parse the simple token
  !! @authors      NT, WN
  !! @param[in]    type : node type
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_simple_token(type)

    ! formal arguments
    integer,                 intent(in)    :: type

    ! local variables
    type(s_sp_node), pointer :: node


    call create_node(type, node)
    call append_node(node)

    return

  end subroutine parse_simple_token

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    create_node
  !> @brief        allocate new node
  !! @authors      NT, WN
  !! @param[in]    type : node type
  !! @param[inout] node : pointer to node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine create_node(type, node)

    ! formal arguments
    integer,                  intent(in)    :: type
    type(s_sp_node), pointer, intent(inout) :: node
    
    ! local variables
    integer                   :: alloc_stat


    allocate(node, stat = alloc_stat)
    if (alloc_stat /= 0)   call error_msg_alloc

    node%type     = type
    node%priority = Priority(type)
    node%capacity = Capacity(type)

    return

  end subroutine create_node

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    append_node
  !> @brief        append additional node
  !! @authors      NT, WN
  !! @param[inout] node   : pointer to node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine append_node(node)

    ! formal arguments
    type(s_sp_node), pointer, intent(inout) :: node

    ! local variables
    type(s_sp_node), pointer  :: prv_node


    if (g_cur_node%capacity > 0) then

      ! append to current node
      !
       
      if (.not. associated(g_cur_node%child_node1)) then
        g_cur_node%child_node1 => node
      else if (.not. associated(g_cur_node%child_node2)) then
        g_cur_node%child_node2 => node
      else
        call error_msg('Sel_parser> unexpected error (1)')
      end if

      g_cur_node%capacity = g_cur_node%capacity - 1

      node%parent_node => g_cur_node

      g_cur_node => node

    else

      ! insert to ancestor node of current node
      !

      if (node%capacity == 0) then
        g_error = '"'//trim(Names(node%type))//'"'
        return
      end if

      node%capacity = node%capacity - 1

      prv_node => g_cur_node
      g_cur_node => g_cur_node%parent_node

      if (.not.associated(g_cur_node)) &
        call error_msg('Sel_parser> unexpected error (2)')

      do while(node%priority < g_cur_node%priority)

        prv_node => g_cur_node
        g_cur_node => g_cur_node%parent_node

        if (.not.associated(g_cur_node)) &
          call error_msg('sel_parser> unexpected error (3)')
      end do

      node%child_node1 => prv_node
      node%parent_node => g_cur_node
      prv_node%parent_node => node

      if (associated(prv_node, g_cur_node%child_node1)) then
        g_cur_node%child_node1 => node
      else
        g_cur_node%child_node2 => node
      end if

      g_cur_node => node

    end if

    return

  end subroutine append_node

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_node
  !> @brief        deallcate node
  !! @authors      NT, WN
  !! @param[inout] node : pointer to node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine dealloc_node(node)

    ! formal arguments
    type(s_sp_node), pointer :: node

    ! local variables
    integer          :: dealloc_stat


    if (associated(node)) then

      call dealloc_selatoms(node%selatoms)

      call dealloc_node(node%child_node1)
      call dealloc_node(node%child_node2)

      node%child_node1 => null()
      node%child_node2 => null()

      deallocate(node, stat = dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc

    end if

    return

  end subroutine dealloc_node

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_node
  !> @brief        check node
  !! @authors      NT, WN
  !! @param[in]    node : pointer to node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine check_node(node)

    ! formal arguments
    type(s_sp_node), pointer :: node


    if (.not. associated(node)) &
      return

    ! check capacity
    if (node%capacity /= 0) then
      g_error = '"'//trim(Names(node%type))//'"'
      return
    end if

    ! check bracket state
    if (node%type == NT_Bracket) then
      if (node%param_str /= '') then
        g_error = 'missing bracket: ")"'
        return
      end if
    end if

    call check_node(node%child_node1)
    call check_node(node%child_node2)

    return

  end subroutine check_node

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_node
  !> @brief        print node
  !! @authors      NT, WN
  !! @param[in]    tab_cnt : number of tab
  !! @param[in]    node    : pointer to node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine print_node(tab_cnt, node)

    ! formal arguments
    integer,                 intent(in)    :: tab_cnt
    type(s_sp_node), pointer               :: node

    ! local variables
    character(100)   :: f, str


    if (.not. associated(node)) &
      return

    write(f, '("(", i4, "x,a)")') tab_cnt+1

    if (len_trim(node%param_str) == 0) then
      str = trim(Names(node%type))
    else
      str = trim(Names(node%type))//'  ['//trim(node%param_str)//']'
    end if

    write(MsgOut, fmt=trim(f)) trim(str)

    call print_node(tab_cnt+2, node%child_node1)
    call print_node(tab_cnt+2, node%child_node2)

    return

  end subroutine print_node

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_int
  !> @brief        check if interger or not
  !! @authors      NT, WN
  !! @param[in]    str : characters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function is_int(str)

    ! return value
    logical        :: is_int

    ! formal arguments
    character(*),            intent(in) :: str

    ! local variables
    integer                  :: len, i


    len = len_trim(str)
    if (len == 0) then
      is_int = .false.
      return
    end if

    do i = 1, len
      if (str(i:i) < '0' .or. str(i:i) > '9') then
        is_int = .false.
        return
      end if
    end do

    is_int = .true.
    return

  end function is_int

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_real
  !> @brief        check if real or not
  !! @authors      NT, WN
  !! @param[in]    str : characters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function is_real(str)

    ! return value
    logical        :: is_real

    ! formal arguments
    character(*),            intent(in) :: str

    ! local variables
    integer                  :: len, i


    len = len_trim(str)
    if (len == 0) then
      is_real = .false.
      return
    end if

    do i = 1, len
      if ((str(i:i) < '0' .or. str(i:i) > '9') .and. &
           str(i:i) /= '.') then
        is_real = .false.
        return
      end if
    end do

    is_real = .true.
    return

  end function is_real

end module select_parser_mod
