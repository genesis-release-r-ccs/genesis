!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_atoms_mod
!> @brief   atom selection modules
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module select_atoms_mod

  use select_contacts_mod
  use select_parser_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! variables
  real(wp),         private, save              :: g_around_scale = 1.0_wp

  ! subroutines
  public  :: select_atom
  public  :: reselect_atom
  private :: make_selatoms
  private :: make_bracket
  private :: make_and
  private :: make_or
  private :: make_not
  private :: make_around_atom
  private :: make_around_resi
  private :: make_around_mole
  private :: make_atom_name
  private :: make_atom_idx
  private :: make_atom_no
  private :: make_res_name
  private :: make_res_no
  private :: make_mol_name
  private :: make_seg_id
  private :: make_hydrogen
  private :: make_backbone
  private :: make_heavy
  private :: make_all
  private :: make_sel_hyd
  private :: calc_selatoms_union
  private :: calc_selatoms_isect
  private :: calc_selatoms_not
  private :: calc_selatoms_not_isect
  private :: quicksort
  private :: reselect_atom_
  private :: adjust_selection
  private :: adjust_selection_alt

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_atom
  !> @brief        return list of atom indexes
  !! @authors      NT, CK
  !! @param[in]    molecule   : structure of molecule information
  !! @param[in]    expression : characters for atom selection
  !! @param[out]   selatoms   : structure of selatoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_atom(molecule, expression, selatoms)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    character(*),            intent(in)    :: expression
    type(s_selatoms),        intent(inout) :: selatoms

    ! local variables
    logical                  :: ret

    type(s_sp_node), pointer :: parse_tree


    ! selection
    !
    call parser_build_tree(expression, ret)

    if (.not. ret) &
      call error_msg('Select_atom> Syntax error. '//trim(parser_get_error()))

    parse_tree => parser_get_tree()

    g_around_scale = 1.0_wp

    if (associated(parse_tree)) then
      call make_selatoms(molecule, parse_tree)
      call alloc_selatoms(selatoms, size(parse_tree%selatoms%idx))
      selatoms%idx(:) = parse_tree%selatoms%idx(:)
    else
      call alloc_selatoms(selatoms, 0)
    end if

    call parser_dealloc_tree

    return

  end subroutine select_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reselect_atom
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule     : structure of molecule information
  !! @param[in]    expression   : characters for atom selection
  !! @param[in]    coord_new    : new atom coordinates
  !! @param[in]    selatoms_org : structure of selatoms
  !! @param[out]   selatoms_new : structure of selatoms
  !! @param[in]    tool_name    : tool name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reselect_atom(molecule,     &
                           expression,   &
                           coord_new,    &
                           selatoms_org, &
                           selatoms_new, &
                           tool_name)

    ! parameters
    real(wp),                parameter     :: DecScale = 0.90_wp
    real(wp),                parameter     :: IncScale = 1.10_wp

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    character(*),            intent(in)    :: expression
    real(wp),                intent(in)    :: coord_new(:,:)
    type(s_selatoms),        intent(in)    :: selatoms_org
    type(s_selatoms),        intent(inout) :: selatoms_new
    character(len=*), optional, intent(in) :: tool_name

    ! local variables
    type(s_molecule)         :: molecule_new
    type(s_selatoms)         :: selatoms_new0
    integer                  :: i, norg, nnew, nnew0
    character(MaxLine)       :: expr0
    logical                  :: ret

    type(s_sp_node), pointer :: parse_tree


    ! check whether re-selection is necessary.
    !
    expr0 = expression
    do i = 1, len_trim(expr0)
      if (expr0(i:i) .ge. 'A' .and. expr0(i:i) .le. 'Z') then
        expr0(i:i) = char(ichar(expr0(i:i)) + 32)
      end if
    end do

    if (index(expr0, "around") == 0) then
      call alloc_selatoms(selatoms_new, size(selatoms_org%idx))
      selatoms_new%idx(:) = selatoms_org%idx(:)
      return
    end if

    norg = size(selatoms_org%idx)
    if (main_rank) then
      write(MsgOut,'(A,I7)') &
         'Reselect_Atom> number of original seleciton : ', norg
    end if

    ! define molecule by new coordinates
    !
    call alloc_molecules(molecule_new, MoleculeAtom, molecule%num_atoms)
    call alloc_molecules(molecule_new, MoleculeMole, molecule%num_molecules)

    molecule_new%num_atoms           = molecule%num_atoms
    molecule_new%atom_no(:)          = molecule%atom_no(:)
    molecule_new%segment_name(:)     = molecule%segment_name(:)
    molecule_new%segment_no(:)       = molecule%segment_no(:)
    molecule_new%residue_no(:)       = molecule%residue_no(:)
    molecule_new%residue_name(:)     = molecule%residue_name(:)
    molecule_new%atom_name(:)        = molecule%atom_name(:)
    molecule_new%molecule_no(:)      = molecule%molecule_no(:)
    molecule_new%chain_id(:)         = molecule%chain_id(:)
    molecule_new%atom_coord(:,:)     = coord_new(:,:)
    molecule_new%num_molecules       = molecule%num_molecules
    molecule_new%molecule_atom_no(:) = molecule%molecule_atom_no(:)
    molecule_new%molecule_mass(:)    = molecule%molecule_mass(:)
    molecule_new%molecule_name(:)    = molecule%molecule_name(:)
    
    ! selection
    !
    call parser_build_tree(expression, ret)

    if (.not. ret) &
      call error_msg('Reselect_Atom> Syntax error. '// &
                     trim(parser_get_error()))

    parse_tree => parser_get_tree()

    ! select (1)
    !
    g_around_scale = 1.0_wp

    if (main_rank) then
      write(MsgOut,'(A)') 'Reselect_Atom> normal selection: '
    endif

    call reselect_atom_(parse_tree, molecule_new, selatoms_new)

    nnew = size(selatoms_new%idx)

    if (main_rank) then
      write(MsgOut,'(A,I7)') &
            'Reselect_Atom> number of new selection : ', nnew
      write(MsgOut,*) ''
    end if

    ! select (2)
    !
    if (norg < nnew) then

      do while (.true.)

        g_around_scale = g_around_scale * DecScale

        if (main_rank) then
          write(MsgOut,'(A,F6.2)') &
                'Reselect_Atom> around scale:', g_around_scale
        end if

        call reselect_atom_(parse_tree, molecule_new, selatoms_new0)

        nnew0 = size(selatoms_new0%idx)

        if (main_rank) then
          write(MsgOut,'(A,I7)') &
               'Reselect_Atom> number of new selection : ', nnew0
          write(MsgOut,*) ''
        end if

        if (norg == nnew0) then

          call alloc_selatoms(selatoms_new, nnew0)
          selatoms_new%idx(:) = selatoms_new0%idx(:)
          exit

        else if (norg > nnew0) then

          if(present(tool_name)) then
            call adjust_selection_alt(norg, &
                                      molecule_new, &
                                      selatoms_new0, &
                                      selatoms_new)
          else
            call adjust_selection(norg, &
                                  molecule_new, &
                                  selatoms_new0, &
                                  selatoms_new)
          end if
          exit

        else
          
          g_around_scale = g_around_scale * (IncScale+0.1_wp)

        end if
      end do

    else if (norg > nnew) then

      do while (.true.)

        g_around_scale = g_around_scale * IncScale


        if (main_rank) then
          write(MsgOut,'(A,F6.2)') &
               'Reselect_Atom> around scale:', g_around_scale
        end if

        call reselect_atom_(parse_tree, molecule_new, selatoms_new0)

        nnew0 = size(selatoms_new0%idx)

        if (main_rank) then

          write(MsgOut,'(A,I7)') &
                'Reselect_Atom> number of new selection : ', nnew0
          write(MsgOut,*) ''
        end if

        if (norg == nnew0) then

          call alloc_selatoms(selatoms_new, nnew0)
          selatoms_new%idx(:) = selatoms_new0%idx(:)
          exit

        else if (norg < nnew0) then

          if(present(tool_name)) then
            call adjust_selection_alt(norg, &
                                      molecule_new, &
                                      selatoms_new, &
                                      selatoms_new0)
          else
            call adjust_selection(norg, &
                                  molecule_new, &
                                  selatoms_new, &
                                  selatoms_new0)
          end if
          call alloc_selatoms(selatoms_new, size(selatoms_new0%idx))
          selatoms_new%idx(:) = selatoms_new0%idx(:)
          exit

        else

          g_around_scale = g_around_scale * (DecScale-0.1_wp)

        end if
      end do

    end if

    ! deallocate
    !
    call parser_dealloc_tree
    call dealloc_selatoms(selatoms_new0)
    call dealloc_molecules(molecule_new, MoleculeAtom)
    call dealloc_molecules(molecule_new, MoleculeMole)

    if (main_rank) then
      write(MsgOut,'(a,i7)') 'Reselect_Atom> number of new selection '//&
            '(final): ', size(selatoms_new%idx)
      write(MsgOut,*) ''
    end if

    return

  end subroutine reselect_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_selatoms
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_selatoms(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    type(s_sp_node), pointer, intent(inout) :: sp_node


    select case(sp_node%type)

    case (NT_Top)

    case (NT_Bracket)
      call make_bracket(molecule, sp_node)

    case (NT_AND)
      call make_and(molecule, sp_node)

    case (NT_OR)
      call make_or(molecule, sp_node)

    case (NT_NOT)
      call make_not(molecule, sp_node)

    case (NT_AroundAtom)
      call make_around_atom(molecule, sp_node)

    case (NT_AroundRes)
      call make_around_resi(molecule, sp_node)

    case (NT_AroundMol)
      call make_around_mole(molecule, sp_node)

    case (NT_AtomName)
      call make_atom_name(molecule, sp_node)

    case (NT_AtomIdx)
      call make_atom_idx(molecule, sp_node)

    case (NT_AtomNo)
      call make_atom_no(molecule, sp_node)

    case (NT_ResName)
      call make_res_name(molecule, sp_node)

    case (NT_ResNo)
      call make_res_no(molecule, sp_node)

    case (NT_MolName)
      call make_mol_name(molecule, sp_node)

    case (NT_SegId)
      call make_seg_id(molecule, sp_node)

    case (NT_Hydrogen)
      call make_hydrogen(molecule, sp_node)

    case (NT_Backbone)
      call make_backbone(molecule, sp_node)

    case (NT_Heavy)
      call make_heavy(molecule, sp_node)

    case (NT_All)
      call make_all(molecule, sp_node)

    end select

    return

  end subroutine make_selatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_bracket
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_bracket(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    integer,         pointer :: idx(:)


    call make_selatoms(molecule, sp_node%child_node1)

    idx => sp_node%child_node1%selatoms%idx

    call alloc_selatoms(sp_node%selatoms, size(idx))
    sp_node%selatoms%idx(:) = idx(:)


    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

  end subroutine make_bracket

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_and
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_and(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node


    call make_selatoms(molecule, sp_node%child_node1)
    call make_selatoms(molecule, sp_node%child_node2)

    call calc_selatoms_isect(sp_node%child_node1%selatoms, &
                             sp_node%child_node2%selatoms, &
                             sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)
    call dealloc_selatoms(sp_node%child_node2%selatoms)

    return

  end subroutine make_and

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_or
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_or(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node


    call make_selatoms(molecule, sp_node%child_node1)
    call make_selatoms(molecule, sp_node%child_node2)

    call calc_selatoms_union(sp_node%child_node1%selatoms, &
                             sp_node%child_node2%selatoms, &
                             sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)
    call dealloc_selatoms(sp_node%child_node2%selatoms)

    return

  end subroutine make_or

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_not
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_not(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node


    call make_selatoms(molecule, sp_node%child_node1)

    call calc_selatoms_not(molecule%num_atoms, &
                           sp_node%child_node1%selatoms, &
                           sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

  end subroutine make_not

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_atom
  !> @brief        return list of atom indices
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_atom(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    type(s_selatoms)         :: selatoms_not
    real(wp)                 :: distance


    read(sp_node%param_str,*,err=900,end=900) distance

    call make_selatoms(molecule, sp_node%child_node1)

    call calc_selatoms_not(molecule%num_atoms, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not)

    if (size(sp_node%child_node1%selatoms%idx) /= 0 .and. &
        size(selatoms_not%idx) /= 0) then

      call select_contacts(SelectModeAtom, &
                           molecule, &
                           molecule, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not, &
                           distance * g_around_scale, &
                           sp_node%selatoms)

    else

      call alloc_selatoms(sp_node%selatoms, 0)

    end if

    call dealloc_selatoms(selatoms_not)
    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

900 call error_msg('Make_Around_Atom> unexpected error.')

    return

  end subroutine make_around_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_resi
  !> @brief        return list of atom indices
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_resi(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    type(s_selatoms)         :: selatoms_not
    real(wp)                 :: distance


    read(sp_node%param_str,*,err=900,end=900) distance

    call make_selatoms(molecule, sp_node%child_node1)

    call calc_selatoms_not(molecule%num_atoms, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not)

    if (size(sp_node%child_node1%selatoms%idx) /= 0 .and. &
        size(selatoms_not%idx) /= 0) then

      call select_contacts(SelectModeResi, &
                           molecule, &
                           molecule, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not, &
                           distance * g_around_scale, &
                           sp_node%selatoms)

    else

      call alloc_selatoms(sp_node%selatoms, 0)

    end if

    call dealloc_selatoms(selatoms_not)
    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

900 call error_msg('Make_Around_Resi> unexpected error.')

  end subroutine make_around_resi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_mole
  !> @brief        return list of atom indices
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_mole(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    type(s_selatoms)         :: selatoms_not
    real(wp)                 :: distance


    read(sp_node%param_str,*,err=900,end=900) distance

    call make_selatoms(molecule, sp_node%child_node1)

    call calc_selatoms_not(molecule%num_atoms, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not)

    if (size(sp_node%child_node1%selatoms%idx) /= 0 .and. &
        size(selatoms_not%idx) /= 0) then

      call select_contacts(SelectModeMole, &
                           molecule, &
                           molecule, &
                           sp_node%child_node1%selatoms, &
                           selatoms_not, &
                           distance * g_around_scale, &
                           sp_node%selatoms)

    else

      call alloc_selatoms(sp_node%selatoms, 0)

    end if

    call dealloc_selatoms(selatoms_not)
    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

900 call error_msg('Make_Around_Mole> unexpected error.')

  end subroutine make_around_mole

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_name
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_name(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    integer                  :: i, num_atm, num_sel, alloc_stat, dealloc_stat
    character(4)             :: an
    
    integer,     allocatable :: idx(:)

    
    num_atm = molecule%num_atoms

    num_sel = 0
    allocate(idx(num_atm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, num_atm
      read(molecule%atom_name(i), *, err=900) an
      if (an == sp_node%param_str) then
        num_sel = num_sel + 1
        idx(num_sel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, num_sel)
    sp_node%selatoms%idx(:) = idx(1:num_sel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

900 call error_msg('Make_atom_name> read error.')

  end subroutine make_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_idx
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_idx(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                  :: i, natm, nsel, bidx, eidx
    character(100)           :: str
    logical                  :: fromto


    str = sp_node%param_str

    natm = molecule%num_atoms

    fromto = .false.
    do i = 1, len(str)
      if (str(i:i) == '-') then
        str(i:i) = ' '
        fromto = .true.
      end if
    end do

    eidx = natm
    read(str,*,err=100,end=100) bidx, eidx

100 if (.not. fromto .or. bidx > eidx) &
      eidx = bidx

    nsel = min(natm, eidx) - max(1, bidx) + 1
    if (nsel < 0) &
      nsel = 0

    call alloc_selatoms(sp_node%selatoms, nsel)

    nsel = 0
    do i = bidx, eidx
      if (i > 0 .and. i <= natm) then
        nsel = nsel + 1
        sp_node%selatoms%idx(nsel) = i
      end if
    end do

    return

  end subroutine make_atom_idx

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_no
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_no(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: natm, nsel, bno, eno, atmno
    character(100)           :: str
    logical                  :: fromto

    integer,     allocatable :: idx(:)


    str = sp_node%param_str

    natm = molecule%num_atoms

    fromto = .false.
    do i = 1, len(str)
      if (str(i:i) == '-') then
        str(i:i) = ' '
        fromto = .true.
      end if
    end do

    eno = -1
    read(str,*,err=100,end=100) bno, eno

100 if (eno == -1) then
      if (fromto) then
        do i = 1, natm
          eno = max(molecule%atom_no(i), eno)
        end do
      else
        eno = bno
      end if
    end if

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    if (eno < bno) &
      eno = bno

    nsel = 0
    do i = 1, natm
      atmno = molecule%atom_no(i)
      if (bno <= atmno .and. atmno <= eno) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)

    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_atom_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_res_name
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_res_name(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    character(10)            :: resnam

    integer,     allocatable :: idx(:)


    resnam = sp_node%param_str

    natm = molecule%num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if (resnam == molecule%residue_name(i)) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_res_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_res_no
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_res_no(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: natm, nsel, bno, eno, resno
    character(100)           :: str
    logical                  :: fromto
    
    integer,     allocatable :: idx(:)
    

    str = sp_node%param_str

    natm = molecule%num_atoms

    fromto = .false.
    do i = 1, len(str)
      if (str(i:i) == '-') then
        str(i:i) = ' '
        fromto = .true.
      end if
    end do

    eno = -1
    read(str,*,err=100,end=100) bno, eno

100 if (eno == -1) then
      if (fromto) then
        do i = 1, natm
          eno = max(molecule%residue_no(i), eno)
        end do
      else
        eno = bno
      end if
    end if

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    if (eno < bno) &
      eno = bno

    nsel = 0
    do i = 1, natm
      resno = molecule%residue_no(i)
      if (bno <= resno .and. resno <= eno) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_res_no

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_mol_name
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_mol_name(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, j, natm, nsel, nmol, first, last
    character(20)            :: molname

    integer,     allocatable :: idx(:)


    molname = sp_node%param_str

    natm = molecule%num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc

    nmol = size(molecule%molecule_name)

    nsel = 0
    do i = 1, nmol
      if (molecule%molecule_name(i) /= molname) &
        cycle

      first = molecule%molecule_atom_no(i)
      if (i /= nmol) then
        last = molecule%molecule_atom_no(i+1) - 1
      else
        last = natm
      end if

      do j = first, last
        nsel = nsel + 1
        idx(nsel) = j
      end do
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_dealloc

    return

  end subroutine make_mol_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_res_no
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_seg_id(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    character(10)            :: segid

    integer,     allocatable :: idx(:)
    
    
    segid = sp_node%param_str

    natm = molecule%num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if (segid == adjustl(molecule%segment_name(i))) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_seg_id

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_hydrogen
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_hydrogen(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node


    call make_sel_hyd(.true., molecule, sp_node)

    return

  end subroutine make_hydrogen

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_backbone
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_backbone(molecule, sp_node)
    
    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                   :: i, natm, nsel, alloc_stat, dealloc_stat

    integer,      allocatable :: idx(:)
    character(4),     pointer :: atmnam(:)


    atmnam => molecule%atom_name

    natm = molecule%num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if (atmnam(i)(1:4) == 'N   ' .or. &
          atmnam(i)(1:4) == 'CA  ' .or. &
          atmnam(i)(1:4) == 'C   ' .or. &
          atmnam(i)(1:4) == 'O   ') then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_backbone

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_heavy
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_heavy(molecule, sp_node)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node


    call make_sel_hyd(.false., molecule, sp_node)

    return

  end subroutine make_heavy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_all
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_all(molecule, sp_node)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    integer                  :: i


    call alloc_selatoms(sp_node%selatoms, molecule%num_atoms)

    do i = 1, molecule%num_atoms
      sp_node%selatoms%idx(i) = i
    end do

    return

  end subroutine make_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_sel_hyd
  !> @brief        return list of atom indexes
  !! @authors      NT
  !! @param[in]    sel_hyd  : flag for select hydrogen or not
  !! @param[in]    molecule : structure of molecule information
  !! @param[inout] sp_node  : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_sel_hyd(sel_hyd, molecule, sp_node)

    ! formal arguments
    logical,                  intent(in)    :: sel_hyd
    type(s_molecule), target, intent(in)    :: molecule
    type(s_sp_node), pointer                :: sp_node

    ! local variables
    integer                   :: i, natm, nsel, alloc_stat, dealloc_stat

    integer,      allocatable :: idx(:)
    character(4),     pointer :: atmnam(:)


    atmnam => molecule%atom_name

    natm = molecule%num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if ((atmnam(i)(1:1) .eq. 'H' .and. sel_hyd) .or. &
          (atmnam(i)(1:1) .ne. 'H' .and. .not. sel_hyd)) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(1:nsel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_sel_hyd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_union
  !> @brief        operate : c = a or b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_union(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer                  :: alloc_stat, numc, i, ia, ib, v
    integer                  :: numa, numb
    integer                  :: dealloc_stat

    integer,     allocatable :: aidx(:), bidx(:), cidx(:)


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = numa + numb

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)

    if (alloc_stat /= 0) &
      call error_msg_alloc

    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort(aidx, 1, numa)
    if (numb > 1) &
      call quicksort(bidx, 1, numb)

    ia = 1
    ib = 1
    numc = 0

    do while (ia <= numa .and. ib <= numb)

      if (aidx(ia) < bidx(ib)) then
        v = aidx(ia)
        ia = ia + 1
      else if(aidx(ia) > bidx(ib)) then
        v = bidx(ib)
        ib = ib + 1
      else
        v = aidx(ia)
        ia = ia + 1
        ib = ib + 1
      end if

      numc = numc + 1
      cidx(numc) = v

    end do

    do i = ia, numa
      numc = numc + 1
      cidx(numc) = aidx(i)
    end do

    do i = ib, numb
      numc = numc + 1
      cidx(numc) = bidx(i)
    end do

    call alloc_selatoms(selatoms_c, numc)
    selatoms_c%idx(1:numc) = cidx(1:numc)

    deallocate(aidx, bidx, cidx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine calc_selatoms_union

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_isect
  !> @brief        operate : c = a and b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_isect(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer                  :: alloc_stat, dealloc_stat, numc, ia, ib
    integer                  :: numa, numb

    integer,     allocatable :: aidx(:), bidx(:), cidx(:)


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = min(numa, numb)

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)

    if (alloc_stat /= 0) &
      call error_msg_alloc

    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort(aidx, 1, numa)
    if (numb > 1) &
      call quicksort(bidx, 1, numb)

    ia = 1
    ib = 1
    numc = 0

    do while (ia <= numa .and. ib <= numb)

      if (aidx(ia) < bidx(ib)) then
        ia = ia + 1
      else if(aidx(ia) > bidx(ib)) then
        ib = ib + 1
      else
        numc = numc + 1
        cidx(numc) = aidx(ia)
        ia = ia + 1
        ib = ib + 1
      end if

    end do

    call alloc_selatoms(selatoms_c, numc)

    selatoms_c%idx(1:numc) = cidx(1:numc)
    deallocate(aidx, bidx, cidx, stat = dealloc_stat)

    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine calc_selatoms_isect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_not
  !> @brief        operate : b = not a
  !! @authors      NT
  !! @param[in]    num_atoms  : number of atoms
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[inout] selatoms_b : list b of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_not(num_atoms, selatoms_a, selatoms_b)

    ! formal arguments
    integer,                  intent(in)    :: num_atoms
    type(s_selatoms),         intent(in)    :: selatoms_a
    type(s_selatoms), target, intent(inout) :: selatoms_b

    ! local variables
    integer                   :: i, alloc_stat, dealloc_stat
    integer                   :: numb, iatm, natm
    integer                   :: numa

    integer,      allocatable :: aidx(:)
    integer,          pointer :: bidx(:)


    numa = size(selatoms_a%idx)

    if (numa /= 0) then
      do i = 1, numa
        natm = max(selatoms_a%idx(i), num_atoms)
      end do
    else
      natm = num_atoms
    end if

    numb = natm - numa

    allocate(aidx(numa), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    aidx(1:numa) = selatoms_a%idx(1:numa)

    if (numa > 1) &
      call quicksort(aidx, 1, numa)

    call alloc_selatoms(selatoms_b, numb)
    bidx => selatoms_b%idx

    iatm = 1
    numb = 0

    do i = 1, natm

      if (iatm <= numa) then

        if (i == aidx(iatm)) then
          iatm = iatm + 1
        else
          numb = numb + 1
          bidx(numb) = i
        end if

      else
        numb = numb + 1
        bidx(numb) = i
      end if

    end do

    deallocate(aidx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine calc_selatoms_not

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_not_isect
  !> @brief        operate : c = a not b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_not_isect(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer                  :: alloc_stat, dealloc_stat, numc, i, ia, ib
    integer                  :: numa, numb

    integer,     allocatable :: aidx(:), bidx(:), cidx(:)


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = numa

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort(aidx, 1, numa)
    if (numb > 1) &
      call quicksort(bidx, 1, numb)

    ia = 1
    ib = 1
    numc = 0

    do while (ia <= numa .and. ib <= numb)

      if (aidx(ia) < bidx(ib)) then
        numc = numc + 1
        cidx(numc) = aidx(ia)
        ia = ia + 1
      else if (aidx(ia) > bidx(ib)) then
        ib = ib + 1
      else
        ia = ia + 1
        ib = ib + 1
      end if

    end do

    do i = ia, numa
      numc = numc + 1
      cidx(numc) = aidx(i)
    end do

    call alloc_selatoms(selatoms_c, numc)
    selatoms_c%idx(1:numc) = cidx(1:numc)

    deallocate(aidx, bidx, cidx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine calc_selatoms_not_isect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort
  !> @brief        sorting function
  !! @authors      NT
  !! @param[inout] a     : list of integers
  !! @param[in]    start : start index
  !! @param[in]    end   : end index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine quicksort(a, start, end)

    ! formal arguments
    integer,                 intent(inout) :: a(*)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    integer                  :: i, j, t, x


    x = a((start + end) / 2)
    i = start
    j = end

    do
      do while (a(i) < x)
        i = i + 1
      end do

      do while (x < a(j))
        j = j - 1
      end do

      if (i >= j) exit

      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call quicksort(a, start, i - 1)
    if (j + 1 < end) &
      call quicksort(a, j + 1, end)

    return

  end subroutine quicksort

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reselect_atom_
  !> @brief        reselect_atom_
  !! @authors      NT
  !! @param[inout] parse_tree   : structure of parse tree node
  !! @param[in]    molecule_new : structure of molecule information
  !! @param[inout] selatoms_new : structure of selatoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reselect_atom_(parse_tree,   &
                            molecule_new, &
                            selatoms_new)

    ! formal arguments
    type(s_sp_node),  pointer, intent(inout) :: parse_tree
    type(s_molecule),          intent(in)    :: molecule_new
    type(s_selatoms),          intent(inout) :: selatoms_new
    

    if (associated(parse_tree)) then
      call make_selatoms(molecule_new, parse_tree)
      call alloc_selatoms(selatoms_new, size(parse_tree%selatoms%idx))
      selatoms_new%idx(:) = parse_tree%selatoms%idx(:)
    else
      call alloc_selatoms(selatoms_new, 0)
    end if

    return

  end subroutine reselect_atom_

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    adjust_selection
  !> @brief        adjust_selection
  !! @authors      NT
  !! @param[in]    target_count : target count
  !! @param[in]    molecule     : structure of molecule information
  !! @param[in]    sa_base      : structure of base selatoms
  !! @param[inout] sa_target    : structure of target selatoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine adjust_selection(target_count, &
                              molecule, &
                              sa_base, &
                              sa_target)

    ! parameter
    real(wp),                parameter     :: AdjustDist = 10.0_wp

    ! formal arguments
    integer,                 intent(in)    :: target_count
    type(s_molecule),        intent(in)    :: molecule
    type(s_selatoms),        intent(in)    :: sa_base
    type(s_selatoms),        intent(inout) :: sa_target

    ! local variables
    type(s_selatoms)         :: sa_adjust1, sa_adjust2


    ! calculate adjust = target - base
    !
    call calc_selatoms_not_isect(sa_target, sa_base, sa_adjust1)


    ! return value : sa_target
    !
    call alloc_selatoms(sa_target, size(sa_base%idx))
    sa_target%idx(:) = sa_base%idx(:)


    ! by molecule
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Adjust_Selection> adjust selection by molecule '
    end if

    call select_contacts(SelectModeMole, &
                         molecule,   &
                         molecule,   &
                         sa_base,    &
                         sa_adjust1, &
                         AdjustDist, &
                         sa_adjust2, &
                         target_count - size(sa_target%idx), &
                         .true.)

    if (main_rank) then
      write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ',  &
                   size(sa_adjust2%idx)
    end if

    call calc_selatoms_union(sa_target, sa_adjust2, sa_target)

    if (target_count == size(sa_target%idx)) &
      return
       
    call calc_selatoms_not_isect(sa_adjust1, sa_adjust2, sa_adjust1)


    ! by residue
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Adjust_Selection> adjust selection by residue : '
    end if

    call select_contacts(SelectModeResi, &
                         molecule,   &
                         molecule,   &
                         sa_base,    &
                         sa_adjust1, &
                         AdjustDist, &
                         sa_adjust2, &
                         target_count - size(sa_target%idx), &
                         .true.)

    if (main_rank) then
      write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ', &
                   size(sa_adjust2%idx)
    end if

    call calc_selatoms_union(sa_target, sa_adjust2, sa_target)

    if (target_count == size(sa_target%idx)) &
      return

    call calc_selatoms_not_isect(sa_adjust1, sa_adjust2, sa_adjust1)


    ! by atom
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Adjust_Selection> adjust selection by atom    : '
    end if

    call select_contacts(SelectModeAtom, &
                         molecule,   &
                         molecule,   &
                         sa_base,    &
                         sa_adjust1, &
                         AdjustDist, &
                         sa_adjust2, &
                         target_count - size(sa_target%idx), &
                         .true.)

    if (main_rank) then
      write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ', &
           size(sa_adjust2%idx)
    end if

    call calc_selatoms_union(sa_target, sa_adjust2, sa_target)

    if (target_count /= size(sa_target%idx)) &
      call error_msg('Adjust_Selection> ERROR: cannot adjust selection.')

    return

  end subroutine adjust_selection


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    adjust_selection_alt
  !> @brief        alter adjust_selection not to use selection by atom
  !! @authors      NT, KYMD
  !! @param[in]    target_count : target count
  !! @param[in]    molecule     : structure of molecule information
  !! @param[in]    sa_base      : structure of base selatoms
  !! @param[inout] sa_target    : structure of target selatoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine adjust_selection_alt(target_count, &
                                  molecule, &
                                  sa_base, &
                                  sa_target)

    ! parameter
    real(wp),                parameter     :: AdjustDist = 10.0_wp

    ! formal arguments
    integer,                 intent(in)    :: target_count
    type(s_molecule),        intent(in)    :: molecule
    type(s_selatoms),        intent(in)    :: sa_base
    type(s_selatoms),        intent(inout) :: sa_target

    ! local variables
    type(s_selatoms)         :: sa_adjust1, sa_adjust2,       &
                                sa_target_org,                &
                                sa_target_tmp, sa_adjust2_tmp, sa_adjust1_tmp
    integer                  :: i, resno_now, resno_next, nadj2


    ! copy
    ! 
    call alloc_selatoms(sa_target_org, size(sa_target%idx))
    sa_target_org%idx(:) = sa_target%idx(:)                      ! supreme number

    ! calculate adjust = target - base
    !
    call calc_selatoms_not_isect(sa_target, sa_base, sa_adjust1) ! adj1=sup-base

    ! return value : sa_target
    !
    call alloc_selatoms(sa_target, size(sa_base%idx))            ! copy base to tar
    sa_target%idx(:) = sa_base%idx(:)

    ! by molecule
    !
    if (main_rank) then
      write(MsgOut,'(A)')  &
        'Adjust_Selection> adjust selection by molecule (make base):'
    end if

    call select_contacts(SelectModeMole, &
                         molecule,   &
                         molecule,   &
                         sa_base,    &
                         sa_adjust1, &
                         AdjustDist, &
                         sa_adjust2, &
                         target_count - size(sa_target%idx), &
                         .true.)

    call calc_selatoms_union(sa_target, sa_adjust2, sa_target)

    if (target_count == size(sa_target%idx)) then
      if (main_rank) then
        write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ',  &
                     size(sa_adjust2%idx)
      end if

      call dealloc_selatoms(sa_adjust1)
      call dealloc_selatoms(sa_adjust2)
      call dealloc_selatoms(sa_target_org)
      return
    end if
       

    nadj2 = size(sa_adjust2%idx)
    if (nadj2 == 0) &
      call error_msg('Adjust_Selection> ERROR: cannot adjust selection.') 
    do i = 1, nadj2 - 1
      resno_now   = molecule%residue_no(sa_adjust2%idx(i))
      resno_next  = molecule%residue_no(sa_adjust2%idx(i+1))
      if(resno_now == resno_next) cycle
      call alloc_selatoms(sa_adjust2_tmp, i) 
      sa_adjust2_tmp%idx(1: i) = sa_adjust2%idx(1: i)

      call calc_selatoms_union(sa_base, sa_adjust2_tmp, sa_target_tmp)         ! base + adj2

      call calc_selatoms_not_isect(sa_adjust1, sa_adjust2_tmp, sa_adjust1_tmp) ! adj1=adj1-adj2

      ! by residue
      !
      if (main_rank) then
        write(MsgOut,'(A, i0, A)') &
          'Adjust_Selection> adjust selection by residue (base atoms: ', i, '): '
      end if

      call select_contacts(SelectModeResi, &
                           molecule,       &
                           molecule,       &
                           sa_base,        &
                           sa_adjust1_tmp, &
                           AdjustDist,     &
                           sa_adjust2_tmp, &
                           target_count - size(sa_target_tmp%idx), &
                           .true.)


      call calc_selatoms_union(sa_target_tmp, sa_adjust2_tmp, sa_target_tmp)

      if (target_count == size(sa_target_tmp%idx)) then
        if (main_rank) then
          write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ',  &
                       i + size(sa_adjust2_tmp%idx)
          !KYMD  write(MsgOut, *)  sa_adjust2_tmp%idx(1: size(sa_adjust2_tmp%idx))
        end if
      
        call alloc_selatoms(sa_target, size(sa_target_tmp%idx))
        sa_target%idx(:) = sa_target_tmp%idx(:)

        call dealloc_selatoms(sa_adjust1)
        call dealloc_selatoms(sa_adjust2)
        call dealloc_selatoms(sa_target_org)
        call dealloc_selatoms(sa_target_tmp)
        call dealloc_selatoms(sa_adjust1_tmp)
        call dealloc_selatoms(sa_adjust2_tmp)
        return
      end if

      if(i + size(sa_adjust2_tmp%idx) == nadj2) then ! res select = mol select
        if (main_rank) then
          write(MsgOut,'(1x, A, i0, A/, A)')  &
           'Adjust_Selection> Any combination of species near the border cannot be ', target_count, '.', &
           '                  Please try to change around radius, frame number, or reffile.'
        end if
        call error_msg('Adjust_Selection> ERROR: Cannot adjust selection.')
      end if

    end do 


    !
    ! reverse order of selection type
    !
    
    call dealloc_selatoms(sa_adjust1)
    call dealloc_selatoms(sa_adjust2)

    ! calculate adjust = target - base (INITIALIZE)
    !
    call calc_selatoms_not_isect(sa_target_org, sa_base, sa_adjust1)

    ! return value : sa_target         (INITIALIZE)
    !
    call alloc_selatoms(sa_target, size(sa_base%idx))
    sa_target%idx(:) = sa_base%idx(:)

    ! by residue
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
        'Adjust_Selection> adjust selection by residue (Make base): '
    end if

    call select_contacts(SelectModeResi, &
                         molecule,   &
                         molecule,   &
                         sa_base,    &
                         sa_adjust1, &
                         AdjustDist, &
                         sa_adjust2, &
                         target_count - size(sa_target%idx), &
                         .true.)

    call calc_selatoms_union(sa_target, sa_adjust2, sa_target) ! base + adj2

    if (target_count == size(sa_target%idx)) then
      if (main_rank) then
        write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ', &
                     size(sa_adjust2%idx)
        !KYMD write(MsgOut, *) sa_adjust2%idx(1: size(sa_adjust2%idx))
      end if

      call dealloc_selatoms(sa_adjust1)
      call dealloc_selatoms(sa_adjust2)
      call dealloc_selatoms(sa_target_org)
      call dealloc_selatoms(sa_target_tmp)
      call dealloc_selatoms(sa_adjust1_tmp)
      call dealloc_selatoms(sa_adjust2_tmp)
      return
    end if


    nadj2 = size(sa_adjust2%idx)
    if (nadj2 == 0) &
      call error_msg('Adjust_Selection> ERROR: cannot adjust selection. ')
    do i = 1, nadj2 - 1
      resno_now  = molecule%residue_no(sa_adjust2%idx(i))
      resno_next = molecule%residue_no(sa_adjust2%idx(i+1))
      if(resno_now == resno_next) cycle
      call alloc_selatoms(sa_adjust2_tmp, i) 
      sa_adjust2_tmp%idx(1: i) = sa_adjust2%idx(1: i)

      call calc_selatoms_union(sa_base, sa_adjust2_tmp, sa_target_tmp) ! base + adj2

      call calc_selatoms_not_isect(sa_adjust1, sa_adjust2_tmp, sa_adjust1_tmp) ! adj1=adj1-adj2

      ! by molecule
      !
      if (main_rank) then
        write(MsgOut,'(A, i0, A)') &
          'Adjust_Selection> adjust selection by molecule (base atoms: ', i, '):'
      end if

      call select_contacts(SelectModeMole, &
                           molecule,       &
                           molecule,       &
                           sa_base,        &
                           sa_adjust1_tmp, &
                           AdjustDist,     &
                           sa_adjust2_tmp, &
                           target_count - size(sa_target_tmp%idx), &
                           .true.)


      call calc_selatoms_union(sa_target_tmp, sa_adjust2_tmp, sa_target_tmp)

      if (target_count == size(sa_target_tmp%idx)) then
        if (main_rank) then
          write(MsgOut,'(A,I7)') 'Adjust_Selection> adjusted : ',  &
                       i + size(sa_adjust2_tmp%idx)
          !KYMD write(MsgOut, *) sa_adjust2_tmp%idx(1: size(sa_adjust2_tmp%idx))
        end if
      
        call alloc_selatoms(sa_target, size(sa_target_tmp%idx))
        sa_target%idx(:) = sa_target_tmp%idx(:)

        call dealloc_selatoms(sa_adjust1)
        call dealloc_selatoms(sa_adjust2)
        call dealloc_selatoms(sa_target_org)
        call dealloc_selatoms(sa_target_tmp)
        call dealloc_selatoms(sa_adjust1_tmp)
        call dealloc_selatoms(sa_adjust2_tmp)
        return
      end if

    end do 


    call error_msg('Adjust_Selection> ERROR: cannot adjust selection. '//&
                   'Please change around radius and/or reffile.')

    return

  end subroutine adjust_selection_alt

end module select_atoms_mod
