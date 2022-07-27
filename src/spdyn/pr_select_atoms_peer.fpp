!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_select_atoms_peer_mod
!> @brief   peer module of lib/select_atoms_mod
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_select_atoms_peer_mod

  use pr_huge_molecule_mod
  use select_parser_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: select_atom_peer
  private :: make_selatoms_peer
  private :: make_bracket_peer
  private :: make_and_peer
  private :: make_or_peer
  private :: make_not_peer
  private :: make_around_atom_peer
  private :: make_around_resi_peer
  private :: make_around_mole_peer
  private :: make_atom_name_peer
  private :: make_atom_idx_peer
  private :: make_atom_no_peer
  private :: make_res_name_peer
  private :: make_res_no_peer
  private :: make_mol_name_peer
  private :: make_seg_id_peer
  private :: make_hydrogen_peer
  private :: make_backbone_peer
  private :: make_heavy_peer
  private :: make_all_peer
  private :: make_sel_hyd_peer
  private :: calc_selatoms_union_peer
  private :: calc_selatoms_isect_peer
  private :: calc_selatoms_not_peer
  private :: calc_selatoms_not_isect_peer
  private :: quicksort_peer

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_atom_peer
  !> @brief        peer subroutine of lib/select_atoms.fpp::select_atom()
  !! @authors      NT
  !! @param[in]    expression : characters for atom selection
  !! @param[inout] selatoms   : structure of selatoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_atom_peer(expression, selatoms)

    ! formal arguments
    character(*),            intent(in)    :: expression
    type(s_selatoms),        intent(inout) :: selatoms

    ! local variables
    logical                  :: ret
    type(s_sp_node), pointer :: parse_tree


    ! selection
    !
    call parser_build_tree(expression, ret)

    if (.not. ret) &
      call error_msg('Select_atom_Peer> Syntax error. '//&
                     trim(parser_get_error()))

    parse_tree => parser_get_tree()

    if (associated(parse_tree)) then
      call make_selatoms_peer(parse_tree)
      call alloc_selatoms(selatoms, size(parse_tree%selatoms%idx))
      selatoms%idx(:) = parse_tree%selatoms%idx(:)
    else
      call alloc_selatoms(selatoms, 0)
    end if

    call parser_dealloc_tree

    return

  end subroutine select_atom_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_selatoms_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_selatoms()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_selatoms_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer, intent(inout) :: sp_node


    select case(sp_node%type)

    case(NT_Top)

    case(NT_Bracket)
      call make_bracket_peer(sp_node)

    case(NT_AND)
      call make_and_peer(sp_node)

    case(NT_OR)
      call make_or_peer(sp_node)

    case(NT_NOT)
      call make_not_peer(sp_node)

    case (NT_AroundAtom)
       call make_around_atom_peer(sp_node)

    case (NT_AroundRes)
       call make_around_resi_peer(sp_node)

    case (NT_AroundMol)
       call make_around_mole_peer(sp_node)

    case(NT_AtomName)
      call make_atom_name_peer(sp_node)

    case(NT_AtomIdx)
      call make_atom_idx_peer(sp_node)

    case(NT_AtomNo)
      call make_atom_no_peer(sp_node)

    case(NT_ResName)
      call make_res_name_peer(sp_node)

    case(NT_ResNo)
      call make_res_no_peer(sp_node)

    case (NT_MolName)
       call make_mol_name_peer(sp_node)

    case(NT_SegId)
      call make_seg_id_peer(sp_node)

    case(NT_Hydrogen)
      call make_hydrogen_peer(sp_node)

    case(NT_Backbone)
      call make_backbone_peer(sp_node)

    case(NT_Heavy)
      call make_heavy_peer(sp_node)

    case(NT_All)
      call make_all_peer(sp_node)

    end select

    return

  end subroutine make_selatoms_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_bracket_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_bracket()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_bracket_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    integer,         pointer :: idx(:)


    call make_selatoms_peer(sp_node%child_node1)

    idx => sp_node%child_node1%selatoms%idx

    call alloc_selatoms(sp_node%selatoms, size(idx))
    sp_node%selatoms%idx(:) = idx(:)


    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

  end subroutine make_bracket_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_and_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_and()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_and_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    call make_selatoms_peer(sp_node%child_node1)
    call make_selatoms_peer(sp_node%child_node2)

    call calc_selatoms_isect_peer(sp_node%child_node1%selatoms, &
                                  sp_node%child_node2%selatoms, &
                                  sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)
    call dealloc_selatoms(sp_node%child_node2%selatoms)

    return

  end subroutine make_and_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_or_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_or()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_or_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    call make_selatoms_peer(sp_node%child_node1)
    call make_selatoms_peer(sp_node%child_node2)

    call calc_selatoms_union_peer(sp_node%child_node1%selatoms, &
                                  sp_node%child_node2%selatoms, &
                                  sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)
    call dealloc_selatoms(sp_node%child_node2%selatoms)

    return

  end subroutine make_or_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_not_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_not()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_not_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    call make_selatoms_peer(sp_node%child_node1)

    call calc_selatoms_not_peer(hm_num_atoms, &
                                sp_node%child_node1%selatoms, &
                                sp_node%selatoms)

    call dealloc_selatoms(sp_node%child_node1%selatoms)

    return

  end subroutine make_not_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_atom_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_around_atom()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_atom_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    !NOT supported yet.

    return

  end subroutine make_around_atom_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_resi_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_around_resi()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_resi_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    !NOT supported yet.

    return

  end subroutine make_around_resi_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_around_mole_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_around_mole()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_around_mole_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    !NOT suppprted yet.

    return

  end subroutine make_around_mole_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_name_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_atom_name()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_name_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    character(4)             :: raw_an, an
    integer,     allocatable :: idx(:)
    integer                  :: i, num_atm, num_sel, alloc_stat, dealloc_stat
    
    
    num_atm = hm_num_atoms

    num_sel = 0
    allocate(idx(num_atm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, num_atm
      raw_an = hm_atom_name(i)
      read(raw_an, *, err=900) an
      if (an == sp_node%param_str) then
        num_sel = num_sel + 1
        idx(num_sel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, num_sel)
    sp_node%selatoms%idx(1:num_sel) = idx(1:num_sel)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

900 call error_msg('Make_atom_name_Peer> read error.')

  end subroutine make_atom_name_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_idx_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_atom_idx()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_idx_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    logical                  :: fromto
    character(MaxLine)       :: str
    integer                  :: i, natm, nsel, bidx, eidx


    str = sp_node%param_str

    natm = hm_num_atoms

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

  end subroutine make_atom_idx_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_atom_no_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_atom_no()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_atom_no_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer  :: sp_node

    ! local variables
    logical                   :: fromto
    character(MaxLine)        :: str
    integer                   :: i, alloc_stat, dealloc_stat
    integer                   :: natm, nsel, bno, eno, atmno
    integer,      allocatable :: idx(:)


    str = sp_node%param_str

    natm = hm_num_atoms

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
          eno = max(hm_atom_no(i), eno)
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
      atmno = hm_atom_no(i)
      if (bno <= atmno .and. atmno <= eno) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_atom_no_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_res_name_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_res_name()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_res_name_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    character(10)            :: resnam
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    integer,     allocatable :: idx(:)


    resnam = sp_node%param_str

    natm = hm_num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if (resnam == hm_residue_name(i)) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_res_name_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_res_no_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_res_no()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_res_no_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    logical                  :: fromto
    character(MaxLine)       :: str
    integer                  :: i, alloc_stat, dealloc_stat
    integer                  :: natm, nsel, bno, eno, resno
    integer,     allocatable :: idx(:)
    
    
    str = sp_node%param_str

    natm = hm_num_atoms

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
          eno = max(hm_residue_no(i), eno)
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
      resno = hm_residue_no(i)
      if (bno <= resno .and. resno <= eno) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_res_no_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_mol_name_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_mol_name()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_mol_name_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    !NOT supported yet.

    return

  end subroutine make_mol_name_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_seg_id_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_seg_id()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_seg_id_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    character(10)            :: segid
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    integer,     allocatable :: idx(:)
    
    
    segid = sp_node%param_str

    natm = hm_num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      if (segid == adjustl(hm_segment_name(i))) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_seg_id_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_hydrogen_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_hydrogen()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_hydrogen_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    call make_sel_hyd_peer(.true., sp_node)

    return

  end subroutine make_hydrogen_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_backbone_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_backbone()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_backbone_peer(sp_node)
    
    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    integer,     allocatable :: idx(:)
    character(4)             :: an


    natm = hm_num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc

    nsel = 0
    do i = 1, natm
      an = hm_atom_name(i)
      if (an(1:4) == 'N   ' .or. &
          an(1:4) == 'CA  ' .or. &
          an(1:4) == 'C   ' .or. &
          an(1:4) == 'O   ') then

        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    do i = 1, nsel
      sp_node%selatoms%idx(i) = idx(i)
    end do
!   sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_backbone_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_heavy_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_heavy()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_heavy_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node


    call make_sel_hyd_peer(.false., sp_node)

    return

  end subroutine make_heavy_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_all_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_all()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_all_peer(sp_node)

    ! formal arguments
    type(s_sp_node), pointer :: sp_node

    ! local variables
    integer                  :: i


    call alloc_selatoms(sp_node%selatoms, hm_num_atoms)

    do i = 1, hm_num_atoms
      sp_node%selatoms%idx(i) = i
    end do

    return

  end subroutine make_all_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    make_sel_hyd_peer
  !> @brief        peer routine of lib/select_atoms_mod::make_sel_hyd()
  !! @authors      NT
  !! @param[inout] sp_node : structure of parse tree node
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine make_sel_hyd_peer(sel_hyd, sp_node)

    ! formal arguments
    logical,                 intent(in)    :: sel_hyd
    type(s_sp_node), pointer               :: sp_node

    ! local variables
    integer                  :: i, natm, nsel, alloc_stat, dealloc_stat
    integer,     allocatable :: idx(:)
    character(4)             :: an


    natm = hm_num_atoms

    allocate(idx(natm), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    nsel = 0
    do i = 1, natm
      an = hm_atom_name(i)
      call toupper(an)
      if ((an(1:1) == 'H' .and.       sel_hyd) .or. &
          (an(1:1) /= 'H' .and. .not. sel_hyd)) then
        nsel = nsel + 1
        idx(nsel) = i
      end if
    end do

    call alloc_selatoms(sp_node%selatoms, nsel)
    sp_node%selatoms%idx(:) = idx(:)

    deallocate(idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    return

  end subroutine make_sel_hyd_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_union_peer
  !> @brief        operate : c = a or b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_union_peer(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer,     allocatable :: aidx(:), bidx(:), cidx(:)
    integer                  :: alloc_stat, numc, i, ia, ib, v
    integer                  :: numa, numb
    integer                  :: dealloc_stat


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = numa + numb

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort_peer(aidx, 1, numa)
    if (numb > 1) &
      call quicksort_peer(bidx, 1, numb)

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
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine calc_selatoms_union_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_isect_peer
  !> @brief        operate : c = a and b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_isect_peer(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer,     allocatable :: aidx(:), bidx(:), cidx(:)
    integer                  :: alloc_stat, dealloc_stat, numc, ia, ib
    integer                  :: numa, numb


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = min(numa, numb)

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort_peer(aidx, 1, numa)
    if (numb > 1) &
      call quicksort_peer(bidx, 1, numb)

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

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine calc_selatoms_isect_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_not_peer
  !> @brief        operate : b = not a
  !! @authors      NT
  !! @param[in]    num_atoms  : number of atoms
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[inout] selatoms_b : list b of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_not_peer(num_atoms, selatoms_a, selatoms_b)

    ! formal arguments
    integer,                  intent(in)    :: num_atoms
    type(s_selatoms),         intent(in)    :: selatoms_a
    type(s_selatoms), target, intent(inout) :: selatoms_b

    ! local variables
    integer,      allocatable :: aidx(:)
    integer,      pointer     :: bidx(:)
    integer                   :: i, alloc_stat, dealloc_stat
    integer                   :: numb, iatm, natm
    integer                   :: numa

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
    if (alloc_stat /= 0)   call error_msg_alloc
    aidx(1:numa) = selatoms_a%idx(1:numa)

    if (numa > 1) &
      call quicksort_peer(aidx, 1, numa)

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
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine calc_selatoms_not_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_selatoms_not_isect_peer
  !> @brief        operate : c = a not b
  !! @authors      NT
  !! @param[in]    selatoms_a : list a of selected atom index
  !! @param[in]    selatoms_b : list b of selected atom index
  !! @param[inout] selatoms_c : list c of selected atom index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_selatoms_not_isect_peer(selatoms_a, selatoms_b, selatoms_c)

    ! formal arguments
    type(s_selatoms),        intent(in)    :: selatoms_a
    type(s_selatoms),        intent(in)    :: selatoms_b
    type(s_selatoms),        intent(inout) :: selatoms_c

    ! local variables
    integer,     allocatable :: aidx(:), bidx(:), cidx(:)
    integer                  :: alloc_stat, dealloc_stat, numc, i, ia, ib
    integer                  :: numa, numb


    numa = size(selatoms_a%idx)
    numb = size(selatoms_b%idx)
    numc = numa

    allocate(aidx(numa), bidx(numb), cidx(numc), stat = alloc_stat)
    if (alloc_stat /= 0)   call error_msg_alloc
    aidx(1:numa) = selatoms_a%idx(1:numa)
    bidx(1:numb) = selatoms_b%idx(1:numb)

    if (numa > 1) &
      call quicksort_peer(aidx, 1, numa)
    if (numb > 1) &
      call quicksort_peer(bidx, 1, numb)

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
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine calc_selatoms_not_isect_peer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort_peer
  !> @brief        sorting function
  !! @authors      NT
  !! @param[inout] a     : list of integers
  !! @param[in]    start : start position
  !! @param[in]    end   : end position
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine quicksort_peer(a, start, end)

    ! formal arguments
    integer,                 intent(inout) ::  a(*)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    integer                :: i, j, t, x


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

    if (start < i - 1) call quicksort_peer(a, start, i - 1)
    if (j + 1 < end)  call quicksort_peer(a, j + 1, end)

    return

  end subroutine quicksort_peer

end module pr_select_atoms_peer_mod
