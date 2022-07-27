!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   select_molecules_mod
!> @brief   subroutines for selecting molecules
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

! in older gfortran versions move_alloc is not pure
! delete pure keyword in such cases to allow compilation
#ifdef MOVE_ALLOC_IS_NOT_PURE
#define pure
#endif /* MOVE_ALLOC_IS_NOT_PURE */

module select_molecules_mod

  use constants_mod
  use one_molecule_str_mod
  use molecule_manipulate_mod
  use molecules_str_mod
  use messages_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use fileio_control_mod

  implicit none

  private
  public :: s_one_molecule
  public :: show_ctrl_molselect
  public :: read_ctrl_molselect
  public :: setup_molselect

  type, public :: s_molsel_info
    integer,              dimension(:), allocatable :: selections
    character(MaxLine),   dimension(:), allocatable :: modes
  end type s_molsel_info

  type, public :: s_selmols
    type(s_one_molecule), dimension(:), allocatable :: mols
  end type s_selmols

  character(3),           dimension(3), parameter   :: accepted_modes = &
                                                        [ "ANY", "ALL", "SET" ]

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_molselect
  !> @brief        show control parameters in MOLECULE_SELECT section
  !! @authors      DS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_molselect

    write(MsgOut,'(A)') '[MOLECULE_SELECTION]'
    write(MsgOut,'(A)') 'selection1     = 1       # group number'
    write(MsgOut,'(A)') 'mode1          = ALL     # [ALL,SET,ANY] molecule selection mode for "selection1"'
    write(MsgOut,'(A)') '                         # ALL: select all molecules in "selection1"'
    write(MsgOut,'(A)') '                         # SET: define "selection1" as a single molecule'
    write(MsgOut,'(A)') '                         # ANY: select molecules containing at least one atom in "selection1"'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_molselect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_molselect
  !> @brief        read MOLECULE_SELECT section in the control file
  !! @authors      DS
  !! @param[in]    handle      : unit number of control file
  !! @param[out]   molsel_info : control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_molselect(handle, molsel_info)

    ! parameters
    character(*),            parameter     :: Section      = 'MOLECULE_SELECTION'
    character(*),            parameter     :: default_mode = 'ALL'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_molsel_info),     intent(out)   :: molsel_info

    ! local variables
    integer                                       :: nsel, sel_num, isel
    character(MaxLine)                            :: field_name
    character(MaxLine), dimension(:), allocatable :: work_carray
    integer,            dimension(:), allocatable :: work_iarray

    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    allocate(molsel_info%selections(0), molsel_info%modes(0))
    nsel = 0
    do
      sel_num = -1
      write(field_name, '("selection", i0)') nsel+1

      call read_ctrlfile_integer(handle, Section, field_name, sel_num)

      if (sel_num < 0) then
        exit
      end if

      nsel = nsel + 1

      allocate(work_carray(nsel), work_iarray(nsel))
      work_iarray(1:nsel-1) = molsel_info%selections
      work_carray(1:nsel-1) = molsel_info%modes
      work_iarray(nsel) = sel_num
      call move_alloc(work_iarray, molsel_info%selections)

      work_carray(nsel) = default_mode
      write(field_name, '("mode", i0)') nsel
      call read_ctrlfile_string(handle, Section, field_name, work_carray(nsel))
      call move_alloc(work_carray, molsel_info%modes)

    end do

    if (nsel == 0) then
      nsel = 1
      deallocate(molsel_info%selections, molsel_info%modes)
      allocate(molsel_info%selections(1), molsel_info%modes(1))
      molsel_info%selections(1) = 0
      molsel_info%modes(1) = default_mode
    end if

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !

    write(MsgOut,'(A)') 'Read_Ctrl_MolSelect> Parameters of Molecule Select'

    do isel = 1, nsel

      write(MsgOut,'(A,I0,A,I0)') '  selection group', isel, '= ', &
        molsel_info%selections(isel)
      write(MsgOut,'(A,I0,A,A)')  '  selection mode', isel, ' = ', &
        trim(molsel_info%modes(isel))
      if (all(molsel_info%modes(isel) /= accepted_modes)) then
        call error_msg('Read_Ctrl_MolSelect> Unknown selection mode: ' &
          // trim(molsel_info%modes(isel)))
      end if

    end do

    return

  end subroutine read_ctrl_molselect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molselect
  !> @brief        setup molecule select
  !! @authors      DS
  !! @param[in]    molsel_info : molecule select information
  !! @param[in]    sel_info    : atom select information
  !! @param[in]    molecule    : molecule information
  !! @param[out]   selmols     : selected molecule array
  !! @param[out]   allmols     : selected molecule array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molselect(molsel_info, sel_info, molecule, selmols, &
    allmols)

    type(s_molsel_info),                             intent(in)  :: molsel_info
    type(s_sel_info),                                intent(in)  :: sel_info
    type(s_molecule),                                intent(in)  :: molecule
    type(s_selmols),      dimension(:), allocatable, intent(out) :: selmols
    type(s_one_molecule), dimension(:), allocatable, intent(out), &
                                                        optional :: allmols

    logical                                         :: select_when_any
    character(MaxLine)                              :: expression, mode
    integer                                         :: i, isel
    type(s_one_molecule), dimension(:), allocatable :: all_molecules

    call check_molecule_structure(molecule)

    ! a list containing all molecules to prevent redundant recalculation
    call populate_molecule_array(molecule, all_molecules)

    allocate(selmols(size(molsel_info%selections)))

    do i = 1, size(molsel_info%selections)

      isel = molsel_info%selections(i)
      mode = molsel_info%modes(i)

      if (isel == 0) then
        expression = "all"
      else
        expression = sel_info%groups(isel)
      end if

      if (mode == "ALL") then
        select_when_any = .false.
      else if (mode == "ANY") then
        select_when_any = .true.
      else
        if (mode == "SET") then
          call set_molecule(molecule, expression, selmols(i)%mols)
          cycle
        else
          call error_msg('Setup_MolSelect> Unknown selection mode: ' &
            // trim(mode))
        end if
      end if

      if (expression == "all" ) then
        allocate(selmols(i)%mols(size(all_molecules)))
        selmols(i)%mols(:) = all_molecules
      else
        call select_molecules(molecule, all_molecules, expression, &
          selmols(i)%mols, select_when_any)
      end if

    end do

    if (present(allmols)) then
      allocate(allmols(size(all_molecules)))
      allmols(:) = all_molecules
    end if

    return

  end subroutine setup_molselect

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_molecules
  !> @brief        select molecules according to atom select information
  !! @authors      DS
  !! @param[in]    molecule        : system molecule information
  !! @param[in]    all_molecules   : array containing all system molecules
  !! @param[in]    expression      : atom select expression
  !! @param[out]   selmols         : selected molecule array
  !! @param[in]    select_when_any : select molecules when any atom matches?
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_molecules(molecule, all_molecules, expression, selmols, &
    select_when_any)

    type(s_molecule),                                intent(in) :: molecule
    type(s_one_molecule), dimension(:),              intent(in) :: all_molecules
    character(*),                                    intent(in) :: expression
    type(s_one_molecule), dimension(:), allocatable, intent(out):: selmols

    logical,                             intent(in), optional :: select_when_any

    logical                                    :: select_any
    logical, dimension(:), allocatable         :: atom_is_selected
    logical, dimension(:), allocatable, target :: all_selected, &
                                                  any_selected
    logical, dimension(:), pointer             :: molecule_is_selected

    integer                                    :: imol, nmol, &
                                                  imol_select, nmol_select

    type(s_selatoms)                           :: selatoms

    if (present(select_when_any)) then
      select_any = select_when_any
    else
      select_any = .false.
    end if

    call select_atom(molecule, expression, selatoms)

    allocate(atom_is_selected(molecule%num_atoms))
    atom_is_selected = .false.
    atom_is_selected(selatoms%idx(:)) = .true.

    nmol = size(all_molecules)
    allocate(all_selected(nmol))
    allocate(any_selected(nmol))

    if (size(all_molecules) == 0) then
      call error_msg("select_molecules > empty all_molecules")
    end if

    do imol = 1, nmol
      all_selected(imol) &
        = all(atom_is_selected(all_molecules(imol)%atom_idx(:)))
      if (all_selected(imol)) then
        any_selected(imol) = .true.
      else
        any_selected(imol) &
          = any(atom_is_selected(all_molecules(imol)%atom_idx(:)))
      end if
    end do

    if (select_any) then
      molecule_is_selected => any_selected
    else
      molecule_is_selected => all_selected
    end if

    nmol_select = count(molecule_is_selected)
    allocate(selmols(nmol_select))

    imol_select = 0
    do imol = 1, nmol
      if (molecule_is_selected(imol)) then
        imol_select = imol_select + 1
        selmols(imol_select) = all_molecules(imol)
      end if
    end do

    deallocate(atom_is_selected, any_selected, all_selected)

    return

  end subroutine select_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_molecule
  !> @brief        create a single molecule structure from selection
  !! @authors      DS
  !! @param[in]    all_molecules   : system molecule information
  !! @param[in]    expression      : atom select expression
  !! @param[out]   selmols         : selected molecule array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_molecule(all_molecules, expression, selmols)

    type(s_molecule),                                intent(in) :: all_molecules
    character(*),                                    intent(in) :: expression
    type(s_one_molecule), dimension(:), allocatable, intent(out):: selmols

    type(s_selatoms)                   :: selatoms

    call select_atom(all_molecules, expression, selatoms)

    call construct_molecule_array_from_atom_array(all_molecules, &
      selatoms%idx, [size(selatoms%idx)], selmols)

    if (size(selmols) /= 1) then
      call error_msg("set_molecule > incorrect molecule_array size")
    end if

    call write_substructure_information(selmols(1))

    return

  end subroutine set_molecule

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    populate_molecule_array
  !> @brief        split molecule structure into separate molecules in a list
  !! @authors      DS
  !! @param[in]    all_molecules : system molecule information
  !! @param[inout] molecule_array : molecule list
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine populate_molecule_array(all_molecules, molecule_array)

    type(s_molecule),     intent(in) :: all_molecules
    type(s_one_molecule), dimension(:), allocatable, intent(out) &
                                     :: molecule_array

    integer, dimension(:),    allocatable :: atoms, atom_counts, bond_counts

    integer, dimension(:, :), allocatable :: bonds

    call cluster_bonds(all_molecules%bond_list, atoms=atoms, &
      atom_counts=atom_counts, bonds=bonds, bond_counts=bond_counts)

    call construct_molecule_array_from_atom_array(all_molecules, &
      atoms, atom_counts, molecule_array, bonds, bond_counts)

    return

  end subroutine populate_molecule_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_molecule_structure
  !> @brief        check if molecule structure is correct
  !! @authors      DS
  !! @param[in]    all_molecules : system molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_molecule_structure(all_molecules)

    type(s_molecule), intent(in) :: all_molecules

    integer :: natm, nbond

    natm = size(all_molecules%atom_no)
    nbond = size(all_molecules%bond_list, 2)

    if (.not. allocated(all_molecules%atom_no)) then
      call error_msg(&
        "check_molecule_structure > s_molecule%atom_no is not allocated")
    end if

    if (.not. allocated(all_molecules%bond_list)) then
      call error_msg(&
        "check_molecule_structure > s_molecule%bond_list is not allocated")
    end if

    if (natm /= size(all_molecules%atom_no)) then
      write(ErrOut, *) &
  "check_molecule_structure > num_atoms and atom_no mismatch in s_molecule"
    end if

    if (nbond /= size(all_molecules%bond_list, 2)) then
      write(ErrOut, *) &
  "check_molecule_structure > num_bonds and bond_list mismatch in s_molecule"
    end if


    return

  end subroutine check_molecule_structure

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    construct_molecule_array_from_atom_array
  !> @brief        constructs moleucle list from atom lists
  !! @authors      DS
  !! @param[in]    all_molecules  : system molecule information
  !! @param[in]    atoms          : array of atom numers
  !! @param[in]    atom_counts    : array atom count in each molecule
  !! @param[out]   molecule_array : list of molecules
  !! @param[in]    bonds          : 2d array of bonds
  !! @param[in]    atom_counts    : array of bond count in each molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine construct_molecule_array_from_atom_array(all_molecules,&
    atoms, atom_counts, molecule_array, bonds, bond_counts)

    type(s_molecule),              intent(in)                  :: all_molecules
    integer,         dimension(:), intent(in)                  :: atoms, &
                                                                  atom_counts
    type(s_one_molecule), dimension(:), allocatable, intent(out) &
                                                               :: molecule_array
    integer,         dimension(:, :),   intent(in),  optional  :: bonds
    integer,         dimension(:),      intent(in),  optional  :: bond_counts


    logical                              :: recalc_bonds
    logical, dimension(:), allocatable   :: bond_done, angl_done, dihe_done, &
                                            impr_done, cmap_done

    integer                              :: natm, iatm, &
                                            nbond, natm_in_mol, &
                                            atom_start, atom_end, &
                                            bond_start, bond_end, &
                                            i, j, cnt

    integer, dimension(:), allocatable   :: atom_idx_array

    if (present(bonds) .and. present(bond_counts)) then
      recalc_bonds = .false.
    else
      recalc_bonds = .true.
    end if

    natm = size(all_molecules%atom_no)
    nbond = size(all_molecules%bond_list, 2)

    allocate(molecule_array(size(atom_counts)))

    allocate(bond_done(nbond))
    bond_done = .false.

    if (allocated(all_molecules%angl_list)) then
      allocate(angl_done(size(all_molecules%angl_list, 2)))
      angl_done = .false.
    end if

    if (allocated(all_molecules%dihe_list)) then
      allocate(dihe_done(size(all_molecules%dihe_list, 2)))
      dihe_done = .false.
    end if

    if (allocated(all_molecules%impr_list)) then
      allocate(impr_done(size(all_molecules%impr_list, 2)))
      impr_done = .false.
    end if

    if (allocated(all_molecules%cmap_list)) then
      allocate(cmap_done(size(all_molecules%cmap_list, 2)))
      cmap_done = .false.
    end if

    atom_end = 0
    bond_end = 0
    do i = 1, size(atom_counts)

      natm_in_mol = atom_counts(i)
      atom_start = atom_end + 1
      atom_end = atom_end + natm_in_mol
      allocate(atom_idx_array(natm_in_mol))
      atom_idx_array(:) = atoms(atom_start:atom_end)

      cnt = 0
      j = 0

      call assign_integer_by_select(atom_idx_array, all_molecules%atom_no,&
        molecule_array(i)%atom_no, molecule_array(i)%num_atoms)

      call assign_character_by_select(atom_idx_array, &
        all_molecules%segment_name, molecule_array(i)%segment_name)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%segment_no, molecule_array(i)%segment_no)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%residue_no, molecule_array(i)%residue_no)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%residue_c_no, molecule_array(i)%residue_c_no)

      call assign_character_by_select(atom_idx_array, &
        all_molecules%residue_name, molecule_array(i)%residue_name)

      call assign_character_by_select(atom_idx_array, &
        all_molecules%atom_name, molecule_array(i)%atom_name)

      call assign_character_by_select(atom_idx_array, &
        all_molecules%atom_cls_name, molecule_array(i)%atom_cls_name)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%atom_cls_no, molecule_array(i)%atom_cls_no)

      call assign_real_by_select(atom_idx_array, &
        all_molecules%charge, molecule_array(i)%charge)

      call assign_real_by_select(atom_idx_array, &
        all_molecules%mass, molecule_array(i)%mass)

      call assign_real_by_select(atom_idx_array, &
        all_molecules%inv_mass, molecule_array(i)%inv_mass)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%imove, molecule_array(i)%imove)

      call assign_character_by_select(atom_idx_array, &
        all_molecules%chain_id, molecule_array(i)%chain_id)

      call assign_real_2d_by_select(atom_idx_array, &
        all_molecules%atom_coord, molecule_array(i)%atom_coord)

      call assign_real_by_select(atom_idx_array, &
        all_molecules%atom_occupancy, molecule_array(i)%atom_occupancy)

      call assign_real_by_select(atom_idx_array, &
        all_molecules%atom_temp_factor, molecule_array(i)%atom_temp_factor)

      call assign_real_2d_by_select(atom_idx_array, &
        all_molecules%atom_velocity, molecule_array(i)%atom_velocity)

      call assign_logical_by_select(atom_idx_array, &
        all_molecules%light_atom_name, molecule_array(i)%light_atom_name)

      call assign_logical_by_select(atom_idx_array, &
        all_molecules%light_atom_mass, molecule_array(i)%light_atom_mass)

      call assign_integer_by_select(atom_idx_array, &
        all_molecules%molecule_no, molecule_array(i)%molecule_no)

      if (recalc_bonds) then

        call assign_by_list(atom_idx_array, all_molecules%bond_list,&
          molecule_array(i)%bond_list, bond_done, &
          molecule_array(i)%num_bonds)

      else

        bond_start = bond_end + 1
        bond_end = bond_end + bond_counts(i)

        allocate(molecule_array(i)%bond_list(2, bond_counts(i)))
        molecule_array(i)%bond_list(:, :) = bonds(:, bond_start:bond_end)
        molecule_array(i)%num_bonds = bond_counts(i)

      end if

      ! disable assigning these lists, because it takes too much time
      ! and they are not used anywhere
      !call assign_by_list(atom_idx_array, all_molecules%angl_list,&
        !molecule_array(i)%angl_list, angl_done, &
        !molecule_array(i)%num_angles)

      !call assign_by_list(atom_idx_array, all_molecules%dihe_list,&
        !molecule_array(i)%dihe_list, dihe_done, &
        !molecule_array(i)%num_dihedrals)

      !call assign_by_list(atom_idx_array, all_molecules%impr_list,&
        !molecule_array(i)%impr_list, impr_done, &
        !molecule_array(i)%num_impropers)

      !call assign_by_list(atom_idx_array, all_molecules%cmap_list,&
        !molecule_array(i)%cmap_list, cmap_done, &
        !molecule_array(i)%num_cmaps)

      if (allocated(all_molecules%molecule_atom_no)) then

        allocate(molecule_array(i)%molecule_atom_no(1), &
                 molecule_array(i)%molecule_mass(1), &
                 molecule_array(i)%molecule_name(1))

        do j = 1, size(all_molecules%molecule_atom_no)

          iatm = all_molecules%molecule_atom_no(j)

          if (any(atom_idx_array == iatm)) then

            molecule_array(i)%molecule_atom_no(1) &
              = all_molecules%molecule_atom_no(j)
            molecule_array(i)%molecule_mass(1) = all_molecules%molecule_mass(j)
            molecule_array(i)%molecule_name(1) = all_molecules%molecule_name(j)

            exit

          end if

        end do

      end if

      molecule_array(i)%num_residues = 1
      molecule_array(i)%num_molecules = 1
      molecule_array(i)%num_segments = 1
      allocate(molecule_array(i)%atom_idx(molecule_array(i)%num_atoms))
      molecule_array(i)%atom_idx(:) = atom_idx_array
      allocate(molecule_array(i)%atom_idx2loc(&
        minval(molecule_array(i)%atom_idx):maxval(molecule_array(i)%atom_idx)))
      molecule_array(i)%atom_idx2loc = 0
      molecule_array(i)%atom_idx2loc(molecule_array(i)%atom_idx) &
        = [ (j,j=1,molecule_array(i)%num_atoms) ]
      allocate(molecule_array(i)%sub_idx(molecule_array(i)%num_atoms))
      molecule_array(i)%sub_idx = 1

      deallocate(atom_idx_array)

    end do

  end subroutine construct_molecule_array_from_atom_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_substrucure_information
  !> @brief        detect substructures inside molecule and write it
  !! @authors      DS
  !! @param[inout] molecule : s_one_molecule structure
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine write_substructure_information(molecule)

    type(s_one_molecule), intent(inout) :: molecule

    integer                                            :: isub, iatm, &
                                                          atom_idx, nsub
    integer,              dimension(:),    allocatable :: atoms, atom_counts, &
                                                          bond_counts
    integer,              dimension(:, :), allocatable :: bonds

    type(s_one_molecule), dimension(:),    allocatable :: molecule_array


    call cluster_bonds(molecule%bond_list, atoms=atoms, &
      atom_counts=atom_counts, bonds=bonds, bond_counts=bond_counts)

    call construct_molecule_array_from_atom_array(molecule%s_molecule, &
      atoms, atom_counts, molecule_array, bonds, bond_counts)

    nsub = size(molecule_array)

    do isub = 1, nsub

      do iatm = 1, molecule_array(isub)%num_atoms
        atom_idx = molecule_array(isub)%atom_idx(iatm)
        molecule%sub_idx(molecule%atom_idx2loc(atom_idx)) = isub
      end do

    end do

    molecule%num_sub = nsub

    return

  end subroutine write_substructure_information

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_real_by_select
  !> @brief        subroutine to assign subsets of 1d real arrays
  !! @authors      DS
  !! @param[in]    select_atoms : integer array with indexes of selected atoms
  !! @param[in]    old_array    : 1d real array
  !! @param[out]   new_array    : 1d real array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_real_by_select(select_atoms, old_array, new_array)

    integer,  dimension(:),              intent(in)  :: select_atoms
    real(wp), dimension(:), allocatable, intent(in)  :: old_array
    real(wp), dimension(:), allocatable, intent(out) :: new_array

    if (.not. allocated(old_array)) return

    allocate(new_array(size(select_atoms)))
    new_array(:) = old_array(select_atoms)

    return

  end subroutine assign_real_by_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_real_2d_by_select
  !> @brief        subroutine to assign subsets of 2d real arrays
  !! @authors      DS
  !! @param[in]    select_atoms : integer array with indexes of selected atoms
  !! @param[in]    old_array    : 2d real array
  !! @param[out]   new_array    : 2d real array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_real_2d_by_select(select_atoms, old_array, new_array)

    integer,  dimension(:),                 intent(in)  :: select_atoms
    real(wp), dimension(:, :), allocatable, intent(in)  :: old_array
    real(wp), dimension(:, :), allocatable, intent(out) :: new_array

    if (.not. allocated(old_array)) return

    allocate(new_array(size(old_array, 1), size(select_atoms)))
    new_array(:, :) = old_array(:, select_atoms)

    return

  end subroutine assign_real_2d_by_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_integer_by_select
  !> @brief        subroutine to assign subsets of 1d integer arrays
  !! @authors      DS
  !! @param[in]    select_atoms : integer array with indexes of selected atoms
  !! @param[in]    old_array    : 1d integer array
  !! @param[out]   new_array    : 1d integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_integer_by_select(select_atoms, old_array, new_array, &
    narray)

    integer, dimension(:),              intent(in)            :: select_atoms
    integer, dimension(:), allocatable, intent(in)            :: old_array
    integer, dimension(:), allocatable, intent(out)           :: new_array
    integer,                            intent(out), optional :: narray

    if (.not. allocated(old_array)) return

    allocate(new_array(size(select_atoms)))
    new_array(:) = old_array(select_atoms)

    if (present(narray)) narray = size(new_array)

    return

  end subroutine assign_integer_by_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_logical_by_select
  !> @brief        subroutine to assign subsets of 1d logical arrays
  !! @authors      DS
  !! @param[in]    select_atoms : integer array with indexes of selected atoms
  !! @param[in]    old_array    : 1d logical array
  !! @param[out]   new_array    : 1d logical array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_logical_by_select(select_atoms, old_array, new_array)

    integer, dimension(:),              intent(in)  :: select_atoms
    logical, dimension(:), allocatable, intent(in)  :: old_array
    logical, dimension(:), allocatable, intent(out) :: new_array

    if (.not. allocated(old_array)) return

    allocate(new_array(size(select_atoms)))
    new_array(:) = old_array(select_atoms)

    return

  end subroutine assign_logical_by_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_character_by_select
  !> @brief        subroutine to assign subsets of 1d character arrays
  !! @authors      DS
  !! @param[in]    select_atoms : integer array with indexes of selected atoms
  !! @param[in]    old_array    : 1d character array
  !! @param[out]   new_array    : 1d character array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_character_by_select(select_atoms, old_array, new_array)

    integer,      dimension(:),              intent(in)  :: select_atoms
    character(*), dimension(:), allocatable, intent(in)  :: old_array
    character(*), dimension(:), allocatable, intent(out) :: new_array

    if (.not. allocated(old_array)) return

    allocate(new_array(size(select_atoms)))
    new_array(:) = old_array(select_atoms)

    return

  end subroutine assign_character_by_select

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_by_list
  !> @brief        subroutine to assign subsets of bond, angle, etc. lists
  !! @authors      DS
  !! @param[in]    atom_array   : array of atom numbers
  !! @param[in]    old_array    : 2d integer array
  !! @param[out]   new_array    : 2d integer array
  !! @param[inout] done_array   : logical array that shows finished bonds, etc.
  !! @param[out]   narray       : number of bonds, etc. in the new array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  pure subroutine assign_by_list(atom_array, old_array, new_array, done_array, &
    narray)

    integer, dimension(:),                 intent(in)             :: atom_array
    integer, dimension(:, :), allocatable, intent(in)             :: old_array
    integer, dimension(:, :), allocatable, intent(out)            :: new_array
    logical, dimension(:),    allocatable, intent(inout)          :: done_array
    integer,                               intent(out),  optional :: narray

    logical                                               :: all_match
    integer                                               :: i, j, d1, d2, iatm
    integer, dimension(:, :), allocatable                 :: temp_i2d

    if (.not. allocated(old_array)) return

    d1 = size(old_array, 1)
    d2 = 0
    allocate(temp_i2d(d1, size(old_array, 2)))
    do i = 1, size(old_array, 2)
      if (done_array(i)) cycle
      all_match = .true.
      do j = 1, d1
        iatm = old_array(j, i)
        if (.not. any(atom_array == iatm)) then
          all_match = .false.
          exit
        end if
      end do
      if (all_match) then
        d2 = d2 + 1
        temp_i2d(:, d2) = old_array(:, i)
        done_array(i) = .true.
      end if
    end do

    allocate(new_array(d1, d2))
    new_array(:, :) = temp_i2d(:, 1:d2)

    if (present(narray)) narray = size(new_array, 2)

    return

  end subroutine assign_by_list

end module select_molecules_mod
