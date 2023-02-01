!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   molecules_str_mod
!> @brief   molecule structure
!! @authors Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module molecules_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_molecule

    integer(iintegers):: num_deg_freedom = 0

    integer           :: num_atoms     = 0
    integer           :: num_bonds     = 0
    integer           :: num_enm_bonds = 0
    integer           :: num_angles    = 0
    integer           :: num_dihedrals = 0
    integer           :: num_impropers = 0
    integer           :: num_cmaps     = 0

    integer           :: num_residues  = 0
    integer           :: num_molecules = 0
    integer           :: num_segments  = 0

    ! shift or not shift
    logical           :: shift_origin   = .false.
    logical           :: special_hydrogen = .false.

    real(wp)          :: total_charge  = 0.0_wp

    ! atom
    integer,          allocatable :: atom_no(:)
    character(4),     allocatable :: segment_name(:)
    integer,          allocatable :: segment_no(:)
    integer,          allocatable :: residue_no(:)
    integer,          allocatable :: residue_c_no(:)
    character(6),     allocatable :: residue_name(:)
    character(4),     allocatable :: atom_name(:)
    character(6),     allocatable :: atom_cls_name(:)
    integer,          allocatable :: atom_cls_no(:)
    real(wp),         allocatable :: charge(:)
    real(wp),         allocatable :: mass(:)
    real(wp),         allocatable :: inv_mass(:)
    integer,          allocatable :: imove(:)
    real(wp),         allocatable :: stokes_radius(:)
    real(wp),         allocatable :: inv_stokes_radius(:)

    character(1),     allocatable :: chain_id(:)
    real(wp),         allocatable :: atom_coord(:,:)
    real(wp),         allocatable :: atom_occupancy(:)
    real(wp),         allocatable :: atom_temp_factor(:)
    real(wp),         allocatable :: atom_velocity(:,:)
    logical,          allocatable :: light_atom_name(:)
    logical,          allocatable :: light_atom_mass(:)

    integer,          allocatable :: molecule_no(:)

    ! bond
    integer,          allocatable :: bond_list(:,:)

    ! enm
    integer,          allocatable :: enm_list(:,:)

    ! angle
    integer,          allocatable :: angl_list(:,:)

    ! dihedral
    integer,          allocatable :: dihe_list(:,:)

    ! improper
    integer,          allocatable :: impr_list(:,:)

    ! cmap
    integer,          allocatable :: cmap_list(:,:)

    ! molecule
    integer,          allocatable :: molecule_atom_no(:)
    real(wp),         allocatable :: molecule_mass(:)
    character(10),    allocatable :: molecule_name(:)

    ! ref coord
    real(wp),         allocatable :: atom_refcoord(:,:)

    ! fit coord
    real(wp),         allocatable :: atom_fitcoord(:,:)

    ! principal component mode
    integer                       :: num_pc_modes
    real(wp),         allocatable :: pc_mode(:)

    ! FEP
    integer                       :: fep_topology       = 1
    integer                       :: num_hbonds_singleA = 0
    integer                       :: num_hbonds_singleB = 0
    integer                       :: num_atoms_fep(5)
    integer                       :: num_bonds_fep(5)
    integer                       :: num_angles_fep(5)
    integer                       :: num_dihedrals_fep(5)
    integer                       :: num_impropers_fep(5)
    integer                       :: num_cmaps_fep(5)
    integer,          allocatable :: bond_list_fep(:,:,:)
    integer,          allocatable :: angl_list_fep(:,:,:)
    integer,          allocatable :: dihe_list_fep(:,:,:)
    integer,          allocatable :: impr_list_fep(:,:,:)
    integer,          allocatable :: cmap_list_fep(:,:,:)
    integer,          allocatable :: id_singleA(:)
    integer,          allocatable :: id_singleB(:)
    integer,          allocatable :: fepgrp(:)
    integer                       :: fepgrp_bond(5,5)
    integer                       :: fepgrp_angl(5,5,5)
    integer                       :: fepgrp_dihe(5,5,5,5)
    integer                       :: fepgrp_cmap(5*5*5*5*5*5*5*5)

  end type s_molecule

  ! parameters for allocatable variables
  integer,      public, parameter :: MoleculeAtom = 1
  integer,      public, parameter :: MoleculeBond = 2
  integer,      public, parameter :: MoleculeAngl = 3
  integer,      public, parameter :: MoleculeDihe = 4
  integer,      public, parameter :: MoleculeImpr = 5
  integer,      public, parameter :: MoleculeCmap = 6
  integer,      public, parameter :: MoleculeMole = 7
  integer,      public, parameter :: MoleculeRefc = 8
  integer,      public, parameter :: MoleculeMode = 9
  integer,      public, parameter :: MoleculeFitc = 10
  integer,      public, parameter :: MoleculeEnm  = 11

  ! subroutines
  public :: init_molecules
  public :: alloc_molecules
  public :: dealloc_molecules
  public :: dealloc_molecules_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_molecules
  !> @brief        initialize molecule information
  !! @authors      YS
  !! @param[out]   molecule : structure of molecule
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_molecules(molecule)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule


    molecule%num_deg_freedom  = 0
    molecule%num_atoms        = 0
    molecule%num_bonds        = 0
    molecule%num_angles       = 0
    molecule%num_dihedrals    = 0
    molecule%num_impropers    = 0
    molecule%num_cmaps        = 0
    molecule%num_residues     = 0
    molecule%num_molecules    = 0
    molecule%num_segments     = 0
    molecule%total_charge     = 0.0_wp
    molecule%special_hydrogen = .false.

    return

  end subroutine init_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_molecules
  !> @brief        allocate molecule information
  !! @authors      YS, CK
  !! @param[inout] molecule : structure of molecule information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_molecules(molecule, variable, var_size)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables

    select case (variable)

    case(MoleculeAtom)

      if (allocated(molecule%atom_no)) then
        if (size(molecule%atom_no) == var_size) return
        deallocate(molecule%atom_no,          &
                   molecule%segment_name,     &
                   molecule%segment_no,       &
                   molecule%residue_no,       &
                   molecule%residue_c_no,     &
                   molecule%residue_name,     &
                   molecule%atom_name,        &
                   molecule%atom_cls_name,    &
                   molecule%atom_cls_no,      &
                   molecule%charge,           &
                   molecule%mass,             &
                   molecule%inv_mass,         &
                   molecule%imove,            &
                   molecule%stokes_radius,    &
                   molecule%inv_stokes_radius,&
                   molecule%chain_id,         &
                   molecule%atom_coord,       &
                   molecule%atom_occupancy,   &
                   molecule%atom_temp_factor, &
                   molecule%atom_velocity,    &
                   molecule%molecule_no,      &
                   molecule%light_atom_name,  &
                   molecule%light_atom_mass,  &
                   stat = dealloc_stat)
      end if

      allocate(molecule%atom_no(var_size),           &
               molecule%segment_name(var_size),      &
               molecule%segment_no(var_size),        &
               molecule%residue_no(var_size),        &
               molecule%residue_c_no(var_size),      &
               molecule%residue_name(var_size),      &
               molecule%atom_name(var_size),         &
               molecule%atom_cls_name(var_size),     &
               molecule%atom_cls_no(var_size),       &
               molecule%charge(var_size),            &
               molecule%mass(var_size),              &
               molecule%inv_mass(var_size),          &
               molecule%imove(var_size),             &
               molecule%stokes_radius(var_size),     &
               molecule%inv_stokes_radius(var_size), &
               molecule%chain_id(var_size),          &
               molecule%atom_coord(3,var_size),      &
               molecule%atom_occupancy(var_size),    &
               molecule%atom_temp_factor(var_size),  &
               molecule%atom_velocity(3,var_size),   &
               molecule%molecule_no(var_size),       &
               molecule%light_atom_mass(var_size),   &
               molecule%light_atom_name(var_size),   &
               stat = alloc_stat)

      molecule%atom_no           (1:var_size) = 0
      molecule%segment_name      (1:var_size) = ''
      molecule%segment_no        (1:var_size) = 0
      molecule%residue_no        (1:var_size) = 0
      molecule%residue_c_no      (1:var_size) = 0
      molecule%residue_name      (1:var_size) = ''
      molecule%atom_name         (1:var_size) = ''
      molecule%atom_cls_name     (1:var_size) = ''
      molecule%atom_cls_no       (1:var_size) = 0
      molecule%charge            (1:var_size) = 0.0_wp
      molecule%mass              (1:var_size) = 0.0_wp
      molecule%inv_mass          (1:var_size) = 0.0_wp
      molecule%imove             (1:var_size) = 0
      molecule%stokes_radius     (1:var_size) = 0.0_wp
      molecule%inv_stokes_radius (1:var_size) = 0.0_wp
      molecule%chain_id          (1:var_size) = ''
      molecule%atom_coord    (1:3,1:var_size) = 0.0_wp
      molecule%atom_occupancy    (1:var_size) = 0.0_wp
      molecule%atom_temp_factor  (1:var_size) = 0.0_wp
      molecule%atom_velocity (1:3,1:var_size) = 0.0_wp
      molecule%molecule_no       (1:var_size) = 0
      molecule%light_atom_name   (1:var_size) = .false.
      molecule%light_atom_mass   (1:var_size) = .false.

    case (MoleculeBond)

      if (allocated(molecule%bond_list)) then
        if (size(molecule%bond_list) == var_size) &
          return
        deallocate(molecule%bond_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%bond_list(2, var_size), &
               stat = alloc_stat)

      molecule%bond_list(1:2,1:var_size) = 0

    case (MoleculeEnm)

      if (allocated(molecule%enm_list)) then
        if (size(molecule%enm_list) == var_size) &
          return
        deallocate(molecule%enm_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%enm_list(2, var_size), &
               stat = alloc_stat)

      molecule%enm_list(1:2,1:var_size) = 0

    case (MoleculeAngl)

      if (allocated(molecule%angl_list)) then
        if (size(molecule%angl_list) == var_size) &
          return
        deallocate(molecule%angl_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%angl_list(3, var_size), &
               stat = alloc_stat)

      molecule%angl_list(1:3,1:var_size) = 0

    case (MoleculeDihe)

      if (allocated(molecule%Dihe_list)) then
        if (size(molecule%dihe_list) == var_size) &
          return
        deallocate(molecule%dihe_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%dihe_list(4, var_size), &
               stat = alloc_stat)

      molecule%dihe_list(1:4,1:var_size) = 0

    case (MoleculeImpr)

      if (allocated(molecule%impr_list)) then
        if (size(molecule%impr_list) == var_size) &
          return
        deallocate(molecule%impr_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%impr_list(4, var_size), &
               stat = alloc_stat)

      molecule%impr_list(1:4,1:var_size) = 0

    case (MoleculeCmap)

      if (allocated(molecule%cmap_list)) then
        if (size(molecule%cmap_list(1,:)) == var_size) &
          return
        deallocate(molecule%cmap_list, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%cmap_list(8, var_size), &
              stat = alloc_stat)

      molecule%cmap_list(1:8,1:var_size) = 0

    case (MoleculeMole)

      if (allocated(molecule%molecule_atom_no)) then
        if (size(molecule%molecule_atom_no) == var_size) &
          return
        deallocate(molecule%molecule_atom_no, &
                   molecule%molecule_mass,    &
                   molecule%molecule_name,    &
                   stat = dealloc_stat)
      end if

      allocate(molecule%molecule_atom_no(var_size), &
               molecule%molecule_mass(var_size),    &
               molecule%molecule_name(var_size),    &
               stat = alloc_stat)

      molecule%molecule_atom_no(1:var_size) = 0
      molecule%molecule_mass   (1:var_size) = 0.0_wp
      molecule%molecule_name   (1:var_size) = ''

    case (MoleculeRefc)

      if (allocated(molecule%atom_refcoord)) then
        if (size(molecule%atom_refcoord(1,:)) == var_size) &
          return
        deallocate(molecule%atom_refcoord, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%atom_refcoord(3,var_size), &
               stat = alloc_stat)

      molecule%atom_refcoord(1:3,1:var_size)  = 0.0_wp

    case (MoleculeFitc)

      if (allocated(molecule%atom_fitcoord)) then
        if (size(molecule%atom_fitcoord(1,:)) == var_size) &
          return
        deallocate(molecule%atom_fitcoord, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%atom_fitcoord(3,var_size), &
               stat = alloc_stat)

      molecule%atom_fitcoord(1:3,1:var_size)  = 0.0_wp

    case (MoleculeMode)

      if (allocated(molecule%pc_mode)) then
        if (size(molecule%pc_mode(:)) == var_size) &
          return
        deallocate(molecule%pc_mode, &
                   stat = dealloc_stat)
      end if

      allocate(molecule%pc_mode(1:var_size), &
               stat = alloc_stat)

      molecule%pc_mode(1:var_size)  = 0.0_wp

    case default

      call error_msg('Alloc_Molecule> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_molecules
  !> @brief        deallocate molecule information
  !! @authors      YS, CK
  !! @param[inout] molecule : structure of molecule information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_molecules(molecule, variable)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case(MoleculeAtom)

      if (allocated(molecule%atom_no)) then
        deallocate (molecule%atom_no,          &
                    molecule%segment_name,     &
                    molecule%segment_no,       &
                    molecule%residue_no,       &
                    molecule%residue_c_no,     &
                    molecule%residue_name,     &
                    molecule%atom_name,        &
                    molecule%atom_cls_name,    &
                    molecule%atom_cls_no,      &
                    molecule%charge,           &
                    molecule%mass,             &
                    molecule%inv_mass,         &
                    molecule%imove,            &
                    molecule%stokes_radius,    &
                    molecule%inv_stokes_radius,&
                    molecule%chain_id,         &
                    molecule%atom_coord,       &
                    molecule%atom_occupancy,   &
                    molecule%atom_temp_factor, &
                    molecule%atom_velocity,    &
                    molecule%molecule_no,      &
                    molecule%light_atom_name,  &
                    molecule%light_atom_mass,  &
                    stat = dealloc_stat)
      end if

    case (MoleculeBond)

      if (allocated(molecule%bond_list)) then
        deallocate (molecule%bond_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeEnm)

      if (allocated(molecule%enm_list)) then
        deallocate (molecule%enm_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeAngl)

      if (allocated(molecule%angl_list)) then
        deallocate (molecule%angl_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeDihe)

      if (allocated(molecule%Dihe_list)) then
        deallocate (molecule%dihe_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeImpr)

      if (allocated(molecule%impr_list)) then
        deallocate (molecule%impr_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeCmap)

      if (allocated(molecule%cmap_list)) then
        deallocate (molecule%cmap_list, &
                    stat = dealloc_stat)
      end if

    case (MoleculeMole)

      if (allocated(molecule%molecule_atom_no)) then
        deallocate (molecule%molecule_atom_no, &
                    molecule%molecule_mass,    &
                    molecule%molecule_name,    &
                    stat = dealloc_stat)
      end if

    case (MoleculeRefc)

      if (allocated(molecule%atom_refcoord)) then
        deallocate (molecule%atom_refcoord, &
                    stat = dealloc_stat)
       end if

    case (MoleculeFitc)

      if (allocated(molecule%atom_fitcoord)) then
        deallocate (molecule%atom_fitcoord, &
                    stat = dealloc_stat)
       end if

    case (MoleculeMode)

      if (allocated(molecule%pc_mode)) then
        deallocate (molecule%pc_mode, &
                    stat = dealloc_stat)
       end if

    case default

      call error_msg('Dealloc_Molecule> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_molecules_all
  !> @brief        deallocate all molecule information
  !! @authors      YS
  !! @param[inout] molecule : structure of molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_molecules_all(molecule)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule


    call dealloc_molecules(molecule, MoleculeAtom)
    call dealloc_molecules(molecule, MoleculeBond)
    call dealloc_molecules(molecule, MoleculeEnm)
    call dealloc_molecules(molecule, MoleculeAngl)
    call dealloc_molecules(molecule, MoleculeDihe)
    call dealloc_molecules(molecule, MoleculeImpr)
    call dealloc_molecules(molecule, MoleculeCmap)
    call dealloc_molecules(molecule, MoleculeMole)
    call dealloc_molecules(molecule, MoleculeRefc)
    call dealloc_molecules(molecule, MoleculeFitc)
    call dealloc_molecules(molecule, MoleculeMode)

    return

  end subroutine dealloc_molecules_all

end module molecules_str_mod
