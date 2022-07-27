!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   molecules_mod
!> @brief   define molecules information
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK), 
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module molecules_mod

  use fileio_top_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use fileio_crd_mod
  use fileio_prmtop_mod
  use fileio_ambcrd_mod
  use fileio_grotop_mod
  use fileio_grocrd_mod
  use fileio_mode_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! local variables
  logical, private :: vervose = .true.   

  ! subroutines
  public  :: define_molecules
  public  :: update_num_deg_freedom
  public  :: export_molecules
  private :: setup_molecule_charmm
  private :: setup_molecule_amber
  private :: setup_molecule_gromacs
  private :: setup_molecule_crd
  private :: setup_molecule_ref
  private :: setup_molecule_fit
  private :: setup_molecule_mode
  private :: setup_molecule_other
  private :: initialize_velocity
  private :: export_molecule_sel
  private :: export_molecule_all
  public  :: check_light_atom_name
  private :: export_molecule_sel_crd
  private :: export_molecule_all_crd
  private :: export_molecule_sel_psf

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_molecules
  !> @brief        a driver subroutine for defining moleucles
  !! @authors      YS, CK, TM, NT, YK
  !! @param[inout] molecule : structure of molecule
  !! @param[in]    pdb      : PDB information (optional)
  !! @param[in]    crd      : CRD information (optional)
  !! @param[in]    top      : TOP information (optional)
  !! @param[in]    par      : PAR information (optional)
  !! @param[in]    gpr      : GPR information (optional)
  !! @param[in]    psf      : PSF information (optional)
  !! @param[in]    ref      : REF information (optional)
  !! @param[in]    fit      : FIT information (optional)
  !! @param[in]    mode     : MODE information (optional)
  !! @param[in]    prmtop   : PRMTOP information (optional)
  !! @param[in]    ambdrd   : AMBCRD information (optional)
  !! @param[in]    ambref   : AMBREF information (optional)
  !! @param[in]    grotop   : GROTOP information (optional)
  !! @param[in]    grocrd   : GROCRD information (optional)
  !! @param[in]    groref   : GROREF information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit,&
                              mode, prmtop, ambcrd, ambref,                    &
                              grotop, grocrd, groref)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_pdb),    optional, intent(in)    :: pdb
    type(s_crd),    optional, intent(in)    :: crd
    type(s_top),    optional, intent(in)    :: top
    type(s_par),    optional, intent(in)    :: par
    type(s_gpr),    optional, intent(in)    :: gpr
    type(s_psf),    optional, intent(in)    :: psf
    type(s_pdb),    optional, intent(in)    :: ref
    type(s_pdb),    optional, intent(in)    :: fit
    type(s_mode),   optional, intent(in)    :: mode
    type(s_prmtop), optional, intent(in)    :: prmtop
    type(s_ambcrd), optional, intent(in)    :: ambcrd
    type(s_ambcrd), optional, intent(in)    :: ambref
    type(s_grotop), optional, intent(in)    :: grotop
    type(s_grocrd), optional, intent(in)    :: grocrd
    type(s_grocrd), optional, intent(in)    :: groref

    ! local variables
    integer         :: i, j
    integer         :: num_atoms_pdb, num_atoms_psf, num_atoms_crd
    integer         :: num_atoms_prmtop, num_atoms_ambcrd
    integer         :: num_atoms_grotop, num_atoms_grocrd
    logical         :: present_pdb, present_crd, present_ref
    logical         :: present_fit
    logical         :: present_mode
    logical         :: present_top, present_par, present_psf, present_gpr
    logical         :: present_prmtop, present_ambcrd, present_ambref
    logical         :: present_grotop, present_grocrd, present_groref


    present_pdb    = .false.
    present_crd    = .false.
    present_psf    = .false.
    present_ref    = .false.
    present_fit    = .false.
    present_top    = .false.
    present_par    = .false.
    present_gpr    = .false.
    present_mode   = .false.
    present_prmtop = .false.
    present_ambcrd = .false.
    present_ambref = .false.
    present_grotop = .false.
    present_grocrd = .false.
    present_groref = .false.

    if (present(pdb)) then
      if (allocated(pdb%atom_no)) present_pdb = .true.
    end if
    if (present(crd)) then
      if (allocated(crd%atom_no)) present_crd = .true.
    end if
    if (present(psf)) then
      if (allocated(psf%atom_no)) present_psf = .true.
    end if
    if (present(ref)) then
      if (allocated(ref%atom_no)) present_ref = .true.
    end if
    if (present(fit)) then
      if (allocated(fit%atom_no)) present_fit = .true.
    end if
    if (present(top)) then
      if (allocated(top%atom_cls_no))   present_top  = .true.
    end if
    if (present(par)) then
      if (allocated(par%bond_atom_cls)) present_par  = .true.
    end if
    if (present(gpr)) then
      if (allocated(gpr%bond_list))     present_gpr  = .true.
    end if
    if (present(mode)) then
      if (allocated(mode%pc_mode))      present_mode = .true.
    end if
    if (present(prmtop)) then
      if (allocated(prmtop%atom_name))  present_prmtop = .true.
    end if
    if (present(ambcrd)) then
      if (allocated(ambcrd%atom_coord)) present_ambcrd = .true.
    end if
    if (present(ambref)) then
      if (allocated(ambref%atom_coord)) present_ambref = .true.
    end if
    if (present(grotop)) then
      if (allocated(grotop%molss))      present_grotop = .true.
    end if
    if (present(grocrd)) then
      if (allocated(grocrd%atom_no))    present_grocrd = .true.
    end if
    if (present(groref)) then
      if (allocated(groref%atom_no))    present_groref = .true.
    end if


    ! Check about the number of atoms in PSF, PDB, CRD
    !
    if (present_pdb) num_atoms_pdb = size(pdb%atom_no)
    if (present_crd) num_atoms_crd = size(crd%atom_no)
    if (present_psf) num_atoms_psf = size(psf%atom_no)

    if (present_psf .and. present_pdb) then
      if (num_atoms_psf /= num_atoms_pdb) then
        write(MsgOut,*) 'Define_Molecules> ATOM in PSF = ', num_atoms_psf
        write(MsgOut,*) 'Define_Molecules> ATOM in PDB = ', num_atoms_pdb
        call error_msg('Define_Molecules> Number of atoms differs.')
      end if
    end if

    if (present_psf .and. present_crd) then
      if (num_atoms_psf /= num_atoms_crd) then
        write(MsgOut,*) 'Define_Molecules> ATOM in PSF = ', num_atoms_psf
        write(MsgOut,*) 'Define_Molecules> ATOM in CRD = ', num_atoms_crd
        call error_msg('Define_Molecules> Number of atoms differs.')
      end if
    end if

    if (present_pdb .and. present_crd) then
      if (num_atoms_pdb /= num_atoms_crd) then
        write(MsgOut,*) 'Define_Molecules> ATOM in PDB = ', num_atoms_pdb
        write(MsgOut,*) 'Define_Molecules> ATOM in CRD = ', num_atoms_crd
        call error_msg('Define_Molecules> Number of atoms differs.')
      end if
    end if


    ! Check about the number of atoms in PRMTOP, AMBCRD
    !
    if (present_prmtop) num_atoms_prmtop = size(prmtop%atom_name)
    if (present_ambcrd) num_atoms_ambcrd = size(ambcrd%atom_coord(1,:))

    if (present_prmtop .and. present_ambcrd) then
      if (num_atoms_prmtop /= num_atoms_ambcrd) then
        write(MsgOut,*) 'Define_Molecules> ATOM in PRMTOP = ', num_atoms_prmtop
        write(MsgOut,*) 'Define_Molecules> ATOM in AMBCRD = ', num_atoms_ambcrd
        call error_msg('Define_Molecules> Number of atoms differs.')
      end if
    end if


    ! Check about the number of atoms in GROTOP, GROCRD
    !
    if (present_grotop) then
      num_atoms_grotop = 0
      do i = 1, grotop%num_molss
        do j = 1, grotop%molss(i)%count
          num_atoms_grotop = &
               num_atoms_grotop + grotop%molss(i)%moltype%mol%num_atoms
        end do
      end do
    end if
    if (present_grocrd) num_atoms_grocrd = size(grocrd%atom_no)

    if (present_grotop .and. present_grocrd) then
      if (num_atoms_grotop /= num_atoms_grocrd) then
        write(MsgOut,*) 'Define_Molecules> ATOM in GROTOP = ', num_atoms_grotop
        write(MsgOut,*) 'Define_Molecules> ATOM in GROCRD = ', num_atoms_grocrd
        call error_msg('Define_Molecules> Number of atoms differs.')
      end if
    end if


    ! Setup molecule
    !
    if (present_top .and. present_par .and. present_psf) then  !< CHARMM

      if (.not. present_pdb .and. .not. present_crd) &
        call error_msg('Define_Molecules> coordinates are not present.')

      call setup_molecule_charmm(molecule, psf, pdb, crd, top, par)
      call setup_molecule_other (molecule)

    else if (present_top .and. present_gpr .and. present_psf) then  !< Go model

      if (.not. present_pdb .and. .not. present_crd) &
        call error_msg('Define_Molecules> coordinates are not present.')

      call setup_molecule_charmm(molecule, psf, pdb, crd)
      call setup_molecule_other (molecule)

    else if (present_prmtop) then !< AMBER

      if (.not. present_ambcrd .and. .not. present_pdb) &
        call error_msg('Define_Molecules> coordinates are not present.')

      call setup_molecule_amber(molecule, prmtop, ambcrd, pdb)
      call setup_molecule_other(molecule)

    else if (present_grotop) then !< GROMACS

      if (.not. present_grocrd .and. .not. present_pdb) &
        call error_msg('Define_Molecules> coordinates are not present.')
           
      call setup_molecule_gromacs(molecule, grotop, grocrd, pdb)
      call setup_molecule_other  (molecule)

    else if (present_psf) then

      if (.not. present_pdb .and. .not. present_crd) then
        call error_msg('Define_Molecules> coordinates are not present.')
      end if

      call setup_molecule_charmm(molecule, psf, pdb, crd)
      call setup_molecule_other (molecule)

    else

      if (.not. present_pdb    .and. &
          .not. present_crd    .and. &
          .not. present_grocrd) then
        call error_msg('Define_Molecules> coordinates are not present.')
      end if

      call setup_molecule_crd(molecule, pdb, crd, grocrd)

    end if

    !
    ! Check light atom name and mass
    !
!   call check_light_atom_name(molecule)


    ! Initialize velocities
    !
    call initialize_velocity(molecule)


    ! Setup reference coordinates in molecule
    !
    if (present_ref .or. present_ambref .or. present_groref) then

      call setup_molecule_ref(molecule, ref, ambref, groref)

    end if


    ! Setup fit coordinates in molecule
    !
    if (present_fit) then

      call setup_molecule_fit(molecule, fit)

    end if


    ! Setup principal component mode in molecule
    !
    if (present_mode) then

      call setup_molecule_mode(molecule, mode)

    end if


    ! write summary of molecule information
    !
    if (main_rank .and. vervose) then

      write(MsgOut,'(A)') 'Define_Molecule> Summary of molecules'

      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_atoms       = ', molecule%num_atoms,     &
           '  num_bonds       = ', molecule%num_bonds
      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_angles      = ', molecule%num_angles,    &
           '  num_dihedrals   = ', molecule%num_dihedrals
      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_impropers   = ', molecule%num_impropers, &
           '  num_cmap_terms  = ', molecule%num_cmaps
      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_residues    = ', molecule%num_residues,  &
           '  num_molecules   = ', molecule%num_molecules
      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_segments    = ', molecule%num_segments,  &
           '  num_deg_freedom = ', molecule%num_deg_freedom
      write(MsgOut,'(A20,F10.3)')                          &
           '  total_charge    = ', molecule%total_charge
      write(MsgOut,'(A)') ' '
      vervose = .false.
    end if

    return

  end subroutine define_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_num_deg_freedom
  !> @brief        update number of degree of freedom
  !! @authors      TM
  !! @param[in]    message         : update message
  !! @param[in]    add_ndegf       : added number of degree of freedom
  !! @param[inout] num_deg_freedom : number of degree of freedom
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_num_deg_freedom(message, add_ndegf, num_deg_freedom)

    ! formal arguments
    character(*),            intent(in)    :: message
    integer,                 intent(in)    :: add_ndegf
    integer(iintegers),      intent(inout) :: num_deg_freedom


    num_deg_freedom = num_deg_freedom + add_ndegf

    if (main_rank) then
      write(MsgOut,'(A)') &
           'Update_Num_Deg_Freedom> Number of degrees of freedom was updated'
      write(MsgOut,'(A20,I10,3A)') &
           '  num_deg_freedom = ', num_deg_freedom, &
           ' (',trim(message), ')'
      write(MsgOut,'(A)') &
           ' '
    end if

    return

  end subroutine update_num_deg_freedom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecules
  !> @brief        a driver subroutine for export moleucles
  !! @authors      NT, KYMD
  !! @param[in]    molecule : molecule information
  !! @param[in]    selatoms : selection information  (optional)
  !! @param[out]   pdb      : PDB information        (optional)
  !! @param[out]   crd      : CRD information        (optional)
  !! @param[out]   top      : TOP information        (optional)
  !! @param[out]   par      : par information        (optional)
  !! @param[out]   gpr      : GPR information        (optional)
  !! @param[out]   psf      : PSF information        (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecules(molecule, selatoms, pdb, crd, top, par, gpr, psf)

    ! formal arguments
    type(s_molecule),           intent(in)    :: molecule
    type(s_selatoms), optional, intent(in)    :: selatoms
    type(s_pdb),      optional, intent(inout) :: pdb
    type(s_crd),      optional, intent(inout) :: crd
    type(s_top),      optional, intent(inout) :: top
    type(s_par),      optional, intent(inout) :: par
    type(s_gpr),      optional, intent(inout) :: gpr
    type(s_psf),      optional, intent(inout) :: psf


    if (present(pdb)) then

      if (present(selatoms)) then
        call export_molecule_sel(molecule, selatoms, pdb)
      else
        call export_molecule_all(molecule, pdb)
      end if

    end if

    if (present(crd)) then

      if (present(selatoms)) then
        call export_molecule_sel_crd(molecule, selatoms, crd)
      else
        call export_molecule_all_crd(molecule, crd)
      end if

    end if

    if (present(top)) then
      ! Not supported yet.
    end if

    if (present(par)) then
      ! Not supported yet.
    end if

    if (present(gpr)) then
      ! Not supported yet.
    end if

    if (present(psf)) then

      if (present(selatoms)) then
        call export_molecule_sel_psf(molecule, selatoms, psf)
      end if

    end if

    return

  end subroutine export_molecules

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_charmm
  !> @brief        setup molecules by CHARMM input informations
  !! @authors      YS, NT
  !! @param[inout] molecule : molecule information
  !! @param[in]    psf      : PSF information
  !! @param[in]    pdb      : PDB information (optional)
  !! @param[in]    crd      : CRD information (optional)
  !! @param[in]    top      : TOP information (optional)
  !! @param[in]    par      : PAR information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_charmm(molecule, psf, pdb, crd, top, par)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_psf),             intent(in)    :: psf
    type(s_pdb),   optional, intent(in)    :: pdb
    type(s_crd),   optional, intent(in)    :: crd
    type(s_top),   optional, intent(in)    :: top
    type(s_par),   optional, intent(in)    :: par

    ! local variables
    integer                  :: natom
    integer                  :: nbond
    integer                  :: nenm 
    integer                  :: nangl
    integer                  :: ndihe
    integer                  :: nimpr
    integer                  :: ncmap
    integer                  :: i, j, k, icls
    logical                  :: present_pdb
    character(6)             :: acls
    character(6)             :: cres


    ! initialize
    !
    call init_molecules(molecule)

    ! allocation
    !
    natom = psf%num_atoms
    nbond = psf%num_bonds
    nenm  = psf%num_enm_bonds
    nangl = psf%num_angles
    ndihe = psf%num_dihedrals
    nimpr = psf%num_impropers
    ncmap = psf%num_cross_terms

    call alloc_molecules(molecule, MoleculeAtom, natom)
    call alloc_molecules(molecule, MoleculeBond, nbond)
    call alloc_molecules(molecule, MoleculeEnm , nenm)
    call alloc_molecules(molecule, MoleculeAngl, nangl)
    call alloc_molecules(molecule, MoleculeDihe, ndihe)
    call alloc_molecules(molecule, MoleculeImpr, nimpr)
    call alloc_molecules(molecule, MoleculeCmap, ncmap)

    ! copy information from PSF
    !
    molecule%atom_no(1:natom)          = psf%atom_no(1:natom)
    molecule%segment_name(1:natom)     = psf%segment_name(1:natom)
    molecule%residue_no(1:natom)       = psf%residue_no(1:natom)
    molecule%residue_name(1:natom)     = psf%residue_name(1:natom)
    molecule%molecule_no(1:natom)      = psf%molecule_no(1:natom)
    molecule%atom_name(1:natom)        = psf%atom_name(1:natom)
    molecule%atom_cls_name(1:natom)    = psf%atom_cls_name(1:natom)
    molecule%charge(1:natom)           = psf%charge(1:natom)
    molecule%mass(1:natom)             = psf%mass(1:natom)
    do i = 1, natom
      if (psf%mass(i) > EPS) then
        molecule%inv_mass(i)           = 1.0_wp/psf%mass(i)
      else
        molecule%inv_mass(i)           = 0.0_wp
      end if
    end do
    molecule%imove(1:natom)            = psf%imove(1:natom)

    if (nbond >= 1) &
    molecule%bond_list(1:2,1:nbond)    = psf%bond_list(1:2,1:nbond)
    if (nenm >= 1) &
    molecule%enm_list(1:2,1:nenm)      = psf%enm_list (1:2,1:nenm)
    if (nangl >= 1) &
    molecule%angl_list(1:3,1:nangl)    = psf%angl_list(1:3,1:nangl)
    if (ndihe >= 1) &
    molecule%dihe_list(1:4,1:ndihe)    = psf%dihe_list(1:4,1:ndihe)
    if (nimpr >= 1) &
    molecule%impr_list(1:4,1:nimpr)    = psf%impr_list(1:4,1:nimpr)
    if (ncmap >= 1) &
    molecule%cmap_list(1:8,1:ncmap)    = psf%cmap_list(1:8,1:ncmap)

    present_pdb = .false.
    if (present(pdb)) then
      if (allocated(pdb%atom_no)) present_pdb = .true.
    end if

    if (present_pdb) then

      ! copy information from PDB
      !
      molecule%chain_id(1:natom)         = pdb%chain_id(1:natom)
      molecule%atom_coord(1:3,1:natom)   = pdb%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = pdb%atom_occupancy(1:natom)
      molecule%atom_temp_factor(1:natom) = pdb%atom_temp_factor(1:natom)
       
    else

      ! copy information from CRD
      !
      molecule%chain_id(1:natom)         = 'A'
      molecule%atom_coord(1:3,1:natom)   = crd%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = 0.0_wp
      molecule%atom_temp_factor(1:natom) = 0.0_wp
       
    end if
 
    ! number of variables
    !
    molecule%num_atoms       = natom
    molecule%num_bonds       = nbond
    molecule%num_enm_bonds   = nenm 
    molecule%num_angles      = nangl
    molecule%num_dihedrals   = ndihe
    molecule%num_impropers   = nimpr
    molecule%num_cmaps       = ncmap
    molecule%num_deg_freedom = 3*natom

    ! redefine atom class name
    !
    if (present(top) .and. psf%type == PsfTypeCHARMM) then

      do i = 1, molecule%num_atoms

        cres = molecule%residue_name(i)
        icls = psf%atom_cls_no(i)

        molecule%atom_cls_name(i) = ''

        do j = 1, top%num_res_type
          if (cres == top%res_name(j)) then
            do k = top%res_atom_cls_bgn(j), top%res_atom_cls_end(j)
              if (icls == top%atom_cls_no(k)) then
                molecule%atom_cls_name(i) = top%atom_cls_name(k)
                goto 10
              end if
            end do
          end if
        end do

10      if (molecule%atom_cls_name(i) == '') &
          call error_msg('Setup_Molecule_Charmm> Unknown atom class no')

      end do
    end if

    ! define atom class number
    !
    if (present(par)) then

      do i = 1, molecule%num_atoms
        acls = molecule%atom_cls_name(i)
        do j = 1, par%num_atom_cls
          if (acls == par%nonb_atom_cls(j)) then
            molecule%atom_cls_no(i) = j
          end if
        end do
      end do

    end if

    return

  end subroutine setup_molecule_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_amber
  !> @brief        setup molecules by AMBER input informations
  !! @authors      NT
  !! @param[inout] molecule : molecule information
  !! @param[in]    prmtop   : PRMTOP information
  !! @param[in]    ambcrd   : AMBCRD information
  !! @param[in]    pdb      : PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_amber(molecule, prmtop, ambcrd, pdb)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_prmtop),           intent(in)    :: prmtop
    type(s_ambcrd), optional, intent(in)    :: ambcrd
    type(s_pdb),    optional, intent(in)    :: pdb

    ! local variables
    real(wp)                  :: center(3)
    integer                   :: i, j, icnt, idih, iimp, natom
    logical                   :: present_pdb


    ! initialize
    !
    call init_molecules(molecule)


    ! atoms
    !
    molecule%num_atoms = prmtop%num_atoms
    molecule%num_deg_freedom = 3 * molecule%num_atoms

    call alloc_molecules(molecule, MoleculeAtom, molecule%num_atoms)

    natom = molecule%num_atoms
    
    ! atom_no
    do i = 1, natom
      molecule%atom_no(i) = i
    end do

    ! atom_name
    molecule%atom_name    (1:natom) = prmtop%atom_name    (1:natom)
    molecule%atom_cls_name(1:natom) = prmtop%atom_cls_name(1:natom)
    molecule%atom_cls_no  (1:natom) = prmtop%atom_cls_no  (1:natom)
    molecule%charge       (1:natom) = prmtop%charge       (1:natom)/18.2223_wp
    molecule%mass         (1:natom) = prmtop%mass         (1:natom)
    molecule%inv_mass     (1:natom) = 1.0_wp/prmtop%mass  (1:natom)

    ! residue_no, residue_name
    if (prmtop%lresidue_label .and. prmtop%lresidue_pointer) then
      do i = 1, prmtop%num_residues-1
        do j = prmtop%res_point(i), prmtop%res_point(i+1)-1
          molecule%residue_no  (j) = i
          molecule%residue_name(j) = prmtop%res_label(i)
        end do
      end do
      do j = prmtop%res_point(i), natom
        molecule%residue_no  (j) = i
        molecule%residue_name(j) = prmtop%res_label(i)
      end do
    end if

    present_pdb = .false.
    if (present(pdb)) then
      if (allocated(pdb%atom_no)) present_pdb = .true.
    end if

    if (present_pdb) then
      molecule%chain_id(1:natom)         = pdb%chain_id(1:natom)
      molecule%atom_coord(1:3,1:natom)   = pdb%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = pdb%atom_occupancy(1:natom)
      molecule%atom_temp_factor(1:natom) = pdb%atom_temp_factor(1:natom)
    else
      molecule%chain_id(1:natom)         = 'A'
      molecule%atom_coord(1:3,1:natom)   = ambcrd%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = 0.0_wp
      molecule%atom_temp_factor(1:natom) = 0.0_wp
    end if

    ! bonds
    !
    molecule%num_bonds = prmtop%num_bondh + prmtop%num_mbonda

    call alloc_molecules(molecule, MoleculeBond, molecule%num_bonds)

    icnt = 0
    do i = 1, prmtop%num_bondh
      icnt = icnt + 1
      molecule%bond_list(1,icnt) = prmtop%bond_inc_hy(1,i) / 3 + 1
      molecule%bond_list(2,icnt) = prmtop%bond_inc_hy(2,i) / 3 + 1
    end do

    do i = 1, prmtop%num_mbonda
      icnt = icnt + 1
      molecule%bond_list(1,icnt) = prmtop%bond_wo_hy(1,i) / 3 + 1
      molecule%bond_list(2,icnt) = prmtop%bond_wo_hy(2,i) / 3 + 1
    end do

    ! angles
    !
    molecule%num_angles = prmtop%num_anglh + prmtop%num_mangla

    call alloc_molecules(molecule, MoleculeAngl, molecule%num_angles)

    icnt = 0
    do i = 1, prmtop%num_anglh
      icnt = icnt + 1
      molecule%angl_list(1,icnt) = prmtop%angl_inc_hy(1,i) / 3 + 1
      molecule%angl_list(2,icnt) = prmtop%angl_inc_hy(2,i) / 3 + 1
      molecule%angl_list(3,icnt) = prmtop%angl_inc_hy(3,i) / 3 + 1
    end do

    do i = 1, prmtop%num_mangla
      icnt = icnt + 1
      molecule%angl_list(1,icnt) = prmtop%angl_wo_hy(1,i) / 3 + 1
      molecule%angl_list(2,icnt) = prmtop%angl_wo_hy(2,i) / 3 + 1
      molecule%angl_list(3,icnt) = prmtop%angl_wo_hy(3,i) / 3 + 1
    end do


    ! dihedrals/impropers
    !
    molecule%num_dihedrals = 0
    molecule%num_impropers = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) then
        molecule%num_impropers = molecule%num_impropers + 1
      else
        molecule%num_dihedrals = molecule%num_dihedrals + 1
      end if
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) then
        molecule%num_impropers = molecule%num_impropers + 1
      else
        molecule%num_dihedrals = molecule%num_dihedrals + 1
      end if
    end do

    call alloc_molecules(molecule, MoleculeDihe, molecule%num_dihedrals)
    call alloc_molecules(molecule, MoleculeImpr, molecule%num_impropers)

    idih = 0
    iimp = 0

    ! check dihe_inc_hy array
    do i = 1, prmtop%num_diheh
      
      ! set improper 
      if (prmtop%dihe_inc_hy(4,i) < 0) then

        iimp = iimp + 1
        molecule%impr_list(1,iimp) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        molecule%impr_list(2,iimp) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        molecule%impr_list(3,iimp) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        molecule%impr_list(4,iimp) = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

      ! set dihedral 
      else

        idih = idih + 1
        molecule%dihe_list(1,idih) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        molecule%dihe_list(2,idih) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        molecule%dihe_list(3,idih) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        molecule%dihe_list(4,idih) =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

      end if

    end do

    ! check dihe_wo_hy array
    !
    do i = 1, prmtop%num_mdihea

      ! set improper 
      if (prmtop%dihe_wo_hy(4,i) < 0) then

        iimp = iimp + 1
        molecule%impr_list(1,iimp) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        molecule%impr_list(2,iimp) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        molecule%impr_list(3,iimp) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        molecule%impr_list(4,iimp) = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

      ! set dihedral 
      else

        idih = idih + 1
        molecule%dihe_list(1,idih) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        molecule%dihe_list(2,idih) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        molecule%dihe_list(3,idih) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        molecule%dihe_list(4,idih) =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

      end if

    end do

    ! cmap 
    !
    molecule%num_cmaps = prmtop%num_cmap 

    call alloc_molecules(molecule, MoleculeCmap, molecule%num_cmaps)

    do i = 1, prmtop%num_cmap

      molecule%cmap_list(1:4,i) = prmtop%cmap_list(1:4,i)
      molecule%cmap_list(5:8,i) = prmtop%cmap_list(2:5,i)

    end do

    ! assigen center of the system and shift coordinates according to center
    !
!   center(1:3) = 0.0_wp
!   do i = 1, natom
!     center(1:3) = center(1:3) + molecule%atom_coord(1:3,i)
!   end do
!   center(1:3) = center(1:3) / real(natom,wp)
!   do i = 1, natom
!     molecule%atom_coord(1:3,i) = molecule%atom_coord(1:3,i) - center(1:3)
!   end do

    return

  end subroutine setup_molecule_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_gromacs
  !> @brief        setup molecules by GROMACS input informations
  !! @authors      NT
  !! @param[inout] molecule : molecule information
  !! @param[in]    grotop   : GROTOP information
  !! @param[in]    grocrd   : GROCRD information
  !! @param[in]    pdb      : PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_gromacs(molecule, grotop, grocrd, pdb)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_grotop),           intent(in)    :: grotop
    type(s_grocrd), optional, intent(in)    :: grocrd
    type(s_pdb),    optional, intent(in)    :: pdb

    ! local variables
    integer                   :: i, j, k, m, ioffset
    integer                   :: natom, nbond, nangl, ndihe, nimpr, nmol
    character(20)             :: name
    character(6)              :: acls
    logical                   :: present_pdb
!    real(wp)                  :: center(1:3)

    type(s_grotop_mol), pointer :: gromol


!TODO (cons_list)       <= s_grotop_mol%s_constr(:)%atom_idx1,atom_idx2
!TODO (cons_leng)       <= s_grotop_mol%s_constr(:)%b0
!TODO (num_constraints) <= s_grotop_mol%num_constrs

    natom = 0
    nbond = 0
    nangl = 0
    ndihe = 0
    nimpr = 0

    call init_molecules(molecule)

    ! setup from GROTOP
    !
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmol   =  grotop%molss(i)%count

      do j = 1, nmol

        natom = natom + gromol%num_atoms

        if (gromol%settles%func == 0) then
          nbond = nbond + gromol%num_bonds
        else
          nbond = nbond + 3
        end if

        if (gromol%settles%func == 0) then
          nangl = nangl + gromol%num_angls
        else
          nangl = nangl + 1
        end if

        do k = 1, gromol%num_dihes
          select case (gromol%dihes(k)%func)
         !shinobu-edited
          case (1,3,4,9,21,22,31,32,33,41,43,52)
            ndihe = ndihe + 1
          case (2)
            nimpr = nimpr + 1
          case default
            call error_msg('Setup_Molecule_Gromacs> Unknown dihedral func type.')
          end select
        end do

        molecule%num_residues = molecule%num_residues + &
                                gromol%atoms(gromol%num_atoms)%residue_no

      end do

    end do

    molecule%num_deg_freedom = 3 * natom

    call alloc_molecules(molecule, MoleculeAtom, natom)
    call alloc_molecules(molecule, MoleculeBond, nbond)
    call alloc_molecules(molecule, MoleculeAngl, nangl)
    call alloc_molecules(molecule, MoleculeDihe, ndihe)
    call alloc_molecules(molecule, MoleculeImpr, nimpr)

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmol   =  grotop%molss(i)%count
      name   =  grotop%molss(i)%name

      do j = 1, nmol

        ioffset = molecule%num_atoms

        ! atom
        do k = 1, gromol%num_atoms
          m                             = molecule%num_atoms + 1
          molecule%atom_no(m)           = m
          molecule%segment_name(m)      = name
          molecule%residue_no(m)        = gromol%atoms(k)%residue_no
          molecule%residue_name(m)      = gromol%atoms(k)%residue_name
          molecule%atom_name(m)         = gromol%atoms(k)%atom_name
          molecule%atom_cls_name(m)     = gromol%atoms(k)%atom_type
          molecule%charge(m)            = gromol%atoms(k)%charge
          molecule%mass(m)              = gromol%atoms(k)%mass
          if (molecule%mass(m) > EPS) then 
            molecule%inv_mass(m)        = 1.0_wp / molecule%mass(m)
          else
            molecule%inv_mass(m)        = 0.0_wp
          end if
          molecule%stokes_radius(m)     = gromol%atoms(k)%stokes_r * 10.0_wp
          !shinobu-edited
          if (molecule%stokes_radius(m) > EPS) then
            molecule%inv_stokes_radius(m) = 0.1_wp / molecule%stokes_radius(m)
          else
            molecule%inv_stokes_radius(m) = 0.0_wp
          end if
          molecule%num_atoms            = m
        end do

        ! bond
        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds
            m                        = molecule%num_bonds + 1
            molecule%bond_list(1, m) = gromol%bonds(k)%atom_idx1 + ioffset
            molecule%bond_list(2, m) = gromol%bonds(k)%atom_idx2 + ioffset
            molecule%num_bonds       = m
          end do

        else

          ! settle O-H (1)
          m                        = molecule%num_bonds + 1
          molecule%bond_list(1, m) = 1 + ioffset
          molecule%bond_list(2, m) = 2 + ioffset
          molecule%num_bonds       = m

          ! settle O-H (2)
          m                        = molecule%num_bonds + 1
          molecule%bond_list(1, m) = 1 + ioffset
          molecule%bond_list(2, m) = 3 + ioffset
          molecule%num_bonds       = m

          ! settle H-H
          m                        = molecule%num_bonds + 1
          molecule%bond_list(1, m) = 2 + ioffset
          molecule%bond_list(2, m) = 3 + ioffset
          molecule%num_bonds       = m

        end if
        
        ! angle
        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls
            m                        = molecule%num_angles + 1
            molecule%angl_list(1, m) = gromol%angls(k)%atom_idx1 + ioffset
            molecule%angl_list(2, m) = gromol%angls(k)%atom_idx2 + ioffset
            molecule%angl_list(3, m) = gromol%angls(k)%atom_idx3 + ioffset
            molecule%num_angles      = m
          end do

        else

          ! settle H-O-H
          m                        = molecule%num_angles + 1
          molecule%angl_list(1, m) = 2 + ioffset
          molecule%angl_list(2, m) = 1 + ioffset
          molecule%angl_list(3, m) = 3 + ioffset
          molecule%num_angles      = m

        end if

        ! dihedral/improper dihedral
        do k = 1, gromol%num_dihes
          select case (gromol%dihes(k)%func)
          !shinobu-edited
          case (1,3,4,9,21,22,31,32,33,41,43,52)
            m                        = molecule%num_dihedrals + 1
            molecule%dihe_list(1, m) = gromol%dihes(k)%atom_idx1 + ioffset
            molecule%dihe_list(2, m) = gromol%dihes(k)%atom_idx2 + ioffset
            molecule%dihe_list(3, m) = gromol%dihes(k)%atom_idx3 + ioffset
            molecule%dihe_list(4, m) = gromol%dihes(k)%atom_idx4 + ioffset
            molecule%num_dihedrals   = m
          case (2)
            m                        = molecule%num_impropers + 1
            molecule%impr_list(1, m) = gromol%dihes(k)%atom_idx1 + ioffset
            molecule%impr_list(2, m) = gromol%dihes(k)%atom_idx2 + ioffset
            molecule%impr_list(3, m) = gromol%dihes(k)%atom_idx3 + ioffset
            molecule%impr_list(4, m) = gromol%dihes(k)%atom_idx4 + ioffset
            molecule%num_impropers   = m
          end select
        end do

      end do

    end do


    natom = molecule%num_atoms

    ! setup atom class number
    !
    do i = 1, natom
      acls = molecule%atom_cls_name(i)
      do j = 1, grotop%num_atomtypes
        if (acls == grotop%atomtypes(j)%type_name) then
          molecule%atom_cls_no(i) = j
          exit
        end if
      end do
    end do

    ! setup from GRCRD / PDB
    !
    present_pdb = .false.
    if (present(pdb)) then
      if (allocated(pdb%atom_no)) present_pdb = .true.
    end if

    if (present_pdb) then
      molecule%chain_id(1:natom)         = pdb%chain_id(1:natom)
      molecule%atom_coord(1:3,1:natom)   = pdb%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = pdb%atom_occupancy(1:natom)
      molecule%atom_temp_factor(1:natom) = pdb%atom_temp_factor(1:natom)
    else
      molecule%chain_id(1:natom)         = 'A'
      molecule%atom_coord(1:3,1:natom)   =grocrd%atom_coord(1:3,1:natom)*10.0_wp
      molecule%atom_occupancy(1:natom)   = 0.0_wp
      molecule%atom_temp_factor(1:natom) = 0.0_wp
    end if

    ! assigen center of the system and shift coordinates according to center
    !
!    center(1:3) = 0.0_wp
!    do i = 1, natom
!      center(1:3) = center(1:3) + molecule%atom_coord(1:3,i)
!    end do
!    center(1:3) = center(1:3) / real(natom,wp)
!    do i = 1, natom
!      molecule%atom_coord(1:3,i) = molecule%atom_coord(1:3,i) - center(1:3)
!    end do

    return

  end subroutine setup_molecule_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_crd
  !> @brief        setup molecules by PDB/CRD information
  !! @authors      YS, NT
  !! @param[inout] molecule : molecule information
  !! @param[in]    pdb      : PDB information (optional)
  !! @param[in]    crd      : CRD information (optional)
  !! @param[in]    grocrd   : GROCRD information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_crd(molecule, pdb, crd, grocrd)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_pdb),   optional, intent(in)    :: pdb
    type(s_crd),   optional, intent(in)    :: crd
    type(s_grocrd),optional, intent(in)    :: grocrd

    ! local variables
    integer                  :: natom
    logical                  :: present_pdb
    logical                  :: present_crd


    ! initialize
    !
    call init_molecules(molecule)

    present_pdb = .false.
    if (present(pdb)) then
      if (allocated(pdb%atom_no)) present_pdb = .true.
    end if

    present_crd = .false.
    if (present(crd)) then
      if (allocated(crd%atom_no)) present_crd = .true.
    end if

    if (present_pdb) then

      ! copy information from PDB
      !

      natom = pdb%num_atoms
      call alloc_molecules(molecule, MoleculeAtom, natom)

      molecule%atom_no(1:natom)          = pdb%atom_no(1:natom)
      molecule%segment_name(1:natom)     = pdb%segment_name(1:natom)
      molecule%segment_no(1:natom)       = 0
      molecule%residue_no(1:natom)       = pdb%residue_no(1:natom)
      molecule%residue_name(1:natom)     = pdb%residue_name(1:natom)
      molecule%atom_name(1:natom)        = pdb%atom_name(1:natom)
      molecule%atom_cls_name(1:natom)    = ''
      molecule%atom_cls_no(1:natom)      = 0
      molecule%charge(1:natom)           = 0.0_wp
      molecule%mass(1:natom)             = 0.0_wp
      molecule%inv_mass(1:natom)         = 0.0_wp
      molecule%imove(1:natom)            = 0
      molecule%chain_id(1:natom)         = pdb%chain_id(1:natom)
      molecule%atom_coord(1:3,1:natom)   = pdb%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = pdb%atom_occupancy(1:natom)
      molecule%atom_temp_factor(1:natom) = pdb%atom_temp_factor(1:natom)

    else if (present_crd) then

      ! copy information from CRD
      !

      natom = crd%num_atoms
      call alloc_molecules(molecule, MoleculeAtom, natom)

      molecule%atom_no(1:natom)          = crd%atom_no(1:natom)
      molecule%segment_name(1:natom)     = crd%segment_name(1:natom)
      molecule%segment_no(1:natom)       = 0
      molecule%residue_no(1:natom)       = crd%residue_no(1:natom)
      molecule%residue_name(1:natom)     = crd%residue_name(1:natom)
      molecule%atom_name(1:natom)        = crd%atom_name(1:natom)
      molecule%atom_cls_name(1:natom)    = ''
      molecule%atom_cls_no(1:natom)      = 0
      molecule%charge(1:natom)           = 0.0_wp
      molecule%mass(1:natom)             = 0.0_wp
      molecule%inv_mass(1:natom)         = 0.0_wp
      molecule%imove(1:natom)            = 0
      molecule%chain_id(1:natom)         = 'A'
      molecule%atom_coord(1:3,1:natom)   = crd%atom_coord(1:3,1:natom)
      molecule%atom_occupancy(1:natom)   = 0.0_wp
      molecule%atom_temp_factor(1:natom) = 0.0_wp

    else

      ! copy information from GROCRD
      !

      natom = grocrd%num_atoms
      call alloc_molecules(molecule, MoleculeAtom, natom)

      molecule%atom_no(1:natom)          = grocrd%atom_no(1:natom)
      molecule%segment_name(1:natom)     = ''
      molecule%segment_no(1:natom)       = 0
      molecule%residue_no(1:natom)       = grocrd%residue_no(1:natom)
      molecule%residue_name(1:natom)     = grocrd%residue_name(1:natom)
      molecule%atom_name(1:natom)        = grocrd%atom_name(1:natom)
      molecule%atom_cls_name(1:natom)    = ''
      molecule%atom_cls_no(1:natom)      = 0
      molecule%charge(1:natom)           = 0.0_wp
      molecule%mass(1:natom)             = 0.0_wp
      molecule%inv_mass(1:natom)         = 0.0_wp
      molecule%imove(1:natom)            = 0
      molecule%chain_id(1:natom)         = 'A'
      molecule%atom_coord(1:3,1:natom)   = grocrd%atom_coord(1:3,1:natom)*10.0_wp
      molecule%atom_occupancy(1:natom)   = 0.0_wp
      molecule%atom_temp_factor(1:natom) = 0.0_wp

    end if


    ! number of variables
    !
    molecule%num_atoms       = natom
    molecule%num_deg_freedom = 3*natom

    return

  end subroutine setup_molecule_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_ref
  !> @brief        setup molecules by Ref information
  !! @authors      CK
  !! @param[inout] molecule : molecule information
  !! @param[in]    ref      : ref information
  !! @param[in]    ambref   : AMBCRD ref information
  !! @param[in]    groref   : GROCRD ref information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_ref(molecule, ref, ambref, groref)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_pdb),    optional, intent(in)    :: ref
    type(s_ambcrd), optional, intent(in)    :: ambref
    type(s_grocrd), optional, intent(in)    :: groref

    ! local variables
    integer                  :: natom


    natom = molecule%num_atoms

    call alloc_molecules(molecule, MoleculeRefc, natom)

    if (present(ref)) then
      if (allocated(ref%atom_coord)) then
        if (natom /= size(ref%atom_coord(1,:))) then
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in Molecule = ', natom
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in PDB REF  = ', &
                                                       size(ref%atom_coord(1,:))
          call error_msg('Setup_Molecule_Ref> Number of atoms differs.')
        end if
        molecule%atom_refcoord(1:3,1:natom) = ref%atom_coord(1:3,1:natom)
      end if
    end if

    if (present(ambref)) then
      if (allocated(ambref%atom_coord)) then
        if (natom /= size(ambref%atom_coord(1,:))) then
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in Molecule = ', natom
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in AMBREF   = ', &
                                                    size(ambref%atom_coord(1,:))
          call error_msg('Setup_Molecule_Ref> Number of atoms differs.')
        end if
        molecule%atom_refcoord(1:3,1:natom) = ambref%atom_coord(1:3,1:natom)
      end if
    end if
        
    if (present(groref)) then
      if (allocated(groref%atom_coord)) then
        if (natom /= size(groref%atom_coord(1,:))) then
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in Molecule = ', natom
          write(MsgOut,*) 'Setup_Molecule_Ref> ATOM in GROREF   = ', &
                                                    size(groref%atom_coord(1,:))
          call error_msg('Setup_Molecule_Ref> Number of atoms differs.')
        end if
        molecule%atom_refcoord(1:3,1:natom) = &
                                  groref%atom_coord(1:3,1:natom)*10.0_wp
      end if
    end if

    return

  end subroutine setup_molecule_ref

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_fit
  !> @brief        setup molecules by Fit information
  !! @authors      KT
  !! @param[inout] molecule : molecule information
  !! @param[in]    fit      : fit information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_fit(molecule, fit)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_pdb),    optional, intent(in)    :: fit

    ! local variables
    integer                  :: natom


    natom = molecule%num_atoms

    if (present(fit)) then
      if (allocated(fit%atom_coord)) then
        if (natom /= size(fit%atom_coord(1,:))) then
          write(MsgOut,*) 'Setup_Molecule_Fit> ATOM in Molecule = ', natom
          write(MsgOut,*) 'Setup_Molecule_Fit> ATOM in PDB FIT  = ', &
                                                       size(fit%atom_coord(1,:))
          call error_msg('Setup_Molecule_Fit> Number of atoms differs.')
        end if

        call alloc_molecules(molecule, MoleculeFitc, natom)

        molecule%atom_fitcoord(1:3,1:natom) = fit%atom_coord(1:3,1:natom)
      end if
    end if

    return

  end subroutine setup_molecule_fit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_mode
  !> @brief        setup molecules by Mode information
  !! @authors      YK
  !! @param[inout] molecule : molecule information
  !! @param[in]    mode     : mode information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_mode(molecule, mode)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_mode),   optional, intent(in)   :: mode

    ! local variables
    integer                  :: nmode


    nmode = mode%num_pc_modes
  
    call alloc_molecules(molecule, MoleculeMode, nmode)

    molecule%pc_mode(:)   = mode%pc_mode(:)
    molecule%num_pc_modes = mode%num_pc_modes

    return

  end subroutine setup_molecule_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_molecule_other
  !> @brief        setup other molecule informations
  !! @authors      NT
  !! @param[inout] molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_molecule_other(molecule)

    ! formal arguments
    !
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    !
    real(wp)                 :: mass
    real(dp)                 :: charge
    integer                  :: i, j
    integer                  :: nmol, imole, imole_pre, initial, final
    integer                  :: col, pcol
    integer                  :: src_col, dst_col, cur_col
    integer                  :: i1, i2
    integer                  :: resno_pre, resno_next
    character(4)             :: segment_pre, segment_next
    character(6)             :: residue_pre, residue_next

    integer, allocatable     :: color(:), color_range(:,:)


    ! total charge
    !
    charge = 0.0_dp
    do i = 1, molecule%num_atoms
      charge = charge + molecule%charge(i)
    end do
    molecule%total_charge = charge


    ! count the number of segments
    !
    do i = 1, molecule%num_atoms
      if (i == 1) then
        molecule%segment_no(1) = 1
        segment_pre = molecule%segment_name(1)
      else
        if (molecule%segment_name(i) /= segment_pre) then
          molecule%segment_no(i) = molecule%segment_no(i-1) + 1
          segment_pre            = molecule%segment_name(i)
        else
          molecule%segment_no(i) = molecule%segment_no(i-1)
        end if
      end if
    end do
    molecule%num_segments = molecule%segment_no(molecule%num_atoms)


    ! count the number of residues
    !
    do i = 1, molecule%num_atoms
      if (i == 1) then
        segment_pre = molecule%segment_name(i)
        residue_pre = molecule%residue_name(i)
        resno_pre   = molecule%residue_no(i)
        molecule%residue_c_no(i) = 1
      else
        segment_pre  = molecule%segment_name(i-1)
        residue_pre  = molecule%residue_name(i-1)
        resno_pre    = molecule%residue_no(i-1)
        segment_next = molecule%segment_name(i)
        residue_next = molecule%residue_name(i)
        resno_next   = molecule%residue_no(i)
        
        if ((segment_pre /= segment_next) .or. &
            (residue_pre /= residue_next) .or. &
            (resno_pre   /= resno_next)) then
          molecule%residue_c_no(i) = molecule%residue_c_no(i-1) + 1
        else
          molecule%residue_c_no(i) = molecule%residue_c_no(i-1)
        end if
      end if
    end do
    molecule%num_residues = molecule%residue_c_no(molecule%num_atoms)


    ! count the number of molecule
    !

    allocate(color(molecule%num_atoms), color_range(2,molecule%num_atoms/2))
    color(:)         = 0
    color_range(:,:) = 0
    
    ! coloring 
    cur_col = 0
    do i = 1, molecule%num_bonds

      i1 = molecule%bond_list(1,i)
      i2 = molecule%bond_list(2,i)

      if (color(i1) == 0 .and. color(i2) == 0) then

        cur_col = cur_col + 1
        color(i1) = cur_col
        color(i2) = cur_col

        color_range(1,cur_col) = min(i1, i2)
        color_range(2,cur_col) = max(i1, i2)

      else if (color(i1) == 0 .or. color(i2) == 0) then

        if (color(i1) /= 0) then
          
          col = color(i1)
          color(i2) = col
          color_range(1,col) = min(color_range(1,col),i2)
          color_range(2,col) = max(color_range(2,col),i2)

        else ! color(i2) /= 0

          col = color(i2)
          color(i1) = col
          color_range(1,col) = min(color_range(1,col),i1)
          color_range(2,col) = max(color_range(2,col),i1)

        end if

      else if (color(i1) /= color(i2)) then

        if (color(i1) < color(i2)) then
          src_col = color(i1)
          dst_col = color(i2)
        else
          src_col = color(i2)
          dst_col = color(i1)
        end if

        color_range(1,src_col) = min(color_range(1,src_col), &
                                     color_range(1,dst_col))
        color_range(2,src_col) = max(color_range(2,src_col), &
                                     color_range(2,dst_col))

        do j = color_range(1,dst_col), color_range(2,dst_col)
          if (color(j) == dst_col) &
            color(j) = src_col
        end do

        color_range(1,dst_col) = 0
        color_range(2,dst_col) = 0

      end if

    end do

    ! define molecule number
    nmol = 0
    pcol = -1

    do i = 1, molecule%num_atoms

      col = color(i)

      if (pcol /= col .or. col == 0) &
           nmol = nmol + 1

      molecule%molecule_no(i) = nmol

      pcol = col
    end do

    deallocate(color, color_range)

    molecule%num_molecules = nmol


    ! allocate molecule information
    !

    call alloc_molecules(molecule, MoleculeMole, nmol)


    ! define information per molecule
    !

    do i = 1, molecule%num_atoms
      if (i == 1) then
        imole = molecule%molecule_no(i)
        molecule%molecule_atom_no(imole) = i
      else
        imole_pre = molecule%molecule_no(i-1)  
        imole     = molecule%molecule_no(i)
        if (imole /= imole_pre) then
          molecule%molecule_atom_no(imole) = i
        end if
      end if
      molecule%molecule_name(imole) = molecule%segment_name(i)
    end do

    do i = 1, molecule%num_molecules
      initial  = molecule%molecule_atom_no(i)
      if (i /= molecule%num_molecules) then
        final = molecule%molecule_atom_no(i+1) - 1
      else
        final = molecule%num_atoms
      end if

      mass = 0.0_wp
      do j = initial, final
        mass = mass + molecule%mass(j)
      end do
      molecule%molecule_mass(i) = mass

    end do

    return

  end subroutine setup_molecule_other

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initialize_velocity
  !> @brief        initialize velocity of atoms
  !! @authors      TM
  !! @param[inout] molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initialize_velocity(molecule)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: natom


    natom = molecule%num_atoms
    molecule%atom_velocity(1:3,1:natom) = 0.0_wp

    return

  end subroutine initialize_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecule_sel
  !> @brief        export molecule to pdb
  !! @authors      NT
  !! @param[in]    molecule : structure of molecule
  !! @param[in]    selatoms : information of atom selection
  !! @param[out]   pdb      : PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecule_sel(molecule, selatoms, pdb)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_selatoms),        intent(in)    :: selatoms
    type(s_pdb),             intent(inout) :: pdb

    ! local variables
    integer                  :: i, idx, natom, nresi, nter, mole1, mole2
    integer,     allocatable :: ter_no(:)


    natom = size(selatoms%idx)
    nresi = 0
    nter  = 0

    allocate(ter_no(natom + 1))

    ! init pdb
    call init_pdb(pdb)

    pdb%hetatm_rec = .false.
    pdb%cryst_rec  = .false.
    pdb%model_rec  = .true.
    pdb%ter_rec    = .true.
    pdb%end_rec    = .true.
    pdb%segment    = .true.

    ! setup atom data
    call alloc_pdb(pdb, PdbAtom, natom)

    do i = 1, natom

      idx = selatoms%idx(i)

      pdb%hetatm(i)           = .false.
      pdb%atom_name(i)        = molecule%atom_name(idx)
      pdb%residue_name(i)     = molecule%residue_name(idx)(1:4)
      pdb%segment_name(i)     = molecule%segment_name(idx)
      pdb%chain_id(i)         = molecule%chain_id(idx)
      pdb%atom_no(i)          = i
      pdb%residue_no(i)       = molecule%residue_no(idx)
      pdb%atom_coord(:,i)     = molecule%atom_coord(:,idx)
      pdb%atom_occupancy(i)   = 1.0_wp
      pdb%atom_temp_factor(i) = 0.0_wp

      nresi = max(nresi, pdb%residue_no(i))

      if (i /= natom) then

        mole1 = molecule%molecule_no(idx)
        mole2 = molecule%molecule_no(selatoms%idx(i + 1))

        if (mole1 /= mole2) then
          nter = nter + 1
          ter_no(nter) = i + 1
        end if

      end if

    end do

    nter = nter + 1
    ter_no(nter) = natom + 1

    ! setup ter data
    call alloc_pdb(pdb, PdbTerLine, nter)

    if (nter >= 1) &
    pdb%ter_line_no(1:nter) = ter_no(1:nter)


    if (natom > 999999) &
      pdb%atom_col7 = .true.

    if (nresi > 9999) then
      pdb%res_col5 = .true.
      if (nresi > 99999) &
        pdb%res_col6 = .true.
    end if

    deallocate(ter_no)

    return

  end subroutine export_molecule_sel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecule_all
  !> @brief        export molecule to pdb
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[out]   pdb      : PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecule_all(molecule, pdb)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_pdb),             intent(inout) :: pdb

    ! local variables
    integer                  :: i, natom, nresi, nter, mole1, mole2
    integer,     allocatable :: ter_no(:)


    natom = size(molecule%atom_name)
    nresi = 0
    nter  = 0

    allocate(ter_no(natom + 1))

    ! init pdb
    call init_pdb(pdb)

    pdb%hetatm_rec = .false.
    pdb%cryst_rec  = .false.
    pdb%model_rec  = .true.
    pdb%ter_rec    = .true.
    pdb%end_rec    = .true.
    pdb%segment    = .true.

    ! setup atom data
    call alloc_pdb(pdb, PdbAtom, natom)

    do i = 1, natom

      pdb%hetatm(i)           = .false.
      pdb%atom_name(i)        = molecule%atom_name(i)
      pdb%residue_name(i)     = molecule%residue_name(i)(1:4)
      pdb%segment_name(i)     = molecule%segment_name(i)
      pdb%chain_id(i)         = molecule%chain_id(i)
      pdb%atom_no(i)          = i
      pdb%residue_no(i)       = molecule%residue_no(i)
      pdb%atom_coord(:,i)     = molecule%atom_coord(:,i)
      pdb%atom_occupancy(i)   = 1.0_wp
      pdb%atom_temp_factor(i) = 0.0_wp

      nresi = max(nresi, pdb%residue_no(i))

      if (i /= natom) then

        mole1 = molecule%molecule_no(i)
        mole2 = molecule%molecule_no(i + 1)

        if (mole1 /= mole2) then
          nter = nter + 1
          ter_no(nter) = i + 1
        end if

      end if

    end do

    nter = nter + 1
    ter_no(nter) = natom + 1

    ! setup ter data
    call alloc_pdb(pdb, PdbTerLine, nter)

    if (nter >= 1) &
    pdb%ter_line_no(1:nter) = ter_no(1:nter)


    if (natom > 999999) &
      pdb%atom_col7 = .true.

    if (nresi > 9999) then
      pdb%res_col5 = .true.
      if (nresi > 99999) &
        pdb%res_col6 = .true.
    end if

    deallocate(ter_no)

    return

  end subroutine export_molecule_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_light_atom_name
  !> @brief        check name of light atoms
  !! @authors      CK
  !! @param[inout] molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_light_atom_name(mass_upper_bound, molecule)

    ! formal arguments
    real(wp),                intent(in   ) :: mass_upper_bound
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i
    integer                  :: natom


    natom = molecule%num_atoms
    do i = 1, natom
      if (molecule%mass(i) .gt. mass_upper_bound) cycle
      molecule%light_atom_mass(i) = .true. 
      if (molecule%atom_name(i)(1:1) .eq. "H" .or. &
          molecule%atom_name(i)(1:1) .eq. "D" .or. &
          molecule%atom_name(i)(1:1) .eq. "h" .or. &
          molecule%atom_name(i)(1:1) .eq. "d") then
        molecule%light_atom_name(i) = .true. 
      endif
      if (molecule%atom_cls_name(i)(1:1) .eq. "H" .or. &
          molecule%atom_cls_name(i)(1:1) .eq. "D" .or. &
          molecule%atom_cls_name(i)(1:1) .eq. "h" .or. &
          molecule%atom_cls_name(i)(1:1) .eq. "d") then
        molecule%light_atom_name(i) = .true. 
      endif

      if (.not. molecule%light_atom_name(i)) then
!        if (main_rank) &
!        write(MsgOut,'(A,i10,A,A4,A6)') &
!           'Check_Light_Atom_Name> ATOM ',i,' has non ordinary name,', &
!           molecule%atom_name(i),molecule%atom_cls_name(i)
        molecule%special_hydrogen = .true.
      endif

    end do

    return

  end subroutine check_light_atom_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecule_sel_crd
  !> @brief        export molecule to crd
  !! @authors      KYMD
  !! @param[in]    molecule : structure of molecule
  !! @param[in]    selatoms : information of atom selection
  !! @param[out]   crd      : CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecule_sel_crd(molecule, selatoms, crd)

    ! formal arguments
    type(s_molecule),        intent(in)  :: molecule
    type(s_selatoms),        intent(in)  :: selatoms
    type(s_crd),             intent(out) :: crd

    ! local variables
    integer                  :: i, idx, natom, residue_no_base, last_residue_no


    natom = size(selatoms%idx)

    ! init crd
    call init_crd(crd)

    ! setup atom data
    call alloc_crd(crd, CrdAtom, natom)

    crd%num_atoms   = natom
    residue_no_base = molecule%residue_no(selatoms%idx(1)) - 1
    last_residue_no = residue_no_base + 1
    do i = 1, natom

      idx = selatoms%idx(i)
      if(molecule%residue_no(idx) /= last_residue_no) then
        residue_no_base = residue_no_base + (molecule%residue_no(idx) - last_residue_no) - 1
        last_residue_no = molecule%residue_no(idx)
      end if

      crd%atom_name(i)        = molecule%atom_name(idx)
      crd%residue_name(i)     = molecule%residue_name(idx) 
      crd%segment_name(i)     = molecule%segment_name(idx)
      crd%atom_no(i)          = i 
      crd%residue_no(i)       = molecule%residue_no(idx) - residue_no_base
      crd%residue_id(i)       = molecule%residue_no(idx)
      crd%atom_coord(:,i)     = molecule%atom_coord(:,idx)

    end do

    return

  end subroutine export_molecule_sel_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecule_all_crd
  !> @brief        export molecule to crd
  !! @authors      KYMD
  !! @param[in]    molecule : molecule information
  !! @param[out]   crd      : crd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecule_all_crd(molecule, crd)

    ! formal arguments
    type(s_molecule),        intent(in)  :: molecule
    type(s_crd),             intent(out) :: crd

    ! local variables
    integer                  :: i, natom, residue_no_base


    natom = size(molecule%atom_name)

    ! init pdb
    call init_crd(crd)

    ! setup atom data
    call alloc_crd(crd, CrdAtom, natom)

    crd%num_atoms   = natom
    residue_no_base = molecule%residue_no(1) - 1
    do i = 1, natom

      crd%atom_name(i)        = molecule%atom_name(i)
      crd%residue_name(i)     = molecule%residue_name(i)(1:4)
      crd%segment_name(i)     = molecule%segment_name(i)
      crd%atom_no(i)          = i
      crd%residue_no(i)       = molecule%residue_no(i) - residue_no_base
      crd%residue_id(i)       = molecule%residue_no(i)
      crd%atom_coord(:,i)     = molecule%atom_coord(:,i)

    end do

    return

  end subroutine export_molecule_all_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    export_molecule_psf
  !> @brief        export molecule to crd
  !! @authors      NT, KY
  !! @param[in]    molecule : structure of molecule
  !! @param[in]    selatoms : information of atom selection
  !! @param[inout] psf      : PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine export_molecule_sel_psf(molecule, selatoms, psf)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_selatoms),        intent(in)    :: selatoms
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i, j, jfirst, jlast
    integer                  :: natom, nbond, ntheta, nphi, nimphi,  &
                                ndon, nacc, nnb, ngrp, nst2, molnt,  &
                                numlp, numlph, ncrterm
    integer                  :: idx
    integer                  :: iatom, ibond, itheta, iphi, iimphi,  &
                                idon, iacc, inb, igrp, ist2, imolnt, &
                                inumlp, inumlph, icrterm
    type(s_psf)              :: tmp_psf
    real(wp)                 :: total_charge
    integer                  :: itmp(8), item 
    integer, allocatable     :: sel_idx(:), inv_sel_idx(:)
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: isGRP


    natom = size(selatoms%idx)
    iatom = psf%num_atoms
    if (natom == iatom .and. psf%total_charge == 0.0_wp) return

    allocate(sel_idx(1: natom), inv_sel_idx(0: iatom), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    sel_idx(1: natom)     = selatoms%idx(1: natom)
    inv_sel_idx(0: iatom) = 0

    ! init psf
    !
    call init_psf(tmp_psf)

    ! setup atom data
    !
    tmp_psf%type = psf%type
    ibond        = psf%num_bonds
    itheta       = psf%num_angles
    iphi         = psf%num_dihedrals
    iimphi       = psf%num_impropers
    igrp         = psf%num_groups
    icrterm      = psf%num_cross_terms
    idon         = psf%num_HB_donors
    iacc         = psf%num_HB_acceptors
    inb          = psf%num_NB_exclusions
    call alloc_psf(tmp_psf, PsfAtom, natom)
    call alloc_psf(tmp_psf, PsfBond, ibond)
    call alloc_psf(tmp_psf, PsfAngl, itheta)
    call alloc_psf(tmp_psf, PsfDihe, iphi)
    call alloc_psf(tmp_psf, PsfImpr, iimphi)
    call alloc_psf(tmp_psf, PsfDonr, idon)
    call alloc_psf(tmp_psf, PsfAcce, iacc)
    call alloc_psf(tmp_psf, PsfNb, inb)
    call alloc_psf(tmp_psf, PsfGrp, igrp)
    call alloc_psf(tmp_psf, PsfCmap, icrterm)

    !  Atomic information
    !
    total_charge             = 0.0_wp
    do i = 1, natom 
      idx = sel_idx(i)

      tmp_psf%atom_no(i)        = i ! psf%atom_no(idx)
      tmp_psf%segment_name(i)   = psf%segment_name(idx)
      tmp_psf%residue_no(i)     = psf%residue_no(idx) 
      tmp_psf%residue_name(i)   = psf%residue_name(idx)
      tmp_psf%atom_name(i)      = psf%atom_name(idx)
      tmp_psf%atom_cls_name(i)  = psf%atom_cls_name(idx)
      tmp_psf%atom_cls_no(i)    = psf%atom_cls_no(idx)
      tmp_psf%charge(i)         = psf%charge(idx)
      tmp_psf%mass(i)           = psf%mass(idx)
      tmp_psf%imove(i)          = psf%imove(idx)
      tmp_psf%molecule_no(i)    = psf%molecule_no(idx) 

      total_charge              = total_charge + psf%charge(idx)
      inv_sel_idx(idx)          = i
    end do

    ! Bond
    !
    nbond = 0
    item  = 2
    do i = 1, ibond
      itmp(1: item) = psf%bond_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      nbond = nbond + 1
      do j = 1, item
        tmp_psf%bond_list(j, nbond) = inv_sel_idx(itmp(j))
      end do
    end do

    ! Angle
    !
    ntheta = 0
    item   = 3
    do i = 1, itheta
      itmp(1: item) = psf%angl_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      ntheta = ntheta + 1
      do j = 1, item
        tmp_psf%angl_list(j, ntheta) = inv_sel_idx(itmp(j))
      end do
    end do

    ! Dihedral
    !
    nphi = 0
    item = 4
    do i = 1, iphi
      itmp(1: item) = psf%dihe_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      nphi = nphi + 1
      do j = 1, item
        tmp_psf%dihe_list(j, nphi) = inv_sel_idx(itmp(j))
      end do
    end do

    ! Improper
    !
    nimphi = 0
    item   = 4
    do i = 1, iimphi
      itmp(1: item) = psf%impr_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      nimphi = nimphi + 1
      do j = 1, item
        tmp_psf%impr_list(j, nimphi) = inv_sel_idx(itmp(j))
      end do
    end do

    ! Donor
    !
    ndon = 0
    item = 2
    do i = 1, idon
      itmp(1: item) = psf%donr_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      ndon = ndon + 1
      do j = 1, item
        tmp_psf%donr_list(j, ndon) = inv_sel_idx(itmp(j))
      end do
    end do

    ! Acceptor
    !
    nacc = 0
    item = 2
    do i = 1, iacc
      itmp(1: item) = psf%acce_list(1: item, i)
      if (.not.any(sel_idx(1: natom) == itmp(1))) cycle
      if (.not.any(sel_idx(1: natom) == itmp(2)) .and. itmp(2) /= 0) cycle
      nacc = nacc + 1
      do j = 1, item
        tmp_psf%acce_list(j, nacc) = inv_sel_idx(itmp(j))
      end do
    end do

    ! NB
    !
    nnb  = 0
    item = 1
    do i = 1, inb
      itmp(1) = psf%nb_list(i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if(j /= item + 1) cycle
      nnb = nnb + 1
      tmp_psf%nb_list(nnb) = inv_sel_idx(itmp(j))
    end do

    ! Group
    !
    ngrp = 0
    item = 3
    do i = 1, igrp 
      itmp(1: item) = psf%grp_list(1: item, i)
      jfirst = itmp(1) + 1
      if (i /= igrp) then
        jlast = psf%grp_list(1, i + 1) 
      else
        jlast = iatom
      end if
      isGRP = .true.
      do j = jfirst, jlast
        if (inv_sel_idx(j) == 0) then
          isGRP = .false.
          exit
        end if
      end do
      if (isGRP) then
        ngrp = ngrp + 1
        tmp_psf%grp_list(1, ngrp)       = inv_sel_idx(jfirst)-1
        tmp_psf%grp_list(2: item, ngrp) = itmp(2: item)
      end if
    end do

    ! Cmap
    !
    ncrterm = 0
    item    = 8
    do i = 1, icrterm
      itmp(1: item) = psf%cmap_list(1: item, i)
      do j = 1, item
        if (.not.any(sel_idx(1: natom) == itmp(j))) exit
      end do
      if (j /= item + 1) cycle
      ncrterm = ncrterm + 1
      do j = 1, item
        tmp_psf%cmap_list(j, ncrterm) = inv_sel_idx(itmp(j))
      end do
    end do


    ! substitution with array size changed
    !
    call dealloc_psf_all(psf)
    call init_psf(psf)
    psf%type              = tmp_psf%type
    psf%num_atoms         = natom
    psf%num_bonds         = nbond
    psf%num_angles        = ntheta
    psf%num_dihedrals     = nphi
    psf%num_impropers     = nimphi
    psf%num_groups        = ngrp
    psf%num_cross_terms   = ncrterm
    psf%num_HB_donors     = ndon
    psf%num_HB_acceptors  = nacc
    psf%num_NB_exclusions = nnb
    psf%total_charge      = total_charge
    call alloc_psf(psf, PsfAtom, natom)
    call alloc_psf(psf, PsfBond, nbond)
    call alloc_psf(psf, PsfAngl, ntheta)
    call alloc_psf(psf, PsfDihe, nphi)
    call alloc_psf(psf, PsfImpr, nimphi)
    call alloc_psf(psf, PsfDonr, ndon)
    call alloc_psf(psf, PsfAcce, nacc)
    call alloc_psf(psf, PsfNb, nnb)
    call alloc_psf(psf, PsfGrp, ngrp)
    call alloc_psf(psf, PsfCmap, ncrterm)

    do i = 1, natom
      psf%atom_no(i)        = i ! tmp_psf%atom_no(i)
      psf%segment_name(i)   = tmp_psf%segment_name(i)
      psf%residue_no(i)     = tmp_psf%residue_no(i) 
      psf%residue_name(i)   = tmp_psf%residue_name(i)
      psf%atom_name(i)      = tmp_psf%atom_name(i)
      psf%atom_cls_name(i)  = tmp_psf%atom_cls_name(i)
      psf%atom_cls_no(i)    = tmp_psf%atom_cls_no(i)
      psf%charge(i)         = tmp_psf%charge(i)
      psf%mass(i)           = tmp_psf%mass(i)
      psf%imove(i)          = tmp_psf%imove(i)
      psf%molecule_no(i)    = tmp_psf%molecule_no(i) 
    end do

    do i = 1, nbond
      psf%bond_list(1: 2, i) = tmp_psf%bond_list(1: 2, i)
    end do

    do i = 1, ntheta
      psf%angl_list(1: 3, i) = tmp_psf%angl_list(1: 3, i)
    end do

    do i = 1, nphi
      psf%dihe_list(1: 4, i) = tmp_psf%dihe_list(1: 4, i)
    end do

    do i = 1, nimphi
      psf%impr_list(1: 4, i) = tmp_psf%impr_list(1: 4, i)
    end do

    do i = 1, ndon
      psf%donr_list(1: 2, i) = tmp_psf%donr_list(1: 2, i)
    end do

    do i = 1, nacc
      psf%acce_list(1: 2, i) = tmp_psf%acce_list(1: 2, i)
    end do

    do i = 1, nnb
      psf%nb_list(i) = tmp_psf%nb_list(i)
    end do

    do i = 1, ngrp
      psf%grp_list(1: 3, i) = tmp_psf%grp_list(1: 3, i)
    end do

    do i = 1, ncrterm
      psf%cmap_list(1: 8, i) = tmp_psf%cmap_list(1: 8, i)
    end do

    call dealloc_psf_all(tmp_psf)

    deallocate(sel_idx, inv_sel_idx, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine export_molecule_sel_psf

end module molecules_mod
