!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_amber2hm_mod
!> @brief   create huge molecule from AMBER data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_amber2hm_mod

  use pr_pdb2hm_mod
  use pr_huge_molecule_mod
  use fileio_prmtop_mod
  use fileio_ambcrd_mod
  use fileio_pdb_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: hm_create_amber
  private :: hm_create
  private :: hm_create_atom
  private :: hm_create_bond
  private :: hm_create_angl
  private :: hm_create_dihe
  private :: hm_create_impr
  private :: hm_create_ambcrd_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_amber
  !> @brief        create huge molecule from AMBER parameter/topology file
  !! @authors      NT
  !! @param[in]    prmtop     : AMBER parameter/topology information
  !! @param[in]    ambcrdfile : AMBER coordiante file
  !! @param[in]    ambreffile : Reference AMBER coordinate file
  !! @param[in]    pdbfile    : PDB file name (optional)
  !! @param[in]    reffile    : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_amber(prmtop, ambcrdfile, ambreffile, pdbfile, reffile)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    character(*),            intent(in)    :: ambcrdfile
    character(*),            intent(in)    :: ambreffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: ambcrd_in, ambref_in, pdb_in, ref_in


    ambcrd_in = 0
    ambref_in = 0
    pdb_in    = 0
    ref_in    = 0

    if (ambcrdfile /= '') &
      call open_file(ambcrd_in, ambcrdfile, IOFileInput)

    if (ambreffile /= '') &
      call open_file(ambref_in, ambreffile, IOFileInput)

    if (pdbfile /= '') &
      call open_file(pdb_in, pdbfile, IOFileInput)

    if (reffile /= '') &
      call open_file(ref_in, reffile, IOFileInput)


    call hm_create(prmtop, ambcrd_in, ambref_in, pdb_in, ref_in)


    if (ref_in /= 0) &
      call close_file(ref_in)

    if (pdb_in /= 0) &
      call close_file(pdb_in)

    if (ambref_in /= 0) &
      call close_file(ambref_in)

    if (ambcrd_in /= 0) &
      call close_file(ambcrd_in)


    if (main_rank) then
      write(MsgOut,'(a)')    'Hm_Create_Amber> Molecule information: '
      write(MsgOut,'(a,i9)') '  num_atoms       = ', hm_num_atoms
      write(MsgOut,'(a,i9)') '  num_bonds       = ', hm_num_bonds
      write(MsgOut,'(a,i9)') '  num_angles      = ', hm_num_angles
      write(MsgOut,'(a,i9)') '  num_dihedrals   = ', hm_num_dihedrals
      write(MsgOut,'(a,i9)') '  num_impropers   = ', hm_num_impropers
      write(MsgOut,'(a)')    ' '
    end if

    return

  end subroutine hm_create_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create
  !> @brief        create huge molecule files
  !! @authors      NT
  !! @param[in]    prmtop    : AMBER parameter/topology information
  !! @param[in]    ambcrd_in : AMBCRD file unit number
  !! @param[in]    ambref_in : AMBREF file unit number
  !! @param[in]    pdb_in    : PDB file unit number
  !! @param[in]    ref_in    : REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create(prmtop, ambcrd_in, ambref_in, pdb_in, ref_in)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    integer,                 intent(in)    :: ambcrd_in
    integer,                 intent(in)    :: ambref_in
    integer,                 intent(in)    :: pdb_in
    integer,                 intent(in)    :: ref_in


    ! create huge molecule from grotop information
    !

    call hm_create_atom(prmtop)

    call hm_create_bond(prmtop)

    call hm_create_angl(prmtop)

    call hm_create_dihe(prmtop)

    call hm_create_impr(prmtop)


    ! read ambcrd / pdb
    !

    if (ambcrd_in /= 0) then

      call hm_create_ambcrd_coord(ambcrd_in, .false.)

    else if (pdb_in /= 0) then

      call hm_create_pdb_coord(pdb_in, .false.)

    end if


    ! read reference pdb
    !

    if (ambref_in /= 0) then

      call hm_create_ambcrd_coord(ambref_in, .true.)

    else if (ref_in /= 0) then

      call hm_create_pdb_coord(ref_in, .true.)

    end if


    if (main_rank) &
      write(MsgOut,'(a)') ''


    return

  end subroutine hm_create

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_atom
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_atom(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, iatom, i, is, ie
    character(4)             :: atom_name
    character(6)             :: atom_cls_name


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Atom> '

    call hm_alloc_atoms(prmtop%num_atoms)

    do iatom = 1, prmtop%num_atoms

      call hm_set_atom_no      (iatom, iatom)
      call hm_set_atom_cls_no  (iatom, prmtop%atom_cls_no  (iatom))
      call hm_set_atom_cls_name(iatom, prmtop%atom_cls_name(iatom))
      call hm_set_atom_name    (iatom, prmtop%atom_name    (iatom))
      call hm_set_charge       (iatom, prmtop%charge       (iatom)/18.2223_wp)
      call hm_set_mass         (iatom, real(prmtop%mass    (iatom),wip))
      call hm_set_segment_name (iatom, '')

    end do

    do i = 1, prmtop%num_residues

      is = prmtop%res_point(i)
      if (i < prmtop%num_residues) then
        ie = prmtop%res_point(i+1)-1
      else
        ie = prmtop%num_atoms
      end if

      do iatom = is, ie
        call hm_set_residue_no  (iatom, i)
        call hm_set_residue_name(iatom, prmtop%res_label(i))
      end do

    end do

    ! check light atom
    !
    do iatom = 1, prmtop%num_atoms
      if (hm_mass(iatom) > LIGHT_ATOM_MASS_LIMIT) cycle
      call hm_set_light_atom_mass(iatom, 1)
      atom_name = hm_atom_name(iatom)
      atom_cls_name = hm_atom_cls_name(iatom)
      if (atom_name(1:1) == 'H' .or. atom_name(1:1) == 'D' .or. &
          atom_name(1:1) == 'h' .or. atom_name(1:1) == 'd')     &
        call hm_set_light_atom_name(iatom, 1)
      if (atom_cls_name(1:1) == 'H' .or. atom_cls_name(1:1) == 'D' .or. &
          atom_cls_name(1:1) == 'h' .or. atom_cls_name(1:1) == 'd')     &
        call hm_set_light_atom_name(iatom, 1)
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_bond
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_bond(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: ibond, i


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Bond> '

    call hm_alloc_bonds(prmtop%num_bondh + prmtop%num_mbonda)

    ibond = 1

    do i = 1, prmtop%num_bondh

      call hm_set_bond_list(1, ibond, prmtop%bond_inc_hy(1,i) / 3 + 1)
      call hm_set_bond_list(2, ibond, prmtop%bond_inc_hy(2,i) / 3 + 1)

      ibond = ibond + 1

    end do

    do i = 1, prmtop%num_mbonda

      call hm_set_bond_list(1, ibond, prmtop%bond_wo_hy(1,i) / 3 + 1)
      call hm_set_bond_list(2, ibond, prmtop%bond_wo_hy(2,i) / 3 + 1)

      ibond = ibond + 1

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_angl
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_angl(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: iangl, i


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Angl> '

    call hm_alloc_angls(prmtop%num_anglh + prmtop%num_mangla)

    iangl = 1

    do i = 1, prmtop%num_anglh

      call hm_set_angl_list(1, iangl, prmtop%angl_inc_hy(1,i) / 3 + 1)
      call hm_set_angl_list(2, iangl, prmtop%angl_inc_hy(2,i) / 3 + 1)
      call hm_set_angl_list(3, iangl, prmtop%angl_inc_hy(3,i) / 3 + 1)

      iangl = iangl + 1

    end do

    do i = 1, prmtop%num_mangla

      call hm_set_angl_list(1, iangl, prmtop%angl_wo_hy(1,i) / 3 + 1)
      call hm_set_angl_list(2, iangl, prmtop%angl_wo_hy(2,i) / 3 + 1)
      call hm_set_angl_list(3, iangl, prmtop%angl_wo_hy(3,i) / 3 + 1)

      iangl = iangl + 1

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_dihe
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_dihe(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: ndihe, idihe, i


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Dihe> '

    ndihe = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) >= 0) &
        ndihe = ndihe + 1
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) >= 0) &
        ndihe = ndihe + 1
    end do

    call hm_alloc_dihes(ndihe)

    idihe = 1

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      call hm_set_dihe_list(1, idihe,      prmtop%dihe_inc_hy(1,i)  / 3 + 1)
      call hm_set_dihe_list(2, idihe,      prmtop%dihe_inc_hy(2,i)  / 3 + 1)
      call hm_set_dihe_list(3, idihe, iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1)
      call hm_set_dihe_list(4, idihe,      prmtop%dihe_inc_hy(4,i)  / 3 + 1)

      idihe = idihe + 1

    end do

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      call hm_set_dihe_list(1, idihe,      prmtop%dihe_wo_hy(1,i)  / 3 + 1)
      call hm_set_dihe_list(2, idihe,      prmtop%dihe_wo_hy(2,i)  / 3 + 1)
      call hm_set_dihe_list(3, idihe, iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1)
      call hm_set_dihe_list(4, idihe,      prmtop%dihe_wo_hy(4,i)  / 3 + 1)

      idihe = idihe + 1

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_impr
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_impr(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: nimpr, iimpr, i


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Impr> '

    nimpr = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        nimpr = nimpr + 1
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        nimpr = nimpr + 1
    end do

    call hm_alloc_imprs(nimpr)

    iimpr = 1

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) >= 0) &
        cycle

      call hm_set_impr_list(1, iimpr,      prmtop%dihe_inc_hy(1,i)  / 3 + 1)
      call hm_set_impr_list(2, iimpr,      prmtop%dihe_inc_hy(2,i)  / 3 + 1)
      call hm_set_impr_list(3, iimpr, iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1)
      call hm_set_impr_list(4, iimpr, iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1)

      iimpr = iimpr + 1

    end do

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) >= 0) &
        cycle

      call hm_set_impr_list(1, iimpr,      prmtop%dihe_wo_hy(1,i)  / 3 + 1)
      call hm_set_impr_list(2, iimpr,      prmtop%dihe_wo_hy(2,i)  / 3 + 1)
      call hm_set_impr_list(3, iimpr, iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1)
      call hm_set_impr_list(4, iimpr, iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1)

      iimpr = iimpr + 1

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_ambcrd_coord
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_ambcrd_coord(crd_in, is_ref)

    ! formal arguments
    integer,                 intent(in)    :: crd_in
    logical,                 intent(in)    :: is_ref

    ! local variables
    real(wp)                 :: center(3), crd(3,2), atom_coord(3,2)
    integer                  :: natom, iatom, i, nr
    character(100)           :: line


    if (.not. is_ref) then
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Ambcrd_Coord> '
    else
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Ambcrd_Coord (REF)> '
    end if

    ! compute center of system
    read(crd_in,'(a)') line
    read(crd_in,*) natom

    center(1:3) = 0.0_wp
    do iatom = 1, natom, 2

      if (natom - iatom < 1) then
        nr = 1
      else
        nr = 2
      end if

      crd(1:3,1:2) = 0.0_wp
      read(crd_in,*) (crd(1:3,i),i=1,nr)
      center(1:3) = center(1:3) + crd(1:3,1)
      center(1:3) = center(1:3) + crd(1:3,2)

    end do
!   center(1:3) = center(1:3) / real(natom,wp)
    center(1:3) = 0.0_wp

    ! read crd file
    rewind(crd_in)

    read(crd_in,'(a)') line
    read(crd_in,*) natom

    do iatom = 1, natom, 2

      if (natom - iatom < 1) then
        nr = 1
      else
        nr = 2
      end if

      read(crd_in,*) (atom_coord(1:3,i),i=1,nr)

      do i = 1, nr
        if (.not. is_ref) then
          call hm_set_atom_coord &
               (1, iatom+i-1, real(atom_coord(1,i) - center(1),wip))
          call hm_set_atom_coord &
               (2, iatom+i-1, real(atom_coord(2,i) - center(2),wip))
          call hm_set_atom_coord &
               (3, iatom+i-1, real(atom_coord(3,i) - center(3),wip))
        else
          call hm_set_atom_refcoord &
               (1, iatom+i-1, real(atom_coord(1,i) - center(1),wip))
          call hm_set_atom_refcoord &
               (2, iatom+i-1, real(atom_coord(2,i) - center(2),wip))
          call hm_set_atom_refcoord &
               (3, iatom+i-1, real(atom_coord(3,i) - center(3),wip))
        end if
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Hm_Create_Ambcrd_Coord> read ERROR. ')

  end subroutine hm_create_ambcrd_coord

end module pr_amber2hm_mod
