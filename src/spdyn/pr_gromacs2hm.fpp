!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_gromacs2hm_mod
!> @brief   create huge molecule from GROMACS data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_gromacs2hm_mod

  use pr_pdb2hm_mod
  use pr_huge_molecule_mod
  use fileio_grotop_mod
  use fileio_grocrd_mod
  use fileio_pdb_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: hm_create_gromacs
  private :: hm_create
  private :: hm_create_atom
  private :: hm_create_bond
  private :: hm_create_angl
  private :: hm_create_dihe
  private :: hm_create_impr
  private :: hm_create_grocrd_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_gromacs
  !> @brief        create huge molecule from GROMACS structure/coordinate file
  !! @authors      NT
  !! @param[in]    grotop     : GROMACS topology information
  !! @param[in]    grocrdfile : GROMACS coordiante file
  !! @param[in]    groreffile : Reference GROMACS coordinate file
  !! @param[in]    pdbfile    : PDB file name (optional)
  !! @param[in]    reffile    : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_gromacs(grotop, grocrdfile, groreffile, pdbfile, reffile)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    character(*),            intent(in)    :: grocrdfile
    character(*),            intent(in)    :: groreffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: grocrd_in, groref_in, pdb_in, ref_in


    grocrd_in = 0
    groref_in = 0
    pdb_in    = 0
    ref_in    = 0

    if (grocrdfile /= '') &
      call open_file(grocrd_in, grocrdfile, IOFileInput)

    if (groreffile /= '') &
      call open_file(groref_in, groreffile, IOFileInput)

    if (pdbfile /= '') &
      call open_file(pdb_in, pdbfile, IOFileInput)

    if (reffile /= '') &
      call open_file(ref_in, reffile, IOFileInput)


    call hm_create(grotop, grocrd_in, groref_in, pdb_in, ref_in)


    if (ref_in /= 0) &
      call close_file(ref_in)

    if (pdb_in /= 0) &
      call close_file(pdb_in)

    if (groref_in /= 0) &
      call close_file(groref_in)

    if (grocrd_in /= 0) &
      call close_file(grocrd_in)


    if (main_rank) then
      write(MsgOut,'(a)')    'Hm_Create_Gromacs> Molecule information: '
      write(MsgOut,'(a,i9)') '  num_atoms       = ', hm_num_atoms
      write(MsgOut,'(a,i9)') '  num_bonds       = ', hm_num_bonds
      write(MsgOut,'(a,i9)') '  num_angles      = ', hm_num_angles
      write(MsgOut,'(a,i9)') '  num_dihedrals   = ', hm_num_dihedrals
      write(MsgOut,'(a,i9)') '  num_impropers   = ', hm_num_impropers
      write(MsgOut,'(a)')    ' '
    end if

    return

  end subroutine hm_create_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create
  !> @brief        create huge molecule data
  !! @authors      NT
  !! @param[in]    grotop    : GROMACS topology information
  !! @param[in]    grocrd_in : GROCRD file unit number
  !! @param[in]    groref_in : GROREF file unit number
  !! @param[in]    pdb_in    : PDB file unit number
  !! @param[in]    ref_in    : REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create(grotop, grocrd_in, groref_in, pdb_in, ref_in)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    integer,                 intent(in)    :: grocrd_in
    integer,                 intent(in)    :: groref_in
    integer,                 intent(in)    :: pdb_in
    integer,                 intent(in)    :: ref_in


    ! create huge molecule from grotop information
    !

    call hm_create_atom(grotop)

    call hm_create_bond(grotop)

    call hm_create_angl(grotop)

    call hm_create_dihe(grotop)

    call hm_create_impr(grotop)


    ! read grocrd / pdb
    !

    if (grocrd_in /= 0) then

      call hm_create_grocrd_coord(grocrd_in, .false.)

    else if (pdb_in /= 0) then

      call hm_create_pdb_coord(pdb_in, .false.)

    end if


    ! read reference pdb
    !

    if (groref_in /= 0) then

      call hm_create_grocrd_coord(grocrd_in, .true.)

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
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_atom(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: iatom, i, j, k, l, cls_no
    integer                  :: natom, nmole
    character(6)             :: cls_name, name

    type(s_grotop_mol),      pointer     :: gromol


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Atom> '

    ! count number of atoms
    natom = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        natom = natom + gromol%num_atoms
      end do
    end do

    call hm_alloc_atoms(natom)


    iatom = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count
      name   =  grotop%molss(i)%name

      do j = 1, nmole

        do k = 1, gromol%num_atoms

          cls_name = gromol%atoms(k)%atom_type
          do l = 1, grotop%num_atomtypes
            if (cls_name == grotop%atomtypes(l)%type_name) then
              cls_no = l
              exit
            end if
          end do
          if (l > grotop%num_atomtypes) &
            write(MsgOut,*) &
            '  WARNING: atom class name is undefined. '//trim(cls_name)

          call hm_set_atom_no      (iatom, iatom)
          call hm_set_atom_cls_no  (iatom, cls_no)
          call hm_set_atom_cls_name(iatom, cls_name)
          call hm_set_atom_name    (iatom, gromol%atoms(k)%atom_name)
          call hm_set_charge       (iatom, gromol%atoms(k)%charge)
          call hm_set_mass         (iatom, real(gromol%atoms(k)%mass,wip))
          call hm_set_residue_no   (iatom, gromol%atoms(k)%residue_no)
          call hm_set_residue_name (iatom, gromol%atoms(k)%residue_name)
          call hm_set_segment_name (iatom, name)

          iatom = iatom + 1

        end do

      end do

    end do

    ! check light atom
    !
    do iatom = 1, natom
      if (hm_mass(iatom) > LIGHT_ATOM_MASS_LIMIT) cycle
      call hm_set_light_atom_mass(iatom, 1)
      name = hm_atom_name(iatom)
      cls_name = hm_atom_cls_name(iatom)
      if (name(1:1) == 'H' .or. name(1:1) == 'D' .or. &
          name(1:1) == 'h' .or. name(1:1) == 'd')     &
        call hm_set_light_atom_name(iatom, 1)
      if (cls_name(1:1) == 'H' .or. cls_name(1:1) == 'D' .or. &
          cls_name(1:1) == 'h' .or. cls_name(1:1) == 'd')     &
        call hm_set_light_atom_name(iatom, 1)
    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_bond
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_bond(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: ibond, ioffset, i, j, k
    integer                  :: nbond, nmole

    type(s_grotop_mol),      pointer     :: gromol


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Bond> '

    ! count number of bonds
    nbond = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        nbond = nbond + gromol%num_bonds
      end do
    end do

    call hm_alloc_bonds(nbond)

    ioffset = 0
    ibond   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_bonds

          call hm_set_bond_list(1,ibond, gromol%bonds(k)%atom_idx1 + ioffset)
          call hm_set_bond_list(2,ibond, gromol%bonds(k)%atom_idx2 + ioffset)

          ibond = ibond + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_angl
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_angl(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: iangl, ioffset, i, j, k
    integer                  :: nangl, nmole

    type(s_grotop_mol),      pointer     :: gromol


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Angl> '

    ! count number of bonds
    nangl = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        nangl = nangl + gromol%num_angls
      end do
    end do

    call hm_alloc_angls(nangl)

    ioffset = 0
    iangl   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_angls

          call hm_set_angl_list(1, iangl, gromol%angls(k)%atom_idx1 + ioffset)
          call hm_set_angl_list(2, iangl, gromol%angls(k)%atom_idx2 + ioffset)
          call hm_set_angl_list(3, iangl, gromol%angls(k)%atom_idx3 + ioffset)

          iangl = iangl + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_dihe
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_dihe(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idihe, ioffset, i, j, k
    integer                  :: ndihe, nmole

    type(s_grotop_mol),      pointer     :: gromol


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Dihe> '

    ! count number of dihedrals
    ndihe = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func == 1 .or. gromol%dihes(k)%func == 3) then
            ndihe = ndihe + 1
          end if
        end do
      end do

    end do

    call hm_alloc_dihes(ndihe)

    ioffset = 0
    idihe   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 3) &
            cycle

          call hm_set_dihe_list(1, idihe, gromol%dihes(k)%atom_idx1 + ioffset)
          call hm_set_dihe_list(2, idihe, gromol%dihes(k)%atom_idx2 + ioffset)
          call hm_set_dihe_list(3, idihe, gromol%dihes(k)%atom_idx3 + ioffset)
          call hm_set_dihe_list(4, idihe, gromol%dihes(k)%atom_idx4 + ioffset)

          idihe = idihe + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_impr
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_impr(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: iimpr, ioffset, i, j, k
    integer                  :: nimpr, nmole

    type(s_grotop_mol),      pointer     :: gromol


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Impr> '

    ! count number of impropers
    nimpr = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func == 2) then
            nimpr = nimpr + 1
          end if
        end do
      end do

    end do

    call hm_alloc_imprs(nimpr)

    ioffset = 0
    iimpr   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 2) &
            cycle

          call hm_set_impr_list(1, iimpr, gromol%dihes(k)%atom_idx1 + ioffset)
          call hm_set_impr_list(2, iimpr, gromol%dihes(k)%atom_idx2 + ioffset)
          call hm_set_impr_list(3, iimpr, gromol%dihes(k)%atom_idx3 + ioffset)
          call hm_set_impr_list(4, iimpr, gromol%dihes(k)%atom_idx4 + ioffset)

          iimpr = iimpr + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_grocrd_coord
  !> @brief        create huge molecule data
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_grocrd_coord(crd_in, is_ref)

    ! formal arguments
    integer,                 intent(in) :: crd_in
    logical,                 intent(in) :: is_ref

    ! local variables
    real(wp)                 :: vel(3), crd(3), center(3), atom_coord(3)
    integer                  :: iatom, natom, atom_no, resi_no
    character(200)           :: line, name
    character(8)             :: resi_name, atom_name


    if (.not. is_ref) then
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Grocrd_Coord> '
    else
      if (main_rank) &
        write(MsgOut,'(a)') 'Hm_Create_Grocrd_Coord (REF)> '
    end if


    ! compute center of system
    read(crd_in,'(a)') line
    read(crd_in,*) natom

    center(1:3) = 0.0_wp
    do iatom = 1, natom
      read(crd_in,'(i5,a4,a6,i5,3f8.3,3f8.3)') &
           resi_no,   &
           resi_name, &
           atom_name, &
           atom_no,   &
           crd(1:3),  &
           vel(1:3)
      center(1:3) = center(1:3) + crd(1:3) * 10.0_wp
    end do
    center(1:3) = center(1:3) / real(natom,wp)

    ! read crd file
    rewind(crd_in)
    read(crd_in,'(a)') line
    read(crd_in,*) natom
    
    do iatom = 1, natom

      read(crd_in,'(i5,a4,a6,i5,3f8.3,3f8.3)') &
           resi_no,   &
           resi_name, &
           atom_name, &
           atom_no,   &
           atom_coord(1:3), &
           vel(1:3)
           
!     atom_coord(1:3) = atom_coord(1:3) * 10.0_wp - center(1:3)
      atom_coord(1:3) = atom_coord(1:3) * 10.0_wp

      if (.not. is_ref) then
        call hm_set_atom_coord   (1, iatom, real(atom_coord(1),wip))
        call hm_set_atom_coord   (2, iatom, real(atom_coord(2),wip))
        call hm_set_atom_coord   (3, iatom, real(atom_coord(3),wip))
      else
        call hm_set_atom_refcoord(1, iatom, real(atom_coord(1),wip))
        call hm_set_atom_refcoord(2, iatom, real(atom_coord(2),wip))
        call hm_set_atom_refcoord(3, iatom, real(atom_coord(3),wip))
      end if

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Hm_Create_Grocrd_Coord> read ERROR. ')

  end subroutine hm_create_grocrd_coord

end module pr_gromacs2hm_mod

