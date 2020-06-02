!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_gromacs2cache_mod
!> @brief   create cached molecule from GROMACS data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_gromacs2cache_mod

  use pr_pdb2cache_mod
  use pr_cached_molecule_mod
  use fileio_grotop_mod
  use fileio_grocrd_mod
  use fileio_pdb_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: cm_create_gromacs
  private :: cm_create_cache
  private :: cm_create_cache_atom
  private :: cm_create_cache_bond
  private :: cm_create_cache_angl
  private :: cm_create_cache_dihe
  private :: cm_create_cache_impr
  private :: cm_create_cache_grocrd_coord
  private :: cm_create_cache_others

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_gromacs
  !> @brief        create cached molecule from GROMACS structure/coordinate file
  !! @authors      NT
  !! @param[in]    grotop     : GROMACS topology information
  !! @param[in]    grocrdfile : GROMACS coordiante file
  !! @param[in]    groreffile : Reference GROMACS coordinate file
  !! @param[in]    pdbfile    : PDB file name (optional)
  !! @param[in]    reffile    : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_gromacs(grotop, grocrdfile, groreffile, pdbfile, reffile)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    character(*),            intent(in)    :: grocrdfile
    character(*),            intent(in)    :: groreffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: grocrd_in, groref_in, pdb_in, ref_in


    if (cm_create_molecules()) then

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


      call cm_create_cache(grotop, grocrd_in, groref_in, pdb_in, ref_in)


      if (ref_in /= 0) &
        call close_file(ref_in)

      if (pdb_in /= 0) &
        call close_file(pdb_in)

      if (groref_in /= 0) &
        call close_file(groref_in)

      if (grocrd_in /= 0) &
        call close_file(grocrd_in)

      call cm_regist_molecules

    end if

    write(MsgOut,'(a)')    'Cm_Create_Gromacs> Molecule information: '
    write(MsgOut,'(a,i9)') '  num_atoms       = ', cm_num_atoms
    write(MsgOut,'(a,i9)') '  num_bonds       = ', cm_num_bonds
    write(MsgOut,'(a,i9)') '  num_angles      = ', cm_num_angles
    write(MsgOut,'(a,i9)') '  num_dihedrals   = ', cm_num_dihedrals
    write(MsgOut,'(a,i9)') '  num_impropers   = ', cm_num_impropers
    write(MsgOut,'(a)')    ' '

    return

  end subroutine cm_create_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache
  !> @brief        create cached files
  !! @authors      NT
  !! @param[in]    grotop    : GROMACS topology information
  !! @param[in]    grocrd_in : GROCRD file unit number
  !! @param[in]    groref_in : GROREF file unit number
  !! @param[in]    pdb_in    : PDB file unit number
  !! @param[in]    ref_in    : REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache(grotop, grocrd_in, groref_in, pdb_in, ref_in)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    integer,                 intent(in)    :: grocrd_in
    integer,                 intent(in)    :: groref_in
    integer,                 intent(in)    :: pdb_in
    integer,                 intent(in)    :: ref_in


    ! create cache from grotop information
    !

    call cm_create_cache_atom(grotop)

    call cm_create_cache_bond(grotop)

    call cm_create_cache_angl(grotop)

    call cm_create_cache_dihe(grotop)

    call cm_create_cache_impr(grotop)


    ! read grocrd / pdb
    !

    if (grocrd_in /= 0) then

      call cm_create_cache_grocrd_coord(grocrd_in, .false.)

    else if (pdb_in /= 0) then

      call cm_create_cache_pdb_coord(pdb_in, .false.)

    end if


    ! read reference pdb
    !

    if (groref_in /= 0) then

      call cm_create_cache_grocrd_coord(grocrd_in, .true.)

    else if (ref_in /= 0) then

      call cm_create_cache_pdb_coord(ref_in, .true.)

    end if


    ! create other cache files (empties)
    !

    call cm_create_cache_others


    write(MsgOut,'(a)') ''


    return

  end subroutine cm_create_cache

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_atom
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_atom(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idata, iatom, i, j, k, l, cls_no
    integer                  :: natom, nmole, num, file_out
    character(6)             :: cls_name, name

    type(s_grotop_mol),      pointer     :: gromol
    integer,                 allocatable :: atom_no      (:)
    integer,                 allocatable :: atom_cls_no  (:)
    character(6),            allocatable :: atom_cls_name(:)
    character(4),            allocatable :: atom_name    (:)
    real(wp),                allocatable :: charge       (:)
    real(wp),                allocatable :: mass         (:)
    integer,                 allocatable :: residue_no   (:)
    character(6),            allocatable :: residue_name (:)
    character(4),            allocatable :: segment_name (:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Atom> '

    ! allocate working memory
    allocate( &
         atom_no       (CMPageSize), &
         atom_cls_no   (CMPageSize), &
         atom_cls_name (CMPageSize), &
         atom_name     (CMPageSize), &
         charge        (CMPageSize), &
         mass          (CMPageSize), &
         residue_no    (CMPageSize), &
         residue_name  (CMPageSize), &
         segment_name  (CMPageSize))


    ! count number of atoms
    natom = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        natom = natom + gromol%num_atoms
      end do
    end do

    cm_num_atoms       = natom


    ! create cache files
    idata = 1
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

          atom_no      (idata) = iatom
          atom_cls_no  (idata) = cls_no
          atom_cls_name(idata) = cls_name
          atom_name    (idata) = gromol%atoms(k)%atom_name
          charge       (idata) = gromol%atoms(k)%charge
          mass         (idata) = gromol%atoms(k)%mass
          residue_no   (idata) = gromol%atoms(k)%residue_no
          residue_name (idata) = gromol%atoms(k)%residue_name
          segment_name (idata) = name

          idata = idata + 1
          
          if (idata > CMPageSize .or. iatom + 1 > natom) then
          
            if (mod(iatom, CMPageSize) == 0) then
              num = iatom - CMPageSize + 1
            else
              num = natom - mod(iatom, CMPageSize) + 1
            end if
          
            ! atom_no
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'atom_no'), &
                                  IOFileOutputNew)
            write(file_out) atom_no(1:CMPageSize)
            call close_file(file_out)
          
            ! atom_cls_no
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'atom_cls_no'), &
                                  IOFileOutputNew)
            write(file_out) atom_cls_no(1:CMPageSize)
            call close_file(file_out)
          
            ! atom_cls_name
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'atom_cls_name'), &
                                  IOFileOutputNew)
            write(file_out) atom_cls_name(1:CMPageSize)
            call close_file(file_out)
          
            ! atom_name
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'atom_name'), &
                                  IOFileOutputNew)
            write(file_out) atom_name(1:CMPageSize)
            call close_file(file_out)
          
            ! charge
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'charge'), &
                                  IOFileOutputNew)
            write(file_out) charge(1:CMPageSize)
            call close_file(file_out)
          
            ! mass
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'mass'), &
                                  IOFileOutputNew)
            write(file_out) mass(1:CMPageSize)
            call close_file(file_out)
          
            ! residue_no
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'residue_no'), &
                                  IOFileOutputNew)
            write(file_out) residue_no(1:CMPageSize)
            call close_file(file_out)
          
            ! residue_name
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'residue_name'), &
                                  IOFileOutputNew)
            write(file_out) residue_name(1:CMPageSize)
            call close_file(file_out)
          
            ! segment_name
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'segment_name'), &
                                  IOFileOutputNew)
            write(file_out) segment_name(1:CMPageSize)
            call close_file(file_out)
          
            idata = 1
          
          end if

          iatom = iatom + 1

        end do

      end do

    end do

    ! deallocate working memory
    deallocate(         &
         atom_no,       &
         atom_cls_no,   &
         atom_cls_name, &
         atom_name,     &
         charge,        &
         mass,          &
         residue_no,    &
         residue_name,  &
         segment_name)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_bond
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_bond(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idata, ibond, ioffset, i, j, k
    integer                  :: nbond, nmole, num, file_out

    type(s_grotop_mol),      pointer     :: gromol
    integer,                 allocatable :: bond_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Bond> '

    ! allocate working memory
    allocate(bond_list(2, CMPageSize))

    ! count number of bonds
    nbond = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        nbond = nbond + gromol%num_bonds
      end do
    end do

    cm_num_bonds = nbond

    ! create cache files
    ioffset = 0
    idata   = 1
    ibond   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_bonds

          bond_list(1,idata) = gromol%bonds(k)%atom_idx1 + ioffset
          bond_list(2,idata) = gromol%bonds(k)%atom_idx2 + ioffset

          idata = idata + 1

          if (idata > CMPageSize .or. ibond + 1 > nbond) then
          
            if (mod(ibond, CMPageSize) == 0) then
              num = ibond - CMPageSize + 1
            else
              num = nbond - mod(ibond, CMPageSize) + 1
            end if
          
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'bond_list'), &
                                  IOFileOutputNew)
            write(file_out) bond_list(1:2,1:CMPageSize)
            call close_file(file_out)

            idata = 1

          end if

          ibond = ibond + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    ! deallocate working memory
    deallocate(bond_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_angl
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_angl(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idata, iangl, ioffset, i, j, k
    integer                  :: nangl, nmole, num, file_out

    type(s_grotop_mol),      pointer     :: gromol
    integer,                 allocatable :: angl_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Angl> '

    ! allocate working memory
    allocate(angl_list(3, CMAnglSize))

    ! count number of bonds
    nangl = 0
    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole
        nangl = nangl + gromol%num_angls
      end do
    end do

    cm_num_angles = nangl

    ! create cache files
    ioffset = 0
    idata   = 1
    iangl   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_angls

          angl_list(1,idata) = gromol%angls(k)%atom_idx1 + ioffset
          angl_list(2,idata) = gromol%angls(k)%atom_idx2 + ioffset
          angl_list(3,idata) = gromol%angls(k)%atom_idx3 + ioffset

          idata = idata + 1

          if (idata > CMAnglSize .or. iangl + 1 > nangl) then
          
            if (mod(iangl, CMAnglSize) == 0) then
              num = iangl - CMAnglSize + 1
            else
              num = nangl - mod(iangl, CMAnglSize) + 1
            end if
          
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'angl_list'), &
                                  IOFileOutputNew)
            write(file_out) angl_list(1:3,1:CMAnglSize)
            call close_file(file_out)

            idata = 1

          end if

          iangl = iangl + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    ! deallocate working memory
    deallocate(angl_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_dihe
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_dihe(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idata, idihe, ioffset, i, j, k
    integer                  :: ndihe, nmole, num, file_out

    type(s_grotop_mol),      pointer     :: gromol
    integer,                 allocatable :: dihe_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Dihe> '

    ! allocate working memory
    allocate(dihe_list(4, CMPageSize))

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

    cm_num_dihedrals = ndihe

    ! create cache files
    ioffset = 0
    idata   = 1
    idihe   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 3) &
            cycle

          dihe_list(1,idata) = gromol%dihes(k)%atom_idx1 + ioffset
          dihe_list(2,idata) = gromol%dihes(k)%atom_idx2 + ioffset
          dihe_list(3,idata) = gromol%dihes(k)%atom_idx3 + ioffset
          dihe_list(4,idata) = gromol%dihes(k)%atom_idx4 + ioffset

          idata = idata + 1

          if (idata > CMPageSize .or. idihe + 1 > ndihe) then
          
            if (mod(idihe, CMPageSize) == 0) then
              num = idihe - CMPageSize + 1
            else
              num = ndihe - mod(idihe, CMPageSize) + 1
            end if
          
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'dihe_list'), &
                                  IOFileOutputNew)
            write(file_out) dihe_list(1:4,1:CMPageSize)
            call close_file(file_out)

            idata = 1

          end if

          idihe = idihe + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    ! deallocate working memory
    deallocate(dihe_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_impr
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_impr(grotop)

    ! formal arguments
    type(s_grotop),          intent(in) :: grotop

    ! local variables
    integer                  :: idata, iimpr, ioffset, i, j, k
    integer                  :: nimpr, nmole, num, file_out

    type(s_grotop_mol),      pointer     :: gromol
    integer,                 allocatable :: impr_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Impr> '

    ! allocate working memory
    allocate(impr_list(4, CMPageSize))

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

    cm_num_impropers = nimpr

    ! create cache files
    ioffset = 0
    idata   = 1
    iimpr   = 1

    do i = 1, grotop%num_molss

      gromol => grotop%molss(i)%moltype%mol
      nmole  =  grotop%molss(i)%count

      do j = 1, nmole

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 2) &
            cycle

          impr_list(1,idata) = gromol%dihes(k)%atom_idx1 + ioffset
          impr_list(2,idata) = gromol%dihes(k)%atom_idx2 + ioffset
          impr_list(3,idata) = gromol%dihes(k)%atom_idx3 + ioffset
          impr_list(4,idata) = gromol%dihes(k)%atom_idx4 + ioffset

          idata = idata + 1

          if (idata > CMPageSize .or. iimpr + 1 > nimpr) then
          
            if (mod(iimpr, CMPageSize) == 0) then
              num = iimpr - CMPageSize + 1
            else
              num = nimpr - mod(iimpr, CMPageSize) + 1
            end if
          
            call open_binary_file(file_out, &
                                  cm_get_cache_name(num, 'impr_list'), &
                                  IOFileOutputNew)
            write(file_out) impr_list(1:4,1:CMPageSize)
            call close_file(file_out)

            idata = 1

          end if

          iimpr = iimpr + 1

        end do

        ioffset = ioffset + gromol%num_atoms

      end do

    end do

    ! deallocate working memory
    deallocate(impr_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_grocrd_coord
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_grocrd_coord(crd_in, is_ref)

    ! formal arguments
    integer,                 intent(in) :: crd_in
    logical,                 intent(in) :: is_ref

    ! local variables
    real(wp)                 :: vel(3), crd(3), center(3)
    integer                  :: iatom, idata, file_out, num
    integer                  :: natom, atom_no, resi_no
    character(200)           :: line, name
    character(8)             :: resi_name, atom_name

    real(wp),    allocatable :: atom_coord(:,:)


    if (.not. is_ref) then
      write(MsgOut,'(a)') 'Cm_Create_Cache_Grocrd_Coord> '
      name = 'atom_coord'
    else
      write(MsgOut,'(a)') 'Cm_Create_Cache_Grocrd_Coord (REF)> '
      name = 'atom_refcoord'
    end if


    allocate(atom_coord(3, CMPageSize))

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
    
    idata = 1

    do iatom = 1, natom

      read(crd_in,'(i5,a4,a6,i5,3f8.3,3f8.3)') &
           resi_no,   &
           resi_name, &
           atom_name, &
           atom_no,   &
           atom_coord(1:3,idata), &
           vel(1:3)
           
      atom_coord(1:3,idata) = atom_coord(1:3,idata) * 10.0_wp - center(1:3)

      idata = idata + 1

      if (idata > CMPageSize .or. iatom + 1 > cm_num_atoms) then

        if (idata - 1 == CMPageSize) then
          num = iatom - CMPageSize + 1
        else
          num = iatom / CMPageSize * CMPageSize + 1
        end if

        ! atom_coord
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, name), &
                              IOFileOutputNew)
        write(file_out) atom_coord(1:3,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(atom_coord)

    write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Cm_Create_Cache_Grocrd_Coord> read ERROR. ')

  end subroutine cm_create_cache_grocrd_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_others
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_others

    ! local variables
    integer                  :: file_out, i, num

    real(wp),    allocatable :: r3_array(:,:)
    integer,     allocatable :: i1_array(:)

    logical                  :: ref_exist


    write(MsgOut,'(a)') 'Cm_Create_Cache_Others> '

    inquire(file=cm_get_cache_name(1, "atom_refcoord"), exist=ref_exist)

    allocate(r3_array(3, CMPageSize), i1_array(CMPageSize))
    r3_array(1:3,1:CMPageSize) = 0
    i1_array(    1:CMPageSize) = 0
    
    do i = 0, (cm_num_atoms / CMPageSize)

      num = i * CMPageSize + 1

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'atom_velocity'), &
                            IOFileOutputNew)
      write(file_out) r3_array(1:3,1:CMPageSize)
      call close_file(file_out)

      call open_binary_file(file_out, &
                            cm_get_cache_name(num, 'molecule_no'), &
                            IOFileOutputNew)
      write(file_out) i1_array(1:CMPageSize)
      call close_file(file_out)

      if (.not. ref_exist) then
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'atom_refcoord'), &
                              IOFileOutputNew)
        write(file_out) r3_array(1:3,1:CMPageSize)
        call close_file(file_out)
      end if

    end do

    deallocate(r3_array, i1_array)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_others

end module pr_gromacs2cache_mod

