!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_amber2cache_mod
!> @brief   create cached molecule from AMBER data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_amber2cache_mod

  use pr_pdb2cache_mod
  use pr_cached_molecule_mod
  use fileio_prmtop_mod
  use fileio_ambcrd_mod
  use fileio_pdb_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: cm_create_amber
  private :: cm_create_cache
  private :: cm_create_cache_atom
  private :: cm_create_cache_bond
  private :: cm_create_cache_angl
  private :: cm_create_cache_dihe
  private :: cm_create_cache_impr
  private :: cm_create_cache_ambcrd_coord
  private :: cm_create_cache_others

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_amber
  !> @brief        create cached molecule from AMBER parameter/topology file
  !! @authors      NT
  !! @param[in]    prmtop     : AMBER parameter/topology information
  !! @param[in]    ambcrdfile : AMBER coordiante file
  !! @param[in]    ambreffile : Reference AMBER coordinate file
  !! @param[in]    pdbfile    : PDB file name (optional)
  !! @param[in]    reffile    : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_amber(prmtop, ambcrdfile, ambreffile, pdbfile, reffile)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    character(*),            intent(in)    :: ambcrdfile
    character(*),            intent(in)    :: ambreffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: ambcrd_in, ambref_in, pdb_in, ref_in


    if (cm_create_molecules()) then

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


      call cm_create_cache(prmtop, ambcrd_in, ambref_in, pdb_in, ref_in)


      if (ref_in /= 0) &
        call close_file(ref_in)

      if (pdb_in /= 0) &
        call close_file(pdb_in)

      if (ambref_in /= 0) &
        call close_file(ambref_in)

      if (ambcrd_in /= 0) &
        call close_file(ambcrd_in)

      call cm_regist_molecules

    end if

    write(MsgOut,'(a)')    'Cm_Create_Amber> Molecule information: '
    write(MsgOut,'(a,i9)') '  num_atoms       = ', cm_num_atoms
    write(MsgOut,'(a,i9)') '  num_bonds       = ', cm_num_bonds
    write(MsgOut,'(a,i9)') '  num_angles      = ', cm_num_angles
    write(MsgOut,'(a,i9)') '  num_dihedrals   = ', cm_num_dihedrals
    write(MsgOut,'(a,i9)') '  num_impropers   = ', cm_num_impropers
    write(MsgOut,'(a)')    ' '

    return

  end subroutine cm_create_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache
  !> @brief        create cached files
  !! @authors      NT
  !! @param[in]    prmtop    : AMBER parameter/topology information
  !! @param[in]    ambcrd_in : AMBCRD file unit number
  !! @param[in]    ambref_in : AMBREF file unit number
  !! @param[in]    pdb_in    : PDB file unit number
  !! @param[in]    ref_in    : REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache(prmtop, ambcrd_in, ambref_in, pdb_in, ref_in)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    integer,                 intent(in)    :: ambcrd_in
    integer,                 intent(in)    :: ambref_in
    integer,                 intent(in)    :: pdb_in
    integer,                 intent(in)    :: ref_in


    ! create cache from grotop information
    !

    call cm_create_cache_atom(prmtop)

    call cm_create_cache_bond(prmtop)

    call cm_create_cache_angl(prmtop)

    call cm_create_cache_dihe(prmtop)

    call cm_create_cache_impr(prmtop)


    ! read ambcrd / pdb
    !

    if (ambcrd_in /= 0) then

      call cm_create_cache_ambcrd_coord(ambcrd_in, .false.)

    else if (pdb_in /= 0) then

      call cm_create_cache_pdb_coord(pdb_in, .false.)

    end if


    ! read reference pdb
    !

    if (ambref_in /= 0) then

      call cm_create_cache_ambcrd_coord(ambref_in, .true.)

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

  subroutine cm_create_cache_atom(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, iatom, i, is, ie
    integer                  :: natom, num, file_out

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

    cm_num_atoms       = prmtop%num_atoms

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

    natom = cm_num_atoms
    idata = 1

    do iatom = 1, natom

      atom_no      (idata) = iatom
      atom_cls_no  (idata) = prmtop%atom_cls_no  (iatom)
      atom_cls_name(idata) = prmtop%atom_cls_name(iatom)
      atom_name    (idata) = prmtop%atom_name    (iatom)
      charge       (idata) = prmtop%charge       (iatom)/18.2223_wp
      mass         (idata) = prmtop%mass         (iatom)
      segment_name (idata) = ''

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

        ! segment_name
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'segment_name'), &
                              IOFileOutputNew)
        write(file_out) segment_name(1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    idata = 1

    do i = 1, prmtop%num_residues

      is = prmtop%res_point(i)
      if (i < prmtop%num_residues) then
        ie = prmtop%res_point(i+1)-1
      else
        ie = natom
      end if

      do iatom = is, ie

        residue_no  (idata) = i
        residue_name(idata) = prmtop%res_label(i)

        idata = idata + 1

        if (idata > CMPageSize .or. iatom + 1 > natom) then

          if (mod(iatom, CMPageSize) == 0) then
            num = iatom - CMPageSize + 1
          else
            num = natom - mod(iatom, CMPageSize) + 1
          end if

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

          idata = 1

        end if
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

  subroutine cm_create_cache_bond(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, ibond, i
    integer                  :: nbond, num, file_out

    integer,                 allocatable :: bond_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Bond> '

    cm_num_bonds = prmtop%num_bondh + prmtop%num_mbonda

    allocate(bond_list(2, CMPageSize))

    idata = 1
    ibond = 1
    nbond = cm_num_bonds

    do i = 1, prmtop%num_bondh

      bond_list(1, idata) = prmtop%bond_inc_hy(1,i) / 3 + 1
      bond_list(2, idata) = prmtop%bond_inc_hy(2,i) / 3 + 1

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

    do i = 1, prmtop%num_mbonda

      bond_list(1, idata) = prmtop%bond_wo_hy(1,i) / 3 + 1
      bond_list(2, idata) = prmtop%bond_wo_hy(2,i) / 3 + 1

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

  subroutine cm_create_cache_angl(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, iangl, i
    integer                  :: nangl, num, file_out

    integer,                 allocatable :: angl_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Angl> '

    cm_num_angles = prmtop%num_anglh + prmtop%num_mangla

    ! allocate working memory
    allocate(angl_list(3, CMAnglSize))


    idata = 1
    iangl = 1
    nangl = cm_num_angles

    do i = 1, prmtop%num_anglh

      angl_list(1, idata) = prmtop%angl_inc_hy(1,i) / 3 + 1
      angl_list(2, idata) = prmtop%angl_inc_hy(2,i) / 3 + 1
      angl_list(3, idata) = prmtop%angl_inc_hy(3,i) / 3 + 1

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

    do i = 1, prmtop%num_mangla

      angl_list(1, idata) = prmtop%angl_wo_hy(1,i) / 3 + 1
      angl_list(2, idata) = prmtop%angl_wo_hy(2,i) / 3 + 1
      angl_list(3, idata) = prmtop%angl_wo_hy(3,i) / 3 + 1

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

  subroutine cm_create_cache_dihe(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, idihe, i
    integer                  :: ndihe, num, file_out

    integer,                 allocatable :: dihe_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Dihe> '

    cm_num_dihedrals = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) >= 0) &
        cm_num_dihedrals = cm_num_dihedrals + 1
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) >= 0) &
        cm_num_dihedrals = cm_num_dihedrals + 1
    end do

    ! allocate working memory
    allocate(dihe_list(4, CMPageSize))

    idata = 1
    idihe = 1
    ndihe = cm_num_dihedrals

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cycle

      dihe_list(1, idata) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
      dihe_list(2, idata) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
      dihe_list(3, idata) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
      dihe_list(4, idata) =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

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

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cycle

      dihe_list(1, idata) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
      dihe_list(2, idata) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
      dihe_list(3, idata) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
      dihe_list(4, idata) =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

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

  subroutine cm_create_cache_impr(prmtop)

    ! formal arguments
    type(s_prmtop),          intent(in) :: prmtop

    ! local variables
    integer                  :: idata, iimpr, i
    integer                  :: nimpr, num, file_out

    integer,                 allocatable :: impr_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Impr> '

    cm_num_impropers = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) &
        cm_num_impropers = cm_num_impropers + 1
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) &
        cm_num_impropers = cm_num_impropers + 1
    end do

    ! allocate working memory
    allocate(impr_list(4, CMPageSize))

    idata = 1
    iimpr = 1
    nimpr = cm_num_impropers

    do i = 1, prmtop%num_diheh

      if (prmtop%dihe_inc_hy(4,i) >= 0) &
        cycle

      impr_list(1, idata) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
      impr_list(2, idata) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
      impr_list(3, idata) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
      impr_list(4, idata) = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

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

    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) >= 0) &
        cycle

      impr_list(1, idata) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
      impr_list(2, idata) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
      impr_list(3, idata) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
      impr_list(4, idata) = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

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

    ! deallocate working memory
    deallocate(impr_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_ambcrd_coord
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_ambcrd_coord(crd_in, is_ref)

    ! formal arguments
    integer,                 intent(in)    :: crd_in
    logical,                 intent(in)    :: is_ref

    ! local variables
    real(wp)                 :: center(3), crd(3,2)
    integer                  :: natom, iatom, idata, i, nr
    integer                  :: file_out, num
    character(100)           :: line, name

    real(wp),    allocatable :: atom_coord(:,:)


    if (.not. is_ref) then
      write(MsgOut,'(a)') 'Cm_Create_Cache_Ambcrd_Coord> '
      name = 'atom_coord'
    else
      write(MsgOut,'(a)') 'Cm_Create_Cache_Ambcrd_Coord (REF)> '
      name = 'atom_refcoord'
    end if


    allocate(atom_coord(3, CMPageSize))

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
    center(1:3) = center(1:3) / real(natom,wp)

    ! read crd file
    rewind(crd_in)

    read(crd_in,'(a)') line
    read(crd_in,*) natom

    idata = 1

    do iatom = 1, natom, 2

      if (natom - iatom < 1) then
        nr = 0
      else
        nr = 1
      end if

      read(crd_in,*) (atom_coord(1, idata + i), &
                      atom_coord(2, idata + i), &
                      atom_coord(3, idata + i), i = 0, nr)

      atom_coord(1:3,idata)    = atom_coord(1:3,idata)    - center(1:3)
      if (nr > 0) &
      atom_coord(1:3,idata+nr) = atom_coord(1:3,idata+nr) - center(1:3)

      idata = idata + 2

      if (idata > CMPageSize .or. iatom + 2 > natom) then

        if (idata - 1 == CMPageSize) then
          num = iatom - CMPageSize + 2
        else
          num = iatom / CMPageSize * CMPageSize + 1
        end if

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

900 call error_msg('Cm_Create_Cache_Ambcrd_Coord> read ERROR. ')

  end subroutine cm_create_cache_ambcrd_coord

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

end module pr_amber2cache_mod
