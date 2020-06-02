!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_charmm2cache_mod
!> @brief   create cached molecule from CHARMM data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_charmm2cache_mod

  use pr_pdb2cache_mod
  use pr_cached_molecule_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_psf_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: cm_create_charmm
  private :: cm_create_cache
  private :: cm_read_psf_header
  private :: cm_create_cache_psf_atom
  private :: cm_create_cache_psf_bond
  private :: cm_create_cache_psf_angl
  private :: cm_create_cache_psf_dihe
  private :: cm_create_cache_psf_impr
  private :: cm_create_cache_psf_cmap
  private :: cm_create_cache_crd_coord
  private :: cm_create_cache_others

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_charmm
  !> @brief        create cached molecule from CHARMM structure/coordinate file
  !! @authors      NT
  !! @param[in]    top     : CHARMM top information
  !! @param[in]    par     : CHARMM par information
  !! @param[in]    psffile : PSF file name
  !! @param[in]    pdbfile : PDB file name (optional)
  !! @param[in]    crdfile : CHARMM CRD file name (optional)
  !! @param[in]    reffile : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_charmm(top, par, psffile, pdbfile, crdfile, reffile)

    ! formal arguments
    type(s_top),             intent(in)    :: top
    type(s_par),             intent(in)    :: par
    character(*),            intent(in)    :: psffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: crdfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: psf_in, pdb_in, crd_in, ref_in


    if (cm_create_molecules()) then

      psf_in = 0
      pdb_in = 0
      crd_in = 0
      ref_in = 0

      call open_file(psf_in, psffile, IOFileInput)

      if (pdbfile /= '') &
        call open_file(pdb_in, pdbfile, IOFileInput)

      if (crdfile /= '') &
        call open_file(crd_in, crdfile, IOFileInput)

      if (reffile /= '') &
        call open_file(ref_in, reffile, IOFileInput)


      call cm_create_cache(top, par, psf_in, pdb_in, crd_in, ref_in)


      if (ref_in /= 0) &
        call close_file(ref_in)

      if (crd_in /= 0) &
        call close_file(crd_in)

      if (pdb_in /= 0) &
        call close_file(pdb_in)

      call close_file(psf_in)


      call cm_regist_molecules

    end if

    write(MsgOut,'(a)')    'Cm_Create_Charmm> Molecule information: '
    write(MsgOut,'(a,i9)') '  num_atoms       = ', cm_num_atoms
    write(MsgOut,'(a,i9)') '  num_bonds       = ', cm_num_bonds
    write(MsgOut,'(a,i9)') '  num_angles      = ', cm_num_angles
    write(MsgOut,'(a,i9)') '  num_dihedrals   = ', cm_num_dihedrals
    write(MsgOut,'(a,i9)') '  num_impropers   = ', cm_num_impropers
    write(MsgOut,'(a,i9)') '  num_cmaps       = ', cm_num_cmaps
    write(MsgOut,'(a)')    ' '

    return

  end subroutine cm_create_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache
  !> @brief        create cached files
  !! @authors      NT
  !! @param[in]    top     CHARMM top information
  !! @param[in]    par     CHARMM par information
  !! @param[in]    psf_in  PSF file unit number
  !! @param[in]    pdb_in  PDB file unit number
  !! @param[in]    crd_in  CHARMM CRD file unit number
  !! @param[in]    ref_in  REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache(top, par, psf_in, pdb_in, crd_in, ref_in)

    ! formal arguments
    type(s_top),             intent(in)    :: top
    type(s_par),             intent(in)    :: par
    integer,                 intent(in)    :: psf_in
    integer,                 intent(in)    :: pdb_in
    integer,                 intent(in)    :: crd_in
    integer,                 intent(in)    :: ref_in

    ! local variables
    integer                  :: i, n1, n2, tmp
    character(100)           :: key
    logical                  :: can_read


    ! read psf
    !

    call cm_read_psf_header(psf_in)

    do while(.true.)

       can_read = .false.
       read(psf_in,*,err=10,end=10) n1, n2, key
       can_read = .true.

10     if (.not. can_read) then
          backspace(psf_in)
          read(psf_in, *, err=100,end=100) n1, key
       end if

       if (key(1:1) == '!') &
         key = key(2:)

       if (key(1:5) == 'NATOM') then
         call cm_create_cache_psf_atom(psf_in, n1, par, top)

       else if (key(1:5) == 'NBOND') then
         call cm_create_cache_psf_bond(psf_in, n1)
 
       else if (key(1:6) == 'NTHETA') then
         call cm_create_cache_psf_angl(psf_in, n1)
 
       else if (key(1:4) == 'NPHI') then
         call cm_create_cache_psf_dihe(psf_in, n1)
 
       else if (key(1:6) =='NIMPHI') then
         call cm_create_cache_psf_impr(psf_in, n1)
 
       else if (key(1:4) == 'NDON') then
         read(psf_in,*) (tmp, tmp, i = 1, n1)
 
       else if (key(1:4) == 'NACC') then
         read(psf_in,*) (tmp, tmp, i = 1, n1)
 
       else if (key(1:3) == 'NNB') then
         read(psf_in,*) (tmp, i = 1, cm_num_atoms)

       else if (key(1:4) == 'NGRP') then
         read(psf_in,*) (tmp, tmp, tmp, i = 1, n1)
 
       else if (key(1:7) == 'NCRTERM') then
         call cm_create_cache_psf_cmap(psf_in, n1)

       end if

    end do
100 continue    


    ! read pdb / crd
    !

    if (pdb_in /= 0) then

      call cm_create_cache_pdb_coord(pdb_in, .false.)

    else if (crd_in /= 0) then

      call cm_create_cache_crd_coord(crd_in)

    end if


    ! read reference pdb
    !

    if (ref_in /= 0) then

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
  !  Subroutine    cm_read_psf_header
  !! @authors      NT
  !! @param[in]    file : unit number of PSF file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_read_psf_header(file)

    ! parameter
    integer,                 parameter     :: MAXROW = 80

    ! formal arguments
    integer,                 intent(in)    :: file

    ! local variables
    integer                  :: i, ntitle
    character(100)           :: line, title
    logical                  :: ipsf, iext, icmap, icheq


    read(file,'(a80)') line

    ipsf  = .false.
    iext  = .false.
    icmap = .false.
    icheq = .false.

    do i = 1, MAXROW
      if (line(i:i+3) == 'PSF ') ipsf  = .true.
      if (line(i:i+3) == 'EXT ') iext  = .true.
      if (line(i:i+3) == 'CMAP') icmap = .true.
      if (line(i:i+3) == 'CHEQ') icheq = .true.
    end do

    if (.not. ipsf) then
      call error_msg('Cm_Read_Psf_Header> ERROR: Format is not correct')
    end if

    read(file,'(a80)') line
    read(file, *) ntitle
    do i = 1, ntitle
      read(file, '(a80)') title
    end do

    return

  end subroutine cm_read_psf_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_atom
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_atom(psf_in, natom, par, top)

    ! parameters
    integer,   parameter    :: PSFTypeCHARMM = 1
    integer,   parameter    :: PSFTypeXPLOR  = 2

    ! formal arguments
    integer,                intent(in)    :: psf_in
    integer,                intent(in)    :: natom
    type(s_par),            intent(in)    :: par
    type(s_top),            intent(in)    :: top

    ! local variables
    integer                 :: iatom, idata, j, icls, psf_type, imove
    integer                 :: num, file_out
    character(8)            :: cstr, cseg_nm, cres_nm, catm_nm
    character(6)            :: ccls
    character(20)           :: ctmp

    integer,                allocatable :: atom_no      (:)
    integer,                allocatable :: atom_cls_no  (:)
    character(6),           allocatable :: atom_cls_name(:)
    character(4),           allocatable :: atom_name    (:)
    real(wp),               allocatable :: charge       (:)
    real(wp),               allocatable :: mass         (:)
    integer,                allocatable :: residue_no   (:)
    character(6),           allocatable :: residue_name (:)
    character(4),           allocatable :: segment_name (:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Atom> '

    cm_num_atoms       = natom

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

    ! check psf file type
    read(psf_in,*) ctmp, ctmp, ctmp, ctmp, ctmp, cstr
    if (cstr(1:1) >= 'A' .and. cstr(1:1) <= 'Z' .or. &
        cstr(1:1) >= 'a' .and. cstr(1:1) <= 'z') then
      psf_type = PsfTypeXPLOR
    else
      psf_type = PsfTypeCHARMM
    end if

    backspace(psf_in)

    ! read psf file
    idata = 1

    do iatom = 1, natom

      select case(psf_type)

      case (PsfTypeCHARMM)

        read(psf_in,*)           &
             atom_no    (idata), &
             cseg_nm,            &
             residue_no (idata), &
             cres_nm,            &
             catm_nm,            &
             atom_cls_no(idata), &
             charge     (idata), &
             mass       (idata), &
             imove

        segment_name(idata) = cseg_nm(1:4)
        residue_name(idata) = cres_nm(1:6)
        atom_name   (idata) = catm_nm(1:4)

        icls = atom_cls_no(idata)
        do j = 1, top%num_atom_cls
          if (icls == top%atom_cls_no(j)) then
            atom_cls_name(idata) = top%atom_cls_name(j)
            exit
          end if
        end do
        if (j > top%num_atom_cls) then
          write(MsgOut,*) '  WARNING: atom class number is not found. ', icls
        end if

        ccls = atom_cls_name(idata)
        do j = 1, par%num_atom_cls
          if (ccls == par%nonb_atom_cls(j)) then
            atom_cls_no(idata) = j
            exit
          end if
        end do
        if (j == par%num_atom_cls + 1) then
          write(MsgOut,*) '  WARNING: atom class name is undefined. '//trim(ccls)
        end if

      case (PsfTypeXPLOR)

        read(psf_in, *)            &
             atom_no      (idata), &
             segment_name (idata), &
             residue_no   (idata), &
             residue_name (idata), &
             atom_name    (idata), &
             atom_cls_name(idata), &
             charge       (idata), &
             mass         (idata), &
             imove

        ccls = atom_cls_name(idata)
        do j = 1, par%num_atom_cls
          if (ccls == par%nonb_atom_cls(j)) then
            atom_cls_no(idata) = j
            exit
          end if
        end do
        if (j == par%num_atom_cls + 1) then
          write(MsgOut,*) '  WARNING: atom class name is undefined. '//trim(ccls)
        end if

      end select

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

    end do

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

  end subroutine cm_create_cache_psf_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_bond
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_bond(psf_in, nbond)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nbond

    ! local variables
    integer                  :: ibond, idata, i, nr
    integer                  :: file_out, num

    integer,     allocatable :: bond_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Bond> '

    cm_num_bonds = nbond

    allocate(bond_list(2, CMPageSize))

    ! read psf file
    idata = 1

    do ibond = 1, nbond, 4

      if (nbond - ibond < 4) then
        nr = nbond - ibond
      else
        nr = 3
      end if
      read(psf_in,*) (bond_list(1, idata + i), &
                      bond_list(2, idata + i), i = 0,nr)

      idata = idata + 4

      if (idata > CMPageSize .or. ibond + 4 > nbond) then

        if (idata - 1 == CMPageSize) then
          num = ibond - CMPageSize + 4
        else
          num = ibond / CMPageSize * CMPageSize + 1
        end if

        ! bond_list
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'bond_list'), &
                              IOFileOutputNew)
        write(file_out) bond_list(1:2,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(bond_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_psf_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_angl
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_angl(psf_in, nangl)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nangl

    ! local variables
    integer                  :: iangl, idata, i, nr
    integer                  :: file_out, num

    integer,    allocatable  :: angl_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Angl> '

    cm_num_angles = nangl

    allocate(angl_list(3, CMAnglSize))

    ! read psf file
    idata = 1

    do iangl = 1, nangl, 3

      if (nangl - iangl < 3) then
        nr = nangl - iangl
      else
        nr = 2
      end if
      read(psf_in,*) (angl_list(1, idata + i), &
                      angl_list(2, idata + i), &
                      angl_list(3, idata + i), i = 0,nr)

      idata = idata + 3

      if (idata > CMAnglSize .or. iangl + 3 > nangl) then

        if (idata - 1 == CMAnglSize) then
          num = iangl - CMAnglSize + 3
        else
          num = iangl / CMAnglSize * CMAnglSize + 1
        end if

        ! angl_list
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'angl_list'), &
                              IOFileOutputNew)
        write(file_out) angl_list(1:3,1:CMAnglSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(angl_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_psf_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_dihe
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_dihe(psf_in, ndihe)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: ndihe

    ! local variables
    integer                  :: idihe, idata, i, nr
    integer                  :: file_out, num

    integer,     allocatable :: dihe_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Dihe> '

    cm_num_dihedrals = ndihe

    allocate(dihe_list(4, CMPageSize))

    ! read psf file
    idata = 1

    do idihe = 1, ndihe, 2

      if (ndihe - idihe < 2) then
        nr = ndihe - idihe
      else
        nr = 1
      end if
      read(psf_in,*) (dihe_list(1, idata + i), &
                      dihe_list(2, idata + i), &
                      dihe_list(3, idata + i), &
                      dihe_list(4, idata + i), i = 0,nr)

      idata = idata + 2

      if (idata > CMPageSize .or. idihe + 2 > ndihe) then

        if (idata - 1 == CMPageSize) then
          num = idihe - CMPageSize + 2
        else
          num = idihe / CMPageSize * CMPageSize + 1
        end if

        ! dihe_list
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'dihe_list'), &
                              IOFileOutputNew)
        write(file_out) dihe_list(1:4,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(dihe_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_psf_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_impr
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_impr(psf_in, nimpr)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nimpr

    ! local variables
    integer                  :: iimpr, idata, i, nr
    integer                  :: file_out, num

    integer,     allocatable :: impr_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Impr> '

    cm_num_impropers = nimpr

    allocate(impr_list(4, CMPageSize))

    ! read psf file
    idata = 1

    do iimpr = 1, nimpr, 2

      if (nimpr - iimpr < 2) then
        nr = nimpr - iimpr
      else
        nr = 1
      end if
      read(psf_in,*) (impr_list(1, idata + i), &
                      impr_list(2, idata + i), &
                      impr_list(3, idata + i), &
                      impr_list(4, idata + i), i = 0,nr)

      idata = idata + 2

      if (idata > CMPageSize .or. iimpr + 2 > nimpr) then

        if (idata - 1 == CMPageSize) then
          num = iimpr - CMPageSize + 2
        else
          num = iimpr / CMPageSize * CMPageSize + 1
        end if

        ! impr_list
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'impr_list'), &
                              IOFileOutputNew)
        write(file_out) impr_list(1:4,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(impr_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_psf_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_psf_cmap
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_psf_cmap(psf_in, ncmap)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: ncmap

    ! local variables
    integer                  :: icmap, idata, i
    integer                  :: file_out, num

    integer,     allocatable :: cmap_list(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Psf_Cmap> '

    cm_num_cmaps = ncmap

    allocate(cmap_list(8, CMPageSize))

    ! read psf file
    idata = 1

    do icmap = 1, ncmap

      read(psf_in,*) cmap_list(1, idata), &
                     cmap_list(2, idata), &
                     cmap_list(3, idata), &
                     cmap_list(4, idata), &
                     cmap_list(5, idata), &
                     cmap_list(6, idata), &
                     cmap_list(7, idata), &
                     cmap_list(8, idata)

      idata = idata + 1

      if (idata > CMPageSize .or. icmap + 1 > ncmap) then

        if (idata - 1 == CMPageSize) then
          num = icmap - CMPageSize + 1
        else
          num = icmap / CMPageSize * CMPageSize + 1
        end if

        ! cmap_list
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'cmap_list'), &
                              IOFileOutputNew)
        write(file_out) cmap_list(1:8,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(cmap_list)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_psf_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cm_create_cache_crd_coord
  !> @brief        create cached files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cm_create_cache_crd_coord(crd_in)

    ! formal arguments
    integer,                 intent(in) :: crd_in

    ! local variables
    integer                  :: iatom, idata, file_out, num
    integer                  :: natom, atom_no, resi_no
    character(200)           :: line
    character(8)             :: flag, resi_name, atom_name, segm_name

    real(wp),    allocatable :: atom_coord(:,:)


    write(MsgOut,'(a)') 'Cm_Create_Cache_Crd_Coord> '


    allocate(atom_coord(3, CMPageSize))

    ! read pdb file
    idata = 1

100 read(crd_in,'(a80)',err=900) line
    if (line(1:1) == '*') then
       goto 100
    else
       read(line,'(i10,a5)') natom, flag
    end if

    do iatom = 1, natom

      read(crd_in,'(a120)') line

      if (natom >= 100000 .or. flag == '  EXT') then
        read(line,'(i10,i10,2x,a8,2x,a8,3f20.10,2x,a8)') &
                   atom_no,             &
                   resi_no,             &
                   resi_name,           &
                   atom_name,           &
                   atom_coord(1,idata), &
                   atom_coord(2,idata), &
                   atom_coord(3,idata), &
                   segm_name
      else
        read(line,'(i5,i5,1x,a4,1x,a4,3f10.5,1x,a4)') &
                   atom_no,             &
                   resi_no,             &
                   resi_name,           &
                   atom_name,           &
                   atom_coord(1,idata), &
                   atom_coord(2,idata), &
                   atom_coord(3,idata), &
                   segm_name
      end if

      idata = idata + 1

      if (idata > CMPageSize .or. iatom + 1 > cm_num_atoms) then

        if (idata - 1 == CMPageSize) then
          num = iatom - CMPageSize + 1
        else
          num = iatom / CMPageSize * CMPageSize + 1
        end if

        ! atom_coord
        call open_binary_file(file_out, &
                              cm_get_cache_name(num, 'atom_coord'), &
                              IOFileOutputNew)
        write(file_out) atom_coord(1:3,1:CMPageSize)
        call close_file(file_out)

        idata = 1

      end if

    end do

    deallocate(atom_coord)

    write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Cm_Create_Cache_Crd_Coord> read ERROR. ')

  end subroutine cm_create_cache_crd_coord

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


    write(MsgOut,'(a)') 'Cm_Create_Cache_Others> '


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

    end do

    deallocate(r3_array, i1_array)

    write(MsgOut,'(a)') '  done.'

    return

  end subroutine cm_create_cache_others

end module pr_charmm2cache_mod
