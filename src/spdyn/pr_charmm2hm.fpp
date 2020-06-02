!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_charmm2hm_mod
!> @brief   create huge molecule from CHARMM data
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_charmm2hm_mod

  use pr_pdb2hm_mod
  use pr_huge_molecule_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_psf_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use timers_mod

  implicit none
  private

  ! subroutines
  public  :: hm_create_charmm
  private :: hm_create
  private :: hm_read_psf_header
  private :: hm_create_psf_atom
  private :: hm_create_psf_bond
  private :: hm_create_psf_angl
  private :: hm_create_psf_dihe
  private :: hm_create_psf_impr
  private :: hm_create_psf_cmap
  private :: hm_create_crd_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_charmm
  !> @brief        create huge molecule from CHARMM structure/coordinate file
  !! @authors      NT
  !! @param[in]    top     : CHARMM top information
  !! @param[in]    par     : CHARMM par information
  !! @param[in]    psffile : PSF file name
  !! @param[in]    pdbfile : PDB file name (optional)
  !! @param[in]    crdfile : CHARMM CRD file name (optional)
  !! @param[in]    reffile : Reference coordinate PDB file name (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_charmm(top, par, psffile, pdbfile, crdfile, reffile)

    ! formal arguments
    type(s_top),             intent(in)    :: top
    type(s_par),             intent(in)    :: par
    character(*),            intent(in)    :: psffile
    character(*),            intent(in)    :: pdbfile
    character(*),            intent(in)    :: crdfile
    character(*),            intent(in)    :: reffile

    ! local variables
    integer                  :: psf_in, pdb_in, crd_in, ref_in


    psf_in = 0
    pdb_in = 0
    crd_in = 0
    ref_in = 0

    call open_file(psf_in, psffile, IOFileInput)

    if (pdbfile /= '') &
      call open_file(pdb_in, pdbfile, IOFileInput)

    if (crdfile /= '') &
      call open_file(crd_in, crdfile, IOFileInput)

    if (reffile /= '' .and. reffile /= pdbfile) &
      call open_file(ref_in, reffile, IOFileInput)

    if (reffile /= '' .and. reffile == pdbfile)  &
      ref_in = pdb_in

    call hm_create(top, par, psf_in, pdb_in, crd_in, ref_in)


    if (ref_in /= 0) &
      call close_file(ref_in)

    if (crd_in /= 0) &
      call close_file(crd_in)

    if (pdb_in /= 0) &
      call close_file(pdb_in)

    call close_file(psf_in)


    if (main_rank) then
      write(MsgOut,'(a)')    'Hm_Create_Charmm> Molecule information: '
      write(MsgOut,'(a,i9)') '  num_atoms       = ', hm_num_atoms
      write(MsgOut,'(a,i9)') '  num_bonds       = ', hm_num_bonds
      write(MsgOut,'(a,i9)') '  num_angles      = ', hm_num_angles
      write(MsgOut,'(a,i9)') '  num_dihedrals   = ', hm_num_dihedrals
      write(MsgOut,'(a,i9)') '  num_impropers   = ', hm_num_impropers
      write(MsgOut,'(a,i9)') '  num_cmaps       = ', hm_num_cmaps
      write(MsgOut,'(a)')    ' '
    end if

    return

  end subroutine hm_create_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create
  !> @brief        create huge molecule files
  !! @authors      NT
  !! @param[in]    top     CHARMM top information
  !! @param[in]    par     CHARMM par information
  !! @param[in]    psf_in  PSF file unit number
  !! @param[in]    pdb_in  PDB file unit number
  !! @param[in]    crd_in  CHARMM CRD file unit number
  !! @param[in]    ref_in  REFERENCE PDB file unit number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create(top, par, psf_in, pdb_in, crd_in, ref_in)

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

    call hm_read_psf_header(psf_in)

    do while(.true.)

      can_read = .false.
      read(psf_in,*,err=10,end=10) n1, n2, key
      can_read = .true.

10    if (.not. can_read) then
        backspace(psf_in)
        read(psf_in, *, err=100,end=100) n1, key
      end if

      if (key(1:1) == '!') &
           key = key(2:)

      if (key(1:5) == 'NATOM') then
        call timer(TimerTest1, TimerOn)
        call hm_create_psf_atom(psf_in, n1, par, top)
        call timer(TimerTest1, TimerOff)

      else if (key(1:5) == 'NBOND') then
        call timer(TimerTest2, TimerOn)
        call hm_create_psf_bond(psf_in, n1)
        call timer(TimerTest2, TimerOff)

      else if (key(1:6) == 'NTHETA') then
        call timer(TimerTest3, TimerOn)
        call hm_create_psf_angl(psf_in, n1)
        call timer(TimerTest3, TimerOff)

      else if (key(1:4) == 'NPHI') then
        call timer(TimerTest4, TimerOn)
        call hm_create_psf_dihe(psf_in, n1)
        call timer(TimerTest4, TimerOff)

      else if (key(1:6) =='NIMPHI') then
        call timer(TimerTest5, TimerOn)
        call hm_create_psf_impr(psf_in, n1)
        call timer(TimerTest5, TimerOff)

      else if (key(1:4) == 'NDON') then
        read(psf_in,*) (tmp, tmp, i = 1, n1)

      else if (key(1:4) == 'NACC') then
        read(psf_in,*) (tmp, tmp, i = 1, n1)

      else if (key(1:3) == 'NNB') then
        read(psf_in,*) (tmp, i = 1, hm_num_atoms)

      else if (key(1:4) == 'NGRP') then
        read(psf_in,*) (tmp, tmp, tmp, i = 1, n1)

      else if (key(1:7) == 'NCRTERM') then
        call hm_create_psf_cmap(psf_in, n1)

      end if

    end do
100 continue


    ! read pdb / crd
    !

    call timer(TimerTest6, TimerOn)
    if (pdb_in /= 0) then

      call hm_create_pdb_coord(pdb_in, .false.)

    else if (crd_in /= 0) then

      call hm_create_crd_coord(crd_in)

    end if
    call timer(TimerTest6, TimerOff)


    ! read reference pdb
    !

    call timer(TimerTest7, TimerOn)
    if (ref_in /= 0) then

      if (ref_in == pdb_in) then
        rewind(ref_in)
      endif
      call hm_create_pdb_coord(ref_in, .true.)

    end if
    call timer(TimerTest7, TimerOff)


    if (main_rank) &
      write(MsgOut,'(a)') ''


    return

  end subroutine hm_create

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_read_psf_header
  !! @authors      NT
  !! @param[in]    file : unit number of PSF file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_read_psf_header(file)

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
      call error_msg('Hm_Read_Psf_Header> ERROR: Format is not correct')
    end if

    read(file,'(a80)') line
    read(file, *) ntitle
    do i = 1, ntitle
      read(file, '(a80)') title
    end do

    return

  end subroutine hm_read_psf_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_atom
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_atom(psf_in, natom, par, top)

    ! parameters
    integer,   parameter    :: PSFTypeCHARMM = 1
    integer,   parameter    :: PSFTypeXPLOR  = 2

    ! formal arguments
    integer,                intent(in)    :: psf_in
    integer,                intent(in)    :: natom
    type(s_par),            intent(in)    :: par
    type(s_top),            intent(in)    :: top

    ! local variables
    integer                 :: iatom, j, icls, psf_type, imove
    character(8)            :: cstr, cseg_nm, cres_nm, catm_nm
    character(6)            :: ccls
    character(20)           :: ctmp

    integer                 :: atom_no
    integer                 :: atom_cls_no
    character(6)            :: atom_cls_name
    character(4)            :: atom_name
    real(wp)                :: charge
    real(wp)                :: mass
    integer                 :: residue_no
    character(6)            :: residue_name
    character(4)            :: segment_name


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Atom> '

    call hm_alloc_atoms(natom)

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
    do iatom = 1, natom

      select case(psf_type)

      case (PsfTypeCHARMM)

        read(psf_in,*)    &
             atom_no,     &
             cseg_nm,     &
             residue_no,  &
             cres_nm,     &
             catm_nm,     &
             atom_cls_no, &
             charge,      &
             mass,        &
             imove

        segment_name = cseg_nm(1:4)
        residue_name = cres_nm(1:6)
        atom_name    = catm_nm(1:4)

        icls = atom_cls_no
        do j = 1, top%num_atom_cls
          if (icls == top%atom_cls_no(j)) then
            atom_cls_name = top%atom_cls_name(j)
            exit
          end if
        end do
        if (j > top%num_atom_cls) then
          write(MsgOut,*) '  WARNING: atom class number is not found. ', icls
        end if

        ccls = atom_cls_name
        do j = 1, par%num_atom_cls
          if (ccls == par%nonb_atom_cls(j)) then
            atom_cls_no = j
            exit
          end if
        end do
        if (j == par%num_atom_cls + 1) then
          write(MsgOut,*)'  WARNING: atom class name is undefined. '//trim(ccls)
        end if

      case (PsfTypeXPLOR)

        read(psf_in, *)     &
             atom_no,       &
             segment_name,  &
             residue_no,    &
             residue_name,  &
             atom_name,     &
             atom_cls_name, &
             charge,        &
             mass,          &
             imove

        ccls = atom_cls_name
        do j = 1, par%num_atom_cls
          if (ccls == par%nonb_atom_cls(j)) then
            atom_cls_no = j
            exit
          end if
        end do
        if (j == par%num_atom_cls + 1) then
          write(MsgOut,*)'  WARNING: atom class name is undefined. '//trim(ccls)
        end if

      end select

      call hm_set_atom_no      (iatom, atom_no)
      call hm_set_atom_cls_no  (iatom, atom_cls_no)
      call hm_set_atom_cls_name(iatom, atom_cls_name)
      call hm_set_atom_name    (iatom, atom_name)
      call hm_set_charge       (iatom, charge)
      call hm_set_mass         (iatom, real(mass,wip))
      call hm_set_residue_no   (iatom, residue_no)
      call hm_set_residue_name (iatom, residue_name)
      call hm_set_segment_name (iatom, segment_name)

    end do

!   ! check light atom
!   !
!   do iatom = 1, natom
!     if (hm_mass(iatom) > LIGHT_ATOM_MASS_LIMIT) cycle
!     call hm_set_light_atom_mass(iatom, 1)
!     atom_name = hm_atom_name(iatom)
!     atom_cls_name = hm_atom_cls_name(iatom)
!     if (atom_name(1:1) == 'H' .or. atom_name(1:1) == 'D' .or. &
!         atom_name(1:1) == 'h' .or. atom_name(1:1) == 'd')     &
!       call hm_set_light_atom_name(iatom, 1)
!     if (atom_cls_name(1:1) == 'H' .or. atom_cls_name(1:1) == 'D' .or. &
!         atom_cls_name(1:1) == 'h' .or. atom_cls_name(1:1) == 'd')     &
!       call hm_set_light_atom_name(iatom, 1)
!   end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_bond
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_bond(psf_in, nbond)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nbond

    ! local variables
    integer                  :: ibond, i, nr
    integer                  :: bond_list(2,4)


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Bond> '

    call hm_alloc_bonds(nbond)

    ! read psf file
    do ibond = 1, nbond, 4

      if (nbond - ibond < 4) then
        nr = nbond - ibond + 1
      else
        nr = 4
      end if
      read(psf_in,*) (bond_list(1,i), &
                      bond_list(2,i), i = 1,nr)

      do i = 1, nr
        call hm_set_bond_list(1, ibond+i-1, bond_list(1,i))
        call hm_set_bond_list(2, ibond+i-1, bond_list(2,i))
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_angl
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_angl(psf_in, nangl)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nangl

    ! local variables
    integer                  :: iangl, i, nr
    integer                  :: angl_list(3,3)


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Angl> '

    call hm_alloc_angls(nangl)

    ! read psf file
    do iangl = 1, nangl, 3

      if (nangl - iangl < 3) then
        nr = nangl - iangl + 1
      else
        nr = 3
      end if
      read(psf_in,*) (angl_list(1,i), &
                      angl_list(2,i), &
                      angl_list(3,i), i = 1,nr)

      do i = 1, nr
        call hm_set_angl_list(1, iangl+i-1, angl_list(1,i))
        call hm_set_angl_list(2, iangl+i-1, angl_list(2,i))
        call hm_set_angl_list(3, iangl+i-1, angl_list(3,i))
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_dihe
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_dihe(psf_in, ndihe)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: ndihe

    ! local variables
    integer                  :: idihe, i, nr
    integer                  :: dihe_list(4,2)


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Dihe> '

    call hm_alloc_dihes(ndihe)

    ! read psf file
    do idihe = 1, ndihe, 2

      if (ndihe - idihe < 2) then
        nr = ndihe - idihe + 1
      else
        nr = 2
      end if
      read(psf_in,*) (dihe_list(1,i), &
                      dihe_list(2,i), &
                      dihe_list(3,i), &
                      dihe_list(4,i), i = 1,nr)

      do i = 1, nr
        call hm_set_dihe_list(1, idihe+i-1, dihe_list(1,i))
        call hm_set_dihe_list(2, idihe+i-1, dihe_list(2,i))
        call hm_set_dihe_list(3, idihe+i-1, dihe_list(3,i))
        call hm_set_dihe_list(4, idihe+i-1, dihe_list(4,i))
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_impr
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_impr(psf_in, nimpr)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: nimpr

    ! local variables
    integer                  :: iimpr, i, nr
    integer                  :: impr_list(4,2)


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Impr> '

    call hm_alloc_imprs(nimpr)

    ! read psf file
    do iimpr = 1, nimpr, 2

      if (nimpr - iimpr < 2) then
        nr = nimpr - iimpr + 1
      else
        nr = 2
      end if
      read(psf_in,*) (impr_list(1,i), &
                      impr_list(2,i), &
                      impr_list(3,i), &
                      impr_list(4,i), i = 1,nr)

      do i = 1, nr
        call hm_set_impr_list(1, iimpr+i-1, impr_list(1,i))
        call hm_set_impr_list(2, iimpr+i-1, impr_list(2,i))
        call hm_set_impr_list(3, iimpr+i-1, impr_list(3,i))
        call hm_set_impr_list(4, iimpr+i-1, impr_list(4,i))
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_psf_cmap
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_psf_cmap(psf_in, ncmap)

    ! formal arguments
    integer,                 intent(in) :: psf_in
    integer,                 intent(in) :: ncmap

    ! local variables
    integer                  :: icmap, i
    integer                  :: cmap_list(8)


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Psf_Cmap> '

    call hm_alloc_cmaps(ncmap)

    ! read psf file
    do icmap = 1, ncmap

      read(psf_in,*) (cmap_list(i), i=1,8)

      do i = 1, 8
        call hm_set_cmap_list(i, icmap, cmap_list(i))
      end do

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

  end subroutine hm_create_psf_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hm_create_crd_coord
  !> @brief        create huge molecule files
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hm_create_crd_coord(crd_in)

    ! formal arguments
    integer,                 intent(in) :: crd_in

    ! local variables
    real(wp)                 :: atom_coord(3)
    integer                  :: iatom, natom, atom_no, resi_no
    character(200)           :: line
    character(8)             :: flag, resi_name, atom_name, segm_name


    if (main_rank) &
      write(MsgOut,'(a)') 'Hm_Create_Crd_Coord> '

    ! read pdb file
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
                   atom_no,       &
                   resi_no,       &
                   resi_name,     &
                   atom_name,     &
                   atom_coord(1), &
                   atom_coord(2), &
                   atom_coord(3), &
                   segm_name
      else
        read(line,'(i5,i5,1x,a4,1x,a4,3f10.5,1x,a4)') &
                   atom_no,       &
                   resi_no,       &
                   resi_name,     &
                   atom_name,     &
                   atom_coord(1), &
                   atom_coord(2), &
                   atom_coord(3), &
                   segm_name
      end if

      call hm_set_atom_coord(1, iatom, real(atom_coord(1),wip))
      call hm_set_atom_coord(2, iatom, real(atom_coord(2),wip))
      call hm_set_atom_coord(3, iatom, real(atom_coord(3),wip))

    end do

    if (main_rank) &
      write(MsgOut,'(a)') '  done.'

    return

900 call error_msg('Hm_Create_Crd_Coord> read ERROR. ')

  end subroutine hm_create_crd_coord

end module pr_charmm2hm_mod
