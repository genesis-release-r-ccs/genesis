!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_psf_mod
!> @brief   read CHARMM PSF file
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ), 
!!          Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_psf_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_psf

    integer                       :: type              = 0

    integer                       :: num_atoms         = 0
    integer                       :: num_bonds         = 0
    integer                       :: num_angles        = 0
    integer                       :: num_dihedrals     = 0
    integer                       :: num_impropers     = 0
    integer                       :: num_HB_donors     = 0
    integer                       :: num_HB_acceptors  = 0
    integer                       :: num_NB_exclusions = 0
    integer                       :: num_groups        = 0
    integer                       :: num_cross_terms   = 0
    real(wp)                      :: total_charge      = 0.0_wp

    ! atom
    integer,          allocatable :: atom_no(:)
    character(4),     allocatable :: segment_name(:)
    integer,          allocatable :: residue_no(:)
    character(6),     allocatable :: residue_name(:)
    character(4),     allocatable :: atom_name(:)
    character(6),     allocatable :: atom_cls_name(:)
    integer,          allocatable :: atom_cls_no(:)
    real(wp),         allocatable :: charge(:)
    real(wp),         allocatable :: mass(:)
    integer,          allocatable :: imove(:)

    integer,          allocatable :: molecule_no(:)

    ! bonds
    integer,          allocatable :: bond_list(:,:)

    ! angles
    integer,          allocatable :: angl_list(:,:)

    ! dihedrals
    integer,          allocatable :: dihe_list(:,:)

    ! impropers
    integer,          allocatable :: impr_list(:,:)

    ! ndonors
    integer,          allocatable :: donr_list(:,:)

    ! acceptors
    integer,          allocatable :: acce_list(:,:)

    ! nb
    integer,          allocatable :: nb_list(:)

    ! grps
    integer,          allocatable :: grp_list(:,:)

    ! cmaps
    integer,          allocatable :: cmap_list(:,:)

    ! numlp

  end type s_psf

  ! parameters for allocatable variables
  integer,      public, parameter :: PsfAtom = 1
  integer,      public, parameter :: PsfBond = 2
  integer,      public, parameter :: PsfAngl = 3
  integer,      public, parameter :: PsfDihe = 4
  integer,      public, parameter :: PsfImpr = 5
  integer,      public, parameter :: PsfDonr = 6
  integer,      public, parameter :: PsfAcce = 7
  integer,      public, parameter :: PsfNb   = 8
  integer,      public, parameter :: PsfGrp  = 9
  integer,      public, parameter :: PsfCmap = 10

  ! parameters
  integer,      public, parameter :: PsfTypeCHARMM = 1
  integer,      public, parameter :: PsfTypeXPLOR  = 2

  ! subroutines
  public  :: input_psf
  public  :: output_psf
  public  :: init_psf
  public  :: alloc_psf
  public  :: dealloc_psf
  public  :: dealloc_psf_all

  private :: read_psf
  private :: read_psf_header
  private :: read_psf_atom
  private :: read_psf_bond
  private :: read_psf_angl
  private :: read_psf_dihe
  private :: read_psf_impr
  private :: read_psf_donr
  private :: read_psf_acce
  private :: read_psf_nb
  private :: read_psf_grp
  private :: read_psf_cmap

  private :: write_psf
  private :: write_psf_header
  private :: write_psf_atom
  private :: write_psf_bond
  private :: write_psf_angl
  private :: write_psf_dihe
  private :: write_psf_impr
  private :: write_psf_donr
  private :: write_psf_acce
  private :: write_psf_nb
  private :: write_psf_grp
  private :: write_psf_cmap
  private :: get_digit

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_psf
  !> @brief        a driver subroutine for reading PSF file
  !! @authors      YS
  !! @param[in]    psf_filename : filename of PSF file
  !! @param[out]   psf          : structure of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_psf(psf_filename, psf)

    ! formal arguments
    character(*),            intent(in)    :: psf_filename
    type(s_psf),             intent(inout) :: psf
    
    ! local variables
    integer                  :: file


    ! open PSF file
    !
    call open_file(file, psf_filename, IOFileInput)

    ! read coordinate data from PSF file
    !
    call read_psf(file, psf)

    ! close PSF file
    !
    call close_file(file)

    return

  end subroutine input_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_psf
  !> @brief        a driver subroutine for writing PSF file
  !! @authors      YS
  !! @param[in]    psf_filename : filename of PSF file
  !! @param[in]    psf          : structure of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_psf(psf_filename, psf)

    ! formal arguments
    character(*),            intent(in)    :: psf_filename
    type(s_psf),             intent(in)    :: psf
    
    ! local variables
    integer                  :: file


    ! open PSF file
    !
    call open_file(file, psf_filename, IOFileOutputNew)

    ! write coordinate data from PSF file
    !
    call write_psf(file, psf)

    ! close PSF file
    !
    call close_file(file)

    return

  end subroutine output_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_psf
  !> @brief        initialize PSF information
  !! @authors      YS
  !! @param[inout] psf : structure of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_psf(psf)

    ! formal arguments
    type(s_psf),             intent(inout) :: psf


    psf%num_atoms         = 0 
    psf%num_bonds         = 0 
    psf%num_angles        = 0
    psf%num_dihedrals     = 0
    psf%num_impropers     = 0
    psf%num_HB_donors     = 0
    psf%num_HB_acceptors  = 0
    psf%num_NB_exclusions = 0
    psf%num_groups        = 0
    psf%num_cross_terms   = 0
    psf%total_charge      = 0.0_wp

    return

  end subroutine init_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_psf
  !> @brief        allocate PSF information
  !! @authors      YS
  !! @param[inout] psf      : structure of PSF information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine alloc_psf(psf, variable, var_size)

    ! formal arguments
    type(s_psf),             intent(inout) :: psf
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    ! initialize
    !
    alloc_stat   = 0
    dealloc_stat = 0


    ! allocate selected variables
    select case (variable)

    case(PsfAtom)

      if (allocated(psf%atom_no)) then
        if (size(psf%atom_no) == var_size) return
        deallocate(psf%atom_no,         &
                   psf%segment_name,    &
                   psf%residue_no,      &
                   psf%residue_name,    &
                   psf%atom_name,       &
                   psf%atom_cls_name,   &
                   psf%atom_cls_no,     &
                   psf%charge,          &
                   psf%mass,            &
                   psf%imove,           &
                   psf%molecule_no,     &
                   stat = dealloc_stat)
      end if

      allocate(psf%atom_no(var_size),         &
               psf%segment_name(var_size),    &
               psf%residue_no(var_size),      &
               psf%residue_name(var_size),    &
               psf%atom_name(var_size),       &
               psf%atom_cls_name(var_size),   &
               psf%atom_cls_no(var_size),     &
               psf%charge(var_size),          &
               psf%mass(var_size),            &
               psf%imove(var_size),           &
               psf%molecule_no(var_size),     &
               stat = alloc_stat)

      psf%atom_no         (1:var_size) = 0
      psf%segment_name    (1:var_size) = ''
      psf%residue_no      (1:var_size) = 0
      psf%residue_name    (1:var_size) = ''
      psf%atom_name       (1:var_size) = ''
      psf%atom_cls_name   (1:var_size) = ''
      psf%atom_cls_no     (1:var_size) = 0
      psf%charge          (1:var_size) = 0.0_wp
      psf%mass            (1:var_size) = 0.0_wp
      psf%imove           (1:var_size) = 0
      psf%molecule_no     (1:var_size) = 0

    case (PsfBond)

      if (allocated(psf%bond_list)) then
        if (size(psf%bond_list(1,:)) == var_size) return
        deallocate(psf%bond_list, stat = dealloc_stat)
      end if

      allocate(psf%bond_list(2, var_size), stat = alloc_stat)

      psf%bond_list(1:2,1:var_size) = 0

    case (PsfAngl)

      if (allocated(psf%angl_list)) then
        if (size(psf%angl_list(1,:)) == var_size) return
        deallocate(psf%angl_list, stat = dealloc_stat)
      end if

      allocate(psf%angl_list(3, var_size), stat = alloc_stat)

      psf%angl_list(1:3,1:var_size) = 0

    case (PsfDihe)
      
      if (allocated(psf%Dihe_list)) then
        if (size(psf%dihe_list(1,:)) == var_size) return
        deallocate(psf%dihe_list, stat = dealloc_stat)
      end if

      allocate(psf%dihe_list(4, var_size), stat = alloc_stat)

      psf%dihe_list(1:4,1:var_size) = 0

    case (PsfImpr)

      if (allocated(psf%impr_list)) then
        if (size(psf%impr_list(1,:)) == var_size) return
        deallocate(psf%impr_list, stat = dealloc_stat)
      end if

      allocate(psf%impr_list(4, var_size), stat = alloc_stat)

      psf%impr_list(1:4,1:var_size) = 0

    case (PsfDonr)
      
      if (allocated(psf%donr_list)) then
        if (size(psf%donr_list(1,:)) == var_size) return
        deallocate(psf%donr_list, stat = dealloc_stat)
      end if

      allocate(psf%donr_list(2, var_size), stat = alloc_stat)

      psf%donr_list(1:2,1:var_size) = 0

    case (PsfAcce)
      
      if (allocated(psf%acce_list)) then
        if (size(psf%acce_list(1,:)) == var_size) return
        deallocate(psf%acce_list, stat = dealloc_stat)
      end if

      allocate(psf%acce_list(2, var_size), stat = alloc_stat)

      psf%acce_list(1:2,1:var_size) = 0

    case (PsfNb)

      if (allocated(psf%nb_list)) then
        if (size(psf%nb_list) == var_size) return
        deallocate(psf%nb_list, stat = dealloc_stat)
      end if

      allocate(psf%nb_list(var_size), stat = alloc_stat)

      psf%nb_list(1:var_size) = 0

    case (PsfGrp)
      
      if (allocated(psf%grp_list)) then
        if (size(psf%grp_list(1,:)) == var_size) return
        deallocate(psf%grp_list, stat = dealloc_stat)
      end if

      allocate(psf%grp_list(3, var_size), stat = alloc_stat)

      psf%grp_list(1:3,1:var_size) = 0

    case (PsfCmap)
      
      if (allocated(psf%cmap_list)) then
        if (size(psf%cmap_list) == var_size) return
        deallocate(psf%cmap_list, stat = dealloc_stat)
      end if

      allocate(psf%cmap_list(8, var_size), stat = alloc_stat)

      psf%cmap_list(1:8,1:var_size) = 0

    case default

      call error_msg('Alloc_Psf> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_psf
  !> @brief        deallocate PSF information
  !! @authors      YS
  !! @param[inout] psf      : structure of PSF information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_psf(psf, variable)

    ! formal arguments
    type(s_psf),             intent(inout) :: psf
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat
   

    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case(PsfAtom)

      if (allocated(psf%atom_no)) then
        deallocate(psf%atom_no,         &
                   psf%segment_name,    &
                   psf%residue_no,      &
                   psf%residue_name,    &
                   psf%atom_name,       &
                   psf%atom_cls_name,   &
                   psf%atom_cls_no,     &
                   psf%charge,          &
                   psf%mass,            &
                   psf%imove,           &
                   psf%molecule_no,     &
                   stat = dealloc_stat)
      end if

    case (PsfBond)

      if (allocated(psf%bond_list)) then
        deallocate(psf%bond_list, stat = dealloc_stat)
      end if

    case (PsfAngl)

      if (allocated(psf%angl_list)) then
         deallocate(psf%angl_list, stat = dealloc_stat)
       end if

    case (PsfDihe)
      
      if (allocated(psf%dihe_list)) then
        deallocate(psf%dihe_list, stat = dealloc_stat)
      end if


    case (PsfImpr)

      if (allocated(psf%impr_list)) then
        deallocate(psf%impr_list, stat = dealloc_stat)
      end if


    case (PsfDonr)

      if (allocated(psf%donr_list)) then
        deallocate(psf%donr_list, stat = dealloc_stat)
      end if


    case (PsfAcce)

      if (allocated(psf%acce_list)) then
        deallocate(psf%acce_list, stat = dealloc_stat)
      end if


    case (PsfNb)

      if (allocated(psf%nb_list)) then
        deallocate(psf%nb_list, stat = dealloc_stat)
      end if


    case (PsfGrp)

      if (allocated(psf%grp_list)) then
        deallocate(psf%grp_list, stat = dealloc_stat)
      end if


    case (PsfCmap)

      if (allocated(psf%cmap_list)) then
        deallocate(psf%cmap_list, stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Psf> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_psf_all
  !> @brief        deallocate all PSF information
  !! @authors      YS
  !! @param[inout] psf : structure of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_psf_all(psf)

    ! formal arguments
    type(s_psf),             intent(inout) :: psf


    call dealloc_psf(psf, PsfAtom)
    call dealloc_psf(psf, PsfBond)
    call dealloc_psf(psf, PsfAngl)
    call dealloc_psf(psf, PsfDihe)
    call dealloc_psf(psf, PsfImpr)
    call dealloc_psf(psf, PsfDonr)
    call dealloc_psf(psf, PsfAcce)
    call dealloc_psf(psf, PsfNb  )
    call dealloc_psf(psf, PsfGrp )
    call dealloc_psf(psf, PsfCmap)

    return

  end subroutine dealloc_psf_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf
  !> @brief        read data from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[out]   psf  : psf data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf(file, psf)
       
    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(out)   :: psf
 
    ! local variables
    integer                  :: n1, n2
    character(100)           :: key
    logical                  :: can_read


    ! deallocate old data
    !
    call dealloc_psf_all(psf)
    call init_psf(psf)

    ! read header section
    !
    call read_psf_header(file)

    do while(.true.)

      can_read = .false.
      read(file,*,err=10,end=10) n1, n2, key
      can_read = .true.

10    if (.not. can_read) then
        backspace(file)
        read(file,*,err=100,end=100) n1, key
      end if

      if (key(1:1) == '!') &
        key = key(2:)

      if (key(1:5) == 'NATOM') then

        ! read atom section
        !
        call read_psf_atom(file, n1, psf)

      else if (key(1:5) == 'NBOND') then

        ! read bond section
        !
        call read_psf_bond(file, n1, psf)

      else if (key(1:6) == 'NTHETA') then

        ! read angle section
        !
        call read_psf_angl(file, n1, psf)

      else if (key(1:4) == 'NPHI') then

        ! read dihedral section
        !
        call read_psf_dihe(file, n1, psf)

      else if (key(1:6) =='NIMPHI') then

        ! read improper section
        !
        call read_psf_impr(file, n1, psf)

      else if (key(1:4) == 'NDON') then

        ! read donor
        !
        call read_psf_donr(file, n1, psf)

      else if (key(1:4) == 'NACC') then

        ! read acceptor
        !
        call read_psf_acce(file, n1, psf)

      else if (key(1:3) == 'NNB') then

        ! read nb
        !
        call read_psf_nb(file, n1, psf)

      else if (key(1:4) == 'NGRP') then

        ! read grp
        !
        call read_psf_grp(file, n1, psf)

      else if (key(1:7) == 'NCRTERM') then

        ! read cmap section
        !
        call read_psf_cmap(file, n1, psf)

      end if

    end do
100 continue

    ! write summary of PSF information
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Psf> Summary of Psffile'
      if (psf%type == PsfTypeXPLOR) then
        write(MsgOut,'(A)') '  psftype         =      xplor'
      else
        write(MsgOut,'(A)') '  psftype         =     charmm'
      endif

      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_atoms       = ', psf%num_atoms,         &
           '  num_bonds       = ', psf%num_bonds
      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_angles      = ', psf%num_angles,        &
           '  num_dihedrals   = ', psf%num_dihedrals
      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_impropers   = ', psf%num_impropers,     &
           '  num_cmap_terms  = ', psf%num_cross_terms
      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_HB_donors   = ', psf%num_HB_donors,     &
           '  num_HB_acceptors= ', psf%num_HB_acceptors
      write(MsgOut,'(A20,I10,A20,I10)')                   &
           '  num_NB_exclusion= ', psf%num_NB_exclusions, &
           '  num_groups      = ', psf%num_groups
      write(MsgOut,'(A20,F10.3)')                         &
           '  total_charge    = ', psf%total_charge
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_header
  !> @brief        read header information of PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_header(file)

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
      call error_msg('Read_Psf_Header> ERROR: Format is not correct')
    end if

    read(file,'(a80)') line
    read(file, *) ntitle
    do i = 1, ntitle
      read(file, '(a80)') title
    end do

    return

  end subroutine read_psf_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_atom
  !> @brief        read atom information from PSF file
  !! @authors      YS, TM
  !! @param[in]    file  : unit number of PSF file
  !! @param[in]    natom : total number of atoms
  !! @param[inout] psf   : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_atom(file, natom, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: natom
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i
    character(8)             :: cstr, cseg_nm, cres_nm, catm_nm
    character(20)            :: ctmp

      
    ! allocate Atom variables
    !
    call alloc_psf(psf, PsfAtom, natom)

    ! Check psf file type
    !
    read(file,*) ctmp, ctmp, ctmp, ctmp, ctmp, cstr
    if (cstr(1:1) >= 'A' .and. cstr(1:1) <= 'Z' .or. &
        cstr(1:1) >= 'a' .and. cstr(1:1) <= 'z') then
      psf%type = PsfTypeXPLOR
    else
       psf%type = PsfTypeCHARMM
    end if

    backspace(file)

    ! Read atom
    !
    psf%num_atoms = natom

    if (psf%type == PsfTypeCHARMM) then

      do i = 1, natom
        read(file, *) psf%atom_no(i),     &
                      cseg_nm,            &
                      psf%residue_no(i),  &
                      cres_nm,            &
                      catm_nm,            &
                      psf%atom_cls_no(i), &
                      psf%charge(i),      &
                      psf%mass(i),        &
                      psf%imove(i)

        psf%segment_name(i) = cseg_nm(1:4)
        psf%residue_name(i) = cres_nm(1:6)
        psf%atom_name(i)    = catm_nm(1:4)

      end do

    else  ! PsfTypeXPLOR 

      do i = 1, natom
        read(file, *) psf%atom_no(i),       &
                      psf%segment_name(i),  &
                      psf%residue_no(i),    &
                      psf%residue_name(i),  &
                      psf%atom_name(i),     &
                      psf%atom_cls_name(i), &
                      psf%charge(i),        &
                      psf%mass(i),          &
                      psf%imove(i)
      end do

    end if

    return

  end subroutine read_psf_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_bond
  !> @brief        read bond information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    nbnd : total number of bonds
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_bond(file, nbnd, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nbnd
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Bond variables
    !
    call alloc_psf(psf, PsfBond, nbnd)


    ! Read Bond
    !
    psf%num_bonds = nbnd

    read(file,*) (psf%bond_list(1, i), &
                  psf%bond_list(2, i), i = 1, nbnd)

    return

  end subroutine read_psf_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_angl
  !> @brief        read angle information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    nang : total number of angles
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_angl(file, nang, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nang
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Angle variables
    !
    call alloc_psf(psf, PsfAngl, nang)


    ! Read Angl
    !
    psf%num_angles = nang

    read(file,*) (psf%angl_list(1, i), &
                  psf%angl_list(2, i), &
                  psf%angl_list(3, i), i = 1, nang)
    
    return

  end subroutine read_psf_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_dihe
  !> @brief        read dihedral angle information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    ndih : total number of dihedral angles
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_dihe(file, ndih, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: ndih
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Dihedral variables
    !
    call alloc_psf(psf, PsfDihe, ndih)


    ! read Dihedral
    !
    psf%num_dihedrals = ndih

    read(file,*) (psf%dihe_list(1, i), &
                  psf%dihe_list(2, i), &
                  psf%dihe_list(3, i), &
                  psf%dihe_list(4, i), i = 1, ndih)

    return

  end subroutine read_psf_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_impr
  !> @brief        read improper angle information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    nimp : total number of improper angles
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_impr(file, nimp, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nimp
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Improper Dihedral variables
    !
    call alloc_psf(psf, PsfImpr, nimp)


    ! read Improper
    !
    psf%num_impropers = nimp

    read(file,*) (psf%impr_list(1, i), &
                  psf%impr_list(2, i), &
                  psf%impr_list(3, i), &
                  psf%impr_list(4, i), i = 1, nimp)

    return

  end subroutine read_psf_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_donr
  !> @brief        read donor information from PSF file
  !! @authors      YS
  !! @param[in]    file  : unit number of PSF file
  !! @param[in]    ndonr : total number of donors
  !! @param[inout] psf   : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_donr(file, ndonr, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: ndonr
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Donor variables
    !
    call alloc_psf(psf, PsfDonr, ndonr)


    ! read donor
    !
    psf%num_HB_donors = ndonr

    read(file,*) (psf%donr_list(1, i), &
                  psf%donr_list(2, i), i = 1, ndonr)

    return

  end subroutine read_psf_donr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_acce
  !> @brief        read acceptor information from PSF file
  !! @authors      YS
  !! @param[in]    file  : unit number of PSF file
  !! @param[in]    nacce : total number of acceptors
  !! @param[inout] psf   : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_acce(file, nacce, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nacce
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate Acceptor variables
    !
    call alloc_psf(psf, PsfAcce, nacce)


    ! read acceptor
    !
    psf%num_HB_acceptors = nacce

    read(file,*) (psf%acce_list(1, i), &
                  psf%acce_list(2, i), i = 1, nacce)

    return

  end subroutine read_psf_acce

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_nb
  !> @brief        read nb information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    nnb  : total number of nbs
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_nb(file, nnb, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nnb
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    ! Allocate NB variables
    !
    call alloc_psf(psf, PsfNb, psf%num_atoms)


    ! read nb
    !
    psf%num_NB_exclusions = nnb

    read(file,*) (psf%nb_list(i), i = 1, psf%num_atoms)

    return

  end subroutine read_psf_nb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_grp
  !> @brief        read group information from PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    ngrp : total number of groups
  !! @param[inout] psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_grp(file, ngrp, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: ngrp
    type(s_psf),             intent(inout) :: psf

    ! integer
    integer                  :: i


    ! Allocate Grp variables
    !
    call alloc_psf(psf, PsfGrp, ngrp)


    ! read nb
    !
    psf%num_groups = ngrp

    read(file,*) (psf%grp_list(1, i), &
                  psf%grp_list(2, i), &
                  psf%grp_list(3, i), i = 1, ngrp)

    return

  end subroutine read_psf_grp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_psf_cmap
  !> @brief        read cmap information from PSF file
  !! @authors      YS
  !! @param[in]    file  : unit number of PSF file
  !! @param[in]    ncmap : total number of cmaps
  !! @param[inout] psf   : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_psf_cmap(file, ncmap, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: ncmap
    type(s_psf),             intent(inout) :: psf

    ! local variables
    integer                  :: i


    !  Allocate cross-term variables
    !
    call alloc_psf(psf, PsfCmap, ncmap)


    ! read cross-terms
    !
    psf%num_cross_terms = ncmap

    read(file,*) (psf%cmap_list(1, i), &
                  psf%cmap_list(2, i), &
                  psf%cmap_list(3, i), &
                  psf%cmap_list(4, i), &
                  psf%cmap_list(5, i), &
                  psf%cmap_list(6, i), &
                  psf%cmap_list(7, i), &
                  psf%cmap_list(8, i), i = 1, ncmap)

    return

  end subroutine read_psf_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf
  !> @brief        write data to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structure of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf


    ! write summary of PSF information
    !
    if (psf%type == PsfTypeXPLOR) then
      write(MsgOut,*) 'Write_Psf> TYPE = XPLOR'
    else
      write(MsgOut,*) 'Write_Psf> TYPE = CHARMM'
    end if

    write(MsgOut,*)'Write_Psf> Number of Atoms         = ',psf%num_atoms
    write(MsgOut,*)'Write_Psf> Number of Bonds         = ',psf%num_bonds
    write(MsgOut,*)'Write_Psf> Number of Angles        = ',psf%num_angles
    write(MsgOut,*)'Write_Psf> Number of Dihedrals     = ',psf%num_dihedrals
    write(MsgOut,*)'Write_Psf> Number of Impropers     = ',psf%num_impropers
    write(MsgOut,*)'Write_Psf> Number of Cross-terms   = ',psf%num_cross_terms
    write(MsgOut,*)'Write_Psf> Number of HB acceptors  = ',psf%num_HB_acceptors
    write(MsgOut,*)'Write_Psf> Number of HB donors     = ',psf%num_HB_donors
    write(MsgOut,*)'Write_Psf> Number of NB exclusions = ',psf%num_NB_exclusions
    write(MsgOut,*)'Write_Psf> Number of Groups        = ',psf%num_groups
    write(MsgOut,*)'Write_Psf> Total charge            = ',psf%total_charge
    write(MsgOut,*)' '

    !  write header section
    !
    call write_psf_header(file)

    !  write atom section
    !
    call write_psf_atom(file, psf)

    !  write bond section
    !
    call write_psf_bond(file, psf)

    !  write angle section
    !
    call write_psf_angl(file, psf)

    !  write dihedral section
    !
    call write_psf_dihe(file, psf)

    !  write improper section
    !
    call write_psf_impr(file, psf)

    !  write donor section
    !
    call write_psf_donr(file, psf)

    !  write acceptor section
    !
    call write_psf_acce(file, psf)

    !  write nb section
    !
    call write_psf_nb(file, psf)

    !  write grp section
    !
    call write_psf_grp(file, psf)

    !  write cmap section
    !
    call write_psf_cmap(file, psf)

    return

  end subroutine write_psf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_header
  !> @brief        write header information of PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_header(file)

    ! parameter
    integer,                 parameter     :: MAXROW = 80

    ! formal arguments
    integer,                 intent(in)    :: file

    ! local variables
    integer                  :: date_time(8)
    character(10)            :: time
    character(8)             :: date
    character(5)             :: zone


    write(file,'(a)') 'PSF CMAP'

    write(file,'(/I8,1X,"!NTITLE")') 1

    call date_and_time(date, time, zone, date_time)
    write(file,'("*  DATE:",i6,"/",i2"/",i2,i7,":",i2,":",i2)') &
             date_time(2), date_time(3), MOD(date_time(1),100), &
             date_time(5), date_time(6), date_time(7)

    return

  end subroutine write_psf_header

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_atom
  !> @brief        write atom information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_atom(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, natom
    character(8)             :: cres_no
    character(2)             :: d

      
    ! write the number of atoms
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    d = get_digit(natom)

    write(file,'(/I'//trim(d)//',1X,"!NATOM")') natom
    

    ! write atom
    !
    if (psf%type == PsfTypeCHARMM) then

      do i = 1, natom
        write(cres_no, '(i4)') psf%residue_no(i)
        write(file, &
            '(I'//trim(d)//',1X,A4,1X,A4,1X,A6,1X,A4,1X,I4,2X,G13.6,G14.6,I8)')&
              psf%atom_no(i),        &
              psf%segment_name(i),   &
              adjustl(cres_no),      &
              psf%residue_name(i),   &
              psf%atom_name(i),      &
              psf%atom_cls_no(i),    &
              psf%charge(i),         &
              psf%mass(i),           &
              psf%imove(i)
      end do

    else  ! PsfTypeXPLOR 

      do i = 1, natom
        write(cres_no, '(i4)') psf%residue_no(i)
        write(file, &
           '(I'//trim(d)//',1X,A4,1X,A4,1X,A6,1X,A4,1X,A6,1X,G13.6,G14.6,I8)') &
              psf%atom_no(i),        &
              psf%segment_name(i),   &
              adjustl(cres_no),      &
              psf%residue_name(i),   &
              psf%atom_name(i),      &
              psf%atom_cls_name(i),  &
              psf%charge(i),         &
              psf%mass(i),           &
              psf%imove(i)
      end do
    end if

    return

  end subroutine write_psf_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_bond
  !> @brief        write bond information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_bond(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, nbnd, natom


    ! write the number of bonds
    !
    if (allocated(psf%bond_list)) then
      nbnd = size(psf%bond_list(1,:))
    else
      nbnd = 0
    end if

    write(file, '(/I'//trim(get_digit(nbnd))//',1X,"!NBOND: bonds")') nbnd
    

    ! write bond
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(8I'//trim(get_digit(natom))//')') &
                         (psf%bond_list(1, i), &
                          psf%bond_list(2, i), i = 1, nbnd)

    return

  end subroutine write_psf_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_angl
  !> @brief        write angle information in PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_angl(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, nang, natom


    ! write the number of angles
    !
    if (allocated(psf%angl_list)) then
      nang = size(psf%angl_list(1,:))
    else
      nang = 0
    end if

    write(file, '(/I'//trim(get_digit(nang))//',1X,"!NTHETA: angles")') nang
    

    ! write angle
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(9I'//trim(get_digit(natom))//')') &
                         (psf%angl_list(1, i), &
                          psf%angl_list(2, i), &
                          psf%angl_list(3, i), i = 1, nang)

    return

  end subroutine write_psf_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_dihe
  !> @brief        write dihedral angle information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_dihe(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, ndih, natom


    ! write the number of angles
    !
    if (allocated(psf%dihe_list)) then
      ndih = size(psf%dihe_list(1,:))
    else
      ndih = 0
    end if

    write(file, '(/I'//trim(get_digit(ndih))//',1X,"!NPHI: dihedrals")') ndih
    

    ! write dihe
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(8I'//trim(get_digit(natom))//')') &
                         (psf%dihe_list(1, i), &
                          psf%dihe_list(2, i), &
                          psf%dihe_list(3, i), &
                          psf%dihe_list(4, i), i = 1, ndih)

    return

  end subroutine write_psf_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_impr
  !> @brief        write improper angle information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_impr(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, nimp, natom


    ! write the number of impropers
    !
    if (allocated(psf%impr_list)) then
      nimp = size(psf%impr_list(1,:))
    else
      nimp = 0
    end if

    write(file, '(/I'//trim(get_digit(nimp))//',1X,"!NIMPHI: impropers")') nimp
    

    ! Write impr
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file,'(8I'//trim(get_digit(natom))//')' ) &
                         (psf%impr_list(1, i), &
                          psf%impr_list(2, i), &
                          psf%impr_list(3, i), &
                          psf%impr_list(4, i), i = 1, nimp)

    return

  end subroutine write_psf_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_donr
  !> @brief        write donor information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_donr(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, ndon, natom


    ! write the number of donors
    !
    if (allocated(psf%donr_list)) then
      ndon = size(psf%donr_list(1,:))
    else
      ndon = 0
    end if

    write(file, '(/I'//trim(get_digit(ndon))//',1X,"!NDON: donors")') ndon

    ! write donors
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(8I'//trim(get_digit(natom))//')') &
                         (psf%donr_list(1, i), &
                          psf%donr_list(2, i), i = 1, ndon)

    return

  end subroutine write_psf_donr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_acce
  !> @brief        write acceptor information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_acce(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, nacc, natom


    ! write the number of acceptors
    !
    if (allocated(psf%acce_list)) then
      nacc = size(psf%acce_list(1,:))
    else
      nacc = 0
    end if

    write(file, '(/I'//trim(get_digit(nacc))//',1X,"!NACC: acceptors")') nacc

    ! write acceptors
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(8I'//trim(get_digit(natom))//')') &
                         (psf%acce_list(1, i), &
                          psf%acce_list(2, i), i = 1, nacc)

    return

  end subroutine write_psf_acce

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_nb
  !> @brief        write nb information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_nb(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, nnb


    nnb = psf%num_atoms

    ! write the number of nb
    !
    write(file,'(/I'//trim(get_digit(psf%num_NB_exclusions))//',1X,"!NNB")') &
         psf%num_NB_exclusions

    ! write nb
    !
    if (allocated(psf%nb_list)) then
      write(file, '(8I'//trim(get_digit(nnb))//')') (psf%nb_list(i), i = 1, nnb)
    else
      write(file, '(8I'//trim(get_digit(nnb))//')') (0, i = 1, nnb)
    end if

    return

  end subroutine write_psf_nb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_grp
  !> @brief        write group information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_grp(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, ngrp, natom


    ! write the number of groups
    !
    if (allocated(psf%grp_list)) then
      ngrp = size(psf%grp_list(1,:))
    else
      ngrp = 0
    end if

    write(file, '(/2I'//trim(get_digit(ngrp))//',1X,"!NGRP")') ngrp, 0

    ! write group
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(9I'//trim(get_digit(natom))//')') &
                         (psf%grp_list(1, i), &
                          psf%grp_list(2, i), &
                          psf%grp_list(3, i), i = 1, ngrp)

    return

  end subroutine write_psf_grp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_psf_cmap
  !> @brief        write cmap information to PSF file
  !! @authors      YS
  !! @param[in]    file : unit number of PSF file
  !! @param[in]    psf  : structureo of PSF information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_psf_cmap(file, psf)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_psf),             intent(in)    :: psf

    ! local variables
    integer                  :: i, ncmap, natom


    ! write the number of cross-term 
    !
    if (allocated(psf%cmap_list)) then
      ncmap = size(psf%cmap_list(1,:))
    else
      ncmap = 0
    end if

    write(file, '(/I'//trim(get_digit(ncmap))//',1X,"!NCRTERM: cross-terms")') &
         ncmap
    

    ! write impr
    !
    if (allocated(psf%atom_no)) then
      natom = size(psf%atom_no)
    else
      natom = 0
    end if

    write(file, '(8I'//trim(get_digit(natom))//')') &
                         (psf%cmap_list(1, i), &
                          psf%cmap_list(2, i), &
                          psf%cmap_list(3, i), &
                          psf%cmap_list(4, i), &
                          psf%cmap_list(5, i), &
                          psf%cmap_list(6, i), &
                          psf%cmap_list(7, i), &
                          psf%cmap_list(8, i), i = 1, ncmap)

    return

  end subroutine write_psf_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_digit
  !> @brief        get the digit number for number of atoms
  !! @authors      NT
  !! @param[in]    n : number of elements
  !! @return       digit number format string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_digit(n)

    ! return values
    character(2)             :: get_digit

    ! formal arguments
    integer,                 intent(in)    :: n

    ! local variables
    integer                  :: cnt
    

    if (n < 10000000) then
      get_digit = '8'
      return
    end if

    cnt = 9

    if (n > 99999999) &
      cnt = cnt + 1

    if (n > 999999999) &
      cnt = cnt + 1

    write(get_digit,'(i0)') cnt
    return

  end function get_digit

end module fileio_psf_mod
