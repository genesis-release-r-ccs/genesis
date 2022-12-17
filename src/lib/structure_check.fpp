!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   structure_check_mod
!> @brief   check protein structure
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2020 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module structure_check_mod

  use molecules_str_mod
  use mpi_parallel_mod
  use constants_mod
  use messages_mod
  use string_mod

  implicit none
  private

  ! for ring check
  logical,                   save  :: do_check
  integer,                   save  :: total_nring_groups
  integer,                   save  :: nring_types
  integer,                   save  :: max_nrings_in_residue
  real(wp),     allocatable, save  :: criteria(:,:)
  integer,      allocatable, save  :: num_rings(:)
  integer,      allocatable, save  :: head(:)
  integer,      allocatable, save  :: NA(:,:,:)
  character(4), allocatable, save  :: RA(:,:,:) 
  character(4), allocatable, save  :: EA(:,:,:)
  character(6), allocatable, save  :: ring_resname(:)
  integer,      allocatable, save  :: ring_list(:)
  integer,      allocatable, save  :: grp2resid_ring(:)
  integer,      allocatable, save  :: grp2resno_ring(:)
  real(wp),     allocatable, save  :: grp2criteria_ring(:)
  integer,                   save  :: nexclusion_ring
  integer,      allocatable, save  :: exclist_ring(:)

  ! for chilarity check
  integer,                   save  :: max_chiral_members
  integer,                   save  :: max_chirals_in_residue
  integer,                   save  :: total_nchiral_groups
  integer,                   save  :: nchiral_types
  integer,      allocatable, save  :: num_chirals(:)
  integer,      allocatable, save  :: nchiral_members(:)
  character(4), allocatable, save  :: CA(:,:,:)
  character(6), allocatable, save  :: chiral_resname(:,:)
  integer,      allocatable, save  :: chiral_list(:)
  character(6), allocatable, save  :: grp2resname_chiral(:)
  integer,      allocatable, save  :: grp2resno_chiral(:)
  integer,                   save  :: nexclusion_chiral
  integer,      allocatable, save  :: exclist_chiral(:)

  ! subroutines
  public   :: setup_structure_check
  private  :: setup_charmm_amber
  private  :: setup_charmm19
  private  :: setup_ring_atomlists
  private  :: setup_chiral_atomlists
  public   :: perform_structure_check
  private  :: check_ring_structure
  private  :: check_chirality

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_structure_check
  !> @brief        setup for structure check
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_structure_check(molecule, forcefield, coord, check, &
                 fix_ring, fix_chiral, exclude_ring_grpid, exclude_chiral_grpid)

    ! formal arguments
    type(s_molecule),   intent(inout) :: molecule
    character(20),      intent(in)    :: forcefield
    real(wp),           intent(inout) :: coord(:,:)
    logical,            intent(in)    :: check
    logical,            intent(in)    :: fix_ring
    logical,            intent(in)    :: fix_chiral
    character(MaxLine), intent(in)    :: exclude_ring_grpid
    character(MaxLine), intent(in)    :: exclude_chiral_grpid

    ! local variables
    logical           :: warning


    if (.not. check) return

    if (main_rank) then
      write(MsgOut,'(a)') 'Setup_Structure_Check> Setup for checking the ring size and chirality errors of proteins and DNA/RNA'
    end if

    ! setup structure check
    !
    select case (forcefield)

    case ('CHARMM','AMBER','GROAMBER')

      call setup_charmm_amber(forcefield)
      do_check = .true.

    case ('CHARMM19')

      call setup_charmm19
      do_check = .true.

    case default

      if (main_rank) then
        write(MsgOut,'(a)') '  WARNING! "check_structure = YES" in [MINIMIZE] is not available in this force field'
        write(MsgOut,'(A)') ' '
      end if
      do_check = .false.

      return

    end select


    ! make ring atom lists
    !
    call setup_ring_atomlists(molecule, exclude_ring_grpid)

    ! setup chilarity check
    !
    call setup_chiral_atomlists(molecule, exclude_chiral_grpid)


    ! write summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A20,I10,A20,I10)')                    &
           '  num_ring_grps   = ', total_nring_groups,     &
           '  num_chiral_grps = ', total_nchiral_groups
    end if


    ! check structure for the initial structure
    !
    warning = .false.
    call perform_structure_check(coord, check, fix_ring, fix_chiral, warning)

    if (main_rank) then
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_structure_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_charmm
  !> @brief        setup ring atoms in CHARMM
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_charmm_amber(forcefield)

    ! formal argument
    character(20),      intent(in)    :: forcefield

    ! local variables
    integer           :: i, ista


    ! define ring residues
    !
    if (forcefield == "CHARMM") then
      nring_types = 12
    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then
      nring_types = 31
    end if
    max_nrings_in_residue = 2

    allocate(num_rings   (nring_types))
    allocate(ring_resname(nring_types))
    allocate(criteria    (nring_types,max_nrings_in_residue))
    allocate(head  (max_nrings_in_residue))
    allocate(NA(nring_types,max_nrings_in_residue,2))
    allocate(RA(nring_types,max_nrings_in_residue,9))
    allocate(EA(nring_types,max_nrings_in_residue,9))

    num_rings(:)    = 0
    ring_resname(:) = "      "
    NA(:,:,:)       = 0
    RA(:,:,:)       = "    "
    EA(:,:,:)       = "    "

    ! PHE
    ring_resname(1) = "PHE"
    num_rings(1)   = 1
    NA(1,1,1:2)    = (/6,5/)
    RA(1,1,1:9)    = (/"CG  ", "CD1 ", "CE1 ", "CZ  ", "CE2 ", "CD2 ", "    ", "    ", "    "/)
    EA(1,1,1:9)    = (/"HD1 ", "HE1 ", "HZ  ", "HE2 ", "HD2 ", "    ", "    ", "    ", "    "/)
    criteria(1,1)  = 1.46_wp

    ! TYR
    ring_resname(2) = "TYR"
    num_rings(2)   = 1
    NA(2,1,1:2)    = (/6,6/)
    RA(2,1,1:9)    = (/"CG  ", "CD1 ", "CE1 ", "CZ  ", "CE2 ", "CD2 ", "    ", "    ", "    "/)
    EA(2,1,1:9)    = (/"HD1 ", "HE1 ", "HD2 ", "HE2 ", "OH  ", "HH  ", "    ", "    ", "    "/)
    criteria(2,1)  = 1.47_wp

    ! HSD
    if (forcefield == "CHARMM") then
      ring_resname(3) = "HSD"
    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then
      ring_resname(3) = "HID"
    end if
    num_rings(3)   = 1
    NA(3,1,1:2)    = (/5,3/)
    RA(3,1,1:9)    = (/"CG  ", "CD2 ", "NE2 ", "CE1 ", "ND1 ", "    ", "    ", "    ", "    "/)
    EA(3,1,1:9)    = (/"HD2 ", "HE1 ", "HD1 ", "    ", "    ", "    ", "    ", "    ", "    "/)
    criteria(3,1)  = 1.46_wp

    ! HSE
    if (forcefield == "CHARMM") then
      ring_resname(4) = "HSE"
    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then
      ring_resname(4) = "HIE"
    end if
    num_rings(4)   = 1
    NA(4,1,1:2)    = (/5,3/)
    RA(4,1,1:9)    = (/"CG  ", "CD2 ", "NE2 ", "CE1 ", "ND1 ", "    ", "    ", "    ", "    "/)
    EA(4,1,1:9)    = (/"HD2 ", "HE2 ", "HE1 ", "    ", "    ", "    ", "    ", "    ", "    "/)
    criteria(4,1)  = 1.46_wp

    ! HSP
    if (forcefield == "CHARMM") then
      ring_resname(5) = "HSP"
    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then
      ring_resname(5) = "HIP"
    end if
    num_rings(5)   = 1
    NA(5,1,1:2)    = (/5,4/)
    RA(5,1,1:9)    = (/"CG  ", "CD2 ", "NE2 ", "CE1 ", "ND1 ", "    ", "    ", "    ", "    "/)
    EA(5,1,1:9)    = (/"HD2 ", "HE2 ", "HD1 ", "HE1 ", "    ", "    ", "    ", "    ", "    "/)
    criteria(5,1)  = 1.46_wp

    ! TRP
    ring_resname(6) = "TRP"
    num_rings(6)   = 1
    NA(6,1,1:2)    = (/9,6/)
    RA(6,1,1:9)    = (/"CG  ", "CD1 ", "NE1 ", "CE2 ", "CD2 ", "CE3 ", "CZ3 ", "CH2 ", "CZ2 "/)
    EA(6,1,1:9)    = (/"HD1 ", "HE1 ", "HZ2 ", "HH2 ", "HZ3 ", "HE3 ", "    ", "    ", "    "/)
    criteria(6,1)  = 1.48_wp

    ! PRO
    ring_resname(7) = "PRO"
    num_rings(7)   = 1
    NA(7,1,1:2)    = (/5,7/)
    RA(7,1,1:9)    = (/"N   ", "CA  ", "CB  ", "CG  ", "CD  ", "    ", "    ", "    ", "    "/)
    EA(7,1,1:9)    = (/"HA  ", "HB1 ", "HB2 ", "HG1 ", "HG2 ", "HD1 ", "HD2 ", "    ", "    "/)
    criteria(7,1)  = 1.60_wp


    ! Nucleic acids
    !
    if (forcefield == "CHARMM") then

      ! GUA
      ring_resname(8) = "GUA"
      num_rings(8)   = 2
      NA(8,1,1:2)    = (/5,6/)
      RA(8,1,1:9)    = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
      EA(8,1,1:9)    = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
      criteria(8,1)  = 1.60_wp

      NA(8,2,1:2)    = (/9,6/)
      RA(8,2,1:9)    = (/"N9  ", "C8  ", "N7  ", "C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
      EA(8,2,1:9)    = (/"H8  ", "O6  ", "H1  ", "N2  ", "H21 ", "H22 ", "    ", "    ", "    "/)
      criteria(8,2)  = 1.50_wp

      ! ADE
      ring_resname(9) = "ADE"
      num_rings(9)   = 2
      NA(9,1,1:2)    = (/5,6/)
      RA(9,1,1:9)    = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
      EA(9,1,1:9)    = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
      criteria(9,1)  = 1.60_wp

      NA(9,2,1:2)    = (/9,4/)
      RA(9,2,1:9)    = (/"N9  ", "C8  ", "N7  ", "C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
      EA(9,2,1:9)    = (/"H8  ", "H61 ", "H62 ", "H2  ", "    ", "    ", "    ", "    ", "    "/)
      criteria(9,2)  = 1.47_wp

      ! CYT
      ring_resname(10) = "CYT"
      num_rings(10)  = 2
      NA(10,1,1:2)   = (/5,6/)
      RA(10,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
      EA(10,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
      criteria(10,1) = 1.60_wp

      NA(10,2,1:2)   = (/6,6/)
      RA(10,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
      EA(10,2,1:9)   = (/"O2  ", "N4  ", "H41 ", "H42 ", "H5  ", "H6  ", "    ", "    ", "    "/)
      criteria(10,2) = 1.48_wp

      ! THY
      ring_resname(11) = "THY"
      num_rings(11)  = 2
      NA(11,1,1:2)   = (/5,6/)
      RA(11,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
      EA(11,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
      criteria(11,1) = 1.60_wp

      NA(11,2,1:2)   = (/6,8/)
      RA(11,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
      EA(11,2,1:9)   = (/"O2  ", "H3  ", "O4  ", "C5M ", "H51 ", "H52 ", "H53 ", "H6  ", "    "/)
      criteria(11,2) = 1.50_wp

      ! URA
      ring_resname(12) = "URA"
      num_rings(12)  = 2
      NA(12,1,1:2)   = (/5,6/)
      RA(12,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
      EA(12,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
      criteria(12,1) = 1.60_wp

      NA(12,2,1:2)   = (/6,5/)
      RA(12,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
      EA(12,2,1:9)   = (/"O2  ", "H3  ", "O4  ", "H5  ", "H6  ", "    ", "    ", "    ", "    "/)
      criteria(12,2) = 1.50_wp

    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then

      ! DG, DG5, DG3, G, G5, G3
      ista = 8
      do i = 0, 5
        if      (i == 0) then
          ring_resname(ista+i) = "DG"
        else if (i == 1) then
          ring_resname(ista+i) = "DG5"
        else if (i == 2) then
          ring_resname(ista+i) = "DG3"
        else if (i == 3) then
          ring_resname(ista+i) = "G"
        else if (i == 4) then
          ring_resname(ista+i) = "G5"
        else if (i == 5) then
          ring_resname(ista+i) = "G3"
        end if

        num_rings(ista+i)   = 2
        NA(ista+i,1,1:2)    = (/5,6/)
        RA(ista+i,1,1:9)    = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
        EA(ista+i,1,1:9)    = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
        criteria(ista+i,1)  = 1.60_wp

        NA(ista+i,2,1:2)    = (/9,6/)
        RA(ista+i,2,1:9)    = (/"N9  ", "C8  ", "N7  ", "C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
        EA(ista+i,2,1:9)    = (/"H8  ", "O6  ", "H1  ", "N2  ", "H21 ", "H22 ", "    ", "    ", "    "/)
        criteria(ista+i,2)  = 1.50_wp
      end do

      ! DA, DA5, DA3, A, A5, A3
      ista = 14
      do i = 0, 5
        if      (i == 0) then
          ring_resname(ista+i) = "DA"
        else if (i == 1) then
          ring_resname(ista+i) = "DA5"
        else if (i == 2) then
          ring_resname(ista+i) = "DA3"
        else if (i == 3) then
          ring_resname(ista+i) = "A"
        else if (i == 4) then
          ring_resname(ista+i) = "A5"
        else if (i == 5) then
          ring_resname(ista+i) = "A3"
        end if

        num_rings(ista+i)   = 2
        NA(ista+i,1,1:2)    = (/5,6/)
        RA(ista+i,1,1:9)    = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
        EA(ista+i,1,1:9)    = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
        criteria(ista+i,1)  = 1.60_wp

        NA(ista+i,2,1:2)    = (/9,4/)
        RA(ista+i,2,1:9)    = (/"N9  ", "C8  ", "N7  ", "C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
        EA(ista+i,2,1:9)    = (/"H8  ", "H61 ", "H62 ", "H2  ", "    ", "    ", "    ", "    ", "    "/)
        criteria(ista+i,2)  = 1.47_wp
      end do

      ! DC, DC5, DC3, C, C5, C3
      ista = 20
      do i = 0, 5
        if      (i == 0) then
          ring_resname(ista+i) = "DC"
        else if (i == 1) then
          ring_resname(ista+i) = "DC5"
        else if (i == 2) then
          ring_resname(ista+i) = "DC3"
        else if (i == 3) then
          ring_resname(ista+i) = "C"
        else if (i == 4) then
          ring_resname(ista+i) = "C5"
        else if (i == 5) then
          ring_resname(ista+i) = "C3"
        end if

        num_rings(ista+i)  = 2
        NA(ista+i,1,1:2)   = (/5,6/)
        RA(ista+i,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
        EA(ista+i,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
        criteria(ista+i,1) = 1.60_wp

        NA(ista+i,2,1:2)   = (/6,6/)
        RA(ista+i,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
        EA(ista+i,2,1:9)   = (/"O2  ", "N4  ", "H41 ", "H42 ", "H5  ", "H6  ", "    ", "    ", "    "/)
        criteria(ista+i,2) = 1.48_wp
      end do

      ! DT, DT5, DT3
      ista = 26
      do i = 0, 2
        if      (i == 0) then
          ring_resname(ista+i) = "DT"
        else if (i == 1) then
          ring_resname(ista+i) = "DT5"
        else if (i == 2) then
          ring_resname(ista+i) = "DT3"
        end if

        num_rings(ista+i)  = 2
        NA(ista+i,1,1:2)   = (/5,6/)
        RA(ista+i,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
        EA(ista+i,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
        criteria(ista+i,1) = 1.60_wp

        NA(ista+i,2,1:2)   = (/6,8/)
        RA(ista+i,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
        EA(ista+i,2,1:9)   = (/"O2  ", "H3  ", "O4  ", "C7  ", "H51 ", "H52 ", "H53 ", "H6  ", "    "/)
        criteria(ista+i,2) = 1.50_wp
      end do

      ! U, U3, U5
      ista = 29
      do i = 0, 2
        if      (i == 0) then
          ring_resname(ista+i) = "U"
        else if (i == 1) then
          ring_resname(ista+i) = "U5"
        else if (i == 2) then
          ring_resname(ista+i) = "U3"
        end if

        num_rings(ista+i)  = 2
        NA(ista+i,1,1:2)   = (/5,6/)
        RA(ista+i,1,1:9)   = (/"C4' ", "C3' ", "C2' ", "C1' ", "O4' ", "    ", "    ", "    ", "    "/)
        EA(ista+i,1,1:9)   = (/"H4' ", "H3' ", "O2' ", "H2' ", "H2''", "H1' ", "    ", "    ", "    "/)
        criteria(ista+i,1) = 1.60_wp

        NA(ista+i,2,1:2)   = (/6,5/)
        RA(ista+i,2,1:9)   = (/"N1  ", "C2  ", "N3  ", "C4  ", "C5  ", "C6  ", "    ", "    ", "    "/)
        EA(ista+i,2,1:9)   = (/"O2  ", "H3  ", "O4  ", "H5  ", "H6  ", "    ", "    ", "    ", "    "/)
        criteria(ista+i,2) = 1.50_wp
      end do

    end if


    ! define chiral groups
    !
    nchiral_types          = 6
    max_chiral_members     = 13
    max_chirals_in_residue = 4

    allocate(num_chirals    (nchiral_types                         ))
    allocate(nchiral_members(nchiral_types                         ))
    allocate(chiral_resname (nchiral_types,max_chiral_members      ))
    allocate(CA             (nchiral_types,max_chirals_in_residue,5))

    nchiral_members(:)        = 0
    num_chirals    (:)        = 0
    chiral_resname (:,:)      = "      "
    CA             (:,:,:)    = "    "

    ! amino acids Group1 (common to CHARMM and AMBER)
    nchiral_members(1)        = 13
    num_chirals    (1)        = 1
    chiral_resname (1,  1:13) = (/"ALA","LEU","VAL","SER","MET","PRO","GLN", &
                                  "TRP","ASN","LYS","ARG","PHE","TYR"/)
    CA             (1,1,1:5 ) = (/"CA  ", "C   ", "CB  ", "N   ", "HA  "/)


    ! amino acids Group2
    if (forcefield == "CHARMM") then

      nchiral_members(2)        = 7
      num_chirals    (2)        = 1
      chiral_resname (2,  1:13) = (/"ASP","GLU","HSD","HSE","HSP","CYS","CYM", &
                                    "   ","   ","   ","   ","   ","   "/)
      CA             (2,1,1:5 ) = (/"CA  ", "C   ", "CB  ", "N   ", "HA  "/)

    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then

      nchiral_members(2)        = 10
      num_chirals    (2)        = 1
      chiral_resname (2,  1:13) = (/"ASP","ASH","GLU","GLH","HID","HIE","HIP", &
                                    "CYS","CYM","CYX","   ","   ","   "/)
      CA             (2,1,1:5 ) = (/"CA  ", "C   ", "CB  ", "N   ", "HA  "/)

    end if

    ! THR
    nchiral_members(3)        = 1
    num_chirals    (3)        = 2
    chiral_resname (3,  1:13) = (/"THR","   ","   ","   ","   ","   ","   ", &
                                  "   ","   ","   ","   ","   ","   "/)
    CA             (3,1,1:5 ) = (/"CA  ", "C   ", "CB  ", "N   ", "HA  "/)
    CA             (3,2,1:5 ) = (/"CB  ", "OG1 ", "CG2 ", "CA  ", "HB  "/)


    ! ILE
    nchiral_members(4)        = 1
    num_chirals    (4)        = 2
    chiral_resname (4,  1:13) = (/"ILE","   ","   ","   ","   ","   ","   ", &
                                  "   ","   ","   ","   ","   ","   "/)
    CA             (4,1,1:5 ) = (/"CA  ", "C   ", "CB  ", "N   ", "HA  "/)
    CA             (4,2,1:5 ) = (/"CB  ", "CG1 ", "CG2 ", "CA  ", "HB  "/)


    ! ribose in DNA/RNA
    if (forcefield == "CHARMM") then

      nchiral_members(5)        = 2
      num_chirals    (5)        = 4
      chiral_resname (5,  1:13) = (/"ADE","GUA","   ","   ","   ","   ", &
                                    "   ","   ","   ","   ","   ","   ","   "/)
      CA             (5,1,1:5 ) = (/"C1' ", "C2' ", "N9  ", "O4' ", "H1' "/)
      CA             (5,2,1:5 ) = (/"C2' ", "C3' ", "C1' ", "O2' ", "H2''"/) ! RNA
      CA             (5,3,1:5 ) = (/"C3' ", "C4' ", "C2' ", "O3' ", "H3' "/)
      CA             (5,4,1:5 ) = (/"C4' ", "C5' ", "C3' ", "O4' ", "H4' "/)


      nchiral_members(6)        = 3
      num_chirals    (6)        = 4
      chiral_resname (6,  1:13) = (/"CYT","URA","THY","   ","   ","   ", &
                                    "   ","   ","   ","   ","   ","   ","   "/)
      CA             (6,1,1:5 ) = (/"C1' ", "C2' ", "N1  ", "O4' ", "H1' "/)
      CA             (6,2,1:5 ) = (/"C2' ", "C3' ", "C1' ", "O2' ", "H2''"/) ! RNA
      CA             (6,3,1:5 ) = (/"C3' ", "C4' ", "C2' ", "O3' ", "H3' "/)
      CA             (6,4,1:5 ) = (/"C4' ", "C5' ", "C3' ", "O4' ", "H4' "/)

    else if (forcefield == "AMBER" .or. forcefield == "GROAMBER") then

      nchiral_members(5)        = 12
      num_chirals    (5)        = 4
      chiral_resname (5,  1:13) = (/"DA ","DA5","DA3","A  ","A5 ","A3 ", &
                                    "DG ","DG5","DG3","G  ","G5 ","G3 ","   "/)
      CA             (5,1,1:5 ) = (/"C1' ", "C2' ", "N9  ", "O4' ", "H1' "/)
      CA             (5,2,1:5 ) = (/"C2' ", "C3' ", "C1' ", "O2' ", "H2''"/) ! RNA
      CA             (5,3,1:5 ) = (/"C3' ", "C4' ", "C2' ", "O3' ", "H3' "/)
      CA             (5,4,1:5 ) = (/"C4' ", "C5' ", "C3' ", "O4' ", "H4' "/)


      nchiral_members(6)        = 12
      num_chirals    (6)        = 4
      chiral_resname (6,  1:13) = (/"DC ","DC5","DC3","C  ","C5 ","C3 ", &
                                    "DT ","DT5","DT3","U  ","U5 ","U3 ","   "/)
      CA             (6,1,1:5 ) = (/"C1' ", "C2' ", "N1  ", "O4' ", "H1' "/)
      CA             (6,2,1:5 ) = (/"C2' ", "C3' ", "C1' ", "O2' ", "H2''"/) ! RNA
      CA             (6,3,1:5 ) = (/"C3' ", "C4' ", "C2' ", "O3' ", "H3' "/)
      CA             (6,4,1:5 ) = (/"C4' ", "C5' ", "C3' ", "O4' ", "H4' "/)

    end if

    return

  end subroutine setup_charmm_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_charmm19
  !> @brief        setup ring atoms in CHARMM19
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_charmm19

    ! setup ring residues
    !
    nring_types           = 5
    max_nrings_in_residue = 1

    allocate(num_rings   (nring_types))
    allocate(ring_resname(nring_types))
    allocate(criteria    (nring_types,max_nrings_in_residue))
    allocate(head  (max_nrings_in_residue))
    allocate(NA(nring_types,max_nrings_in_residue,2))
    allocate(RA(nring_types,max_nrings_in_residue,9))
    allocate(EA(nring_types,max_nrings_in_residue,9))

    num_rings(:)    = 0
    ring_resname(:) = "      "
    NA(:,:,:)       = 0
    RA(:,:,:)       = "    "
    EA(:,:,:)       = "    "

    ! PHE
    ring_resname(1) = "PHE"
    num_rings(1)  = 1
    NA(1,1,1:2)   = (/6,0/)
    RA(1,1,1:9)   = (/"CG  ", "CD1 ", "CE1 ", "CZ  ", "CE2 ", "CD2 ", "    ", "    ", "    "/)
    criteria(1,1) = 1.40_wp

    ! TYR
    ring_resname(2) = "TYR"
    num_rings(2)  = 1
    NA(2,1,1:2)   = (/6,2/)
    RA(2,1,1:9)   = (/"CG  ", "CD1 ", "CE1 ", "CZ  ", "CE2 ", "CD2 ", "    ", "    ", "    "/)
    EA(2,1,1:9)   = (/"OH  ", "HH  ", "    ", "    ", "    ", "    ", "    ", "    ", "    "/)
    criteria(2,1) = 1.40_wp

    ! HIS
    ring_resname(3) = "HIS"
    num_rings(3)  = 1
    NA(3,1,1:2)   = (/5,1/)
    RA(3,1,1:9)   = (/"CG  ", "CD2 ", "NE2 ", "CE1 ", "ND1 ", "    ", "    ", "    ", "    "/)
    EA(3,1,1:9)   = (/"HD1 ", "    ", "    ", "    ", "    ", "    ", "    ", "    ", "    "/)
    criteria(3,1) = 1.40_wp

    ! TRP
    ring_resname(4) = "TRP"
    num_rings(4)  = 1
    NA(4,1,1:2)   = (/9,1/)
    RA(4,1,1:9)   = (/"CG  ", "CD1 ", "NE1 ", "CE2 ", "CD2 ", "CE3 ", "CZ3 ", "CH2 ", "CZ2 "/)
    EA(4,1,1:9)   = (/"HE1 ", "    ", "    ", "    ", "    ", "    ", "    ", "    ", "    "/)
    criteria(4,1) = 1.40_wp

    ! PRO
    ring_resname(5) = "PRO"
    num_rings(5)  = 1
    NA(5,1,1:2)   = (/5,0/)
    RA(5,1,1:9)   = (/"N   ", "CA  ", "CB  ", "CG  ", "CD  ", "    ", "    ", "    ", "    "/) 
    criteria(5,1) = 1.62_wp


    ! define chiral groups
    !
    nchiral_types          = 0
    max_chiral_members     = 0
    max_chirals_in_residue = 0


    return

  end subroutine setup_charmm19

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_ring_atomlists
  !> @brief        make atom lists used in bond length calculation
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    exclude_ring_gripid : exclude ring gripid
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_ring_atomlists(molecule, exclude_ring_grpid)

    ! formal arguments
    type(s_molecule),     target, intent(in) :: molecule
    character(MaxLine),           intent(in) :: exclude_ring_grpid

    ! local variables
    integer                :: i, j, k, l, m, jres
    integer                :: prev_head, prev_jres, prev_resno, curr_resno
    integer                :: icycle, ifound_grp, ifound_atm
    logical                :: residue_changed, is_ring_residue
    character(4)           :: curr_segid, prev_segid
    integer                :: natoms

    integer,       pointer :: resno(:)
    character(4),  pointer :: segid(:)
    character(4),  pointer :: atomname(:)
    character(6),  pointer :: resname(:)


    natoms   =  molecule%num_atoms
    resno    => molecule%residue_no
    segid    => molecule%segment_name
    atomname => molecule%atom_name
    resname  => molecule%residue_name

    do icycle = 1, 2

      ! count the number of ring residues in the system
      !
      ifound_grp = 0
      ifound_atm = 0
      do i = 1, natoms
        do j = 1, nring_types
          do k = 1, num_rings(j)
            if (resname(i) .eq. ring_resname(j) .and. &
                atomname(i) .eq. RA(j,k,1)) then
              ifound_grp = ifound_grp + 1
              ifound_atm = ifound_atm + NA(j,k,1) + NA(j,k,2) + 2
              if (icycle == 2) then
                grp2resid_ring   (ifound_grp) = j
                grp2resno_ring   (ifound_grp) = resno(i)
                grp2criteria_ring(ifound_grp) = criteria(j,k)
              end if
            end if
          end do
        end do
      end do

      if (icycle == 1) then
        allocate(ring_list(ifound_atm), grp2resid_ring(ifound_grp), &
                 grp2resno_ring(ifound_grp), grp2criteria_ring(ifound_grp))
      end if

    end do
    total_nring_groups = ifound_grp

    ! make a list of ring atoms
    !
    ring_list(:) = -1
    jres = 0

    do i = 1, natoms

      ! check residue was changed or not
      !
      curr_resno = resno(i)
      curr_segid = segid(i)
      residue_changed = .true.
      if (i > 1) then
        if (curr_resno == prev_resno) then
          if (curr_segid == prev_segid) then
            residue_changed = .false.
          end if
        end if
      end if
      prev_resno = resno(i)
      prev_segid = segid(i)

      ! check ring residue or not
      !
      if (residue_changed) then
        is_ring_residue = .false.
        do j = 1, nring_types
          if (resname(i) .eq. ring_resname(j)) then
            is_ring_residue = .true.
            prev_jres = jres
            jres = j
            exit
          end if
        end do
      end if

      if (.not. is_ring_residue) cycle

      ! define header information in ring_list
      !
      if (residue_changed) then
        do k = 1, num_rings(jres)
          if (k == 1) then
            if (prev_jres /= 0) then
              head(1) = prev_head
              do l = 1, num_rings(prev_jres)
                head(1) = head(1) + NA(prev_jres,l,1) + NA(prev_jres,l,2) + 2
              end do
            else
              head(1) = 1
            end if
          else
            head(k) = head(k-1) + NA(jres,k-1,1) + NA(jres,k-1,2) + 2
          end if
          ring_list(head(k)  ) = NA(jres,k,1)
          ring_list(head(k)+1) = NA(jres,k,2)
        end do
        residue_changed = .false.
        prev_head = head(1)
      end if

      ! make ring_list
      !
      do k = 1, num_rings(jres)
        ! ring atoms
        do m = 1, NA(jres,k,1)
          if (atomname(i) .eq. RA(jres,k,m)) then
            ring_list(head(k)+m+1) = i
            exit
          end if
        end do
        ! extra atoms
        do m = 1, NA(jres,k,2)
          if (atomname(i) .eq. EA(jres,k,m)) then
            ring_list(head(k)+NA(jres,k,1)+m+1) = i
            exit
          end if
        end do
      end do
    end do

    ! make exclusion list
    !
    nexclusion_ring = split_num(trim(exclude_ring_grpid))
    if (nexclusion_ring /= 0) then
      allocate(exclist_ring(nexclusion_ring))
      call split(nexclusion_ring, nexclusion_ring, exclude_ring_grpid, exclist_ring)
    end if


    return

  end subroutine setup_ring_atomlists

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_chiral_atomlists
  !> @brief        make atom lists used in bond length calculation
  !! @authors      TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    exclude_ring_gripid : exclude ring gripid
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_chiral_atomlists(molecule, exclude_chiral_grpid)

    ! formal arguments
    type(s_molecule), target, intent(in) :: molecule
    character(MaxLine),       intent(in) :: exclude_chiral_grpid

    ! local variables
    integer                :: i, j, k, l, m
    integer                :: icycle, ifound_grp, jtype, lcount
    integer                :: ihead0, prev_ncatoms
    integer                :: prev_jres, prev_resno, curr_resno
    logical                :: residue_changed, is_chiral_residue
    character(4)           :: curr_segid, prev_segid

    integer                :: natoms

    integer,       pointer :: resno(:)
    character(4),  pointer :: segid(:)
    character(4),  pointer :: atomname(:)
    character(6),  pointer :: resname(:)


    natoms   =  molecule%num_atoms
    resno    => molecule%residue_no
    segid    => molecule%segment_name
    atomname => molecule%atom_name
    resname  => molecule%residue_name

    DO icycle = 1, 2

      ! count the number of chiral groups in the system
      !
      ifound_grp = 0
      do i = 1, natoms
        do j = 1, nchiral_types
          do k = 1, num_chirals(j)
            if (atomname(i) == CA(j,k,1)) then
              do m = 1, nchiral_members(j)
                if (resname(i) .eq. chiral_resname(j,m)) then
                  ifound_grp = ifound_grp + 1
                  if (icycle == 2) then
                    grp2resname_chiral(ifound_grp) = resname(i)
                    grp2resno_chiral  (ifound_grp) = resno  (i)
                  end if
                end if
              end do
            end if
          end do
        end do
      end do

      if (icycle == 1) then
        allocate(chiral_list(ifound_grp*5), grp2resname_chiral(ifound_grp), &
                 grp2resno_chiral(ifound_grp))
      end if

    END DO
    total_nchiral_groups = ifound_grp


    ! make a list of chiral atoms
    !
    chiral_list(:) = -1
    jtype = 0
    ihead0 = 1
    prev_ncatoms = 0

    do i = 1, natoms

      ! check residue was changed or not
      !
      curr_resno = resno(i)
      curr_segid = segid(i)
      residue_changed = .true.
      if (i > 1) then
        if (curr_resno == prev_resno) then
          if (curr_segid == prev_segid) then
            residue_changed = .false.
          end if
        end if
      end if
      prev_resno = resno(i)
      prev_segid = segid(i)

      ! check chiral residue or not
      !
      if (residue_changed) then
        is_chiral_residue = .false.
        do j = 1, nchiral_types
          do k = 1, nchiral_members(j)
            if (resname(i) .eq. chiral_resname(j,k)) then
              is_chiral_residue = .true.
              jtype = j
              exit
            end if
          end do
          if (is_chiral_residue) exit
        end do
      end if

      if (.not. is_chiral_residue) cycle

      if (residue_changed) then
        ihead0 = ihead0 + prev_ncatoms
        prev_ncatoms = num_chirals(jtype)*5
        residue_changed = .false.
      end if

      ! make chiral_list
      !
      lcount = 0
      do k = 1, num_chirals(jtype)
        do l = 1, 5
          lcount = lcount + 1
          if (atomname(i) .eq. CA(jtype,k,l)) then
            chiral_list(ihead0 + lcount - 1) = i
          end if
        end do
      end do

    end do

    ! make exclusion list
    !
    nexclusion_chiral = split_num(trim(exclude_chiral_grpid))
    if (nexclusion_chiral /= 0) then
      allocate(exclist_chiral(nexclusion_chiral))
      call split(nexclusion_chiral, nexclusion_chiral, exclude_chiral_grpid, exclist_chiral)
    end if


    return

  end subroutine setup_chiral_atomlists

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    perform_structure_check
  !> @brief        check protein structure
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine perform_structure_check(coord, check, fix_ring, fix_chiral, warning)

    ! formal arguments
    real(wp), intent(inout) :: coord(:,:)
    logical,  intent(in)    :: check
    logical,  intent(in)    :: fix_ring
    logical,  intent(in)    :: fix_chiral
    logical,  intent(in)    :: warning


    ! check aromatic ring size
    !
    call check_ring_structure(coord, check, fix_ring, warning)

    ! check chilarity
    !
    call check_chirality(coord, check, fix_chiral, warning)


    return

  end subroutine perform_structure_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_ring_structure
  !> @brief        check ring bond length and fix ring penetration error
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_ring_structure(coord, check, fix, warning)

    ! formal arguments
    real(wp), intent(inout) :: coord(:,:)
    logical,  intent(in)    :: check
    logical,  intent(in)    :: fix
    logical,  intent(in)    :: warning

    ! local variables
    integer  :: i, j
    integer  :: i1, i2, nr, ne
    integer  :: ihead, icount, jcount
    real(wp) :: length, max_length
    logical  :: sus_found, def_error, do_fix


    if (.not. check) return
    if (.not. do_check) return

    if (main_rank .and. warning) then
      write(MsgOut,'(A)') 'Check_Ring_Structure> Check ring structure'
    end if

    ! Compute bond length
    !
    sus_found = .false.
    icount    = 0
    do i = 1, total_nring_groups

      ! number of ring atoms
      icount = icount + 1
      nr = ring_list(icount)

      ! number of extra atoms
      icount = icount + 1
      ne = ring_list(icount)

      ! compute bond length in the ring
      !
      ihead = icount + 1
      max_length = 0.0_wp
      do j = 1, nr - 1
        icount = icount + 1
        i1 = ring_list(icount    )
        i2 = ring_list(icount + 1)
        if (i1 == -1 .or. i2 == -1) cycle
        length = sqrt( (coord(1,i1) - coord(1,i2))**2 &
                     + (coord(2,i1) - coord(2,i2))**2 &
                     + (coord(3,i1) - coord(3,i2))**2 )
        if (length > max_length) max_length = length
      end do
      icount = icount + 1

      do j = 1, ne
        icount = icount + 1
      end do


      ! error detected
      !
      if (max_length > grp2criteria_ring(i)) then

        sus_found = .true.

        if (main_rank) then
          if (warning) then
            write(MsgOut,'(a31,i10,a3,a7,i5,a9,i10,a1,a19,f8.3)')                           &
              '  suspicious ring group id   = ', i, ' : ', ring_resname(grp2resid_ring(i)), &
              grp2resno_ring(i), ' (atom = ', ring_list(ihead), ')',                        &
              ' max_bond_length = ', max_length
          end if
        end if

        ! Reduce ring size
        !
        if (fix) then

          ! exclusion
          do_fix = .true.
          do j = 1, nexclusion_ring
            if (i == exclist_ring(j)) then
              do_fix = .false.
              exit
            end if
          end do

          if (do_fix) then
            if (main_rank) then
              write(MsgOut,'(a24,i10,a3,a7,i5,a9,i10,a1,a19,f8.3)')                    &
                '  fix ring group id   = ', i, ' : ', ring_resname(grp2resid_ring(i)), &
                grp2resno_ring(i), ' (atom = ', ring_list(ihead), ')',                 &
                ' max_bond_length = ', max_length
            end if

            jcount = 0
            i1 = ring_list(ihead)
            do j = 1, nr + ne - 1
              jcount = jcount + 1
              i2     = ring_list(ihead + jcount)
              if (i1 == -1 .or. i2 == -1) cycle
              coord(1,i2) = (coord(1,i1)*8.0_wp + coord(1,i2)*2.0_wp)/10.0_wp
              coord(2,i2) = (coord(2,i1)*8.0_wp + coord(2,i2)*2.0_wp)/10.0_wp
              coord(3,i2) = (coord(3,i1)*8.0_wp + coord(3,i2)*2.0_wp)/10.0_wp
            end do
          end if

        end if

      end if

    end do

    if (main_rank) then
      if (warning) then
        if (sus_found) then
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '  WARNING!'
          write(MsgOut,'(A)') '  Some suspicious residues were detected. Minimization might be too short, or "ring penetration"'
          write(MsgOut,'(A)') '  might happen in the above residues. Check the structure of those residues very carefully.'
          write(MsgOut,'(A)') '  If you found a ring penetration, try to perform an energy minimization again '
          write(MsgOut,'(A)') '  with the options "check_structure = YES" and "fix_ring_error = YES" in [MINIMIZE].'
          write(MsgOut,'(A)') '  The energy minimization should be restarted from the restart file obtained in "this" run.'
          write(MsgOut,'(A)') '  For more information, see the chapter on [MINIMIZE] in the user manual.'
          write(MsgOut,'(A)') ' '
        else
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '  No suspicious residue was detected.'
          write(MsgOut,'(A)') ' '
        end if
      end if
    end if

    return

  end subroutine check_ring_structure

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_chirality
  !> @brief        check chirality in amino acids and RNA/DNA
  !! @authors      TM
  !! @param[in]    coord : coordinates of atoms in the system
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_chirality(coord, check, fix, warning)

    ! formal arguments
    real(wp), intent(inout) :: coord(:,:)
    logical,  intent(in)    :: check
    logical,  intent(in)    :: fix
    logical,  intent(in)    :: warning

    ! local variables
    integer  :: i, j, i1, i2, i3, i4, i5, icount
    real(wp) :: vec42(3), vec43(3), vec15(3), nvec(3), tmp
    real(wp) :: norm15, normn
    logical  :: sus_found, do_fix

    if (.not. check) return
    if (.not. do_check) return

    if (main_rank .and. warning) then
      write(MsgOut,'(A)') 'Check_Chirality> Check chirality'
    end if

    sus_found = .false.

    icount = 0
    do i = 1, total_nchiral_groups

      icount = icount + 1
      i1     = chiral_list(icount)
      icount = icount + 1
      i2     = chiral_list(icount)
      icount = icount + 1
      i3     = chiral_list(icount)
      icount = icount + 1
      i4     = chiral_list(icount)
      icount = icount + 1
      i5     = chiral_list(icount)

      ! skip definition error
      !
      if (i1 == -1 .or. i2 == -1 .or. i3 == -1 .or. i4 == -1 .or. i5 == -1) then
        cycle
      end if

      ! compute vectors around the chirality center
      !
      vec43(1:3) = coord(1:3,i3) - coord(1:3,i4)
      vec42(1:3) = coord(1:3,i2) - coord(1:3,i4)
      vec15(1:3) = coord(1:3,i5) - coord(1:3,i1)

      nvec(1) = vec43(2)*vec42(3) - vec43(3)*vec42(2)
      nvec(2) = vec43(3)*vec42(1) - vec43(1)*vec42(3)
      nvec(3) = vec43(1)*vec42(2) - vec43(2)*vec42(1)

      norm15 = sqrt(vec15(1)**2 + vec15(2)**2 + vec15(3)**2)
      normn  = sqrt(nvec (1)**2 + nvec (2)**2 + nvec (3)**2)
      tmp    = dot_product(vec15,nvec)/(norm15*normn)

      if (tmp < 0.7_wp) then

        sus_found = .true.

        if (main_rank) then
          if (warning) then
            write(MsgOut,'(a31,i10,a3,a7,i5,a9,i10,a1,a9,f8.3)') &
              '  suspicious chiral group id = ', i, ' : ',       &
              grp2resname_chiral(i), grp2resno_chiral(i),        &
              ' (atom = ', i1, ') ',                             &
              ' angle = ', acos(tmp)/RAD
          end if
        end if

        ! change hydrogen atom position
        !
        if (fix) then

          ! exclusion
          do_fix = .true.
          do j = 1, nexclusion_chiral
            if (i == exclist_chiral(j)) then
              do_fix = .false.
              exit
            end if
          end do

          if (do_fix) then
            if (main_rank) then
              write(MsgOut,'(a24,i10,a3,a7,i5,a9,i10,a1,a9,f8.3)') &
                '  fix chiral group id = ', i, ' : ',              &
                grp2resname_chiral(i), grp2resno_chiral(i),        &
                ' (atom = ', i1, ')',                              &
                ' angle = ', acos(tmp)/RAD
            end if
            coord(1:3,i5) = coord(1:3,i1) - vec15(1:3)
          end if

        end if
      end if

    end do

    if (main_rank) then
      if (warning) then
        if (sus_found) then
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '  WARNING!'
          write(MsgOut,'(A)') '  Some suspicious residues were detected. Minimization might be too short, or "chirality error"'
          write(MsgOut,'(A)') '  might happen in the above residues. Check the structure of those residues very carefully.'
          write(MsgOut,'(A)') '  If you found a chirality error, try to perform an energy minimization again'
          write(MsgOut,'(A)') '  with the options "check_structure = YES" and "fix_chirality_error = YES" in [MINIMIZE].'
          write(MsgOut,'(A)') '  The energy minimization should be restarted from the restart file obtained in "this" run.'
          write(MsgOut,'(A)') '  For more information, see the chapter on [MINIMIZE] in the user manual.'
          write(MsgOut,'(A)') ' '
        else
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '  No suspicious residue was detected.'
          write(MsgOut,'(A)') ' '
        end if
      end if
    end if

    return

  end subroutine check_chirality

end module structure_check_mod
