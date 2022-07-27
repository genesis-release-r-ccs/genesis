!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_rtp_mod
!> @brief   read/write GROMACS residue topology file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module fileio_rtp_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rtp_atom
    character(6)          :: name   = ''
    character(6)          :: type   = ''
    real(wp)              :: charge = 0.0_wp
    integer               :: charge_group = 0
  end type s_rtp_atom

  type, public :: s_rtp_bond
    character(6)          :: atom1  = ''
    character(6)          :: atom2  = ''
    real(wp)              :: b0     = 0.0_wp
    real(wp)              :: kb     = 0.0_wp
  end type s_rtp_bond

  type, public :: s_rtp_excl
    character(6)          :: atom1  = ''
    character(6)          :: atom2  = ''
  end type s_rtp_excl

  type, public :: s_rtp_angl
    character(6)          :: atom1  = ''
    character(6)          :: atom2  = ''
    character(6)          :: atom3  = ''
    real(wp)              :: th0    = 0.0_wp
    real(wp)              :: cth    = 0.0_wp
  end type s_rtp_angl

  type, public :: s_rtp_dihe
    character(6)          :: atom1  = ''
    character(6)          :: atom2  = ''
    character(6)          :: atom3  = ''
    character(6)          :: atom4  = ''
    real(wp)              :: phi0   = 0.0_wp
    real(wp)              :: cp     = 0.0_wp
    integer               :: multi  = 0
  end type s_rtp_dihe

  type, public :: s_rtp_impr
    character(6)          :: atom1  = ''
    character(6)          :: atom2  = ''
    character(6)          :: atom3  = ''
    character(6)          :: atom4  = ''
    real(wp)              :: q0     = 0.0_wp
    real(wp)              :: cq     = 0.0_wp
  end type s_rtp_impr

  type, public :: s_rtp_res
    character(6)                  :: name = ''
    integer                       :: attribute = 0
    type(s_rtp_atom), allocatable :: atoms(:)
    type(s_rtp_bond), allocatable :: bonds(:)
    type(s_rtp_excl), allocatable :: excls(:)
    type(s_rtp_angl), allocatable :: angls(:)
    type(s_rtp_dihe), allocatable :: dihes(:)
    type(s_rtp_impr), allocatable :: imprs(:)
  end type s_rtp_res

  type,  public :: s_rtp_backbone
    integer                       :: attribute = 0
    character(6),     allocatable :: atoms(:)
  end type s_rtp_backbone

  type, public :: s_rtp
    integer                           :: bondedtypes(9) = 0
    type(s_rtp_res),      allocatable :: ress(:)
    type(s_rtp_backbone), allocatable :: backbones(:)
  end type s_rtp

  ! constants for bonded types
  integer,      public, parameter :: RTPBondTypeBonds        = 1
  integer,      public, parameter :: RTPBondTypeAngles       = 2
  integer,      public, parameter :: RTPBondTypeDihedrals    = 3
  integer,      public, parameter :: RTPBondTypeImpropers    = 4
  integer,      public, parameter :: RTPBondTypeAllDihedrals = 5
  integer,      public, parameter :: RTPBondTypeNrexcl       = 6
  integer,      public, parameter :: RTPBondTypeNrexclRes    = 7
  integer,      public, parameter :: RTPBondTypeHH14         = 8
  integer,      public, parameter :: RTPBondTypeRemoveDih    = 9

  ! constants for attributes
  integer,      public, parameter :: RTPAttributeProtein     = 1
  integer,      public, parameter :: RTPAttributeDNA         = 2
  integer,      public, parameter :: RTPAttributeRNA         = 3
  integer,      public, parameter :: RTPAttributeLigand      = 4

  ! constants for RTP structure allocatable variables
  integer,      public, parameter :: RTPRes          = 1
  integer,      public, parameter :: RTPBackbone     = 2

  ! constants for RES structure allocatable variables
  integer,      public, parameter :: RTPResAtom      = 1
  integer,      public, parameter :: RTPResBond      = 2
  integer,      public, parameter :: RTPResExcl      = 3
  integer,      public, parameter :: RTPResAngl      = 4
  integer,      public, parameter :: RTPResDihe      = 5
  integer,      public, parameter :: RTPResImpr      = 6

  ! constants for BAKBONE structure allocatable variables
  integer,      public, parameter :: RTPBackboneAtom = 1

  ! constants for directive
  integer,     private, parameter :: DBondedTypes    = 1
  integer,     private, parameter :: D_residue       = 2
  !    res directive
  integer,     private, parameter :: DAtoms          = 100
  integer,     private, parameter :: DBonds          = 101
  integer,     private, parameter :: DExclusions     = 102
  integer,     private, parameter :: DAngles         = 103
  integer,     private, parameter :: DDihedrals      = 104
  integer,     private, parameter :: DImpropers      = 105
  !    not regular directive in GRO-RTP file
  integer,     private, parameter :: DProtein        = 1000
  integer,     private, parameter :: DDNA            = 1001
  integer,     private, parameter :: DRNA            = 1002
  integer,     private, parameter :: DLigand         = 1003
  integer,     private, parameter :: DBackbones      = 1004

  ! subroutines
  public  :: input_rtp
  public  :: alloc_rtp
  public  :: alloc_rtp_res
  public  :: alloc_rtp_backbone
  public  :: dealloc_rtp
  public  :: dealloc_rtp_res
  public  :: dealloc_rtp_backbone
  public  :: dealloc_rtp_all
  public  :: merge_rtp
  private :: read_rtp
  private :: read_bondedtypes
  private :: read_protein
  private :: read_dna
  private :: read_rna
  private :: read_ligand
  private :: read_residue
  private :: read_res_atom
  private :: read_res_bond
  private :: read_res_excl
  private :: read_res_angl
  private :: read_res_dihe
  private :: read_res_impr
  private :: read_backbone
  private :: read_section_lines
  private :: check_directive

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rtp
  !> @brief        a driver subroutine for reading GROMACS RTP file
  !! @authors      NT
  !! @param[in]    rtp_filename : filename of RTP file
  !! @param[inout] rtp          : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rtp(rtp_filename, rtp)

    ! formal arguments
    character(*),            intent(in)    :: rtp_filename
    type(s_rtp),             intent(inout) :: rtp

    ! local variables
    type(s_rtp)              :: rtp0
    integer                  :: file, i, nfile

    character(MaxFilename), allocatable :: files(:)


    ! parse filename string
    !

    nfile = split_num(rtp_filename)
    allocate(files(nfile))

    call split(nfile, nfile, rtp_filename, files)


    do i = 1, nfile

      ! open RTP file
      !
      call open_file(file, files(i), IOFileInput)


      ! read coordinate data from RTP file
      !
      call read_rtp(file, rtp0)


      ! close RTP file
      !
      call close_file(file)


      ! merge RTP data
      !
      call merge_rtp(rtp, rtp0)

      call dealloc_rtp_all(rtp0)

    end do

    deallocate(files)

    return

  end subroutine input_rtp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rtp
  !> @brief        allocate GROMACS RTP data
  !! @authors      NT
  !! @param[inout] rtp      : structure of GROMACS RTP data
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rtp(rtp, variable, var_size)

    ! formal arguments
    type(s_rtp),             intent(inout) :: rtp
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RtpRes)

      if (allocated(rtp%ress)) then
        if (size(rtp%ress) == var_size) return
        deallocate(rtp%ress, &
                   stat = dealloc_stat)
      end if

      allocate(rtp%ress(var_size),&
               stat = alloc_stat)

    case(RtpBackbone)

      if (allocated(rtp%backbones)) then
        if (size(rtp%backbones) == var_size) return
        deallocate(rtp%backbones, &
                   stat = dealloc_stat)
      end if

      allocate(rtp%backbones(var_size),&
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Rtp> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rtp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rtp_res
  !> @brief        allocate GROMACS RTP residue data
  !! @authors      NT
  !! @param[inout] res      : structure of GROMACS RTP residue data
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rtp_res(res, variable, var_size)

    ! formal arguments
    type(s_rtp_res),         intent(inout) :: res
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RtpResAtom)

      if (allocated(res%atoms)) then
        if (size(res%atoms) == var_size) return
        deallocate(res%atoms, stat = dealloc_stat)
      end if

      allocate(res%atoms(var_size), stat = alloc_stat)

    case(RtpResBond)

      if (allocated(res%bonds)) then
        if (size(res%bonds) == var_size) return
        deallocate(res%bonds, stat = dealloc_stat)
      end if

      allocate(res%bonds(var_size), stat = alloc_stat)

    case(RtpResExcl)

      if (allocated(res%excls)) then
        if (size(res%excls) == var_size) return
        deallocate(res%excls, stat = dealloc_stat)
      end if

      allocate(res%excls(var_size), stat = alloc_stat)

    case(RtpResAngl)

      if (allocated(res%angls)) then
        if (size(res%angls) == var_size) return
        deallocate(res%angls, stat = dealloc_stat)
      end if

      allocate(res%angls(var_size), stat = alloc_stat)

    case(RtpResDihe)

      if (allocated(res%dihes)) then
        if (size(res%dihes) == var_size) return
        deallocate(res%dihes, stat = dealloc_stat)
      end if

      allocate(res%dihes(var_size), stat = alloc_stat)

    case(RtpResImpr)

      if (allocated(res%imprs)) then
        if (size(res%imprs) == var_size) return
        deallocate(res%imprs, stat = dealloc_stat)
      end if

      allocate(res%imprs(var_size), stat = alloc_stat)

    case default

      call error_msg('Alloc_Rtp_Res> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rtp_res

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rtp_backbone
  !> @brief        allocate GROMACS RTP backbone data
  !! @authors      NT
  !! @param[inout] backbone : structure of GROMACS RTP backbone data
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rtp_backbone(backbone, variable, var_size)

    ! formal arguments
    type(s_rtp_backbone),    intent(inout) :: backbone
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RtpBackboneAtom)

      if (allocated(backbone%atoms)) then
        if (size(backbone%atoms) == var_size) return
        deallocate(backbone%atoms, stat = dealloc_stat)
      end if

      allocate(backbone%atoms(var_size), stat = alloc_stat)

    case default

      call error_msg('Alloc_Rtp_Backbone> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rtp_backbone

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rtp
  !> @brief        deallocate GROMACS RTP data
  !! @authors      NT
  !! @param[inout] rtp      : structure of GROMACS RTP data
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rtp(rtp, variable)

    ! formal arguments
    type(s_rtp),             intent(inout) :: rtp
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(RtpRes)

      if (allocated(rtp%ress)) then
        deallocate(rtp%ress, &
                   stat = dealloc_stat)
      end if

    case(RtpBackbone)

      if (allocated(rtp%backbones)) then
        deallocate(rtp%backbones, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Rtp> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rtp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rtp_res
  !> @brief        deallocate GROMACS RTP residue data
  !! @authors      NT
  !! @param[inout] res      : structure of GROMACS RTP residue data
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rtp_res(res, variable)

    ! formal arguments
    type(s_rtp_res),         intent(inout) :: res
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RtpResAtom)

      if (allocated(res%atoms)) then
        deallocate(res%atoms, stat = dealloc_stat)
      end if

    case(RtpResBond)

      if (allocated(res%bonds)) then
        deallocate(res%bonds, stat = dealloc_stat)
      end if

    case(RtpResExcl)

      if (allocated(res%excls)) then
        deallocate(res%excls, stat = dealloc_stat)
      end if

    case(RtpResAngl)

      if (allocated(res%angls)) then
        deallocate(res%angls, stat = dealloc_stat)
      end if

    case(RtpResDihe)

      if (allocated(res%dihes)) then
        deallocate(res%dihes, stat = dealloc_stat)
      end if

    case(RtpResImpr)

      if (allocated(res%imprs)) then
        deallocate(res%imprs, stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Rtp_Res> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rtp_res

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rtp_backbone
  !> @brief        deallocate GROMACS RTP backbone data
  !! @authors      NT
  !! @param[inout] backbone : structure of GROMACS RTP backbone data
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rtp_backbone(backbone, variable)

    ! formal arguments
    type(s_rtp_backbone),    intent(inout) :: backbone
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RtpBackboneAtom)

      if (allocated(backbone%atoms)) then
        deallocate(backbone%atoms, stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Rtp_Backbone> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rtp_backbone

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rtp_all
  !> @brief        deallocate all GROMACS RTP data
  !! @authors      NT
  !! @param[inout] rtp : structure of GROMACS RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rtp_all(rtp)

    ! formal arguments
    type(s_rtp),             intent(inout) :: rtp

    ! local variables
    integer                  :: i


    do i = 1, size(rtp%backbones)
      call dealloc_rtp_backbone(rtp%backbones(i), RTPBackboneAtom)
    end do

    do i = 1, size(rtp%ress)
      call dealloc_rtp_res(rtp%ress(i), RTPResAtom)
      call dealloc_rtp_res(rtp%ress(i), RTPResBond)
      call dealloc_rtp_res(rtp%ress(i), RTPResExcl)
      call dealloc_rtp_res(rtp%ress(i), RTPResAngl)
      call dealloc_rtp_res(rtp%ress(i), RTPResDihe)
      call dealloc_rtp_res(rtp%ress(i), RTPResImpr)
    end do

    call dealloc_rtp(rtp, RtpRes)

    return

  end subroutine dealloc_rtp_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_rtp
  !> @brief        merge two RTP data
  !! @authors      NT
  !! @param[inout] rtp0 : structure of GROMACS RTP data
  !! @param[in]    rtp1 : structure of GROMACS RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine merge_rtp(rtp0, rtp1)

    ! formal arguments
    type(s_rtp),             intent(inout) :: rtp0
    type(s_rtp),             intent(in)    :: rtp1

    ! local variables
    type(s_rtp)              :: rtp
    integer                  :: i, n, nn


    ! overwrite bonded types
    !
    n = size(rtp0%bondedtypes)

    do i = 1, n
      if (rtp0%bondedtypes(i) /= 0 .and. &
          rtp0%bondedtypes(i) /= rtp1%bondedtypes(i)) then
        write(MsgOut,'(A)') &
             'Merge_Rtp> WARNING: bondedtypes were overwrote by rtp1.'
        exit
      end if
    end do

    rtp0%bondedtypes(1:n) = rtp1%bondedtypes(1:n)

    ! res
    !
    if (.not. allocated(rtp0%ress)) then

      n = size(rtp1%ress)
      call alloc_rtp(rtp0, RTPRes, n)
      rtp0%ress(1:n) = rtp1%ress(1:n)

    else

      ! copy rtp0 -> rtp
      n = size(rtp0%ress)
      call alloc_rtp(rtp, RTPRes, n)
      rtp%ress(1:n) = rtp0%ress(1:n)

      ! re-allocate rtp0
      call dealloc_rtp(rtp0, RTPRes)
      call alloc_rtp  (rtp0, RTPRes, size(rtp%ress) + size(rtp1%ress))

      ! copy rtp  -> rtp0
      n = size(rtp%ress)
      rtp0%ress(1:n) = rtp%ress(1:n)

      ! copy rtp1 -> rtp0
      nn = size(rtp%ress)
      n  = size(rtp1%ress)
      rtp0%ress(nn+1:nn+n) = rtp1%ress(1:n)

      ! deallocate
      call dealloc_rtp(rtp, RTPRes)

    end if

    ! backbone
    !
    if (.not. allocated(rtp0%backbones)) then

      n = size(rtp1%backbones)
      call alloc_rtp(rtp0, RTPBackbone, n)
      rtp0%backbones(1:n) = rtp1%backbones(1:n)

    else

      ! copy rtp0 -> rtp
      n = size(rtp0%backbones)
      call alloc_rtp(rtp, RTPBackbone, n)
      rtp%backbones(1:n) = rtp0%backbones(1:n)

      ! re-allocate rtp0
      call dealloc_rtp(rtp0, RTPBackbone)
      call alloc_rtp  (rtp0, RTPBackbone, size(rtp%backbones) + &
                                          size(rtp1%backbones))

      ! copy rtp  -> rtp0
      n = size(rtp%backbones)
      rtp0%backbones(1:n) = rtp%backbones(1:n)

      ! copy rtp1 -> rtp0
      nn = size(rtp%backbones)
      n  = size(rtp1%backbones)
      rtp0%backbones(nn+1:nn+n) = rtp1%backbones(1:n)

      ! deallocate
      call dealloc_rtp(rtp, RTPBackbone)

    end if

    return

  end subroutine merge_rtp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rtp
  !> @brief        read data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rtp(file, rtp)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp),             intent(inout) :: rtp

    ! local variables
    integer                  :: directive, nres, nbb, i
    integer                  :: ni 
    character(200)           :: line


    ! pre-read
    !
    nres = 0
    nbb  = 0

    do while(.true.)
      read(file,'(A)',end=100) line

      ! check directive
      if (.not. check_directive(line, directive)) &
        cycle

      ! read per directive
      select case(directive)

      case (D_residue)
        nres = nres + 1

      case (DBackbones)
        nbb  = nbb  + 1

      end select

    end do

    ! actual-read
    !

100 call alloc_rtp(rtp, RtpRes,     nres)
    call alloc_rtp(rtp, RtpBackbone, nbb)

    nres = 0
    nbb  = 0

    rewind(file)

    do while(.true.)
      read(file,'(A)',end=200) line

      ! check directive
      if (.not. check_directive(line, directive)) &
        cycle

      ! read per directive
      select case(directive)

      case (DBondedTypes)
        call read_bondedtypes(file, rtp)

      case (DProtein)
        call read_protein(line, rtp)

      case (DDNA)
        call read_dna(line, rtp)

      case (DRNA)
        call read_rna(line, rtp)

      case (DLigand)
        call read_ligand(line, rtp)

      case (D_residue)
        nres = nres + 1
        call read_residue(line, rtp%ress(nres))

      case (DAtoms)
        call read_res_atom(file, rtp%ress(nres))

      case (DBonds)
        call read_res_bond(file, rtp%ress(nres))

      case (DExclusions)
        call read_res_excl(file, rtp%ress(nres))

      case (DAngles)
        call read_res_angl(file, rtp%ress(nres))

      case (DDihedrals)
        call read_res_dihe(file, rtp%ress(nres))

      case (DImpropers)
        call read_res_impr(file, rtp%ress(nres))

      case (DBackbones)
        nbb = nbb + 1
        call read_backbone(file, rtp%backbones(nbb))

      end select

    end do

200 continue


    ! write summary of RTP information
    !
    if (main_rank) then

      write(MsgOut,'(a)') 'Read_Rtp> Summary of Rtp file'
      write(MsgOut,'(a20,i10)') &
           '  num_residues    = ', size(rtp%ress)
      do i = 1, size(rtp%ress)
        write(MsgOut,'(4x,i2,1x,a6,":",$)') i, rtp%ress(i)%name

        ni = 0
        if (allocated(rtp%ress(i)%atoms)) then
          ni = size(rtp%ress(i)%atoms)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "natom:", ni

        ni = 0
        if (allocated(rtp%ress(i)%bonds)) then
          ni = size(rtp%ress(i)%bonds)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "nbond:", ni

        ni = 0
        if (allocated(rtp%ress(i)%excls)) then
          ni = size(rtp%ress(i)%excls)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "nexcl:", ni

        ni = 0
        if (allocated(rtp%ress(i)%angls)) then
          ni = size(rtp%ress(i)%angls)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "angls:", ni

        ni = 0
        if (allocated(rtp%ress(i)%dihes)) then
          ni = size(rtp%ress(i)%dihes)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "ndihe:", ni

        ni = 0
        if (allocated(rtp%ress(i)%imprs)) then
          ni = size(rtp%ress(i)%imprs)
        endif
        write(MsgOut,'(a6,i3,"|",$)') "nimpr:", ni
        write(MsgOut,'(a)') ''
       
      end do
      write(MsgOut,'(a20,i10)') &
           '  num_backbone def= ', size(rtp%backbones)

      write(MsgOut,'(a)') ''

    end if

    return

  end subroutine read_rtp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_bondedtypes
  !> @brief        read bonded types data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_bondedtypes(file, rtp)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp),             intent(inout) :: rtp

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    if (nline /= 1) &
      call error_msg('Read_Bondedtypes> bad format')

    read(lines(1),*) rtp%bondedtypes(RTPBondTypeBonds:RTPBondTypeRemoveDih)

    return

  end subroutine read_bondedtypes

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_protein
  !> @brief        read protein data from RTP file
  !! @authors      NT
  !! @param[in]    line : directive line string
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_protein(line, rtp)

    ! formal arguments
    character(*),            intent(in)    :: line
    type(s_rtp),             intent(inout) :: rtp


    rtp%ress     (1:size(rtp%ress))     %attribute = RTPAttributeProtein
    rtp%backbones(1:size(rtp%backbones))%attribute = RTPAttributeProtein

    return

  end subroutine read_protein

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_dna
  !> @brief        read DNA data from RTP file
  !! @authors      NT
  !! @param[in]    line : directive line string
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_dna(line, rtp)

    ! formal arguments
    character(*),            intent(in)    :: line
    type(s_rtp),             intent(inout) :: rtp


    rtp%ress     (1:size(rtp%ress))     %attribute = RTPAttributeDNA
    rtp%backbones(1:size(rtp%backbones))%attribute = RTPAttributeDNA

    return

  end subroutine read_dna

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rna
  !> @brief        read RNA data from RTP file
  !! @authors      NT
  !! @param[in]    line : directive line string
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rna(line, rtp)

    ! formal arguments
    character(*),            intent(in)    :: line
    type(s_rtp),             intent(inout) :: rtp


    rtp%ress     (1:size(rtp%ress))     %attribute = RTPAttributeRNA
    rtp%backbones(1:size(rtp%backbones))%attribute = RTPAttributeRNA

    return

  end subroutine read_rna

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ligand
  !> @brief        read ligand data from RTP file
  !! @authors      NT
  !! @param[in]    line : directive line string
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ligand(line, rtp)

    ! formal arguments
    character(*),            intent(in)    :: line
    type(s_rtp),             intent(inout) :: rtp


    rtp%ress     (1:size(rtp%ress))     %attribute = RTPAttributeLigand
    rtp%backbones(1:size(rtp%backbones))%attribute = RTPAttributeLigand

    return

  end subroutine read_ligand

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_residue
  !> @brief        read residue data from RTP file
  !! @authors      NT
  !! @param[in]    line : directive line string
  !! @param[inout] res  : structure of RTP-RES data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_residue(line, res)

    ! formal arguments
    character(*),            intent(in)    :: line
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    character(10)           :: bl, br


    read(line,*) bl, res%name, br

    return

  end subroutine read_residue

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_atom
  !> @brief        read residue atom data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_atom(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResAtom, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%atoms(i)%name,   &
                       res%atoms(i)%type,   &
                       res%atoms(i)%charge, &
                       res%atoms(i)%charge_group
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_bond
  !> @brief        read residue bond data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_bond(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResBond, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%bonds(i)%atom1, &
                       res%bonds(i)%atom2, &
                       res%bonds(i)%b0,    &
                       res%bonds(i)%kb
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_excl
  !> @brief        read residue exclusion data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_excl(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResExcl, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%excls(i)%atom1, &
                       res%excls(i)%atom2
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_excl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_angl
  !> @brief        read residue angle data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_angl(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResAngl, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%angls(i)%atom1, &
                       res%angls(i)%atom2, &
                       res%angls(i)%atom3, &
                       res%angls(i)%th0,   &
                       res%angls(i)%cth
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_dihe
  !> @brief        read residue dihedral data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_dihe(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResDihe, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%dihes(i)%atom1, &
                       res%dihes(i)%atom2, &
                       res%dihes(i)%atom3, &
                       res%dihes(i)%atom4, &
                       res%dihes(i)%phi0,  &
                       res%dihes(i)%cp,    &
                       res%dihes(i)%multi
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_res_impr
  !> @brief        read residue improper data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] rtp  : structure of RTP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_res_impr(file, res)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_res),         intent(inout) :: res

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_res(res, RTPResImpr, nline)

    do i = 1, nline
      read(lines(i),*,end=1) &
                       res%imprs(i)%atom1, &
                       res%imprs(i)%atom2, &
                       res%imprs(i)%atom3, &
                       res%imprs(i)%atom4, &
                       res%imprs(i)%q0,    &
                       res%imprs(i)%cq
1     continue
    end do

    deallocate(lines)

    return

  end subroutine read_res_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_backbone
  !> @brief        read backbone data from RTP file
  !! @authors      NT
  !! @param[in]    file : file unit number
  !! @param[inout] backbone : structure of RTP_BACKBONE data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_backbone(file, backbone)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_rtp_backbone),    intent(inout) :: backbone

    ! local variables
    integer                  :: i, nline
    character(100), allocatable :: lines(:)


    call read_section_lines(file, lines)

    nline = size(lines)

    call alloc_rtp_backbone(backbone, RTPBackboneAtom, nline)

    do i = 1, nline
      read(lines(i),*) backbone%atoms(i)
    end do

    deallocate(lines)

    return

  end subroutine read_backbone

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_section_lines
  !> @brief        read section lines from current file pointer
  !! @authors      NT
  !! @param[in]    file  : file unit number
  !! @param[inout] lines : section lines
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_section_lines(file, lines)

    ! formal arguments
    integer,                   intent(in)    :: file
    character(*), allocatable, intent(inout) :: lines(:)

    ! local variables
    integer                  :: i, nread, nline, directive, alloc_stat
    character(200)           :: line
    logical                  :: alloc_ok


    alloc_ok = .false.

10  nread = 0
    nline = 0

    do while(.true.)
      read(file,'(A)',end=20) line
      if (check_directive(line, directive)) &
        exit
      nread = nread + 1

      i = index(line,';')
      if (i /= 0) &
        line = line(:i-1)
      if (trim(line) /= '') then
        nline = nline + 1
        if (alloc_ok) &
          lines(nline) = line
      end if
    end do

20  do i = 1, nread + 1
      backspace(file)
    end do

    if (.not. alloc_ok) then

      alloc_stat = 0
      allocate(lines(nline), stat=alloc_stat)
      if (alloc_stat /= 0) &
        call error_msg('Read_Section_Lines> memory allocation error.')

      alloc_ok = .true.
      goto 10
    end if

    return

  end subroutine read_section_lines

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      check_directive
  !> @brief        check the directive
  !! @authors      NT
  !! @param[in]    line      : read line
  !! @param[out]   directive : directive ID
  !! @return       flag for directive was found or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function check_directive(line, directive)

    ! return value
    logical  :: check_directive

    ! formal arguments
    character(*),            intent(in)    :: line
    integer,                 intent(out)   :: directive

    ! local variables
    integer                  :: idxb, idxe
    character(20)            :: dir_str


    check_directive = .false.

    idxb = index(line,'[')
    if (idxb == 0) &
       return

    idxe = index(line,']')
    if (idxe == 0) &
       idxe = len(line)

    read(line(idxb+1:idxe-1),*,err=100,end=100) dir_str
100 continue

    select case (dir_str)

    case ('bondedtypes')
      directive = DBondedTypes

    case default
      directive = D_residue

    case ('atoms')
      directive = DAtoms

    case ('bonds')
      directive = DBonds

    case ('exclusions')
      directive = DExclusions

    case ('angles')
      directive = DAngles

    case ('dihedrals')
      directive = DDihedrals

    case ('impropers')
      directive = DImpropers

    case ('protein')
      directive = DProtein

    case ('dna')
      directive = DDNA

    case ('rna')
      directive = DRNA

    case ('ligand')
      directive = DLigand

    case ('backbones')
      directive = DBackbones

    end select

    check_directive = .true.
    return

  end function check_directive

end module fileio_rtp_mod
