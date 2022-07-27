!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_crd_mod
!> @brief   read coordinate data from CHARMM CRD file
!! @authors Yuji Sugita (YS), Jaewoon Jung (JJ), Takaharu Mori (TM), 
!!          Kenta YAMADA (KYMD), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_crd_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_crd

     integer                       :: num_atoms = 0
     character(4),     allocatable :: atom_name(:)    
     character(6),     allocatable :: residue_name(:) 
     character(4),     allocatable :: segment_name(:)
     integer,          allocatable :: atom_no(:)
     integer,          allocatable :: residue_no(:)
     integer,          allocatable :: residue_id(:)
     real(wp),         allocatable :: atom_coord(:,:)

  end type s_crd

  ! parameters for allocatable variables
  integer,     public, parameter   :: CrdAtom    = 1

  ! subroutines
  public  :: input_crd
  public  :: output_crd
  public  :: init_crd
  public  :: alloc_crd  
  public  :: dealloc_crd  
  public  :: dealloc_crd_all
  private :: read_crd
  private :: write_crd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_crd
  !> @brief        a driver subroutine for reading CHARMM CRD file
  !! @authors      YS
  !! @param[in]    crd_filename : filename of CRD file
  !! @param[out]   crd          : structure of CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_crd(crd_filename, crd)

    ! formal arguments
    character(*),            intent(in)    :: crd_filename
    type(s_crd),             intent(inout) :: crd
   
    ! local variables
    integer                  :: file
    

    ! open CRD file
    !
    call open_file(file, crd_filename, IOFileInput)

    ! read coordinate data from CRD file
    !
    call read_crd(file, crd)

    ! close CRD file
    !
    call close_file(file)

    return

  end subroutine input_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_crd
  !> @brief        a driver subroutine for writing CHARMM CRD file
  !! @authors      YS, KY
  !! @param[in]    crd_filename : filename of CRD file
  !! @param[in]    crd          : structure of CRD information
  !! @param[optional]   ioext   : write in an extended form, if true
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_crd(crd_filename, crd, ioext)

    ! formal arguments
    character(*),            intent(in)    :: crd_filename
    type(s_crd),             intent(in)    :: crd
    logical, optional                      :: ioext
   
    ! local variables
    integer                  :: file
    logical                  :: lio
    

    if (present(ioext)) then
       lio = ioext
    else
       lio = .false.
    end if

    ! open CRD file
    !
    call open_file(file, crd_filename, IOFileOutputNew)

    ! write coordinate data from CRD file
    !
    call write_crd(file, crd, lio)

    ! close CRD file
    !
    call close_file(file)

    return

  end subroutine output_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_crd
  !> @brief        initialize CHARMM CRD information
  !! @authors      YS
  !! @param[out]   crd : structure of CHARMM CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_crd(crd)

    ! formal arguments
    type(s_crd),             intent(inout) :: crd


    crd%num_atoms = 0

    return

  end subroutine init_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_crd
  !> @brief        allocate CHARMM CRD information
  !! @authors      YS, KYMD
  !! @param[inout] crd      : structure of CHARMM CRD information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_crd(crd, variable, var_size)

    ! formal arguments
    type(s_crd),             intent(inout) :: crd
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
    !
    select case (variable)

    case(CrdAtom)

      if (allocated(crd%atom_name)) then
        if (size(crd%atom_name) == var_size) return
        deallocate(crd%atom_name,    &
                   crd%residue_name, &
                   crd%segment_name, &
                   crd%atom_no,      &
                   crd%residue_no,   &
                   crd%residue_id,   &
                   crd%atom_coord,   &
                   stat = dealloc_stat)
      end if

      allocate(crd%atom_name(var_size),    &
               crd%residue_name(var_size), &
               crd%segment_name(var_size), &
               crd%atom_no(var_size),      &
               crd%residue_no(var_size),   &
               crd%residue_id(var_size),   &
               crd%atom_coord(3,var_size), &
               stat = alloc_stat)

      crd%atom_name     (1:var_size) = ''  
      crd%residue_name  (1:var_size) = ''
      crd%segment_name  (1:var_size) = ''
      crd%atom_no       (1:var_size) = 0
      crd%residue_no    (1:var_size) = 0
      crd%residue_id    (1:var_size) = 0
      crd%atom_coord(1:3,1:var_size) = 0.0_wp

    case default

      call error_msg('Alloc_Crd> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_crd
  !> @brief        deallocate CHARMM CRD information
  !! @authors      YS, KYMD
  !! @param[inout] crd      : structure of CHARMM CRD information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_crd(crd, variable)

    ! formal arguments
    type(s_crd),             intent(inout) :: crd
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    select case (variable)

    case(CrdAtom)

      if (allocated(crd%atom_name)) then
        deallocate(crd%atom_name,    &
                   crd%residue_name, &
                   crd%segment_name, &
                   crd%atom_no,      &
                   crd%residue_no,   &
                   crd%residue_id,   &
                   crd%atom_coord,   &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Crd> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_crd_all
  !> @brief        deallocate all CHARMM CRD information
  !! @authors      YS
  !! @param[inout] crd : structure of CHARMM CRD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_crd_all(crd)

    ! formal arguments
    type(s_crd),             intent(inout) :: crd


    call dealloc_crd(crd, CrdAtom)

    return

  end subroutine dealloc_crd_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_crds
  !> @brief        read data from CHARMM CRD file
  !! @authors      YS, JJ, TM, KYMD
  !! @param[in]    file : unit number of CHARMM CRD file
  !! @param[out]   crd  : CHARMM CRD data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_crd(file, crd)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_crd),             intent(inout) :: crd

    ! local variables
    character(MaxLine)       :: line
    integer                  :: natm, i, j
    character(8)             :: cres_nm, catm_nm, cseg_nm, cres_id
    character(5)             :: flag


100 read(file,'(a80)',err=900) line
    if (line(1:1) == '*') then
      goto 100
    else
      read(line,'(i10,a5)') natm, flag
    end if

    call alloc_crd(crd, CrdAtom, natm)

    do i = 1, natm

      read(file,'(a120)') line

      if (natm >= 100000 .or. flag == '  EXT') then
        read(line,'(i10,i10,2x,a8,2x,a8,3f20.10,2x,a8,2x,a8)') &      
             crd%atom_no(i),              &
             crd%residue_no(i),           &
             cres_nm,                     &
             catm_nm,                     &
             (crd%atom_coord(j,i),j=1,3), &
             cseg_nm,                     &
             cres_id

        crd%residue_name(i) = cres_nm(1:6)
        crd%atom_name(i)    = catm_nm(1:4)
        crd%segment_name(i) = cseg_nm(1:4)
        read(cres_id, '(i8)') crd%residue_id(i)
      else
        read(line,'(i5,i5,1x,a4,1x,a4,3f10.5,1x,a4,1x,a4)') &
             crd%atom_no(i),              &
             crd%residue_no(i),           &
             crd%residue_name(i),         &
             crd%atom_name(i),            &
             (crd%atom_coord(j,i),j=1,3), &
             crd%segment_name(i),         &
             cres_id

        read(cres_id(1: 4), '(i4)') crd%residue_id(i)
      end if

    end do

    crd%num_atoms = natm

    ! write summary of CRD information
    !
    if (main_rank) then
      write(MsgOut,'(A)')       'Read_Crd> Summary of CRD file'
      write(MsgOut,'(A20,I10)') '  num_atoms       = ', crd%num_atoms
      write(MsgOut,'(A)')       ' '
    end if

    return

900 call error_msg('Read_Crd> read error ')
    
  end subroutine read_crd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_crd
  !> @brief        write data from CHARMM CRD file
  !! @authors      YS, KYMD, KY
  !! @param[in]    file  : unit number of CHARMM CRD file
  !! @param[in]    crd   : CHARMM CRD data
  !! @param[in]    ioext : write in an extended form, if true
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_crd(file, crd, ioext)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_crd),             intent(in)    :: crd
    logical,                 intent(in)    :: ioext

    ! local variables
    integer                  :: natm, i, j


    natm = crd%num_atoms

    if (.not.allocated(crd%atom_coord)) &
      call error_msg('Write_Crd> not allocated: crd%atom_coord')

    if (ioext) then
      write(file,'(i10,''  EXT'')') natm
      do i = 1, natm
          ! SI: Original source code
          !write(file,'(i10,i10,1x,a9,1x,a9,3f20.10,1x,a9,1x,i9)') &
          ! SI: Correct CHARMM EXT Format
          !write(line,'(i10,i10,2x,a8,2x,a8,3f20.10,2x,a8,2x,a8)') &
          ! SI: Pseudo CHARMM EXT Format for GENESIS
          write(file, '(i10,i10,2x,a4,6x,a4,4x,3f20.10,2x,a4,6x,i0)') &
                crd%atom_no(i),              &
                crd%residue_no(i),           &
                crd%residue_name(i),         &
                crd%atom_name(i),            &
                (crd%atom_coord(j,i),j=1,3), &
                crd%segment_name(i),         &
                crd%residue_id(i)
      end do

    else
      write(file,'(i5)') natm
      do i = 1, natm
          write(file,'(i5,i5,1x,a4,1x,a4,3f10.5,1x,a4,1x,i4)') &
                crd%atom_no(i),              &
                crd%residue_no(i),           &
                crd%residue_name(i),         &
                crd%atom_name(i),            &
                (crd%atom_coord(j,i),j=1,3), &
                crd%segment_name(i),         &
                crd%residue_id(i)
      end do

    end if

    return

  end subroutine write_crd
    
end module fileio_crd_mod
