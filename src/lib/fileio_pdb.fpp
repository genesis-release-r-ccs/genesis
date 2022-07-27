!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_pdb_mod
!> @brief   read coordinate data from PDB file
!! @authors Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_pdb_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pdb

    logical           :: hetatm_rec   = .false. !< flag for HETATM record
    logical           :: cryst_rec    = .false. !< flag for CRYSTAL record
    logical           :: model_rec    = .false. !< flag for MODEL record
    logical           :: ter_rec      = .false. !< flag for TER record
    logical           :: end_rec      = .false. !< flag for END record

    logical           :: atom_col7    = .false. !< flag for 7-column atom number
    logical           :: res_col5     = .false. !< flag for 5-column resi number
    logical           :: res_col6     = .false. !< flag for 6-column resi number
    logical           :: segment      = .false. !< flag for segment number

    integer           :: num_atoms    = 0
    real(wp)          :: pbc_box(3,3) = 0.0_wp
    integer           :: model_no     = 0       !< model number

    ! atom
    logical,          allocatable :: hetatm(:)  !< hetero or atom 
    character(4),     allocatable :: atom_name(:)
    character(4),     allocatable :: residue_name(:)
    character(4),     allocatable :: segment_name(:)
    character(1),     allocatable :: chain_id(:)
    integer,          allocatable :: atom_no(:)
    integer,          allocatable :: residue_no(:)
    real(wp),         allocatable :: atom_coord(:,:)
    real(wp),         allocatable :: atom_occupancy(:)
    real(wp),         allocatable :: atom_temp_factor(:)

    ! ter
    integer,          allocatable :: ter_line_no(:) !< TER record

  end type s_pdb

  ! parameters for allocatable variables
  integer,      public, parameter :: PdbAtom    = 1
  integer,      public, parameter :: PdbTerLine = 2

  ! local variables
  logical,                private :: vervose = .true.

  ! subroutines
  public  :: input_pdb
  public  :: output_pdb
  public  :: init_pdb
  public  :: alloc_pdb
  public  :: dealloc_pdb
  public  :: dealloc_pdb_all
  private :: read_pdb
  private :: write_pdb
  private :: check_pdb

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_pdb
  !> @brief        a driver subroutine for reading PDB file
  !! @authors      YS
  !! @param[in]    pdb_filename : filename of PDB file
  !! @param[out]   pdb          : structure of PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_pdb(pdb_filename, pdb)

    ! formal arguments
    character(*),            intent(in)    :: pdb_filename
    type(s_pdb),             intent(inout) :: pdb
   
    ! local variables
    integer                  :: file
    

    ! open PDB file
    !
    call open_file(file, pdb_filename, IOFileInput)

    ! read coordinate data from PDB file
    !
    call read_pdb(file, pdb)

    ! close PDB file
    !
    call close_file(file)

    return

  end subroutine input_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_pdb
  !> @brief        a driver subroutine for writing PDB file
  !! @authors      YS
  !! @param[in]    pdb_filename : filename of PDB file
  !! @param[in]    pdb          : structure of PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_pdb(pdb_filename, pdb)

    ! formal arguments
    character(*),            intent(in) :: pdb_filename
    type(s_pdb),             intent(in) :: pdb
   
    ! local variables
    integer                  :: file
    

    ! open PDB file
    !
    call open_file(file, pdb_filename, IOFileOutputNew)

    ! write coordinate data from PDB file
    !
    call write_pdb(file, pdb)

    ! close PDB file
    !
    call close_file(file)

    return

  end subroutine output_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_pdb
  !> @brief        initialize PDB information
  !! @authors      YS
  !! @param[out]   pdb : structure of PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_pdb(pdb)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb


    pdb%hetatm_rec       = .false.
    pdb%cryst_rec        = .false.
    pdb%model_rec        = .false.
    pdb%ter_rec          = .false.
    pdb%end_rec          = .false.

    pdb%atom_col7        = .false.
    pdb%res_col5         = .false.
    pdb%res_col6         = .false.
    pdb%segment          = .false.

    pdb%num_atoms        = 0
    pdb%model_no         = 0

    pdb%pbc_box(1:3,1:3) = 0.0_wp

    return

  end subroutine init_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pdb
  !> @brief        allocate PDB information
  !! @authors      YS
  !! @param[inout] pdb      : structure of PDB information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pdb(pdb, variable, var_size)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb
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

    case(PdbAtom)

      if (allocated(pdb%hetatm)) then
        if (size(pdb%hetatm) == var_size) return
        deallocate(pdb%hetatm,           &
                   pdb%atom_name,        &
                   pdb%residue_name,     &
                   pdb%segment_name,     &
                   pdb%chain_id,         &
                   pdb%atom_no,          &
                   pdb%residue_no,       &
                   pdb%atom_coord,       &
                   pdb%atom_occupancy,   &
                   pdb%atom_temp_factor, &
                   stat = dealloc_stat)
      end if

      allocate(pdb%hetatm(var_size),           & 
               pdb%atom_name(var_size),        &
               pdb%residue_name(var_size),     &
               pdb%segment_name(var_size),     &
               pdb%chain_id(var_size),         &
               pdb%atom_no(var_size),          &
               pdb%residue_no(var_size),       &
               pdb%atom_coord(3,var_size),     &
               pdb%atom_occupancy(var_size),   &
               pdb%atom_temp_factor(var_size), &
               stat = alloc_stat)

      pdb%hetatm          (1:var_size) = .false.
      pdb%atom_name       (1:var_size) = ''
      pdb%residue_name    (1:var_size) = ''
      pdb%segment_name    (1:var_size) = ''
      pdb%chain_id        (1:var_size) = ''
      pdb%atom_no         (1:var_size) = 0
      pdb%residue_no      (1:var_size) = 0
      pdb%atom_coord  (1:3,1:var_size) = 0.0_wp
      pdb%atom_occupancy  (1:var_size) = 0.0_wp
      pdb%atom_temp_factor(1:var_size) = 0.0_wp

    case(PdbTerLine)

      if (allocated(pdb%ter_line_no)) then
        if (size(pdb%ter_line_no) == var_size) return
        deallocate(pdb%ter_line_no, stat = dealloc_stat)
      end if

      allocate(pdb%ter_line_no(var_size), stat = alloc_stat)

      pdb%ter_line_no(1:var_size) = 0

    case default

      call error_msg('Alloc_Pdb> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pdb
  !> @brief        deallocate PDB information
  !! @authors      YS
  !! @param[inout] pdb      : structure of PDB information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pdb(pdb, variable)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !  
    select case (variable)

    case(PdbAtom)

      if (allocated(pdb%hetatm)) then
        deallocate (pdb%hetatm,           &
                    pdb%atom_name,        &
                    pdb%residue_name,     &
                    pdb%segment_name,     &
                    pdb%chain_id,         &
                    pdb%atom_no,          &
                    pdb%residue_no,       &
                    pdb%atom_coord,       &
                    pdb%atom_occupancy,   &
                    pdb%atom_temp_factor, &
                    stat = dealloc_stat)
      end if

    case(PdbTerLine)

      if (allocated(pdb%ter_line_no)) then
        deallocate (pdb%ter_line_no, stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Pdb> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pdb_all
  !> @brief        deallocate all PDB information
  !! @authors      YS
  !! @param[inout] pdb : structure of PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pdb_all(pdb)

    ! formal arguments
    type(s_pdb),             intent(inout) :: pdb


    call dealloc_pdb(pdb, PdbAtom)
    call dealloc_pdb(pdb, PdbTerLine)

    return

  end subroutine dealloc_pdb_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pdb
  !> @brief        read data from PDB file
  !! @authors      YS
  !! @param[in]    file : unit number of PDB file
  !! @param[out]   pdb  : pdb data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_pdb(file, pdb)

    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_pdb),             intent(inout) :: pdb

    ! local variables
    integer                  :: i, iatm, iter, num_read_atoms, num_read_ters
    character(80)            :: line
    logical                  :: sk_seg, no_atm
    logical                  :: a_col7, r_col5, r_col6
    logical                  :: no_seg


    ! deallocate old data
    !
    call dealloc_pdb_all(pdb)
    call init_pdb(pdb)


    ! check line count, atom number, segment
    !
    num_read_atoms = 0
    num_read_ters  = 0

    call check_pdb(file, num_read_atoms, num_read_ters, &
                   sk_seg, no_atm, a_col7, r_col5, r_col6)

    ! allocate buffer
    !
    call alloc_pdb(pdb, PdbAtom,    num_read_atoms)
    call alloc_pdb(pdb, PdbTerLine, num_read_ters)

    ! read data from PDB file 
    !
    iatm   = 0
    iter   = 0
    no_seg = .false.

    do while(.true.)

      read(file, fmt='(a80)', end=200, err=910) line

      if (line(1:4) == 'ATOM' .or. line(1:4) == 'HETA') then

        iatm = iatm + 1

        ! atom or hetatm
        pdb%hetatm(iatm)       = (line(1:1) == 'H')

        ! atom name
        if ((line(13:13) .eq. ' ') .or. &
            (line(13:13) > '1' .and. line(13:13) < '9')) then
          write(pdb%atom_name(iatm),'(a3,a1)') line(14:16), line(13:13)
        else
          pdb%atom_name(iatm)  = line(13:16)
        endif

        ! residue name
        pdb%residue_name(iatm) = line(18:21)

        ! chain Id
        if (r_col5 .or. r_col6) then
          pdb%chain_id(iatm)   = ' '
        else
          pdb%chain_id(iatm)   = line(22:22)
        end if

        ! segment name
        if (sk_seg .or. line(73:76) == '    ') then
          pdb%segment_name(iatm) = 'A   '
          no_seg = .true.
        else
          pdb%segment_name(iatm) = line(73:76)
        end if

        ! atom serial no
        if (no_atm) then
          pdb%atom_no(iatm) = iatm
        else if (a_col7) then
          read(line, fmt='(4x,i7)',err=911) pdb%atom_no(iatm)
        else
          read(line, fmt='(6x,i5)',err=911) pdb%atom_no(iatm)
        end if

        ! residue sequence no
        if (r_col6) then
          read(line, fmt='(21x,i6)',err=911) pdb%residue_no(iatm)
        else if (r_col5) then
          read(line, fmt='(21x,i5)',err=911) pdb%residue_no(iatm)
        else
          read(line, fmt='(22x,i4)',err=911) pdb%residue_no(iatm)
        end if

        ! atom coordinate
        ! atom occupancy
        ! atom coordinate
        pdb%atom_occupancy(iatm)   = 0.0_wp
        pdb%atom_temp_factor(iatm) = 0.0_wp

        read(line, fmt='(30x,3f8.3,f6.2,f6.2)', err=911)  &
             (pdb%atom_coord(i,iatm), i=1,3),             &
              pdb%atom_occupancy(iatm),                   &
              pdb%atom_temp_factor(iatm)

      else if (line(1:4) == 'TER ') then

        iter = iter + 1

        ! atom serial no
        pdb%ter_line_no(iter) = 0

        if (.not. no_atm) then
          if (a_col7) then
            read(line, fmt='(4x,i7)', err=300) pdb%ter_line_no(iter)
          else
            read(line, fmt='(6x,i5)', err=300) pdb%ter_line_no(iter)
          end if
        end if

300     continue

        if (pdb%ter_line_no(iter) == 0 .and. iatm /= 0) &
          pdb%ter_line_no(iter) = pdb%atom_no(iatm) + 1

        if (pdb%ter_line_no(iter) == 0) &
          goto 912

        pdb%ter_rec = .true.

      else if (line(1:4) == 'CRYS') then

        ! crystal info
        read(line, fmt='(6x,3f9.3)', err=913) &
             (pdb%pbc_box(i,i), i=1,3)

        pdb%cryst_rec = .true.

      else if (line(1:4) == 'END ') then

        ! end
        pdb%end_rec = .true.

      else if (line(1:4) == 'MODE') then

        ! model
        pdb%model_rec = .true.
        read(line, fmt='(10x,i4)',err=914) pdb%model_no

      end if

    end do

200 continue

    if (no_seg) then
      if (main_rank .and. vervose) then
        write(MsgOut,'(A)') 'Read_Pdb> there are no segment (renamed "A   ")'
        write(MsgOut,'(A)') ' '
      end if
    end if

    pdb%hetatm_rec = .true.
    pdb%atom_col7  = a_col7
    pdb%res_col5   = r_col5
    pdb%res_col6   = r_col6
    pdb%segment    = .not. no_seg
    pdb%num_atoms  = num_read_atoms

    ! write summary of PDB information
    !
    if (main_rank .and. vervose) then
      write(MsgOut,'(A)') 'Read_Pdb> Summary of Data in PDB file'
      write(MsgOut,'(A20,I10)') '  num_atoms       = ', pdb%num_atoms
      write(MsgOut,'(A)') ' '
      vervose = .false.
    end if

    return

910 call error_msg('Read_Pdb> read error (read line) :'//line)
911 call error_msg('Read_Pdb> read error (ATOM/HETATM record) :'//line)
912 call error_msg('Read_Pdb> read error (TER record) :'//line)
913 call error_msg('Read_Pdb> read error (CRYST record) :'//line)
914 call error_msg('Read_Pdb> read error (MODEL record) :'//line)

  end subroutine read_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_pdb
  !> @brief        write data to PDB file
  !! @authors      YS
  !! @param[in]    file : unit number of PDB file
  !! @param[in]    pdb  : structure of PDB information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_pdb(file, pdb)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_pdb),             intent(in)    :: pdb

    ! local variables
    integer                  :: i, j, len
    integer                  :: iter, num_atoms, num_ters
    character(80)            :: fmt_a, fmt_b, fmt_t
    character(6)             :: crec
    character(4)             :: catm, cres, cstr, cseg
    character(7)             :: nchar_atom
    logical                  :: use_cid
    

    if (.not.allocated(pdb%hetatm)) &
      call error_msg('Out_Pdb> not allocated: pdb%hetatm')
    if (.not.allocated(pdb%atom_name)) &
      call error_msg('Out_Pdb> not allocated: pdb%atom_name')
    if (.not.allocated(pdb%residue_name)) &
      call error_msg('Out_Pdb> not allocated: pdb%residue_name')
    if (.not.allocated(pdb%segment_name)) &
      call error_msg('Out_Pdb> not allocated: pdb%segment_name')
    if (.not.allocated(pdb%chain_id)) &
      call error_msg('Out_Pdb> not allocated: pdb%chain_id')
    if (.not.allocated(pdb%atom_no)) &
      call error_msg('Out_Pdb> not allocated: pdb%atom_no')
    if (.not.allocated(pdb%residue_no)) &
      call error_msg('Out_Pdb> not allocated: pdb%residue_no')
    if (.not.allocated(pdb%atom_coord)) &
      call error_msg('Out_Pdb> not allocated: pdb%atom_coord')
    if (.not.allocated(pdb%atom_occupancy)) &
      call error_msg('Out_Pdb> not allocated: pdb%atom_occupancy')
    if (.not.allocated(pdb%atom_temp_factor)) &
      call error_msg('Out_Pdb> not allocated: pdb%atom_temp_factor')

    num_atoms = size(pdb%hetatm)
    num_ters  = size(pdb%ter_line_no)

    if (pdb%atom_col7) then
      if (pdb%res_col6) then
        use_cid = .false.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,i6,48x,a1)'
      else if (pdb%res_col5) then
        use_cid = .false.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_b   = '(a6,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,a1,i4,49x,a1)'
      end if
    else
      if (pdb%res_col6) then
        use_cid = .false.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,i6,48x,a1)'
      else if (pdb%res_col5) then
        use_cid = .false.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_b   = '(a6,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,a1,i4,49x,a1)'
      end if
    end if

    ! write CRYST record
    if (pdb%cryst_rec) then
      write(file, fmt='(a6,3(f9.3),a42)') 'CRYST1', &
            (pdb%pbc_box(j,j), j=1,3), &
            '  90.00  90.00  90.00 P 1           1       '
    end if

    ! write MODEL record
    if (pdb%model_rec) then
      write(file, fmt='(a5,5x,i4,61x,a1)') 'MODEL', &
            pdb%model_no, ' '
    end if


    iter = 1

    do i = 1, num_atoms

      ! write ATOM/HETATM record
      if (pdb%hetatm_rec) then
         if (.not. pdb%hetatm(i) .or. num_atoms > 99999) then
            crec = 'ATOM  '
         else
            crec = 'HETATM'
         end if
      else
         crec = 'ATOM  '
      end if

      if (pdb%atom_col7) then
        write(nchar_atom,'(i7)') i
      else
        if (i <= 99999) then
          write(nchar_atom,'(i5)') i
        else
          write(nchar_atom,'(a5)') '*****'
        end if
      end if

      read(pdb%atom_name(i), *) cstr
      len = len_trim(cstr)
      if (len < 4) then
        write(catm, fmt='(1x,a3)') cstr
      else
        catm = pdb%atom_name(i)
      end if

      read(pdb%residue_name(i),*) cstr
      len = len_trim(cstr)
      if (len == 2) then
        write(cres, fmt='(1x,a2)') cstr
      else if (len == 1) then
        write(cres, fmt='(2x,a1)') cstr
      else
        cres = cstr
      end if

      if (pdb%segment) then
        cseg = pdb%segment_name(i)
      else
        cseg = '    '
      end if

      if (use_cid) then
        if (pdb%residue_no(i) < 10000) then
          write(file, fmt=fmt_a)              &
                crec,                         &
                nchar_atom,                   &
                catm,                         &
                cres,                         &
                pdb%chain_id(i),              &
                pdb%residue_no(i),            &
                (pdb%atom_coord(j,i), j=1,3), &
                pdb%atom_occupancy(i),        &
                pdb%atom_temp_factor(i),      &
                cseg
        else
          write(file, fmt=fmt_b)              &
                crec,                         &
                nchar_atom,                   &
                catm,                         &
                cres,                         &
                pdb%chain_id(i),              &
                pdb%residue_no(i),            &
                (pdb%atom_coord(j,i), j=1,3), &
                pdb%atom_occupancy(i),        &
                pdb%atom_temp_factor(i),      &
                cseg
        end if
      else
        write(file, fmt=fmt_a)              &
              crec,                         &
              nchar_atom,                   &
              catm,                         &
              cres,                         &
              pdb%residue_no(i),            &
              (pdb%atom_coord(j,i), j=1,3), &
              pdb%atom_occupancy(i),        &
              pdb%atom_temp_factor(i),      &
              cseg
      end if

      ! write TER record
      if (.not. pdb%ter_rec .or. iter > num_ters) then
        cycle
      end if
      
      if (pdb%ter_line_no(iter) <= pdb%atom_no(i)) then
        iter = iter + 1
        cycle
      end if

      if (i < num_atoms) then
        if (pdb%ter_line_no(iter) > pdb%atom_no(i+1)) &
          cycle
      end if

      if (pdb%atom_col7) then
        write(nchar_atom,'(i7)') pdb%atom_no(i)+1
      else
        if (i <= 99999) then
          write(nchar_atom,'(i5)') pdb%atom_no(i)+1
        else
          write(nchar_atom,'(a5)') '*****'
        end if
      end if
      if (use_cid) then

        write(file, fmt=fmt_t)   &
             'TER   '          , &
             nchar_atom        , &
             cres              , &
             pdb%chain_id(i)   , &
             pdb%residue_no(i) , &
             ' '

      else

        write(file, fmt=fmt_t)   &
             'TER   '          , &
             nchar_atom        , &
             cres              , &
             pdb%residue_no(i) , &
             ' '

      end if

      iter = iter + 1

    end do

    ! write ENDMDL record
    if (pdb%model_rec) then
      write(file, fmt='(a6,69x,a1)') 'ENDMDL', ' '
    end if

    ! write END record
    if (pdb%end_rec) then
      write(file, fmt='(a6,69x,a1)') 'END   ', ' '
    end if

    return 

  end subroutine write_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_pdb
  !> @brief        check format of PDB file
  !! @authors      YS
  !! @param[in]    file           : unit number of PDB file
  !! @param[out]   num_read_atoms : number of read atoms
  !! @param[out]   num_read_ters  : number of read ters
  !! @param[out]   sk_seg         : flag for skip segment read
  !! @param[out]   no_atm         : flag for atom serial No. is not defined
  !! @param[out]   a_col7         : atom no column size is 7 or not
  !! @param[out]   r_col5         : residue no column size is 5 or not
  !! @param[out]   r_col6         : residue no column size is 6 or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_pdb(file, num_read_atoms, num_read_ters, &
                       sk_seg, no_atm, a_col7, r_col5, r_col6)

    ! formal arugments
    integer,                 intent(in)    :: file
    integer,                 intent(out)   :: num_read_atoms
    integer,                 intent(out)   :: num_read_ters
    logical,                 intent(out)   :: sk_seg
    logical,                 intent(out)   :: no_atm
    logical,                 intent(out)   :: a_col7
    logical,                 intent(out)   :: r_col5
    logical,                 intent(out)   :: r_col6
 
    ! local variables
    character(80)            :: line
    character                :: c1, c2, c3
    logical                  :: ck_seg


    ! initialize
    !
    num_read_atoms = 0
    num_read_ters  = 0
    sk_seg         = .false.
    no_atm         = .false.
    a_col7         = .false.
    r_col5         = .false.
    r_col6         = .false.

    ck_seg         = .true.


    ! check format of PDB data
    !
    do while(.true.)

      read(file, fmt='(a80)', end=100, err=910) line

      if (line(1:4) .eq. 'ATOM' .or. line(1:4) .eq. 'HETA') then

        num_read_atoms = num_read_atoms + 1

        ! check column for segment
        if (ck_seg) then
          if (len_trim(line) < 76) sk_seg = .true.
          ck_seg = .false.
        end if

        ! check atom number overflow
        c3 = line(11:11)
        if (.not. no_atm .and. line(7:11) .eq. '*****') then
          no_atm = .true.
        else if (.not. no_atm .and. (c3<'0' .or. c3>'9')) then
          no_atm = .true.
        end if

        ! check 7-column over atom number
        c1 = line(5:5)
        c2 = line(6:6)
        if ((c1>='0' .and. c1<='9') .or. &
             (c2>='0' .and. c2<='9')) then
          a_col7 = .true.
        end if

        ! check 4-column over residue number
        c1 = line(22:22)
        if (c1>='0' .and. c1<='9') then
          r_col5 = .true.
        end if

        c1 = line(27:27)
        if (c1>='0' .and. c1<='9') then
          r_col5 = .true.
          r_col6 = .true.
        end if

      else if (line(1:4) .eq. 'TER ') then

        num_read_ters = num_read_ters + 1

      end if

    end do

100 rewind(file)

    if (main_rank .and. vervose) then
      if (no_atm) then
        write(MsgOut,'(A)') 'Check_Pdb> No atom serial number (Renumbered)'
        write(MsgOut,'(A)') ' '
      else if (a_col7) then
        write(MsgOut,'(A)') 'Check_Pdb> 7-column atom serial number'
        write(MsgOut,'(A)') ' '
      end if

      if (r_col6) then
        write(MsgOut,'(A)') 'Check_Pdb> 6-column residue sequence number'
        write(MsgOut,'(A)') ' '
      else if (r_col5) then
        write(MsgOut,'(A)') 'Check_Pdb> 5-column residue sequence number'
        write(MsgOut,'(A)') ' '
      end if
    end if

    return

910 call error_msg('Check_Pdb> read error (read line)')

  end subroutine check_pdb

end module fileio_pdb_mod
