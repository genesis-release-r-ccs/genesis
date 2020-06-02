!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_eef1_mod
!> @brief   read EEF1/IMM1 solvation parameter file
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_eef1_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_eef1

    integer                       :: num_atoms
    character(6),     allocatable :: atom_name(:)
    real(wp),         allocatable :: volume(:,:)
    real(wp),         allocatable :: gref(:,:)
    real(wp),         allocatable :: gfree(:,:)
    real(wp),         allocatable :: href(:,:)
    real(wp),         allocatable :: cpref(:,:)
    real(wp),         allocatable :: sigw(:,:)

  end type s_eef1

  ! parameters
  integer,     private, parameter :: MAXROW = 80

  ! subroutines
  public  :: input_eef1
  public  :: alloc_eef1
  public  :: dealloc_eef1
  public  :: read_eef1

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_eef1
  !> @brief        open, read, and close eef1file 
  !! @authors      TM
  !! @param[in]    eef1_filename : filename of eef1file
  !! @param[out]   eef1          : CHARMM EEF1 information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_eef1(eef1_filename, eef1)

    ! formal arguments
    character(*),            intent(in)    :: eef1_filename
    type(s_eef1),            intent(inout) :: eef1

    ! local variables
    integer                  :: unit_no

    ! open eef1file
    !
    call open_file(unit_no, eef1_filename, IOFileInput)

    ! read eef1file
    !
    call read_eef1(unit_no, eef1)

    ! close eef1file
    !
    call close_file(unit_no)


    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Eef1> Summary of eef1file'
      write(MsgOut,'(A20,I10)') '  num_atoms       = ', eef1%num_atoms
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine input_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_eef1
  !> @brief        allocate CHARMM EEF1 information
  !! @authors      TM
  !! @param[inout] eef1     : CHARMM EEF1 information
  !! @param[in]    var_size : allocation size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_eef1(eef1, var_size)

    ! formal arguments
    type(s_eef1),            intent(inout) :: eef1
    integer,                 intent(in)    :: var_size
    
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    character(6),     allocatable :: atom_name(:)
    real(wp),         allocatable :: volume(:,:)
    real(wp),         allocatable :: gref(:,:)
    real(wp),         allocatable :: gfree(:,:)
    real(wp),         allocatable :: href(:,:)
    real(wp),         allocatable :: cpref(:,:)
    real(wp),         allocatable :: sigw(:,:)

    
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    if (allocated(eef1%atom_name)) then
      if (size(eef1%atom_name) == var_size) &
        return
      deallocate(eef1%atom_name,      &
                 eef1%volume,         &
                 eef1%gref,           &
                 eef1%gfree,          &
                 eef1%href,           &
                 eef1%cpref,          &
                 eef1%sigw,           &
                 stat = dealloc_stat)
    end if

    allocate(eef1%atom_name(var_size),   &
             eef1%volume   (2,var_size), &
             eef1%gref     (2,var_size), &
             eef1%gfree    (2,var_size), &
             eef1%href     (2,var_size), &
             eef1%cpref    (2,var_size), &
             eef1%sigw     (2,var_size), &
             stat = alloc_stat)

    eef1%atom_name(1:var_size)     = ''
    eef1%volume   (1:2,1:var_size) = 0.0_wp
    eef1%gref     (1:2,1:var_size) = 0.0_wp
    eef1%gfree    (1:2,1:var_size) = 0.0_wp
    eef1%href     (1:2,1:var_size) = 0.0_wp
    eef1%cpref    (1:2,1:var_size) = 0.0_wp
    eef1%sigw     (1:2,1:var_size) = 0.0_wp

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_eef1
  !> @brief        deallocate CHARMM EEF1 information
  !! @authors      TM
  !! @param[inout] eef1     : CHARMM EEF1 information
  !! @param[in]    variable : an variable to be allocated 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_eef1(eef1)

    ! formal arguments
    type(s_eef1),            intent(inout) :: eef1

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    if (allocated(eef1%atom_name)) then
      deallocate(eef1%atom_name,      &
                 eef1%volume,         &
                 eef1%gref,           &
                 eef1%gfree,          &
                 eef1%href,           &
                 eef1%cpref,          &
                 eef1%sigw,           &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_eef1
  !> @brief        read data from eef1file
  !! @authors      TM
  !! @param[in]    unit_no : unit number of eef1file
  !! @param[out]   eef1    : CHARMM EEF1 inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_eef1(unit_no, eef1)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_eef1),            intent(inout) :: eef1

    ! local variables
    integer       :: i, nchar, nsta, nend, ndata
    integer       :: ichex, iwater
    real(wp)      :: volume, gref, gfree, href, cpref, sigw
    character(80) :: line
    character(6)  :: an
    logical       :: chex_on, water_on, first_end

    ichex      = 0
    iwater     = 0
    chex_on    = .false.
    water_on   = .false.
    first_end  = .false.

    do while(.true.)

      read(unit_no, '(a80)') line

      call char_line(MAXROW, line, nchar)

      if (nchar <= 0) cycle

      nsta = 1
      nend = nchar

      call read_ndata(line, nsta, nend, ndata)

      if (ndata == 1 .and. line(1:4) == 'CHEX') then
        chex_on    = .true.
        water_on   = .false.
      end if
      if (ndata == 1 .and. line(1:5) == 'WATER') then
        chex_on    = .false.
        water_on   = .true.
      end if
      if (ndata == 1 .and. line(1:3) == 'END' .and. first_end) then
        exit
      end if
      if (ndata == 1 .and. line(1:3) == 'END') then
        first_end  = .true.
      end if

      if (chex_on .and. ndata == 7) then
        ichex  = ichex  + 1
      end if
      if (water_on .and. ndata == 7) then
        iwater = iwater + 1
      end if

    end do

    if (ichex /= iwater) then
      call error_msg('Read_Eef1> Inconsistent atom number between CHEX and WATER in eef1file')
    else if (ichex == iwater) then
      eef1%num_atoms = ichex
    else
      call error_msg('Read_Eef1> Error in eef1file')
    end if

    call alloc_eef1(eef1,eef1%num_atoms)

    ! restart reading
    !
    rewind(unit_no)


    ichex      = 0
    iwater     = 0
    chex_on    = .false.
    water_on   = .false.
    first_end  = .false.

    do while(.true.)

      read(unit_no, '(a80)') line

      call char_line(MAXROW, line, nchar)

      if (nchar <= 0) cycle

      nsta = 1
      nend = nchar

      call read_ndata(line, nsta, nend, ndata)

      if (ndata == 1 .and. line(1:4) == 'CHEX') then
        chex_on    = .true.
        water_on   = .false.
      end if
      if (ndata == 1 .and. line(1:5) == 'WATER') then
        chex_on    = .false.
        water_on   = .true.
      end if
      if (ndata == 1 .and. line(1:3) == 'END' .and. first_end) then
        exit
      end if
      if (ndata == 1 .and. line(1:3) == 'END') then
        first_end  = .true.
      end if

      if (chex_on .and. ndata == 7) then

        ichex  = ichex  + 1

        !  Read information
        !
        call read_word(line, nsta, nend, 6, an)
        call read_real(line, nsta, nend, volume)
        call read_real(line, nsta, nend, gref)
        call read_real(line, nsta, nend, gfree)
        call read_real(line, nsta, nend, href)
        call read_real(line, nsta, nend, cpref)
        call read_real(line, nsta, nend, sigw)

        !  Store parameters
        !
        eef1%atom_name(ichex) = an
        eef1%volume(1,ichex)  = volume
        eef1%gref  (1,ichex)  = gref
        eef1%gfree (1,ichex)  = gfree
        eef1%href  (1,ichex)  = href
        eef1%cpref (1,ichex)  = cpref
        eef1%sigw  (1,ichex)  = sigw

      end if

      if (water_on .and. ndata == 7) then

        iwater = iwater + 1

        !  Read information
        !
        call read_word(line, nsta, nend, 6, an)
        call read_real(line, nsta, nend, volume)
        call read_real(line, nsta, nend, gref)
        call read_real(line, nsta, nend, gfree)
        call read_real(line, nsta, nend, href)
        call read_real(line, nsta, nend, cpref)
        call read_real(line, nsta, nend, sigw)

        if (an /= eef1%atom_name(iwater)) then
          call error_msg('Read_Eef1> CHEX and WATER do not have the same order of atom lists')
        end if

        !  Store parameters
        !
        eef1%volume(2,iwater) = volume
        eef1%gref  (2,iwater) = gref
        eef1%gfree (2,iwater) = gfree
        eef1%href  (2,iwater) = href
        eef1%cpref (2,iwater) = cpref
        eef1%sigw  (2,iwater) = sigw

      end if

    end do

    return

  end subroutine read_eef1

end module fileio_eef1_mod
