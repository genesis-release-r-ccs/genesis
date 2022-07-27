!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_rstmep_mod
!> @brief   read / write information for rpath-MEP
!! @authors Yoshinobu Akinaga (YA) and Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_rstmep_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rstmep

    logical               :: restart     = .false.
    integer               :: mep_natoms  = 0
    real(wp), allocatable :: mep_coord(:)
    integer               :: qm_natoms   = 0
    real(wp), allocatable :: qm_charge(:)
    integer               :: num_fep     = 0
    real(wp)              :: qm_energy   = 0.0_wp
    real(wp), allocatable :: pfunc(:)

  end type s_rstmep

  logical :: first = .true.

  ! subroutines
  public  :: input_rstmep
  public  :: output_rstmep
  private :: read_rstmep
  private :: write_rstmep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rstmep
  !> @brief        a driver subroutine for reading rstmep file
  !! @authors      YA, KY
  !! @param[in]    filename : filename of rstmep file
  !! @param[out]   rstmep   : structure of rstmep information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rstmep(filename, rstmep)

    ! formal arguments
    character(*),    intent(in)    :: filename
    type(s_rstmep),  intent(inout) :: rstmep
   
    ! local variables
    integer   :: iounit
    

    ! open rstmep file
    !
    call open_file(iounit, filename, IOFileInput)

    ! read coordinate data from mep file
    !
    call read_rstmep(iounit, rstmep)

    ! close rstmep file
    !
    call close_file(iounit)

    return

  end subroutine input_rstmep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_rstmep
  !> @brief        a driver subroutine for writing rstmep file
  !! @authors      YA, KY
  !! @param[in]    filename : filename of rstmep file
  !! @param[in]    rstmep   : structure of rstmep information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_rstmep(filename, rstmep)

    ! formal arguments
    character(*),     intent(in) :: filename
    type(s_rstmep),   intent(in) :: rstmep

    ! local
    integer     :: iounit 


    ! open rstmep file
    !
    call open_file(iounit, filename, IOFileOutputReplace)

    ! write coordinate data from rstmep file
    !
    call write_rstmep(iounit, rstmep)

    ! close rstmep file
    !
    call close_file(iounit)

    return

  end subroutine output_rstmep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rstmep
  !> @brief        read data from rstmep file
  !! @authors      YA, KY
  !! @param[in]    mepunit : unit number of rstmep file
  !! @param[out]   rstmep  : rstmep data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rstmep(mepunit, rstmep)

    ! formal arugments
    integer,          intent(in)    :: mepunit
    type(s_rstmep),   intent(inout) :: rstmep

    ! local variables
    integer :: i, ii
    character(20) :: line


    read(mepunit,*) 
    read(mepunit,*) rstmep%mep_natoms

    allocate(rstmep%mep_coord(3*rstmep%mep_natoms))

    ii = 0
    do i = 1, rstmep%mep_natoms
      read(mepunit, *) rstmep%mep_coord(ii+1:ii+3)
      ii = ii + 3
    end do
    read(mepunit,'(a)') line
    
    if (index(line,'CHARGE') > 0) then
      read(mepunit,*) rstmep%qm_natoms
      allocate(rstmep%qm_charge(rstmep%qm_natoms))
      read(mepunit,*) rstmep%qm_charge
      read(mepunit,*)
      read(mepunit,*) rstmep%qm_energy
      read(mepunit,'(a)') line
    end if
    if (index(line,'FEP') > 0) then
      read(mepunit,*) rstmep%num_fep
      allocate(rstmep%pfunc(2))
      read(mepunit,*) rstmep%pfunc
    end if

    rstmep%restart = .true.

    return

  end subroutine read_rstmep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_rstmep
  !> @brief        write data into rstmep file
  !! @authors      YA, KY
  !! @param[in]    mepunit : unit number of rstmep file
  !! @param[in]    rstmep  : MEP data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_rstmep(mepunit, rstmep)

    ! formal arugments
    integer,          intent(in) :: mepunit
    type(s_rstmep),   intent(in) :: rstmep

    ! local variables
    integer :: i, ii


    write(mepunit, '("MEP COORD")')
    write(mepunit, '(I8)') rstmep%mep_natoms
    ii = 0
    do i = 1, rstmep%mep_natoms
      write(mepunit, '(3E30.20)') rstmep%mep_coord(ii+1:ii+3)
      ii = ii + 3
    end do
    if (allocated(rstmep%qm_charge)) then
      write(mepunit, '("QM CHARGE")')
      write(mepunit, '(I8)') rstmep%qm_natoms
      write(mepunit, '(5E30.20)') rstmep%qm_charge
      write(mepunit, '("QM INTERNAL ENERGY")')
      write(mepunit, '(E30.20)') rstmep%qm_energy
    end if
    if (allocated(rstmep%pfunc)) then
      write(mepunit, '("FEP INFO")')
      write(mepunit, '(I8)') rstmep%num_fep
      write(mepunit, '(2E30.20)') rstmep%pfunc
    end if
    write(mepunit, '("---")')

    return

  end subroutine write_rstmep
    
end module fileio_rstmep_mod
