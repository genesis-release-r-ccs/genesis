!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   trajectory_str_mod
!> @brief   structure of trajectory information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module trajectory_str_mod

  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_trajectory
    real(wp),         allocatable :: coord(:,:)
    real(wp)                      :: pbc_box(3,3)
  end type s_trajectory

  type, public :: s_trj_list
    character(MaxFilename), allocatable :: filenames(:)
    integer,                allocatable :: md_steps(:)
    integer,                allocatable :: mdout_periods(:)
    integer,                allocatable :: ana_periods(:)
    integer,                allocatable :: start_steps(:)
    integer                             :: trj_format
    integer                             :: trj_type
  end type s_trj_list

  type, public :: s_trj_file
    character(MaxFilename)              :: name
    integer                             :: trj_format
    integer                             :: trj_type
    integer                             :: in_out
                                  
    integer                             :: unit_no
    logical                             :: byte_swap
    logical                             :: rw_header
    integer                             :: step_no
  end type s_trj_file

  ! parameters
  integer,      public, parameter :: TrjFormatPDB       = 1
  integer,      public, parameter :: TrjFormatAmber     = 2
  integer,      public, parameter :: TrjFormatDCD       = 3
  integer,      public, parameter :: TrjFormatGromacs   = 4
  integer,      public, parameter :: TrjFormatCharmmRst = 5
  integer,      public, parameter :: TrjFormatNamdRst   = 6

  integer,      public, parameter :: TrjTypeCoor        = 1
  integer,      public, parameter :: TrjTypeCoorBox     = 2

  character(*), public, parameter :: TrjFormatTypes(6)  = (/'PDB       ', &
                                                            'AMBER     ', &
                                                            'DCD       ', &
                                                            'GROMACS   ', &
                                                            'CHARMM_RST', &
                                                            'NAMD_RST  '/)

  character(*), public, parameter :: TrjTypeTypes(2)    = (/'COOR    ', &
                                                            'COOR+BOX'/)

  ! subroutines
  public  :: alloc_trj_list
  public  :: dealloc_trj_list
  public  :: alloc_trajectory
  public  :: dealloc_trajectory

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_trj_list
  !> @brief        allocate trajectory files
  !! @authors      NT
  !! @param[inout] trj_list : trajectory files information
  !! @param[in]    ntrj     : number of trajectory files
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine alloc_trj_list(trj_list, ntrj)

    ! formal argments
    type(s_trj_list),        intent(inout) :: trj_list
    integer,                 intent(in)    :: ntrj

    ! local variables
    integer                 :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    if (allocated(trj_list%filenames)) then
      if (size(trj_list%filenames) /= ntrj) &
        deallocate(trj_list%filenames,      &
                   trj_list%md_steps,       &
                   trj_list%mdout_periods,  &
                   trj_list%ana_periods,    &
                   stat = dealloc_stat)
    end if

    if (.not. allocated(trj_list%filenames)) &
      allocate(trj_list%filenames    (ntrj), &
               trj_list%md_steps     (ntrj), &
               trj_list%mdout_periods(ntrj), &
               trj_list%ana_periods  (ntrj), &
               trj_list%start_steps  (ntrj), &
               stat = alloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    trj_list%filenames    (1:ntrj) = ''
    trj_list%md_steps     (1:ntrj) = 0
    trj_list%mdout_periods(1:ntrj) = 0
    trj_list%ana_periods  (1:ntrj) = 0
    trj_list%start_steps  (1:ntrj) = 0

    return

  end subroutine alloc_trj_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_trj_list
  !> @brief        deallocate trajectory files
  !! @authors      NT
  !! @param[inout] trj_list : trajectory files information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_trj_list(trj_list)

    ! formal argments
    type(s_trj_list),        intent(inout) :: trj_list

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    if (allocated(trj_list%filenames)) then
      deallocate(trj_list%filenames,     &
                 trj_list%md_steps,      &
                 trj_list%mdout_periods, &
                 trj_list%ana_periods,   &
                 trj_list%start_steps,   &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_trj_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_trajectory
  !> @brief        allocate trajectory
  !! @authors      NT
  !! @param[inout] trajectory : trajectory information
  !! @param[in]    natom      : number of trajectory atom count
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine alloc_trajectory(trajectory, natom)

    ! formal argments
    type(s_trajectory),      intent(inout) :: trajectory
    integer,                 intent(in)    :: natom

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    if (allocated(trajectory%coord)) then
      if (size(trajectory%coord(1,:)) /= natom) &
        deallocate(trajectory%coord, stat = dealloc_stat)
    end if

    if (.not. allocated(trajectory%coord)) &
      allocate(trajectory%coord(3,natom), stat = alloc_stat)

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    trajectory%coord(1:3,1:natom) = 0.0_wp

    return

  end subroutine alloc_trajectory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_trajectory
  !> @brief        deallocate trajectory
  !! @authors      NT
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_trajectory(trajectory)

    ! formal argments
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    if (allocated(trajectory%coord)) then
      deallocate(trajectory%coord, stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_trajectory

end module trajectory_str_mod
