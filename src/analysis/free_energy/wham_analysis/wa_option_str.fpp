!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   wa_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module wa_option_str_mod

  use select_atoms_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only

    ! wham variables
    integer                       :: dimension
    integer                       :: nblocks
    real(wp)                      :: temperature
    real(wp)                      :: tolerance
    integer,          allocatable :: rest_func_num(:)
    integer,          allocatable :: rest_func_no(:)
    real(wp),         allocatable :: grid_min(:)
    real(wp),         allocatable :: grid_max(:)
    integer,          allocatable :: num_grids(:)
    logical,          allocatable :: is_periodic(:)
    real(wp),         allocatable :: box_size(:)

    ! selection variables
    integer                       :: num_atoms
    type(s_selatoms), allocatable :: selatoms(:)

    ! restraints variables
    integer,          allocatable :: rest_funcs(:)
    integer,          allocatable :: rest_sel_index(:,:)
    integer,          allocatable :: rest_nreplica(:)
    real(wp),         allocatable :: rest_constants(:,:)
    real(wp),         allocatable :: rest_references(:,:)

  end type s_option

  ! parameters for RESTRAINTS functions
  integer,      public, parameter :: RestraintsFuncPOSI      = 1
  integer,      public, parameter :: RestraintsFuncDIST      = 2
  integer,      public, parameter :: RestraintsFuncDISTCOM   = 3
  integer,      public, parameter :: RestraintsFuncRMSD      = 4
  integer,      public, parameter :: RestraintsFuncRMSDCOM   = 5
  integer,      public, parameter :: RestraintsFuncANGLE     = 6
  integer,      public, parameter :: RestraintsFuncANGLECOM  = 7
  integer,      public, parameter :: RestraintsFuncDIHED     = 8
  integer,      public, parameter :: RestraintsFuncDIHEDCOM  = 9
  
  character(*), public, parameter :: RestraintsFuncs(9) = (/'POSI     ', &
                                                            'DIST     ', &
                                                            'DISTMASS ', &
                                                            'RMSD     ', &
                                                            'RMSDMASS ', &
                                                            'ANGLE    ', &
                                                            'ANGLEMASS', &
                                                            'DIHED    ', &
                                                            'DIHEDMASS'/)

  ! subroutines
  public  :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      NT
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, dealloc_stat


    if (allocated(option%rest_func_num)) then
      deallocate(option%rest_func_num,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%rest_func_no)) then
      deallocate(option%rest_func_no,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%grid_min)) then
      deallocate(option%grid_min,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%grid_max)) then
      deallocate(option%grid_max,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%num_grids)) then
      deallocate(option%num_grids,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%is_periodic)) then
      deallocate(option%is_periodic,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%box_size)) then
      deallocate(option%box_size,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%selatoms)) then
      do i = 1, size(option%selatoms)
        call dealloc_selatoms(option%selatoms(i))
      end do
      deallocate(option%selatoms,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%rest_funcs)) then
      deallocate(option%rest_funcs,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%rest_sel_index)) then
      deallocate(option%rest_sel_index,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%rest_constants)) then
      deallocate(option%rest_constants,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    if (allocated(option%rest_references)) then
      deallocate(option%rest_references,stat=dealloc_stat)
      if (dealloc_stat /= 0) &
        call error_msg_dealloc
    end if

    return

  end subroutine dealloc_option

end module wa_option_str_mod
