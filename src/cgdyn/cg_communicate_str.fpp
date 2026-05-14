!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_communicate_mod
!> @brief   Provides utilities for MPI communication.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_communicate_str_mod

  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_comm
    integer,          allocatable :: num_cell_coord(:)
    integer,          allocatable :: num_cell_force(:)
    integer,          allocatable :: send_force_size(:)
    integer,          allocatable :: recv_force_size(:)
    integer,          allocatable :: send_coord_size(:)
    integer,          allocatable :: recv_coord_size(:)
    integer,          allocatable :: irequest(:,:)
    integer,          allocatable :: send_size(:,:)
    integer,          allocatable :: recv_size(:,:)
    integer,          allocatable :: send_address(:,:)
    integer,          allocatable :: recv_address(:,:)

    integer,          allocatable :: ic_send(:,:)
    integer,          allocatable :: ic_recv(:,:)
    integer,          allocatable :: if_send(:,:)
    integer,          allocatable :: if_recv(:,:)
    integer,          allocatable :: ic_send_list(:,:)
    integer,          allocatable :: ic_recv_list(:,:)
    integer,          allocatable :: if_send_list(:,:)
    integer,          allocatable :: if_recv_list(:,:)
    integer,          allocatable :: int_send(:)
    integer,          allocatable :: int_recv(:)
    real(wip),        allocatable :: buf_send(:)
    real(wip),        allocatable :: buf_recv(:)
  end type s_comm

  ! parameters for allocatable variables
  integer,      public, parameter :: CommProc     = 1
  integer,      public, parameter :: CommCell     = 2
  integer,      public, parameter :: CommBuffer   = 3

  ! subroutines
  public  :: alloc_comm
  public  :: dealloc_comm
  public  :: dealloc_comm_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_comm
  !> @brief        Allocates the communication structure.
  !! @authors      JJ
  !! @param[inout] comm      : communicate structure
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_comm(comm, variable, var_size, var_size1, var_size2)

    ! formal arguments
    type(s_comm),            intent(inout) :: comm
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (CommProc)

      if (allocated(comm%num_cell_coord)) then
        if (size(comm%num_cell_coord(:)) /= var_size) &
          deallocate(comm%num_cell_coord,   &
                     comm%num_cell_force,   &
                     comm%recv_force_size,  &
                     comm%send_force_size,  &
                     comm%recv_coord_size,  &
                     comm%send_coord_size,  &
                     comm%irequest,         &
                     comm%send_size,        &
                     comm%recv_size,        &
                     comm%send_address,     &
                     comm%recv_address,     &
                     stat = dealloc_stat)
      end if

      if (.not.allocated(comm%num_cell_coord))    &
        allocate(comm%num_cell_coord (var_size),  &
                 comm%num_cell_force (var_size),  &
                 comm%recv_force_size(var_size),  &
                 comm%send_force_size(var_size),  &
                 comm%recv_coord_size(var_size),  &
                 comm%send_coord_size(var_size),  &
                 comm%irequest     (4,var_size),  &
                 comm%send_size    (2,var_size),  &
                 comm%recv_size    (2,var_size),  &
                 comm%send_address(2,var_size+1), &
                 comm%recv_address(2,var_size+1), &
                 stat = alloc_stat)

    case (CommCell)

      if (allocated(comm%ic_send)) then
        deallocate(comm%ic_send,      &
                   comm%ic_recv,      &
                   comm%if_send,      &
                   comm%if_recv,      &
                   comm%ic_send_list, &
                   comm%ic_recv_list, &
                   comm%if_send_list, &
                   comm%if_recv_list, &
                   stat = dealloc_stat)
      end if

      if (.not.allocated(comm%ic_send))    &
        allocate(comm%ic_send     (var_size1, var_size), &
                 comm%ic_recv     (var_size1, var_size), &
                 comm%if_send     (var_size1, var_size), &
                 comm%if_recv     (var_size1, var_size), &
                 comm%ic_send_list(var_size2, var_size), &
                 comm%ic_recv_list(var_size2, var_size), &
                 comm%if_send_list(var_size2, var_size), &
                 comm%if_recv_list(var_size2, var_size), &
                 stat = alloc_stat)

    case (CommBuffer)

      if (allocated(comm%buf_send)) then
        if (size(comm%buf_send(:)) /= var_size) &
          deallocate(comm%buf_send,             &
                     comm%buf_recv,             &
                     comm%int_send,             &
                     comm%int_recv,             &
                     stat = dealloc_stat)
      end if

      if (.not.allocated(comm%buf_send))    &
        allocate(comm%buf_send(var_size),   &
                 comm%int_send(var_size),   &
                 comm%buf_recv(var_size1),  &
                 comm%int_recv(var_size1),  &
                 stat = alloc_stat)

    case default

      call error_msg('Error: Alloc_Comm - Invalid variable specified.')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_comm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_comm
  !> @brief        Deallocates the communication structure.
  !! @authors      JJ
  !! @param[inout] comm      : communicate structure
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_comm(comm, variable)

    ! formal arguments
    type(s_comm),            intent(inout) :: comm
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat

    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (CommProc)

      if (allocated(comm%num_cell_coord)) then
        deallocate(comm%num_cell_coord,   &
                   comm%num_cell_force,   &
                   comm%recv_force_size,  &
                   comm%send_force_size,  &
                   comm%recv_coord_size,  &
                   comm%send_coord_size,  &
                   comm%irequest,         &
                   comm%send_size,        &
                   comm%recv_size,        &
                   comm%send_address,     &
                   comm%recv_address,     &
                   stat = dealloc_stat)
    end if

    case (CommCell)

      if (allocated(comm%ic_send)) then
        deallocate(comm%ic_send,      &
                   comm%ic_recv,      &
                   comm%if_send,      &
                   comm%if_recv,      &
                   comm%ic_send_list, &
                   comm%ic_recv_list, &
                   comm%if_send_list, &
                   comm%if_recv_list, &
                   stat = dealloc_stat)
      end if

    case (CommBuffer)

      if (allocated(comm%buf_send)) then
        deallocate(comm%buf_send,             &
                   comm%buf_recv,             &
                   comm%int_send,             &
                   comm%int_recv,             &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Error: Dealloc_Comm - Invalid variable specified.')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_comm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
    !  Subroutine    dealloc_comm_all
  !> @brief        Deallocates all communication structures.
  !! @authors      JJ
  !! @param[inout] comm : communication structure  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_comm_all(comm)

    ! formal arguments
    type(s_comm),          intent(inout) :: comm

    call dealloc_comm(comm, CommProc)
    call dealloc_comm(comm, CommCell)
    call dealloc_comm(comm, CommBuffer)

    return

  end subroutine dealloc_comm_all

end module cg_communicate_str_mod


