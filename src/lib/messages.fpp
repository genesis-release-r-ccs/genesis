!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   messages_mod
!> @brief   utilities for taking messages to Std or Err output
!! @authors Yuji Sugita (YS)
!! @note    MsgOut and ErrOut are defined here.
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module messages_mod

  use mpi_parallel_mod

  implicit none
  private

  ! parameters
  integer, public, parameter :: MsgOut  = 6
  integer, public, parameter :: ErrOut  = 6

  ! subroutines and functions
  public :: error_msg
  public :: error_msg_alloc
  public :: error_msg_dealloc
  public :: error_msg_fileio

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg
  !> @brief        write error message and stop execution
  !! @param[in]    message     : Error message (optional)
  !! @param[in]    file_name   : File name calling subroutine (optional)
  !! @param[in]    line_number : Line number where calling happens (optional)
  !! @authors      YS,DU
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg(message, file_name, line_number)

#if defined(INTEL)
    use ifcore, only:tracebackqq

#elif defined(KCOMP)
    use service_routines, only:errtra

#endif

    ! formal arguments
    character(*), optional, intent(in) :: message
    character(*), optional, intent(in) :: file_name
    integer,      optional, intent(in) :: line_number

    ! local variables
    character(1024) :: msg
    character(1024) :: loc
    character(1024) :: buf

    ! Check optional variables
    msg = ""
    if (present(message)) then
      write(msg, "(A)") message
    end if

    loc = ""
    if (present(file_name) .and. present(line_number)) then
      write(loc, "(A,A1,I05)") trim(file_name), ":", line_number
    end if

    ! Write buffer
    if (len(trim(loc)) == 0) then
      write(buf, "(A,A13,I5)") trim(msg), " at rank_no: ", my_world_rank
    else
      write(buf, "(A1,A,A2,A,A13,I5)") "(", trim(loc), ") ", trim(msg), " at rank_no: ", my_world_rank
    end if

#if defined(INTEL)
    call TracebackQQ(trim(msg),-1)
    call exit(1)

#elif defined(KCOMP)
    write(ErrOut,*) trim(msg)
    call errtra
    call setrcd(1)
    call exit

#else
    write(ErrOut,*) trim(msg)
    call exit(1)

#endif

  end subroutine error_msg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_alloc
  !> @brief        write error message in allocation and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_alloc

    call error_msg('Memory allocation error')

  end subroutine error_msg_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_dealloc
  !> @brief        write error message in deallocation and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_dealloc

    call error_msg('Memory deallocation error')

  end subroutine error_msg_dealloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_fileio
  !> @brief        write error message in fileio and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_fileio

    call error_msg('File I/O error')

  end subroutine error_msg_fileio

end module messages_mod
