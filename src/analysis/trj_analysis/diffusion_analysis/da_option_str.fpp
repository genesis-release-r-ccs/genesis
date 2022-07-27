!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_option_str_mod
!> @brief   structure of option information
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_option_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option
    logical                            :: check_only
    real(wp)                           :: time_step
    real(wp)                           :: distance_unit
    integer                            :: start_step    = -1
    real(wp)                           :: start_time    = -1_wp
    real(wp)                           :: start_percent = -1_wp
    integer                            :: stop_step     = -1
    real(wp)                           :: stop_time     = -1_wp
    real(wp)                           :: stop_percent  = -1_wp
    integer, dimension(:), allocatable :: ndofs
  end type s_option

  public :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      DS, TM
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option

    return

  end subroutine dealloc_option

end module da_option_str_mod
