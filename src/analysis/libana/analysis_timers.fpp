!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   analysis_timers_mod
!> @brief   timers_mod routines adapted to analysis tools
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module analysis_timers_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  ! import only variables and routines that are needed
  use timers_mod, only : TimerOn
  use timers_mod, only : TimerOff
  use timers_mod, only : TimerTotal
  use timers_mod, only : TimerDynamics
  use timers_mod, only : timer
!TM(191230)
!  use timers_mod, only : calculate_time
!TM(191230)

  implicit none
  private

  ! parameters from timers_mod to be accessible by using analysis_timers_mod
  public                     :: TimerOn
  public                     :: TimerOff
  public                     :: TimerTotal
  public                     :: TimerDynamics

  ! timers_mod timer routine will be accessible by using analysis_timers_mod
  public                     :: timer

  integer, public, parameter :: TimerAnalysis = TimerDynamics


  ! output_time routine from this module will be used instead of timers_mods
  public :: output_time

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_time
  !> @brief        output analysis timer profile
  !! @authors      DS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_time

    ! local variables
    real(dp)                      :: avetime_setup
    real(dp)        , allocatable :: avetime(:)
    real(dp)        , allocatable :: maxtime(:), mintime(:)


!TM(191230)
!    call calculate_time(avetime, mintime, maxtime)
    return
!TM(191230)

    if (main_rank) then
      avetime_setup   = avetime(TimerTotal) - avetime(TimerAnalysis)

      write(MsgOut,'(a)') 'Output_Time> Analysis timer profile'
      write(MsgOut,'(a,f12.3)') '  total time      =', avetime(TimerTotal)
      write(MsgOut,'(a,f12.3)') '    setup         =', avetime_setup
      write(MsgOut,'(a,f12.3)') '    analysis      =', avetime(TimerAnalysis)
      write(MsgOut,'(a)') ''

    end if

    return

  end subroutine output_time

end module analysis_timers_mod
