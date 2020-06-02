!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_minimize_str_mod
!> @brief   structure of energy minimization
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_minimize_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_minimize
    integer             :: method 
    integer             :: nsteps
    integer             :: eneout_period
    integer             :: crdout_period
    integer             :: rstout_period
    integer             :: nbupdate_period 
    real(wp)            :: force_scale_init
    real(wp)            :: force_scale_max 
    logical             :: verbose
  end type s_minimize

  ! parameters
  integer,      public, parameter :: MinimizeMethodSD = 1
  
  character(*), public, parameter :: MinimizeMethodTypes(1) = (/'SD'/)


  ! subroutines
  public  ::  init_minimize

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_minimize
  !> @brief        initialize minimize information
  !! @authors      TM
  !! @param[out]   minimize : minimize information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_minimize(minimize)

    ! formal arguments
    type(s_minimize),  intent(inout) :: minimize


    minimize%method          = 0
    minimize%nsteps          = 0
    minimize%eneout_period   = 0
    minimize%crdout_period   = 0
    minimize%rstout_period   = 0
    minimize%nbupdate_period = 0
    minimize%verbose         = .false.

    return

  end subroutine init_minimize

end module sp_minimize_str_mod
