!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_dynamics_str_mod
!> @brief   structure of dynamics information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_dynamics_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_dynamics
    logical             :: restart
    logical             :: multiple_random_seed
    integer             :: integrator
    integer             :: nsteps
    real(wp)            :: timestep
    integer             :: istart_step
    integer             :: iend_step
    integer             :: eneout_period
    integer             :: crdout_period
    integer             :: velout_period
    integer             :: rstout_period
    integer             :: stoptr_period
    integer             :: nbupdate_period
    integer             :: qmsave_period
    integer             :: iseed
    real(wp)            :: initial_time
    logical             :: stop_com_translation
    logical             :: stop_com_rotation
    logical             :: annealing
    integer             :: anneal_period
    real(wp)            :: dtemperature
    logical             :: target_md
    logical             :: steered_md
    real(wp)            :: initial_value
    real(wp)            :: final_value
    logical             :: verbose
    logical             :: random_restart
    logical             :: shrink_box
    integer             :: shrink_period
    real(wp)            :: dbox_x
    real(wp)            :: dbox_y
    real(wp)            :: dbox_z
    logical             :: esp_mm
    integer             :: calc_qm_period
    logical             :: avg_qm_charge
  end type s_dynamics

  ! parameters
  integer,      public, parameter :: IntegratorLEAP = 1
  integer,      public, parameter :: IntegratorVVER = 2
  integer,      public, parameter :: IntegratorBDEM = 3
  integer,      public, parameter :: IntegratorBD2N = 4
  integer,      public, parameter :: IntegratorSDMP = 5
  integer,      public, parameter :: IntegratorVVER_CG = 6

  character(*), public, parameter :: IntegratorTypes(6)  = &
                                      (/'LEAP   ','VVER   ','BDEM   ',&
                                        'BD2N   ','SDMP   ','VVER_CG'/)


  ! subroutines
  public  ::  init_dynamics

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_dynamics
  !> @brief        initialize dynamics information
  !! @authors      TM
  !! @param[out]   dynamics : dynamics information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_dynamics(dynamics)

    ! formal arguments
    type(s_dynamics), intent(inout) :: dynamics


    dynamics%restart              = .false.
    dynamics%integrator           = 0
    dynamics%nsteps               = 0
    dynamics%timestep             = 0.0_wp
    dynamics%eneout_period        = 0
    dynamics%crdout_period        = 0
    dynamics%velout_period        = 0
    dynamics%rstout_period        = 0
    dynamics%stoptr_period        = 0
    dynamics%nbupdate_period      = 0
    dynamics%iseed                = 0
    dynamics%initial_time         = 0.0_wp
    dynamics%stop_com_translation = .false.
    dynamics%stop_com_rotation    = .false.
    dynamics%annealing            = .false.
    dynamics%anneal_period        = 0
    dynamics%dtemperature         = 0.0_wp
    dynamics%istart_step          = 0
    dynamics%iend_step            = 0
    dynamics%verbose              = .false.
    dynamics%target_md            = .false.
    dynamics%steered_md           = .false.
    dynamics%initial_value        = 0.0_wp
    dynamics%final_value          = 0.0_wp
    dynamics%random_restart       = .true.

    return

  end subroutine init_dynamics

end module at_dynamics_str_mod
