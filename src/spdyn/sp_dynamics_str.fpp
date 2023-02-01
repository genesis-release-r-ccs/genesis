!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_dynamics_str_mod
!> @brief   structure of dynamics information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_dynamics_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_dynamics
    logical             :: restart
    integer             :: integrator
    integer             :: nsteps
    real(wip)           :: timestep
    integer             :: istart_step
    integer             :: iend_step
    integer             :: crdout_period
    integer             :: velout_period
    integer             :: eneout_period
    integer             :: rstout_period
    integer             :: stoptr_period
    integer             :: nbupdate_period
    integer             :: elec_long_period
    integer             :: iseed
    integer             :: iseed_init_velocity
    logical             :: iseed_read
    real(wip)           :: initial_time
    logical             :: stop_com_translation
    logical             :: stop_com_rotation
    logical             :: annealing
    integer             :: anneal_period
    real(wip)           :: dtemperature
    logical             :: xi_respa
    logical             :: xo_respa
    integer             :: thermo_period
    integer             :: baro_period
    logical             :: verbose
    logical             :: target_md
    logical             :: steered_md
    logical             :: hydrogen_mr
    integer             :: hmr_target
    real(wp)            :: hmr_ratio
    real(wp)            :: hmr_ratio_xh1
    real(wp)            :: hmr_ratio_xh2
    real(wp)            :: hmr_ratio_xh3
    real(wp)            :: hmr_ratio_ring
    real(wp)            :: initial_rmsd
    real(wp)            :: final_rmsd
    ! FEP
    integer             :: fepout_period
    integer             :: equilsteps
  end type s_dynamics

  ! parameters
  integer,      public, parameter :: IntegratorLEAP = 1
  integer,      public, parameter :: IntegratorVVER = 2
  integer,      public, parameter :: IntegratorVRES = 3
  integer,      public, parameter :: IntegratorPMTS = 4

  character(*), public, parameter :: IntegratorTypes(4)  = (/'LEAP', &
                                                             'VVER', &
                                                             'VRES', &
                                                             'PTMS'/)

  integer,      public, parameter :: HmrTargetAll     = 1
  integer,      public, parameter :: HmrTargetSolvent = 2
  integer,      public, parameter :: HmrTargetSolute  = 3

  character(*), public, parameter :: HmrTargetTypes(3)  = (/'ALL    ', &
                                                            'SOLVENT', &
                                                            'SOLUTE '/)

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
    type(s_dynamics),        intent(inout) :: dynamics


    dynamics%restart              = .false.
    dynamics%integrator           = 0
    dynamics%timestep             = 0.0_dp 
    dynamics%nsteps               = 0
    dynamics%initial_time         = 0.0_wip
    dynamics%crdout_period        = 0
    dynamics%velout_period        = 0
    dynamics%eneout_period        = 0
    dynamics%rstout_period        = 0
    dynamics%stoptr_period        = 0
    dynamics%nbupdate_period      = 0
    dynamics%elec_long_period     = 0
    dynamics%iseed                = 0
    dynamics%iseed_init_velocity  = 0
    dynamics%stop_com_translation = .false.
    dynamics%stop_com_rotation    = .false.
    dynamics%annealing            = .false.
    dynamics%anneal_period        = 0
    dynamics%dtemperature         = 0.0_wip
    dynamics%xi_respa             = .true.
    dynamics%xo_respa             = .false.
    dynamics%thermo_period        = 0
    dynamics%baro_period          = 0
    dynamics%istart_step          = 0
    dynamics%iend_step            = 0
    dynamics%verbose              = .false.
    dynamics%iseed_read           = .false.
    dynamics%hydrogen_mr          = .false.
    dynamics%hmr_ratio            = 3.0_wp
    ! FEP
    dynamics%fepout_period        = 0
    dynamics%equilsteps           = 0

    return

  end subroutine init_dynamics

end module sp_dynamics_str_mod
