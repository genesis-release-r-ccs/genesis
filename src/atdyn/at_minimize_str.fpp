!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_minimize_str_mod
!> @brief   structure of energy minimization
!! @authors Takaharu Mori (TM), Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_minimize_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_minimize
    integer                :: method 
    integer                :: nsteps
    integer                :: eneout_period
    integer                :: crdout_period
    integer                :: rstout_period
    integer                :: nbupdate_period 
    logical                :: verbose
    integer                :: qmsave_period
    integer, allocatable   :: fixedatom_id(:)
    integer                :: num_fixatoms
    integer, allocatable   :: optatom_id(:)
    integer                :: num_optatoms
    real(wp)               :: tol_rmsg
    real(wp)               :: tol_maxg
    logical                :: eneout
    logical                :: crdout
    logical                :: rstout
    logical                :: eneout_short

    ! For SD
    real(wp)               :: force_scale_init
    real(wp)               :: force_scale_max 

    ! For LBFGS
    integer                :: ncorrection
    logical                :: lbfgs_bnd
    logical                :: lbfgs_bnd_qmonly
    real(wp)               :: lbfgs_bnd_maxmove

    ! For LBFGS - micro-iteration
    logical                :: macro
    integer                :: start_micro
    integer                :: nsteps_micro
    real(wp)               :: tol_rmsg_micro
    real(wp)               :: tol_maxg_micro
    integer, allocatable   :: optatom_macro_id(:)
    integer                :: num_optatoms_macro
    integer, allocatable   :: optatom_micro_id(:)
    integer                :: num_optatoms_micro

    character(len=60)      :: csave_micro
    logical                :: lsave_micro(4)
    integer                :: isave_micro(44)
    real(wp)               :: dsave_micro(29)
    real(wp), allocatable  :: vec_micro(:)
    real(wp), allocatable  :: upper_micro(:)
    real(wp), allocatable  :: lower_micro(:)
    real(wp), allocatable  :: gradient_micro(:)
    real(wp), allocatable  :: work_lbfgs_micro(:)
    integer,  allocatable  :: list_bound_micro(:)
    integer,  allocatable  :: iwork_lbfgs_micro(:)

    logical                :: check_structure

!ky what is this for?
    integer                :: pairalloc_period

  end type s_minimize

  ! parameters
  integer,      public, parameter :: MinimizeMethodSD    = 1
  integer,      public, parameter :: MinimizeMethodLBFGS = 2
  
  character(*), public, parameter :: MinimizeMethodTypes(2) = (/'SD   ',&
                                                                'LBFGS'/)

  integer, public, parameter :: MinimizeLBFGS      = 1

  ! subroutines
  public  ::  init_minimize
  public  ::  alloc_minimize
  public  ::  dealloc_minimize

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


    minimize%method           = 0
    minimize%nsteps           = 0
    minimize%eneout_period    = 0
    minimize%crdout_period    = 0
    minimize%rstout_period    = 0
    minimize%nbupdate_period  = 0
    minimize%pairalloc_period = 0
    minimize%verbose          = .false.
    minimize%macro            = .false.
    minimize%nsteps_micro     = 0
    minimize%lbfgs_bnd         = .true.
    minimize%lbfgs_bnd_qmonly  = .true.
    minimize%lbfgs_bnd_maxmove = 0.1D+00
    minimize%check_structure   = .true.

    return

  end subroutine init_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_minimize
  !> @brief        allocate minimize information
  !! @authors      KY
  !! @param[inout] minimize  : structure of restraints information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size2 : size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_minimize(minimize, variable, var_size, var_size2)

    ! formal arguments
    type(s_minimize),        intent(inout) :: minimize
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,      optional,  intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (MinimizeLBFGS)

      if (allocated(minimize%vec_micro)) then
        if (size(minimize%vec_micro(:)) == var_size) return
        deallocate( &
                   minimize%vec_micro,            &
                   minimize%upper_micro,          &
                   minimize%lower_micro,          &
                   minimize%gradient_micro,       &
                   minimize%work_lbfgs_micro,     &
                   minimize%list_bound_micro,     &
                   minimize%iwork_lbfgs_micro,    &
                   stat = dealloc_stat)
      end if

      allocate( &
               minimize%vec_micro(var_size),           &
               minimize%upper_micro(var_size),         &
               minimize%lower_micro(var_size),         &
               minimize%gradient_micro(var_size),      &
               minimize%work_lbfgs_micro(var_size2),   &
               minimize%iwork_lbfgs_micro(var_size*3), &
               minimize%list_bound_micro(var_size),    &
               stat = alloc_stat)

      minimize%vec_micro        (1:var_size) = 0.0_wp
      minimize%upper_micro      (1:var_size) = 0.0_wp
      minimize%lower_micro      (1:var_size) = 0.0_wp
      minimize%gradient_micro   (1:var_size) = 0.0_wp
      minimize%work_lbfgs_micro (1:var_size2) = 0.0_wp
      minimize%iwork_lbfgs_micro(1:var_size*3) = 0
      minimize%list_bound_micro (1:var_size) = 0

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_minimize
  !> @brief        deallocate minimize information
  !! @authors      KY
  !! @param[inout] minimize : structure of restraints information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_minimize(minimize, variable)

    ! formal arguments
    type(s_minimize),        intent(inout) :: minimize
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    select case (variable)

    case (MinimizeLBFGS)

      if (allocated(minimize%vec_micro)) then
        deallocate( &
                   minimize%vec_micro,            &
                   minimize%upper_micro,          &
                   minimize%lower_micro,          &
                   minimize%gradient_micro,       &
                   minimize%work_lbfgs_micro,     &
                   minimize%list_bound_micro,     &
                   minimize%iwork_lbfgs_micro,    &
                   stat = dealloc_stat)
      end if

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_minimize

end module at_minimize_str_mod
