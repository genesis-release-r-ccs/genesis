!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_rpath_str_mod
!> @brief   structure of replica path information
!! @authors Yasuaki Komuro (YK), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_rpath_str_mod

  use sp_output_str_mod
  use molecules_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rpath
    integer                           :: ncycle
    integer                           :: dimension
    integer                           :: nreplica
    integer                           :: rpath_period
    real(wp)                          :: delta
    real(wp)                          :: smooth
    real(wp)                          :: sum_distance
    real(wp)                          :: distance_prev
    real(wp)                          :: distance_init
    logical                           :: fix_terminal
    logical                           :: use_restart
    logical                           :: equilibration_only
    logical                           :: avoid_shrinkage
    integer                           :: fitting_method
    character(256)                    :: fitting_atom
    integer,              allocatable :: rest_function(:)
    real(wp),             allocatable :: rest_constants(:,:,:)
    real(wp),             allocatable :: rest_reference(:,:,:)
    real(wp),             allocatable :: rest_reference_prev(:,:)
    real(wp),             allocatable :: rest_reference_init(:,:)
    real(dp),             allocatable :: force(:)
    real(dp),             allocatable :: metric(:,:)
    real(dp),             allocatable :: before_gather(:)
    real(dp),             allocatable :: after_gather(:)
    type(s_output)                    :: output
  end type s_rpath

  ! parameters for allocatable variables
  integer, public, parameter      :: RpathReplicas  = 1
  integer, public, parameter      :: RpathUmbrellas = 2

  ! subroutines
  public :: alloc_rpath
  public :: dealloc_rpath
  public :: dealloc_rpath_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rpath
  !> @brief        allocate rpath information
  !! @authors      TM
  !! @param[inout] rpath     : information of rpath
  !! @param[in]    variable  : allocatable variable
  !! @param[in]    var_size1 : size of variables (REMD dimension)
  !! @param[in]    var_size2 : size of variables (Total number of replicas)
  !! @param[in]    var_size3 : size of variables (max number of replica)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rpath(rpath, variable, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,      optional,  intent(in)    :: var_size2
    integer,      optional,  intent(in)    :: var_size3

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RpathReplicas)
      if (allocated(rpath%rest_function)) then
        if (size(rpath%rest_function(:)) == var_size1) return
        deallocate(rpath%rest_function,       &
                   rpath%force,               &
                   rpath%metric,              &
                   rpath%before_gather,       &
                   rpath%after_gather,        &
                   stat = dealloc_stat)
      end if

      allocate(rpath%rest_function(var_size1),             &
               rpath%force(var_size1),                     &
               rpath%metric(var_size1,var_size1),          &
               rpath%before_gather(var_size1),             &
               rpath%after_gather(var_size1*var_size2),    &
               stat = alloc_stat)

      rpath%rest_function(1:var_size1)                 = 0
      rpath%force(1:var_size1)                         = 0.0_dp
      rpath%metric(1:var_size1,1:var_size1)            = 0.0_dp
      rpath%before_gather(1:var_size1)                 = 0.0_dp
      rpath%after_gather(1:var_size1*var_size2)        = 0.0_dp

    case(RpathUmbrellas)

      if (allocated(rpath%rest_constants)) then
        if (size(rpath%rest_constants(:,:,:)) == (4*var_size1*var_size2)) return
        deallocate(rpath%rest_constants,          &
                   rpath%rest_reference,          &
                   rpath%rest_reference_prev,     &
                   rpath%rest_reference_init,     &
                   stat = dealloc_stat)
      end if

      allocate(rpath%rest_constants(4,var_size1,var_size2), &
               rpath%rest_reference(2,var_size1,var_size2), &
               rpath%rest_reference_prev(var_size1,var_size2), &
               rpath%rest_reference_init(var_size1,var_size2), &
               stat = alloc_stat)

      rpath%rest_constants(1:4,1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference(1:2,1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference_prev(1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference_init(1:var_size1,1:var_size2) = 0.0_wp

    case default

      call error_msg('Alloc_Rpath> bad variable')

    end select

    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath
  !> @brief        deallocate rpath information
  !! @authors      TM
  !! @param[inout] rpath     : rpath information
  !! @param[in]    variable : allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath(rpath, variable)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (RpathReplicas)

      if (allocated(rpath%rest_function)) then
        deallocate(rpath%rest_function,       &
                   rpath%force,               &
                   rpath%metric,              &
                   rpath%before_gather,       &
                   rpath%after_gather,        &
                   stat = dealloc_stat)
      end if

    case (RpathUmbrellas)

      if (allocated(rpath%rest_constants)) then
        deallocate(rpath%rest_constants,       &
                   rpath%rest_reference,      &
                   rpath%rest_reference_prev, &
                   rpath%rest_reference_init, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Rpath> bad variable')

    end select 


    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath_all
  !> @brief        deallocate all rpath information
  !! @authors      TM
  !! @param[inout] rpath : information of rpath
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath_all(rpath)

    ! formal arguments
    type(s_rpath),         intent(inout) :: rpath


    call dealloc_rpath(rpath, RpathReplicas)
    call dealloc_rpath(rpath, RpathUmbrellas)

    return

  end subroutine dealloc_rpath_all

end module sp_rpath_str_mod
