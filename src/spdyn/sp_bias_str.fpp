!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_bias_str_mod
!> @brief   structure of bias-exchange replica information
!! @authors Shingo Ito (SI)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_bias_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures

  type, public :: s_bias
    integer                           :: dimension
    integer                           :: nreplicas        
    integer                           :: nwalkers
    integer                           :: ndocks 
    integer                           :: type            
    integer                           :: ncycles         
    integer                           :: exchange_period 
    integer                           :: iseed           
    logical                           :: equilibration_only 
    logical                           :: reduce_bias_potential
    logical                           :: force_exchange ! For debugging

    integer,              allocatable :: dockid2repsetid(:)
    integer,              allocatable :: dockid2repsetid_ref(:)
    integer,              allocatable :: repsetid2dockid(:)
    integer,              allocatable :: repsetid2paramid(:,:)
    integer,              allocatable :: num_criteria_rep(:,:)
    integer,              allocatable :: num_exchanges_rep(:,:)  
    integer,              allocatable :: reac_function(:)
    integer,              allocatable :: umbrid2funcid(:)
    integer,              allocatable :: ifunc(:,:)
    integer,              allocatable :: umbrid2numfuncs(:)
    integer,              allocatable :: umbrid2funclist(:,:)
    real(wp),             allocatable :: dparameters(:)
    real(wp),             allocatable :: potential_energy(:)
    real(wp),             allocatable :: area(:)
    real(wp),             allocatable :: after_gather (:)
    real(wp),             allocatable :: random_table (:,:)
    logical,              allocatable :: cyclic_params(:)
  end type s_bias

  ! parameters for allocatable variables
  integer, public, parameter      :: ABMDBias    = 1

  ! parameters for REMD type
  integer, public, parameter      :: BiasExchange       = 1
  !integer, public, parameter      :: MWBiasExchage      = 2
  !integer, public, parameter      :: PTBiasExchage      = 3

  character(5), public, parameter :: DummyWaterName      = 'ZZZZZ'
  character(*), public, parameter :: BiasTypes(2)  = (/'BIAS     ',&
                                                       'MWBIAS   '/)

  ! subroutines
  public :: alloc_bias
  public :: dealloc_bias
  public :: dealloc_bias_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_bias
  !> @brief        allocate BEABMD information
  !! @authors      SI
  !! @param[inout] bias      : information of bias-exchange
  !! @param[in]    variable  : allocatable variable
  !! @param[in]    var_size1 : size of variables (Total number of docks)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_bias(bias, variable, var_size1)

    ! formal arguments
    type(s_bias),            intent(inout) :: bias
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(ABMDBias)

      if (allocated(bias%cyclic_params)) then
        deallocate(bias%cyclic_params,          &
                   bias%dockid2repsetid,        &
                   bias%dockid2repsetid_ref,    &
                   bias%repsetid2dockid,        &
                   bias%repsetid2paramid,       &
                   bias%num_criteria_rep,       &
                   bias%num_exchanges_rep,      &                 
                   bias%ifunc,                  &
                   bias%dparameters,            &
                   bias%reac_function,          &
                   bias%area,                   &
                   bias%after_gather,           &
                   bias%random_table,           &
                   stat = dealloc_stat)
      end if

      allocate(bias%cyclic_params(var_size1),                    &
               bias%dockid2repsetid(var_size1),                  &
               bias%dockid2repsetid_ref(var_size1),              &
               bias%repsetid2dockid(var_size1),                  &
               bias%repsetid2paramid(var_size1,2),               &
               bias%num_criteria_rep(var_size1,2),               &
               bias%num_exchanges_rep(var_size1,2),              &          
               bias%ifunc(var_size1,4),                          &
               bias%dparameters(var_size1),                      &
               bias%reac_function(var_size1),                    &
               bias%potential_energy(var_size1),                 &
               bias%area(var_size1),                             &
               bias%after_gather(var_size1),                     &
               bias%random_table(var_size1,var_size1),           &
               stat = alloc_stat)

      bias%cyclic_params(1:var_size1)                        = .false.
      bias%dockid2repsetid(1:var_size1)                      = 0
      bias%dockid2repsetid_ref(1:var_size1)                  = 0
      bias%repsetid2dockid(1:var_size1)                      = 0
      bias%repsetid2paramid(1:var_size1,1:2)                 = 0
      bias%num_criteria_rep(1:var_size1,1:2)                 = 0
      bias%num_exchanges_rep(1:var_size1,1:2)                = 0
      bias%ifunc(1:var_size1,1:4)                            = 0
      bias%dparameters(1:var_size1)                          = 0.0_wp
      bias%reac_function(1:var_size1)                        = 0

    case default

      call error_msg('Alloc_Bias> Invalid variable.')

    end select


    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_bias

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_bias
  !> @brief        deallocate bias-exchange information
  !! @authors      SI
  !! @param[inout] bias     : bias-exchange information
  !! @param[in]    variable : allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_bias(bias, variable)

    ! formal arguments
    type(s_bias),            intent(inout) :: bias
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (ABMDBias)

      if (allocated(bias%cyclic_params)) then
        deallocate(bias%cyclic_params,          &
                   bias%dockid2repsetid,        &
                   bias%dockid2repsetid_ref,    &
                   bias%repsetid2dockid,        &
                   bias%repsetid2paramid,       &
                   bias%num_criteria_rep,       &
                   bias%num_exchanges_rep,      &
                   bias%ifunc,                  &
                   bias%dparameters,            &
                   bias%reac_function,          &
                   bias%area,                   &
                   bias%after_gather,           &
                   bias%random_table,           &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Bias> Invalid variable.')

    end select 


    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_bias

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_bias_all
  !> @brief        deallocate all bias-exchange information
  !! @authors      SI
  !! @param[inout] bias : information of bias-exchange
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_bias_all(bias)

    ! formal arguments
    type(s_bias),         intent(inout) :: bias


    call dealloc_bias(bias, ABMDBias)

    return

  end subroutine dealloc_bias_all

end module sp_bias_str_mod
