!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_alchemy_str
!> @brief   structure of alchemy information
!! @authors Hiraku Oshima (HO), Nobuhiko Kato (NK)
!
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_alchemy_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_alchemy
    integer               :: equilsteps
    integer               :: fepout_period
    integer               :: num_fep_windows
    real(wp)              :: sc_alpha
    real(wp)              :: sc_beta
    real(wp), allocatable :: lambljA(:)
    real(wp), allocatable :: lambljB(:)
    real(wp), allocatable :: lambelA(:)
    real(wp), allocatable :: lambelB(:)
    real(wp), allocatable :: lambbondA(:)
    real(wp), allocatable :: lambbondB(:)
    real(wp), allocatable :: lambrest(:)
    integer               :: num_fep_neighbor
    integer               :: fep_direction
    integer               :: fep_topology
    integer               :: fep_md_type
    integer               :: ref_lambid
    integer               :: lambid

    ! For REST2-like
    integer, allocatable  :: atom_cls_no_org_singleA(:)
    integer, allocatable  :: atom_cls_no_org_singleB(:)
    integer, allocatable  :: atom_cls_no_org_dualA(:)
    integer, allocatable  :: atom_cls_no_org_dualB(:)
    integer               :: num_atom_cls_single
    integer               :: num_atom_cls_dualA
    integer               :: num_atom_cls_dualB
    integer               :: istart_atom_cls_single
    integer               :: istart_atom_cls_dualA
    integer               :: istart_atom_cls_dualB
  end type s_alchemy

  ! parameters
  integer,      public, parameter :: FEP_Bothsides   = 1
  integer,      public, parameter :: FEP_Forward     = 2
  integer,      public, parameter :: FEP_Reverse     = 3
  integer,      public, parameter :: FEP_Nodirection = 4
  character(*), public, parameter :: FEPDirectionTypes(4) = (/'Bothsides', &
                                                              'Forward  ', &
                                                              'Reverse  ', &
                                                              'NONE     '/)

  integer,      public, parameter :: FEP_Hybrid = 1
  integer,      public, parameter :: FEP_Dual   = 2
  character(*), public, parameter :: FEPTopologyTypes(2) = (/'Hybrid', &
                                                             'Dual  '/)

  integer,      public, parameter :: FEP_Serial   = 1
  integer,      public, parameter :: FEP_Parallel = 2
  integer,      public, parameter :: FEP_Single   = 3
  character(*), public, parameter :: FEPMDTypes(3) = (/'Serial  ', &
                                                       'Parallel', &
                                                       'Single  '/)

  ! subroutines
  public  ::  init_alchemy
  public  ::  alloc_alchemy
  public  ::  dealloc_alchemy

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_alchemy
  !> @brief        initialize alchemy information
  !! @authors      HO
  !! @param[out]   alchemy : alchemy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_alchemy(alchemy)

    ! formal arguments
    type(s_alchemy),        intent(inout) :: alchemy

    alchemy%fepout_period    = 0
    alchemy%equilsteps       = 0
    alchemy%num_fep_windows  = 0
    alchemy%sc_alpha         = 5.0
    alchemy%sc_beta          = 0.5
    alchemy%fep_direction    = FEP_Bothsides
    alchemy%num_fep_neighbor = 2
    alchemy%lambid           = 1

    return

  end subroutine init_alchemy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_alchemy
  !> @brief        allocate md information for FEP
  !! @authors      HO
  !! @param[inout] alchemy  : information of alchemy
  !! @param[in]    var_size1 : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_alchemy(alchemy, var_size1)

    ! formal arguments
    type(s_alchemy),        intent(inout) :: alchemy
    integer,                 intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat, dealloc_stat

    alloc_stat = 0
    dealloc_stat = 0

    if (allocated(alchemy%lambljA)) then
        deallocate(alchemy%lambljA, &
            alchemy%lambljB,        &
            alchemy%lambelA,        &
            alchemy%lambelB,        &
            alchemy%lambbondA,        &
            alchemy%lambbondB,        &
            alchemy%lambrest,       &
            stat = dealloc_stat)
    end if

    allocate(alchemy%lambljA(var_size1), &
        alchemy%lambljB(var_size1),      &
        alchemy%lambelA(var_size1),      &
        alchemy%lambelB(var_size1),      &
        alchemy%lambbondA(var_size1),      &
        alchemy%lambbondB(var_size1),      &
        alchemy%lambrest(var_size1),     &
        stat = alloc_stat)

    alchemy%lambljA(1:var_size1) = 0.0_wp
    alchemy%lambljB(1:var_size1) = 0.0_wp
    alchemy%lambelA(1:var_size1) = 0.0_wp
    alchemy%lambelB(1:var_size1) = 0.0_wp
    alchemy%lambbondA(1:var_size1) = 0.0_wp
    alchemy%lambbondB(1:var_size1) = 0.0_wp
    alchemy%lambrest(1:var_size1) = 0.0_wp

    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_alchemy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_alchemy
  !> @brief        deallocate alchemy information for FEP
  !! @authors      HO
  !! @param[inout] alchemy : alchemy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_alchemy(alchemy)

    ! formal arguments
    type(s_alchemy), intent(inout) :: alchemy

    ! local variables
    integer                  :: dealloc_stat

    dealloc_stat = 0

    deallocate(alchemy%lambljA, &
        alchemy%lambljB,        &
        alchemy%lambelA,        &
        alchemy%lambelB,        &
        alchemy%lambbondA,        &
        alchemy%lambbondB,        &
        alchemy%lambrest,        &
        stat = dealloc_stat)

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_alchemy


end module sp_alchemy_str_mod
