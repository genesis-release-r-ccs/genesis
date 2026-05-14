!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_pairlist_str_mod
!> @brief   Defines the structure of pairlist information.
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_pairlist_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pairlist
    real(wp)                        :: pairlistdist
    real(wp)                        :: cg_pairlistdist_vdw
    real(wp)                        :: cg_pairlistdist_ele
    real(wp)                        :: cg_pairlistdist_126
    real(wp)                        :: cg_pairlistdist_WF
    real(wp)                        :: cg_pairlistdist_PWMcos
    real(wp)                        :: cg_pairlistdist_DNAbp
    real(wp)                        :: cg_pairlistdist_exv
    integer                         :: num_exv
    integer                         :: num_dna_exv
    integer                         :: num_dna_base
    integer,            allocatable :: num_nb15_calc(:,:)
    integer,            allocatable :: num_nb15_calc1(:,:)
    integer,            allocatable :: num_nb15_nobc(:)
    integer,            allocatable :: num_cg_exv_calc(:)
    integer,            allocatable :: exv_list(:)
    integer,            allocatable :: num_cg_DNA_exv_calc(:)
    integer,            allocatable :: dna_exv_list(:)
    integer,            allocatable :: num_cg_DNA_base_calc(:)
    integer,            allocatable :: dna_base_list(:)
    integer,            allocatable :: num_cg_ele_calc(:)
    integer,            allocatable :: num_cg_pwmcos_calc(:)
    integer,            allocatable :: num_cg_pwmcosns_calc(:)
    integer,            allocatable :: num_cg_kh_calc(:)
    integer,            allocatable :: num_cg_idr_kh_calc(:)
    integer,            allocatable :: num_cg_idr_hps_calc(:)
    integer,            allocatable :: num_cg_idr_mpipi_calc(:)
    integer,            allocatable :: nb15_calc_list(:,:)
    integer,            allocatable :: nb15_calc_list1(:,:)
    integer,            allocatable :: nb15_calc_list_nobc(:,:)
    integer,            allocatable :: cg_exv_list(:,:)
    integer,            allocatable :: cg_DNA_exv_list(:,:)
    integer,            allocatable :: cg_DNA_base_list(:,:)
    integer,            allocatable :: cg_ele_list(:,:)
    integer,            allocatable :: cg_pwmcos_list(:,:)
    integer,            allocatable :: cg_pwmcosns_list(:,:)
    integer,            allocatable :: cg_kh_list(:,:)
    integer,            allocatable :: cg_idr_kh_list(:,:)
    integer,            allocatable :: cg_idr_hps_list(:,:)
    integer,            allocatable :: cg_idr_mpipi_list(:,:)
    integer,            allocatable :: nb15_cell(:)
    integer,            allocatable :: nb15_list(:,:)
    ! for GPU
    integer,            allocatable :: univ_ij_load(:)         !  (ij)
#ifndef PGICUDA
    integer(1),         allocatable :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer(1),         allocatable :: pack_univ_mask2(:)      !  (ix*iy*ij/8)
    integer,            allocatable :: univ_ij_sort_list(:)    !  (ij)
    integer,            allocatable :: univ_ix_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_ix_list(:,:)       !  (iix, ij)
    integer,            allocatable :: univ_iy_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_iy_list(:,:)       !  (iiy, ij)
#else
    integer(1), allocatable, pinned :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer,    allocatable, pinned :: univ_ij_sort_list(:)    !  (ij)
    integer,    allocatable, pinned :: univ_ix_natom(:)        !  (ij)
    integer(1), allocatable, pinned :: univ_ix_list(:,:)       !  (iix, ij)
    integer,    allocatable, pinned :: univ_iy_natom(:)        !  (ij)
    integer(1), allocatable, pinned :: univ_iy_list(:,:)       !  (iiy, ij)
#endif
    integer                 :: univ_ncell_nonzero
    integer                 :: univ_update = 0
    integer                 :: univ_mask2_size = 0
    integer                 :: pack_univ_mask2_size = 0

    integer                 :: realloc_exv
    integer                 :: realloc_DNA_exv
    integer                 :: realloc_DNA_base
    integer                 :: realloc_elec
    integer                 :: realloc_KH
    integer                 :: realloc_IDR_KH
    integer                 :: realloc_IDR_HPS
    integer                 :: realloc_IDR_MPIPI
    integer                 :: realloc_PWMcos
    integer                 :: realloc_PWMcosns

  end type s_pairlist

  ! parameters for allocatable variables
  integer,        public, parameter :: PairListCGExvNum      = 1
  integer,        public, parameter :: PairListCGEleNum      = 2
  integer,        public, parameter :: PairListCGDNAExvNum   = 3
  integer,        public, parameter :: PairListCGDNABaseNum  = 4
  integer,        public, parameter :: PairListCGKHNum       = 5
  integer,        public, parameter :: PairListCGIDRKHNum    = 6
  integer,        public, parameter :: PairListCGIDRHPSNum   = 7
  integer,        public, parameter :: PairListPWMCOSNum     = 8
  integer,        public, parameter :: PairListPWMCOSnsNum   = 9 
  integer,        public, parameter :: PairListCGIDRMPIPINum = 10
  ! 
  integer,        public, parameter :: PairListCGExvList      = 111
  integer,        public, parameter :: PairListCGEleList      = 112
  integer,        public, parameter :: PairListCGDNAExvList   = 113
  integer,        public, parameter :: PairListCGDNABaseList  = 114
  integer,        public, parameter :: PairListCGKHList       = 115
  integer,        public, parameter :: PairListCGIDRKHList    = 116
  integer,        public, parameter :: PairListCGIDRHPSList   = 117
  integer,        public, parameter :: PairListPWMCOSList     = 118
  integer,        public, parameter :: PairListPWMCOSnsList   = 119
  integer,        public, parameter :: PairListCGIDRMPIPIList = 120
  
  ! variables for maximum numbers in one cell
  integer,        public            :: Max_exv_nb15
  integer,        public            :: Max_dna_exv_nb15
  integer,        public            :: Max_dna_base_nb15
  integer,        public            :: Max_elec_nb15
  integer,        public            :: Max_KH_nb15
  integer,        public            :: Max_IDR_KH_nb15
  integer,        public            :: Max_IDR_HPS_nb15
  integer,        public            :: Max_IDR_MPIPI_nb15
  integer,        public            :: Max_pwmcos_nb15
  integer,        public            :: Max_pwmcosns_nb15

  ! subroutines
  public  :: init_pairlist
  public  :: alloc_pairlist
  public  :: dealloc_pairlist
  public  :: dealloc_pairlist_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_pairlist
  !> @brief        Initializes pairlist information.
  !! @authors      YS
  !! @param[out]   pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_pairlist(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    pairlist%pairlistdist = 0.0_wp

    return

  end subroutine init_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist
  !> @brief        Allocates pairlist information.
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist(pairlist, variable, var_size, var_size1)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,    optional,    intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (PairListCGExvNum)

      if (allocated(pairlist%num_cg_exv_calc)) then
        if (size(pairlist%num_cg_exv_calc) /= var_size)              &
          deallocate(pairlist%num_cg_exv_calc,                       &
                     pairlist%exv_list,                              &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_exv_calc))                 &
        allocate(pairlist%num_cg_exv_calc(var_size),                 &
                 pairlist%exv_list(var_size),                        &
                 stat = alloc_stat)
      pairlist%num_cg_exv_calc(1:var_size)  = 0

    case (PairListCGEleNum)

      if (allocated(pairlist%num_cg_ele_calc)) then
        if (size(pairlist%num_cg_ele_calc) /= var_size)              &
          deallocate(pairlist%num_cg_ele_calc,                       & 
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_ele_calc))                 &
        allocate(pairlist%num_cg_ele_calc(var_size),                 &
                 stat = alloc_stat)
      pairlist%num_cg_ele_calc(1:var_size) = 0

    case (PairListCGDNAExvNum)

      if (allocated(pairlist%num_cg_DNA_exv_calc)) then
        if (size(pairlist%num_cg_DNA_exv_calc) /= var_size)          &
          deallocate(pairlist%num_cg_DNA_exv_calc,                   &
                     pairlist%dna_exv_list,                          &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_DNA_exv_calc))             &
        allocate(pairlist%num_cg_DNA_exv_calc(var_size),             &
                 pairlist%dna_exv_list(var_size),                    &
                 stat = alloc_stat)
      pairlist%num_cg_DNA_exv_calc(1:var_size)  = 0

    case (PairListCGDNABaseNum)

      if (allocated(pairlist%num_cg_DNA_base_calc)) then
        if (size(pairlist%num_cg_DNA_base_calc) /= var_size)         &
          deallocate(pairlist%num_cg_DNA_base_calc,                  &
                     pairlist%dna_base_list,                         &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_DNA_base_calc))            &
        allocate(pairlist%num_cg_DNA_base_calc(var_size),            &
                 pairlist%dna_base_list(var_size),                   &
                 stat = alloc_stat)
      pairlist%num_cg_DNA_base_calc(1:var_size)  = 0


    case (PairListPWMCOSNum)

      if (allocated(pairlist%num_cg_pwmcos_calc)) then
        if (size(pairlist%num_cg_pwmcos_calc(:)) /= var_size)        &
          deallocate(pairlist%num_cg_pwmcos_calc,                    &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_pwmcos_calc))              &
        allocate(pairlist%num_cg_pwmcos_calc(var_size),              &
                 stat = alloc_stat)
      pairlist%num_cg_pwmcos_calc(1:var_size) = 0

    case (PairListPWMCOSnsNum)

      if (allocated(pairlist%num_cg_pwmcosns_calc)) then
        if (size(pairlist%num_cg_pwmcosns_calc(:)) /= var_size)      &
          deallocate(pairlist%num_cg_pwmcosns_calc,                  &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_pwmcosns_calc))            &
        allocate(pairlist%num_cg_pwmcosns_calc(var_size),            &
                 stat = alloc_stat)
      pairlist%num_cg_pwmcosns_calc(1:var_size) = 0

    case (PairListCGKHNum)

      if (allocated(pairlist%num_cg_kh_calc)) then
        if (size(pairlist%num_cg_kh_calc) /= var_size)               &
          deallocate(pairlist%num_cg_kh_calc,                        &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_kh_calc))                  &
        allocate(pairlist%num_cg_kh_calc(var_size),                  &
                 stat = alloc_stat)
      pairlist%num_cg_kh_calc  (1:var_size) = 0

    case (PairListCGIDRKHNum)

      if (allocated(pairlist%num_cg_idr_kh_calc)) then
        if (size(pairlist%num_cg_idr_kh_calc) /= var_size)           &
          deallocate(pairlist%num_cg_idr_kh_calc,                    &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_idr_kh_calc))              &
        allocate(pairlist%num_cg_idr_kh_calc(var_size),              &
                 stat = alloc_stat)
      pairlist%num_cg_idr_kh_calc  (1:var_size) = 0

    case (PairListCGIDRHPSNum)

      if (allocated(pairlist%num_cg_idr_hps_calc)) then
        if (size(pairlist%num_cg_idr_hps_calc) /= var_size)          &
          deallocate(pairlist%num_cg_idr_hps_calc,                   &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_idr_hps_calc))             &
        allocate(pairlist%num_cg_idr_hps_calc(var_size),             &
                 stat = alloc_stat)
      pairlist%num_cg_idr_hps_calc  (1:var_size) = 0

    case (PairListCGIDRMPIPINum)

      if (allocated(pairlist%num_cg_idr_mpipi_calc)) then
        if (size(pairlist%num_cg_idr_mpipi_calc) /= var_size)          &
            deallocate(pairlist%num_cg_idr_mpipi_calc,                   &
            stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%num_cg_idr_mpipi_calc))             &
          allocate(pairlist%num_cg_idr_mpipi_calc(var_size),             &
          stat = alloc_stat)
      pairlist%num_cg_idr_mpipi_calc  (1:var_size) = 0

    case (PairListCGExvList)

      if (allocated(pairlist%cg_exv_list)) then
        if (size(pairlist%cg_exv_list) /= var_size*var_size1)        &
          deallocate(pairlist%cg_exv_list,                           &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_exv_list))                     &
        allocate(pairlist%cg_exv_list(var_size, var_size1),          &
                 stat = alloc_stat)

    case (PairListCGEleList)

      if (allocated(pairlist%cg_ele_list)) then
        if (size(pairlist%cg_ele_list) /= var_size*var_size1)        &
          deallocate(pairlist%cg_ele_list,                           &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_ele_list))                     &
        allocate(pairlist%cg_ele_list(var_size, var_size1),          &
                 stat = alloc_stat)

    case (PairListCGDNAExvList)

      if (allocated(pairlist%cg_DNA_exv_list)) then
        if (size(pairlist%cg_DNA_exv_list) /= var_size*var_size1)    &
          deallocate(pairlist%cg_DNA_exv_list,                       &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_DNA_exv_list))                 &
        allocate(pairlist%cg_DNA_exv_list(var_size, var_size1),      &
                 stat = alloc_stat)

    case (PairListCGDNABaseList)

      if (allocated(pairlist%cg_DNA_base_list)) then
        if (size(pairlist%cg_DNA_base_list) /= var_size*var_size1)   &
          deallocate(pairlist%cg_DNA_base_list,                      &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_DNA_base_list))                &
        allocate(pairlist%cg_DNA_base_list(var_size, var_size1),     &
                 stat = alloc_stat)


    case (PairListPwmCosList)

      if (allocated(pairlist%cg_pwmcos_list)) then
        if (size(pairlist%cg_pwmcos_list) /= var_size*var_size1)     &
          deallocate(pairlist%cg_pwmcos_list,                        &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_pwmcos_list))                  &
        allocate(pairlist%cg_pwmcos_list(var_size, var_size1),       &
                 stat = alloc_stat)

    case (PairListPWMCOSnsList)

      if (allocated(pairlist%cg_pwmcosns_list)) then
        if (size(pairlist%cg_pwmcosns_list) /= var_size*var_size1)   &
          deallocate(pairlist%cg_pwmcosns_list,                      &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_pwmcosns_list))                &
        allocate(pairlist%cg_pwmcosns_list(var_size,var_size1),      &
                 stat = alloc_stat)

    case (PairListCGKHList)

      if (allocated(pairlist%cg_kh_list)) then
        if (size(pairlist%cg_kh_list) /= var_size*var_size1)         &
          deallocate(pairlist%cg_kh_list,                            &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_kh_list))                      &
        allocate(pairlist%cg_kh_list(var_size, var_size1),           &
                 stat = alloc_stat)

    case (PairListCGIDRKHList)

      if (allocated(pairlist%cg_idr_kh_list)) then
        if (size(pairlist%cg_idr_kh_list) /= var_size*var_size1)     &
          deallocate(pairlist%cg_idr_kh_list,                        &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_idr_kh_list))                  &
        allocate(pairlist%cg_idr_kh_list(var_size, var_size1),       &
                 stat = alloc_stat)

    case (PairListCGIDRHPSList)

      if (allocated(pairlist%cg_idr_hps_list)) then
        if (size(pairlist%cg_idr_hps_list) /= var_size*var_size1)    &
          deallocate(pairlist%cg_idr_hps_list,                       &
                     stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_idr_hps_list))                 &
        allocate(pairlist%cg_idr_hps_list(var_size, var_size1),      &
                 stat = alloc_stat)

    case (PairListCGIDRMPIPIList)

      if (allocated(pairlist%cg_idr_mpipi_list)) then
        if (size(pairlist%cg_idr_mpipi_list) /= var_size*var_size1)    &
            deallocate(pairlist%cg_idr_mpipi_list,                       &
            stat = dealloc_stat)
      end if
      if (.not. allocated(pairlist%cg_idr_mpipi_list))                 &
          allocate(pairlist%cg_idr_mpipi_list(var_size, var_size1),      &
          stat = alloc_stat)

    end select

    if (alloc_stat /=0)   &
      call error_msg('Alloc_Pairlist: Failed to allocate memory.')
    if (dealloc_stat /=0) &
      call error_msg('Alloc_Pairlist: Failed to deallocate memory.')

    return

  end subroutine alloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist
  !> @brief        Deallocates pairlist information.
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist(pairlist, variable)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (PairListCGExvNum)

      if (allocated(pairlist%num_cg_exv_calc)) then
          deallocate(pairlist%num_cg_exv_calc,                     &
                     stat = dealloc_stat)
      end if

    case (PairListCGEleNum)

      if (allocated(pairlist%num_cg_ele_calc)) then
          deallocate(pairlist%num_cg_ele_calc,                     & 
                     stat = dealloc_stat)
      end if

    case (PairListCGDNAExvNum)

      if (allocated(pairlist%num_cg_DNA_exv_calc)) then
          deallocate(pairlist%num_cg_DNA_exv_calc,                 &
                     stat = dealloc_stat)
      end if

    case (PairListCGDNABaseNum)

      if (allocated(pairlist%num_cg_DNA_base_calc)) then
          deallocate(pairlist%num_cg_DNA_base_calc,                &
                     stat = dealloc_stat)
      end if

    case (PairListPWMCOSNum)

      if (allocated(pairlist%num_cg_pwmcos_calc)) then
          deallocate(pairlist%num_cg_pwmcos_calc,                  &
                     stat = dealloc_stat)
      end if

    case (PairListPWMCOSnsNum)

      if (allocated(pairlist%num_cg_pwmcosns_calc)) then
          deallocate(pairlist%num_cg_pwmcosns_calc,                &
                     stat = dealloc_stat)
      end if

    case (PairListCGKHNum)

      if (allocated(pairlist%num_cg_kh_calc)) then
          deallocate(pairlist%num_cg_kh_calc,                      &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRKHNum)

      if (allocated(pairlist%num_cg_idr_kh_calc)) then
          deallocate(pairlist%num_cg_idr_kh_calc,                  &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRHPSNum)

      if (allocated(pairlist%num_cg_idr_hps_calc)) then
          deallocate(pairlist%num_cg_idr_hps_calc,                 &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRMPIPINum)

      if (allocated(pairlist%num_cg_idr_mpipi_calc)) then
        deallocate(pairlist%num_cg_idr_mpipi_calc,                 &
            stat = dealloc_stat)
      end if

    case (PairListCGExvList)

      if (allocated(pairlist%cg_exv_list)) then
          deallocate(pairlist%cg_exv_list,                         &
                     stat = dealloc_stat)
      end if

    case (PairListCGEleList)

      if (allocated(pairlist%cg_ele_list)) then
          deallocate(pairlist%cg_ele_list,                         &
                     stat = dealloc_stat)
      end if

    case (PairListCGDNAExvList)

      if (allocated(pairlist%cg_DNA_exv_list)) then
          deallocate(pairlist%cg_DNA_exv_list,                     &
                     stat = dealloc_stat)
      end if

    case (PairListCGDNABaseList)

      if (allocated(pairlist%cg_DNA_base_list)) then
          deallocate(pairlist%cg_DNA_base_list,                     &
                     stat = dealloc_stat)
      end if

    case (PairListPWMCOSList)

      if (allocated(pairlist%cg_pwmcos_list)) then
          deallocate(pairlist%cg_pwmcos_list,                       &
                     stat = dealloc_stat)
      end if

    case (PairListPWMCOSnsList)

      if (allocated(pairlist%cg_pwmcosns_list)) then
          deallocate(pairlist%num_cg_pwmcosns_calc,                 &
                     stat = dealloc_stat)
      end if

    case (PairListCGKHList)

      if (allocated(pairlist%cg_kh_list)) then
          deallocate(pairlist%cg_kh_list,                           &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRKHList)

      if (allocated(pairlist%cg_idr_kh_list)) then
          deallocate(pairlist%cg_idr_kh_list,                       &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRHPSList)

      if (allocated(pairlist%cg_idr_hps_list)) then
          deallocate(pairlist%cg_idr_hps_list,                      &
                     stat = dealloc_stat)
      end if

    case (PairListCGIDRMPIPIList)

      if (allocated(pairlist%cg_idr_mpipi_list)) then
        deallocate(pairlist%cg_idr_mpipi_list,                      &
            stat = dealloc_stat)
      end if

    end select

    if (dealloc_stat /=0) &
      call error_msg('Dealloc_Pairlist: Failed to deallocate memory.')

    return

  end subroutine dealloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist_all
  !> @brief        Deallocates all pairlist information.
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist_all(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


!   call dealloc_pairlist(pairlist, PairListNoTable)
!   call dealloc_pairlist(pairlist, PairListTable)

    return

  end subroutine dealloc_pairlist_all

end module cg_pairlist_str_mod
