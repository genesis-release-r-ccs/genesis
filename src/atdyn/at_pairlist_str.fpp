!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_pairlist_str_mod
!> @brief   structure of pairlist information
!! @authors Yuji Sugita (YS), Takashi Imai (TI), Jaewoon Jung (JJ), 
!!          Takaharu Mori (TM), Motoshi Kamiya (MK), Kiyoshi Yagi (KY)
!!          Cheng Tan (CT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_pairlist_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_table_pair
    integer                       :: num_nb15_max
    integer                       :: water_nb15_max
    integer,          allocatable :: num_nb15_calc(:)
    integer,          allocatable :: num_nb15_calcw(:)
    integer,          allocatable :: nb15_calc_list(:,:)
    integer,          allocatable :: nb15_calc_listw(:,:)
    integer,          allocatable :: cell_linked_list(:)
    integer,          allocatable :: cell_linked_listw(:)
    integer,          allocatable :: atom_cell_index(:)
    integer,          allocatable :: num_nb15_calc_water(:)
    integer,          allocatable :: nb15_calc_list_water(:,:)
    integer,          allocatable :: cell_linked_list_water(:)
    integer,          allocatable :: atom_cell_index_water(:)
    integer,          allocatable :: num_list_water(:)
  end type s_table_pair

  type, public :: s_pairlist
    type(s_table_pair)            :: table
    logical                       :: allocation
    logical                       :: allocate_nobc
    logical                       :: allocate_pbc
    logical                       :: allocate_solsol
    logical                       :: allocate_solwat
    logical                       :: allocate_watwat
    integer                       :: num_nb15_max
    integer,          allocatable :: num_nb15_pre(:)
    integer,          allocatable :: num_nb15(:)
    integer,          allocatable :: num_nb15_calc(:,:)
    integer,          allocatable :: nb15_calc_list(:,:)
    real(wp)                      :: pairlistdist
    integer                       :: ecqm_num_nb15
    integer, allocatable          :: ecqm_nb15_list(:,:)
    ! for GBSA
    integer                       :: num_all_max
    integer,          allocatable :: num_all_pre(:)
    integer,          allocatable :: num_all(:)
    integer,          allocatable :: num_all_calc(:,:)
    integer,          allocatable :: all_calc_list(:,:)
    ! for CG models: AICG2P + 3SPN.2C
    logical                       :: allocate_nobc_cg
    logical                       :: allocate_nobc_cg_pwmcos
    logical                       :: allocate_nobc_cg_pwmcosns
    logical                       :: allocate_nobc_cg_IDR_HPS
    logical                       :: allocate_nobc_cg_IDR_KH
    ! 
    integer,          allocatable :: cg_DNA_basepair_list(:,:)
    integer,          allocatable :: cg_DNA_exv_list(:,:)
    integer,          allocatable :: cg_ele_list(:,:)
    real(wp),         allocatable :: cg_ele_scaling_list(:,:)
    integer,          allocatable :: cg_exv_list(:,:)
    integer,          allocatable :: cg_pwmcos_list(:,:)
    integer,          allocatable :: cg_pwmcosns_list(:,:)
    integer,          allocatable :: cg_IDR_HPS_list(:,:)
    integer,          allocatable :: cg_IDR_KH_list(:,:)
    integer,          allocatable :: cg_KH_list(:,:)
    integer,          allocatable :: cg_KH_model_list(:,:)
    ! 
    integer,          allocatable :: num_cg_DNA_basepair_calc(:,:)
    integer,          allocatable :: num_cg_DNA_exv_calc(:,:)
    integer,          allocatable :: num_cg_ele_calc(:,:)
    integer,          allocatable :: num_cg_exv_calc(:,:)
    integer,          allocatable :: num_cg_pwmcos_calc(:,:)
    integer,          allocatable :: num_cg_pwmcosns_calc(:,:)
    integer,          allocatable :: num_cg_IDR_HPS_calc(:,:)
    integer,          allocatable :: num_cg_IDR_KH_calc(:,:)
    integer,          allocatable :: num_cg_KH_calc(:,:)
    ! 
    integer,          allocatable :: num_cg_DNA_basepair_pre(:)
    integer,          allocatable :: num_cg_DNA_exv_pre(:)
    integer,          allocatable :: num_cg_ele_pre(:)
    integer,          allocatable :: num_cg_exv_pre(:)
    integer,          allocatable :: num_cg_pwmcos_pre(:)
    integer,          allocatable :: num_cg_pwmcosns_pre(:)
    integer,          allocatable :: num_cg_IDR_HPS_pre(:)
    integer,          allocatable :: num_cg_IDR_KH_pre(:)
    integer,          allocatable :: num_cg_KH_pre(:)
    ! 
    integer,          allocatable :: num_cg_DNA_basepair(:)
    integer,          allocatable :: num_cg_DNA_exv(:)
    integer,          allocatable :: num_cg_ele(:)
    integer,          allocatable :: num_cg_exv(:)
    integer,          allocatable :: num_cg_pwmcos(:)
    integer,          allocatable :: num_cg_pwmcosns(:)
    integer,          allocatable :: num_cg_IDR_HPS(:)
    integer,          allocatable :: num_cg_IDR_KH(:)
    integer,          allocatable :: num_cg_KH(:)
    ! 
    integer                       :: num_cg_DNA_basepair_max
    integer                       :: num_cg_DNA_exv_max
    integer                       :: num_cg_ele_max
    integer                       :: num_cg_exv_max
    integer                       :: num_cg_pwmcos_max
    integer                       :: num_cg_pwmcosns_max
    integer                       :: num_cg_IDR_HPS_max
    integer                       :: num_cg_IDR_KH_max
    integer                       :: num_cg_KH_max
    ! 
    integer,          allocatable :: cell_index_cg_all(:)
    integer,          allocatable :: cell_linked_list_cg_all(:)
    integer,          allocatable :: cell_linked_list_cg_DNA(:)
    integer,          allocatable :: cell_linked_list_cg_DNA_phos(:)
    integer,          allocatable :: cell_linked_list_cg_DNA_base(:)
    integer,          allocatable :: cell_linked_list_cg_IDR_HPS(:)
    integer,          allocatable :: cell_linked_list_cg_IDR_KH(:)
    integer,          allocatable :: cell_linked_list_cg_KH(:)
    integer,          allocatable :: cell_linked_list_cg_charged(:)
    integer,          allocatable :: cell_head_index_cg_all(:)
    integer,          allocatable :: cell_head_index_cg_DNA(:)
    integer,          allocatable :: cell_head_index_cg_DNA_phos(:)
    integer,          allocatable :: cell_head_index_cg_DNA_base(:)
    integer,          allocatable :: cell_head_index_cg_IDR_HPS(:)
    integer,          allocatable :: cell_head_index_cg_IDR_KH(:)
    integer,          allocatable :: cell_head_index_cg_KH(:)
    integer,          allocatable :: cell_head_index_cg_charged(:)
    !
    logical                       :: allocate_pbc_cg_exv
    logical                       :: allocate_pbc_cg_ele
    logical                       :: allocate_pbc_cg_DNA_bp
    logical                       :: allocate_pbc_cg_DNA_exv
    logical                       :: allocate_pbc_cg_pwmcos
    logical                       :: allocate_pbc_cg_pwmcosns
    logical                       :: allocate_pbc_cg_IDR_HPS
    logical                       :: allocate_pbc_cg_IDR_KH
    logical                       :: allocate_pbc_cg_KH

    ! 
    real(wp)                      :: cg_pairlistdist_ele
    real(wp)                      :: cg_pairlistdist_126
    real(wp)                      :: cg_pairlistdist_PWMcos
    real(wp)                      :: cg_pairlistdist_DNAbp
    real(wp)                      :: cg_pairlistdist_exv
  end type s_pairlist

  ! parameters for allocatable variables
  integer,      public, parameter :: PairListAtomNobc        = 1
  integer,      public, parameter :: PairListIntNobc         = 2
  integer,      public, parameter :: PairListPbcSolute       = 3
  integer,      public, parameter :: PairListPbcWater        = 4
  integer,      public, parameter :: PairListPbcSoluteSolute = 5
  integer,      public, parameter :: PairListPbcSoluteWater  = 6
  integer,      public, parameter :: PairListPbcWaterWater   = 7
  integer,      public, parameter :: PairListNthreads        = 8
  integer,      public, parameter :: PairListEcqm            = 9
  integer,      public, parameter :: PairListAtomNobcGbsa    = 10
  integer,      public, parameter :: PairListIntNobcGbsa     = 11
  integer,      public, parameter :: PairListNthreadsGbsa    = 12
  integer,      public, parameter :: PairListAtomNobcCG      = 13
  integer,      public, parameter :: PairListNthreadsCG      = 14
  integer,      public, parameter :: PairListCGDNABP         = 15
  integer,      public, parameter :: PairListCGDNAexv        = 16
  integer,      public, parameter :: PairListCGele           = 17
  integer,      public, parameter :: PairListCGexv           = 18
  integer,      public, parameter :: PairListCGPWMcos        = 19
  integer,      public, parameter :: PairListCGPWMcosns      = 20
  integer,      public, parameter :: PairListCGIDRHPS        = 21
  integer,      public, parameter :: PairListCGIDRKH         = 22
  integer,      public, parameter :: PairListCGKH            = 23
  integer,      public, parameter :: PairListNthCGPWMcos     = 24
  integer,      public, parameter :: PairListNthCGPWMcosns   = 25
  ! 
  integer,      public, parameter :: PairListNthreadsPbcCG   = 30
  integer,      public, parameter :: PairListAtomPbcCGexv    = 31
  integer,      public, parameter :: PairListAtomPbcCGele    = 32
  integer,      public, parameter :: PairListAtomPbcCGDNAexv = 33
  integer,      public, parameter :: PairListAtomPbcCGDNAbp  = 34
  integer,      public, parameter :: PairListAtomPbcCGPWMcos = 35
  integer,      public, parameter :: PairListAtomPbcCGPWMcosns = 36
  integer,      public, parameter :: PairListAtomPbcCGIDRHPS = 37
  integer,      public, parameter :: PairListAtomPbcCGIDRKH  = 38
  integer,      public, parameter :: PairListAtomPbcCGKH     = 39
  ! 
  integer,      public, parameter :: PairListCellsPbcCG      = 50
  integer,      public, parameter :: PairListPbcCGDNAbp      = 51
  integer,      public, parameter :: PairListPbcCGDNAexv     = 52
  integer,      public, parameter :: PairListPbcCGele        = 53
  integer,      public, parameter :: PairListPbcCGexv        = 54
  integer,      public, parameter :: PairListPbcCGPWMcos     = 55
  integer,      public, parameter :: PairListPbcCGPWMcosns   = 56
  integer,      public, parameter :: PairListPbcCGIDRHPS     = 57
  integer,      public, parameter :: PairListPbcCGIDRKH      = 58
  integer,      public, parameter :: PairListPbcCGKH         = 59

  ! subroutines
  public :: init_pairlist
  public :: alloc_pairlist
  public :: alloc_pairlist2
  public :: dealloc_pairlist
  public :: dealloc_pairlist_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_pairlist
  !> @brief        initialize pairlist information
  !! @authors      YS
  !! @param[out]   pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_pairlist(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    pairlist%table%num_nb15_max          = 0
    pairlist%table%water_nb15_max        = 0

    pairlist%allocation                  = .true.
    pairlist%allocate_nobc               = .true.
    pairlist%allocate_pbc                = .true.
    pairlist%allocate_solsol             = .true.
    pairlist%allocate_solwat             = .true.
    pairlist%allocate_watwat             = .true.
    pairlist%allocate_nobc_cg            = .true.
    pairlist%allocate_nobc_cg_pwmcos     = .true.
    pairlist%allocate_nobc_cg_pwmcosns   = .true.
    pairlist%allocate_nobc_cg_IDR_HPS    = .true.
    pairlist%allocate_nobc_cg_IDR_KH     = .true.
    pairlist%allocate_pbc_cg_exv         = .true.
    pairlist%allocate_pbc_cg_ele         = .true.
    pairlist%allocate_pbc_cg_DNA_bp      = .true.
    pairlist%allocate_pbc_cg_DNA_exv     = .true.
    pairlist%allocate_pbc_cg_pwmcos      = .true.
    pairlist%allocate_pbc_cg_pwmcosns    = .true.
    pairlist%allocate_pbc_cg_IDR_HPS     = .true.
    pairlist%allocate_pbc_cg_IDR_KH      = .true.
    pairlist%allocate_pbc_cg_KH          = .true.
    pairlist%num_cg_DNA_basepair_max     = 0
    pairlist%num_cg_DNA_exv_max          = 0
    pairlist%num_cg_ele_max              = 0
    pairlist%num_cg_exv_max              = 0
    pairlist%num_cg_pwmcos_max           = 0
    pairlist%num_cg_pwmcosns_max         = 0
    pairlist%num_cg_IDR_HPS_max          = 0
    pairlist%num_cg_IDR_KH_max           = 0
    pairlist%num_cg_KH_max               = 0
    pairlist%num_nb15_max                = 0
    pairlist%pairlistdist                = 0.0_wp
    pairlist%cg_pairlistdist_ele         = 0.0_wp
    pairlist%cg_pairlistdist_126         = 0.0_wp
    pairlist%cg_pairlistdist_PWMcos      = 0.0_wp
    pairlist%cg_pairlistdist_DNAbp       = 0.0_wp
    pairlist%cg_pairlistdist_exv         = 0.0_wp

    return

  end subroutine init_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist
  !> @brief        allocate pairlist information
  !! @authors      YS, TI, JJ, TM, CT
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist(pairlist, variable, var_size, var_size2)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,       optional, intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: nthread, omp_get_max_threads


    alloc_stat   = 0
    dealloc_stat = 0

    ! get number of threads
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! allocate selected variables
    !
    select case (variable)

    case (PairListAtomNobc)

      if (allocated(pairlist%num_nb15_calc)) then
        if (size(pairlist%num_nb15_calc) == var_size*nthread) return
        deallocate(pairlist%num_nb15_calc, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%num_nb15_calc(var_size,nthread), &
               stat = alloc_stat)

    case (PairListIntNobc)

      if (allocated(pairlist%nb15_calc_list)) then
        if (size(pairlist%nb15_calc_list) == var_size*nthread) return
        deallocate(pairlist%nb15_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%nb15_calc_list(var_size,nthread), &
               stat = alloc_stat)
      
    case (PairListAtomNobcCG)

      if (allocated(pairlist%num_cg_DNA_basepair_calc)) then
        if (size(pairlist%num_cg_DNA_basepair_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_DNA_basepair_calc, &
            pairlist%num_cg_DNA_exv_calc,             &
            pairlist%num_cg_ele_calc,                 &
            pairlist%num_cg_exv_calc,                 &
            pairlist%num_cg_IDR_HPS_calc,             &
            pairlist%num_cg_IDR_KH_calc,              &
            pairlist%num_cg_KH_calc,                  &
            stat = dealloc_stat)
      end if

      allocate(pairlist%num_cg_DNA_basepair_calc(var_size,nthread), &
          pairlist%num_cg_DNA_exv_calc(var_size,nthread),           &
          pairlist%num_cg_ele_calc(var_size,nthread),               &
          pairlist%num_cg_exv_calc(var_size,nthread),               &
          pairlist%num_cg_IDR_HPS_calc(var_size,nthread),           &
          pairlist%num_cg_IDR_KH_calc(var_size,nthread),            &
          pairlist%num_cg_KH_calc(var_size,nthread),                &
          stat = alloc_stat)

    case (PairListCGDNABP)

      if (allocated(pairlist%cg_DNA_basepair_list)) then
        if (size(pairlist%cg_DNA_basepair_list) == var_size*nthread) return
        deallocate(pairlist%cg_DNA_basepair_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_DNA_basepair_list(var_size,nthread), &
          stat = alloc_stat)

    case (PairListCGDNAexv)
      if (allocated(pairlist%cg_DNA_exv_list)) then
        if (size(pairlist%cg_DNA_exv_list) == var_size*nthread) return
        deallocate(pairlist%cg_DNA_exv_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_DNA_exv_list(var_size,nthread), &
          stat = alloc_stat)

    case (PairListCGele)
      if (allocated(pairlist%cg_ele_list)) then
        if (size(pairlist%cg_ele_list) == var_size*nthread) return
        deallocate (pairlist%cg_ele_list, &
            pairlist%cg_ele_scaling_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_ele_list(var_size,nthread), &
          stat = alloc_stat)
      allocate(pairlist%cg_ele_scaling_list(var_size,nthread), &
          stat = alloc_stat)
      pairlist%cg_ele_scaling_list(:, :) = 1.0_wp

    case (PairListCGKH)
      if (allocated(pairlist%cg_KH_list)) then
        if (size(pairlist%cg_KH_list) == var_size*nthread) return
        deallocate (pairlist%cg_KH_list, &
            pairlist%cg_KH_model_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_KH_list(var_size,nthread), &
          stat = alloc_stat)
      allocate(pairlist%cg_KH_model_list(var_size,nthread), &
          stat = alloc_stat)
      pairlist%cg_KH_model_list(:, :) = 0

    case (PairListCGIDRHPS)
      if (allocated(pairlist%cg_IDR_HPS_list)) then
        if (size(pairlist%cg_IDR_HPS_list) == var_size*nthread) return
        deallocate (pairlist%cg_IDR_HPS_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_IDR_HPS_list(var_size,nthread), &
          stat = alloc_stat)

    case (PairListCGIDRKH)
      if (allocated(pairlist%cg_IDR_KH_list)) then
        if (size(pairlist%cg_IDR_KH_list) == var_size*nthread) return
        deallocate (pairlist%cg_IDR_KH_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_IDR_KH_list(var_size,nthread), &
          stat = alloc_stat)

    case (PairListCGexv)
      if (allocated(pairlist%cg_exv_list)) then
        if (size(pairlist%cg_exv_list) == var_size*nthread) return
        deallocate(pairlist%cg_exv_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_exv_list(var_size,nthread), &
          stat = alloc_stat)

    case (PairListCGPWMcos)
      if (allocated(pairlist%cg_pwmcos_list)) then
        if (size(pairlist%cg_pwmcos_list) == var_size*var_size2) return
        deallocate(pairlist%cg_pwmcos_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_pwmcos_list(var_size,var_size2), &
          stat = alloc_stat)

    case (PairListCGPWMcosns)
      if (allocated(pairlist%cg_pwmcosns_list)) then
        if (size(pairlist%cg_pwmcosns_list) == var_size*var_size2) return
        deallocate(pairlist%cg_pwmcosns_list, &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_pwmcosns_list(var_size,var_size2), &
          stat = alloc_stat)

    case (PairListPbcSolute)

      if (allocated(pairlist%table%num_nb15_calc)) then
        if (size(pairlist%table%num_nb15_calc) == var_size) return
        deallocate(pairlist%table%num_nb15_calc,    &
                   pairlist%table%num_nb15_calcw,   &
                   pairlist%table%atom_cell_index,  &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%num_nb15_calc(var_size),    &
               pairlist%table%num_nb15_calcw(var_size),   &
               pairlist%table%cell_linked_list(var_size), &
               pairlist%table%atom_cell_index(var_size),  &
               stat = alloc_stat)

    case (PairListPbcWater)

      if (allocated(pairlist%table%num_nb15_calc_water)) then
        if (size(pairlist%table%num_nb15_calc_water) == var_size) return
        deallocate(pairlist%table%num_nb15_calc_water,    &
                   pairlist%table%cell_linked_list_water, &
                   pairlist%table%atom_cell_index_water,  &
                   pairlist%table%cell_linked_listw,      &
                   pairlist%table%num_list_water,         &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%num_nb15_calc_water(var_size),    &
               pairlist%table%cell_linked_list_water(var_size), &
               pairlist%table%atom_cell_index_water(var_size),  &
               pairlist%table%cell_linked_listw(3*var_size),    &
               pairlist%table%num_list_water(var_size),         &
               stat = alloc_stat)

    case (PairListNthreads)

      if (allocated(pairlist%num_nb15_pre)) then
        if (size(pairlist%num_nb15_pre) == nthread) return
        deallocate(pairlist%num_nb15_pre, &
                   pairlist%num_nb15,     &
                   stat = dealloc_stat)
      end if
      allocate(pairlist%num_nb15_pre(nthread), &
               pairlist%num_nb15    (nthread), &
               stat = alloc_stat)
      
    case (PairListNthreadsCG)

      if (allocated(pairlist%num_cg_DNA_basepair_pre)) then
        if (size(pairlist%num_cg_DNA_basepair_pre) == nthread) return
        deallocate(pairlist%num_cg_DNA_basepair_pre,       &
            pairlist%num_cg_DNA_basepair,                  &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_DNA_basepair_pre (nthread), &
          pairlist%num_cg_DNA_basepair          (nthread), &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_DNA_exv_pre)) then
        if (size(pairlist%num_cg_DNA_exv_pre) == nthread) return
        deallocate(pairlist%num_cg_DNA_exv_pre,            &
            pairlist%num_cg_DNA_exv,                       &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_DNA_exv_pre (nthread),      &
          pairlist%num_cg_DNA_exv          (nthread),      &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_ele_pre)) then
        if (size(pairlist%num_cg_ele_pre) == nthread) return
        deallocate(pairlist%num_cg_ele_pre,                &
            pairlist%num_cg_ele,                           &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_ele_pre (nthread),          &
          pairlist%num_cg_ele          (nthread),          &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_KH_pre)) then
        if (size(pairlist%num_cg_KH_pre) == nthread) return
        deallocate(pairlist%num_cg_KH_pre,                 &
            pairlist%num_cg_KH,                            &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_KH_pre (nthread),           &
          pairlist%num_cg_KH          (nthread),           &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_IDR_HPS_pre)) then
        if (size(pairlist%num_cg_IDR_HPS_pre) == nthread) return
        deallocate(pairlist%num_cg_IDR_HPS_pre,            &
            pairlist%num_cg_IDR_HPS,                       &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_IDR_HPS_pre (nthread),      &
          pairlist%num_cg_IDR_HPS          (nthread),      &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_IDR_KH_pre)) then
        if (size(pairlist%num_cg_IDR_KH_pre) == nthread) return
        deallocate(pairlist%num_cg_IDR_KH_pre,            &
            pairlist%num_cg_IDR_KH,                       &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_IDR_KH_pre (nthread),      &
          pairlist%num_cg_IDR_KH          (nthread),      &
          stat = alloc_stat)

      if (allocated(pairlist%num_cg_exv_pre)) then
        if (size(pairlist%num_cg_exv_pre) == nthread) return
        deallocate(pairlist%num_cg_exv_pre,                &
            pairlist%num_cg_exv,                           &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_exv_pre (nthread),          &
          pairlist%num_cg_exv          (nthread),          &
          stat = alloc_stat)

    case (PairListNthCGPWMcos)
      
      if (allocated(pairlist%num_cg_pwmcos_pre)) then
        if (size(pairlist%num_cg_pwmcos_pre) == var_size) return
        deallocate(pairlist%num_cg_pwmcos_pre,       &
            pairlist%num_cg_pwmcos,                  &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_pwmcos_pre(var_size), &
          pairlist%num_cg_pwmcos(var_size),          &
          stat = alloc_stat)
      
    case (PairListNthCGPWMcosns)

      if (allocated(pairlist%num_cg_pwmcosns_pre)) then
        if (size(pairlist%num_cg_pwmcosns_pre) == var_size) return
        deallocate(pairlist%num_cg_pwmcosns_pre,       &
            pairlist%num_cg_pwmcosns,                  &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_pwmcosns_pre(var_size), &
          pairlist%num_cg_pwmcosns(var_size),          &
          stat = alloc_stat)

    case (PairListEcqm)

      if (allocated(pairlist%ecqm_nb15_list)) then
        if (size(pairlist%ecqm_nb15_list) == var_size*2) return
        deallocate(pairlist%ecqm_nb15_list, stat = dealloc_stat)
      end if
      allocate(pairlist%ecqm_nb15_list(2, var_size), stat = alloc_stat)

    case (PairListAtomNobcGbsa)

      if (allocated(pairlist%num_all_calc)) then
        if (size(pairlist%num_all_calc) == var_size*nthread) return
        deallocate(pairlist%num_all_calc, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%num_all_calc(var_size,nthread), &
               stat = alloc_stat)

    case (PairListIntNobcGbsa)

      if (allocated(pairlist%all_calc_list)) then
        if (size(pairlist%all_calc_list) == var_size*nthread) return
        deallocate(pairlist%all_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%all_calc_list(var_size,nthread), &
               stat = alloc_stat)

    case (PairListNthreadsGbsa)

      if (allocated(pairlist%num_all_pre)) then
        if (size(pairlist%num_all_pre) == nthread) return
        deallocate(pairlist%num_all_pre, &
                   pairlist%num_all,     &
                   stat = dealloc_stat)
      end if
      allocate(pairlist%num_all_pre(nthread), &
               pairlist%num_all    (nthread), &
               stat = alloc_stat)

    ! =============
    ! CG model: PBC
    ! =============
    ! 
    case (PairListAtomPbcCGDNAexv)

      if (allocated(pairlist%num_cg_DNA_exv_calc)) then
        if (size(pairlist%num_cg_DNA_exv_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_DNA_exv_calc,                    &
            pairlist%cell_linked_list_cg_DNA,                       &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_DNA_exv_calc(var_size,nthread),      &
          pairlist%cell_linked_list_cg_DNA(var_size),               &
          stat = alloc_stat)

    case (PairListAtomPbcCGDNAbp)

      if (allocated(pairlist%num_cg_DNA_basepair_calc)) then
        if (size(pairlist%num_cg_DNA_basepair_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_DNA_basepair_calc,               &
            pairlist%cell_linked_list_cg_DNA_base,                  &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_DNA_basepair_calc(var_size,nthread), &
          pairlist%cell_linked_list_cg_DNA_base(var_size),          &
          stat = alloc_stat)

    case (PairListAtomPbcCGele)

      if (allocated(pairlist%num_cg_ele_calc)) then
        if (size(pairlist%num_cg_ele_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_ele_calc,                           &
            pairlist%cell_linked_list_cg_charged,                      &
            pairlist%cell_linked_list_cg_DNA_phos,                     &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_ele_calc(var_size,nthread),             &
          pairlist%cell_linked_list_cg_charged(var_size),              &
          pairlist%cell_linked_list_cg_DNA_phos(var_size2),            &
          stat = alloc_stat)

    case (PairListAtomPbcCGexv)

      if (allocated(pairlist%num_cg_exv_calc)) then
        if (size(pairlist%num_cg_exv_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_exv_calc,                        &
            pairlist%cell_linked_list_cg_all,                       &
            pairlist%cell_index_cg_all,                             &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_exv_calc(var_size,nthread),          &
          pairlist%cell_linked_list_cg_all(var_size),               &
          pairlist%cell_index_cg_all(var_size),                     &
          stat = alloc_stat)

    case (PairListAtomPbcCGPWMcos)

      if (allocated(pairlist%num_cg_pwmcos_calc)) then
        if (size(pairlist%num_cg_pwmcos_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_pwmcos_calc,                     &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_pwmcos_calc(var_size,nthread),       &
          stat = alloc_stat)
      
    case (PairListAtomPbcCGPWMcosns)

      if (allocated(pairlist%num_cg_pwmcosns_calc)) then
        if (size(pairlist%num_cg_pwmcosns_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_pwmcosns_calc,                     &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_pwmcosns_calc(var_size,nthread),       &
          stat = alloc_stat)

    case (PairListAtomPbcCGIDRHPS)

      if (allocated(pairlist%num_cg_IDR_HPS_calc)) then
        if (size(pairlist%num_cg_IDR_HPS_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_IDR_HPS_calc,                    &
            pairlist%cell_linked_list_cg_IDR_HPS,                   &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_IDR_HPS_calc(var_size,nthread),      &
          pairlist%cell_linked_list_cg_IDR_HPS(var_size),           &
          stat = alloc_stat)

    case (PairListAtomPbcCGIDRKH)

      if (allocated(pairlist%num_cg_IDR_KH_calc)) then
        if (size(pairlist%num_cg_IDR_KH_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_IDR_KH_calc,                    &
            pairlist%cell_linked_list_cg_IDR_KH,                   &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_IDR_KH_calc(var_size,nthread),      &
          pairlist%cell_linked_list_cg_IDR_KH(var_size),           &
          stat = alloc_stat)

    case (PairListAtomPbcCGKH)

      if (allocated(pairlist%num_cg_KH_calc)) then
        if (size(pairlist%num_cg_KH_calc) == var_size*nthread) return
        deallocate(pairlist%num_cg_KH_calc,                        &
            pairlist%cell_linked_list_cg_KH,                       &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_KH_calc(var_size,nthread),          &
          pairlist%cell_linked_list_cg_KH(var_size),               &
          stat = alloc_stat)

    ! ---------------------------
    ! numbering in each thread...
    ! ---------------------------
    ! 
    case (PairListNthreadsPbcCG)

      if (allocated(pairlist%num_cg_DNA_basepair_pre)) then
        if (size(pairlist%num_cg_DNA_basepair_pre) == nthread) return
        deallocate(pairlist%num_cg_DNA_basepair_pre,       &
            pairlist%num_cg_DNA_basepair,                  &
            pairlist%num_cg_DNA_exv_pre,                   &
            pairlist%num_cg_DNA_exv,                       &
            pairlist%num_cg_ele_pre,                       &
            pairlist%num_cg_ele,                           &
            pairlist%num_cg_IDR_HPS_pre,                   &
            pairlist%num_cg_IDR_HPS,                       &
            pairlist%num_cg_IDR_KH_pre,                    &
            pairlist%num_cg_IDR_KH,                        &
            pairlist%num_cg_KH_pre,                        &
            pairlist%num_cg_KH,                            &
            pairlist%num_cg_exv_pre,                       &
            pairlist%num_cg_exv,                           &
            pairlist%num_cg_pwmcos_pre,                    &
            pairlist%num_cg_pwmcos,                        &
            pairlist%num_cg_pwmcosns_pre,                  &
            pairlist%num_cg_pwmcosns,                      &
            stat = dealloc_stat)
      end if
      allocate(pairlist%num_cg_DNA_basepair_pre (nthread), &
          pairlist%num_cg_DNA_basepair          (nthread), &
          pairlist%num_cg_DNA_exv_pre           (nthread), &
          pairlist%num_cg_DNA_exv               (nthread), &
          pairlist%num_cg_ele_pre               (nthread), &
          pairlist%num_cg_ele                   (nthread), &
          pairlist%num_cg_IDR_HPS_pre           (nthread), &
          pairlist%num_cg_IDR_HPS               (nthread), &
          pairlist%num_cg_IDR_KH_pre            (nthread), &
          pairlist%num_cg_IDR_KH                (nthread), &
          pairlist%num_cg_KH_pre                (nthread), &
          pairlist%num_cg_KH                    (nthread), &
          pairlist%num_cg_exv_pre               (nthread), &
          pairlist%num_cg_exv                   (nthread), &
          pairlist%num_cg_pwmcos_pre            (nthread), &
          pairlist%num_cg_pwmcos                (nthread), &
          pairlist%num_cg_pwmcosns_pre          (nthread), &
          pairlist%num_cg_pwmcosns              (nthread), &
          stat = alloc_stat)

    ! ----------------
    ! the real list...
    ! ----------------
    ! 
    case (PairListPbcCGDNAbp)

      if (allocated(pairlist%cg_DNA_basepair_list)) then
        if (size(pairlist%cg_DNA_basepair_list) == var_size*nthread) return
        deallocate(pairlist%cg_DNA_basepair_list, stat = dealloc_stat)
      end if
      allocate(pairlist%cg_DNA_basepair_list(var_size,nthread), stat = alloc_stat)

    case (PairListPbcCGDNAexv)

      if (allocated(pairlist%cg_DNA_exv_list)) then
        if (size(pairlist%cg_DNA_exv_list) == var_size*nthread) return
        deallocate(pairlist%cg_DNA_exv_list, stat = dealloc_stat)
      end if
      allocate(pairlist%cg_DNA_exv_list(var_size,nthread), stat = alloc_stat)

    case (PairListPbcCGele)

      if (allocated(pairlist%cg_ele_list)) then
        if (size(pairlist%cg_ele_list) == var_size*nthread) return
        deallocate (pairlist%cg_ele_list,                 &
            pairlist%cg_ele_scaling_list,                 &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_ele_list(var_size,nthread),    &
          pairlist%cg_ele_scaling_list(var_size,nthread), &
          stat = alloc_stat)
      pairlist%cg_ele_scaling_list(:, :) = 1.0_wp

    case (PairListPbcCGIDRHPS)

      if (allocated(pairlist%cg_IDR_HPS_list)) then
        if (size(pairlist%cg_IDR_HPS_list) == var_size*nthread) return
        deallocate (pairlist%cg_IDR_HPS_list, stat = dealloc_stat)
      end if

      allocate(pairlist%cg_IDR_HPS_list(var_size,nthread), stat = alloc_stat)

    case (PairListPbcCGIDRKH)

      if (allocated(pairlist%cg_IDR_KH_list)) then
        if (size(pairlist%cg_IDR_KH_list) == var_size*nthread) return
        deallocate (pairlist%cg_IDR_KH_list, stat = dealloc_stat)
      end if

      allocate(pairlist%cg_IDR_KH_list(var_size,nthread), stat = alloc_stat)

    case (PairListPbcCGKH)

      if (allocated(pairlist%cg_KH_list)) then
        if (size(pairlist%cg_KH_list) == var_size*nthread) return
        deallocate (pairlist%cg_KH_list, &
            pairlist%cg_KH_model_list,   &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cg_KH_list(var_size,nthread), stat = alloc_stat)
      allocate(pairlist%cg_KH_model_list(var_size,nthread), stat = alloc_stat)

      pairlist%cg_KH_model_list(:, :) = 0

    case (PairListPbcCGexv)
      if (allocated(pairlist%cg_exv_list)) then
        if (size(pairlist%cg_exv_list) == var_size*nthread) return
        deallocate(pairlist%cg_exv_list, stat = dealloc_stat)
      end if

      allocate(pairlist%cg_exv_list(var_size,nthread), stat = alloc_stat)

    case (PairListPbcCGPWMcos)
      if (allocated(pairlist%cg_pwmcos_list)) then
        if (size(pairlist%cg_pwmcos_list) == var_size*nthread) return
        deallocate(pairlist%cg_pwmcos_list, stat = dealloc_stat)
      end if

      allocate(pairlist%cg_pwmcos_list(var_size,nthread), stat = alloc_stat)
      
    case (PairListPbcCGPWMcosns)
      if (allocated(pairlist%cg_pwmcosns_list)) then
        if (size(pairlist%cg_pwmcosns_list) == var_size*nthread) return
        deallocate(pairlist%cg_pwmcosns_list, stat = dealloc_stat)
      end if

      allocate(pairlist%cg_pwmcosns_list(var_size,nthread), stat = alloc_stat)

    ! ---------------------------
    ! cell head of linked list...
    ! ---------------------------
    !
    case (PairListCellsPbcCG)
      if (allocated(pairlist%cell_head_index_cg_all)) then
        deallocate(pairlist%cell_head_index_cg_all,       &
            pairlist%cell_head_index_cg_DNA,              &
            pairlist%cell_head_index_cg_DNA_phos,         &
            pairlist%cell_head_index_cg_DNA_base,         &
            pairlist%cell_head_index_cg_IDR_HPS,          &
            pairlist%cell_head_index_cg_IDR_KH,           &
            pairlist%cell_head_index_cg_KH,               &
            pairlist%cell_head_index_cg_charged,          &
            stat = dealloc_stat)
      end if

      allocate(pairlist%cell_head_index_cg_all(var_size), &
          pairlist%cell_head_index_cg_DNA(var_size),      &
          pairlist%cell_head_index_cg_DNA_phos(var_size), &
          pairlist%cell_head_index_cg_DNA_base(var_size), &
          pairlist%cell_head_index_cg_IDR_HPS(var_size),  &
          pairlist%cell_head_index_cg_IDR_KH(var_size),   &
          pairlist%cell_head_index_cg_KH(var_size),       &
          pairlist%cell_head_index_cg_charged(var_size),  &
          stat = alloc_stat)

    case  default

      call error_msg('Alloc_PairList> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist2
  !> @brief        allocate pairlist information
  !! @authors      MK
  !! @param[inout] pairlist  : pairlist information
  !! @param[in]    variable  : allocatable variables
  !! @param[in]    var_size1 : size of variables (1st dimension)
  !! @param[in]    var_size2 : size of variables (2nd dimension)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist2(pairlist, variable, var_size1, var_size2)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: nthread, omp_get_max_threads


    alloc_stat   = 0
    dealloc_stat = 0

    ! get number of threads
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! allocate selected variables
    !
    select case (variable)

    case (PairListPbcSoluteSolute)

      if (allocated(pairlist%table%nb15_calc_list)) then
        if (size(pairlist%table%nb15_calc_list) == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_list, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_list(var_size1,var_size2), &
               stat = alloc_stat)

    case (PairListPbcSoluteWater)

      if (allocated(pairlist%table%nb15_calc_listw)) then
        if (size(pairlist%table%nb15_calc_listw) == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_listw, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_listw(var_size1,var_size2), &
               stat = alloc_stat)

    case (PairListPbcWaterWater)

      if (allocated(pairlist%table%nb15_calc_list_water)) then
        if (size(pairlist%table%nb15_calc_list_water)               &
            == var_size1*var_size2) return
        deallocate(pairlist%table%nb15_calc_list_water, &
                   stat = dealloc_stat)
      end if

      allocate(pairlist%table%nb15_calc_list_water(var_size1,var_size2), &
               stat = alloc_stat)

    case  default

      call error_msg('Alloc_PairList2> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_pairlist2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist
  !> @brief        deallocate pairlist information
  !! @authors      YS, TI, JJ
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

    select case (variable)

    case (PairListAtomNobc)

      if (allocated(pairlist%num_nb15_calc)) then
        deallocate (pairlist%num_nb15_calc,   &
                    stat = dealloc_stat)
      end if

    case (PairListIntNobc)

      if (allocated(pairlist%nb15_calc_list)) then
        deallocate (pairlist%nb15_calc_list, &
                    stat = dealloc_stat)
      end if
      
    case (PairListCGDNABP)

      if (allocated(pairlist%cg_DNA_basepair_list)) then
        deallocate (pairlist%cg_DNA_basepair_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGDNAexv)

      if (allocated(pairlist%cg_DNA_exv_list)) then
        deallocate (pairlist%cg_DNA_exv_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGele)

      if (allocated(pairlist%cg_ele_list)) then
        deallocate (pairlist%cg_ele_list, &
            pairlist%cg_ele_scaling_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGIDRHPS)

      if (allocated(pairlist%cg_IDR_HPS_list)) then
        deallocate (pairlist%cg_IDR_HPS_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGIDRKH)

      if (allocated(pairlist%cg_IDR_KH_list)) then
        deallocate (pairlist%cg_IDR_KH_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGKH)

      if (allocated(pairlist%cg_KH_list)) then
        deallocate (pairlist%cg_KH_list, &
            pairlist%cg_KH_model_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGexv)

      if (allocated(pairlist%cg_exv_list)) then
        deallocate (pairlist%cg_exv_list, &
            stat = dealloc_stat)
      end if

    case (PairListCGPWMcos)

      if (allocated(pairlist%cg_pwmcos_list)) then
        deallocate (pairlist%cg_pwmcos_list, &
            stat = dealloc_stat)
      end if
      
    case (PairListCGPWMcosns)

      if (allocated(pairlist%cg_pwmcosns_list)) then
        deallocate (pairlist%cg_pwmcosns_list, &
            stat = dealloc_stat)
      end if

    case (PairListPbcSolute)

      if (allocated(pairlist%table%num_nb15_calc)) then
        deallocate (pairlist%table%num_nb15_calc,    &
                    pairlist%table%num_nb15_calcw,   &
                    pairlist%table%cell_linked_list, &
                    pairlist%table%atom_cell_index,  &
                    stat = dealloc_stat)
      end if

    case (PairListPbcWater)

      if (allocated(pairlist%table%num_nb15_calc_water)) then
        deallocate (pairlist%table%num_nb15_calc_water,    &
                    pairlist%table%cell_linked_list_water, &
                    pairlist%table%atom_cell_index_water,  &
                    stat = dealloc_stat)
      end if

    case (PairListPbcSoluteSolute)

      if (allocated(pairlist%table%nb15_calc_list)) then
        deallocate( pairlist%table%nb15_calc_list, &
                    stat = dealloc_stat)
      end if

    case (PairListPbcSoluteWater)

      if (allocated(pairlist%table%nb15_calc_listw)) then
        deallocate (pairlist%table%nb15_calc_listw, &
                    stat = dealloc_stat)
      end if

    case (PairListPbcWaterWater)

      if (allocated(pairlist%table%nb15_calc_list_water)) then
        deallocate (pairlist%table%nb15_calc_list_water, &
                    stat = dealloc_stat)
      end if

    case (PairListNthreads)

      if (allocated(pairlist%num_nb15_pre)) then
        deallocate (pairlist%num_nb15_pre,     &
                    pairlist%num_nb15,         &
                    stat = dealloc_stat)
      end if

    case (PairListEcqm)

      if (allocated(pairlist%ecqm_nb15_list)) then
        deallocate (pairlist%ecqm_nb15_list, stat = dealloc_stat)
      end if

    case (PairListAtomNobcGbsa)

      if (allocated(pairlist%num_all_calc)) then
        deallocate (pairlist%num_all_calc,   &
                    stat = dealloc_stat)
      end if

    case (PairListIntNobcGbsa)

      if (allocated(pairlist%all_calc_list)) then
        deallocate (pairlist%all_calc_list, &
                    stat = dealloc_stat)
      end if

    case (PairListNthreadsGbsa)

      if (allocated(pairlist%num_all_pre)) then
        deallocate (pairlist%num_all_pre,     &
                    pairlist%num_all,         &
                    stat = dealloc_stat)
      end if
      
    case (PairListAtomNobcCG)

      if (allocated(pairlist%num_cg_DNA_basepair_calc)) then
        deallocate(pairlist%num_cg_DNA_basepair_calc, &
            pairlist%num_cg_DNA_exv_calc,             &
            pairlist%num_cg_ele_calc,                 &
            pairlist%num_cg_exv_calc,                 &
            pairlist%num_cg_IDR_HPS_calc,             &
            pairlist%num_cg_IDR_KH_calc,              &
            pairlist%num_cg_KH_calc,                  &
            stat = dealloc_stat)
      end if
      
    case (PairListNthreadsCG)

      if (allocated(pairlist%num_cg_DNA_basepair_pre)) then
        deallocate (pairlist%num_cg_DNA_basepair_pre, &
            pairlist%num_cg_DNA_basepair,             &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_DNA_exv_pre))      then
        deallocate (pairlist%num_cg_DNA_exv_pre,      &
            pairlist%num_cg_DNA_exv,                  &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_ele_pre))          then
        deallocate (pairlist%num_cg_ele_pre,          &
            pairlist%num_cg_ele,                      &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_IDR_HPS_pre))      then
        deallocate (pairlist%num_cg_IDR_HPS_pre,      &
            pairlist%num_cg_IDR_HPS,                  &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_IDR_KH_pre))      then
        deallocate (pairlist%num_cg_IDR_KH_pre,      &
            pairlist%num_cg_IDR_KH,                  &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_KH_pre))          then
        deallocate (pairlist%num_cg_KH_pre,          &
            pairlist%num_cg_KH,                      &
            stat = dealloc_stat)
      end if
      if (allocated(pairlist%num_cg_exv_pre))          then
        deallocate (pairlist%num_cg_exv_pre,          &
            pairlist%num_cg_exv,                      &
            stat = dealloc_stat)
      end if
    case (PairListNthCGPWMcos)
      if (allocated(pairlist%num_cg_pwmcos))       then
        deallocate (pairlist%num_cg_pwmcos,       &
            stat = dealloc_stat)
      end if
    case (PairListNthCGPWMcosns)
      if (allocated(pairlist%num_cg_pwmcosns))       then
        deallocate (pairlist%num_cg_pwmcosns,       &
            stat = dealloc_stat)
      end if

    ! =============
    ! CG model: PBC
    ! =============
    ! 
    case (PairListAtomPbcCGDNAexv)

      if (allocated(pairlist%num_cg_DNA_exv_calc)) then
        deallocate(pairlist%num_cg_DNA_exv_calc,      &
            pairlist%cell_linked_list_cg_DNA,         &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGDNAbp)

      if (allocated(pairlist%num_cg_DNA_basepair_calc)) then
        deallocate(pairlist%num_cg_DNA_basepair_calc, &
            pairlist%cell_linked_list_cg_DNA_base,    &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGele)

      if (allocated(pairlist%num_cg_ele_calc)) then
        deallocate(pairlist%num_cg_ele_calc,          &
            pairlist%cell_linked_list_cg_charged,     &
            pairlist%cell_linked_list_cg_DNA_phos,    &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGexv)

      if (allocated(pairlist%num_cg_exv_calc)) then
        deallocate(pairlist%num_cg_exv_calc,          &
            pairlist%cell_linked_list_cg_all,         &
            pairlist%cell_index_cg_all,               &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGPWMcos)

      if (allocated(pairlist%num_cg_pwmcos_calc)) then
        deallocate(pairlist%num_cg_pwmcos_calc,       &
            stat = dealloc_stat)
      end if
      
    case (PairListAtomPbcCGPWMcosns)

      if (allocated(pairlist%num_cg_pwmcosns_calc)) then
        deallocate(pairlist%num_cg_pwmcosns_calc,       &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGIDRHPS)

      if (allocated(pairlist%num_cg_IDR_HPS_calc)) then
        deallocate(pairlist%num_cg_IDR_HPS_calc,      &
            pairlist%cell_linked_list_cg_IDR_HPS,     &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGIDRKH)

      if (allocated(pairlist%num_cg_IDR_KH_calc)) then
        deallocate(pairlist%num_cg_IDR_KH_calc,      &
            pairlist%cell_linked_list_cg_IDR_KH,     &
            stat = dealloc_stat)
      end if

    case (PairListAtomPbcCGKH)

      if (allocated(pairlist%num_cg_KH_calc)) then
        deallocate(pairlist%num_cg_KH_calc,      &
            pairlist%cell_linked_list_cg_KH,     &
            stat = dealloc_stat)
      end if

    case (PairListNthreadsPbcCG)

      if (allocated(pairlist%num_cg_DNA_basepair_pre)) then
        deallocate(pairlist%num_cg_DNA_basepair_pre,  &
            pairlist%num_cg_DNA_basepair,             &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_DNA_exv_pre)) then
        deallocate(pairlist%num_cg_DNA_exv_pre,       &
            pairlist%num_cg_DNA_exv,                  &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_ele_pre)) then
        deallocate(pairlist%num_cg_ele_pre,           &
            pairlist%num_cg_ele,                      &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_IDR_HPS_pre)) then
        deallocate(pairlist%num_cg_IDR_HPS_pre,       &
            pairlist%num_cg_IDR_HPS,                  &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_IDR_KH_pre)) then
        deallocate(pairlist%num_cg_IDR_KH_pre,       &
            pairlist%num_cg_IDR_KH,                  &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_KH_pre)) then
        deallocate(pairlist%num_cg_KH_pre,       &
            pairlist%num_cg_KH,                  &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_exv_pre)) then
        deallocate(pairlist%num_cg_exv_pre,           &
            pairlist%num_cg_exv,                      &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_pwmcos)) then
        deallocate(pairlist%num_cg_pwmcos_pre,        &
            pairlist%num_cg_pwmcos,                   &
            stat = dealloc_stat)
      end if

      if (allocated(pairlist%num_cg_pwmcosns)) then
        deallocate(pairlist%num_cg_pwmcosns_pre,        &
            pairlist%num_cg_pwmcosns,                   &
            stat = dealloc_stat)
      end if

    ! ----------------
    ! the real list...
    ! ----------------
    ! 
    case (PairListPbcCGDNAbp)

      if (allocated(pairlist%cg_DNA_basepair_list)) then
        deallocate(pairlist%cg_DNA_basepair_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGDNAexv)

      if (allocated(pairlist%cg_DNA_exv_list)) then
        deallocate(pairlist%cg_DNA_exv_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGele)

      if (allocated(pairlist%cg_ele_list)) then
        deallocate (pairlist%cg_ele_list, &
            pairlist%cg_ele_scaling_list, &
            stat = dealloc_stat)
      end if

    case (PairListPbcCGIDRHPS)

      if (allocated(pairlist%cg_IDR_HPS_list)) then
        deallocate (pairlist%cg_IDR_HPS_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGIDRKH)

      if (allocated(pairlist%cg_IDR_KH_list)) then
        deallocate (pairlist%cg_IDR_KH_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGKH)

      if (allocated(pairlist%cg_KH_list)) then
        deallocate (pairlist%cg_KH_list, &
            pairlist%cg_KH_model_list,   &
            stat = dealloc_stat)
      end if

    case (PairListPbcCGexv)
      if (allocated(pairlist%cg_exv_list)) then
        deallocate(pairlist%cg_exv_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGPWMcos)
      if (allocated(pairlist%cg_pwmcos_list)) then
        deallocate(pairlist%cg_pwmcos_list, stat = dealloc_stat)
      end if

    case (PairListPbcCGPWMcosns)
      if (allocated(pairlist%cg_pwmcosns_list)) then
        deallocate(pairlist%cg_pwmcosns_list, stat = dealloc_stat)
      end if

    ! ---------------------------
    ! cell head of linked list...
    ! ---------------------------
    !
    case (PairListCellsPbcCG)
      if (allocated(pairlist%cell_head_index_cg_all)) then
        deallocate(pairlist%cell_head_index_cg_all, &
            pairlist%cell_head_index_cg_DNA,        &
            pairlist%cell_head_index_cg_DNA_phos,   &
            pairlist%cell_head_index_cg_DNA_base,   &
            pairlist%cell_head_index_cg_IDR_HPS,    &
            pairlist%cell_head_index_cg_IDR_KH,     &
            pairlist%cell_head_index_cg_KH,         &
            pairlist%cell_head_index_cg_charged,    &
            stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_PairList> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist_all
  !> @brief        deallocate all pairlist information
  !! @authors      YS, TI
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist_all(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    call dealloc_pairlist(pairlist, PairListAtomNobc)
    call dealloc_pairlist(pairlist, PairListIntNobc)
    call dealloc_pairlist(pairlist, PairListPbcSolute)
    call dealloc_pairlist(pairlist, PairListPbcWater)
    call dealloc_pairlist(pairlist, PairListPbcSoluteSolute)
    call dealloc_pairlist(pairlist, PairListPbcSoluteWater)
    call dealloc_pairlist(pairlist, PairListPbcWaterWater)
    call dealloc_pairlist(pairlist, PairListNthreads)
    call dealloc_pairlist(pairlist, PairListEcqm)
    call dealloc_pairlist(pairlist, PairListAtomNobcGbsa)
    call dealloc_pairlist(pairlist, PairListIntNobcGbsa)
    call dealloc_pairlist(pairlist, PairListNthreadsGbsa)
    call dealloc_pairlist(pairlist, PairListNthreadsCG)
    call dealloc_pairlist(pairlist, PairListAtomNobcCG)
    call dealloc_pairlist(pairlist, PairListCGDNABP)
    call dealloc_pairlist(pairlist, PairListCGDNAexv)
    call dealloc_pairlist(pairlist, PairListCGele)
    call dealloc_pairlist(pairlist, PairListCGexv)
    call dealloc_pairlist(pairlist, PairListCGPWMcos)
    call dealloc_pairlist(pairlist, PairListCGPWMcosns)
    call dealloc_pairlist(pairlist, PairListCGIDRHPS)
    call dealloc_pairlist(pairlist, PairListCGIDRKH)
    call dealloc_pairlist(pairlist, PairListCGKH)
    call dealloc_pairlist(pairlist, PairListNthCGPWMcos)
    call dealloc_pairlist(pairlist, PairListNthCGPWMcosns)
    ! 
    call dealloc_pairlist(pairlist, PairListAtomPbcCGexv)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGele)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGDNAexv)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGDNAbp)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGPWMcos)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGPWMcosns)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGIDRHPS)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGIDRKH)
    call dealloc_pairlist(pairlist, PairListAtomPbcCGKH)
    call dealloc_pairlist(pairlist, PairListNthreadsPbcCG)
    call dealloc_pairlist(pairlist, PairListPbcCGDNAbp)
    call dealloc_pairlist(pairlist, PairListPbcCGDNAexv)
    call dealloc_pairlist(pairlist, PairListPbcCGexv)
    call dealloc_pairlist(pairlist, PairListPbcCGele)
    call dealloc_pairlist(pairlist, PairListPbcCGPWMcos)
    call dealloc_pairlist(pairlist, PairListPbcCGPWMcosns)
    call dealloc_pairlist(pairlist, PairListPbcCGIDRHPS)
    call dealloc_pairlist(pairlist, PairListPbcCGIDRKH)
    call dealloc_pairlist(pairlist, PairListPbcCGKH)
    call dealloc_pairlist(pairlist, PairListCellsPbcCG)

    return

  end subroutine dealloc_pairlist_all

end module at_pairlist_str_mod
