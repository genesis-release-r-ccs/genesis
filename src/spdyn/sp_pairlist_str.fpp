!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_str_mod
!> @brief   structure of pairlist information
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_str_mod

  use iso_c_binding

  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pairlist
    real(wp)                        :: pairlistdist
    integer,            allocatable :: num_nb15_calc(:,:)
    integer,            allocatable :: num_nb15_calc_solute(:)
    integer,            allocatable :: num_nb15_calc_water(:)
    integer,            allocatable :: num_nb15_calc_water_p(:)
    integer,            allocatable :: num_nb15_calc1(:,:)
    integer,            allocatable :: nb15_calc_list(:,:)
    integer,            allocatable :: nb15_calc_list_fugaku(:,:,:)
    integer,            allocatable :: nb15_calc_list_solute(:,:)
    integer,            allocatable :: nb15_calc_list_water(:,:)
    integer,            allocatable :: nb15_calc_list_water_p(:,:)
    integer,            allocatable :: nb15_calc_list1(:,:)
    integer,            allocatable :: nb15_cell(:)
    integer,            allocatable :: nb15_list(:,:)
    ! for GPU
    integer,            allocatable :: univ_ij_load(:)         !  (ij)
#ifndef PGICUDA
    integer(1),         allocatable :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer(1),         allocatable :: pack_univ_mask2(:)      !  (ix*iy, ij)
    integer,            allocatable :: univ_ij_sort_list(:)    !  (ij)
    integer,            allocatable :: univ_ix_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_ix_list(:,:)       !  (iix, ij)
    integer,            allocatable :: univ_iy_natom(:)        !  (ij)
    integer(1),         allocatable :: univ_iy_list(:,:)       !  (iiy, ij)
#else
    integer(1), allocatable, pinned :: univ_mask2(:,:)         !  (ix*iy, ij)
    integer(1), allocatable, pinned :: pack_univ_mask2(:)      !  (ix*iy, ij)
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
    ! for FEP
    integer,            allocatable :: num_nb15_calc_fep(:,:)
    integer,            allocatable :: num_nb15_calc1_fep(:,:)
    integer,            allocatable :: nb15_calc_list_fep(:,:)
    integer,            allocatable :: nb15_calc_list1_fep(:,:)
    integer,            allocatable :: nb15_cell_fep(:)
    integer,            allocatable :: nb15_list_fep(:,:)
  end type s_pairlist

  ! parameters for allocatable variables
  integer,        public, parameter :: PairListGeneric     = 1
  integer,        public, parameter :: PairListGPU         = 2
  integer,        public, parameter :: PairListIntelGeneric= 3
  integer,        public, parameter :: PairListFugaku      = 4
  integer,        public, parameter :: PairListFEP         = 6

  ! variables for maximum numbers in one cell
  integer,        public            :: MaxNb15        = 15000
  integer,        public            :: MaxNb15_Fugaku = 3000
  integer,        public            :: MaxNb15Water   = 1500
  integer,        public            :: MaxNb15_S      = 1000
  integer,        public            :: MaxNb15_W      = 1000
  integer,        public            :: MaxNb15_W_P    = 1000
  ! FEP: MaxNb15_fep will be changed depending on the numbe of perturbed atoms.
  integer,        public            :: MaxNb15_fep    = 15000

  ! subroutines
  public  :: init_pairlist
  public  :: alloc_pairlist
  public  :: dealloc_pairlist
  public  :: dealloc_pairlist_all

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


    pairlist%pairlistdist = 0.0_wp

    return

  end subroutine init_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist
  !> @brief        allocate pairlist information
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !! @param[in]    variable : allocatable variables
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist(pairlist, variable, var_size)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (PairListGeneric)

      if (allocated(pairlist%num_nb15_calc1)) then
        if (size(pairlist%num_nb15_calc1) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1,  &
                     pairlist%num_nb15_calc,   &
                     pairlist%nb15_calc_list1, &
                     pairlist%nb15_calc_list,  &
                     pairlist%nb15_cell,       &
                     pairlist%nb15_list,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1)) &
        allocate(pairlist%num_nb15_calc1 (MaxAtom, var_size), &
                 pairlist%num_nb15_calc  (MaxAtom, maxcell),  &
                 pairlist%nb15_calc_list1(MaxNb15, var_size), &
                 pairlist%nb15_calc_list (MaxNb15, maxcell),  &
                 pairlist%nb15_cell      (maxcell),           &
                 pairlist%nb15_list      (MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%num_nb15_calc1 (1:MaxAtom, 1:var_size)  = 0
      pairlist%num_nb15_calc  (1:MaxAtom, 1:maxcell)   = 0

    case (PairListFugaku)

      if (allocated(pairlist%num_nb15_calc)) then
        if (size(pairlist%num_nb15_calc) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc,          &
                     pairlist%nb15_calc_list_fugaku,  &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1))           &
        allocate(pairlist%num_nb15_calc  (MaxAtom*var_size,1),&
                 pairlist%nb15_calc_list_fugaku               &
                      (MaxNb15_Fugaku, MaxAtom*var_size,1),   &
                 stat = alloc_stat)

      pairlist%num_nb15_calc  (1:MaxAtom*var_size,1)   = 0

    case (PairListIntelGeneric)

      if (allocated(pairlist%nb15_cell)) then
        deallocate(pairlist%nb15_cell,       &
                   pairlist%nb15_list,       &
                   stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%nb15_cell))                  &
        allocate(pairlist%nb15_cell      (maxcell),             &
                 pairlist%nb15_list      (2*MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%nb15_cell(1:maxcell)  = 0
      pairlist%nb15_list(1:2*MaxAtom, 1:maxcell)   = 0

    case (PairListGPU)

      if (allocated(pairlist%univ_ij_load)) then
        if (size(pairlist%univ_ij_load) /= univ_maxcell1) then
#ifdef USE_GPU
          call unset_pinned_memory(pairlist%univ_ij_load)
          call unset_pinned_memory(pairlist%univ_ij_sort_list)
          call unset_pinned_memory(pairlist%univ_ix_natom)   
          call unset_pinned_memory(pairlist%univ_ix_list)   
          call unset_pinned_memory(pairlist%univ_iy_natom)   
          call unset_pinned_memory(pairlist%univ_iy_list)   
#endif
          deallocate(pairlist%univ_ij_load,      &
                     pairlist%univ_ij_sort_list, &
                     pairlist%univ_ix_natom,     &
                     pairlist%univ_ix_list,      &
                     pairlist%univ_iy_natom,     &
                     pairlist%univ_iy_list,      &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(pairlist%univ_ij_load)) then
        allocate(pairlist%univ_ij_load(univ_maxcell1),            &
                 pairlist%univ_ij_sort_list(univ_maxcell1),       &
                 pairlist%univ_ix_natom(univ_maxcell1),           &
                 pairlist%univ_ix_list(MaxAtom, univ_maxcell1),   &
                 pairlist%univ_iy_natom(univ_maxcell1),           &
                 pairlist%univ_iy_list(MaxAtom, univ_maxcell1),   &
                 stat = alloc_stat)
#ifdef USE_GPU
        call set_pinned_memory(pairlist%univ_ij_sort_list, univ_maxcell1*4)
        call set_pinned_memory(pairlist%univ_ix_natom, univ_maxcell1*4)
        call set_pinned_memory(pairlist%univ_iy_natom, univ_maxcell1*4)
        call set_pinned_memory(pairlist%univ_ix_list, MaxAtom*univ_maxcell1)
        call set_pinned_memory(pairlist%univ_iy_list, MaxAtom*univ_maxcell1)
#endif
      end if

      ! initialization
      pairlist%univ_ij_load(:)      = 0
      pairlist%univ_ij_sort_list(:) = 0
      pairlist%univ_ix_natom(:)     = 0
      pairlist%univ_ix_list(:,:)    = 0
      pairlist%univ_iy_natom(:)     = 0
      pairlist%univ_iy_list(:,:)    = 0

    case (PairListFEP)

      if (allocated(pairlist%num_nb15_calc1_fep)) then
        if (size(pairlist%num_nb15_calc1_fep) /= var_size*MaxAtom) &
          deallocate(pairlist%num_nb15_calc1_fep,  &
                     pairlist%num_nb15_calc_fep,   &
                     pairlist%nb15_calc_list1_fep, &
                     pairlist%nb15_calc_list_fep,  &
                     pairlist%nb15_cell_fep,       &
                     pairlist%nb15_list_fep,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(pairlist%num_nb15_calc1_fep)) &
        allocate(pairlist%num_nb15_calc1_fep (MaxAtom, var_size), &
                 pairlist%num_nb15_calc_fep  (MaxAtom, maxcell),  &
                 pairlist%nb15_calc_list1_fep(MaxNb15_fep, var_size), &
                 pairlist%nb15_calc_list_fep (MaxNb15_fep, maxcell),  &
                 pairlist%nb15_cell_fep      (maxcell),           &
                 pairlist%nb15_list_fep      (MaxAtom, maxcell),  &
                 stat = alloc_stat)

      pairlist%num_nb15_calc1_fep (1:MaxAtom, 1:var_size)  = 0
      pairlist%num_nb15_calc_fep  (1:MaxAtom, 1:maxcell)   = 0

    end select

    if (alloc_stat /=0)   call error_msg_alloc
    if (dealloc_stat /=0) call error_msg_dealloc

    return

  end subroutine alloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_pairlist
  !> @brief        deallocate pairlist information
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

    case (PairListGeneric)

      if (allocated(pairlist%num_nb15_calc1)) then
        deallocate(pairlist%num_nb15_calc1,  &
                   pairlist%num_nb15_calc,   &
                   pairlist%nb15_calc_list1, &
                   pairlist%nb15_calc_list,  &
                   pairlist%nb15_cell,       &
                   pairlist%nb15_list,       &
                   stat = dealloc_stat)
      end if

    case (PairListGPU)

      if (allocated(pairlist%univ_ij_load)) then
        deallocate(pairlist%univ_ij_load,      &
                   pairlist%univ_ij_sort_list, &
                   pairlist%univ_ix_natom,     &
                   pairlist%univ_ix_list,      &
                   pairlist%univ_iy_natom,     &
                   pairlist%univ_iy_list,      &
                   stat = dealloc_stat)
      end if

    case (PairListIntelGeneric)

      if (allocated(pairlist%nb15_cell)) then
        deallocate(pairlist%nb15_cell,  &
                   pairlist%nb15_list,  &
                   stat = dealloc_stat)
      end if

    case (PairListFugaku)

      if (allocated(pairlist%num_nb15_calc)) then
        deallocate(pairlist%num_nb15_calc,          &
                   pairlist%nb15_calc_list_fugaku,  &
                   stat = dealloc_stat)
      end if

    case (PairListFEP)

      if (allocated(pairlist%num_nb15_calc1_fep)) then
        deallocate(pairlist%num_nb15_calc1_fep,  &
                   pairlist%num_nb15_calc_fep,   &
                   pairlist%nb15_calc_list1_fep, &
                   pairlist%nb15_calc_list_fep,  &
                   pairlist%nb15_cell_fep,       &
                   pairlist%nb15_list_fep,       &
                   stat = dealloc_stat)
      end if
 
    end select

    if (dealloc_stat /=0) call error_msg_dealloc

    return

  end subroutine dealloc_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pairlist_all
  !> @brief        deallocate all pairlist information
  !! @authors      JJ
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pairlist_all(pairlist)

    ! formal arguments
    type(s_pairlist),        intent(inout) :: pairlist


    call dealloc_pairlist(pairlist, PairListGeneric)
    call dealloc_pairlist(pairlist, PairListIntelGeneric)
    call dealloc_pairlist(pairlist, PairListGPU)
    call dealloc_pairlist(pairlist, PairListFEP)

    return

  end subroutine dealloc_pairlist_all

end module sp_pairlist_str_mod
