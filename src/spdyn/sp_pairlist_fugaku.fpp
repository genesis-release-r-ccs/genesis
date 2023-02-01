!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_fugaku_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_fugaku_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: update_pairlist_pbc_fugaku
  ! FEP
  public  :: update_pairlist_pbc_fugaku_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_fugaku
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_fugaku(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2_work(1:MaxAtom)
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: num_nb15
    integer                   :: id, omp_get_thread_num

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)


    ncell             => domain%num_cell_local
    nboundary         => domain%num_cell_boundary
    natom             => domain%num_atom
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    coord_pbc         => domain%coord_pbc
    cell_pairlist1    => domain%cell_pairlist1
    ncell_pairlist1   => domain%ncell_pairlist1
    istcell_pairlist1 => domain%istcell_pairlist1
    ncell_pairlist2   => domain%ncell_pairlist2
    istcell_pairlist2 => domain%istcell_pairlist2

    num_nb15_calc     => pairlist%num_nb15_calc
    nb15_calc_list    => pairlist%nb15_calc_list_fugaku
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    !$omp parallel default(shared)                                     &
    !$omp private(num_add_prefetch, k, start_i, start_j, ixx, iyy)     &
    !$omp private(id, i, j, ij, ix, iy, num_nb15, rtmp, dij, rij2_work)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell+nboundary, nthread
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = start_i + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    k = domain%num_atom_domain
    num_nb15_calc(1:k,1)  = 0

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      start_i = start_atom(i)
#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif
      do ix = 1, natom(i)-1

        num_nb15 = 0
        ixx = start_i + ix

        do iy = ix+1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            num_nb15 = num_nb15 + 1
            iyy = start_i + iy
#ifdef DEBUG
            if (num_nb15 > MaxNb15_Fugaku) &
              call error_msg(              &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
            nb15_calc_list(num_nb15,ixx,1) = iyy

          end if

        end do

        num_nb15_calc(ixx,1) = num_nb15

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = id+1, ncell+nboundary, nthread

      if (ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif
       
        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15 = num_nb15_calc(ixx,1)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2 .and. &
                exclusion_mask(iy,ix,ij) == 1) then
              iyy = start_j + iy
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 exceeds MaxNb15_Fugaku')
#endif
              nb15_calc_list(num_nb15,ixx,1) = iyy
            end if

          end do

          num_nb15_calc(ixx,1) = num_nb15

        end do
      end do

    end do

    !$omp barrier

    do i = id+1, ncell+nboundary, nthread

      if (ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15 = num_nb15_calc(ixx,1)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2) then
              iyy = start_j + iy
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ixx,1) = iyy
            end if

          end do

          num_nb15_calc(ixx,1) = num_nb15

        end do
      end do

    end do

    !$omp barrier

!   do i = id+1, ncell, nthread
!     start_i = start_atom(i)
!     do ix = 1, natom(i)-1

!       ixx = start_i + ix
!       ! prefetch size:64
!       !
!       num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ixx,1), &
!                                 num_nb15_calc(ixx+1,1))
!       do k = 1, num_add_prefetch
!         nb15_calc_list(num_nb15_calc(ixx,1)+k,ixx,1) = nb15_calc_list(k,ixx+1,1)
!       end do
!     end do
!     ix = natom(i)
!     ixx = start_i + ix
!     num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ixx,1))
!     if (num_nb15_calc(ixx,1) > 0) then
!     do k = 1, num_add_prefetch
!       nb15_calc_list(num_nb15_calc(ixx,1)+k,ixx,1) &
!             = nb15_calc_list(num_nb15_calc(ixx,1),ixx,1)
!     end do
!     end if
!   end do

    !$omp end parallel

    return

  end subroutine update_pairlist_pbc_fugaku

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_fugaku_fep
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition for FEP
  !! @authors      HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_fugaku_fep(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2_work(1:MaxAtom)
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: num_nb15
    integer                   :: id, omp_get_thread_num

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary
    integer,          pointer,contiguous :: nb15_calc_list(:,:,:)
    integer,          pointer,contiguous :: num_nb15_calc(:,:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)

    ! FEP
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_cell_fep(:), nb15_list_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer                   :: num_nb15_pre_fep, num_nb15_fep
    integer                   :: num_nb15_cell_fep

    ncell             => domain%num_cell_local
    nboundary         => domain%num_cell_boundary
    natom             => domain%num_atom
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    coord_pbc         => domain%coord_pbc
    cell_pairlist1    => domain%cell_pairlist1
    ncell_pairlist1   => domain%ncell_pairlist1
    istcell_pairlist1 => domain%istcell_pairlist1
    ncell_pairlist2   => domain%ncell_pairlist2
    istcell_pairlist2 => domain%istcell_pairlist2

    num_nb15_calc     => pairlist%num_nb15_calc
    nb15_calc_list    => pairlist%nb15_calc_list_fugaku
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    ! FEP
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep

    !$omp parallel default(shared)                                      &
    !$omp private(num_add_prefetch, k, start_i, start_j, ixx, iyy)      &
    !$omp private(id, i, j, ij, ix, iy, num_nb15, rtmp, dij, rij2_work, &
    !$omp         fg1, fg2, num_nb15_fep, num_nb15_pre_fep,             &
    !$omp         num_nb15_cell_fep)                                    
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell+nboundary, nthread
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = start_i + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    k = domain%num_atom_domain
    num_nb15_calc(1:k,1)  = 0

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      start_i = start_atom(i)
#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0

      do ix = 1, natom(i)-1

        num_nb15 = 0
        ixx = start_i + ix

        do iy = ix+1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            ! FEP
            fg1 = domain%fepgrp(ix,i)
            fg2 = domain%fepgrp(iy,i)
            if (enefunc%fepgrp_nonb(fg1,fg2) == 0) then
              ! FEP: Exclude partA-partB interactions
              cycle
            else if (enefunc%fepgrp_nonb(fg1,fg2) /= 5) then
              ! FEP: perturbed
              num_nb15_fep = num_nb15_fep + 1
              nb15_calc_list1_fep(num_nb15_fep,i) = iy
              cycle
            end if

            num_nb15 = num_nb15 + 1
            iyy = start_i + iy
#ifdef DEBUG
            if (num_nb15 > MaxNb15_Fugaku) &
              call error_msg(              &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
            nb15_calc_list(num_nb15,ixx,1) = iyy

          end if

        end do

        num_nb15_calc(ixx,1) = num_nb15

        ! FEP
        num_nb15_calc1_fep(ix,i) = num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = id+1, ncell+nboundary, nthread

      if (ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif

        ! FEP
        num_nb15_fep = 0
        num_nb15_pre_fep = 0
        num_nb15_cell_fep = 0
       
        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15 = num_nb15_calc(ixx,1)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2 .and. &
                exclusion_mask(iy,ix,ij) == 1) then

              ! FEP
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if (enefunc%fepgrp_nonb(fg1,fg2) == 0) then
                ! FEP: Exclude partA-partB interactions
                cycle
              else if (enefunc%fepgrp_nonb(fg1,fg2) /= 5) then
                ! FEP: perturbed
                num_nb15_fep = num_nb15_fep + 1
                nb15_calc_list_fep(num_nb15_fep,ij) = iy
                cycle
              end if

              iyy = start_j + iy
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 exceeds MaxNb15_Fugaku')
#endif
              nb15_calc_list(num_nb15,ixx,1) = iyy
            end if

          end do

          num_nb15_calc(ixx,1) = num_nb15

          ! FEP
          num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
          if (num_nb15_fep /= num_nb15_pre_fep) then
            num_nb15_cell_fep = num_nb15_cell_fep + 1
            nb15_list_fep(num_nb15_cell_fep,ij) = ix
          end if
          num_nb15_pre_fep = num_nb15_fep

        end do

        ! FEP
        nb15_cell_fep(ij) = num_nb15_cell_fep

      end do

    end do

    !$omp barrier

    do i = id+1, ncell+nboundary, nthread

      if (ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif

        ! FEP
        num_nb15_fep = 0
        num_nb15_pre_fep = 0
        num_nb15_cell_fep = 0
 
        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15 = num_nb15_calc(ixx,1)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2) then

              ! FEP
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if (enefunc%fepgrp_nonb(fg1,fg2) == 0) then
                ! FEP: Exclude partA-partB interactions
                cycle
              else if (enefunc%fepgrp_nonb(fg1,fg2) /= 5) then
                ! FEP: perturbed
                num_nb15_fep = num_nb15_fep + 1
                nb15_calc_list_fep(num_nb15_fep,ij) = iy
                cycle
              end if

              iyy = start_j + iy
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ixx,1) = iyy
            end if

          end do

          num_nb15_calc(ixx,1) = num_nb15

          ! FEP
          num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
          if (num_nb15_fep /= num_nb15_pre_fep) then
            num_nb15_cell_fep = num_nb15_cell_fep + 1
            nb15_list_fep(num_nb15_cell_fep,ij) = ix
          end if
          num_nb15_pre_fep = num_nb15_fep

        end do

        ! FEP
        nb15_cell_fep(ij) = num_nb15_cell_fep

      end do

    end do

    !$omp barrier

!   do i = id+1, ncell, nthread
!     start_i = start_atom(i)
!     do ix = 1, natom(i)-1

!       ixx = start_i + ix
!       ! prefetch size:64
!       !
!       num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ixx,1), &
!                                 num_nb15_calc(ixx+1,1))
!       do k = 1, num_add_prefetch
!         nb15_calc_list(num_nb15_calc(ixx,1)+k,ixx,1) = nb15_calc_list(k,ixx+1,1)
!       end do
!     end do
!     ix = natom(i)
!     ixx = start_i + ix
!     num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ixx,1))
!     if (num_nb15_calc(ixx,1) > 0) then
!     do k = 1, num_add_prefetch
!       nb15_calc_list(num_nb15_calc(ixx,1)+k,ixx,1) &
!             = nb15_calc_list(num_nb15_calc(ixx,1),ixx,1)
!     end do
!     end if
!   end do

    !$omp end parallel

    return

  end subroutine update_pairlist_pbc_fugaku_fep

end module sp_pairlist_fugaku_mod
