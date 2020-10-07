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
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: update_pairlist_pbc_fugaku

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_fugaku
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_fugaku(enefunc, domain, pairlist)
  
#ifdef PKTIMER
    use Ctim
#endif
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2_work(1:MaxAtom)
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy
    integer                   :: num_nb15
    integer                   :: id, omp_get_thread_num

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    integer,          pointer,contiguous :: natom(:)
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

#ifdef PKTIMER
    call gettod(sas)
#ifdef FJ_PROF_FAPP
    call fapp_start("PairList",13,0)
#endif
    call timer_sta(13)
#endif

    !$omp parallel default(none)                                       &
    !$omp shared(nthread, ncell, nboundary, natom, coord_pbc, coord,   &
    !$omp        trans1, num_nb15_calc, MaxAtom, exclusion_mask1,      &
    !$omp        nb15_calc_list, maxcell_near, cell_pairlist1,         &
    !$omp        exclusion_mask, pairdist2, maxcell)                   &
    !$omp shared(MaxNb15_Fugaku)                                       &
    !$omp shared(ncell_pairlist1, istcell_pairlist1, ncell_pairlist2,  &
    !$omp       istcell_pairlist2 )                                    &
    !$omp private(num_add_prefetch, k) &
    !$omp private(id, i, j, ij, ix, iy, num_nb15, rtmp, dij, rij2_work)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    num_nb15_calc(1:MaxAtom, 1:ncell+nboundary)  = 0

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif
      do ix = 1, natom(i)-1

        num_nb15 = 0

        do iy = ix+1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            num_nb15 = num_nb15 + 1
#ifdef DEBUG
            if (num_nb15 > MaxNb15_Fugaku) &
              call error_msg(              &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
            nb15_calc_list(num_nb15,ix,i) = (i-1)*MaxAtom+iy

          end if

        end do

        num_nb15_calc(ix,i) = num_nb15

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = id+1, ncell+nboundary, nthread

      if(ncell_pairlist1(i) == 0) cycle

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1


        j = cell_pairlist1(2,ij)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif
       
        do ix = 1, natom(i)

          num_nb15 = num_nb15_calc(ix,i)
          rtmp(1) = coord_pbc(ix,1,i)
          rtmp(2) = coord_pbc(ix,2,i)
          rtmp(3) = coord_pbc(ix,3,i)

          do iy = 1, natom(j)

            dij(1) = rtmp(1) - coord_pbc(iy,1,j)
            dij(2) = rtmp(2) - coord_pbc(iy,2,j)
            dij(3) = rtmp(3) - coord_pbc(iy,3,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2 .and. &
                exclusion_mask(iy,ix,ij) == 1) then
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 exceeds MaxNb15_Fugaku')
#endif
              nb15_calc_list(num_nb15,ix,i) = (j-1)*MaxAtom+iy
            end if

          end do

          num_nb15_calc(ix,i) = num_nb15

        end do
      end do

    end do

    !$omp barrier

    do i = id+1, ncell+nboundary, nthread

      if( ncell_pairlist2(i) == 0 ) cycle

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1


        j = cell_pairlist1(2,ij)

#ifdef DEBUG
        if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
          call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) exceeds MaxAtom')
#endif

        do ix = 1, natom(i)

          num_nb15 = num_nb15_calc(ix,i)
          rtmp(1) = coord_pbc(ix,1,i)
          rtmp(2) = coord_pbc(ix,2,i)
          rtmp(3) = coord_pbc(ix,3,i)

          do iy = 1, natom(j)

            dij(1) = rtmp(1) - coord_pbc(iy,1,j)
            dij(2) = rtmp(2) - coord_pbc(iy,2,j)
            dij(3) = rtmp(3) - coord_pbc(iy,3,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          end do

!ocl loop_nofusion
          do iy = 1, natom(j)

            if (rij2_work(iy) < pairdist2) then
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15_Fugaku) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ix,i) = (j-1)*MaxAtom+iy
            end if

          end do

          num_nb15_calc(ix,i) = num_nb15

        end do
      end do

    end do

    !$omp barrier

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)-1

        ! prefetch size:64
        !
        num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ix,i), &
                                  num_nb15_calc(ix+1,i))
        do k = 1, num_add_prefetch
          nb15_calc_list(num_nb15_calc(ix,i)+k,ix,i) = nb15_calc_list(k,ix+1,i)
        end do
      end do
      ix = natom(i)
      num_add_prefetch = min(64,MaxNb15_Fugaku-num_nb15_calc(ix,i))
      if (num_nb15_calc(ix,i) > 0) then
      do k = 1, num_add_prefetch
        nb15_calc_list(num_nb15_calc(ix,i)+k,ix,i) &
              = nb15_calc_list(num_nb15_calc(ix,i),ix,i)
      end do
      end if
    end do

    !$omp end parallel

#ifdef PKTIMER
    call timer_end(13)
#ifdef FJ_PROF_FAPP
    call fapp_stop("PairList",13,0)
#endif
    call gettod(eae)
    Timc(5)=Timc(5)+(eae-sas)
#endif

    return

  end subroutine update_pairlist_pbc_fugaku

end module sp_pairlist_fugaku_mod
