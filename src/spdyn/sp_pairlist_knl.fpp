!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_knl_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_knl_mod

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
  public  :: update_pairlist_pbc_knl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_knl
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_knl(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3), trans(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15
    integer                   :: num_nb15_total
    integer                   :: id, my_id, my_id_real, omp_get_thread_num
    integer                   :: num_nb15_cell

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0

    do ij = 1, maxcell
      nb15_cell(ij) = 0
    end do

    !$omp parallel                                            &
    !$omp private(id, my_id, my_id_real, i, ix, num_nb15, iy, &
    !$omp         k, ij, j, rtmp, dij, rij2, num_nb15_cell, trans)                              
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell, nthread
      nb15_cell(ij) = 0
    end do

    !$omp barrier

    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      do ix = 1, natom(i)

        num_nb15 = 0
        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          if (exclusion_mask(iy,ix,ij) == 1) then

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 <= pairdist2) then
              num_nb15 = iy
              exit
            end if

          end if
        end do

        nb15_list(2*ix-1,ij) = num_nb15

        if (num_nb15 > 0) then

          do iy = natom(j), 1, -1

            if (exclusion_mask(iy,ix,ij) == 1) then

              dij(1) = rtmp(1) - trans2(iy,1,j)
              dij(2) = rtmp(2) - trans2(iy,2,j)
              dij(3) = rtmp(3) - trans2(iy,3,j)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

              if (rij2 <= pairdist2) then
                num_nb15 = iy
                exit
              end if

            end if
          end do

          nb15_list(2*ix,ij) = num_nb15

        end if

      end do

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      num_nb15_cell = 0

      do ix = 1, natom(i)

        num_nb15 = 0
        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          dij(1) = rtmp(1) - trans2(iy,1,j)
          dij(2) = rtmp(2) - trans2(iy,2,j)
          dij(3) = rtmp(3) - trans2(iy,3,j)
          rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 <= pairdist2) then
            num_nb15 = iy
            num_nb15_cell = 1
            exit
          end if

        end do

        nb15_list(2*ix-1,ij) = num_nb15
        if (num_nb15_cell == 1) nb15_cell(ij) = 1

        if (num_nb15 > 0) then

          do iy = natom(j), 1, -1

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 <= pairdist2) then
              num_nb15 = iy
              num_nb15_cell = 1
              exit
            end if

          end do

        end if

        nb15_list(2*ix,ij) = num_nb15

      end do

    end do

    !$omp end parallel

    return

  end subroutine update_pairlist_pbc_knl

end module sp_pairlist_knl_mod
