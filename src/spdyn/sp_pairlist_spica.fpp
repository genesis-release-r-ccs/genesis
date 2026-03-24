!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_spica_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_spica_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: update_pairlist_spica
  public  :: alloc_pairlist_spica
  public  :: update_pairlist_pspica
  public  :: alloc_pairlist_pspica

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist_spica
  !> @brief        alloc pairlist size in each domain 
  !!               condition 
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist_spica(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy, iatmcls, jatmcls
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: ncharge, ic
    integer                   :: num_nb15, num_nb15w, num_nb15w_p
    integer                   :: num_nb15c

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary, num_charge
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w_p(:,:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w_p(:)
    integer,          pointer,contiguous :: num_nb15_calc_c(:)
    integer,          pointer,contiguous :: charge_list(:)
    integer,          pointer,contiguous :: charge_list_inv(:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)
    integer,          pointer,contiguous :: atmcls(:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)
    integer,          pointer            :: nonb_ljtype(:,:)


    ncell             => domain%num_cell_local
    nboundary         => domain%num_cell_boundary
    natom             => domain%num_atom
    num_charge        => domain%num_charge_domain
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    coord_pbc         => domain%coord_pbc
    cell_pairlist1    => domain%cell_pairlist1
    ncell_pairlist1   => domain%ncell_pairlist1
    istcell_pairlist1 => domain%istcell_pairlist1
    ncell_pairlist2   => domain%ncell_pairlist2
    istcell_pairlist2 => domain%istcell_pairlist2
    atmcls            => domain%atom_cls_no
    charge            => domain%charge
    nonb_ljtype       => enefunc%nonb_ljtype
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    ncharge = 0
    do i = 1, ncell+nboundary
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = start_i + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
        if (abs(charge(ix,i)) > EPS) then
          ncharge = ncharge + 1
        end if
      end do
    end do
    num_charge = ncharge

    k = domain%num_atom_domain
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, ncharge, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxAtom_domain = k * 2
    MaxCharge_domain = ncharge * 2

    call alloc_pairlist(pairlist, PairListSPICA1, MaxAtom_domain, &
                        MaxCharge_domain)
    call alloc_domain(domain, DomainCHARGE, MaxAtom_domain, MaxCharge_domain)
    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_w_p => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    charge_list       => domain%charge_list
    charge_list_inv   => domain%charge_list_inv

    ncharge = 0
    do i = 1, ncell+nboundary
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = start_i + ix
        if (abs(charge(ix,i)) > EPS) then
          ncharge = ncharge + 1
          charge_list(ncharge) = ixx
          charge_list_inv(ixx) = ncharge
        end if
      end do
    end do
    
    k = domain%num_atom_domain
    num_nb15_calc_s(1:k)  = 0
    num_nb15_calc_w(1:k)  = 0
    num_nb15_calc_w_p(1:k)  = 0
    num_nb15_calc_c  (1:num_charge)  = 0

    ! make a pairlist in the same cell
    !
    do i = 1, ncell

      start_i = start_atom(i)

      do ix = 1, natom(i)-1

        ixx = start_i + ix
        num_nb15  = 0
        num_nb15w = 0
        num_nb15w_p = 0
        num_nb15c = 0
        rtmp(1) = coord_pbc(ixx,1,1)
        rtmp(2) = coord_pbc(ixx,2,1)
        rtmp(3) = coord_pbc(ixx,3,1)
        iatmcls = atmcls(ix,i)
        ixx = start_i + ix
        if (abs(charge(ix,i)) > EPS) ic = charge_list_inv(ixx)

        do iy = ix+1, natom(i)

          iyy = start_i + iy
          dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
          dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
          dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
          jatmcls = atmcls(iy,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < pairdist2 .and. exclusion_mask1(iy,ix,i) == 1) then
            if (nonb_ljtype(jatmcls,iatmcls) == 9) then
              num_nb15 = num_nb15 + 1
            else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
              num_nb15w = num_nb15w + 1
            else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
              num_nb15w_p = num_nb15w_p + 1
            end if
            if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,i)) > EPS) then
              num_nb15c = num_nb15c + 1
            end if
          end if
        end do

        num_nb15_calc_s(ixx) = num_nb15
        num_nb15_calc_w(ixx) = num_nb15w
        num_nb15_calc_w_p(ixx) = num_nb15w_p
        if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = 1, ncell+nboundary

      if(ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          if (abs(charge(ix,i)) > EPS) then
            ic = charge_list_inv(ixx)
            num_nb15c = num_nb15_calc_c(ic)
          end if
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)  

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2 .and. exclusion_mask(iy,ix,ij) == 1) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
              end if
              if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,j)) > EPS) then
                num_nb15c = num_nb15c + 1
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p
          if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

        end do
      end do

    end do

    do i = 1, ncell+nboundary

      if(ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          if (abs(charge(ix,i)) > EPS) then
            ic = charge_list_inv(ixx)
            num_nb15c = num_nb15_calc_c(ic)
          end if
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
              end if
              if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,j)) > EPS) then
                num_nb15c = num_nb15c + 1
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p
          if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

        end do
      end do

    end do

    num_nb15    = 0
    num_nb15w   = 0
    num_nb15w_p = 0
    num_nb15c   = 0
    do ixx = 1, domain%num_atom_domain
      num_nb15    = max(num_nb15,   num_nb15_calc_s(ixx))
      num_nb15w   = max(num_nb15w,  num_nb15_calc_w(ixx))
      num_nb15w_p = max(num_nb15w_p,num_nb15_calc_w_p(ixx))
      num_nb15c   = max(num_nb15c  ,num_nb15_calc_c(ixx))
    end do
    do ic = 1, num_charge
      num_nb15c   = max(num_nb15c  ,num_nb15_calc_c(ic))
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_nb15, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, num_nb15w, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, num_nb15w_p, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, num_nb15c  , 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxNb15_S = num_nb15 * 2
    MaxNb15_W = num_nb15w * 2
    MaxNb15_W_P = num_nb15w_p * 2
    MaxNb15_C = num_nb15c * 2

!#ifdef DEBUG
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_S    : ', MaxNb15_S
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_W    : ', MaxNb15_W
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_W_P  : ', MaxNb15_W_P
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_C    : ', MaxNb15_C  
      write(MsgOut,*) ' '
    end if
!#endif

    call alloc_pairlist(pairlist, PairListSPICA2, MaxAtom_domain, &
                        MaxCharge_domain)


    return

  end subroutine alloc_pairlist_spica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_pairlist_pspica
  !> @brief        alloc pairlist size in each domain 
  !!               condition 
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_pairlist_pspica(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy, iatmcls, jatmcls
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: num_nb15, num_nb15w, num_nb15w_p

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w_p(:,:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w_p(:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)
    integer,          pointer,contiguous :: atmcls(:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)
    integer,          pointer            :: nonb_ljtype(:,:)


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
    atmcls            => domain%atom_cls_no
    nonb_ljtype       => enefunc%nonb_ljtype
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    do i = 1, ncell+nboundary
      start_i = start_atom(i)
      do ix = 1, natom(i)
        ixx = start_i + ix
        coord_pbc(ixx,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ixx,2,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ixx,3,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    k = domain%num_atom_domain
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, k, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxAtom_domain = k * 2
    call alloc_pairlist(pairlist, PairListPSPICA1, MaxAtom_domain)
    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_w_p => pairlist%num_nb15_calc_water_p
    
    k = domain%num_atom_domain
    num_nb15_calc_s(1:k)  = 0
    num_nb15_calc_w(1:k)  = 0
    num_nb15_calc_w_p(1:k)  = 0

    ! make a pairlist in the same cell
    !
    do i = 1, ncell

      start_i = start_atom(i)

      do ix = 1, natom(i)-1

        ixx = start_i + ix
        num_nb15  = 0
        num_nb15w = 0
        num_nb15w_p = 0
        rtmp(1) = coord_pbc(ixx,1,1)
        rtmp(2) = coord_pbc(ixx,2,1)
        rtmp(3) = coord_pbc(ixx,3,1)
        iatmcls = atmcls(ix,i)
        ixx = start_i + ix

        do iy = ix+1, natom(i)

          iyy = start_i + iy
          dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
          dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
          dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
          jatmcls = atmcls(iy,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < pairdist2 .and. exclusion_mask1(iy,ix,i) == 1) then
            if (nonb_ljtype(jatmcls,iatmcls) == 9) then
              num_nb15 = num_nb15 + 1
            else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
              num_nb15w = num_nb15w + 1
            else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
              num_nb15w_p = num_nb15w_p + 1
            end if
          end if
        end do

        num_nb15_calc_s(ixx) = num_nb15
        num_nb15_calc_w(ixx) = num_nb15w
        num_nb15_calc_w_p(ixx) = num_nb15w_p

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = 1, ncell+nboundary

      if(ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)  

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2 .and. exclusion_mask(iy,ix,ij) == 1) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p

        end do
      end do

    end do

    do i = 1, ncell+nboundary

      if(ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1

        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p

        end do
      end do

    end do

    num_nb15    = 0
    num_nb15w   = 0
    num_nb15w_p = 0
    do ixx = 1, domain%num_atom_domain
      num_nb15    = max(num_nb15,   num_nb15_calc_s(ixx))
      num_nb15w   = max(num_nb15w,  num_nb15_calc_w(ixx))
      num_nb15w_p = max(num_nb15w_p,num_nb15_calc_w_p(ixx))
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, num_nb15, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, num_nb15w, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, num_nb15w_p, 1, mpi_integer, &
                       mpi_max, mpi_comm_country, ierror)
#endif
    MaxNb15_S = num_nb15 * 2
    MaxNb15_W = num_nb15w * 2
    MaxNb15_W_P = num_nb15w_p * 2
    MaxNb15_C = 1

!#ifdef DEBUG
    if (main_rank) then
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_S    : ', MaxNb15_S
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_W    : ', MaxNb15_W
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_W_P  : ', MaxNb15_W_P
      write(MsgOut,*) ' '
    end if
!#endif

    call alloc_pairlist(pairlist, PairListPSPICA2, MaxAtom_domain, 1)


    return

  end subroutine alloc_pairlist_pspica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_spica
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_spica(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy, iatmcls, jatmcls, ic
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: num_nb15, num_nb15w, num_nb15w_p
    integer                   :: num_nb15c, ncharge
    integer                   :: id 

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w_p(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_c(:,:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w_p(:)
    integer,          pointer,contiguous :: num_nb15_calc_c(:)
    integer,          pointer,contiguous :: charge_list(:), charge_list_inv(:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)
    integer,          pointer,contiguous :: atmcls(:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)
    integer,          pointer            :: nonb_ljtype(:,:)


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
    atmcls            => domain%atom_cls_no
    charge            => domain%charge
    charge_list       => domain%charge_list
    charge_list_inv   => domain%charge_list_inv
    nonb_ljtype       => enefunc%nonb_ljtype

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_w_p => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_w_p=> pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge 
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    if (domain%nonbond_kernel == NBK_Spica) then
      ncharge = 0
      do i = 1, ncell
        start_i = start_atom(i)
        do ix = 1, natom(i)
          ixx = start_i + ix
          if (abs(charge(ix,i)) > EPS) then
            ncharge = ncharge + 1
            charge_list(ncharge) = ixx
            charge_list_inv(ixx) = ncharge
          end if
        end do
      end do
      domain%num_charge_domain_inside = ncharge
      do i = ncell+1, ncell+nboundary
        start_i = start_atom(i)
        do ix = 1, natom(i)
          ixx = start_i + ix
          if (abs(charge(ix,i)) > EPS) then
            ncharge = ncharge + 1
            charge_list(ncharge) = ixx
            charge_list_inv(ixx) = ncharge
          end if
        end do
      end do
      domain%num_charge_domain = ncharge
    end if

    !$omp parallel default(shared)                                     &
    !$omp private(num_add_prefetch, k, start_i, start_j, ixx, iyy)     &
    !$omp private(id, i, j, ij, ix, iy, num_nb15, num_nb15w, rtmp,     &
    !$omp         dij, rij2, iatmcls, jatmcls, num_nb15w_p, num_nb15c, ic)
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
    num_nb15_calc_s(1:k)  = 0
    num_nb15_calc_w(1:k)  = 0
    num_nb15_calc_w_p(1:k)  = 0
    k = domain%num_charge_domain
    num_nb15_calc_c(1:k)  = 0

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      start_i = start_atom(i)

      do ix = 1, natom(i)-1

        ixx = start_i + ix
        num_nb15  = 0
        num_nb15w = 0
        num_nb15w_p = 0
        num_nb15c = 0
        rtmp(1) = coord_pbc(ixx,1,1)
        rtmp(2) = coord_pbc(ixx,2,1)
        rtmp(3) = coord_pbc(ixx,3,1)
        iatmcls = atmcls(ix,i)
        ixx = start_i + ix
        if (abs(charge(ix,i)) > EPS) ic = charge_list_inv(ixx)

        do iy = ix+1, natom(i)

          iyy = start_i + iy
          dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
          dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
          dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
          jatmcls = atmcls(iy,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < pairdist2 .and. exclusion_mask1(iy,ix,i) == 1) then
            if (nonb_ljtype(jatmcls,iatmcls) == 9) then
              num_nb15 = num_nb15 + 1
              nb15_calc_list_s(num_nb15,ixx) = iyy
            else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
              num_nb15w = num_nb15w + 1
              nb15_calc_list_w(num_nb15w,ixx) = iyy
            else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
              num_nb15w_p = num_nb15w_p + 1
              nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
            end if
            if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,i)) > EPS) then
              num_nb15c = num_nb15c + 1
              nb15_calc_list_c(num_nb15c,ic) = iyy
            end if
          end if
        end do

        num_nb15_calc_s(ixx) = num_nb15
        num_nb15_calc_w(ixx) = num_nb15w
        num_nb15_calc_w_p(ixx) = num_nb15w_p
        if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = id+1, ncell+nboundary, nthread

      if(ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          if (abs(charge(ix,i)) > EPS) then
            ic = charge_list_inv(ixx)
            num_nb15c = num_nb15_calc_c(ic)
          end if
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)  

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2 .and. exclusion_mask(iy,ix,ij) == 1) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then

                num_nb15 = num_nb15 + 1
                nb15_calc_list_s(num_nb15,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
                nb15_calc_list_w(num_nb15w,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
                nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
              end if
              if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,j)) > EPS) then
                num_nb15c = num_nb15c + 1
                nb15_calc_list_c(num_nb15c,ic) = iyy
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p
          if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

        end do
      end do

    end do

    do i = id+1, ncell+nboundary, nthread

      if(ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          if (abs(charge(ix,i)) > EPS) then
            ic = charge_list_inv(ixx)
            num_nb15c = num_nb15_calc_c(ic)
          end if
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
                nb15_calc_list_s(num_nb15,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
                nb15_calc_list_w(num_nb15w,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
                nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
              end if
              if (abs(charge(ix,i)) > EPS .and. abs(charge(iy,j)) > EPS) then
                num_nb15c = num_nb15c + 1
                nb15_calc_list_c(num_nb15c,ic) = iyy
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p
          if (abs(charge(ix,i)) > EPS) num_nb15_calc_c(ic) = num_nb15c

        end do
      end do

    end do

    !$omp barrier

    !$omp end parallel

    return

  end subroutine update_pairlist_spica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pspica
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pspica(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, ix, iy, iatmcls, jatmcls
    integer                   :: start_i, start_j, ixx, iyy, num_atom
    integer                   :: num_nb15, num_nb15w, num_nb15w_p
    integer                   :: id

    integer :: num_add_prefetch, k

    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:), coord_pbc(:,:,:)
    integer,          pointer,contiguous :: natom(:), start_atom(:)
    integer(int2),    pointer,contiguous :: cell_pairlist1(:,:)
    integer,          pointer            :: ncell, nboundary
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w_p(:,:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w_p(:)
    integer(1),       pointer,contiguous :: exclusion_mask(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:)
    integer,          pointer,contiguous :: atmcls(:,:)

    real(dp)                             :: sas,eae

    integer,          pointer,contiguous :: ncell_pairlist1(:)
    integer,          pointer,contiguous :: istcell_pairlist1(:)
    integer,          pointer,contiguous :: ncell_pairlist2(:)
    integer,          pointer,contiguous :: istcell_pairlist2(:)
    integer,          pointer            :: nonb_ljtype(:,:)


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
    atmcls            => domain%atom_cls_no
    nonb_ljtype       => enefunc%nonb_ljtype

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_w_p => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_w_p  => pairlist%nb15_calc_list_water_p
    exclusion_mask    => enefunc%exclusion_mask
    exclusion_mask1   => enefunc%exclusion_mask1

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    !$omp parallel default(shared)                                     &
    !$omp private(num_add_prefetch, k, start_i, start_j, ixx, iyy)     &
    !$omp private(id, i, j, ij, ix, iy, num_nb15, num_nb15w, rtmp,     &
    !$omp         dij, rij2, iatmcls, jatmcls, num_nb15w_p)
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
    num_nb15_calc_s(1:k)  = 0
    num_nb15_calc_w(1:k)  = 0
    num_nb15_calc_w_p(1:k)  = 0

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      start_i = start_atom(i)

      do ix = 1, natom(i)-1

        ixx = start_i + ix
        num_nb15  = 0
        num_nb15w = 0
        num_nb15w_p = 0
        rtmp(1) = coord_pbc(ixx,1,1)
        rtmp(2) = coord_pbc(ixx,2,1)
        rtmp(3) = coord_pbc(ixx,3,1)
        iatmcls = atmcls(ix,i)
        ixx = start_i + ix

        do iy = ix+1, natom(i)

          iyy = start_i + iy
          dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
          dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
          dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
          jatmcls = atmcls(iy,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < pairdist2 .and. exclusion_mask1(iy,ix,i) == 1) then
            if (nonb_ljtype(jatmcls,iatmcls) == 9) then
              num_nb15 = num_nb15 + 1
              nb15_calc_list_s(num_nb15,ixx) = iyy
            else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
              num_nb15w = num_nb15w + 1
              nb15_calc_list_w(num_nb15w,ixx) = iyy
            else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
              num_nb15w_p = num_nb15w_p + 1
              nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
            end if
          end if
        end do

        num_nb15_calc_s(ixx) = num_nb15
        num_nb15_calc_w(ixx) = num_nb15w
        num_nb15_calc_w_p(ixx) = num_nb15w_p

      end do
    end do

    !$omp barrier

    ! Make a pairlist between different cells
    !
    do i = id+1, ncell+nboundary, nthread

      if(ncell_pairlist1(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist1(i), istcell_pairlist1(i)+ncell_pairlist1(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)  

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2 .and. exclusion_mask(iy,ix,ij) == 1) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then

                num_nb15 = num_nb15 + 1
                nb15_calc_list_s(num_nb15,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
                nb15_calc_list_w(num_nb15w,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
                nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p

        end do
      end do

    end do

    do i = id+1, ncell+nboundary, nthread

      if(ncell_pairlist2(i) == 0) cycle
      start_i = start_atom(i)

      do ij = istcell_pairlist2(i), istcell_pairlist2(i)+ncell_pairlist2(i)-1


        j = cell_pairlist1(2,ij)
        start_j = start_atom(j)

        do ix = 1, natom(i)

          ixx = start_i + ix
          num_nb15  = num_nb15_calc_s(ixx)
          num_nb15w = num_nb15_calc_w(ixx)
          num_nb15w_p = num_nb15_calc_w_p(ixx)
          rtmp(1) = coord_pbc(ixx,1,1)
          rtmp(2) = coord_pbc(ixx,2,1)
          rtmp(3) = coord_pbc(ixx,3,1)
          iatmcls = atmcls(ix,i)

          do iy = 1, natom(j)

            iyy = start_j + iy
            dij(1) = rtmp(1) - coord_pbc(iyy,1,1)
            dij(2) = rtmp(2) - coord_pbc(iyy,2,1)
            dij(3) = rtmp(3) - coord_pbc(iyy,3,1)
            jatmcls = atmcls(iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then
              if (nonb_ljtype(jatmcls,iatmcls) == 9) then
                num_nb15 = num_nb15 + 1
                nb15_calc_list_s(num_nb15,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 12) then
                num_nb15w = num_nb15w + 1
                nb15_calc_list_w(num_nb15w,ixx) = iyy
              else if (nonb_ljtype(jatmcls,iatmcls) == 5) then
                num_nb15w_p = num_nb15w_p + 1
                nb15_calc_list_w_p(num_nb15w_p,ixx) = iyy
              end if
            end if

          end do

          num_nb15_calc_s(ixx) = num_nb15
          num_nb15_calc_w(ixx) = num_nb15w
          num_nb15_calc_w_p(ixx) = num_nb15w_p

        end do
      end do

    end do

    !$omp barrier

    !$omp end parallel

    return

  end subroutine update_pairlist_pspica

end module sp_pairlist_spica_mod
