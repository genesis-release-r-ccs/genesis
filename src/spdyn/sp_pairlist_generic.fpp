!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_generic_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_generic_mod

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
  public  :: update_pairlist_pbc_generic
  public  :: update_pairlist_pbc_check_generic
  ! FEP
  public  :: update_pairlist_pbc_generic_fep
  public  :: update_pairlist_pbc_check_generic_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_generic
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_generic(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3), trans(3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table


    real(dp)                             :: sas,eae


    table           => enefunc%table%table
    table_order     => enefunc%table%table_order
    nonb_limiter    => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell, trans)                       &
    !$omp reduction(+:num_nb15_total)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      num_nb15 = 0
      num_nb15_pre = 0

#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      do ix = 1, natom(i) - 1

        rtmp(1) = trans2(ix,1,i) 
        rtmp(2) = trans2(ix,2,i) 
        rtmp(3) = trans2(ix,3,i) 

        do iy = ix + 1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            dij(1) = rtmp(1) - trans2(iy,1,i) 
            dij(2) = rtmp(2) - trans2(iy,2,i) 
            dij(3) = rtmp(3) - trans2(iy,3,i) 
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                dij(1) = trans2(ix,1,i) - trans2(iy,1,i)
                dij(2) = trans2(ix,2,i) - trans2(iy,2,i)
                dij(3) = trans2(ix,3,i) - trans2(iy,3,i)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < err_minimum_contact) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
              end if
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list1(num_nb15,i) = iy
            end if
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)
#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          if (exclusion_mask(iy,ix,ij)==1) then

            dij(1) = rtmp(1) - trans2(iy,1,j) 
            dij(2) = rtmp(2) - trans2(iy,2,j) 
            dij(3) = rtmp(3) - trans2(iy,3,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                if (rij2 < err_minimum_contact) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
              end if
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    !$omp end parallel

!#ifdef HAVE_MPI_GENESIS
!    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
!         0, mpi_comm_country, ierror)
!#endif

    return

  end subroutine update_pairlist_pbc_generic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check_generic
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ, CK, NT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check_generic(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3), trans(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: small_contact
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table


    table           => enefunc%table%table
    table_order     => enefunc%table%table_order
    nonb_limiter    => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0
    small_contact   =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell, trans)                       &
    !$omp reduction(+:num_nb15_total) reduction(+:small_contact)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      num_nb15 = 0
      num_nb15_pre = 0

      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      do ix = 1, natom(i) - 1

        do iy = ix + 1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            num_nb15 = num_nb15 + 1

            if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
              dij(1) = trans2(ix,1,i) - trans2(iy,1,i)
              dij(2) = trans2(ix,2,i) - trans2(iy,2,i)
              dij(3) = trans2(ix,3,i) - trans2(iy,3,i)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < err_minimum_contact) then
                if (.not. nonb_limiter) &
                  call error_msg( &
                  'Debug: Update_Pairlist_Pbc> too small contact')
                small_contact = small_contact + 1
              end if
            endif
            if (num_nb15 > MaxNb15) &
              call error_msg( &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')

            nb15_calc_list1(num_nb15,i) = iy
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          if (exclusion_mask(iy,ix,ij)==1) then

            dij(1) = rtmp(1) - trans2(iy,1,j) 
            dij(2) = rtmp(2) - trans2(iy,2,j) 
            dij(3) = rtmp(3) - trans2(iy,3,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1

              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
                  small_contact = small_contact + 1
                end if
              endif
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              num_nb15 = num_nb15 + 1

              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    !$omp end parallel

    if (small_contact > 0) then
      write(MsgOut, *) "Warning: small contacts exist in inner loop"
    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc_check_generic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_generic_fep
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition for FEP
  !! @authors      HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_generic_fep(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3), trans(3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table


    real(dp)                             :: sas,eae

    ! FEP
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_cell_fep(:), nb15_list_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer                   :: num_nb15_pre_fep, num_nb15_fep
    integer                   :: num_nb15_total_fep, num_nb15_total1_fep
    integer                   :: num_nb15_cell_fep

    table           => enefunc%table%table
    table_order     => enefunc%table%table_order
    nonb_limiter    => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0

    ! FEP
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell, trans, fg1, fg2,             &
    !$omp         num_nb15_fep, num_nb15_pre_fep, num_nb15_cell_fep)           &
    !$omp reduction(+:num_nb15_total), reduction(+:num_nb15_total_fep)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      num_nb15 = 0
      num_nb15_pre = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0

#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      do ix = 1, natom(i) - 1

        rtmp(1) = trans2(ix,1,i) 
        rtmp(2) = trans2(ix,2,i) 
        rtmp(3) = trans2(ix,3,i) 

        do iy = ix + 1, natom(i)

          if (exclusion_mask1(iy,ix,i) == 1) then

            dij(1) = rtmp(1) - trans2(iy,1,i) 
            dij(2) = rtmp(2) - trans2(iy,2,i) 
            dij(3) = rtmp(3) - trans2(iy,3,i) 
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < pairdist2) then

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
#ifdef DEBUG
              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                dij(1) = trans2(ix,1,i) - trans2(iy,1,i)
                dij(2) = trans2(ix,2,i) - trans2(iy,2,i)
                dij(3) = trans2(ix,3,i) - trans2(iy,3,i)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < err_minimum_contact) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
              end if
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list1(num_nb15,i) = iy
            end if
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc1_fep(ix,i) = num_nb15_fep - num_nb15_pre_fep
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)
#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          if (exclusion_mask(iy,ix,ij)==1) then

            dij(1) = rtmp(1) - trans2(iy,1,j) 
            dij(2) = rtmp(2) - trans2(iy,2,j) 
            dij(3) = rtmp(3) - trans2(iy,3,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

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

              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                if (rij2 < err_minimum_contact) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
              end if
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

      nb15_cell(ij) = num_nb15_cell

      ! FEP
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

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
 
              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

      nb15_cell(ij) = num_nb15_cell

      ! FEP
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    !$omp end parallel

!#ifdef HAVE_MPI_GENESIS
!    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
!         0, mpi_comm_country, ierror)
!#endif

    return

  end subroutine update_pairlist_pbc_generic_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check_generic_fep
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition for FEP
  !! @authors      HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check_generic_fep(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3), trans(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: small_contact
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: natom(:)
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table

    ! FEP
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_cell_fep(:), nb15_list_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer                   :: num_nb15_pre_fep, num_nb15_fep
    integer                   :: num_nb15_total_fep, num_nb15_total1_fep
    integer                   :: num_nb15_cell_fep

    table           => enefunc%table%table
    table_order     => enefunc%table%table_order
    nonb_limiter    => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0
    small_contact   =  0

    ! FEP
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell, trans, fg1, fg2,             &
    !$omp         num_nb15_fep, num_nb15_pre_fep, num_nb15_cell_fep)           &
    !$omp reduction(+:num_nb15_total), reduction(+:num_nb15_total_fep)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      num_nb15 = 0
      num_nb15_pre = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0

      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      do ix = 1, natom(i) - 1

        do iy = ix + 1, natom(i)

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
              if (table .and. table_order == 1 .and. (.not. nonb_limiter))then
                dij(1) = trans2(ix,1,i) - trans2(iy,1,i)
                dij(2) = trans2(ix,2,i) - trans2(iy,2,i)
                dij(3) = trans2(ix,3,i) - trans2(iy,3,i)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                    call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
                  small_contact = small_contact + 1
                end if
              endif
              if (num_nb15_fep > MaxNb15) &
                call error_msg( &
                'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
              nb15_calc_list1_fep(num_nb15_fep,i) = iy
              cycle
            end if
 
            num_nb15 = num_nb15 + 1

            if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
              dij(1) = trans2(ix,1,i) - trans2(iy,1,i)
              dij(2) = trans2(ix,2,i) - trans2(iy,2,i)
              dij(3) = trans2(ix,3,i) - trans2(iy,3,i)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < err_minimum_contact) then
                if (.not. nonb_limiter) &
                  call error_msg( &
                  'Debug: Update_Pairlist_Pbc> too small contact')
                small_contact = small_contact + 1
              end if
            endif
            if (num_nb15 > MaxNb15) &
              call error_msg( &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')

            nb15_calc_list1(num_nb15,i) = iy
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc1_fep(ix,i) = num_nb15_fep - num_nb15_pre_fep
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          if (exclusion_mask(iy,ix,ij)==1) then

            dij(1) = rtmp(1) - trans2(iy,1,j) 
            dij(2) = rtmp(2) - trans2(iy,2,j) 
            dij(3) = rtmp(3) - trans2(iy,3,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              ! FEP
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if (enefunc%fepgrp_nonb(fg1,fg2) == 0) then
                ! FEP: Exclude partA-partB interactions
                cycle
              else if (enefunc%fepgrp_nonb(fg1,fg2) /= 5) then
                ! FEP: perturbed
                num_nb15_fep = num_nb15_fep + 1
                if (table .and. table_order == 1 .and. (.not. nonb_limiter))then
                  if (rij2 < err_minimum_contact) then
                    if (.not. nonb_limiter) &
                      call error_msg( &
                      'Debug: Update_Pairlist_Pbc> too small contact')
                    small_contact = small_contact + 1
                  end if
                endif
                if (num_nb15_fep > MaxNb15) &
                  call error_msg( &
                  'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
                nb15_calc_list_fep(num_nb15_fep,ij) = iy
                cycle
              end if

              num_nb15 = num_nb15 + 1
              if (table .and. table_order == 1 .and. (.not. nonb_limiter)) then
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                  call error_msg( &
                    'Debug: Update_Pairlist_Pbc> too small contact')
                  small_contact = small_contact + 1
                end if
              endif
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

      nb15_cell(ij) = num_nb15_cell

      ! FEP
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      trans(1:3) = real(cell_move(1:3,j,i),wp)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      ! FEP
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        rtmp(1) = trans2(ix,1,i) + trans(1)*system_size(1)
        rtmp(2) = trans2(ix,2,i) + trans(2)*system_size(2)
        rtmp(3) = trans2(ix,3,i) + trans(3)*system_size(3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1) = rtmp(1) - trans2(iy,1,j)
            dij(2) = rtmp(2) - trans2(iy,2,j)
            dij(3) = rtmp(3) - trans2(iy,3,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              ! FEP
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if (enefunc%fepgrp_nonb(fg1,fg2) == 0) then
                ! FEP: Exclude partA-partB interactions
                cycle
              else if (enefunc%fepgrp_nonb(fg1,fg2) /= 5) then
                ! FEP: perturbed
                num_nb15_fep = num_nb15_fep + 1
                if (num_nb15_fep > MaxNb15) &
                  call error_msg( &
                  'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
                nb15_calc_list_fep(num_nb15_fep,ij) = iy
                cycle
              end if
 
              num_nb15 = num_nb15 + 1

              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

      nb15_cell(ij) = num_nb15_cell

      ! FEP
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    !$omp end parallel

    if (small_contact > 0) then
      write(MsgOut, *) "Warning: small contacts exist in inner loop"
    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc_check_generic_fep

end module sp_pairlist_generic_mod
