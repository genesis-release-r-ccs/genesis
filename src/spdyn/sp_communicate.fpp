!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_communicate_mod
!> @brief   utilities for mpi communication
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_communicate_mod

  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_comm
    integer                       :: num_cell(3)
    integer                       :: send_force_lower_size(3)
    integer                       :: send_force_upper_size(3)
    integer                       :: recv_force_lower_size(3)
    integer                       :: recv_force_upper_size(3)
    integer                       :: send_coord_lower_size(3)
    integer                       :: send_coord_upper_size(3)
    integer                       :: recv_coord_lower_size(3)
    integer                       :: recv_coord_upper_size(3)
    integer                       :: send_enm_lower_size(3)
    integer                       :: send_enm_upper_size(3)
    integer                       :: recv_enm_lower_size(3)
    integer                       :: recv_enm_upper_size(3)
    integer,          allocatable :: ic_lower_send(:,:)
    integer,          allocatable :: ic_upper_send(:,:)
    integer,          allocatable :: ic_lower_recv(:,:)
    integer,          allocatable :: ic_upper_recv(:,:)
    integer,          allocatable :: if_lower_send(:,:)
    integer,          allocatable :: if_upper_send(:,:)
    integer,          allocatable :: if_lower_recv(:,:)
    integer,          allocatable :: if_upper_recv(:,:)
    integer,          allocatable :: int_send(:,:)
    integer,          allocatable :: int_recv(:,:)
    real(wip),        allocatable :: buf_send(:,:)
    real(wip),        allocatable :: buf_recv(:,:)
  end type s_comm

  ! subroutines
  public  :: setup_communicate
  public  :: setup_communicate_size
  public  :: communicate_force
  public  :: communicate_coor
  public  :: communicate_ptl
  public  :: communicate_constraints
  public  :: communicate_bond
  public  :: update_cell_boundary
  public  :: update_cell_boundary_enm
  public  :: update_cell_size_constraints
  public  :: update_cell_size_enm
  public  :: update_random_number
  ! FEP
  public  :: communicate_constraints_fep
  public  :: communicate_bond_fep
  public  :: update_cell_boundary_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate
  !> @brief        Initial setup of communicate among processors
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !! @param[out]   comm     : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate(boundary, domain, comm)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variable
    integer                  :: cell(3), max_cell
    integer                  :: i, j, k
    integer                  :: index_upper, index_lower, local_index, bd_index
    integer                  :: ix, iy, iz
    integer                  :: alloc_stat, size
    integer,         pointer :: ncell_local
    integer,         pointer :: cell_start(:), cell_end(:), cell_length(:)
    integer(int2),   pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: num_domain(:)


    num_domain  => boundary%num_domain

    ncell_local => domain%num_cell_local
    cell_start  => domain%cell_start
    cell_end    => domain%cell_end
    cell_length => domain%cell_length
    cell_g2b    => domain%cell_g2b
    cell_g2l    => domain%cell_g2l

    cell(1)     =  boundary%num_cells_x
    cell(2)     =  boundary%num_cells_y
    cell(3)     =  boundary%num_cells_z

    ! Memory allocation
    !
    alloc_stat = 0
    max_cell = max(cell_length(2)*cell_length(3),    &
                  (cell_length(1)+2)*cell_length(3), &
                  (cell_length(1)+2)*(cell_length(2)+2))

    size = 8 * MaxAtom * max_cell
    allocate(comm%ic_upper_send(max_cell,3), comm%ic_lower_send(max_cell,3), &
             comm%ic_upper_recv(max_cell,3), comm%ic_lower_recv(max_cell,3), &
             comm%if_upper_send(max_cell,3), comm%if_lower_send(max_cell,3), &
             comm%if_upper_recv(max_cell,3), comm%if_lower_recv(max_cell,3), &
             comm%buf_send(size,2), comm%buf_recv(size,2),                   &
             comm%int_send(size,2), comm%int_recv(size,2),                   &
             stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! Communicate in x direction
    !
    k = 0
    if (num_domain(1) > 1) then
      do i = cell_start(2), cell_end(2)
        do j = cell_start(3), cell_end(3)

          k = k + 1

          ! sending coordinate, receving force
          !
          index_upper = cell_end(1) + (i-1)*cell(1) + (j-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_upper)

          if (local_index /= 0) then
            comm%ic_upper_send(k,1) = local_index
            comm%if_upper_recv(k,1) = local_index
          end if

          index_lower = cell_start(1) + (i-1)*cell(1) + (j-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_lower)

          if (local_index /= 0) then
            comm%ic_lower_send(k,1) = local_index
            comm%if_lower_recv(k,1) = local_index
          end if

          ! receving coordinate, sending force
          !
          if (cell_end(1) == cell(1)) then
            ix = 1
          else
            ix = cell_end(1) + 1
          end if

          index_upper = ix + (i-1)*cell(1) + (j-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_upper)

          if (bd_index /= 0) then
            comm%ic_upper_recv(k,1) = bd_index + ncell_local
            comm%if_upper_send(k,1) = bd_index + ncell_local
          end if

          if (cell_start(1) == 1) then
            ix = cell(1)
          else
            ix = cell_start(1) - 1
          end if

          index_lower = ix + (i-1)*cell(1) + (j-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_lower)

          if (bd_index /= 0) then
            comm%ic_lower_recv(k,1) = bd_index + ncell_local
            comm%if_lower_send(k,1) = bd_index + ncell_local
          end if

        end do
      end do
    end if

    comm%num_cell(1) = k

    ! communicate in y direction
    !
    k = 0
    if (num_domain(2) > 1) then
      do i = cell_start(1)-1, cell_end(1)+1
        if (i == 0) then
          ix = cell(1)
        else if (i == (cell(1)+1)) then
          ix = 1
        else
          ix = i
        end if

        do j = cell_start(3), cell_end(3)

          k = k + 1

          ! sending coordinate, receving force
          !
          index_upper = ix + (cell_end(2)-1)*cell(1) + (j-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_upper)
          bd_index = cell_g2b(index_upper)

          if (local_index /= 0) then
            comm%ic_upper_send(k,2) = local_index
            comm%if_upper_recv(k,2) = local_index
          else if (bd_index /= 0) then
            comm%ic_upper_send(k,2) = bd_index + ncell_local
            comm%if_upper_recv(k,2) = bd_index + ncell_local
          end if

          index_lower = ix + (cell_start(2)-1)*cell(1) + (j-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_lower)
          bd_index = cell_g2b(index_lower)

          if (local_index /= 0) then
            comm%ic_lower_send(k,2) = local_index
            comm%if_lower_recv(k,2) = local_index
          else if (bd_index /= 0) then
            comm%ic_lower_send(k,2) = bd_index + ncell_local
            comm%if_lower_recv(k,2) = bd_index + ncell_local
          end if

          ! receving coordinate, sending force
          !
          if (cell_end(2) == cell(2)) then
            iy = 1
          else
            iy = cell_end(2) + 1
          end if

          index_upper = ix + (iy-1)*cell(1) + (j-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_upper)

          if (bd_index /= 0) then
            comm%ic_upper_recv(k,2) = bd_index + ncell_local
            comm%if_upper_send(k,2) = bd_index + ncell_local
          end if

          if (cell_start(2) == 1) then
            iy = cell(2)
          else
            iy = cell_start(2) - 1
          end if

          index_lower = ix + (iy-1)*cell(1) + (j-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_lower)

          if (bd_index /= 0) then
            comm%ic_lower_recv(k,2) = bd_index + ncell_local
            comm%if_lower_send(k,2) = bd_index + ncell_local
          end if

        end do
      end do
    end if

    comm%num_cell(2) = k

    ! sending cell index (coordinate, z direction)
    !
    k = 0

    if (num_domain(3) > 1) then

      do i = cell_start(1)-1, cell_end(1)+1
        if (i == 0) then
          ix = cell(1)
        else if (i == (cell(1)+1)) then
          ix = 1
        else
          ix = i
        end if

        do j = cell_start(2)-1, cell_end(2)+1

          if (j == 0) then
            iy = cell(2)
          else if (j == (cell(2)+1)) then
            iy = 1
          else
            iy = j
          end if

          k = k + 1

          ! sending coordinates, receving forces
          !
          index_upper = ix + (iy-1)*cell(1) + (cell_end(3)-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_upper)
          bd_index = cell_g2b(index_upper)

          if (local_index /= 0) then
            comm%ic_upper_send(k,3) = local_index
            comm%if_upper_recv(k,3) = local_index
          else if (bd_index /= 0) then
            comm%ic_upper_send(k,3) = bd_index + ncell_local
            comm%if_upper_recv(k,3) = bd_index + ncell_local
          end if

          index_lower = ix + (iy-1)*cell(1) + (cell_start(3)-1)*cell(1)*cell(2)
          local_index = cell_g2l(index_lower)
          bd_index = cell_g2b(index_lower)

          if (local_index /= 0) then
            comm%ic_lower_send(k,3) = local_index
            comm%if_lower_recv(k,3) = local_index
          else if (bd_index /= 0) then
            comm%ic_lower_send(k,3) = bd_index + ncell_local
            comm%if_lower_recv(k,3) = bd_index + ncell_local
          end if

          ! receving coordinates, sending forces
          !
          if (cell_end(3) == cell(3)) then
            iz = 1
          else
            iz = cell_end(3) + 1
          end if

          index_upper = ix + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_upper)

          if (bd_index /= 0) then
            comm%ic_upper_recv(k,3) = bd_index + ncell_local
            comm%if_upper_send(k,3) = bd_index + ncell_local
          end if

          if (cell_start(3) == 1) then
            iz = cell(3)
          else
            iz = cell_start(3) - 1
          end if

          index_lower = ix + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          bd_index = cell_g2b(index_lower)

          if (bd_index /= 0) then
            comm%ic_lower_recv(k,3) = bd_index + ncell_local
            comm%if_lower_send(k,3) = bd_index + ncell_local
          end if

        end do
      end do
    end if

    comm%num_cell(3) = k

    return

  end subroutine setup_communicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate_size
  !> @brief        Check the size of the transfer data
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[out]   comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate_size(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: send_size(4), recv_size(4)
    integer                  :: i, k, ii, ic
    integer                  :: irequest0, irequest1, irequest2, irequest3
    integer                  :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    integer,         pointer :: num_cell(:), natom(:)
    integer,         pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,         pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)


    natom         => domain%num_atom
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    num_cell      => comm%num_cell

    if_lower_send => comm%if_lower_send
    if_upper_send => comm%if_upper_send
    if_lower_recv => comm%if_lower_recv
    if_upper_recv => comm%if_upper_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv

    do ii = 1, 3

      ! check the send size of forces (lower)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = if_lower_send(i,ii)
        k = k + natom(ic)
      end do
      send_size(1) = k

      ! check the send size of forces (upper)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = if_upper_send(i,ii)
        k = k + natom(ic)
      end do
      send_size(2) = k

      ! check the send size of coordinates (lower)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_lower_send(i,ii)
        k = k + natom(ic)
      end do
      send_size(3) = k

      ! check the send size of coordinates (upper)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_upper_send(i,ii)
        k = k + natom(ic)
      end do
      send_size(4) = k

#ifdef HAVE_MPI_GENESIS
      ! send the size of the data (lower)
      !
      call mpi_irecv(recv_size(1), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii),     &
                     mpi_comm_city, irequest0, ierror)
      call mpi_isend(send_size(1), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank,     &
                     mpi_comm_city, irequest1, ierror)
      ! send the size of the data (upper)
      !
      call mpi_irecv(recv_size(2), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank+nproc_city,     &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(2), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii)+nproc_city,     &
                     mpi_comm_city, irequest3, ierror)
      ! send the size of the data (lower)
      !
      call mpi_irecv(recv_size(3), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii)+2*nproc_city,     &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(send_size(3), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank+2*nproc_city,     &
                     mpi_comm_city, irequest5, ierror)
      ! send the size of the data (upper)
      !
      call mpi_irecv(recv_size(4), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank+3*nproc_city,     &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(send_size(4), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii)+3*nproc_city,     &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest0, istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)
      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      comm%send_force_lower_size(ii) = send_size(1)
      comm%recv_force_upper_size(ii) = recv_size(1)
      comm%send_force_upper_size(ii) = send_size(2)
      comm%recv_force_lower_size(ii) = recv_size(2)
      comm%send_coord_lower_size(ii) = send_size(3)
      comm%recv_coord_upper_size(ii) = recv_size(3)
      comm%send_coord_upper_size(ii) = send_size(4)
      comm%recv_coord_lower_size(ii) = recv_size(4)

      k = 0
      do i = 1, num_cell(ii)
        ic = if_upper_recv(i,ii)
        k = k + natom(ic)
      end do
       if (recv_size(1) /= k) &
         call error_msg('Setup_Communicate_Size> Disagreement between'//&
                        ' sending and receving processors')

       k = 0
       do i = 1, num_cell(ii)
         ic = if_lower_recv(i,ii)
         k = k + natom(ic)
       end do
       if (recv_size(2) /= k) &
         call error_msg('Setup_Communicate_Size> Disagreement between'//&
                        ' sending and receving processors')

       k = 0
       do i = 1, num_cell(ii)
         ic = ic_upper_recv(i,ii)
         k = k + natom(ic)
       end do
       if (recv_size(3) /= k) &
         call error_msg('Setup_Communicate_Size> Disagreement between'//&
                         'sending and receving processors')

       k = 0
       do i = 1, num_cell(ii)
         ic = ic_lower_recv(i,ii)
         k = k + natom(ic)
       end do
       if (recv_size(4) /= k) &
         call error_msg('Setup_Communicate_Size> Disagreement between'//&
                        'sending and receving processors')

    end do

    return

  end subroutine setup_communicate_size

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_force
  !> @brief        Pack the force date and trasfer it
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[inout] comm   : communication information
  !! @param[inout] force  : forces of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_force(domain, comm, force)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm
    real(wip),               intent(inout) :: force(:,:,:)

    ! local variable
    integer                  :: i, k, ii, ic, ix, id, omp_get_thread_num
#ifdef HAVE_MPI_GENESIS
    integer                  :: num_req
    integer                  :: ireq(4)
    integer                  :: istat(mpi_status_size,4)
#endif
!    integer                  :: irequest, irequest1, irequest2, irequest3
!#ifdef HAVE_MPI_GENESIS
!    integer                  :: istatus(mpi_status_size)
!#endif

    real(wip),       pointer :: buf_send(:,:), buf_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:), natom(:)
    integer,         pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,         pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,         pointer :: upper_send(:), upper_recv(:)
    integer,         pointer :: lower_send(:), lower_recv(:)


    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    natom         => domain%num_atom

    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    if_lower_send => comm%if_lower_send
    if_upper_send => comm%if_upper_send
    if_lower_recv => comm%if_lower_recv
    if_upper_recv => comm%if_upper_recv
    upper_send    => comm%send_force_upper_size
    upper_recv    => comm%recv_force_upper_size
    lower_send    => comm%send_force_lower_size
    lower_recv    => comm%recv_force_lower_size

    !$omp parallel private(id, ii, k, i, ic, ix)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do ii = 3, 1, -1

      ! Pack the force data (lower)
      !
      if (mod(id,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = if_lower_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(3*k-2:3*k,1) = force(1:3,ix,ic)
          end do
        end do
      end if

      ! pack the force data (upper)
      !
      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = if_upper_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(3*k-2:3*k,2) = force(1:3,ix,ic)
          end do
        end do
      end if

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      ! send the data (lower)
      !
      num_req = 1
      call mpi_irecv(buf_recv(1,1), 3*upper_recv(ii), mpi_wip_real, &
                     iproc_upper(ii),                               &
                     iproc_upper(ii),  &
                     mpi_comm_city, ireq(num_req), ierror)

      num_req = num_req + 1
      call mpi_isend(buf_send(1,1), 3*lower_send(ii), mpi_wip_real, &
                     iproc_lower(ii),                               &
                     my_city_rank,  &
                     mpi_comm_city, ireq(num_req), ierror)

      ! send the data (upper)
      !
      num_req = num_req + 1
      call mpi_irecv(buf_recv(1,2), 3*lower_recv(ii), mpi_wip_real, &
                     iproc_lower(ii),                               &
                     my_city_rank+nproc_city,  &
                     mpi_comm_city, ireq(num_req), ierror)

      num_req = num_req + 1
      call mpi_isend(buf_send(1,2), 3*upper_send(ii), mpi_wip_real, &
                     iproc_upper(ii),                               &
                     iproc_upper(ii)+nproc_city,  &
                     mpi_comm_city, ireq(num_req), ierror)

      call mpi_waitall(num_req, ireq, istat, ierror)
!     call mpi_wait(irequest,  istatus, ierror)
!     call mpi_wait(irequest1, istatus, ierror)
!     call mpi_wait(irequest2, istatus, ierror)
!     call mpi_wait(irequest3, istatus, ierror)
#endif
      !$omp end master
      !$omp barrier

      ! get the force
      !
      if (mod(id,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = if_upper_recv(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            force(1:3,ix,ic) = force(1:3,ix,ic) + buf_recv(3*k-2:3*k,1)
          end do
        end do
      end if

      ! get the force
      !
      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = if_lower_recv(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            force(1:3,ix,ic) = force(1:3,ix,ic) + buf_recv(3*k-2:3*k,2)
          end do
        end do
      end if

      !$omp barrier

    end do

    !$omp end parallel

    return

  end subroutine communicate_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_coor
  !> @brief        Pack the coordinate data and trasfer it
  !! @authors      JJ
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_coor(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, k, ii, ic, ix, id, omp_get_thread_num
#ifdef HAVE_MPI_GENESIS
    integer                  :: num_req
    integer                  :: ireq(4)
    integer                  :: istat(mpi_status_size,4)
#endif
!    integer                  :: irequest, irequest1, irequest2, irequest3
!#ifdef HAVE_MPI_GENESIS
!    integer                  :: istatus(mpi_status_size)
!#endif

    real(wip),       pointer :: coord(:,:,:), buf_send(:,:), buf_recv(:,:)
    integer,         pointer :: natom(:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: upper_send(:), upper_recv(:)
    integer,         pointer :: lower_send(:), lower_recv(:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)


    coord         => domain%coord
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    natom         => domain%num_atom

    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv
    upper_send    => comm%send_coord_upper_size
    upper_recv    => comm%recv_coord_upper_size
    lower_send    => comm%send_coord_lower_size
    lower_recv    => comm%recv_coord_lower_size

    !$omp parallel private(id, ii, k, i, ic, ix)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do ii = 1, 3

      ! Pack the coordinate data (lower)
      !
      if (mod(id,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = ic_lower_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(3*k-2:3*k,1) = coord(1:3,ix,ic)
          end do
        end do
      end if

      ! pack the coordinate data (upper)
      !
      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, comm%num_cell(ii)
          ic = ic_upper_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(3*k-2:3*k,2) = coord(1:3,ix,ic)
          end do
        end do
      end if

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      ! send the data(lower)
      !
      num_req = 1
      call mpi_irecv(buf_recv(1,1), 3*upper_recv(ii), mpi_wip_real, &
                     iproc_upper(ii),                               &
                     iproc_upper(ii),  &
                     mpi_comm_city, ireq(num_req), ierror)

      num_req = num_req + 1
      call mpi_isend(buf_send(1,1), 3*lower_send(ii), mpi_wip_real, &
                     iproc_lower(ii),                               &
                     my_city_rank,  &
                     mpi_comm_city, ireq(num_req), ierror)

      ! send the data(upper)
      !
      num_req = num_req + 1
      call mpi_irecv(buf_recv(1,2), 3*lower_recv(ii), mpi_wip_real, &
                     iproc_lower(ii),                               &
                     my_city_rank+nproc_city,  &
                     mpi_comm_city, ireq(num_req), ierror)

      num_req = num_req + 1
      call mpi_isend(buf_send(1,2), 3*upper_send(ii), mpi_wip_real, &
                     iproc_upper(ii),                               &
                     iproc_upper(ii)+nproc_city,  &
                     mpi_comm_city, ireq(num_req), ierror)
      
      call mpi_waitall(num_req, ireq, istat, ierror)
!     call mpi_wait(irequest,  istatus, ierror)
!     call mpi_wait(irequest1, istatus, ierror)
!     call mpi_wait(irequest2, istatus, ierror)
!     call mpi_wait(irequest3, istatus, ierror)
#endif
      !$omp end master
      !$omp barrier

      ! get the coordinate data (from upper)
      !
      if (mod(id,nthread) == 0) then
      k = 0
      do i = 1, comm%num_cell(ii)
        ic = ic_upper_recv(i,ii)
        do ix = 1, natom(ic)
          k = k + 1
          coord(1:3,ix,ic) = buf_recv(3*k-2:3*k,1)
        end do
      end do
      end if

      ! get the coordinate data (from lower)
      !
      if (mod(id+1,nthread) == 0) then
      k = 0
      do i = 1, comm%num_cell(ii)
        ic = ic_lower_recv(i,ii)
        do ix = 1, natom(ic)
          k = k + 1
          coord(1:3,ix,ic) = buf_recv(3*k-2:3*k,2)
        end do
      end do
      end if

      !$omp barrier

    end do

    !$omp end parallel

    return

  end subroutine communicate_coor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_ptl
  !> @brief        Pack the incoming particles data of each boundary cell
  !! @authors      JJ
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_ptl(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: send_size(2,2), recv_size(2,2)
    integer                  :: i, j, k, ii, ic, ix, nadd
    integer                  :: irequest, irequest1, irequest2, irequest3
    integer                  :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: buf_send(:,:), buf_recv(:,:), buf_real(:,:,:)
    integer,         pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,         pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)
    integer,         pointer :: ptl_add(:)
    integer,         pointer :: int_send(:,:), int_recv(:,:), buf_int(:,:,:)


    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    ptl_add       => domain%ptl_add
    buf_int       => domain%buf_integer
    buf_real      => domain%buf_real

    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    int_send      => comm%int_send
    int_recv      => comm%int_recv
    if_lower_send => comm%if_lower_send
    if_upper_send => comm%if_upper_send
    if_lower_recv => comm%if_lower_recv
    if_upper_recv => comm%if_upper_recv

    do ii = 3, 1, -1

      ! Pack outgoing data (lower)
      !
      k = 0

      do i = 1, comm%num_cell(ii)

        ic = if_lower_send(i,ii)
        int_send(i,1) = ptl_add(ic)

        do ix = 1, ptl_add(ic)

          k = k + 1
          do j = 1, 8
            buf_send(8*(k-1)+j,1) = buf_real(j,ix,ic)
          end do

          int_send(2*k-1+comm%num_cell(ii),1) = buf_int(1,ix,ic)
          int_send(2*k  +comm%num_cell(ii),1) = buf_int(2,ix,ic)

        end do
      end do

      send_size(1,1) = 2*k + comm%num_cell(ii)
      send_size(2,1) = 8*k

      ! pack the outgoing date (upper)
      !
      k = 0

      do i = 1, comm%num_cell(ii)

        ic = if_upper_send(i,ii)
        int_send(i,2) = ptl_add(ic)

        do ix = 1, ptl_add(ic)

          k = k + 1
          do j = 1, 8
            buf_send(8*(k-1)+j,2) = buf_real(j,ix,ic)
          end do

          int_send(2*k-1+comm%num_cell(ii),2) = buf_int(1,ix,ic)
          int_send(2*k  +comm%num_cell(ii),2) = buf_int(2,ix,ic)

        end do
      end do

      send_size(1,2) = 2*k + comm%num_cell(ii)
      send_size(2,2) = 8*k

#ifdef HAVE_MPI_GENESIS
      ! send the size of data
      !
      call mpi_irecv(recv_size(1,1), 2, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii),       &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(send_size(1,1), 2, mpi_integer, iproc_lower(ii), &
                     my_city_rank,       &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(recv_size(1,2), 2, mpi_integer, iproc_lower(ii), &
                     my_city_rank+nproc_city,       &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(1,2), 2, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii)+nproc_city,       &
                     mpi_comm_city,  irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,1), recv_size(1,1), mpi_integer,      &
                     iproc_upper(ii),                                 &
                     iproc_upper(ii),       &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), send_size(1,1), mpi_integer,      &
                     iproc_lower(ii),                                 &
                     my_city_rank,       &
                     mpi_comm_city, irequest1, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,2), recv_size(1,2), mpi_integer,      &
                     iproc_lower(ii),                                 &
                     my_city_rank+nproc_city,       &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), send_size(1,2), mpi_integer,      &
                     iproc_upper(ii),                                 &
                     iproc_upper(ii)+nproc_city,       &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)


      call mpi_irecv(buf_recv(1,1), recv_size(2,1), mpi_wip_real,   &
                     iproc_upper(ii),                               &
                     iproc_upper(ii)+2*nproc_city,     &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), send_size(2,1), mpi_wip_real,   &
                     iproc_lower(ii),                               &
                     my_city_rank+2*nproc_city,     &
                     mpi_comm_city, irequest5, ierror)

      call mpi_irecv(buf_recv(1,2), recv_size(2,2), mpi_wip_real,   &
                     iproc_lower(ii),                               &
                     my_city_rank+3*nproc_city,     &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), send_size(2,2), mpi_wip_real,   &
                     iproc_upper(ii),                               &
                     iproc_upper(ii)+3*nproc_city,     &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the imcoming data
      !
      k = 0

      do i = 1, comm%num_cell(ii)

        ic = if_upper_recv(i,ii)
        nadd = int_recv(i,1)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do j = 1, 8
            buf_real(j,ix,ic) = buf_recv(8*(k-1)+j,1)
          end do

          buf_int(1,ix,ic) = int_recv(2*k-1+comm%num_cell(ii),1)
          buf_int(2,ix,ic) = int_recv(2*k  +comm%num_cell(ii),1)

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

      ! get the imcoming data
      !
      k = 0
      do i = 1, comm%num_cell(ii)

        ic = if_lower_recv(i,ii)
        nadd = int_recv(i,2)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do j = 1, 8
            buf_real(j,ix,ic) = buf_recv(8*(k-1)+j,2)
          end do

          buf_int(1,ix,ic) = int_recv(2*k-1+comm%num_cell(ii),2)
          buf_int(2,ix,ic) = int_recv(2*k  +comm%num_cell(ii),2)

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

    end do

    return

  end subroutine communicate_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_constraints
  !> @brief        Pack the incoming particles data of each boundary cell
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_constraints(domain, comm, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_comm),        target, intent(inout) :: comm
    type(s_constraints), target, intent(inout) :: constraints

    ! local variable
    integer                      :: send_size(2,2), recv_size(2,2)
    integer                      :: real_size, int_size
    integer                      :: i, j, k, l, ii, ic, ix, nadd, water_atom
    integer                      :: list, list1, size
    integer                      :: irequest, irequest1, irequest2, irequest3
    integer                      :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                      :: istatus(mpi_status_size)
#endif
    real(wip),           pointer :: buf_send(:,:), buf_recv(:,:)
    real(wip),           pointer :: buf_real(:,:,:), water_move_real(:,:,:)
    real(wip),           pointer :: HGr_move_real(:,:,:,:)
    integer,             pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,             pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,             pointer :: iproc_upper(:), iproc_lower(:), num_cell(:)
    integer,             pointer :: ptl_add(:), water_move(:)
    integer,             pointer :: buf_int(:,:,:), water_move_int(:,:,:)
    integer,             pointer :: HGr_move(:,:), HGr_move_int(:,:,:,:)
    integer,             pointer :: int_send(:,:), int_recv(:,:)


    iproc_upper     => domain%iproc_upper
    iproc_lower     => domain%iproc_lower
    ptl_add         => domain%ptl_add
    buf_real        => domain%buf_real
    buf_int         => domain%buf_integer
    water_move      => domain%water%move
    water_move_real => domain%water%move_real
    water_move_int  => domain%water%move_integer

    buf_send        => comm%buf_send
    buf_recv        => comm%buf_recv
    int_send        => comm%int_send
    int_recv        => comm%int_recv
    if_lower_send   => comm%if_lower_send
    if_upper_send   => comm%if_upper_send
    if_lower_recv   => comm%if_lower_recv
    if_upper_recv   => comm%if_upper_recv
    num_cell        => comm%num_cell

    HGr_move        => constraints%HGr_move
    HGr_move_int    => constraints%HGr_move_int
    HGr_move_real   => constraints%HGr_move_real

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    size = constraints%connect + 2
    do ii = 3, 1, -1

      ! Pack outgoing data (lower)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)

        int_send(size*(i-1)+1,1) = ptl_add(ic)
        int_send(size*(i-1)+2,1) = water_move(ic)

        do j = 1, constraints%connect
          int_send(size*(i-1)+j+2,1) = HGr_move(j,ic)
        end do

        do ix = 1, ptl_add(ic)

          k = k + 1
          do list = 1, 8
            buf_send(8*(k-1)+list,1) = buf_real(list,ix,ic)
          end do

          do list = 1, 3
            list1 = 3*(k-1) + list + size*num_cell(ii)
            int_send(list1,1) = buf_int(list,ix,ic)
          end do

        end do

      end do

      real_size = 8*k
      int_size = 3*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)
        do j = 1, constraints%connect
          do ix = 1, HGr_move(j,ic)
            do l = 1, j+1

              k = k + 1
              do list = 1, 9
                buf_send(real_size+9*(k-1)+list,1) = &
                     HGr_move_real(9*(l-1)+list,ix,j,ic)
              end do

              do list = 1, 3
                list1 = 3*(k-1) + list + int_size
                int_send(list1,1) = HGr_move_int(3*(l-1)+list,ix,j,ic)
              end do

            end do
          end do
        end do

      end do

      real_size = real_size + 9*k
      int_size = int_size + 3*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)
        do ix = 1, water_move(ic)

          k = k + 1
          do j = 1, 6*water_atom
            buf_send(real_size+6*water_atom*(k-1)+j,1) &
                            = water_move_real(j,ix,ic)
          end do

          do j = 1, water_atom
            int_send(int_size+water_atom*(k-1)+j,1)    &
                          = water_move_int(j,ix,ic)
          end do

        end do
      end do

      send_size(1,1) = int_size + water_atom*k
      send_size(2,1) = real_size + 6*water_atom*k

      ! pack the outgoing date (upper)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)

        int_send(size*(i-1)+1,2) = ptl_add(ic)
        int_send(size*(i-1)+2,2) = water_move(ic)

        do j = 1, constraints%connect
          int_send(size*(i-1)+j+2,2) = HGr_move(j,ic)
        end do

        do ix = 1, ptl_add(ic)

          k = k + 1
          do list = 1, 8
            buf_send(8*(k-1)+list,2) = buf_real(list,ix,ic)
          end do

          do list = 1, 3
            list1 = 3*(k-1) + list + size*num_cell(ii)
            int_send(list1,2) = buf_int(list,ix,ic)
          end do

        end do

      end do

      real_size = 8*k
      int_size = 3*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)
        do j = 1, constraints%connect
          do ix = 1, HGr_move(j,ic)

            do l = 1, j+1

              k = k + 1
              do list = 1, 9
                buf_send(real_size+9*(k-1)+list,2) = &
                     HGr_move_real(9*(l-1)+list,ix,j,ic)
              end do

              do list = 1, 3
                list1 = 3*(k-1) + list + int_size
                int_send(list1,2) = HGr_move_int(3*(l-1)+list,ix,j,ic)
              end do

            end do
          end do
        end do

      end do

      real_size = real_size + 9*k
      int_size  = int_size + 3*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)
        do ix = 1, water_move(ic)

          k = k + 1
          do j = 1, 6*water_atom
            buf_send(real_size+6*water_atom*(k-1)+j,2) &
                           =  water_move_real(j,ix,ic)
          end do

          do j = 1, water_atom
            int_send(int_size+water_atom*(k-1)+j,2)   &
                          = water_move_int(j,ix,ic)
          end do

        end do
      end do

      send_size(1,2) = int_size + water_atom*k
      send_size(2,2) = real_size + 6*water_atom*k
#ifdef HAVE_MPI_GENESIS

      ! send the size of data
      !
      call mpi_irecv(recv_size(1,1), 2, mpi_integer, iproc_upper(ii),     &
                     iproc_upper(ii),                                     &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(send_size(1,1), 2, mpi_integer, iproc_lower(ii),     &
                     my_city_rank,                                        &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(recv_size(1,2), 2, mpi_integer, iproc_lower(ii),     &
                     my_city_rank+nproc_city,                             &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(1,2), 2, mpi_integer, iproc_upper(ii),     &
                     iproc_upper(ii)+nproc_city,                          &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,1), recv_size(1,1), mpi_integer,          &
                     iproc_upper(ii),                                     &
                     iproc_upper(ii),                                     &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), send_size(1,1), mpi_integer,          &
                     iproc_lower(ii),                                     &
                     my_city_rank,                                        &
                     mpi_comm_city, irequest1, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,2), recv_size(1,2), mpi_integer,          &
                     iproc_lower(ii),                                     &
                     my_city_rank+nproc_city,                             &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), send_size(1,2), mpi_integer,          &
                     iproc_upper(ii),                                     &
                     iproc_upper(ii)+nproc_city,                          &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      call mpi_irecv(buf_recv(1,1), recv_size(2,1), mpi_wip_real,       &
                     iproc_upper(ii),                                   &
                     iproc_upper(ii)+2*nproc_city,                      &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), send_size(2,1), mpi_wip_real,       &
                     iproc_lower(ii),                                   &
                     my_city_rank+2*nproc_city,                         &
                     mpi_comm_city, irequest5, ierror)

      call mpi_irecv(buf_recv(1,2), recv_size(2,2), mpi_wip_real,       &
                     iproc_lower(ii),                                   &
                     my_city_rank+3*nproc_city,                         &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), send_size(2,2), mpi_wip_real,       &
                     iproc_upper(ii),                                   &
                     iproc_upper(ii)+3*nproc_city,                      &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the imcoming data (from upper)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        nadd = int_recv(size*(i-1)+1,1)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do list = 1, 8
            buf_real(list,ix,ic) = buf_recv(8*(k-1)+list,1)
          end do

          do list = 1, 3
            list1 = 3*(k-1) + list + size*num_cell(ii)
            buf_int(list,ix,ic) = int_recv(list1,1)
          end do

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

      real_size = 8*k
      int_size = 3*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        do j = 1, constraints%connect

          nadd = int_recv(size*(i-1)+j+2,1)
          do ix = HGr_move(j,ic)+1, HGr_move(j,ic)+nadd

            do l = 1, j+1

              k = k + 1
              do list = 1, 9
                HGr_move_real(9*(l-1)+list,ix,j,ic) = &
                     buf_recv(real_size+9*(k-1)+list,1)
              end do

              do list = 1, 3
                list1 = 3*(k-1) + list + int_size
                HGr_move_int(3*(l-1)+list,ix,j,ic) = int_recv(list1,1)
              end do

            end do
          end do

          HGr_move(j,ic) = HGr_move(j,ic) + nadd

        end do
      end do

      real_size = real_size + 9*k
      int_size  = int_size + 3*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        nadd = int_recv(size*(i-1)+2,1)

        do ix = water_move(ic)+1, water_move(ic)+nadd

          k = k + 1
          do j = 1, 6*water_atom
            water_move_real(j,ix,ic) &
              = buf_recv(real_size+6*water_atom*(k-1)+j,1)
          end do

          do j = 1, water_atom
            water_move_int(j,ix,ic) &
            = int_recv(int_size+water_atom*(k-1)+j,1)
          end do

        end do

        water_move(ic) = water_move(ic) + nadd

      end do

      ! get the imcoming data (from lower)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        nadd = int_recv(size*(i-1)+1,2)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do list = 1, 8
            buf_real(list,ix,ic) = buf_recv(8*(k-1)+list,2)
          end do

          do list = 1, 3
            list1 = 3*(k-1) + list + size*num_cell(ii)
            buf_int(list,ix,ic) = int_recv(list1,2)
          end do

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

      real_size = 8*k
      int_size = 3*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        do j = 1, constraints%connect

          nadd = int_recv(size*(i-1)+j+2,2)
          do ix = HGr_move(j,ic)+1, HGr_move(j,ic)+nadd
            do l = 1, j+1

              k = k + 1
              do list = 1, 9
                HGr_move_real(9*(l-1)+list,ix,j,ic) = &
                     buf_recv(real_size+9*(k-1)+list,2)
              end do

              do list = 1, 3
                list1 = 3*(k-1) + list + int_size
                HGr_move_int(3*(l-1)+list,ix,j,ic) = int_recv(list1,2)
              end do

            end do
          end do

          HGr_move(j,ic) = HGr_move(j,ic) + nadd

        end do

      end do

      real_size = real_size + 9*k
      int_size = int_size + 3*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        nadd = int_recv(size*(i-1)+2,2)

        do ix = water_move(ic)+1, water_move(ic)+nadd

          k = k + 1
          do j = 1, 6*water_atom
            water_move_real(j,ix,ic) &
              = buf_recv(real_size+6*water_atom*(k-1)+j,2)
          end do

          do j = 1, water_atom
            water_move_int(j,ix,ic) &
               = int_recv(int_size+water_atom*(k-1)+j,2)
          end do

        end do

        water_move(ic) = water_move(ic) + nadd

      end do

    end do

    return

  end subroutine communicate_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_bond
  !> @brief        Pack the incoming bonding data of each boundary cell
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_bond(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: send_size(2,2),recv_size(2,2)
    integer                   :: i, j, k, ii, ic, ix, nadd
    integer                   :: irequest, irequest1, irequest2, irequest3
    integer                   :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:,:), buf_recv(:,:)
    real(wp),         pointer :: buf_bond_real(:,:,:), buf_angle_real(:,:,:)
    real(wp),         pointer :: buf_dihed_real(:,:,:), buf_impr_real(:,:,:)
    real(wp),         pointer :: buf_rb_dihed_real(:,:,:)
    real(wp),         pointer :: buf_restraint_real(:,:,:)
    real(wp),         pointer :: buf_fitting_real(:,:,:)
    integer,          pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,          pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,          pointer :: iproc_upper(:), iproc_lower(:)
    integer,          pointer :: int_send(:,:), int_recv(:,:)
    integer,          pointer :: bond_add(:), buf_bond_integer(:,:,:)
    integer,          pointer :: angle_add(:), buf_angle_integer(:,:,:)
    integer,          pointer :: dihed_add(:), buf_dihed_integer(:,:,:)
    integer,          pointer :: rb_dihed_add(:), buf_rb_dihed_integer(:,:,:)
    integer,          pointer :: impr_add(:), buf_impr_integer(:,:,:)
    integer,          pointer :: cmap_add(:), buf_cmap_integer(:,:,:)
    integer,          pointer :: restraint_add(:), buf_restraint_integer(:,:)
    integer,          pointer :: fitting_add(:), buf_fitting_integer(:,:)


    iproc_upper           => domain%iproc_upper
    iproc_lower           => domain%iproc_lower

    buf_send              => comm%buf_send
    buf_recv              => comm%buf_recv
    if_lower_send         => comm%if_lower_send
    if_upper_send         => comm%if_upper_send
    if_lower_recv         => comm%if_lower_recv
    if_upper_recv         => comm%if_upper_recv
    int_send              => comm%int_send
    int_recv              => comm%int_recv

    bond_add              => enefunc%bond_add
    buf_bond_real         => enefunc%buf_bond_real
    buf_bond_integer      => enefunc%buf_bond_integer
    angle_add             => enefunc%angle_add
    buf_angle_real        => enefunc%buf_angle_real
    buf_angle_integer     => enefunc%buf_angle_integer
    dihed_add             => enefunc%dihed_add
    buf_dihed_real        => enefunc%buf_dihed_real
    buf_dihed_integer     => enefunc%buf_dihed_integer
    rb_dihed_add          => enefunc%rb_dihed_add
    buf_rb_dihed_real     => enefunc%buf_rb_dihed_real
    buf_rb_dihed_integer  => enefunc%buf_rb_dihed_integer
    impr_add              => enefunc%impr_add
    buf_impr_real         => enefunc%buf_impr_real
    buf_impr_integer      => enefunc%buf_impr_integer
    cmap_add              => enefunc%cmap_add
    buf_cmap_integer      => enefunc%buf_cmap_integer
    restraint_add         => enefunc%restraint_add
    buf_restraint_real    => enefunc%buf_restraint_real
    buf_restraint_integer => enefunc%buf_restraint_integer
    fitting_add           => enefunc%fitting_add
    buf_fitting_real      => enefunc%buf_fitting_real
    buf_fitting_integer   => enefunc%buf_fitting_integer


    do ii = 3, 1, -1

      ! Pack outgoing data (lower)
      !

      k = 8*comm%num_cell(ii)
      j = 0

      do i = 1, comm%num_cell(ii)

        ic = if_lower_send(i,ii)
        int_send(8*i-7,1) = bond_add(ic)
        int_send(8*i-6,1) = angle_add(ic)
        int_send(8*i-5,1) = dihed_add(ic)
        int_send(8*i-4,1) = rb_dihed_add(ic)
        int_send(8*i-3,1) = impr_add(ic)
        int_send(8*i-2,1) = cmap_add(ic)
        int_send(8*i-1,1) = restraint_add(ic)
        int_send(8*i  ,1) = fitting_add(ic)

        do ix = 1, bond_add(ic)
          int_send(k+1:k+4,1) = buf_bond_integer(1:4,ix,ic)
          k = k + 4
          buf_send(j+1:j+2,1) = buf_bond_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, angle_add(ic)
          int_send(k+1:k+7,1) = buf_angle_integer(1:7,ix,ic)
          k = k + 7
          buf_send(j+1:j+4,1) = buf_angle_real(1:4,ix,ic)
          j = j + 4
        end do

        do ix = 1, dihed_add(ic)
          int_send(k+1:k+9,1) = buf_dihed_integer(1:9,ix,ic)
          k = k + 9
          buf_send(j+1:j+2,1) = buf_dihed_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, rb_dihed_add(ic)
          int_send(k+1:k+7,1) = buf_rb_dihed_integer(1:7,ix,ic)
          k = k + 7
          buf_send(j+1:j+6,1) = buf_rb_dihed_real(1:6,ix,ic)
          j = j + 6
        end do

        do ix = 1, impr_add(ic)
          int_send(k+1:k+8,1) = buf_impr_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+2,1) = buf_impr_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, cmap_add(ic)
          int_send(k+1:k+15,1) = buf_cmap_integer(1:15,ix,ic)
          k = k + 15
        end do

        do ix = 1, restraint_add(ic)
          int_send(k+1,1) = buf_restraint_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+7,1) = buf_restraint_real(1:7,ix,ic)
          j = j + 7
        end do

        do ix = 1, fitting_add(ic)
          int_send(k+1,1) = buf_fitting_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+3,1) = buf_fitting_real(1:3,ix,ic)
          j = j + 3
        end do

      end do

      send_size(1,1) = k
      send_size(2,1) = j

      ! pack the outgoing date (upper)
      !
      k = 8*comm%num_cell(ii)
      j = 0

      do i = 1, comm%num_cell(ii)

        ic = if_upper_send(i,ii)
        int_send(8*i-7,2) = bond_add(ic)
        int_send(8*i-6,2) = angle_add(ic)
        int_send(8*i-5,2) = dihed_add(ic)
        int_send(8*i-4,2) = rb_dihed_add(ic)
        int_send(8*i-3,2) = impr_add(ic)
        int_send(8*i-2,2) = cmap_add(ic)
        int_send(8*i-1,2) = restraint_add(ic)
        int_send(8*i  ,2) = fitting_add(ic)

        do ix = 1, bond_add(ic)
          int_send(k+1:k+4,2) = buf_bond_integer(1:4,ix,ic)
          k = k + 4
          buf_send(j+1:j+2,2) = buf_bond_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, angle_add(ic)
          int_send(k+1:k+7,2) = buf_angle_integer(1:7,ix,ic)
          k = k + 7
          buf_send(j+1:j+4,2) = buf_angle_real(1:4,ix,ic)
          j = j + 4
        end do

        do ix = 1, dihed_add(ic)
          int_send(k+1:k+9,2) = buf_dihed_integer(1:9,ix,ic)
          k = k + 9
          buf_send(j+1:j+2,2) = buf_dihed_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, rb_dihed_add(ic)
          int_send(k+1:k+7,2) = buf_rb_dihed_integer(1:7,ix,ic)
          k = k + 7
          buf_send(j+1:j+6,2) = buf_rb_dihed_real(1:6,ix,ic)
          j = j + 6
        end do

        do ix = 1, impr_add(ic)
          int_send(k+1:k+8,2) = buf_impr_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+2,2) = buf_impr_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, cmap_add(ic)
          int_send(k+1:k+15,2) = buf_cmap_integer(1:15,ix,ic)
          k = k + 15
        end do

        do ix = 1, restraint_add(ic)
          int_send(k+1,2) = buf_restraint_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+7,2) = buf_restraint_real(1:7,ix,ic)
          j = j + 7
        end do

        do ix = 1, fitting_add(ic)
          int_send(k+1,2) = buf_fitting_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+3,2) = buf_fitting_real(1:3,ix,ic)
          j = j + 3
        end do

      end do

      send_size(1,2) = k
      send_size(2,2) = j

#ifdef HAVE_MPI_GENESIS

       ! send the size of data
       !
       call mpi_irecv(recv_size(1,1), 2, mpi_integer, iproc_upper(ii),     &
                      iproc_upper(ii),                                     &
                      mpi_comm_city, irequest, ierror)
       call mpi_isend(send_size(1,1), 2, mpi_integer, iproc_lower(ii),     &
                      my_city_rank,                                        &
                      mpi_comm_city, irequest1, ierror)
       call mpi_irecv(recv_size(1,2), 2, mpi_integer, iproc_lower(ii),     &
                      my_city_rank+nproc_city,                             &
                      mpi_comm_city, irequest2, ierror)
       call mpi_isend(send_size(1,2), 2, mpi_integer, iproc_upper(ii),     &
                      iproc_upper(ii)+nproc_city,                          &
                      mpi_comm_city, irequest3, ierror)

       call mpi_wait(irequest,  istatus, ierror)
       call mpi_wait(irequest1, istatus, ierror)
       call mpi_wait(irequest2, istatus, ierror)
       call mpi_wait(irequest3, istatus, ierror)

       ! send the integer data (lower)
       !
       call mpi_irecv(int_recv(1,1), recv_size(1,1), mpi_integer,          &
                      iproc_upper(ii),                                     &
                      iproc_upper(ii),                                     &
                      mpi_comm_city, irequest, ierror)
       call mpi_isend(int_send(1,1), send_size(1,1), mpi_integer,          &
                      iproc_lower(ii),                                     &
                      my_city_rank,                                        &
                      mpi_comm_city, irequest1, ierror)

       ! send the integer data (upper)
       !
       call mpi_irecv(int_recv(1,2), recv_size(1,2), mpi_integer,          &
                      iproc_lower(ii),                                     &
                      my_city_rank+nproc_city,                             &
                      mpi_comm_city, irequest2, ierror)
       call mpi_isend(int_send(1,2), send_size(1,2), mpi_integer,          &
                      iproc_upper(ii),                                     &
                      iproc_upper(ii)+nproc_city,                          &
                      mpi_comm_city, irequest3, ierror)

       call mpi_wait(irequest,  istatus, ierror)
       call mpi_wait(irequest1, istatus, ierror)
       call mpi_wait(irequest2, istatus, ierror)
       call mpi_wait(irequest3, istatus, ierror)

       ! send the real data (lower)
       !
       call mpi_irecv(buf_recv(1,1), recv_size(2,1), mpi_wip_real,      &
                      iproc_upper(ii),                                  &
                      iproc_upper(ii),                                  &
                      mpi_comm_city, irequest4, ierror)
       call mpi_isend(buf_send(1,1), send_size(2,1), mpi_wip_real,      &
                      iproc_lower(ii),                                  &
                      my_city_rank,                                     &
                      mpi_comm_city, irequest5, ierror)

       ! send the real data
       !
       call mpi_irecv(buf_recv(1,2), recv_size(2,2), mpi_wip_real,       &
                      iproc_lower(ii),                                   &
                      my_city_rank+nproc_city,                           &
                      mpi_comm_city, irequest6, ierror)
       call mpi_isend(buf_send(1,2), send_size(2,2), mpi_wip_real,       &
                      iproc_upper(ii),                                   &
                      iproc_upper(ii)+nproc_city,                        &
                      mpi_comm_city, irequest7, ierror)

       call mpi_wait(irequest4, istatus, ierror)
       call mpi_wait(irequest5, istatus, ierror)
       call mpi_wait(irequest6, istatus, ierror)
       call mpi_wait(irequest7, istatus, ierror)
#endif

       ! get the imcoming data
       !

       k = 8*comm%num_cell(ii)
       j = 0

       do i = 1, comm%num_cell(ii)
          ic = if_upper_recv(i,ii)

          nadd = int_recv(8*i-7,1)
          do ix = bond_add(ic)+1, bond_add(ic)+nadd
             buf_bond_integer(1:4,ix,ic) = int_recv(k+1:k+4,1)
             k = k + 4
             buf_bond_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          bond_add(ic) = bond_add(ic) + nadd

          nadd = int_recv(8*i-6,1)
          do ix = angle_add(ic)+1, angle_add(ic)+nadd
             buf_angle_integer(1:7,ix,ic) = int_recv(k+1:k+7,1)
             k = k + 7
             buf_angle_real(1:4,ix,ic) = buf_recv(j+1:j+4,1)
             j = j + 4
          end do
          angle_add(ic) = angle_add(ic) + nadd

          nadd = int_recv(8*i-5,1)
          do ix = dihed_add(ic)+1, dihed_add(ic)+nadd
             buf_dihed_integer(1:9,ix,ic) = int_recv(k+1:k+9,1)
             k = k + 9
             buf_dihed_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          dihed_add(ic) = dihed_add(ic) + nadd

          nadd = int_recv(8*i-4,1)
          do ix = rb_dihed_add(ic)+1, rb_dihed_add(ic)+nadd
             buf_rb_dihed_integer(1:7,ix,ic) = int_recv(k+1:k+7,1)
             k = k + 7
             buf_rb_dihed_real(1:6,ix,ic) = buf_recv(j+1:j+6,1)
             j = j + 6
          end do
          rb_dihed_add(ic) = rb_dihed_add(ic) + nadd

          nadd = int_recv(8*i-3,1)
          do ix = impr_add(ic)+1, impr_add(ic)+nadd
             buf_impr_integer(1:8,ix,ic) = int_recv(k+1:k+8,1)
             k = k + 8
             buf_impr_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          impr_add(ic) = impr_add(ic) + nadd

          nadd = int_recv(8*i-2,1)
          do ix = cmap_add(ic)+1, cmap_add(ic)+nadd
             buf_cmap_integer(1:15,ix,ic) = int_recv(k+1:k+15,1)
             k = k + 15
          end do
          cmap_add(ic) = cmap_add(ic) + nadd

          nadd = int_recv(8*i-1,1)
          do ix = restraint_add(ic)+1, restraint_add(ic)+nadd
             buf_restraint_integer(ix,ic) = int_recv(k+1,1)
             k = k + 1
             buf_restraint_real(1:7,ix,ic) = buf_recv(j+1:j+7,1)
             j = j + 7
          end do
          restraint_add(ic) = restraint_add(ic) + nadd

          nadd = int_recv(8*i,1)
          do ix = fitting_add(ic)+1, fitting_add(ic)+nadd
             buf_fitting_integer(ix,ic) = int_recv(k+1,1)
             k = k + 1
             buf_fitting_real(1:3,ix,ic) = buf_recv(j+1:j+3,1)
             j = j + 3
          end do
          fitting_add(ic) = fitting_add(ic) + nadd

       end do

       ! get the incoming data
       !

       k = 8*comm%num_cell(ii)
       j = 0

       do i = 1, comm%num_cell(ii)

          ic = if_lower_recv(i,ii)

          nadd = int_recv(8*i-7,2)
          do ix = bond_add(ic)+1, bond_add(ic)+nadd
             buf_bond_integer(1:4,ix,ic) = int_recv(k+1:k+4,2)
             k = k + 4
             buf_bond_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          bond_add(ic) = bond_add(ic) + nadd

          nadd = int_recv(8*i-6,2)
          do ix = angle_add(ic)+1, angle_add(ic)+nadd
             buf_angle_integer(1:7,ix,ic) = int_recv(k+1:k+7,2)
             k = k + 7
             buf_angle_real(1:4,ix,ic) = buf_recv(j+1:j+4,2)
             j = j + 4
          end do
          angle_add(ic) = angle_add(ic) + nadd

          nadd = int_recv(8*i-5,2)
          do ix = dihed_add(ic)+1, dihed_add(ic)+nadd
             buf_dihed_integer(1:9,ix,ic) = int_recv(k+1:k+9,2)
             k = k + 9
             buf_dihed_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          dihed_add(ic) = dihed_add(ic) + nadd

          nadd = int_recv(8*i-4,2)
          do ix = rb_dihed_add(ic)+1, rb_dihed_add(ic)+nadd
             buf_rb_dihed_integer(1:7,ix,ic) = int_recv(k+1:k+7,2)
             k = k + 7
             buf_rb_dihed_real(1:6,ix,ic) = buf_recv(j+1:j+6,2)
             j = j + 6
          end do
          rb_dihed_add(ic) = rb_dihed_add(ic) + nadd

          nadd = int_recv(8*i-3,2)
          do ix = impr_add(ic)+1, impr_add(ic)+nadd
             buf_impr_integer(1:8,ix,ic) = int_recv(k+1:k+8,2)
             k = k + 8
             buf_impr_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          impr_add(ic) = impr_add(ic) + nadd

          nadd = int_recv(8*i-2,2)
          do ix = cmap_add(ic)+1, cmap_add(ic)+nadd
             buf_cmap_integer(1:15,ix,ic) = int_recv(k+1:k+15,2)
             k = k + 15
          end do
          cmap_add(ic) = cmap_add(ic) + nadd

          nadd = int_recv(8*i-1,2)
          do ix = restraint_add(ic)+1, restraint_add(ic)+nadd
             k = k + 1
             buf_restraint_integer(ix,ic) = int_recv(k,2)
             buf_restraint_real(1:7,ix,ic) = buf_recv(j+1:j+7,2)
             j = j + 7
          end do
          restraint_add(ic) = restraint_add(ic) + nadd

          nadd = int_recv(8*i,2)
          do ix = fitting_add(ic)+1, fitting_add(ic)+nadd
             k = k + 1
             buf_fitting_integer(ix,ic) = int_recv(k,2)
             buf_fitting_real(1:3,ix,ic) = buf_recv(j+1:j+3,2)
             j = j + 3
          end do
          fitting_add(ic) = fitting_add(ic) + nadd

       end do

    end do

    return

  end subroutine communicate_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    updae_cell_boundary
  !> @brief        Update coordinates of the boundary cell with constraint
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_boundary(domain, comm, boundary, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_comm),        target, intent(inout) :: comm
    type(s_boundary),    target, intent(inout) :: boundary
    type(s_constraints), target, intent(inout) :: constraints

    ! local variable
    real(wip)                :: x_shift, y_shift, z_shift
    real(wip)                :: move(3)
    integer                  :: i, k, j, ii, ic, ix, iwater, list
    integer                  :: irequest, irequest1, irequest2, irequest3
    integer                  :: irequest4, irequest5, irequest6, irequest7
    integer                  :: id, omp_get_thread_num
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),        pointer :: trans(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: cell_pbc_move(:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: buf_send(:,:), buf_recv(:,:)
    real(wip),       pointer :: bsize_x, bsize_y, bsize_z
    integer,         pointer :: num_cell(:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)
    integer,         pointer :: upper_recv(:), lower_recv(:)
    integer,         pointer :: upper_send(:), lower_send(:)
    integer,         pointer :: int_send(:,:), int_recv(:,:)
    integer,         pointer :: ncell_local, ncell_bd, natom(:), nsolute(:)
    integer,         pointer :: nwater(:), water_list(:,:,:), atmcls(:,:)
    integer,         pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: No_HGr(:), HGr_local(:,:),HGr_bond_list(:,:,:,:)


    ncell_local   => domain%num_cell_local
    ncell_bd      => domain%num_cell_boundary
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water
    water_list    => domain%water_list
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    inv_mass      => domain%inv_mass
    atmcls        => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    id_g2l        => domain%id_g2l
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    trans         => domain%trans_vec
    cell_pbc_move => domain%cell_pbc_move

    num_cell      => comm%num_cell
    int_send      => comm%int_send
    int_recv      => comm%int_recv
    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv
    upper_recv    => comm%recv_coord_upper_size
    lower_recv    => comm%recv_coord_lower_size
    upper_send    => comm%send_coord_upper_size
    lower_send    => comm%send_coord_lower_size

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list

    do ii = 1, 3

      !$omp parallel private(k, i, ic, ix, id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      ! Pack the data (lower)
      !
      if (id == 0) then
        k = 0
        do i = 1, num_cell(ii)

          ic = ic_lower_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(8*k-7,1) = coord(1,ix,ic)
            buf_send(8*k-6,1) = coord(2,ix,ic)
            buf_send(8*k-5,1) = coord(3,ix,ic)
            buf_send(8*k-4,1) = velocity(1,ix,ic)
            buf_send(8*k-3,1) = velocity(2,ix,ic)
            buf_send(8*k-2,1) = velocity(3,ix,ic)
            buf_send(8*k-1,1) = charge(ix,ic)
            buf_send(8*k,1)   = mass(ix,ic)
            int_send(3*k-2,1) = atmcls(ix,ic)
            int_send(3*k-1,1) = id_l2g(ix,ic)
            int_send(3*k  ,1) = id_l2g_sol(ix,ic)
          end do

        end do
      end if

      ! pack the data (upper)
      !
      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, num_cell(ii)

          ic = ic_upper_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(8*k-7,2) = coord(1,ix,ic)
            buf_send(8*k-6,2) = coord(2,ix,ic)
            buf_send(8*k-5,2) = coord(3,ix,ic)
            buf_send(8*k-4,2) = velocity(1,ix,ic)
            buf_send(8*k-3,2) = velocity(2,ix,ic)
            buf_send(8*k-2,2) = velocity(3,ix,ic)
            buf_send(8*k-1,2) = charge(ix,ic)
            buf_send(8*k,2)   = mass(ix,ic)
            int_send(3*k-2,2) = atmcls(ix,ic)
            int_send(3*k-1,2) = id_l2g(ix,ic)
            int_send(3*k  ,2) = id_l2g_sol(ix,ic)
          end do

        end do
      end if

      !$omp end parallel

#ifdef HAVE_MPI_GENESIS

      ! send the data
      !
      call mpi_irecv(int_recv(1,1), 3*upper_recv(ii), mpi_integer,           &
                     iproc_upper(ii),                                        &
                     iproc_upper(ii),                                        &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), 3*lower_send(ii), mpi_integer,           &
                     iproc_lower(ii),                                        &
                     my_city_rank,                                           &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(int_recv(1,2), 3*lower_recv(ii), mpi_integer,           &
                     iproc_lower(ii),                                        &
                     my_city_rank+nproc_city,                                &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), 3*upper_send(ii), mpi_integer,           &
                     iproc_upper(ii),                                        &
                     iproc_upper(ii)+nproc_city,                             &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      call mpi_irecv(buf_recv(1,1), 8*upper_recv(ii), mpi_wip_real,        &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii),                                      &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), 8*lower_send(ii), mpi_wip_real,        &
                     iproc_lower(ii),                                      &
                     my_city_rank,                                         &
                     mpi_comm_city, irequest5, ierror)

      call mpi_irecv(buf_recv(1,2), 8*lower_recv(ii), mpi_wip_real,        &
                     iproc_lower(ii),                                      &
                     my_city_rank+nproc_city,                              &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), 8*upper_send(ii), mpi_wip_real,        &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii)+nproc_city,                           &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the data of the boundary cells
      !

      !$omp parallel private(k, i, ic, ix, id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (id == 0) then
      k = 0
      do i = 1, num_cell(ii)

        ic = ic_upper_recv(i,ii)
        do ix = 1, natom(ic)
          k = k + 1
          coord(1,ix,ic)    = buf_recv(8*k-7,1)
          coord(2,ix,ic)    = buf_recv(8*k-6,1)
          coord(3,ix,ic)    = buf_recv(8*k-5,1)
          velocity(1,ix,ic) = buf_recv(8*k-4,1)
          velocity(2,ix,ic) = buf_recv(8*k-3,1)
          velocity(3,ix,ic) = buf_recv(8*k-2,1)
          charge(ix,ic)     = buf_recv(8*k-1,1)
          mass(ix,ic)       = buf_recv(8*k,1)
          atmcls(ix,ic)     = int_recv(3*k-2,1)
          id_l2g(ix,ic)     = int_recv(3*k-1,1)
          id_l2g_sol(ix,ic) = int_recv(3*k  ,1)
        end do

      end do
      end if

      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, num_cell(ii)

          ic = ic_lower_recv(i,ii)

          do ix = 1, natom(ic)
            k = k + 1
            coord(1,ix,ic)    = buf_recv(8*k-7,2)
            coord(2,ix,ic)    = buf_recv(8*k-6,2)
            coord(3,ix,ic)    = buf_recv(8*k-5,2)
            velocity(1,ix,ic) = buf_recv(8*k-4,2)
            velocity(2,ix,ic) = buf_recv(8*k-3,2)
            velocity(3,ix,ic) = buf_recv(8*k-2,2)
            charge(ix,ic)     = buf_recv(8*k-1,2)
            mass(ix,ic)       = buf_recv(8*k,2)
            atmcls(ix,ic)     = int_recv(3*k-2,2)
            id_l2g(ix,ic)     = int_recv(3*k-1,2)
            id_l2g_sol(ix,ic) = int_recv(3*k  ,2)
          end do

        end do
      end if

      !$omp end parallel

    end do

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ix, x_shift, y_shift, z_shift, move, ic, iwater,&
    !$omp         list, ii)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! inverse of mass
    !
    do i = id+1, ncell_local, nthread
      do ix = 1, natom(i)
        inv_mass(ix,i) = 0.0_wip
        if (mass(ix,i) > EPS) inv_mass(ix,i) = 1.0_wip / mass(ix,i)
      end do
    end do

    ! check the translation for each particle
    !
    do i = id+1, ncell_local, nthread

!ocl nosimd
      do ix = 1, No_HGr(i)

        !coordinate shifted against the origin
        !
        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)
        trans(1,ix,i) = move(1) - boundary%origin_x
        trans(2,ix,i) = move(2) - boundary%origin_y
        trans(3,ix,i) = move(3) - boundary%origin_z
        ii = id_l2g_sol(ix,i)
        id_g2l(1,ii) = i
        id_g2l(2,ii) = ix

      end do

      k = No_HGr(i)
      do j = 1, constraints%connect

!ocl nosimd
        do ix = 1, HGr_local(j,i)
          list = HGr_bond_list(1,ix,j,i)
          x_shift = coord(1,list,i) - boundary%origin_x
          y_shift = coord(2,list,i) - boundary%origin_y
          z_shift = coord(3,list,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do list = 1, j+1
            k = k + 1
            trans(1,k,i) = move(1) - boundary%origin_x
            trans(2,k,i) = move(2) - boundary%origin_y
            trans(3,k,i) = move(3) - boundary%origin_z
            ii = id_l2g_sol(k,i)
            id_g2l(1,ii) = i
            id_g2l(2,ii) = k
          end do

        end do
      end do
      nsolute(i) = k

      if (constraints%water_type == TIP4) then

        do ix = nsolute(i)+1, natom(i)-3, 4

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          do j = 1,4
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do

      else if (constraints%water_type == TIP3) then

        do ix = nsolute(i)+1, natom(i)-2, 3

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do j = 1,3
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do
      else if (constraints%water_type == TIP2) then

        do ix = nsolute(i)+1, natom(i)-1, 2

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do j = 1,2
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do

      else if (constraints%water_type == TIP1) then

        do ix = nsolute(i)+1, natom(i)

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          trans(1,ix,i) = move(1) - boundary%origin_x
          trans(2,ix,i) = move(2) - boundary%origin_y
          trans(3,ix,i) = move(3) - boundary%origin_z

        end do

      end if

    end do

    do i = id+1, ncell_bd, nthread
      ic = i + ncell_local

!ocl nosimd
      do ix = 1, No_HGr(ic)

        !coordinate shifted against the origin
        !
        x_shift = coord(1,ix,ic) - boundary%origin_x
        y_shift = coord(2,ix,ic) - boundary%origin_y
        z_shift = coord(3,ix,ic) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_wip-bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip-bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip-bsize_z*anint(z_shift/bsize_z)
        if (domain%nonbond_kernel == NBK_Fugaku .or. &
            domain%nonbond_kernel == NBK_Intel) then
          move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
          move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
          move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
        end if
        trans(1,ix,ic) = move(1) - boundary%origin_x
        trans(2,ix,ic) = move(2) - boundary%origin_y
        trans(3,ix,ic) = move(3) - boundary%origin_z
        ii = id_l2g_sol(ix,ic)
        id_g2l(1,ii) = ic
        id_g2l(2,ii) = ix
      end do

      k = No_HGr(ic)

      do j = 1, constraints%connect

!ocl nosimd
        do ix = 1, HGr_local(j,ic)

          k = k + 1
          HGr_bond_list(1,ix,j,ic) = k
          x_shift = coord(1,k,ic) - boundary%origin_x
          y_shift = coord(2,k,ic) - boundary%origin_y
          z_shift = coord(3,k,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if
          trans(1,k,ic) = move(1) - boundary%origin_x
          trans(2,k,ic) = move(2) - boundary%origin_y
          trans(3,k,ic) = move(3) - boundary%origin_z
          ii = id_l2g_sol(k,ic)
          id_g2l(1,ii) = ic
          id_g2l(2,ii) = k

          do list = 2, j+1
            k = k + 1
            HGr_bond_list(list,ix,j,ic) = k
            trans(1,k,ic) = move(1) - boundary%origin_x
            trans(2,k,ic) = move(2) - boundary%origin_y
            trans(3,k,ic) = move(3) - boundary%origin_z
            ii = id_l2g_sol(k,ic)
            id_g2l(1,ii) = ic
            id_g2l(2,ii) = k
          end do

        end do
      end do
      nsolute(ic) = k

      iwater = 0
      if (constraints%water_type == TIP4) then

        do ix = nsolute(ic)+1, natom(ic)-3, 4

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1, 4
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do

      else if (constraints%water_type == TIP3) then

        do ix = nsolute(ic)+1, natom(ic)-2, 3

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1,3
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do
      else if (constraints%water_type == TIP2) then

        do ix = nsolute(ic)+1, natom(ic)-1, 2

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1,2
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do

      else if (constraints%water_type == TIP1) then

        do ix = nsolute(ic)+1, natom(ic)

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if
          trans(1,ix,ic) = move(1) - boundary%origin_x
          trans(2,ix,ic) = move(2) - boundary%origin_y
          trans(3,ix,ic) = move(3) - boundary%origin_z
        end do

      end if

      nwater(ic) = iwater

    end do
    !$omp end parallel

    return

  end subroutine update_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    updae_cell_boundary_enm
  !> @brief        Update coordinates of the boundary cell
  !! @authors      JJ
  !! @param[inout] domain   : domain information
  !! @param[inout] comm     : communication information
  !! @param[inout] boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_boundary_enm(domain, enefunc, comm)

    ! formal arguments
    type(s_domain),   target, intent(in   ) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_comm),     target, intent(inout) :: comm

    ! local variable
    integer                   :: i, k, ii, ic, ix
    integer                   :: irequest, irequest1, irequest2, irequest3
    integer                   :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wp),         pointer :: force(:,:), dist(:,:)
    real(wip),        pointer :: buf_send(:,:), buf_recv(:,:)
    integer,          pointer :: enm(:), list(:,:,:), func(:,:)
    integer,          pointer :: num_cell(:)
    integer,          pointer :: int_send(:,:), int_recv(:,:)
    integer,          pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,          pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,          pointer :: upper_recv(:), lower_recv(:)
    integer,          pointer :: upper_send(:), lower_send(:)
    integer,          pointer :: iproc_upper(:), iproc_lower(:)


    enm                 => enefunc%num_enm
    list                => enefunc%enm_list
    force               => enefunc%enm_force_const
    dist                => enefunc%enm_dist_min

    iproc_upper         => domain%iproc_upper
    iproc_lower         => domain%iproc_lower

    num_cell            => comm%num_cell
    int_send            => comm%int_send
    int_recv            => comm%int_recv
    buf_send            => comm%buf_send
    buf_recv            => comm%buf_recv
    ic_lower_send       => comm%ic_lower_send
    ic_upper_send       => comm%ic_upper_send
    ic_lower_recv       => comm%ic_lower_recv
    ic_upper_recv       => comm%ic_upper_recv
    upper_recv          => comm%recv_enm_upper_size
    lower_recv          => comm%recv_enm_lower_size
    upper_send          => comm%send_enm_upper_size
    lower_send          => comm%send_enm_lower_size

    do ii = 1, 3

      ! Pack the data (lower)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_lower_send(i,ii)
        do ix = 1, enm(ic)
          k = k + 1
          buf_send(2*k-1,1) = force(ix,ic)
          buf_send(2*k  ,1) = dist (ix,ic)
          int_send(2*k-1,1) = list(1,ix,ic)
          int_send(2*k  ,1) = list(2,ix,ic)
        end do
      end do

      ! pack the data (upper)
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_upper_send(i,ii)
        do ix = 1, enm(ic)
          k = k + 1
          buf_send(2*k-1,2) = force(ix,ic)
          buf_send(2*k  ,2) = dist (ix,ic)
          int_send(2*k-1,2) = list(1,ix,ic)
          int_send(2*k  ,2) = list(2,ix,ic)
        end do
      end do

#ifdef HAVE_MPI_GENESIS
      ! send the data
      !
      call mpi_irecv(int_recv(1,1), 2*upper_recv(ii), mpi_integer,  &
                     iproc_upper(ii), iproc_upper(ii),              &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), 2*lower_send(ii), mpi_integer,  &
                     iproc_lower(ii), my_city_rank,                 &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(int_recv(1,2), 2*lower_recv(ii), mpi_integer,  &
                     iproc_lower(ii), my_city_rank+nproc_city,      &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), 2*upper_send(ii), mpi_integer,  &
                     iproc_upper(ii), iproc_upper(ii)+nproc_city,   &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      call mpi_irecv(buf_recv(1,1), 2*upper_recv(ii), mpi_wip_real, &
                     iproc_upper(ii), iproc_upper(ii),              &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), 2*lower_send(ii),  mpi_wip_real,&
                     iproc_lower(ii), my_city_rank,                 &
                     mpi_comm_city, irequest5, ierror)
      call mpi_irecv(buf_recv(1,2), 2*lower_recv(ii), mpi_wip_real, &
                     iproc_lower(ii), my_city_rank+nproc_city,      &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), 2*upper_send(ii), mpi_wip_real, &
                     iproc_upper(ii), iproc_upper(ii)+nproc_city,   &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the data of the boundary cells
      !
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_upper_recv(i,ii)
        do ix = 1, enm(ic)
          k = k + 1
          force(  ix,ic) = buf_recv(2*k-1,1)
          dist (  ix,ic) = buf_recv(2*k  ,1)
          list (1,ix,ic) = int_recv(2*k-1,1)
          list (2,ix,ic) = int_recv(2*k  ,1)
        end do
      end do

      k = 0
      do i = 1, num_cell(ii)
        ic = ic_lower_recv(i,ii)
        do ix = 1, enm(ic)
          k = k + 1
          force(  ix,ic) = buf_recv(2*k-1,2)
          dist (  ix,ic) = buf_recv(2*k  ,2)
          list (1,ix,ic) = int_recv(2*k-1,2)
          list (2,ix,ic) = int_recv(2*k  ,2)
        end do
      end do
    end do

    return

  end subroutine update_cell_boundary_enm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_cell_size_enm
  !> @brief        Update cell size of enm function
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !! @param[inout] comm    : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_size_enm(domain, enefunc, comm)

    ! formal arguments
    type(s_domain),  target, intent(in   ) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, k, ii, ic, size, ncell, ncell_local
    integer                  :: irequest, irequest1, irequest2, irequest3
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    integer,         pointer :: enm(:)
    integer,         pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,         pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)
    integer,         pointer :: upper_send(:), lower_send(:)
    integer,         pointer :: upper_recv(:), lower_recv(:)
    integer,         pointer :: num_cell(:)
    integer,     allocatable :: send_size(:,:), recv_size(:,:)


    enm           => enefunc%num_enm
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower

    upper_send    => comm%send_enm_upper_size
    lower_send    => comm%send_enm_lower_size
    upper_recv    => comm%recv_enm_upper_size
    lower_recv    => comm%recv_enm_lower_size

    num_cell      => comm%num_cell
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv


    size = max(num_cell(1), num_cell(2), num_cell(3))
    allocate(recv_size(size,2))
    allocate(send_size(size,2))

    do ii = 1, 3

      ! check the send size of coordinates
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_lower_send(i,ii)
        send_size(i,1) = enm(ic)
        k = k + enm(ic)
      end do
      lower_send(ii) = k

      ! check the send size of coordinates
      k = 0
      do i = 1, num_cell(ii)
        ic = ic_upper_send(i,ii)
        send_size(i,2) = enm(ic)
        k = k + enm(ic)
      end do
      upper_send(ii) = k

#ifdef HAVE_MPI_GENESIS
      ! send the size of the data
      !
      call mpi_irecv(recv_size(1,1), num_cell(ii), mpi_integer,   &
                     iproc_upper(ii), iproc_upper(ii),            &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(send_size(1,1), num_cell(ii), mpi_integer,   &
                     iproc_lower(ii), my_city_rank,               &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(recv_size(1,2), num_cell(ii), mpi_integer,   &
                     iproc_lower(ii), my_city_rank+nproc_city,    &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(1,2), num_cell(ii), mpi_integer,   &
                     iproc_upper(ii), iproc_upper(ii)+nproc_city, &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)
#endif

      k = 0
      do i = 1, num_cell(ii)
        ic = ic_upper_recv(i,ii)
        enm(ic) = recv_size(i,1)
        k = k + enm(ic)
      end do
      upper_recv(ii) = k

      k = 0
      do i = 1, num_cell(ii)
        ic = ic_lower_recv(i,ii)
        enm(ic) = recv_size(i,2)
        k = k + enm(ic)
      end do
      lower_recv(ii) = k

    end do

    deallocate(recv_size)
    deallocate(send_size)

    return

  end subroutine update_cell_size_enm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_cell_size_constraints
  !> @brief        Update cell size with constraints
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_size_constraints(domain, comm, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_comm),        target, intent(inout) :: comm
    type(s_constraints), target, intent(inout) :: constraints

    ! local variable
    integer                  :: i, j, k, ii, ic, size, size1
    integer                  :: irequest, irequest1, irequest2, irequest3
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    integer,         pointer :: num_cell(:), natom(:)
    integer,         pointer :: No_HGr(:), HGr_local(:,:)
    integer,         pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,         pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)
    integer,         pointer :: upper_send(:), lower_send(:)
    integer,         pointer :: upper_recv(:), lower_recv(:)
    integer,         pointer :: upper_fsend(:), lower_fsend(:)
    integer,         pointer :: upper_frecv(:), lower_frecv(:)
    integer,     allocatable :: send_size(:,:), recv_size(:,:)


    natom         => domain%num_atom
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower

    num_cell      => comm%num_cell
    if_lower_send => comm%if_lower_send
    if_upper_send => comm%if_upper_send
    if_lower_recv => comm%if_lower_recv
    if_upper_recv => comm%if_upper_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv
    upper_send    => comm%send_coord_upper_size
    lower_send    => comm%send_coord_lower_size
    upper_recv    => comm%recv_coord_upper_size
    lower_recv    => comm%recv_coord_lower_size
    upper_fsend   => comm%send_force_upper_size
    lower_fsend   => comm%send_force_lower_size
    upper_frecv   => comm%recv_force_upper_size
    lower_frecv   => comm%recv_force_lower_size

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local


    size = max(num_cell(1), num_cell(2), num_cell(3), 1)
    size1 = constraints%connect + 2
    allocate(recv_size(size1*size,2))
    allocate(send_size(size1*size,2))

    do ii = 1, 3

      ! check the send size of coordinates
      k = 0
      do i = 1, num_cell(ii)

        ic = ic_lower_send(i,ii)
        send_size(size1*(i-1)+1,1) = natom(ic)
        send_size(size1*(i-1)+2,1) = No_HGr(ic)

        do j = 1, constraints%connect
          send_size(size1*(i-1)+j+2,1) = HGr_local(j,ic)
        end do

        k = k + natom(ic)

      end do
      lower_send(ii) = k

      k = 0
      do i = 1, comm%num_cell(ii)

        ic = ic_upper_send(i,ii)
        send_size(size1*(i-1)+1,2) = natom(ic)
        send_size(size1*(i-1)+2,2) = No_HGr(ic)

        do j = 1, constraints%connect
          send_size(size1*(i-1)+j+2,2) = HGr_local(j,ic)
        end do

        k = k + natom(ic)

      end do
      upper_send(ii) = k

#ifdef HAVE_MPI_GENESIS

      ! send the size of the data
      !
      call mpi_irecv(recv_size(1,1), size1*num_cell(ii), mpi_integer,  &
                     iproc_upper(ii),                                  &
                     iproc_upper(ii),                                  &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(send_size(1,1), size1*num_cell(ii), mpi_integer,  &
                     iproc_lower(ii),                                  &
                     my_city_rank,                                     &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(recv_size(1,2), size1*num_cell(ii), mpi_integer,  &
                     iproc_lower(ii),                                  &
                     my_city_rank+nproc_city,                          &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(1,2), size1*num_cell(ii), mpi_integer,  &
                     iproc_upper(ii),                                  &
                     iproc_upper(ii)+nproc_city,                       &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)
#endif

      k = 0
      do i = 1, num_cell(ii)

        ic = ic_upper_recv(i,ii)
        natom(ic)  = recv_size(size1*(i-1)+1,1)
        No_HGr(ic) = recv_size(size1*(i-1)+2,1)

        do j = 1, constraints%connect
          HGr_local(j,ic) = recv_size(size1*(i-1)+j+2,1)
        end do

        k = k + natom(ic)

      end do
      upper_recv(ii) = k

      k = 0
      do i = 1, num_cell(ii)

        ic = ic_lower_recv(i,ii)
        natom(ic)  = recv_size(size1*(i-1)+1,2)
        No_HGr(ic) = recv_size(size1*(i-1)+2,2)

        do j = 1, constraints%connect
          HGr_local(j,ic) = recv_size(size1*(i-1)+j+2,2)
        end do

        k = k + natom(ic)

      end do
      lower_recv(ii) = k

    end do

    do ii = 3, 1, -1

      ! check the send size of forces
      k = 0
      do i = 1, num_cell(ii)
        ic = if_lower_send(i,ii)
        k = k + natom(ic)
      end do
      lower_fsend(ii) = k

      ! check the send size of forces
      k = 0
      do i = 1, num_cell(ii)
        ic = if_upper_send(i,ii)
        k = k + natom(ic)
      end do
      upper_fsend(ii) = k

#ifdef HAVE_MPI_GENESIS
      ! send the size of the data
      !
      call mpi_irecv(upper_frecv(ii), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii),                                  &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(lower_fsend(ii), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank,                                     &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(lower_frecv(ii), 1, mpi_integer, iproc_lower(ii), &
                     my_city_rank+nproc_city,                          &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(upper_fsend(ii), 1, mpi_integer, iproc_upper(ii), &
                     iproc_upper(ii)+nproc_city,                       &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)
#endif

      k = 0
      do i = 1, num_cell(ii)
        ic = if_upper_recv(i,ii)
        k = k + natom(ic)
      end do

      if (upper_frecv(ii) /= k) &
        call error_msg('Update_Cell_Size_Constraints> Disagreement between'//&
                       ' sending and receving processors')

      k = 0
      do i = 1, num_cell(ii)
        ic = if_lower_recv(i,ii)
        k = k + natom(ic)
      end do

      if (lower_frecv(ii) /= k) &
        call error_msg('Update_Cell_Size_Constraints> Disagreement between'//&
                       ' sending and receving processors')

    end do

    deallocate(recv_size)
    deallocate(send_size)

    return

  end subroutine update_cell_size_constraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_random_number
  !> @brief        Transfer random number information
  !! @authors      JJ
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_random_number(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, ii, ic
    integer                  :: irequest, irequest1, irequest2, irequest3
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: buf_send(:,:), buf_recv(:,:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)


    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower

    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv


    do ii = 1, 3

      ! Pack the data (lower)
      !
      do i = 1, comm%num_cell(ii)
        ic = ic_lower_send(i,ii)
        buf_send(i,1) = domain%random(ic)
      end do

      ! Pack the data (upper)
      !
      do i = 1, comm%num_cell(ii)
        ic = ic_upper_send(i,ii)
        buf_send(i,2) = domain%random(ic)
      end do

#ifdef HAVE_MPI_GENESIS

      ! send the data
      !
      call mpi_irecv(buf_recv(1,1), comm%num_cell(ii), mpi_wip_real,       &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii),                                      &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(buf_send(1,1), comm%num_cell(ii), mpi_wip_real,       &
                     iproc_lower(ii),                                      &
                     my_city_rank,                                         &
                     mpi_comm_city, irequest1, ierror)

      ! send the data(upper)
      !
      call mpi_irecv(buf_recv(1,2), comm%num_cell(ii), mpi_wip_real,       &
                     iproc_lower(ii),                                      &
                     my_city_rank+nproc_city,                              &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(buf_send(1,2), comm%num_cell(ii), mpi_wip_real,       &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii)+nproc_city,                           &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)
#endif

      ! Get the data (from upper)
      do i = 1, comm%num_cell(ii)
        ic = ic_upper_recv(i,ii)
        domain%random(ic) = buf_recv(i,1)
      end do

      ! Get the data (from lower)
      do i = 1, comm%num_cell(ii)
        ic = ic_lower_recv(i,ii)
        domain%random(ic) = buf_recv(i,2)
      end do

    end do

    return

  end subroutine update_random_number


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_constraints_fep
  !> @brief        Pack the incoming particles data of each boundary cell
  !                for FEP
  !! @authors      HO
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_constraints_fep(domain, comm, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_comm),        target, intent(inout) :: comm
    type(s_constraints), target, intent(inout) :: constraints

    ! local variable
    integer                      :: send_size(2,2), recv_size(2,2)
    integer                      :: real_size, int_size
    integer                      :: i, j, k, l, ii, ic, ix, nadd, water_atom
    integer                      :: list, list1, size
    integer                      :: irequest, irequest1, irequest2, irequest3
    integer                      :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                      :: istatus(mpi_status_size)
#endif
    real(wip),           pointer :: buf_send(:,:), buf_recv(:,:)
    real(wip),           pointer :: buf_real(:,:,:), water_move_real(:,:,:)
    real(wip),           pointer :: HGr_move_real(:,:,:,:)
    integer,             pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,             pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,             pointer :: iproc_upper(:), iproc_lower(:), num_cell(:)
    integer,             pointer :: ptl_add(:), water_move(:)
    integer,             pointer :: buf_int(:,:,:), water_move_int(:,:,:)
    integer,             pointer :: HGr_move(:,:), HGr_move_int(:,:,:,:)
    integer,             pointer :: int_send(:,:), int_recv(:,:)

    ! FEP
    integer                      :: num_buf_real, num_buf_int
    integer                      :: num_HGr_real, num_HGr_int


    iproc_upper     => domain%iproc_upper
    iproc_lower     => domain%iproc_lower
    ptl_add         => domain%ptl_add
    buf_real        => domain%buf_real
    buf_int         => domain%buf_integer
    water_move      => domain%water%move
    water_move_real => domain%water%move_real
    water_move_int  => domain%water%move_integer

    buf_send        => comm%buf_send
    buf_recv        => comm%buf_recv
    int_send        => comm%int_send
    int_recv        => comm%int_recv
    if_lower_send   => comm%if_lower_send
    if_upper_send   => comm%if_upper_send
    if_lower_recv   => comm%if_lower_recv
    if_upper_recv   => comm%if_upper_recv
    num_cell        => comm%num_cell

    HGr_move        => constraints%HGr_move
    HGr_move_int    => constraints%HGr_move_int
    HGr_move_real   => constraints%HGr_move_real

    if (constraints%water_type == TIP4) then
      water_atom = 4
    else if (constraints%water_type == TIP3) then
      water_atom = 3
    else if (constraints%water_type == TIP2) then
      water_atom = 2
    else if (constraints%water_type == TIP1) then
      water_atom = 1
    end if

    num_buf_real = 10
    num_buf_int  = 5
    num_HGr_real = 11
    num_HGr_int  = 5

    size = constraints%connect + 2
    do ii = 3, 1, -1

      ! Pack outgoing data (lower)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)

        int_send(size*(i-1)+1,1) = ptl_add(ic)
        int_send(size*(i-1)+2,1) = water_move(ic)

        do j = 1, constraints%connect
          int_send(size*(i-1)+j+2,1) = HGr_move(j,ic)
        end do

        do ix = 1, ptl_add(ic)

          k = k + 1
          do list = 1, num_buf_real
            buf_send(num_buf_real*(k-1)+list,1) = buf_real(list,ix,ic)
          end do

          do list = 1, num_buf_int
            list1 = num_buf_int*(k-1) + list + size*num_cell(ii)
            int_send(list1,1) = buf_int(list,ix,ic)
          end do

        end do

      end do

      real_size = num_buf_real*k
      int_size = num_buf_int*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)
        do j = 1, constraints%connect
          do ix = 1, HGr_move(j,ic)
            do l = 1, j+1

              k = k + 1
              do list = 1, num_HGr_real
                buf_send(real_size+num_HGr_real*(k-1)+list,1) = &
                     HGr_move_real(num_HGr_real*(l-1)+list,ix,j,ic)
              end do

              do list = 1, num_HGr_int
                list1 = num_HGr_int*(k-1) + list + int_size
                int_send(list1,1) = HGr_move_int(num_HGr_int*(l-1)+list,ix,j,ic)
              end do

            end do
          end do
        end do

      end do

      real_size = real_size + num_HGr_real*k
      int_size = int_size + num_HGr_int*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_send(i,ii)
        do ix = 1, water_move(ic)

          k = k + 1
          do j = 1, 6*water_atom
            buf_send(real_size+6*water_atom*(k-1)+j,1) &
                            = water_move_real(j,ix,ic)
          end do

          do j = 1, 2*water_atom
            int_send(int_size+2*water_atom*(k-1)+j,1)    &
                          = water_move_int(j,ix,ic)
          end do

        end do
      end do

      send_size(1,1) = int_size + 2*water_atom*k
      send_size(2,1) = real_size + 6*water_atom*k

      ! pack the outgoing date (upper)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)

        int_send(size*(i-1)+1,2) = ptl_add(ic)
        int_send(size*(i-1)+2,2) = water_move(ic)

        do j = 1, constraints%connect
          int_send(size*(i-1)+j+2,2) = HGr_move(j,ic)
        end do

        do ix = 1, ptl_add(ic)

          k = k + 1
          do list = 1, num_buf_real
            buf_send(num_buf_real*(k-1)+list,2) = buf_real(list,ix,ic)
          end do

          do list = 1, num_buf_int
            list1 = num_buf_int*(k-1) + list + size*num_cell(ii)
            int_send(list1,2) = buf_int(list,ix,ic)
          end do

        end do

      end do

      real_size = num_buf_real*k
      int_size = num_buf_int*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)
        do j = 1, constraints%connect
          do ix = 1, HGr_move(j,ic)

            do l = 1, j+1

              k = k + 1
              do list = 1, num_HGr_real
                buf_send(real_size+num_HGr_real*(k-1)+list,2) = &
                     HGr_move_real(num_HGr_real*(l-1)+list,ix,j,ic)
              end do

              do list = 1, num_HGr_int
                list1 = num_HGr_int*(k-1) + list + int_size
                int_send(list1,2) = HGr_move_int(num_HGr_int*(l-1)+list,ix,j,ic)
              end do

            end do
          end do
        end do

      end do

      real_size = real_size + num_HGr_real*k
      int_size  = int_size + num_HGr_int*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_send(i,ii)
        do ix = 1, water_move(ic)

          k = k + 1
          do j = 1, 6*water_atom
            buf_send(real_size+6*water_atom*(k-1)+j,2) &
                           =  water_move_real(j,ix,ic)
          end do

          do j = 1, 2*water_atom
            int_send(int_size+2*water_atom*(k-1)+j,2)   &
                          = water_move_int(j,ix,ic)
          end do

        end do
      end do

      send_size(1,2) = int_size + 2*water_atom*k
      send_size(2,2) = real_size + 6*water_atom*k
#ifdef HAVE_MPI_GENESIS

      ! send the size of data
      !
      call mpi_irecv(recv_size(1,1), 2, mpi_integer, iproc_upper(ii),     &
                     iproc_upper(ii),           &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(send_size(1,1), 2, mpi_integer, iproc_lower(ii),     &
                     my_city_rank,           &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(recv_size(1,2), 2, mpi_integer, iproc_lower(ii),     &
                     my_city_rank+nproc_city,           &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(send_size(1,2), 2, mpi_integer, iproc_upper(ii),     &
                     iproc_upper(ii)+nproc_city,           &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,1), recv_size(1,1), mpi_integer,          &
                     iproc_upper(ii),                                     &
                     iproc_upper(ii),           &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), send_size(1,1), mpi_integer,          &
                     iproc_lower(ii),                                     &
                     my_city_rank,           &
                     mpi_comm_city, irequest1, ierror)

      ! send the data
      !
      call mpi_irecv(int_recv(1,2), recv_size(1,2), mpi_integer,          &
                     iproc_lower(ii),                                     &
                     my_city_rank+nproc_city,           &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), send_size(1,2), mpi_integer,          &
                     iproc_upper(ii),                                     &
                     iproc_upper(ii)+nproc_city,           &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      call mpi_irecv(buf_recv(1,1), recv_size(2,1), mpi_wip_real,       &
                     iproc_upper(ii),                                   &
                     iproc_upper(ii)+2*nproc_city,         &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), send_size(2,1), mpi_wip_real,       &
                     iproc_lower(ii),                                   &
                     my_city_rank+2*nproc_city,         &
                     mpi_comm_city, irequest5, ierror)

      call mpi_irecv(buf_recv(1,2), recv_size(2,2), mpi_wip_real,       &
                     iproc_lower(ii),                                   &
                     my_city_rank+3*nproc_city,         &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), send_size(2,2), mpi_wip_real,       &
                     iproc_upper(ii),                                   &
                     iproc_upper(ii)+3*nproc_city,         &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the imcoming data (from upper)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        nadd = int_recv(size*(i-1)+1,1)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do list = 1, num_buf_real
            buf_real(list,ix,ic) = buf_recv(num_buf_real*(k-1)+list,1)
          end do

          do list = 1, num_buf_int
            list1 = num_buf_int*(k-1) + list + size*num_cell(ii)
            buf_int(list,ix,ic) = int_recv(list1,1)
          end do

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

      real_size = num_buf_real*k
      int_size = num_buf_int*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        do j = 1, constraints%connect

          nadd = int_recv(size*(i-1)+j+2,1)
          do ix = HGr_move(j,ic)+1, HGr_move(j,ic)+nadd

            do l = 1, j+1

              k = k + 1
              do list = 1, num_HGr_real 
                HGr_move_real(num_HGr_real*(l-1)+list,ix,j,ic) = &
                     buf_recv(real_size+num_HGr_real*(k-1)+list,1)
              end do

              do list = 1, num_HGr_int
                list1 = num_HGr_int*(k-1) + list + int_size
                HGr_move_int(num_HGr_int*(l-1)+list,ix,j,ic) = int_recv(list1,1)
              end do

            end do
          end do

          HGr_move(j,ic) = HGr_move(j,ic) + nadd

        end do
      end do

      real_size = real_size + num_HGr_real*k
      int_size  = int_size + num_HGr_int*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_upper_recv(i,ii)
        nadd = int_recv(size*(i-1)+2,1)

        do ix = water_move(ic)+1, water_move(ic)+nadd

          k = k + 1
          do j = 1, 6*water_atom
            water_move_real(j,ix,ic) &
              = buf_recv(real_size+6*water_atom*(k-1)+j,1)
          end do

          do j = 1, 2*water_atom
            water_move_int(j,ix,ic) &
            = int_recv(int_size+2*water_atom*(k-1)+j,1)
          end do

        end do

        water_move(ic) = water_move(ic) + nadd

      end do

      ! get the imcoming data (from lower)
      !
      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        nadd = int_recv(size*(i-1)+1,2)

        do ix = ptl_add(ic)+1, ptl_add(ic)+nadd

          k = k + 1
          do list = 1, num_buf_real
            buf_real(list,ix,ic) = buf_recv(num_buf_real*(k-1)+list,2)
          end do

          do list = 1, num_buf_int
            list1 = num_buf_int*(k-1) + list + size*num_cell(ii)
            buf_int(list,ix,ic) = int_recv(list1,2)
          end do

        end do

        ptl_add(ic) = ptl_add(ic) + nadd

      end do

      real_size = num_buf_real*k
      int_size = num_buf_int*k + size*num_cell(ii)

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        do j = 1, constraints%connect

          nadd = int_recv(size*(i-1)+j+2,2)
          do ix = HGr_move(j,ic)+1, HGr_move(j,ic)+nadd
            do l = 1, j+1

              k = k + 1
              do list = 1, num_HGr_real 
                HGr_move_real(num_HGr_real*(l-1)+list,ix,j,ic) = &
                     buf_recv(real_size+num_HGr_real*(k-1)+list,2)
              end do

              do list = 1, num_HGr_int
                list1 = num_HGr_int*(k-1) + list + int_size
                HGr_move_int(num_HGr_int*(l-1)+list,ix,j,ic) = int_recv(list1,2)
              end do

            end do
          end do

          HGr_move(j,ic) = HGr_move(j,ic) + nadd

        end do

      end do

      real_size = real_size + num_HGr_real*k
      int_size = int_size + num_HGr_int*k

      k = 0
      do i = 1, num_cell(ii)

        ic = if_lower_recv(i,ii)
        nadd = int_recv(size*(i-1)+2,2)

        do ix = water_move(ic)+1, water_move(ic)+nadd

          k = k + 1
          do j = 1, 6*water_atom
            water_move_real(j,ix,ic) &
              = buf_recv(real_size+6*water_atom*(k-1)+j,2)
          end do

          do j = 1, 2*water_atom
            water_move_int(j,ix,ic) &
               = int_recv(int_size+2*water_atom*(k-1)+j,2)
          end do

        end do

        water_move(ic) = water_move(ic) + nadd

      end do

    end do

    return

  end subroutine communicate_constraints_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_bond_fep
  !> @brief        Pack the incoming bonding data of each boundary cell
  !                for FEP
  !! @authors      HO
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_bond_fep(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: send_size(2,2),recv_size(2,2)
    integer                   :: i, j, k, ii, ic, ix, nadd
    integer                   :: irequest, irequest1, irequest2, irequest3
    integer                   :: irequest4, irequest5, irequest6, irequest7
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:,:), buf_recv(:,:)
    real(wp),         pointer :: buf_bond_real(:,:,:), buf_angle_real(:,:,:)
    real(wp),         pointer :: buf_dihed_real(:,:,:), buf_impr_real(:,:,:)
    real(wp),         pointer :: buf_rb_dihed_real(:,:,:)
    real(wp),         pointer :: buf_restraint_real(:,:,:)
    real(wp),         pointer :: buf_fitting_real(:,:,:)
    integer,          pointer :: if_lower_send(:,:), if_upper_send(:,:)
    integer,          pointer :: if_lower_recv(:,:), if_upper_recv(:,:)
    integer,          pointer :: iproc_upper(:), iproc_lower(:)
    integer,          pointer :: int_send(:,:), int_recv(:,:)
    integer,          pointer :: bond_add(:), buf_bond_integer(:,:,:)
    integer,          pointer :: angle_add(:), buf_angle_integer(:,:,:)
    integer,          pointer :: dihed_add(:), buf_dihed_integer(:,:,:)
    integer,          pointer :: rb_dihed_add(:), buf_rb_dihed_integer(:,:,:)
    integer,          pointer :: impr_add(:), buf_impr_integer(:,:,:)
    integer,          pointer :: cmap_add(:), buf_cmap_integer(:,:,:)
    integer,          pointer :: restraint_add(:), buf_restraint_integer(:,:)
    integer,          pointer :: fitting_add(:), buf_fitting_integer(:,:)


    iproc_upper           => domain%iproc_upper
    iproc_lower           => domain%iproc_lower

    buf_send              => comm%buf_send
    buf_recv              => comm%buf_recv
    if_lower_send         => comm%if_lower_send
    if_upper_send         => comm%if_upper_send
    if_lower_recv         => comm%if_lower_recv
    if_upper_recv         => comm%if_upper_recv
    int_send              => comm%int_send
    int_recv              => comm%int_recv

    bond_add              => enefunc%bond_add
    buf_bond_real         => enefunc%buf_bond_real
    buf_bond_integer      => enefunc%buf_bond_integer
    angle_add             => enefunc%angle_add
    buf_angle_real        => enefunc%buf_angle_real
    buf_angle_integer     => enefunc%buf_angle_integer
    dihed_add             => enefunc%dihed_add
    buf_dihed_real        => enefunc%buf_dihed_real
    buf_dihed_integer     => enefunc%buf_dihed_integer
    rb_dihed_add          => enefunc%rb_dihed_add
    buf_rb_dihed_real     => enefunc%buf_rb_dihed_real
    buf_rb_dihed_integer  => enefunc%buf_rb_dihed_integer
    impr_add              => enefunc%impr_add
    buf_impr_real         => enefunc%buf_impr_real
    buf_impr_integer      => enefunc%buf_impr_integer
    cmap_add              => enefunc%cmap_add
    buf_cmap_integer      => enefunc%buf_cmap_integer
    restraint_add         => enefunc%restraint_add
    buf_restraint_real    => enefunc%buf_restraint_real
    buf_restraint_integer => enefunc%buf_restraint_integer
    fitting_add           => enefunc%fitting_add
    buf_fitting_real      => enefunc%buf_fitting_real
    buf_fitting_integer   => enefunc%buf_fitting_integer


    do ii = 3, 1, -1

      ! Pack outgoing data (lower)
      !

      k = 8*comm%num_cell(ii)
      j = 0

      do i = 1, comm%num_cell(ii)

        ic = if_lower_send(i,ii)
        int_send(8*i-7,1) = bond_add(ic)
        int_send(8*i-6,1) = angle_add(ic)
        int_send(8*i-5,1) = dihed_add(ic)
        int_send(8*i-4,1) = rb_dihed_add(ic)
        int_send(8*i-3,1) = impr_add(ic)
        int_send(8*i-2,1) = cmap_add(ic)
        int_send(8*i-1,1) = restraint_add(ic)
        int_send(8*i  ,1) = fitting_add(ic)

        do ix = 1, bond_add(ic)
          int_send(k+1:k+5,1) = buf_bond_integer(1:5,ix,ic)
          k = k + 5
          buf_send(j+1:j+2,1) = buf_bond_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, angle_add(ic)
          int_send(k+1:k+8,1) = buf_angle_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+4,1) = buf_angle_real(1:4,ix,ic)
          j = j + 4
        end do

        do ix = 1, dihed_add(ic)
          int_send(k+1:k+10,1) = buf_dihed_integer(1:10,ix,ic)
          k = k + 10
          buf_send(j+1:j+2,1) = buf_dihed_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, rb_dihed_add(ic)
          int_send(k+1:k+8,1) = buf_rb_dihed_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+6,1) = buf_rb_dihed_real(1:6,ix,ic)
          j = j + 6
        end do

        do ix = 1, impr_add(ic)
          int_send(k+1:k+9,1) = buf_impr_integer(1:9,ix,ic)
          k = k + 9
          buf_send(j+1:j+2,1) = buf_impr_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, cmap_add(ic)
          int_send(k+1:k+16,1) = buf_cmap_integer(1:16,ix,ic)
          k = k + 16
        end do

        do ix = 1, restraint_add(ic)
          int_send(k+1,1) = buf_restraint_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+7,1) = buf_restraint_real(1:7,ix,ic)
          j = j + 7
        end do

        do ix = 1, fitting_add(ic)
          int_send(k+1,1) = buf_fitting_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+3,1) = buf_fitting_real(1:3,ix,ic)
          j = j + 3
        end do

      end do

      send_size(1,1) = k
      send_size(2,1) = j

      ! pack the outgoing date (upper)
      !
      k = 8*comm%num_cell(ii)
      j = 0

      do i = 1, comm%num_cell(ii)

        ic = if_upper_send(i,ii)
        int_send(8*i-7,2) = bond_add(ic)
        int_send(8*i-6,2) = angle_add(ic)
        int_send(8*i-5,2) = dihed_add(ic)
        int_send(8*i-4,2) = rb_dihed_add(ic)
        int_send(8*i-3,2) = impr_add(ic)
        int_send(8*i-2,2) = cmap_add(ic)
        int_send(8*i-1,2) = restraint_add(ic)
        int_send(8*i  ,2) = fitting_add(ic)

        do ix = 1, bond_add(ic)
          int_send(k+1:k+5,2) = buf_bond_integer(1:5,ix,ic)
          k = k + 5
          buf_send(j+1:j+2,2) = buf_bond_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, angle_add(ic)
          int_send(k+1:k+8,2) = buf_angle_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+4,2) = buf_angle_real(1:4,ix,ic)
          j = j + 4
        end do

        do ix = 1, dihed_add(ic)
          int_send(k+1:k+10,2) = buf_dihed_integer(1:10,ix,ic)
          k = k + 10
          buf_send(j+1:j+2,2) = buf_dihed_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, rb_dihed_add(ic)
          int_send(k+1:k+8,2) = buf_rb_dihed_integer(1:8,ix,ic)
          k = k + 8
          buf_send(j+1:j+6,2) = buf_rb_dihed_real(1:6,ix,ic)
          j = j + 6
        end do

        do ix = 1, impr_add(ic)
          int_send(k+1:k+9,2) = buf_impr_integer(1:9,ix,ic)
          k = k + 9
          buf_send(j+1:j+2,2) = buf_impr_real(1:2,ix,ic)
          j = j + 2
        end do

        do ix = 1, cmap_add(ic)
          int_send(k+1:k+16,2) = buf_cmap_integer(1:16,ix,ic)
          k = k + 16
        end do

        do ix = 1, restraint_add(ic)
          int_send(k+1,2) = buf_restraint_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+7,2) = buf_restraint_real(1:7,ix,ic)
          j = j + 7
        end do

        do ix = 1, fitting_add(ic)
          int_send(k+1,2) = buf_fitting_integer(ix,ic)
          k = k + 1
          buf_send(j+1:j+3,2) = buf_fitting_real(1:3,ix,ic)
          j = j + 3
        end do

      end do

      send_size(1,2) = k
      send_size(2,2) = j

#ifdef HAVE_MPI_GENESIS

       ! send the size of data
       !
       call mpi_irecv(recv_size(1,1), 2, mpi_integer, iproc_upper(ii),     &
                      iproc_upper(ii),           &
                      mpi_comm_city, irequest, ierror)
       call mpi_isend(send_size(1,1), 2, mpi_integer, iproc_lower(ii),     &
                      my_city_rank,           &
                      mpi_comm_city, irequest1, ierror)
       call mpi_irecv(recv_size(1,2), 2, mpi_integer, iproc_lower(ii),     &
                      my_city_rank+nproc_city,           &
                      mpi_comm_city, irequest2, ierror)
       call mpi_isend(send_size(1,2), 2, mpi_integer, iproc_upper(ii),     &
                      iproc_upper(ii)+nproc_city,           &
                      mpi_comm_city, irequest3, ierror)

       call mpi_wait(irequest,  istatus, ierror)
       call mpi_wait(irequest1, istatus, ierror)
       call mpi_wait(irequest2, istatus, ierror)
       call mpi_wait(irequest3, istatus, ierror)

       ! send the integer data (lower)
       !
       call mpi_irecv(int_recv(1,1), recv_size(1,1), mpi_integer,          &
                      iproc_upper(ii),                                     &
                      iproc_upper(ii),           &
                      mpi_comm_city, irequest, ierror)
       call mpi_isend(int_send(1,1), send_size(1,1), mpi_integer,          &
                      iproc_lower(ii),                                     &
                      my_city_rank,           &
                      mpi_comm_city, irequest1, ierror)

       ! send the integer data (upper)
       !
       call mpi_irecv(int_recv(1,2), recv_size(1,2), mpi_integer,          &
                      iproc_lower(ii),                                     &
                      my_city_rank+nproc_city,           &
                      mpi_comm_city, irequest2, ierror)
       call mpi_isend(int_send(1,2), send_size(1,2), mpi_integer,          &
                      iproc_upper(ii),                                     &
                      iproc_upper(ii)+nproc_city,           &
                      mpi_comm_city, irequest3, ierror)

       call mpi_wait(irequest,  istatus, ierror)
       call mpi_wait(irequest1, istatus, ierror)
       call mpi_wait(irequest2, istatus, ierror)
       call mpi_wait(irequest3, istatus, ierror)

       ! send the real data (lower)
       !
       call mpi_irecv(buf_recv(1,1), recv_size(2,1), mpi_wip_real,      &
                      iproc_upper(ii),                                  &
                      iproc_upper(ii),        &
                      mpi_comm_city, irequest4, ierror)
       call mpi_isend(buf_send(1,1), send_size(2,1), mpi_wip_real,      &
                      iproc_lower(ii),                                  &
                      my_city_rank,        &
                      mpi_comm_city, irequest5, ierror)

       ! send the real data
       !
       call mpi_irecv(buf_recv(1,2), recv_size(2,2), mpi_wip_real,       &
                      iproc_lower(ii),                                   &
                      my_city_rank+nproc_city,         &
                      mpi_comm_city, irequest6, ierror)
       call mpi_isend(buf_send(1,2), send_size(2,2), mpi_wip_real,       &
                      iproc_upper(ii),                                   &
                      iproc_upper(ii)+nproc_city,         &
                      mpi_comm_city, irequest7, ierror)

       call mpi_wait(irequest4, istatus, ierror)
       call mpi_wait(irequest5, istatus, ierror)
       call mpi_wait(irequest6, istatus, ierror)
       call mpi_wait(irequest7, istatus, ierror)
#endif

       ! get the imcoming data
       !

       k = 8*comm%num_cell(ii)
       j = 0

       do i = 1, comm%num_cell(ii)
          ic = if_upper_recv(i,ii)

          nadd = int_recv(8*i-7,1)
          do ix = bond_add(ic)+1, bond_add(ic)+nadd
             buf_bond_integer(1:5,ix,ic) = int_recv(k+1:k+5,1)
             k = k + 5
             buf_bond_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          bond_add(ic) = bond_add(ic) + nadd

          nadd = int_recv(8*i-6,1)
          do ix = angle_add(ic)+1, angle_add(ic)+nadd
             buf_angle_integer(1:8,ix,ic) = int_recv(k+1:k+8,1)
             k = k + 8
             buf_angle_real(1:4,ix,ic) = buf_recv(j+1:j+4,1)
             j = j + 4
          end do
          angle_add(ic) = angle_add(ic) + nadd

          nadd = int_recv(8*i-5,1)
          do ix = dihed_add(ic)+1, dihed_add(ic)+nadd
             buf_dihed_integer(1:10,ix,ic) = int_recv(k+1:k+10,1)
             k = k + 10
             buf_dihed_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          dihed_add(ic) = dihed_add(ic) + nadd

          nadd = int_recv(8*i-4,1)
          do ix = rb_dihed_add(ic)+1, rb_dihed_add(ic)+nadd
             buf_rb_dihed_integer(1:8,ix,ic) = int_recv(k+1:k+8,1)
             k = k + 8
             buf_rb_dihed_real(1:6,ix,ic) = buf_recv(j+1:j+6,1)
             j = j + 6
          end do
          rb_dihed_add(ic) = rb_dihed_add(ic) + nadd

          nadd = int_recv(8*i-3,1)
          do ix = impr_add(ic)+1, impr_add(ic)+nadd
             buf_impr_integer(1:9,ix,ic) = int_recv(k+1:k+9,1)
             k = k + 9
             buf_impr_real(1:2,ix,ic) = buf_recv(j+1:j+2,1)
             j = j + 2
          end do
          impr_add(ic) = impr_add(ic) + nadd

          nadd = int_recv(8*i-2,1)
          do ix = cmap_add(ic)+1, cmap_add(ic)+nadd
             buf_cmap_integer(1:16,ix,ic) = int_recv(k+1:k+16,1)
             k = k + 16
          end do
          cmap_add(ic) = cmap_add(ic) + nadd

          nadd = int_recv(8*i-1,1)
          do ix = restraint_add(ic)+1, restraint_add(ic)+nadd
             buf_restraint_integer(ix,ic) = int_recv(k+1,1)
             k = k + 1
             buf_restraint_real(1:7,ix,ic) = buf_recv(j+1:j+7,1)
             j = j + 7
          end do
          restraint_add(ic) = restraint_add(ic) + nadd

          nadd = int_recv(8*i,1)
          do ix = fitting_add(ic)+1, fitting_add(ic)+nadd
             buf_fitting_integer(ix,ic) = int_recv(k+1,1)
             k = k + 1
             buf_fitting_real(1:3,ix,ic) = buf_recv(j+1:j+3,1)
             j = j + 3
          end do
          fitting_add(ic) = fitting_add(ic) + nadd

       end do

       ! get the incoming data
       !

       k = 8*comm%num_cell(ii)
       j = 0

       do i = 1, comm%num_cell(ii)

          ic = if_lower_recv(i,ii)

          nadd = int_recv(8*i-7,2)
          do ix = bond_add(ic)+1, bond_add(ic)+nadd
             buf_bond_integer(1:5,ix,ic) = int_recv(k+1:k+5,2)
             k = k + 5
             buf_bond_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          bond_add(ic) = bond_add(ic) + nadd

          nadd = int_recv(8*i-6,2)
          do ix = angle_add(ic)+1, angle_add(ic)+nadd
             buf_angle_integer(1:8,ix,ic) = int_recv(k+1:k+8,2)
             k = k + 8
             buf_angle_real(1:4,ix,ic) = buf_recv(j+1:j+4,2)
             j = j + 4
          end do
          angle_add(ic) = angle_add(ic) + nadd

          nadd = int_recv(8*i-5,2)
          do ix = dihed_add(ic)+1, dihed_add(ic)+nadd
             buf_dihed_integer(1:10,ix,ic) = int_recv(k+1:k+10,2)
             k = k + 10
             buf_dihed_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          dihed_add(ic) = dihed_add(ic) + nadd

          nadd = int_recv(8*i-4,2)
          do ix = rb_dihed_add(ic)+1, rb_dihed_add(ic)+nadd
             buf_rb_dihed_integer(1:8,ix,ic) = int_recv(k+1:k+8,2)
             k = k + 8
             buf_rb_dihed_real(1:6,ix,ic) = buf_recv(j+1:j+6,2)
             j = j + 6
          end do
          rb_dihed_add(ic) = rb_dihed_add(ic) + nadd

          nadd = int_recv(8*i-3,2)
          do ix = impr_add(ic)+1, impr_add(ic)+nadd
             buf_impr_integer(1:9,ix,ic) = int_recv(k+1:k+9,2)
             k = k + 9
             buf_impr_real(1:2,ix,ic) = buf_recv(j+1:j+2,2)
             j = j + 2
          end do
          impr_add(ic) = impr_add(ic) + nadd

          nadd = int_recv(8*i-2,2)
          do ix = cmap_add(ic)+1, cmap_add(ic)+nadd
             buf_cmap_integer(1:16,ix,ic) = int_recv(k+1:k+16,2)
             k = k + 16
          end do
          cmap_add(ic) = cmap_add(ic) + nadd

          nadd = int_recv(8*i-1,2)
          do ix = restraint_add(ic)+1, restraint_add(ic)+nadd
             k = k + 1
             buf_restraint_integer(ix,ic) = int_recv(k,2)
             buf_restraint_real(1:7,ix,ic) = buf_recv(j+1:j+7,2)
             j = j + 7
          end do
          restraint_add(ic) = restraint_add(ic) + nadd

          nadd = int_recv(8*i,2)
          do ix = fitting_add(ic)+1, fitting_add(ic)+nadd
             k = k + 1
             buf_fitting_integer(ix,ic) = int_recv(k,2)
             buf_fitting_real(1:3,ix,ic) = buf_recv(j+1:j+3,2)
             j = j + 3
          end do
          fitting_add(ic) = fitting_add(ic) + nadd

       end do

    end do

    return

  end subroutine communicate_bond_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    updae_cell_boundary_fep
  !> @brief        Update coordinates of the boundary cell with constraint
  !                for FEP
  !! @authors      HO
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_boundary_fep(domain, comm, boundary, constraints)

    ! formal arguments
    type(s_domain),      target, intent(inout) :: domain
    type(s_comm),        target, intent(inout) :: comm
    type(s_boundary),    target, intent(inout) :: boundary
    type(s_constraints), target, intent(inout) :: constraints

    ! local variable
    real(wip)                :: x_shift, y_shift, z_shift
    real(wip)                :: move(3)
    integer                  :: i, k, j, ii, ic, ix, iwater, list
    integer                  :: irequest, irequest1, irequest2, irequest3
    integer                  :: irequest4, irequest5, irequest6, irequest7
    integer                  :: id, omp_get_thread_num
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),        pointer :: trans(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: cell_pbc_move(:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: buf_send(:,:), buf_recv(:,:)
    real(wip),       pointer :: bsize_x, bsize_y, bsize_z
    integer,         pointer :: num_cell(:)
    integer,         pointer :: ic_lower_send(:,:), ic_upper_send(:,:)
    integer,         pointer :: ic_lower_recv(:,:), ic_upper_recv(:,:)
    integer,         pointer :: iproc_upper(:), iproc_lower(:)
    integer,         pointer :: upper_recv(:), lower_recv(:)
    integer,         pointer :: upper_send(:), lower_send(:)
    integer,         pointer :: int_send(:,:), int_recv(:,:)
    integer,         pointer :: ncell_local, ncell_bd, natom(:), nsolute(:)
    integer,         pointer :: nwater(:), water_list(:,:,:), atmcls(:,:)
    integer,         pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: No_HGr(:), HGr_local(:,:),HGr_bond_list(:,:,:,:)

    ! FEP
    integer,         pointer :: fepgrp(:,:)
    integer,         pointer :: fep_atmcls_singleB(:,:)
    real(wp),        pointer :: fep_chargeA(:,:)
    real(wp),        pointer :: fep_chargeB(:,:)
    integer                  :: num_buf_real, num_buf_int

    ncell_local   => domain%num_cell_local
    ncell_bd      => domain%num_cell_boundary
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water
    water_list    => domain%water_list
    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    inv_mass      => domain%inv_mass
    atmcls        => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    id_l2g_sol    => domain%id_l2g_solute
    id_g2l        => domain%id_g2l
    iproc_upper   => domain%iproc_upper
    iproc_lower   => domain%iproc_lower
    trans         => domain%trans_vec
    cell_pbc_move => domain%cell_pbc_move

    num_cell      => comm%num_cell
    int_send      => comm%int_send
    int_recv      => comm%int_recv
    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    ic_lower_send => comm%ic_lower_send
    ic_upper_send => comm%ic_upper_send
    ic_lower_recv => comm%ic_lower_recv
    ic_upper_recv => comm%ic_upper_recv
    upper_recv    => comm%recv_coord_upper_size
    lower_recv    => comm%recv_coord_lower_size
    upper_send    => comm%send_coord_upper_size
    lower_send    => comm%send_coord_lower_size

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list

    ! FEP
    fepgrp        => domain%fepgrp
    fep_chargeA   => domain%fep_chargeA
    fep_chargeB   => domain%fep_chargeB
    fep_atmcls_singleB => domain%fep_atmcls_singleB

    num_buf_real = 10
    num_buf_int  = 5

    do ii = 1, 3

      !$omp parallel private(k, i, ic, ix, id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      ! Pack the data (lower)
      !
      if (id == 0) then
      k = 0
      do i = 1, num_cell(ii)

        ic = ic_lower_send(i,ii)
        do ix = 1, natom(ic)
          k = k + 1
          buf_send(num_buf_real*k-num_buf_real+ 1,1) = coord(1,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 2,1) = coord(2,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 3,1) = coord(3,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 4,1) = velocity(1,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 5,1) = velocity(2,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 6,1) = velocity(3,ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 7,1) = charge(ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 8,1) = mass(ix,ic)
          buf_send(num_buf_real*k-num_buf_real+ 9,1) = fep_chargeA(ix,ic)
          buf_send(num_buf_real*k-num_buf_real+10,1) = fep_chargeB(ix,ic)
          int_send(num_buf_int*k-num_buf_int+1,1) = atmcls(ix,ic)
          int_send(num_buf_int*k-num_buf_int+2,1) = id_l2g(ix,ic)
          int_send(num_buf_int*k-num_buf_int+3,1) = id_l2g_sol(ix,ic)
          int_send(num_buf_int*k-num_buf_int+4,1) = fepgrp(ix,ic)
          int_send(num_buf_int*k-num_buf_int+5,1) = fep_atmcls_singleB(ix,ic)
        end do

      end do
      end if

      ! pack the data (upper)
      !
      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, num_cell(ii)

          ic = ic_upper_send(i,ii)
          do ix = 1, natom(ic)
            k = k + 1
            buf_send(num_buf_real*k-num_buf_real+ 1,2) = coord(1,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 2,2) = coord(2,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 3,2) = coord(3,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 4,2) = velocity(1,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 5,2) = velocity(2,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 6,2) = velocity(3,ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 7,2) = charge(ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 8,2) = mass(ix,ic)
            buf_send(num_buf_real*k-num_buf_real+ 9,2) = fep_chargeA(ix,ic)
            buf_send(num_buf_real*k-num_buf_real+10,2) = fep_chargeB(ix,ic)
            int_send(num_buf_int*k-num_buf_int+1,2) = atmcls(ix,ic)
            int_send(num_buf_int*k-num_buf_int+2,2) = id_l2g(ix,ic)
            int_send(num_buf_int*k-num_buf_int+3,2) = id_l2g_sol(ix,ic)
            int_send(num_buf_int*k-num_buf_int+4,2) = fepgrp(ix,ic)
            int_send(num_buf_int*k-num_buf_int+5,2) = fep_atmcls_singleB(ix,ic)
          end do

        end do
      end if

      !$omp end parallel

#ifdef HAVE_MPI_GENESIS

      ! send the data
      !
      call mpi_irecv(int_recv(1,1), num_buf_int*upper_recv(ii), mpi_integer, &
                     iproc_upper(ii),                                        &
                     iproc_upper(ii),              &
                     mpi_comm_city, irequest, ierror)
      call mpi_isend(int_send(1,1), num_buf_int*lower_send(ii), mpi_integer, &
                     iproc_lower(ii),                                        &
                     my_city_rank,              &
                     mpi_comm_city, irequest1, ierror)
      call mpi_irecv(int_recv(1,2), num_buf_int*lower_recv(ii), mpi_integer, &
                     iproc_lower(ii),                                        &
                     my_city_rank+nproc_city,              &
                     mpi_comm_city, irequest2, ierror)
      call mpi_isend(int_send(1,2), num_buf_int*upper_send(ii), mpi_integer, &
                     iproc_upper(ii),                                        &
                     iproc_upper(ii)+nproc_city,              &
                     mpi_comm_city, irequest3, ierror)

      call mpi_wait(irequest,  istatus, ierror)
      call mpi_wait(irequest1, istatus, ierror)
      call mpi_wait(irequest2, istatus, ierror)
      call mpi_wait(irequest3, istatus, ierror)

      call mpi_irecv(buf_recv(1,1), num_buf_real*upper_recv(ii), mpi_wip_real, &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii),            &
                     mpi_comm_city, irequest4, ierror)
      call mpi_isend(buf_send(1,1), num_buf_real*lower_send(ii), mpi_wip_real, &
                     iproc_lower(ii),                                      &
                     my_city_rank,            &
                     mpi_comm_city, irequest5, ierror)

      call mpi_irecv(buf_recv(1,2), num_buf_real*lower_recv(ii), mpi_wip_real, &
                     iproc_lower(ii),                                      &
                     my_city_rank+nproc_city,            &
                     mpi_comm_city, irequest6, ierror)
      call mpi_isend(buf_send(1,2), num_buf_real*upper_send(ii), mpi_wip_real, &
                     iproc_upper(ii),                                      &
                     iproc_upper(ii)+nproc_city,            &
                     mpi_comm_city, irequest7, ierror)

      call mpi_wait(irequest4, istatus, ierror)
      call mpi_wait(irequest5, istatus, ierror)
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      ! get the data of the boundary cells
      !

      !$omp parallel private(k, i, ic, ix, id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (id == 0) then
      k = 0
      do i = 1, num_cell(ii)

        ic = ic_upper_recv(i,ii)
        do ix = 1, natom(ic)
          k = k + 1
          coord(1,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 1,1)
          coord(2,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 2,1)
          coord(3,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 3,1)
          velocity(1,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 4,1)
          velocity(2,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 5,1)
          velocity(3,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 6,1)
          charge(ix,ic)      = buf_recv(num_buf_real*k-num_buf_real+ 7,1)
          mass(ix,ic)        = buf_recv(num_buf_real*k-num_buf_real+ 8,1)
          fep_chargeA(ix,ic) = buf_recv(num_buf_real*k-num_buf_real+ 9,1)
          fep_chargeB(ix,ic) = buf_recv(num_buf_real*k-num_buf_real+10,1)
          atmcls(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+1,1)
          id_l2g(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+2,1)
          id_l2g_sol(ix,ic)         = int_recv(num_buf_int*k-num_buf_int+3,1)
          fepgrp(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+4,1)
          fep_atmcls_singleB(ix,ic) = int_recv(num_buf_int*k-num_buf_int+5,1)
        end do

      end do
      end if

      if (mod(id+1,nthread) == 0) then
        k = 0
        do i = 1, num_cell(ii)

          ic = ic_lower_recv(i,ii)

          do ix = 1, natom(ic)
            k = k + 1
            coord(1,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 1,2)
            coord(2,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 2,2)
            coord(3,ix,ic)     = buf_recv(num_buf_real*k-num_buf_real+ 3,2)
            velocity(1,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 4,2)
            velocity(2,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 5,2)
            velocity(3,ix,ic)  = buf_recv(num_buf_real*k-num_buf_real+ 6,2)
            charge(ix,ic)      = buf_recv(num_buf_real*k-num_buf_real+ 7,2)
            mass(ix,ic)        = buf_recv(num_buf_real*k-num_buf_real+ 8,2)
            fep_chargeA(ix,ic) = buf_recv(num_buf_real*k-num_buf_real+ 9,2)
            fep_chargeB(ix,ic) = buf_recv(num_buf_real*k-num_buf_real+10,2)
            atmcls(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+1,2)
            id_l2g(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+2,2)
            id_l2g_sol(ix,ic)         = int_recv(num_buf_int*k-num_buf_int+3,2)
            fepgrp(ix,ic)             = int_recv(num_buf_int*k-num_buf_int+4,2)
            fep_atmcls_singleB(ix,ic) = int_recv(num_buf_int*k-num_buf_int+5,2)
          end do

        end do
      end if

      !$omp end parallel

    end do

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ix, x_shift, y_shift, z_shift, move, ic, iwater,&
    !$omp         list, ii)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! inverse of mass
    !
    do i = id+1, ncell_local, nthread
      do ix = 1, natom(i)
        inv_mass(ix,i) = 0.0_wip
        if (mass(ix,i) > EPS) inv_mass(ix,i) = 1.0_wip / mass(ix,i)
      end do
    end do

    ! check the translation for each particle
    !
    do i = id+1, ncell_local, nthread

!ocl nosimd
      do ix = 1, No_HGr(i)

        !coordinate shifted against the origin
        !
        x_shift = coord(1,ix,i) - boundary%origin_x
        y_shift = coord(2,ix,i) - boundary%origin_y
        z_shift = coord(3,ix,i) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)
        trans(1,ix,i) = move(1) - boundary%origin_x
        trans(2,ix,i) = move(2) - boundary%origin_y
        trans(3,ix,i) = move(3) - boundary%origin_z
        ii = id_l2g_sol(ix,i)
        id_g2l(1,ii) = i
        id_g2l(2,ii) = ix

      end do

      k = No_HGr(i)
      do j = 1, constraints%connect

!ocl nosimd
        do ix = 1, HGr_local(j,i)
          list = HGr_bond_list(1,ix,j,i)
          x_shift = coord(1,list,i) - boundary%origin_x
          y_shift = coord(2,list,i) - boundary%origin_y
          z_shift = coord(3,list,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do list = 1, j+1
            k = k + 1
            trans(1,k,i) = move(1) - boundary%origin_x
            trans(2,k,i) = move(2) - boundary%origin_y
            trans(3,k,i) = move(3) - boundary%origin_z
            ii = id_l2g_sol(k,i)
            id_g2l(1,ii) = i
            id_g2l(2,ii) = k
          end do

        end do
      end do
      nsolute(i) = k

      if (constraints%water_type == TIP4) then

        do ix = nsolute(i)+1, natom(i)-3, 4

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          do j = 1,4
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do

      else if (constraints%water_type == TIP3) then

        do ix = nsolute(i)+1, natom(i)-2, 3

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do j = 1,3
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do
      else if (constraints%water_type == TIP2) then

        do ix = nsolute(i)+1, natom(i)-1, 2

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          do j = 1,2
            trans(1,ix-1+j,i) = move(1) - boundary%origin_x
            trans(2,ix-1+j,i) = move(2) - boundary%origin_y
            trans(3,ix-1+j,i) = move(3) - boundary%origin_z
          end do

        end do

      else if (constraints%water_type == TIP1) then

        do ix = nsolute(i)+1, natom(i)

          x_shift = coord(1,ix,i) - boundary%origin_x
          y_shift = coord(2,ix,i) - boundary%origin_y
          z_shift = coord(3,ix,i) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          trans(1,ix,i) = move(1) - boundary%origin_x
          trans(2,ix,i) = move(2) - boundary%origin_y
          trans(3,ix,i) = move(3) - boundary%origin_z

        end do

      end if

    end do

    do i = id+1, ncell_bd, nthread
      ic = i + ncell_local

!ocl nosimd
      do ix = 1, No_HGr(ic)

        !coordinate shifted against the origin
        !
        x_shift = coord(1,ix,ic) - boundary%origin_x
        y_shift = coord(2,ix,ic) - boundary%origin_y
        z_shift = coord(3,ix,ic) - boundary%origin_z

        !coordinate shifted to the first quadrant and set into the boundary box
        !
        move(1) = bsize_x*0.5_wip-bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip-bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip-bsize_z*anint(z_shift/bsize_z)
        if (domain%nonbond_kernel == NBK_Fugaku .or. &
            domain%nonbond_kernel == NBK_Intel) then
          move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
          move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
          move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
        end if
        trans(1,ix,ic) = move(1) - boundary%origin_x
        trans(2,ix,ic) = move(2) - boundary%origin_y
        trans(3,ix,ic) = move(3) - boundary%origin_z
        ii = id_l2g_sol(ix,ic)
        id_g2l(1,ii) = ic
        id_g2l(2,ii) = ix
      end do

      k = No_HGr(ic)

      do j = 1, constraints%connect

!ocl nosimd
        do ix = 1, HGr_local(j,ic)

          k = k + 1
          HGr_bond_list(1,ix,j,ic) = k
          x_shift = coord(1,k,ic) - boundary%origin_x
          y_shift = coord(2,k,ic) - boundary%origin_y
          z_shift = coord(3,k,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if
          trans(1,k,ic) = move(1) - boundary%origin_x
          trans(2,k,ic) = move(2) - boundary%origin_y
          trans(3,k,ic) = move(3) - boundary%origin_z
          ii = id_l2g_sol(k,ic)
          id_g2l(1,ii) = ic
          id_g2l(2,ii) = k

          do list = 2, j+1
            k = k + 1
            HGr_bond_list(list,ix,j,ic) = k
            trans(1,k,ic) = move(1) - boundary%origin_x
            trans(2,k,ic) = move(2) - boundary%origin_y
            trans(3,k,ic) = move(3) - boundary%origin_z
            ii = id_l2g_sol(k,ic)
            id_g2l(1,ii) = ic
            id_g2l(2,ii) = k
          end do

        end do
      end do
      nsolute(ic) = k

      iwater = 0
      if (constraints%water_type == TIP4) then

        do ix = nsolute(ic)+1, natom(ic)-3, 4

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1, 4
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do

      else if (constraints%water_type == TIP3) then

        do ix = nsolute(ic)+1, natom(ic)-2, 3

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1,3
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do
      else if (constraints%water_type == TIP2) then

        do ix = nsolute(ic)+1, natom(ic)-1, 2

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if

          do j = 1,2
            trans(1,ix-1+j,ic) = move(1) - boundary%origin_x
            trans(2,ix-1+j,ic) = move(2) - boundary%origin_y
            trans(3,ix-1+j,ic) = move(3) - boundary%origin_z
            water_list(j,iwater,ic) = ix-1+j
          end do
        end do

      else if (constraints%water_type == TIP1) then

        do ix = nsolute(ic)+1, natom(ic)

          iwater = iwater + 1
          x_shift = coord(1,ix,ic) - boundary%origin_x
          y_shift = coord(2,ix,ic) - boundary%origin_y
          z_shift = coord(3,ix,ic) - boundary%origin_z
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
          end if
          trans(1,ix,ic) = move(1) - boundary%origin_x
          trans(2,ix,ic) = move(2) - boundary%origin_y
          trans(3,ix,ic) = move(3) - boundary%origin_z
        end do

      end if

      nwater(ic) = iwater

    end do
    !$omp end parallel

    return

  end subroutine update_cell_boundary_fep

end module sp_communicate_mod
