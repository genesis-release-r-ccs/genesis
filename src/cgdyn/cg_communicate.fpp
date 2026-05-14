!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   cg_communicate_mod
!> @brief   Provides utilities for MPI communication.                 
!! @authors Jaewoon Jung (JJ) 
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_communicate_mod

  use cg_boundary_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use cg_communicate_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only:  omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_communicate
  public  :: update_communicate_size
  public  :: communicate_coor
  public  :: communicate_force
  public  :: communicate_ptl
  public  :: communicate_ptl_martini
  public  :: update_cell_size
  public  :: update_cell_boundary
  public  :: communicate_bond
  public  :: communicate_angl
  public  :: communicate_dihe
  public  :: communicate_stack
  public  :: communicate_contact
  public  :: communicate_pwmcos
  public  :: communicate_pwmcosns
  public  :: communicate_restraint

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate
  !> @brief        Initializes communication among processors.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !! @param[out]   comm     : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate(boundary, domain, comm)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    call setup_communicate_alloc(boundary, domain, comm)
    call setup_communicate_coord_cell(boundary, domain, comm)
    call setup_communicate_force_cell(boundary, domain, comm)

    return

  end subroutine setup_communicate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate_alloc
  !> @brief        Determines the allocation size.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !! @param[out]   comm     : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate_alloc(boundary, domain, comm)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm  

    ! local variable
    integer                  :: cell(3), max_cell
    integer                  :: k
    integer                  :: icell_start(3), icell_end(3)
    integer                  :: jcell_start(3), jcell_end(3)
    integer                  :: size1, size2, ip
    integer,     allocatable :: rank_g2l(:)
    integer,         pointer :: ncell_local, ncell_boundary
    integer,         pointer :: cell_start(:), cell_end(:), cell_to_rank(:,:,:)
    integer,         pointer :: num_proc
    integer,         pointer :: iproc(:)
    integer,         pointer :: cell_start_proc(:,:), cell_end_proc(:,:)
    integer,         pointer :: cell_g2l(:), cell_g2b(:), num_domain(:)

    num_domain  => boundary%num_domain

    ncell_local     => domain%num_cell_local
    ncell_boundary  => domain%num_cell_boundary
    num_proc        => domain%num_comm_proc
    cell_start      => domain%cell_start
    cell_end        => domain%cell_end
    iproc           => domain%iproc      
    cell_start_proc => domain%cell_start_proc
    cell_end_proc   => domain%cell_end_proc
    cell_g2b        => domain%cell_g2b
    cell_g2l        => domain%cell_g2l
    cell_to_rank    => domain%natom_global

    cell(1)     =  boundary%num_cells_x
    cell(2)     =  boundary%num_cells_y
    cell(3)     =  boundary%num_cells_z

    call alloc_comm(comm, CommProc, num_proc, 1, 1)

    ! cell index for communication
    !
    max_cell = 0
    do k = 1, num_proc

      comm%num_cell_coord(k) = 0
      ip = iproc(k)
      jcell_start(1:3) = cell_start_proc(ip+1,1:3)
      jcell_end  (1:3) = cell_end_proc  (ip+1,1:3)

      call setup_cellc_face_alloc(k, cell, cell_to_rank,      &
                                  jcell_start, jcell_end,     &
                                  comm%num_cell_coord)
      max_cell = max(max_cell, comm%num_cell_coord(k))

    end do

    icell_start(1:3) = cell_start(1:3)
    icell_end  (1:3) = cell_end  (1:3)

    allocate(rank_g2l(nproc_city))

    do k = 1, num_proc
      ip = iproc(k)
      rank_g2l(ip+1) = k
      comm%num_cell_force(k) = 0
    end do

    call setup_cellf_face_alloc(cell, cell_to_rank,           &
                                icell_start, icell_end,       &
                                rank_g2l, comm%num_cell_force)

    do k = 1, num_proc
      max_cell = max(max_cell, comm%num_cell_force(k))
    end do
    call alloc_comm(comm, CommCell, num_proc, max_cell, MaxAtom_domain)


    Max_alloc_size1 = max(8*MaxAtom_Domain, 2*MaxBond, 4*MaxAngl, &
                          3*MaxDihe, 3*MaxStack, 60*MaxPwmCos,    &
                          24*MaxPwmCosns, 3*MaxContact)
    Max_alloc_size2 = max(6*MaxAtom_Domain, 3*MaxBond, 4*MaxAngl, &
                          6*MaxDihe, 3*MaxStack, 5*MaxPwmCos,     &
                          10*MaxPwmCosns, 4*MaxContact)
    call alloc_domain(domain, DomainPtlMove, ncell_local+ncell_boundary, &
                      num_proc, 1)
    call alloc_domain(domain, DomainPtlArray, Max_alloc_size1, &
                      Max_alloc_size2, 1)
    
    !allocation 
    !
    size1 = max(8*MaxAtom_domain*num_proc, 3*MaxContact, 1)
    size2 = max(6*MaxAtom_domain*num_proc, 4*MaxContact, 1)
    call alloc_comm(comm, CommBuffer, size1, size2, 1)

    deallocate(rank_g2l)

    return

  end subroutine setup_communicate_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate_coord_cell
  !> @brief        Determines the communication cell (coordinates).
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !! @param[out]   comm     : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate_coord_cell(boundary, domain, comm)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: cell(3)
    integer                  :: k, ip, i, ic
    integer                  :: jcell_start(3), jcell_end(3)
    integer,         pointer :: ncell_local
    integer,         pointer :: cell_start(:), cell_end(:), cell_to_rank(:,:,:)
    integer,         pointer :: num     
    integer,         pointer :: iproc(:)
    integer,         pointer :: num_proc(:), proc_list(:,:)
    integer,         pointer :: cell_start_proc(:,:), cell_end_proc(:,:)
    integer,         pointer :: cell_g2l(:), cell_g2b(:)
    integer,         pointer :: num_cellc(:)
    integer,         pointer :: ic_send(:,:)
    integer,         pointer :: if_recv(:,:)

    ncell_local     => domain%num_cell_local
    num             => domain%num_comm_proc
    cell_start      => domain%cell_start
    cell_end        => domain%cell_end
    iproc           => domain%iproc
    num_proc        => domain%num_proc
    proc_list       => domain%proc_list
    cell_start_proc => domain%cell_start_proc
    cell_end_proc   => domain%cell_end_proc
    cell_g2b        => domain%cell_g2b
    cell_g2l        => domain%cell_g2l
    cell_to_rank    => domain%natom_global
    num_cellc       => comm%num_cell_coord
    ic_send         => comm%ic_send
    if_recv         => comm%if_recv

    cell(1)     =  boundary%num_cells_x
    cell(2)     =  boundary%num_cells_y
    cell(3)     =  boundary%num_cells_z

    ! cell index for communication
    !
    do k = 1, num

      num_cellc(k) = 0
      ip = iproc(k)
      jcell_start(1:3) = cell_start_proc(ip+1,1:3)
      jcell_end  (1:3) = cell_end_proc  (ip+1,1:3)

      call setup_cellc_face(k, cell, cell_to_rank,    &
                            jcell_start, jcell_end,   &
                            cell_g2l, num_cellc,      &
                            ic_send, if_recv)

    end do

    do ip = 1, num
      do i = 1, num_cellc(ip)
        ic = ic_send(i,ip)
        num_proc(ic) = num_proc(ic) + 1
        k = num_proc(ic)
        proc_list(k,ic) = ip
      end do
    end do

    return

  end subroutine setup_communicate_coord_cell

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_communicate_force_cell
  !> @brief        Determines the communication cell (force).
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    domain   : domain information
  !! @param[out]   comm     : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_communicate_force_cell(boundary, domain, comm)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: cell(3)
    integer                  :: k
    integer                  :: icell_start(3), icell_end(3)
    integer                  :: ip
    integer,     allocatable :: rank_g2l(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: cell_start(:), cell_end(:), cell_to_rank(:,:,:)
    integer,         pointer :: num_proc
    integer,         pointer :: iproc(:)
    integer,         pointer :: cell_g2l(:), cell_g2b(:), num_domain(:)
    integer,         pointer :: num_cellf(:)
    integer,         pointer :: ic_recv(:,:)
    integer,         pointer :: if_send(:,:)

    num_domain  => boundary%num_domain

    ncell_local     => domain%num_cell_local
    num_proc        => domain%num_comm_proc
    cell_start      => domain%cell_start
    cell_end        => domain%cell_end
    iproc           => domain%iproc
    cell_g2b        => domain%cell_g2b
    cell_g2l        => domain%cell_g2l
    cell_to_rank    => domain%natom_global
    num_cellf       => comm%num_cell_force
    ic_recv         => comm%ic_recv
    if_send         => comm%if_send

    cell(1)     =  boundary%num_cells_x
    cell(2)     =  boundary%num_cells_y
    cell(3)     =  boundary%num_cells_z

    icell_start(1:3) = cell_start(1:3)
    icell_end  (1:3) = cell_end  (1:3)

    allocate(rank_g2l(nproc_city))

    do k = 1, num_proc
      ip = iproc(k)
      rank_g2l(ip+1) = k
      num_cellf(k) = 0
    end do

    call setup_cellf_face(ncell_local,                   &
                          cell, cell_to_rank,            &
                          icell_start, icell_end,        &
                          cell_g2b, rank_g2l, num_cellf, &
                          ic_recv, if_send)

    deallocate(rank_g2l)

    return

  end subroutine setup_communicate_force_cell

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_communicate_size
  !> @brief        Checks the size of the transfer data.
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[out]   comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_communicate_size(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm  

    integer                  :: i, k, ic, iproc, ip, start_ic, ix
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    integer,         pointer :: natom(:), start_atom(:), num
    integer,         pointer :: num_cellf(:)
    integer,         pointer :: num_cellc(:)
    integer,         pointer :: if_send(:,:)
    integer,         pointer :: if_recv(:,:)
    integer,         pointer :: ic_send(:,:)
    integer,         pointer :: ic_recv(:,:)
    integer,         pointer :: if_send_list(:,:)
    integer,         pointer :: if_recv_list(:,:)
    integer,         pointer :: ic_send_list(:,:)
    integer,         pointer :: ic_recv_list(:,:)
    integer,         pointer :: iproc_proc(:)
    integer,         pointer :: send_size(:,:), recv_size(:,:)
    integer,         pointer :: irequest(:,:)

    natom           => domain%num_atom
    start_atom      => domain%start_atom
    iproc_proc      => domain%iproc
    num             => domain%num_comm_proc
    num_cellc       => comm%num_cell_coord
    num_cellf       => comm%num_cell_force
    if_send         => comm%if_send
    if_recv         => comm%if_recv
    ic_send         => comm%ic_send
    ic_recv         => comm%ic_recv
    if_send_list    => comm%if_send_list
    if_recv_list    => comm%if_recv_list
    ic_send_list    => comm%ic_send_list
    ic_recv_list    => comm%ic_recv_list
    send_size       => comm%send_size
    recv_size       => comm%recv_size
    irequest        => comm%irequest 

    ! check the send size of forces 
    !
    do iproc = 1, num
      k = 0
      do i = 1, num_cellf(iproc)
        ic = if_send(i,iproc)
        k  = k + natom(ic)
      end do
      send_size(1,iproc) = k
    end do

    ! check the send size of coordinates
    !
    do iproc = 1, num
      k = 0
      do i = 1, num_cellc(iproc)
        ic = ic_send(i,iproc)
        k  = k + natom(ic)
      end do
      send_size(2,iproc) = k
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the size of the data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc) 
      call mpi_irecv(recv_size(1,iproc), 1, mpi_integer, ip,            &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,     &
                     irequest(1,iproc), ierror)
      call mpi_irecv(recv_size(2,iproc), 1, mpi_integer, ip,            &
                     2*(my_city_rank+1)*nproc_city+ip, mpi_comm_city,   &
                     irequest(2,iproc), ierror)
    end do
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_isend(send_size(1,iproc), 1, mpi_integer, ip,            &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,     &
                     irequest(3,iproc), ierror)
      call mpi_isend(send_size(2,iproc), 1, mpi_integer, ip,            &
                     2*(ip+1)*nproc_city+my_city_rank, mpi_comm_city,   &
                     irequest(4,iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror) 
      call mpi_wait(irequest(2,iproc), istatus, ierror) 
      call mpi_wait(irequest(3,iproc), istatus, ierror) 
      call mpi_wait(irequest(4,iproc), istatus, ierror) 
    end do
#endif

    do iproc = 1, num
      comm%send_force_size(iproc) = send_size(1,iproc)
      comm%send_coord_size(iproc) = send_size(2,iproc)
      comm%recv_force_size(iproc) = recv_size(1,iproc) 
      comm%recv_coord_size(iproc) = recv_size(2,iproc)
    end do

    do iproc = 1, num
      k = 0
      do i = 1, num_cellc(iproc)
        ic = if_recv(i,iproc)
        k = k + natom(ic)
      end do
      if (recv_size(1,iproc) /= k) then
        call error_msg('Error: Setup_Communicate_Size - Disagreement '//&
                       'between sending and receiving processors.')
      end if
      k = 0
      do i = 1, num_cellf(iproc)
        ic = ic_recv(i,iproc)
        k = k + natom(ic)
      end do
      if (recv_size(2,iproc) /= k) then
        call error_msg('Error: Setup_Communicate_Size - Disagreement '//&
                       'between sending and receiving processors.')
      end if
    end do     

    ! Atom lists
    !
    do iproc = 1, num
      k = 0
      do i = 1, num_cellf(iproc)
        ic = if_send(i,iproc)
        start_ic = start_atom(ic)
        do ix = 1, natom(ic)
          k = k + 1
          if_send_list(k,iproc) = start_ic + ix
          ic_recv_list(k,iproc) = start_ic + ix
        end do
      end do
      k = 0
      do i = 1, num_cellc(iproc)
        ic = ic_send(i,iproc)
        start_ic = start_atom(ic)
        do ix = 1, natom(ic)
          k = k + 1
          ic_send_list(k,iproc) = start_ic + ix
          if_recv_list(k,iproc) = start_ic + ix
        end do
      end do
    end do
 
    return

  end subroutine update_communicate_size

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_coor
  !> @brief        Packs and transfers the coordinate data.
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
    integer                  :: iproc, i, k, ix, ip
    integer                  :: start_iproc
    integer,     allocatable :: irequest1(:), irequest2(:)
    integer,     allocatable :: send1(:), recv1(:)
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: coord(:,:)
    real(wip),       pointer :: buf_send(:), buf_recv(:)
    integer,         pointer :: natom(:), start_atom(:)
    integer,         pointer :: num
    integer,         pointer :: ic_send(:,:), ic_recv(:,:)
    integer,         pointer :: ic_send_list(:,:), ic_recv_list(:,:)
    integer,         pointer :: send(:), recv(:)
    integer,         pointer :: iproc_proc(:)

    coord           => domain%coord
    iproc_proc      => domain%iproc
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    num             => domain%num_comm_proc

    buf_send        => comm%buf_send
    buf_recv        => comm%buf_recv
    ic_send         => comm%ic_send
    ic_recv         => comm%ic_recv
    ic_send_list    => comm%ic_send_list
    ic_recv_list    => comm%ic_recv_list
    send            => comm%send_coord_size
    recv            => comm%recv_coord_size

    k = num
    allocate(irequest1(k), irequest2(k), send1(k), recv1(k))

    send1(1) = 0
    recv1(1) = 0
    i = 0
    k = 0
    do iproc = 2, num
      i = i + 3*send(iproc-1)
      k = k + 3*recv(iproc-1)
      send1(iproc) = i
      recv1(iproc) = k
    end do

    ! Pack the coordinate data 
    !
    do iproc = 1, num
      start_iproc = send1(iproc)
      !$omp parallel do private(i,ix)
      do i = 1, send(iproc)
        ix = ic_send_list(i,iproc)
        buf_send(3*i-2+start_iproc) = coord(ix,1)
        buf_send(3*i-1+start_iproc) = coord(ix,2)
        buf_send(3*i  +start_iproc) = coord(ix,3)
      end do
      !$omp end parallel do
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(iproc)+1), 3*recv(iproc),           &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,  &
                     mpi_comm_city, irequest1(iproc), ierror)
    end do
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_isend(buf_send(send1(iproc)+1), 3*send(iproc),           &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,  &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

#endif

    ! get the coordinate data 
    !
    do iproc = 1, num
      start_iproc = recv1(iproc)
      !$omp parallel do private(i,ix)
      do i = 1, recv(iproc)
        ix = ic_recv_list(i,iproc)
        coord(ix,1) = buf_recv(3*i-2+start_iproc)
        coord(ix,2) = buf_recv(3*i-1+start_iproc)
        coord(ix,3) = buf_recv(3*i  +start_iproc)
      end do
      !$omp end parallel do
    end do

    deallocate(irequest1, irequest2, send1, recv1)

    return

  end subroutine communicate_coor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_force
  !> @brief        Packs and transfers the force data. 
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[inout] comm   : communication information
  !! @param[inout] force  : forces of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_force(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, k, ix, iproc, ip
    integer                  :: start_iproc
    integer,     allocatable :: irequest1(:), irequest2(:)
    integer,     allocatable :: send1(:), recv1(:)
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: buf_send(:), buf_recv(:)
    integer,         pointer :: natom(:), start_atom(:)
    integer,         pointer :: num
    integer,         pointer :: if_send(:,:)
    integer,         pointer :: if_recv(:,:)
    integer,         pointer :: if_send_list(:,:)
    integer,         pointer :: if_recv_list(:,:)
    integer,         pointer :: send(:), recv(:)
    integer,         pointer :: iproc_proc(:)

    iproc_proc      => domain%iproc
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    num             => domain%num_comm_proc
    force           => domain%force

    buf_send        => comm%buf_send
    buf_recv        => comm%buf_recv
    if_send         => comm%if_send
    if_recv         => comm%if_recv
    if_send_list    => comm%if_send_list
    if_recv_list    => comm%if_recv_list
    send            => comm%send_force_size
    recv            => comm%recv_force_size

    k = num
    allocate(irequest1(k), irequest2(k), send1(k), recv1(k))

    send1(1) = 0
    recv1(1) = 0
    i = 0
    k = 0
    do iproc = 2, num
      i = i + 3*send(iproc-1)
      k = k + 3*recv(iproc-1)
      send1(iproc) = i
      recv1(iproc) = k
    end do

    ! Pack the force data 
    !
    do iproc = 1, num
      start_iproc = send1(iproc)
      !$omp parallel do private(i, ix)
      do i = 1, send(iproc)
        ix = if_send_list(i,iproc)
        buf_send(3*i-2+start_iproc) = force(ix,1)
        buf_send(3*i-1+start_iproc) = force(ix,2)
        buf_send(3*i  +start_iproc) = force(ix,3)
      end do
      !$omp end parallel do
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(iproc)+1), 3*recv(iproc),           &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,  &
                     mpi_comm_city, irequest1(iproc), ierror)
    end do
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_isend(buf_send(send1(iproc)+1), 3*send(iproc),           &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,  &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the force
    !
    do iproc = 1, num
      start_iproc = recv1(iproc)
      !$omp parallel do private(i, ix)
      do i = 1, recv(iproc)
        ix = if_recv_list(i,iproc)
        force(ix,1) = force(ix,1) + buf_recv(3*i-2+start_iproc)
        force(ix,2) = force(ix,2) + buf_recv(3*i-1+start_iproc)
        force(ix,3) = force(ix,3) + buf_recv(3*i  +start_iproc)
      end do
      !$omp end parallel do
    end do

    deallocate(irequest1, irequest2, send1, recv1)

    return

  end subroutine communicate_force

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_ptl
  !> @brief        Handles particle migration between processes.
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
    integer                  :: i, j, k, ic, list, iproc, ip
    integer                  :: k1, j1
    integer                  :: start_iproc(2)
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: buf_send(:), buf_recv(:)
    real(wip),       pointer :: ptl_comm_real(:)
    integer,         pointer :: num
    integer,         pointer :: send_size(:,:), recv_size(:,:)
    integer,         pointer :: send1(:,:), recv1(:,:)
    integer,         pointer :: irequest(:,:)
    integer,         pointer :: iproc_proc(:)
    integer,         pointer :: charge_move(:), nocharge_move(:)
    integer,         pointer :: charge_comm(:), nocharge_comm(:)
    integer,         pointer :: charge_comm_move(:), nocharge_comm_move(:)
    integer,         pointer :: int_send(:), int_recv(:)
    integer,         pointer :: cell_g2l(:)
    integer,         pointer :: ptl_comm_int(:)

    num                => domain%num_comm_proc
    iproc_proc         => domain%iproc
    cell_g2l           => domain%cell_g2l
    charge_move        => domain%type1_move
    nocharge_move      => domain%type2_move
    charge_comm        => domain%type1_comm
    nocharge_comm      => domain%type2_comm
    charge_comm_move   => domain%type1_comm_move
    nocharge_comm_move => domain%type2_comm_move
    ptl_comm_real      => domain%buf_var0_comm_real
    ptl_comm_int       => domain%buf_var0_comm_int

    send_size          => comm%send_size
    recv_size          => comm%recv_size
    send1              => comm%send_address
    recv1              => comm%recv_address
    irequest           => comm%irequest
    buf_send           => comm%buf_send 
    buf_recv           => comm%buf_recv 
    int_recv           => comm%int_recv 
    int_send           => comm%int_send 

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 2 + 5*charge_comm(iproc) + 5*nocharge_comm(iproc)
      j =     8*charge_comm(iproc) + 8*nocharge_comm(iproc)
      int_send(1+start_iproc(1)) = charge_comm(iproc)
      int_send(2+start_iproc(1)) = nocharge_comm(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    charge_comm(1:num) = 0
    nocharge_comm(1:num) = 0
    do i = 1, domain%charge_comm_domain
      iproc = ptl_comm_int(6*i-5)
      start_iproc(1:2) = send1(1:2,iproc)
      charge_comm(iproc) = charge_comm(iproc) + 1
      k1 = start_iproc(1) + 2 + 5*(charge_comm(iproc)-1)
      int_send(k1+1:k1+5) = ptl_comm_int(6*i-4:6*i)
      j1 = start_iproc(2) + 8*(charge_comm(iproc)-1)
      buf_send(j1+1:j1+8) = ptl_comm_real(8*i-7:8*i)
    end do
    do i = domain%charge_comm_domain+1, domain%charge_comm_domain+domain%nocharge_comm_domain
      iproc = ptl_comm_int(6*i-5)
      start_iproc(1:2) = send1(1:2,iproc)
      nocharge_comm(iproc) = nocharge_comm(iproc) + 1
      k1 = start_iproc(1) + 2 + 5*(charge_comm(iproc)+nocharge_comm(iproc)-1)
      int_send(k1+1:k1+5) = ptl_comm_int(6*i-4:6*i)
      j1 = start_iproc(2) + 8*(charge_comm(iproc)+nocharge_comm(iproc)-1)
      buf_send(j1+1:j1+8) = ptl_comm_real(8*i-7:8*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,           &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,    &
                     irequest(1,iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,           &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,    &
                     irequest(2,iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! send the data (integer values)
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),  &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest(1,iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),  &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest(2,iproc), ierror)
    end do

    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do

    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),  &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,&
                     mpi_comm_city, irequest(1,iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),  &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,&
                     mpi_comm_city, irequest(2,iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    charge_comm(1:num) = 0
    nocharge_comm(1:num) = 0
    list = 0

    do iproc = 1, num

      start_iproc(1:2) = recv1(1:2,iproc)
      charge_comm(iproc) = int_recv(1+start_iproc(1))
      nocharge_comm(iproc) = int_recv(2+start_iproc(1))

      k = 2
      j = 0
 
      do i = 1, charge_comm(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 5*(i-1)
        ptl_comm_int (5*list-4:5*list) = int_recv(k1+1:k1+5) 
        ic = cell_g2l(ptl_comm_int(5*list-4))
        charge_move(ic) = charge_move(ic) + 1
        j1 = start_iproc(2) + 8*(i-1)
        ptl_comm_real(8*list-7:8*list) = buf_recv(j1+1:j1+8)
      end do
      do i = charge_comm(iproc)+1, charge_comm(iproc)+nocharge_comm(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 5*(i-1)
        ptl_comm_int (5*list-4:5*list) = int_recv(k1+1:k1+5) 
        ic = cell_g2l(ptl_comm_int(5*list-4))
        nocharge_move(ic) = nocharge_move(ic) + 1
        j1 = start_iproc(2) + 8*(i-1)
        ptl_comm_real(8*list-7:8*list) = buf_recv(j1+1:j1+8)
      end do
     
    end do
    domain%charge_comm_domain = list
 
    return

  end subroutine communicate_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_ptl_martini
  !> @brief        Handles particle migration between processes.
  !! @authors      JJ  
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_ptl_martini(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, j, k, ic, list, iproc, ip
    integer                  :: k1, j1
    integer                  :: start_iproc(2)
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    real(wip),       pointer :: buf_send(:), buf_recv(:)
    real(wip),       pointer :: ptl_comm_real(:)
    integer,         pointer :: num
    integer,         pointer :: send_size(:,:), recv_size(:,:)
    integer,         pointer :: send1(:,:), recv1(:,:)
    integer,         pointer :: irequest(:,:)
    integer,         pointer :: iproc_proc(:)
    integer,         pointer :: charge_move(:), nocharge_move(:)
    integer,         pointer :: charge_comm(:), nocharge_comm(:)
    integer,         pointer :: charge_comm_move(:), nocharge_comm_move(:)
    integer,         pointer :: int_send(:), int_recv(:)
    integer,         pointer :: cell_g2l(:)
    integer,         pointer :: ptl_comm_int(:)

    num                => domain%num_comm_proc
    iproc_proc         => domain%iproc
    cell_g2l           => domain%cell_g2l
    charge_move        => domain%type1_move
    nocharge_move      => domain%type2_move
    charge_comm        => domain%type1_comm
    nocharge_comm      => domain%type2_comm
    charge_comm_move   => domain%type1_comm_move
    nocharge_comm_move => domain%type2_comm_move
    ptl_comm_real      => domain%buf_var0_comm_real
    ptl_comm_int       => domain%buf_var0_comm_int

    send_size          => comm%send_size
    recv_size          => comm%recv_size
    send1              => comm%send_address
    recv1              => comm%recv_address
    irequest           => comm%irequest
    buf_send           => comm%buf_send 
    buf_recv           => comm%buf_recv 
    int_recv           => comm%int_recv 
    int_send           => comm%int_send 

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 2 + 3*charge_comm(iproc) + 3*nocharge_comm(iproc)
      j =     8*charge_comm(iproc) + 8*nocharge_comm(iproc)
      int_send(1+start_iproc(1)) = charge_comm(iproc)
      int_send(2+start_iproc(1)) = nocharge_comm(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    charge_comm(1:num) = 0
    nocharge_comm(1:num) = 0
    do i = 1, domain%charge_comm_domain
      iproc = ptl_comm_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      charge_comm(iproc) = charge_comm(iproc) + 1
      k1 = start_iproc(1) + 2 + 3*(charge_comm(iproc)-1)
      int_send(k1+1:k1+3) = ptl_comm_int(4*i-2:4*i)
      j1 = start_iproc(2) + 8*(charge_comm(iproc)-1)
      buf_send(j1+1:j1+8) = ptl_comm_real(8*i-7:8*i)
    end do
    do i = domain%charge_comm_domain+1, domain%charge_comm_domain+domain%nocharge_comm_domain
      iproc = ptl_comm_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      nocharge_comm(iproc) = nocharge_comm(iproc) + 1
      k1 = start_iproc(1) + 2 + 3*(charge_comm(iproc)+nocharge_comm(iproc)-1)
      int_send(k1+1:k1+3) = ptl_comm_int(4*i-2:4*i)
      j1 = start_iproc(2) + 8*(charge_comm(iproc)+nocharge_comm(iproc)-1)
      buf_send(j1+1:j1+8) = ptl_comm_real(8*i-7:8*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,           &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,    &
                     irequest(1,iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,           &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,    &
                     irequest(2,iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! send the data (integer values)
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),  &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest(1,iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),  &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest(2,iproc), ierror)
    end do

    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do

    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),  &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,&
                     mpi_comm_city, irequest(1,iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),  &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,&
                     mpi_comm_city, irequest(2,iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest(1,iproc), istatus, ierror)
      call mpi_wait(irequest(2,iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    charge_comm(1:num) = 0
    nocharge_comm(1:num) = 0
    list = 0

    do iproc = 1, num

      start_iproc(1:2) = recv1(1:2,iproc)
      charge_comm(iproc) = int_recv(1+start_iproc(1))
      nocharge_comm(iproc) = int_recv(2+start_iproc(1))

      k = 2
      j = 0
 
      do i = 1, charge_comm(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 3*(i-1)
        ptl_comm_int (3*list-2:3*list) = int_recv(k1+1:k1+3) 
        ic = cell_g2l(ptl_comm_int(3*list-2))
        charge_move(ic) = charge_move(ic) + 1
        j1 = start_iproc(2) + 8*(i-1)
        ptl_comm_real(8*list-7:8*list) = buf_recv(j1+1:j1+8)
      end do
      do i = charge_comm(iproc)+1, charge_comm(iproc)+nocharge_comm(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 3*(i-1)
        ptl_comm_int (3*list-2:3*list) = int_recv(k1+1:k1+3) 
        ic = cell_g2l(ptl_comm_int(3*list-2))
        nocharge_move(ic) = nocharge_move(ic) + 1
        j1 = start_iproc(2) + 8*(i-1)
        ptl_comm_real(8*list-7:8*list) = buf_recv(j1+1:j1+8)
      end do
     
    end do
    domain%charge_comm_domain = list
 
    return

  end subroutine communicate_ptl_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_cell_size
  !> @brief        Updates the cell size.
  !! @authors      JJ  
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_size(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    call update_cell_size1(domain, comm)
    call update_communicate_size(domain, comm)

    return
  
  end subroutine update_cell_size

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_cell_size1
  !> @brief        Updates the cell size.
  !! @authors      JJ  
  !! @param[inout] domain : domain information
  !! @param[inout] comm   : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_size1(domain, comm)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_comm),    target, intent(inout) :: comm

    ! local variable
    integer                  :: i, k, ic, ncell, ncell_local
    integer                  :: max_cell, iproc, ip
    integer,     allocatable :: irequest1(:), irequest2(:)
    integer,     allocatable :: irequest3(:), irequest4(:)
    integer,     allocatable :: send_size(:,:), recv_size(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                  :: istatus(mpi_status_size)
#endif

    integer,         pointer :: num
    integer,         pointer :: num_cellc(:)
    integer,         pointer :: num_cellf(:)
    integer,         pointer :: natom(:), start_atom(:)
    integer,         pointer :: ncharge(:)
    integer,         pointer :: if_send(:,:)
    integer,         pointer :: if_recv(:,:)
    integer,         pointer :: ic_send(:,:)
    integer,         pointer :: ic_recv(:,:)
    integer,         pointer :: iproc_proc(:)
    integer,         pointer :: send(:)
    integer,         pointer :: recv(:)
    integer,         pointer :: fsend(:)
    integer,         pointer :: frecv(:)

    num             => domain%num_comm_proc
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    ncharge         => domain%num_charge
    iproc_proc      => domain%iproc

    send            => comm%send_coord_size
    recv            => comm%recv_coord_size
    fsend           => comm%send_force_size
    frecv           => comm%recv_force_size

    num_cellc       => comm%num_cell_coord
    num_cellf       => comm%num_cell_force
    if_send         => comm%if_send 
    if_recv         => comm%if_recv 
    ic_send         => comm%ic_send 
    ic_recv         => comm%ic_recv 

    k = num
    max_cell = 0
    do iproc = 1, num
      max_cell = max(max_cell, num_cellc(iproc), num_cellf(iproc))
    end do
    allocate(send_size(2*max_cell,k), recv_size(2*max_cell,k), &
             irequest1(k), irequest2(k), irequest3(k), irequest4(k))

    ! check the send size of coordinates
    !
    do iproc = 1, num
      k = 0
      do i = 1, num_cellc(iproc)
        ic = ic_send(i,iproc)
        send_size(2*i-1,iproc) = ncharge(ic)
        send_size(2*i  ,iproc) = natom  (ic)
        k = k + natom(ic)
      end do
      send(iproc) = k
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the size of the data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2*num_cellf(iproc),            &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,   &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2*num_cellc(iproc),            &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,   &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif
    do iproc = 1, num
      k = 0
      do i = 1, num_cellf(iproc)
        ic = ic_recv(i,iproc)
        ncharge(ic) = recv_size(2*i-1,iproc)
        natom  (ic) = recv_size(2*i  ,iproc)
        k = k + natom(ic)
      end do
      recv(iproc) = k
    end do
  
    deallocate(recv_size, send_size, irequest1, irequest2, irequest3, &
               irequest4)

    ncell_local = domain%num_cell_local
    ncell       = ncell_local + domain%num_cell_boundary
    k = domain%num_atom_domain
    do i = ncell_local+1, ncell
      domain%start_atom(i) = k
      k = k + natom(i)
    end do

    return

  end subroutine update_cell_size1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    updae_cell_boundary
  !> @brief        Updates coordinates of the boundary cell.
  !! @authors      JJ  
  !! @param[inout] domain   : domain information
  !! @param[inout] comm     : communication information
  !! @param[inout] boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_boundary(domain, enefunc, comm, boundary)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_comm),     target, intent(inout) :: comm
    type(s_boundary), target, intent(inout) :: boundary

    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3)
    integer                   :: i, k, ic, ix, ixx, start_i, list
    integer                   :: num_dna, num_base, num_elec
    integer                   :: num_idr_kh, num_idr_hps, num_kh
    integer                   :: num_idr_mpipi
    integer                   :: kx, jx, base_type
    integer                   :: id
    integer                   :: ixx1, start_i1

    real(wip),        pointer :: coord(:,:), velocity(:,:)
    real(wp),         pointer :: trans(:,:)
    real(wp),         pointer :: charge(:), cell_pbc_move(:,:)
    real(wip),        pointer :: mass(:)
    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    integer,          pointer :: start_atom(:)
    integer,          pointer :: atom_2_cell(:)
    integer,          pointer :: ncell_local, ncell_bd
    integer,          pointer :: natom(:), ncharge(:), atmcls(:)
    integer,          pointer :: nbase(:), nphos(:)
    integer(1),       pointer :: charge_type(:), cg_pro_use_KH(:)
    integer(1),       pointer :: cg_IDR_KH(:), cg_IDR_HPS(:)
    integer(1),       pointer :: cg_IDR_MPIPI(:)
    integer(1),       pointer :: dna_check(:)
    integer,          pointer :: base_list(:), phos_list(:)
    integer,          pointer :: id_l2g(:), id_g2l(:), chain_id(:)
    integer,          pointer :: atom_type(:)
    integer,          pointer :: cg_elec_list(:), cg_elec_list_inv(:)
    integer,          pointer :: cg_dna_list(:), cg_dna_list_inv(:)
    integer,          pointer :: cg_base_list(:), cg_base_list_inv(:)
    integer,          pointer :: cg_kh_list(:), cg_kh_list_inv(:)
    integer,          pointer :: cg_idr_kh_list(:), cg_idr_kh_list_inv(:)
    integer,          pointer :: cg_idr_hps_list(:), cg_idr_hps_list_inv(:)
    integer,          pointer :: cg_idr_mpipi_list(:), cg_idr_mpipi_list_inv(:)

    ncell_local         => domain%num_cell_local
    ncell_bd            => domain%num_cell_boundary
    natom               => domain%num_atom
    ncharge             => domain%num_charge
    start_atom          => domain%start_atom
    atom_2_cell         => domain%atom_2_cell
    charge_type         => domain%charge_type
    coord               => domain%coord
    velocity            => domain%velocity
    charge              => domain%charge
    mass                => domain%mass
    atmcls              => domain%atom_cls_no
    atom_type           => domain%NA_base_type
    id_l2g              => domain%id_l2g
    id_g2l              => domain%id_g2l
    chain_id            => domain%mol_chain_id
    trans               => domain%trans_vec
    cell_pbc_move       => domain%cell_pbc_move
    cg_pro_use_KH       => domain%cg_pro_use_KH
    cg_IDR_KH           => domain%cg_IDR_KH     
    cg_IDR_HPS          => domain%cg_IDR_HPS    
    cg_IDR_MPIPI        => domain%cg_IDR_MPIPI
    dna_check           => domain%dna_check
    nbase               => domain%num_base
    nphos               => domain%num_phos
    base_list           => domain%base_list
    phos_list           => domain%phos_list

    cg_elec_list        => enefunc%cg_elec_list
    cg_elec_list_inv    => enefunc%cg_elec_list_inv
    cg_dna_list         => enefunc%cg_dna_list
    cg_dna_list_inv     => enefunc%cg_dna_list_inv
    cg_base_list        => enefunc%cg_base_list
    cg_base_list_inv    => enefunc%cg_base_list_inv
    cg_kh_list          => enefunc%cg_kh_list
    cg_kh_list_inv      => enefunc%cg_kh_list_inv
    cg_idr_kh_list      => enefunc%cg_idr_kh_list
    cg_idr_kh_list_inv  => enefunc%cg_idr_kh_list_inv
    cg_idr_hps_list     => enefunc%cg_idr_hps_list
    cg_idr_hps_list_inv => enefunc%cg_idr_hps_list_inv
    cg_idr_mpipi_list     => enefunc%cg_idr_mpipi_list
    cg_idr_mpipi_list_inv => enefunc%cg_idr_mpipi_list_inv

    bsize_x       => boundary%box_size_x
    bsize_y       => boundary%box_size_y
    bsize_z       => boundary%box_size_z

    call update_cell_boundary1(domain, comm)

    k = 0
    do i = ncell_local+1, ncell_local+ncell_bd
      k = k + natom(i)
    end do
    domain%num_atom_boundary = k

    do i = 1, ncell_local+ncell_bd
      start_i1 = start_atom(i)
      do ix = 1, natom(i)
        ixx1 = ix + start_i1
        atom_2_cell(ixx1) = i
      end do
    end do

    !$omp parallel                                                &
    !$omp private(id, i, ix, x_shift, y_shift, z_shift, move, ic, &
    !$omp         start_i, list)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! check the translation for each particle
    !
    if (boundary%type == BoundaryTypePBC) then

      do i = id+1, domain%num_atom_domain, nthread
        x_shift = coord(i,1)-boundary%origin_x
        y_shift = coord(i,2)-boundary%origin_y
        z_shift = coord(i,3)-boundary%origin_z
        move(1) = bsize_x*0.5_wip-bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip-bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip-bsize_z*anint(z_shift/bsize_z)
        x_shift = x_shift + move(1)
        y_shift = y_shift + move(2)
        z_shift = z_shift + move(3)
        trans(i,1) = move(1) - boundary%origin_x
        trans(i,2) = move(2) - boundary%origin_y
        trans(i,3) = move(3) - boundary%origin_z
      end do

      do i = id+domain%num_atom_domain+1, &
             domain%num_atom_domain+domain%num_atom_boundary, nthread
        ic = atom_2_cell(i)
        x_shift = coord(i,1)-boundary%origin_x
        y_shift = coord(i,2)-boundary%origin_y
        z_shift = coord(i,3)-boundary%origin_z
        move(1) = bsize_x*0.5_wip-bsize_x*anint(x_shift/bsize_x)
        move(2) = bsize_y*0.5_wip-bsize_y*anint(y_shift/bsize_y)
        move(3) = bsize_z*0.5_wip-bsize_z*anint(z_shift/bsize_z)
        move(1) = move(1) + cell_pbc_move(1,ic)*bsize_x
        move(2) = move(2) + cell_pbc_move(2,ic)*bsize_y
        move(3) = move(3) + cell_pbc_move(3,ic)*bsize_z
        trans(i,1) = move(1) - boundary%origin_x
        trans(i,2) = move(2) - boundary%origin_y
        trans(i,3) = move(3) - boundary%origin_z
      end do

    end if

    !$omp end parallel

    if (enefunc%forcefield == ForcefieldResidCG) then

      num_kh = 0
      num_dna = 0
      num_base = 0
      num_idr_kh = 0
      num_idr_hps = 0
      num_idr_mpipi = 0
      do i = 1, ncell_local+ncell_bd
        kx = 0
        jx = 0
        start_i = start_atom(i)
        do ix = 1, natom(i)
          ixx = ix + start_i
          cg_pro_use_KH(ixx) = 0
          cg_IDR_KH(ixx) = 0
          cg_IDR_HPS(ixx) = 0
          cg_IDR_MPIPI(ixx) = 0
          if (enefunc%cg_ele_calc) cg_elec_list_inv(ixx) = 0
          if (enefunc%cg_DNA_exv_calc) cg_dna_list_inv(ixx) = 0
          if (enefunc%cg_DNA_base_pair_calc) cg_base_list_inv(ixx) = 0
          if (enefunc%cg_KH_calc) cg_kh_list_inv(ixx) = 0
          if (enefunc%cg_IDR_KH_calc) cg_idr_kh_list_inv(ixx) = 0
          if (enefunc%cg_IDR_HPS_calc) cg_idr_hps_list_inv(ixx) = 0
          if (enefunc%cg_IDR_MPIPI_calc) cg_idr_mpipi_list_inv(ixx) = 0
          k = atom_type(ixx)
          if (k <= NABaseTypeDBMax .or. k == NABaseTypeDP .or. &
              k == NABaseTypeDS) then
            num_dna = num_dna + 1
            cg_dna_list(num_dna) = ixx
            cg_dna_list_inv(ixx) = num_dna
            if (k <= NABaseTypeDBMax) then
              num_base = num_base + 1
              cg_base_list(num_base) = ixx
              cg_base_list_inv(ixx) = num_base
              dna_check(ixx) = 1
              kx = kx + 1
              base_list(kx+start_i) = ixx
            else if (k == NABaseTypeDP) then
              dna_check(ix+start_i) = 2
              jx = jx + 1
              phos_list(jx+start_i) = ix + start_i
            else if (k == NABaseTypeDS) then
              dna_check(ix+start_i) = 3
            end if
          else
            dna_check(ix+start_i) = 0
            if (atom_type(ixx) == NABaseTypeKH .or.       &
                atom_type(ixx) == NABaseTypeBothKH) then
              cg_pro_use_KH(ixx) = 1
              num_kh = num_kh + 1
              cg_kh_list(num_kh) = ixx
              cg_kh_list_inv(ixx) = num_kh
            end if
            if (atom_type(ixx) == NABaseTypeIDRKH .or.    &
                atom_type(ixx) == NABaseTypeBothKH) then
              cg_IDR_KH(ixx) = 1
              num_idr_kh = num_idr_kh + 1
              cg_idr_kh_list(num_idr_kh) = ixx
              cg_idr_kh_list_inv(ixx) = num_idr_kh
            end if
            if (atom_type(ix+start_i) == NABaseTypeIDRHPS) then
              cg_IDR_HPS(ix+start_i) = 1
              num_idr_hps = num_idr_hps + 1
              cg_idr_hps_list(num_idr_hps) = ixx
              cg_idr_hps_list_inv(ixx) = num_idr_hps
            end if
            if (atom_type(ix+start_i) == NABaseTypeIDRMPIPI) then
              cg_IDR_MPIPI(ix+start_i) = 1
              num_idr_mpipi = num_idr_mpipi + 1
              cg_idr_mpipi_list(num_idr_mpipi) = ixx
              cg_idr_mpipi_list_inv(ixx) = num_idr_mpipi
            end if
          end if
        end do
        nbase(i) = kx
        nphos(i) = jx
      end do
      enefunc%num_cg_DNA = num_dna
      enefunc%num_cg_base = num_base
      enefunc%num_cg_KH = num_kh
      enefunc%num_cg_IDR_KH = num_idr_kh
      enefunc%num_cg_IDR_HPS = num_idr_hps
      enefunc%num_cg_IDR_MPIPI = num_idr_mpipi

    end if

    ! molecule type, coordinate, charge in elec cell
    !
    num_elec = 0
    if (enefunc%forcefield /= ForcefieldGroMartini) then
      do i = 1, ncell_local+ncell_bd
        start_i = start_atom(i)
        do ix = 1, ncharge(i)
          ixx = ix + start_i
          base_type = atom_type(ixx)
          if (base_type == NABaseTypeDP) then
            charge_type(ixx) = 1
          else if (base_type > NABaseTypeNAMax) then
            charge_type(ixx) = 2
          else
            charge_type(ixx) = 3
          end if
          num_elec = num_elec + 1
          cg_elec_list(num_elec) = ixx
          cg_elec_list_inv(ixx) = num_elec
        end do
      end do
    else
      do i = 1, ncell_local+ncell_bd
        start_i = start_atom(i)
        do ix = 1, ncharge(i)
          ixx = ix + start_i
          num_elec = num_elec + 1
          cg_elec_list(num_elec) = ixx
          cg_elec_list_inv(ixx) = num_elec
        end do
      end do
    end if  
    enefunc%num_cg_elec = num_elec

    return

  end subroutine update_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    updae_cell_boundary1
  !> @brief        Updates coordinates of the boundary cell.
  !! @authors      JJ  
  !! @param[inout] domain   : domain information
  !! @param[inout] comm     : communication information
  !! @param[inout] boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_cell_boundary1(domain, comm)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm

    ! local variable
    integer                  :: i, k, list, ip, iproc
    integer                  :: start_iproc(2)
    integer,     allocatable :: irequest1(:), irequest2(:)
    integer,     allocatable :: irequest3(:), irequest4(:)
    integer,     allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: coord(:,:), velocity(:,:)
    real(wp),         pointer :: charge(:)
    real(wip),        pointer :: mass(:)
    real(wip),        pointer :: buf_send(:), buf_recv(:)
    integer,          pointer :: num
    integer,          pointer :: start_atom(:)
    integer,          pointer :: ic_send(:,:)
    integer,          pointer :: ic_recv(:,:)
    integer,          pointer :: ic_send_list(:,:)
    integer,          pointer :: ic_recv_list(:,:)
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: recv(:)
    integer,          pointer :: send(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: ncell_local, ncell_bd
    integer,          pointer :: natom(:), ncharge(:), atmcls(:)
    integer,          pointer :: id_l2g(:), id_g2l(:), chain_id(:)
    integer,          pointer :: atom_type(:)

    num              => domain%num_comm_proc
    ncell_local      => domain%num_cell_local
    ncell_bd         => domain%num_cell_boundary
    natom            => domain%num_atom
    ncharge          => domain%num_charge
    start_atom       => domain%start_atom
    coord            => domain%coord
    velocity         => domain%velocity
    charge           => domain%charge
    mass             => domain%mass
    atmcls           => domain%atom_cls_no
    atom_type        => domain%NA_base_type
    id_l2g           => domain%id_l2g
    id_g2l           => domain%id_g2l
    chain_id         => domain%mol_chain_id
    iproc_proc       => domain%iproc

    buf_send         => comm%buf_send
    buf_recv         => comm%buf_recv
    int_send         => comm%int_send
    int_recv         => comm%int_recv
    ic_send          => comm%ic_send 
    ic_recv          => comm%ic_recv 
    ic_send_list     => comm%ic_send_list
    ic_recv_list     => comm%ic_recv_list
    recv             => comm%recv_coord_size
    send             => comm%send_coord_size

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send1(2,k), recv1(2,k))

    send1(1:2,1) = 0
    recv1(1:2,1) = 0
    i = 0
    k = 0
    do iproc = 2, num
      i = i + send(iproc-1)
      k = k + recv(iproc-1)
      send1(1,iproc) = 4*i
      recv1(1,iproc) = 4*k
      send1(2,iproc) = 8*i
      recv1(2,iproc) = 8*k
    end do

    ! Pack the data 
    !
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      do i = 1, send(iproc)
        list = ic_send_list(i,iproc)
        buf_send(8*i-7+start_iproc(2)) = coord    (list,1)
        buf_send(8*i-6+start_iproc(2)) = coord    (list,2)
        buf_send(8*i-5+start_iproc(2)) = coord    (list,3)
        buf_send(8*i-4+start_iproc(2)) = velocity (list,1)
        buf_send(8*i-3+start_iproc(2)) = velocity (list,2)
        buf_send(8*i-2+start_iproc(2)) = velocity (list,3)
        buf_send(8*i-1+start_iproc(2)) = charge   (list)
        buf_send(8*i  +start_iproc(2)) = mass     (list)
        int_send(4*i-3+start_iproc(1)) = atmcls   (list)
        int_send(4*i-2+start_iproc(1)) = id_l2g   (list)
        int_send(4*i-1+start_iproc(1)) = chain_id (list)
        int_send(4*i  +start_iproc(1)) = atom_type(list)
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    ! send the data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), 4*recv(iproc),          &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,    &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), 8*recv(iproc),          &
                     mpi_wip_real, ip, 2*(my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest2(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), 4*send(iproc),          &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,    &
                     mpi_comm_city, irequest3(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), 8*send(iproc),          &
                     mpi_wip_real, ip, 2*(ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest4(iproc), ierror)
    end do

    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
      call mpi_wait(irequest3(iproc), istatus, ierror)
      call mpi_wait(irequest4(iproc), istatus, ierror)
    end do
#endif
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, recv(iproc)
        list = ic_recv_list(i,iproc)
        coord    (list,1) = buf_recv(8*i-7+start_iproc(2))
        coord    (list,2) = buf_recv(8*i-6+start_iproc(2))
        coord    (list,3) = buf_recv(8*i-5+start_iproc(2))
        velocity (list,1) = buf_recv(8*i-4+start_iproc(2))
        velocity (list,2) = buf_recv(8*i-3+start_iproc(2))
        velocity (list,3) = buf_recv(8*i-2+start_iproc(2))
        charge   (list  ) = buf_recv(8*i-1+start_iproc(2))
        mass     (list  ) = buf_recv(8*i  +start_iproc(2))
        atmcls   (list  ) = int_recv(4*i-3+start_iproc(1))
        id_l2g   (list  ) = int_recv(4*i-2+start_iproc(1))
        chain_id (list  ) = int_recv(4*i-1+start_iproc(1))
        atom_type(list  ) = int_recv(4*i  +start_iproc(1))
        id_g2l(id_l2g(list)) = list
      end do
    end do

    deallocate(irequest1, irequest2, irequest3, irequest4, send1, recv1)

    return

  end subroutine update_cell_boundary1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_bond
  !> @brief        Packs incoming bonding data for bonds.
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
    integer                   :: start_iproc(2), list
    integer                   :: i, j, k, iproc, ip, j1, k1
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wip),        pointer :: buf_bond_move_real(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: bondsq_move(:), bond_move(:)
    integer,          pointer :: buf_bond_move_int(:)

    num                     => domain%num_comm_proc      
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv
    int_send                => comm%int_send
    int_recv                => comm%int_recv

    bondsq_move             => domain%type1_comm
    bond_move               => domain%type2_comm
    buf_bond_move_real      => domain%buf_var0_comm_real
    buf_bond_move_int       => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 2 + 3*bondsq_move(iproc) + 3*bond_move(iproc)
      j =     2*bondsq_move(iproc) + 2*bond_move(iproc)
      int_send(1+start_iproc(1)) = bondsq_move(iproc)
      int_send(2+start_iproc(1)) = bond_move(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    bondsq_move(1:num) = 0
    bond_move  (1:num) = 0
    do i = 1, enefunc%bonds_comm_domain
      iproc = buf_bond_move_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      bondsq_move(iproc) = bondsq_move(iproc) + 1
      k1 = start_iproc(1) + 2 + 3*(bondsq_move(iproc)-1)
      int_send(k1+1:k1+3) = buf_bond_move_int(4*i-2:4*i)
      j1 = start_iproc(2) + 2*(bondsq_move(iproc)-1)
      buf_send(j1+1:j1+2) = buf_bond_move_real(2*i-1:2*i)
    end do
    do i = enefunc%bonds_comm_domain+1, enefunc%bonds_comm_domain+enefunc%bondq_comm_domain
      iproc = buf_bond_move_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      bond_move(iproc) = bond_move(iproc) + 1
      k1 = start_iproc(1) + 2 + 3*(bondsq_move(iproc)+bond_move(iproc)-1)
      int_send(k1+1:k1+3) = buf_bond_move_int(4*i-2:4*i)
      j1 = start_iproc(2) + 2*(bondsq_move(iproc)+bond_move(iproc)-1)
      buf_send(j1+1:j1+2) = buf_bond_move_real(2*i-1:2*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,             &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,      &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,             &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,      &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),      &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,     &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),      &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,     &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),      &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,    &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),      &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,    &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    bondsq_move(1:num) = 0
    bond_move  (1:num) = 0
    list = 0

    do iproc = 1, num

      start_iproc(1:2) = recv1(1:2,iproc)
      bondsq_move(iproc) = int_recv(1+start_iproc(1))
      bond_move(iproc) = int_recv(2+start_iproc(1))

      do i = 1, bondsq_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 3*(i-1)
        buf_bond_move_int (3*list-2:3*list) = int_recv(k1+1:k1+3)
        j1 = start_iproc(2) + 2*(i-1)
        buf_bond_move_real(2*list-1:2*list) = buf_recv(j1+1:j1+2)
        j = j + 2
      end do
      do i = bondsq_move(iproc)+1, bondsq_move(iproc)+bond_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + 3*(i-1)
        buf_bond_move_int (3*list-2:3*list) = int_recv(k1+1:k1+3)
        j1 = start_iproc(2) + 2*(i-1)
        buf_bond_move_real(2*list-1:2*list) = buf_recv(j1+1:j1+2)
        j = j + 2
      end do

    end do

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)
    return

  end subroutine communicate_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_angl
  !> @brief        Packs incoming bonding data for angles.
  !! @authors      JJ  
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_angl(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: start_iproc(2)
    integer                   :: i, j, k, iproc, ip, j1, k1
    integer                   :: nadd1, nadd2, list
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    real(wip),        pointer :: buf_angl_move_real(:)
    integer,          pointer :: anglf_move(:), angll_move(:)
    integer,          pointer :: angl_move(:)
    integer,          pointer :: buf_angl_move_int(:)

    num                     => domain%num_comm_proc      
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv
    int_send                => comm%int_send
    int_recv                => comm%int_recv

    anglf_move              => domain%type1_comm
    angll_move              => domain%type2_comm
    angl_move               => domain%type3_comm
    buf_angl_move_int       => domain%buf_var0_comm_int
    buf_angl_move_real      => domain%buf_var0_comm_real

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 3 + 4*anglf_move(iproc) + 3*angll_move(iproc) + 4*angl_move(iproc)
      j =                           3*angll_move(iproc) + 4*angl_move(iproc)
      int_send(1+start_iproc(1)) = anglf_move(iproc)
      int_send(2+start_iproc(1)) = angll_move(iproc)
      int_send(3+start_iproc(1)) = angl_move(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    anglf_move(1:num) = 0
    angll_move(1:num) = 0
    angl_move (1:num) = 0

    do i = 1, enefunc%anglf_comm_domain
      iproc = buf_angl_move_int(5*i-4)
      start_iproc(1:2) = send1(1:2,iproc)
      anglf_move(iproc) = anglf_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 4*(anglf_move(iproc)-1)
      int_send(k1+1:k1+4) = buf_angl_move_int(5*i-3:5*i)
    end do
    nadd1 = 5*enefunc%anglf_comm_domain
    do i = 1, enefunc%angll_comm_domain
      iproc = buf_angl_move_int(4*i-3+nadd1)
      start_iproc(1:2) = send1(1:2,iproc)
      angll_move(iproc) = angll_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 4*anglf_move(iproc) + 3*(angll_move(iproc)-1)
      int_send(k1+1:k1+3) = buf_angl_move_int(4*i-2+nadd1:4*i+nadd1)
      j1 = start_iproc(2) + 3*(angll_move(iproc)-1)
      buf_send(j1+1:j1+3) = buf_angl_move_real(3*i-2:3*i)
    end do
    nadd1 = 5*enefunc%anglf_comm_domain + 4*enefunc%angll_comm_domain
    nadd2 = 3*enefunc%angll_comm_domain
    do i = 1, enefunc%angl_comm_domain
      iproc = buf_angl_move_int(5*i-4+nadd1)
      start_iproc(1:2) = send1(1:2,iproc)
      angl_move(iproc) = angl_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 4*anglf_move(iproc) + 3*angll_move(iproc) &
                              + 4*(angl_move(iproc)-1)
      int_send(k1+1:k1+4) = buf_angl_move_int(5*i-3+nadd1:5*i+nadd1)
      j1 = start_iproc(2) + 3*angll_move(iproc) + 4*(angl_move(iproc)-1)
      buf_send(j1+1:j1+4) = buf_angl_move_real(4*i-3+nadd2:4*i+nadd2)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,             &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,      &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,             &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,      &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),      &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,     &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),      &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,     &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),      &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,    &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),      &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,    &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    anglf_move(1:num) = 0
    angll_move(1:num) = 0
    angl_move(1:num) = 0

    nadd1 = 0
    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      anglf_move(iproc) = int_recv(1+start_iproc(1))
      angll_move(iproc) = int_recv(2+start_iproc(1))
      angl_move(iproc)  = int_recv(3+start_iproc(1))
      do i = 1, anglf_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 4*(i-1)
        buf_angl_move_int(4*list-3:4*list) = int_recv(k1+1:k1+4)
      end do
    end do
    enefunc%anglf_comm_domain = list
    nadd1 = list * 4 

    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, angll_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 4*anglf_move(iproc) + 3*(i-1)
        buf_angl_move_int(3*list-2+nadd1:3*list+nadd1) = int_recv(k1+1:k1+3)
        j1 = start_iproc(2) + 3*(i-1)
        buf_angl_move_real(3*list-2:3*list) = buf_recv(j1+1:j1+3)
      end do
    end do
    enefunc%angll_comm_domain = list
    nadd1 = nadd1 + list * 3 
    nadd2 = list * 3 

    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, angl_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 4*anglf_move(iproc) + 3*angll_move(iproc) &
                                + 4*(i-1)
        buf_angl_move_int(4*list-3+nadd1:4*list+nadd1) = int_recv(k1+1:k1+4)
        j1 = start_iproc(2) + 3*angll_move(iproc) + 4*(i-1)
        buf_angl_move_real(4*list-3+nadd2:4*list+nadd2) = buf_recv(j1+1:j1+4)
      end do
    end do
    enefunc%angl_comm_domain = list
        
    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)
    return

  end subroutine communicate_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_dihe
  !> @brief        Packs incoming bonding data for dihedrals.
  !! @authors      JJ  
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_dihe(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: start_iproc(2)
    integer                   :: i, j, k, iproc, ip, j1, k1
    integer                   :: nadd1, nadd2, list
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    real(wip),        pointer :: buf_dihe_move_real(:)
    integer,          pointer :: dihef_move(:), dihel_move(:)
    integer,          pointer :: dihe_move(:)
    integer,          pointer :: buf_dihe_move_int(:)

    num                     => domain%num_comm_proc      
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv
    int_send                => comm%int_send
    int_recv                => comm%int_recv

    dihef_move              => domain%type1_comm
    dihel_move              => domain%type2_comm
    dihe_move               => domain%type3_comm
    buf_dihe_move_int       => domain%buf_var0_comm_int
    buf_dihe_move_real      => domain%buf_var0_comm_real

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 3 + 6*dihef_move(iproc) + 5*dihel_move(iproc) + 6*dihe_move(iproc)
      j =                           3*dihel_move(iproc) + 2*dihe_move(iproc)
      int_send(1+start_iproc(1)) = dihef_move(iproc)
      int_send(2+start_iproc(1)) = dihel_move(iproc)
      int_send(3+start_iproc(1)) = dihe_move(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    dihef_move(1:num) = 0
    dihel_move(1:num) = 0
    dihe_move (1:num) = 0

    do i = 1, enefunc%dihef_comm_domain
      iproc = buf_dihe_move_int(7*i-6)
      start_iproc(1:2) = send1(1:2,iproc)
      dihef_move(iproc) = dihef_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 6*(dihef_move(iproc)-1)
      int_send(k1+1:k1+6) = buf_dihe_move_int(7*i-5:7*i)
    end do
    nadd1 = 7*enefunc%dihef_comm_domain
    do i = 1, enefunc%dihel_comm_domain
      iproc = buf_dihe_move_int(6*i-5+nadd1)
      start_iproc(1:2) = send1(1:2,iproc)
      dihel_move(iproc) = dihel_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 6*dihef_move(iproc) + 5*(dihel_move(iproc)-1)
      int_send(k1+1:k1+5) = buf_dihe_move_int(6*i-4+nadd1:6*i+nadd1)
      j1 = start_iproc(2) + 3*(dihel_move(iproc)-1)
      buf_send(j1+1:j1+3) = buf_dihe_move_real(3*i-2:3*i)
    end do
    nadd1 = 7*enefunc%dihef_comm_domain + 6*enefunc%dihel_comm_domain
    nadd2 = 3*enefunc%dihel_comm_domain
    do i = 1, enefunc%dihe_comm_domain
      iproc = buf_dihe_move_int(7*i-6+nadd1)
      start_iproc(1:2) = send1(1:2,iproc)
      dihe_move(iproc) = dihe_move(iproc) + 1
      k1 = start_iproc(1) + 3 + 6*dihef_move(iproc) + 5*dihel_move(iproc) &
                              + 6*(dihe_move(iproc)-1)
      int_send(k1+1:k1+6) = buf_dihe_move_int(7*i-5+nadd1:7*i+nadd1)
      j1 = start_iproc(2) + 3*dihel_move(iproc) + 2*(dihe_move(iproc)-1)
      buf_send(j1+1:j1+2) = buf_dihe_move_real(2*i-1+nadd2:2*i+nadd2)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,             &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,      &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,             &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,      &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),      &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,     &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),      &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,     &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),      &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip,    &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),      &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank,    &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    dihef_move(1:num) = 0
    dihel_move(1:num) = 0
    dihe_move(1:num) = 0

    nadd1 = 0
    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      dihef_move(iproc) = int_recv(1+start_iproc(1))
      dihel_move(iproc) = int_recv(2+start_iproc(1))
      dihe_move(iproc)  = int_recv(3+start_iproc(1))
      do i = 1, dihef_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 6*(i-1)
        buf_dihe_move_int(6*list-5:6*list) = int_recv(k1+1:k1+6)
      end do
    end do
    enefunc%dihef_comm_domain = list
    nadd1 = list * 6 

    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, dihel_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 6*dihef_move(iproc) + 5*(i-1)
        buf_dihe_move_int(5*list-4+nadd1:5*list+nadd1) = int_recv(k1+1:k1+5)
        j1 = start_iproc(2) + 3*(i-1)
        buf_dihe_move_real(3*list-2:3*list) = buf_recv(j1+1:j1+3)
      end do
    end do
    enefunc%dihel_comm_domain = list
    nadd1 = nadd1 + list * 5 
    nadd2 = list * 3 

    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, dihe_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 3 + 6*dihef_move(iproc) + 5*dihel_move(iproc) &
                                + 6*(i-1)
        buf_dihe_move_int(6*list-5+nadd1:6*list+nadd1) = int_recv(k1+1:k1+6)
        j1 = start_iproc(2) + 3*dihel_move(iproc) + 2*(i-1)
        buf_dihe_move_real(2*list-1+nadd2:2*list+nadd2) = buf_recv(j1+1:j1+2)
      end do
    end do
    enefunc%dihe_comm_domain = list
        
    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)
    return

  end subroutine communicate_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_stack
  !> @brief        Communicates stack terms. 
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_stack(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, k, j1, k1, iproc, ip, list
    integer                   :: start_iproc(2)
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wip),        pointer :: buf_stack_move_real(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: stack_move(:)
    integer,          pointer :: buf_stack_move_int(:)

    num                     => domain%num_comm_proc
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv
    int_send                => comm%int_send
    int_recv                => comm%int_recv

    stack_move              => domain%type1_comm
    buf_stack_move_real     => domain%buf_var0_comm_real
    buf_stack_move_int      => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      int_send(1+start_iproc(1)) = stack_move(iproc)
      send_size(1,iproc) = 1 + 3*stack_move(iproc)
      send_size(2,iproc) = 3*stack_move(iproc)
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    stack_move(1:num) = 0

    do i = 1, enefunc%stack_comm_domain
      iproc = buf_stack_move_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      stack_move(iproc) = stack_move(iproc) + 1
      k1 = start_iproc(1) + 1 + 3*(stack_move(iproc)-1)
      int_send(k1+1:k1+3) = buf_stack_move_int(4*i-2:4*i)
      j1 = start_iproc(2) + 3*(stack_move(iproc)-1)
      buf_send(j1+1:j1+3) = buf_stack_move_real(3*i-2:3*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,           &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,    &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,           &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,    &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),   &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,  & 
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),   &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,  &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),   &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),   &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    stack_move(1:num) = 0
    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      stack_move(iproc) = int_recv(1+start_iproc(1))
      do i = 1, stack_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 1 + 3*(i-1)
        buf_stack_move_int(3*list-2:3*list) = int_recv(k1+1:k1+3)
        j1 = start_iproc(2) + 3*(i-1)
        buf_stack_move_real(3*list-2:3*list) = buf_recv(j1+1:j1+3)
      end do
    end do
    enefunc%stack_comm_domain = list

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)

    return

  end subroutine communicate_stack

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pwmcos
  !> @brief        Communicates PWMCOS terms. 
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pwmcos(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, j, k, iproc, ip
    integer                   :: j1, k1, list
    integer                   :: start_iproc(2)
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: recv1(:,:), send1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wip),        pointer :: buf_pwmcos_move_real(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: pwmcos_move(:)
    integer,          pointer :: buf_pwmcos_move_int(:)

    num                     => domain%num_comm_proc      
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_send                => comm%buf_send
    int_send                => comm%int_send
    int_recv                => comm%int_recv
    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv

    pwmcos_move             => domain%type1_comm
    buf_pwmcos_move_real    => domain%buf_var0_comm_real
    buf_pwmcos_move_int     => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 1 + 4*pwmcos_move(iproc)
      j =    60*pwmcos_move(iproc)
      int_send(1+start_iproc(1)) = pwmcos_move(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    pwmcos_move(1:num) = 0

    do i = 1, enefunc%pwmcos_comm_domain
      iproc = buf_pwmcos_move_int(5*i-4)
      start_iproc(1:2) = send1(1:2,iproc)
      pwmcos_move(iproc) = pwmcos_move(iproc) + 1
      k1 = start_iproc(1) + 1 + 4*(pwmcos_move(iproc)-1)
      int_send(k1+1:k1+4) = buf_pwmcos_move_int(5*i-3:5*i)
      j1 = start_iproc(2) + 60*(pwmcos_move(iproc)-1)
      buf_send(j1+1:j1+60) = buf_pwmcos_move_real(60*i-59:60*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,         &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,  &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,         &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,  &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),  &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),  &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),   &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),   &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    pwmcos_move(1:num) = 0
    list = 0

    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      pwmcos_move(iproc) = int_recv(1+start_iproc(1))
      do i = 1, pwmcos_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 1 + 4*(i-1)
        buf_pwmcos_move_int(4*list-3:4*list) = int_recv(k1+1:k1+4)
        j1 = start_iproc(2) + 60*(i-1)
        buf_pwmcos_move_real(60*list-59:60*list) = buf_recv(j1+1:j1+60)
      end do
    end do
    enefunc%pwmcos_comm_domain = list

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)

    return

  end subroutine communicate_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_pwmcosns
  !> @brief        Communicates PWMCOSNS terms. 
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_pwmcosns(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, j, k, iproc, ip
    integer                   :: j1, k1, list
    integer                   :: start_iproc(2)
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: recv1(:,:), send1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wip),        pointer :: buf_pwmcosns_move_real(:)
    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: pwmcosns_move(:)
    integer,          pointer :: buf_pwmcosns_move_int(:)

    num                     => domain%num_comm_proc      
    iproc_proc              => domain%iproc

    buf_send                => comm%buf_send
    buf_recv                => comm%buf_recv
    int_send                => comm%int_send
    int_recv                => comm%int_recv

    pwmcosns_move           => domain%type1_comm
    buf_pwmcosns_move_real  => domain%buf_var0_comm_real
    buf_pwmcosns_move_int   => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data 
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      k = 1 +10*pwmcosns_move(iproc)
      j =    24*pwmcosns_move(iproc)
      int_send(1+start_iproc(1)) = pwmcosns_move(iproc)
      send_size(1,iproc) = k
      send_size(2,iproc) = j
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    pwmcosns_move(1:num) = 0

    do i = 1, enefunc%pwmcosns_comm_domain
      iproc = buf_pwmcosns_move_int(11*i-10)
      start_iproc(1:2) = send1(1:2,iproc)
      pwmcosns_move(iproc) = pwmcosns_move(iproc) + 1
      k1 = start_iproc(1) + 1 + 10*(pwmcosns_move(iproc)-1)
      int_send(k1+1:k1+10) = buf_pwmcosns_move_int(11*i-9:11*i)
      j1 = start_iproc(2) + 24*(pwmcosns_move(iproc)-1)
      buf_send(j1+1:j1+24) = buf_pwmcosns_move_real(24*i-23:24*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,          &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,   &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,          &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,   &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),  &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),  &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data 
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),   &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),   &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    pwmcosns_move(1:num) = 0
    list = 0

    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      pwmcosns_move(iproc) = int_recv(1+start_iproc(1))
      do i = 1, pwmcosns_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 1 + 10*(i-1)
        buf_pwmcosns_move_int (10*list-9:10*list) = int_recv(k1+1:k1+10)
        j1 = start_iproc(2) + 24*(i-1)
        buf_pwmcosns_move_real(24*list-23:24*list) = buf_recv(j1+1:j1+24)
      end do
    end do
    enefunc%pwmcosns_comm_domain = list

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)

    return

  end subroutine communicate_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_contact
  !> @brief        Sends contact information to neighboring processes.
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_contact(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, k, j1, k1, iproc, ip, list
    integer                   :: start_iproc(2)
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: contact_move(:)
    integer,          pointer :: contactlist(:,:), contactfunc(:)
    integer,          pointer :: buf_comm_int(:)
    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wp),         pointer :: lj12(:), lj10(:), lj6(:)
    real(wip),        pointer :: buf_comm_real(:)

    num           => domain%num_comm_proc
    iproc_proc    => domain%iproc

    buf_send      => comm%buf_send
    buf_recv      => comm%buf_recv
    int_send      => comm%int_send
    int_recv      => comm%int_recv
    contactlist   => enefunc%contact_list
    contactfunc   => enefunc%contact_func
    lj12          => enefunc%contact_lj12
    lj10          => enefunc%contact_lj10
    lj6           => enefunc%contact_lj6

    contact_move  => domain%type1_comm
    buf_comm_real => domain%buf_var0_comm_real
    buf_comm_int  => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      int_send(1+start_iproc(1)) = contact_move(iproc)
      send_size(1,iproc) = 1 + 3*contact_move(iproc)
      send_size(2,iproc) = 3*contact_move(iproc)
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    contact_move(1:num) = 0
    do i = 1, enefunc%num_contact_boundary
      iproc = buf_comm_int(4*i-3)
      start_iproc(1:2) = send1(1:2,iproc)
      contact_move(iproc) = contact_move(iproc) + 1
      k1 = start_iproc(1) + 1 + 3*(contact_move(iproc)-1)
      int_send(k1+1:k1+3) = buf_comm_int(4*i-2:4*i)
      j1 = start_iproc(2) + 3*(contact_move(iproc)-1)
      buf_send(j1+1:j1+3) = buf_comm_real(3*i-2:3*i)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,           &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,    &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,           &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,    &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),   &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,  &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),   &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,  &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),   &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),   &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    contact_move(1:num) = 0
    list = enefunc%num_contact_domain
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      contact_move(iproc) = int_recv(1+start_iproc(1))
      do i = 1, contact_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 1 + 3*(i-1)
        contactlist(1,list) = int_recv(k1+1)
        contactlist(2,list) = int_recv(k1+2)
        contactfunc(  list) = int_recv(k1+3)
        j1 = start_iproc(2) + 3*(i-1)
        lj12(list) = buf_recv(j1+1)
        lj10(list) = buf_recv(j1+2)
        lj6 (list) = buf_recv(j1+3)
      end do
    end do
    enefunc%num_contact_boundary = list - enefunc%num_contact_domain

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)

    return

  end subroutine communicate_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    communicate_restraint
  !> @brief        Sends restraint information to neighboring processes.
  !! @authors      JJ
  !! @param[inout] domain  : domain information
  !! @param[inout] comm    : communication information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine communicate_restraint(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_comm),     target, intent(inout) :: comm
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variable
    integer                   :: i, k, j1, k1, iproc, ip, list
    integer                   :: start_iproc(2), nadd1, nadd2
    integer,      allocatable :: irequest1(:), irequest2(:)
    integer,      allocatable :: irequest3(:), irequest4(:)
    integer,      allocatable :: send_size(:,:)
    integer,      allocatable :: recv_size(:,:)
    integer,      allocatable :: send1(:,:), recv1(:,:)
#ifdef HAVE_MPI_GENESIS
    integer                   :: istatus(mpi_status_size)
#endif

    integer,          pointer :: num
    integer,          pointer :: iproc_proc(:)
    integer,          pointer :: int_send(:), int_recv(:)
    integer,          pointer :: rest_move(:), fit_move(:)
    integer,          pointer :: buf_comm_int(:)
    real(wip),        pointer :: buf_send(:), buf_recv(:)
    real(wip),        pointer :: buf_comm_real(:)

    num             => domain%num_comm_proc
    iproc_proc      => domain%iproc

    buf_send        => comm%buf_send
    buf_recv        => comm%buf_recv
    int_send        => comm%int_send
    int_recv        => comm%int_recv

    rest_move       => domain%type1_comm
    fit_move        => domain%type2_comm
    buf_comm_real   => domain%buf_var0_comm_real
    buf_comm_int    => domain%buf_var0_comm_int

    k = num
    allocate(irequest1(k), irequest2(k), irequest3(k), irequest4(k), &
             send_size(2,k+1), recv_size(2,k+1), send1(2,k+1), recv1(2,k+1))

    ! Pack outgoing data
    !
    send1(1:2,1) = 0
    do iproc = 1, num
      start_iproc(1:2) = send1(1:2,iproc)
      int_send(1+start_iproc(1)) = rest_move(iproc)
      int_send(2+start_iproc(1)) = fit_move(iproc)
      send_size(1,iproc) = 2 + rest_move(iproc) + fit_move(iproc)
      send_size(2,iproc) = 7*rest_move(iproc) + 3*fit_move(iproc)
      send1(1:2,iproc+1) = send1(1:2,iproc) + send_size(1:2,iproc)
    end do

    rest_move(1:num) = 0
    fit_move(1:num) = 0
    do i = 1, enefunc%rest_comm_domain
      iproc = buf_comm_int(2*i-1)
      start_iproc(1:2) = send1(1:2,iproc)
      rest_move(iproc) = rest_move(iproc) + 1
      k1 = start_iproc(1) + 2 + rest_move(iproc)
      int_send(k1) = buf_comm_int(2*i)
      j1 = start_iproc(2) + 7*(rest_move(iproc)-1)
      buf_send(j1+1:j1+7) = buf_comm_real(7*i-6:7*i)
    end do
    nadd1 = 2*enefunc%rest_comm_domain
    nadd2 = 7*enefunc%rest_comm_domain
    do i = 1, enefunc%fit_comm_domain
      iproc = buf_comm_int(2*i-1+nadd1)
      start_iproc(1:2) = send1(1:2,iproc)
      fit_move(iproc) = fit_move(iproc) + 1
      k1 = start_iproc(1) + 2 + rest_move(iproc) + fit_move(iproc)
      int_send(k1) = buf_comm_int(2*i+nadd1)
      j1 = start_iproc(2) + 7*rest_move(iproc) + 3*(fit_move(iproc)-1)
      buf_send(j1+1:j1+3) = buf_comm_real(3*i-2+nadd2:3*i+nadd2)
    end do

#ifdef HAVE_MPI_GENESIS
    ! communication of the size of data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(recv_size(1,iproc), 2, mpi_integer, ip,           &
                     (my_city_rank+1)*nproc_city+ip, mpi_comm_city,    &
                     irequest1(iproc), ierror)
      call mpi_isend(send_size(1,iproc), 2, mpi_integer, ip,           &
                     (ip+1)*nproc_city+my_city_rank, mpi_comm_city,    &
                     irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    recv1(1:2,1) = 0
    start_iproc(1:2) = 0
    do iproc = 2, num
      start_iproc(1:2) = start_iproc(1:2) + recv_size(1:2,iproc-1)
      recv1(1:2,iproc) = start_iproc(1:2)
    end do

    ! communication of the integer data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(int_recv(recv1(1,iproc)+1), recv_size(1,iproc),   &
                     mpi_integer, ip, (my_city_rank+1)*nproc_city+ip,  &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(int_send(send1(1,iproc)+1), send_size(1,iproc),   &
                     mpi_integer, ip, (ip+1)*nproc_city+my_city_rank,  &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do

    ! communication of the real data
    !
    do iproc = 1, num
      ip = iproc_proc(iproc)
      call mpi_irecv(buf_recv(recv1(2,iproc)+1), recv_size(2,iproc),   &
                     mpi_wip_real, ip, (my_city_rank+1)*nproc_city+ip, &
                     mpi_comm_city, irequest1(iproc), ierror)
      call mpi_isend(buf_send(send1(2,iproc)+1), send_size(2,iproc),   &
                     mpi_wip_real, ip, (ip+1)*nproc_city+my_city_rank, &
                     mpi_comm_city, irequest2(iproc), ierror)
    end do
    do iproc = 1, num
      call mpi_wait(irequest1(iproc), istatus, ierror)
      call mpi_wait(irequest2(iproc), istatus, ierror)
    end do
#endif

    ! get the imcoming data
    !
    rest_move(1:num) = 0
    fit_move (1:num) = 0

    list = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      rest_move(iproc) = int_recv(1+start_iproc(1))
      fit_move (iproc) = int_recv(2+start_iproc(1))
      do i = 1, rest_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + i
        buf_comm_int (list) = int_recv(k1)
        j1 = start_iproc(2) + 7*(i-1)
        buf_comm_real(7*list-6:7*list) = buf_recv(j1+1:j1+7)
      end do
    end do
    enefunc%rest_comm_domain = list

    nadd1 = list
    nadd2 = list*7
    list  = 0
    do iproc = 1, num
      start_iproc(1:2) = recv1(1:2,iproc)
      do i = 1, fit_move(iproc)
        list = list + 1
        k1 = start_iproc(1) + 2 + rest_move(iproc) + i
        buf_comm_int (list+nadd1) = int_recv(k1)
        j1 = start_iproc(2) + 7*rest_move(iproc) + 3*(i-1)
        buf_comm_real(3*list-2+nadd2:3*list+nadd2) = buf_recv(j1+1:j1+3)
      end do
    end do
    enefunc%fit_comm_domain = list

    deallocate(irequest1, irequest2, irequest3, irequest4, send_size, &
               recv_size, send1, recv1)
    return

  end subroutine communicate_restraint

 
  subroutine setup_cellc_face_alloc(iproc, cell, cell_to_rank,   &
                                    jcell_start, jcell_end,      &
                                    num_cellc)

    ! formal arguments
    integer,                  intent(in)    :: iproc
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(in)    :: jcell_start(:), jcell_end(:)
    integer,                  intent(inout) :: num_cellc(:)

    integer                   :: a, ix, ixx, iy, iyy, iz, k

    ! x face
    !
    a = jcell_end(1) + 1
    if (a == cell(1)+1) a = 1
    do iz = jcell_start(3), jcell_end(3)
      do iy = jcell_start(2), jcell_end(2)
        if (my_city_rank == cell_to_rank(a,iy,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
        end if
      end do
    end do

    a = jcell_start(1) - 1
    if (a == 0) a = cell(1)
    do iz = jcell_start(3), jcell_end(3)
      do iy = jcell_start(2), jcell_end(2)
        if (my_city_rank == cell_to_rank(a,iy,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
        end if
      end do
    end do

    ! y face
    !
    a = jcell_end(2) + 1
    if (a == cell(2)+1) a = 1
    do iz = jcell_start(3), jcell_end(3)
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,a,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
        end if
      end do
    end do

    a = jcell_start(2) - 1
    if (a == 0) a = cell(2)
    do iz = jcell_start(3), jcell_end(3)
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,a,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
        end if
      end do
    end do

    ! z face
    !
    a = jcell_end(3) + 1
    if (a == cell(3)+1) a = 1
    do iy = jcell_start(2)-1, jcell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,iyy,a)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
        end if
      end do
    end do

    ! z face
    a = jcell_start(3) - 1
    if (a == 0) a = cell(3)
    do iy = jcell_start(2)-1, jcell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,iyy,a)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
        end if
      end do
    end do

    return

  end subroutine setup_cellc_face_alloc


  subroutine setup_cellf_face_alloc(cell, cell_to_rank,          &
                                    icell_start, icell_end,      &
                                    rank_g2l, num_cellf)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(in)    :: icell_start(:), icell_end(:)
    integer,                  intent(in)    :: rank_g2l(:)
    integer,                  intent(inout) :: num_cellf(:)

    integer                   :: a, ix, ixx, iy, iyy, iz, iproc, ip

    ! x face
    !
    a = icell_start(1) - 1
    if (a == 0) a = cell(1)
    do iz = icell_start(3), icell_end(3)
      do iy = icell_start(2), icell_end(2)
        iproc = cell_to_rank(a,iy,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    a = icell_end(1) + 1
    if (a == cell(1)+1) a = 1
    do iz = icell_start(3), icell_end(3)
      do iy = icell_start(2), icell_end(2)
        iproc = cell_to_rank(a,iy,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    ! y face
    !
    a = icell_start(2) - 1
    if (a == 0)  a = cell(2)
    do iz = icell_start(3), icell_end(3)
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,a,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    a = icell_end(2) + 1
    if (a == cell(2)+1)  a = 1
    do iz = icell_start(3), icell_end(3)
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,a,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    ! z face
    !
    a = icell_start(3) - 1
    if (a == 0) a = cell(3)
    do iy = icell_start(2)-1, icell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,iyy,a)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    a = icell_end(3) + 1
    if (a == cell(3)+1) a = 1
    do iy = icell_start(2)-1, icell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,iyy,a)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
      end do
    end do

    return

   end subroutine setup_cellf_face_alloc
        
  subroutine setup_cellc_face(iproc, cell, cell_to_rank,   &
                              jcell_start, jcell_end,      &
                              cell_g2l, num_cellc,         &
                              ic_send, if_recv)

    ! formal arguments
    integer,                  intent(in)    :: iproc
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(in)    :: jcell_start(:), jcell_end(:)
    integer,                  intent(in)    :: cell_g2l(:)
    integer,                  intent(inout) :: num_cellc(:)
    integer,                  intent(inout) :: ic_send(:,:)
    integer,                  intent(inout) :: if_recv(:,:)

    integer                   :: a, ix, ixx, iy, iyy, iz, k, ic_global, ic_local

    ! x face
    !
    a = jcell_start(1) - 1
    if (a == 0) a = cell(1)
    do iz = jcell_start(3), jcell_end(3)
      do iy = jcell_start(2), jcell_end(2)
        if (my_city_rank == cell_to_rank(a,iy,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = a + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    a = jcell_end(1) + 1
    if (a == cell(1)+1) a = 1
    do iz = jcell_start(3), jcell_end(3)
      do iy = jcell_start(2), jcell_end(2)
        if (my_city_rank == cell_to_rank(a,iy,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = a + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    ! y face
    !
    a = jcell_start(2) - 1
    if (a == 0) a = cell(2)
    do iz = jcell_start(3), jcell_end(3)
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,a,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = ixx + (a-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    a = jcell_end(2) + 1
    if (a == cell(2)+1) a = 1
    do iz = jcell_start(3), jcell_end(3)
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,a,iz)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = ixx + (a-1)*cell(1) + (iz-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    ! z face
    !
    a = jcell_start(3) - 1
    if (a == 0) a = cell(3)
    do iy = jcell_start(2)-1, jcell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,iyy,a)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = ixx + (iyy-1)*cell(1) + (a-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    a = jcell_end(3) + 1
    if (a == cell(3)+1) a = 1
    do iy = jcell_start(2)-1, jcell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = jcell_start(1)-1, jcell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        if (my_city_rank == cell_to_rank(ixx,iyy,a)) then
          num_cellc(iproc) = num_cellc(iproc) + 1
          k = num_cellc(iproc)
          ic_global = ixx + (iyy-1)*cell(1) + (a-1)*cell(1)*cell(2)
          ic_local  = cell_g2l(ic_global)
          ic_send(k,iproc) = ic_local
          if_recv(k,iproc) = ic_local
        end if
      end do
    end do

    return

  end subroutine setup_cellc_face

  subroutine setup_cellf_face(num_cell_local,                &
                              cell, cell_to_rank,            &
                              icell_start, icell_end,        &
                              cell_g2b, rank_g2l, num_cellf, &
                              ic_recv, if_send)

    ! formal arguments
    integer,                  intent(in)    :: num_cell_local
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(in)    :: icell_start(:), icell_end(:)
    integer,                  intent(in)    :: cell_g2b(:), rank_g2l(:)
    integer,                  intent(inout) :: num_cellf(:)
    integer,                  intent(inout) :: ic_recv(:,:)
    integer,                  intent(inout) :: if_send(:,:)

    integer                   :: a, ix, ixx, iy, iyy, iz, k, iproc, ip
    integer                   :: ic_global, ic_local

    ! x face
    !
    a = icell_start(1) - 1
    if (a == 0) a = cell(1)
    do iz = icell_start(3), icell_end(3)
      do iy = icell_start(2), icell_end(2)
        iproc = cell_to_rank(a,iy,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = a + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    a = icell_end(1) + 1
    if (a == cell(1)+1) a = 1
    do iz = icell_start(3), icell_end(3)
      do iy = icell_start(2), icell_end(2)
        iproc = cell_to_rank(a,iy,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = a + (iy-1)*cell(1) + (iz-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    ! y face
    !
    a = icell_start(2) - 1
    if (a == 0)  a = cell(2)
    do iz = icell_start(3), icell_end(3)
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,a,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = ixx + (a-1)*cell(1) + (iz-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    a = icell_end(2) + 1
    if (a == cell(2)+1)  a = 1
    do iz = icell_start(3), icell_end(3)
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,a,iz)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = ixx + (a-1)*cell(1) + (iz-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    ! z face
    !
    a = icell_start(3) - 1
    if (a == 0) a = cell(3)
    do iy = icell_start(2)-1, icell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,iyy,a)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = ixx + (iyy-1)*cell(1) + (a-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    a = icell_end(3) + 1
    if (a == cell(3)+1) a = 1
    do iy = icell_start(2)-1, icell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = icell_start(1)-1, icell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        iproc = cell_to_rank(ixx,iyy,a)
        ip    = rank_g2l(iproc+1)
        num_cellf(ip) = num_cellf(ip) + 1
        k = num_cellf(ip)
        ic_global = ixx + (iyy-1)*cell(1) + (a-1)*cell(1)*cell(2)
        ic_local  = cell_g2b(ic_global) + num_cell_local
        ic_recv(k,ip) = ic_local
        if_send(k,ip) = ic_local
      end do
    end do

    return

   end subroutine setup_cellf_face
        
end module cg_communicate_mod
