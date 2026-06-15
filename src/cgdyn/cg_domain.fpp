!--------1---------2---------3---------4---------5---------6---------7---------8 ! 
!  Module   cg_domain_mod
!> @brief   Provides utilities for domain decomposition.                 
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_domain_mod

  use cg_boundary_mod
  use cg_energy_mod
  use cg_communicate_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use molecules_str_mod
  use molecules_mod
  use timers_mod
  use random_mod
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

  ! parameters for setup cell capacity
  real(wip),        private, parameter :: VolumeBox8      = 5120.0_wip
  real(wip),        private, parameter :: VboxRate        = 1.4_wip
  real(wip),        private, parameter :: ShrinkRate      = 1.2_wip
  integer,          private, parameter :: NAtomBox8       =  90
  integer,          private, parameter :: NBondBox8       =  90
  integer,          private, parameter :: NAnglBox8       = 100
  integer,          private, parameter :: NDiheBox8       = 200
  integer,          private, parameter :: NImprBox8       = 50
  integer,          private, parameter :: NHGrpBox8       = 40
  integer,          private, parameter :: NHMovBox8       = 8
  integer,          private, parameter :: NContBox8       = 100
#ifdef USE_GPU
  integer,          private, parameter :: NAtomMax_in_CUDA = 255
#endif

  ! subroutines
  public  :: setup_domain
  public  :: setup_domain_interaction
  public  :: setup_cell_to_rank
  public  :: setup_cell_boundary
  public  :: update_nobc_boundary
  public  :: copy_domain_information
  public  :: molecule_accum
  public  :: setup_domain_lb
  private :: setup_atom_pbc
  private :: setup_atom_nobc
  private :: assign_neighbor_cells
  private :: assign_cell_atoms
  private :: assign_cell_interaction
  private :: molecule_to_domain
  private :: check_atom_coord

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain
  !> @brief        Sets up domain information.
  !! @authors      JJ  
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain(ene_info, boundary, molecule, enefunc, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: cell(3), nc(3)
    integer                  :: ncel_local, ncel_bound, ncel_all

    ! initialize structure informations
    !
    call init_domain(domain)
    call init_enefunc(enefunc)

    domain%num_atom_all   = molecule%num_atoms
    domain%system_size(1) = real(boundary%box_size_x,wp)
    domain%system_size(2) = real(boundary%box_size_y,wp)
    domain%system_size(3) = real(boundary%box_size_z,wp)
    domain%cell_size(1)   = real(boundary%cell_size_x,wp)
    domain%cell_size(2)   = real(boundary%cell_size_y,wp)
    domain%cell_size(3)   = real(boundary%cell_size_z,wp)

    cell(1) = boundary%num_cells_x
    cell(2) = boundary%num_cells_y
    cell(3) = boundary%num_cells_z
    nc(1:3) = boundary%num_domain(1:3)
    domain%num_cell(1:3) = cell(1:3)

    call alloc_domain(domain, DomainCellGlobal, cell(1), cell(2), cell(3))

    ! decision of domains
    !
    if (boundary%load_balance) then
      call setup_domain_range(cell, molecule, boundary, domain)
    else
      call setup_domain_range_nolb(cell, boundary, domain)
    end if

    ! assign the rank of each dimension from my_rank
    !
    call setup_cell_to_rank(cell, boundary, domain)

    ! check the MPI rank that should be communicated
    !
    call alloc_domain(domain, DomainNeighborRank, nproc_country, 1, 1)
    call setup_comm_rank(cell, domain)

    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local      = domain%num_cell_local
    ncel_bound      = domain%num_cell_boundary
    ncel_all        = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)

    ! assign global<->local mapping of cell index
    !
    call setup_cell_local(cell, domain)

    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, nc, domain)

    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all,            1, 1)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)

    if (boundary%type == BoundaryTypePBC) then
      call setup_atom_pbc(molecule, boundary, enefunc, domain)
    else
      call setup_atom_nobc(molecule, boundary, enefunc, domain)
    end if

    ! assign the interaction cell for each interaction
    !
    enefunc%forcefield = ene_info%forcefield
    call assign_neighbor_cells(enefunc, boundary, domain)

    call setup_domain_interaction(boundary, domain)

    call check_atom_coord(boundary, domain)

    return

  end subroutine setup_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_range
  !> @brief        Determines the domain range.
  !! @authors      JJ  
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_range(cell, molecule, boundary, domain)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3), origin(3)
    integer                   :: i, icx, icy, icz
    integer                   :: nproc_xyz, nproc_x, nproc_y, nproc_z
    integer                   :: nc(3), iproc(3), dime
    integer                   :: my_xx_rank, my_yy_rank, my_zz_rank
    integer                   :: index_x, index_y, index_z
    integer                   :: natom_all

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wp),         pointer :: coord(:,:)
    integer,          pointer :: min_cell
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: natom_global(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)

    coord           => molecule%atom_coord

    bsize_x         => boundary%box_size_x
    bsize_y         => boundary%box_size_y
    bsize_z         => boundary%box_size_z
    csize_x         => boundary%cell_size_x
    csize_y         => boundary%cell_size_y
    csize_z         => boundary%cell_size_z
    min_cell        => boundary%min_domain_cell
    ncel_x          => boundary%num_cells_x
    ncel_y          => boundary%num_cells_y
    ncel_z          => boundary%num_cells_z
 
    natom_global    => domain%natom_global
    cell_start      => domain%cell_start
    cell_end        => domain%cell_end  
 
    origin(1)       = boundary%origin_x
    origin(2)       = boundary%origin_y
    origin(3)       = boundary%origin_z
    nc(1)           = boundarY%num_domain(1)
    nc(2)           = boundarY%num_domain(2)
    nc(3)           = boundarY%num_domain(3)

    natom_all       = domain%num_atom_all

    ! count the number of particles in each cell
    !
    do i = 1, natom_all

      x_shift = coord(1,i) - boundary%origin_x
      y_shift = coord(2,i) - boundary%origin_y
      z_shift = coord(3,i) - boundary%origin_z

      move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
      move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
      move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

      x_shift = x_shift + move(1)
      y_shift = y_shift + move(2)
      z_shift = z_shift + move(3)

      icx = int(x_shift/csize_x)
      icy = int(y_shift/csize_y)
      icz = int(z_shift/csize_z)
      if (icx == ncel_x) icx = icx - 1
      if (icy == ncel_y) icy = icy - 1
      if (icz == ncel_z) icz = icz - 1
      icx = icx + 1
      icy = icy + 1
      icz = icz + 1

      natom_global(icx,icy,icz) = natom_global(icx,icy,icz) + 1

    end do

    nproc_xyz = nc(1) * nc(2) * nc(3)
    nproc_x   = nc(1)
    nproc_y   = nc(2)
    nproc_z   = nc(3)
    cell_start(1) = 1
    cell_start(2) = 1
    cell_start(3) = 1
    cell_end(1)   = boundary%num_cells_x
    cell_end(2)   = boundary%num_cells_y
    cell_end(3)   = boundary%num_cells_z

    ! define mpi_rank in each dimension
    !
    iproc(1) = mod(my_country_rank, nc(1))
    iproc(2) = mod(my_country_rank/nc(1),nc(2))
    iproc(3) = my_country_rank/(nc(1)*nc(2))
    index_x  = iproc(2)*nproc_country + iproc(3)
    index_y  = iproc(3)*nproc_country + iproc(1)
    index_z  = iproc(1)*nproc_country + iproc(2)
    call mpi_comm_split(mpi_comm_country, index_x, my_country_rank, &
                        grid_commx,ierror)
    call mpi_comm_rank (grid_commx, my_x_rank, ierror)
    call mpi_comm_split(mpi_comm_country, index_y, my_country_rank, &
                        grid_commy,ierror)
    call mpi_comm_rank (grid_commy, my_y_rank, ierror)
    call mpi_comm_split(mpi_comm_country, index_z, my_country_rank, &
                        grid_commz,ierror)
    call mpi_comm_rank (grid_commz, my_z_rank, ierror)
    my_xx_rank = my_x_rank
    my_yy_rank = my_y_rank
    my_zz_rank = my_z_rank

    i = 0

    do while (nproc_xyz > 1)

      i = i + 1
      ! decide the dimension that should be divided
      !
      call decide_separate_dim(nproc_x, nproc_y, nproc_z, dime)

      ! decide the boundary
      !
      if (dime == 1) call select_boundary_x(cell, natom_global, min_cell, nproc_x, &
                              nproc_xyz, my_xx_rank, cell_start, cell_end)
      if (dime == 2) call select_boundary_y(cell, natom_global, min_cell, nproc_y, &
                              nproc_xyz, my_yy_rank, cell_start, cell_end)
      if (dime == 3) call select_boundary_z(cell, natom_global, min_cell, nproc_z, &
                              nproc_xyz, my_zz_rank, cell_start, cell_end)

    end do

    return

  end subroutine setup_domain_range

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_range_nolb
  !> @brief        Determines the domain range without load balancing.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_range_nolb(cell, boundary, domain)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: nc(3), iproc(3)
    integer                   :: i, j, k, quotient, remainder

    integer,          pointer :: cell_start(:), cell_end(:)

    cell_start  => domain%cell_start
    cell_end    => domain%cell_end  
    nc(1)       = boundarY%num_domain(1)
    nc(2)       = boundarY%num_domain(2)
    nc(3)       = boundarY%num_domain(3)

    ! Assign the rank of each dimension from my_rank
    !
    iproc(1) = mod(my_city_rank, nc(1))
    iproc(2) = mod(my_city_rank/nc(1), nc(2))
    iproc(3) = my_city_rank/(nc(1)*nc(2))

    ! Assign the cell index for each processor and the total number of cells
    !
    do i = 1, 3

      quotient = cell(i) / nc(i)
      remainder = mod(cell(i), nc(i))

      if (iproc(i) <= (remainder -1)) then
        quotient       = quotient + 1
        cell_start(i)  = quotient * iproc(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
      else
        cell_start(i)  = (quotient+1)*remainder + quotient*(iproc(i)-remainder)
        cell_start(i)  = cell_start(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
      end if

    end do

    return

  end subroutine setup_domain_range_nolb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    decide_separate_dim
  !> @brief        Determines the dimension to be separated in a KD-tree.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine decide_separate_dim(nproc_x, nproc_y, nproc_z, dime)

    ! formal arguments
    integer,                  intent(in)    :: nproc_x
    integer,                  intent(in)    :: nproc_y
    integer,                  intent(in)    :: nproc_z
    integer,                  intent(inout) :: dime

    integer                   :: dime_max

    if (nproc_x > 1 .and. nproc_y > 1 .and. nproc_z > 1) then

      dime_max = max(nproc_x, nproc_y, nproc_z)
      if (nproc_x == dime_max) then
        dime = 1
      else if (nproc_y == dime_max) then
        dime = 2
      else
        dime = 3
      end if

    else if (nproc_x > 1 .and. nproc_y > 1) then

      dime_max = max(nproc_x, nproc_y)
      if (nproc_x == dime_max) then
        dime = 1
      else
        dime = 2
      end if

    else if (nproc_x > 1 .and. nproc_z > 1) then

      dime_max = max(nproc_x, nproc_z)
      if (nproc_x == dime_max) then
        dime = 1
      else
        dime = 3
      end if

    else if (nproc_y > 1 .and. nproc_z > 1) then

      dime_max = max(nproc_y, nproc_z)
      if (nproc_y == dime_max) then
        dime = 2
      else
        dime = 3
      end if

    else if (nproc_x > 1) then

      dime = 1

    else if (nproc_y > 1) then

      dime = 2

    else if (nproc_z > 1) then

      dime = 3

    end if

    return

  end subroutine decide_separate_dim

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_boundary_x
  !> @brief        Determines the boundary that separates domains in the X-dimension.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_boundary_x(nc, natom_global, min_cell, nproc_x, nproc_xyz, &
                               my_xx_rank, cell_start, cell_end)

    integer,                  intent(in   ) :: nc(:)
    integer,                  intent(in   ) :: natom_global(:,:,:)
    integer,                  intent(in   ) :: min_cell
    integer,                  intent(inout) :: nproc_x, nproc_xyz
    integer,                  intent(inout) :: my_xx_rank
    integer,                  intent(inout) :: cell_start(:)
    integer,                  intent(inout) :: cell_end(:)

    integer                   :: sum_total, sum_half, i, j, k, ia, ib
    integer                   :: alloc_stat, half_point
    integer, allocatable      :: sum_natom(:), sum_accum(:)

    allocate(sum_natom(nc(1)), sum_accum(nc(1)), &
             stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    sum_total = 0
    do i = cell_start(1), cell_end(1)
      sum_natom(i) = 0
      do j = cell_start(3), cell_end(3)
        do k = cell_start(2), cell_end(2)
          sum_natom(i) = sum_natom(i) + natom_global(i,k,j)
        end do
      end do
      sum_total = sum_total + sum_natom(i)
      sum_accum(i) = sum_total
    end do

    sum_half = sum_total / 2

    if (min_cell == 2) then
      ia = cell_start(1) + nproc_x - 1
      ib = cell_end(1) - nproc_x 
    else if (min_cell == 1) then
      ia = cell_start(1) + nproc_x/2 - 1
      ib = cell_end(1) - nproc_x/2 
    end if
    if (sum_accum(ia) >= sum_half) then
      half_point = ia
    else if (sum_accum(ib) <= sum_half) then
      half_point = ib
    else
      do i = ia, ib-1
        if (sum_accum(i) < sum_half .and. sum_accum(i+1) >= sum_half) &
          half_point = i
      end do
    end if

    if (my_xx_rank < nproc_x/2) then
      cell_end(1) = half_point
    else
      cell_start(1) = half_point + 1
      my_xx_rank = my_xx_rank - nproc_x/2
    end if
    nproc_x = nproc_x / 2
    nproc_xyz = nproc_xyz / 2

    deallocate(sum_natom, sum_accum, stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine select_boundary_x

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_boundary_y
  !> @brief        Determines the boundary that separates domains in the Y-dimension.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_boundary_y(nc, natom_global, min_cell, nproc_y, nproc_xyz, &
                               my_yy_rank, cell_start, cell_end)

    integer,                  intent(in   ) :: nc(:)
    integer,                  intent(in   ) :: natom_global(:,:,:)
    integer,                  intent(in   ) :: min_cell
    integer,                  intent(inout) :: nproc_y, nproc_xyz
    integer,                  intent(inout) :: my_yy_rank
    integer,                  intent(inout) :: cell_start(:)
    integer,                  intent(inout) :: cell_end(:)

    integer                   :: sum_total, sum_half, i, j, k, ia, ib
    integer                   :: alloc_stat, half_point
    integer, allocatable      :: sum_natom(:), sum_accum(:)

    allocate(sum_natom(nc(2)), sum_accum(nc(2)), &
             stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    sum_total = 0
    do i = cell_start(2), cell_end(2)
      sum_natom(i) = 0
      do j = cell_start(3), cell_end(3)
        do k = cell_start(1), cell_end(1)
          sum_natom(i) = sum_natom(i) + natom_global(k,i,j)
        end do
      end do
      sum_total = sum_total + sum_natom(i)
      sum_accum(i) = sum_total
    end do

    sum_half = sum_total / 2

    if (min_cell == 2) then
      ia = cell_start(2) + nproc_y - 1
      ib = cell_end(2) - nproc_y  
    else if (min_cell == 1) then
      ia = cell_start(2) + nproc_y/2 - 1
      ib = cell_end(2) - nproc_y/2
    end if
    if (sum_accum(ia) >= sum_half) then
      half_point = ia
    else if (sum_accum(ib) <= sum_half) then
      half_point = ib
    else
      do i = ia, ib-1
        if (sum_accum(i) < sum_half .and. sum_accum(i+1) >= sum_half) &
          half_point = i
      end do
    end if

    if (my_yy_rank < nproc_y/2) then
      cell_end(2) = half_point
    else
      cell_start(2) = half_point + 1
      my_yy_rank = my_yy_rank - nproc_y/2
    end if
    nproc_y = nproc_y / 2
    nproc_xyz = nproc_xyz / 2

    deallocate(sum_natom, sum_accum, stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine select_boundary_y

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_boundary_z
  !> @brief        Determines the boundary that separates domains in the Z-dimension.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_boundary_z(nc, natom_global, min_cell, nproc_z, nproc_xyz, &
                               my_zz_rank, cell_start, cell_end)

    integer,                  intent(in   ) :: nc(:)
    integer,                  intent(in   ) :: natom_global(:,:,:)
    integer,                  intent(in   ) :: min_cell
    integer,                  intent(inout) :: nproc_z, nproc_xyz
    integer,                  intent(inout) :: my_zz_rank
    integer,                  intent(inout) :: cell_start(:)
    integer,                  intent(inout) :: cell_end(:)

    integer                   :: sum_total, sum_half, i, j, k, ia, ib
    integer                   :: alloc_stat, half_point
    integer, allocatable      :: sum_natom(:), sum_accum(:)

    allocate(sum_natom(nc(3)), sum_accum(nc(3)), &
             stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    sum_total = 0
    do i = cell_start(3), cell_end(3)
      sum_natom(i) = 0
      do j = cell_start(2), cell_end(2)
        do k = cell_start(1), cell_end(1)
          sum_natom(i) = sum_natom(i) + natom_global(k,j,i)
        end do
      end do
      sum_total = sum_total + sum_natom(i)
      sum_accum(i) = sum_total
    end do

    sum_half = sum_total / 2

    if (min_cell == 2) then
      ia = cell_start(3) + nproc_z - 1
      ib = cell_end(3) - nproc_z  
    else if (min_cell == 1) then
      ia = cell_start(3) + nproc_z/2 - 1
      ib = cell_end(3) - nproc_z/2
    end if

    if (sum_accum(ia) >= sum_half) then
      half_point = ia
    else if (sum_accum(ib) <= sum_half) then
      half_point = ib
    else
      do i = ia, ib-1
        if (sum_accum(i) < sum_half .and. sum_accum(i+1) >= sum_half) &
          half_point = i
      end do
    end if

    if (my_zz_rank < nproc_z/2) then
      cell_end(3) = half_point
    else
      cell_start(3) = half_point + 1
      my_zz_rank = my_zz_rank - nproc_z/2
    end if
    nproc_z = nproc_z / 2
    nproc_xyz = nproc_xyz / 2

    deallocate(sum_natom, sum_accum, stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine select_boundary_z

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_local
  !> @brief        Defines the relationship between local and global cell indices.
  !! @authors      JJ
  !! @param[in]    cell     : number of cell
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_local(cell, domain)

    integer,                  intent(in)    :: cell(:)
    type(s_domain),   target, intent(inout) :: domain

    integer                   :: i, j, k, icel_local, icel
    integer,          pointer :: cell_start(:), cell_end(:)
    integer,          pointer :: cell_g2l(:), cell_l2g(:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx_orig(:), cell_l2gy_orig(:)
    integer,          pointer :: cell_l2gz_orig(:), cell_gxyz2l(:,:,:)

    cell_start        => domain%cell_start
    cell_end          => domain%cell_end  
    cell_g2l          => domain%cell_g2l
    cell_l2g          => domain%cell_l2g
    cell_l2gx         => domain%cell_l2gx
    cell_l2gy         => domain%cell_l2gy
    cell_l2gz         => domain%cell_l2gz
    cell_l2gx_orig    => domain%cell_l2gx_orig
    cell_l2gy_orig    => domain%cell_l2gy_orig
    cell_l2gz_orig    => domain%cell_l2gz_orig
    cell_gxyz2l       => domain%cell_gxyz2l    

    icel_local = 0
    do i = cell_start(3), cell_end(3)
      do j = cell_start(2), cell_end(2)
        do k = cell_start(1), cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          cell_g2l(icel) = icel_local
          cell_l2g(icel_local) = icel
          cell_l2gx(icel_local) = k
          cell_l2gy(icel_local) = j
          cell_l2gz(icel_local) = i
          cell_l2gx_orig(icel_local) = k
          cell_l2gy_orig(icel_local) = j
          cell_l2gz_orig(icel_local) = i
          cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do

    return

  end subroutine setup_cell_local

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_interaction
  !> @brief        Defines the pairwise interaction between cells.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_interaction(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, j, ij, inbc
    integer                   :: ncel_local, nboundary, cell(3)

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wp),         pointer :: natom(:)
    integer,          pointer :: cell_start(:), cell_end(:)
    integer,          pointer :: cell_gxyz2l(:,:,:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx1(:), cell_l2gy1(:), cell_l2gz1(:)
    integer,          pointer :: cell_g2l(:)
    integer,          pointer :: num_domain(:)
    integer,          pointer :: natom_t0(:), ncharge_t0(:)
    integer,          pointer :: cell_to_rank(:,:,:)
    integer,          pointer :: near_cells_count(:), near_cells(:,:)
    integer,          pointer :: far_cells_count(:), far_cells(:,:)
  
    real(wp),         allocatable :: ncharge(:)

    bsize_x           => boundary%box_size_x
    bsize_y           => boundary%box_size_y
    bsize_z           => boundary%box_size_z
    num_domain        => boundary%num_domain

    cell_start        => domain%cell_start
    cell_end          => domain%cell_end
    cell_gxyz2l       => domain%cell_gxyz2l
    cell_l2gx         => domain%cell_l2gx
    cell_l2gy         => domain%cell_l2gy
    cell_l2gz         => domain%cell_l2gz
    cell_l2gx1        => domain%cell_l2gx_orig
    cell_l2gy1        => domain%cell_l2gy_orig
    cell_l2gz1        => domain%cell_l2gz_orig
    cell_g2l          => domain%cell_g2l
    natom_t0          => domain%num_atom
    ncharge_t0        => domain%num_charge
    cell_to_rank      => domain%natom_global
    near_cells_count  => domain%near_cells_count
    near_cells        => domain%near_cells
    far_cells_count   => domain%far_cells_count
    far_cells         => domain%far_cells
    natom             => domain%num_atom_t0

    cell(1)      = boundary%num_cells_x
    cell(2)      = boundary%num_cells_y
    cell(3)      = boundary%num_cells_z
    ncel_local   = domain%num_cell_local
    nboundary    = domain%num_cell_boundary

    allocate(ncharge(1:ncel_local+nboundary))

    ! number of atoms in each cell (proceesor number is also considered)
    !
    call assign_cell_atoms(natom, ncharge, natom_t0, ncharge_t0,  &
                           cell_l2gx1, cell_l2gy1, cell_l2gz1,    &
                           cell_to_rank, ncel_local, nboundary)

    ! assign the interaction cell for each interaction
    !
    call assign_cell_interaction(natom, ncharge, cell,                  &
                                 cell_l2gx,  cell_l2gy,  cell_l2gz,     &
                                 cell_gxyz2l,                &
                                 ncel_local, nboundary,                 &
                                 near_cells_count, near_cells,          &
                                 far_cells_count, far_cells)


    if (boundary%type == BoundaryTypeNOBC) then

      do i = 1, ncel_local+nboundary
        ij = 0
        do inbc = 1, domain%near_cells_count(i)
          j = domain%near_cells(inbc,i)
          if (abs(cell_l2gx(i)-cell_l2gx(j)) < cell(1)/2 .and. &
              abs(cell_l2gy(i)-cell_l2gy(j)) < cell(2)/2 .and. &
              abs(cell_l2gz(i)-cell_l2gz(j)) < cell(3)/2) then
            ij = ij + 1
            domain%near_cells(ij,i) = j
          end if
        end do
        domain%near_cells_count(i) = ij
        ij = 0
        do inbc = 1, domain%far_cells_count(i)
          j = domain%far_cells(inbc,i)
          if (abs(cell_l2gx(i)-cell_l2gx(j)) < cell(1)/2 .and. &
              abs(cell_l2gy(i)-cell_l2gy(j)) < cell(2)/2 .and. &
              abs(cell_l2gz(i)-cell_l2gz(j)) < cell(3)/2) then
            ij = ij + 1
            domain%far_cells(ij,i) = j
          end if
        end do
        domain%far_cells_count(i) = ij
      end do

    end if
    
    deallocate(ncharge)

    return

  end subroutine setup_domain_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_to_rank
  !> @brief        Defines the processor rank and the number of cells &
  !!               in each domain.
  !! @authors      JJ
  !! @param[in]    cell     : cell number in each dimension
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_to_rank(cell, boundary, domain)

    ! formal arguments
    integer,                 intent(in)    :: cell(3)
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k

    integer,         pointer :: ncell_local, ncell_boundary, num_domain(:)
    integer,         pointer :: cell_start(:), cell_end(:), cell_length(:)
    integer,         pointer :: cell_to_rank(:,:,:)

    num_domain          => boundary%num_domain

    ncell_local         => domain%num_cell_local
    ncell_boundary      => domain%num_cell_boundary
    cell_start          => domain%cell_start
    cell_end            => domain%cell_end  
    cell_length         => domain%cell_length
    cell_to_rank        => domain%natom_global
    cell_start          => domain%cell_start
    cell_end            => domain%cell_end  

    cell_to_rank(1:cell(1),1:cell(2),1:cell(3)) = 0
      
    do k = cell_start(3), cell_end(3)
      do j = cell_start(2), cell_end(2)
        do i = cell_start(1), cell_end(1)
          cell_to_rank(i,j,k) = my_country_rank
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, cell_to_rank, cell(1)*cell(2)*cell(3), &
                       mpi_integer, mpi_sum, mpi_comm_country, ierror)
  
    ! Assign the cell length in each dimension
    !
    cell_length(1:3) = cell_end(1:3) - cell_start(1:3) + 1

    ! define the number of cells in each domain
    !
    ncell_local = cell_length(1)*cell_length(2)*cell_length(3)
   
    ! maximum number of boundary cells
    !
    ncell_boundary = 2*(cell_length(1) * cell_length(2)                     &
                       +cell_length(2) * cell_length(3)                     &
                       +cell_length(1) * cell_length(3))                    &
                   + 4*(cell_length(1) + cell_length(2) + cell_length(3))   &
                   + 8

    return

  end subroutine setup_cell_to_rank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_comm_rank
  !> @brief        Defines the processor rank that should be communicated.
  !! @authors      JJ
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_comm_rank(cell, domain)

    ! formal arguments
    integer,         target, intent(in   ) :: cell(:)
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: num_proc

    integer,         pointer :: cell_start(:), cell_end(:)
    integer,         pointer :: cell_to_rank(:,:,:)
    integer,         pointer :: iproc(:)
    integer,         pointer :: cell_start_proc(:,:), cell_end_proc(:,:)

    cell_start      => domain%cell_start
    cell_end        => domain%cell_end  
    cell_to_rank    => domain%natom_global
    iproc           => domain%iproc      
    cell_start_proc => domain%cell_start_proc
    cell_end_proc   => domain%cell_end_proc

    num_proc = 0

    !cell_start and cell_end of all processor
    !
    cell_start_proc(my_country_rank+1,1:3) = cell_start(1:3)
    cell_end_proc  (my_country_rank+1,1:3) = cell_end  (1:3)

    call mpi_allreduce(mpi_in_place, cell_start_proc, 3*nproc_country, &
                       mpi_integer, mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, cell_end_proc  , 3*nproc_country, &
                       mpi_integer, mpi_sum, mpi_comm_country, ierror)

    ! face
    !
    call neighbor_rank_x(cell, cell_start, cell_end, cell_to_rank,  &
                         num_proc, iproc)

    call neighbor_rank_y(cell, cell_start, cell_end, cell_to_rank,  &
                         num_proc, iproc)

    call neighbor_rank_z(cell, cell_start, cell_end, cell_to_rank,  &
                         num_proc, iproc)

    !check the maximum number of processor to be communicated
    !
    domain%num_comm_proc = num_proc

    call alloc_domain(domain, DomainComm, num_proc, 1, 1)

    return

  end subroutine setup_comm_rank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_boundary
  !> @brief        Defines boundary cells in each domain.
  !! @authors      JJ
  !! @param[in]    cell     : cell count for each axis
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_boundary(cell, ndom, domain)

    ! formal arguments
    integer,                 intent(in)    :: cell(3)
    integer,                 intent(in)    :: ndom(3)
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ic, jc, kc
    integer                  :: icel, icel_local

    integer,         pointer :: ncell, ncell_boundary
    integer,         pointer :: cell_start(:), cell_end(:)
    integer,         pointer :: cell_g2b(:), cell_b2g(:)
    integer,         pointer :: cell_gxyz2l(:,:,:)
    integer,         pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,         pointer :: cell_l2gx_orig(:)
    integer,         pointer :: cell_l2gy_orig(:)
    integer,         pointer :: cell_l2gz_orig(:)
    integer,         pointer :: cell_to_rank(:,:,:)
    integer,         pointer :: iproc(:)
    integer,         pointer :: domain_cell_rank(:)
    real(wp),        pointer :: cell_pbc_move(:,:)

    ncell                => domain%num_cell_local
    ncell_boundary       => domain%num_cell_boundary
    cell_start           => domain%cell_start
    cell_end             => domain%cell_end  
    cell_g2b             => domain%cell_g2b
    cell_b2g             => domain%cell_b2g
    cell_gxyz2l          => domain%cell_gxyz2l
    cell_l2gx            => domain%cell_l2gx
    cell_l2gy            => domain%cell_l2gy
    cell_l2gz            => domain%cell_l2gz
    cell_l2gx_orig       => domain%cell_l2gx_orig
    cell_l2gy_orig       => domain%cell_l2gy_orig
    cell_l2gz_orig       => domain%cell_l2gz_orig
    cell_to_rank         => domain%natom_global
    domain_cell_rank     => domain%domain_cell_rank
    cell_pbc_move        => domain%cell_pbc_move
    iproc                => domain%iproc

    icel_local = 0

    ! face boundary (x direction, upper)
    !
    do j = cell_start(2), cell_end(2)
      do k = cell_start(3), cell_end(3)
        icel_local = icel_local + 1
        if (cell_end(1) == cell(1)) then
          i = 1
          ic = cell(1) + 1
        else
          i = cell_end(1) + 1
          ic = i
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_b2g(icel_local) = icel
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = j
        cell_l2gz(icel_local+ncell) = k
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,j,k) = icel_local+ncell
      end do
    end do

    ! face boundary (x direction, lower)
    !
    do j = cell_start(2), cell_end(2)
      do k = cell_start(3), cell_end(3)
        icel_local = icel_local + 1
        if (cell_start(1) == 1) then
          i = cell(1)
          ic = 0
        else
          i = cell_start(1) - 1
          ic = i
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_b2g(icel_local) = icel
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = j
        cell_l2gz(icel_local+ncell) = k
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,j,k) = icel_local+ncell
      end do
    end do

    ! face boundary (y direction, upper)
    !
    do ic = cell_start(1)-1, cell_end(1)+1
      if (ic == 0) then
        i = cell(1)
      else if (ic == (cell(1)+1)) then
        i = 1
      else
        i = ic
      end if
      do k = cell_start(3), cell_end(3)
        icel_local = icel_local + 1
        if (cell_end(2) == cell(2)) then
          j = 1
          jc = cell(2) + 1
        else
          j = cell_end(2) + 1
          jc = j
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        cell_b2g(icel_local) = icel
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = jc
        cell_l2gz(icel_local+ncell) = k
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,jc,k) = icel_local+ncell
      end do
    end do

    ! face boundary (y direction, lower)
    !
    do ic = cell_start(1)-1, cell_end(1)+1
      if (ic == 0) then
        i = cell(1)
      else if (ic == (cell(1)+1)) then
        i = 1
      else
        i = ic
      end if
      do k = cell_start(3), cell_end(3)
        icel_local = icel_local + 1
        if (cell_start(2) == 1) then
          j = cell(2)
          jc = 0
        else
          j = cell_start(2) - 1
          jc = j
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        cell_b2g(icel_local) = icel
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = jc
        cell_l2gz(icel_local+ncell) = k
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,jc,k) = icel_local+ncell
      end do
    end do

    ! face boundary (z direction, upper)
    !
    do ic = cell_start(1)-1, cell_end(1)+1
      if (ic == 0) then
        i = cell(1)
      else if (ic == (cell(1)+1)) then
        i = 1
      else
        i = ic
      end if
      do jc = cell_start(2)-1, cell_end(2)+1

        if (jc == 0) then
          j = cell(2)
        else if (jc == (cell(2)+1)) then
          j = 1
        else
          j = jc
        end if
        icel_local = icel_local + 1
        if (cell_end(3) == cell(3)) then
          k = 1
          kc = cell(3) + 1
        else
          k = cell_end(3) + 1
          kc = k
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        cell_b2g(icel_local) = icel
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = jc
        cell_l2gz(icel_local+ncell) = kc
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,jc,kc) = icel_local+ncell
      end do
    end do

    ! face boundary (z direction, lower)
    !
    do ic = cell_start(1)-1, cell_end(1)+1
      if (ic == 0) then
        i = cell(1)
      else if (ic == (cell(1)+1)) then
        i = 1
      else
        i = ic
      end if
      do jc = cell_start(2)-1, cell_end(2)+1
        if (jc == 0) then
          j = cell(2)
        else if (jc == (cell(2)+1)) then
          j = 1
        else
          j = jc
        end if
        icel_local = icel_local + 1
        if (cell_start(3) == 1) then
          k = cell(3)
          kc = 0
        else
          k = cell_start(3) - 1
          kc = k
        end if
        icel = i + (j-1)*cell(1) + (k-1)*cell(1)*cell(2)
        cell_g2b(icel) = icel_local
        cell_b2g(icel_local) = icel
        domain_cell_rank(icel_local+ncell) = cell_to_rank(i,j,k)
        cell_l2gx(icel_local+ncell) = ic
        cell_l2gy(icel_local+ncell) = jc
        cell_l2gz(icel_local+ncell) = kc
        cell_l2gx_orig(icel_local+ncell) = i
        cell_l2gy_orig(icel_local+ncell) = j
        cell_l2gz_orig(icel_local+ncell) = k
        cell_gxyz2l(ic,jc,kc) = icel_local+ncell
      end do
    end do

    ! total number of boundary cells
    !
    ncell_boundary = icel_local

    do i = 1, ncell_boundary
      icel = i + ncell
      if (cell_l2gx(icel) .eq. 0) cell_pbc_move(1,icel) = -1.0_wp
      if (cell_l2gy(icel) .eq. 0) cell_pbc_move(2,icel) = -1.0_wp
      if (cell_l2gz(icel) .eq. 0) cell_pbc_move(3,icel) = -1.0_wp
      if (cell_l2gx(icel) .eq. (cell(1)+1)) cell_pbc_move(1,icel) = 1.0_wp
      if (cell_l2gy(icel) .eq. (cell(2)+1)) cell_pbc_move(2,icel) = 1.0_wp
      if (cell_l2gz(icel) .eq. (cell(3)+1)) cell_pbc_move(3,icel) = 1.0_wp
    end do

    ! decide the processor of the boundary cell
    !
    do j = ncell+1, ncell+ncell_boundary
      do i = 1, domain%num_comm_proc
        k = iproc(i)
        if (domain_cell_rank(j) == k) then
          domain_cell_rank(j) = i
          exit
        end if
      end do
    end do
 
    return

  end subroutine setup_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_nobc_boundary
  !> @brief        Updates the system size under NOBC conditions.
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[inout] boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_nobc_boundary(domain, boundary)

    ! formal arguments
    type(s_domain),           intent(in)    :: domain
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    integer         :: i
    real(wip)       :: coord_min(1:3), coord_max(1:3), box_size(1:3)

    coord_min(1:3)          =  1000000000000.0_wip
    coord_max(1:3)          = -1000000000000.0_wip
    do i = 1, domain%num_atom_domain
      coord_min(1:3) = min(coord_min(1:3), domain%coord(i,1:3))
      coord_max(1:3) = max(coord_max(1:3), domain%coord(i,1:3))
    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, coord_min, 3, mpi_real8, mpi_min, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, coord_max, 3, mpi_real8, mpi_max, &
                       mpi_comm_country, ierror)
#endif
    box_size(1:3) = max(-coord_min(1:3), coord_max(1:3)) + 0.1_wip
    boundary%box_size_x = box_size(1)*2.0_wip
    boundary%box_size_y = box_size(2)*2.0_wip
    boundary%box_size_z = box_size(3)*2.0_wip
    if (boundary%box_size_x > boundary%box_size_x_max) &
      call error_msg('Error: Calculated box_size_x is greater than '//&
                     'box_size_x_max.')
    if (boundary%box_size_y > boundary%box_size_y_max) &
      call error_msg('Error: Calculated box_size_y is greater than '//&
                     'box_size_y_max.')
    if (boundary%box_size_z > boundary%box_size_z_max) &
      call error_msg('Error: Calculated box_size_z is greater than '//&
                     'box_size_z_max.')
    if (boundary%box_size_x < boundary%box_size_x_min) &
      write(Msgout, '(A)') 'Calculated box_size_x is less than box_size_x_min.'
    if (boundary%box_size_y < boundary%box_size_y_min) &
      write(Msgout, '(A)') 'Calculated box_size_y is less than box_size_y_min.'
    if (boundary%box_size_z < boundary%box_size_z_min) &
      write(Msgout, '(A)') 'Calculated box_size_z is less than box_size_z_min.'

    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,wp)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,wp)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,wp)

    return

  end subroutine update_nobc_boundary


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_pbc
  !> @brief        Sets up atom maps with all atoms.
  !! @authors      JJ  
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_pbc(molecule, boundary, enefunc, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(inout) :: enefunc 
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3), origin(3)
    integer                   :: step
    integer                   :: i, icx, icy, icz, icel, k, l1, l2
    integer                   :: id
    integer                   :: num_atom, num_charge
    integer                   :: icel_local, ncel_local
    integer                   :: ncel, natom_all
    integer                   :: patm, pchg
    integer                   :: start_i1
    integer                   :: ixx1, ix

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wp),         pointer :: cell_pbc_move(:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), ncharge(:), start_atom(:)

    bsize_x         => boundary%box_size_x
    bsize_y         => boundary%box_size_y
    bsize_z         => boundary%box_size_z
    csize_x         => boundary%cell_size_x
    csize_y         => boundary%cell_size_y
    csize_z         => boundary%cell_size_z
    ncel_x          => boundary%num_cells_x
    ncel_y          => boundary%num_cells_y
    ncel_z          => boundary%num_cells_z
 
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    cell_pbc_move   => domain%cell_pbc_move
    ncharge         => domain%num_charge
 
    origin(1)       = boundary%origin_x
    origin(2)       = boundary%origin_y
    origin(3)       = boundary%origin_z

    ncel            = domain%num_cell_local + domain%num_cell_boundary
    ncel_local      = domain%num_cell_local
    natom_all       = domain%num_atom_all


    ! charged particles
    !
    do step = 1, 2

      natom(1:ncel) = 0
      ncharge(1:ncel) = 0

      do i = 1, natom_all

        if (abs(molecule%charge(i)) > EPS .or. &
            molecule%residue_name(i)(1:3) == 'HIS') then

          !coordinate shifted against the origin
          !
          x_shift = molecule%atom_coord(1,i) - boundary%origin_x
          y_shift = molecule%atom_coord(2,i) - boundary%origin_y
          z_shift = molecule%atom_coord(3,i) - boundary%origin_z

          !coordinate shifted to the first quadrant and set into the boundary box
          !
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)

          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)

          !assign which cell
          !
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

          ! atoms inside the domain
          !
          if (cell_g2l(icel) /= 0) then

            ! local cell index
            !
            icel_local = cell_g2l(icel)
  
            patm = natom(icel_local)
            pchg = ncharge(icel_local)
            patm = patm + 1
            pchg = pchg + 1 
 
            if (step == 2) & 
              call molecule_to_domain(molecule, move, origin, i, &
                                      domain, domain%start_atom,        &
                                      icel_local, patm)
       
            natom(icel_local) = patm
            ncharge(icel_local) = pchg
  
          ! atoms in a boundary
          !
          else if (cell_g2b(icel) /= 0) then
  
            ! local cell index
            !
            icel_local = cell_g2b(icel) + ncel_local
  
            patm = natom(icel_local)
            pchg = ncharge(icel_local)
            patm = patm + 1
            pchg = pchg + 1

            move(1) = move(1) + cell_pbc_move(1,icel_local)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,icel_local)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,icel_local)*bsize_z
 
            if (step == 2) & 
              call molecule_to_domain(molecule, move, origin, i, &
                                      domain, domain%start_atom,        &
                                      icel_local, patm)
  
            natom(icel_local) = patm
            ncharge(icel_local) = pchg
  
          end if
  
        end if 
  
      end do
  
      do i = 1, natom_all

        if (abs(molecule%charge(i)) <= EPS .and. &
            molecule%residue_name(i)(1:3) /= 'HIS') then

          !coordinate shifted against the origin
          !
          x_shift = molecule%atom_coord(1,i) - boundary%origin_x
          y_shift = molecule%atom_coord(2,i) - boundary%origin_y
          z_shift = molecule%atom_coord(3,i) - boundary%origin_z

          !coordinate shifted to the first quadrant and set into the boundary box
          !
          move(1) = bsize_x*0.5_wip - bsize_x*anint(x_shift/bsize_x)
          move(2) = bsize_y*0.5_wip - bsize_y*anint(y_shift/bsize_y)
          move(3) = bsize_z*0.5_wip - bsize_z*anint(z_shift/bsize_z)
  
          x_shift = x_shift + move(1)
          y_shift = y_shift + move(2)
          z_shift = z_shift + move(3)
  
          !assign which cell
          !
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
  
          ! atoms inside the domain
          !
          if (cell_g2l(icel) /= 0) then
  
            ! local cell index
            !
            icel_local = cell_g2l(icel)
            patm = natom(icel_local)
            patm = patm + 1
 
            if (step == 2) & 
             call molecule_to_domain(molecule, move, origin, i, &
                                     domain, domain%start_atom,        &
                                     icel_local, patm)
  
            natom(icel_local) = patm
  
          ! atoms in a boundary
          !
          else if (cell_g2b(icel) /= 0) then
  
            ! local cell index
            !
            icel_local = cell_g2b(icel) + ncel_local
            patm = natom(icel_local)
            patm = patm + 1 
  
            move(1) = move(1) + cell_pbc_move(1,icel_local)*bsize_x
            move(2) = move(2) + cell_pbc_move(2,icel_local)*bsize_y
            move(3) = move(3) + cell_pbc_move(3,icel_local)*bsize_z
 
            if (step == 2) & 
              call molecule_to_domain(molecule, move, origin, i, &
                                      domain, domain%start_atom,        &
                                      icel_local, patm)
         
            natom(icel_local) = patm
  
          end if
  
        end if 

      end do

      if (step == 1) then
        domain%start_atom(1:ncel) = 0
        do i = 1, ncel-1
          domain%start_atom(i+1) = domain%start_atom(i) + natom(i)
        end do
        k = 0
        l1 = 0
        l2 = 0
        do i = 1, ncel_local
          k = k + natom(i)
          l1 = max(l1,ncharge(i))
          l2 = max(l2,natom(i)-ncharge(i))
        end do
        domain%num_atom_domain = k

        !$omp parallel private(id)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
        call get_para_range(1, k, nthread, id, domain%start_index(id), &
                            domain%end_index(id))
        !$omp end parallel

        k = 0
        do i = ncel_local+1, ncel
          k = k + natom(i)
          l1 = max(l1,ncharge(i))
          l2 = max(l2,natom(i)-ncharge(i))
        end do
        domain%num_atom_boundary = k
        num_atom = domain%num_atom_domain + domain%num_atom_boundary
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, num_atom, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
#endif
        MaxAtom_domain = max(num_atom*2, min_atom)
        num_atom = MaxAtom_domain 
        call alloc_domain(domain, DomainDynvar_Atom, num_atom, 1, 1)
        num_charge = 0
        do i = 1, ncel
          num_charge = num_charge + ncharge(i)
        end do
        enefunc%num_cg_elec = num_charge
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, num_charge, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
#endif
        Max_cg_elec = max(num_charge*2, min_atom)
        num_charge = Max_cg_elec
        call alloc_enefunc(enefunc, EneFuncCGElecList, num_charge) 
        call alloc_enefunc(enefunc, EneFuncCGElecInvList, MaxAtom_domain) 
      end if 
    end do

    do i = 1, domain%num_atom_domain
      domain%coord_prev(i,1) = domain%coord(i,1)
      domain%coord_prev(i,2) = domain%coord(i,2)
      domain%coord_prev(i,3) = domain%coord(i,3)
    end do

    k = 0
    do i = 1, ncel
      start_i1 = domain%start_atom(i)
      do ix = 1, natom(i)
        ixx1 = ix + start_i1
        domain%atom_2_cell(ixx1) = i
      end do
      do ix = 1, ncharge(i)
        ixx1 = ix + start_i1
        k = k + 1
        enefunc%cg_elec_list(k) = ixx1
        enefunc%cg_elec_list_inv(ixx1) = k
      end do 
    end do

    return

  end subroutine setup_atom_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_nobc
  !> @brief        Sets up atom maps with all atoms.
  !! @authors      JJ  
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_nobc(molecule, boundary, enefunc, domain)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: origin(3)
    integer                   :: num_atom, num_charge
    integer                   :: id
    integer                   :: i, icx, icy, icz, icel, k, step, l1, l2
    integer                   :: icel_local
    integer                   :: ncel_local
    integer                   :: ncel, natom_all
    integer                   :: patm, pchg
    integer                   :: start_i1
    integer                   :: ixx1, ix

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wp),         pointer :: cell_pbc_move(:,:)
    real(wip),        pointer :: vel(:,:)
    real(wp),         pointer :: charge(:)
    real(wip),        pointer :: mass(:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: natom(:), ncharge(:), start_atom(:)
    integer,          pointer :: atom_cls(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: atom_2_cel(:)

    bsize_x         => boundary%box_size_x
    bsize_y         => boundary%box_size_y
    bsize_z         => boundary%box_size_z
    csize_x         => boundary%cell_size_x
    csize_y         => boundary%cell_size_y
    csize_z         => boundary%cell_size_z
    ncel_x          => boundary%num_cells_x
    ncel_y          => boundary%num_cells_y
    ncel_z          => boundary%num_cells_z
 
    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    cell_pbc_move   => domain%cell_pbc_move
    ncharge         => domain%num_charge
    start_atom      => domain%start_atom
    vel             => domain%velocity
    charge          => domain%charge
    mass            => domain%mass   
    atom_cls        => domain%atom_cls_no
    chain_id        => domain%mol_chain_id
    id_l2g          => domain%id_l2g
    id_g2l          => domain%id_g2l
    atom_2_cel      => domain%atom_2_cell
 
    origin(1)       = boundary%origin_x
    origin(2)       = boundary%origin_y
    origin(3)       = boundary%origin_z

    ncel            = domain%num_cell_local + domain%num_cell_boundary
    ncel_local      = domain%num_cell_local
    natom_all       = domain%num_atom_all


    ! charged particles
    !
    do step = 1, 2

      natom(1:ncel) = 0
      ncharge(1:ncel) = 0

      do i = 1, natom_all

        if (abs(molecule%charge(i)) > EPS .or. &
            molecule%residue_name(i)(1:3) == 'HIS') then

          !coordinate shifted against the origin
          !
          x_shift = molecule%atom_coord(1,i) - boundary%origin_x
          y_shift = molecule%atom_coord(2,i) - boundary%origin_y
          z_shift = molecule%atom_coord(3,i) - boundary%origin_z
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip

          !assign which cell
          !
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
  
          ! atoms inside the domain
          !
          if (cell_g2l(icel) /= 0) then

            ! local cell index
            !
            icel_local = cell_g2l(icel)

            patm = natom(icel_local)
            pchg = ncharge(icel_local)
            patm = patm + 1
            pchg = pchg + 1 

            if (step == 2) &
              call molecule_to_domain_nobc(molecule, i, domain, start_atom, &
                                           icel_local, patm)
     
            natom(icel_local) = patm
            ncharge(icel_local) = pchg

          ! atoms in a boundary
          !
          else if (cell_g2b(icel) /= 0) then

            ! local cell index
            !
            icel_local = cell_g2b(icel) + ncel_local

            patm = natom(icel_local)
            pchg = ncharge(icel_local)
            patm = patm + 1
            pchg = pchg + 1

            if (step == 2) &
              call molecule_to_domain_nobc(molecule, i, domain, start_atom, &
                                           icel_local, patm)

            natom(icel_local) = patm
            ncharge(icel_local) = pchg

          end if

        end if 

      end do

      do i = 1, natom_all

        if (abs(molecule%charge(i)) <= EPS .and. &
            molecule%residue_name(i)(1:3) /= 'HIS') then

          !coordinate shifted against the origin
          !
          x_shift = molecule%atom_coord(1,i) - boundary%origin_x
          y_shift = molecule%atom_coord(2,i) - boundary%origin_y
          z_shift = molecule%atom_coord(3,i) - boundary%origin_z
          x_shift = x_shift + bsize_x*0.5_wip
          y_shift = y_shift + bsize_y*0.5_wip
          z_shift = z_shift + bsize_z*0.5_wip
  
          !assign which cell
          !
          icx = int(x_shift/csize_x)
          icy = int(y_shift/csize_y)
          icz = int(z_shift/csize_z)
          if (icx == ncel_x) icx = icx - 1
          if (icy == ncel_y) icy = icy - 1
          if (icz == ncel_z) icz = icz - 1
          icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y
  
          ! atoms inside the domain
          !
          if (cell_g2l(icel) /= 0) then
  
            ! local cell index
            !
            icel_local = cell_g2l(icel)
            patm = natom(icel_local)
            patm = patm + 1

            if (step == 2) &
              call molecule_to_domain_nobc(molecule, i, domain, start_atom, &
                                           icel_local, patm)

            natom(icel_local) = patm

          ! atoms in a boundary
          !
          else if (cell_g2b(icel) /= 0) then

            ! local cell index
            !
            icel_local = cell_g2b(icel) + ncel_local
            patm = natom(icel_local)
            patm = patm + 1 

            if (step == 2) &
              call molecule_to_domain_nobc(molecule, i, domain, start_atom, &
                                           icel_local, patm)
       
            natom(icel_local) = patm

          end if

        end if 

      end do

      if (step == 1) then
        start_atom(1:ncel) = 0
        do i = 1, ncel-1
          start_atom(i+1) = start_atom(i) + natom(i)
        end do
        k = 0
        l1 = 0
        l2 = 0
        do i = 1, ncel_local
          k = k + natom(i)
          l1 = max(l1,ncharge(i))
          l2 = max(l2,natom(i)-ncharge(i))
        end do
        domain%num_atom_domain = k
        k = 0
        do i = ncel_local+1, ncel
          k = k + natom(i)
          l1 = max(l1,ncharge(i))
          l2 = max(l2,natom(i)-ncharge(i))
        end do
        domain%num_atom_boundary = k

        !$omp parallel private(id)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        call get_para_range(1, k, nthread, id, domain%start_index(id), &
                            domain%end_index(id))
        !$omp end parallel

        num_atom = domain%num_atom_domain + domain%num_atom_boundary
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, num_atom, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
#endif
        MaxAtom_domain = num_atom * 2
        num_atom = MaxAtom_domain
        call alloc_domain(domain, DomainDynvar_Atom, num_atom, 1, 1)
        num_charge = 0
        do i = 1, ncel
          num_charge = num_charge + ncharge(i)
        end do
        enefunc%num_cg_elec = num_charge
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, num_charge, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
#endif
        Max_cg_elec = num_charge * 2
        num_charge = Max_cg_elec
        call alloc_enefunc(enefunc, EneFuncCGElecList, num_charge)
        call alloc_enefunc(enefunc, EneFuncCGElecInvList, MaxAtom_domain)
      end if
    end do

    k = 0
    do i = 1, ncel
      start_i1 = start_atom(i)
      do ix = 1, natom(i)
        ixx1 = ix + start_i1
        domain%atom_2_cell(ixx1) = i
      end do
      do ix = 1, ncharge(i)
        ixx1 = ix + start_i1
        k = k + 1
        enefunc%cg_elec_list(k) = ixx1
        enefunc%cg_elec_list_inv(ixx1) = k
      end do
    end do

    return

  end subroutine setup_atom_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_neighboring_cells
  !> @brief        Assigns neighboring cells.
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_neighbor_cells(enefunc, boundary, domain)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, j, ii, ic, jc, inbc, k
    integer                   :: ncel_local, nboundary

    integer,          pointer :: cell_g2l(:), cell_l2g(:)
    integer,          pointer :: cell_b2g(:), cell_g2b(:)

    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    cell_b2g        => domain%cell_b2g
    cell_l2g        => domain%cell_l2g

    ncel_local      =  domain%num_cell_local
    nboundary       =  domain%num_cell_boundary

    call alloc_domain(domain, DomainNeighbourCell, ncel_local+nboundary, 1, 1)

    do i = 1, ncel_local

      k  = 0
      ic = cell_l2g(i)

      do inbc = 1, 27

        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2l(jc)
        if (j /= 0) then
          domain%near_neighbor_cells(inbc,i) = j
        else if (cell_g2b(jc) /= 0) then
          domain%near_neighbor_cells(inbc,i) = cell_g2b(jc)+ncel_local
        end if
        if (j > i) then
          k  = k  + 1
          domain%neighbor_cells(k,i) = j
        else if (cell_g2b(jc) /= 0) then
          j  = cell_g2b(jc) + ncel_local
          k  = k  + 1
          domain%neighbor_cells(k,i) = j
        end if
      end do
      domain%near_neighbor_cell_count(i)  = k

      if (enefunc%forcefield /= ForcefieldGroMartini) then

        do inbc = 28, 125
          jc = boundary%neighbor_cells(inbc,ic)
          j  = cell_g2l(jc)
          if (j > i) then
            k  = k  + 1
            domain%neighbor_cells(k,i) = j
          else if (cell_g2b(jc) /= 0) then
            j  = cell_g2b(jc) + ncel_local
            k  = k  + 1
            domain%neighbor_cells(k,i) = j
          end if
        end do
        domain%neighbor_cell_count(i) = k

      end if

    end do

    do ii = 1, nboundary

      k  = 0
      i  = ii + ncel_local
      ic = cell_b2g(ii)

      do inbc = 1, 27
        jc = boundary%neighbor_cells(inbc,ic)
        j  = cell_g2b(jc) + ncel_local
        if (j > i) then
          k  = k  + 1
          domain%neighbor_cells(k,i) = j
        end if
      end do

      if (enefunc%forcefield /= ForcefieldGroMartini) then

        do inbc = 28, 125
          jc = boundary%neighbor_cells(inbc,ic)
          j  = cell_g2b(jc) + ncel_local 
          if (j > i) then
            k  = k  + 1
            domain%neighbor_cells(k,i) = j
          end if
        end do
        domain%neighbor_cell_count(i) = k

      end if

    end do

    return

  end subroutine assign_neighbor_cells
    
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_atoms
  !> @brief        Compares the total number of atoms in each cell.
  !! @authors      JJ  
  !! @param[in]    
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_atoms(natom, ncharge, num_atom, num_charge, &
                               cell_l2gx,  cell_l2gy, cell_l2gz,     &
                               neighbor,  ncel_local, nboundary)

    ! formal arguments
    real(wp),                intent(inout) :: natom(:)
    real(wp),                intent(inout) :: ncharge(:)
    integer,                 intent(in)    :: num_atom(:) 
    integer,                 intent(in)    :: num_charge(:) 
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: neighbor(:,:,:)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary

    ! local variables
    integer                  :: i, ii, i1, i2, i3 

    do i = 1, ncel_local
      natom(i)   = real((num_atom(i)+1)*nproc_city+my_city_rank, wp)
      ncharge(i) = real((num_charge(i)+1)*nproc_city+my_city_rank, wp)
    end do

    do ii = 1, nboundary

      i  = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)
      natom(i)   = real((num_atom(i)+1)*nproc_city+neighbor(i1,i2,i3),wp)
      ncharge(i) = real((num_charge(i)+1)*nproc_city+neighbor(i1,i2,i3),wp)

    end do

  end subroutine assign_cell_atoms


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_interactions
  !> @brief        Assigns the cell index for given cell-cell interaction.
  !! @authors      JJ  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_interaction(natom, ncharge, cell,               &
                                     cell_l2gx,  cell_l2gy,  cell_l2gz,  &
                                     cell_gxyz2l,               &
                                     ncel_local, nboundary,              &
                                     near_cells_count, near_cells,       &
                                     far_cells_count, far_cells)

    ! formal arguments
    real(wp),                intent(in)    :: natom(:)
    real(wp),                intent(in)    :: ncharge(:)
    integer,                 intent(in)    :: cell(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_gxyz2l(:,:,:)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary
    integer,                 intent(inout) :: near_cells_count(:)
    integer,                 intent(inout) :: near_cells(:,:)
    integer,                 intent(inout) :: far_cells_count(:)
    integer,                 intent(inout) :: far_cells(:,:)

    ! local variables
    real(wip)                :: ic1, ic2
    integer                  :: i, i1, i2, i3, j, ii, jj, ij
    integer                  :: icx1, icy1, icz1, icx2, icy2, icz2
    integer                  :: icx, icy, icz, k_near, k_far

    ij = 0

    ! assign the interaction cell for each interaction
    !
    do i = 1, ncel_local-1

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)
      k_near = 0
      k_far  = 0

      do j = i+1, ncel_local

        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (icx <= 1 .and. icy <= 1 .and. icz <= 1) then
          k_near = k_near + 1
          near_cells(k_near,i) = j
        else if (icx <= 2 .and. icy <= 2 .and. icz <= 2) then
          k_far = k_far + 1
          far_cells(k_far,i) = j
        end if

      end do

      near_cells_count(i)     = k_near
      far_cells_count(i)      = k_far

    end do

    do i = 1, ncel_local

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)
      k_near = near_cells_count(i)
      k_far  = far_cells_count(i)

      do jj = 1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)
        icx = abs(icx1-icx2)
        icy = abs(icy1-icy2)
        icz = abs(icz1-icz2)

        if (icx <= 1 .and. icy <= 1 .and. icz <= 1) then

          ic1 = natom(i)
          ic2 = natom(j)
          if (ic1 < ic2) then
            k_near = k_near + 1
            near_cells(k_near,i) = j
          end if

        else if (icx <= 2 .and. icy <= 2 .and. icz <= 2) then

          i1 = -1
          if (icx == 2) then
            i1 = (icx1+icx2)/2
          end if
          i2 = -1
          if (icy == 2) then
            i2 = (icy1+icy2)/2
          end if
          i3 = -1
          if (icz == 2) then
            i3 = (icz1+icz2)/2
          end if

          if (i1 >= 0 .and. i2 >= 0 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(i1+1,i2+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,i2+1,icz2+1))
            if (ic1 < ic2) then
              i3 = icz1
            else
              i3 = icz2
            end if
          else if (i1 >= 0 .and. i2 == -1 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(i1+1,icy1+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,icy2+1,i3+1))
            if (ic1 < ic2) then
              i2 = icy1
            else
              i2 = icy2
            end if
          else if (i1 == -1 .and. i2 >= 0 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,i2+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,i2+1,i3+1))
            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else if (i1 == -1 .and. i2 == -1 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,icy1+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,icy2+1,i3+1))
            if (ic1 < ic2) then
              i1 = icx1
              i2 = icy1
            else
              i1 = icx2
              i2 = icy2
            end if
          else if (i1 == -1 .and. i2 >= 0 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,i2+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,i2+1,icz2+1))
            if (ic1 < ic2) then
              i1 = icx1
              i3 = icz1
            else
              i1 = icx2
              i3 = icz2
            end if
          else if (i1 >= 0 .and. i2 == -1 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(i1+1,icy1+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,icy2+1,icz2+1))
            if (ic1 < ic2) then
              i2 = icy1
              i3 = icz1 
            else
              i2 = icy2
              i3 = icz2
            end if
          end if
          ij = cell_gxyz2l(i1+1,i2+1,i3+1)
          if (ij > 0 .and. ij <= ncel_local) then
            k_far = k_far + 1
            far_cells(k_far,i) = j
          end if

        end if

      end do

      near_cells_count(i)     = k_near
      far_cells_count(i)      = k_far

    end do

    do ii = 1, nboundary-1

      i = ii + ncel_local
      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)
      k_far  = 0

      do jj = ii+1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = abs(icx1-icx2)
        icy = abs(icy1-icy2)
        icz = abs(icz1-icz2)

        if (icx <= 1 .and. icy <= 1 .and. icz <= 1) then

          ic1 = natom(i)
          ic2 = natom(j)

        else if (icx <= 2 .and. icy <= 2 .and. icz <= 2) then

          i1 = -1
          if (icx == 2) then
            i1 = (icx1+icx2)/2
          end if
          i2 = -1
          if (icy == 2) then
            i2 = (icy1+icy2)/2
          end if
          i3 = -1
          if (icz == 2) then
            i3 = (icz1+icz2)/2
          end if

          if (i1 >= 0 .and. i2 >= 0 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(i1+1,i2+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,i2+1,icz2+1))
            if (ic1 < ic2) then
              i3 = icz1
            else
              i3 = icz2
            end if
          else if (i1 >= 0 .and. i2 == -1 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(i1+1,icy1+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,icy2+1,i3+1))
            if (ic1 < ic2) then
              i2 = icy1
            else
              i2 = icy2
            end if
          else if (i1 == -1 .and. i2 >= 0 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,i2+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,i2+1,i3+1))
            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else if (i1 == -1 .and. i2 == -1 .and. i3 >= 0) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,icy1+1,i3+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,icy2+1,i3+1))
            if (ic1 < ic2) then
              i1 = icx1
              i2 = icy1
            else
              i1 = icx2
              i2 = icy2
            end if
          else if (i1 == -1 .and. i2 >= 0 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(icx1+1,i2+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(icx2+1,i2+1,icz2+1))
            if (ic1 < ic2) then
              i1 = icx1
              i3 = icz1
            else
              i1 = icx2
              i3 = icz2
            end if
          else if (i1 >= 0 .and. i2 == -1 .and. i3 == -1) then
            ic1 = ncharge(cell_gxyz2l(i1+1,icy1+1,icz1+1))
            ic2 = ncharge(cell_gxyz2l(i1+1,icy2+1,icz2+1))
            if (ic1 < ic2) then
              i2 = icy1
              i3 = icz1
            else
              i2 = icy2
              i3 = icz2
            end if
          end if
          ij = cell_gxyz2l(i1+1,i2+1,i3+1)
          if (ij > 0 .and. ij <= ncel_local) then
            k_far = k_far + 1
            far_cells(k_far,i) = j
          end if
        end if

      end do

      far_cells_count(i)      = k_far

    end do

    return
  
  end subroutine assign_cell_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain
  !> @brief        Copies molecule information to the domain.
  !! @authors      JJ  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain(molecule, move, origin, iatom, domain, &
                                start_atom, icel, icel_atom)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    real(wip),                intent(in)    :: move(3)
    real(wip),                intent(in)    :: origin(3)
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: start_atom(:)
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)
    integer,          pointer :: molecule_no(:)
    integer                   :: start_icel

    coord        => molecule%atom_coord
    vel          => molecule%atom_velocity
    charge       => molecule%charge
    mass         => molecule%mass
    atom_class   => molecule%atom_cls_no
    molecule_no  => molecule%molecule_no

    start_icel = start_atom(icel)

    domain%id_l2g(icel_atom+start_icel) = iatom
    domain%id_g2l(iatom) = icel_atom + start_icel
    domain%mol_chain_id(icel_atom+start_icel) = molecule_no(iatom)

    domain%coord      (icel_atom+start_icel,1:3) = coord (1:3,iatom) - origin(1:3)
    domain%coord_ref  (icel_atom+start_icel,1:3) = coord (1:3,iatom) - origin(1:3)
    domain%velocity   (icel_atom+start_icel,1:3) = vel   (1:3,iatom)
    domain%charge     (icel_atom+start_icel)     = charge    (iatom)
    domain%mass       (icel_atom+start_icel)     = mass      (iatom)
    domain%atom_cls_no(icel_atom+start_icel)     = atom_class(iatom)
    domain%trans_vec  (icel_atom+start_icel,1:3) = real(move(1:3),wp)

    return

  end subroutine molecule_to_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain_nobc
  !> @brief        Copies molecule information to the domain.
  !! @authors      JJ  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain_nobc(molecule, iatom, domain, start_atom, &
                                     icel, icel_atom)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: start_atom(:)
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom

    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)
    integer,          pointer :: molecule_no(:)
    integer                   :: start_icel

    coord        => molecule%atom_coord
    vel          => molecule%atom_velocity
    charge       => molecule%charge
    mass         => molecule%mass
    atom_class   => molecule%atom_cls_no
    molecule_no  => molecule%molecule_no

    start_icel = start_atom(icel)

    domain%id_l2g(icel_atom+start_icel) = iatom
    domain%id_g2l(iatom) = icel_atom + start_icel
    domain%mol_chain_id(icel_atom+start_icel) = molecule_no(iatom)

    domain%coord      (icel_atom+start_icel,1:3) = coord (1:3,iatom) 
    domain%coord_ref  (icel_atom+start_icel,1:3) = coord (1:3,iatom)
    domain%velocity   (icel_atom+start_icel,1:3) = vel   (1:3,iatom)
    domain%charge     (icel_atom+start_icel)     = charge    (iatom)
    domain%mass       (icel_atom+start_icel)     = mass      (iatom)
    domain%atom_cls_no(icel_atom+start_icel)     = atom_class(iatom)

    return

  end subroutine molecule_to_domain_nobc


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord
  !> @brief        Checks atom coordinates.
  !! @authors      CK  
  !! @param[in]    ene_info      : ENERGY section control parameters information
  !! @param[in]    boundary      : boundary condition information
  !! @param[in]    contact_check : flag for contact_check
  !! @param[inout] domain        : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    integer                   :: i, ix, iy, ij, j, start_i, start_j
    integer                   :: id
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:), start_atom(:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: id_l2g(:)
    integer,          pointer :: near_cells_count(:)
    integer,          pointer :: near_cells(:,:)
    real(wp),         pointer :: trans1(:,:), trans2(:,:)
    real(wip),        pointer :: mass(:), coord(:,:) 
    real(wp),         pointer :: system_size(:)

    mass             => domain%mass 
    coord            => domain%coord
    trans1           => domain%trans_vec
    trans2           => domain%translated
    system_size      => domain%system_size
    natom            => domain%num_atom
    start_atom       => domain%start_atom
    ncell            => domain%num_cell_local
    nboundary        => domain%num_cell_boundary
    id_l2g           => domain%id_l2g
    near_cells_count => domain%near_cells_count
    near_cells       => domain%near_cells

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr, start_i, start_j)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, domain%num_atom_domain + domain%num_atom_boundary
        trans2(i,1:3) = real(coord(i,1:3),wp) + trans1(i,1:3)
      end do
    else
      do i = id+1, domain%num_atom_domain + domain%num_atom_boundary
        trans2(i,1:3) = real(coord(i,1:3),wp)
      end do
    endif
    !$omp barrier

    do i = id+1, ncell, nthread
      start_i = start_atom(i)
      do ix = 1, natom(i)
        dir(1:3) = trans2(ix+start_i,1:3)
        do iy = ix + 1, natom(i)
          dr(1:3) = dir(1:3) - trans2(iy+start_i,1:3)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          if (rr < EPS) then
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
            'Warning: Too small distances:',  &
            id_l2g(ix+start_i), id_l2g(iy+start_i), sqrt(rr)
            !$omp end critical
          endif
        end do
      end do
    end do

    do i = id+1, ncell+nboundary, nthread

      start_i = start_atom(i)

      do ix = 1, natom(i)

        dir(1:3) = trans2(ix+start_i,1:3) 

        do ij = 1, near_cells_count(i)
   
          j = near_cells(ij,i)
          start_j = start_atom(j)

          do iy = 1, natom(j)

            dr(1:3) = dir(1:3) - real(coord(iy+start_j,1:3),wp) 
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)')      &
                     'Warning: Too small distances:',  &
                      id_l2g(ix+start_i), id_l2g(iy+start_j), sqrt(rr)
              !$omp end critical
            endif
     
          end do
        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine check_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_domain_information
  !> @brief        Saves/loads data using temporary data before/after reallocation.
  !! @authors      JJ  
  !! @param[in]    direction   : 1 to temporary, 2 from temporary
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_domain_information(direction, domain)

    ! formal arguments
    integer,                 intent(in   )  :: direction
    type(s_domain),  target, intent(inout)  :: domain

    integer                  :: i
    real(wip),       pointer :: coord(:,:), velocity(:,:)        
    real(wip),       pointer :: mass(:)
    real(wp),        pointer :: charge(:)
    real(wip),       pointer :: buf_real(:)
    integer,         pointer :: atmcls(:), chain_id(:), id_l2g(:), atom_type(:)
    integer,         pointer :: buf_int(:)

    coord         => domain%coord
    velocity      => domain%velocity
    charge        => domain%charge
    mass          => domain%mass
    atmcls        => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    chain_id      => domain%mol_chain_id
    atom_type     => domain%NA_base_type
    buf_real      => domain%buf_var0_stay_real
    buf_int       => domain%buf_var0_stay_int

    if (direction == 1) then

      do i = 1, domain%num_atom_domain
        buf_real(8*i-7) = coord    (i,1)
        buf_real(8*i-6) = coord    (i,2)
        buf_real(8*i-5) = coord    (i,3)
        buf_real(8*i-4) = velocity (i,1)
        buf_real(8*i-3) = velocity (i,2)
        buf_real(8*i-2) = velocity (i,3)
        buf_real(8*i-1) = mass     (i  )
        buf_real(8*i  ) = charge   (i  )
        buf_int (4*i-3) = atmcls   (i  )
        buf_int (4*i-2) = id_l2g   (i  )
        buf_int (4*i-1) = chain_id (i  )
        buf_int (4*i  ) = atom_type(i  )
      end do

    else if (direction == 2) then

      do i = 1, domain%num_atom_domain
        coord    (i,1) = buf_real(8*i-7)
        coord    (i,2) = buf_real(8*i-6)
        coord    (i,3) = buf_real(8*i-5)
        velocity (i,1) = buf_real(8*i-4)
        velocity (i,2) = buf_real(8*i-3)
        velocity (i,3) = buf_real(8*i-2)
        mass     (i  ) = buf_real(8*i-1)
        charge   (i  ) = buf_real(8*i  )
        atmcls   (i  ) = buf_int (4*i-3)
        id_l2g   (i  ) = buf_int (4*i-2)
        chain_id (i  ) = buf_int (4*i-1)
        atom_type(i  ) = buf_int (4*i  )
      end do

    end if

    return

  end subroutine copy_domain_information

  subroutine neighbor_rank_x(cell, cell_start, cell_end, cell_to_rank, &
                             num_proc, iproc)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_start(:)
    integer,                  intent(in)    :: cell_end  (:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(inout) :: num_proc
    integer,                  intent(inout) :: iproc(:)

    integer                   :: a, iy, iz, key, rank, k

    a = cell_end(1) + 1
    if (a == cell(1)+1) a = 1
    do iz = cell_start(3), cell_end(3)
      do iy = cell_start(2), cell_end(2)
        rank = cell_to_rank(a,iy,iz)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    a = cell_start(1) - 1
    if (a == 0) a = cell(1)
    do iz = cell_start(3), cell_end(3)
      do iy = cell_start(2), cell_end(2)
        rank = cell_to_rank(a,iy,iz)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    return

  end subroutine neighbor_rank_x

  subroutine neighbor_rank_y(cell, cell_start, cell_end, cell_to_rank, &
                             num_proc, iproc)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_start(:)
    integer,                  intent(in)    :: cell_end  (:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(inout) :: num_proc
    integer,                  intent(inout) :: iproc(:)

    integer                   :: a, ix, ixx, iz, key, rank, k

    a = cell_end(2) + 1
    if (a == cell(2)+1) a = 1
    do iz = cell_start(3), cell_end(3)
      do ix = cell_start(1)-1, cell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        rank = cell_to_rank(ixx,a,iz)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    a = cell_start(2) - 1
    if (a == 0) a = cell(2)
    do iz = cell_start(3), cell_end(3)
      do ix = cell_start(1)-1, cell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        rank = cell_to_rank(ixx,a,iz)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    return

  end subroutine neighbor_rank_y

  subroutine neighbor_rank_z(cell, cell_start, cell_end, cell_to_rank, &
                             num_proc, iproc)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    integer,                  intent(in)    :: cell_start(:)
    integer,                  intent(in)    :: cell_end  (:)
    integer,                  intent(in)    :: cell_to_rank(:,:,:)
    integer,                  intent(inout) :: num_proc
    integer,                  intent(inout) :: iproc(:)

    integer                   :: a, ix, ixx, iy, iyy, key, rank, k

    a = cell_end(3) + 1
    if (a == cell(3)+1) a = 1
    do iy = cell_start(2)-1, cell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = cell_start(1)-1, cell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        rank = cell_to_rank(ixx,iyy,a)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    a = cell_start(3) - 1
    if (a == 0) a = cell(3)
    do iy = cell_start(2)-1, cell_end(2)+1
      if (iy == cell(2)+1) then
        iyy = 1
      else if (iy == 0) then
        iyy = cell(2)
      else
        iyy = iy
      end if
      do ix = cell_start(1)-1, cell_end(1)+1
        if (ix == cell(1)+1) then
          ixx = 1
        else if (ix == 0) then
          ixx = cell(1)
        else
          ixx = ix
        end if
        rank = cell_to_rank(ixx,iyy,a)
        key = 1
        do k = 1, num_proc
          if (rank == iproc(k)) then
            key = 0
            exit
          end if
        end do
        if (key == 1) then
          num_proc = num_proc + 1
          iproc(num_proc) = rank
        end if
      end do
    end do

    return

  end subroutine neighbor_rank_z

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_accum
  !> @brief        Copies domain data to the molecule.
  !! @authors      JJ  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_accum(boundary, domain, molecule)

    ! formal arguments
    type(s_boundary),         intent(in   ) :: boundary
    type(s_domain),   target, intent(inout) :: domain
    type(s_molecule), target, intent(inout) :: molecule

    integer                   :: i, ig, cell(3), ic(3), natom
    real(wp)                  :: shift(3), origin(3), bsize(3), move(3)
    real(wp)                  :: csize(3)

    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)
    integer,          pointer :: molecule_no(:)
    real(wip),        pointer :: coord_domain(:,:), vel_domain(:,:)
    real(wip),        pointer :: mass_domain(:)
    real(wp),         pointer :: charge_domain(:)
    integer,          pointer :: chain_id(:), class_domain(:), id_l2g(:)
    integer,          pointer :: natom_global(:,:,:)
    integer                   :: start_icel

    coord         => molecule%atom_coord
    vel           => molecule%atom_velocity
    charge        => molecule%charge
    mass          => molecule%mass
    atom_class    => molecule%atom_cls_no
    molecule_no   => molecule%molecule_no
    
    coord_domain  => domain%coord
    vel_domain    => domain%velocity
    charge_domain => domain%charge
    mass_domain   => domain%mass
    chain_id      => domain%mol_chain_id
    class_domain  => domain%atom_cls_no
    id_l2g        => domain%id_l2g
    natom_global  => domain%natom_global

    cell(1:3)     = domain%num_cell(1:3) 
    origin(1)     = boundary%origin_x
    origin(2)     = boundary%origin_y
    origin(3)     = boundary%origin_z
    bsize(1)      = boundary%box_size_x
    bsize(2)      = boundary%box_size_y
    bsize(3)      = boundary%box_size_z
    csize(1)      = boundary%cell_size_x
    csize(2)      = boundary%cell_size_y
    csize(3)      = boundary%cell_size_z
    natom         = molecule%num_atoms

    call alloc_domain(domain, DomainCellGlobal, cell(1), cell(2), cell(3))

    do i = 1, natom
      coord(1:3,i) = 0.0_wp
      vel  (1:3,i) = 0.0_wp
      charge(   i) = 0.0_wp
      mass  (   i) = 0.0_wp
      atom_class(i) = 0
      molecule_no(i) = 0
    end do
    natom_global(1:cell(1),1:cell(2),1:cell(3)) = 0

    do i = 1, domain%num_atom_domain

      ig = id_l2g(i)
      coord      (1:3,ig) = coord_domain (i,1:3)
      vel        (1:3,ig) = vel_domain   (i,1:3)
      charge     (    ig) = charge_domain(i)
      mass       (    ig) = mass_domain  (i)
      molecule_no(    ig) = chain_id     (i)
      atom_class (    ig) = class_domain (i)

      shift(1:3) = coord_domain(i,1:3) - origin(1:3)
      if (boundary%type == BoundaryTypePBC) then
        move(1:3) = bsize(1:3)*0.5_wip - bsize(1:3)*anint(shift(1:3)/bsize(1:3))
      else
        move(1:3) = bsize(1:3)*0.5_wip
      end if
      shift(1:3) = shift(1:3) + move(1:3)

      ic(1:3) = int(shift(1:3)/csize(1:3))
      if (ic(1) == cell(1)) ic(1) = ic(1) - 1
      if (ic(2) == cell(2)) ic(2) = ic(2) - 1
      if (ic(3) == cell(3)) ic(3) = ic(3) - 1
      ic(1:3) = ic(1:3) + 1

      natom_global(ic(1),ic(2),ic(3)) = natom_global(ic(1),ic(2),ic(3)) + 1

    end do

    call mpi_allreduce(mpi_in_place, coord, 3*natom, mpi_wp_real, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, vel  , 3*natom, mpi_wp_real, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, charge  , natom, mpi_wp_real, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, mass  , natom, mpi_wp_real, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, molecule_no, natom, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atom_class, natom, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, natom_global, cell(1)*cell(2)*cell(3), &
                       mpi_integer, mpi_sum, mpi_comm_country, ierror)

    return

  end subroutine molecule_accum

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_lb
  !> @brief        Sets up domain information with load balancing.
  !! @authors      JJ  
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_lb(boundary, molecule, enefunc, domain)

    ! formal arguments
    type(s_boundary),        intent(in   ) :: boundary
    type(s_molecule),        intent(in   ) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: cell(3), nc(3)
    integer                  :: ncel_local, ncel_bound, ncel_all

    cell(1:3) = domain%num_cell(1:3) 
    nc(1:3) = boundary%num_domain(1:3)

    ! decision of domains
    !
    call setup_domain_range_lb(cell, molecule, boundary, domain)

    ! assign the rank of each dimension from my_rank
    !
    call setup_cell_to_rank(cell, boundary, domain)

    ! check the MPI rank that should be communicated
    !
    call alloc_domain(domain, DomainNeighborRank, nproc_country, 1, 1)
    call setup_comm_rank(cell, domain)

    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local      = domain%num_cell_local
    ncel_bound      = domain%num_cell_boundary
    ncel_all        = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)

    ! assign global<->local mapping of cell index
    !
    call setup_cell_local(cell, domain)

    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, nc, domain)

    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all,            1, 1)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)

    if (boundary%type == BoundaryTypePBC) then
      call setup_atom_pbc(molecule, boundary, enefunc, domain)
    else
      call setup_atom_nobc(molecule, boundary, enefunc, domain)
    end if

    ! assign the interaction cell for each interaction
    !
    call assign_neighbor_cells(enefunc, boundary, domain)

    call setup_domain_interaction(boundary, domain)

    call check_atom_coord(boundary, domain)

    return

  end subroutine setup_domain_lb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_range_lb
  !> @brief        Determines the domain range for load balancing.
  !! @authors      JJ  
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_range_lb(cell, molecule, boundary, domain)

    ! formal arguments
    integer,                  intent(in)    :: cell(:)
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain


    ! local variable
    integer                   :: i, j, k, l
    integer                   :: nproc_xyz, nproc_x, nproc_y, nproc_z
    integer                   :: nc(3), iproc(3), dime
    integer                   :: my_xx_rank, my_yy_rank, my_zz_rank
    integer                   :: index_x, index_y, index_z
    integer                   :: natom_all

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wp),         pointer :: coord(:,:)
    integer,          pointer :: min_cell
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: natom_global(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)

    natom_global    => domain%natom_global
    cell_start      => domain%cell_start
    cell_end        => domain%cell_end

    nc(1)           = boundarY%num_domain(1)
    nc(2)           = boundarY%num_domain(2)
    nc(3)           = boundarY%num_domain(3)
    min_cell        => boundary%min_domain_cell

    natom_all       = domain%num_atom_all

    nproc_xyz = nc(1) * nc(2) * nc(3)
    nproc_x   = nc(1)
    nproc_y   = nc(2)
    nproc_z   = nc(3)
    cell_start(1) = 1
    cell_start(2) = 1
    cell_start(3) = 1
    cell_end(1)   = boundary%num_cells_x
    cell_end(2)   = boundary%num_cells_y
    cell_end(3)   = boundary%num_cells_z

    ! define mpi_rank in each dimension
    !
    iproc(1) = mod(my_country_rank, nc(1))
    iproc(2) = mod(my_country_rank/nc(1),nc(2))
    iproc(3) = my_country_rank/(nc(1)*nc(2))
    index_x  = iproc(2)*nproc_country + iproc(3)
    index_y  = iproc(3)*nproc_country + iproc(1)
    index_z  = iproc(1)*nproc_country + iproc(2)
    call mpi_comm_split(mpi_comm_country, index_x, my_country_rank, &
                        grid_commx,ierror)
    call mpi_comm_rank (grid_commx, my_x_rank, ierror)
    call mpi_comm_split(mpi_comm_country, index_y, my_country_rank, &
                        grid_commy,ierror)
    call mpi_comm_rank (grid_commy, my_y_rank, ierror)
    call mpi_comm_split(mpi_comm_country, index_z, my_country_rank, &
                        grid_commz,ierror)
    call mpi_comm_rank (grid_commz, my_z_rank, ierror)
    my_xx_rank = my_x_rank
    my_yy_rank = my_y_rank
    my_zz_rank = my_z_rank

    i = 0

    do while (nproc_xyz > 1)

      i = i + 1
      ! decide the dimension that should be divided
      !
      call decide_separate_dim(nproc_x, nproc_y, nproc_z, dime)

      ! decide the boundary
      !
      if (dime == 1) call select_boundary_x(cell, natom_global, min_cell, nproc_x, &
                              nproc_xyz, my_xx_rank, cell_start, cell_end)
      if (dime == 2) call select_boundary_y(cell, natom_global, min_cell, nproc_y, &
                              nproc_xyz, my_yy_rank, cell_start, cell_end)
      if (dime == 3) call select_boundary_z(cell, natom_global, min_cell, nproc_z, &
                              nproc_xyz, my_zz_rank, cell_start, cell_end)

    end do


    return

  end subroutine setup_domain_range_lb

end module cg_domain_mod
