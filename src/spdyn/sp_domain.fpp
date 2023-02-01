!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_domain_mod
!> @brief   utilities for domain decomposition
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_domain_mod

  use sp_constraints_mod
  use sp_boundary_mod
  use sp_energy_mod
  use sp_communicate_mod
  use sp_energy_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use hardwareinfo_mod
  use molecules_mod
  use timers_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters for setup cell capacity
  real(wip),        private, parameter :: VolumeBox8 = 512.0_wip
  real(wip),        private, parameter :: VboxRate   = 1.4_wip
  real(wip),        private, parameter :: ShrinkRate = 1.2_wip
  integer,          private, parameter :: NAtomBox8  = 120
  integer,          private, parameter :: NBondBox8  = 120
  integer,          private, parameter :: NAnglBox8  = 300
  integer,          private, parameter :: NDiheBox8  = 700
  integer,          private, parameter :: NImprBox8  = 50
  integer,          private, parameter :: NCmapBox8  = 12
  integer,          private, parameter :: NHGrpBox8  = 40
  integer,          private, parameter :: NHMovBox8  = 8
  integer,          private, parameter :: NAtomMax_in_CUDA = 255

  ! subroutines
  public  :: setup_domain
  public  :: setup_domain_pio
  public  :: setup_domain_interaction
  public  :: update_domain_interaction
  public  :: update_outgoing_ptl
  public  :: update_incoming_ptl
  public  :: setup_processor_rank
  public  :: setup_cell_capacity
  public  :: setup_cell_boundary
  private :: setup_solute_and_water
  private :: setup_hbond_group
  private :: setup_atom_by_HBond
  private :: assign_neighbor_cells
  private :: assign_cell_atoms
  private :: assign_cell_cpu
  private :: assign_cell_interaction
  private :: molecule_to_domain
  private :: check_atom_coord
  private :: check_atom_coord_pio
  private :: select_kernel
  ! FEP
  public  :: setup_domain_fep
  private :: setup_fep_correspondence

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain
  !> @brief        setup domain information
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain(ene_info, con_info, &
                          boundary, molecule, enefunc, constraints, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all
    integer                  :: num_solute_all


    ! initialize structure informations
    !
    call init_domain(domain)
    call init_enefunc(enefunc)
    call init_constraints(constraints)

    enefunc%table%num_all     = molecule%num_atoms

    domain%num_atom_all       = molecule%num_atoms            &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)
    domain%num_duplicate      = boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)

    molecule%num_deg_freedom  = molecule%num_deg_freedom      &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)

    if (main_rank .and. domain%num_atom_all /= molecule%num_atoms) &
      write(MsgOut,*) 'Total number of atoms   : ', &
                      domain%num_atom_all

    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z
    domain%system_size_ini(1) = boundary%box_size_x
    domain%system_size_ini(2) = boundary%box_size_y
    domain%system_size_ini(3) = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type

    ! select kernel
    !
    call select_kernel(ene_info, boundary, domain)

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)

    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity(boundary, domain, molecule)


    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)

    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do

    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)

    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all, 1, 1)

    ! decide hydrogen atom from mass
    !
    call check_light_atom_name(con_info%hydrogen_mass_upper_bound, molecule)

    call setup_solute_and_water(molecule, enefunc,       &
                                constraints%water_model, &
                                constraints%water_type,  &
                                constraints%num_water)

    num_solute_all            = enefunc%table%num_solute      &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)

    if (constraints%water_type == TIP4) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
    else if (constraints%water_type == TIP3) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
    else if (constraints%water_type == TIP2) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 2, 1)
    else if (constraints%water_type == TIP1) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 1, 1)
    end if

    call setup_hbond_group     (molecule, domain, enefunc, constraints)
    call setup_atom_by_HBond   (molecule, boundary, enefunc, constraints, &
                                domain)
    call setup_global_to_local_atom_index(enefunc, domain)

!   call setup_ring_check      (molecule, enefunc, constraints, domain)

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    call check_atom_coord(ene_info, boundary, domain)

    ! setup water molecule information
    !
    if (constraints%water_type >= TIP3) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
      domain%water%charge(1)      = enefunc%table%charge_O
      domain%water%charge(2)      = enefunc%table%charge_H
      domain%water%charge(3)      = enefunc%table%charge_H
      domain%water%charge(4)      = enefunc%table%charge_D
      domain%water%mass(1)        = enefunc%table%mass_O
      domain%water%mass(2)        = enefunc%table%mass_H
      domain%water%mass(3)        = enefunc%table%mass_H
      domain%water%mass(4)        = enefunc%table%mass_D
    else if (constraints%water_type == TIP2) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_1
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_2
      domain%water%charge(1)      = enefunc%table%charge_1
      domain%water%charge(2)      = enefunc%table%charge_2
      domain%water%mass(1)        = enefunc%table%mass_1
      domain%water%mass(2)        = enefunc%table%mass_2
    else
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_1
      domain%water%charge(1)      = enefunc%table%charge_1
      domain%water%mass(1)        = enefunc%table%mass_1
    end if

    return

  end subroutine setup_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_pio
  !> @brief        setup domain information
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_pio(ene_info, con_info, &
                              boundary, enefunc, constraints, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all


    ! initialize structure informations
    !
    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    enefunc%contact_check     = ene_info%contact_check
    enefunc%minimum_contact   = ene_info%minimum_contact

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type

    ! select kernel
    !
    call select_kernel(ene_info, boundary, domain)

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)

    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity(boundary, domain)

    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)

    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do

    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)

    ! allocation
    !
    call alloc_domain(domain, DomainDynvar, ncel_all, 1, 1)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)

    if (constraints%water_type == TIP4) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
    else if (constraints%water_type == TIP3) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
    else if (constraints%water_type == TIP2) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 2, 1)
    else if (constraints%water_type == TIP1) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 1, 1)
    end if

    ! assign particles in each domain
    !
    call setup_atom_by_HBond_pio(boundary, enefunc, constraints, domain)

    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    if (enefunc%contact_check) &
      call check_atom_coord_pio(enefunc, boundary, domain)

    ! setup water molecule information
    !
!    if (constraints%water_type >= TIP3) then
    if (constraints%water_type >= TIP2) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
      domain%water%charge(1)      = enefunc%table%charge_O
      domain%water%charge(2)      = enefunc%table%charge_H
      domain%water%charge(3)      = enefunc%table%charge_H
      domain%water%charge(4)      = enefunc%table%charge_D
      domain%water%mass(1)        = enefunc%table%mass_O
      domain%water%mass(2)        = enefunc%table%mass_H
      domain%water%mass(3)        = enefunc%table%mass_H
      domain%water%mass(4)        = enefunc%table%mass_D
    else
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_1
      domain%water%charge(1)      = enefunc%table%charge_1
      domain%water%mass(1)        = enefunc%table%mass_1
    end if

    return

  end subroutine setup_domain_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_interaction
  !> @brief        define the pairwise interaction between cells
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
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: cell_gxyz2l(:,:,:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx1(:), cell_l2gy1(:), cell_l2gz1(:)
    integer,          pointer :: neighbor(:,:,:), num_domain(:)
    integer,          pointer :: nsolute(:), nwater(:), natom_t0(:)
    integer(int1),    pointer :: virial_check(:,:)

    real(wp),         allocatable :: natom(:)


    bsize_x     => boundary%box_size_x
    bsize_y     => boundary%box_size_y
    bsize_z     => boundary%box_size_z
    num_domain  => boundary%num_domain

    cell_start   => domain%cell_start
    cell_end     => domain%cell_end
    cell_pair    => domain%cell_pair
    cell_move    => domain%cell_move
    virial_check => domain%virial_check
    cell_gxyz2l  => domain%cell_gxyz2l
    cell_l2gx    => domain%cell_l2gx
    cell_l2gy    => domain%cell_l2gy
    cell_l2gz    => domain%cell_l2gz
    cell_l2gx1   => domain%cell_l2gx_orig
    cell_l2gy1   => domain%cell_l2gy_orig
    cell_l2gz1   => domain%cell_l2gz_orig
    nsolute      => domain%num_solute
    nwater       => domain%num_water
    natom_t0     => domain%num_atom
    neighbor     => domain%neighbor

    cell(1)      = boundary%num_cells_x
    cell(2)      = boundary%num_cells_y
    cell(3)      = boundary%num_cells_z
    ncel_local   = domain%num_cell_local
    nboundary    = domain%num_cell_boundary

    allocate(natom(1:ncel_local+nboundary))

    ! check neighboring cells of each local cell
    !
    call assign_neighbor_cells(boundary, domain)

    ! number of atoms in each cell (proceesor number is also considered)
    !
    call assign_cell_atoms(natom, natom_t0, nsolute, nwater, &
                           cell_l2gx, cell_l2gy, cell_l2gz,  &
                           cell_start, cell_end,             &
                           neighbor, ncel_local, nboundary)

    ! assign the interaction cell for each interaction
    !
    call assign_cell_interaction(natom, num_domain, cell,             &
                                 cell_l2gx,  cell_l2gy,  cell_l2gz,   &
                                 cell_l2gx1, cell_l2gy1, cell_l2gz1,  &
                                 cell_gxyz2l, cell_start, cell_end,   &
                                 ncel_local, nboundary, cell_move,    &
                                 cell_pair, virial_check)

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%near_neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
        end if
      end do
    end do
    maxcell_near = ij
    do i = 1, ncel_local+nboundary
      do inbc = domain%near_neighbor_cell_count(i)+1,                    &
                domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
        end if
      end do
    end do
    maxcell = ij
    univ_maxcell = ij + ncel_local

    call alloc_domain(domain, DomainCellPairList, maxcell, &
                      ncel_local+nboundary, 1)

    if (domain%pairlist_kernel == PLK_GPU) &
    call alloc_domain(domain, DomainUnivCellPairList, univ_maxcell, &
                      ncel_local+nboundary, 1)

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%near_neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
          domain%cell_pairlist1(1,ij) = i
          domain%cell_pairlist1(2,ij) = j
          domain%cell_pairlist2(j,i) = ij
          domain%cell_pairlist2(i,j) = ij
          if (domain%ncell_pairlist1(i) == 0) &
            domain%istcell_pairlist1(i) = ij
          domain%ncell_pairlist1(i) = domain%ncell_pairlist1(i) + 1
        end if
      end do
    end do
    maxcell_near = ij
    do i = 1, ncel_local+nboundary
      do inbc = domain%near_neighbor_cell_count(i)+1,                       &
                domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
          domain%cell_pairlist1(1,ij) = i
          domain%cell_pairlist1(2,ij) = j
          domain%cell_pairlist2(j,i) = ij
          domain%cell_pairlist2(i,j) = ij
          if (domain%ncell_pairlist2(i) == 0) &
            domain%istcell_pairlist2(i) = ij
          domain%ncell_pairlist2(i) = domain%ncell_pairlist2(i) + 1
        end if
      end do
    end do
    maxcell = ij

    if (domain%pairlist_kernel == PLK_GPU) then

      ij = 0

      ! cell-pairs (same cell)
      !
      do i = 1, ncel_local
        ij = ij + 1
        domain%univ_cell_pairlist1(1,ij) = i
        domain%univ_cell_pairlist1(2,ij) = i
        domain%univ_cell_pairlist2(i,i) = ij
      end do

      ! cell-pairs (different cell)
      !
      do i = 1, ncel_local+nboundary
        do inbc = 1, domain%near_neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%univ_cell_pairlist1(1,ij) = i
            domain%univ_cell_pairlist1(2,ij) = j
            domain%univ_cell_pairlist2(j,i) = ij
            domain%univ_cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      do i = 1, ncel_local+nboundary
        do inbc = domain%near_neighbor_cell_count(i)+1,                      &
             domain%neighbor_cell_count(i)
          j = domain%neighbor_cells(inbc,i)
          if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
            ij = ij + 1
            domain%univ_cell_pairlist1(1,ij) = i
            domain%univ_cell_pairlist1(2,ij) = j
            domain%univ_cell_pairlist2(j,i) = ij
            domain%univ_cell_pairlist2(i,j) = ij
          end if
        end do
      end do
      univ_maxcell = ij
      univ_ncell_near = ncel_local + maxcell_near
    end if

    deallocate(natom)

    return

  end subroutine setup_domain_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_domain_interaction
  !> @brief        update the pairwise interaction between cells
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_domain_interaction(boundary, enefunc, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wip)                 :: calc_time
    integer                   :: i, j, ij, inbc
    integer                   :: ncel_local, nboundary, cell(3)

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: cell_start(:), cell_end(:)
    integer(int2),    pointer :: cell_pair(:,:)
    integer(int2),    pointer :: cell_gxyz2l(:,:,:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,          pointer :: cell_l2gx1(:), cell_l2gy1(:), cell_l2gz1(:)
    integer,          pointer :: neighbor(:,:,:), num_domain(:)
    integer(int1),    pointer :: virial_check(:,:)

    real(wp),         allocatable :: cpu_time(:), time_proc(:)


    bsize_x      => boundary%box_size_x
    bsize_y      => boundary%box_size_y
    bsize_z      => boundary%box_size_z
    num_domain   => boundary%num_domain

    cell_start   => domain%cell_start
    cell_end     => domain%cell_end
    cell_pair    => domain%cell_pair
    cell_move    => domain%cell_move
    cell_gxyz2l  => domain%cell_gxyz2l
    cell_l2gx    => domain%cell_l2gx
    cell_l2gy    => domain%cell_l2gy
    cell_l2gz    => domain%cell_l2gz
    cell_l2gx1   => domain%cell_l2gx_orig
    cell_l2gy1   => domain%cell_l2gy_orig
    cell_l2gz1   => domain%cell_l2gz_orig
    neighbor     => domain%neighbor
    virial_check => domain%virial_check

    cell(1)      = boundary%num_cells_x
    cell(2)      = boundary%num_cells_y
    cell(3)      = boundary%num_cells_z
    ncel_local   = domain%num_cell_local
    nboundary    = domain%num_cell_boundary


    allocate(cpu_time(1:ncel_local+nboundary), time_proc(1:nproc_city))

    ! cpu time(real) of each processor
    !
    if (enefunc%pme_use) then
      calc_time = total_time(3)+total_time(4)+total_time(5)+total_time(8)
    else
      calc_time = total_time(3)+total_time(4)+total_time(5)+total_time(7)
    end if

    calc_time      = calc_time - calc_time_prev
    calc_time_prev = calc_time_prev + calc_time

#ifdef HAVE_MPI_GENESIS
    call mpi_allgather(calc_time, 1, mpi_wip_real, &
                       time_proc, 1, mpi_wp_real, mpi_comm_country, ierror)
#else
    time_proc(1) = calc_time
#endif

    ! cpu time (cpu time of domain * random number of each cell)
    !
    call assign_cell_cpu(time_proc, domain%random,        &
                         cell_l2gx, cell_l2gy, cell_l2gz, &
                         cell_start, cell_end, neighbor,  &
                         ncel_local, nboundary, cpu_time)

    ! assign the interaction cell for each interaction
    !
    call assign_cell_interaction(cpu_time, num_domain, cell,          &
                                 cell_l2gx,  cell_l2gy,  cell_l2gz,   &
                                 cell_l2gx1, cell_l2gy1, cell_l2gz1,  &
                                 cell_gxyz2l, cell_start, cell_end,   &
                                 ncel_local, nboundary, cell_move,    &
                                 cell_pair, virial_check)

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
        end if
      end do
    end do
    maxcell = ij

    !!call alloc_domain(domain, DomainCellPairList, maxcell, &
    !!                  ncel_local+nboundary,1)
    domain%cell_pairlist1(1:2,1:maxcell) = 0
    domain%cell_pairlist2(1:ncel_local+nboundary,1:ncel_local+nboundary) = 0

    ij = 0
    do i = 1, ncel_local+nboundary
      do inbc = 1, domain%neighbor_cell_count(i)
        j = domain%neighbor_cells(inbc,i)
        if (cell_pair(j,i) > 0 .and. cell_pair(j,i) <= ncel_local) then
          ij = ij + 1
          domain%cell_pairlist1(1,ij) = i
          domain%cell_pairlist1(2,ij) = j
          domain%cell_pairlist2(j,i) = ij
          domain%cell_pairlist2(i,j) = ij
        end if
      end do
    end do
    maxcell = ij

    deallocate(cpu_time, time_proc)

    return

  end subroutine update_domain_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_outgoing_ptl
  !> @brief        update particles in each cell
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_outgoing_ptl(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    real(wip)                 :: x_shift, y_shift, z_shift
    real(wip)                 :: move(3)
    integer                   :: i, k, ix, icx, icy, icz, icel, ncel
    integer                   :: icel_local, icel_bd

    real(wip),        pointer :: bsize_x, bsize_y, bsize_z
    real(wip),        pointer :: csize_x, csize_y, csize_z
    real(wip),        pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wip),        pointer :: mass(:,:)
    real(wip),        pointer :: buf_real(:,:,:)
    integer,          pointer :: ncel_x, ncel_y, ncel_z
    integer,          pointer :: atmcls(:,:), id_l2g(:,:)
    integer,          pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,          pointer :: buf_int(:,:,:)


    bsize_x  => boundary%box_size_x
    bsize_y  => boundary%box_size_y
    bsize_z  => boundary%box_size_z
    ncel_x   => boundary%num_cells_x
    ncel_y   => boundary%num_cells_y
    ncel_z   => boundary%num_cells_z
    csize_x  => boundary%cell_size_x
    csize_y  => boundary%cell_size_y
    csize_z  => boundary%cell_size_z

    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real


    ! initializaiton
    !
    ncel = domain%num_cell_local + domain%num_cell_boundary
    ptl_exit(1:domain%num_cell_local) = 0
    ptl_add (1:ncel) = 0

    ! Check outgoing particles
    !
    do i = 1, domain%num_cell_local

      k = 0
      do ix = 1, domain%num_atom(i)

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

        !assign which cell
        !
        icx = int(x_shift/csize_x)
        icy = int(y_shift/csize_y)
        icz = int(z_shift/csize_z)
        if (icx == ncel_x) icx = icx - 1
        if (icy == ncel_y) icy = icy - 1
        if (icz == ncel_z) icz = icz - 1

        icel = 1 + icx + icy*ncel_x + icz*ncel_x*ncel_y

        icel_local = domain%cell_g2l(icel)
        icel_bd    = domain%cell_g2b(icel)

        if (icel_local /= i) then

          ptl_exit(i) = ptl_exit(i) + 1
          ptlindex(ptl_exit(i),i) = ix

          if (icel_local /= 0) then

            ptl_add(icel_local) = ptl_add(icel_local) + 1
            buf_real(1,ptl_add(icel_local),icel_local) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_local),icel_local) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_local),icel_local) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_local),icel_local) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_local),icel_local) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_local),icel_local) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_local),icel_local) = charge(ix,i)
            buf_real(8,ptl_add(icel_local),icel_local) = mass(ix,i)
            buf_int (1,ptl_add(icel_local),icel_local) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_local),icel_local) = id_l2g(ix,i)

          else if (icel_bd /= 0) then

            icel_bd = icel_bd + domain%num_cell_local
            ptl_add(icel_bd) = ptl_add(icel_bd) + 1
            buf_real(1,ptl_add(icel_bd),icel_bd) = coord(1,ix,i)
            buf_real(2,ptl_add(icel_bd),icel_bd) = coord(2,ix,i)
            buf_real(3,ptl_add(icel_bd),icel_bd) = coord(3,ix,i)
            buf_real(4,ptl_add(icel_bd),icel_bd) = velocity(1,ix,i)
            buf_real(5,ptl_add(icel_bd),icel_bd) = velocity(2,ix,i)
            buf_real(6,ptl_add(icel_bd),icel_bd) = velocity(3,ix,i)
            buf_real(7,ptl_add(icel_bd),icel_bd) = charge(ix,i)
            buf_real(8,ptl_add(icel_bd),icel_bd) = mass(ix,i)
            buf_int (1,ptl_add(icel_bd),icel_bd) = atmcls(ix,i)
            buf_int (2,ptl_add(icel_bd),icel_bd) = id_l2g(ix,i)

          end if
        end if
      end do
    end do

    return

  end subroutine update_outgoing_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_incoming_ptl
  !> @brief        update particles in each cell
  !! @authors      JJ
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_incoming_ptl(domain)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ix, kx
    logical                  :: insert

    real(wip),       pointer :: coord(:,:,:), velocity(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wip),       pointer :: mass(:,:)
    real(wip),       pointer :: buf_real(:,:,:)
    integer,         pointer :: atmcls(:,:), id_l2g(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: ptl_add(:), ptl_exit(:), ptlindex(:,:)
    integer,         pointer :: buf_int(:,:,:)


    coord    => domain%coord
    velocity => domain%velocity
    charge   => domain%charge
    mass     => domain%mass
    atmcls   => domain%atom_cls_no
    id_l2g   => domain%id_l2g
    id_g2l   => domain%id_g2l
    ptl_add  => domain%ptl_add
    ptl_exit => domain%ptl_exit
    ptlindex => domain%ptl_exit_index
    buf_int  => domain%buf_integer
    buf_real => domain%buf_real


    ! Incoming particles
    !
    do i = 1, domain%num_cell_local

      ! When the number of coming particles is larger than that of outgoing ones
      !
      if (ptl_add(i) >= ptl_exit(i)) then

        do k = 1, ptl_exit(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        do k = ptl_exit(i)+1, ptl_add(i)
          ix = k + domain%num_atom(i) - ptl_exit(i)
          coord(1,ix,i)               = buf_real(1,k,i)
          coord(2,ix,i)               = buf_real(2,k,i)
          coord(3,ix,i)               = buf_real(3,k,i)
          velocity(1,ix,i)            = buf_real(4,k,i)
          velocity(2,ix,i)            = buf_real(5,k,i)
          velocity(3,ix,i)            = buf_real(6,k,i)
          charge(ix,i)                = buf_real(7,k,i)
          mass(ix,i)                  = buf_real(8,k,i)
          atmcls(ix,i)                = buf_int(1,k,i)
          id_l2g(ix,i)                = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ix
        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)

      ! When the number of coming particles is less than that of outgoing ones
      !
      else

        do k = 1, ptl_add(i)
          coord(1,ptlindex(k,i),i)    = buf_real(1,k,i)
          coord(2,ptlindex(k,i),i)    = buf_real(2,k,i)
          coord(3,ptlindex(k,i),i)    = buf_real(3,k,i)
          velocity(1,ptlindex(k,i),i) = buf_real(4,k,i)
          velocity(2,ptlindex(k,i),i) = buf_real(5,k,i)
          velocity(3,ptlindex(k,i),i) = buf_real(6,k,i)
          charge(ptlindex(k,i),i)     = buf_real(7,k,i)
          mass(ptlindex(k,i),i)       = buf_real(8,k,i)
          atmcls(ptlindex(k,i),i)     = buf_int(1,k,i)
          id_l2g(ptlindex(k,i),i)     = buf_int(2,k,i)
          id_g2l(1,buf_int(2,k,i))    = i
          id_g2l(2,buf_int(2,k,i))    = ptlindex(k,i)
        end do

        j  = 0
        ix = domain%num_atom(i)
        k  = ptl_add(i) + 1

        do while (j < (ptl_exit(i)-ptl_add(i)))

          insert = .true.
          do kx = k, ptl_exit(i)
            if (ix == ptlindex(kx,i)) then
              insert = .false.
              j = j + 1
              exit
            end if
          end do

          if (insert) then
            kx = ptlindex(k,i)
            coord(1,kx,i)          = coord(1,ix,i)
            coord(2,kx,i)          = coord(2,ix,i)
            coord(3,kx,i)          = coord(3,ix,i)
            velocity(1,kx,i)       = velocity(1,ix,i)
            velocity(2,kx,i)       = velocity(2,ix,i)
            velocity(3,kx,i)       = velocity(3,ix,i)
            charge(kx,i)           = charge(ix,i)
            mass(kx,i)             = mass(ix,i)
            atmcls(kx,i)           = atmcls(ix,i)
            id_l2g(kx,i)           = id_l2g(ix,i)
            id_g2l(1,id_l2g(kx,i)) = i
            id_g2l(2,id_l2g(kx,i)) = kx

            j = j + 1
            k = k + 1

          end if

          ix = ix - 1

        end do

        domain%num_atom(i) = domain%num_atom(i) + ptl_add(i) - ptl_exit(i)

      end if

      domain%num_solute(i) = domain%num_atom(i)

    end do

    return

  end subroutine update_incoming_ptl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_rank
  !> @brief        define the processor rank in each dimension and decide the
  !!               number of cells in each domain
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !! @@aram[out]   cell     : cells in boundary
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_rank(boundary, domain, cell)

    ! formal arguments
    type(s_boundary),target, intent(in)    :: boundary
    type(s_domain),  target, intent(inout) :: domain
    integer,                 intent(inout) :: cell(3)

    ! local variables
    integer                  :: iproc(3)
    integer                  :: i, j, k, quotient, remainder

    integer, pointer         :: ncell_local, ncell_boundary, num_domain(:)
    integer, pointer         :: cell_start(:), cell_end(:), cell_length(:)
    integer, pointer         :: iproc_lower(:), iproc_upper(:), neighbor(:,:,:)


    num_domain     => boundary%num_domain

    ncell_local    => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end
    cell_length    => domain%cell_length
    iproc_lower    => domain%iproc_lower
    iproc_upper    => domain%iproc_upper
    neighbor       => domain%neighbor

    ! Assign the rank of each dimension from my_rank
    !
    iproc(1) = mod(my_city_rank, num_domain(1))
    iproc(2) = mod(my_city_rank/num_domain(1), num_domain(2))
    iproc(3) = my_city_rank/(num_domain(1)*num_domain(2))

    ! Cell number of each dimension
    !
    cell(1) = boundary%num_cells_x
    cell(2) = boundary%num_cells_y
    cell(3) = boundary%num_cells_z

    ! Default value of the number of cell in each domain
    !
    ncell_local = 1

    ! Assign the cell index for each processor and the total number of cells
    !
    do i = 1, 3

      quotient = cell(i) / num_domain(i)
      remainder = mod(cell(i), num_domain(i))

      if (iproc(i) <= (remainder -1)) then
        quotient       = quotient + 1
        cell_start(i)  = quotient * iproc(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      else
        cell_start(i)  = (quotient+1)*remainder + quotient*(iproc(i)-remainder)
        cell_start(i)  = cell_start(i) + 1
        cell_end(i)    = cell_start(i) + quotient - 1
        cell_length(i) = quotient
        ncell_local    = ncell_local * quotient
      end if

    end do

    ! Assign the lower processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = num_domain(i) - 1

      iproc_lower(i) = iproc(1) + &
                       iproc(2)*num_domain(1) + &
                       iproc(3)*num_domain(1)*num_domain(2)

      iproc(i) = iproc(i) + 1
      if (iproc(i) == num_domain(i)) &
        iproc(i) = 0

    end do

    ! Assign the upper processor index
    !
    do i = 1, 3

      iproc(i) = iproc(i) + 1
      if (iproc(i) == num_domain(i)) &
        iproc(i) = 0

      iproc_upper(i) = iproc(1) + &
                       iproc(2)*num_domain(1) + &
                       iproc(3)*num_domain(1)*num_domain(2)

      iproc(i) = iproc(i) - 1
      if (iproc(i) == -1) &
        iproc(i) = num_domain(i) - 1

    end do

    ! Assign the neighboring process index
    !
    do i = -1, 1

      if (i == -1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = num_domain(1) - 1

      else if (i == 1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == num_domain(1)) &
          iproc(1) = 0

      end if

      do j = -1, 1

        if (j == -1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = num_domain(2) - 1

        else if (j == 1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == num_domain(2)) &
            iproc(2) = 0

        end if

        do k = -1, 1

          if (k == -1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = num_domain(3) - 1

          else if (k == 1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == num_domain(3)) &
              iproc(3) = 0
          end if

          neighbor(i,j,k) = iproc(1) + &
                            iproc(2)*num_domain(1) + &
                            iproc(3)*num_domain(1)*num_domain(2)

          if (k == -1) then

            iproc(3) = iproc(3) + 1
            if (iproc(3) == num_domain(3)) &
              iproc(3) = 0

          else if (k == 1) then

            iproc(3) = iproc(3) - 1
            if (iproc(3) == -1) &
              iproc(3) = num_domain(3) - 1

          end if
        end do

        if (j == -1) then

          iproc(2) = iproc(2) + 1
          if (iproc(2) == num_domain(2)) &
            iproc(2) = 0

        else if (j == 1) then

          iproc(2) = iproc(2) - 1
          if (iproc(2) == -1) &
            iproc(2) = num_domain(2) - 1
        end if

      end do

      if (i == -1) then

        iproc(1) = iproc(1) + 1
        if (iproc(1) == num_domain(1)) &
          iproc(1) = 0

      else if (i == 1) then

        iproc(1) = iproc(1) - 1
        if (iproc(1) == -1) &
          iproc(1) = num_domain(1) - 1

      end if
    end do

    ! maximum number of boundary cells
    !
    ncell_boundary = 2*(cell_length(1) * cell_length(2) +                   &
                        cell_length(2) * cell_length(3) +                   &
                        cell_length(1) * cell_length(3))+                   &
                     4*(cell_length(1) + cell_length(2) + cell_length(3)) + &
                     8
    return

  end subroutine setup_processor_rank

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_capacity
  !> @brief        setup cell capacity for memory allocation
  !! @authors      JJ, HO
  !! @param[in]    boundary : boundary information
  !! @param[in]    domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_capacity(boundary, domain, molecule)

    ! formal arguments
    type(s_boundary),           intent(in) :: boundary
    type(s_domain),             intent(in) :: domain
    type(s_molecule), optional, intent(in) :: molecule

    ! local variables
    real(wip)                :: v, v_rate
    integer                  :: natom, natom2, nbond, nangle
    integer                  :: ndihedral, nimproper, ncmap

    ! FEP
    integer                  :: natom_partA, natom_partB

    if (present(molecule)) then
      natom = molecule%num_atoms
      nbond = molecule%num_bonds
      nangle = molecule%num_angles
      ndihedral = molecule%num_dihedrals
      nimproper = molecule%num_impropers
      ncmap = molecule%num_cmaps
    end if

    v = boundary%box_size_x * &
        boundary%box_size_y * &
        boundary%box_size_z

    v = v / real(boundary%num_domain(1)* &
                 boundary%num_domain(2)* &
                 boundary%num_domain(3),wip)
    v = v / real(domain%num_cell_local,wip)

    v_rate = VboxRate*v / VolumeBox8


    ! sp_domain_str
    !

    MaxAtom      = int(v_rate * real(NAtomBox8,wip))
    if (domain%pairlist_kernel == PLK_GPU) &
    MaxAtom      = min(MaxAtom,NAtomMax_in_CUDA)
    inv_MaxAtom  = 1.0_wp / real(MaxAtom,wp)

    MaxWater     = int(v_rate * real(NAtomBox8,wip) / 3.0_wip)
    MaxMove      = int(v_rate * real(NAtomBox8,wip) / 5.0_wip)
    MaxWaterMove = int(v_rate * real(NAtomBox8,wip) / 7.0_wip)
    MaxWater = MaxWater * 3

    ! If the system is vacuum, MaxAtom, MaxWater, etc. are set to 
    ! the number of atoms in the system.
    if (present(molecule)) then
      MaxAtom      = min(natom, MaxAtom)
      MaxWater     = min(natom, MaxWater)
      MaxMove      = min(natom, MaxMove)
      MaxWaterMove = min(natom, MaxWaterMove)
    end if


    ! sp_enefunc_str
    !

    MaxBond      = int(v_rate * real(NBondBox8,wip) * ShrinkRate)
    MaxAngle     = int(v_rate * real(NAnglBox8,wip) * ShrinkRate)
    MaxDihe      = int(v_rate * real(NDiheBox8,wip) * ShrinkRate)
    MaxImpr      = int(v_rate * real(NImprBox8,wip) * ShrinkRate)
    MaxCmap      = int(v_rate * real(NCmapBox8,wip) * ShrinkRate)

    BondMove     = MaxBond  / 2
    AngleMove    = MaxAngle / 2
    DiheMove     = MaxDihe  / 2
    ImprMove     = MaxImpr  / 2
    CmapMove     = MaxCmap  / 2

    ! If vacuum, MaxBond, MaxAngle, etc. are set to the number of
    ! bonds, angles, etc. of the molecule.
    if (present(molecule)) then
      MaxBond      = min(nbond, MaxBond)
      MaxAngle     = min(nangle, MaxAngle)
      MaxDihe      = min(10*ndihedral, MaxDihe)
      MaxImpr      = min(nimproper, MaxImpr)
      MaxCmap      = min(ncmap, MaxCmap)
      BondMove     = min(nbond, BondMove)
      AngleMove    = min(nangle, AngleMove)
      DiheMove     = min(10*ndihedral, DiheMove)
      ImprMove     = min(nimproper, ImprMove)
      CmapMove     = min(ncmap, CmapMove)
    end if


    ! sp_pairlist_str
    !

    MaxNb15      = int(v_rate * real(NAtomBox8,wip))
    if (present(molecule)) MaxNb15 = min(natom, MaxNb15)
    MaxNb15      = MaxNb15 * MaxNb15
    MaxNb15      = MaxNb15 * domain%scale_pairlist_generic

    MaxNb15_Fugaku = int(v_rate * real(NAtomBox8,wip)) * 12
    if (present(molecule)) MaxNb15_Fugaku = min(natom, MaxNb15_Fugaku)
    MaxNb15_Fugaku = MaxNb15_Fugaku * domain%scale_pairlist_fugaku

    MaxNb15Water = int(v_rate * real(NAtomBox8,wip))
    if (present(molecule)) MaxNb15Water = min(natom, MaxNb15Water)
    MaxNb15Water = MaxNb15Water ** 2
    MaxNb15water = MaxNb15Water / 10

    ! FEP: Number of pairs between perturbed atoms and others
    ! is, at most, MaxNb15 * (Number of perturbed atoms).
    if (domain%fep_use) then
      MaxNb15_fep = int(v_rate * real(NAtomBox8,wip))
      natom_partA = molecule%num_atoms_fep(1)+molecule%num_atoms_fep(2)
      natom_partB = molecule%num_atoms_fep(3)+molecule%num_atoms_fep(4)
      MaxNb15_fep = min(natom, MaxNb15_fep) * max(natom_partA, natom_partB)
      MaxNb15_fep = MaxNb15_fep * domain%scale_pairlist_generic
    end if

    ! If vacuum, MaxNb15 and MaxNb15Water are set to natom**2.
    ! MaxContact and ContactMove are also set to natom**2.
    ! If the number of contact in the system is given, change natom**2 to it.
    if (present(molecule)) then
      MaxContact   = min(MaxNb15, MaxContact)
      ContactMove  = min(MaxNb15, ContactMove)
    end if


    ! sp_constraints_str
    !

    HGroupMax    = int(v_rate * real(NHGrpBox8,wip) * ShrinkRate)
    HGrpMaxMove  = int(v_rate * real(NHMovBox8,wip) * ShrinkRate)

    ! If vacuum, HGroupMax and HGrpMaxMove are set to natom.
    if (present(molecule)) then
      HGroupMax    = min(natom, HGroupMax)
      HGrpMaxMove  = min(natom, HGrpMaxMove)
    end if


#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) 'Cell volume                      : ', v
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxBond      : ', MaxBond
      write(MsgOut,*) 'sp_enefunc_str    ::MaxAngle     : ', MaxAngle
      write(MsgOut,*) 'sp_enefunc_str    ::MaxDihe      : ', MaxDihe
      write(MsgOut,*) 'sp_enefunc_str    ::MaxImpr      : ', MaxImpr
      write(MsgOut,*) 'sp_enefunc_str    ::MaxCmap      : ', MaxCmap
      write(MsgOut,*) 'sp_enefunc_str    ::BondMove     : ', BondMove
      write(MsgOut,*) 'sp_enefunc_str    ::AngleMove    : ', AngleMove
      write(MsgOut,*) 'sp_enefunc_str    ::DiheMove     : ', DiheMove
      write(MsgOut,*) 'sp_enefunc_str    ::ImprMove     : ', ImprMove
      write(MsgOut,*) 'sp_enefunc_str    ::CmapMove     : ', CmapMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      if (domain%fep_use) &
        write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15_fep  : ', MaxNb15_fep
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    return

  end subroutine setup_cell_capacity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cell_boundary
  !> @brief        define boundary cells in each domain
  !! @authors      JJ
  !! @param[in]    cell   : cell count for each axis
  !! @param[in]    ndom   : domain count for each axis
  !! @param[inout] domain : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cell_boundary(cell, ndom, domain)

    ! formal arguments
    integer,                 intent(in)    :: cell(3)
    integer,                 intent(in)    :: ndom(3)
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, ic, jc, kc, icel, icel_local

    integer,         pointer :: ncell, ncell_boundary
    integer,         pointer :: cell_start(:), cell_end(:)
    integer,         pointer :: cell_b2g(:)
    integer(int2),   pointer :: cell_g2b(:)
    integer(int2),   pointer :: cell_gxyz2l(:,:,:)
    integer,         pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)
    integer,         pointer :: cell_l2gx_orig(:)
    integer,         pointer :: cell_l2gy_orig(:)
    integer,         pointer :: cell_l2gz_orig(:)
    real(wp),        pointer :: cell_pbc_move(:,:)


    ncell          => domain%num_cell_local
    ncell_boundary => domain%num_cell_boundary
    cell_start     => domain%cell_start
    cell_end       => domain%cell_end
    cell_g2b       => domain%cell_g2b
    cell_b2g       => domain%cell_b2g
    cell_gxyz2l    => domain%cell_gxyz2l
    cell_l2gx      => domain%cell_l2gx
    cell_l2gy      => domain%cell_l2gy
    cell_l2gz      => domain%cell_l2gz
    cell_l2gx_orig => domain%cell_l2gx_orig
    cell_l2gy_orig => domain%cell_l2gy_orig
    cell_l2gz_orig => domain%cell_l2gz_orig
    cell_pbc_move  => domain%cell_pbc_move

    icel_local = 0

    if (ndom(1) > 1) then
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

      if (ndom(2) > 1) then
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
            cell_l2gx(icel_local+ncell) = ic
            cell_l2gy(icel_local+ncell) = jc
            cell_l2gz(icel_local+ncell) = k
            cell_l2gx_orig(icel_local+ncell) = i
            cell_l2gy_orig(icel_local+ncell) = j
            cell_l2gz_orig(icel_local+ncell) = k
            cell_gxyz2l(ic,jc,k) = icel_local+ncell

          end do
        end do

        if (ndom(3) > 1) then

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
              cell_l2gx(icel_local+ncell) = ic
              cell_l2gy(icel_local+ncell) = jc
              cell_l2gz(icel_local+ncell) = kc
              cell_l2gx_orig(icel_local+ncell) = i
              cell_l2gy_orig(icel_local+ncell) = j
              cell_l2gz_orig(icel_local+ncell) = k
              cell_gxyz2l(ic,jc,kc) = icel_local+ncell

            end do
          end do
        end if
      end if
    end if

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

    return

  end subroutine setup_cell_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_solute_and_water
  !> @brief        setup solute and water list
  !! @authors      JJ
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] water_model : water model
  !! @param[inout] water_type  : water type
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_solute_and_water(molecule, enefunc, water_model, &
                                    water_type, num_water)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    character(*),            intent(inout) :: water_model
    integer,                 intent(inout) :: water_type
    integer,                 intent(inout) :: num_water

    ! local variables
    integer                  :: i, k, natom, nwater, nsolute, water_atom
    integer                  :: io, ih, id(4)
    real(wp)                 :: mass(4), mass_max, mass_min


    natom   = molecule%num_atoms
    nwater  = 0
    nsolute = 0

    do i = 1, natom
      if ((molecule%residue_name(i)(1:3) .eq. 'TIP' .or.  &
           molecule%residue_name(i)(1:3) .eq. 'WAT' .or.  &
           molecule%residue_name(i)(1:3) .eq. 'MIZ' .or.  &
           molecule%residue_name(i)(1:3) .eq. 'SOL').and. &
          .not.(molecule%light_atom_mass(i))      .and. &
          .not.(molecule%light_atom_name(i))      .and. &
           molecule%atom_name(i)(1:1) .ne. 'H') then
        nwater = nwater + 1
      else if (molecule%residue_name(i)(1:3) .ne. 'TIP' .and. &
               molecule%residue_name(i)(1:3) .ne. 'WAT' .and. &
               molecule%residue_name(i)(1:3) .ne. 'MIZ' .and. &
               molecule%residue_name(i)(1:3) .ne. 'SOL') then
        nsolute = nsolute + 1
      end if
    end do

    if ((natom-nsolute)/4 == nwater .and. nwater > 0) then
      water_type = TIP4
    else if ((natom-nsolute)/3 == nwater) then
      water_type = TIP3
    else if ((natom-nsolute)/2 == nwater .and. nwater > 0) then
      water_type = TIP2
    else if (natom-nsolute == nwater .and. nwater > 0) then
      water_type = TIP1
    else
      call error_msg('Setup_Solute_And_Water> # of water is incorrect.')
    end if

    enefunc%table%num_water  = nwater
    enefunc%table%num_solute = nsolute
    num_water = nwater
    call alloc_enefunc(enefunc, EneFuncTableWat, nwater,  1)
    call alloc_enefunc(enefunc, EneFuncTableSol, nsolute, natom)

    if (water_type == TIP4) then
      water_atom = 4
    else if (water_type == TIP3) then
      water_atom = 3
    else if (water_type == TIP2) then
      water_atom = 2
    else if (water_type == TIP1) then
      water_atom = 1
    end if

    i       = 1
    nwater  = 0
    nsolute = 0

    do while (.true.)

      if (i > natom) exit

      if (molecule%residue_name(i)(1:3) .eq. 'TIP' .or. &
          molecule%residue_name(i)(1:3) .eq. 'WAT' .or. &
          molecule%residue_name(i)(1:3) .eq. 'MIZ' .or. &
          molecule%residue_name(i)(1:3) .eq. 'SOL') then

        do k = 1, water_atom
          mass(k) = molecule%mass(i-1+k)
        end do
        mass_max = -1000.0_wp
        mass_min =  1000.0_wp
        do k = 1, water_atom
          mass_max = max(mass_max,mass(k))
          mass_min = min(mass_min,mass(k))
        end do

        id(1:water_atom) = 0

        if (water_atom == 4) then
          do k = 1, water_atom
            if (mass(k) == mass_max) id(1) = i-1+k
            if (mass(k) == mass_min) id(water_atom) = i-1+k
            if (mass(k) > mass_min .and. mass(k) < mass_max) then
              if (id(2) == 0) id(2) = i-1+k
              if (id(2) /= 0) id(3) = i-1+k
            end if
          end do
        else if (water_atom == 3) then
          do k = 1, water_atom
            if (mass(k) == mass_max) id(1) = i-1+k
            if (mass(k) == mass_min) then
              if (id(2) == 0) id(2) = i-1+k
              if (id(2) /= 0) id(3) = i-1+k
            end if
          end do
        else if (water_atom == 2) then
          id(1)=i
          id(2)=i+1
        else if (water_atom == 1) then
          id(1) = i
        end if
        nwater = nwater + 1
        enefunc%table%water_list(1:water_atom,nwater) = id(1:water_atom)
        i = i + water_atom

      else

        nsolute = nsolute + 1
        enefunc%table%solute_list(nsolute) = i
        enefunc%table%solute_list_inv(i) = nsolute
        i = i + 1

      end if

    end do

    if (nwater /= enefunc%table%num_water) &
      call error_msg('Setup_Solute_And_Water> number of water is incorrect')

    if (nsolute /= enefunc%table%num_solute) &
      call error_msg('Setup_Solute_And_Water> number of solute is incorrect')


    ! set water oxygen and hydrogen
    !
    ! if (size(enefunc%table%water_list(1,:)) /= 0 .and. &
    !     water_atom >= 3) then
    if (size(enefunc%table%water_list(1,:)) /= 0 .and. &
        water_atom >= 3) then

      io = enefunc%table%water_list(1,1)
      ih = enefunc%table%water_list(2,1)

      enefunc%table%atom_cls_no_O = molecule%atom_cls_no(io)
      enefunc%table%atom_cls_no_H = molecule%atom_cls_no(ih)
      enefunc%table%charge_O      = molecule%charge(io)
      enefunc%table%charge_H      = molecule%charge(ih)
      enefunc%table%mass_O        = molecule%mass(io)
      enefunc%table%mass_H        = molecule%mass(ih)

      ! dummy atom
      !
      if (water_atom == 4) then
        io = enefunc%table%water_list(4,1)
        enefunc%table%atom_cls_no_D = molecule%atom_cls_no(io)
        enefunc%table%charge_D      = molecule%charge(io)
        enefunc%table%mass_D        = molecule%mass(io)
      end if

    else if (size(enefunc%table%water_list(1,:)) /= 0 .and. &
             water_atom == 2) then

      io = enefunc%table%water_list(1,1)
      ih = enefunc%table%water_list(2,1)
      enefunc%table%atom_cls_no_1  = molecule%atom_cls_no(io)
      enefunc%table%atom_cls_no_2  = molecule%atom_cls_no(ih)
      enefunc%table%charge_1       = molecule%charge(io)
      enefunc%table%charge_2       = molecule%charge(ih)
      enefunc%table%mass_1         = molecule%mass(io)
      enefunc%table%mass_2         = molecule%mass(ih)

    else if (size(enefunc%table%water_list(1,:)) /= 0 .and. &
             water_atom == 1) then

      io = enefunc%table%water_list(1,1)
      enefunc%table%atom_cls_no_1 = molecule%atom_cls_no(io)
      enefunc%table%charge_1      = molecule%charge(io)
      enefunc%table%mass_1        = molecule%mass(io)

    end if

    return

  end subroutine setup_solute_and_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_ring_check
  !> @brief        check rings
  !! @authors      JJ
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_ring_check(molecule, enefunc, constraints, domain)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain      

    ! local variables
    integer                  :: i, i1, i2, nhmax, connect1, connect2
    integer                  :: j, k, l, m
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2

    ! count the number of bonds including hydrogen
    !
    call alloc_constraints(constraints, ConstraintsRing, molecule%num_atoms)

    do i = 1, molecule%num_atoms
      constraints%ring(i) = 0
      constraints%num_connect(1:3,i) = 0
    end do

    ! first connectivity
    !
    do i = 1, molecule%num_bonds

      i1 = molecule%bond_list(1,i)
      i2 = molecule%bond_list(2,i)
      mi1 = molecule%light_atom_mass(i1)
      mi2 = molecule%light_atom_mass(i2)
      cl1 = molecule%light_atom_name(i1)
      cl2 = molecule%light_atom_name(i2)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
        cl2 = mi2
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1)
        cl2 = (cl2 .or. mi2)
      end if

      if (enefunc%table%solute_list_inv(i1) /= 0 .and. &
          enefunc%table%solute_list_inv(i2) /= 0) then
        if (.not.cl1 .and. .not.cl2) then
          connect1 = constraints%num_connect(1,i1)
          connect2 = constraints%num_connect(1,i2)
          connect1 = connect1 + 1
          connect2 = connect2 + 1
          constraints%connectivity(connect1,1,i1) = i2
          constraints%connectivity(connect2,1,i2) = i1
          constraints%num_connect(1,i1) = connect1
          constraints%num_connect(1,i2) = connect2
        end if
      end if

    end do

    ! second connectivity
    !
    do i = 1, molecule%num_atoms
      if (enefunc%table%solute_list_inv(i) /= 0) then
        do j = 1, constraints%num_connect(1,i)
          i1 = constraints%connectivity(j,1,i)
          do k = 1, constraints%num_connect(1,i1)
            i2 = constraints%connectivity(k,1,i1)
            if (i /= i2) then
              connect1 = constraints%num_connect(2,i)
              connect1 = connect1 + 1
              constraints%connectivity(connect1,2,i) = i2
              constraints%num_connect(2,i) = connect1 
            end if
          end do
        end do
      end if
    end do

    ! third connectivity
    !
    do i = 1, molecule%num_atoms
      if (enefunc%table%solute_list_inv(i) /= 0) then
        do j = 1, constraints%num_connect(2,i)
          i1 = constraints%connectivity(j,2,i)
          do k = 1, constraints%num_connect(1,i1)
            i2 = constraints%connectivity(k,1,i1)
            cl1 = .true.
            if (i == i2 ) cl1 = .false.
            do l = 1, constraints%num_connect(1,i)
              m = constraints%connectivity(l,1,i)
              if (i2 == m) then
                cl1 = .false.
                exit
              end if
            end do
            if (cl1) then
              connect1 = constraints%num_connect(3,i)
              connect1 = connect1 + 1
              constraints%connectivity(connect1,3,i) = i2
              constraints%num_connect(3,i) = connect1
            end if
          end do
        end do
      end if
    end do

    ! check ring
    !
    do i = 1, molecule%num_atoms
      if (enefunc%table%solute_list_inv(i) /= 0) then
        do j = 1, 3
          do k = 2, 3
            do l = 1, constraints%num_connect(j,i)
              i1 = constraints%connectivity(l,j,i) 
              do m = 1, constraints%num_connect(k,i)
                i2 = constraints%connectivity(m,k,i)
                if (i1 == i2 .and. ((j == k .and. l /= m) .or. (j /= k))) then
                  constraints%ring(i) = 1
                  exit
                end if
              end do
              if (constraints%ring(i) == 1) exit
            end do
            if (constraints%ring(i) == 1) exit
          end do
          if (constraints%ring(i) == 1) exit
        end do
      end if
    end do

    ! ring in each domain
    !
    do i = 1, enefunc%table%num_all
      i1 = domain%id_g2l(1,i)
      if (enefunc%table%solute_list_inv(i) /= 0 .and. i1 > 0) then
        i2 = domain%id_g2l(2,i)
        domain%ring(i2,i1) = constraints%ring(i)
      end if
    end do

    call dealloc_constraints(constraints, ConstraintsRing)
         
    return

  end subroutine setup_ring_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hbond_group
  !> @brief        setup bonds including hydrogen atom
  !! @authors      JJ, HO
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hbond_group(molecule, domain, enefunc, constraints)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints

    ! local variables
    integer                  :: i, i1, i2, nhmax, k
    logical                  :: mi1, mi2
    logical                  :: cl1, cl2


    ! count the number of bonds including hydrogen
    !
    call alloc_constraints(constraints, ConstraintsHBond, molecule%num_atoms)

    do i = 1, molecule%num_bonds

      i1 = molecule%bond_list(1,i)
      i2 = molecule%bond_list(2,i)
      mi1 = molecule%light_atom_mass(i1)
      mi2 = molecule%light_atom_mass(i2)
      cl1 = molecule%light_atom_name(i1)
      cl2 = molecule%light_atom_name(i2)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
        cl2 = mi2
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1)
        cl2 = (cl2 .or. mi2)
      end if

      if (enefunc%table%solute_list_inv(i1) /= 0 .and. &
          enefunc%table%solute_list_inv(i2) /= 0) then

        if (cl1 .or. cl2) then

          if (domain%fep_use) then
            ! Hydrogen rewiring in FEP:
            ! When singleA-dualA or singleB-dualB bonds have hydrogen atoms,
            ! SHAKE problems occur. If coordinate and velocity are synchonized
            ! after SHAKE, coordinate and velocity of atoms involved in the bonds
            ! are moved and then the constraints are sometimes broken. To avoid
            ! this problem, hydrogen atoms beloging to dualB is rewired to connect
            ! with atoms of singleA. The constraints of the hydrogen atoms are
            ! considered in performing SHAKE for molecule A. After SHAKE, singleB
            ! atoms are synchonized to singleA.
            !
            ! singleA--dualA
            !                 ==>   singleA--dualA
            ! singleB--dualB               \_dualB
            ! 
            ! To rewire the bonds, for singleA-dualB bond including hydrogen,
            ! the atom index of singleB is replaced with the corresponding atom
            ! index of singleA.
            if ((int(molecule%fepgrp(i1)) == 2) .and. &
                (int(molecule%fepgrp(i2)) == 4)) then
              do k = 1, molecule%num_atoms_fep(1)
                if (molecule%id_singleB(k) == i1) then
                  i1 = molecule%id_singleA(k)
                  exit
                end if
              end do
            else if ((int(molecule%fepgrp(i1)) == 4) .and. &
                     (int(molecule%fepgrp(i2)) == 2)) then
              do k = 1, molecule%num_atoms_fep(2)
                if (molecule%id_singleB(k) == i2) then
                  i2 = molecule%id_singleA(k)
                  exit
                end if
              end do
            end if
          end if

          if (cl1) then
            constraints%duplicate(i2) = constraints%duplicate(i2) + 1
            constraints%H_index(constraints%duplicate(i2),i2) = i1
          else
            constraints%duplicate(i1) = constraints%duplicate(i1) + 1
            constraints%H_index(constraints%duplicate(i1),i1) = i2
          end if

        end if

      end if

    end do


    ! count XHn group for each number n
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then
        constraints%nh(1) = constraints%nh(1) + 1

      else if (constraints%duplicate(i1) == 2) then
        constraints%nh(2) = constraints%nh(2) + 1

      else if (constraints%duplicate(i1) == 3) then
        constraints%nh(3) = constraints%nh(3) + 1

      else if (constraints%duplicate(i1) == 4) then
        constraints%nh(4) = constraints%nh(4) + 1

      else if (constraints%duplicate(i1) == 5) then
        constraints%nh(5) = constraints%nh(5) + 1

      else if (constraints%duplicate(i1) == 6) then
        constraints%nh(6) = constraints%nh(6) + 1

      else if (constraints%duplicate(i1) == 7) then
        constraints%nh(7) = constraints%nh(7) + 1

      else if (constraints%duplicate(i1) == 8) then
        constraints%nh(8) = constraints%nh(8) + 1

      else if (constraints%duplicate(i1) >= 8) then
        call error_msg( &
             'Setup_HBond_Group> Bond(>8) for one atom is not considered')

      end if

    end do

    constraints%connect  = 0
    if (constraints%nh(1) /= 0) constraints%connect = 1
    if (constraints%nh(2) /= 0) constraints%connect = 2
    if (constraints%nh(3) /= 0) constraints%connect = 3
    if (constraints%nh(4) /= 0) constraints%connect = 4
    if (constraints%nh(5) /= 0) constraints%connect = 5
    if (constraints%nh(6) /= 0) constraints%connect = 6
    if (constraints%nh(7) /= 0) constraints%connect = 7
    if (constraints%nh(8) /= 0) constraints%connect = 8

    nhmax = max(constraints%nh(1),constraints%nh(2),constraints%nh(3), &
                constraints%nh(4),constraints%nh(5),constraints%nh(6), &
                constraints%nh(7),constraints%nh(8))

    call alloc_constraints(constraints, ConstraintsBondGroup, nhmax)


    ! Make a list of XHn
    !
    constraints%nh(1:8) = 0

    do i = 1, enefunc%table%num_solute

      i1 = enefunc%table%solute_list(i)

      if (constraints%duplicate(i1) == 1) then

        constraints%nh(1) = constraints%nh(1) + 1
        constraints%H_Group(1,constraints%nh(1),1) = i1
        constraints%H_Group(2,constraints%nh(1),1) = constraints%H_index(1,i1)

      else if (constraints%duplicate(i1) == 2) then

        constraints%nh(2) = constraints%nh(2) + 1
        constraints%H_Group(1,constraints%nh(2),2) = i1
        constraints%H_Group(2,constraints%nh(2),2) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(2),2) = constraints%H_index(2,i1)

      else if (constraints%duplicate(i1) == 3) then

        constraints%nh(3) = constraints%nh(3) + 1
        constraints%H_Group(1,constraints%nh(3),3) = i1
        constraints%H_Group(2,constraints%nh(3),3) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(3),3) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(3),3) = constraints%H_index(3,i1)

      else if (constraints%duplicate(i1) == 4) then

        constraints%nh(4) = constraints%nh(4) + 1
        constraints%H_Group(1,constraints%nh(4),4) = i1
        constraints%H_Group(2,constraints%nh(4),4) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(4),4) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(4),4) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(4),4) = constraints%H_index(4,i1)

      else if (constraints%duplicate(i1) == 5) then

        constraints%nh(5) = constraints%nh(5) + 1
        constraints%H_Group(1,constraints%nh(5),5) = i1
        constraints%H_Group(2,constraints%nh(5),5) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(5),5) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(5),5) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(5),5) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(5),5) = constraints%H_index(5,i1)

      else if (constraints%duplicate(i1) == 6) then

        constraints%nh(6) = constraints%nh(6) + 1
        constraints%H_Group(1,constraints%nh(6),6) = i1
        constraints%H_Group(2,constraints%nh(6),6) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(6),6) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(6),6) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(6),6) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(6),6) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(6),6) = constraints%H_index(6,i1)

      else if (constraints%duplicate(i1) == 7) then

        constraints%nh(7) = constraints%nh(7) + 1
        constraints%H_Group(1,constraints%nh(7),7) = i1
        constraints%H_Group(2,constraints%nh(7),7) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(7),7) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(7),7) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(7),7) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(7),7) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(7),7) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(7),7) = constraints%H_index(7,i1)

      else if (constraints%duplicate(i1) == 8) then

        constraints%nh(8) = constraints%nh(8) + 1
        constraints%H_Group(1,constraints%nh(8),8) = i1
        constraints%H_Group(2,constraints%nh(8),8) = constraints%H_index(1,i1)
        constraints%H_Group(3,constraints%nh(8),8) = constraints%H_index(2,i1)
        constraints%H_Group(4,constraints%nh(8),8) = constraints%H_index(3,i1)
        constraints%H_Group(5,constraints%nh(8),8) = constraints%H_index(4,i1)
        constraints%H_Group(6,constraints%nh(8),8) = constraints%H_index(5,i1)
        constraints%H_Group(7,constraints%nh(8),8) = constraints%H_index(6,i1)
        constraints%H_Group(8,constraints%nh(8),8) = constraints%H_index(7,i1)
        constraints%H_Group(9,constraints%nh(8),8) = constraints%H_index(8,i1)

      end if

    end do

    return

  end subroutine setup_hbond_group

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_HBond
  !> @brief        setup atom maps with H-bond connection groups
  !! @authors      JJ
  !! @param[inout] molecule    : molecule information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_HBond(molecule, boundary, enefunc, constraints, &
                                 domain)

    ! formal arguments
    type(s_molecule),    target, intent(inout) :: molecule
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_enefunc),     target, intent(in)    :: enefunc
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: shift(3), bsize(3), bsize_orig(3), csize(3)
    real(wip)                    :: move(3), origin(3), coord_temp(3)
    integer                      :: ioffset
    integer                      :: i, j, ic(1:3), icel, i1, j1, k1, k
    integer                      :: dupl_x, dupl_y, dupl_z
    integer                      :: dupl(3), num_cell(3), num_cell_block(3)
    integer                      :: duplicate_cell_startx, duplicate_cell_endx
    integer                      :: duplicate_cell_starty, duplicate_cell_endy
    integer                      :: duplicate_cell_startz, duplicate_cell_endz
    integer                      :: isolute, iwater, ih1, ih2, id
    integer                      :: icel_local
    integer                      :: ncel_local, ncel
    integer                      :: patm, psol, pwat, pnoh, phgl
    character(4)                 :: ci1
    logical                      :: mi1, cl1

    real(wip),           pointer :: bsize_x, bsize_y, bsize_z
    real(wip),           pointer :: csize_x, csize_y, csize_z
    real(wp),            pointer :: cell_pbc_move(:,:)
    integer,             pointer :: ncel_x, ncel_y, ncel_z
    integer(int2),       pointer :: cell_g2l(:), cell_g2b(:)
    integer,             pointer :: natom(:), nsolute(:), nwater(:)
    integer,             pointer :: solute_list(:,:), water_list(:,:,:)
    integer,             pointer :: No_HGr(:), HGr_local(:,:)
    integer,             pointer :: HGr_bond_list(:,:,:,:)


    ncel          = domain%num_cell_local + domain%num_cell_boundary
    ncel_local    = domain%num_cell_local

    call alloc_constraints(constraints, ConstraintsDomainBond, ncel, &
                           constraints%connect)

    cell_g2l      => domain%cell_g2l
    cell_g2b      => domain%cell_g2b
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water
    solute_list   => domain%solute_list
    water_list    => domain%water_list
    cell_pbc_move => domain%cell_pbc_move

    No_HGr        => constraints%No_HGr
    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list

    bsize(1)      = boundary%box_size_x
    bsize(2)      = boundary%box_size_y
    bsize(3)      = boundary%box_size_z
    bsize_orig(1) = boundary%box_size_orig(1)
    bsize_orig(2) = boundary%box_size_orig(2)
    bsize_orig(3) = boundary%box_size_orig(3)
    csize(1)      = boundary%cell_size_x
    csize(2)      = boundary%cell_size_y
    csize(3)      = boundary%cell_size_z
    num_cell(1)   = boundary%num_cells_x
    num_cell(2)   = boundary%num_cells_y
    num_cell(3)   = boundary%num_cells_z

    origin(1)     = boundary%origin_x
    origin(2)     = boundary%origin_y
    origin(3)     = boundary%origin_z

    natom  (1:ncel) = 0
    nsolute(1:ncel) = 0
    nwater (1:ncel) = 0
    No_HGr (1:ncel) = 0
    do i = 1, constraints%connect
      HGr_local(i,1:ncel) = 0
    end do

    ! FEP: make id list for single part
    if (domain%fep_use) then
      call setup_fep_correspondence(molecule, domain)
    end if

    ! solute atoms (not bonded to hydrogen) in each domain
    !
    do dupl_x = 1, boundary%num_duplicate(1)
      dupl(1) = dupl_x
      do dupl_y = 1, boundary%num_duplicate(2)
        dupl(2) = dupl_y
        do dupl_z = 1, boundary%num_duplicate(3)
          dupl(3) = dupl_z

          ioffset = dupl(1)-1 + (dupl(2)-1)*boundary%num_duplicate(1) &
                  + (dupl(3)-1)*boundary%num_duplicate(1)*boundary%num_duplicate(2)

          ioffset = ioffset * molecule%num_atoms

          do isolute = 1, enefunc%table%num_solute

            i   = enefunc%table%solute_list(isolute)

            ! FEP: skip singleB atoms
            if (domain%fep_use) then
              if (molecule%fepgrp(i) == 2) cycle
            end if

            ci1 = molecule%atom_name(i)

            mi1 = molecule%light_atom_mass(i)
            cl1 = molecule%light_atom_name(i)
            if (constraints%hydrogen_type == ConstraintAtomMass) then
              cl1 = mi1
            else if (constraints%hydrogen_type == ConstraintAtomBoth) then
              cl1 = (cl1 .or. mi1)
            end if

            if (constraints%duplicate(i) == 0 .and. .not. cl1) then

              !coordinate shifted against the origin
              !
              shift(1:3) = molecule%atom_coord(1:3,i) &
                         + bsize_orig(1:3)*real(dupl(1:3)-1,wip) - origin(1:3)

              !coordinate shifted to the first quadrant and set into the boundary box
              !
              move(1:3)  = bsize(1:3)*0.5_wip &
                         - bsize(1:3)*anint(shift(1:3)/bsize(1:3))
              shift(1:3) = shift(1:3) + move(1:3)

              !assign which cell
              !
              ic(1:3) = int(shift(1:3)/csize(1:3))
              if (ic(1) == num_cell(1)) ic(1) = ic(1) - 1
              if (ic(2) == num_cell(2)) ic(2) = ic(2) - 1
              if (ic(3) == num_cell(3)) ic(3) = ic(3) - 1
              icel = 1 + ic(1) + ic(2)*num_cell(1) + ic(3)*num_cell(1)*num_cell(2)

              ! atoms inside the domain
              !
              if (cell_g2l(icel) /= 0) then

                ! local cell index
                !
                icel_local = cell_g2l(icel)

                patm =  natom  (icel_local)
                psol =  nsolute(icel_local)
                pnoh =  No_HGr (icel_local)

                ! local_count : total number of atoms in each cell
                !
                patm = patm + 1
                psol = psol + 1
                pnoh = pnoh + 1

                solute_list(psol,icel_local) = patm

                call molecule_to_domain(molecule, move, origin, i,   &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)

                natom  (icel_local) = patm
                nsolute(icel_local) = psol
                No_HGr (icel_local) = pnoh

              ! atoms in the boundary
              !
              else if (cell_g2b(icel) /= 0) then

                ! local cell index
                !
                icel_local = cell_g2b(icel) + ncel_local

                patm =  natom  (icel_local)
                psol =  nsolute(icel_local)
                pnoh =  No_HGr (icel_local)

                ! local_count : total number of atoms in each cell
                !
                patm = patm + 1
                psol = psol + 1
                pnoh = pnoh + 1

                solute_list(psol,icel_local) = patm

                if (domain%nonbond_kernel == NBK_Fugaku .or.           &
                    domain%nonbond_kernel == NBK_Intel)                &
                  move(1:3) = move(1:3)                                &
                            + cell_pbc_move(1:3,icel_local)*bsize(1:3)
                call molecule_to_domain(molecule, move, origin, i,     &
                                        domain, icel_local, patm,      &
                                        dupl, bsize_orig, ioffset)

                natom  (icel_local) = patm
                nsolute(icel_local) = psol
                No_HGr (icel_local) = pnoh

              end if

            end if

          end do

        end do
      end do
    end do

    ! Hydrogen bonding group
    !
    do dupl_x = 1, boundary%num_duplicate(1)
      dupl(1) = dupl_x
      do dupl_y = 1, boundary%num_duplicate(2)
        dupl(2) = dupl_y
        do dupl_z = 1, boundary%num_duplicate(3)
          dupl(3) = dupl_z

          ioffset = dupl(1)-1 + (dupl(2)-1)*boundary%num_duplicate(1) &
                  + (dupl(3)-1)*boundary%num_duplicate(1)*boundary%num_duplicate(2)

          ioffset = ioffset * molecule%num_atoms

          do isolute = 1, constraints%connect

            do j = 1, constraints%nh(isolute)

              i = constraints%H_Group(1,j,isolute)

              ! FEP: skip singleB atoms
              if (domain%fep_use) then
                if (molecule%fepgrp(i) == 2) cycle
              end if

              coord_temp(1:3) = molecule%atom_coord(1:3,i)
              shift(1:3) = molecule%atom_coord(1:3,i) &
                         + bsize_orig(1:3)*real(dupl(1:3)-1,wip) - origin(1:3)

              move(1:3)  = bsize(1:3)*0.5_wip &
                         - bsize(1:3)*anint(shift(1:3)/bsize(1:3))
              shift(1:3) = shift(1:3) + move(1:3)

              !assign which cell
              !
              ic(1:3) = int(shift(1:3)/csize(1:3))
              if (ic(1) == num_cell(1)) ic(1) = ic(1) - 1
              if (ic(2) == num_cell(2)) ic(2) = ic(2) - 1
              if (ic(3) == num_cell(3)) ic(3) = ic(3) - 1
              icel = 1 + ic(1) + ic(2)*num_cell(1) + ic(3)*num_cell(1)*num_cell(2)

              ! atoms inside the domain
              !
              if (cell_g2l(icel) /= 0) then

                icel_local = cell_g2l(icel)

                patm =  natom  (icel_local)
                psol =  nsolute(icel_local)
                phgl =  HGr_local(isolute,icel_local)

                ! hydrogen atom
                !
                patm = patm + 1
                psol = psol + 1
                phgl = phgl + 1

                solute_list(psol,icel_local) = patm
                HGr_bond_list(1,phgl,isolute,icel_local) = patm

                call molecule_to_domain(molecule, move, origin, i,   &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)

                HGr_local(isolute,icel_local) = phgl

                do ih1 = 1, isolute

                  i = constraints%H_Group(ih1+1,j,isolute)

                  ! FEP: skip singleB atoms
                  if (domain%fep_use) then
                    if (molecule%fepgrp(i) == 2) cycle
                  end if

                  do k = 1, 3
                    if (molecule%atom_coord(k,i) < &
                           coord_temp(k) - 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,i) = &
                             molecule%atom_coord(k,i) + bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,i) =           &
                             molecule%atom_refcoord(k,i) + bsize(k)
                    else if (molecule%atom_coord(k,i) > &
                           coord_temp(k) + 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,i) = &
                             molecule%atom_coord(k,i) - bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,i) =           &
                             molecule%atom_refcoord(k,i) - bsize(k)
                    end if
                  end do

                  patm = patm + 1
                  psol = psol + 1

                  solute_list(psol,icel_local) = patm
                  HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm
      
                  call molecule_to_domain(molecule, move, origin, i,   &
                                          domain, icel_local, patm,    &
                                          dupl, bsize_orig, ioffset)

                end do

                natom  (icel_local) = patm
                nsolute(icel_local) = psol

              else if (cell_g2b(icel) /= 0) then

                icel_local = cell_g2b(icel) + ncel_local

                patm =  natom  (icel_local)
                psol =  nsolute(icel_local)
                phgl =  HGr_local(isolute,icel_local)

                ! hydrogen atoms
                !
                patm = patm + 1
                psol = psol + 1
                phgl = phgl + 1

                solute_list(psol,icel_local) = patm
                HGr_bond_list(1,phgl,isolute,icel_local) = patm
    
                if (domain%nonbond_kernel == NBK_Fugaku .or.           &
                    domain%nonbond_kernel == NBK_Intel)                &
                  move(1:3) = move(1:3)                                &
                            + cell_pbc_move(1:3,icel_local)*bsize(1:3)
 
                call molecule_to_domain(molecule, move, origin, i,   &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)
 
                HGr_local(isolute,icel_local) = phgl

                do ih1 = 1, isolute

                  i = constraints%H_Group(ih1+1,j,isolute)

                  ! FEP: skip singleB atoms
                  if (domain%fep_use) then
                    if (molecule%fepgrp(i) == 2) cycle
                  end if

                  do k = 1, 3
                    if (molecule%atom_coord(k,i) < &
                           coord_temp(k) - 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,i) = &
                             molecule%atom_coord(k,i) + bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,i) =           &
                             molecule%atom_refcoord(k,i) + bsize(k)
                    else if (molecule%atom_coord(k,i) > &
                           coord_temp(k) + 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,i) = &
                             molecule%atom_coord(k,i) - bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,i) =           &
                             molecule%atom_refcoord(k,i) - bsize(k)
                    end if
                  end do

                  patm = patm + 1
                  psol = psol + 1

                  solute_list(psol,icel_local) = patm
                  HGr_bond_list(ih1+1,phgl,isolute,icel_local) = patm

                  call molecule_to_domain(molecule, move, origin, i,   &
                                          domain, icel_local, patm,    &
                                          dupl, bsize_orig, ioffset)

                end do

                natom  (icel_local) = patm
                nsolute(icel_local) = psol

              end if

            end do

          end do

        end do
      end do
    end do

    ! water atoms in each domain
    !
    do dupl_x = 1, boundary%num_duplicate(1)
      dupl(1) = dupl_x
      do dupl_y = 1, boundary%num_duplicate(2)
        dupl(2) = dupl_y
        do dupl_z = 1, boundary%num_duplicate(3)
          dupl(3) = dupl_z

          ioffset = dupl(1)-1 + (dupl(2)-1)*boundary%num_duplicate(1) &
                  + (dupl(3)-1)*boundary%num_duplicate(1)*boundary%num_duplicate(2)

          ioffset = ioffset * molecule%num_atoms

          do iwater = 1, enefunc%table%num_water

            i   = enefunc%table%water_list(1,iwater)
            if (constraints%water_type == TIP2) then
              ih1 = enefunc%table%water_list(2,iwater)
            end if
            if (constraints%water_type == TIP3) then
              ih1 = enefunc%table%water_list(2,iwater)
              ih2 = enefunc%table%water_list(3,iwater)
            end if
            if (constraints%water_type == TIP4) then
              ih1 = enefunc%table%water_list(2,iwater)
              ih2 = enefunc%table%water_list(3,iwater)
              id = enefunc%table%water_list(4,iwater)
            end if

            coord_temp(1:3) = molecule%atom_coord(1:3,i)

            !coordinate shifted against the origin
            !
            shift(1:3) = molecule%atom_coord(1:3,i) &
                       + bsize_orig(1:3)*real(dupl(1:3)-1,wip) - origin(1:3)

            !coordinate shifted to the first quadrant and set into the boundary box
            !
            move(1:3)  = bsize(1:3)*0.5_wip &
                       - bsize(1:3)*anint(shift(1:3)/bsize(1:3))
            shift(1:3) = shift(1:3) + move(1:3)

            !assign which cell
            !
            ic(1:3) = int(shift(1:3)/csize(1:3))
            if (ic(1) == num_cell(1)) ic(1) = ic(1) - 1
            if (ic(2) == num_cell(2)) ic(2) = ic(2) - 1
            if (ic(3) == num_cell(3)) ic(3) = ic(3) - 1
            icel = 1 + ic(1) + ic(2)*num_cell(1) + ic(3)*num_cell(1)*num_cell(2)

            ! atoms inside the domain
            !
            if (cell_g2l(icel) /= 0) then

              ! local cell index
              !
              icel_local = cell_g2l(icel)
              patm =  natom (icel_local)
              pwat =  nwater(icel_local)
              pwat = pwat + 1

              ! oxygen atoms
              !
              patm = patm + 1
              water_list(1,pwat,icel_local) = patm
              call molecule_to_domain(molecule, move, origin, i, &
                                      domain, icel_local, patm,  &
                                      dupl, bsize_orig, ioffset)

              if (constraints%water_type == TIP2 .or. &
                  constraints%water_type == TIP3 .or. &
                  constraints%water_type == TIP4) then

                ! first hydrogen atoms
                !
                patm = patm + 1
                water_list(2,pwat,icel_local) = patm
                do k = 1, 3
                  if (molecule%atom_coord(k,ih1) < &
                         coord_temp(k) - 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih1) = &
                           molecule%atom_coord(k,ih1) + bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih1) =           &
                           molecule%atom_refcoord(k,ih1) + bsize(k)
                  else if (molecule%atom_coord(k,ih1) > &
                         coord_temp(k) + 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih1) = &
                           molecule%atom_coord(k,ih1) - bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih1) =           &
                           molecule%atom_refcoord(k,ih1) - bsize(k)
                  end if
                end do
                call molecule_to_domain(molecule, move, origin, ih1, &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)

                if (constraints%water_type == TIP3 .or. &
                    constraints%water_type == TIP4) then

                  ! second hydrogen atoms
                  !
                  patm = patm + 1
                  water_list(3,pwat,icel_local) = patm
                  do k = 1, 3
                    if (molecule%atom_coord(k,ih2) < &
                           coord_temp(k) - 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,ih2) = &
                             molecule%atom_coord(k,ih2) + bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,ih2) =           &
                             molecule%atom_refcoord(k,ih2) + bsize(k)
                    else if (molecule%atom_coord(k,ih2) > &
                           coord_temp(k) + 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,ih2) = &
                             molecule%atom_coord(k,ih2) - bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,ih2) =           &
                             molecule%atom_refcoord(k,ih2) - bsize(k)
                    end if
                  end do
                  call molecule_to_domain(molecule, move, origin, ih2, &
                                          domain, icel_local, patm,    &
                                          dupl, bsize_orig, ioffset)

                end if

                ! dummy atoms
                !
                if (constraints%water_type == TIP4) then
                  patm = patm + 1
                  water_list(4,pwat,icel_local) = patm
                  do k = 1, 3
                    if (molecule%atom_coord(k,id) < &
                           coord_temp(k) - 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,id) = &
                             molecule%atom_coord(k,id) + bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,id) =           &
                             molecule%atom_refcoord(k,id) + bsize(k)
                    else if (molecule%atom_coord(k,id) > &
                           coord_temp(k) + 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,id) = &
                             molecule%atom_coord(k,id) - bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,id) =           &
                             molecule%atom_refcoord(k,id) - bsize(k)
                    end if
                  end do
                  call molecule_to_domain(molecule, move, origin, id, &
                                          domain, icel_local, patm,   &
                                          dupl, bsize_orig, ioffset)
                end if

              end if
              natom(icel_local)  = patm
              nwater(icel_local) = pwat

            ! atoms in the boundary
            !
            else if (cell_g2b(icel) /= 0) then

              ! local cell index
              !
              icel_local = cell_g2b(icel) + ncel_local

              patm =  natom (icel_local)
              pwat =  nwater(icel_local)
              pwat = pwat + 1

              ! oxygen atoms
              !
              patm = patm + 1
              water_list(1,pwat,icel_local) = patm

              if (domain%nonbond_kernel == NBK_Fugaku .or.           &
                  domain%nonbond_kernel == NBK_Intel)                &
                move(1:3) = move(1:3)                                &
                          + cell_pbc_move(1:3,icel_local)*bsize(1:3)

              call molecule_to_domain(molecule, move, origin, i, &
                                      domain, icel_local, patm,  &
                                      dupl, bsize_orig, ioffset)


              if (constraints%water_type == TIP2 .or. &
                  constraints%water_type == TIP3 .or. &
                  constraints%water_type == TIP4) then

                ! first hydrogen atoms
                !
                patm = patm + 1
                water_list(2,pwat,icel_local) = patm
                do k = 1, 3
                  if (molecule%atom_coord(k,ih1) < &
                         coord_temp(k) - 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih1) = &
                           molecule%atom_coord(k,ih1) + bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih1) =           &
                           molecule%atom_refcoord(k,ih1) + bsize(k)
                  else if (molecule%atom_coord(k,ih1) > &
                         coord_temp(k) + 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih1) = &
                           molecule%atom_coord(k,ih1) - bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih1) =           &
                           molecule%atom_refcoord(k,ih1) - bsize(k)
                  end if
                end do
                call molecule_to_domain(molecule, move, origin, ih1, &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)
  
              if (constraints%water_type == TIP3 .or. &
                  constraints%water_type == TIP4) then

                ! second hydrogen atoms
                !
                patm = patm + 1
                water_list(3,pwat,icel_local) = patm
                do k = 1, 3
                  if (molecule%atom_coord(k,ih2) < &
                         coord_temp(k) - 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih2) = &
                           molecule%atom_coord(k,ih2) + bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih2) =           &
                           molecule%atom_refcoord(k,ih2) + bsize(k)
                  else if (molecule%atom_coord(k,ih2) > &
                       coord_temp(k) + 0.5_wp*bsize(k)) then
                    molecule%atom_coord(k,ih2) = &
                           molecule%atom_coord(k,ih2) - bsize(k)
                    if (allocated(molecule%atom_refcoord)) &
                      molecule%atom_refcoord(k,ih2) =           &
                           molecule%atom_refcoord(k,ih2) - bsize(k)
                  end if
                end do
                call molecule_to_domain(molecule, move, origin, ih2, &
                                        domain, icel_local, patm,    &
                                        dupl, bsize_orig, ioffset)
             end if
  
                ! dummy atoms
                !
                if (constraints%water_type == TIP4) then
                  patm = patm + 1
                  water_list(4,pwat,icel_local) = patm
                  do k = 1, 3
                    if (molecule%atom_coord(k,id) < &
                           coord_temp(k) - 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,id) = &
                             molecule%atom_coord(k,id) + bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,id) =           &
                             molecule%atom_refcoord(k,id) + bsize(k)
                    else if (molecule%atom_coord(k,id) > &
                           coord_temp(k) + 0.5_wp*bsize(k)) then
                      molecule%atom_coord(k,id) = &
                             molecule%atom_coord(k,id) - bsize(k)
                      if (allocated(molecule%atom_refcoord)) &
                        molecule%atom_refcoord(k,id) =           &
                             molecule%atom_refcoord(k,id) - bsize(k)
                    end if
                  end do
                  call molecule_to_domain(molecule, move, origin, id, &
                                          domain, icel_local, patm,   &
                                          dupl, bsize_orig, ioffset)
                end if

              end if

              natom(icel_local)  = patm
              nwater(icel_local) = pwat

            end if

          end do

        end do
      end do
    end do

    return

  end subroutine setup_atom_by_HBond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atom_by_HBond_pio
  !> @brief        setup atom maps with H-bond connection groups
  !! @authors      JJ
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atom_by_HBond_pio(boundary, enefunc, constraints, domain)

    ! formal arguments
    type(s_boundary),    target, intent(in)    :: boundary
    type(s_enefunc),     target, intent(in)    :: enefunc
    type(s_constraints), target, intent(inout) :: constraints
    type(s_domain),      target, intent(inout) :: domain

    ! local variable
    real(wip)                    :: shift(3), bsize(3)
    real(wip)                    :: move(3)
    integer                      :: file_num, file_tot_num
    integer                      :: i, j, k, ix, list, icel, ic
    integer                      :: icel_local
    integer                      :: ncel_local, ncel
    integer                      :: ncell_local_pio, ncell_pio
    integer                      :: water_atom

    integer,             pointer :: natom(:), nsolute(:), nwater(:)
    integer,             pointer :: water_list(:,:,:), water_list_pio(:,:,:,:)
    integer,             pointer :: cell_l2g_pio(:,:)
    integer,             pointer :: natom_pio(:,:)
    integer,             pointer :: nsolute_pio(:,:), nwater_pio(:,:)
    integer,             pointer :: solute_list(:,:)
    integer,             pointer :: id_l2g(:,:), id_l2g_sol(:,:)
    integer(int2),       pointer :: id_g2l(:,:), cell_g2l(:), cell_g2b(:)
    integer,             pointer :: id_l2g_pio(:,:,:), id_l2g_sol_pio(:,:,:)
    integer,             pointer :: No_HGr(:), HGr_local(:,:)
    integer,             pointer :: HGr_bond_list(:,:,:,:)
    integer,             pointer :: No_HGr_pio(:,:), HGr_local_pio(:,:,:)
    integer,             pointer :: HGr_bond_list_pio(:,:,:,:,:)
    integer,             pointer :: aclass(:,:), aclass_pio(:,:,:)
    integer,             pointer :: ring(:,:), ring_pio(:,:,:)
    real(wp),            pointer :: trans(:,:,:)
    real(wp),            pointer :: charge(:,:), charge_pio(:,:,:)
    real(wp),            pointer :: cell_pbc_move(:,:)
    real(wip),           pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),           pointer :: coord_pio(:,:,:,:)
    real(wip),           pointer :: vel(:,:,:), vel_pio(:,:,:,:)
    real(wip),           pointer :: mass(:,:), mass_pio(:,:,:), inv_mass(:,:)
    real(wip),           pointer :: HGr_bond_dist(:,:,:,:)
    real(wip),           pointer :: HGr_bond_dist_pio(:,:,:,:,:)


    ncel            = domain%num_cell_local + domain%num_cell_boundary
    ncel_local      = domain%num_cell_local
    ncell_pio       = domain%ncell_pio
    ncell_local_pio = domain%ncell_local_pio

    call alloc_constraints(constraints, ConstraintsDomainBond, ncel, &
                           constraints%connect)

    cell_g2l        => domain%cell_g2l
    cell_g2b        => domain%cell_g2b
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    solute_list     => domain%solute_list
    water_list      => domain%water_list
    cell_l2g_pio    => domain%cell_l2g_pio
    natom_pio       => domain%num_atom_pio
    nsolute_pio     => domain%num_solute_pio
    nwater_pio      => domain%num_water_pio
    water_list_pio  => domain%water_list_pio
    id_l2g          => domain%id_l2g
    id_l2g_sol      => domain%id_l2g_solute
    id_g2l          => domain%id_g2l
    id_l2g_pio      => domain%id_l2g_pio
    id_l2g_sol_pio  => domain%id_l2g_solute_pio
    trans           => domain%trans_vec
    coord           => domain%coord
    vel             => domain%velocity
    coord_ref       => domain%coord_ref
    mass            => domain%mass
    charge          => domain%charge
    aclass          => domain%atom_cls_no
    aclass_pio      => domain%atom_cls_no_pio
    ring            => domain%ring
    ring_pio        => domain%ring_pio
    coord_pio       => domain%coord_pio
    vel_pio         => domain%velocity_pio
    mass_pio        => domain%mass_pio
    charge_pio      => domain%charge_pio
    inv_mass        => domain%inv_mass
    cell_pbc_move   => domain%cell_pbc_move

    No_HGr          => constraints%No_HGr
    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list
    HGr_bond_dist   => constraints%HGr_bond_dist
    No_HGr_pio      => constraints%No_HGr_pio
    HGr_local_pio   => constraints%HGr_local_pio
    HGr_bond_list_pio => constraints%HGr_bond_list_pio
    HGr_bond_dist_pio => constraints%HGr_bond_dist_pio

    bsize(1)      = boundary%box_size_x
    bsize(2)      = boundary%box_size_y
    bsize(3)      = boundary%box_size_z

    natom  (1:ncel) = 0
    nsolute(1:ncel) = 0
    nwater (1:ncel) = 0
    No_HGr (1:ncel) = 0
    do i = 1, constraints%connect
      HGr_local(i,1:ncel) = 0
    end do

    file_num = 0
    file_tot_num = boundary%multiple_file(1) &
                  *boundary%multiple_file(2) &
                  *boundary%multiple_file(3)
    domain%file_tot_num = file_tot_num

    do file_num = 1, file_tot_num

      do icel = 1, ncell_pio

        ic = cell_l2g_pio(icel,file_num)

        if (ic /= 0) then
        if (cell_g2l(ic) /= 0) then

          i = cell_g2l(ic)

          natom  (i) = natom_pio  (icel,file_num)
          nsolute(i) = nsolute_pio(icel,file_num)
          nwater (i) = nwater_pio (icel,file_num)
          No_HGr (i) = No_HGr_pio (icel,file_num)

          do ix = 1, natom(i)
            coord    (1:3,ix,i) = coord_pio (1:3,ix,icel,file_num)
            coord_ref(1:3,ix,i) = coord_pio (1:3,ix,icel,file_num)
            vel      (1:3,ix,i) = vel_pio   (1:3,ix,icel,file_num)
            charge   (    ix,i) = charge_pio(    ix,icel,file_num)
            aclass   (    ix,i) = aclass_pio(    ix,icel,file_num)
            ring     (    ix,i) = ring_pio  (    ix,icel,file_num)
            mass     (    ix,i) = mass_pio  (    ix,icel,file_num)
            id_l2g   (    ix,i) = id_l2g_pio(    ix,icel,file_num)
            inv_mass (    ix,i) = 0.0_wip
            if (mass(ix,i) > EPS) inv_mass(ix,i) = 1.0_wip / mass(ix,i)
          end do

          do ix = 1, nsolute(i)
            id_l2g_sol(ix,i) = id_l2g_sol_pio(ix,icel,file_num)
            id_g2l(1,id_l2g_sol(ix,i)) = i
            id_g2l(2,id_l2g_sol(ix,i)) = ix
          end do

          do ix = 1, nwater(i)
            if (constraints%water_type == TIP4) then
              water_list(1:4,ix,i) = water_list_pio(1:4,ix,icel,file_num)
            else if (constraints%water_type == TIP3) then
              water_list(1:3,ix,i) = water_list_pio(1:3,ix,icel,file_num)
            else if (constraints%water_type == TIP2) then
              water_list(1:2,ix,i) = water_list_pio(1:2,ix,icel,file_num)
            else if (constraints%water_type == TIP1) then
              water_list(1:1,ix,i) = water_list_pio(1:1,ix,icel,file_num)
            end if
          end do

          do ix = 1, No_HGr(i)
            trans(1:3,ix,i) = bsize(1:3)*0.5_wip &
                            - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
          end do

          k = No_HGr(i)
          do j = 1, constraints%connect
            HGr_local(j,i) = HGr_local_pio(j,icel,file_num)
            do ix = 1, HGr_local(j,i)
              HGr_bond_list(1:j+1,ix,j,i) = &
                            HGr_bond_list_pio(1:j+1,ix,j,icel,file_num)
              if (icel <= ncell_local_pio)     & 
              HGr_bond_dist(1:j+1, ix, j, i) = &
                            HGr_bond_dist_pio(1:j+1,ix,j,icel,file_num)
              list = HGr_bond_list(1,ix,j,i)
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,list,i)/bsize(1:3))
              do list = 1, j+1
                k = k + 1
                trans(1:3,k,i) = move(1:3)
              end do
            end do
          end do

          if (constraints%water_type == TIP4) then
            do ix = nsolute(i)+1, natom(i)-3, 4
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              do j = 1, 4
                trans(1:3,ix-1+j,i) = move(1:3)
              end do
            end do
          else if (constraints%water_type == TIP3) then
            do ix = nsolute(i)+1, natom(i)-2, 3
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              do j = 1, 3
                trans(1:3,ix-1+j,i) = move(1:3)
              end do
            end do
          else if (constraints%water_type == TIP2) then
               do ix = nsolute(i)+1, natom(i)-2, 2
                  move(1:3) = bsize(1:3)*0.5_wip &
                       - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
                  do j = 1, 2
                     trans(1:3,ix-1+j,i) = move(1:3)
                  end do
               end do
          else if (constraints%water_type == TIP1) then
            do ix = nsolute(i)+1, natom(i)
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              trans(1:3,ix,i) = move(1:3)
            end do
          end if

        else if (cell_g2b(ic) /= 0) then

          i = cell_g2b(ic) + ncel_local

          natom  (i) = natom_pio  (icel,file_num)
          nsolute(i) = nsolute_pio(icel,file_num)
          nwater (i) = nwater_pio (icel,file_num)
          No_HGr (i) = No_HGr_pio (icel,file_num)

          do ix = 1, natom(i)
            coord    (1:3,ix,i) = coord_pio (1:3,ix,icel,file_num)
            coord_ref(1:3,ix,i) = coord_pio (1:3,ix,icel,file_num)
            vel      (1:3,ix,i) = vel_pio   (1:3,ix,icel,file_num)
            charge   (    ix,i) = charge_pio(    ix,icel,file_num)
            aclass   (    ix,i) = aclass_pio(    ix,icel,file_num)
            ring     (    ix,i) = ring_pio  (    ix,icel,file_num)
            mass     (    ix,i) = mass_pio  (    ix,icel,file_num)
            id_l2g   (    ix,i) = id_l2g_pio(    ix,icel,file_num)
            inv_mass (    ix,i) = 0.0_wip
            if (mass(ix,i) > EPS) inv_mass(ix,i) = 1.0_wip / mass(ix,i)
          end do

          do ix = 1, nsolute(i)
            id_l2g_sol(ix,i) = id_l2g_sol_pio(ix,icel,file_num)
            id_g2l(1,id_l2g_sol(ix,i)) = i
            id_g2l(2,id_l2g_sol(ix,i)) = ix
          end do

          do ix = 1, nwater(i)
            if (constraints%water_type == TIP4) then
              water_list(1:4,ix,i) = water_list_pio(1:4,ix,icel,file_num)
            else if (constraints%water_type == TIP3) then
              water_list(1:3,ix,i) = water_list_pio(1:3,ix,icel,file_num)
            else if (constraints%water_type == TIP2) then
              water_list(1:2,ix,i) = water_list_pio(1:2,ix,icel,file_num)
            else if (constraints%water_type == TIP1) then
              water_list(1:1,ix,i) = water_list_pio(1:1,ix,icel,file_num)
            end if
          end do

          if (domain%nonbond_kernel == NBK_Fugaku .or. &
              domain%nonbond_kernel == NBK_Intel) then
            do ix = 1, No_HGr(i)
              trans(1:3,ix,i) = bsize(1:3)*0.5_wip                           &
                              - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3)) &
                              + bsize(1:3)*cell_pbc_move(1:3,i)
            end do
          else
            do ix = 1, No_HGr(i)
              trans(1:3,ix,i) = bsize(1:3)*0.5_wip &
                              - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
            end do
          end if

          k = No_HGr(i)
          do j = 1, constraints%connect
            HGr_local(j,i) = HGr_local_pio(j,icel,file_num)
            do ix = 1, HGr_local(j,i)
              HGr_bond_list(1:j+1,ix,j,i) = &
                            HGr_bond_list_pio(1:j+1,ix,j,icel,file_num)
              if (icel <= ncell_local_pio)     & 
              HGr_bond_dist(1:j+1, ix, j, i) = &
                            HGr_bond_dist_pio(1:j+1,ix,j,icel,file_num)
              do list = 2,j+1
              end do
              list = HGr_bond_list(1,ix,j,i)
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,list,i)/bsize(1:3))
              if (domain%nonbond_kernel == NBK_Fugaku .or. &
                  domain%nonbond_kernel == NBK_Intel)      &
                move(1:3) = move(1:3) + cell_pbc_move(1:3,i)*bsize(1:3)
              do list = 1, j+1
                k = k + 1
                trans(1:3,k,i) = move(1:3)
              end do
            end do
          end do


          if (constraints%water_type == TIP4) then
            do ix = nsolute(i)+1, natom(i)-3, 4
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              if (domain%nonbond_kernel == NBK_Fugaku) &
                move(1:3) = move(1:3) + cell_pbc_move(1:3,i)*bsize(1:3)
              do j = 1, 4
                trans(1:3,ix-1+j,i) = move(1:3)
              end do
            end do
          else if (constraints%water_type == TIP3) then
            do ix = nsolute(i)+1, natom(i)-2, 3
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              if (domain%nonbond_kernel == NBK_Fugaku) &
                move(1:3) = move(1:3) + cell_pbc_move(1:3,i)*bsize(1:3)
              do j = 1, 3
                trans(1:3,ix-1+j,i) = move(1:3)
              end do
            end do
          else if (constraints%water_type == TIP2) then
            do ix = nsolute(i)+1, natom(i)-1, 2
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              if (domain%nonbond_kernel == NBK_Fugaku) &
                move(1:3) = move(1:3) + cell_pbc_move(1:3,i)*bsize(1:3)
              do j = 1, 2
                trans(1:3,ix-1+j,i) = move(1:3)
              end do
            end do
          else if (constraints%water_type == TIP1) then
            do ix = nsolute(i)+1, natom(i)
              move(1:3) = bsize(1:3)*0.5_wip &
                        - bsize(1:3)*anint(coord(1:3,ix,i)/bsize(1:3))
              if (domain%nonbond_kernel == NBK_Fugaku) &
                move(1:3) = move(1:3) + cell_pbc_move(1:3,i)*bsize(1:3)
              trans(1:3,ix,i) = move(1:3)
            end do
          end if

        end if
        end if

      end do

    end do

    call dealloc_domain(domain, DomainDynvar_Atom_pio)     
    
    return

  end subroutine setup_atom_by_HBond_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_global_to_local_atom_index
  !> @brief        relationship between global and local indices
  !! @authors      JJ
  !! @param[in]    enefunc  : energy function information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_global_to_local_atom_index(enefunc, domain)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, ix, ig, ig1, ig2, is
    integer                   :: ncel, natom_dupl

    integer(2),       pointer :: id_g2l(:,:)
    integer,          pointer :: id_l2g(:,:), id_l2g_solute(:,:)
    integer,          pointer :: natom(:), nsolute(:)
    integer,          pointer :: sollist(:)


    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g
    id_l2g_solute => domain%id_l2g_solute
    natom         => domain%num_atom
    nsolute       => domain%num_solute
    sollist       => enefunc%table%solute_list_inv

    ncel         = domain%num_cell_local + domain%num_cell_boundary
    natom_dupl   = domain%num_atom_all / domain%num_duplicate

    ! local variable
    do i = 1, ncel
      do ix = 1, nsolute(i)
        ig  = id_l2g(ix,i)
        id_l2g_solute(ix,i) = ig
        id_g2l(1,ig) = i
        id_g2l(2,ig) = ix
      end do
    end do

    return

  end subroutine setup_global_to_local_atom_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_neighboring_cells
  !> @brief        check the neighboring cells of each cell
  !! @authors      JJ
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_neighbor_cells(boundary, domain)

    ! formal arguments
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: ii, jj, kk, ic, jc
    integer                   :: i, j, k
    integer                   :: inb, jnb, knb
    integer                   :: inbs, jnbs, knbs
    integer                   :: ncell_x, ncell_y, ncell_z
    integer                   :: ncel_local, nboundary

    integer,          pointer :: cell_l2g(:), cell_b2g(:)
    integer(int2),    pointer :: cell_g2l(:), cell_g2b(:)
    integer,          pointer :: cell_l2gx(:), cell_l2gy(:), cell_l2gz(:)


    cell_g2l   => domain%cell_g2l
    cell_g2b   => domain%cell_g2b
    cell_b2g   => domain%cell_b2g
    cell_l2g   => domain%cell_l2g
    cell_l2gx  => domain%cell_l2gx
    cell_l2gy  => domain%cell_l2gy
    cell_l2gz  => domain%cell_l2gz

    ncel_local = domain%num_cell_local
    nboundary  = domain%num_cell_boundary

    ncell_x    = boundary%num_cells_x
    ncell_y    = boundary%num_cells_y
    ncell_z    = boundary%num_cells_z

    call alloc_domain(domain, DomainNeighbourCell, ncel_local+nboundary, 1, 1)

    do ii = 1, ncel_local

      kk = 0
      ic = cell_l2g(ii)
      i  = cell_l2gx(ii)-1
      j  = cell_l2gy(ii)-1
      k  = cell_l2gz(ii)-1

     do knb = k-1, k+1

        knbs = mod(knb+ncell_z,ncell_z)

        do jnb = j-1, j+1

          jnbs = mod(jnb+ncell_y,ncell_y)

          do inb = i-1, i+1

            inbs = mod(inb+ncell_x,ncell_x)

            jc = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            jj  = cell_g2l(jc)

            if (jj > ii) then
              kk = kk + 1
              domain%neighbor_cells(kk,ii) = jj
            else if (cell_g2b(jc) /= 0) then
              jj = cell_g2b(jc) + ncel_local
              kk = kk + 1
              domain%neighbor_cells(kk,ii) = jj
            end if

          end do
        end do
      end do

      domain%near_neighbor_cell_count(ii) = kk

      do knb = k-2, k+2

        knbs = mod(knb+ncell_z,ncell_z)

        do jnb = j-2, j+2

          jnbs = mod(jnb+ncell_y,ncell_y)

          do inb = i-2, i+2

            inbs = mod(inb+ncell_x,ncell_x)

            if (abs(inb-i)>=2 .or. abs(jnb-j)>=2 .or. abs(knb-k)>=2) then

              jc = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
              jj  = cell_g2l(jc)

              if (jj > ii) then
                kk = kk + 1
                domain%neighbor_cells(kk,ii) = jj
              else if (cell_g2b(jc) /= 0) then
                jj = cell_g2b(jc) + ncel_local
                kk = kk + 1
                domain%neighbor_cells(kk,ii) = jj
              end if

            end if

          end do
        end do
      end do

      domain%neighbor_cell_count(ii) = kk

    end do

    do ii = ncel_local+1, ncel_local+nboundary

      kk = 0
      ic = ii - ncel_local
      ic = cell_b2g(ic)
      i  = cell_l2gx(ii)-1
      i  = mod(i+ncell_x,ncell_x)
      j  = cell_l2gy(ii)-1
      j  = mod(j+ncell_y,ncell_y)
      k  = cell_l2gz(ii)-1
      k  = mod(k+ncell_z,ncell_z)

      do knb = k-1, k+1

        knbs = mod(knb+ncell_z,ncell_z)

        do jnb = j-1, j+1

          jnbs = mod(jnb+ncell_y,ncell_y)

          do inb = i-1, i+1

            inbs = mod(inb+ncell_x,ncell_x)

            jc = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            jj = cell_g2b(jc) + ncel_local

            if (jj > ii) then

              kk = kk + 1
              domain%neighbor_cells(kk,ii) = jj

            end if

          end do
        end do
      end do

      domain%near_neighbor_cell_count(ii) = kk

      do knb = k-2, k+2

        knbs = mod(knb+ncell_z,ncell_z)

        do jnb = j-2, j+2

          jnbs = mod(jnb+ncell_y,ncell_y)

          do inb = i-2, i+2

            inbs = mod(inb+ncell_x,ncell_x)

            if (abs(inb-i)>=2 .or. abs(jnb-j)>=2 .or. abs(knb-k)>=2) then

              jc = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
              jj = cell_g2b(jc) + ncel_local

              if (jj > ii) then

                kk = kk + 1
                domain%neighbor_cells(kk,ii) = jj

              end if

            end if

          end do
        end do
      end do

      domain%neighbor_cell_count(ii) = kk

    end do

    return

  end subroutine assign_neighbor_cells

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_atoms
  !> @brief        compate the total number of atoms in each cell
  !! @authors      JJ
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_atoms(natom, num_atom, nsolute, nwater, &
                               cell_l2gx,  cell_l2gy, cell_l2gz, &
                               cell_start, cell_end,  neighbor,  &
                               ncel_local, nboundary)

    ! formal arguments
    real(wp),                intent(inout) :: natom(:)
    integer,                 intent(in)    :: num_atom(:)
    integer,                 intent(in)    :: nsolute(:)
    integer,                 intent(in)    :: nwater(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: neighbor(-1:1,-1:1,-1:1)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary

    ! local variables
    integer                  :: i, ii, i1, i2, i3


    do i = 1, ncel_local
      natom(i) = real(num_atom(i)*nproc_city + my_city_rank, wp)
    end do

    do ii = 1, nboundary

      i  = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)

      if (i1 == cell_start(1)-1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(-1,0,0),wp)
          end if
        end if
      else if (i1 == cell_end(1)+1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(1,0,0),wp)
          end if
        end if
      else if (i1 >= cell_start(1) .and. i1 <= cell_end(1)) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,-1,0),wp)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,1,0),wp)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,-1),wp)
          else if (i3 == cell_end(3)+1) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,1),wp)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            natom(i) = real(num_atom(i)*nproc_city+neighbor(0,0,0),wp)
          end if
        end if
      end if
    end do

    return

  end subroutine assign_cell_atoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_cpu
  !> @brief        compute the total cput time of each cell (randomized)
  !! @authors      JJ
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_cpu(time_proc, random,               &
                             cell_l2gx, cell_l2gy, cell_l2gz, &
                             cell_start, cell_end, neighbor,  &
                             ncel_local, nboundary, cpu_time)

    ! formal arguments
    real(wp),                intent(in)    :: time_proc(:)
    real(wip),               intent(in)    :: random(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: neighbor(-1:1,-1:1,-1:1)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary
    real(wp),                intent(inout) :: cpu_time(:)

    ! local variables
    integer                  :: i, ii, i1, i2, i3, j1, j2, j3


    do i = 1, ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)
      j1 = i1 + 1
      j2 = i2 + 1
      j3 = i3 + 1
      cpu_time(i) = time_proc(my_city_rank+1)*random(i)
    end do

    do ii = 1, nboundary

      i = ii + ncel_local
      i1 = cell_l2gx(i)
      i2 = cell_l2gy(i)
      i3 = cell_l2gz(i)
      j1 = i1 + 1
      j2 = i2 + 1
      j3 = i3 + 1

      if (i1 == cell_start(1)-1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(-1,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(-1,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(-1,0,0)+1)*random(i)
          end if
        end if
      else if (i1 == cell_end(1)+1) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(1,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(1,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(1,0,0)+1)*random(i)
          end if
        end if
      else if (i1 >= cell_start(1) .and. i1 <= cell_end(1)) then
        if (i2 == cell_start(2)-1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,-1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,-1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,-1,0)+1)*random(i)
          end if
        else if (i2 == cell_end(2)+1) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,1,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,1,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,1,0)+1)*random(i)
          end if
        else if (i2 >= cell_start(2) .and. i2 <= cell_end(2)) then
          if (i3 == cell_start(3)-1) then
            cpu_time(i) = time_proc(neighbor(0,0,-1)+1)*random(i)
          else if (i3 == cell_end(3)+1) then
            cpu_time(i) = time_proc(neighbor(0,0,1)+1)*random(i)
          else if (i3 >= cell_start(3) .and. i3 <= cell_end(3)) then
            cpu_time(i) = time_proc(neighbor(0,0,0)+1)*random(i)
          end if
        end if
      end if
    end do

    return

  end subroutine assign_cell_cpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_cell_interactions
  !> @brief        assign the cell index for given cell-cell interaction
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_cell_interaction(natom, num_domain, cell,            &
                                     cell_l2gx,  cell_l2gy,  cell_l2gz,  &
                                     cell_l2gx1, cell_l2gy1, cell_l2gz1, &
                                     cell_gxyz2l,cell_start, cell_end,   &
                                     ncel_local, nboundary,              &
                                     cell_move, cell_pair, virial_check)

    ! formal arguments
    real(wp),                intent(in)    :: natom(:)
    integer,                 intent(in)    :: num_domain(:)
    integer,                 intent(in)    :: cell(:)
    integer,                 intent(in)    :: cell_l2gx(:)
    integer,                 intent(in)    :: cell_l2gy(:)
    integer,                 intent(in)    :: cell_l2gz(:)
    integer,                 intent(in)    :: cell_l2gx1(:)
    integer,                 intent(in)    :: cell_l2gy1(:)
    integer,                 intent(in)    :: cell_l2gz1(:)
    integer(int2),           intent(in)    :: cell_gxyz2l(:,:,:)
    integer,                 intent(in)    :: cell_start(:)
    integer,                 intent(in)    :: cell_end(:)
    integer,                 intent(in)    :: ncel_local
    integer,                 intent(in)    :: nboundary
    integer(1),              intent(inout) :: cell_move(:,:,:)
    integer(int2),           intent(inout) :: cell_pair(:,:)
    integer(int1),           intent(inout) :: virial_check(:,:)

    ! local variables
    real(wip)                :: ic1, ic2
    integer                  :: i, i1, i2, i3, j, ii, jj, ij
    integer                  :: icx1, icy1, icz1, icx2, icy2, icz2
    integer                  :: icx, icy, icz, movex, movey, movez


    ij = 0

    ! assign the interaction cell for each interaction
    !
    do i = 1, ncel_local-1

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do j = i+1, ncel_local

        cell_pair(i,j) = 0
        cell_pair(j,i) = 0

        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then
          icx = (icx1 + icx2)/2
          icy = (icy1 + icy2)/2
          icz = (icz1 + icz2)/2
          cell_pair(i,j) = i
          cell_pair(j,i) = i

          icx =min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
          icy =min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
          icz =min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

          if (icx == abs(icx1-icx2-cell(1))) cell_move(1,j,i) = -1
          if (icx == abs(icx1-icx2+cell(1))) cell_move(1,j,i) =  1
          if (icx == abs(icx1-icx2))         cell_move(1,j,i) =  0
          if (icy == abs(icy1-icy2-cell(2))) cell_move(2,j,i) = -1
          if (icy == abs(icy1-icy2+cell(2))) cell_move(2,j,i) =  1
          if (icy == abs(icy1-icy2))         cell_move(2,j,i) =  0
          if (icz == abs(icz1-icz2-cell(3))) cell_move(3,j,i) = -1
          if (icz == abs(icz1-icz2+cell(3))) cell_move(3,j,i) =  1
          if (icz == abs(icz1-icz2))         cell_move(3,j,i) =  0
          if (cell_move(1,j,i) == 0 .and. &
              cell_move(2,j,i) == 0 .and. &
              cell_move(3,j,i) == 0) then
            virial_check(j,i) = 0
            virial_check(i,j) = 0
          end if

          cell_move(1:3,i,j) = - cell_move(1:3,j,i)

          ij = ij + 1
        end if

      end do
    end do

    do i = 1, ncel_local
      cell_pair(i,i) = i
      cell_move(1:3,i,i) = 0
      ij = ij + 1
    end do

    do i = 1, ncel_local

      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do jj = 1, nboundary

        j = jj + ncel_local
        cell_pair(j,i) = 0
        cell_pair(i,j) = 0

        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        movex = 0
        movey = 0
        movez = 0

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (icx == abs(icx1-icx2-cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = -1
          movex = -cell(1)
        else if (icx == abs(icx1-icx2+cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = 1
          movex = cell(1)
        else
          cell_move(1,j,i) = 0
          movex = 0
        end if

        if (icy == abs(icy1-icy2-cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = -1
          movey = -cell(2)
        else if (icy == abs(icy1-icy2+cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = 1
          movey = cell(2)
        else
          cell_move(2,j,i) = 0
          movey = 0
        end if

        if (icz == abs(icz1-icz2-cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = -1
          movez = -cell(3)
        else if (icz == abs(icz1-icz2+cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = 1
          movez = cell(3)
        else
          cell_move(3,j,i) = 0
          movez = 0
        end if

        cell_move(1:3,i,j) = - cell_move(1:3,j,i)

        icx = icx1 - icx2 + movex
        icy = icy1 - icy2 + movey
        icz = icz1 - icz2 + movez

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

          ij = ij + 1
          icx = icx1 + icx2 + movex
          icy = icy1 + icy2 + movey
          icz = icz1 + icz2 + movez

          if (icx == (2*cell_start(1)-1) .or. icx == (2*cell_end(1)+1)) then

            ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
            ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))

            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else
            i1 = icx / 2
          end if

          if ((icy == (2*cell_start(2)-1) .or. icy == (2*cell_end(2)+1))) then
            if (num_domain(2) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i2 = icy1
              else
                i2 = icy2
              end if
            else
              i2 = (icy1 + icy2) / 2
            end if
          else
            i2 = (icy1 + icy2) / 2
          end if

          if ((icz == (2*cell_start(3)-1) .or. icz == (2*cell_end(3)+1))) then
            if (num_domain(3) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i3 = icz1
              else
                i3 = icz2
              end if
            else
              i3 = (icz1 + icz2) / 2
            end if
          else
            i3 = (icz1 + icz2) / 2
          end if

          cell_pair(i,j) = cell_gxyz2l(i1+1,i2+1,i3+1)
          cell_pair(j,i) = cell_gxyz2l(i1+1,i2+1,i3+1)

        end if

      end do
    end do

    do i = 1, ncel_local

      icx1 = cell_l2gx1(i)
      icy1 = cell_l2gy1(i)
      icz1 = cell_l2gz1(i)

      do jj = 1, nboundary

        j = jj + ncel_local

        icx2 = cell_l2gx1(j)
        icy2 = cell_l2gy1(j)
        icz2 = cell_l2gz1(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

        if (icx == abs(icx1-icx2-cell(1))) then
          cell_move(1,j,i) = -1
        else if (icx == abs(icx1-icx2+cell(1))) then
          cell_move(1,j,i) = 1
        else
          cell_move(1,j,i) = 0
        end if

        if (icy == abs(icy1-icy2-cell(2))) then
          cell_move(2,j,i) = -1
        else if (icy == abs(icy1-icy2+cell(2))) then
          cell_move(2,j,i) = 1
        else
          cell_move(2,j,i) = 0
        end if

        if (icz == abs(icz1-icz2-cell(3))) then
          cell_move(3,j,i) = -1
        else if (icz == abs(icz1-icz2+cell(3))) then
          cell_move(3,j,i) = 1
        else
          cell_move(3,j,i) = 0
        end if
        cell_move(1:3,i,j) = -cell_move(1:3,j,i)

        if (cell_move(1,j,i) == 0 .and. &
            cell_move(2,j,i) == 0 .and. &
            cell_move(3,j,i) == 0) then
          virial_check(j,i) = 0
          virial_check(i,j) = 0
        end if

        end if

      end do
    end do

    do ii = 1, nboundary
      ij = ij + 1
      i = ii + ncel_local
      cell_pair(i,i) = i
      cell_move(1:3,i,i) = 0
    end do

    do ii = 1, nboundary-1

      i = ii + ncel_local
      icx1 = cell_l2gx(i)
      icy1 = cell_l2gy(i)
      icz1 = cell_l2gz(i)

      do jj = ii+1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx(j)
        icy2 = cell_l2gy(j)
        icz2 = cell_l2gz(j)
        cell_pair(j,i) = 0
        cell_pair(i,j) = 0

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        movex = 0
        movey = 0
        movez = 0
        icx = min(abs(icx1-icx2),abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2),abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2),abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (icx == abs(icx1-icx2-cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = -1
          movex = -cell(1)
        else if (icx == abs(icx1-icx2+cell(1)) .and. num_domain(1) == 1) then
          cell_move(1,j,i) = 1
          movex = cell(1)
        else
          cell_move(1,j,i) = 0
          movex = 0
        end if

        if (icy == abs(icy1-icy2-cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = -1
          movey = -cell(2)
        else if (icy == abs(icy1-icy2+cell(2)) .and. num_domain(2) == 1) then
          cell_move(2,j,i) = 1
          movey = cell(2)
        else
          cell_move(2,j,i) = 0
          movey = 0
        end if

        if (icz == abs(icz1-icz2-cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = -1
          movez = -cell(3)
        else if (icz == abs(icz1-icz2+cell(3)) .and. num_domain(3) == 1) then
          cell_move(3,j,i) = 1
          movez = cell(3)
        else
          cell_move(3,j,i) = 0
          movez = 0
        end if

        cell_move(1:3,i,j) = - cell_move(1:3,j,i)

        icx = icx1 - icx2 + movex
        icy = icy1 - icy2 + movey
        icz = icz1 - icz2 + movez

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

          ij = ij + 1
          icx = icx1 + icx2
          icy = icy1 + icy2 + movey
          icz = icz1 + icz2 + movez

          if (icx == (2*cell_start(1)-1) .or. icx == (2*cell_end(1)+1)) then
            ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
            ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
            if (ic1 < ic2) then
              i1 = icx1
            else
              i1 = icx2
            end if
          else
            i1 = icx / 2
          end if

          if ((icy == (2*cell_start(2)-1) .or. icy == (2*cell_end(2)+1))) then
            if (num_domain(2) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i2 = icy1
              else
                i2 = icy2
              end if
            else
              i2 = (icy1 + icy2) / 2
            end if
          else
            i2 = (icy1 + icy2) / 2
          end if

          if ((icz == (2*cell_start(3)-1) .or. icz == (2*cell_end(3)+1))) then
            if (num_domain(3) > 1) then
              ic1 = natom(cell_gxyz2l(icx1+1,icy1+1,icz1+1))
              ic2 = natom(cell_gxyz2l(icx2+1,icy2+1,icz2+1))
              if (ic1 < ic2) then
                i3 = icz1
              else
                i3 = icz2
              end if
            else
              i3 = (icz1 + icz2) / 2
            end if
          else
            i3 = (icz1 + icz2) / 2
          end if

          cell_pair(i,j) = cell_gxyz2l(i1+1,i2+1,i3+1)
          cell_pair(j,i) = cell_gxyz2l(i1+1,i2+1,i3+1)

        end if
      end do
    end do

    do ii = 1, nboundary-1

      i = ii + ncel_local
      icx1 = cell_l2gx1(i)
      icy1 = cell_l2gy1(i)
      icz1 = cell_l2gz1(i)

      do jj = ii+1, nboundary

        j = jj + ncel_local
        icx2 = cell_l2gx1(j)
        icy2 = cell_l2gy1(j)
        icz2 = cell_l2gz1(j)

        icx = icx1 - icx2
        icy = icy1 - icy2
        icz = icz1 - icz2

        icx = min(abs(icx1-icx2), abs(icx1-icx2-cell(1)),abs(icx1-icx2+cell(1)))
        icy = min(abs(icy1-icy2), abs(icy1-icy2-cell(2)),abs(icy1-icy2+cell(2)))
        icz = min(abs(icz1-icz2), abs(icz1-icz2-cell(3)),abs(icz1-icz2+cell(3)))

        if (abs(icx) <= 2 .and. abs(icy) <= 2 .and. abs(icz) <= 2) then

        if (icx == abs(icx1-icx2-cell(1))) then
          cell_move(1,j,i) = -1
        else if (icx == abs(icx1-icx2+cell(1))) then
          cell_move(1,j,i) = 1
        else
          cell_move(1,j,i) = 0
        end if

        if (icy == abs(icy1-icy2-cell(2))) then
          cell_move(2,j,i) = -1
        else if (icy == abs(icy1-icy2+cell(2))) then
          cell_move(2,j,i) = 1
        else
          cell_move(2,j,i) = 0
        end if

        if (icz == abs(icz1-icz2-cell(3))) then
          cell_move(3,j,i) = -1
        else if (icz == abs(icz1-icz2+cell(3))) then
          cell_move(3,j,i) = 1
        else
          cell_move(3,j,i) = 0
        end if

        cell_move(1,i,j) = -cell_move(1,j,i)
        cell_move(2,i,j) = -cell_move(2,j,i)
        cell_move(3,i,j) = -cell_move(3,j,i)

        if (abs(cell_move(1,j,i)) <= EPS .and. &
            abs(cell_move(2,j,i)) <= EPS .and. &
            abs(cell_move(3,j,i)) <= EPS) then
          virial_check(j,i) = 0
          virial_check(i,j) = 0
        end if

        end if

      end do
    end do

    univ_maxcell1 = ij + ncel_local

    return

  end subroutine assign_cell_interaction

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    molecule_to_domain
  !> @brief        copy molecule information to domain
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine molecule_to_domain(molecule, move, origin, iatom, domain, &
                                icel, icel_atom, dupl, box, ioffset)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    real(wip),                intent(in)    :: move(3)
    real(wip),                intent(in)    :: origin(3)
    integer,                  intent(in)    :: iatom
    type(s_domain),   target, intent(inout) :: domain
    integer,                  intent(in)    :: icel
    integer,                  intent(in)    :: icel_atom
    integer,                  intent(in)    :: dupl(3)
    real(wip),                intent(in)    :: box(3)
    integer,                  intent(in)    :: ioffset

    ! local variables
    real(wp),         pointer :: coord(:,:), vel(:,:)
    real(wp),         pointer :: charge(:), mass(:)
    integer,          pointer :: atom_class(:)

    real(wip),        pointer :: coord_local (:,:,:)
    real(wip),        pointer :: ref_local   (:,:,:)
    real(wip),        pointer :: vel_local   (:,:,:)
    real(wp),         pointer :: charge_local(:,:)
    real(wip),        pointer :: mass_local  (:,:), inv_mass(:,:)
    real(wp),         pointer :: trans       (:,:,:)
    integer,          pointer :: class_local (:,:)
    integer,          pointer :: id_l2g      (:,:)
    integer(int2),    pointer :: id_g2l      (:,:)

    ! FEP
    integer                   :: iatom2, k

    coord        => molecule%atom_coord
    vel          => molecule%atom_velocity
    charge       => molecule%charge
    mass         => molecule%mass
    atom_class   => molecule%atom_cls_no

    id_g2l       => domain%id_g2l
    id_l2g       => domain%id_l2g
    coord_local  => domain%coord
    ref_local    => domain%coord_ref
    vel_local    => domain%velocity
    charge_local => domain%charge
    mass_local   => domain%mass
    inv_mass     => domain%inv_mass
    class_local  => domain%atom_cls_no
    trans        => domain%trans_vec


    id_l2g(icel_atom,icel)    = iatom + ioffset

    coord_local(1:3,icel_atom,icel) = coord (1:3,iatom)              &
                                    + real(dupl(1:3)-1,wip)*box(1:3) &
                                    - origin(1:3)
    ref_local  (1:3,icel_atom,icel) = coord_local(1:3,icel_atom,icel)
    vel_local  (1:3,icel_atom,icel) = vel   (1:3,iatom)
    charge_local   (icel_atom,icel) = charge    (iatom)
    mass_local     (icel_atom,icel) = mass      (iatom)
    if (mass_local(icel_atom,icel) > EPS) &
      inv_mass(icel_atom,icel) = 1.0_wip / mass_local(icel_atom, icel)
    class_local    (icel_atom,icel) = atom_class(iatom)
    trans      (1:3,icel_atom,icel) = move(1:3)

    ! FEP
    if (domain%fep_use) then
      if (molecule%fepgrp(iatom)==1) then
        domain%fep_chargeA(icel_atom,icel) = charge(iatom)
        if (domain%num_atom_single_all > 0) then
          do k = 1, domain%num_atom_single_all
            if (domain%id_singleA(k) == iatom) then
              iatom2 = domain%id_singleB(k)
              exit
            end if
          end do
          domain%fep_chargeB(icel_atom,icel) = charge(iatom2)
          domain%fep_atmcls_singleB(icel_atom,icel) = atom_class(iatom2)
        else
          domain%fep_chargeB(icel_atom,icel) = charge(iatom)
          domain%fep_atmcls_singleB(icel_atom,icel) = atom_class(iatom)
        end if
      else if (molecule%fepgrp(iatom)==3) then
        domain%fep_chargeA(icel_atom,icel) = charge(iatom)
      else if (molecule%fepgrp(iatom)==4) then
        domain%fep_chargeB(icel_atom,icel) = charge(iatom)
      end if
      domain%fepgrp(icel_atom,icel) = int(molecule%fepgrp(iatom))
    end if

    return

  end subroutine molecule_to_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord
  !> @brief        check_atom_coordinate
  !! @authors      CK, HO
  !! @param[in]    ene_info      : ENERGY section control parameters information
  !! @param[in]    boundary      : boundary condition information
  !! @param[inout] domain        : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord(ene_info, boundary, domain)

    ! formal arguments
    type(s_ene_info), target, intent(in)    :: ene_info
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, ix, iy, ij, j
    integer                   :: id, omp_get_thread_num
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:)
    integer,          pointer :: ncell, nboundary
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: id_l2g(:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    real(wip),        pointer :: mass(:,:)
    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)

    ! FEP
    integer                   :: fg1, fg2

    mass           => domain%mass
    coord          => domain%coord
    cell_move      => domain%cell_move
    system_size    => domain%system_size
    natom          => domain%num_atom
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    cell_pairlist1 => domain%cell_pairlist1
    id_l2g         => domain%id_l2g

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        dir(1) = coord(1,ix,i)
        dir(2) = coord(2,ix,i)
        dir(3) = coord(3,ix,i)
        do iy = ix + 1, natom(i)

          if (domain%fep_use) then
            ! FEP: Skip perturbed region
            fg1 = domain%fepgrp(ix,i)
            fg2 = domain%fepgrp(iy,i)
            if((fg1 == 1) .and. (fg2 == 2)) cycle
            if((fg1 == 2) .and. (fg2 == 1)) cycle
            if((fg1 == 1) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 2)) cycle
            if((fg1 == 3) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 3)) cycle
            if((fg1 == 5) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 5)) cycle
            if((fg1 == 5) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 5)) cycle
            if((fg1 == 3) .and. (fg2 == 3)) cycle
            if((fg1 == 4) .and. (fg2 == 4)) cycle
            if((fg1 == 1) .and. (fg2 == 3)) cycle
            if((fg1 == 3) .and. (fg2 == 1)) cycle
            if((fg1 == 2) .and. (fg2 == 4)) cycle
            if((fg1 == 4) .and. (fg2 == 2)) cycle
          end if

          dr(1) = dir(1) - coord(1,iy,i)
          dr(2) = dir(2) - coord(2,iy,i)
          dr(3) = dir(3) - coord(3,iy,i)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          if (rr < EPS) &
            call error_msg('Check_Atom_Coord> coordinates are too close')
          if (rr < ene_info%minimum_contact .and. &
              mass(ix,i)*mass(iy,i) > EPS) then
            !$omp critical
            write(MsgOut, '(A,I10,I10,F10.5)') &
            'WARNING: too small distances: ', id_l2g(ix,i), id_l2g(iy,i), &
                                              sqrt(rr)
            !$omp end critical
            if (rr < ene_info%err_minimum_contact .and. &
                mass(ix,i)*mass(iy,i) > EPS      .and. &
                .not. ene_info%contact_check) &
              call error_msg('Check_Atom_Coord> too small contact: please use nonb_limiter')
          end if
        end do
      end do
    end do

    if (boundary%type == BoundaryTypePBC) then
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)

        do ix = 1, natom(i)

          dir(1:3) = coord(1:3,ix,i)

          do iy = 1, natom(j)

            if (domain%fep_use) then
              ! FEP: Skip perturbed region
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if((fg1 == 1) .and. (fg2 == 2)) cycle
              if((fg1 == 2) .and. (fg2 == 1)) cycle
              if((fg1 == 1) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 1)) cycle
              if((fg1 == 2) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 2)) cycle
              if((fg1 == 3) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 3)) cycle
              if((fg1 == 5) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 5)) cycle
              if((fg1 == 5) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 5)) cycle
              if((fg1 == 3) .and. (fg2 == 3)) cycle
              if((fg1 == 4) .and. (fg2 == 4)) cycle
              if((fg1 == 1) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 1)) cycle
              if((fg1 == 2) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 2)) cycle
            end if

            dr(1) = dir(1) - coord(1,iy,j)
            dr(2) = dir(2) - coord(2,iy,j)
            dr(3) = dir(3) - coord(3,iy,j)
            dr(1) = dr(1) - system_size(1)*anint(dr(1)/system_size(1))
            dr(2) = dr(2) - system_size(2)*anint(dr(2)/system_size(2))
            dr(3) = dr(3) - system_size(3)*anint(dr(3)/system_size(3))
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) &
              call error_msg('Check_Atom_Coord> coordinates are too close')
            if (rr < ene_info%minimum_contact) then
              !$omp critical
              write(MsgOut, '(A,I10,I10,F10.5)') &
              'WARNING: too small distances: ', id_l2g(ix,i), id_l2g(iy,j), &
                                                sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and. &
                  .not. ene_info%contact_check) &
                call error_msg('Check_Atom_Coord> too small contact: please use nonb_limiter')
            end if
          end do
        end do
      end do
    else
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)
        do ix = 1, natom(i)

          dir(1:3) = coord(1:3,ix,i)

          do iy = 1, natom(j)

            if (domain%fep_use) then
              ! FEP: Skip perturbed region
              fg1 = domain%fepgrp(ix,i)
              fg2 = domain%fepgrp(iy,j)
              if((fg1 == 1) .and. (fg2 == 2)) cycle
              if((fg1 == 2) .and. (fg2 == 1)) cycle
              if((fg1 == 1) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 1)) cycle
              if((fg1 == 2) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 2)) cycle
              if((fg1 == 3) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 3)) cycle
              if((fg1 == 5) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 5)) cycle
              if((fg1 == 5) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 5)) cycle
              if((fg1 == 3) .and. (fg2 == 3)) cycle
              if((fg1 == 4) .and. (fg2 == 4)) cycle
              if((fg1 == 1) .and. (fg2 == 3)) cycle
              if((fg1 == 3) .and. (fg2 == 1)) cycle
              if((fg1 == 2) .and. (fg2 == 4)) cycle
              if((fg1 == 4) .and. (fg2 == 2)) cycle
            end if

            dr(1) = dir(1) - coord(1,iy,j)
            dr(2) = dir(2) - coord(2,iy,j)
            dr(3) = dir(3) - coord(3,iy,j)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) &
              call error_msg('Check_Atom_Coord> coordinates are too close')
            if (rr < ene_info%minimum_contact) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), &
                                              sqrt(rr)
              !$omp end critical
              if (rr < ene_info%err_minimum_contact .and.  &
                 .not. ene_info%contact_check) &
                call error_msg('Check_Atom_Coord> too small contact: please use nonb_limiter')
            end if
          end do
        end do
      end do
    end if

    !$omp end parallel

    return

  end subroutine check_atom_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_atom_coord_pio
  !> @brief        check_atom_coordinate
  !! @authors      CK
  !! @param[in]    enefunc  : energy potential function information
  !! @param[in]    boundary : boundary condition information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_atom_coord_pio(enefunc, boundary, domain)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_domain),   target, intent(inout) :: domain

    ! local variable
    integer                   :: i, ix, iy, ij, j
    integer                   :: id, omp_get_thread_num
    real(wp)                  :: rr, dr(3), dir(3)
    real(wp)                  :: tr(3)

    integer,          pointer :: natom(:)
    integer,          pointer :: ncell, nboundary
    integer(int2),    pointer :: cell_pairlist1(:,:)
    integer,          pointer :: id_l2g(:,:)
    integer(1),       pointer :: cell_move(:,:,:)
    real(wip),        pointer :: mass(:,:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: system_size(:)


    mass           => domain%mass
    coord          => domain%coord
    trans1         => domain%trans_vec
    trans2         => domain%translated
    cell_move      => domain%cell_move
    system_size    => domain%system_size
    natom          => domain%num_atom
    ncell          => domain%num_cell_local
    nboundary      => domain%num_cell_boundary
    cell_pairlist1 => domain%cell_pairlist1
    id_l2g         => domain%id_l2g

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy,  ij,  dir, dr, rr, tr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
          trans2(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
          trans2(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
        end do
      end do
    else
      do i = id+1, ncell+nboundary, nthread
        do ix = 1, natom(i)
          trans2(ix,1,i) = coord(1,ix,i)
          trans2(ix,2,i) = coord(2,ix,i)
          trans2(ix,3,i) = coord(3,ix,i)
        end do
      end do
    end if
    !$omp barrier

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        dir(1:3) = trans2(ix,1:3,i)
        do iy = ix + 1, natom(i)
          dr(1) = dir(1) - trans2(iy,1,i)
          dr(2) = dir(2) - trans2(iy,2,i)
          dr(3) = dir(3) - trans2(iy,3,i)
          rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          if (rr < EPS) &
            call error_msg('Check_Atom_Coord_Pio> coordinates are too close')

          if (rr < enefunc%minimum_contact .and. &
              mass(ix,i)*mass(iy,i) > EPS) then
            !$omp critical
            write(MsgOut,'(A,I10,I10,F10.5)') &
            'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, i), &
                                            sqrt(rr)
            !$omp end critical
          end if
        end do
      end do
    end do

    if (boundary%type == BoundaryTypePBC) then
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)

        tr(1:3) = real(cell_move(1:3,j,i),wp)*system_size(1:3)
        do ix = 1, natom(i)

          dir(1:3) = trans2(ix,1:3,i)

          do iy = 1, natom(j)
            dr(1) = dir(1) - trans2(iy,1,j)+ tr(1)
            dr(2) = dir(2) - trans2(iy,2,j)+ tr(2)
            dr(3) = dir(3) - trans2(iy,3,j)+ tr(3)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) &
            call error_msg('Check_Atom_Coord_Pio> coordinates are too '//&
            'close near')

            if (rr < enefunc%minimum_contact) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), &
                                              sqrt(rr)
              !$omp end critical
            end if

          end do
        end do
      end do
    else
      do ij = id+1, maxcell_near, nthread
        i = cell_pairlist1(1,ij)
        j = cell_pairlist1(2,ij)
        do ix = 1, natom(i)

          dir(1:3) = trans2(ix,1:3,i)

          do iy = 1, natom(j)
            dr(1) = dir(1) - trans2(iy,1,j)
            dr(2) = dir(2) - trans2(iy,2,j)
            dr(3) = dir(3) - trans2(iy,3,j)
            rr = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            if (rr < EPS) &
              call error_msg('Check_Atom_Coord_Pio> coordinates are too close')

            if (rr < enefunc%minimum_contact) then
              !$omp critical
              write(MsgOut,'(A,I10,I10,F10.5)') &
              'WARNING: too small distances:',id_l2g(ix,i), id_l2g(iy, j), &
                                              sqrt(rr)
              !$omp end critical
            end if

          end do
        end do
      end do
    end if

    !$omp end parallel

    return

  end subroutine check_atom_coord_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    select_kernel
  !> @brief        kernel selection
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    boundary : boundary information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine select_kernel(ene_info, boundary, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i
    logical                  :: has_avx2, has_avx512
    character(20), allocatable :: cpuflags(:)


    select case(ene_info%nonbond_kernel)

    case (NBK_Autoselect)

#if defined(FUGAKU)

      domain%pairlist_kernel = PLK_Fugaku
      domain%nonbond_kernel  = NBK_Fugaku

#elif defined(KCOMP)

      domain%pairlist_kernel = PLK_Generic
      domain%nonbond_kernel  = NBK_Fugakusmall

#elif defined(USE_GPU)

      domain%pairlist_kernel = PLK_GPU
      domain%nonbond_kernel  = NBK_GPU

#else

      call get_cpu_flags(cpuflags)

      has_avx2 = .false.
      has_avx512 = .false.
      do i = 1, size(cpuflags)

        if (cpuflags(i)(1:4) .eq. 'avx2') then
          has_avx2 = .true.
        else if (cpuflags(i)(1:6) .eq. 'avx512') then
          has_avx512 = .true.
        end if

      end do

      if (has_avx512 .or. has_avx2) then

        domain%pairlist_kernel = PLK_Intel
        domain%nonbond_kernel  = NBK_Intel

      else

        domain%pairlist_kernel = PLK_Generic
        domain%nonbond_kernel  = NBK_Generic

      end if

      deallocate(cpuflags)

#endif

    case (NBK_Generic)

      domain%pairlist_kernel = PLK_Generic
      domain%nonbond_kernel  = NBK_Generic

    case (NBK_Intel)

      domain%pairlist_kernel = PLK_Intel
      domain%nonbond_kernel  = NBK_Intel

    case (NBK_Fugaku)

      domain%pairlist_kernel = PLK_Fugaku
      domain%nonbond_kernel  = NBK_Fugaku

    end select

    if (ene_info%nonb_limiter) then

      domain%pairlist_kernel = PLK_Generic
      domain%nonbond_kernel  = NBK_Generic

    end if

    if (ene_info%electrostatic == ElectrostaticCutoff) then
      domain%pairlist_kernel = PLK_Generic
      domain%nonbond_kernel  = NBK_Generic
    end if

    if ((boundary%num_domain(1) == 1 .or. &
         boundary%num_domain(2) == 1 .or. &
         boundary%num_domain(3) == 1) .and. &
         domain%nonbond_kernel /= NBK_GPU) then
      domain%pairlist_kernel = PLK_Generic
      domain%nonbond_kernel  = NBK_Generic
    end if

    if (ene_info%nonbond_kernel == NBK_GPU) then

      if (ene_info%electrostatic == ElectrostaticCUTOFF) &
        call error_msg( &
           'Read_Ctrl_Energy> Electrostatic cutoff is not available with GPU')

      if (ene_info%nonb_limiter) &
        call error_msg( &
           'Read_Ctrl_Energy> nonb_limiter is not available with GPU')

      if (ene_info%structure_check /= StructureCheckNone) &
        call error_msg( &
           'Read_Ctrl_Energy> structure_check is not available with GPU')

    end if

    domain%scale_pairlist_fugaku  = ene_info%scale_pairlist_fugaku
    domain%scale_pairlist_generic = ene_info%scale_pairlist_generic

    if (main_rank) then

      write(MsgOut,'(a)')   'Select_kernel> '
      write(MsgOut,'(a,a)') '  Pairlist        = ', PLK_Strings(domain%pairlist_kernel)
      write(MsgOut,'(a,a)') '  Nonbond         = ', NBK_Types(domain%nonbond_kernel)
      write(MsgOut,'(a)')   ''

    end if

    return

  end subroutine select_kernel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_domain_fep
  !> @brief        setup domain information for FEP
  !! @authors      HO
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    con_info : CONSTRAINTS section control parameters information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] enefunc     : energy potential function information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_domain_fep(ene_info, con_info, &
                          boundary, molecule, enefunc, constraints, domain)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_cons_info),       intent(in)    :: con_info
    type(s_boundary),        intent(inout) :: boundary
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: i, j, k, cell(3)
    integer                  :: icel_local, icel
    integer                  :: ncel_local, ncel_bound, ncel_all
    integer                  :: num_solute_all


    ! initialize structure informations
    !
    call init_domain(domain)
    call init_enefunc(enefunc)
    call init_constraints(constraints)

    enefunc%table%num_all     = molecule%num_atoms

    domain%num_atom_all       = molecule%num_atoms            &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)
    domain%num_duplicate      = boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)

    molecule%num_deg_freedom  = molecule%num_deg_freedom      &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)

    if (main_rank .and. domain%num_atom_all /= molecule%num_atoms) &
      write(MsgOut,*) 'Total number of atoms   : ', &
                      domain%num_atom_all

    domain%system_size(1)     = boundary%box_size_x
    domain%system_size(2)     = boundary%box_size_y
    domain%system_size(3)     = boundary%box_size_z
    domain%system_size_ini(1) = boundary%box_size_x
    domain%system_size_ini(2) = boundary%box_size_y
    domain%system_size_ini(3) = boundary%box_size_z

    domain%cell_size(1)       = real(boundary%cell_size_x,wp)
    domain%cell_size(2)       = real(boundary%cell_size_y,wp)
    domain%cell_size(3)       = real(boundary%cell_size_z,wp)

    enefunc%table%table       = ene_info%table
    enefunc%table%water_model = ene_info%water_model

    constraints%rigid_bond    = con_info%rigid_bond
    constraints%fast_water    = con_info%fast_water
    constraints%water_model   = con_info%water_model
    constraints%hydrogen_type = con_info%hydrogen_type

    ! FEP
    domain%num_atom_single_all  = molecule%num_atoms_fep(2)
    domain%fep_use     = .true.

    ! select kernel
    !
    call select_kernel(ene_info, boundary, domain)

    ! assign the rank of each dimension from my_rank
    !
    call setup_processor_rank(boundary, domain, cell)

    ! decide cell capacity (max**) for memory allocation
    !
    call setup_cell_capacity(boundary, domain, molecule)


    ! memory allocaltion of maps connecting local to global cell indices
    !
    ncel_local = domain%num_cell_local
    ncel_bound = domain%num_cell_boundary
    ncel_all   = ncel_local + ncel_bound

    call alloc_domain(domain, DomainCellGlobal, cell(1),cell(2),cell(3))
    call alloc_domain(domain, DomainCellLocal,    ncel_local, 1, 1)
    call alloc_domain(domain, DomainCellLocBou,   ncel_all,   1, 1)
    call alloc_domain(domain, DomainCellBoundary, ncel_bound, 1, 1)
    call alloc_domain(domain, DomainCellPair,     ncel_all,   1, 1)

    ! assign global<->local mapping of cell indexa
    !
    icel_local = 0
    do i = domain%cell_start(3), domain%cell_end(3)
      do j = domain%cell_start(2), domain%cell_end(2)
        do k = domain%cell_start(1), domain%cell_end(1)
          icel_local = icel_local + 1
          icel = k + (j-1)*cell(1) + (i-1)*cell(1)*cell(2)
          domain%cell_g2l(icel) = icel_local
          domain%cell_l2g(icel_local) = icel
          domain%cell_l2gx(icel_local) = k
          domain%cell_l2gy(icel_local) = j
          domain%cell_l2gz(icel_local) = i
          domain%cell_l2gx_orig(icel_local) = k
          domain%cell_l2gy_orig(icel_local) = j
          domain%cell_l2gz_orig(icel_local) = i
          domain%cell_gxyz2l(k,j,i) = icel_local
        end do
      end do
    end do

    ! assigin each boundary cell
    !
    call setup_cell_boundary(cell, boundary%num_domain, domain)

    ! assign of atom maps connecting global local to global atom indices
    !
    call alloc_domain(domain, DomainDynvar, ncel_all, 1, 1)

    ! decide hydrogen atom from mass
    !
    call check_light_atom_name(con_info%hydrogen_mass_upper_bound, molecule)

    call setup_solute_and_water(molecule, enefunc,       &
                                constraints%water_model, &
                                constraints%water_type,  &
                                constraints%num_water)

    num_solute_all            = enefunc%table%num_solute      &
                               *boundary%num_duplicate(1)     &
                               *boundary%num_duplicate(2)     &
                               *boundary%num_duplicate(3)
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)

    if (constraints%water_type == TIP4) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 4, 1)
    else if (constraints%water_type == TIP3) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 3, 1)
    else if (constraints%water_type == TIP2) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 2, 1)
    else if (constraints%water_type == TIP1) then
      call alloc_domain(domain, DomainDynvar_Atom, ncel_all, 1, 1)
    end if

    ! FEP
    call alloc_domain(domain, DomainFEP, ncel_all, 1, 1)

    call setup_hbond_group     (molecule, domain, enefunc, constraints)
    call setup_atom_by_HBond   (molecule, boundary, enefunc, constraints, &
                                domain)
    call setup_global_to_local_atom_index(enefunc, domain)

!   call setup_ring_check      (molecule, enefunc, constraints, domain)

#ifdef DEBUG
    ! debug
    !

    if (main_rank) then
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_domain_str     ::MaxAtom      : ', MaxAtom
      write(MsgOut,*) 'sp_domain_str     ::MaxWater     : ', MaxWater
      write(MsgOut,*) 'sp_domain_str     ::MaxMove      : ', MaxMove
      write(MsgOut,*) 'sp_domain_str     ::MaxWaterMove : ', MaxWaterMove
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15      : ', MaxNb15
      write(MsgOut,*) 'sp_enefunc_str    ::MaxNb15Water : ', MaxNb15Water
      write(MsgOut,*) ''
      write(MsgOut,*) 'sp_constraints_str::HGroupMax    : ', HGroupMax
      write(MsgOut,*) 'sp_constraints_str::HGrpMaxMove  : ', HGrpMaxMove
      write(MsgOut,*) ''

    end if
#endif

    ! assign the interaction cell for each interaction
    !
    call setup_domain_interaction(boundary, domain)

    call check_atom_coord(ene_info, boundary, domain)

    ! setup water molecule information
    !
    if (constraints%water_type >= TIP3) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_O
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(3) = enefunc%table%atom_cls_no_H
      domain%water%atom_cls_no(4) = enefunc%table%atom_cls_no_D
      domain%water%charge(1)      = enefunc%table%charge_O
      domain%water%charge(2)      = enefunc%table%charge_H
      domain%water%charge(3)      = enefunc%table%charge_H
      domain%water%charge(4)      = enefunc%table%charge_D
      domain%water%mass(1)        = enefunc%table%mass_O
      domain%water%mass(2)        = enefunc%table%mass_H
      domain%water%mass(3)        = enefunc%table%mass_H
      domain%water%mass(4)        = enefunc%table%mass_D
    else if (constraints%water_type == TIP2) then
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_1
      domain%water%atom_cls_no(2) = enefunc%table%atom_cls_no_2
      domain%water%charge(1)      = enefunc%table%charge_1
      domain%water%charge(2)      = enefunc%table%charge_2
      domain%water%mass(1)        = enefunc%table%mass_1
      domain%water%mass(2)        = enefunc%table%mass_2
    else
      domain%water%atom_cls_no(1) = enefunc%table%atom_cls_no_1
      domain%water%charge(1)      = enefunc%table%charge_1
      domain%water%mass(1)        = enefunc%table%mass_1
    end if

    return

  end subroutine setup_domain_fep



  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_fep_correspondence(molecule, domain)

    ! formal arguments
    type(s_molecule), intent(inout) :: molecule
    type(s_domain),   intent(inout) :: domain

    ! local variables
    integer                         :: iA, iB, k

    if (domain%num_atom_single_all > 0) then

      allocate(domain%id_singleA(domain%num_atom_single_all))
      allocate(domain%id_singleB(domain%num_atom_single_all))

      iA = 1
      iB = 1
      do k = 1, molecule%num_atoms
        if (molecule%fepgrp(k) == 1) then
          domain%id_singleA(iA) = k
          iA = iA + 1
        else if (molecule%fepgrp(k) == 2) then
          domain%id_singleB(iB) = k
          iB = iB + 1
        end if
      end do

      if ((molecule%num_atoms_fep(2) > 0) .and. (iA /= iB)) then
        call error_msg('The number of atoms in singleB is different from that in singleA.')
      end if

    end if

    return

  end subroutine setup_fep_correspondence

end module sp_domain_mod
