!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   cg_domain_str
!> @brief   Structure of domain information.
!! @authors Jaewoon Jung (JJ) 
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_domain_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  type, public :: s_domain

    integer                       :: num_atom_all
    integer                       :: num_atom_domain
    integer                       :: num_atom_boundary
    integer                       :: num_deg_freedom
    integer                       :: num_cell_local
    integer                       :: num_cell_boundary
    integer                       :: max_num_atom
    integer                       :: num_comm_proc
    integer                       :: start_index(0:48), end_index(0:48)

    integer                       :: cell_start(3)
    integer                       :: cell_end(3)
    integer                       :: num_upper(3)
    integer                       :: num_lower(3)
    integer                       :: num_cell(3)
    integer                       :: cell_length(3)
    integer,          allocatable :: iproc(:)
    integer,          allocatable :: cell_start_proc(:,:)
    integer,          allocatable :: cell_end_proc  (:,:)

    ! system size
    real(wp)                      :: system_size(3)
    real(wp)                      :: cell_size(3)

    !  DomainCellGlobal
    integer,          allocatable :: cell_g2l(:)
    integer,          allocatable :: cell_g2b(:)
    integer,          allocatable :: cell_gxyz2l(:,:,:)
    integer,          allocatable :: natom_global(:,:,:)
    !  DomainCellLocal
    integer,          allocatable :: cell_l2g(:)
    !  DomainCellLocalBoundary
    integer,          allocatable :: cell_l2gx(:)
    integer,          allocatable :: cell_l2gy(:)
    integer,          allocatable :: cell_l2gz(:)
    integer,          allocatable :: cell_l2gx_orig(:)
    integer,          allocatable :: cell_l2gy_orig(:)
    integer,          allocatable :: cell_l2gz_orig(:)
    !  DomainCellBoundary
    integer,          allocatable :: cell_b2g(:)
    !  DomainCellPair
    real(wp),         allocatable :: cell_pbc_move(:,:)
    !  DomainNeighbourCell
    integer,          allocatable :: near_neighbor_cell_count(:)
    integer,          allocatable :: neighbor_cell_count(:)
    integer,          allocatable :: neighbor_cells(:,:)
    integer,          allocatable :: near_neighbor_cells(:,:)
    integer,          allocatable :: near_cells_count(:)
    integer,          allocatable :: near_cells(:,:)
    integer,          allocatable :: far_cells_count(:)
    integer,          allocatable :: far_cells(:,:)
    !  DomainDynvar
    integer,          allocatable :: id_l2g(:)
    integer,          allocatable :: base_list(:)
    integer,          allocatable :: phos_list(:)
    integer,          allocatable :: num_atom(:)
    integer,          allocatable :: domain_cell_rank(:)
    real(wp),         allocatable :: num_atom_t0(:)
    integer,          allocatable :: start_atom(:)
    integer,          allocatable :: num_charge(:)
    integer,          allocatable :: num_nocharge(:)
    integer,          allocatable :: num_base(:)
    integer,          allocatable :: num_phos(:)
    integer,          allocatable :: atom_cls_no(:)
    integer,          allocatable :: atom_2_cell(:)
    integer,          allocatable :: mol_chain_id(:)
    integer(1),       allocatable :: charge_type(:)
    integer,          allocatable :: na_base_type(:)
    integer(1),       allocatable :: dna_check(:)
    integer(1),       allocatable :: pwmcos_protein(:)
    integer(1),       allocatable :: pwmcosns_protein(:)
    integer(1),       allocatable :: cg_pro_use_KH(:)
    integer(1),       allocatable :: cg_IDR_KH(:)
    integer(1),       allocatable :: cg_IDR_HPS(:)
    integer(1),       allocatable :: cg_IDR_MPIPI(:)
    real(wp),         allocatable :: trans_vec(:,:)
    real(wp),         allocatable :: translated(:,:)
    real(wp),         allocatable :: charge(:)
    real(wp),         allocatable :: force_omp(:,:,:)
    real(wp),         allocatable :: force_pbc(:,:,:)
    real(wip),        allocatable :: coord(:,:)
    real(wip),        allocatable :: coord_ref(:,:)
    real(wip),        allocatable :: coord_old(:,:)
    real(wip),        allocatable :: coord_prev(:,:)
    real(wip),        allocatable :: velocity(:,:)
    real(wip),        allocatable :: velocity_ref(:,:)
    real(wip),        allocatable :: velocity_full(:,:)
    real(wip),        allocatable :: velocity_half(:,:)
    real(wip),        allocatable :: force(:,:)
    real(wip),        allocatable :: force_short(:,:)
    real(wip),        allocatable :: force_long(:,:)
    real(wip),        allocatable :: mass(:)
    real(wip),        allocatable :: random(:)
    real(dp),         allocatable :: virial_cellpair(:,:)
    !  DomainGlobal
    integer,          allocatable :: id_g2l(:)
    !  DomainPtlMove
    integer                       :: n_stay1
    integer                       :: n_stay2
    integer                       :: charge_move_domain
    integer                       :: charge_comm_domain
    integer                       :: nocharge_move_domain
    integer                       :: nocharge_comm_domain
    integer,          allocatable :: num_proc(:)
    integer,          allocatable :: proc_list(:,:)
    real(wip),        allocatable :: buf_var0_move_real(:)
    real(wip),        allocatable :: buf_var0_comm_real(:)
    real(wip),        allocatable :: buf_var0_stay_real(:)
    integer,          allocatable :: buf_var0_move_int (:)
    integer,          allocatable :: buf_var0_comm_int (:)
    integer,          allocatable :: buf_var0_stay_int (:)
    integer,          allocatable :: type1_stay(:)
    integer,          allocatable :: type2_stay(:)
    integer,          allocatable :: type3_stay(:)
    integer,          allocatable :: type1_move(:)
    integer,          allocatable :: type2_move(:)
    integer,          allocatable :: type1_comm(:)
    integer,          allocatable :: type2_comm(:)
    integer,          allocatable :: type3_comm(:)
    integer,          allocatable :: type1_comm_move(:)
    integer,          allocatable :: type2_comm_move(:)

  end type s_domain

  ! parameters for allocatable variables
  integer,      public, parameter :: DomainCellGlobal       = 1
  integer,      public, parameter :: DomainCellLocal        = 2
  integer,      public, parameter :: DomainCellLocBou       = 3
  integer,      public, parameter :: DomainCellBoundary     = 4
  integer,      public, parameter :: DomainCellPair         = 5
  integer,      public, parameter :: DomainCellPairList     = 6
  integer,      public, parameter :: DomainCellPairListELEC = 7
  integer,      public, parameter :: DomainNeighbourCell    = 8
  integer,      public, parameter :: DomainDynvar           = 9
  integer,      public, parameter :: DomainDynvar_Atom      = 10
  integer,      public, parameter :: DomainGlobal           = 11
  integer,      public, parameter :: DomainPtlMove          = 12
  integer,      public, parameter :: DomainPtlArray         = 13
  integer,      public, parameter :: DomainUnivCellPairList = 14
  integer,      public, parameter :: DomainNeighborRank     = 15
  integer,      public, parameter :: DomainComm             = 16

  ! variables for maximum numbers in one cell
  integer,      public            :: MaxAtom_domain         = 0
  integer,      public            :: Max_cg_Base
  integer,      public            :: Max_cg_Dna 
  integer,      public            :: Max_cg_KH
  integer,      public            :: Max_cg_IDR_KH
  integer,      public            :: Max_cg_IDR_HPS
  integer,      public            :: Max_cg_IDR_MPIPI
  integer,      public            :: Max_cg_Elec
  integer,      public            :: Max_cg_Pwmcos
  integer,      public            :: Max_cg_Pwmcosns
  integer,      public            :: Max_alloc_size1
  integer,      public            :: Max_alloc_size2
  integer,      public            :: Max_alloc_size3
  integer,      public            :: Max_alloc_size4
  integer,      public            :: Max_alloc_size5
  integer,      public            :: Max_alloc_size6
  integer,      public            :: Max_alloc_size7
  integer,      public            :: Max_alloc_size8
  integer,      public            :: min_atom               = 160

  ! variables for maximum cells
  integer,      public            :: maxcell, maxcell_near

  ! variables : cpu time of each processor
  real(dp),     public            :: calc_time
  real(dp),     public            :: calc_time_prev

  ! subroutines
  public  :: init_domain
  public  :: alloc_domain
  public  :: dealloc_domain
  public  :: dealloc_domain_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_domain
  !> @brief        Initializes domain information.
  !! @authors      JJ
  !! @param[out]   domain  : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_domain(domain)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain


    domain%num_deg_freedom          = 0
    domain%num_cell_local           = 0
    domain%num_cell_boundary        = 0
    domain%max_num_atom             = 0
    domain%num_atom_all             = 0
    domain%cell_start(1:3)          = 0
    domain%cell_end(1:3)            = 0
    domain%cell_length(1:3)         = 0

    return

  end subroutine init_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_domain
  !> @brief        Allocates domain information.
  !! @authors      JJ    
  !! @param[inout] domain    : domain information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable
  !! @param[in]    var_size2 : 3rd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_domain(domain, variable, var_size, var_size1, var_size2)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (DomainCellGlobal)

      if (allocated(domain%cell_g2l)) then
        if (size(domain%cell_g2l(:)) /= var_size*var_size1*var_size2) &
          deallocate(domain%cell_g2l,           &
                     domain%cell_g2b,           &
                     domain%cell_gxyz2l,        &
                     domain%natom_global,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_g2l)) &
        allocate(domain%cell_g2l     (var_size*var_size1*var_size2),     &
                 domain%cell_g2b     (var_size*var_size1*var_size2),     &
                 domain%cell_gxyz2l  (0:var_size+1,0:var_size1+1,        &
                                                       0:var_size2+1),   &
                 domain%natom_global (var_size,var_size1,var_size2),     &
                 stat = alloc_stat)

      domain%cell_g2l     (1:var_size*var_size1*var_size2)   = 0
      domain%cell_g2b     (1:var_size*var_size1*var_size2)   = 0
      domain%cell_gxyz2l  (0:var_size+1, 0:var_size1+1, 0:var_size2+1) = 0
      domain%natom_global (1:var_size, 1:var_size1, 1:var_size2) = 0

    case (DomainCellLocal)

      if (allocated(domain%cell_l2g)) then
        if (size(domain%cell_l2g(:)) /= var_size) &
          deallocate(domain%cell_l2g, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_l2g)) &
        allocate(domain%cell_l2g(var_size), stat = alloc_stat)

      domain%cell_l2g(1:var_size) = 0

    case (DomainCellLocBou)

      if (allocated(domain%cell_l2gx)) then
        if (size(domain%cell_l2gx(:)) /= var_size) &
          deallocate(domain%cell_l2gx,           &
                     domain%cell_l2gy,           &
                     domain%cell_l2gz,           &
                     domain%cell_l2gx_orig,      &
                     domain%cell_l2gy_orig,      &
                     domain%cell_l2gz_orig,      &
                     domain%domain_cell_rank,    &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_l2gx))            &
        allocate(domain%cell_l2gx           (var_size), &
                 domain%cell_l2gy           (var_size), &
                 domain%cell_l2gz           (var_size), &
                 domain%cell_l2gx_orig      (var_size), &
                 domain%cell_l2gy_orig      (var_size), &
                 domain%cell_l2gz_orig      (var_size), &
                 domain%domain_cell_rank    (var_size), &
                 stat = alloc_stat)

      domain%cell_l2gx           (1:var_size) = 0
      domain%cell_l2gy           (1:var_size) = 0
      domain%cell_l2gz           (1:var_size) = 0
      domain%cell_l2gx_orig      (1:var_size) = 0
      domain%cell_l2gy_orig      (1:var_size) = 0
      domain%cell_l2gz_orig      (1:var_size) = 0
      domain%domain_cell_rank    (1:var_size) = 0

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        if (size(domain%cell_b2g(:)) /= var_size) &
          deallocate(domain%cell_b2g,             &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_b2g))       &
        allocate(domain%cell_b2g(var_size),       &
                 stat = alloc_stat)

      domain%cell_b2g     (1:var_size) = 0

    case (DomainCellPair)

      if (allocated(domain%cell_pbc_move)) then
        if (size(domain%cell_pbc_move(1,:)) /= var_size) then
          deallocate(domain%cell_pbc_move,     &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%cell_pbc_move)) then
        allocate(domain%cell_pbc_move (3, var_size),  &
                 stat = alloc_stat)
      end if

      domain%cell_pbc_move (1:3, 1:var_size) = 0.0_wp

    case (DomainNeighborRank)

      if (allocated(domain%iproc)) then
        if (size(domain%iproc(:)) /= var_size)   &
          deallocate(domain%iproc,               &
                     domain%cell_start_proc,     &
                     domain%cell_end_proc,       &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%iproc))               &
        allocate(domain%iproc          (var_size   ),  &
                 domain%cell_start_proc(var_size, 3),  &
                 domain%cell_end_proc  (var_size, 3),  &
                 stat = alloc_stat)

      domain%iproc          (1:var_size     ) = 0
      domain%cell_start_proc(1:var_size, 1:3) = 0
      domain%cell_end_proc  (1:var_size, 1:3) = 0

    case (DomainNeighbourCell)

      if (allocated(domain%neighbor_cell_count)) then
        if (size(domain%neighbor_cell_count(:)) /= var_size) &
          deallocate(domain%neighbor_cell_count,      &
                     domain%near_neighbor_cell_count, &
                     domain%neighbor_cells,           &
                     domain%near_neighbor_cells,      &
                     domain%near_cells_count,         &
                     domain%near_cells,               &
                     domain%far_cells_count,          &
                     domain%far_cells,                &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%neighbor_cell_count))      &
        allocate(domain%neighbor_cell_count(var_size),      &
                 domain%near_neighbor_cell_count(var_size), &
                 domain%neighbor_cells(125, var_size),      &
                 domain%near_neighbor_cells(27, var_size),  &
                 domain%near_cells_count(var_size),         &
                 domain%near_cells(27, var_size),           &
                 domain%far_cells_count(var_size),          &
                 domain%far_cells(98, var_size),            &
                 stat = alloc_stat)

      domain%neighbor_cell_count  (1:var_size)      = 0
      domain%near_neighbor_cell_count  (1:var_size) = 0
      domain%neighbor_cells(1:125, 1:var_size)      = 0
      domain%near_neighbor_cells(1:27, 1:var_size)  = 0
      domain%near_cells_count(1:var_size)           = 0
      domain%near_cells(1:27, 1:var_size)           = 0
      domain%far_cells_count(1:var_size)            = 0
      domain%far_cells(1:98, 1:var_size)            = 0

    case (DomainDynvar)

      if (allocated(domain%num_atom)) then
        if (size(domain%num_atom) /= var_size) then
          deallocate(domain%num_atom,               &
                     domain%num_atom_t0,            &
                     domain%start_atom,             &
                     domain%num_charge,             &
                     domain%num_nocharge,           &
                     domain%random,                 &
                     domain%num_base,               &
                     domain%num_phos,               &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%num_atom)) then
          allocate(domain%num_atom(var_size),         &
                   domain%num_atom_t0(var_size),      &
                   domain%start_atom(var_size),       &
                   domain%num_charge(var_size),       &
                   domain%num_nocharge(var_size),     &
                   domain%random(var_size),           &
                   domain%num_base(var_size),         &
                   domain%num_phos(var_size),         &
                   stat = alloc_stat)
      end if
      domain%num_atom        (1:var_size) = 0
      domain%num_atom_t0     (1:var_size) = 0
      domain%start_atom      (1:var_size) = 0
      domain%num_charge      (1:var_size) = 0
      domain%num_nocharge    (1:var_size) = 0
      domain%random          (1:var_size) = 0.0_wp
      domain%num_base        (1:var_size) = 0
      domain%num_phos        (1:var_size) = 0

    case (DomainDynvar_Atom)

      if (allocated(domain%id_l2g)) then
        if (size(domain%id_l2g(:)) /= var_size) then
          deallocate(domain%id_l2g,             &
                     domain%base_list,          &
                     domain%phos_list,          &
                     domain%atom_cls_no,        &
                     domain%atom_2_cell,        &
                     domain%mol_chain_id,       &
                     domain%NA_base_type,       &
                     domain%charge_type,        &
                     domain%dna_check,          &
                     domain%cg_pro_use_KH,      &
                     domain%cg_IDR_KH,          &
                     domain%cg_IDR_HPS,         &
                     domain%cg_IDR_MPIPI,       &
                     domain%pwmcos_protein,     &
                     domain%pwmcosns_protein,   &
                     domain%coord,              &
                     domain%coord_ref,          &
                     domain%coord_old,          &
                     domain%coord_prev,         &
                     domain%velocity,           &
                     domain%velocity_ref,       &
                     domain%velocity_full,      &
                     domain%velocity_half,      &
                     domain%trans_vec,          &
                     domain%translated,         &
                     domain%force,              &
                     domain%force_short,        &
                     domain%force_long,         &
                     domain%force_omp,          &
                     domain%force_pbc,          &
                     domain%charge,             &
                     domain%mass,               &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(domain%id_l2g)) then
        allocate(domain%id_l2g        (var_size),             &
                 domain%base_list     (var_size),             &
                 domain%phos_list     (var_size),             &
                 domain%atom_cls_no   (var_size),             &
                 domain%atom_2_cell   (var_size),             &
                 domain%mol_chain_id  (var_size),             &
                 domain%NA_base_type  (var_size),             &
                 domain%charge_type   (var_size),             & 
                 domain%dna_check     (var_size),             & 
                 domain%cg_pro_use_KH (var_size),             & 
                 domain%cg_IDR_KH     (var_size),             & 
                 domain%cg_IDR_HPS    (var_size),             & 
                 domain%cg_IDR_MPIPI  (var_size),             & 
                 domain%pwmcos_protein(var_size),             & 
                 domain%pwmcosns_protein(var_size),           & 
                 domain%coord         (var_size,3),           &
                 domain%coord_ref     (var_size,3),           &
                 domain%coord_old     (var_size,3),           &
                 domain%coord_prev    (var_size,3),           &
                 domain%velocity      (var_size,3),           &
                 domain%velocity_ref  (var_size,3),           &
                 domain%velocity_full (var_size,3),           &
                 domain%velocity_half (var_size,3),           &
                 domain%trans_vec     (var_size,3),           &
                 domain%translated    (var_size,3),           &
                 domain%force         (var_size,3),           &
                 domain%force_short   (var_size,3),           &
                 domain%force_long    (var_size,3),           &
                 domain%force_omp     (var_size,3, nthread),  &
                 domain%force_pbc     (var_size,3, nthread),  &
                 domain%charge        (var_size),             &
                 domain%mass          (var_size),             &
                 stat = alloc_stat)
      end if

      domain%id_l2g        (1:var_size)               = 0
      domain%base_list     (1:var_size)               = 0
      domain%phos_list     (1:var_size)               = 0
      domain%atom_cls_no   (1:var_size)               = 0
      domain%atom_2_cell   (1:var_size)               = 0
      domain%mol_chain_id  (1:var_size)               = 0
      domain%NA_base_type  (1:var_size)               = 0
      domain%charge_type   (1:var_size)               = 0
      domain%dna_check     (1:var_size)               = 0
      domain%cg_pro_use_KH (1:var_size)               = 0
      domain%cg_IDR_KH     (1:var_size)               = 0
      domain%cg_IDR_HPS    (1:var_size)               = 0
      domain%cg_IDR_MPIPI  (1:var_size)               = 0
      domain%pwmcos_protein(1:var_size)               = 0
      domain%pwmcosns_protein(1:var_size)             = 0
      domain%coord         (1:var_size,1:3)           = 0.0_wip
      domain%coord_ref     (1:var_size,1:3)           = 0.0_wip
      domain%coord_old     (1:var_size,1:3)           = 0.0_wip
      domain%coord_prev    (1:var_size,1:3)           = 0.0_wip
      domain%velocity      (1:var_size,1:3)           = 0.0_wip
      domain%velocity_ref  (1:var_size,1:3)           = 0.0_wip
      domain%velocity_full (1:var_size,1:3)           = 0.0_wip
      domain%trans_vec     (1:var_size,1:3)           = 0.0_wp
      domain%translated    (1:var_size,1:3)           = 0.0_wp
      domain%force         (1:var_size,1:3)           = 0.0_wip
      domain%force_short   (1:var_size,1:3)           = 0.0_wp
      domain%force_long    (1:var_size,1:3)           = 0.0_wp
      domain%force_omp     (1:var_size,1:3,1:nthread) = 0.0_wp
      domain%force_pbc     (1:var_size,1:3,1:nthread) = 0.0_wp
      domain%charge        (1:var_size)               = 0.0_wp
      domain%mass          (1:var_size)               = 0.0_wp

    case (DomainGlobal)

      if (allocated(domain%id_g2l)) then
        if (size(domain%id_g2l(:)) /= var_size) &
          deallocate(domain%id_g2l,               &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_g2l))         &
        allocate(domain%id_g2l(var_size),       &
                 stat = alloc_stat)

      domain%id_g2l     (1:var_size) = 0

    case (DomainPtlMove)

      if (allocated(domain%type1_stay)) then
        deallocate(domain%type1_stay,         &
                   domain%type2_stay,         &
                   domain%type1_move,         &
                   domain%type2_move,         &
                   domain%type1_comm_move,    &
                   domain%type2_comm_move,    &
                   domain%num_proc,           &
                   domain%proc_list,          &
                   stat = dealloc_stat)
      end if

      if (.not. allocated(domain%type1_stay)) &
        allocate(domain%type1_move        (           var_size), &
                 domain%type2_move        (           var_size), &
                 domain%type1_stay        (           var_size), &
                 domain%type2_stay        (           var_size), &
                 domain%type1_comm_move   (           var_size), &
                 domain%type2_comm_move   (           var_size), &
                 domain%num_proc          (           var_size), &
                 domain%proc_list         (var_size1, var_size), &
                 stat = alloc_stat)

      domain%type1_stay        (            1:var_size) = 0
      domain%type2_stay        (            1:var_size) = 0
      domain%type1_move        (            1:var_size) = 0
      domain%type2_move        (            1:var_size) = 0
      domain%type1_comm_move   (            1:var_size) = 0
      domain%type2_comm_move   (            1:var_size) = 0
      domain%num_proc          (            1:var_size) = 0
      domain%proc_list         (1:var_size1,1:var_size) = 0

    case (DomainPtlArray)

      if (allocated(domain%buf_var0_stay_real)) then
        if (size(domain%buf_var0_stay_real(:)) /= var_size) &
          deallocate(domain%buf_var0_stay_real, &
                     domain%buf_var0_stay_int , &
                     domain%buf_var0_move_real, &
                     domain%buf_var0_move_int , &
                     domain%buf_var0_comm_real, &
                     domain%buf_var0_comm_int , &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%buf_var0_stay_real))  &
        allocate(domain%buf_var0_stay_real(var_size ), &
                 domain%buf_var0_stay_int (var_size1), &
                 domain%buf_var0_move_real(var_size ), &
                 domain%buf_var0_move_int (var_size1), &
                 domain%buf_var0_comm_real(var_size ), &
                 domain%buf_var0_comm_int (var_size1), &
                 stat = alloc_stat)

      domain%buf_var0_stay_real(1:var_size ) = 0.0_wp
      domain%buf_var0_stay_int (1:var_size1) = 0
      domain%buf_var0_move_real(1:var_size ) = 0.0_wp
      domain%buf_var0_move_int (1:var_size1) = 0
      domain%buf_var0_comm_real(1:var_size ) = 0.0_wp
      domain%buf_var0_comm_int (1:var_size1) = 0

    case (DomainComm)

      if (allocated(domain%type1_comm)) then
        if (size(domain%type1_comm(:)) /= var_size) &
          deallocate(domain%type1_comm,         &
                     domain%type2_comm,         &
                     domain%type3_comm,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%type1_comm))  &
        allocate(domain%type1_comm (var_size), &
                 domain%type2_comm (var_size), &
                 domain%type3_comm (var_size), &
                 stat = alloc_stat)

      domain%type1_comm(1:var_size) = 0
      domain%type2_comm(1:var_size) = 0
      domain%type3_comm(1:var_size) = 0

    case default

      call error_msg('Error: Alloc_Domain - Invalid variable selected.')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain
  !> @brief        Deallocates domain information.
  !! @authors      JJ    
  !! @@aram[inout] domain   : domain information
  !! @param[in]    variable : selected variable  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_domain(domain, variable)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    select case (variable)


    case (DomainCellGlobal)

      if (allocated(domain%cell_g2l)) then
        deallocate(domain%cell_g2l,     &
                   domain%cell_g2b,     &
                   domain%cell_gxyz2l,  &
                   domain%natom_global, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocal)

      if (allocated(domain%cell_l2g)) then
        deallocate(domain%cell_l2g, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocBou)

      if (allocated(domain%cell_l2gx)) then
        deallocate(domain%cell_l2gx,        &
                   domain%cell_l2gy,        &
                   domain%cell_l2gz,        &
                   domain%cell_l2gx_orig,   &
                   domain%cell_l2gy_orig,   &
                   domain%cell_l2gz_orig,   &
                   domain%domain_cell_rank, &
                   stat = dealloc_stat)
      end if

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        deallocate(domain%cell_b2g, &
                   stat = dealloc_stat)
      end if

    case (DomainCellPair)

      if (allocated(domain%cell_pbc_move)) then
        deallocate(domain%cell_pbc_move,&
                   stat = dealloc_stat)
      end if

    case (DomainNeighborRank)

      if (allocated(domain%iproc)) then
        deallocate(domain%iproc,               &
                   domain%cell_start_proc,     &
                   domain%cell_end_proc,       &
                   stat = dealloc_stat)
      end if

    case (DomainNeighbourCell)

      if (allocated(domain%neighbor_cell_count)) then
        deallocate(domain%neighbor_cell_count,      &
                   domain%near_neighbor_cell_count, &
                   domain%neighbor_cells,           &
                   domain%near_neighbor_cells,      &
                   domain%near_cells_count,         &
                   domain%near_cells,               &
                   domain%far_cells_count,          &
                   domain%far_cells,                &
                   stat = dealloc_stat)

      end if

    case (DomainDynvar)

      if (allocated(domain%num_atom)) then
        deallocate(domain%num_atom,       &
                   domain%num_atom_t0,    &
                   domain%start_atom,     &
                   domain%num_charge,     &
                   domain%num_nocharge,   &
                   domain%random,         &
                   domain%num_base,       &
                   domain%num_phos,       &
                   stat = dealloc_stat)
      end if

    case (DomainDynvar_Atom)

      if (allocated(domain%id_l2g)) then
        deallocate(domain%id_l2g,             &
                   domain%base_list,          &
                   domain%phos_list,          &
                   domain%atom_cls_no,        &
                   domain%atom_2_cell,        &
                   domain%mol_chain_id,       &
                   domain%NA_base_type,       &
                   domain%charge_type,        &
                   domain%dna_check,          &
                   domain%cg_pro_use_KH,      &
                   domain%cg_IDR_KH,          &
                   domain%cg_IDR_HPS,         &
                   domain%cg_IDR_MPIPI,       &
                   domain%pwmcos_protein,     &
                   domain%pwmcosns_protein,   &
                   domain%coord,              &
                   domain%coord_ref,          &
                   domain%coord_old,          &
                   domain%coord_prev,         &
                   domain%velocity,           &
                   domain%velocity_ref,       &
                   domain%velocity_full,      &
                   domain%velocity_half,      &
                   domain%trans_vec,          &
                   domain%translated,         &
                   domain%force,              &
                   domain%force_short,        &
                   domain%force_long,         &
                   domain%force_omp,          &
                   domain%force_pbc,          &
                   domain%charge,             &
                   domain%mass,               &
                   stat = dealloc_stat)
      end if

    case (DomainGlobal)

      if (allocated(domain%id_g2l)) then
        deallocate(domain%id_g2l,        &
                   stat = dealloc_stat)
      end if

    case (DomainPtlMove) 

      if (allocated(domain%type1_stay)) then
        deallocate(domain%type1_stay,         &
                   domain%type2_stay,         &
                   domain%type1_move,         &
                   domain%type2_move,         &
                   domain%type1_comm_move,    &
                   domain%type2_comm_move,    &
                   domain%num_proc,           &
                   domain%proc_list,          &
                   stat = dealloc_stat)
      end if

    case (DomainPtlArray)

      if (allocated(domain%buf_var0_stay_real)) then
        deallocate(domain%buf_var0_stay_real, &
                   domain%buf_var0_stay_int , &
                   domain%buf_var0_move_real, &
                   domain%buf_var0_move_int , &
                   domain%buf_var0_comm_real, &
                   domain%buf_var0_comm_int , &
                   stat = dealloc_stat)
      end if

    case (DomainComm)

      if (allocated(domain%type1_comm)) then
        deallocate(domain%type1_comm,         &
                   domain%type2_comm,         &
                   domain%type3_comm,         &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Error: Dealloc_Domain - Invalid variable selected.')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_domain


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain_all
  !> @brief        Deallocates all domain information.
  !! @authors      JJ
  !! @param[inout] domain : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_domain_all(domain)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain


    call dealloc_domain(domain, DomainCellGlobal)
    call dealloc_domain(domain, DomainCellLocal)
    call dealloc_domain(domain, DomainCellLocBou)
    call dealloc_domain(domain, DomainCellBoundary)
    call dealloc_domain(domain, DomainCellPair)
    call dealloc_domain(domain, DomainNeighborRank)
    call dealloc_domain(domain, DomainNeighbourCell)
    call dealloc_domain(domain, DomainDynvar)
    call dealloc_domain(domain, DomainDynvar_Atom)
    call dealloc_domain(domain, DomainGlobal)
    call dealloc_domain(domain, DomainPtlMove)
    call dealloc_domain(domain, DomainPtlArray)
    call dealloc_domain(domain, DomainComm)

    return

  end subroutine dealloc_domain_all

end module cg_domain_str_mod
