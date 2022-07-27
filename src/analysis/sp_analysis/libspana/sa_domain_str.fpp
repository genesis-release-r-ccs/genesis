!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_domain_str
!> @brief   structure of domain
!! @authors Jaewoon Jung (JJ), Isseki Yu(IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_domain_str_mod

  use sa_boundary_str_mod
  use sa_option_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_domain

    integer                       :: num_atom_all
    integer                       :: num_atom_selec
    integer                       :: num_group
    integer                       :: num_cell_local
    integer                       :: num_cell_boundary
    integer                       :: max_boundary_layer

    integer                       :: cell_start(3)
    integer                       :: cell_end(3)
    integer                       :: cell_length(3)
    integer                       :: iproc_lower(3)
    integer                       :: iproc_upper(3)
    integer                       :: neighbor(-1:1,-1:1,-1:1)

    !  DomainCellGlobal
    !  cell_gxyz2l(:,:,:) gives index of the cell in local domain from
    !  the global position of the cell (IY)
    integer,          allocatable :: cell_g2l(:)
    integer,          allocatable :: cell_g2b(:)
    integer,          allocatable :: cell_gxyz2l(:,:,:)
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
    !  DomainDynvar
    integer,          allocatable :: id_l2s(:,:)
    integer,          allocatable :: num_atom(:)
    integer,          allocatable :: num_atom_group(:,:)
    integer,          allocatable :: num_atom_t0(:)
    integer,          allocatable :: num_solute(:)
    integer,          allocatable :: num_water(:)
    real(wp),         allocatable :: coord(:,:,:)
    !  DomainGlobal
    integer,          allocatable :: id_global2selec(:)
    !  DomainSelec
    integer,          allocatable :: id_s2l(:,:)
    integer,          allocatable :: selec_atom2cell(:)
    integer,          allocatable :: id_selec2global(:)

  end type s_domain

  ! parameters for allocatable variables
  integer,      public, parameter :: DomainCellGlobal       = 1
  integer,      public, parameter :: DomainCellLocal        = 2
  integer,      public, parameter :: DomainCellLocBou       = 3
  integer,      public, parameter :: DomainCellBoundary     = 4
  integer,      public, parameter :: DomainDynvar           = 5
  integer,      public, parameter :: DomainDynvar_Atom      = 6
  integer,      public, parameter :: DomainGlobal           = 7
  integer,      public, parameter :: DomainSelec            = 8

  ! variables for maximum numbers in one cell
  integer,      public            :: MaxAtom                = 0
  integer,      public            :: MaxSolute              = 0
  integer,      public            :: MaxWater               = 0
  integer,      public            :: MaxMove                = 30
  integer,      public            :: MaxWaterMove           = 20

  ! subroutines
  public  :: init_domain
  public  :: alloc_domain
  public  :: dealloc_domain
  public  :: dealloc_domain_all


  !---(IY) these constants are imported from sp_energy_str
  ! variables for maximum numbers in one cell (these number will be updated)
  integer,      public            :: MaxBond   = 0
  integer,      public            :: MaxAngle  = 0
  integer,      public            :: MaxDihe   = 0
  integer,      public            :: MaxDihe_rb= 0
  integer,      public            :: MaxImpr   = 0
  integer,      public            :: MaxCmap   = 0
  integer,      public            :: MaxRest   = 100
  integer,      public            :: MaxExcl   = 36
  integer,      public            :: MaxContact= 0

  integer,      public            :: BondMove  = 0
  integer,      public            :: AngleMove = 0
  integer,      public            :: DiheMove  = 0
  integer,      public            :: ImprMove  = 0
  integer,      public            :: CmapMove  = 10
  integer,      public            :: RestMove  = 50
  integer,      public            :: ContactMove = 0

  integer,      public, parameter :: MaxAtomCls = 1000
  real(wp),     public            :: lj_coef(2,MaxAtomCls)

  !(IY) these constants are imported from sp_pairlist_str
  ! variables for maximum numbers in one cell
  integer,      public            :: MaxNb15      = 15000
  integer,      public            :: MaxNb15_NOBC = 3000
  integer,      public            :: MaxNb15Water = 1500
  !(IY) these constants are imported from sp_constraint_str
  integer,      public            :: HGroupMax   = 0
  integer,      public            :: HGrpMaxMove = 20
  !------

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_domain
  !> @brief        initialize domain information
  !! @authors      JJ
  !! @param[out]   domain  : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_domain(boundary, domain, option)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_boundary),        intent(in) :: boundary
    type(s_option),          intent(in) :: option

    real(wp)                 :: csize_xyz(3)
    real(wp)                 :: min_csize


    domain%num_cell_local           = 0
    domain%num_cell_boundary        = 0
    domain%num_atom_all             = 0
    domain%cell_start(3)            = 0
    domain%cell_end(3)              = 0
    domain%cell_length(3)           = 0
    domain%iproc_lower(3)           = 0
    domain%iproc_upper(3)           = 0
    domain%neighbor(-1:1,-1:1,-1:1) = 0

    !determin the max boundary layer
    csize_xyz(1)   =  boundary%cell_size_x
    csize_xyz(2)   =  boundary%cell_size_y
    csize_xyz(3)   =  boundary%cell_size_z
    min_csize      =  minval(csize_xyz)
    domain%max_boundary_layer = int(option%buffer/min_csize)+1

    return

  end subroutine init_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_domain
  !> @brief        allocate domain information
  !! @authors      JJ , IY
  !! @param[inout] domain    : domain information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable
  !! @param[in]    var_size1 : 2nd size of the selected variable
  !! @param[in]    var_size2 : 3rd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_domain(domain, option, variable, var_size,var_size1,var_size2)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2
    type(s_option),          intent(in)    :: option

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: mbl, bound_min, bound_max


    mbl          = domain%max_boundary_layer
    bound_min    = 1-mbl
    bound_max    = mbl

    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (DomainCellGlobal)

      if (allocated(domain%cell_g2l)) then
        if (size(domain%cell_g2l(:)) /= var_size*var_size1*var_size2) &
          deallocate(domain%cell_g2l,    &
                     domain%cell_g2b,    &
                     domain%cell_gxyz2l, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_g2l)) &
        allocate(domain%cell_g2l   (var_size*var_size1*var_size2), &
                 domain%cell_g2b   (var_size*var_size1*var_size2), &
                 ! size depends on boundary layer thickness (IY)
                 domain%cell_gxyz2l(bound_min:var_size +bound_max, &
                                    bound_min:var_size1+bound_max, &
                                    bound_min:var_size2+bound_max),&
                                    stat = alloc_stat)

      domain%cell_g2l   (1:var_size*var_size1*var_size2) = 0
      domain%cell_g2b   (1:var_size*var_size1*var_size2) = 0
      domain%cell_gxyz2l(bound_min:var_size  + bound_max,&
                         bound_min:var_size1 + bound_max,&
                         bound_min:var_size2 + bound_max) = 0

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
          deallocate(domain%cell_l2gx,      &
                     domain%cell_l2gy,      &
                     domain%cell_l2gz,      &
                     domain%cell_l2gx_orig, &
                     domain%cell_l2gy_orig, &
                     domain%cell_l2gz_orig, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_l2gx))      &
        allocate(domain%cell_l2gx     (var_size), &
                 domain%cell_l2gy     (var_size), &
                 domain%cell_l2gz     (var_size), &
                 domain%cell_l2gx_orig(var_size), &
                 domain%cell_l2gy_orig(var_size), &
                 domain%cell_l2gz_orig(var_size), &
                 stat = alloc_stat)

      domain%cell_l2gx     (1:var_size) = 0
      domain%cell_l2gy     (1:var_size) = 0
      domain%cell_l2gz     (1:var_size) = 0
      domain%cell_l2gx_orig(1:var_size) = 0
      domain%cell_l2gy_orig(1:var_size) = 0
      domain%cell_l2gz_orig(1:var_size) = 0

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        if (size(domain%cell_b2g(:)) /= var_size) &
          deallocate(domain%cell_b2g, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%cell_b2g)) &
        allocate(domain%cell_b2g(var_size), stat = alloc_stat)

      domain%cell_b2g(1:var_size) = 0

    case (DomainDynvar)

      if (allocated(domain%coord)) then
        if (size(domain%num_atom) /= var_size) &
          deallocate(domain%id_l2s,        &
                     domain%num_atom,      &
                     domain%num_atom_t0,   &
                     domain%num_atom_group, &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%coord)) &
        allocate(domain%num_atom     (var_size),              &
                 domain%num_atom_t0  (var_size),              &
                 domain%num_atom_group (var_size1, var_size), &
                 stat = alloc_stat)

      domain%num_atom     (1:var_size)                = 0
      domain%num_atom_t0  (1:var_size)                = 0
      domain%num_atom_group (1:var_size1, 1:var_size) = 0

    case (DomainDynvar_Atom)

      if (allocated(domain%id_l2s)) then
        if (size(domain%id_l2s(1,:)) /= var_size .or. &
            size(domain%id_l2s(:,1)) /= MaxAtom) &
          deallocate(domain%id_l2s,        &
                     domain%coord,         &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_l2s)) then
        allocate(domain%id_l2s(MaxAtom, var_size),     &
                 domain%coord (3, MaxAtom, var_size),  &
                 stat = alloc_stat)
      end if

      domain%id_l2s(1:MaxAtom, 1:var_size)        = 0
      domain%coord (1:3, 1:MaxAtom, 1:var_size)   = 0.0_wp

    case (DomainGlobal)

      if (allocated(domain%id_global2selec)) then
        if (size(domain%id_global2selec(:)) /= var_size) &
          deallocate(domain%id_global2selec, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_global2selec)) &
        allocate(domain%id_global2selec(var_size), stat = alloc_stat)
        domain%id_global2selec(1:var_size) = 0

    case (DomainSelec)

      if (allocated(domain%id_s2l)) then
        if (size(domain%id_s2l(2,:)) /= var_size) &
          deallocate(domain%id_s2l, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_s2l)) &
        allocate(domain%id_s2l(2,var_size), stat = alloc_stat)

      domain%id_s2l(1:2,1:var_size) = 0

      ! allocation for selected atoms (IY)
      if (allocated(domain%selec_atom2cell)) then
        if (size(domain%selec_atom2cell(:)) /= var_size) &
          deallocate(domain%selec_atom2cell, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%selec_atom2cell)) &
        allocate(domain%selec_atom2cell(var_size), stat = alloc_stat)

      domain%selec_atom2cell(1:var_size) = 0

      if (allocated(domain%id_selec2global)) then
        if (size(domain%id_selec2global(:)) /= var_size) &
          deallocate(domain%id_selec2global, stat = dealloc_stat)
      end if

      if (.not. allocated(domain%id_selec2global)) &
        allocate(domain%id_selec2global(var_size), stat = alloc_stat)
        domain%id_selec2global(1:var_size) = 0

    case default

      call error_msg('Alloc_Domain> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain
  !> @brief        deallocate domain information
  !! @authors      JJ
  !! @@aram[inout] domain   : domain information
  !! @param[in]    variable : selected variable
  !
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
        deallocate(domain%cell_g2l,    &
                   domain%cell_g2b,    &
                   domain%cell_gxyz2l, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocal)

      if (allocated(domain%cell_l2g)) then
        deallocate(domain%cell_l2g, &
                   stat = dealloc_stat)
      end if

    case (DomainCellLocBou)

      if (allocated(domain%cell_l2gx)) then
        deallocate(domain%cell_l2gx,      &
                   domain%cell_l2gy,      &
                   domain%cell_l2gz,      &
                   domain%cell_l2gx_orig, &
                   domain%cell_l2gy_orig, &
                   domain%cell_l2gz_orig, &
                   stat = dealloc_stat)
      end if

    case (DomainCellBoundary)

      if (allocated(domain%cell_b2g)) then
        deallocate(domain%cell_b2g, &
                   stat = dealloc_stat)
      end if

    case (DomainDynvar)

      if (allocated(domain%coord)) then
        deallocate(domain%id_l2s,        &
                   domain%num_atom,      &
                   domain%num_atom_t0,   &
                   domain%coord,         &
                   stat = dealloc_stat)
      end if

    case (DomainDynvar_Atom)

      if (allocated(domain%id_l2s)) then
        deallocate(domain%id_l2s,        &
                   domain%coord,         &
                   stat = dealloc_stat)
      end if

    case (DomainGlobal)

      if (allocated(domain%id_global2selec)) then
        deallocate(domain%id_global2selec, &
                   stat = dealloc_stat)
      end if

    case (DomainSelec)

      if (allocated(domain%id_s2l)) then
        deallocate(domain%id_s2l, &
                   domain%selec_atom2cell, &
                   domain%id_selec2global, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Domain> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_domain

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_domain_all
  !> @brief        deallocate all domain information
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
    call dealloc_domain(domain, DomainDynvar)
    call dealloc_domain(domain, DomainDynvar_Atom)
    call dealloc_domain(domain, DomainGlobal)
    call dealloc_domain(domain, DomainSelec)

    return

  end subroutine dealloc_domain_all

end module sa_domain_str_mod
