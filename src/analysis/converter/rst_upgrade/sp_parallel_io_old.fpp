!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_parallel_io_md
!> @brief   Parallel I/O module
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_parallel_io_old_mod

  use sp_parallel_io_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use fileio_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! parameters
  integer,       parameter      :: PioFileHeaderMarker = 12345

  ! global variables
  type(s_constraints),  public, save, target :: pio_constraints
  type(s_dynvars),      public, save, target :: pio_dynvars
  type(s_dynamics),     public, save, target :: pio_dynamics

  ! subroutines
  public  :: pio_read_domain_rst_old

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_domain_rst
  !> @brief        read restart information for domain
  !! @authors      NT
  !! @param[in]    filename    : restart filename 
  !! @param[inout] boundary    : boundary information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] constraints : constraints information (optional)
  !! @param[inout] dynvars     : dynvars information (optional)
  !! @param[inout] dynamics    : dynamics information (optional)
  !! @param[inout] t0_info     : t=0 information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_domain_rst_old(filename,    &
                                     boundary,    &
                                     domain,      &
                                     enefunc,     &
                                     constraints, &
                                     dynvars,     &
                                     dynamics,    &
                                     t0_info)
                             
    ! formal arguments
    character(*),                   intent(in)    :: filename
    type(s_boundary),               intent(inout) :: boundary
    type(s_domain),                 intent(inout) :: domain
    type(s_enefunc),                intent(inout) :: enefunc
    type(s_constraints),  optional, intent(inout) :: constraints
    type(s_dynvars),      optional, intent(inout) :: dynvars
    type(s_dynamics),     optional, intent(inout) :: dynamics
    type(s_pio_t0_info),  optional, intent(inout) :: t0_info

    ! local variables
    integer               :: file
    integer               :: i, j, ncell, nvar, nvar1, nvar2, nvar3
    integer               :: ix, icel, ig, nbyte
    integer               :: iseed_tmp
    character,allocatable :: bytes(:)
    real(wp)              :: lj_coef(2,MaxAtomCls)


    ! initialize
    !
    call init_boundary(boundary)
    call init_domain  (domain)
    call init_enefunc (enefunc)
    if (present(constraints)) call init_constraints(constraints)
    if (present(dynvars))     call init_dynvars    (dynvars)
    if (present(dynamics))    call init_dynamics   (dynamics)


    call pio_open_file(file, filename, IOFileInput)


    ! read constant values
    !
    read(file) ELECOEF


    ! read maximum value for allocatable variables
    !

    ! sp_domain_str
    read(file) MaxAtom
    read(file) MaxWater
    read(file) MaxMove
    read(file) MaxWaterMove

    ! sp_enefunc_str
    read(file) MaxBond
    read(file) MaxAngle
    read(file) MaxDihe
    read(file) MaxImpr
    read(file) MaxCmap
    read(file) BondMove
    read(file) AngleMove
    read(file) DiheMove
    read(file) ImprMove
    read(file) CmapMove

    ! sp_pairlist_str
    read(file) MaxNb15
    read(file) MaxNb15Water

    ! sp_constraints_str
    read(file) HGroupMax
    read(file) HGrpMaxMove


    ! read boundary information
    !

    read(file) boundary%box_size_x
    read(file) boundary%box_size_y
    read(file) boundary%box_size_z
    read(file) boundary%origin_x
    read(file) boundary%origin_y
    read(file) boundary%origin_z
    read(file) boundary%box_size_x_ref
    read(file) boundary%box_size_y_ref
    read(file) boundary%box_size_z_ref
    read(file) boundary%num_domain(1:3)
    read(file) boundary%num_cells_x
    read(file) boundary%num_cells_y
    read(file) boundary%num_cells_z


    ! read s_domain information
    !

    read(file) domain%num_deg_freedom
    read(file) domain%num_atom_all

    ! domain_str::Dynvar variables
    read(file) ncell
    call alloc_domain(domain, DomainDynvar, ncell, 1, 1)
    read(file) domain%num_atom   (1:ncell)
    read(file) domain%num_atom_t0(1:ncell)
    read(file) domain%num_solute (1:ncell)
    read(file) domain%num_water  (1:ncell)

    do i = 1, ncell

      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      read(file) domain%coord      (1:3, 1:nvar, i)
      read(file) domain%velocity   (1:3, 1:nvar, i)
      read(file) domain%trans_vec  (1:3, 1:nvar, i)
      read(file) domain%charge     (     1:nvar, i)
      read(file) domain%mass       (     1:nvar, i)
      read(file) domain%id_l2g     (     1:nvar, i)
      read(file) domain%atom_cls_no(     1:nvar, i)
      read(file) domain%solute_list(     1:nvar, i)

      ! 'MaxWater' arrays
      nvar = domain%num_water(i)
      read(file) domain%water_list (1:3, 1:nvar, i)

    end do

    ! domain_str::Global variables
    read(file) nvar
    call alloc_domain(domain, DomainGlobal, nvar, 1, 1)
    !! Re-compute id_g2l by id_l2g
    do icel = 1, ncell
      do ix = 1, domain%num_atom(icel)
        ig = domain%id_l2g(ix, icel)
        domain%id_g2l(1, ig) = icel
        domain%id_g2l(2, ig) = ix
      end do
    end do
    

    ! read s_enefunc information
    !

    ! enefunc_str::Bond variables
    ! enefunc_str::Angl variables
    ! enefunc_str::Dihe variables
    ! enefunc_str::RBDihe variables
    ! enefunc_str::Impr variables

    read(file) ncell
    call alloc_enefunc(enefunc, EneFuncBase,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncBond,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncAngl,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncDihe,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncRBDihe, ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncImpr,   ncell, ncell)
    read(file) enefunc%num_bond       (1:ncell)
    read(file) enefunc%num_angle      (1:ncell)
    read(file) enefunc%num_dihedral   (1:ncell)
    read(file) enefunc%num_rb_dihedral(1:ncell)
    read(file) enefunc%num_improper   (1:ncell)
    read(file) enefunc%notation_14types

    do i = 1, ncell

      nvar = enefunc%num_bond(i)
      read(file) enefunc%bond_list       (1:2, 1:nvar, i)
      read(file) enefunc%bond_force_const(     1:nvar, i)
      read(file) enefunc%bond_dist_min   (     1:nvar, i)

      nvar = enefunc%num_angle(i)
      read(file) enefunc%angle_list       (1:3, 1:nvar, i)
      read(file) enefunc%angle_force_const(     1:nvar, i)
      read(file) enefunc%angle_theta_min  (     1:nvar, i)
      read(file) enefunc%urey_force_const (     1:nvar, i)
      read(file) enefunc%urey_rmin        (     1:nvar, i)

      nvar = enefunc%num_dihedral(i)
      read(file) enefunc%dihe_list       (1:4, 1:nvar, i)
      read(file) enefunc%dihe_force_const(     1:nvar, i)
      read(file) enefunc%dihe_periodicity(     1:nvar, i)
      read(file) enefunc%dihe_phase      (     1:nvar, i)

      nvar = enefunc%num_rb_dihedral(i)
      read(file) enefunc%rb_dihe_list    (1:4, 1:nvar, i)
      read(file) enefunc%rb_dihe_c       (1:6, 1:nvar, i)

      nvar = enefunc%num_improper(i)
      read(file) enefunc%impr_list       (1:4, 1:nvar, i)
      read(file) enefunc%impr_force_const(     1:nvar, i)
      read(file) enefunc%impr_periodicity(     1:nvar, i)
      read(file) enefunc%impr_phase      (     1:nvar, i)

    end do

    ! enefunc_str::EneFuncRest variables

    read(file) enefunc%restraint

    read(file) ncell
    call alloc_enefunc(enefunc, EneFuncRest, ncell, ncell)
    read(file) enefunc%num_restraint(1:ncell)

    do i = 1, ncell

      nvar = enefunc%num_restraint(i)
      read(file) enefunc%restraint_atom (     1:nvar, i)
      read(file) enefunc%restraint_force(1:4, 1:nvar, i)
      read(file) enefunc%restraint_coord(1:3, 1:nvar, i)

    end do

    ! enefunc_str::Cmap variables

    read(file) nvar1
    read(file) nvar2
    call alloc_enefunc(enefunc, EneFuncCmap, ncell, nvar1, nvar2)
    read(file) enefunc%num_cmap(1:ncell)
    
    do i = 1, ncell

      nvar = enefunc%num_cmap(i)
      read(file) enefunc%cmap_list  (1:8, 1:nvar, i)
      read(file) enefunc%cmap_type  (     1:nvar, i)

    end do

    read(file) enefunc%cmap_resolution(1:nvar2)
    read(file) enefunc%cmap_coef      (1:4, 1:4, 1:nvar1, 1:nvar1, 1:nvar2)


    ! enefunc_str::nonbond parameters
    
    read(file) nvar
    call alloc_enefunc(enefunc, EneFuncNbon, nvar)

    read(file) enefunc%nb14_lj12(1:nvar,1:nvar)
    read(file) enefunc%nb14_lj6 (1:nvar,1:nvar)
    read(file) enefunc%nonb_lj12(1:nvar,1:nvar)
    read(file) enefunc%nonb_lj6 (1:nvar,1:nvar)
    read(file) lj_coef          (1:2,   1:nvar)

    enefunc%num_atom_cls = nvar

    ! enefunc_str:: AMBERScale
    read(file) nvar
    call alloc_enefunc(enefunc, EneFuncAMBERScale, nvar)

    read(file) enefunc%dihe_scnb(1:nvar)
    read(file) enefunc%dihe_scee(1:nvar)

    ! enefunc_str::other variables

    read(file) enefunc%table%atom_cls_no_O
    read(file) enefunc%table%atom_cls_no_H
    read(file) enefunc%table%charge_O
    read(file) enefunc%table%charge_H
    read(file) enefunc%table%mass_O
    read(file) enefunc%table%mass_H
    read(file) enefunc%table%water_model
    read(file) enefunc%fudge_lj
    read(file) enefunc%fudge_qq
    read(file) enefunc%excl_level


    ! read s_constraints information
    !
    if (present(constraints)) then

      read(file) constraints%connect
      read(file) constraints%water_rHH
      read(file) constraints%water_rOH
      read(file) constraints%water_massO
      read(file) constraints%water_massH

      ! domain_constraint_str::ConstraintsDomainBond
      read(file) nvar
      read(file) nvar1
      read(file) nvar2  ! == nvar1 + 1

      call alloc_constraints(constraints, ConstraintsDomainBond, nvar, nvar1)
      read(file) constraints%No_HGr   (         1:nvar)
      read(file) constraints%HGr_local(1:nvar1, 1:nvar)
  
      do i = 1, nvar
        do j = 1, nvar1
          read(file) nvar3
          read(file) constraints%HGr_bond_list(1:nvar2, 1:nvar3, j, i)
          read(file) constraints%HGr_bond_dist(1:nvar2, 1:nvar3, j, i)
        end do
      end do

    else

      read(file) pio_constraints%connect
      read(file) pio_constraints%water_rHH
      read(file) pio_constraints%water_rOH
      read(file) pio_constraints%water_massO
      read(file) pio_constraints%water_massH

      ! domain_constraint_str::ConstraintsDomainBond
      read(file) nvar
      read(file) nvar1
      read(file) nvar2  ! == nvar1 + 1

      call alloc_constraints(pio_constraints,ConstraintsDomainBond,nvar,nvar1)
      read(file) pio_constraints%No_HGr   (         1:nvar)
      read(file) pio_constraints%HGr_local(1:nvar1, 1:nvar)
  
      do i = 1, nvar
        do j = 1, nvar1
          read(file) nvar3
          read(file) pio_constraints%HGr_bond_list(1:nvar2, 1:nvar3, j, i)
          read(file) pio_constraints%HGr_bond_dist(1:nvar2, 1:nvar3, j, i)
        end do
      end do

    end if


    ! dynvars
    !

    if (present(dynvars)) then

      read(file) dynvars%thermostat_momentum
      read(file) dynvars%barostat_momentum(1:3)

    else

      read(file) pio_dynvars%thermostat_momentum
      read(file) pio_dynvars%barostat_momentum(1:3)

    end if


    ! dynamics
    !

    if (present(dynamics)) then

      read(file) dynamics%iseed
      iseed_tmp = dynamics%iseed

    else

      read(file) pio_dynamics%iseed
      iseed_tmp = pio_dynamics%iseed

    end if


    ! random internal states
    !

    read(file) nbyte
    allocate(bytes(nbyte))
    read(file) bytes(1:nbyte)
    
    call random_init(iseed_tmp)
    call random_stock_frombyte(bytes)
    call random_pull_stock
    
    deallocate(bytes)


    ! t0 information
    !

    if (present(t0_info)) then

      read(file) t0_info%input_topfile
      read(file) t0_info%input_parfile
      read(file) t0_info%energy_pairlistdist
      read(file) t0_info%energy_table
      read(file) t0_info%energy_watermodel
      read(file) t0_info%constraint_rigidbond
      read(file) t0_info%constraint_fastwater
      read(file) t0_info%constraint_watermodel
      read(file) t0_info%ensemble_type
      read(file) t0_info%boundary_boxsizex
      read(file) t0_info%boundary_boxsizey
      read(file) t0_info%boundary_boxsizez
      read(file) t0_info%boundary_originx
      read(file) t0_info%boundary_originy
      read(file) t0_info%boundary_originz
      read(file) t0_info%boundary_domain_xyz

      read(file) nvar
      if (allocated(t0_info%selection_group)) &
        deallocate(t0_info%selection_group)
      allocate(t0_info%selection_group(nvar))
      read(file) t0_info%selection_group(1:nvar)

      read(file) nvar
      if (allocated(t0_info%restraint_func))  &
        deallocate(t0_info%restraint_func,    &
                   t0_info%restraint_const,   &
                   t0_info%restraint_index)
      allocate(t0_info%restraint_func(nvar),  &
               t0_info%restraint_const(nvar), &
               t0_info%restraint_index(nvar))
      read(file) t0_info%restraint_func(1:nvar)
      read(file) t0_info%restraint_const(1:nvar)
      read(file) t0_info%restraint_index(1:nvar)

    else

      read(file) pio_t0_info%input_topfile
      read(file) pio_t0_info%input_parfile
      read(file) pio_t0_info%energy_pairlistdist
      read(file) pio_t0_info%energy_table
      read(file) pio_t0_info%energy_watermodel
      read(file) pio_t0_info%constraint_rigidbond
      read(file) pio_t0_info%constraint_fastwater
      read(file) pio_t0_info%constraint_watermodel
      read(file) pio_t0_info%ensemble_type
      read(file) pio_t0_info%boundary_boxsizex
      read(file) pio_t0_info%boundary_boxsizey
      read(file) pio_t0_info%boundary_boxsizez
      read(file) pio_t0_info%boundary_originx
      read(file) pio_t0_info%boundary_originy
      read(file) pio_t0_info%boundary_originz
      read(file) pio_t0_info%boundary_domain_xyz

      read(file) nvar
      if (allocated(pio_t0_info%selection_group)) &
        deallocate(pio_t0_info%selection_group)
      allocate(pio_t0_info%selection_group(nvar))
      read(file) pio_t0_info%selection_group(1:nvar)

      read(file) nvar
      if (allocated(pio_t0_info%restraint_func))  &
        deallocate(pio_t0_info%restraint_func,    &
                   pio_t0_info%restraint_const,   &
                   pio_t0_info%restraint_index)
      allocate(pio_t0_info%restraint_func(nvar),  &
               pio_t0_info%restraint_const(nvar), &
               pio_t0_info%restraint_index(nvar))
      read(file) pio_t0_info%restraint_func(1:nvar)
      read(file) pio_t0_info%restraint_const(1:nvar)
      read(file) pio_t0_info%restraint_index(1:nvar)

    end if


    ! read restart information
    !       (Domain restart file is START or RESTART)
    read(file) pio_restart


    call close_file(file)

    return

  end subroutine pio_read_domain_rst_old

end module sp_parallel_io_old_mod
