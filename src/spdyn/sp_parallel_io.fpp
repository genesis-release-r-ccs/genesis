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

module sp_parallel_io_mod

  use sp_remd_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_restraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use fileio_data_mod
  use fileio_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters
  integer,       parameter      :: PioFileHeaderMarker = 12345

  ! structures
  type, public :: s_pio_t0_info

    ! input files
    character(200)              :: input_topfile
    character(200)              :: input_parfile

    ! energy settings
    real(wp)                    :: energy_pairlistdist
    logical                     :: energy_table
    character(5)                :: energy_watermodel
    
    ! constraint settings
    logical                     :: constraint_rigidbond
    logical                     :: constraint_fastwater
    character(5)                :: constraint_watermodel

    ! ensemble settings
    integer                     :: ensemble_type

    ! boundary settings
    real(wp)                    :: boundary_boxsizex
    real(wp)                    :: boundary_boxsizey
    real(wp)                    :: boundary_boxsizez
    real(wp)                    :: boundary_originx
    real(wp)                    :: boundary_originy
    real(wp)                    :: boundary_originz
    integer                     :: boundary_domain_xyz

    ! selection settings    
    character(1000),allocatable :: selection_group(:)

    ! restraint settings
    integer,        allocatable :: restraint_func(:)
    character(256), allocatable :: restraint_const(:)
    character(256), allocatable :: restraint_index(:)

  end type s_pio_t0_info

  ! global variables
  logical,              public  :: pio_restart  = .false.

  type(s_pio_t0_info),  public, save, target :: pio_t0_info

  logical                       :: pio_crdtrj_hdr = .true.
  logical                       :: pio_veltrj_hdr = .true.

  integer,          allocatable :: tmp_atom_idx(:)
  real(sp),         allocatable :: tmp_atom_crd(:,:)


  ! subroutines
  public  :: pio_write_domain_trj
  public  :: pio_write_selection
  public  :: pio_write_domain_rst
  public  :: pio_read_selection
  public  :: pio_read_selection_blank
  public  :: pio_read_domain_rst
  public  :: pio_read_domain_str
  public  :: pio_check_ranked_file
  public  :: pio_get_ranked_filename
  public  :: pio_open_file
  private :: pio_get_file_endian

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_write_domain_trj
  !> @brief        write trajectory for domain
  !! @authors      NT
  !! @param[in]    unit_no    : file unit number
  !! @param[in]    boundary   : boundary information
  !! @param[in]    domain     : domain information
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    snapshot   : number of snapshot
  !! @param[in]    crd_trj    : flag for coordinate or trajectory
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_write_domain_trj(unit_no,    &
                                  boundary,   &
                                  domain,     &
                                  trajectory, &
                                  snapshot,   &
                                  crd_trj)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    real(wip),               intent(in)    :: trajectory(:,:,:)
    integer,                 intent(in)    :: snapshot
    logical,                 intent(in)    :: crd_trj

    ! local variables
    real(sp)                 :: box(6)
    integer                  :: i, ix, natom_domain


    ! write header
    !
    if (      crd_trj .and. pio_crdtrj_hdr .or. &
        .not. crd_trj .and. pio_veltrj_hdr) then

      ! my rank number
      write(unit_no) my_world_rank
      ! # of all atoms
      write(unit_no) domain%num_atom_all
      ! # of snapshot
      write(unit_no) snapshot

      if (crd_trj) then
        pio_crdtrj_hdr = .false.
      else
        pio_veltrj_hdr = .false.
      end if

    end if

    if (.not. allocated(tmp_atom_idx)) &
      allocate(tmp_atom_idx(domain%num_cell_local * MaxAtom))

    if (.not. allocated(tmp_atom_crd)) &
      allocate(tmp_atom_crd(3,domain%num_cell_local * MaxAtom))

    ! write body
    !
    natom_domain = 1
    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        tmp_atom_idx(    natom_domain) = domain%id_l2g(ix,i)
        tmp_atom_crd(1:3,natom_domain) = trajectory(1:3,ix,i)
        natom_domain = natom_domain + 1
      end do
    end do

    natom_domain = natom_domain - 1

    box = 0.0_sp
    box(1) = boundary%box_size_x_ref
    box(3) = boundary%box_size_y_ref
    box(6) = boundary%box_size_z_ref

    ! box information
    write(unit_no) box
    ! # of domain atoms
    write(unit_no) natom_domain
    ! global atom indices
    write(unit_no) tmp_atom_idx(    1:natom_domain)
    ! atom coordinates
    write(unit_no) tmp_atom_crd(1:3,1:natom_domain)

    return

  end subroutine pio_write_domain_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_write_selection
  !> @brief        write selection list
  !! @authors      JJ
  !! @param[in]    filename    : selection filename 
  !! @param[in]    restraints  : selection list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_write_selection(filename, restraints)
                             
    ! formal arguments
    character(*),                   intent(in) :: filename
    type(s_restraints),             intent(in) :: restraints 

    ! local variables
    integer               :: file
    integer               :: i, j, ngroups, max_natom, nvar

    !$omp critical
    call open_data(filename, IOFileDataWriteAscii, file)
    !$omp end critical

    ! write selection indices
    !
    ngroups = restraints%num_groups
    max_natom = restraints%max_atoms
    call write_data_integer &
         (file, 'restraints:num_groups', ngroups)
    call write_data_integer &
         (file, 'restraints:max_atoms', max_natom)
    call write_data_integer_array        &
         (file, 'restraints:num_atoms',  &
             (/ngroups/), restraints%num_atoms(1:ngroups))
    do i = 1, ngroups
      nvar = restraints%num_atoms(i)
      call write_data_integer_array      &
           (file, 'restraints:atomlist', &
               (/nvar/), restraints%atomlist(1:nvar,i))
    end do

    !$omp critical
    call close_data(file)
    !$omp end critical

    return

  end subroutine pio_write_selection

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_write_domain_rst
  !> @brief        write restart information for domain
  !! @authors      NT
  !! @param[in]    filename    : restart filename 
  !! @param[in]    boundary    : boundary information
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    constraints : constraints information (optional)
  !! @param[in]    dynvars     : dynvars information (optional)
  !! @param[in]    dynamics    : dynamics information (optional)
  !! @param[in]    t0_info     : t=0 information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_write_domain_rst(filename,        &
                                  boundary,        &
                                  domain,          &
                                  enefunc,         &
                                  constraints,     &
                                  dynvars,         &
                                  dynamics,        &
                                  remd)
                             
    ! formal arguments
    character(*),                   intent(in) :: filename
    type(s_boundary),               intent(in) :: boundary
    type(s_domain),                 intent(in) :: domain
    type(s_enefunc),                intent(in) :: enefunc
    type(s_constraints),            intent(in) :: constraints
    type(s_dynvars),      optional, intent(in) :: dynvars
    type(s_dynamics),     optional, intent(in) :: dynamics
    type(s_remd),         optional, intent(in) :: remd

    ! local variables
    integer               :: file
    integer               :: i, j, ncell, ncell_local, nvar, nvar1, nvar2, nvar3
    integer               :: nrep, ndim
    logical               :: write_file
    character,allocatable :: bytes(:)


    !$omp critical
    if (get_extension(filename) /= 'rsa') then
      call open_data(filename, IOFileDataWrite, file)
    else
      call open_data(filename, IOFileDataWriteAscii, file)
    end if
    !$omp end critical


    ! write constant values
    !
    call write_data_real_wp &
         (file, 'ELECOEF', ELECOEF)

    ! write maximum value for allocatable variables
    !

    ! sp_domain_str
    call write_data_integer &
         (file, 'MaxAtom', MaxAtom)
    call write_data_integer &
         (file, 'MaxWater', MaxWater)
    call write_data_integer &
         (file, 'MaxMove', MaxMove)
    call write_data_integer &
         (file, 'MaxWaterMove', MaxWaterMove)

    ! sp_enefunc_str
    call write_data_integer &
         (file, 'MaxBond', MaxBond)
    call write_data_integer &
         (file, 'MaxAngle', MaxAngle)
    call write_data_integer &
         (file, 'MaxDihe', MaxDihe)
    call write_data_integer &
         (file, 'MaxImpr', MaxImpr)
    call write_data_integer &
         (file, 'MaxCmap', MaxCmap)
    call write_data_integer &
         (file, 'BondMove', BondMove)
    call write_data_integer &
         (file, 'AngleMove', AngleMove)
    call write_data_integer &
         (file, 'DiheMove', DiheMove)
    call write_data_integer &
         (file, 'ImprMove', ImprMove)
    call write_data_integer &
         (file, 'CmapMove', CmapMove)

    ! sp_pairlist_str
    call write_data_integer &
         (file, 'MaxNb15', MaxNb15)
    call write_data_integer &
         (file, 'MaxNb15Water', MaxNb15Water)

    ! sp_constraints_str
    call write_data_integer &
         (file, 'HGroupMax', HGroupMax)
    call write_data_integer &
         (file, 'HGrpMaxMove', HGrpMaxMove)

    ! write boundary information
    !
    call write_data_real_wip &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call write_data_real_wip &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call write_data_real_wip &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call write_data_real_wp &
         (file, 'boundary:origin_x', boundary%origin_x)
    call write_data_real_wp &
         (file, 'boundary:origin_y', boundary%origin_y)
    call write_data_real_wp &
         (file, 'boundary:origin_z', boundary%origin_z)
    call write_data_real_wip &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call write_data_real_wip &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call write_data_real_wip &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call write_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain(1:3))
    call write_data_real_wip &
         (file, 'boundary:cell_size_x', boundary%cell_size_x)
    call write_data_real_wip &
         (file, 'boundary:cell_size_y', boundary%cell_size_y)
    call write_data_real_wip &
         (file, 'boundary:cell_size_z', boundary%cell_size_z)
    call write_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call write_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call write_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)

    ! write tip4 information
    !
    call write_data_integer   &
         (file, 'constraints:tip4', constraints%water_type)

    ! write s_domain information
    !
    call write_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call write_data_integer &
         (file, 'domain:num_cell_local', domain%num_cell_local)

    ! domain_str::Dynvar variables
    ncell = size(domain%num_atom(:))
    ncell_local = domain%num_cell_local

    call write_data_integer &
         (file, 'domain:ncell', ncell)

    call write_data_integer_array &
         (file, 'domain:num_atom', &
             (/ncell/), domain%num_atom  (1:ncell))
    call write_data_integer_array &
         (file, 'domain:num_solute', &
             (/ncell/), domain%num_solute(1:ncell))
    call write_data_integer_array &
         (file, 'domain:num_water', &
             (/ncell/), domain%num_water (1:ncell))
    call write_data_integer_array &
         (file, 'domain:cell_l2g', &
             (/ncell_local/), domain%cell_l2g (1:ncell_local))
    call write_data_integer_array &
         (file, 'domain:cell_b2g', &
             (/ncell-ncell_local/), domain%cell_b2g (1:ncell-ncell_local))

    do i = 1, ncell

      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      call write_data_real_wip_array &
         (file, 'domain:coord', &
             (/3,nvar,1/), domain%coord      (1:3, 1:nvar, i))
      call write_data_real_wip_array &
         (file, 'domain:velocity', &
             (/3,nvar,1/), domain%velocity   (1:3, 1:nvar, i))
      call write_data_real_wp_array &
         (file, 'domain:charge', &
             (/nvar,1/),   domain%charge     (1:nvar, i))
      call write_data_real_wip_array &
         (file, 'domain:mass', &
             (/nvar,1/),   domain%mass       (1:nvar, i))
      call write_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g     (1:nvar, i))
      call write_data_integer_array &
         (file, 'domain:id_l2g_solute', &
             (/nvar,1/),   domain%id_l2g_solute(1:nvar, i))
      call write_data_integer_array &
         (file, 'domain:atom_cls_no', &
             (/nvar,1/),   domain%atom_cls_no(1:nvar, i))

      nvar = domain%num_water(i)
      if (nvar == 0) &
        cycle
      if (constraints%water_type == TIP4) then
        nvar1 = 4
      else if (constraints%water_type == TIP3) then
        nvar1 = 3
      else if (constraints%water_type == TIP1) then
        nvar1 = 1
      end if
      call write_data_integer_array &
         (file, 'domain:water_list', &
             (/nvar1,nvar,1/), domain%water_list (1:nvar1,1:nvar,i))

    end do

    ! write s_enefunc information
    !
    call write_data_integer &
         (file, 'enefunc:num_solute', enefunc%table%num_solute)
    call write_data_real_wp &
         (file, 'enefunc:OH_bond', enefunc%table%OH_bond)
    call write_data_real_wp &
         (file, 'enefunc:OH_force', enefunc%table%OH_force)
    call write_data_real_wp &
         (file, 'enefunc:HOH_angle', enefunc%table%HOH_angle)
    call write_data_real_wp &
         (file, 'enefunc:HOH_force', enefunc%table%HOH_force)

    ncell = domain%num_cell_local

    ! enefunc_str::Bond variables
    ! enefunc_str::Angl variables
    ! enefunc_str::Dihe variables
    ! enefunc_str::Impr variables

    call write_data_integer_array &
         (file, 'enefunc:num_bond', &
             (/ncell/), enefunc%num_bond       (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_angle', &
             (/ncell/), enefunc%num_angle      (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_dihedral', &
             (/ncell/), enefunc%num_dihedral   (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_rb_dihedral', &
             (/ncell/), enefunc%num_rb_dihedral(1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_improper', &
             (/ncell/), enefunc%num_improper   (1:ncell))
    call write_data_integer &
         (file, 'enefunc:notation_14types', enefunc%notation_14types)

    do i = 1, ncell

      nvar = enefunc%num_bond(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:bond_list', &
                 (/2,nvar,1/), enefunc%bond_list     (1:2, 1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:bond_force_const', &
                 (/nvar,1/), enefunc%bond_force_const(     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:bond_dist_min', &
                 (/nvar,1/), enefunc%bond_dist_min   (     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:bond_pbc', &
                 (/nvar,1/), enefunc%bond_pbc        (     1:nvar, i))
      end if

      nvar = enefunc%num_angle(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:angle_list', &
                 (/3,nvar,1/), enefunc%angle_list     (1:3, 1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:angle_force_const', &
                 (/nvar,1/), enefunc%angle_force_const(     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:angle_theta_min', &
                 (/nvar,1/), enefunc%angle_theta_min  (     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:urey_force_const', &
                 (/nvar,1/), enefunc%urey_force_const (     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:urey_rmin', &
                 (/nvar,1/), enefunc%urey_rmin        (     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:angle_pbc', &
                 (/3,nvar,1/), enefunc%angle_pbc      (1:3, 1:nvar, i))
      end if

      nvar = enefunc%num_dihedral(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:dihe_list', &
                 (/4,nvar,1/), enefunc%dihe_list     (1:4, 1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:dihe_force_const', &
                 (/nvar,1/), enefunc%dihe_force_const(     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:dihe_periodicity', &
                 (/nvar,1/), enefunc%dihe_periodicity(     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:dihe_phase', &
                 (/nvar,1/), enefunc%dihe_phase      (     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:dihe_pbc', &
                 (/3,nvar,1/), enefunc%dihe_pbc      (1:3, 1:nvar, i))
      end if

      nvar = enefunc%num_rb_dihedral(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:rb_dihe_list', &
                 (/4,nvar,1/), enefunc%rb_dihe_list(1:4, 1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:rb_dihe_c', &
                 (/6,nvar,1/), enefunc%rb_dihe_c   (1:6, 1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:rb_dihe_pbc', &
                 (/3,nvar,1/), enefunc%rb_dihe_pbc (1:3, 1:nvar, i))
      end if

      nvar = enefunc%num_improper(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:impr_list', &
                 (/4,nvar,1/), enefunc%impr_list     (1:4, 1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:impr_force_const', &
                 (/nvar,1/), enefunc%impr_force_const(     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:impr_periodicity', &
                 (/nvar,1/), enefunc%impr_periodicity(     1:nvar, i))
        call write_data_real_wp_array &
             (file, 'enefunc:impr_phase', &
                 (/nvar,1/), enefunc%impr_phase      (     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:impr_pbc', &
                 (/3,nvar,1/), enefunc%impr_pbc      (1:3, 1:nvar, i))
      end if

    end do

    ! enefunc_str::EneFuncRest variables

    call write_data_logical &
         (file, 'enefunc:restraint', enefunc%restraint)

    call write_data_integer_array &
         (file, 'enefunc:num_restraint', &
             (/ncell/), enefunc%num_restraint(1:ncell))

    do i = 1, ncell

      nvar = enefunc%num_restraint(i)
      if (nvar == 0) &
        cycle
      call write_data_integer_array &
           (file, 'enefunc:restraint_atom', &
               (/nvar,1/),    enefunc%restraint_atom (     1:nvar, i))
      call write_data_real_wp_array &
           (file, 'enefunc:restraint_force', &
               (/4,nvar,1/),  enefunc%restraint_force(1:4, 1:nvar, i))
      call write_data_real_wp_array &
           (file, 'enefunc:restraint_coord', &
               (/3,nvar,1/),  enefunc%restraint_coord(1:3, 1:nvar, i))

    end do

    ! enefunc_str::Cmap variables

    if (allocated(enefunc%cmap_coef)) then
      nvar1 = size(enefunc%cmap_coef(1,1,:,1,1))
      nvar2 = size(enefunc%cmap_coef(1,1,1,1,:))
    else
      nvar1 = 0
      nvar2 = 0
    end if

    call write_data_integer &
         (file, 'enefunc:cmap_ngrid0', nvar1)
    call write_data_integer &
         (file, 'enefunc:cmap_ncmap_p', nvar2)

    call write_data_integer_array &
         (file, 'enefunc:num_cmap', &
             (/ncell/), enefunc%num_cmap(1:ncell))

    do i = 1, ncell

      nvar = enefunc%num_cmap(i)
      if (nvar == 0) &
        cycle
      call write_data_integer_array &
           (file, 'enefunc:cmap_list', &
               (/8,nvar,1/), enefunc%cmap_list(1:8, 1:nvar, i))
      call write_data_integer_array &
           (file, 'enefunc:cmap_type', &
               (/nvar,1/),   enefunc%cmap_type(     1:nvar, i))
      call write_data_integer_array &
           (file, 'enefunc:cmap_pbc', &
               (/6,nvar,1/), enefunc%cmap_pbc (1:6, 1:nvar, i))

    end do

    if (nvar1 > 0 .and. nvar2 > 0) then
      call write_data_integer_array &
           (file, 'enefunc:cmap_resolution', &
               (/nvar2/), enefunc%cmap_resolution(1:nvar2))
      call write_data_real_wp_array &
           (file, 'enefunc:cmap_coef', &
               (/4,4,nvar1,nvar1,nvar2/), &
                   enefunc%cmap_coef(1:4, 1:4, 1:nvar1, 1:nvar1, 1:nvar2))
    end if

    ! enefunc_str::nonbond parameters
    
    nvar = enefunc%num_atom_cls

    call write_data_integer &
         (file, 'enefunc:num_atom_cls', nvar)

    if (nvar > 0) then
      call write_data_real_wp_array &
           (file, 'enefunc:nb14_lj12', &
               (/nvar,nvar/), enefunc%nb14_lj12(1:nvar,1:nvar))
      call write_data_real_wp_array &
           (file, 'enefunc:nb14_lj6', &
               (/nvar,nvar/), enefunc%nb14_lj6 (1:nvar,1:nvar))
      call write_data_real_wp_array &
           (file, 'enefunc:nonb_lj12', &
               (/nvar,nvar/), enefunc%nonb_lj12(1:nvar,1:nvar))
      call write_data_real_wp_array &
           (file, 'enefunc:nonb_lj6', &
               (/nvar,nvar/), enefunc%nonb_lj6 (1:nvar,1:nvar))
      call write_data_real_wp_array &
           (file, 'enefunc:nonb_lj6_factor', &
               (/     nvar/), enefunc%nonb_lj6_factor(1:nvar))
    end if

    ! enefunc_str:: AMBERScale
    if (allocated(enefunc%dihe_scnb)) then
      nvar = size(enefunc%dihe_scnb)
    else
      nvar = 0
    end if

    if (nvar > 0) then
      call write_data_real_wp_array &
           (file, 'enefunc:dihe_scnb', &
               (/nvar/), enefunc%dihe_scnb(0:nvar-1))
      call write_data_real_wp_array &
           (file, 'enefunc:dihe_scee', &
               (/nvar/), enefunc%dihe_scee(0:nvar-1))

    end if

    ! enefunc_str::other variables

    call write_data_integer &
         (file, 'enefunc:table:num_water', enefunc%table%num_water)
    call write_data_integer &
         (file, 'enefunc:table:atom_cls_no_O', enefunc%table%atom_cls_no_O)
    call write_data_integer &
         (file, 'enefunc:table:atom_cls_no_H', enefunc%table%atom_cls_no_H)
    if (constraints%water_type == TIP4) &
      call write_data_integer &
           (file, 'enefunc:table:atom_cls_no_D', enefunc%table%atom_cls_no_D)
    call write_data_real_wp &
         (file, 'enefunc:table:charge_O', enefunc%table%charge_O)
    call write_data_real_wp &
         (file, 'enefunc:table:charge_H', enefunc%table%charge_H)
    if (constraints%water_type == TIP4) &
      call write_data_real_wp &
           (file, 'enefunc:table:charge_D', enefunc%table%charge_D)
    call write_data_real_wp &
         (file, 'enefunc:table:mass_O', enefunc%table%mass_O)
    call write_data_real_wp &
         (file, 'enefunc:table:mass_H', enefunc%table%mass_H)
    if (constraints%water_type == TIP4) &
      call write_data_real_wp &
           (file, 'enefunc:table:mass_D', enefunc%table%mass_D)
    call write_data_byte_array &
         (file, 'enefunc:table:water_model', (/5/), enefunc%table%water_model)
    call write_data_real_wp &
         (file, 'enefunc:table:fudge_lj', enefunc%fudge_lj)
    call write_data_real_wp &
         (file, 'enefunc:table:fudge_qq', enefunc%fudge_qq)
    call write_data_integer &
         (file, 'enefunc:table:excl_level', enefunc%excl_level)


    ! write s_constraints information
    !
    call write_data_integer &
         (file, 'constraints:connect', constraints%connect)
    call write_data_real_wip &
         (file, 'constraints:water_rHH', constraints%water_rHH)
    call write_data_real_wip &
         (file, 'constraints:water_rOH', constraints%water_rOH)
    if (constraints%water_type == TIP4) &
      call write_data_real_wip &
           (file, 'constraints:water_rOD', constraints%water_rOD)
    call write_data_real_wip &
         (file, 'constraints:water_massO', constraints%water_massO)
    call write_data_real_wip &
         (file, 'constraints:water_massH', constraints%water_massH)

    ! domain_constraint_str::ConstraintsDomainBond
    nvar  = size(constraints%No_HGr(:))
    nvar1 = size(constraints%HGr_local(:,1))
    nvar2 = size(constraints%HGr_bond_list(:,1,1,1))

    call write_data_integer &
         (file, 'constraints:No_HGr_size', nvar)
    call write_data_integer &
         (file, 'constraints:HGr_local_size', nvar1)
    call write_data_integer &
         (file, 'constraints:HGr_bond_list_size', nvar2)

    if (nvar > 0 .and. nvar1 > 0) then

      call write_data_integer_array &
           (file, 'constraints:No_HGr', &
               (/nvar/),       constraints%No_HGr   (         1:nvar))
      call write_data_integer_array &
           (file, 'constraints:HGr_local', &
               (/nvar1,nvar/), constraints%HGr_local(1:nvar1, 1:nvar))

    end if

    do i = 1, nvar
      do j = 1, nvar1

        nvar3 = constraints%HGr_local(j,i)
        if (nvar2 == 0 .or. nvar3 == 0) &
          cycle

        call write_data_integer_array &
             (file, 'constraints:HGr_bond_list', &
       (/nvar2,nvar3,1,1/), constraints%HGr_bond_list(1:nvar2, 1:nvar3, j, i))
        call write_data_real_wip_array &
             (file, 'constraints:HGr_bond_dist', &
       (/nvar2,nvar3,1,1/), constraints%HGr_bond_dist(1:nvar2, 1:nvar3, j, i))
      end do
    end do

    ! dynvars information
    !
    if (present(dynvars)) then

      write_file = .true.
      call write_data_logical &
         (file, 'dynvar:read_dynvar', write_file)
      call write_data_real_wip &
           (file, 'dynvars:thermostat_momentum', dynvars%thermostat_momentum)
      call write_data_real_wip_array &
           (file, 'dynvars:barostat_momentum', &
               (/3/), dynvars%barostat_momentum(1:3))

    else
      write_file = .false.
      call write_data_logical &
         (file, 'dynvar:read_dynvar', write_file)

    end if

    ! dynamics information
    !
    if (present(dynamics)) then

      write_file = .true.
      call write_data_logical &
           (file, 'dynamics:read_iseed', write_file)
      call write_data_integer &
           (file, 'dynamics:iseed', dynamics%iseed)
      call write_data_logical &
           (file, 'pio_restart', pio_restart)

    else
      write_file = .false.
      call write_data_logical &
           (file, 'dynamics:read_iseed', write_file)
    end if

    ! random internal states
    !
    call random_stock_tobyte(bytes)

    call write_data_byte_array &
         (file, 'random', (/size(bytes)/), bytes)

    deallocate(bytes)

    ! remd information
    !
    if (present(remd)) then
      write_file = .true.
      call write_data_logical &
         (file, 'remd:read_remd', write_file)
      nrep = remd%total_nreplicas
      ndim = remd%dimension
      call write_data_integer &
           (file, 'remd:iseed', remd%iseed)
      call write_data_integer &
           (file, 'remd:total_nreplicas', nrep)
      call write_data_integer &
           (file, 'remd:dimension', ndim)
      call write_data_integer_array          &
           (file, 'remd:repid2parmsetid',    &
                  (/nrep/), remd%repid2parmsetid(1:nrep))
      call write_data_integer_array          &
           (file, 'remd:num_criteria',       &
                  (/nrep,ndim,2/), remd%num_criteria(1:nrep,1:ndim,1:2))
      call write_data_integer_array          &
           (file, 'remd:num_exchanges',      &
                  (/nrep,ndim,2/), remd%num_exchanges(1:nrep,1:ndim,1:2))
    else
      write_file = .false.
      call write_data_logical &
         (file, 'remd:read_remd', write_file)
    end if

    !$omp critical
    call close_data(file)
    !$omp end critical

    return

  end subroutine pio_write_domain_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_selection
  !> @brief        read selection list
  !! @authors      JJ
  !! @param[in]    filename    : selection filename 
  !! @param[in]    restraints  : selection list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_selection(filename, restraints)

    ! formal arguments
    character(*),                   intent(in)    :: filename
    type(s_restraints),             intent(inout) :: restraints

    ! local variables
    integer               :: file
    integer               :: i, j, ngroup, max_natom, nvar

    call open_data(filename, IOFileDataRead, file)

    call read_data_integer &
         (file, 'restraints:num_groups', ngroup)
    call read_data_integer &
         (file, 'restraints:max_atoms', max_natom)

    restraints%num_groups = ngroup
    restraints%max_atoms  = max_natom

    call alloc_restraints(restraints, RestraintsList, ngroup, max_natom)

    call read_data_integer_array        &
         (file, 'restraints:num_atoms',  &
             (/ngroup/), restraints%num_atoms(1:ngroup))
    do i = 1, ngroup
      nvar = restraints%num_atoms(i)
      call read_data_integer_array      &
           (file, 'restraints:atomlist', &
               (/nvar/), restraints%atomlist(1:nvar,i))
    end do

    call close_data(file)

    return

  end subroutine pio_read_selection

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_selection_blank
  !> @brief        read selection list
  !! @authors      JJ
  !! @param[in]    restraints  : selection list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_selection_blank(restraints)

    ! formal arguments
    type(s_restraints),             intent(inout) :: restraints

    ! local variables
    integer               :: i, j, ngroup, max_natom, nvar

    ngroup = 0
    restraints%num_groups = ngroup

    call alloc_restraints(restraints, RestraintsList, 1, 1)

    restraints%num_atoms(1) = 0
    restraints%atomlist(1,1) = 0

    return

  end subroutine pio_read_selection_blank

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
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_domain_rst(filename,     &
                                 filenum,      &
                                 file_tot_num, &
                                 boundary,     &
                                 domain,       &
                                 enefunc,      &
                                 constraints,  &
                                 dynvars,      &
                                 dynamics,     &
                                 remd)
                             
    ! formal arguments
    character(*),                   intent(in)    :: filename
    integer,                        intent(in)    :: filenum
    integer,                        intent(in)    :: file_tot_num
    type(s_boundary),               intent(inout) :: boundary
    type(s_domain),                 intent(inout) :: domain
    type(s_enefunc),                intent(inout) :: enefunc
    type(s_constraints),            intent(inout) :: constraints
    type(s_dynvars),     optional,  intent(inout) :: dynvars
    type(s_dynamics),    optional,  intent(inout) :: dynamics
    type(s_remd),        optional,  intent(inout) :: remd

    ! local variables
    integer               :: file
    integer               :: i, j, ncell, ncell_local 
    integer               :: nvar, nvar1, nvar2, nvar3
    integer               :: nrep, ndim
    integer               :: ix, icel, ig, nbyte
    integer               :: iseed_tmp = 0
    integer               :: temp(6)
    logical               :: read_file
    character,allocatable :: bytes(:)


    ! initialize
    !
    call init_domain  (domain)
    call init_enefunc (enefunc)
    call init_constraints(constraints)
    if (present(dynvars)) call init_dynvars    (dynvars)
    if (present(dynamics)) call init_dynamics(dynamics)


    call open_data(filename, IOFileDataRead, file)


    ! read constant values
    !
    call read_data_real_wp &
         (file, 'ELECOEF', ELECOEF)


    ! read maximum value for allocatable variables
    !

    ! sp_domain_str
    call read_data_integer &
         (file, 'MaxAtom', MaxAtom)
    inv_MaxAtom = 1.0_wp / real(MaxAtom,wp)
    call read_data_integer &
         (file, 'MaxWater', MaxWater)
    call read_data_integer &
         (file, 'MaxMove', MaxMove)
    call read_data_integer &
         (file, 'MaxWaterMove', MaxWaterMove)

    ! sp_enefunc_str
    call read_data_integer &
         (file, 'MaxBond', MaxBond)
    call read_data_integer &
         (file, 'MaxAngle', MaxAngle)
    call read_data_integer &
         (file, 'MaxDihe', MaxDihe)
    call read_data_integer &
         (file, 'MaxImpr', MaxImpr)
    call read_data_integer &
         (file, 'MaxCmap', MaxCmap)
    call read_data_integer &
         (file, 'BondMove', BondMove)
    call read_data_integer &
         (file, 'AngleMove', AngleMove)
    call read_data_integer &
         (file, 'DiheMove', DiheMove)
    call read_data_integer &
         (file, 'ImprMove', ImprMove)
    call read_data_integer &
         (file, 'CmapMove', CmapMove)

    ! sp_pairlist_str
    call read_data_integer &
         (file, 'MaxNb15', MaxNb15)
    call read_data_integer &
         (file, 'MaxNb15Water', MaxNb15Water)

    ! sp_constraints_str
    call read_data_integer &
         (file, 'HGroupMax', HGroupMax)
    call read_data_integer &
         (file, 'HGrpMaxMove', HGrpMaxMove)

    ! read boundary information
    !
    call read_data_real_wip &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call read_data_real_wip &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call read_data_real_wip &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call read_data_real_wp &
         (file, 'boundary:origin_x', boundary%origin_x)
    call read_data_real_wp &
         (file, 'boundary:origin_y', boundary%origin_y)
    call read_data_real_wp &
         (file, 'boundary:origin_z', boundary%origin_z)
    call read_data_real_wip &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call read_data_real_wip &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call read_data_real_wip &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call read_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain_pio(1:3))
    call read_data_real_wip &
         (file, 'boundary:cell_size_x', boundary%cell_size_x)
    call read_data_real_wip &
         (file, 'boundary:cell_size_y', boundary%cell_size_y)
    call read_data_real_wip &
         (file, 'boundary:cell_size_z', boundary%cell_size_z)
    call read_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call read_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call read_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)

    call pio_boundary_check(boundary)

    ! read tip4 information
    !
    call read_data_integer &
         (file, 'constraints:tip4', constraints%water_type)

    ! read s_domain information
    !
    call read_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call read_data_integer &
         (file, 'domain:num_cell_local', ncell_local)
    call read_data_integer &
         (file, 'domain:ncell', ncell)

    domain%ncell_pio = ncell
    domain%ncell_local_pio = ncell_local

    ! domain_str::Dynvar variables
    if (filenum == 1) &
    call alloc_domain(domain, DomainDynvar_pio, ncell, file_tot_num, 1)

    call read_data_integer_array       &
         (file, 'domain:num_atom',     &
             (/ncell/), domain%num_atom_pio(1:ncell,filenum))
    call read_data_integer_array       &
         (file, 'domain:num_solute',   &
             (/ncell/), domain%num_solute_pio(1:ncell,filenum))
    call read_data_integer_array       &
         (file, 'domain:num_water',    &
             (/ncell/), domain%num_water_pio(1:ncell,filenum))
    call read_data_integer_array       &
         (file, 'domain:cell_l2g',     &
             (/ncell_local/), domain%cell_l2g_pio(1:ncell_local,filenum))
    call read_data_integer_array       &
         (file, 'domain:cell_b2g',     &
             (/ncell-ncell_local/), domain%cell_b2g_pio(1:ncell-ncell_local,filenum))
    do i = ncell_local+1,ncell
      domain%cell_l2g_pio(i,filenum) = &
                   domain%cell_b2g_pio(i-ncell_local,filenum)
    end do

    if (constraints%water_type == TIP4) then
      nvar1 = 4
    else if (constraints%water_type == TIP3) then
      nvar1 = 3
    else if (constraints%water_type == TIP1) then
      nvar1 = 1
    end if
    if (filenum == 1) &
      call alloc_domain(domain, DomainDynvar_Atom_pio, ncell, &
                        nvar1, file_tot_num)

    do i = 1, ncell

      ! 'MaxAtom' arrays
      nvar = domain%num_atom_pio(i,filenum)
      call read_data_real_wip_array &
         (file, 'domain:coord', &
             (/3,nvar,1/), domain%coord_pio   (1:3, 1:nvar, i,filenum))
      call read_data_real_wip_array &
         (file, 'domain:velocity', &
             (/3,nvar,1/), domain%velocity_pio(1:3, 1:nvar, i,filenum))
      call read_data_real_wp_array &
         (file, 'domain:charge', &
             (/nvar,1/),   domain%charge_pio       (1:nvar, i,filenum))
      call read_data_real_wip_array &
         (file, 'domain:mass', &
             (/nvar,1/),   domain%mass_pio         (1:nvar, i,filenum))
      call read_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g_pio       (1:nvar, i,filenum))
      call read_data_integer_array &
         (file, 'domain:id_l2g_solute', &
             (/nvar,1/),   domain%id_l2g_solute_pio(1:nvar, i,filenum))
      call read_data_integer_array &
         (file, 'domain:atom_cls_no', &
             (/nvar,1/),   domain%atom_cls_no_pio  (1:nvar, i,filenum))

      nvar = domain%num_water_pio(i,filenum)
      call read_data_integer_array &
         (file, 'domain:water_list', &
             (/nvar1,nvar,1/), domain%water_list_pio(1:nvar1,1:nvar,i,filenum))
    end do

    ! read s_enefunc information
    !
    call read_data_integer &
         (file, 'enefunc:num_solute', enefunc%table%num_solute)
    call read_data_real_wp &
         (file, 'enefunc:OH_bond', enefunc%table%OH_bond)
    call read_data_real_wp &
         (file, 'enefunc:OH_force', enefunc%table%OH_force)
    call read_data_real_wp &
         (file, 'enefunc:HOH_angle', enefunc%table%HOH_angle)
    call read_data_real_wp &
         (file, 'enefunc:HOH_force', enefunc%table%HOH_force)

    ! enefunc_str::Bond variables
    ! enefunc_str::Angl variables
    ! enefunc_str::Dihe variables
    ! enefunc_str::RBDihe variables
    ! enefunc_str::Impr variables

    if (filenum == 1) then
      call alloc_enefunc(enefunc, EneFuncBase_pio,   ncell_local, file_tot_num)
      call alloc_enefunc(enefunc, EneFuncBond_pio,   ncell_local, file_tot_num)
      call alloc_enefunc(enefunc, EneFuncAngl_pio,   ncell_local, file_tot_num)
      call alloc_enefunc(enefunc, EneFuncDihe_pio,   ncell_local, file_tot_num)
      call alloc_enefunc(enefunc, EneFuncRBDihe_pio, ncell_local, file_tot_num)
      call alloc_enefunc(enefunc, EneFuncImpr_pio,   ncell_local, file_tot_num)
    end if

    call read_data_integer_array &
         (file, 'enefunc:num_bond', &
             (/ncell_local/), enefunc%num_bond_pio      (1:ncell_local,filenum))
    call read_data_integer_array &
         (file, 'enefunc:num_angle', &
             (/ncell_local/), enefunc%num_angle_pio     (1:ncell_local,filenum))
    call read_data_integer_array &
         (file, 'enefunc:num_dihedral', &
             (/ncell_local/), enefunc%num_dihedral_pio  (1:ncell_local,filenum))
    call read_data_integer_array &
         (file, 'enefunc:num_rb_dihedral', &
            (/ncell_local/), enefunc%num_rb_dihedral_pio(1:ncell_local,filenum))
    call read_data_integer_array &
         (file, 'enefunc:num_improper', &
             (/ncell_local/), enefunc%num_improper_pio  (1:ncell_local,filenum))
    call read_data_integer &
         (file, 'enefunc:notation_14types', enefunc%notation_14types)

    do i = 1, ncell_local

      nvar = enefunc%num_bond_pio(i,filenum)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:bond_list', &
              (/2,nvar,1/), enefunc%bond_list_pio     (1:2, 1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:bond_force_const', &
              (/nvar,1/), enefunc%bond_force_const_pio(     1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:bond_dist_min', &
              (/nvar,1/), enefunc%bond_dist_min_pio   (     1:nvar, i,filenum))
        call read_data_integer_array &
             (file, 'enefunc:bond_pbc', &
              (/nvar,1/), enefunc%bond_pbc_pio        (     1:nvar, i,filenum))
      end if

      nvar = enefunc%num_angle_pio(i,filenum)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:angle_list', &
              (/3,nvar,1/), enefunc%angle_list_pio     (1:3, 1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:angle_force_const', &
              (/nvar,1/), enefunc%angle_force_const_pio(     1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:angle_theta_min', &
              (/nvar,1/), enefunc%angle_theta_min_pio  (     1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:urey_force_const', &
              (/nvar,1/), enefunc%urey_force_const_pio (     1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:urey_rmin', &
              (/nvar,1/), enefunc%urey_rmin_pio        (     1:nvar, i,filenum))
        call read_data_integer_array &
             (file, 'enefunc:angle_pbc', &
              (/3,nvar,1/), enefunc%angle_pbc_pio      (1:3, 1:nvar, i,filenum))
      end if

      nvar = enefunc%num_dihedral_pio(i,filenum)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:dihe_list', &
              (/4,nvar,1/), enefunc%dihe_list_pio     (1:4, 1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:dihe_force_const', &
              (/nvar,1/), enefunc%dihe_force_const_pio(     1:nvar, i,filenum))
        call read_data_integer_array &
             (file, 'enefunc:dihe_periodicity', &
              (/nvar,1/), enefunc%dihe_periodicity_pio(     1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:dihe_phase', &
              (/nvar,1/), enefunc%dihe_phase_pio      (     1:nvar, i,filenum))
        call read_data_integer_array &
             (file, 'enefunc:dihe_pbc', &
              (/3,nvar,1/), enefunc%dihe_pbc_pio      (1:3, 1:nvar, i,filenum))
      end if

      nvar = enefunc%num_rb_dihedral_pio(i,filenum)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:rb_dihe_list', &
                 (/4,nvar,1/), enefunc%rb_dihe_list_pio(1:4, 1:nvar, i,filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:rb_dihe_c', &
                 (/6,nvar,1/), enefunc%rb_dihe_c_pio   (1:6, 1:nvar, i,filenum))
        call read_data_integer_array &
             (file, 'enefunc:rb_dihe_pbc', &
                 (/3,nvar,1/), enefunc%rb_dihe_pbc_pio (1:3, 1:nvar, i,filenum))
      end if

      nvar = enefunc%num_improper_pio(i, filenum)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:impr_list', &
              (/4,nvar,1/), enefunc%impr_list_pio     (1:4, 1:nvar, i, filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:impr_force_const', &
              (/nvar,1/), enefunc%impr_force_const_pio(     1:nvar, i, filenum))
        call read_data_integer_array &
             (file, 'enefunc:impr_periodicity', &
              (/nvar,1/), enefunc%impr_periodicity_pio(     1:nvar, i, filenum))
        call read_data_real_wp_array &
             (file, 'enefunc:impr_phase', &
              (/nvar,1/), enefunc%impr_phase_pio      (     1:nvar, i, filenum))
        call read_data_integer_array &
             (file, 'enefunc:impr_pbc', &
              (/3,nvar,1/), enefunc%impr_pbc_pio      (1:3, 1:nvar, i, filenum))
      end if

    end do


    ! enefunc_str::EneFuncRest variables

    call read_data_logical &
         (file, 'enefunc:restraint', enefunc%restraint)

    call read_data_integer_array &
         (file, 'enefunc:num_restraint', &
             (/ncell_local/), enefunc%num_restraint_pio(1:ncell_local, filenum))

    if (filenum == 1) &
      call alloc_enefunc(enefunc, EneFuncRest_pio, ncell_local, file_tot_num)

    do i = 1, ncell_local

      nvar = enefunc%num_restraint_pio(i,filenum)
      if (nvar == 0) &
        cycle
      call read_data_integer_array &
           (file, 'enefunc:restraint_atom', &
            (/nvar,1/),    enefunc%restraint_atom_pio (     1:nvar, i, filenum))
      call read_data_real_wp_array &
           (file, 'enefunc:restraint_force', &
            (/4,nvar,1/),  enefunc%restraint_force_pio(1:4, 1:nvar, i, filenum))
      call read_data_real_wp_array &
           (file, 'enefunc:restraint_coord', &
            (/3,nvar,1/),  enefunc%restraint_coord_pio(1:3, 1:nvar, i, filenum))

    end do

    ! enefunc_str::Cmap variables

    call read_data_integer &
         (file, 'enefunc:cmap_ngrid0', nvar1)
    call read_data_integer &
         (file, 'enefunc:cmap_ncmap_p', nvar2)

    enefunc%cmap_ngrid0  = nvar1
    enefunc%cmap_ncmap_p = nvar2

    call read_data_integer_array &
         (file, 'enefunc:num_cmap', &
             (/ncell_local/), enefunc%num_cmap_pio(1:ncell_local, filenum))

    if (filenum == 1) then
      call alloc_enefunc(enefunc, EneFuncCmapcoef_pio, nvar1, nvar2)
      call alloc_enefunc(enefunc, EneFuncCmap_pio, ncell_local, file_tot_num)
    end if

    do i = 1, ncell_local

      nvar = enefunc%num_cmap_pio(i,filenum)
      if (nvar == 0) &
        cycle
      call read_data_integer_array &
           (file, 'enefunc:cmap_list', &
               (/8,nvar,1/), enefunc%cmap_list_pio(1:8, 1:nvar, i,filenum))
      call read_data_integer_array &
           (file, 'enefunc:cmap_type', &
               (/nvar,1/),   enefunc%cmap_type_pio(     1:nvar, i,filenum))
      call read_data_integer_array &
           (file, 'enefunc:cmap_pbc', &
               (/6,nvar,1/), enefunc%cmap_pbc_pio (1:6, 1:nvar, i,filenum))

    end do

    if (nvar1 > 0 .and. nvar2 > 0) then
      call read_data_integer_array &
           (file, 'enefunc:cmap_resolution', &
               (/nvar2/), enefunc%cmap_resolution_pio(1:nvar2))
      call read_data_real_wp_array &
           (file, 'enefunc:cmap_coef', &
               (/4,4,nvar1,nvar1,nvar2/), &
                   enefunc%cmap_coef_pio(1:4,1:4,1:nvar1,1:nvar1,1:nvar2))
    end if

    ! enefunc_str::nonbond parameters
    
    call read_data_integer &
         (file, 'enefunc:num_atom_cls', enefunc%num_atom_cls)

    nvar = enefunc%num_atom_cls

    if (filenum == 1) &
      call alloc_enefunc(enefunc, EneFuncNbon, nvar)

    if (nvar > 0) then
      call read_data_real_wp_array &
           (file, 'enefunc:nb14_lj12', &
               (/nvar,nvar/), enefunc%nb14_lj12(1:nvar,1:nvar))
      call read_data_real_wp_array &
           (file, 'enefunc:nb14_lj6', &
               (/nvar,nvar/), enefunc%nb14_lj6 (1:nvar,1:nvar))
      call read_data_real_wp_array &
           (file, 'enefunc:nonb_lj12', &
               (/nvar,nvar/), enefunc%nonb_lj12(1:nvar,1:nvar))
      call read_data_real_wp_array &
           (file, 'enefunc:nonb_lj6', &
               (/nvar,nvar/), enefunc%nonb_lj6 (1:nvar,1:nvar))
      call read_data_real_wp_array &
           (file, 'enefunc:nonb_lj6_factor', &
               (/     nvar/), enefunc%nonb_lj6_factor(1:nvar))
    end if

    ! enefunc_str:: AMBERScale
    call get_data_size &
         (file, 'enefunc:dihe_scnb', nvar, .false.)

    if (filenum == 1) &
      call alloc_enefunc(enefunc, EneFuncAMBERScale, nvar)

    if (nvar > 0) then
      call read_data_real_wp_array &
           (file, 'enefunc:dihe_scnb', &
               (/nvar/), enefunc%dihe_scnb(0:nvar-1))
      call read_data_real_wp_array &
           (file, 'enefunc:dihe_scee', &
               (/nvar/), enefunc%dihe_scee(0:nvar-1))
    end if

    ! enefunc_str::other variables

    call read_data_integer &
         (file, 'enefunc:table:num_water', enefunc%table%num_water)
    call read_data_integer &
         (file, 'enefunc:table:atom_cls_no_O', enefunc%table%atom_cls_no_O)
    call read_data_integer &
         (file, 'enefunc:table:atom_cls_no_H', enefunc%table%atom_cls_no_H)
    constraints%num_water = enefunc%table%num_water
    if (constraints%water_type == TIP4) &
      call read_data_integer &
           (file, 'enefunc:table:atom_cls_no_D', enefunc%table%atom_cls_no_D)
    call read_data_real_wp &
         (file, 'enefunc:table:charge_O', enefunc%table%charge_O)
    call read_data_real_wp &
         (file, 'enefunc:table:charge_H', enefunc%table%charge_H)
    if (constraints%water_type == TIP4) &
      call read_data_real_wp &
           (file, 'enefunc:table:charge_D', enefunc%table%charge_D)
    call read_data_real_wp &
         (file, 'enefunc:table:mass_O', enefunc%table%mass_O)
    call read_data_real_wp &
         (file, 'enefunc:table:mass_H', enefunc%table%mass_H)
    if (constraints%water_type == TIP4) &
      call read_data_real_wp &
           (file, 'enefunc:table:mass_D', enefunc%table%mass_D)
    call read_data_byte_array &
         (file, 'enefunc:table:water_model', (/5/), enefunc%table%water_model)
    call read_data_real_wp &
         (file, 'enefunc:table:fudge_lj', enefunc%fudge_lj)
    call read_data_real_wp &
         (file, 'enefunc:table:fudge_qq', enefunc%fudge_qq)
    call read_data_integer &
         (file, 'enefunc:table:excl_level', enefunc%excl_level)


    ! read s_constraints information
    !
    call read_data_integer &
         (file, 'constraints:connect', constraints%connect)
    call read_data_real_wip &
         (file, 'constraints:water_rHH', constraints%water_rHH)
    call read_data_real_wip &
         (file, 'constraints:water_rOH', constraints%water_rOH)
    if (constraints%water_type == TIP4) &
      call read_data_real_wip &
           (file, 'constraints:water_rOD', constraints%water_rOD)
    call read_data_real_wip &
         (file, 'constraints:water_massO', constraints%water_massO)
    call read_data_real_wip &
         (file, 'constraints:water_massH', constraints%water_massH)

    ! domain_constraint_str::ConstraintsDomainBond
    call read_data_integer &
         (file, 'constraints:No_HGr_size', nvar)
    call read_data_integer &
         (file, 'constraints:HGr_local_size', nvar1)
    call read_data_integer &
         (file, 'constraints:HGr_bond_list_size', nvar2) ! == nvar1 + 1

    if (filenum == 1) &
      call alloc_constraints(constraints, ConstraintsDomainBond_pio, nvar, &
                             nvar1, file_tot_num)

    if (nvar > 0 .and. nvar1 > 0) then

      call read_data_integer_array &
           (file, 'constraints:No_HGr', &
            (/nvar/),       constraints%No_HGr_pio   (         1:nvar, filenum))
      call read_data_integer_array &
           (file, 'constraints:HGr_local', &
            (/nvar1,nvar/), constraints%HGr_local_pio(1:nvar1, 1:nvar, filenum))

    end if

    do i = 1, nvar
      do j = 1, nvar1

        nvar3 = constraints%HGr_local_pio(j,i, filenum)
        if (nvar2 == 0 .or. nvar3 == 0) &
          cycle

        call read_data_integer_array             &
             (file, 'constraints:HGr_bond_list', &
              (/nvar2,nvar3,1,1/),               &
              constraints%HGr_bond_list_pio(1:nvar2, 1:nvar3, j, i, filenum))
        call read_data_real_wip_array            &
             (file, 'constraints:HGr_bond_dist', &
              (/nvar2,nvar3,1,1/),               &
               constraints%HGr_bond_dist_pio(1:nvar2, 1:nvar3, j, i, filenum))

      end do
    end do

    ! dynvars
    !

    call read_data_logical &
         (file, 'dynvar:read_dynvar', read_file)

    if (present(dynvars) .and. read_file) then
      call read_data_real_wip &
         (file, 'dynvars:thermostat_momentum', dynvars%thermostat_momentum)
      call read_data_real_wip_array &
         (file, 'dynvars:barostat_momentum', &
             (/3/), dynvars%barostat_momentum(1:3))
    end if

    ! dynamics
    !
 
    call read_data_logical &
         (file, 'dynamics:read_iseed', read_file)

    if (present(dynamics) .and. read_file) then
      call read_data_integer &
           (file, 'dynamics:iseed', dynamics%iseed)
      iseed_tmp = dynamics%iseed
      call read_data_logical &
           (file, 'pio_restart', pio_restart)
    end if

    ! random internal states
    !
    call get_data_size &
         (file, 'random', nbyte, .false.)
    allocate(bytes(nbyte))

    call read_data_byte_array &
         (file, 'random', (/nbyte/), bytes(1:nbyte))

    call random_init(iseed_tmp)
    call random_stock_frombyte(bytes)
    call random_pull_stock
    
    deallocate(bytes)

    call read_data_logical &
         (file, 'remd:read_remd', read_file)
    if (read_file .and. present(remd)) then
      ndim = 0
      nrep = 0
      call read_data_integer &
           (file, 'remd:iseed', remd%iseed)
      call read_data_integer &
           (file, 'remd:total_nreplicas', remd%total_nreplicas)
      nrep = remd%total_nreplicas
      call read_data_integer &
           (file, 'remd:dimension', remd%dimension)
      ndim = remd%dimension
      if (ndim > 0 .and. nrep > 0 .and. filenum == 1) &
        call alloc_remd(remd, RemdReplicas_rst, ndim, nrep)

      if (ndim > 0 .and. nrep > 0) then
        call read_data_integer_array                      &
             (file, 'remd:repid2parmsetid', (/nrep/),     &
              remd%repid2parmsetid(1:nrep)) 
        call read_data_integer_array                      &
             (file, 'remd:num_criteria', (/nrep,ndim,2/), &
              remd%num_criteria(1:nrep,1:ndim,1:2)) 
        call read_data_integer_array                      &
             (file, 'remd:num_exchanges', (/nrep,ndim,2/),&
              remd%num_exchanges(1:nrep,1:ndim,1:2))
      end if
    end if

    call close_data(file)

    return

  end subroutine pio_read_domain_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_domain_str
  !> @brief        read domain information
  !! @authors      NT
  !! @param[in]    filename : restart filename 
  !! @param[inout] boundary : boundary information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_domain_str(filename, boundary, domain, convert)
                             
    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    logical,                 intent(in)    :: convert

    ! local variables
    integer                  :: file
    integer                  :: i, ncell, nvar


    call open_data(filename, IOFileDataRead, file)


    ! read constant values
    !
    if (.not. convert) &
      call read_data_real_wp &
         (file, 'ELECOEF', ELECOEF)


    ! read boundary information
    !

    call read_data_real_wip &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call read_data_real_wip &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call read_data_real_wip &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call read_data_real_wp &
         (file, 'boundary:origin_x', boundary%origin_x)
    call read_data_real_wp &
         (file, 'boundary:origin_y', boundary%origin_y)
    call read_data_real_wp &
         (file, 'boundary:origin_z', boundary%origin_z)
    call read_data_real_wip &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call read_data_real_wip &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call read_data_real_wip &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call read_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain(1:3))
    call read_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call read_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call read_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)


    ! read s_domain information
    !

    call read_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call read_data_integer &
         (file, 'domain:num_cell_local', domain%num_cell_local)

    ! domain_str::Dynvar variables
    call read_data_integer &
         (file, 'domain:ncell', ncell)
    call alloc_domain(domain, DomainDynvar, ncell, 1, 1)
    call alloc_domain(domain, DomainDynvar_Atom, ncell, 1, 1)

    call read_data_integer_array &
         (file, 'domain:num_atom', &
             (/ncell/), domain%num_atom   (1:ncell))

    do i = 1, ncell

      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      call read_data_real_wip_array &
         (file, 'domain:coord', &
             (/3,nvar,1/), domain%coord   (1:3, 1:nvar, i))
      call read_data_real_wip_array &
         (file, 'domain:velocity', &
             (/3,nvar,1/), domain%velocity(1:3, 1:nvar, i))
      call read_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g  (1:nvar, i))

    end do

    call close_data(file)

    return

  end subroutine pio_read_domain_str

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    pio_boundary_check
  !> @brief        check the compatibility of boundary condition
  !! @authors      JJ
  !! @param[in]    boundary : boundary structure
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_boundary_check(boundary)

    ! formal arguments
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    integer                  :: i, j

    ! check if the number of cells are multple of 2*domain
    !
    if (mod(boundary%num_cells_x, boundary%num_domain(1)) /= 0)            &
      call error_msg('Pio_Boundary_Check> Cell number in x dimesion is '// &
                     'not multiple of domain_x')
    if (boundary%num_cells_x <= boundary%num_domain(1))                    &
      call error_msg('Pio_Boundary_Check> Cell number in x dimesion is '// &
                     'not greater than domain_x')
    if (mod(boundary%num_cells_y, boundary%num_domain(2)) /= 0)            &
      call error_msg('Pio_Boundary_Check> Cell number in y dimesion is '// &
                     'not multiple of domain_y')
    if (boundary%num_cells_y <= boundary%num_domain(2))                    &
      call error_msg('Pio_Boundary_Check> Cell number in y dimesion is '// &
                     'not greater than domain_y')
    if (mod(boundary%num_cells_z, boundary%num_domain(3)) /= 0)            &
      call error_msg('Pio_Boundary_Check> Cell number in z dimesion is '// &
                     'not multiple of domain_z')
    if (boundary%num_cells_z <= boundary%num_domain(3))                    &
      call error_msg('Pio_Boundary_Check> Cell number in z dimesion is '// &
                     'not greater than domain_z')

    ! check the domain number in rst and control input
    !
    if (boundary%num_pio_domain(1) /= 0 .or. &
        boundary%num_pio_domain(2) /= 0 .or. &
        boundary%num_pio_domain(3) /= 0) then
      if (boundary%num_domain_pio(1) /= boundary%num_pio_domain(1)) &
        call error_msg('Pio_Boundary_Check> domain_x in rst and '// &
                       'pio_domain_x in control input are different')
      if (boundary%num_domain_pio(2) /= boundary%num_pio_domain(2)) &
        call error_msg('Pio_Boundary_Check> domain_y in rst and '// &
                       'pio_domain_y in control input are different')
      if (boundary%num_domain_pio(3) /= boundary%num_pio_domain(3)) &
        call error_msg('Pio_Boundary_Check> domain_z in rst and '// &
                       'pio_domain_z in control input are different')
    end if

    return

  end subroutine pio_boundary_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_check_ranked_file
  !> @brief        check ranked file name or not
  !! @authors      NT
  !! @param[in]    filename : file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_check_ranked_file(filename)

    ! function
    logical                  :: pio_check_ranked_file

    ! formal arguments
    character(100),          intent(in) :: filename

    ! local variables
    integer                  :: ci1, ci2


    ! check filename
    !
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    pio_check_ranked_file = (ci1 /= 0 .and. ci2 /= 0 .and. ci1 < ci2)

    return

  end function pio_check_ranked_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_get_ranked_filename
  !> @brief        get ranked file name
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[out]   nplace   : number of places
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_get_ranked_filename(filename, my_local_rank, nplace)

    ! function
    character(100)           :: pio_get_ranked_filename

    ! formal arguments
    character(*),            intent(in) :: filename
    integer,                 intent(in) :: my_local_rank
    integer,       optional, intent(in) :: nplace

    ! local variables
    integer                  :: ci1, ci2, ip
    character(50)            :: fmt_str

    ! constants
    integer,       parameter :: Places (7) = &
         (/10, 100, 1000, 10000, 100000, 1000000, 10000000/)


    if (present(nplace)) then
      ip = nplace
    else
      do ip = 5, size(Places)
        if (nproc_world < Places(ip)) &
          exit
      end do
    end if

    ! check filename
    !
    pio_get_ranked_filename = filename
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    if (ci1 == 0 .or. ci2 ==0 .or. ci1 > ci2) then
      call error_msg( &
      'Pio_Get_Ranked_Filename> Filename is not correctly ranked in [OUTPUT]')
    end if

    write(fmt_str,'(A,I0,A,I0,A)') '(A,I',ip,'.',ip,',A)'

    write(pio_get_ranked_filename,fmt=fmt_str)   &
          trim(filename(1:ci1-1)), my_local_rank, &
          trim(filename(ci2+1:))

    return

  end function pio_get_ranked_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_open_file
  !> @brief        open domain decomposition restart file
  !! @authors      NT
  !! @param[inout] file     : file unit number
  !! @param[in]    filename : file name
  !! @param[in]    in_out   : input or output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_open_file(file, filename, in_out)

    ! formal arguments
    integer,                 intent(inout) :: file
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out

    ! local variables
    integer                  :: imark, ist


    select case(in_out)

    case (IOFileInput)

      call open_binary_file(file, filename, in_out, &
                            pio_get_file_endian(filename))
      read(file, iostat=ist) imark

      if (ist /= 0 .or. imark /= PioFileHeaderMarker) &
        call error_msg('Pio_Open_File> Unknown file format :'//&
        trim(filename))

    case (IOFileOutputNew, IOFileOutputReplace)

      call open_binary_file(file, filename, in_out)
      write(file) PioFileHeaderMarker

    case default

      call error_msg('Pio_Open_File> Unknown I/O mode')

    end select

    return

  end subroutine pio_open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_get_file_endian
  !> @brief      
  !! @authors      NT
  !! @param[in]  
  !! @param[out] 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_get_file_endian(filename)

    ! function
    integer          :: pio_get_file_endian

    ! formal arguments
    character(*),    intent(in) :: filename

    ! local variables
    integer          :: file, ival


    file = get_unit_no()

    open(file, &
         file    = filename,        &
         status  = 'old',           &
         form    = 'unformatted',   &
         convert = 'little_endian', &
         access  = 'direct',        &
         recl    = 4)
    read(file,rec=1) ival

    if (ival /= 4) then

      close(file)

      open(file, &
           file    = filename,      &
           status  = 'old',         &
           form    = 'unformatted', &
           convert = 'big_endian',  &
           access  = 'direct',      &
           recl    = 4)
      read(file,rec=1) ival

      if (ival /= 4) &
        call error_msg('Pio_Get_File_Endian> unknown file format.')

      pio_get_file_endian = IOFileBigEndian

    else

      pio_get_file_endian = IOFileLittleEndian

    end if

    close(file)
    call free_unit_no(file)

    return

  end function pio_get_file_endian

end module sp_parallel_io_mod
