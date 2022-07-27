!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_gpu_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_gpu_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_tbl_lnr_gpu
  public  :: compute_energy_nonbond_tbl_ljpme_gpu
  public  :: compute_force_nonbond_tbl_lnr_gpu
  public  :: compute_force_nonbond_tbl_ljpme_gpu
  public  :: compute_energy_nonbond_notbl_gpu
  public  :: compute_force_nonbond_notbl_gpu

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_lnr_gpu
  !> @brief        calculate nonbonded energy with gpu
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_lnr_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_tbl_lnr_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_tbl_ljpme_gpu
  !> @brief        calculate nonbonded energy with gpu (LJ-PME)
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_tbl_ljpme_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_ljpme_univ(  &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_tbl_ljpme_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_gpu
  !> @brief        calculate nonbonded energy with gpu
  !! @authors      JJ, NT
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_gpu( &
                                                 domain, enefunc, pairlist, &
                                                 npt, coord_pbc, force,     &
                                                 virial, eelec, evdw,       &
                                                 ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j, ix, ixx, start_i
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,ixx,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_notable_univ(      &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         cell_move, natom, start_atom,                        &
         nonb_lj12, nonb_lj6, nonb_lj6_factor, table_ene,     &
         table_grad, univ_cell_pairlist1, univ_mask2,         &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         domain%num_atom_domain, MaxAtom, MaxAtomCls,         &
         num_atom_cls, ncell_local, ncell_bound, ncell_max,   &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2, system_size(1), system_size(2), system_size(3) )

#endif

    return

  end subroutine compute_energy_nonbond_notbl_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_lnr_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_lnr_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,              &            
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &            
         ncell_local, ncell_bound, ncell_max, cutoff_int,    &           
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_tbl_lnr_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_tbl_ljpme_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_tbl_ljpme_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_table_ljpme_univ(  &
         coord_pbc, force(:,:,:,1), ene_virial,              &  
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &  
         ncell_local, ncell_bound, ncell_max, cutoff_int,    & 
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_tbl_ljpme_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : logical to check virial calculation is required
  !! @param[in]    cpu_calc   : flag for cpu calculation or not
  !! @param[inout] coord_pbc  : pbc oriented coordinates
  !! @param[inout] force_omp  : temporary forces of target system
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_gpu( &
                                                domain, enefunc, pairlist, &
                                                npt, cpu_calc, coord_pbc,  &
                                                force_omp, force, virial,  &
                                                ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

#ifdef USE_GPU

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: nonb_lj6_factor(:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer(1),       pointer :: cell_move(:,:,:)
    integer,          pointer :: atmcls(:,:)
    integer(int2),    pointer :: cell_pairlist(:,:)
    integer,          pointer :: natom(:), start_atom(:)
    integer(int1),    pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx
    integer                   :: ij, start_i, ix, i


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    start_atom      => domain%start_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    nonb_lj6_factor => enefunc%nonb_lj6_factor
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ij = domain%num_atom_domain
    !$omp parallel do private(i,ix,start_i)
    do i = 1, ncell
      start_i = start_atom(i)
      do ix = 1, natom(i)
        coord_pbc(     start_i+ix,1,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(  ij+start_i+ix,1,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(2*ij+start_i+ix,1,1) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do
    !$omp end parallel do

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    univ_gpu_start = ncell_local
    call gpu_launch_compute_force_nonbond_notable_univ(      &
         coord_pbc, force(:,:,:,1), ene_virial,              &  
         cell_move, natom, start_atom, nonb_lj12, nonb_lj6,  &
         nonb_lj6_factor, table_grad, univ_cell_pairlist1,   &
         univ_mask2, univ_ix_natom, univ_ix_list,            &
         univ_iy_natom, univ_iy_list, univ_ij_sort_list,     &
         virial_check, domain%num_atom_domain, MaxAtom,      &
         MaxAtomCls, num_atom_cls,                           &  
         ncell_local, ncell_bound, ncell_max, cutoff_int,    & 
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

#endif

    return

  end subroutine compute_force_nonbond_notbl_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cpu_compute_force_tbl_univ    
  !> @brief        calculate nonbonded force within each cell
  !! @authors      JJ 
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cpu_compute_force_tbl_univ(coord, trans1, cell_move,          &
             system_size, charge, atmcls, natom, nonb_lj12, nonb_lj6,      &
             table_grad, univ_cell_pairlist1, univ_mask2,                  &
             univ_ix_natom, univ_ix_list, univ_iy_natom, univ_iy_list,     &
             univ_ij_sort_list, ncell_local, ncell_bound,                  &
             univ_ncell_nonzero, univ_update, density, cutoff2,            &
             coord_pbc, force)

    ! formal arguments
    real(wip),                intent(in)    :: coord(:,:,:)
    real(wp),                 intent(in)    :: trans1(:,:,:)
    integer(1),               intent(in)    :: cell_move(:,:,:)
    real(wp),                 intent(in)    :: system_size(:)
    real(wp),                 intent(in)    :: charge(:,:)
    integer,                  intent(in)    :: atmcls(:,:)
    integer,                  intent(in)    :: natom(:)
    real(wp),                 intent(in)    :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),                 intent(in)    :: table_grad(:)
    integer,                  intent(in)    :: univ_cell_pairlist1(:,:)
    integer(1),               intent(in)    :: univ_mask2(:,:)
    integer,                  intent(in)    :: univ_ix_natom(:)
    integer(1),               intent(in)    :: univ_ix_list(:,:)
    integer,                  intent(in)    :: univ_iy_natom(:)
    integer(1),               intent(in)    :: univ_iy_list(:,:)
    integer,                  intent(in)    :: univ_ij_sort_list(:)
    integer,                  intent(in)    :: ncell_local, ncell_bound
    integer,                  intent(in)    :: univ_ncell_nonzero
    integer,                  intent(in)    :: univ_update 
    real(wp),                 intent(in)    :: cutoff2, density
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij_list(1:3,MaxAtom), rij2_list(MaxAtom)
    real(wp)                  :: cell_move_local(1:3), dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: work(1:3), grad_coef
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), factor
    integer                   :: j_list(MaxAtom)
    integer                   :: index, univ_ij
    integer                   :: i, j, ix, iy, iix, iiy, k, L, L1, idx
    integer                   :: ix_natom, iy_natom
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: num_count


    !$omp parallel default(shared)                                            &
    !$omp private(id, i, j, ix, iy, iix, iiy, k, univ_ij, index, idx, j_list, &
    !$omp         ix_natom, iy_natom, num_count, L, L1, iatmcls, jatmcls,     &
    !$omp         cell_move_local, rtmp, qtmp, rij2, R, term_lj12, term_lj6,  &
    !$omp         term_elec, grad_coef, work, dij, inv_r2, lj12, lj6,         &
    !$omp         dij_list, rij2_list, force_local, factor)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell_local+ncell_bound, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(ix,1:3,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        force_local(1:3) = 0.0_wp

        do iy = ix + 1, natom(i)

          idx = iy + (ix-1)*univ_natom_max

          factor = real(univ_mask2(idx,i),wp)

          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (univ_mask2(idx,i) == 1 .and. rij2 <= cutoff2) then

          rij2    = cutoff2*density/rij2
          jatmcls = atmcls(iy,i)
          lj6     = nonb_lj6(iatmcls,jatmcls)
          lj12    = nonb_lj12(iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2

          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,i,id+1) = force(1,iy,i,id+1) + work(1)
          force(2,iy,i,id+1) = force(2,iy,i,id+1) + work(2)
          force(3,iy,i,id+1) = force(3,iy,i,id+1) + work(3)

          end if

        end do

        force(1,ix,i,id+1) = force(1,ix,i,id+1) + force_local(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) + force_local(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) + force_local(3)

      end do

    end do

    !$omp end parallel 

    return

  end subroutine cpu_compute_force_tbl_univ

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cpu_compute_force_notbl_univ    
  !> @brief        calculate nonbonded force within each cell
  !! @authors      JJ 
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cpu_compute_force_notbl_univ(coord, trans1, cell_move,        &
             system_size, charge, atmcls, natom, nonb_lj12, nonb_lj6,      &
             table_grad, univ_cell_pairlist1, univ_mask2,                  &
             univ_ix_natom, univ_ix_list, univ_iy_natom, univ_iy_list,     &
             univ_ij_sort_list, ncell_local, ncell_bound,                  &
             univ_ncell_nonzero, univ_update, density, cutoff2,            &
             coord_pbc, force)

    ! formal arguments
    real(wip),                intent(in)    :: coord(:,:,:)
    real(wp),                 intent(in)    :: trans1(:,:,:)
    integer(1),               intent(in)    :: cell_move(:,:,:)
    real(wp),                 intent(in)    :: system_size(:)
    real(wp),                 intent(in)    :: charge(:,:)
    integer,                  intent(in)    :: atmcls(:,:)
    integer,                  intent(in)    :: natom(:)
    real(wp),                 intent(in)    :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),                 intent(in)    :: table_grad(:)
    integer,                  intent(in)    :: univ_cell_pairlist1(:,:)
    integer(1),               intent(in)    :: univ_mask2(:,:)
    integer,                  intent(in)    :: univ_ix_natom(:)
    integer(1),               intent(in)    :: univ_ix_list(:,:)
    integer,                  intent(in)    :: univ_iy_natom(:)
    integer(1),               intent(in)    :: univ_iy_list(:,:)
    integer,                  intent(in)    :: univ_ij_sort_list(:)
    integer,                  intent(in)    :: ncell_local, ncell_bound
    integer,                  intent(in)    :: univ_ncell_nonzero
    integer,                  intent(in)    :: univ_update 
    real(wp),                 intent(in)    :: cutoff2, density
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij_list(1:3,MaxAtom), rij2_list(MaxAtom)
    real(wp)                  :: cell_move_local(1:3), dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: work(1:3), grad_coef
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), factor
    integer                   :: j_list(MaxAtom)
    integer                   :: index, univ_ij
    integer                   :: i, j, ix, iy, iix, iiy, k, L, L1, idx
    integer                   :: ix_natom, iy_natom
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: num_count


    !$omp parallel default(shared)                                            &
    !$omp private(id, i, j, ix, iy, iix, iiy, k, univ_ij, index, idx, j_list, &
    !$omp         ix_natom, iy_natom, num_count, L, L1, iatmcls, jatmcls,     &
    !$omp         cell_move_local, rtmp, qtmp, rij2, R, term_lj12, term_lj6,  &
    !$omp         term_elec, grad_coef, work, dij, inv_r2, lj12, lj6,         &
    !$omp         dij_list, rij2_list, force_local, factor)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell_local+ncell_bound, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(ix,1:3,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        force_local(1:3) = 0.0_wp

        do iy = ix + 1, natom(i)

          idx = iy + (ix-1)*univ_natom_max

          factor = real(univ_mask2(idx,i),wp)

          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (univ_mask2(idx,i) == 1 .and. rij2 <= cutoff2) then

          inv_r2 = 1.0_wp / rij2
          rij2    = cutoff2*density*inv_r2
          jatmcls = atmcls(iy,i)
          lj6     = nonb_lj6(iatmcls,jatmcls)
          lj12    = nonb_lj12(iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L

          term_lj6  = inv_r2 * inv_r2 * inv_r2
          term_lj12 = term_lj6 * term_lj6
          term_lj12 = -12.0_wp * term_lj12 * inv_r2
          term_lj6  = - 6.0_wp * term_lj6  * inv_r2
          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,i,id+1) = force(1,iy,i,id+1) + work(1)
          force(2,iy,i,id+1) = force(2,iy,i,id+1) + work(2)
          force(3,iy,i,id+1) = force(3,iy,i,id+1) + work(3)

          end if

        end do

        force(1,ix,i,id+1) = force(1,ix,i,id+1) + force_local(1)
        force(2,ix,i,id+1) = force(2,ix,i,id+1) + force_local(2)
        force(3,ix,i,id+1) = force(3,ix,i,id+1) + force_local(3)

      end do

    end do

    !$omp end parallel 

    return

  end subroutine cpu_compute_force_notbl_univ

end module sp_energy_nonbond_gpu_mod
