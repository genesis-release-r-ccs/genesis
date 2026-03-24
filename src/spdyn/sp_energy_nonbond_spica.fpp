!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_spica_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_spica_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod
#ifdef OMP
  use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_notbl_spica
  public  :: compute_energy_nonbond_notbl_pspica
  public  :: compute_energy_nonbond_notbl_spica_check
  public  :: compute_force_nonbond_notbl_spica_check
  public  :: compute_energy_nonbond_notbl_pspica_check
  public  :: compute_force_nonbond_notbl_pspica_check
  public  :: compute_force_nonbond_notbl_spica
  public  :: compute_force_nonbond_notbl_pspica

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_spica_check
  !> @brief        calculate nonbonded energy with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_spica_check( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial,  &
                                                eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB
    real(wp)                  :: R, minimum_contact
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:), charge_list(:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: num_nb15_calc_c(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)
    integer,          pointer :: nb15_calc_list_c(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge
    charge_list       => domain%charge_list

    table_ene         => enefunc%table%table_ene
    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff
    minimum_contact   =  enefunc%minimum_contact

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         elec_temp, evdw_temp, dij, viri, within_cutoff,          &
    !$omp         rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, domain%num_charge_domain, nthread

      kki = charge_list(i)
      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_c(i)

        ij = nb15_calc_list_c(k,i)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_spica_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_spica_check
  !> @brief        calculate nonbonded force with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_spica_check( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB
    real(wp)                  :: R, minimum_contact
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: atmcls(:,:), charge_list(:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: num_nb15_calc_c(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)
    integer,          pointer :: nb15_calc_list_c(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge
    charge_list       => domain%charge_list

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff
    minimum_contact   =  enefunc%minimum_contact

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         dij, viri, within_cutoff, rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, domain%num_charge_domain, nthread

      kki = charge_list(i)
      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_c(i)

        ij = nb15_calc_list_c(k,i)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_spica_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_pspica_check
  !> @brief        calculate nonbonded energy with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_pspica_check( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial,  &
                                                eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB
    real(wp)                  :: R, minimum_contact
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge

    table_ene         => enefunc%table%table_ene
    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff
    minimum_contact   =  enefunc%minimum_contact

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         elec_temp, evdw_temp, dij, viri, within_cutoff,          &
    !$omp         rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_pspica_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_pspica_check
  !> @brief        calculate nonbonded force with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_pspica_check( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB
    real(wp)                  :: R, minimum_contact
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff
    minimum_contact   =  enefunc%minimum_contact

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         dij, viri, within_cutoff, rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2   = max(rij2, minimum_contact)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_pspica_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_spica
  !> @brief        calculate nonbonded energy with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_spica( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial,  &
                                                eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB 
    real(wp)                  :: R
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:), charge_list(:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: num_nb15_calc_c(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)
    integer,          pointer :: nb15_calc_list_c(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge
    charge_list       => domain%charge_list

    table_ene         => enefunc%table%table_ene
    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water 
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         elec_temp, evdw_temp, dij, viri, within_cutoff,          &
    !$omp         rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, domain%num_charge_domain, nthread

      kki = charge_list(i)
      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_c(i)

        ij = nb15_calc_list_c(k,i)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp

    end do

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB 
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_spica 

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_notbl_pspica
  !> @brief        calculate nonbonded energy with lookup table for elec
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] eelec      : electrostatic energy of target systems
  !! @param[inout] evdw       : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_notbl_pspica( &
                                                domain, enefunc, pairlist, &
                                                atmcls_pbc,                &
                                                coord_pbc, force, virial,  &
                                                eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, rij2_inv, rij3_inv, rij_inv
    real(wp)                  :: ljA, ljB 
    real(wp)                  :: R
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    real(wip),        pointer :: coord(:,:,:)
    real(wp),         pointer :: charge(:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer :: num_nb15_calc_s(:)
    integer,          pointer :: num_nb15_calc_w(:)
    integer,          pointer :: num_nb15_calc_ww(:)
    integer,          pointer :: nb15_calc_list_w(:,:)
    integer,          pointer :: nb15_calc_list_ww(:,:)
    integer,          pointer :: nb15_calc_list_s(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge

    table_ene         => enefunc%table%table_ene
    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water 
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, k, ki, kki, iy, ij, j, rij2, L, L1, R, rtmp,  &
    !$omp         qtmp, jqtmp, iatmcls, jatmcls, term_ljA, term_ljB,       &
    !$omp         term_elec, grad_coef, work, force_local, ljA, ljB,       &
    !$omp         elec_temp, evdw_temp, dij, viri, within_cutoff,          &
    !$omp         rij2_inv, rij3_inv, rij_inv)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      iatmcls = atmcls_pbc(kki)

      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      elec_temp      = 0.0_wp
      evdw_temp      = 0.0_wp
      viri(1)        = 0.0_wp
      viri(2)        = 0.0_wp
      viri(3)        = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 4.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij2_inv * rij2_inv * within_cutoff
        term_ljA  = term_ljB * term_ljB * term_ljB
        term_ljB  = term_ljB * rij_inv
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -12.0_wp * term_ljA * rij2_inv
        term_ljB  = - 5.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        dij(1) = rtmp(1) - coord_pbc(1,ij,1)
        dij(2) = rtmp(2) - coord_pbc(2,ij,1)
        dij(3) = rtmp(3) - coord_pbc(3,ij,1)
        rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)
        jqtmp   = coord_pbc(4,ij,1)

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
        term_elec = term_elec * within_cutoff
        evdw_temp = evdw_temp + term_ljA*ljA - term_ljB*ljB
        elec_temp = elec_temp + qtmp*jqtmp*term_elec
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eelec(id+1) = eelec(id+1) + elec_temp
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_notbl_pspica 

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_spica
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_spica( &
                                                domain, enefunc, pairlist,   &
                                                npt, atmcls_pbc,             &
                                                coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)


    if (npt) then
      call compute_force_nonbond_notbl_spica_npt( &  
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force, virial)

    else
      call compute_force_nonbond_notbl_spica_nonpt( & 
           domain, enefunc, pairlist, atmcls_pbc,      &
           coord_pbc, force)
    end if

    return

  end subroutine compute_force_nonbond_notbl_spica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_pspica
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_pspica( &
                                                domain, enefunc, pairlist,   &
                                                npt, atmcls_pbc,             &
                                                coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)


    if (npt) then
      call compute_force_nonbond_notbl_pspica_npt( &
           domain, enefunc, pairlist, atmcls_pbc,  &
           coord_pbc, force, virial)

    else
      call compute_force_nonbond_notbl_pspica_nonpt( &
           domain, enefunc, pairlist, atmcls_pbc,    &
           coord_pbc, force)
    end if

    return

  end subroutine compute_force_nonbond_notbl_pspica

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_spica_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_spica_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

#ifdef PKTIMER
    use Ctim
#endif
    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, rij2_inv
    real(wp)                  :: rij3_inv, rij_inv
    real(wp)                  :: rij3_invd
    real(wp)                  :: R, ljA, ljB
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    
    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer,contiguous :: charge_list(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_ww(:)
    integer,          pointer,contiguous :: num_nb15_calc_c(:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_ww(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_c(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge
    charge_list       => domain%charge_list

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge 

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1, R,    &
    !$omp         term_ljA, term_ljB, grad_coef, work, term_elec,           &
    !$omp         force_local, ljA, ljB, iatmcls, jatmcls, dij1, dij2,      &
    !$omp         dij3, jqtmp, viri, within_cutoff, rij2_inv, rij3_inv,     &
    !$omp         rij_inv)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, domain%num_charge_domain, nthread

      kki = charge_list(i)
      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_c(i)

        ij = nb15_calc_list_c(k,i)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      iatmcls = atmcls_pbc(kki)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 4.0_wp * term_ljB * rij2_inv 

        grad_coef = work(1)*ljA - work(2)*ljB 
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        term_ljB = term_ljB * rij_inv 
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 5.0_wp * term_ljB * rij2_inv

        grad_coef = work(1)*ljA - work(2)*ljB
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

#if defined (FUGAKU)
      !ocl norecurrence
#endif
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1 = rtmp(1) - coord_pbc(1,ij,1)
        dij2 = rtmp(2) - coord_pbc(2,ij,1)
        dij3 = rtmp(3) - coord_pbc(3,ij,1)

        rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB 
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_spica_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_pspica_npt
  !> @brief        calculate nonbonded force with virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  !! @param[inout] virial     : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_pspica_npt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force, virial)

#ifdef PKTIMER
    use Ctim
#endif
    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, rij2_inv
    real(wp)                  :: rij3_inv, rij_inv
    real(wp)                  :: rij3_invd
    real(wp)                  :: R, ljA, ljB
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), viri(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    
    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_ww(:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_ww(:,:)


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1, R,    &
    !$omp         term_ljA, term_ljB, grad_coef, work, term_elec,           &
    !$omp         force_local, ljA, ljB, iatmcls, jatmcls, dij1, dij2,      &
    !$omp         dij3, jqtmp, viri, within_cutoff, rij2_inv, rij3_inv,     &
    !$omp         rij_inv)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      iatmcls = atmcls_pbc(kki)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 4.0_wp * term_ljB * rij2_inv 
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*ljA - work(2)*ljB + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        term_ljB = term_ljB * rij_inv 
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 5.0_wp * term_ljB * rij2_inv
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*ljA - work(2)*ljB + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

      end do

#if defined (FUGAKU)
      !ocl norecurrence
#else
      !$omp simd
#endif
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1 = rtmp(1) - coord_pbc(1,ij,1)
        dij2 = rtmp(2) - coord_pbc(2,ij,1)
        dij3 = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        viri(1) = viri(1) + dij1*work(1)
        viri(2) = viri(2) + dij2*work(2)
        viri(3) = viri(3) + dij3*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_pspica_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_spica_nonpt
  !> @brief        calculate nonbonded force without virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_spica_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

#ifdef PKTIMER
    use Ctim
#endif
    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, rij2_inv
    real(wp)                  :: rij3_inv, rij_inv
    real(wp)                  :: R, ljA, ljB
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    
    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer,contiguous :: charge_list(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_ww(:)
    integer,          pointer,contiguous :: num_nb15_calc_c(:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_ww(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_c(:,:)

    real(dp)                  :: sas,eae


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge
    charge_list       => domain%charge_list

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    num_nb15_calc_c   => pairlist%num_nb15_calc_charge
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p
    nb15_calc_list_c  => pairlist%nb15_calc_list_charge 

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1, R,    &
    !$omp         term_ljA, term_ljB, grad_coef, work, term_elec,           &
    !$omp         force_local, ljA, ljB, iatmcls, jatmcls, dij1, dij2,      &
    !$omp         dij3, jqtmp, within_cutoff, rij2_inv, rij3_inv, rij_inv)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, domain%num_charge_domain, nthread

      kki = charge_list(i)
      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_c(i)

        ij = nb15_calc_list_c(k,i)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)

    end do

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      iatmcls = atmcls_pbc(kki)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 4.0_wp * term_ljB * rij2_inv 

        grad_coef = work(1)*ljA - work(2)*ljB 
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        term_ljB = term_ljB * rij_inv
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 5.0_wp * term_ljB * rij2_inv

        grad_coef = work(1)*ljA - work(2)*ljB 
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1 = rtmp(1) - coord_pbc(1,ij,1)
        dij2 = rtmp(2) - coord_pbc(2,ij,1)
        dij3 = rtmp(3) - coord_pbc(3,ij,1)

        rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        grad_coef = term_ljA*ljA - term_ljB*ljB
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_spica_nonpt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_notbl_pspica_nonpt
  !> @brief        calculate nonbonded force without virial
  !! @authors      JJ 
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions
  !! @param[in]    pairlist   : interaction list in each domain
  !! @param[inout] atmcls_pbc : atom class number
  !! @param[inout] coord_pbc  : coordinates for each cell
  !! @param[inout] force      : forces for each cell
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_notbl_pspica_nonpt( &
                                                 domain, enefunc, pairlist, &
                                                 atmcls_pbc,                &
                                                 coord_pbc, force)

#ifdef PKTIMER
    use Ctim
#endif
    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij1, dij2, dij3, rij2, rij2_inv
    real(wp)                  :: rij3_inv, rij_inv
    real(wp)                  :: R, ljA, ljB
    real(wp)                  :: term_ljA, term_ljB, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3), within_cutoff
    integer                   :: i, ix, iy, j, k, ki, kki, ij, L, L1
    integer                   :: id
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    
    real(wip),        pointer,contiguous :: coord(:,:,:)
    real(wp),         pointer,contiguous :: charge(:,:)
    real(wp),         pointer,contiguous :: trans1(:,:,:)
    real(wp),         pointer,contiguous :: nonb_ljA(:,:), nonb_ljB(:,:)
    real(wp),         pointer            :: density, cutoff
    real(wp),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: atmcls(:,:)
    integer,          pointer,contiguous :: natom(:), nsolute(:), start_atom(:)
    integer,          pointer,contiguous :: num_nb15_calc_s(:)
    integer,          pointer,contiguous :: num_nb15_calc_w(:)
    integer,          pointer,contiguous :: num_nb15_calc_ww(:)
    integer,          pointer,contiguous :: nb15_calc_list_s(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_w(:,:)
    integer,          pointer,contiguous :: nb15_calc_list_ww(:,:)

    real(dp)                  :: sas,eae


    atmcls            => domain%atom_cls_no
    natom             => domain%num_atom
    nsolute           => domain%num_solute
    start_atom        => domain%start_atom
    coord             => domain%coord
    trans1            => domain%trans_vec
    charge            => domain%charge

    table_grad        => enefunc%table%table_grad
    density           => enefunc%table%density
    cutoff            => enefunc%cutoffdist
    nonb_ljA          => enefunc%nonb_lj12
    nonb_ljB          => enefunc%nonb_lj6

    num_nb15_calc_s   => pairlist%num_nb15_calc_solute
    num_nb15_calc_w   => pairlist%num_nb15_calc_water
    num_nb15_calc_ww  => pairlist%num_nb15_calc_water_p
    nb15_calc_list_s  => pairlist%nb15_calc_list_solute
    nb15_calc_list_w  => pairlist%nb15_calc_list_water
    nb15_calc_list_ww => pairlist%nb15_calc_list_water_p

    ncell             =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local       =  domain%num_cell_local
    cutoff2           =  cutoff * cutoff

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, ki, kki, ij, rij2, L, L1, R,    &
    !$omp         term_ljA, term_ljB, grad_coef, work, term_elec,           &
    !$omp         force_local, ljA, ljB, iatmcls, jatmcls, dij1, dij2,      &
    !$omp         dij3, jqtmp, within_cutoff, rij2_inv, rij3_inv, rij_inv)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      ki = start_atom(i)
      do ix = 1, natom(i)
        kki = ki + ix
        coord_pbc(1,kki,1) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,kki,1) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,kki,1) = coord(3,ix,i) + trans1(3,ix,i)
        coord_pbc(4,kki,1) = charge(ix,i)
        atmcls_pbc(kki)    = atmcls(ix,i)
      end do
    end do

    !$omp barrier

    do kki = id+1, domain%num_atom_domain, nthread

      rtmp(1) = coord_pbc(1,kki,1)
      rtmp(2) = coord_pbc(2,kki,1)
      rtmp(3) = coord_pbc(3,kki,1)
      qtmp    = coord_pbc(4,kki,1)
      iatmcls = atmcls_pbc(kki)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp

!ocl norecurrence
      do k = 1, num_nb15_calc_w(kki)

        ij = nb15_calc_list_w(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 4.0_wp * term_ljB * rij2_inv 
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*ljA - work(2)*ljB + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_ww(kki)

        ij = nb15_calc_list_ww(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1    = rtmp(1) - coord_pbc(1,ij,1)
        dij2    = rtmp(2) - coord_pbc(2,ij,1)
        dij3    = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2 = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = sqrt(rij2_inv)
        rij2   = cutoff2 * density * rij2_inv
        L      = int(rij2)
        R      = rij2 - L

        term_ljB = rij2_inv * rij2_inv * within_cutoff
        term_ljA = term_ljB * term_ljB * term_ljB
        term_ljB = term_ljB * rij_inv
        work(1)  = -12.0_wp * term_ljA * rij2_inv
        work(2)  = - 5.0_wp * term_ljB * rij2_inv
        work(3)  = table_grad(L)
        work(3)  = work(3) + R*(table_grad(L+1)-work(3))
        work(3)  = work(3) * within_cutoff

        grad_coef = work(1)*ljA - work(2)*ljB + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

!ocl norecurrence
      do k = 1, num_nb15_calc_s(kki)

        ij = nb15_calc_list_s(k,kki)

        jatmcls = atmcls_pbc(ij)
        ljA     = nonb_ljA(jatmcls,iatmcls)
        ljB     = nonb_ljB(jatmcls,iatmcls)

        dij1 = rtmp(1) - coord_pbc(1,ij,1)
        dij2 = rtmp(2) - coord_pbc(2,ij,1)
        dij3 = rtmp(3) - coord_pbc(3,ij,1)

        jqtmp   = coord_pbc(4,ij,1)

        rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
        within_cutoff = merge(0.0_wp,1.0_wp,rij2>cutoff2)
        rij2_inv = 1.0_wp / rij2
        rij3_inv = rij2_inv * sqrt(rij2_inv)
        rij2    = cutoff2*density*rij2_inv

        L  = int(rij2)
        R  = rij2 - L

        term_ljB  = rij3_inv * rij3_inv * within_cutoff
        term_ljA  = term_ljB * rij3_inv
        term_ljA  = -9.0_wp * term_ljA * rij2_inv
        term_ljB  = -6.0_wp * term_ljB * rij2_inv
        term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
        term_elec = term_elec * within_cutoff
        grad_coef = term_ljA*ljA - term_ljB*ljB + qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(1,ij,1,id+1) = force(1,ij,1,id+1) + work(1)
        force(2,ij,1,id+1) = force(2,ij,1,id+1) + work(2)
        force(3,ij,1,id+1) = force(3,ij,1,id+1) + work(3)

      end do

      force(1,kki,1,id+1) = force(1,kki,1,id+1) + force_local(1)
      force(2,kki,1,id+1) = force(2,kki,1,id+1) + force_local(2)
      force(3,kki,1,id+1) = force(3,kki,1,id+1) + force_local(3)

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_notbl_pspica_nonpt

end module sp_energy_nonbond_spica_mod
