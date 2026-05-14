!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_martini_mod
!> @brief   Calculates Martini nonlocal energy.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_martini_mod

  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public :: compute_energy_martini_vdw
  public :: compute_energy_martini_ele
  public :: compute_energy_martini_rf
  public :: compute_force_martini_vdw
  public :: compute_force_martini_ele
  public :: compute_force_martini_rf

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_martini_vdw
  !> @brief        Calculates van der Waals energy.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_martini_vdw(domain, enefunc, pairlist, &
                                        coord, force, evdw, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: Coef_A12(2), Coef_B12(2), Coef_C12
    real(wp)                  :: Coef_A6(2), Coef_B6(2), Coef_C6
    real(wp)                  :: Coef_Aelec, Coef_Belec, Coef_Celec
    real(wp)                  :: cutoff, switch
    integer                   :: i, ix, iy, j, k, ki, ij, L, num_atom
    integer                   :: id, num_nb15
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGexv, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cg_cutoffdist_vdw
    switch          =  enefunc%cg_switchdist_vdw
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_cg_exv_calc
    nb15_calc_list  => pairlist%cg_exv_list

    num_atom        =  domain%num_atom_domain + domain%num_atom_boundary
    cutoff2         =  cutoff * cutoff
    switch2         =  switch * switch

    ! coefficient for shift
    !
    diff = cutoff - switch
    Coef_A12(1) = 13.0_wp*switch - 16.0_wp*cutoff
    Coef_A12(1) = 12.0_wp*Coef_A12(1) / (cutoff**14) / (diff**2)
    Coef_A12(2) = 0.0_wp
    Coef_B12(1) = 13.0_wp*switch - 15.0_wp*cutoff
    Coef_B12(1) = -12.0_wp*Coef_B12(1) / (cutoff**14) / (diff**3)
    Coef_B12(2) = 0.0_wp
    Coef_C12    = 1.0_wp/(cutoff**12) - Coef_A12(1)/3.0_wp*(diff**3) &
                                      - Coef_B12(1)/4.0_wp*(diff**4)
    Coef_A6 (1) = 7.0_wp*switch - 10.0_wp*cutoff
    Coef_A6 (1) = 6.0_wp*Coef_A6(1) / (cutoff**8) / (diff**2)
    Coef_A6 (2) = 0.0_wp
    Coef_B6 (1) = 7.0_wp*switch - 9.0_wp*cutoff
    Coef_B6 (1) = -6.0_wp*Coef_B6(1) / (cutoff**8) / (diff**3)
    Coef_B6 (2) = 0.0_wp
    Coef_C6     = 1.0_wp/(cutoff**6) - Coef_A6(1)/3.0_wp*(diff**3) &
                                     - Coef_B6(1)/4.0_wp*(diff**4)
    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, j, k, rtmp, force_local, iatmcls, num_nb15,   &
    !$omp         evdw_temp, dij, rij2,  within_cutoff, rij, rij2_inv,     &
    !$omp         rij_inv, rij_diff, rij_diff2, rij_diff3, rij_diff4,      &
    !$omp         jatmcls, lj12, lj6, rij6_inv, rij12_inv, L, term_lj12,   &
    !$omp         term_lj6, grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ix = id+1, num_atom, nthread

      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp
      iatmcls  = atmcls(ix)
      evdw_temp = 0.0_wp
      num_nb15 = num_nb15_calc(ix)

      do k = 1, num_nb15

        j  = nb15_calc_list(k,ix)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2 > cutoff2)
        L = merge(1, 2, rij2 > switch2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = 1.0_wp / rij
        rij_diff = rij - switch
        rij_diff2 = rij_diff*rij_diff
        rij_diff3 = rij_diff2*rij_diff
        rij_diff4 = rij_diff3*rij_diff

        jatmcls = atmcls(j)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        rij6_inv  = rij2_inv*rij2_inv*rij2_inv
        rij12_inv = rij6_inv*rij6_inv
        term_lj12 = rij12_inv - Coef_A12(L)/3.0_wp*rij_diff3 &
                              - Coef_B12(L)/4.0_wp*rij_diff4 - Coef_C12
        term_lj6  = rij6_inv  - Coef_A6(L)/3.0_wp*rij_diff3 &
                              - Coef_B6(L)/4.0_wp*rij_diff4 - Coef_C6

        evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)*within_cutoff

        term_lj12 = -12.0_wp*rij12_inv*rij2_inv &
                    - (Coef_A12(L)*rij_diff2+Coef_B12(L)*rij_diff3)*rij_inv
        term_lj6  = - 6.0_wp*rij6_inv*rij2_inv  &
                    - (Coef_A6(L)*rij_diff2+Coef_B6(L)*rij_diff3)*rij_inv
        grad_coef = term_lj12*lj12 - term_lj6*lj6
        grad_coef = grad_coef * within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      evdw(id+1) = evdw(id+1) + evdw_temp

    end do

    !$omp end parallel

    call timer(TimerCGexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_martini_vdw


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_martini_ele
  !> @brief        Calculates electrostatic energy.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_martini_ele(domain, enefunc, pairlist, &
                                        coord, force, elec, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: elec(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff, cutoff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: Coef_A12(2), Coef_B12(2), Coef_C12
    real(wp)                  :: Coef_A6(2), Coef_B6(2), Coef_C6
    real(wp)                  :: Coef_Aelec, Coef_Belec, Coef_Celec
    integer                   :: i, ix, iy, j, k, ki, ij, num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: cg_elec_list(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: charge(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    start_atom      => domain%start_atom
    charge          => domain%charge

    cg_elec_list    => enefunc%cg_elec_list
    num_nb15_pbc    => pairlist%num_cg_ele_calc
    nb15_calc_list  => pairlist%cg_ele_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff          =  enefunc%cg_cutoffdist_ele
    cutoff2         =  cutoff*cutoff

    ! coefficient for shift
    !
    Coef_Aelec  = -5.0_wp/(cutoff**4)
    Coef_Belec  =  4.0_wp/(cutoff**5)
    Coef_Celec  =  1.0_wp/cutoff - Coef_Aelec/3.0_wp*(cutoff**3) &
                                 - Coef_Belec/4.0_wp*(cutoff**4)
    el_fact     = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, rtmp, qtmp, num_nb15, elec_temp, force_local, &
    !$omp         k, j, dij, rij2, within_cutoff, rij, rij2_inv, rij_inv,  &
    !$omp         jqtmp, term_elec, grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, enefunc%num_cg_elec, nthread

      ix = cg_elec_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      qtmp    = charge(ix)

      num_nb15 = num_nb15_pbc(i)
      elec_temp   = 0.0_wp
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      !ocl norecurrence
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = rtmp(1) - coord(j,1)
        dij(2) = rtmp(2) - coord(j,2)
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        within_cutoff = merge(0.0_wp, 1.0_wp, rij2>cutoff2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv  = 1.0_wp / rij

        jqtmp   = charge(j)

        term_elec = rij_inv   - Coef_Aelec/3.0_wp*rij*rij2 &
                              - Coef_Belec/4.0_wp*rij2*rij2 - Coef_Celec 
        term_elec = term_elec * el_fact

        elec_temp = elec_temp + qtmp*jqtmp*term_elec*within_cutoff

        term_elec = - rij_inv*rij2_inv - Coef_Aelec*rij - Coef_Belec*rij2
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec*within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      elec(id+1) = elec(id+1) + elec_temp

    end do

    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_martini_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_martini_rf
  !> @brief        Calculates electrostatic energy.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_martini_rf(domain, enefunc, pairlist, &
                                       coord, force, elec, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: elec(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff, cutoff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: elec_temp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: epsilon_rf, epsilon
    real(wp)                  :: krf, crf
    integer                   :: i, ix, iy, j, k, ki, ij, num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: cg_elec_list(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: charge(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    start_atom      => domain%start_atom
    charge          => domain%charge

    cg_elec_list    => enefunc%cg_elec_list
    num_nb15_pbc    => pairlist%num_cg_ele_calc
    nb15_calc_list  => pairlist%cg_ele_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff          =  enefunc%cg_cutoffdist_ele
    cutoff2         =  cutoff*cutoff

    ! coefficient for shift
    !
    epsilon_rf = enefunc%epsilon_rf
    epsilon    = enefunc%dielec_const
    if (epsilon_rf < EPS) then
      krf = 1.0_wp/(2.0_wp*cutoff2*cutoff)
    else
      krf = (epsilon_rf-epsilon) / (2.0_wp*epsilon_rf+epsilon) 
      krf = krf / (cutoff2*cutoff)
    end if  
    crf = 1.0_wp/cutoff + krf*cutoff2
    el_fact     = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, rtmp, qtmp, num_nb15, elec_temp, force_local, &
    !$omp         k, j, dij, rij2, within_cutoff, rij, rij2_inv, rij_inv,  &
    !$omp         jqtmp, term_elec, grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    viri(1) = 0.0_wp
    viri(2) = 0.0_wp
    viri(3) = 0.0_wp

    do i = id+1, enefunc%num_excl, nthread

      ix = enefunc%excl_list(1,i)
      j  = enefunc%excl_list(2,i)
      qtmp    = charge(ix)
      jqtmp   = charge(j)

      if (abs(qtmp) > EPS .and. abs(jqtmp) > EPS) then
        dij(1) = coord(ix,1) - coord(j,1)
        dij(2) = coord(ix,2) - coord(j,2)
        dij(3) = coord(ix,3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        term_elec = krf*rij2 - crf
        term_elec = term_elec * el_fact
        elec(id+1) = elec(id+1) + qtmp*jqtmp*term_elec
        term_elec = 2.0_wp*krf
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)
        force(ix,1,id+1) = force(ix,1,id+1) - work(1)
        force(ix,2,id+1) = force(ix,2,id+1) - work(2)
        force(ix,3,id+1) = force(ix,3,id+1) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)
      end if

    end do

    virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
    virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
    virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    do i = id+1, enefunc%num_cg_elec, nthread

      ix = cg_elec_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      qtmp    = charge(ix)

      num_nb15 = num_nb15_pbc(i)
      elec_temp   = 0.0_wp
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      !ocl norecurrence
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = rtmp(1) - coord(j,1)
        dij(2) = rtmp(2) - coord(j,2)
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        within_cutoff = merge(0.0_wp, 1.0_wp, rij2>cutoff2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv  = 1.0_wp / rij

        jqtmp   = charge(j)

        term_elec = rij_inv + krf*rij2 - crf
        term_elec = term_elec * el_fact

        elec_temp = elec_temp + qtmp*jqtmp*term_elec*within_cutoff

        term_elec = - rij_inv*rij2_inv + 2.0_wp*krf
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec*within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      elec(id+1) = elec(id+1) + elec_temp

    end do

    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_martini_rf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_martini_vdw
  !> @brief        Calculates van der Waals force.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_martini_vdw(domain, enefunc, pairlist, &
                                       coord, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: Coef_A12(2), Coef_B12(2), Coef_C12
    real(wp)                  :: Coef_A6(2), Coef_B6(2), Coef_C6
    real(wp)                  :: Coef_Aelec, Coef_Belec, Coef_Celec
    real(wp)                  :: cutoff, switch
    integer                   :: i, ix, iy, j, k, ki, ij, L, num_atom
    integer                   :: id, num_nb15
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGexv, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cg_cutoffdist_vdw
    switch          =  enefunc%cg_switchdist_vdw
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    num_nb15_calc   => pairlist%num_cg_exv_calc
    nb15_calc_list  => pairlist%cg_exv_list

    num_atom        =  domain%num_atom_domain + domain%num_atom_boundary
    cutoff2         =  cutoff * cutoff
    switch2         =  switch * switch

    ! coefficient for shift
    !
    diff = cutoff - switch
    Coef_A12(1) = 13.0_wp*switch - 16.0_wp*cutoff
    Coef_A12(1) = 12.0_wp*Coef_A12(1) / (cutoff**14) / (diff**2)
    Coef_A12(2) = 0.0_wp
    Coef_B12(1) = 13.0_wp*switch - 15.0_wp*cutoff
    Coef_B12(1) = -12.0_wp*Coef_B12(1) / (cutoff**14) / (diff**3)
    Coef_B12(2) = 0.0_wp
    Coef_A6 (1) = 7.0_wp*switch - 10.0_wp*cutoff
    Coef_A6 (1) = 6.0_wp*Coef_A6(1) / (cutoff**8) / (diff**2)
    Coef_A6 (2) = 0.0_wp
    Coef_B6 (1) = 7.0_wp*switch - 9.0_wp*cutoff
    Coef_B6 (1) = -6.0_wp*Coef_B6(1) / (cutoff**8) / (diff**3)
    Coef_B6 (2) = 0.0_wp

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, j, k, rtmp, force_local, iatmcls, num_nb15,   &
    !$omp         dij, rij2,  within_cutoff, rij, rij2_inv, rij_inv,       &
    !$omp         rij_diff, rij_diff2, rij_diff3, rij_diff4, jatmcls,      &
    !$omp         lj12, lj6, rij6_inv, rij12_inv, L, term_lj12, term_lj6,  &
    !$omp         grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ix = id+1, num_atom, nthread

      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp
      iatmcls  = atmcls(ix)
      num_nb15 = num_nb15_calc(ix)

      do k = 1, num_nb15

        j  = nb15_calc_list(k,ix)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        within_cutoff = merge(0.0_wp,1.0_wp,rij2 > cutoff2)
        L = merge(1, 2, rij2 > switch2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv = 1.0_wp / rij
        rij_diff = rij - switch
        rij_diff2 = rij_diff*rij_diff
        rij_diff3 = rij_diff2*rij_diff

        jatmcls = atmcls(j)
        lj12    = nonb_lj12(jatmcls,iatmcls)
        lj6     = nonb_lj6 (jatmcls,iatmcls)

        rij6_inv  = rij2_inv*rij2_inv*rij2_inv
        rij12_inv = rij6_inv*rij6_inv

        term_lj12 = -12.0_wp*rij12_inv*rij2_inv &
                    - (Coef_A12(L)*rij_diff2+Coef_B12(L)*rij_diff3)*rij_inv
        term_lj6  = - 6.0_wp*rij6_inv*rij2_inv  &
                    - (Coef_A6(L)*rij_diff2+Coef_B6(L)*rij_diff3)*rij_inv
        grad_coef = term_lj12*lj12 - term_lj6*lj6
        grad_coef = grad_coef * within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    call timer(TimerCGexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_force_martini_vdw


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_martini_ele
  !> @brief        Calculates electrostatic force.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_martini_ele(domain, enefunc, pairlist, &
                                        coord, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff, cutoff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: Coef_A12(2), Coef_B12(2), Coef_C12
    real(wp)                  :: Coef_A6(2), Coef_B6(2), Coef_C6
    real(wp)                  :: Coef_Aelec, Coef_Belec, Coef_Celec
    integer                   :: i, ix, iy, j, k, ki, ij, num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: cg_elec_list(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: charge(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    start_atom      => domain%start_atom
    charge          => domain%charge

    cg_elec_list    => enefunc%cg_elec_list
    num_nb15_pbc    => pairlist%num_cg_ele_calc
    nb15_calc_list  => pairlist%cg_ele_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff          =  enefunc%cg_cutoffdist_ele
    cutoff2         =  cutoff*cutoff

    ! coefficient for shift
    !
    Coef_Aelec  = -5.0_wp/(cutoff**4)
    Coef_Belec  =  4.0_wp/(cutoff**5)
    el_fact     = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, rtmp, qtmp, num_nb15, force_local, k, j, dij, &
    !$omp         rij2, within_cutoff, rij, rij2_inv, rij_inv, jqtmp,      &
    !$omp         term_elec, grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, enefunc%num_cg_elec, nthread

      ix = cg_elec_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      qtmp    = charge(ix)

      num_nb15 = num_nb15_pbc(i)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      !ocl norecurrence
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = rtmp(1) - coord(j,1)
        dij(2) = rtmp(2) - coord(j,2)
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        within_cutoff = merge(0.0_wp, 1.0_wp, rij2>cutoff2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv  = 1.0_wp / rij

        jqtmp   = charge(j)

        term_elec = - rij_inv*rij2_inv - Coef_Aelec*rij - Coef_Belec*rij2
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec*within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_force_martini_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_martini_rf
  !> @brief        Calculates electrostatic force.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_martini_rf(domain, enefunc, pairlist, &
                                      coord, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, rij2_inv, rij_inv
    real(wp)                  :: rij6_inv, rij12_inv
    real(wp)                  :: lj12, lj6, el_fact
    real(wp)                  :: diff, cutoff
    real(wp)                  :: rij_diff, rij_diff2, rij_diff3, rij_diff4
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, switch2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: viri(3), within_cutoff
    real(wp)                  :: epsilon_rf, epsilon
    real(wp)                  :: krf, crf
    integer                   :: i, ix, iy, j, k, ki, ij, num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: ncell, ncell_local

    integer,          pointer :: cg_elec_list(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: charge(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    start_atom      => domain%start_atom
    charge          => domain%charge

    cg_elec_list    => enefunc%cg_elec_list
    num_nb15_pbc    => pairlist%num_cg_ele_calc
    nb15_calc_list  => pairlist%cg_ele_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff          =  enefunc%cg_cutoffdist_ele
    cutoff2         =  cutoff*cutoff

    ! coefficient for shift
    !
    epsilon_rf = enefunc%epsilon_rf
    epsilon    = enefunc%dielec_const
    if (epsilon_rf < EPS) then
      krf = 1.0_wp/(2.0_wp*cutoff2*cutoff)
    else
      krf = (epsilon_rf-epsilon) / (2.0_wp*epsilon_rf+epsilon) 
      krf = krf / (cutoff2*cutoff)
    end if  
    crf = 1.0_wp/cutoff + krf*cutoff2
    el_fact     = ELECOEF / enefunc%dielec_const

    !$omp parallel default(shared)                                         &
    !$omp private(id, i, ix, rtmp, qtmp, num_nb15, force_local, k, j, dij, &
    !$omp         rij2, within_cutoff, rij, rij2_inv, rij_inv, jqtmp,      &
    !$omp         term_elec, grad_coef, work, viri)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    viri(1) = 0.0_wp
    viri(2) = 0.0_wp
    viri(3) = 0.0_wp

    do i = id+1, enefunc%num_excl, nthread

      ix = enefunc%excl_list(1,i)
      j  = enefunc%excl_list(2,i)
      qtmp    = charge(ix)
      jqtmp   = charge(j)

      if (abs(qtmp) > EPS .and. abs(jqtmp) > EPS) then
        dij(1) = coord(ix,1) - coord(j,1)
        dij(2) = coord(ix,2) - coord(j,2)
        dij(3) = coord(ix,3) - coord(j,3)
        term_elec = 2.0_wp*krf
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)
        force(ix,1,id+1) = force(ix,1,id+1) - work(1)
        force(ix,2,id+1) = force(ix,2,id+1) - work(2)
        force(ix,3,id+1) = force(ix,3,id+1) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)
      end if

    end do

    virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
    virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
    virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    do i = id+1, enefunc%num_cg_elec, nthread

      ix = cg_elec_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)
      qtmp    = charge(ix)

      num_nb15 = num_nb15_pbc(i)
      force_local(1) = 0.0_wp
      force_local(2) = 0.0_wp
      force_local(3) = 0.0_wp
      viri(1) = 0.0_wp
      viri(2) = 0.0_wp
      viri(3) = 0.0_wp

      !ocl norecurrence
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = rtmp(1) - coord(j,1)
        dij(2) = rtmp(2) - coord(j,2)
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        within_cutoff = merge(0.0_wp, 1.0_wp, rij2>cutoff2)
        rij = sqrt(rij2)
        rij2_inv = 1.0_wp / rij2
        rij_inv  = 1.0_wp / rij

        jqtmp   = charge(j)

        term_elec = - rij_inv*rij2_inv + 2.0_wp*krf
        term_elec = term_elec * el_fact
        grad_coef = qtmp*jqtmp*term_elec*within_cutoff
        work(1) = grad_coef*dij(1)
        work(2) = grad_coef*dij(2)
        work(3) = grad_coef*dij(3)
        viri(1) = viri(1) + dij(1)*work(1)
        viri(2) = viri(2) + dij(2)*work(2)
        viri(3) = viri(3) + dij(3)*work(3)

        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)

    end do

    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_force_martini_rf

end module cg_energy_martini_mod

