!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_nonlocal_mod
!> @brief   Calculates nonlocal energy.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_nonlocal_mod

  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
#endif

  implicit none
  private

  ! subroutines
  public :: compute_energy_PWMcos
  public :: compute_energy_PWMcosns
  public :: compute_energy_general_exv_AICG2P
  public :: compute_energy_DNA_exv
  public :: compute_energy_dna_base_pairing
  public :: compute_energy_CG_ele
  public :: compute_energy_CG_KH
  public :: compute_energy_CG_IDR_KH
  public :: compute_energy_CG_IDR_HPS
  public :: compute_force_CG_IDR_HPS
  public :: compute_energy_CG_IDR_MPIPI
  public :: compute_force_CG_IDR_MPIPI

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_general_exv_AICG2P
  !> @brief        Calculates excluded volume energy with pairlist (PBC).
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_general_exv_AICG2P(domain, enefunc, pairlist, &
                                               coord, force, eexv, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: eexv(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12, term_lj12
    real(wp)                  :: lj12
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: epsilon
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3), factor
    real(wp)                  :: exv_temp
    integer                   :: i, ix, iy, j, k, start_i
    integer                   :: num_atom, num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp),         pointer :: param_sigma(:,:), param_epsilon(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp), parameter :: exv_cutoff_factor = 4.0_wp
    real(wp), parameter :: exv_energy_shift = 1.0_wp / 4096.0_wp

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGexv, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cutoffdist
    param_epsilon   => enefunc%nonb_aicg_eps
    param_sigma     => enefunc%nonb_aicg_sig

    num_nb15_calc   => pairlist%num_cg_exv_calc
    nb15_calc_list  => pairlist%cg_exv_list

    num_atom        =  domain%num_atom_domain + domain%num_atom_boundary
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy, k, rtmp, num_nb15, exv_temp, &
    !$omp         dij, rij2, inv_rij2, inv_rij6, inv_rij12,      &
    !$omp         term_lj12, grad_coef, work, force_local,       &
    !$omp         iatmcls, jatmcls, lj12, epsilon, sigma,        &
    !$omp         sigma_sqr, sigma_6th, start_i, factor)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ix = id+1, num_atom, nthread

      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_calc(ix)
      iatmcls  = atmcls(ix)
      exv_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,ix)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        jatmcls = atmcls(j)
        epsilon = param_epsilon(iatmcls,jatmcls)
        sigma   = param_sigma(iatmcls,jatmcls)
        sigma_sqr = sigma * sigma

        if (rij2 > sigma_sqr * exv_cutoff_factor) cycle

        sigma_6th = sigma_sqr * sigma_sqr * sigma_sqr

        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij6  = sigma_6th * inv_rij6
        term_lj12 = inv_rij6 * inv_rij6

        exv_temp  = exv_temp + epsilon*(term_lj12 - exv_energy_shift)

        grad_coef   = - inv_rij2*12.0_wp*epsilon*term_lj12
        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      eexv(id+1) = eexv(id+1) + exv_temp

    end do

    !$omp end parallel

    call timer(TimerCGexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_general_exv_AICG2P

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_HPS
  !> @brief        Calculates IDR interaction in the HPS model (PBC).
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] ehps     : energy of IDR HPS interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_HPS(domain, enefunc, pairlist, &
                                       coord, force, ehps, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: ehps(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, term_lj6, term_lj12
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: lambda, sigma, epsilon
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: ehps_temp, ehps_ene
    integer                   :: i, ix, j, k
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp),         pointer :: param_sigma(:,:), param_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_cg_IDR_HPS, IDR_list(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRHPS, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cg_cutoffdist_126
    cutoff2         =  cutoff * cutoff
    epsilon         =  enefunc%cg_IDR_HPS_epsilon
    num_cg_IDR_HPS  => enefunc%num_cg_IDR_HPS
    IDR_list        => enefunc%cg_IDR_HPS_list
    param_lambda    => enefunc%cg_IDR_HPS_lambda
    param_sigma     => enefunc%cg_IDR_HPS_sigma

    num_nb15_pbc    => pairlist%num_cg_idr_hps_calc
    nb15_calc_list  => pairlist%cg_idr_hps_list

    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                              &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15, ehps_temp,   &
    !$omp         dij, rij2, inv_rij2, inv_rij6, term_lj6,      &
    !$omp         term_lj12, grad_coef, work, force_local,      &
    !$omp         iatmcls, jatmcls, lambda, sigma, ehps_ene)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_IDR_HPS, nthread

      ix = IDR_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_pbc(i)
      iatmcls  = atmcls(ix)
      ehps_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 > cutoff2) cycle

        jatmcls = atmcls(j)
        lambda  = param_lambda(iatmcls,jatmcls)
        sigma   = param_sigma (iatmcls,jatmcls)

        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        term_lj6  = sigma * inv_rij6
        term_lj12 = term_lj6 * term_lj6
        ehps_ene  = 4.0_wp*epsilon*(term_lj12 - term_lj6)

        if (term_lj6 >= 0.5_wp) then
          ehps_temp = ehps_temp + ehps_ene + (1.0_wp-lambda)*epsilon
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*epsilon*inv_rij2*grad_coef
        else
          ehps_temp = ehps_temp + lambda*ehps_ene
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*lambda*epsilon*inv_rij2*grad_coef
        end if

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      ehps(id+1) = ehps(id+1) + ehps_temp

    end do

    !$omp end parallel

    call timer(TimerCGIDRHPS, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_HPS

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_CG_IDR_HPS
  !> @brief        Calculates IDR interaction in the HPS model (PBC).
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] ehps     : energy of IDR HPS interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_CG_IDR_HPS(domain, enefunc, pairlist, &
                                      coord, force, ehps, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: ehps(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, term_lj6, term_lj12
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: lambda, sigma, epsilon
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: ehps_temp, ehps_ene
    real(wp)                  :: scale1, scale2, factor(2)
    integer                   :: i, ix, j, k
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp),         pointer :: param_sigma(:,:), param_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_cg_IDR_HPS, IDR_list(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRHPS, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cg_cutoffdist_126
    cutoff2         =  cutoff * cutoff
    epsilon         =  enefunc%cg_IDR_HPS_epsilon
    num_cg_IDR_HPS  => enefunc%num_cg_IDR_HPS
    IDR_list        => enefunc%cg_IDR_HPS_list
    param_lambda    => enefunc%cg_IDR_HPS_lambda
    param_sigma     => enefunc%cg_IDR_HPS_sigma

    num_nb15_pbc    => pairlist%num_cg_idr_hps_calc
    nb15_calc_list  => pairlist%cg_idr_hps_list

    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                              &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15, ehps_temp,   &
    !$omp         dij, rij2, inv_rij2, inv_rij6, term_lj6,      &
    !$omp         term_lj12, grad_coef, work, force_local,      &
    !$omp         iatmcls, jatmcls, lambda, sigma, ehps_ene, scale1, scale2)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_IDR_HPS, nthread

      ix = IDR_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_pbc(i)
      iatmcls  = atmcls(ix)
      force_local(1:3) = 0.0_wp

      !!dir$ simd
      !ocl simd
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        scale1 = merge(0.0_wp, 1.0_wp, rij2 > cutoff2)

        jatmcls = atmcls(j)
        lambda  = param_lambda(jatmcls,iatmcls)
        sigma   = param_sigma (jatmcls,iatmcls)

        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        term_lj6  = sigma * inv_rij6
        term_lj12 = term_lj6 * term_lj6

        scale2 = merge(1.0_wp, lambda, term_lj6 >= 0.5_wp)
        grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
        grad_coef = 4.0_wp*scale2*epsilon*inv_rij2*grad_coef
        grad_coef = grad_coef * scale1

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

!       ! virial
!       !
!       virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
!       virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
!       virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)

    end do

    !$omp end parallel

    call timer(TimerCGIDRHPS, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_force_CG_IDR_HPS


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_MPIPI
  !> @brief        Calculates IDR interaction in the MPIPI model (PBC).
  !! @authors      CT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] empipi   : energy of IDR MPIPI interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_MPIPI(domain, enefunc, pairlist, &
                                       coord, force, empipi, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: empipi(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: cutoff_scale
    real(wp)                  :: cutoff_scale_sqr
    real(wp)                  :: sigma, sigma_sqr
    real(wp)                  :: epsilon, alpha
    integer                   :: mu, mu2
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: inv_rij_sqr
    real(wp)                  :: sig_over_rij_sqr
    real(wp)                  :: sig_over_rij_2muth
    real(wp)                  :: three_2muth
    real(wp)                  :: tmp_val1
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: empipi_temp, empipi_ene
    integer                   :: i, ix, j, k
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp), pointer         :: mpipi_epsilon(:, :)
    real(wp), pointer         :: mpipi_sigma(:, :)
    integer,  pointer         :: mpipi_mu(:, :)
    integer,  pointer         :: mpipi_nu(:, :)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_cg_IDR_MPIPI, IDR_list(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRMPIPI, TimerOn)

    atmcls           => domain%atom_cls_no

    cutoff_scale     = enefunc%cg_WangFrenkel_cutoff_scale
    cutoff_scale_sqr = cutoff_scale * cutoff_scale

    num_cg_IDR_MPIPI => enefunc%num_cg_IDR_MPIPI
    IDR_list         => enefunc%cg_IDR_MPIPI_list
    mpipi_sigma      => enefunc%cg_WangFrenkel_sigma
    mpipi_epsilon    => enefunc%cg_WangFrenkel_epsilon
    mpipi_mu         => enefunc%cg_WangFrenkel_mu
    mpipi_nu         => enefunc%cg_WangFrenkel_nu

    num_nb15_pbc     => pairlist%num_cg_idr_mpipi_calc
    nb15_calc_list   => pairlist%cg_idr_mpipi_list


    !$omp parallel default(shared)                      &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15,      &
    !$omp         empipi_temp, dij, rij2,               &
    !$omp         sigma, sigma_sqr, epsilon, alpha,     &
    !$omp         mu, mu2, cutoff, cutoff_sqr,          &
    !$omp         inv_rij_sqr,  three_2muth,            &
    !$omp         sig_over_rij_sqr, sig_over_rij_2muth, &
    !$omp         tmp_val1,                             &
    !$omp         grad_coef, work, force_local,         &
    !$omp         iatmcls, jatmcls, empipi_ene)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_IDR_MPIPI, nthread

      ix = IDR_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_pbc(i)
      iatmcls  = atmcls(ix)
      empipi_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = coord(j,1) - rtmp(1)
        dij(2) = coord(j,2) - rtmp(2)
        dij(3) = coord(j,3) - rtmp(3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        jatmcls = atmcls(j)

        sigma   = mpipi_sigma (iatmcls,jatmcls)
        sigma_sqr  = sigma * sigma
        cutoff_sqr = sigma_sqr * cutoff_scale_sqr

        if (rij2 > cutoff_sqr) cycle

        epsilon     = mpipi_epsilon(iatmcls, jatmcls)
        mu          = mpipi_mu(iatmcls, jatmcls)
        mu2         = mu + mu
        three_2muth = 3.0_wp**mu2
        alpha       = 2.0_wp * three_2muth * (1.5_wp / (three_2muth-1))**3.0_wp

        inv_rij_sqr       = 1.0_wp / rij2
        sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr
        sig_over_rij_2muth = sig_over_rij_sqr**mu

        tmp_val1   = three_2muth * sig_over_rij_2muth - 1.0_wp
        empipi_ene = epsilon*alpha*(sig_over_rij_2muth-1)*tmp_val1*tmp_val1

        empipi_temp = empipi_temp + empipi_ene
        grad_coef = mu2 * epsilon * alpha * &
              sig_over_rij_2muth * inv_rij_sqr * tmp_val1 * &
              (3.0_wp * tmp_val1 + 2.0_wp - 2.0_wp * three_2muth)

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        ! virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        ! virial(2,1,id+1) = virial(2,1,id+1) - dij(2)*work(1)
        ! virial(3,1,id+1) = virial(3,1,id+1) - dij(3)*work(1)
        ! virial(1,2,id+1) = virial(1,2,id+1) - dij(1)*work(2)
        ! virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        ! virial(3,2,id+1) = virial(3,2,id+1) - dij(3)*work(2)
        ! virial(1,3,id+1) = virial(1,3,id+1) - dij(1)*work(3)
        ! virial(2,3,id+1) = virial(2,3,id+1) - dij(2)*work(3)
        ! virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      empipi(id+1) = empipi(id+1) + empipi_temp

    end do

    !$omp end parallel

    call timer(TimerCGIDRMPIPI, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_MPIPI

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_CG_IDR_MPIPI
  !> @brief        Calculates IDR interaction in the MPIPI model (PBC).
  !! @authors      CT
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] empipi   : energy of IDR MPIPI interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_CG_IDR_MPIPI(domain, enefunc, pairlist, &
                                      coord, force, empipi, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: empipi(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: cutoff_scale
    real(wp)                  :: cutoff_scale_sqr
    real(wp)                  :: sigma, sigma_sqr
    real(wp)                  :: epsilon, alpha
    integer                   :: mu, mu2
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: inv_rij_sqr
    real(wp)                  :: sig_over_rij_sqr
    real(wp)                  :: sig_over_rij_2muth
    real(wp)                  :: three_2muth
    real(wp)                  :: tmp_val1
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: empipi_temp, empipi_ene
    real(wp)                  :: scale_simd
    integer                   :: i, ix, j, k
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp), pointer         :: mpipi_epsilon(:, :)
    real(wp), pointer         :: mpipi_sigma(:, :)
    integer,  pointer         :: mpipi_mu(:, :)
    integer,  pointer         :: mpipi_nu(:, :)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_cg_IDR_MPIPI, IDR_list(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRMPIPI, TimerOn)

    atmcls           => domain%atom_cls_no

    cutoff_scale     = enefunc%cg_WangFrenkel_cutoff_scale
    cutoff_scale_sqr = cutoff_scale * cutoff_scale

    num_cg_IDR_MPIPI => enefunc%num_cg_IDR_MPIPI
    IDR_list         => enefunc%cg_IDR_MPIPI_list
    mpipi_sigma      => enefunc%cg_WangFrenkel_sigma
    mpipi_epsilon    => enefunc%cg_WangFrenkel_epsilon
    mpipi_mu         => enefunc%cg_WangFrenkel_mu
    mpipi_nu         => enefunc%cg_WangFrenkel_nu

    num_nb15_pbc     => pairlist%num_cg_idr_mpipi_calc
    nb15_calc_list   => pairlist%cg_idr_mpipi_list


    !$omp parallel default(shared)                      &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15,      &
    !$omp         empipi_temp, dij, rij2,               &
    !$omp         sigma, sigma_sqr, epsilon, alpha,     &
    !$omp         mu, mu2, cutoff, cutoff_sqr,          &
    !$omp         inv_rij_sqr,  three_2muth,            &
    !$omp         sig_over_rij_sqr, sig_over_rij_2muth, &
    !$omp         tmp_val1,                             &
    !$omp         grad_coef, work, force_local,         &
    !$omp         iatmcls, jatmcls, scale_simd,         &
    !$omp         empipi_ene)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_IDR_MPIPI, nthread

      ix = IDR_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_pbc(i)
      iatmcls  = atmcls(ix)
      force_local(1:3) = 0.0_wp

!     !ocl simd
!     !ocl norecurrence(force)
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1) = coord(j,1) - rtmp(1)
        dij(2) = coord(j,2) - rtmp(2)
        dij(3) = coord(j,3) - rtmp(3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        jatmcls = atmcls(j)

        sigma   = mpipi_sigma (iatmcls,jatmcls)
        sigma_sqr  = sigma * sigma
        cutoff_sqr = sigma_sqr * cutoff_scale_sqr

        if (rij2 > cutoff_sqr) cycle

!       scale_simd = merge(0.0_wp, 1.0_wp, rij2 > cutoff_sqr)

        epsilon     = mpipi_epsilon(iatmcls, jatmcls)
        mu          = mpipi_mu(iatmcls, jatmcls)
        mu2         = mu + mu
        three_2muth = 3.0_wp**mu2
        alpha       = 2.0_wp * three_2muth * (1.5_wp / (three_2muth - 1))**3

        inv_rij_sqr       = 1.0_wp / rij2
        sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr
        sig_over_rij_2muth = sig_over_rij_sqr**mu

        tmp_val1   = three_2muth * sig_over_rij_2muth - 1.0_wp

        grad_coef = mu2 * epsilon * alpha * &
              sig_over_rij_2muth * inv_rij_sqr * tmp_val1 * &
              (3.0_wp * tmp_val1 + 2.0_wp - 2.0_wp * three_2muth)

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        ! virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        ! virial(2,1,id+1) = virial(2,1,id+1) - dij(2)*work(1)
        ! virial(3,1,id+1) = virial(3,1,id+1) - dij(3)*work(1)
        ! virial(1,2,id+1) = virial(1,2,id+1) - dij(1)*work(2)
        ! virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        ! virial(3,2,id+1) = virial(3,2,id+1) - dij(3)*work(2)
        ! virial(1,3,id+1) = virial(1,3,id+1) - dij(1)*work(3)
        ! virial(2,3,id+1) = virial(2,3,id+1) - dij(2)*work(3)
        ! virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)

    end do

    !$omp end parallel

    call timer(TimerCGIDRMPIPI, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_force_CG_IDR_MPIPI


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_KH
  !> @brief        Calculates IDR interaction in the KH model (PBC).
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] ekh      : energy of IDR KH interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_KH(domain, enefunc, pairlist, &
                                      coord, force, ekh, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: ekh(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, term_lj6, term_lj12
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: sigma, epsilon, lambda
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: ekh_temp, ekh_ene
    integer                   :: i, ix, j, k
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls

    real(wp),         pointer :: param_sigma(:,:), param_epsilon(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_cg_IDR_KH, IDR_list(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRKH, TimerOn)

    atmcls          => domain%atom_cls_no

    cutoff          =  enefunc%cg_cutoffdist_126
    cutoff2         =  cutoff * cutoff
    num_cg_IDR_KH   => enefunc%num_cg_IDR_KH
    IDR_list        => enefunc%cg_IDR_KH_list
    param_sigma     => enefunc%cg_IDR_KH_sigma
    param_epsilon   => enefunc%cg_IDR_KH_epsilon

    num_nb15_calc   => pairlist%num_cg_idr_kh_calc
    nb15_calc_list  => pairlist%cg_idr_kh_list

    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                           &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15, ekh_temp, &
    !$omp         dij, rij2, inv_rij2, inv_rij6, term_lj6,   &
    !$omp         term_lj12, grad_coef, work, force_local,   &
    !$omp         iatmcls, jatmcls, epsilon, sigma, lambda, ekh_ene)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_IDR_KH, nthread

      ix = IDR_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_calc(i)
      iatmcls  = atmcls(ix)
      ekh_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 > cutoff2) cycle

        jatmcls = atmcls(j)
        epsilon = param_epsilon(iatmcls,jatmcls)
        sigma   = param_sigma  (iatmcls,jatmcls)
        if (epsilon > 0.0_wp) then
          lambda  = -1.0_wp
        else
          lambda  =  1.0_wp
          epsilon = -epsilon
        end if

        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        term_lj6  = sigma * inv_rij6
        term_lj12 = term_lj6 * term_lj6
        ekh_ene  = 4.0_wp*epsilon*(term_lj12 - term_lj6)

        if (term_lj6 >= 0.5_wp) then
          ekh_temp = ekh_temp + ekh_ene + (1.0_wp-lambda)*epsilon
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*epsilon*inv_rij2*grad_coef
        else
          ekh_temp = ekh_temp + lambda*ekh_ene
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*lambda*epsilon*inv_rij2*grad_coef
        end if

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      ekh(id+1) = ekh(id+1) + ekh_temp

    end do

    !$omp end parallel

    call timer(TimerCGIDRKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_KH
  !> @brief        Calculates interaction in the KH model.
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] ekh      : energy of KH interactions
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_KH(domain, enefunc, pairlist, &
                                  coord, force, ekh, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: ekh(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, term_lj6, term_lj12
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: grad_coef, work(3)
    real(wp)                  :: sigma, epsilon, lambda
    real(wp)                  :: alpha, eps0
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: ekh_temp, ekh_ene
    integer                   :: i, ix, j, k
    integer                   :: chain_idx, chain_idy
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: iatmcls, jatmcls
    integer                   :: kh_model

    real(wp)                  :: cg_KH_mod_A_lambda, cg_KH_mod_A_eps_0
    real(wp)                  :: cg_KH_mod_B_lambda, cg_KH_mod_B_eps_0
    real(wp)                  :: cg_KH_mod_C_lambda, cg_KH_mod_C_eps_0
    real(wp)                  :: cg_KH_mod_D_lambda, cg_KH_mod_D_eps_0
    real(wp)                  :: cg_KH_mod_E_lambda, cg_KH_mod_E_eps_0
    real(wp)                  :: cg_KH_mod_F_lambda, cg_KH_mod_F_eps_0

    real(wp),         pointer :: param_sigma(:,:), param_epsilon(:,:)
    integer(1),       pointer :: KH_mol_pair(:,:)
    integer,          pointer :: atmcls(:), chain_id(:)
    integer,          pointer :: num_cg_KH, KH_list(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGKH, TimerOn)

    atmcls          => domain%atom_cls_no
    chain_id        => domain%mol_chain_id

    cutoff          =  enefunc%cg_cutoffdist_126
    num_cg_KH       => enefunc%num_cg_KH
    KH_list         => enefunc%cg_KH_list
    param_sigma     => enefunc%cg_KH_sigma
    param_epsilon   => enefunc%cg_KH_epsilon
    kh_mol_pair     => enefunc%cg_KH_mol_pair

    num_nb15_calc   => pairlist%num_cg_kh_calc
    nb15_calc_list  => pairlist%cg_kh_list

    cutoff2         =  cutoff * cutoff

    cg_KH_mod_A_lambda = enefunc%cg_KH_mod_A_lambda
    cg_KH_mod_A_eps_0  = enefunc%cg_KH_mod_A_eps_0
    cg_KH_mod_B_lambda = enefunc%cg_KH_mod_B_lambda
    cg_KH_mod_B_eps_0  = enefunc%cg_KH_mod_B_eps_0
    cg_KH_mod_C_lambda = enefunc%cg_KH_mod_C_lambda
    cg_KH_mod_C_eps_0  = enefunc%cg_KH_mod_C_eps_0
    cg_KH_mod_D_lambda = enefunc%cg_KH_mod_D_lambda
    cg_KH_mod_D_eps_0  = enefunc%cg_KH_mod_D_eps_0 
    cg_KH_mod_E_lambda = enefunc%cg_KH_mod_E_lambda
    cg_KH_mod_E_eps_0  = enefunc%cg_KH_mod_E_eps_0 
    cg_KH_mod_F_lambda = enefunc%cg_KH_mod_F_lambda
    cg_KH_mod_F_eps_0  = enefunc%cg_KH_mod_F_eps_0 

    !$omp parallel default(shared)                           &
    !$omp private(id, i, ix, j, k, rtmp, num_nb15, ekh_temp, &
    !$omp         dij, rij2, inv_rij2, inv_rij6, term_lj6,   &
    !$omp         term_lj12, grad_coef, work, force_local,   &
    !$omp         iatmcls, jatmcls, chain_idx, chain_idy,    &
    !$omp         epsilon, sigma, lambda, alpha, eps0,       &
    !$omp         ekh_ene, kh_model)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, num_cg_KH, nthread

      ix = KH_list(i)
      rtmp(1) = coord(ix,1)
      rtmp(2) = coord(ix,2)
      rtmp(3) = coord(ix,3)

      num_nb15 = num_nb15_calc(ix)
      iatmcls  = atmcls(ix)
      chain_idx = chain_id(ix)
      ekh_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,ix)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 > cutoff2) cycle

        chain_idy = chain_id(j)

        kh_model = KH_mol_pair(chain_idy, chain_idx)

        select case (kh_model)
        case ( 1 )
          alpha = cg_KH_mod_A_lambda
          eps0  = cg_KH_mod_A_eps_0
        case ( 2 )
          alpha = cg_KH_mod_B_lambda
          eps0  = cg_KH_mod_B_eps_0
        case ( 3 )
          alpha = cg_KH_mod_C_lambda
          eps0  = cg_KH_mod_C_eps_0
        case ( 4 )
          alpha = cg_KH_mod_D_lambda
          eps0  = cg_KH_mod_D_eps_0
        case ( 5 )
          alpha = cg_KH_mod_E_lambda
          eps0  = cg_KH_mod_E_eps_0
        case ( 6 )
          alpha = cg_KH_mod_F_lambda
          eps0  = cg_KH_mod_F_eps_0
        end select

        jatmcls = atmcls(j)
        epsilon = alpha*(param_epsilon(iatmcls,jatmcls)-eps0)*0.593_wp
!       epsilon = alpha*(param_epsilon(iatmcls,jatmcls)-eps0)
        sigma   = param_sigma  (iatmcls,jatmcls)

        if (epsilon > 0.0_wp) then
          lambda  = -1.0_wp
        else
          lambda  =  1.0_wp
          epsilon = -epsilon
        end if

        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        term_lj6  = sigma * inv_rij6
        term_lj12 = term_lj6 * term_lj6
        ekh_ene  = 4.0_wp*epsilon*(term_lj12 - term_lj6)

        if (term_lj6 >= 0.5_wp) then
          ekh_temp = ekh_temp + ekh_ene + (1.0_wp-lambda)*epsilon
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*epsilon*inv_rij2*grad_coef
        else
          ekh_temp = ekh_temp + lambda*ekh_ene
          grad_coef = 6.0_wp*term_lj6 - 12.0_wp*term_lj12
          grad_coef = 4.0_wp*lambda*epsilon*inv_rij2*grad_coef
        end if

        work(1) = grad_coef * dij(1)
        work(2) = grad_coef * dij(2)
        work(3) = grad_coef * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      ekh(id+1) = ekh(id+1) + ekh_temp

    end do

    !$omp end parallel

    call timer(TimerCGKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_exv
  !> @brief        Calculates excluded volume energy with pairlist (PBC).
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_exv(domain, enefunc, pairlist, &
                                    coord, force, eexv, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: eexv(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: grad_coef_exv, work(3)
    real(wp)                  :: rtmp(3)
    real(wp)                  :: force_local(3)
    real(wp)                  :: inv_rij_sqr, inv_rij_6th
    real(wp)                  :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)                  :: exv_temp
    integer                   :: i, ix, iy, j, k, start_i
    integer                   :: i_base, j_base
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: ncell

    real(wp),         pointer :: param_exv_sigma(:,:), param_exv_epsilon
    integer,          pointer :: natom(:)
    integer,          pointer :: base_type(:)
    integer,          pointer :: cg_dna_list(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDNAexv, TimerOn)

    natom              => domain%num_atom
    base_type          => domain%na_base_type

    param_exv_sigma    => enefunc%cgDNA_exv_sigma
    param_exv_epsilon  => enefunc%cgDNA_exv_epsilon
    cg_dna_list        => enefunc%cg_dna_list

    num_nb15_calc      => pairlist%num_cg_DNA_exv_calc
    nb15_calc_list     => pairlist%cg_DNA_exv_list

    ncell              =  domain%num_cell_local + domain%num_cell_boundary

    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, j, iy, k, rtmp, num_nb15, exv_temp, &
    !$omp         dij, rij2, inv_rij_sqr, inv_rij_6th, sigma,    &
    !$omp         sigma_sqr, sigma_6th, sig_over_rij_6th,        &
    !$omp         sig_over_rij_12th, grad_coef_exv,              &
    !$omp         work, force_local, i_base, j_base, start_i)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, enefunc%num_cg_DNA, nthread

      ix = cg_dna_list(i)
      rtmp(1:3) = coord(ix,1:3)

      num_nb15 = num_nb15_calc(i)
      i_base   = base_type(ix)
      exv_temp = 0.0_wp
      force_local(1:3) = 0.0_wp

      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !   
        dij(1) = rtmp(1) - coord(j,1) 
        dij(2) = rtmp(2) - coord(j,2) 
        dij(3) = rtmp(3) - coord(j,3) 
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        j_base  = base_type(j)
        sigma = param_exv_sigma(i_base,j_base)
        sigma_sqr = sigma * sigma

        if (rij2 > sigma_sqr) cycle

        sigma_6th         = sigma_sqr * sigma_sqr * sigma_sqr

        inv_rij_sqr       = 1.0_wp / rij2
        inv_rij_6th       = inv_rij_sqr * inv_rij_sqr * inv_rij_sqr

        sig_over_rij_6th  = sigma_6th * inv_rij_6th
        sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

        exv_temp          = exv_temp &
                + param_exv_epsilon*(sig_over_rij_12th-2*sig_over_rij_6th+1)

        ! gradient
        !
        grad_coef_exv =  inv_rij_sqr*12.0_wp*param_exv_epsilon &
                         *(sig_over_rij_6th - sig_over_rij_12th)

        work(1) = grad_coef_exv * dij(1)
        work(2) = grad_coef_exv * dij(2)
        work(3) = grad_coef_exv * dij(3)

        ! store force
        !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1:3,id+1) = force(ix,1:3,id+1) + force_local(1:3)
      eexv(id+1) = eexv(id+1) + exv_temp

    end do

    !$omp end parallel

    call timer(TimerCGDNAexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_exv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_pairing
  !> @brief        Calculates base-pairing energy with pairlist (PBC).
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epair    : base-pairing energy
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_pairing(domain, enefunc, pairlist, &
                                             coord_pbc, force, epair,   &
                                             virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord_pbc(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: epair(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    integer           :: id
    integer           :: natom_all, num_nb15
    integer           :: i, iix, ix, j, k, l, m, ig_x, ig_y, start_i
    integer           :: ncell
    real(wp)          :: cutoff, cutoff_sqr
    logical           :: cg_infinite_DNA

    integer           :: ig_S1, ig_S3, ig_B2, ig_B4, ig_B5, ig_B6
    integer           :: ix_S1, ix_S3, ix_B2, ix_B4, ix_B5, ix_B6
    integer           :: icel_S1, icel_S3, icel_B2, icel_B4, icel_B5, icel_B6
    integer           :: type_B2, type_B4, type_B5, type_B6
    integer           :: aindex(1:4)
    real(wp)          :: viri(3)
    real(wp)          :: nopbc_S1(3), nopbc_B2(3), nopbc_S3(3)
    real(wp)          :: nopbc_B4(3), nopbc_B5(3), nopbc_B6(3)
    real(wp)          :: pbc_S1(3), pbc_B2(3), pbc_S3(3)
    real(wp)          :: pbc_B4(3), pbc_B5(3), pbc_B6(3)
    real(wp)          :: pbc(3), nopbc(3)
    real(wp)          :: d21(3), r21_sqr, r21, r21_inv, e21(3)
    real(wp)          :: d42(3), r42_sqr, r42, r42_inv, e42(3)
    real(wp)          :: d24(3), e24(3)
    real(wp)          :: d23(3)
    real(wp)          :: d43(3), r43_sqr, r43, r43_inv, e43(3)
    real(wp)          :: d25(3), r25_sqr, r25, r25_inv, e25(3)
    real(wp)          :: d46(3), r46_sqr, r46, r46_inv, e46(3)
    real(wp)          :: cos_124, cos_342, sin_124, sin_342
    real(wp)          :: cos_125, cos_346, sin_125, sin_346
    real(wp)          :: cos_1243, sin_1243
    real(wp)          :: angle_124, angle_342
    real(wp)          :: angle_125, angle_346
    real(wp)          :: angle_1243
    real(wp)          :: delta_angle_124, delta_angle_342
    real(wp)          :: delta_angle_125, delta_angle_346
    real(wp)          :: delta_angle_1243
    real(wp)          :: abs_delta_angle_124, abs_delta_angle_342
    real(wp)          :: abs_delta_angle_125, abs_delta_angle_346
    real(wp)          :: abs_delta_angle_1243
    real(wp)          :: grad_angle_124(1:6)
    real(wp)          :: grad_angle_342(1:6)
    real(wp)          :: grad_angle_125(1:6)
    real(wp)          :: grad_angle_346(1:6)
    real(wp)          :: grad_angle_1243(1:6)
    real(wp)          :: force_angle_124_tmp(1:6)
    real(wp)          :: force_angle_342_tmp(1:6)
    real(wp)          :: force_angle_125_tmp(1:6)
    real(wp)          :: force_angle_346_tmp(1:6)
    real(wp)          :: force_angle_1243_tmp(1:6)
    real(wp)          :: cos_dih_24, sin_dih_24
    real(wp)          :: grad_dih_24(1:9)
    real(wp)          :: force_dih_24_tmp(1:9)
    real(wp)          :: grad_repl(1:3), ene_repl
    real(wp)          :: grad_attr(1:3), ene_attr
    real(wp)          :: force_repl_tmp(1:3)
    real(wp)          :: force_attr_tmp(1:3)
    real(wp)          :: grad_coef_phi, ene_coef_phi
    real(wp)          :: grad_coef_theta_1, ene_coef_theta_1
    real(wp)          :: grad_coef_theta_2, ene_coef_theta_2
    real(wp)          :: grad_coef_theta_3, ene_coef_theta_3
    real(wp)          :: grad_coef_theta_cs, ene_coef_theta_cs
    real(wp)          :: force_coef_tmp
    real(wp)          :: var_tmp
    real(wp)          :: force_S1(3), force_S3(3)
    real(wp)          :: force_B2(3), force_B4(3)
    real(wp)          :: force_B5(3), force_B6(3)

    real(wp)          :: bp_theta_1_0
    real(wp)          :: bp_theta_2_0
    real(wp)          :: bp_theta_3_0
    real(wp)          :: bp_phi_1_0
    real(wp)          :: bp_sigma
    real(wp)          :: bp_epsilon
    real(wp)          :: cs_1_epsilon
    real(wp)          :: cs_1_sigma
    real(wp)          :: cs_1_theta_cs_0
    real(wp)          :: cs_2_epsilon
    real(wp)          :: cs_2_sigma
    real(wp)          :: cs_2_theta_cs_0

    integer,          pointer :: natom(:), nbase(:), base_list(:)
    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: cg_base_list(:)
    integer,          pointer :: num_nb15_pbc(:)
    integer,          pointer :: nb15_calc_list(:,:)

    integer,          pointer :: base_type(:)
    integer,          pointer :: chain_id(:)
    integer,          pointer :: DNA_end(:,:)
    real(wp),         pointer :: param_bp_theta_1(:)
    real(wp),         pointer :: param_bp_theta_2(:)
    real(wp),         pointer :: param_bp_theta_3(:)
    real(wp),         pointer :: param_bp_phi_1(:)
    real(wp),         pointer :: param_bp_sigma(:)
    real(wp),         pointer :: param_bp_epsilon(:)
    real(wp),         pointer :: param_cs_1_epsilon(:,:)
    real(wp),         pointer :: param_cs_1_sigma(:,:)
    real(wp),         pointer :: param_cs_1_theta_cs(:,:)
    real(wp),         pointer :: param_cs_2_epsilon(:,:)
    real(wp),         pointer :: param_cs_2_sigma(:,:)
    real(wp),         pointer :: param_cs_2_theta_cs(:,:)

    real(wp)          :: param_bp_alpha
    real(wp)          :: param_bp_K
    real(wp)          :: param_cs_alpha
    real(wp)          :: param_cs_K
    real(wp)          :: bp_theta_threshold_2
    real(wp)          :: bp_theta_threshold_1
    real(wp)          :: cs_theta_threshold_2
    real(wp)          :: cs_theta_threshold_1

    real(wp)          :: e_tmp_bp
    real(wp)          :: e_tmp_cs

    !==========================================================================
    !                                   func
    ! U_bp = U_repl + U_attr
    !
    ! U_repl = U_morse                          ( r < r0 )
    !        = 0                                ( r > r0 )
    !
    ! U_attr = f_bell * ( - epsilon )           ( r < r0 )
    !        = f_bell * ( U_morse - epsilon)    ( r > r0 )
    !
    ! U_cstk = U_attr
    !
    ! U_attr = f_bell * ( - epsilon )           ( r < r0 )
    !        = f_bell * ( U_morse - epsilon)    ( r > r0 )
    !--------------------------------------------------------------------------
    !                                 topology
    !    d21   d42   d43
    ! S1----B2    B4----S3
    ! |       \  /       |
    ! |        \/        |   phi: dihedral S1-B2=B4-S3
    ! |   d46  /\  d25   |
    ! |       /  \       |
    ! o-----B6    B5-----o
    !
    !==========================================================================


    call timer(TimerNonBond, TimerOn)
    call timer(TimerBasePair, TimerOn)

    cutoff               = enefunc%cg_cutoffdist_DNAbp
    cutoff_sqr           = cutoff * cutoff

    ncell                =  domain%num_cell_local + domain%num_cell_boundary
    natom_all            =  domain%num_atom_all
    natom                => domain%num_atom
    nbase                => domain%num_base
    base_list            => domain%base_list
    base_type            => domain%na_base_type
    start_atom           => domain%start_atom
    id_l2g               => domain%id_l2g
    id_g2l               => domain%id_g2l
    chain_id             => domain%mol_chain_id

    num_nb15_pbc         => pairlist%num_cg_DNA_base_calc
    nb15_calc_list       => pairlist%cg_DNA_base_list

    DNA_end              => enefunc%DNA_end
    cg_base_list         => enefunc%cg_base_list
    param_bp_theta_1     => enefunc%base_pair_theta_1
    param_bp_theta_2     => enefunc%base_pair_theta_2
    param_bp_theta_3     => enefunc%base_pair_theta_3
    param_bp_phi_1       => enefunc%base_pair_phi_1
    param_bp_sigma       => enefunc%base_pair_sigma
    param_bp_epsilon     => enefunc%base_pair_epsilon
    param_cs_1_epsilon   => enefunc%base_cross_1_epsilon
    param_cs_1_sigma     => enefunc%base_cross_1_sigma
    param_cs_1_theta_cs  => enefunc%base_cross_1_theta_cs
    param_cs_2_epsilon   => enefunc%base_cross_2_epsilon
    param_cs_2_sigma     => enefunc%base_cross_2_sigma
    param_cs_2_theta_cs  => enefunc%base_cross_2_theta_cs

    param_bp_alpha       = enefunc%base_pair_alpha
    param_bp_K           = enefunc%base_pair_K
    param_cs_alpha       = enefunc%base_cross_alpha
    param_cs_K           = enefunc%base_cross_K
    cg_infinite_DNA      = enefunc%cg_infinite_DNA

    bp_theta_threshold_2 = PI / param_bp_K
    bp_theta_threshold_1 = bp_theta_threshold_2 / 2.0_wp
    cs_theta_threshold_2 = PI / param_cs_K
    cs_theta_threshold_1 = cs_theta_threshold_2 / 2.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel default(shared)                             &
    !$omp private(id, nopbc, pbc,                                &
    !$omp         nopbc_S1, nopbc_B2, nopbc_S3, nopbc_B4, nopbc_B5, nopbc_B6, &
    !$omp         pbc_S1, pbc_B2, pbc_S3, pbc_B4, pbc_B5, pbc_B6,             &
    !$omp         i, iix, ix, j, k, l, m, ig_x, ig_y,        &
    !$omp         ig_S1, ig_B2, ig_S3, ig_B4, ig_B5, ig_B6,  &
    !$omp         ix_S1, ix_B2, ix_S3, ix_B4, ix_B5, ix_B6,  &
    !$omp         icel_S1, icel_B2, icel_S3, icel_B4, icel_B5, icel_B6,  &
    !$omp         type_B2, type_B4, type_B5, type_B6,        &
    !$omp         aindex,                     &
    !$omp         d21, r21_sqr, r21, r21_inv, e21,           &
    !$omp         d42, r42_sqr, r42, r42_inv, e42, d24, e24, &
    !$omp         d23, d43, r43_sqr, r43, r43_inv, e43,      &
    !$omp         d25, r25_sqr, r25, r25_inv, e25,           &
    !$omp         d46, r46_sqr, r46, r46_inv, e46,           &
    !$omp         cos_124, cos_342, sin_124, sin_342,        &
    !$omp         cos_125, cos_346, sin_125, sin_346,        &
    !$omp         cos_1243, sin_1243,                        &
    !$omp         angle_124, angle_342,                      &
    !$omp         angle_125, angle_346,                      &
    !$omp         angle_1243,                                &
    !$omp         delta_angle_124, delta_angle_342,          &
    !$omp         delta_angle_125, delta_angle_346,          &
    !$omp         delta_angle_1243,                          &
    !$omp         abs_delta_angle_124, abs_delta_angle_342,  &
    !$omp         abs_delta_angle_125, abs_delta_angle_346,  &
    !$omp         abs_delta_angle_1243,                      &
    !$omp         grad_angle_124, force_angle_124_tmp,       &
    !$omp         grad_angle_342, force_angle_342_tmp,       &
    !$omp         grad_angle_125, force_angle_125_tmp,       &
    !$omp         grad_angle_346, force_angle_346_tmp,       &
    !$omp         grad_angle_1243, force_angle_1243_tmp,     &
    !$omp         cos_dih_24, sin_dih_24,                    &
    !$omp         grad_dih_24, force_dih_24_tmp,             &
    !$omp         grad_repl, ene_repl,                       &
    !$omp         grad_attr, ene_attr,                       &
    !$omp         force_repl_tmp, force_attr_tmp,            &
    !$omp         grad_coef_phi, ene_coef_phi,               &
    !$omp         grad_coef_theta_1, ene_coef_theta_1,       &
    !$omp         grad_coef_theta_2, ene_coef_theta_2,       &
    !$omp         grad_coef_theta_3, ene_coef_theta_3,       &
    !$omp         grad_coef_theta_cs, ene_coef_theta_cs,     &
    !$omp         force_coef_tmp, var_tmp,                   &
    !$omp         force_S1, force_S3,                        &
    !$omp         force_B2, force_B4,                        &
    !$omp         force_B5, force_B6,                        &
    !$omp         bp_theta_1_0, bp_theta_2_0,                &
    !$omp         bp_theta_3_0, bp_phi_1_0,                  &
    !$omp         bp_sigma, bp_epsilon,                      &
    !$omp         cs_1_epsilon, cs_2_epsilon,                &
    !$omp         cs_1_sigma, cs_2_sigma,                    &
    !$omp         cs_1_theta_cs_0, cs_2_theta_cs_0,          &
    !$omp         e_tmp_bp, e_tmp_cs, num_nb15, start_i, viri)
    !
#ifdef OMP
    id      = omp_get_thread_num()
    nthread = omp_get_num_threads()
#else
    id      = 0
    nthread = 1
#endif

    do i = id+1, enefunc%num_cg_base, nthread

      ix = cg_base_list(i)

      num_nb15 = num_nb15_pbc(i)

      pbc  (1:3) = coord_pbc(ix,1:3)
      ig_x       = id_l2g(ix)

      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)
        ig_y = id_l2g(j)

        if (ig_y > ig_x) then
          pbc_B2  (1:3) = pbc(1:3)
          pbc_B4  (1:3) = coord_pbc(j,1:3)
          d42(1:3)      = pbc_B2(1:3) - pbc_B4(1:3) 
        else
          pbc_B2  (1:3) = coord_pbc(j,1:3)
          pbc_B4  (1:3) = pbc(1:3)
          d42(1:3)      = pbc_B2(1:3) - pbc_B4(1:3) 
        end if

        r42_sqr = d42(1) * d42(1) + d42(2) * d42(2) + d42(3) * d42(3)

        if ( r42_sqr > cutoff_sqr ) then
          cycle
        end if

        if (ig_y > ig_x) then
          ig_B2   = ig_x
          ig_B4   = ig_y
          ix_B2   = ix 
          ix_B4   = j
          type_B2 = base_type(ix_B2)
          type_B4 = base_type(j)

        else
          ig_B2   = ig_y
          ig_B4   = ig_x
          ix_B2   = j
          ix_B4   = ix 
          type_B2 = base_type(j)
          type_B4 = base_type(ix_B4)
        end if

        ig_S1   = ig_B2 - 1
        ig_S3   = ig_B4 - 1
        ig_B5   = ig_B4 - 3
        ig_B6   = ig_B2 + 3
        ix_S1   = id_g2l(ig_S1) 
        ix_S3   = id_g2l(ig_S3) 
        ix_B5   = id_g2l(ig_B5) 
        ix_B6   = id_g2l(ig_B6) 

        if (cg_infinite_DNA) then
          if (ig_B6 > natom_all .or. chain_id(ix_B6) /= chain_id(ix_B2)) then
            ig_B6 = DNA_end(1,chain_id(ix_B2))
            ix_B6 = id_g2l(ig_B6)
          end if
          if (ig_B5 <= 0 .or. chain_id(ix_B5) /= chain_id(ix_B4) ) then
            ig_B5 = DNA_end(2,chain_id(ix_B4))
            ix_B5 = id_g2l(ig_B5)
          end if
        end if

        pbc_S1  (1:3) = coord_pbc(ix_S1,1:3)
        pbc_S3  (1:3) = coord_pbc(ix_S3,1:3)
        pbc_B5  (1:3) = coord_pbc(ix_B5,1:3)
        pbc_B6  (1:3) = coord_pbc(ix_B6,1:3)

        bp_theta_1_0  = param_bp_theta_1 (type_B2)
        bp_theta_2_0  = param_bp_theta_2 (type_B2)
        bp_theta_3_0  = param_bp_theta_3 (type_B2)
        bp_phi_1_0    = param_bp_phi_1   (type_B2)
        bp_sigma      = param_bp_sigma   (type_B2)
        bp_epsilon    = param_bp_epsilon (type_B2)

        ! information of S1, S3, B5 and B6
        e_tmp_bp = 0.0_wp
        e_tmp_cs = 0.0_wp

        force_S1(1:3) = 0.0_wp
        force_S3(1:3) = 0.0_wp
        force_B2(1:3) = 0.0_wp
        force_B4(1:3) = 0.0_wp
        force_B5(1:3) = 0.0_wp
        force_B6(1:3) = 0.0_wp

        !=====================================================================
        !               Base Pairing Energy / Force Calculation             
        !=====================================================================
        !
        !  ____                   ____       _      _
        ! | __ )  __ _ ___  ___  |  _ \ __ _(_)_ __(_)_ __   __ _
        ! |  _ \ / _` / __|/ _ \ | |_) / _` | | '__| | '_ \ / _` |
        ! | |_) | (_| \__ \  __/ |  __/ (_| | | |  | | | | | (_| |
        ! |____/ \__,_|___/\___| |_|   \__,_|_|_|  |_|_| |_|\__, |
        !                                                   |___/
        !

        ! -----------------
        ! 1 -- 2 <== 4 -- 3
        ! -----------------
        !
        d24(1:3) = -d42(1:3)
        r42      = sqrt(r42_sqr)
        r42_inv  = 1.0_wp / r42
        e42(1:3) = d42(1:3) * r42_inv
        e24(1:3) = d24(1:3) * r42_inv

        ! -----------------
        ! 1 <== 2 -- 4 -- 3
        ! -----------------
        !
        d21(1:3) = pbc_S1(1:3) - pbc_B2(1:3)          
        r21_sqr  = d21(1) * d21(1) + d21(2) * d21(2) + d21(3) * d21(3)
        r21      = sqrt(r21_sqr)
        r21_inv  = 1.0_wp / r21
        e21(1:3) = d21(1:3) * r21_inv

        ! -----------------
        ! 1 -- 2 -- 4 ==> 3
        ! -----------------
        !
        d43(1:3) = pbc_S3(1:3) - pbc_B4(1:3) 
        r43_sqr  = d43(1) * d43(1) + d43(2) * d43(2) + d43(3) * d43(3)
        r43      = sqrt(r43_sqr)
        r43_inv  = 1.0_wp / r43
        e43(1:3) = d43(1:3) * r43_inv

        ! -----------------------
        ! Angle B2: 1 -- ~2~ -- 4
        ! -----------------------
        !
        cos_124 = e21(1) * e24(1) + e21(2) * e24(2) + e21(3) * e24(3)
        if (cos_124 >  1.0_wp) cos_124 =  1.0_wp
        if (cos_124 < -1.0_wp) cos_124 = -1.0_wp
        sin_124 = sqrt(1.0_wp - cos_124 * cos_124)
        if (sin_124 < EPS) sin_124 = EPS
        angle_124 = acos(cos_124)

        ! -----------------------
        ! Angle B4: 2 -- ~4~ -- 3
        ! -----------------------
        !
        cos_342 = e43(1) * e42(1) + e43(2) * e42(2) + e43(3) * e42(3)
        if (cos_342 >  1.0_wp) cos_342 =  1.0_wp
        if (cos_342 < -1.0_wp) cos_342 = -1.0_wp
        sin_342 = sqrt(1.0_wp - cos_342 * cos_342)
        if (sin_342 < EPS) sin_342 = EPS
        angle_342 = acos(cos_342)

        ! -------------------------
        ! Dihedral 1 -- 2 == 4 -- 3
        ! -------------------------
        !
        call calculate_nonlocal_dihedral(d21, d42, d43, cos_dih_24, &
                                         sin_dih_24, grad_dih_24)
        ene_coef_phi  =   0.5_wp * (1.0_wp + cos_dih_24 * cos(bp_phi_1_0) &
                                    + sin_dih_24 * sin(bp_phi_1_0))
        grad_coef_phi = - 0.5_wp * (         sin_dih_24 * cos(bp_phi_1_0) &
                                    - cos_dih_24 * sin(bp_phi_1_0))

        ! ============================================================
        ! basepairing interaction: energy/force calculation @@@@@@@...
        ! ============================================================
        !
        if (r42 < bp_sigma) then

          call calculate_repulsive(d24, r42, bp_sigma, param_bp_alpha, &
                                   bp_epsilon, ene_repl, grad_repl)

          e_tmp_bp = e_tmp_bp + ene_repl

          force_repl_tmp(1:3) = - grad_repl(1:3)
          force_B2(1:3)       = force_B2(1:3) + force_repl_tmp(1:3)
          force_B4(1:3)       = force_B4(1:3) - force_repl_tmp(1:3)

        end if

        delta_angle_124     = angle_124 - bp_theta_1_0
        abs_delta_angle_124 = abs(delta_angle_124)
        delta_angle_342     = angle_342 - bp_theta_2_0
        abs_delta_angle_342 = abs(delta_angle_342)

        call calculate_attractive(d24, r42, bp_sigma, param_bp_alpha, &
                                  bp_epsilon, ene_attr, grad_attr)

        if ( abs_delta_angle_124 <= bp_theta_threshold_1  ) then
          if ( abs_delta_angle_342 <= bp_theta_threshold_1  ) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_attr

            force_dih_24_tmp(1:9) = -grad_coef_phi*ene_attr*grad_dih_24(1:9)
            force_S1(1:3)         = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)         = force_B2(1:3) - force_dih_24_tmp(1:3) &
                                  + force_dih_24_tmp(4:6)
            force_B4(1:3)         = force_B4(1:3) - force_dih_24_tmp(4:6) &
                                  - force_dih_24_tmp(7:9)
            force_S3(1:3)         = force_S3(1:3) + force_dih_24_tmp(7:9)
            force_attr_tmp(1:3)   = - ene_coef_phi * grad_attr(1:3)
            force_B2(1:3)         = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)         = force_B4(1:3) - force_attr_tmp(1:3)

          else if ( abs_delta_angle_342 <= bp_theta_threshold_2  ) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, &
                                 grad_angle_342)
            var_tmp           = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K &
                                                 * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_2 * ene_attr

            force_coef_tmp           = grad_coef_phi*ene_coef_theta_2*ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) &
                                     + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) &
                                     - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)
            force_coef_tmp           = ene_coef_phi*grad_coef_theta_2*ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) &
                                       + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) &
                                     - force_angle_342_tmp(1:3) &
                                     - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) &
                                     + force_angle_342_tmp(4:6)
            force_coef_tmp           = ene_coef_phi * ene_coef_theta_2
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          end if

        else if ( abs_delta_angle_124 <= bp_theta_threshold_2  ) then

          call calculate_angle(d21, d24, r21, r42, cos_124, sin_124, &
                               grad_angle_124)
          var_tmp = sin(param_bp_K * delta_angle_124)
          ene_coef_theta_1  = var_tmp * var_tmp
          grad_coef_theta_1 = param_bp_K * sin(2.0_wp * param_bp_K &
                                               * delta_angle_124)

          if ( abs_delta_angle_342 <= bp_theta_threshold_1  ) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 * ene_attr
            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 &
                                      * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) &
                                     + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) &
                                     - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)
            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 &
                                      * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) &
                                     + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) &
                                     - force_angle_124_tmp(1:3) &
                                     - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) &
                                     + force_angle_124_tmp(4:6)
            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          else if ( abs_delta_angle_342 <= bp_theta_threshold_2  ) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, &
                                 grad_angle_342)
            var_tmp = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K &
                                                 * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 &
                                 * ene_coef_theta_2 * ene_attr
            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 &
                                      * ene_coef_theta_2 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) &
                                     + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) &
                                     - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)
            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 &
                                      * ene_coef_theta_2 * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) &
                                     + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) &
                                     - force_angle_124_tmp(1:3) &
                                     - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) &
                                     + force_angle_124_tmp(4:6)
            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 &
                                      * grad_coef_theta_2 * ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) &
                                     + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) &
                                     - force_angle_342_tmp(1:3) &
                                     - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) &
                                     + force_angle_342_tmp(4:6)
            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 &
                                      * ene_coef_theta_2
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          end if
        end if

        ! ====================================================================
        ! base cross-stacking interaction: energy/force calculation @@@@@@@...
        ! ====================================================================
        !   ____                     ____  _             _    _
        !  / ___|_ __ ___  ___ ___  / ___|| |_ __ _  ___| | _(_)_ __   __ _
        ! | |   | '__/ _ \/ __/ __| \___ \| __/ _` |/ __| |/ / | '_ \ / _` |
        ! | |___| | | (_) \__ \__ \  ___) | || (_| | (__|   <| | | | | (_| |
        !  \____|_|  \___/|___/___/ |____/ \__\__,_|\___|_|\_\_|_| |_|\__, |
        !                                                             |___/
        !

        ! ---------------------------
        ! Angle 1234: 1 -- 2 - 3 -- 4
        ! ---------------------------
        !
        cos_1243 = e21(1) * e43(1) + e21(2) * e43(2) + e21(3) * e43(3)
        if (cos_1243 >  1.0_wp) cos_1243 =  1.0_wp
        if (cos_1243 < -1.0_wp) cos_1243 = -1.0_wp
        sin_1243   = sqrt(1.0_wp - cos_1243 * cos_1243)
        angle_1243 = acos(cos_1243)
        if (sin_1243 < EPS) sin_1243 = EPS

        delta_angle_1243     = angle_1243 - bp_theta_3_0
        abs_delta_angle_1243 = abs(delta_angle_1243)

        ! ================================
        ! Cross-stacking between B2 and B5
        ! ================================
        !
        ! -----------------
        ! 1 -- 2    4 -- 3
        !       \
        !        \
        !         \
        !          \
        !           5
        ! -----------------
        !
        if ( ig_B5 > 0 .and. chain_id(ix_B5) == chain_id(ix_B4) ) then

          type_B5       = base_type           (ix_B5)
          cs_1_epsilon  = param_cs_1_epsilon  (type_B2, type_B5)
          cs_1_sigma    = param_cs_1_sigma    (type_B2, type_B5)
          cs_1_theta_cs_0 = param_cs_1_theta_cs (type_B2, type_B5)

          ! -------
          ! 2 ==> 5
          ! -------
          !
          pbc_B5(1:3) = coord_pbc(ix_B5,1:3) 

          d25(1:3) = pbc_B5(1:3) - pbc_B2(1:3) 
          r25_sqr  = d25(1) * d25(1) + d25(2) * d25(2) + d25(3) * d25(3)
          r25      = sqrt(r25_sqr)
          r25_inv  = 1.0_wp / r25
          e25(1:3) = d25(1:3) * r25_inv

          ! -----------------------
          ! Angle B2: 1 -- ~2~ -- 5
          ! -----------------------
          !
          cos_125 = e21(1) * e25(1) + e21(2) * e25(2) + e21(3) * e25(3)
          if (cos_125 >  1.0_wp) cos_125 =  1.0_wp
          if (cos_125 < -1.0_wp) cos_125 = -1.0_wp
          sin_125   = sqrt(1.0_wp - cos_125 * cos_125)
          angle_125 = acos(cos_125)
          if (sin_125 < EPS) sin_125 = EPS

          delta_angle_125     = angle_125 - cs_1_theta_cs_0
          abs_delta_angle_125 = abs(delta_angle_125)
  
          call calculate_attractive(d25, r25, cs_1_sigma, param_cs_alpha, &
                                    cs_1_epsilon, ene_attr, grad_attr)

          if ( abs_delta_angle_125 < cs_theta_threshold_1 ) then

            if ( abs_delta_angle_1243 <= bp_theta_threshold_1  ) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B2(1:3)  = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)  = force_B5(1:3) - force_attr_tmp(1:3)

            else if ( abs_delta_angle_1243 <= bp_theta_threshold_2  ) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, &
                                   grad_angle_1243)
              var_tmp = sin(param_bp_K * delta_angle_1243)
              ene_coef_theta_3  = var_tmp * var_tmp
              grad_coef_theta_3 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_1243)

              e_tmp_cs = e_tmp_cs + ene_coef_theta_3 * ene_attr

              force_coef_tmp            = grad_coef_theta_3 * ene_attr
              force_angle_1243_tmp(1:6) = - force_coef_tmp * grad_angle_1243(1:6)
              force_S1(1:3)             = force_S1(1:3) + force_angle_1243_tmp(1:3)
              force_B2(1:3)             = force_B2(1:3) - force_angle_1243_tmp(1:3)
              force_B4(1:3)             = force_B4(1:3) - force_angle_1243_tmp(4:6)
              force_S3(1:3)             = force_S3(1:3) + force_angle_1243_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B2(1:3)             = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)             = force_B5(1:3) - force_attr_tmp(1:3)

            end if

          else if ( abs_delta_angle_125 <= cs_theta_threshold_2 ) then

            call calculate_angle(d21, d25, r21, r25, cos_125, sin_125, grad_angle_125)
            var_tmp = sin(param_cs_K * delta_angle_125)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_125)

            if ( abs_delta_angle_1243 <= bp_theta_threshold_1  ) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_125_tmp(1:6) = - force_coef_tmp * grad_angle_125(1:6)
              force_S1(1:3)            = force_S1(1:3) + force_angle_125_tmp(1:3)
              force_B2(1:3)            = force_B2(1:3) - force_angle_125_tmp(1:3) &
                                       - force_angle_125_tmp(4:6)
              force_B5(1:3)            = force_B5(1:3) + force_angle_125_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)            = force_B5(1:3) - force_attr_tmp(1:3)
  
            else if ( abs_delta_angle_1243 <= bp_theta_threshold_2  ) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, &
                                   grad_angle_1243)
              var_tmp = sin(param_bp_K * delta_angle_1243)
              ene_coef_theta_3  = var_tmp * var_tmp
              grad_coef_theta_3 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_1243)
  
              e_tmp_cs = e_tmp_cs + ene_coef_theta_3 * ene_coef_theta_cs * ene_attr

              force_coef_tmp            = grad_coef_theta_3 * ene_coef_theta_cs * ene_attr
              force_angle_1243_tmp(1:6) = - force_coef_tmp * grad_angle_1243(1:6)
              force_S1(1:3)             = force_S1(1:3) + force_angle_1243_tmp(1:3)
              force_B2(1:3)             = force_B2(1:3) - force_angle_1243_tmp(1:3)
              force_B4(1:3)             = force_B4(1:3) - force_angle_1243_tmp(4:6)
              force_S3(1:3)             = force_S3(1:3) + force_angle_1243_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * grad_coef_theta_cs * ene_attr
              force_angle_125_tmp(1:6)  = - force_coef_tmp * grad_angle_125(1:6)
              force_S1(1:3)             = force_S1(1:3) + force_angle_125_tmp(1:3)
              force_B2(1:3)             = force_B2(1:3) - force_angle_125_tmp(1:3) &
                                        - force_angle_125_tmp(4:6)
              force_B5(1:3)             = force_B5(1:3) + force_angle_125_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * ene_coef_theta_cs
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B2(1:3)             = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)             = force_B5(1:3) - force_attr_tmp(1:3)
  
            end if
  
          end if

        end if

        ! ================================
        ! Cross-stacking between B4 and B6
        ! ================================
        !
        ! -----------------
        ! 1 -- 2    4 -- 3
        !          /
        !         /
        !        /
        !       /
        !      6
        ! -----------------
        !
        if ( ig_B6 <= natom_all .and. &
             chain_id(ix_B6) == chain_id(ix_B2) ) then

          type_B6       = base_type           (ix_B6)
          cs_2_epsilon  = param_cs_2_epsilon  (type_B4, type_B6)
          cs_2_sigma    = param_cs_2_sigma    (type_B4, type_B6)
          cs_2_theta_cs_0 = param_cs_2_theta_cs (type_B4, type_B6)

          ! -------
          ! 6 <== 4
          ! -------
          !
          !
          pbc_B6(1:3) = coord_pbc(ix_B6,1:3)

          d46(1:3) = pbc_B6(1:3) - pbc_B4(1:3) 
          r46_sqr  = d46(1) * d46(1) + d46(2) * d46(2) + d46(3) * d46(3)
          r46      = sqrt(r46_sqr)
          r46_inv  = 1.0_wp / r46
          e46(1:3) = d46(1:3) * r46_inv

          ! -----------------------
          ! Angle B2: 6 -- ~4~ -- 3
          ! -----------------------
          !
          cos_346 = e43(1) * e46(1) + e43(2) * e46(2) + e43(3) * e46(3)
          if (cos_346 >  1.0_wp) cos_346 =  1.0_wp
          if (cos_346 < -1.0_wp) cos_346 = -1.0_wp
          sin_346   = sqrt(1.0_wp - cos_346 * cos_346)
          angle_346 = acos(cos_346)
          if (sin_346 < EPS) sin_346 = EPS

          delta_angle_346     = angle_346 - cs_2_theta_cs_0
          abs_delta_angle_346 = abs(delta_angle_346)

          call calculate_attractive(d46, r46, cs_2_sigma, param_cs_alpha, &
                                    cs_2_epsilon, ene_attr, grad_attr)

          if ( abs_delta_angle_346 < cs_theta_threshold_1 ) then

            if ( abs_delta_angle_1243 <= bp_theta_threshold_1  ) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B4(1:3)  = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)  = force_B6(1:3) - force_attr_tmp(1:3)

            else if ( abs_delta_angle_1243 <= bp_theta_threshold_2  ) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, &
                                   grad_angle_1243)
              var_tmp = sin(param_bp_K * delta_angle_1243)
              ene_coef_theta_3  = var_tmp * var_tmp
              grad_coef_theta_3 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_1243)

              e_tmp_cs = e_tmp_cs + ene_coef_theta_3 * ene_attr

              force_coef_tmp            = grad_coef_theta_3 * ene_attr
              force_angle_1243_tmp(1:6) = - force_coef_tmp * grad_angle_1243(1:6)
              force_S1(1:3)             = force_S1(1:3) + force_angle_1243_tmp(1:3)
              force_B2(1:3)             = force_B2(1:3) - force_angle_1243_tmp(1:3)
              force_B4(1:3)             = force_B4(1:3) - force_angle_1243_tmp(4:6)
              force_S3(1:3)             = force_S3(1:3) + force_angle_1243_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B4(1:3)             = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)             = force_B6(1:3) - force_attr_tmp(1:3)

            end if

          else if ( abs_delta_angle_346 <= cs_theta_threshold_2 ) then

            call calculate_angle(d43, d46, r43, r46, cos_346, sin_346, grad_angle_346)
            var_tmp = sin(param_cs_K * delta_angle_346)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_346)

            if ( abs_delta_angle_1243 <= bp_theta_threshold_1  ) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_346_tmp(1:6) = - force_coef_tmp * grad_angle_346(1:6)
              force_S3(1:3)            = force_S3(1:3) + force_angle_346_tmp(1:3)
              force_B4(1:3)            = force_B4(1:3) - force_angle_346_tmp(1:3) &
                                       - force_angle_346_tmp(4:6)
              force_B6(1:3)            = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B4(1:3)            = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)            = force_B6(1:3) - force_attr_tmp(1:3)

            else if ( abs_delta_angle_1243 <= bp_theta_threshold_2  ) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, &
                                   grad_angle_1243)
              var_tmp = sin(param_bp_K * delta_angle_1243)
              ene_coef_theta_3  = var_tmp * var_tmp
              grad_coef_theta_3 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_1243)

              e_tmp_cs = e_tmp_cs + ene_coef_theta_3 * ene_coef_theta_cs * ene_attr

              force_coef_tmp            = grad_coef_theta_3 * ene_coef_theta_cs * ene_attr
              force_angle_1243_tmp(1:6) = - force_coef_tmp * grad_angle_1243(1:6)
              force_S1(1:3)             = force_S1(1:3) + force_angle_1243_tmp(1:3)
              force_B2(1:3)             = force_B2(1:3) - force_angle_1243_tmp(1:3)
              force_B4(1:3)             = force_B4(1:3) - force_angle_1243_tmp(4:6)
              force_S3(1:3)             = force_S3(1:3) + force_angle_1243_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * grad_coef_theta_cs * ene_attr
              force_angle_346_tmp(1:6)  = - force_coef_tmp * grad_angle_346(1:6)
              force_S3(1:3)             = force_S3(1:3) + force_angle_346_tmp(1:3)
              force_B4(1:3)             = force_B4(1:3) - force_angle_346_tmp(1:3) &
                                        - force_angle_346_tmp(4:6)
              force_B6(1:3)             = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * ene_coef_theta_cs
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B4(1:3)             = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)             = force_B6(1:3) - force_attr_tmp(1:3)

            end if

          end if

        end if

        ! Calc energy
        epair(id+1) = epair(id+1) + e_tmp_bp + e_tmp_cs

        ! store force
        !
        force(ix_S1,1:3,id+1) = force(ix_S1,1:3,id+1) + force_S1(1:3)
        force(ix_B2,1:3,id+1) = force(ix_B2,1:3,id+1) + force_B2(1:3)
        force(ix_B4,1:3,id+1) = force(ix_B4,1:3,id+1) + force_B4(1:3)
        force(ix_S3,1:3,id+1) = force(ix_S3,1:3,id+1) + force_S3(1:3)
        force(ix_B5,1:3,id+1) = force(ix_B5,1:3,id+1) + force_B5(1:3)
        force(ix_B6,1:3,id+1) = force(ix_B6,1:3,id+1) + force_B6(1:3)

        ! virial
        !
        viri(1:3) = force_S1(1:3) * coord_pbc(ix_S1,1:3) &
                  + force_B2(1:3) * coord_pbc(ix_B2,1:3) &
                  + force_B4(1:3) * coord_pbc(ix_B4,1:3) &
                  + force_S3(1:3) * coord_pbc(ix_S3,1:3) &
                  + force_B5(1:3) * coord_pbc(ix_B5,1:3) &
                  + force_B6(1:3) * coord_pbc(ix_B6,1:3)

        virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

      end do

    end do
    !$omp end parallel

    call timer(TimerBasePair, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_base_pairing


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_ele
  !> @brief        Calculates electrostatic energy/force 
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_ele(domain, enefunc, pairlist, &
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
    real(wp)                  :: dij(3), rij, rij2
    real(wp)                  :: inv_rij
    real(wp)                  :: debye_length
    real(wp)                  :: inv_debye_length
    real(wp)                  :: ele_coef
    real(wp)                  :: diele_const
    real(wp)                  :: inv_diele_const
    real(wp)                  :: ele_tmp_e_T, ele_tmp_a_C
    real(wp)                  :: ele_tmp_sol_T, ele_tmp_sol_C
    real(wp)                  :: scale_Q(3,3)
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: term_elec, grad_coef, work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, force_local(3), factor
    integer                   :: i, ix, iy, j, k, start_i, ilist
    integer                   :: type1, type2
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: ncell

    integer,          pointer :: cg_elec_list(:)
    integer,          pointer :: start_atom(:)
    integer(1),       pointer :: charge_type(:)
    integer,          pointer :: num_nb15_nobc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    real(wp),         pointer :: charge(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    start_atom      => domain%start_atom
    charge_type     => domain%charge_type
    charge          => domain%charge

    cg_elec_list    => enefunc%cg_elec_list
    num_nb15_nobc   => pairlist%num_cg_ele_calc
    nb15_calc_list  => pairlist%cg_ele_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff          =  enefunc%cg_cutoffdist_ele
    cutoff2         =  cutoff*cutoff

    ! electrostatic parameters
    ele_tmp_sol_T = enefunc%cg_ele_sol_T
    ele_tmp_sol_C = enefunc%cg_ele_sol_IC
    ele_tmp_e_T   = 2.494e2_wp - 7.88e-1_wp * ele_tmp_sol_T &
        + 7.2e-4_wp * ele_tmp_sol_T * ele_tmp_sol_T
    ele_tmp_a_C   = 1.0e0_wp - 2.551e-1_wp * ele_tmp_sol_C  &
        + 5.151e-2_wp * ele_tmp_sol_C * ele_tmp_sol_C       &
        - 6.889e-3_wp * ele_tmp_sol_C * ele_tmp_sol_C * ele_tmp_sol_C
    diele_const = ele_tmp_e_T * ele_tmp_a_C

    debye_length =              &
         sqrt(                                     &
        (CAL2JOU * ELECTRIC_CONST1                  &
        * diele_const * KBOLTZ * ele_tmp_sol_T) &
        /                                           &
        (2.0_dp * AVOGADRO1 * AVOGADRO1               &
        * ELEMENT_CHARGE1 * ELEMENT_CHARGE1 * ele_tmp_sol_C)  &
        )

    inv_debye_length = 1.0_wp / debye_length
    ele_coef         = enefunc%cg_ele_coef
    inv_diele_const  = 1.0_wp / diele_const

    scale_Q(1:3,1:3) = 1.0_wp
    scale_Q(1  ,2  ) = - enefunc%cg_pro_DNA_ele_scale_Q / 0.6
    scale_Q(2  ,1  ) = scale_Q(1,2)

    !$omp parallel default(shared)                                   &
    !$omp private(id, i, ix, j, iy, k, rtmp, type1, type2, num_nb15, &
    !$omp         dij, rij2,  inv_rij,  rij, elec_temp, term_elec,   &
    !$omp         grad_coef, work, force_local, start_i, factor,     &
    !$omp         qtmp, ilist)
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
      type1   = charge_type(ix)
      qtmp    = charge(ix)

      num_nb15 = num_nb15_nobc(i)
      elec_temp   = 0.0_wp
      force_local(1:3) = 0.0_wp

      !!dir$ simd
      !ocl norecurrence
      do k = 1, num_nb15

        j  = nb15_calc_list(k,i)

        ! compute distance
        !
        type2     = charge_type(j)
        dij(1) = rtmp(1) - coord(j,1)
        dij(2) = rtmp(2) - coord(j,2)
        dij(3) = rtmp(3) - coord(j,3)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        factor = merge(0.0_wp,1.0_wp,rij2>cutoff2)

!       if (rij2 > cutoff2) cycle

        rij = sqrt(rij2)
        inv_rij = 1.0_wp / rij

        term_elec = ele_coef * qtmp * charge(j)  &
                  * inv_diele_const * inv_rij    &
                  * exp(-rij*inv_debye_length)
        term_elec = term_elec * scale_Q(type1,type2) * factor

        elec_temp = elec_temp + term_elec

        grad_coef = -term_elec*(inv_debye_length+inv_rij)*inv_rij

        work(1)   = grad_coef * dij(1)
        work(2)   = grad_coef * dij(2)
        work(3)   = grad_coef * dij(3)

        ! store force
          !
        force_local(1) = force_local(1) - work(1)
        force_local(2) = force_local(2) - work(2)
        force_local(3) = force_local(3) - work(3)
        force(j,1,id+1)  = force(j,1,id+1) + work(1)
        force(j,2,id+1)  = force(j,2,id+1) + work(2)
        force(j,3,id+1)  = force(j,3,id+1) + work(3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
        virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
        virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

      end do

      force(ix,1,id+1) = force(ix,1,id+1) + force_local(1)
      force(ix,2,id+1) = force(ix,2,id+1) + force_local(2)
      force(ix,3,id+1) = force(ix,3,id+1) + force_local(3)
      elec(id+1) = elec(id+1) + elec_temp

    end do

    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_PWMcos
  !> @brief        Calculates PWMcos energy with pairlist
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_PWMcos(domain, enefunc, pairlist, &
                                   coord, force, epwmcos, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: epwmcos(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: trans(3)
    real(wp)                  :: sigma, sig_sqr, phi, phi2
    real(wp)                  :: rtmp(3), viri(3)
    real(wp)                  :: epwmcos_temp
    integer                   :: i, k, ix, i1, j, iy, icount, start_i
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: ncell

    ! indicex + params
    integer                   :: i_CA_0, i_CA_N, i_CA_C
    integer                   :: i_CA_N_cell, i_CA_N_atom
    integer                   :: i_CA_C_cell, i_CA_C_atom
    integer                   :: i_DB_0, i_DB_5, i_DB_3
    integer                   :: i_DS_0
    integer                   :: i_DB_5_cell, i_DB_5_atom
    integer                   :: i_DB_3_cell, i_DB_3_atom
    integer                   :: i_DS_0_cell, i_DS_0_atom
    integer                   :: j_base_type
    real(wp)                  :: pwm_shift, pwm_factor
    real(wp)                  :: cutoff, cutoff_sqr
    ! native values
    real(wp)                  :: r00, t10, t20, t30
    ! vectors and distances
    real(wp)                  :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)                  :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)                  :: v53(3), r53, r53_sqr, r53_inv, e53(3)
    real(wp)                  :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)                  :: t1, t2, t3
    real(wp)                  :: cos_t1, cos_t2, cos_t3
    real(wp)                  :: sin_t1, sin_t2, sin_t3
    ! modulating functions
    real(wp)                  :: dt1, dt2, dt3
    real(wp)                  :: abs_dt1, abs_dt2, abs_dt3
    real(wp)                  :: cos_dt1, cos_dt2, cos_dt3
    real(wp)                  :: sin_dt1, sin_dt2, sin_dt3
    real(wp)                  :: ktheta_2, ktheta
    ! energy calc
    real(wp)                  :: pwm_score
    real(wp)                  :: pwm_score_A, pwm_score_C
    real(wp)                  :: pwm_score_G, pwm_score_T
    real(wp)                  :: ene_coef_r0
    real(wp)                  :: grad_coef_t1, ene_coef_t1
    real(wp)                  :: grad_coef_t2, ene_coef_t2
    real(wp)                  :: grad_coef_t3, ene_coef_t3
    real(wp)                  :: etmp0
    real(wp)                  :: e_tmp_pwmcos
    ! force coef
    real(wp)                  :: var_tmp
    real(wp)                  :: grad_r0(3)
    real(wp)                  :: grad_t1(6)
    real(wp)                  :: grad_t2(6)
    real(wp)                  :: grad_t3(6)
    ! force calc
    real(wp)                  :: f_r0_coef_tmp
    real(wp)                  :: f_r0_tmp(6)
    real(wp)                  :: f_t1_coef_tmp
    real(wp)                  :: f_t1_tmp(6)
    real(wp)                  :: f_t2_coef_tmp
    real(wp)                  :: f_t2_tmp(6)
    real(wp)                  :: f_t3_coef_tmp
    real(wp)                  :: f_t3_tmp(6)
    real(wp)                  :: f_CA_0(3)
    real(wp)                  :: f_CA_N(3)
    real(wp)                  :: f_CA_C(3)
    real(wp)                  :: f_DB_0(3)
    real(wp)                  :: f_DB_5(3)
    real(wp)                  :: f_DB_3(3)
    real(wp)                  :: f_DS_0(3)
    real(wp)                  :: force_CA_0(3)
    real(wp)                  :: force_CA_N(3)
    real(wp)                  :: force_CA_C(3)

    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: pwmcos_id(:)
    integer,          pointer :: num_nb15_calc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    integer,          pointer :: pwmcos_id_N(:)
    integer,          pointer :: pwmcos_id_C(:)
    integer,          pointer :: base_type(:)
    real(wp),         pointer :: pwmcos_r0(:,:)
    real(wp),         pointer :: pwmcos_theta1(:,:)
    real(wp),         pointer :: pwmcos_theta2(:,:)
    real(wp),         pointer :: pwmcos_theta3(:,:)
    real(wp),         pointer :: pwmcos_ene_A(:,:)
    real(wp),         pointer :: pwmcos_ene_C(:,:)
    real(wp),         pointer :: pwmcos_ene_G(:,:)
    real(wp),         pointer :: pwmcos_ene_T(:,:)
    real(wp),         pointer :: pwmcos_gamma(:,:)
    real(wp),         pointer :: pwmcos_eps(:,:)
    integer,          pointer :: pwmcos_count(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPWMcos, TimerOn)

    id_l2g             => domain%id_l2g
    id_g2l             => domain%id_g2l
    start_atom         => domain%start_atom
    base_type          => domain%na_base_type

    pwmcos_id          => enefunc%pwmcos_protein_id
    pwmcos_id_N        => enefunc%pwmcos_protein_id_N
    pwmcos_id_C        => enefunc%pwmcos_protein_id_C
    pwmcos_r0          => enefunc%pwmcos_r0
    pwmcos_theta1      => enefunc%pwmcos_theta1
    pwmcos_theta2      => enefunc%pwmcos_theta2
    pwmcos_theta3      => enefunc%pwmcos_theta3
    pwmcos_ene_A       => enefunc%pwmcos_ene_A
    pwmcos_ene_C       => enefunc%pwmcos_ene_C
    pwmcos_ene_G       => enefunc%pwmcos_ene_G
    pwmcos_ene_T       => enefunc%pwmcos_ene_T
    pwmcos_gamma       => enefunc%pwmcos_gamma
    pwmcos_eps         => enefunc%pwmcos_eps
    pwmcos_count       => enefunc%pwmcos_count

    num_nb15_calc      => pairlist%num_cg_pwmcos_calc
    nb15_calc_list     => pairlist%cg_pwmcos_list

    ncell              =  domain%num_cell_local 

    ! set parameters
    !
    sigma     = enefunc%pwmcos_sigma
    sig_sqr   = sigma * sigma
    !
    phi       = enefunc%pwmcos_phi
    phi2      = phi * 2.0_wp
    !
    ktheta_2  = PI / phi
    ktheta    = ktheta_2 / 2.0_wp

    !$omp parallel default(shared)                                        &
    !$omp private(id, i, i1, ix, iy, j, k, i_CA_N_cell, i_CA_N_atom,      &
    !$omp         i_CA_C_cell, i_CA_C_atom, i_DB_5_cell, i_DB_5_atom,     &
    !$omp         i_DB_3_cell, i_DB_3_atom, i_DS_0_cell, i_DS_0_atom,     &
    !$omp         rtmp, force_CA_0, force_CA_C, force_CA_N, epwmcos_temp, &
    !$omp         num_nb15, trans, i_CA_0, i_CA_N, i_CA_C, i_DB_0,        &
    !$omp         i_DB_5, i_DB_3, i_DS_0, j_base_type, pwm_shift,         &
    !$omp         pwm_factor, cutoff, cutoff_sqr, r00, t10, t20, t30,     &
    !$omp         vbc, vbs, v53, vcn, rbc, rbs, r53, rcn, rbc_sqr,        &
    !$omp         rbs_sqr, r53_sqr, rcn_sqr, rbc_inv, rbs_inv, r53_inv,   &
    !$omp         rcn_inv, ebc, ebs, e53, ecn, t1, t2, t3, cos_t1,        &
    !$omp         cos_t2, cos_t3, sin_t1, sin_t2, sin_t3, dt1, dt2, dt3,  &
    !$omp         abs_dt1, abs_dt2, abs_dt3, cos_dt1, cos_dt2, cos_dt3,   &
    !$omp         sin_dt1, sin_dt2, sin_dt3, pwm_score_A, pwm_score_C,    &
    !$omp         pwm_score_G, pwm_score_T, ene_coef_r0,                  &
    !$omp         grad_coef_t1, ene_coef_t1, grad_coef_t2, ene_coef_t2,   &
    !$omp         grad_coef_t3, ene_coef_t3, etmp0, e_tmp_pwmcos,         &
    !$omp         var_tmp, grad_r0, grad_t1, grad_t2, grad_t3,            &
    !$omp         f_r0_coef_tmp, f_r0_tmp, f_t1_coef_tmp, f_t1_tmp,       &
    !$omp         f_t2_coef_tmp, f_t2_tmp, f_t3_coef_tmp, f_t3_tmp,       &
    !$omp         f_CA_0, f_CA_N, f_CA_C,  f_DB_0, f_DB_5, f_DB_3,        &
    !$omp         f_DS_0, start_i, icount, pwm_score, viri) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, enefunc%num_pwmcos_domain, nthread

      ix = id_g2l(pwmcos_id(i))
      rtmp(1:3) = coord(ix,1:3)
 
      i_CA_C      = pwmcos_id_C(i)
      i_CA_C_atom = id_g2l(i_CA_C) 
      i_CA_N      = pwmcos_id_N(i)
      i_CA_N_atom = id_g2l(i_CA_N)
      vcn(1:3)    = coord(i_CA_N_atom,1:3) &
                  - coord(i_CA_C_atom,1:3) 
      rcn_sqr     = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
      rcn         = sqrt(rcn_sqr)
      rcn_inv     = 1.0_wp / rcn
      ecn(1:3)    = vcn(1:3) * rcn_inv

      num_nb15 = num_nb15_calc(i)
      epwmcos_temp = 0.0_wp
      force_CA_0(1:3) = 0.0_wp
      force_CA_N(1:3) = 0.0_wp
      force_CA_C(1:3) = 0.0_wp

      do icount = 1, pwmcos_count(i)

        r00         = pwmcos_r0    (icount, i)
        t10         = pwmcos_theta1(icount, i)
        t20         = pwmcos_theta2(icount, i)
        t30         = pwmcos_theta3(icount, i)
        pwm_score_A = pwmcos_ene_A (icount, i)
        pwm_score_C = pwmcos_ene_C (icount, i)
        pwm_score_G = pwmcos_ene_G (icount, i)
        pwm_score_T = pwmcos_ene_T (icount, i)
        pwm_shift   = pwmcos_eps   (icount, i)
        pwm_factor  = pwmcos_gamma (icount, i)
        cutoff      = r00 + 5.0_wp
        cutoff_sqr  = cutoff * cutoff

        do k = 1, num_nb15

          j  = nb15_calc_list(k,i)

          vbc(1:3) = coord(ix,1:3) - coord(j,1:3) 
          rbc_sqr = vbc(1)*vbc(1) + vbc(2)*vbc(2) + vbc(3)*vbc(3)

          if (rbc_sqr >= cutoff_sqr) cycle

          f_CA_0(1:3) = 0.0_wp
          f_CA_N(1:3) = 0.0_wp
          f_CA_C(1:3) = 0.0_wp
          f_DB_0(1:3) = 0.0_wp
          f_DB_5(1:3) = 0.0_wp
          f_DB_3(1:3) = 0.0_wp
          f_DS_0(1:3) = 0.0_wp

          e_tmp_pwmcos = 0.0_wp

          ! ---------------------------------
          ! PWM energy: sequence specificity!
          ! ---------------------------------
          !
          j_base_type = base_type(j)
          if (j_base_type == NABaseTypeDBA) pwm_score = pwm_score_A
          if (j_base_type == NABaseTypeDBC) pwm_score = pwm_score_C
          if (j_base_type == NABaseTypeDBG) pwm_score = pwm_score_G
          if (j_base_type == NABaseTypeDBT) pwm_score = pwm_score_T
          etmp0 = pwm_factor * ( pwm_score + pwm_shift )
            
          ! ==================================
          ! Distance and angle calculations...
          ! ==================================
          i_DB_0 = id_l2g(j)
          i_DB_5 = i_DB_0 - 3
          i_DB_3 = i_DB_0 + 3
          i_DS_0 = i_DB_0 - 1

          ! --------------
          ! sugar <== base
          ! --------------
          !
          i_DS_0_atom = id_g2l(i_DS_0)
          vbs(1:3)    = coord(i_DS_0_atom,1:3) - coord(j,1:3)  
          rbs_sqr     = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
          rbs         = sqrt(rbs_sqr)
          rbs_inv     = 1.0_wp / rbs
          ebs(1:3)    = vbs(1:3) * rbs_inv

          ! -----------
          ! base ==> Ca
          ! -----------
          !
          rbc      = sqrt(rbc_sqr)
          rbc_inv  = 1.0_wp / rbc
          ebc(1:3) = vbc(1:3) * rbc_inv

          ! -------------------
          ! base 5' ==> base 3'
          ! -------------------
          !
          i_DB_5_atom = id_g2l(i_DB_5) 
          i_DB_3_atom = id_g2l(i_DB_3) 
          v53(1:3) = coord(i_DB_3_atom,1:3) - coord(i_DB_5_atom,1:3) 
          r53_sqr  = v53(1) * v53(1) + v53(2) * v53(2) + v53(3) * v53(3)
          r53      = sqrt(r53_sqr)
          r53_inv  = 1.0_wp / r53
          e53(1:3) = v53(1:3) * r53_inv

          ! -----------------------------
          ! Angle t1: sugar -- base -- Ca
          ! -----------------------------
          !
          cos_t1 = ebc(1) * ebs(1) + ebc(2) * ebs(2) + ebc(3) * ebs(3)
          if (cos_t1 >  1.0_wp) cos_t1 =  1.0_wp
          if (cos_t1 < -1.0_wp) cos_t1 = -1.0_wp
          sin_t1 = sqrt(1.0_wp - cos_t1 * cos_t1)
          if (sin_t1 < EPS) sin_t1 = EPS
          t1     = acos(cos_t1)

          ! ------------------------------------------
          ! Angle t2: base5' - base3' <==> base0 -- Ca
          ! ------------------------------------------
          !
          cos_t2 = ebc(1) * e53(1) + ebc(2) * e53(2) + ebc(3) * e53(3)
          if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
          if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
          sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
          if (sin_t2 < EPS) sin_t2 = EPS
          t2     = acos(cos_t2)

          ! ---------------------------------------
          ! Angle t3: base0 -- Ca <==> Ca_N -- Ca_C
          ! ---------------------------------------
          !
          cos_t3 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
          if (cos_t3 >  1.0_wp) cos_t3 =  1.0_wp
          if (cos_t3 < -1.0_wp) cos_t3 = -1.0_wp
          sin_t3 = sqrt(1.0_wp - cos_t3 * cos_t3)
          if (sin_t3 < EPS) sin_t3 = EPS
          t3     = acos(cos_t3)

          ! ---------------
          ! Delta angles...
          ! ---------------
          !
          dt1     = t1 - t10
          dt2     = t2 - t20
          dt3     = t3 - t30
          abs_dt1 = abs(dt1)
          abs_dt2 = abs(dt2)
          abs_dt3 = abs(dt3)

          ! ==================================================================
          ! Energy/Force Calculation
          ! ==================================================================
          !
          ! -----------------
          ! Simple angle test
          ! -----------------
          !
          if (abs_dt1 >= phi2 .or. &
              abs_dt2 >= phi2 .or. &
              abs_dt3 >= phi2 ) then
            cycle
          end if

          ! -------------------
          ! Distance modulating
          ! -------------------
          !
          call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, &
                                  grad_r0)

          ! ----------------
          ! Angle modulating
          ! ----------------
          if ( abs_dt1 < phi ) then
            ene_coef_t1  = 1.0_wp
            grad_coef_t1 = 0.0_wp
            grad_t1(1:6) = 0.0_wp
          else
            sin_dt1      = sin(ktheta * dt1)
            ene_coef_t1  = sin_dt1 * sin_dt1
            call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
            grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
          end if
          if ( abs_dt2 < phi ) then
            ene_coef_t2  = 1.0_wp
            grad_coef_t2 = 0.0_wp
            grad_t2(1:6) = 0.0_wp
          else
            sin_dt2      = sin(ktheta * dt2)
            ene_coef_t2  = sin_dt2 * sin_dt2
            call calculate_angle(vbc, v53, rbc, r53, cos_t2, sin_t2, grad_t2)
            grad_coef_t2 = ktheta * sin(2.0_wp * ktheta * dt2)
          end if
          if ( abs_dt3 < phi ) then
            ene_coef_t3  = 1.0_wp
            grad_coef_t3 = 0.0_wp
            grad_t3(1:6) = 0.0_wp
          else
            sin_dt3      = sin(ktheta * dt3)
            ene_coef_t3  = sin_dt3 * sin_dt3
            call calculate_angle(vbc, vcn, rbc, rcn, cos_t3, sin_t3, grad_t3)
            grad_coef_t3 = ktheta * sin(2.0_wp * ktheta * dt3)
          end if

          ! ======
          ! Energy
          ! ======
          !
          e_tmp_pwmcos = etmp0 * ene_coef_r0 * ene_coef_t1 &
                       * ene_coef_t2 * ene_coef_t3
          epwmcos_temp = epwmcos_temp + e_tmp_pwmcos

          ! ======
          ! Forces
          ! ======
          !
          ! --
          ! r0
          ! --
          !
          f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2 * ene_coef_t3
          f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_r0_tmp(1:3)
          !
          ! --
          ! t1
          ! --
          !
          f_t1_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t2 * ene_coef_t3 * grad_coef_t1
          f_t1_tmp(1:6) = - f_t1_coef_tmp * grad_t1(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t1_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_t1_tmp(1:3) - f_t1_tmp(4:6)
          f_DS_0(1:3)   = f_DS_0(1:3) + f_t1_tmp(4:6)
          !
          ! --
          ! t2
          ! --
          !
          f_t2_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t3 * grad_coef_t2
          f_t2_tmp(1:6) = - f_t2_coef_tmp * grad_t2(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t2_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_t2_tmp(1:3)
          f_DB_5(1:3)   = f_DB_5(1:3) - f_t2_tmp(4:6)
          f_DB_3(1:3)   = f_DB_3(1:3) + f_t2_tmp(4:6)
          !
          ! --
          ! t3
          ! --
          !
          f_t3_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2 * grad_coef_t3
          f_t3_tmp(1:6) = - f_t3_coef_tmp * grad_t3(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t3_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_t3_tmp(1:3)
          f_CA_C(1:3)   = f_CA_C(1:3) - f_t3_tmp(4:6)
          f_CA_N(1:3)   = f_CA_N(1:3) + f_t3_tmp(4:6)

          ! ------------------------
          ! update all the forces...
          ! ------------------------
          !
          force_CA_0(1:3) = force_CA_0(1:3) + f_CA_0(1:3)
          force_CA_N(1:3) = force_CA_N(1:3) + f_CA_N(1:3)
          force_CA_C(1:3) = force_CA_C(1:3) + f_CA_C(1:3)
          force(j,1:3,id+1) = force(j,1:3,id+1) + f_DB_0(1:3)
          force(i_DB_5_atom,1:3,id+1) = &
                force(i_DB_5_atom,1:3,id+1) + f_DB_5(1:3)
          force(i_DB_3_atom,1:3,id+1) = &
                force(i_DB_3_atom,1:3,id+1) + f_DB_3(1:3)
          force(i_DS_0_atom,1:3,id+1) = &
                force(i_DS_0_atom,1:3,id+1) + f_DS_0(1:3)

          viri(1:3) = f_CA_0(1:3) * coord(ix,1:3)          &
                    + f_CA_N(1:3) * coord(i_CA_N_atom,1:3) &
                    + f_CA_C(1:3) * coord(i_cA_C_atom,1:3) &
                    + f_DB_0(1:3) * coord(j,1:3)           &
                    + f_DB_5(1:3) * coord(i_DB_5_atom,1:3) &
                    + f_DB_3(1:3) * coord(i_DB_3_atom,1:3) &
                    + f_DS_0(1:3) * coord(i_DS_0_atom,1:3)
          virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
          virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
          virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

        end do

      end do

      force(ix,1:3,id+1) = force(ix,1:3,id+1) + force_CA_0(1:3)
      force(i_CA_C_atom,1:3,id+1) = &
          force(i_CA_C_atom,1:3,id+1) + force_CA_C(1:3)
      force(i_CA_N_atom,1:3,id+1) = &
          force(i_CA_N_atom,1:3,id+1) + force_CA_N(1:3)
      epwmcos(id+1) = epwmcos(id+1) + epwmcos_temp

    end do

    !$omp end parallel

    call timer(TimerCGPWMcos, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_PWMcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_PWMcosns
  !> @brief        Calculates PWMcosns energy with pairlist 
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_PWMcosns(domain, enefunc, pairlist, &
                                     coord, force, epwmcosns, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(dp),                 intent(inout) :: epwmcosns(nthread)
    real(dp),                 intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                  :: trans(3)
    real(wp)                  :: sigma, sig_sqr, phi, phi2
    real(wp)                  :: rtmp(3), viri(3)
    real(wp)                  :: epwmcosns_temp
    integer                   :: i, ix, i1, j, iy, k, icount, start_i
    integer                   :: num_nb15
    integer                   :: id
    integer                   :: ncell

    ! indicex + params
    integer                   :: i_CA_N_atom, i_CA_C_atom
    integer                   :: i_DS_0_atom
    real(wp)                  :: cutoff, cutoff_sqr
    ! native values
    real(wp)                  :: r00, t10, t20
    ! vectors and distances
    real(wp)                  :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)                  :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)                  :: v53(3), r53, r53_sqr, r53_inv, e53(3)
    real(wp)                  :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)                  :: t1, t2
    real(wp)                  :: cos_t1, cos_t2
    real(wp)                  :: sin_t1, sin_t2
    ! modulating functions
    real(wp)                  :: dt1, dt2
    real(wp)                  :: abs_dt1, abs_dt2
    real(wp)                  :: cos_dt1, cos_dt2
    real(wp)                  :: sin_dt1, sin_dt2
    real(wp)                  :: ktheta_2, ktheta
    ! energy calc
    real(wp)                  :: pwm_score_A, pwm_score_C
    real(wp)                  :: pwm_score_G, pwm_score_T
    real(wp)                  :: ene_coef_r0
    real(wp)                  :: grad_coef_t1, ene_coef_t1
    real(wp)                  :: grad_coef_t2, ene_coef_t2
    real(wp)                  :: etmp0
    real(wp)                  :: e_tmp_pwmcosns
    ! force coef
    real(wp)                  :: var_tmp
    real(wp)                  :: grad_r0(3)
    real(wp)                  :: grad_t1(6)
    real(wp)                  :: grad_t2(6)
    ! force calc
    real(wp)                  :: f_r0_coef_tmp
    real(wp)                  :: f_r0_tmp(6)
    real(wp)                  :: f_t1_coef_tmp
    real(wp)                  :: f_t1_tmp(6)
    real(wp)                  :: f_t2_coef_tmp
    real(wp)                  :: f_t2_tmp(6)
    real(wp)                  :: f_CA_0(3)
    real(wp)                  :: f_CA_N(3)
    real(wp)                  :: f_CA_C(3)
    real(wp)                  :: f_DP_0(3)
    real(wp)                  :: f_DS_0(3)
    real(wp)                  :: force_CA_0(3)
    real(wp)                  :: force_CA_N(3)
    real(wp)                  :: force_CA_C(3)

    integer,          pointer :: id_l2g(:), id_g2l(:)
    integer,          pointer :: start_atom(:)
    integer,          pointer :: pwmcosns_id(:)
    integer,          pointer :: num_nb15_nobc(:)
    integer,          pointer :: nb15_calc_list(:,:)
    integer,          pointer :: pwmcosns_id_N(:)
    integer,          pointer :: pwmcosns_id_C(:)
    integer,          pointer :: base_type(:)
    real(wp),         pointer :: pwmcosns_r0(:,:)
    real(wp),         pointer :: pwmcosns_theta1(:,:)
    real(wp),         pointer :: pwmcosns_theta2(:,:)
    real(wp),         pointer :: pwmcosns_ene(:,:)
    integer,          pointer :: pwmcosns_count(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPwmcosns, TimerOn)

    id_l2g             => domain%id_l2g
    id_g2l             => domain%id_g2l
    start_atom         => domain%start_atom
    base_type          => domain%na_base_type

    pwmcosns_id        => enefunc%pwmcosns_protein_id
    pwmcosns_id_N      => enefunc%pwmcosns_protein_id_N
    pwmcosns_id_C      => enefunc%pwmcosns_protein_id_C
    pwmcosns_r0        => enefunc%pwmcosns_r0
    pwmcosns_theta1    => enefunc%pwmcosns_theta1
    pwmcosns_theta2    => enefunc%pwmcosns_theta2
    pwmcosns_ene       => enefunc%pwmcosns_ene
    pwmcosns_count     => enefunc%pwmcosns_count

    num_nb15_nobc      => pairlist%num_cg_pwmcosns_calc
    nb15_calc_list     => pairlist%cg_pwmcosns_list

    ncell              =  domain%num_cell_local 

    ! set parameters
    !
    sigma     = enefunc%pwmcosns_sigma
    sig_sqr   = sigma * sigma
    !
    phi       = enefunc%pwmcosns_phi
    phi2      = phi * 2.0_wp
    !
    ktheta_2  = PI / phi
    ktheta    = ktheta_2 / 2.0_wp

    !$omp parallel default(shared)                                        &
    !$omp private(id, i, i1, ix, iy, j, k, i_CA_N_atom, i_CA_C_atom,      &
    !$omp         i_DS_0_atom, rtmp, force_CA_0, force_CA_C, force_CA_N,  &
    !$omp         epwmcosns_temp, num_nb15, trans,         &
    !$omp         cutoff, cutoff_sqr, r00, t10, t20,          &
    !$omp         vbc, vbs, v53, vcn, rbc, rbs, r53, rcn, rbc_sqr,        &
    !$omp         rbs_sqr, r53_sqr, rcn_sqr, rbc_inv, rbs_inv, r53_inv,   &
    !$omp         rcn_inv, ebc, ebs, e53, ecn, t1, t2, cos_t1, cos_t2,    &
    !$omp         sin_t1, sin_t2, dt1, dt2,  &
    !$omp         abs_dt1, abs_dt2, cos_dt1, cos_dt2,    &
    !$omp         sin_dt1, sin_dt2, pwm_score_A, pwm_score_C,    &
    !$omp         pwm_score_G, pwm_score_T, ene_coef_r0,                  &
    !$omp         grad_coef_t1, ene_coef_t1, grad_coef_t2, ene_coef_t2,   &
    !$omp         etmp0, e_tmp_pwmcosns, viri,        &
    !$omp         var_tmp, grad_r0, grad_t1, grad_t2,     &
    !$omp         f_r0_coef_tmp, f_r0_tmp, f_t1_coef_tmp, f_t1_tmp,       &
    !$omp         f_t2_coef_tmp, f_t2_tmp,       &
    !$omp         f_CA_0, f_CA_N, f_CA_C,  f_DP_0, f_DS_0, start_i) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, enefunc%num_pwmcosns_domain, nthread

      ix = id_g2l(pwmcosns_id(i)) ! local index of i_CA_0 in atdyn
      rtmp(1:3) = coord(ix,1:3)
 
      i_CA_C_atom = id_g2l(pwmcosns_id_C(i))
      i_CA_N_atom = id_g2l(pwmcosns_id_N(i))
      vcn(1:3)    = coord(i_CA_N_atom,1:3) - coord(i_CA_C_atom,1:3) 
      rcn_sqr     = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
      rcn         = sqrt(rcn_sqr)
      rcn_inv     = 1.0_wp / rcn
      ecn(1:3)    = vcn(1:3) * rcn_inv

      num_nb15 = num_nb15_nobc(i)
      epwmcosns_temp = 0.0_wp
      force_CA_0(1:3) = 0.0_wp
      force_CA_N(1:3) = 0.0_wp
      force_CA_C(1:3) = 0.0_wp

      do icount = 1, pwmcosns_count(i)

        r00         = pwmcosns_r0    (icount, i)
        t10         = pwmcosns_theta1(icount, i)
        t20         = pwmcosns_theta2(icount, i)
        etmp0       = pwmcosns_ene(icount, i)

        cutoff      = r00 + 5.0_wp
        cutoff_sqr  = cutoff * cutoff

        do k = 1, num_nb15

          j  = nb15_calc_list(k,i) ! local index of i_DP_0 in atdyn

          vbc(1:3) = coord(ix,1:3) - coord(j,1:3) 
          rbc_sqr = vbc(1)*vbc(1) + vbc(2)*vbc(2) + vbc(3)*vbc(3)

          if (rbc_sqr >= cutoff_sqr) cycle

          f_CA_0(1:3) = 0.0_wp
          f_CA_N(1:3) = 0.0_wp
          f_CA_C(1:3) = 0.0_wp
          f_DP_0(1:3) = 0.0_wp
          f_DS_0(1:3) = 0.0_wp

          e_tmp_pwmcosns = 0.0_wp

          ! ==================================
          ! Distance and angle calculations...
          ! ==================================
          i_DS_0_atom = id_g2l(id_l2g(j) + 1)

          ! --------------
          ! phos <== base
          ! --------------
          !
          vbs(1:3)    = coord(i_DS_0_atom,1:3) - coord(j,1:3)  
          rbs_sqr     = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
          rbs         = sqrt(rbs_sqr)
          rbs_inv     = 1.0_wp / rbs
          ebs(1:3)    = vbs(1:3) * rbs_inv

          ! -----------
          ! base ==> Ca
          ! -----------
          !
          rbc      = sqrt(rbc_sqr)
          rbc_inv  = 1.0_wp / rbc
          ebc(1:3) = vbc(1:3) * rbc_inv

          ! -----------------------------
          ! Angle t1: sugar -- phos -- Ca
          ! -----------------------------
          !
          cos_t1 = ebc(1) * ebs(1) + ebc(2) * ebs(2) + ebc(3) * ebs(3)
          if (cos_t1 >  1.0_wp) cos_t1 =  1.0_wp
          if (cos_t1 < -1.0_wp) cos_t1 = -1.0_wp
          sin_t1 = sqrt(1.0_wp - cos_t1 * cos_t1)
          if (sin_t1 < EPS) sin_t1 = EPS
          t1     = acos(cos_t1)

          ! ---------------------------------------
          ! Angle t2: base0 -- Ca <==> Ca_N -- Ca_C
          ! ---------------------------------------
          !
          cos_t2 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
          if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
          if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
          sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
          if (sin_t2 < EPS) sin_t2 = EPS
          t2     = acos(cos_t2)

          ! ---------------
          ! Delta angles...
          ! ---------------
          !
          dt1     = t1 - t10
          dt2     = t2 - t20
          abs_dt1 = abs(dt1)
          abs_dt2 = abs(dt2)

          ! ============================================================================
          ! Energy/Force Calculation: gamma * (pwm_score + pwm_shift) * f * g1 * g2 * g3
          ! ============================================================================
          ! -----------------
          ! Simple angle test
          ! -----------------
          !
          if (abs_dt1 >= phi2 .or. abs_dt2 >= phi2 ) cycle

          ! -------------------
          ! Distance modulating
          ! -------------------
          !
          call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, grad_r0)

          ! ----------------
          ! Angle modulating
          ! ----------------
          if ( abs_dt1 < phi ) then
            ene_coef_t1  = 1.0_wp
            grad_coef_t1 = 0.0_wp
            grad_t1(1:6) = 0.0_wp
          else
            sin_dt1      = sin(ktheta * dt1)
            ene_coef_t1  = sin_dt1 * sin_dt1
            call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
            grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
          end if

          if ( abs_dt2 < phi ) then
            ene_coef_t2  = 1.0_wp
            grad_coef_t2 = 0.0_wp
            grad_t2(1:6) = 0.0_wp
          else
            sin_dt2      = sin(ktheta * dt2)
            ene_coef_t2  = sin_dt2 * sin_dt2
            call calculate_angle(vbc, vcn, rbc, rcn, cos_t2, sin_t2, grad_t2)
            grad_coef_t2 = ktheta * sin(2.0_wp * ktheta * dt2)
          end if

          ! ======
          ! Energy
          ! ======
          !
          e_tmp_pwmcosns = etmp0*ene_coef_r0*ene_coef_t1*ene_coef_t2
          epwmcosns_temp = epwmcosns_temp + e_tmp_pwmcosns

          ! ======
          ! Forces
          ! ======
          !
          ! --
          ! r0
          ! --
          !
          f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2
          f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
          f_DP_0(1:3)   = f_DP_0(1:3) - f_r0_tmp(1:3)
          !
          ! --
          ! t1
          ! --
          !
          f_t1_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t2 * grad_coef_t1
          f_t1_tmp(1:6) = - f_t1_coef_tmp * grad_t1(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t1_tmp(1:3)
          f_DP_0(1:3)   = f_DP_0(1:3) - f_t1_tmp(1:3) - f_t1_tmp(4:6)
          f_DS_0(1:3)   = f_DS_0(1:3) + f_t1_tmp(4:6)
          !
          ! --
          ! t2
          ! --
          !
          f_t2_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t1 * grad_coef_t2
          f_t2_tmp(1:6) = - f_t2_coef_tmp * grad_t2(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t2_tmp(1:3)
          f_DP_0(1:3)   = f_DP_0(1:3) - f_t2_tmp(1:3)
          f_CA_C(1:3)   = f_CA_C(1:3) - f_t2_tmp(4:6)
          f_CA_N(1:3)   = f_CA_N(1:3) + f_t2_tmp(4:6)
          !
          ! ------------------------
          ! update all the forces...
          ! ------------------------
          !
          force_CA_0(1:3) = force_CA_0(1:3) + f_CA_0(1:3)
          force_CA_N(1:3) = force_CA_N(1:3) + f_CA_N(1:3)
          force_CA_C(1:3) = force_CA_C(1:3) + f_CA_C(1:3)
          force(j,1:3,id+1) = force(j,1:3,id+1) + f_DP_0(1:3)
          force(i_DS_0_atom,1:3,id+1) = force(i_DS_0_atom,1:3,id+1) &
                                        + f_DS_0(1:3)

          ! virial
          !
          viri(1:3) = f_CA_0(1:3) * coord(ix,1:3)          &
                    + f_CA_N(1:3) * coord(i_CA_N_atom,1:3) &
                    + f_CA_C(1:3) * coord(i_CA_C_atom,1:3) &
                    + f_DP_0(1:3) * coord(j,1:3)           &
                    + f_DS_0(1:3) * coord(i_DS_0_atom,1:3)
          virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
          virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
          virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

        end do

      end do

      force(ix,1:3,id+1) = force(ix,1:3,id+1) + force_CA_0(1:3)
      force(i_CA_C_atom,1:3,id+1) = force(i_CA_C_atom,1:3,id+1) &
                                  + force_CA_C(1:3)
      force(i_CA_N_atom,1:3,id+1) = force(i_CA_N_atom,1:3,id+1) &
                                  + force_CA_N(1:3)
      epwmcosns(id+1) = epwmcosns(id+1) + epwmcosns_temp

    end do

    !$omp end parallel

    call timer(TimerCGPwmcosns, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_PWMcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_nonlocal_dihedral
  !> @brief        Calculates dihedral angle.
  !! @authors      CK; CT
  !! @param[in]    aindex  : atom indices
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    s_vec   : translational vector
  !! @param[inout] cos_dih : cosine of dihedral angles
  !! @param[inout] sin_dih : sine of dihedral angles
  !! @param[inout] grad    : gradient of dihedral angles
  !
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !! @note         Copied from src/lib/dihedral_libs.fpp
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_nonlocal_dihedral(dij, djk, dlk, cos_dih, sin_dih, grad)

    ! formal arguments
    real(wp),  intent(in)    :: dij(3), djk(3), dlk(3)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(1:9)

    ! local variable
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl

    aijk(1) = dij(2)*djk(3) - dij(3)*djk(2)
    aijk(2) = dij(3)*djk(1) - dij(1)*djk(3)
    aijk(3) = dij(1)*djk(2) - dij(2)*djk(1)

    ajkl(1) = dlk(2)*djk(3) - dlk(3)*djk(2)
    ajkl(2) = dlk(3)*djk(1) - dlk(1)*djk(3)
    ajkl(3) = dlk(1)*djk(2) - dlk(2)*djk(1)

    raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
    rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)

    inv_raijk2  = 1.0_wp / raijk2
    inv_rajkl2  = 1.0_wp / rajkl2

    inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

    cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl

    rjk     = sqrt(djk(1)*djk(1) + djk(2)*djk(2) + djk(3)*djk(3))
    inv_rjk = 1.0_wp/rjk

    tmp(1)  = aijk(1)*dlk(1) + aijk(2)*dlk(2) + aijk(3)*dlk(3)
    sin_dih = tmp(1) * rjk * inv_raijkl

    dotpro_ijk = dij(1)*djk(1) + dij(2)*djk(2) + dij(3)*djk(3)
    dotpro_jkl = djk(1)*dlk(1) + djk(2)*dlk(2) + djk(3)*dlk(3)

    tmp(1) = rjk*inv_raijk2
    tmp(2) = rjk*inv_rajkl2

    tmp(3) =  dotpro_ijk*inv_raijk2*inv_rjk
    tmp(4) =  dotpro_jkl*inv_rajkl2*inv_rjk

    ! This part is different from the one in /lib/dihedral_libs
    !
    grad(1:3) = -tmp(1)*aijk(1:3)
    grad(4:6) =  tmp(3)*aijk(1:3) - tmp(4)*ajkl(1:3)
    grad(7:9) =                   + tmp(2)*ajkl(1:3)

    return

  end subroutine calculate_nonlocal_dihedral

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_repulsive
  !> @brief        Calculates repulsive (potential) gradients.
  !! @authors      CT
  !! @param[in]    vec     : vector connecting two particles
  !! @param[in]    dist    : distance
  !! @param[in]    sigma   : distance in reference struct
  !! @param[in]    alpha   : alpha in U_morse
  !! @param[in]    epsilon : energy coefficient
  !! @param[inout] ene     : repulsive energy
  !! @param[inout] grad    : repulsive gradient
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_repulsive(vec, dist, sigma, alpha, epsilon, ene, grad)

    ! formal arguments
    !
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: dist
    real(wp),                intent(in)    :: sigma
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: epsilon
    real(wp),                intent(inout) :: ene
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    !
    real(wp)                 :: expon, exptmp, etmp
    real(wp)                 :: f_coef

    if ( dist < sigma ) then
      expon  = - alpha * (dist - sigma)
      exptmp = exp(expon)
      etmp   = 1.0_wp - exptmp
      ene    = epsilon * etmp * etmp
      f_coef = -2.0_wp * alpha * epsilon * etmp * exptmp / dist
      grad(1:3) = f_coef * vec(1:3)
    else
      ene    = 0.0_wp
      grad(1:3)   = 0.0_wp
    end if

  end subroutine calculate_repulsive

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_attractive
  !> @brief        Calculates 3SPN.2C Morse (potential) gradients.
  !! @authors      CT
  !! @param[in]    vec     : vector connecting two particles
  !! @param[in]    dist    : distance
  !! @param[in]    sigma   : distance in reference struct
  !! @param[in]    alpha   : alpha in U_morse
  !! @param[in]    epsilon : energy coefficient
  !! @param[inout] ene     : attractive energy
  !! @param[inout] grad    : attractive gradient
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_attractive(vec, dist, sigma, alpha, epsilon, ene, grad)

    ! formal arguments
    !
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: dist
    real(wp),                intent(in)    :: sigma
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: epsilon
    real(wp),                intent(inout) :: ene
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    !
    real(wp)                 :: expon, exptmp, etmp
    real(wp)                 :: f_coef

    if ( dist > sigma ) then
      expon  = - alpha * (dist - sigma)
      exptmp = exp(expon)
      etmp   = 1.0_wp - exptmp
      ene    = epsilon * etmp * etmp - epsilon
      f_coef = -2.0_wp * alpha * epsilon * etmp * exptmp / dist
      grad(1:3) = f_coef * vec(1:3)
    else
      ene    = - epsilon
      grad(1:3)   = 0.0_wp
    end if

  end subroutine calculate_attractive

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_angle
  !> @brief        Calculates angle gradients.
  !! @authors      CT
  !! @param[in]    dr1     : vector connecting particle 1 and 2
  !! @param[in]    dr2     : vector connecting particle 2 and 3
  !! @param[in]    r1      : distance between particle 1 and 2
  !! @param[in]    r2      : distance between particle 2 and 3
  !! @param[in]    cos_ang : cosine of angle
  !! @param[in]    sin_ang : sine of angle
  !! @param[inout] grad    : gradient of angle
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_angle(dr1, dr2, r1, r2, cos_ang, sin_ang, grad)

    ! formal arguments
    !
    real(wp),  intent(in)    :: dr1(:)
    real(wp),  intent(in)    :: dr2(:)
    real(wp),  intent(in)    :: r1
    real(wp),  intent(in)    :: r2
    real(wp),  intent(in)    :: cos_ang
    real(wp),  intent(in)    :: sin_ang
    real(wp),  intent(inout) :: grad(1:6)

    ! local variables
    !
    real(wp) :: tmp1
    real(wp) :: tmp2
    real(wp) :: tmp3

    tmp1 = 1.0_wp  / ( sin_ang * r1 * r2 )
    tmp2 = cos_ang / ( sin_ang * r1 * r1 )
    tmp3 = cos_ang / ( sin_ang * r2 * r2 )

    grad(1:3) = - tmp1 * dr2(1:3) + tmp2 * dr1(1:3)
    grad(4:6) = - tmp1 * dr1(1:3) + tmp3 * dr2(1:3)

    return

  end subroutine calculate_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_gaussian
  !> @brief        Calculates Gaussian (potential) gradients.
  !! @authors      CT
  !! @param[in]    vec     : vector connecting two particles
  !! @param[in]    r : r
  !! @param[in]    r0 : r0
  !! @param[in]    sigma_sqr   : distance in reference struct
  !! @param[inout] enecoef : gaussian energy coefficient
  !! @param[inout] grad    : gaussian gradient
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_gaussian(vec, r, r0, sigma_sqr, enecoef, grad)

    ! formal arguments
    !
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: r
    real(wp),                intent(in)    :: r0
    real(wp),                intent(in)    :: sigma_sqr
    real(wp),                intent(inout) :: enecoef
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    !
    real(wp)                 :: delta_r
    real(wp)                 :: expon
    real(wp)                 :: f_coef

    delta_r = r - r0
    expon  = - 0.5_wp * delta_r * delta_r / sigma_sqr
    enecoef   = exp(expon)
    f_coef = - delta_r * enecoef / sigma_sqr / r
    grad(1:3) = f_coef * vec(1:3)

  end subroutine calculate_gaussian

end module cg_energy_nonlocal_mod

