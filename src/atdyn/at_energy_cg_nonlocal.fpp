!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_cg_nonlocal_mod
!> @brief   calculate base-base interaction energy
!! @authors Cheng Tan (CT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_cg_nonlocal_mod

  use at_pairlist_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_DNA_base_stacking
  public  :: compute_energy_DNA_base_stacking_pbc
  public  :: compute_energy_DNA_base_pairing
  public  :: compute_energy_DNA_base_pairing_pbc
  public  :: compute_energy_DNA_exv
  public  :: compute_energy_DNA_exv_pbc
  public  :: compute_energy_CG_ele
  public  :: compute_energy_CG_ele_pbc
  public  :: compute_energy_CG_pwmcos
  public  :: compute_energy_CG_pwmcos_pbc
  public  :: compute_energy_CG_pwmcosns
  public  :: compute_energy_CG_pwmcosns_pbc
  public  :: compute_energy_CG_IDR_HPS
  public  :: compute_energy_CG_IDR_HPS_pbc
  public  :: compute_energy_CG_IDR_KH
  public  :: compute_energy_CG_IDR_KH_pbc
  public  :: compute_energy_CG_KH
  public  :: compute_energy_CG_KH_pbc

  private :: calculate_repulsive
  private :: calculate_attractive
  private :: calculate_nonlocal_dihedral
  private :: calculate_nonlocal_dihedral_pbc
  private :: calculate_angle
  private :: calculate_gaussian

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_stacking
  !> @brief        calculate base stacking energy
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ebase   : base energy of target systems
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_stacking(enefunc, coord, force, virial, &
                                              ebase)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ebase

    ! local variables
    integer                  :: i, j, k
    integer                  :: istart, iend
    integer                  :: i_sugar, i_base5, i_base3
    real(wp)                 :: v21(1:3), r21_sqr, r21
    real(wp)                 :: v23(1:3), r23_sqr, r23
    real(wp)                 :: dot_v21_v23, cos_t_bs, sin_t_bs
    real(wp)                 :: delta_t_bs, abs_delta_t_bs, delta_r23
    real(wp)                 :: expnt, exp_expnt, m_exp_expnt
    real(wp)                 :: K_theta, alpha
    real(wp)                 :: ene_morse, grad_morse
    real(wp)                 :: ene_repl, grad_repl, ene_attr, grad_attr
    real(wp)                 :: cos_dt, coef_bell
    real(wp)                 :: grad_bell_coef_0, grad_bell_coef_1, grad_bell_coef_3
    real(wp)                 :: vtmp
    real(wp)                 :: g_attr_1(1:3), g_attr_3(1:3)
    real(wp)                 :: theta_bs_threshold_1, theta_bs_threshold_2

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: ene_coef(:), sigma(:), theta_bs(:)
    real(wp),        pointer :: work(:,:)


    !==========================================================================
    !                                   func
    ! U_bs = U_repl + U_attr
    !
    ! U_repl = U_morse                          ( r < r0 )
    !        = 0                                ( r > r0 )
    !
    ! U_attr = f_bell * ( - epsilon )           ( r < r0 )
    !        = f_bell * ( U_morse - epsilon)    ( r > r0 )
    !--------------------------------------------------------------------------
    !                                 topology
    !    r12
    ! S-----B5
    ! |        r23       ;     theta_bs: S-B5-B3
    ! o-----B3
    !
    !==========================================================================

    call timer(TimerBaseStack, TimerOn)

    ! use pointers
    !
    istart   = enefunc%istart_base_stack
    iend     = enefunc%iend_base_stack
    list     => enefunc%base_stack_list
    ene_coef => enefunc%base_stack_epsilon
    sigma    => enefunc%base_stack_sigma
    theta_bs => enefunc%base_stack_theta_bs
    work     => enefunc%work

    K_theta = enefunc%base_stack_K
    alpha = enefunc%base_stack_alpha

    theta_bs_threshold_2 = PI / K_theta
    theta_bs_threshold_1 = theta_bs_threshold_2 / 2.0_wp

    ! calculation of base energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, i_sugar, i_base5, i_base3,        &
    !$omp         v21, r21_sqr, r21, v23, r23_sqr, r23,      &
    !$omp         dot_v21_v23, cos_t_bs, sin_t_bs,           &
    !$omp         delta_t_bs, abs_delta_t_bs, delta_r23,     &
    !$omp         expnt, exp_expnt, m_exp_expnt,             &
    !$omp         ene_morse, grad_morse, ene_repl,           &
    !$omp         grad_repl, ene_attr, grad_attr,            &
    !$omp         cos_dt, coef_bell, grad_bell_coef_0,       &
    !$omp         grad_bell_coef_1, grad_bell_coef_3,        &
    !$omp         g_attr_1, g_attr_3, vtmp)                  &
    !$omp shared(istart, iend, ene_coef, sigma, theta_bs,    &
    !$omp        theta_bs_threshold_2, theta_bs_threshold_1, &
    !$omp        K_theta, alpha,                             &
    !$omp        work, coord, list)                          &
    !$omp reduction(+:ebase) reduction(+:virial)
    !
    do i = istart, iend

      ! base stack energy
      !
      i_sugar = list(1, i)
      i_base5 = list(2, i)
      i_base3 = list(3, i)

      ! set all temporary "work" elements to be zero!
      !
      work(1:9, i) = 0.0_wp

      ! r21: Distance between sugar and base5
      !
      v21(1:3) = coord(1:3, i_sugar) - coord(1:3, i_base5)
      r21_sqr  = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
      r21      = sqrt(r21_sqr)

      ! r23: Distance between base5 and base3
      !
      v23(1:3) = coord(1:3, i_base3) - coord(1:3, i_base5)
      r23_sqr  = v23(1)*v23(1) + v23(2)*v23(2) + v23(3)*v23(3)
      r23      = sqrt(r23_sqr)

      ! cos_t_bs: cos(theta_bs)
      ! sin_t_bs: sin(theta_bs)
      !
      dot_v21_v23 = v21(1) * v23(1) + v21(2) * v23(2) + v21(3) * v23(3)
      cos_t_bs    = dot_v21_v23 / (r21 * r23)
      if (cos_t_bs > 1.0_wp)  cos_t_bs =  1.0_wp
      if (cos_t_bs < -1.0_wp) cos_t_bs = -1.0_wp
      sin_t_bs       = sqrt(1.0_wp - cos_t_bs * cos_t_bs)

      ! Delta theta_bs
      !
      delta_t_bs     = acos(cos_t_bs) - theta_bs(i)
      abs_delta_t_bs = abs(delta_t_bs)

      ! Delta distance: r23 - sigma
      !
      delta_r23      = r23 - sigma(i)

      ! ------------
      ! common terms
      ! ------------
      expnt       = - alpha * delta_r23 ! - alpha (r - r0)
      exp_expnt   = exp(expnt)          ! exp(-alpha (r - r0))
      m_exp_expnt = 1.0_wp - exp_expnt  ! 1 - exp(-alpha (r - r0))
      ene_morse   = ene_coef(i) * m_exp_expnt * m_exp_expnt
      grad_morse  = 2.0_wp * alpha * &
          ene_coef(i) * exp_expnt * m_exp_expnt / r23

      ! ---------
      ! Repulsion
      ! ---------
      if (delta_r23 < 0.0_wp) then
        ene_repl  = ene_morse
        grad_repl = grad_morse
      else
        ene_repl  = 0.0_wp
        grad_repl = 0.0_wp
      end if
      work(4:6, i) =   grad_repl * v23(1:3)
      work(7:9, i) = - grad_repl * v23(1:3)

      ! ----------
      ! Attraction
      ! ----------
      if (delta_r23 < 0.0_wp) then
        ene_attr  = - ene_coef(i)
        grad_attr = 0.0_wp
      else
        ene_attr  = ene_morse - ene_coef(i)
        grad_attr = grad_morse
      end if

      if (abs_delta_t_bs < theta_bs_threshold_1) then
        coef_bell = 1.0_wp
        work(4:6, i) = work(4:6, i) + grad_attr * v23(1:3)
        work(7:9, i) = work(7:9, i) - grad_attr * v23(1:3)
      else if (abs_delta_t_bs < theta_bs_threshold_2) then
        cos_dt    = cos(K_theta*delta_t_bs)
        coef_bell = 1.0_wp - cos_dt * cos_dt
        grad_bell_coef_0 = - ene_attr * K_theta    &
            * sin(2.0_wp * K_theta * delta_t_bs) &
            / (sin_t_bs * r21 * r23)
        grad_bell_coef_1 = - grad_bell_coef_0 * dot_v21_v23 / r21_sqr
        grad_bell_coef_3 = - grad_bell_coef_0 * dot_v21_v23 / r23_sqr
        g_attr_1(1:3) = grad_bell_coef_0 * v23(1:3) + grad_bell_coef_1 * v21(1:3)
        g_attr_3(1:3) = grad_bell_coef_0 * v21(1:3) + grad_bell_coef_3 * v23(1:3)
        work(1:3, i) = work(1:3, i) - g_attr_1(1:3)
        work(4:6, i) = work(4:6, i) + g_attr_1(1:3) + g_attr_3(1:3) + coef_bell * grad_attr * v23(1:3)
        work(7:9, i) = work(7:9, i) - g_attr_3(1:3) - coef_bell * grad_attr * v23(1:3)
      else
        coef_bell = 0.0_wp
      end if

      ebase = ebase + ene_repl + coef_bell * ene_attr

      do j = 1, 3
        do k = j+1, 3
          vtmp = work(j, i) * coord(k, i_sugar)    &
              + work(j + 3, i) * coord(k, i_base5) &
              + work(j + 6, i) * coord(k, i_base3)
          virial(k, j) = virial(k, j) - vtmp
          virial(j, k) = virial(j, k) - vtmp
        end do
        vtmp =    work(j, i) * coord(j, i_sugar) &
            + work(j + 3, i) * coord(j, i_base5) &
            + work(j + 6, i) * coord(j, i_base3)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    do i = istart, iend
      force(1:3, list(1, i),1) = force(1:3, list(1, i),1) + work(1:3, i)
      force(1:3, list(2, i),1) = force(1:3, list(2, i),1) + work(4:6, i)
      force(1:3, list(3, i),1) = force(1:3, list(3, i),1) + work(7:9, i)
    end do

    call timer(TimerBaseStack, TimerOff)

    return

  end subroutine compute_energy_DNA_base_stacking

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_stacking_pbc
  !> @brief        calculate base stacking energy
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ebase    : base energy of target systems
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_stacking_pbc(enefunc, boundary, &
                                                  coord, force, virial, ebase)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ebase

    ! local variables
    integer           :: i, j, k
    integer           :: istart, iend
    integer           :: i_sugar, i_base5, i_base3
    real(wp)          :: v21(1:3), r21_sqr, r21
    real(wp)          :: v23(1:3), r23_sqr, r23
    real(wp)          :: dot_v21_v23, cos_t_bs, sin_t_bs
    real(wp)          :: delta_t_bs, abs_delta_t_bs, delta_r23
    real(wp)          :: expnt, exp_expnt, m_exp_expnt
    real(wp)          :: K_theta, alpha
    real(wp)          :: ene_morse, grad_morse
    real(wp)          :: ene_repl, grad_repl, ene_attr, grad_attr
    real(wp)          :: cos_dt, coef_bell
    real(wp)          :: grad_bell_coef_0, grad_bell_coef_1, grad_bell_coef_3
    real(wp)          :: vtmp
    real(wp)          :: g_attr_1(1:3), g_attr_3(1:3)
    real(wp)          :: theta_bs_threshold_1, theta_bs_threshold_2
    real(wp)          :: bsize(3), half_bsize(3)

    integer,  pointer :: list(:,:)
    real(wp), pointer :: ene_coef(:), sigma(:), theta_bs(:)
    real(wp), pointer :: work(:,:)


    !==========================================================================
    !                                   func
    ! U_bs = U_repl + U_attr
    !
    ! U_repl = U_morse                          ( r < r0 )
    !        = 0                                ( r > r0 )
    !
    ! U_attr = f_bell * ( - epsilon )           ( r < r0 )
    !        = f_bell * ( U_morse - epsilon)    ( r > r0 )
    !--------------------------------------------------------------------------
    !                                 topology
    !    r12
    ! S-----B5
    ! |        r23       ;     theta_bs: S-B5-B3
    ! o-----B3
    !
    !==========================================================================

    call timer(TimerBaseStack, TimerOn)

    ! use pointers
    !
    istart   = enefunc%istart_base_stack
    iend     = enefunc%iend_base_stack
    list     => enefunc%base_stack_list
    ene_coef => enefunc%base_stack_epsilon
    sigma    => enefunc%base_stack_sigma
    theta_bs => enefunc%base_stack_theta_bs
    work     => enefunc%work

    K_theta = enefunc%base_stack_K
    alpha = enefunc%base_stack_alpha

    theta_bs_threshold_2 = PI / K_theta
    theta_bs_threshold_1 = theta_bs_threshold_2 / 2.0_wp

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of base energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, i_sugar, i_base5, i_base3,        &
    !$omp         v21, r21_sqr, r21, v23, r23_sqr, r23,      &
    !$omp         dot_v21_v23, cos_t_bs, sin_t_bs,           &
    !$omp         delta_t_bs, abs_delta_t_bs, delta_r23,     &
    !$omp         expnt, exp_expnt, m_exp_expnt,             &
    !$omp         ene_morse, grad_morse, ene_repl,           &
    !$omp         grad_repl, ene_attr, grad_attr,            &
    !$omp         cos_dt, coef_bell, grad_bell_coef_0,       &
    !$omp         grad_bell_coef_1, grad_bell_coef_3,        &
    !$omp         g_attr_1, g_attr_3, vtmp)                  &
    !$omp shared(istart, iend, ene_coef, sigma, theta_bs,    &
    !$omp        theta_bs_threshold_2, theta_bs_threshold_1, &
    !$omp        K_theta, alpha,                             &
    !$omp        bsize, half_bsize,                          &
    !$omp        work, coord, list)                          &
    !$omp reduction(+:ebase) reduction(+:virial)
    !
    do i = istart, iend

      ! base stack energy
      !
      i_sugar = list(1, i)
      i_base5 = list(2, i)
      i_base3 = list(3, i)

      ! set all temporary "work" elements to be zero!
      !
      work(1:9, i) = 0.0_wp

      ! r21: Distance between sugar and base5
      !
      v21(1:3) = coord(1:3, i_sugar) - coord(1:3, i_base5)
      if (v21(1) > half_bsize(1)) then
        v21(1) = v21(1) - bsize(1)
      else if (v21(1) < -half_bsize(1)) then
        v21(1) = v21(1) + bsize(1)
      end if
      if (v21(2) > half_bsize(2)) then
        v21(2) = v21(2) - bsize(2)
      else if (v21(2) < -half_bsize(2)) then
        v21(2) = v21(2) + bsize(2)
      end if
      if (v21(3) > half_bsize(3)) then
        v21(3) = v21(3) - bsize(3)
      else if (v21(3) < -half_bsize(3)) then
        v21(3) = v21(3) + bsize(3)
      end if

      r21_sqr  = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
      r21      = sqrt(r21_sqr)

      ! r23: Distance between base5 and base3
      !
      v23(1:3) = coord(1:3, i_base3) - coord(1:3, i_base5)
      if (v23(1) > half_bsize(1)) then
        v23(1) = v23(1) - bsize(1)
      else if (v23(1) < -half_bsize(1)) then
        v23(1) = v23(1) + bsize(1)
      end if
      if (v23(2) > half_bsize(2)) then
        v23(2) = v23(2) - bsize(2)
      else if (v23(2) < -half_bsize(2)) then
        v23(2) = v23(2) + bsize(2)
      end if
      if (v23(3) > half_bsize(3)) then
        v23(3) = v23(3) - bsize(3)
      else if (v23(3) < -half_bsize(3)) then
        v23(3) = v23(3) + bsize(3)
      end if

      r23_sqr  = v23(1)*v23(1) + v23(2)*v23(2) + v23(3)*v23(3)
      r23      = sqrt(r23_sqr)

      ! cos_t_bs: cos(theta_bs)
      ! sin_t_bs: sin(theta_bs)
      !
      dot_v21_v23 = v21(1) * v23(1) + v21(2) * v23(2) + v21(3) * v23(3)
      cos_t_bs    = dot_v21_v23 / (r21 * r23)
      if (cos_t_bs > 1.0_wp)  cos_t_bs =  1.0_wp
      if (cos_t_bs < -1.0_wp) cos_t_bs = -1.0_wp
      sin_t_bs       = sqrt(1.0_wp - cos_t_bs * cos_t_bs)

      ! Delta theta_bs
      !
      delta_t_bs     = acos(cos_t_bs) - theta_bs(i)
      abs_delta_t_bs = abs(delta_t_bs)

      ! Delta distance: r23 - sigma
      !
      delta_r23      = r23 - sigma(i)

      ! ------------
      ! common terms
      ! ------------
      expnt       = - alpha * delta_r23 ! - alpha (r - r0)
      exp_expnt   = exp(expnt)          ! exp(-alpha (r - r0))
      m_exp_expnt = 1.0_wp - exp_expnt  ! 1 - exp(-alpha (r - r0))
      ene_morse   = ene_coef(i) * m_exp_expnt * m_exp_expnt
      grad_morse  = 2.0_wp * alpha * &
          ene_coef(i) * exp_expnt * m_exp_expnt / r23

      ! ---------
      ! Repulsion
      ! ---------
      if (delta_r23 < 0.0_wp) then
        ene_repl  = ene_morse
        grad_repl = grad_morse
      else
        ene_repl  = 0.0_wp
        grad_repl = 0.0_wp
      end if
      work(4:6, i) =   grad_repl * v23(1:3)
      work(7:9, i) = - grad_repl * v23(1:3)

      ! ----------
      ! Attraction
      ! ----------
      if (delta_r23 < 0.0_wp) then
        ene_attr  = - ene_coef(i)
        grad_attr = 0.0_wp
      else
        ene_attr  = ene_morse - ene_coef(i)
        grad_attr = grad_morse
      end if

      if (abs_delta_t_bs < theta_bs_threshold_1) then
        coef_bell = 1.0_wp
        work(4:6, i) = work(4:6, i) + grad_attr * v23(1:3)
        work(7:9, i) = work(7:9, i) - grad_attr * v23(1:3)
      else if (abs_delta_t_bs < theta_bs_threshold_2) then
        cos_dt    = cos(K_theta*delta_t_bs)
        coef_bell = 1.0_wp - cos_dt * cos_dt
        grad_bell_coef_0 = - ene_attr * K_theta    &
            * sin(2.0_wp * K_theta * delta_t_bs) &
            / (sin_t_bs * r21 * r23)
        grad_bell_coef_1 = - grad_bell_coef_0 * dot_v21_v23 / r21_sqr
        grad_bell_coef_3 = - grad_bell_coef_0 * dot_v21_v23 / r23_sqr
        g_attr_1(1:3) = grad_bell_coef_0 * v23(1:3) + grad_bell_coef_1 * v21(1:3)
        g_attr_3(1:3) = grad_bell_coef_0 * v21(1:3) + grad_bell_coef_3 * v23(1:3)
        work(1:3, i) = work(1:3, i) - g_attr_1(1:3)
        work(4:6, i) = work(4:6, i) + g_attr_1(1:3) + g_attr_3(1:3) + coef_bell * grad_attr * v23(1:3)
        work(7:9, i) = work(7:9, i) - g_attr_3(1:3) - coef_bell * grad_attr * v23(1:3)
      else
        coef_bell = 0.0_wp
      end if

      ebase = ebase + ene_repl + coef_bell * ene_attr

      do j = 1, 3
        do k = j+1, 3
          vtmp = work(j, i) * coord(k, i_sugar)    &
              + work(j + 3, i) * coord(k, i_base5) &
              + work(j + 6, i) * coord(k, i_base3)
          virial(k, j) = virial(k, j) - vtmp
          virial(j, k) = virial(j, k) - vtmp
        end do
        vtmp =    work(j, i) * coord(j, i_sugar) &
            + work(j + 3, i) * coord(j, i_base5) &
            + work(j + 6, i) * coord(j, i_base3)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    do i = istart, iend
      force(1:3, list(1, i),1) = force(1:3, list(1, i),1) + work(1:3, i)
      force(1:3, list(2, i),1) = force(1:3, list(2, i),1) + work(4:6, i)
      force(1:3, list(3, i),1) = force(1:3, list(3, i),1) + work(7:9, i)
    end do

    call timer(TimerBaseStack, TimerOff)

    return

  end subroutine compute_energy_DNA_base_stacking_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_pairing
  !> @brief        calculate base-pairing energy with pairlist (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epair    : base-pairing energy
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_pairing(enefunc, pairlist, coord, force, &
                                             virial, epair)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epair

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_bp, ini_bp, fin_bp
    integer                   :: i, j, k, l, m
    real(wp)                  :: cutoff, cutoff_sqr

    integer                   :: i_S1, i_S3
    integer                   :: i_B2, i_B4
    integer                   :: i_B5, i_B6
    integer                   :: type_B2, type_B4
    integer                   :: type_B5, type_B6
    integer                   :: aindex(1:4)
    real(wp)                  :: d21(3), r21_sqr, r21, r21_inv, e21(3)
    real(wp)                  :: d42(3), r42_sqr, r42, r42_inv, e42(3)
    real(wp)                  :: d24(3), e24(3)
    real(wp)                  :: d43(3), r43_sqr, r43, r43_inv, e43(3)
    real(wp)                  :: d25(3), r25_sqr, r25, r25_inv, e25(3)
    real(wp)                  :: d46(3), r46_sqr, r46, r46_inv, e46(3)
    real(wp)                  :: cos_124, cos_342, sin_124, sin_342
    real(wp)                  :: cos_125, cos_346, sin_125, sin_346
    real(wp)                  :: cos_1243, sin_1243
    real(wp)                  :: angle_124, angle_342
    real(wp)                  :: angle_125, angle_346
    real(wp)                  :: angle_1243
    real(wp)                  :: delta_angle_124, delta_angle_342
    real(wp)                  :: delta_angle_125, delta_angle_346
    real(wp)                  :: delta_angle_1243
    real(wp)                  :: abs_delta_angle_124, abs_delta_angle_342
    real(wp)                  :: abs_delta_angle_125, abs_delta_angle_346
    real(wp)                  :: abs_delta_angle_1243
    real(wp)                  :: grad_angle_124(1:6)
    real(wp)                  :: grad_angle_342(1:6)
    real(wp)                  :: grad_angle_125(1:6)
    real(wp)                  :: grad_angle_346(1:6)
    real(wp)                  :: grad_angle_1243(1:6)
    real(wp)                  :: force_angle_124_tmp(1:6)
    real(wp)                  :: force_angle_342_tmp(1:6)
    real(wp)                  :: force_angle_125_tmp(1:6)
    real(wp)                  :: force_angle_346_tmp(1:6)
    real(wp)                  :: force_angle_1243_tmp(1:6)
    real(wp)                  :: cos_dih_24, sin_dih_24
    real(wp)                  :: grad_dih_24(1:9)
    real(wp)                  :: force_dih_24_tmp(1:9)
    real(wp)                  :: grad_repl(1:3), ene_repl
    real(wp)                  :: grad_attr(1:3), ene_attr
    real(wp)                  :: force_repl_tmp(1:3)
    real(wp)                  :: force_attr_tmp(1:3)
    real(wp)                  :: grad_coef_phi, ene_coef_phi
    real(wp)                  :: grad_coef_theta_1, ene_coef_theta_1
    real(wp)                  :: grad_coef_theta_2, ene_coef_theta_2
    real(wp)                  :: grad_coef_theta_3, ene_coef_theta_3
    real(wp)                  :: grad_coef_theta_cs, ene_coef_theta_cs
    real(wp)                  :: force_coef_tmp
    real(wp)                  :: var_tmp
    real(wp)                  :: force_S1(3), force_S3(3)
    real(wp)                  :: force_B2(3), force_B4(3)
    real(wp)                  :: force_B5(3), force_B6(3)

    real(wp)                  :: bp_theta_1_0
    real(wp)                  :: bp_theta_2_0
    real(wp)                  :: bp_theta_3_0
    real(wp)                  :: bp_phi_1_0
    real(wp)                  :: bp_sigma
    real(wp)                  :: bp_epsilon
    real(wp)                  :: cs_1_epsilon
    real(wp)                  :: cs_1_sigma
    real(wp)                  :: cs_1_theta_cs_0
    real(wp)                  :: cs_2_epsilon
    real(wp)                  :: cs_2_sigma
    real(wp)                  :: cs_2_theta_cs_0

    integer,  pointer         :: bp_list(:,:)
    integer,  pointer         :: num_bp_calc(:,:)
    integer,  pointer         :: base_type(:)
    integer,  pointer         :: chain_id(:)
    real(wp), pointer         :: param_bp_theta_1(:)
    real(wp), pointer         :: param_bp_theta_2(:)
    real(wp), pointer         :: param_bp_theta_3(:)
    real(wp), pointer         :: param_bp_phi_1(:)
    real(wp), pointer         :: param_bp_sigma(:)
    real(wp), pointer         :: param_bp_epsilon(:)
    real(wp), pointer         :: param_cs_1_epsilon(:,:)
    real(wp), pointer         :: param_cs_1_sigma(:,:)
    real(wp), pointer         :: param_cs_1_theta_cs(:,:)
    real(wp), pointer         :: param_cs_2_epsilon(:,:)
    real(wp), pointer         :: param_cs_2_sigma(:,:)
    real(wp), pointer         :: param_cs_2_theta_cs(:,:)

    real(wp)                  :: param_bp_alpha
    real(wp)                  :: param_bp_K
    real(wp)                  :: param_cs_alpha
    real(wp)                  :: param_cs_K
    real(wp)                  :: bp_theta_threshold_2
    real(wp)                  :: bp_theta_threshold_1
    real(wp)                  :: cs_theta_threshold_2
    real(wp)                  :: cs_theta_threshold_1

    real(wp)                  :: e_tmp_bp
    real(wp)                  :: e_tmp_cs


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

    natom                = size(coord(1,:))
    cutoff               = enefunc%cg_cutoffdist_DNAbp
    cutoff_sqr           = cutoff * cutoff

    bp_list              => pairlist%cg_DNA_basepair_list
    num_bp_calc          => pairlist%num_cg_DNA_basepair_calc
    base_type            => enefunc%NA_base_type
    chain_id             => enefunc%mol_chain_id
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

    bp_theta_threshold_2 = PI / param_bp_K
    bp_theta_threshold_1 = bp_theta_threshold_2 / 2.0_wp
    cs_theta_threshold_2 = PI / param_cs_K
    cs_theta_threshold_1 = cs_theta_threshold_2 / 2.0_wp

    num_bp               = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                             &
    !$omp firstprivate(num_bp)                               &
    !$omp private(my_id, id,                                 &
    !$omp         ini_bp, fin_bp,                            &
    !$omp         i, k, j, l, m,                             &
    !$omp         i_S1, i_S3, i_B2, i_B4, i_B5, i_B6,        &
    !$omp         type_B2, type_B4, type_B5, type_B6,        &
    !$omp         aindex,                                    &
    !$omp         d21, r21_sqr, r21, r21_inv, e21,           &
    !$omp         d42, r42_sqr, r42, r42_inv, e42, d24, e24, &
    !$omp         d43, r43_sqr, r43, r43_inv, e43,           &
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
    !$omp         e_tmp_bp, e_tmp_cs                         &
    !$omp         )                                          &
    !$omp shared(coord, force, my_city_rank, nproc_city,     &
    !$omp        nthread, natom,                             &
    !$omp        cutoff, cutoff_sqr,                         &
    !$omp        bp_list, num_bp_calc, base_type, chain_id,  &
    !$omp        param_bp_theta_1, param_bp_theta_2,         &
    !$omp        param_bp_theta_3, param_bp_phi_1,           &
    !$omp        param_bp_sigma, param_bp_epsilon,           &
    !$omp        param_cs_1_epsilon, param_cs_2_epsilon,     &
    !$omp        param_cs_1_sigma, param_cs_2_sigma,         &
    !$omp        param_cs_1_theta_cs, param_cs_2_theta_cs,   &
    !$omp        param_bp_alpha, param_bp_K,                 &
    !$omp        param_cs_alpha, param_cs_K,                 &
    !$omp        bp_theta_threshold_1, cs_theta_threshold_1, &
    !$omp        bp_theta_threshold_2, cs_theta_threshold_2) &
    !$omp reduction(+:virial) reduction(+:epair)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1

    do i = 1, natom-1

      ini_bp = num_bp + 1
      fin_bp = num_bp + num_bp_calc(i,id)
      num_bp = fin_bp

      ! if (mod(i-1,nproc_city*nthread) .ne. my_id) cycle

      do k = ini_bp, fin_bp

        i_B2 = i
        i_B4 = bp_list(k,id)

        d42(1:3) = coord(1:3, i_B2) - coord(1:3, i_B4)
        r42_sqr = d42(1) * d42(1) + d42(2) * d42(2) + d42(3) * d42(3)

        if (r42_sqr > cutoff_sqr) then
          cycle
        end if

        e_tmp_bp = 0.0_wp
        e_tmp_cs = 0.0_wp

        force_S1(1:3) = 0.0_wp
        force_S3(1:3) = 0.0_wp
        force_B2(1:3) = 0.0_wp
        force_B4(1:3) = 0.0_wp
        force_B5(1:3) = 0.0_wp
        force_B6(1:3) = 0.0_wp

        !==============================================================================
        !               Base Pairing Energy / Force Calculation                       !
        !==============================================================================
        !
        !  ____                   ____       _      _
        ! | __ )  __ _ ___  ___  |  _ \ __ _(_)_ __(_)_ __   __ _
        ! |  _ \ / _` / __|/ _ \ | |_) / _` | | '__| | '_ \ / _` |
        ! | |_) | (_| \__ \  __/ |  __/ (_| | | |  | | | | | (_| |
        ! |____/ \__,_|___/\___| |_|   \__,_|_|_|  |_|_| |_|\__, |
        !                                                   |___/
        !
        i_S1 = i_B2 - 1
        i_S3 = i_B4 - 1

        type_B2      = base_type        (i_B2)
        type_B4      = base_type        (i_B4)
        bp_theta_1_0 = param_bp_theta_1 (type_B2)
        bp_theta_2_0 = param_bp_theta_2 (type_B2)
        bp_theta_3_0 = param_bp_theta_3 (type_B2)
        bp_phi_1_0   = param_bp_phi_1   (type_B2)
        bp_sigma     = param_bp_sigma   (type_B2)
        bp_epsilon   = param_bp_epsilon (type_B2)

        ! write (*,*) " i: ", i_B2, " j: ", i_B4, type_B2, type_B4, bp_theta_1_0 / RAD, bp_theta_2_0 / RAD, bp_phi_1_0 / RAD
        ! write (*,*) " i: ", i_B2, " j: ", i_B4, type_B2, type_B4, bp_sigma, bp_epsilon / JOU2CAL

        ! -----------------
        ! 1 -- 2 <== 4 -- 3
        ! -----------------
        !
        d24(1:3) = coord(1:3, i_B4) - coord(1:3, i_B2)
        r42      = sqrt(r42_sqr)
        r42_inv  = 1.0_wp / r42
        e42(1:3) = d42(1:3) * r42_inv
        e24(1:3) = d24(1:3) * r42_inv

        ! -----------------
        ! 1 <== 2 -- 4 -- 3
        ! -----------------
        !
        d21(1:3) = coord(1:3, i_S1) - coord(1:3, i_B2)
        r21_sqr  = d21(1) * d21(1) + d21(2) * d21(2) + d21(3) * d21(3)
        r21      = sqrt(r21_sqr)
        r21_inv  = 1.0_wp / r21
        e21(1:3) = d21(1:3) * r21_inv

        ! -----------------
        ! 1 -- 2 -- 4 ==> 3
        ! -----------------
        !
        d43(1:3) = coord(1:3, i_S3) - coord(1:3, i_B4)
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
        angle_124 = acos(cos_124)

        ! -----------------------
        ! Angle B4: 2 -- ~4~ -- 3
        ! -----------------------
        !
        cos_342 = e43(1) * e42(1) + e43(2) * e42(2) + e43(3) * e42(3)
        if (cos_342 >  1.0_wp) cos_342 =  1.0_wp
        if (cos_342 < -1.0_wp) cos_342 = -1.0_wp
        sin_342 = sqrt(1.0_wp - cos_342 * cos_342)
        angle_342 = acos(cos_342)

        ! -------------------------
        ! Dihedral 1 -- 2 == 4 -- 3
        ! -------------------------
        !
        aindex(1) = i_S1
        aindex(2) = i_B2
        aindex(3) = i_B4
        aindex(4) = i_S3
        call calculate_nonlocal_dihedral(aindex, coord, cos_dih_24, sin_dih_24, grad_dih_24)
        ene_coef_phi  =   0.5_wp * (1.0_wp + cos_dih_24 * cos(bp_phi_1_0) + sin_dih_24 * sin(bp_phi_1_0))
        grad_coef_phi = - 0.5_wp * (         sin_dih_24 * cos(bp_phi_1_0) - cos_dih_24 * sin(bp_phi_1_0))

        ! write (*,*) "id=", my_id, " i: ", i_B2, " j: ", i_B4, cos_dih_24, sin_dih_24

        ! ============================================================
        ! basepairing interaction: energy/force calculation @@@@@@@...
        ! ============================================================
        !
        if (r42 < bp_sigma) then

          call calculate_repulsive(d24, r42, bp_sigma, param_bp_alpha, bp_epsilon, ene_repl, grad_repl)

          e_tmp_bp = e_tmp_bp + ene_repl

          force_repl_tmp(1:3) = - grad_repl(1:3)
          force_B2(1:3)       = force_B2(1:3) + force_repl_tmp(1:3)
          force_B4(1:3)       = force_B4(1:3) - force_repl_tmp(1:3)

        end if

        delta_angle_124     = angle_124 - bp_theta_1_0
        abs_delta_angle_124 = abs(delta_angle_124)
        delta_angle_342     = angle_342 - bp_theta_2_0
        abs_delta_angle_342 = abs(delta_angle_342)

        call calculate_attractive(d24, r42, bp_sigma, param_bp_alpha, bp_epsilon, ene_attr, grad_attr)

        if (abs_delta_angle_124 <= bp_theta_threshold_1) then
          if (abs_delta_angle_342 <= bp_theta_threshold_1) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_attr

            force_dih_24_tmp(1:9) = - grad_coef_phi * ene_attr * grad_dih_24(1:9)
            force_S1(1:3)         = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)         = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)         = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)         = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_attr_tmp(1:3)   = - ene_coef_phi * grad_attr(1:3)
            force_B2(1:3)         = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)         = force_B4(1:3) - force_attr_tmp(1:3)

          else if (abs_delta_angle_342 <= bp_theta_threshold_2) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, grad_angle_342)
            var_tmp           = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_2 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_2 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_2 * ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_angle_342_tmp(1:3) - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) + force_angle_342_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_2
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          end if
        else if (abs_delta_angle_124 <= bp_theta_threshold_2) then

          call calculate_angle(d21, d24, r21, r42, cos_124, sin_124, grad_angle_124)
          var_tmp = sin(param_bp_K * delta_angle_124)
          ene_coef_theta_1  = var_tmp * var_tmp
          grad_coef_theta_1 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_124)

          if (abs_delta_angle_342 <= bp_theta_threshold_1) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_angle_124_tmp(1:3) - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) + force_angle_124_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          else if (abs_delta_angle_342 <= bp_theta_threshold_2) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, grad_angle_342)
            var_tmp = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 * ene_coef_theta_2 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 * ene_coef_theta_2 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 * ene_coef_theta_2 * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_angle_124_tmp(1:3) - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) + force_angle_124_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 * grad_coef_theta_2 * ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_angle_342_tmp(1:3) - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) + force_angle_342_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 * ene_coef_theta_2
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
        i_B6 = i_B2 + 3
        i_B5 = i_B4 - 3

        ! ---------------------------
        ! Angle 1234: 1 -- 2 - 3 -- 4
        ! ---------------------------
        !
        cos_1243 = e21(1) * e43(1) + e21(2) * e43(2) + e21(3) * e43(3)
        if (cos_1243 >  1.0_wp) cos_1243 =  1.0_wp
        if (cos_1243 < -1.0_wp) cos_1243 = -1.0_wp
        sin_1243   = sqrt(1.0_wp - cos_1243 * cos_1243)
        angle_1243 = acos(cos_1243)

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
        if (i_B5 > 0 .and. chain_id(i_B5) == chain_id(i_B4)) then

          type_B5       = base_type           (i_B5)
          cs_1_epsilon  = param_cs_1_epsilon  (type_B2, type_B5)
          cs_1_sigma    = param_cs_1_sigma    (type_B2, type_B5)
          cs_1_theta_cs_0 = param_cs_1_theta_cs (type_B2, type_B5)

          ! -------
          ! 2 ==> 5
          ! -------
          !
          d25(1:3) = coord(1:3, i_B5) - coord(1:3, i_B2)
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

          delta_angle_125     = angle_125 - cs_1_theta_cs_0
          abs_delta_angle_125 = abs(delta_angle_125)

          call calculate_attractive(d25, r25, cs_1_sigma, param_cs_alpha, cs_1_epsilon, ene_attr, grad_attr)

          if (abs_delta_angle_125 < cs_theta_threshold_1) then

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B2(1:3)  = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)  = force_B5(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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

          else if (abs_delta_angle_125 <= cs_theta_threshold_2) then

            call calculate_angle(d21, d25, r21, r25, cos_125, sin_125, grad_angle_125)
            var_tmp = sin(param_cs_K * delta_angle_125)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_125)

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_125_tmp(1:6) = - force_coef_tmp * grad_angle_125(1:6)
              force_S1(1:3)            = force_S1(1:3) + force_angle_125_tmp(1:3)
              force_B2(1:3)            = force_B2(1:3) - force_angle_125_tmp(1:3) - force_angle_125_tmp(4:6)
              force_B5(1:3)            = force_B5(1:3) + force_angle_125_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)            = force_B5(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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
              force_B2(1:3)             = force_B2(1:3) - force_angle_125_tmp(1:3) - force_angle_125_tmp(4:6)
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
        if (i_B6 <= natom .and. chain_id(i_B6) == chain_id(i_B2)) then

          type_B6       = base_type           (i_B6)
          cs_2_epsilon  = param_cs_2_epsilon  (type_B4, type_B6)
          cs_2_sigma    = param_cs_2_sigma    (type_B4, type_B6)
          cs_2_theta_cs_0 = param_cs_2_theta_cs (type_B4, type_B6)

          ! -------
          ! 6 <== 4
          ! -------
          !
          d46(1:3) = coord(1:3, i_B6) - coord(1:3, i_B4)
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

          delta_angle_346     = angle_346 - cs_2_theta_cs_0
          abs_delta_angle_346 = abs(delta_angle_346)

          call calculate_attractive(d46, r46, cs_2_sigma, param_cs_alpha, cs_2_epsilon, ene_attr, grad_attr)

          if (abs_delta_angle_346 < cs_theta_threshold_1) then

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B4(1:3)  = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)  = force_B6(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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

          else if (abs_delta_angle_346 <= cs_theta_threshold_2) then

            call calculate_angle(d43, d46, r43, r46, cos_346, sin_346, grad_angle_346)
            var_tmp = sin(param_cs_K * delta_angle_346)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_346)

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_346_tmp(1:6) = - force_coef_tmp * grad_angle_346(1:6)
              force_S3(1:3)            = force_S3(1:3) + force_angle_346_tmp(1:3)
              force_B4(1:3)            = force_B4(1:3) - force_angle_346_tmp(1:3) - force_angle_346_tmp(4:6)
              force_B6(1:3)            = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B4(1:3)            = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)            = force_B6(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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
              force_B4(1:3)             = force_B4(1:3) - force_angle_346_tmp(1:3) - force_angle_346_tmp(4:6)
              force_B6(1:3)             = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * ene_coef_theta_cs
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B4(1:3)             = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)             = force_B6(1:3) - force_attr_tmp(1:3)

            end if

          end if

        end if

        ! Calc energy
        epair = epair + e_tmp_bp + e_tmp_cs

        ! store force
        !
        force(1:3,i_S1,id) = force(1:3,i_S1,id) + force_S1(1:3)
        force(1:3,i_B2,id) = force(1:3,i_B2,id) + force_B2(1:3)
        force(1:3,i_B4,id) = force(1:3,i_B4,id) + force_B4(1:3)
        force(1:3,i_S3,id) = force(1:3,i_S3,id) + force_S3(1:3)
        force(1:3,i_B5,id) = force(1:3,i_B5,id) + force_B5(1:3)
        force(1:3,i_B6,id) = force(1:3,i_B6,id) + force_B6(1:3)

        ! virial
        !
        ! do l = 1, 3
        !   virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        ! end do

        do l = 1, 3
          do m = l+1, 3
            var_tmp =force_S1(l) * coord(m, i_S1) &
                + force_B2(l) * coord(m, i_B2) &
                + force_B4(l) * coord(m, i_B4) &
                + force_S3(l) * coord(m, i_S3) &
                + force_B5(l) * coord(m, i_B5) &
                + force_B6(l) * coord(m, i_B6)
            virial(m, l) = virial(m, l) - var_tmp
            virial(l, m) = virial(l, m) - var_tmp
          end do
          var_tmp =force_S1(l) * coord(l, i_S1) &
              + force_B2(l) * coord(l, i_B2) &
              + force_B4(l) * coord(l, i_B4) &
              + force_S3(l) * coord(l, i_S3) &
              + force_B5(l) * coord(l, i_B5) &
              + force_B6(l) * coord(l, i_B6)
          virial(l,l) = virial(l,l) - var_tmp
        end do

      end do

    end do
    !$omp end parallel

    call timer(TimerBasePair, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_base_pairing

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_pairing_pbc
  !> @brief        calculate base-pairing energy with pairlist (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : PBC coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epair    : base-pairing energy
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_pairing_pbc(enefunc, boundary, pairlist, &
                                                 coord, force, virial, epair)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epair

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: natom, nbase
    integer           :: num_bp, ini_bp, fin_bp
    integer           :: ibase
    integer           :: i, k, l, m
    logical           :: cg_infinite_DNA

    integer           :: i1, i2, i3, k_pbc ! for PBC
    integer           :: i_S1, i_S3, i_B2, i_B4, i_B5, i_B6
    integer           :: type_B2, type_B4, type_B5, type_B6
    integer           :: aindex(1:4)
    real(wp)          :: cutoff, cutoff_sqr
    real(wp)          :: dij_pbc(1:3)
    real(wp)          :: d21(3), r21_sqr, r21, r21_inv, e21(3)
    real(wp)          :: d42(3), r42_sqr, r42, r42_inv, e42(3)
    real(wp)          :: d24(3), e24(3)
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

    integer,  pointer :: bp_list(:,:)
    integer,  pointer :: num_bp_calc(:,:)
    integer,  pointer :: base_type(:)
    integer,  pointer :: chain_id(:)
    integer,  pointer :: cg_list_base(:)
    real(wp), pointer :: param_bp_theta_1(:)
    real(wp), pointer :: param_bp_theta_2(:)
    real(wp), pointer :: param_bp_theta_3(:)
    real(wp), pointer :: param_bp_phi_1(:)
    real(wp), pointer :: param_bp_sigma(:)
    real(wp), pointer :: param_bp_epsilon(:)
    real(wp), pointer :: param_cs_1_epsilon(:,:)
    real(wp), pointer :: param_cs_1_sigma(:,:)
    real(wp), pointer :: param_cs_1_theta_cs(:,:)
    real(wp), pointer :: param_cs_2_epsilon(:,:)
    real(wp), pointer :: param_cs_2_sigma(:,:)
    real(wp), pointer :: param_cs_2_theta_cs(:,:)

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

    real(wp)          :: bsize(3), half_bsize(3)

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

    natom                = size(coord(1,:))
    nbase                = enefunc%num_cg_particle_DNA_base
    cutoff               = enefunc%cg_cutoffdist_DNAbp
    cutoff_sqr           = cutoff * cutoff

    bp_list              => pairlist%cg_DNA_basepair_list
    num_bp_calc          => pairlist%num_cg_DNA_basepair_calc
    base_type            => enefunc%NA_base_type
    chain_id             => enefunc%mol_chain_id
    cg_list_base         => enefunc%cg_particle_DNA_base
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

    bsize(1)             = boundary%box_size_x
    bsize(2)             = boundary%box_size_y
    bsize(3)             = boundary%box_size_z
    half_bsize(1:3)      = 0.5_wp * bsize(1:3)

    num_bp               = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                             &
    !$omp firstprivate(num_bp)                               &
    !$omp private(my_id, id,                                 &
    !$omp         ini_bp, fin_bp,                            &
    !$omp         i, k, l, m,                                &
    !$omp         i1, i2, i3, k_pbc,                         &
    !$omp         i_S1, i_S3, i_B2, i_B4, i_B5, i_B6,        &
    !$omp         type_B2, type_B4, type_B5, type_B6,        &
    !$omp         aindex, ibase, dij_pbc,                    &
    !$omp         d21, r21_sqr, r21, r21_inv, e21,           &
    !$omp         d42, r42_sqr, r42, r42_inv, e42, d24, e24, &
    !$omp         d43, r43_sqr, r43, r43_inv, e43,           &
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
    !$omp         e_tmp_bp, e_tmp_cs                         &
    !$omp         )                                          &
    !$omp shared(coord,                                      &
    !$omp        force, my_city_rank, nproc_city,            &
    !$omp        nthread, natom, nbase,                      &
    !$omp        cutoff, cutoff_sqr, cg_infinite_DNA,        &
    !$omp        bp_list, num_bp_calc,                       &
    !$omp        base_type, chain_id, cg_list_base,          &
    !$omp        param_bp_theta_1, param_bp_theta_2,         &
    !$omp        param_bp_theta_3, param_bp_phi_1,           &
    !$omp        param_bp_sigma, param_bp_epsilon,           &
    !$omp        param_cs_1_epsilon, param_cs_2_epsilon,     &
    !$omp        param_cs_1_sigma, param_cs_2_sigma,         &
    !$omp        param_cs_1_theta_cs, param_cs_2_theta_cs,   &
    !$omp        param_bp_alpha, param_bp_K,                 &
    !$omp        param_cs_alpha, param_cs_K,                 &
    !$omp        bp_theta_threshold_1, cs_theta_threshold_1, &
    !$omp        bp_theta_threshold_2, cs_theta_threshold_2, &
    !$omp        bsize, half_bsize)                          &
    !$omp reduction(+:virial) reduction(+:epair)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1

    do ibase = 1, nbase - 1

      ini_bp = num_bp + 1
      fin_bp = num_bp + num_bp_calc(ibase,id)
      num_bp = fin_bp

      i = cg_list_base(ibase)

      do k = ini_bp, fin_bp

        i_B2 = i

        ! i_B4 = bp_list(k,id)
        k_pbc = bp_list(k,id)
        i_B4 = k_pbc / 27
        k_pbc = k_pbc - i_B4*27
        i3 = k_pbc / 9
        k_pbc = k_pbc - i3*9
        i2 = k_pbc / 3
        i1 = k_pbc - i2*3
        i1 = i1 - 1
        i2 = i2 - 1
        i3 = i3 - 1

        d42(1)  = coord(1, i_B2) - coord(1, i_B4) - bsize(1) * real(i1, wp)
        d42(2)  = coord(2, i_B2) - coord(2, i_B4) - bsize(2) * real(i2, wp)
        d42(3)  = coord(3, i_B2) - coord(3, i_B4) - bsize(3) * real(i3, wp)
        dij_pbc(1:3) = coord(1:3, i_B2) - coord(1:3, i_B4) - d42(1:3)

        r42_sqr = d42(1) * d42(1) + d42(2) * d42(2) + d42(3) * d42(3)

        if (r42_sqr > cutoff_sqr) then
          cycle
        end if

        e_tmp_bp = 0.0_wp
        e_tmp_cs = 0.0_wp

        force_S1(1:3) = 0.0_wp
        force_S3(1:3) = 0.0_wp
        force_B2(1:3) = 0.0_wp
        force_B4(1:3) = 0.0_wp
        force_B5(1:3) = 0.0_wp
        force_B6(1:3) = 0.0_wp

        !==============================================================================
        !               Base Pairing Energy / Force Calculation                       !
        !==============================================================================
        !
        !  ____                   ____       _      _
        ! | __ )  __ _ ___  ___  |  _ \ __ _(_)_ __(_)_ __   __ _
        ! |  _ \ / _` / __|/ _ \ | |_) / _` | | '__| | '_ \ / _` |
        ! | |_) | (_| \__ \  __/ |  __/ (_| | | |  | | | | | (_| |
        ! |____/ \__,_|___/\___| |_|   \__,_|_|_|  |_|_| |_|\__, |
        !                                                   |___/
        !
        i_S1 = i_B2 - 1
        i_S3 = i_B4 - 1

        type_B2      = base_type        (i_B2)
        type_B4      = base_type        (i_B4)
        bp_theta_1_0 = param_bp_theta_1 (type_B2)
        bp_theta_2_0 = param_bp_theta_2 (type_B2)
        bp_theta_3_0 = param_bp_theta_3 (type_B2)
        bp_phi_1_0   = param_bp_phi_1   (type_B2)
        bp_sigma     = param_bp_sigma   (type_B2)
        bp_epsilon   = param_bp_epsilon (type_B2)

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
        d21(1:3) = coord(1:3, i_S1) - coord(1:3, i_B2)
        if (d21(1) > half_bsize(1)) then
          d21(1) = d21(1) - bsize(1)
        else if (d21(1) < -half_bsize(1)) then
          d21(1) = d21(1) + bsize(1)
        end if
        if (d21(2) > half_bsize(2)) then
          d21(2) = d21(2) - bsize(2)
        else if (d21(2) < -half_bsize(2)) then
          d21(2) = d21(2) + bsize(2)
        end if
        if (d21(3) > half_bsize(3)) then
          d21(3) = d21(3) - bsize(3)
        else if (d21(3) < -half_bsize(3)) then
          d21(3) = d21(3) + bsize(3)
        end if
        !
        r21_sqr  = d21(1) * d21(1) + d21(2) * d21(2) + d21(3) * d21(3)
        r21      = sqrt(r21_sqr)
        r21_inv  = 1.0_wp / r21
        e21(1:3) = d21(1:3) * r21_inv

        ! -----------------
        ! 1 -- 2 -- 4 ==> 3
        ! -----------------
        !
        d43(1:3) = coord(1:3, i_S3) - coord(1:3, i_B4)
        if (d43(1) > half_bsize(1)) then
          d43(1) = d43(1) - bsize(1)
        else if (d43(1) < -half_bsize(1)) then
          d43(1) = d43(1) + bsize(1)
        end if
        if (d43(2) > half_bsize(2)) then
          d43(2) = d43(2) - bsize(2)
        else if (d43(2) < -half_bsize(2)) then
          d43(2) = d43(2) + bsize(2)
        end if
        if (d43(3) > half_bsize(3)) then
          d43(3) = d43(3) - bsize(3)
        else if (d43(3) < -half_bsize(3)) then
          d43(3) = d43(3) + bsize(3)
        end if
        !
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
        angle_124 = acos(cos_124)

        ! -----------------------
        ! Angle B4: 2 -- ~4~ -- 3
        ! -----------------------
        !
        cos_342 = e43(1) * e42(1) + e43(2) * e42(2) + e43(3) * e42(3)
        if (cos_342 >  1.0_wp) cos_342 =  1.0_wp
        if (cos_342 < -1.0_wp) cos_342 = -1.0_wp
        sin_342 = sqrt(1.0_wp - cos_342 * cos_342)
        angle_342 = acos(cos_342)

        ! -------------------------
        ! Dihedral 1 -- 2 == 4 -- 3
        ! -------------------------
        !
        aindex(1) = i_S1
        aindex(2) = i_B2
        aindex(3) = i_B4
        aindex(4) = i_S3
        call calculate_nonlocal_dihedral_pbc(aindex, coord, bsize, &
            cos_dih_24, sin_dih_24, grad_dih_24)
        ene_coef_phi  =   0.5_wp * (1.0_wp + cos_dih_24 * cos(bp_phi_1_0) &
            + sin_dih_24 * sin(bp_phi_1_0))
        grad_coef_phi = - 0.5_wp * (         sin_dih_24 * cos(bp_phi_1_0) &
            - cos_dih_24 * sin(bp_phi_1_0))

        ! ============================================================
        ! basepairing interaction: energy/force calculation @@@@@@@...
        ! ============================================================
        !
        if (r42 < bp_sigma) then

          call calculate_repulsive(d24, r42, bp_sigma, param_bp_alpha, bp_epsilon, ene_repl, grad_repl)

          e_tmp_bp = e_tmp_bp + ene_repl

          force_repl_tmp(1:3) = - grad_repl(1:3)
          force_B2(1:3)       = force_B2(1:3) + force_repl_tmp(1:3)
          force_B4(1:3)       = force_B4(1:3) - force_repl_tmp(1:3)

        end if

        delta_angle_124     = angle_124 - bp_theta_1_0
        abs_delta_angle_124 = abs(delta_angle_124)
        delta_angle_342     = angle_342 - bp_theta_2_0
        abs_delta_angle_342 = abs(delta_angle_342)

        call calculate_attractive(d24, r42, bp_sigma, param_bp_alpha, bp_epsilon, ene_attr, grad_attr)

        if (abs_delta_angle_124 <= bp_theta_threshold_1) then
          if (abs_delta_angle_342 <= bp_theta_threshold_1) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_attr

            force_dih_24_tmp(1:9) = - grad_coef_phi * ene_attr * grad_dih_24(1:9)
            force_S1(1:3)         = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)         = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)         = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)         = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_attr_tmp(1:3)   = - ene_coef_phi * grad_attr(1:3)
            force_B2(1:3)         = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)         = force_B4(1:3) - force_attr_tmp(1:3)

          else if (abs_delta_angle_342 <= bp_theta_threshold_2) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, grad_angle_342)
            var_tmp           = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_2 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_2 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_2 * ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_angle_342_tmp(1:3) - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) + force_angle_342_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_2
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          end if

        else if (abs_delta_angle_124 <= bp_theta_threshold_2) then

          call calculate_angle(d21, d24, r21, r42, cos_124, sin_124, grad_angle_124)
          var_tmp = sin(param_bp_K * delta_angle_124)
          ene_coef_theta_1  = var_tmp * var_tmp
          grad_coef_theta_1 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_124)

          if (abs_delta_angle_342 <= bp_theta_threshold_1) then

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_angle_124_tmp(1:3) - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) + force_angle_124_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1
            force_attr_tmp(1:3)      = - force_coef_tmp * grad_attr(1:3)
            force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_attr_tmp(1:3)

          else if (abs_delta_angle_342 <= bp_theta_threshold_2) then

            call calculate_angle(d42, d43, r42, r43, cos_342, sin_342, grad_angle_342)
            var_tmp = sin(param_bp_K * delta_angle_342)
            ene_coef_theta_2  = var_tmp * var_tmp
            grad_coef_theta_2 = param_bp_K * sin(2.0_wp * param_bp_K * delta_angle_342)

            e_tmp_bp = e_tmp_bp + ene_coef_phi * ene_coef_theta_1 * ene_coef_theta_2 * ene_attr

            force_coef_tmp           = grad_coef_phi * ene_coef_theta_1 * ene_coef_theta_2 * ene_attr
            force_dih_24_tmp(1:9)    = - force_coef_tmp * grad_dih_24(1:9)
            force_S1(1:3)            = force_S1(1:3) + force_dih_24_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_dih_24_tmp(1:3) + force_dih_24_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) - force_dih_24_tmp(4:6) - force_dih_24_tmp(7:9)
            force_S3(1:3)            = force_S3(1:3) + force_dih_24_tmp(7:9)

            force_coef_tmp           = ene_coef_phi * grad_coef_theta_1 * ene_coef_theta_2 * ene_attr
            force_angle_124_tmp(1:6) = - force_coef_tmp * grad_angle_124(1:6)
            force_S1(1:3)            = force_S1(1:3) + force_angle_124_tmp(1:3)
            force_B2(1:3)            = force_B2(1:3) - force_angle_124_tmp(1:3) - force_angle_124_tmp(4:6)
            force_B4(1:3)            = force_B4(1:3) + force_angle_124_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 * grad_coef_theta_2 * ene_attr
            force_angle_342_tmp(1:6) = - force_coef_tmp * grad_angle_342(1:6)
            force_B2(1:3)            = force_B2(1:3) + force_angle_342_tmp(1:3)
            force_B4(1:3)            = force_B4(1:3) - force_angle_342_tmp(1:3) - force_angle_342_tmp(4:6)
            force_S3(1:3)            = force_S3(1:3) + force_angle_342_tmp(4:6)

            force_coef_tmp           = ene_coef_phi * ene_coef_theta_1 * ene_coef_theta_2
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
        i_B6 = i_B2 + 3
        i_B5 = i_B4 - 3

        if (cg_infinite_DNA) then
          if (i_B6 > natom .or. chain_id(i_B6) /= chain_id(i_B2)) then
            i_B6 = i_B2 - 3
            do while(.true.)
              i_B6 = i_B6 - 3
              if (i_B6 <= 0) exit
              if (chain_id(i_B6) /= chain_id(i_B2)) exit
            end do
            i_B6 = i_B6 + 3
          end if
          if (i_B5 <= 0 .or. chain_id(i_B5) /= chain_id(i_B4)) then
            i_B5 = i_B4 + 3
            do while(.true.)
              i_B5 = i_B5 + 3
              if (i_B5 > natom) exit
              if (chain_id(i_B5) /= chain_id(i_B4)) exit
            end do
            i_B5 = i_B5 - 3
          end if
        end if

        ! ---------------------------
        ! Angle 1234: 1 -- 2 - 3 -- 4
        ! ---------------------------
        !
        cos_1243 = e21(1) * e43(1) + e21(2) * e43(2) + e21(3) * e43(3)
        if (cos_1243 >  1.0_wp) cos_1243 =  1.0_wp
        if (cos_1243 < -1.0_wp) cos_1243 = -1.0_wp
        sin_1243   = sqrt(1.0_wp - cos_1243 * cos_1243)
        angle_1243 = acos(cos_1243)

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
        if (i_B5 > 0 .and. chain_id(i_B5) == chain_id(i_B4)) then

          type_B5       = base_type           (i_B5)
          cs_1_epsilon  = param_cs_1_epsilon  (type_B2, type_B5)
          cs_1_sigma    = param_cs_1_sigma    (type_B2, type_B5)
          cs_1_theta_cs_0 = param_cs_1_theta_cs (type_B2, type_B5)

          ! -------
          ! 2 ==> 5
          ! -------
          !
          d25(1:3) = coord(1:3, i_B5) - coord(1:3, i_B2)
          if (d25(1) > half_bsize(1)) then
            d25(1) = d25(1) - bsize(1)
          else if (d25(1) < -half_bsize(1)) then
            d25(1) = d25(1) + bsize(1)
          end if
          if (d25(2) > half_bsize(2)) then
            d25(2) = d25(2) - bsize(2)
          else if (d25(2) < -half_bsize(2)) then
            d25(2) = d25(2) + bsize(2)
          end if
          if (d25(3) > half_bsize(3)) then
            d25(3) = d25(3) - bsize(3)
          else if (d25(3) < -half_bsize(3)) then
            d25(3) = d25(3) + bsize(3)
          end if
          !
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

          delta_angle_125     = angle_125 - cs_1_theta_cs_0
          abs_delta_angle_125 = abs(delta_angle_125)

          call calculate_attractive(d25, r25, cs_1_sigma, param_cs_alpha, cs_1_epsilon, ene_attr, grad_attr)

          if (abs_delta_angle_125 < cs_theta_threshold_1) then

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B2(1:3)  = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)  = force_B5(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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

          else if (abs_delta_angle_125 <= cs_theta_threshold_2) then

            call calculate_angle(d21, d25, r21, r25, cos_125, sin_125, grad_angle_125)
            var_tmp = sin(param_cs_K * delta_angle_125)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_125)

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_125_tmp(1:6) = - force_coef_tmp * grad_angle_125(1:6)
              force_S1(1:3)            = force_S1(1:3) + force_angle_125_tmp(1:3)
              force_B2(1:3)            = force_B2(1:3) - force_angle_125_tmp(1:3) - force_angle_125_tmp(4:6)
              force_B5(1:3)            = force_B5(1:3) + force_angle_125_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B2(1:3)            = force_B2(1:3) + force_attr_tmp(1:3)
              force_B5(1:3)            = force_B5(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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
              force_B2(1:3)             = force_B2(1:3) - force_angle_125_tmp(1:3) - force_angle_125_tmp(4:6)
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
        if (i_B6 <= natom .and. chain_id(i_B6) == chain_id(i_B2)) then

          type_B6       = base_type           (i_B6)
          cs_2_epsilon  = param_cs_2_epsilon  (type_B4, type_B6)
          cs_2_sigma    = param_cs_2_sigma    (type_B4, type_B6)
          cs_2_theta_cs_0 = param_cs_2_theta_cs (type_B4, type_B6)

          ! -------
          ! 6 <== 4
          ! -------
          !
          d46(1:3) = coord(1:3, i_B6) - coord(1:3, i_B4)
          if (d46(1) > half_bsize(1)) then
            d46(1) = d46(1) - bsize(1)
          else if (d46(1) < -half_bsize(1)) then
            d46(1) = d46(1) + bsize(1)
          end if
          if (d46(2) > half_bsize(2)) then
            d46(2) = d46(2) - bsize(2)
          else if (d46(2) < -half_bsize(2)) then
            d46(2) = d46(2) + bsize(2)
          end if
          if (d46(3) > half_bsize(3)) then
            d46(3) = d46(3) - bsize(3)
          else if (d46(3) < -half_bsize(3)) then
            d46(3) = d46(3) + bsize(3)
          end if
          !
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

          delta_angle_346     = angle_346 - cs_2_theta_cs_0
          abs_delta_angle_346 = abs(delta_angle_346)

          call calculate_attractive(d46, r46, cs_2_sigma, param_cs_alpha, cs_2_epsilon, ene_attr, grad_attr)

          if (abs_delta_angle_346 < cs_theta_threshold_1) then

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_attr

              force_attr_tmp(1:3) = - grad_attr(1:3)
              force_B4(1:3)  = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)  = force_B6(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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

          else if (abs_delta_angle_346 <= cs_theta_threshold_2) then

            call calculate_angle(d43, d46, r43, r46, cos_346, sin_346, grad_angle_346)
            var_tmp = sin(param_cs_K * delta_angle_346)
            ene_coef_theta_cs  = var_tmp * var_tmp
            grad_coef_theta_cs = param_cs_K * sin(2.0_wp * param_cs_K * delta_angle_346)

            if (abs_delta_angle_1243 <= bp_theta_threshold_1) then

              e_tmp_cs = e_tmp_cs + ene_coef_theta_cs * ene_attr

              force_coef_tmp           = grad_coef_theta_cs * ene_attr
              force_angle_346_tmp(1:6) = - force_coef_tmp * grad_angle_346(1:6)
              force_S3(1:3)            = force_S3(1:3) + force_angle_346_tmp(1:3)
              force_B4(1:3)            = force_B4(1:3) - force_angle_346_tmp(1:3) - force_angle_346_tmp(4:6)
              force_B6(1:3)            = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_attr_tmp(1:3)      = - ene_coef_theta_cs * grad_attr(1:3)
              force_B4(1:3)            = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)            = force_B6(1:3) - force_attr_tmp(1:3)

            else if (abs_delta_angle_1243 <= bp_theta_threshold_2) then

              call calculate_angle(d21, d43, r21, r43, cos_1243, sin_1243, grad_angle_1243)
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
              force_B4(1:3)             = force_B4(1:3) - force_angle_346_tmp(1:3) - force_angle_346_tmp(4:6)
              force_B6(1:3)             = force_B6(1:3) + force_angle_346_tmp(4:6)

              force_coef_tmp            = ene_coef_theta_3 * ene_coef_theta_cs
              force_attr_tmp(1:3)       = - force_coef_tmp * grad_attr(1:3)
              force_B4(1:3)             = force_B4(1:3) + force_attr_tmp(1:3)
              force_B6(1:3)             = force_B6(1:3) - force_attr_tmp(1:3)

            end if

          end if

        end if

        ! Calc energy
        epair = epair + e_tmp_bp + e_tmp_cs

        ! store force
        !
        force(1:3,i_S1,id) = force(1:3,i_S1,id) + force_S1(1:3)
        force(1:3,i_B2,id) = force(1:3,i_B2,id) + force_B2(1:3)
        force(1:3,i_B4,id) = force(1:3,i_B4,id) + force_B4(1:3)
        force(1:3,i_S3,id) = force(1:3,i_S3,id) + force_S3(1:3)
        force(1:3,i_B5,id) = force(1:3,i_B5,id) + force_B5(1:3)
        force(1:3,i_B6,id) = force(1:3,i_B6,id) + force_B6(1:3)

        ! virial
        !
        do l = 1, 3
          do m = l+1, 3
            var_tmp =force_S1(l) * coord(m, i_S1) &
                + force_B2(l) * coord(m, i_B2) &
                + force_B4(l) * coord(m, i_B4) &
                + force_S3(l) * coord(m, i_S3) &
                + force_B5(l) * coord(m, i_B5) &
                + force_B6(l) * coord(m, i_B6)
            virial(m, l) = virial(m, l) - var_tmp
            virial(l, m) = virial(l, m) - var_tmp
          end do
          var_tmp =force_S1(l) * coord(l, i_S1) &
              + force_B2(l) * coord(l, i_B2) &
              + force_B4(l) * coord(l, i_B4) &
              + force_S3(l) * coord(l, i_S3) &
              + force_B5(l) * coord(l, i_B5) &
              + force_B6(l) * coord(l, i_B6)
          virial(l,l) = virial(l,l) - var_tmp
        end do

      end do

    end do
    !$omp end parallel

    call timer(TimerBasePair, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_base_pairing_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_exv
  !> @brief        calculate DNA excluded volume energy with pairlist (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eexv     : DNA excluded volume energy
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_exv(enefunc, pairlist, coord, force, virial, eexv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eexv

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_exv, ini_exv, fin_exv
    integer                   :: i, j, k, l, m
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: dij(3), rij_sqr
    real(wp)                  :: inv_rij_sqr, inv_rij_6th
    real(wp)                  :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)                  :: grad_coef_exv, grad(3)
    real(wp)                  :: eexv_tmp

    integer,  pointer         :: exv_list(:,:)
    integer,  pointer         :: num_exv_calc(:,:)

    integer,  pointer         :: base_type(:)
    real(wp), pointer         :: param_exv_sigma(:,:)
    real(wp)                  :: param_exv_epsilon


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDNAexv, TimerOn)

    natom             = size(coord(1,:))

    exv_list          => pairlist%cg_DNA_exv_list
    num_exv_calc      => pairlist%num_cg_DNA_exv_calc

    base_type         => enefunc%NA_base_type
    param_exv_sigma   => enefunc%cgDNA_exv_sigma
    param_exv_epsilon = enefunc%cgDNA_exv_epsilon

    num_exv           = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                         &
    !$omp firstprivate(num_exv)                          &
    !$omp private(my_id, id,                             &
    !$omp         ini_exv, fin_exv,                      &
    !$omp         i, k, j, l, m,                         &
    !$omp         sigma, sigma_sqr, sigma_6th,           &
    !$omp         dij, rij_sqr,                          &
    !$omp         inv_rij_sqr, inv_rij_6th,              &
    !$omp         sig_over_rij_6th, sig_over_rij_12th,   &
    !$omp         grad_coef_exv, grad, eexv_tmp          &
    !$omp         )                                      &
    !$omp shared(coord, force, my_city_rank, nproc_city, &
    !$omp        nthread, natom,                         &
    !$omp        exv_list, num_exv_calc, base_type,      &
    !$omp        param_exv_sigma, param_exv_epsilon      &
    !$omp        )                                       &
    !$omp reduction(+:virial) reduction(+:eexv)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i = 1, natom-1

      ini_exv = num_exv + 1
      fin_exv = num_exv + num_exv_calc(i,id)
      num_exv = fin_exv

      ! if (mod(i-1,nproc_city*nthread) .ne. my_id) cycle

      do k = ini_exv, fin_exv

        j = exv_list(k,id)

        sigma     = param_exv_sigma(base_type(i), base_type(j))
        sigma_sqr = sigma * sigma

        dij(1:3)  = coord(1:3,j) - coord(1:3,i)
        rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij_sqr > sigma_sqr) cycle

        sigma_6th         = sigma_sqr * sigma_sqr * sigma_sqr

        inv_rij_sqr       = 1.0_wp / rij_sqr
        inv_rij_6th       = inv_rij_sqr * inv_rij_sqr * inv_rij_sqr

        sig_over_rij_6th  = sigma_6th * inv_rij_6th
        sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

        eexv_tmp          = param_exv_epsilon * (sig_over_rij_12th - 2 * sig_over_rij_6th + 1)
        eexv              = eexv + eexv_tmp

        ! gradient
        !
        grad_coef_exv = - inv_rij_sqr * 12.0_wp * param_exv_epsilon * (sig_over_rij_6th - sig_over_rij_12th)
        grad(1:3)     = grad_coef_exv * dij(1:3)

        ! store force
        !
        force(1:3,i,id)  = force(1:3,i,id) - grad(1:3)
        force(1:3,j,id)  = force(1:3,j,id) + grad(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
        end do

      end do
    end do
    !$omp end parallel

    call timer(TimerCGDNAexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_exv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_exv_pbc
  !> @brief        calculate DNA excluded volume energy with pairlist (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eexv     : DNA excluded volume energy
  !! @note         3SPN.2C
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_exv_pbc(enefunc, boundary, pairlist, &
      coord, force, virial, eexv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eexv

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: n_dna, i_dna
    integer           :: num_exv, ini_exv, fin_exv
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc
    real(wp)          :: sigma, sigma_sqr
    real(wp)          :: dij(3), rij_sqr
    real(wp)          :: inv_rij_sqr, sig_over_rij_sqr
    real(wp)          :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)          :: grad_coef_exv, grad(3)
    real(wp)          :: eexv_tmp
    !
    real(wp)            :: coord_tmp(3)
    integer             :: basetype_i, basetype_j
    real(wp)            :: force_tmp(3), factor
    real(wp)            :: virial_tmp(3,3)
    real(wp)            :: virial_omp(3,3,nthread)
    real(wp)            :: eexv_omp(nthread), eexv_omp_tmp
    !
    integer,  pointer :: cg_list_dna(:)
    integer,  pointer :: exv_list(:,:)
    integer,  pointer :: num_exv_calc(:,:)
    !
    integer,  pointer :: base_type(:)
    real(wp), pointer :: param_exv_sigma(:,:)
    real(wp)          :: param_exv_epsilon
    !
    real(wp)          :: bsize(3)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDNAexv, TimerOn)

    n_dna              = enefunc%num_cg_particle_DNA_all
    !
    exv_list           => pairlist%cg_DNA_exv_list
    num_exv_calc       => pairlist%num_cg_DNA_exv_calc
    !
    base_type          => enefunc%NA_base_type
    param_exv_sigma    => enefunc%cgDNA_exv_sigma
    param_exv_epsilon  = enefunc%cgDNA_exv_epsilon
    cg_list_dna        => enefunc%cg_particle_DNA_all
    !
    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    !
    num_exv            = 0

    eexv_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp firstprivate(num_exv)                        &
    !$omp private(my_id, id,                           &
    !$omp         ini_exv, fin_exv,                    &
    !$omp         i, k, j, l, m,                       &
    !$omp         i1, i2, i3, k_pbc,                   &
    !$omp         sigma, sigma_sqr,                    &
    !$omp         dij, rij_sqr,                        &
    !$omp         inv_rij_sqr, sig_over_rij_sqr,       &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         coord_tmp, force_tmp, factor,        &
    !$omp         basetype_i, basetype_j,              &
    !$omp         virial_tmp, eexv_omp_tmp,            &
    !$omp         grad_coef_exv, grad, eexv_tmp        &
    !$omp         )                                    &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, n_dna, cg_list_dna,          &
    !$omp        exv_list, num_exv_calc, base_type,    &
    !$omp        param_exv_sigma, param_exv_epsilon,   &
    !$omp        virial_omp, eexv_omp,                 &
    !$omp        bsize,            force               &
    !$omp        )                                     &
    !$omp reduction(+:virial) reduction(+:eexv)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i_dna = 1, n_dna - 1

      ini_exv = num_exv + 1
      fin_exv = num_exv + num_exv_calc(i_dna,id)
      num_exv = fin_exv

      if (fin_exv >= ini_exv) then

        i = cg_list_dna(i_dna)

        coord_tmp(1:3)      = coord(1:3, i)
        basetype_i          = base_type(i)
        force_tmp(1:3)      = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp
        eexv_omp_tmp        = 0.0_wp

        do k = ini_exv, fin_exv

          k_pbc = exv_list(k, id)
          j     = k_pbc / 27
          k_pbc = k_pbc - j*27
          i3    = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2    = k_pbc / 3
          i1    = k_pbc - i2*3
          i1    = i1 - 1
          i2    = i2 - 1
          i3    = i3 - 1

          basetype_j          = base_type(j)

          sigma     = param_exv_sigma(basetype_i, basetype_j)
          sigma_sqr = sigma * sigma

          dij(1)  = coord(1,j) - coord_tmp(1) + bsize(1)*real(i1,wp)
          dij(2)  = coord(2,j) - coord_tmp(2) + bsize(2)*real(i2,wp)
          dij(3)  = coord(3,j) - coord_tmp(3) + bsize(3)*real(i3,wp)
          rij_sqr  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij_sqr > sigma_sqr) cycle

          inv_rij_sqr       = 1.0_wp / rij_sqr
          sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr

          sig_over_rij_6th  = sig_over_rij_sqr * sig_over_rij_sqr * sig_over_rij_sqr
          sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

          eexv_tmp = param_exv_epsilon * (sig_over_rij_12th - 2 * sig_over_rij_6th + 1)
          eexv_omp_tmp = eexv_omp_tmp + eexv_tmp

          ! gradient
          !
          grad_coef_exv = - inv_rij_sqr * 12.0_wp * param_exv_epsilon * (sig_over_rij_6th - sig_over_rij_12th)
          grad(1)       = grad_coef_exv * dij(1)
          grad(2)       = grad_coef_exv * dij(2)
          grad(3)       = grad_coef_exv * dij(3)

          ! store force
          !
          force_tmp(1) = force_tmp(1) - grad(1)
          force_tmp(2) = force_tmp(2) - grad(2)
          force_tmp(3) = force_tmp(3) - grad(3)
          force(1,j,id) = force(1,j,id) + grad(1)
          force(2,j,id) = force(2,j,id) + grad(2)
          force(3,j,id) = force(3,j,id) + grad(3)

          ! virial
          !
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        eexv_omp(id) = eexv_omp(id) + eexv_omp_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)

      end if

    end do
    !$omp end parallel

    do i = 1, nthread
      eexv = eexv + eexv_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGDNAexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_DNA_exv_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_ele
  !> @brief        calculate electrostatic energy with pairlist (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eele     : electrostatic energy
  !! @note         3SPN.2C + protein
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_ele(enefunc, pairlist, coord, force, virial, &
                                   eele)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eele

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_ele, ini_ele, fin_ele
    integer                   :: i, j, k, l, m
    real(wp)                  :: r_tmp
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: dij(3), rij_sqr, rij, inv_rij
    real(wp)                  :: debye_length
    real(wp)                  :: inv_debye_length
    real(wp)                  :: ele_coef
    real(wp)                  :: diele_const
    real(wp)                  :: inv_diele_const
    real(wp)                  :: ele_tmp_e_T, ele_tmp_a_C
    real(wp)                  :: ele_tmp_sol_T, ele_tmp_sol_C
    real(wp)                  :: grad_coef_ele, grad(3)
    real(wp)                  :: e_tmp_ele

    integer,  pointer         :: ele_list(:,:)
    integer,  pointer         :: num_ele_calc(:,:)
    real(wp), pointer         :: ele_scaling(:,:)

    real(wp), pointer         :: charge(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    natom            = size(coord(1,:))

    ele_list         => pairlist%cg_ele_list
    ele_scaling      => pairlist%cg_ele_scaling_list
    num_ele_calc     => pairlist%num_cg_ele_calc

    charge           => enefunc%cg_charge

    ! set parameters
    !
    cutoff           = enefunc%cg_cutoffdist_ele
    cutoff_sqr       = cutoff * cutoff

    ! electrostatic parameters
    ele_tmp_sol_T = enefunc%cg_ele_sol_T
    ele_tmp_sol_C = enefunc%cg_ele_sol_IC
    ele_tmp_e_T   = 2.494e2_wp - 7.88e-1_wp * ele_tmp_sol_T &
        + 7.2e-4_wp * ele_tmp_sol_T * ele_tmp_sol_T
    ele_tmp_a_C   = 1.0e0_wp - 2.551e-1_wp * ele_tmp_sol_C  &
        + 5.151e-2_wp * ele_tmp_sol_C * ele_tmp_sol_C       &
        - 6.889e-3_wp * ele_tmp_sol_C * ele_tmp_sol_C * ele_tmp_sol_C
    diele_const = ele_tmp_e_T * ele_tmp_a_C

    debye_length = 1.0e10_wp                     &
        * sqrt(                                  &
        (CAL2JOU * ELECTRIC_CONST                &
        * diele_const * KBOLTZ * ele_tmp_sol_T)  &
        /                                        &
        (2.0_wp * AVOGADRO * AVOGADRO            &
        * ELEMENT_CHARGE * ELEMENT_CHARGE        &
        * ele_tmp_sol_C))
    inv_debye_length = 1.0_wp / debye_length

    ele_coef         = enefunc%cg_ele_coef
    inv_diele_const  = 1.0_wp / diele_const

    num_ele          = 0

    ! write (*,*) cutoff, inv_debye_length, ele_coef, inv_diele_const

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                  &
    !$omp firstprivate(num_ele)                   &
    !$omp private(my_id, id,                      &
    !$omp         ini_ele, fin_ele,               &
    !$omp         i, k, j, l, m,                  &
    !$omp         r_tmp,                          &
    !$omp         dij, rij_sqr, rij, inv_rij,     &
    !$omp         grad_coef_ele, grad, e_tmp_ele  &
    !$omp         )                               &
    !$omp shared(coord, my_city_rank, nproc_city, &
    !$omp        nthread, natom,                  &
    !$omp        cutoff, cutoff_sqr,              &
    !$omp        inv_debye_length,                &
    !$omp        ele_coef, inv_diele_const,       &
    !$omp        ele_list, num_ele_calc,          &
    !$omp        ele_scaling,                     &
    !$omp        charge, force                    &
    !$omp        )                                &
    !$omp reduction(+:virial) reduction(+:eele)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do i = 1, natom-1

      ini_ele = num_ele + 1
      fin_ele = num_ele + num_ele_calc(i,id)
      num_ele = fin_ele

      ! if (mod(i-1,nproc_city*nthread) .ne. my_id) cycle

      do k = ini_ele, fin_ele

        j = ele_list(k,id)

        dij(1:3)  = coord(1:3,j) - coord(1:3,i)
        rij_sqr = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij_sqr > cutoff_sqr) cycle

        rij   = sqrt(rij_sqr)
        inv_rij = 1.0_wp / rij

        ! energy
        !
        e_tmp_ele = ele_coef * charge(i) * charge(j) &
            * inv_diele_const * inv_rij              &
            * exp(-rij * inv_debye_length)           &
            * ele_scaling(k, id)

        eele = eele + e_tmp_ele

        ! gradient
        !
        grad_coef_ele = e_tmp_ele * (inv_debye_length + inv_rij) * inv_rij
        grad(1:3)     = grad_coef_ele * dij(1:3)

        ! store force
        !
        force(1:3,i,id)  = force(1:3,i,id) - grad(1:3)
        force(1:3,j,id)  = force(1:3,j,id) + grad(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
        end do

      end do
    end do
    !$omp end parallel

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_ele

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_ele_pbc
  !> @brief        calculate CG electrostatic energy with pairlist (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eele     : electrostatic energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_ele_pbc(enefunc, boundary, pairlist, &
                                       coord, force, virial, eele)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eele

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: n_charge
    integer           :: num_ele, ini_ele, fin_ele
    integer           :: i, j, l, i1, i2, i3, k
    integer           :: i_atom, j_atom
    logical           :: proceed
    real(wp)          :: r_tmp(3), qtmp
    real(wp)          :: cutoff, cutoff_sqr
    real(wp)          :: dij(3), rij_sqr, rij, inv_rij
    real(wp)          :: debye_length
    real(wp)          :: inv_debye_length
    real(wp)          :: ele_coef
    real(wp)          :: diele_const
    real(wp)          :: inv_diele_const
    real(wp)          :: ele_tmp_e_T, ele_tmp_a_C
    real(wp)          :: ele_tmp_sol_T, ele_tmp_sol_C
    real(wp)          :: grad_coef_ele, grad(3)
    real(wp)          :: e_tmp_ele
    !
    real(wp)          :: bsize(3), inv_bsize(3)
    real(wp)          :: force_temp(3), factor
    real(wp)          :: virial_omp(3,3,nthread), virial_tmp(3,3)
    real(wp)          :: eele_omp(nthread), eele_tmp
    !
    real(wp), pointer :: charge(:)
    integer,  pointer :: cg_list_charged(:)
    integer,  pointer :: ele_list(:,:)
    integer,  pointer :: num_ele_calc(:,:)
    real(wp), pointer :: ele_scaling(:,:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGDebye, TimerOn)

    n_charge         = enefunc%num_cg_particle_charged

    ele_list         => pairlist%cg_ele_list
    ele_scaling      => pairlist%cg_ele_scaling_list
    num_ele_calc     => pairlist%num_cg_ele_calc

    charge           => enefunc%cg_charge
    cg_list_charged  => enefunc%cg_particle_charged

    bsize(1)         =  boundary%box_size_x
    bsize(2)         =  boundary%box_size_y
    bsize(3)         =  boundary%box_size_z
    inv_bsize(1:3)   = 1.0_wp/bsize(1:3)

    ! set parameters
    !
    cutoff           = enefunc%cg_cutoffdist_ele
    cutoff_sqr       = cutoff * cutoff

    ! electrostatic parameters
    ele_tmp_sol_T = enefunc%cg_ele_sol_T
    ele_tmp_sol_C = enefunc%cg_ele_sol_IC
    ele_tmp_e_T   = 2.494e2_wp - 7.88e-1_wp * ele_tmp_sol_T &
        + 7.2e-4_wp * ele_tmp_sol_T * ele_tmp_sol_T
    ele_tmp_a_C   = 1.0e0_wp - 2.551e-1_wp * ele_tmp_sol_C  &
        + 5.151e-2_wp * ele_tmp_sol_C * ele_tmp_sol_C       &
        - 6.889e-3_wp * ele_tmp_sol_C * ele_tmp_sol_C * ele_tmp_sol_C
    diele_const = ele_tmp_e_T * ele_tmp_a_C

    debye_length = 1.0e10_wp                     &
        * sqrt(                                  &
        (CAL2JOU * ELECTRIC_CONST                &
        * diele_const * KBOLTZ * ele_tmp_sol_T)  &
        /                                        &
        (2.0_wp * AVOGADRO * AVOGADRO            &
        * ELEMENT_CHARGE * ELEMENT_CHARGE        &
        * ele_tmp_sol_C))
    inv_debye_length = 1.0_wp / debye_length

    ele_coef         = enefunc%cg_ele_coef
    inv_diele_const  = 1.0_wp / diele_const

    num_ele          = 0

    eele_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    !$omp parallel default(none)                  &
    !$omp firstprivate(num_ele)                   &
    !$omp private(my_id, id, force_temp,          &
    !$omp         ini_ele, fin_ele,               &
    !$omp         i, j, l, k, i1, i2, i3,         &
    !$omp         i_atom, j_atom,                 &
    !$omp         r_tmp, factor, qtmp,            &
    !$omp         dij, rij_sqr, rij, inv_rij,     &
    !$omp         grad_coef_ele, grad, e_tmp_ele, &
    !$omp         proceed, eele_tmp, virial_tmp)  &
    !$omp shared(coord, my_city_rank, nproc_city, &
    !$omp        nthread, n_charge,               &
    !$omp        cutoff, cutoff_sqr,              &
    !$omp        inv_debye_length,                &
    !$omp        ele_coef, inv_diele_const,       &
    !$omp        ele_list, num_ele_calc,          &
    !$omp        ele_scaling,                     &
    !$omp        charge, cg_list_charged,         &
    !$omp        bsize, inv_bsize, force,         &
    !$omp        virial_omp, eele_omp, my_country_rank)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do i = 1, n_charge - 1

      proceed = .false.

      ini_ele = num_ele + 1
      fin_ele = num_ele + num_ele_calc(i,id)
      num_ele = fin_ele

      if (fin_ele >= ini_ele) proceed = .true.

      if (proceed) then

        i_atom = cg_list_charged(i)

        r_tmp(1:3) = coord(1:3,i_atom)
        qtmp = charge(i_atom)
        force_temp(1:3) = 0.0_wp
        eele_tmp = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp

        !#ifdef IFORT19
        !dir$ simd
        !#endif
        do j = ini_ele, fin_ele

          k = ele_list(j,id)
          j_atom = k / 27
          k = k - j_atom*27
          i3 = k / 9
          k = k - i3*9
          i2 = k / 3
          i1 = k - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          dij(1)  = coord(1,j_atom) - r_tmp(1) + bsize(1)*real(i1,wp)
          dij(2)  = coord(2,j_atom) - r_tmp(2) + bsize(2)*real(i2,wp)
          dij(3)  = coord(3,j_atom) - r_tmp(3) + bsize(3)*real(i3,wp)
          rij_sqr = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          factor = merge(0.0_wp,1.0_wp,rij_sqr>cutoff_sqr)
          ! if (rij_sqr > cutoff_sqr) cycle

          rij = sqrt(rij_sqr)
          inv_rij = 1.0_wp / rij

          ! energy
          !
          e_tmp_ele = ele_coef                  &
              * qtmp * charge(j_atom) &
              * inv_diele_const * inv_rij       &
              * exp(-rij * inv_debye_length)    &
              * ele_scaling(j, id)*factor

          !       eele = eele + e_tmp_ele
          !       eele_omp(id) = eele_omp(id) + e_tmp_ele
          eele_tmp = eele_tmp + e_tmp_ele

          ! gradient
          !
          grad_coef_ele = e_tmp_ele * (inv_debye_length + inv_rij) * inv_rij
          grad(1)     = grad_coef_ele * dij(1)
          grad(2)     = grad_coef_ele * dij(2)
          grad(3)     = grad_coef_ele * dij(3)

          ! store force
          !
          !       force(1:3,i_atom) = force(1:3,i_atom) - grad(1:3)
          force_temp(1) = force_temp(1) - grad(1)
          force_temp(2) = force_temp(2) - grad(2)
          force_temp(3) = force_temp(3) - grad(3)
          force(1,j_atom,id) = force(1,j_atom,id) + grad(1)
          force(2,j_atom,id) = force(2,j_atom,id) + grad(2)
          force(3,j_atom,id) = force(3,j_atom,id) + grad(3)

          ! virial
          !
          !       do l = 1, 3
          !         virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
          !         virial_omp(1:3, l,id) = virial_omp(1:3, l,id) - dij(1:3) * grad(l)
          !         virial_tmp(1:3, l) = virial_tmp(1:3, l) - dij(1:3) * grad(l)
          !       end do
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do

        force(1:3,i_atom,id) = force(1:3,i_atom,id) + force_temp(1:3)
        eele_omp(id) = eele_omp(id) + eele_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)

      end if
    end do
    !$omp end parallel

    do i = 1, nthread
      eele = eele + eele_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGDebye, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_ele_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_HPS
  !> @brief        calculate IDR interaction in the HPS model (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ehps     : energy of IDR HPS interactions
  !! @note         IDR HPS model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_HPS(enefunc, pairlist, coord, force, virial, ehps)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ehps

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_hps, ini_hps, fin_hps
    integer                   :: i, j, k, l, m
    real(wp)                  :: lambda, epsilon
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: dij(3), rij_sqr
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: inv_rij_sqr, inv_rij_6th
    real(wp)                  :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)                  :: grad_coef_hps, grad(3)
    real(wp)                  :: ehps_tmp

    integer,  pointer         :: hps_list(:,:)
    integer,  pointer         :: num_hps_calc(:,:)

    real(wp), pointer         :: hps_sigma_half(:)
    real(wp), pointer         :: hps_lambda_half(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRHPS, TimerOn)

    natom            = size(coord(1,:))

    hps_list         => pairlist%cg_IDR_HPS_list
    num_hps_calc     => pairlist%num_cg_IDR_HPS_calc

    hps_sigma_half   => enefunc%cg_IDR_HPS_sigma_half
    hps_lambda_half  => enefunc%cg_IDR_HPS_lambda_half
    epsilon          =  enefunc%cg_IDR_HPS_epsilon

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

    num_hps          = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp firstprivate(num_hps)                        &
    !$omp private(my_id, id,                           &
    !$omp         ini_hps, fin_hps,                    &
    !$omp         i, k, j, l, m,                       &
    !$omp         lambda,                              &
    !$omp         sigma, sigma_sqr, sigma_6th,         &
    !$omp         dij, rij_sqr,                        &
    !$omp         inv_rij_sqr, inv_rij_6th,            &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         grad_coef_hps, grad, ehps_tmp        &
    !$omp         )                                    &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, natom, epsilon,              &
    !$omp        cutoff, cutoff_sqr,                   &
    !$omp        hps_list, num_hps_calc,               &
    !$omp        hps_sigma_half, hps_lambda_half,      &
    !$omp        force)                                &
    !$omp reduction(+:virial) reduction(+:ehps)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i = 1, natom-1

      ini_hps = num_hps + 1
      fin_hps = num_hps + num_hps_calc(i,id)
      num_hps = fin_hps

      do k = ini_hps, fin_hps

        j = hps_list(k,id)

        lambda    = hps_lambda_half(i) + hps_lambda_half(j)
        sigma     = hps_sigma_half(i) + hps_sigma_half(j)
        sigma_sqr = sigma * sigma

        dij(1:3)  = coord(1:3,j) - coord(1:3,i)
        rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij_sqr > cutoff_sqr) cycle

        sigma_6th         = sigma_sqr * sigma_sqr * sigma_sqr

        inv_rij_sqr       = 1.0_wp / rij_sqr
        inv_rij_6th       = inv_rij_sqr * inv_rij_sqr * inv_rij_sqr

        sig_over_rij_6th  = sigma_6th * inv_rij_6th
        sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

        ehps_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

        if (sig_over_rij_6th >= 0.5_wp) then
          ehps = ehps + ehps_tmp + (1.0_wp - lambda) * epsilon
          grad_coef_hps = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        else
          ehps = ehps + lambda * ehps_tmp
          grad_coef_hps = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        end if

        ! gradient
        !
        grad(1:3)     = grad_coef_hps * dij(1:3)

        ! store force
        !
        force(1:3,i,id)  = force(1:3,i,id) - grad(1:3)
        force(1:3,j,id)  = force(1:3,j,id) + grad(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
        end do

      end do
    end do
    !$omp end parallel

    call timer(TimerCGIDRHPS, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_HPS

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_HPS_pbc
  !> @brief        calculate IDR interaction in the HPS model (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ehps     : energy of IDR HPS interactions
  !! @note         IDR HPS model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_HPS_pbc(enefunc, boundary, pairlist, &
                                           coord, force, virial, ehps)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ehps

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: n_idr, i_idr
    integer           :: num_hps, ini_hps, fin_hps
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc
    real(wp)          :: coor_i_tmp(3)
    real(wp)          :: hps_sigma_half_i_tmp, hps_lambda_half_i_tmp
    real(wp)          :: lambda, epsilon
    real(wp)          :: sigma, sigma_sqr, sigma_6th
    real(wp)          :: dij(3), rij_sqr
    real(wp)          :: cutoff, cutoff_sqr
    real(wp)          :: inv_rij_sqr, inv_rij_6th
    real(wp)          :: sig_over_rij_sqr
    real(wp)          :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)          :: grad_coef_hps, grad(3)
    real(wp)          :: ehps_tmp
    !
    real(wp)          :: bsize(3), inv_bsize(3)
    !
    real(wp)          :: force_tmp(3), factor
    real(wp)          :: virial_tmp(3,3)
    real(wp)          :: virial_omp(3,3,nthread)
    real(wp)          :: ehps_omp(nthread), ehps_omp_tmp
    !
    integer,  pointer :: cg_list_idr(:)
    integer,  pointer :: hps_list(:,:)
    integer,  pointer :: num_hps_calc(:,:)
    !
    real(wp), pointer :: hps_sigma_half(:)
    real(wp), pointer :: hps_lambda_half(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRHPS, TimerOn)

    n_idr           =  enefunc%num_cg_particle_IDR_HPS
    cg_list_idr     => enefunc%cg_particle_IDR_HPS

    hps_list        => pairlist%cg_IDR_HPS_list
    num_hps_calc    => pairlist%num_cg_IDR_HPS_calc

    hps_sigma_half  => enefunc%cg_IDR_HPS_sigma_half
    hps_lambda_half => enefunc%cg_IDR_HPS_lambda_half
    epsilon         =  enefunc%cg_IDR_HPS_epsilon

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    inv_bsize(1:3)  = 1.0_wp/bsize(1:3)

    num_hps         = 0

    ehps_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp firstprivate(num_hps)                        &
    !$omp private(my_id, id,                           &
    !$omp         ini_hps, fin_hps,                    &
    !$omp         i, k, j, l, m, i_idr,                &
    !$omp         i1, i2, i3, k_pbc,                   &
    !$omp         lambda,                              &
    !$omp         coor_i_tmp,                          &
    !$omp         hps_sigma_half_i_tmp,                &
    !$omp         hps_lambda_half_i_tmp,               &
    !$omp         sigma, sigma_sqr, sigma_6th,         &
    !$omp         dij, rij_sqr,                        &
    !$omp         inv_rij_sqr, sig_over_rij_sqr,       &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         force_tmp, factor,                   &
    !$omp         virial_tmp, ehps_omp_tmp,            &
    !$omp         grad_coef_hps, grad, ehps_tmp        &
    !$omp         )                                    &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, n_idr, epsilon,              &
    !$omp        cutoff, cutoff_sqr,                   &
    !$omp        hps_list, num_hps_calc,               &
    !$omp        hps_sigma_half, hps_lambda_half,      &
    !$omp        cg_list_idr, bsize, inv_bsize,        &
    !$omp        virial_omp, ehps_omp,                 &
    !$omp        force)                                &
    !$omp reduction(+:virial) reduction(+:ehps)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i_idr = 1, n_idr - 1

      ini_hps = num_hps + 1
      fin_hps = num_hps + num_hps_calc(i_idr, id)
      num_hps = fin_hps

      if (fin_hps >= ini_hps) then

        i = cg_list_idr(i_idr)

        coor_i_tmp(1:3) = coord(1:3, i)
        hps_sigma_half_i_tmp = hps_sigma_half(i)
        hps_lambda_half_i_tmp = hps_lambda_half(i)
        force_tmp(1:3)      = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp
        ehps_omp_tmp        = 0.0_wp

        do k = ini_hps, fin_hps

          k_pbc = hps_list(k,id)
          j = k_pbc / 27
          k_pbc = k_pbc - j*27
          i3 = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2 = k_pbc / 3
          i1 = k_pbc - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          lambda    = hps_lambda_half_i_tmp + hps_lambda_half(j)
          sigma     = hps_sigma_half_i_tmp + hps_sigma_half(j)
          sigma_sqr = sigma * sigma

          dij(1)  = coord(1,j) - coor_i_tmp(1) + bsize(1) * real(i1, wp)
          dij(2)  = coord(2,j) - coor_i_tmp(2) + bsize(2) * real(i2, wp)
          dij(3)  = coord(3,j) - coor_i_tmp(3) + bsize(3) * real(i3, wp)
          rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij_sqr > cutoff_sqr) cycle

          inv_rij_sqr       = 1.0_wp / rij_sqr
          sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr

          sig_over_rij_6th  = sig_over_rij_sqr * sig_over_rij_sqr * sig_over_rij_sqr
          sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

          ehps_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

          if (sig_over_rij_6th >= 0.5_wp) then
            ehps_omp_tmp = ehps_omp_tmp + ehps_tmp + (1.0_wp - lambda) * epsilon
            grad_coef_hps = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          else
            ehps_omp_tmp = ehps_omp_tmp + lambda * ehps_tmp
            grad_coef_hps = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          end if

          ! gradient
          !
          grad(1) = grad_coef_hps * dij(1)
          grad(2) = grad_coef_hps * dij(2)
          grad(3) = grad_coef_hps * dij(3)

          ! store force
          !
          force_tmp(1) = force_tmp(1) - grad(1)
          force_tmp(2) = force_tmp(2) - grad(2)
          force_tmp(3) = force_tmp(3) - grad(3)
          force(1,j,id) = force(1,j,id) + grad(1)
          force(2,j,id) = force(2,j,id) + grad(2)
          force(3,j,id) = force(3,j,id) + grad(3)

          ! virial
          !
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        ehps_omp(id) = ehps_omp(id) + ehps_omp_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)
      end if

    end do                      ! i_idr
    !$omp end parallel

    do i = 1, nthread
      ehps = ehps + ehps_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGIDRHPS, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_HPS_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_KH
  !> @brief        calculate IDR interaction in the KH model (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ekh      : energy of IDR KH interactions
  !! @note         IDR KH model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_KH(enefunc, pairlist, coord, force, virial, ekh)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ekh

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_kh, ini_kh, fin_kh
    integer                   :: i, j, k, l, m
    real(wp)                  :: lambda, epsilon
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: dij(3), rij_sqr
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: inv_rij_sqr, inv_rij_6th
    real(wp)                  :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)                  :: grad_coef_kh, grad(3)
    real(wp)                  :: ekh_tmp

    integer,  pointer         :: kh_list(:,:)
    integer,  pointer         :: num_kh_calc(:,:)

    integer,  pointer         :: atom_type(:)
    real(wp), pointer         :: kh_sigma_half(:)
    real(wp), pointer         :: kh_epsilon(:, :)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRKH, TimerOn)

    natom            = size(coord(1,:))

    kh_list         => pairlist%cg_IDR_KH_list
    num_kh_calc     => pairlist%num_cg_IDR_KH_calc

    atom_type       => enefunc%atom_cls
    kh_sigma_half   => enefunc%cg_IDR_KH_sigma_half
    kh_epsilon      => enefunc%cg_IDR_KH_epsilon_D

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

    num_kh          = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp firstprivate(num_kh)                         &
    !$omp private(my_id, id,                           &
    !$omp         ini_kh, fin_kh,                      &
    !$omp         i, k, j, l, m,                       &
    !$omp         lambda,  epsilon,                    &
    !$omp         sigma, sigma_sqr, sigma_6th,         &
    !$omp         dij, rij_sqr,                        &
    !$omp         inv_rij_sqr, inv_rij_6th,            &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         grad_coef_kh, grad, ekh_tmp          &
    !$omp         )                                    &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, natom,                       &
    !$omp        cutoff, cutoff_sqr,                   &
    !$omp        kh_list, num_kh_calc, atom_type,      &
    !$omp        kh_sigma_half, kh_epsilon,            &
    !$omp        force)                                &
    !$omp reduction(+:virial) reduction(+:ekh)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i = 1, natom-1

      ini_kh = num_kh + 1
      fin_kh = num_kh + num_kh_calc(i,id)
      num_kh = fin_kh

      do k = ini_kh, fin_kh

        j = kh_list(k,id)

        sigma     = kh_sigma_half(i) + kh_sigma_half(j)
        sigma_sqr = sigma * sigma
        sigma_6th = sigma_sqr * sigma_sqr * sigma_sqr

        epsilon   = kh_epsilon(atom_type(i), atom_type(j))
        if (epsilon > 0) then
          lambda  = -1.0_wp
        else
          lambda  =  1.0_wp
          epsilon = -epsilon
        end if

        dij(1:3)  = coord(1:3,j) - coord(1:3,i)
        rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij_sqr > cutoff_sqr) cycle

        inv_rij_sqr       = 1.0_wp / rij_sqr
        inv_rij_6th       = inv_rij_sqr * inv_rij_sqr * inv_rij_sqr

        sig_over_rij_6th  = sigma_6th * inv_rij_6th
        sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

        ekh_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

        if (sig_over_rij_6th >= 0.5_wp) then
          ekh = ekh + ekh_tmp + (1.0_wp - lambda) * epsilon
          grad_coef_kh = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        else
          ekh = ekh + lambda * ekh_tmp
          grad_coef_kh = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        end if

        ! gradient
        !
        grad(1:3)     = grad_coef_kh * dij(1:3)

        ! store force
        !
        force(1:3,i,id)  = force(1:3,i,id) - grad(1:3)
        force(1:3,j,id)  = force(1:3,j,id) + grad(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
        end do

      end do
    end do
    !$omp end parallel

    call timer(TimerCGIDRKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_IDR_KH_pbc
  !> @brief        calculate IDR interaction in the KH model (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ekh     : energy of IDR KH interactions
  !! @note         IDR KH model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_IDR_KH_pbc(enefunc, boundary, pairlist, &
      coord, force, virial, ekh)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ekh

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: n_idr, i_idr
    integer           :: num_kh, ini_kh, fin_kh
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc
    integer           :: atomtype_i, atomtype_j
    real(wp)          :: coor_i_tmp(3)
    real(wp)          :: lambda, epsilon
    real(wp)          :: si_tmp ! sigma of i
    real(wp)          :: sigma, sigma_sqr, sigma_6th
    real(wp)          :: dij(3), rij_sqr
    real(wp)          :: cutoff, cutoff_sqr
    real(wp)          :: inv_rij_sqr, inv_rij_6th
    real(wp)          :: sig_over_rij_sqr
    real(wp)          :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)          :: grad_coef_kh, grad(3)
    real(wp)          :: ekh_tmp
    !
    real(wp)          :: bsize(3), inv_bsize(3)
    !
    real(wp)          :: force_tmp(3), factor
    real(wp)          :: virial_tmp(3,3)
    real(wp)          :: virial_omp(3,3,nthread)
    real(wp)          :: ekh_omp(nthread), ekh_omp_tmp
    !
    integer,  pointer :: cg_list_idr(:)
    integer,  pointer :: kh_list(:,:)
    integer,  pointer :: num_kh_calc(:,:)
    !
    integer,  pointer :: atom_type(:)
    real(wp), pointer :: kh_sigma_half(:)
    real(wp), pointer :: kh_epsilon(:, :)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGIDRKH, TimerOn)

    n_idr          =  enefunc%num_cg_particle_IDR_KH
    cg_list_idr    => enefunc%cg_particle_IDR_KH

    kh_list        => pairlist%cg_IDR_KH_list
    num_kh_calc    => pairlist%num_cg_IDR_KH_calc

    atom_type      => enefunc%atom_cls
    kh_sigma_half  => enefunc%cg_IDR_KH_sigma_half
    kh_epsilon     => enefunc%cg_IDR_KH_epsilon_D

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    inv_bsize(1:3) = 1.0_wp/bsize(1:3)

    num_kh         = 0

    ekh_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                       &
    !$omp firstprivate(num_kh)                         &
    !$omp private(my_id, id,                           &
    !$omp         ini_kh, fin_kh,                      &
    !$omp         i, k, j, l, m, i_idr,                &
    !$omp         i1, i2, i3, k_pbc,                   &
    !$omp         lambda, epsilon, si_tmp,             &
    !$omp         sigma, sigma_sqr, sigma_6th,         &
    !$omp         atomtype_i, atomtype_j,              &
    !$omp         dij, rij_sqr, coor_i_tmp,            &
    !$omp         inv_rij_sqr, sig_over_rij_sqr,       &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         force_tmp, factor,                   &
    !$omp         virial_tmp, ekh_omp_tmp,             &
    !$omp         grad_coef_kh, grad, ekh_tmp          &
    !$omp         )                                    &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, n_idr,                       &
    !$omp        cutoff, cutoff_sqr,                   &
    !$omp        kh_list, num_kh_calc, atom_type,      &
    !$omp        kh_sigma_half, kh_epsilon,            &
    !$omp        cg_list_idr, bsize, inv_bsize,        &
    !$omp        virial_omp, ekh_omp,                  &
    !$omp        force)                                &
    !$omp reduction(+:virial) reduction(+:ekh)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i_idr = 1, n_idr - 1

      ini_kh = num_kh + 1
      fin_kh = num_kh + num_kh_calc(i_idr, id)
      num_kh = fin_kh

      if (fin_kh >= ini_kh) then

        i = cg_list_idr(i_idr)

        coor_i_tmp(1:3) = coord(1:3, i)
        si_tmp = kh_sigma_half(i)
        atomtype_i = atom_type(i)
        force_tmp(1:3)      = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp
        ekh_omp_tmp        = 0.0_wp

        do k = ini_kh, fin_kh

          k_pbc = kh_list(k,id)
          j = k_pbc / 27
          k_pbc = k_pbc - j*27
          i3 = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2 = k_pbc / 3
          i1 = k_pbc - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          sigma     = si_tmp + kh_sigma_half(j)
          sigma_sqr = sigma * sigma

          atomtype_j = atom_type(j)
          epsilon   = kh_epsilon(atomtype_i, atomtype_j)
          if (epsilon > 0) then
            lambda  = -1.0_wp
          else
            lambda  =  1.0_wp
            epsilon = -epsilon
          end if

          dij(1)  = coord(1,j) - coor_i_tmp(1) + bsize(1) * real(i1, wp)
          dij(2)  = coord(2,j) - coor_i_tmp(2) + bsize(2) * real(i2, wp)
          dij(3)  = coord(3,j) - coor_i_tmp(3) + bsize(3) * real(i3, wp)
          rij_sqr = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij_sqr > cutoff_sqr) cycle

          inv_rij_sqr       = 1.0_wp / rij_sqr
          sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr

          sig_over_rij_6th  = sig_over_rij_sqr * sig_over_rij_sqr * sig_over_rij_sqr
          sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

          ekh_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

          if (sig_over_rij_6th >= 0.5_wp) then
            ekh_omp_tmp = ekh_omp_tmp + ekh_tmp + (1.0_wp - lambda) * epsilon
            grad_coef_kh = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          else
            ekh_omp_tmp = ekh_omp_tmp + lambda * ekh_tmp
            grad_coef_kh = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          end if

          ! gradient
          !
          grad(1) = grad_coef_kh * dij(1)
          grad(2) = grad_coef_kh * dij(2)
          grad(3) = grad_coef_kh * dij(3)

          ! store force
          !
          force_tmp(1) = force_tmp(1) - grad(1)
          force_tmp(2) = force_tmp(2) - grad(2)
          force_tmp(3) = force_tmp(3) - grad(3)
          force(1,j,id) = force(1,j,id) + grad(1)
          force(2,j,id) = force(2,j,id) + grad(2)
          force(3,j,id) = force(3,j,id) + grad(3)

          ! virial
          !
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do                  ! k

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        ekh_omp(id) = ekh_omp(id) + ekh_omp_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)

      end if

    end do
    !$omp end parallel

    do i = 1, nthread
      ekh = ekh + ekh_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGIDRKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_IDR_KH_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_KH
  !> @brief        calculate interaction in the KH model (NOBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ekh      : energy of KH interactions
  !! @note         KH model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_KH(enefunc, pairlist, coord, force, virial, ekh)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ekh

    ! local variables
    integer                   :: my_id, id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: natom
    integer                   :: num_kh, ini_kh, fin_kh
    integer                   :: i, j, k, l, m
    integer                   :: model_set
    real(wp)                  :: alpha, eps0
    real(wp)                  :: lambda, epsilon
    real(wp)                  :: sigma, sigma_sqr, sigma_6th
    real(wp)                  :: dij(3), rij_sqr
    real(wp)                  :: cutoff, cutoff_sqr
    real(wp)                  :: inv_rij_sqr, inv_rij_6th
    real(wp)                  :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)                  :: grad_coef_kh, grad(3)
    real(wp)                  :: ekh_tmp
    real(wp)                  :: cg_KH_mod_A_lambda, cg_KH_mod_A_eps_0
    real(wp)                  :: cg_KH_mod_B_lambda, cg_KH_mod_B_eps_0
    real(wp)                  :: cg_KH_mod_C_lambda, cg_KH_mod_C_eps_0
    real(wp)                  :: cg_KH_mod_D_lambda, cg_KH_mod_D_eps_0
    real(wp)                  :: cg_KH_mod_E_lambda, cg_KH_mod_E_eps_0
    real(wp)                  :: cg_KH_mod_F_lambda, cg_KH_mod_F_eps_0

    integer,  pointer         :: kh_list(:,:)
    integer,  pointer         :: kh_model_list(:,:)
    integer,  pointer         :: num_kh_calc(:,:)

    integer,  pointer         :: atom_type(:)
    real(wp), pointer         :: kh_sigma_half(:)
    real(wp), pointer         :: kh_epsilon(:, :)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGKH, TimerOn)

    natom           = size(coord(1,:))

    kh_list         => pairlist%cg_KH_list
    kh_model_list   => pairlist%cg_KH_model_list
    num_kh_calc     => pairlist%num_cg_KH_calc

    atom_type       => enefunc%atom_cls
    kh_sigma_half   => enefunc%cg_KH_sigma_half
    kh_epsilon      => enefunc%cg_KH_epsilon

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

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

    num_kh          = 0

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                         &
    !$omp firstprivate(num_kh)                           &
    !$omp private(my_id, id,                             &
    !$omp         ini_kh, fin_kh,                        &
    !$omp         i, k, j, l, m,                         &
    !$omp         model_set, lambda,  epsilon,           &
    !$omp         alpha, eps0,                           &
    !$omp         sigma, sigma_sqr, sigma_6th,           &
    !$omp         dij, rij_sqr,                          &
    !$omp         inv_rij_sqr, inv_rij_6th,              &
    !$omp         sig_over_rij_6th, sig_over_rij_12th,   &
    !$omp         grad_coef_kh, grad, ekh_tmp            &
    !$omp         )                                      &
    !$omp shared(coord, my_city_rank, nproc_city,        &
    !$omp        nthread, natom,                         &
    !$omp        cutoff, cutoff_sqr,                     &
    !$omp        kh_list, kh_model_list,                 &
    !$omp        num_kh_calc, atom_type,                 &
    !$omp        cg_KH_mod_A_lambda, cg_KH_mod_A_eps_0,  &
    !$omp        cg_KH_mod_B_lambda, cg_KH_mod_B_eps_0,  &
    !$omp        cg_KH_mod_C_lambda, cg_KH_mod_C_eps_0,  &
    !$omp        cg_KH_mod_D_lambda, cg_KH_mod_D_eps_0,  &
    !$omp        cg_KH_mod_E_lambda, cg_KH_mod_E_eps_0,  &
    !$omp        cg_KH_mod_F_lambda, cg_KH_mod_F_eps_0,  &
    !$omp        kh_sigma_half, kh_epsilon,              &
    !$omp        force)                                  &
    !$omp reduction(+:virial) reduction(+:ekh)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i = 1, natom-1

      ini_kh = num_kh + 1
      fin_kh = num_kh + num_kh_calc(i,id)
      num_kh = fin_kh

      do k = ini_kh, fin_kh

        j = kh_list(k,id)

        select case (kh_model_list(k, id))
        case (1)
          alpha = cg_KH_mod_A_lambda
          eps0  = cg_KH_mod_A_eps_0
        case (2)
          alpha = cg_KH_mod_B_lambda
          eps0  = cg_KH_mod_B_eps_0
        case (3)
          alpha = cg_KH_mod_C_lambda
          eps0  = cg_KH_mod_C_eps_0
        case (4)
          alpha = cg_KH_mod_D_lambda
          eps0  = cg_KH_mod_D_eps_0
        case (5)
          alpha = cg_KH_mod_E_lambda
          eps0  = cg_KH_mod_E_eps_0
        case (6)
          alpha = cg_KH_mod_F_lambda
          eps0  = cg_KH_mod_F_eps_0
        end select

        sigma     = kh_sigma_half(i) + kh_sigma_half(j)
        sigma_sqr = sigma * sigma
        sigma_6th = sigma_sqr * sigma_sqr * sigma_sqr

        epsilon   = alpha * (kh_epsilon(atom_type(i), atom_type(j)) - eps0) * 0.593_wp
        if (epsilon > 0) then
          lambda  = -1.0_wp
        else
          lambda  =  1.0_wp
          epsilon = -epsilon
        end if

        dij(1:3)  = coord(1:3,j) - coord(1:3,i)
        rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij_sqr > cutoff_sqr) cycle

        inv_rij_sqr       = 1.0_wp / rij_sqr
        inv_rij_6th       = inv_rij_sqr * inv_rij_sqr * inv_rij_sqr

        sig_over_rij_6th  = sigma_6th * inv_rij_6th
        sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

        ekh_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

        if (sig_over_rij_6th >= 0.5_wp) then
          ekh = ekh + ekh_tmp + (1.0_wp - lambda) * epsilon
          grad_coef_kh = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        else
          ekh = ekh + lambda * ekh_tmp
          grad_coef_kh = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
        end if

        ! gradient
        !
        grad(1:3)     = grad_coef_kh * dij(1:3)

        ! store force
        !
        force(1:3,i,id)  = force(1:3,i,id) - grad(1:3)
        force(1:3,j,id)  = force(1:3,j,id) + grad(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3, l) = virial(1:3, l) - dij(1:3) * grad(l)
        end do

      end do
    end do
    !$omp end parallel

    call timer(TimerCGKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_KH

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_KH_pbc
  !> @brief        calculate interaction in the KH model (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ekh     : energy of KH interactions
  !! @note         KH model
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_KH_pbc(enefunc, boundary, pairlist, &
                                      coord, force, virial, ekh)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ekh

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: n_idr, i_idr
    integer           :: num_kh, ini_kh, fin_kh
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc ! variables for PBC
    integer           :: atomtype_i, atomtype_j
    integer           :: model_set
    real(wp)          :: coor_i_tmp(3)
    real(wp)          :: alpha, eps0
    real(wp)          :: lambda, epsilon
    real(wp)          :: si_tmp ! sigma of i
    real(wp)          :: sigma, sigma_sqr, sigma_6th
    real(wp)          :: dij(3), rij_sqr
    real(wp)          :: cutoff, cutoff_sqr
    real(wp)          :: inv_rij_sqr, inv_rij_6th
    real(wp)          :: sig_over_rij_sqr
    real(wp)          :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)          :: grad_coef_kh, grad(3)
    real(wp)          :: ekh_tmp
    !
    real(wp)          :: cg_KH_mod_A_lambda, cg_KH_mod_A_eps_0
    real(wp)          :: cg_KH_mod_B_lambda, cg_KH_mod_B_eps_0
    real(wp)          :: cg_KH_mod_C_lambda, cg_KH_mod_C_eps_0
    real(wp)          :: cg_KH_mod_D_lambda, cg_KH_mod_D_eps_0
    real(wp)          :: cg_KH_mod_E_lambda, cg_KH_mod_E_eps_0
    real(wp)          :: cg_KH_mod_F_lambda, cg_KH_mod_F_eps_0
    !
    real(wp)          :: bsize(3), inv_bsize(3)
    !
    real(wp)          :: force_tmp(3), factor
    real(wp)          :: virial_tmp(3,3)
    real(wp)          :: virial_omp(3,3,nthread)
    real(wp)          :: ekh_omp(nthread), ekh_omp_tmp
    !
    integer,  pointer :: cg_list_idr(:)
    integer,  pointer :: kh_list(:,:)
    integer,  pointer :: kh_model_list(:,:)
    integer,  pointer :: num_kh_calc(:,:)
    !
    integer,  pointer :: atom_type(:)
    real(wp), pointer :: kh_sigma_half(:)
    real(wp), pointer :: kh_epsilon(:, :)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGKH, TimerOn)

    n_idr          =  enefunc%num_cg_particle_KH
    cg_list_idr    => enefunc%cg_particle_KH

    kh_list        => pairlist%cg_KH_list
    kh_model_list  => pairlist%cg_KH_model_list
    num_kh_calc    => pairlist%num_cg_KH_calc

    atom_type      => enefunc%atom_cls
    kh_sigma_half  => enefunc%cg_KH_sigma_half
    kh_epsilon     => enefunc%cg_KH_epsilon

    cutoff           = enefunc%cg_cutoffdist_126
    cutoff_sqr       = cutoff * cutoff

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

    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    inv_bsize(1:3) = 1.0_wp/bsize(1:3)

    num_kh         = 0

    ekh_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                        &
    !$omp firstprivate(num_kh)                          &
    !$omp private(my_id, id,                            &
    !$omp         ini_kh, fin_kh,                       &
    !$omp         i, k, j, l, m, i_idr,                 &
    !$omp         i1, i2, i3, k_pbc,                    &
    !$omp         model_set, lambda,  epsilon,          &
    !$omp         alpha, eps0, si_tmp,                  &
    !$omp         atomtype_i, atomtype_j,               &
    !$omp         sigma, sigma_sqr, sigma_6th,          &
    !$omp         dij, rij_sqr, coor_i_tmp,             &
    !$omp         inv_rij_sqr, sig_over_rij_sqr,        &
    !$omp         sig_over_rij_6th, sig_over_rij_12th,  &
    !$omp         force_tmp, factor,                    &
    !$omp         virial_tmp, ekh_omp_tmp,              &
    !$omp         grad_coef_kh, grad, ekh_tmp           &
    !$omp         )                                     &
    !$omp shared(coord, my_city_rank, nproc_city,       &
    !$omp        nthread, n_idr,                        &
    !$omp        cutoff, cutoff_sqr,                    &
    !$omp        kh_list, kh_model_list,                &
    !$omp        num_kh_calc, atom_type,                &
    !$omp        kh_sigma_half, kh_epsilon,             &
    !$omp        cg_KH_mod_A_lambda, cg_KH_mod_A_eps_0, &
    !$omp        cg_KH_mod_B_lambda, cg_KH_mod_B_eps_0, &
    !$omp        cg_KH_mod_C_lambda, cg_KH_mod_C_eps_0, &
    !$omp        cg_KH_mod_D_lambda, cg_KH_mod_D_eps_0, &
    !$omp        cg_KH_mod_E_lambda, cg_KH_mod_E_eps_0, &
    !$omp        cg_KH_mod_F_lambda, cg_KH_mod_F_eps_0, &
    !$omp        cg_list_idr, bsize, inv_bsize,         &
    !$omp        virial_omp, ekh_omp,                   &
    !$omp        force)                                 &
    !$omp reduction(+:virial) reduction(+:ekh)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id = id + 1


    do i_idr = 1, n_idr - 1

      ini_kh = num_kh + 1
      fin_kh = num_kh + num_kh_calc(i_idr, id)
      num_kh = fin_kh

      if (fin_kh >= ini_kh) then

        i = cg_list_idr(i_idr)

        coor_i_tmp(1:3) = coord(1:3, i)
        si_tmp = kh_sigma_half(i)
        atomtype_i = atom_type(i)
        force_tmp(1:3)      = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp
        ekh_omp_tmp         = 0.0_wp

        do k = ini_kh, fin_kh

          k_pbc = kh_list(k,id)
          j = k_pbc / 27
          k_pbc = k_pbc - j*27
          i3 = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2 = k_pbc / 3
          i1 = k_pbc - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          select case (kh_model_list(k, id))
          case (1)
            alpha = cg_KH_mod_A_lambda
            eps0  = cg_KH_mod_A_eps_0
          case (2)
            alpha = cg_KH_mod_B_lambda
            eps0  = cg_KH_mod_B_eps_0
          case (3)
            alpha = cg_KH_mod_C_lambda
            eps0  = cg_KH_mod_C_eps_0
          case (4)
            alpha = cg_KH_mod_D_lambda
            eps0  = cg_KH_mod_D_eps_0
          case (5)
            alpha = cg_KH_mod_E_lambda
            eps0  = cg_KH_mod_E_eps_0
          case (6)
            alpha = cg_KH_mod_F_lambda
            eps0  = cg_KH_mod_F_eps_0
          end select

          sigma     = si_tmp + kh_sigma_half(j)
          sigma_sqr = sigma * sigma

          atomtype_j = atom_type(j)
          epsilon   = alpha * (kh_epsilon(atomtype_i, atomtype_j) - eps0) * 0.593_wp
          if (epsilon > 0) then
            lambda  = -1.0_wp
          else
            lambda  =  1.0_wp
            epsilon = -epsilon
          end if

          dij(1)  = coord(1,j) - coor_i_tmp(1) + bsize(1) * real(i1, wp)
          dij(2)  = coord(2,j) - coor_i_tmp(2) + bsize(2) * real(i2, wp)
          dij(3)  = coord(3,j) - coor_i_tmp(3) + bsize(3) * real(i3, wp)
          rij_sqr = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij_sqr > cutoff_sqr) cycle

          inv_rij_sqr       = 1.0_wp / rij_sqr
          sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr

          sig_over_rij_6th  = sig_over_rij_sqr * sig_over_rij_sqr * sig_over_rij_sqr
          sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

          ekh_tmp          = 4.0_wp * epsilon * (sig_over_rij_12th - sig_over_rij_6th)

          if (sig_over_rij_6th >= 0.5_wp) then
            ekh_omp_tmp = ekh_omp_tmp + ekh_tmp + (1.0_wp - lambda) * epsilon
            grad_coef_kh = - inv_rij_sqr * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          else
            ekh_omp_tmp = ekh_omp_tmp + lambda * ekh_tmp
            grad_coef_kh = - inv_rij_sqr * lambda * 4.0_wp * epsilon * (6.0_wp * sig_over_rij_6th - 12.0_wp * sig_over_rij_12th)
          end if

          ! gradient
          !
          grad(1) = grad_coef_kh * dij(1)
          grad(2) = grad_coef_kh * dij(2)
          grad(3) = grad_coef_kh * dij(3)

          ! store force
          !
          force_tmp(1) = force_tmp(1) - grad(1)
          force_tmp(2) = force_tmp(2) - grad(2)
          force_tmp(3) = force_tmp(3) - grad(3)
          force(1,j,id) = force(1,j,id) + grad(1)
          force(2,j,id) = force(2,j,id) + grad(2)
          force(3,j,id) = force(3,j,id) + grad(3)

          ! virial
          !
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do                  ! k

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        ekh_omp(id) = ekh_omp(id) + ekh_omp_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)
      end if

    end do                      ! i_idr
    !$omp end parallel

    do i = 1, nthread
      ekh = ekh + ekh_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGKH, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_KH_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_repulsive
  !> @brief        calculate repulsive gradients in CG DNA basepairing !
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
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: dist
    real(wp),                intent(in)    :: sigma
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: epsilon
    real(wp),                intent(inout) :: ene
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    real(wp)                 :: expon, exptmp, etmp
    real(wp)                 :: f_coef


    if (dist < sigma) then
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

    return

  end subroutine calculate_repulsive

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_pwmcos
  !> @brief        calculate PWMcos interactions between protein and DNA
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epwmcos  : PWMcos energy
  !! @note         3SPN.2C + protein
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_pwmcos(enefunc, pairlist, coord, force, &
                                      virial, epwmcos)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epwmcos

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: natom, npwmcos
    integer           :: ini_pwmcos, fin_pwmcos
    integer           :: k, l, m
    integer           :: ipwm, ipair
    ! indicex + params
    integer           :: i_CA_0, i_CA_N, i_CA_C
    integer           :: i_DB_0, i_DB_5, i_DB_3
    integer           :: i_DS_0
    integer           :: j_base_type
    real(wp)          :: pwm_shift, pwm_factor
    real(wp)          :: cutoff, cutoff_sqr
    ! native values
    real(wp)          :: r00, t10, t20, t30
    ! vectors and distances
    real(wp)          :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)          :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)          :: v53(3), r53, r53_sqr, r53_inv, e53(3)
    real(wp)          :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)          :: t1, t2, t3
    real(wp)          :: cos_t1, cos_t2, cos_t3
    real(wp)          :: sin_t1, sin_t2, sin_t3
    ! modulating functions
    real(wp)          :: dt1, dt2, dt3
    real(wp)          :: abs_dt1, abs_dt2, abs_dt3
    real(wp)          :: cos_dt1, cos_dt2, cos_dt3
    real(wp)          :: sin_dt1, sin_dt2, sin_dt3
    real(wp)          :: sigma, phi, phi2
    real(wp)          :: sig_sqr
    real(wp)          :: ktheta_2, ktheta
    ! energy calc
    real(wp)          :: pwm_score
    real(wp)          :: ene_coef_r0
    real(wp)          :: grad_coef_t1, ene_coef_t1
    real(wp)          :: grad_coef_t2, ene_coef_t2
    real(wp)          :: grad_coef_t3, ene_coef_t3
    real(wp)          :: etmp0
    real(wp)          :: e_tmp_pwmcos
    ! force coef
    real(wp)          :: var_tmp
    real(wp)          :: grad_r0(3)
    real(wp)          :: grad_t1(6)
    real(wp)          :: grad_t2(6)
    real(wp)          :: grad_t3(6)
    ! force calc
    real(wp)          :: f_r0_coef_tmp
    real(wp)          :: f_r0_tmp(6)
    real(wp)          :: f_t1_coef_tmp
    real(wp)          :: f_t1_tmp(6)
    real(wp)          :: f_t2_coef_tmp
    real(wp)          :: f_t2_tmp(6)
    real(wp)          :: f_t3_coef_tmp
    real(wp)          :: f_t3_tmp(6)
    real(wp)          :: f_CA_0(3)
    real(wp)          :: f_CA_N(3)
    real(wp)          :: f_CA_C(3)
    real(wp)          :: f_DB_0(3)
    real(wp)          :: f_DB_5(3)
    real(wp)          :: f_DB_3(3)
    real(wp)          :: f_DS_0(3)

    integer,  pointer :: pwmcos_list(:,:)
    integer,  pointer :: n_pwmcos_pair(:)
    integer,  pointer :: pwmcos_id(:)
    integer,  pointer :: pwmcos_id_N(:)
    integer,  pointer :: pwmcos_id_C(:)
    real(wp), pointer :: pwmcos_r0(:)
    real(wp), pointer :: pwmcos_theta1(:)
    real(wp), pointer :: pwmcos_theta2(:)
    real(wp), pointer :: pwmcos_theta3(:)
    real(wp), pointer :: pwmcos_ene_A(:)
    real(wp), pointer :: pwmcos_ene_C(:)
    real(wp), pointer :: pwmcos_ene_G(:)
    real(wp), pointer :: pwmcos_ene_T(:)
    real(wp), pointer :: pwmcos_gamma(:)
    real(wp), pointer :: pwmcos_eps(:)
    integer,  pointer :: pwmcos_spec(:)
    integer,  pointer :: pwmcos_to_pair(:)
    integer,  pointer :: base_type(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPWMcos, TimerOn)

    natom   = size(coord(1,:))
    npwmcos = enefunc%num_pwmcos_terms

    pwmcos_list    => pairlist%cg_pwmcos_list
    n_pwmcos_pair  => pairlist%num_cg_pwmcos

    pwmcos_id      => enefunc%pwmcos_protein_id
    pwmcos_id_N    => enefunc%pwmcos_protein_id_N
    pwmcos_id_C    => enefunc%pwmcos_protein_id_C
    pwmcos_r0      => enefunc%pwmcos_r0
    pwmcos_theta1  => enefunc%pwmcos_theta1
    pwmcos_theta2  => enefunc%pwmcos_theta2
    pwmcos_theta3  => enefunc%pwmcos_theta3
    pwmcos_ene_A   => enefunc%pwmcos_ene_A
    pwmcos_ene_C   => enefunc%pwmcos_ene_C
    pwmcos_ene_G   => enefunc%pwmcos_ene_G
    pwmcos_ene_T   => enefunc%pwmcos_ene_T
    pwmcos_gamma   => enefunc%pwmcos_gamma
    pwmcos_eps     => enefunc%pwmcos_eps
    pwmcos_spec    => enefunc%pwmcos_specificity
    pwmcos_to_pair => enefunc%pwmcos_to_pairlist_id
    base_type      => enefunc%NA_base_type

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

    !$omp parallel default(none)                  &
    !$omp private(my_id, id,                      &
    !$omp         ini_pwmcos, fin_pwmcos,         &
    !$omp         k, l, m,                        &
    !$omp         ipwm, ipair,                    &
    !$omp         i_CA_0, i_CA_N, i_CA_C,         &
    !$omp         i_DB_0, i_DB_5, i_DB_3,         &
    !$omp         i_DS_0,                         &
    !$omp         j_base_type,                    &
    !$omp         pwm_shift, pwm_factor,          &
    !$omp         cutoff, cutoff_sqr,             &
    !$omp         r00, t10, t20, t30,             &
    !$omp         vbc, vbs, v53, vcn,             &
    !$omp         rbc, rbs, r53, rcn,             &
    !$omp         rbc_sqr, rbs_sqr,               &
    !$omp         r53_sqr, rcn_sqr,               &
    !$omp         rbc_inv, rbs_inv,               &
    !$omp         r53_inv, rcn_inv,               &
    !$omp         ebc, ebs, e53, ecn,             &
    !$omp         t1, t2, t3,                     &
    !$omp         cos_t1, cos_t2, cos_t3,         &
    !$omp         sin_t1, sin_t2, sin_t3,         &
    !$omp         dt1, dt2, dt3,                  &
    !$omp         abs_dt1, abs_dt2, abs_dt3,      &
    !$omp         cos_dt1, cos_dt2, cos_dt3,      &
    !$omp         sin_dt1, sin_dt2, sin_dt3,      &
    !$omp         pwm_score,                      &
    !$omp         ene_coef_r0,                    &
    !$omp         grad_coef_t1, ene_coef_t1,      &
    !$omp         grad_coef_t2, ene_coef_t2,      &
    !$omp         grad_coef_t3, ene_coef_t3,      &
    !$omp         etmp0, e_tmp_pwmcos,            &
    !$omp         var_tmp, grad_r0,               &
    !$omp         grad_t1, grad_t2, grad_t3,      &
    !$omp         f_r0_coef_tmp, f_r0_tmp,        &
    !$omp         f_t1_coef_tmp, f_t1_tmp,        &
    !$omp         f_t2_coef_tmp, f_t2_tmp,        &
    !$omp         f_t3_coef_tmp, f_t3_tmp,        &
    !$omp         f_CA_0, f_CA_N, f_CA_C,         &
    !$omp         f_DB_0, f_DB_5, f_DB_3,         &
    !$omp         f_DS_0                          &
    !$omp         )                               &
    !$omp shared(coord, my_city_rank, nproc_city, &
    !$omp        nthread, natom, npwmcos,         &
    !$omp        sig_sqr,                         &
    !$omp        sigma, phi, phi2,                &
    !$omp        ktheta_2, ktheta,                &
    !$omp        pwmcos_list, n_pwmcos_pair,      &
    !$omp        pwmcos_id, pwmcos_id_N,          &
    !$omp        pwmcos_id_C, pwmcos_to_pair,     &
    !$omp        pwmcos_r0, pwmcos_theta1,        &
    !$omp        pwmcos_theta2, pwmcos_theta3,    &
    !$omp        pwmcos_ene_A, pwmcos_ene_C,      &
    !$omp        pwmcos_ene_G, pwmcos_ene_T,      &
    !$omp        pwmcos_gamma, pwmcos_eps,        &
    !$omp        pwmcos_spec, base_type,          &
    !$omp        force)                           &
    !$omp reduction(+:virial) reduction(+:epwmcos)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do ipwm = 1, npwmcos

      ! if (mod(ipwm - 1, nproc_city*nthread) .ne. my_id) cycle

      i_CA_0     = pwmcos_id      (ipwm)
      ipair      = pwmcos_to_pair (ipwm)

      if (mod(ipair - 1, nproc_city*nthread) .ne. my_id) cycle

      i_CA_N     = pwmcos_id_N    (ipwm)
      i_CA_C     = pwmcos_id_C    (ipwm)

      r00        = pwmcos_r0      (ipwm)
      t10        = pwmcos_theta1  (ipwm)
      t20        = pwmcos_theta2  (ipwm)
      t30        = pwmcos_theta3  (ipwm)

      pwm_shift  = pwmcos_eps     (ipwm)
      pwm_factor = pwmcos_gamma   (ipwm)

      cutoff     = pwmcos_r0(ipwm) + 5.0
      cutoff_sqr = cutoff * cutoff

      ! -------------
      ! Ca_C ==> Ca_N
      ! -------------
      !
      vcn(1:3)   = coord(1:3, i_CA_N) - coord(1:3, i_CA_C)
      rcn_sqr    = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
      rcn        = sqrt(rcn_sqr)
      rcn_inv    = 1.0_wp / rcn
      ecn(1:3)   = vcn(1:3) * rcn_inv

      ini_pwmcos = 1
      fin_pwmcos = n_pwmcos_pair(ipair)

      do k = ini_pwmcos, fin_pwmcos

        i_DB_0   = pwmcos_list(k, ipair)

        vbc(1:3) = coord(1:3, i_CA_0) - coord(1:3, i_DB_0)
        rbc_sqr  = vbc(1) * vbc(1) + vbc(2) * vbc(2) + vbc(3) * vbc(3)

        if (rbc_sqr >= cutoff_sqr) cycle

        f_CA_0(1:3) = 0.0_wp
        f_CA_N(1:3) = 0.0_wp
        f_CA_C(1:3) = 0.0_wp
        f_DB_0(1:3) = 0.0_wp
        f_DB_5(1:3) = 0.0_wp
        f_DB_3(1:3) = 0.0_wp
        f_DS_0(1:3) = 0.0_wp

        etmp0        = 0.0_wp
        e_tmp_pwmcos = 0.0_wp

        ! ==================================
        ! Distance and angle calculations...
        ! ==================================

        j_base_type = base_type(i_DB_0)

        i_DB_5 = i_DB_0 - 3
        i_DB_3 = i_DB_0 + 3
        i_DS_0 = i_DB_0 - 1

        ! --------------
        ! sugar <== base
        ! --------------
        !
        vbs(1:3) = coord(1:3, i_DS_0) - coord(1:3, i_DB_0)
        rbs_sqr  = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
        rbs      = sqrt(rbs_sqr)
        rbs_inv  = 1.0_wp / rbs
        ebs(1:3) = vbs(1:3) * rbs_inv

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
        v53(1:3) = coord(1:3, i_DB_3) - coord(1:3, i_DB_5)
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
        t1     = acos(cos_t1)

        ! ------------------------------------------
        ! Angle t2: base5' - base3' <==> base0 -- Ca
        ! ------------------------------------------
        !
        cos_t2 = ebc(1) * e53(1) + ebc(2) * e53(2) + ebc(3) * e53(3)
        if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
        if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
        sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
        t2     = acos(cos_t2)

        ! ---------------------------------------
        ! Angle t3: base0 -- Ca <==> Ca_N -- Ca_C
        ! ---------------------------------------
        !
        cos_t3 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
        if (cos_t3 >  1.0_wp) cos_t3 =  1.0_wp
        if (cos_t3 < -1.0_wp) cos_t3 = -1.0_wp
        sin_t3 = sqrt(1.0_wp - cos_t3 * cos_t3)
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


        ! ============================================================================
        ! Energy/Force Calculation: gamma * (pwm_score + pwm_shift) * f * g1 * g2 * g3
        ! ============================================================================
        !
        ! -----------------
        ! Simple angle test
        ! -----------------
        !
        if (abs_dt1 >= phi2 .or. &
            abs_dt2 >= phi2 .or. &
            abs_dt3 >= phi2) then
          cycle
        end if


        ! ---------------------------------
        ! PWM energy: sequence specificity!
        ! ---------------------------------
        !
        if (j_base_type == NABaseTypeDBA) then
          pwm_score = pwmcos_ene_A(ipwm)
        else if (j_base_type == NABaseTypeDBC) then
          pwm_score = pwmcos_ene_C(ipwm)
        else if (j_base_type == NABaseTypeDBG) then
          pwm_score = pwmcos_ene_G(ipwm)
        else if (j_base_type == NABaseTypeDBT) then
          pwm_score = pwmcos_ene_T(ipwm)
        end if
        etmp0 = pwm_factor * (pwm_score + pwm_shift)

        ! -------------------
        ! Distance modulating
        ! -------------------
        !
        call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, grad_r0)

        ! ----------------
        ! Angle modulating
        ! ----------------
        if (abs_dt1 < phi) then
          ene_coef_t1  = 1.0_wp
          grad_coef_t1 = 0.0_wp
          grad_t1(1:6) = 0.0_wp
        else
          sin_dt1      = sin(ktheta * dt1)
          ene_coef_t1  = sin_dt1 * sin_dt1
          call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
          grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
        end if
        if (abs_dt2 < phi) then
          ene_coef_t2  = 1.0_wp
          grad_coef_t2 = 0.0_wp
          grad_t2(1:6) = 0.0_wp
        else
          sin_dt2      = sin(ktheta * dt2)
          ene_coef_t2  = sin_dt2 * sin_dt2
          call calculate_angle(vbc, v53, rbc, r53, cos_t2, sin_t2, grad_t2)
          grad_coef_t2 = ktheta * sin(2.0_wp * ktheta * dt2)
        end if
        if (abs_dt3 < phi) then
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
        e_tmp_pwmcos = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2 * ene_coef_t3
        epwmcos      = epwmcos + e_tmp_pwmcos

        ! ======
        ! Forces
        ! ======
        !
        ! --
        ! r0
        ! --
        !
        ! f_r0_coef_tmp = e_tmp_pwmcos / ene_coef_r0
        f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2 * ene_coef_t3
        f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
        f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
        f_DB_0(1:3)   = f_DB_0(1:3) - f_r0_tmp(1:3)
        !
        ! --
        ! t1
        ! --
        !
        ! f_t1_coef_tmp = e_tmp_pwmcos / ene_coef_t1 * grad_coef_t1
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
        ! f_t2_coef_tmp = e_tmp_pwmcos / ene_coef_t2 * grad_coef_t2
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
        ! f_t3_coef_tmp = e_tmp_pwmcos / ene_coef_t3 * grad_coef_t3
        f_t3_coef_tmp =  etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2 * grad_coef_t3
        f_t3_tmp(1:6) = - f_t3_coef_tmp * grad_t3(1:6)
        f_CA_0(1:3)   = f_CA_0(1:3) + f_t3_tmp(1:3)
        f_DB_0(1:3)   = f_DB_0(1:3) - f_t3_tmp(1:3)
        f_CA_C(1:3)   = f_CA_C(1:3) - f_t3_tmp(4:6)
        f_CA_N(1:3)   = f_CA_N(1:3) + f_t3_tmp(4:6)


        ! ------------------------
        ! update all the forces...
        ! ------------------------
        !
        force(1:3, i_CA_0,id) = force(1:3, i_CA_0,id) + f_CA_0(1:3)
        force(1:3, i_CA_N,id) = force(1:3, i_CA_N,id) + f_CA_N(1:3)
        force(1:3, i_CA_C,id) = force(1:3, i_CA_C,id) + f_CA_C(1:3)
        force(1:3, i_DB_0,id) = force(1:3, i_DB_0,id) + f_DB_0(1:3)
        force(1:3, i_DB_5,id) = force(1:3, i_DB_5,id) + f_DB_5(1:3)
        force(1:3, i_DB_3,id) = force(1:3, i_DB_3,id) + f_DB_3(1:3)
        force(1:3, i_DS_0,id) = force(1:3, i_DS_0,id) + f_DS_0(1:3)

        ! ------
        ! Virial
        ! ------
        do l = 1, 3
          do m = l+1, 3
            var_tmp &
                = f_CA_0(l) * coord(m, i_CA_0) &
                + f_CA_N(l) * coord(m, i_CA_N) &
                + f_CA_C(l) * coord(m, i_CA_C) &
                + f_DB_0(l) * coord(m, i_DB_0) &
                + f_DB_5(l) * coord(m, i_DB_5) &
                + f_DB_3(l) * coord(m, i_DB_3) &
                + f_DS_0(l) * coord(m, i_DS_0)
            virial(m, l) = virial(m, l) - var_tmp
            virial(l, m) = virial(l, m) - var_tmp
          end do
          var_tmp &
              = f_CA_0(l) * coord(l, i_CA_0) &
              + f_CA_N(l) * coord(l, i_CA_N) &
              + f_CA_C(l) * coord(l, i_CA_C) &
              + f_DB_0(l) * coord(l, i_DB_0) &
              + f_DB_5(l) * coord(l, i_DB_5) &
              + f_DB_3(l) * coord(l, i_DB_3) &
              + f_DS_0(l) * coord(l, i_DS_0)
          virial(l,l) = virial(l,l) - var_tmp
        end do

      end do

    end do
    !$omp end parallel

    call timer(TimerCGPWMcos, TimerOff)
    call timer(TimerNonBond, TimerOff)


    return

  end subroutine compute_energy_CG_pwmcos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_pwmcos_pbc
  !> @brief        calculate PWMcos interactions between protein and DNA (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : PBC coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epwmcos  : PWMcos energy
  !! @note         3SPN.2C + protein
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_pwmcos_pbc(enefunc, boundary, pairlist, &
                                          coord, force, virial, epwmcos)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epwmcos

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: npwmcos
    integer           :: num_pwmcos, ini_pwmcos, fin_pwmcos
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc ! for PBC
    integer           :: i_atom, j_atom
    integer           :: ipwm
    ! indicex + params
    integer           :: i_CA_0, i_CA_N, i_CA_C
    integer           :: i_DB_0, i_DB_5, i_DB_3
    integer           :: i_DS_0
    integer           :: j_base_type
    real(wp)          :: pwm_shift, pwm_factor
    real(wp)          :: cutoff, cutoff_sqr
    ! native values
    real(wp)          :: r00, t10, t20, t30
    ! vectors and distances
    real(wp)          :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)          :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)          :: v53(3), r53, r53_sqr, r53_inv, e53(3)
    real(wp)          :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)          :: t1, t2, t3
    real(wp)          :: cos_t1, cos_t2, cos_t3
    real(wp)          :: sin_t1, sin_t2, sin_t3
    ! modulating functions
    real(wp)          :: dt1, dt2, dt3
    real(wp)          :: abs_dt1, abs_dt2, abs_dt3
    real(wp)          :: cos_dt1, cos_dt2, cos_dt3
    real(wp)          :: sin_dt1, sin_dt2, sin_dt3
    real(wp)          :: sigma, phi, phi2
    real(wp)          :: sig_sqr
    real(wp)          :: ktheta_2, ktheta
    ! energy calc
    real(wp)          :: pwm_score
    real(wp)          :: ene_coef_r0
    real(wp)          :: grad_coef_t1, ene_coef_t1
    real(wp)          :: grad_coef_t2, ene_coef_t2
    real(wp)          :: grad_coef_t3, ene_coef_t3
    real(wp)          :: etmp0
    real(wp)          :: e_tmp_pwmcos
    ! force coef
    real(wp)          :: var_tmp
    real(wp)          :: grad_r0(3)
    real(wp)          :: grad_t1(6)
    real(wp)          :: grad_t2(6)
    real(wp)          :: grad_t3(6)
    ! force calc
    real(wp)          :: f_r0_coef_tmp
    real(wp)          :: f_r0_tmp(6)
    real(wp)          :: f_t1_coef_tmp
    real(wp)          :: f_t1_tmp(6)
    real(wp)          :: f_t2_coef_tmp
    real(wp)          :: f_t2_tmp(6)
    real(wp)          :: f_t3_coef_tmp
    real(wp)          :: f_t3_tmp(6)
    real(wp)          :: f_CA_0(3)
    real(wp)          :: f_CA_N(3)
    real(wp)          :: f_CA_C(3)
    real(wp)          :: f_DB_0(3)
    real(wp)          :: f_DB_5(3)
    real(wp)          :: f_DB_3(3)
    real(wp)          :: f_DS_0(3)

    integer,  pointer :: cg_list_base(:)
    integer,  pointer :: pwmcos_list(:,:)
    integer,  pointer :: n_pwmcos_pair(:, :)
    integer,  pointer :: pwmcos_id(:)
    integer,  pointer :: pwmcos_id_N(:)
    integer,  pointer :: pwmcos_id_C(:)
    real(wp), pointer :: pwmcos_r0(:)
    real(wp), pointer :: pwmcos_theta1(:)
    real(wp), pointer :: pwmcos_theta2(:)
    real(wp), pointer :: pwmcos_theta3(:)
    real(wp), pointer :: pwmcos_ene_A(:)
    real(wp), pointer :: pwmcos_ene_C(:)
    real(wp), pointer :: pwmcos_ene_G(:)
    real(wp), pointer :: pwmcos_ene_T(:)
    real(wp), pointer :: pwmcos_gamma(:)
    real(wp), pointer :: pwmcos_eps(:)
    integer,  pointer :: pwmcos_spec(:)
    integer,  pointer :: base_type(:)

    real(wp)          :: bsize(3), half_bsize(3)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPWMcos, TimerOn)

    npwmcos        = enefunc%num_pwmcos_terms

    pwmcos_list    => pairlist%cg_pwmcos_list
    n_pwmcos_pair  => pairlist%num_cg_pwmcos_calc
    cg_list_base   => enefunc%cg_particle_DNA_base

    pwmcos_id      => enefunc%pwmcos_protein_id
    pwmcos_id_N    => enefunc%pwmcos_protein_id_N
    pwmcos_id_C    => enefunc%pwmcos_protein_id_C
    pwmcos_r0      => enefunc%pwmcos_r0
    pwmcos_theta1  => enefunc%pwmcos_theta1
    pwmcos_theta2  => enefunc%pwmcos_theta2
    pwmcos_theta3  => enefunc%pwmcos_theta3
    pwmcos_ene_A   => enefunc%pwmcos_ene_A
    pwmcos_ene_C   => enefunc%pwmcos_ene_C
    pwmcos_ene_G   => enefunc%pwmcos_ene_G
    pwmcos_ene_T   => enefunc%pwmcos_ene_T
    pwmcos_gamma   => enefunc%pwmcos_gamma
    pwmcos_eps     => enefunc%pwmcos_eps
    pwmcos_spec    => enefunc%pwmcos_specificity
    base_type      => enefunc%NA_base_type

    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

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

    num_pwmcos = 0

    !$omp parallel default(none)               &
    !$omp firstprivate(num_pwmcos)             &
    !$omp private(my_id, id,                   &
    !$omp         ini_pwmcos, fin_pwmcos,      &
    !$omp         i, j, k, l, m,               &
    !$omp         i1, i2, i3, k_pbc,           &
    !$omp         i_atom, j_atom,              &
    !$omp         ipwm,                        &
    !$omp         i_CA_0, i_CA_N, i_CA_C,      &
    !$omp         i_DB_0, i_DB_5, i_DB_3,      &
    !$omp         i_DS_0,                      &
    !$omp         j_base_type,                 &
    !$omp         pwm_shift, pwm_factor,       &
    !$omp         cutoff, cutoff_sqr,          &
    !$omp         r00, t10, t20, t30,          &
    !$omp         vbc, vbs, v53, vcn,          &
    !$omp         rbc, rbs, r53, rcn,          &
    !$omp         rbc_sqr, rbs_sqr,            &
    !$omp         r53_sqr, rcn_sqr,            &
    !$omp         rbc_inv, rbs_inv,            &
    !$omp         r53_inv, rcn_inv,            &
    !$omp         ebc, ebs, e53, ecn,          &
    !$omp         t1, t2, t3,                  &
    !$omp         cos_t1, cos_t2, cos_t3,      &
    !$omp         sin_t1, sin_t2, sin_t3,      &
    !$omp         dt1, dt2, dt3,               &
    !$omp         abs_dt1, abs_dt2, abs_dt3,   &
    !$omp         cos_dt1, cos_dt2, cos_dt3,   &
    !$omp         sin_dt1, sin_dt2, sin_dt3,   &
    !$omp         pwm_score,                   &
    !$omp         ene_coef_r0,                 &
    !$omp         grad_coef_t1, ene_coef_t1,   &
    !$omp         grad_coef_t2, ene_coef_t2,   &
    !$omp         grad_coef_t3, ene_coef_t3,   &
    !$omp         etmp0, e_tmp_pwmcos,         &
    !$omp         var_tmp, grad_r0,            &
    !$omp         grad_t1, grad_t2, grad_t3,   &
    !$omp         f_r0_coef_tmp, f_r0_tmp,     &
    !$omp         f_t1_coef_tmp, f_t1_tmp,     &
    !$omp         f_t2_coef_tmp, f_t2_tmp,     &
    !$omp         f_t3_coef_tmp, f_t3_tmp,     &
    !$omp         f_CA_0, f_CA_N, f_CA_C,      &
    !$omp         f_DB_0, f_DB_5, f_DB_3,      &
    !$omp         f_DS_0                       &
    !$omp         )                            &
    !$omp shared(coord,                        &
    !$omp        my_city_rank, nproc_city,     &
    !$omp        nthread, npwmcos,             &
    !$omp        sig_sqr,                      &
    !$omp        sigma, phi, phi2,             &
    !$omp        ktheta_2, ktheta,             &
    !$omp        cg_list_base,                 &
    !$omp        pwmcos_list, n_pwmcos_pair,   &
    !$omp        pwmcos_id, pwmcos_id_N,       &
    !$omp        pwmcos_id_C,                  &
    !$omp        pwmcos_r0, pwmcos_theta1,     &
    !$omp        pwmcos_theta2, pwmcos_theta3, &
    !$omp        pwmcos_ene_A, pwmcos_ene_C,   &
    !$omp        pwmcos_ene_G, pwmcos_ene_T,   &
    !$omp        pwmcos_gamma, pwmcos_eps,     &
    !$omp        pwmcos_spec, base_type,       &
    !$omp        bsize, half_bsize,            &
    !$omp        force)                        &
    !$omp reduction(+:virial) reduction(+:epwmcos)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do ipwm = 1, npwmcos

      ini_pwmcos = num_pwmcos + 1
      fin_pwmcos = num_pwmcos + n_pwmcos_pair(ipwm, id)
      num_pwmcos = fin_pwmcos

      if (fin_pwmcos >= ini_pwmcos) then

        i_CA_0     = pwmcos_id      (ipwm)
        i_CA_N     = pwmcos_id_N    (ipwm)
        i_CA_C     = pwmcos_id_C    (ipwm)

        r00        = pwmcos_r0      (ipwm)
        t10        = pwmcos_theta1  (ipwm)
        t20        = pwmcos_theta2  (ipwm)
        t30        = pwmcos_theta3  (ipwm)

        pwm_shift  = pwmcos_eps     (ipwm)
        pwm_factor = pwmcos_gamma   (ipwm)

        cutoff     = pwmcos_r0(ipwm) + 5.0_wp
        cutoff_sqr = cutoff * cutoff

        ! -------------
        ! Ca_C ==> Ca_N
        ! -------------
        !
        vcn(1:3)   = coord(1:3, i_CA_N) - coord(1:3, i_CA_C)
        if (vcn(1) > half_bsize(1)) then
          vcn(1) = vcn(1) - bsize(1)
        else if (vcn(1) < -half_bsize(1)) then
          vcn(1) = vcn(1) + bsize(1)
        end if
        if (vcn(2) > half_bsize(2)) then
          vcn(2) = vcn(2) - bsize(2)
        else if (vcn(2) < -half_bsize(2)) then
          vcn(2) = vcn(2) + bsize(2)
        end if
        if (vcn(3) > half_bsize(3)) then
          vcn(3) = vcn(3) - bsize(3)
        else if (vcn(3) < -half_bsize(3)) then
          vcn(3) = vcn(3) + bsize(3)
        end if
        !
        rcn_sqr    = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
        rcn        = sqrt(rcn_sqr)
        rcn_inv    = 1.0_wp / rcn
        ecn(1:3)   = vcn(1:3) * rcn_inv

        do k = ini_pwmcos, fin_pwmcos

          ! i_DB_0   = pwmcos_list(k, id)
          k_pbc = pwmcos_list(k,id)
          i_DB_0 = k_pbc / 27
          k_pbc = k_pbc - i_DB_0*27
          i3 = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2 = k_pbc / 3
          i1 = k_pbc - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          vbc(1)  = coord(1, i_CA_0) - coord(1, i_DB_0) - bsize(1) * real(i1, wp)
          vbc(2)  = coord(2, i_CA_0) - coord(2, i_DB_0) - bsize(2) * real(i2, wp)
          vbc(3)  = coord(3, i_CA_0) - coord(3, i_DB_0) - bsize(3) * real(i3, wp)
          rbc_sqr  = vbc(1) * vbc(1) + vbc(2) * vbc(2) + vbc(3) * vbc(3)

          if (rbc_sqr >= cutoff_sqr) cycle

          f_CA_0(1:3) = 0.0_wp
          f_CA_N(1:3) = 0.0_wp
          f_CA_C(1:3) = 0.0_wp
          f_DB_0(1:3) = 0.0_wp
          f_DB_5(1:3) = 0.0_wp
          f_DB_3(1:3) = 0.0_wp
          f_DS_0(1:3) = 0.0_wp

          etmp0        = 0.0_wp
          e_tmp_pwmcos = 0.0_wp

          ! ==================================
          ! Distance and angle calculations...
          ! ==================================

          j_base_type = base_type(i_DB_0)

          i_DB_5 = i_DB_0 - 3
          i_DB_3 = i_DB_0 + 3
          i_DS_0 = i_DB_0 - 1

          ! --------------
          ! sugar <== base
          ! --------------
          !
          vbs(1:3) = coord(1:3, i_DS_0) - coord(1:3, i_DB_0)
          if (vbs(1) > half_bsize(1)) then
            vbs(1) = vbs(1) - bsize(1)
          else if (vbs(1) < -half_bsize(1)) then
            vbs(1) = vbs(1) + bsize(1)
          end if
          if (vbs(2) > half_bsize(2)) then
            vbs(2) = vbs(2) - bsize(2)
          else if (vbs(2) < -half_bsize(2)) then
            vbs(2) = vbs(2) + bsize(2)
          end if
          if (vbs(3) > half_bsize(3)) then
            vbs(3) = vbs(3) - bsize(3)
          else if (vbs(3) < -half_bsize(3)) then
            vbs(3) = vbs(3) + bsize(3)
          end if
          !
          rbs_sqr  = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
          rbs      = sqrt(rbs_sqr)
          rbs_inv  = 1.0_wp / rbs
          ebs(1:3) = vbs(1:3) * rbs_inv

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
          v53(1:3) = coord(1:3, i_DB_3) - coord(1:3, i_DB_5)
          if (v53(1) > half_bsize(1)) then
            v53(1) = v53(1) - bsize(1)
          else if (v53(1) < -half_bsize(1)) then
            v53(1) = v53(1) + bsize(1)
          end if
          if (v53(2) > half_bsize(2)) then
            v53(2) = v53(2) - bsize(2)
          else if (v53(2) < -half_bsize(2)) then
            v53(2) = v53(2) + bsize(2)
          end if
          if (v53(3) > half_bsize(3)) then
            v53(3) = v53(3) - bsize(3)
          else if (v53(3) < -half_bsize(3)) then
            v53(3) = v53(3) + bsize(3)
          end if
          !
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
          t1     = acos(cos_t1)

          ! ------------------------------------------
          ! Angle t2: base5' - base3' <==> base0 -- Ca
          ! ------------------------------------------
          !
          cos_t2 = ebc(1) * e53(1) + ebc(2) * e53(2) + ebc(3) * e53(3)
          if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
          if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
          sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
          t2     = acos(cos_t2)

          ! ---------------------------------------
          ! Angle t3: base0 -- Ca <==> Ca_N -- Ca_C
          ! ---------------------------------------
          !
          cos_t3 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
          if (cos_t3 >  1.0_wp) cos_t3 =  1.0_wp
          if (cos_t3 < -1.0_wp) cos_t3 = -1.0_wp
          sin_t3 = sqrt(1.0_wp - cos_t3 * cos_t3)
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

          ! ============================================================================
          ! Energy/Force Calculation: gamma * (pwm_score + pwm_shift) * f * g1 * g2 * g3
          ! ============================================================================
          !
          ! -----------------
          ! Simple angle test
          ! -----------------
          !
          if (abs_dt1 >= phi2 .or. &
              abs_dt2 >= phi2 .or. &
              abs_dt3 >= phi2) then
            cycle
          end if


          ! ---------------------------------
          ! PWM energy: sequence specificity!
          ! ---------------------------------
          !
          if (j_base_type == NABaseTypeDBA) then
            pwm_score = pwmcos_ene_A(ipwm)
          else if (j_base_type == NABaseTypeDBC) then
            pwm_score = pwmcos_ene_C(ipwm)
          else if (j_base_type == NABaseTypeDBG) then
            pwm_score = pwmcos_ene_G(ipwm)
          else if (j_base_type == NABaseTypeDBT) then
            pwm_score = pwmcos_ene_T(ipwm)
          end if
          etmp0 = pwm_factor * (pwm_score + pwm_shift)

          ! -------------------
          ! Distance modulating
          ! -------------------
          !
          call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, grad_r0)

          ! ----------------
          ! Angle modulating
          ! ----------------
          if (abs_dt1 < phi) then
            ene_coef_t1  = 1.0_wp
            grad_coef_t1 = 0.0_wp
            grad_t1(1:6) = 0.0_wp
          else
            sin_dt1      = sin(ktheta * dt1)
            ene_coef_t1  = sin_dt1 * sin_dt1
            call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
            grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
          end if
          if (abs_dt2 < phi) then
            ene_coef_t2  = 1.0_wp
            grad_coef_t2 = 0.0_wp
            grad_t2(1:6) = 0.0_wp
          else
            sin_dt2      = sin(ktheta * dt2)
            ene_coef_t2  = sin_dt2 * sin_dt2
            call calculate_angle(vbc, v53, rbc, r53, cos_t2, sin_t2, grad_t2)
            grad_coef_t2 = ktheta * sin(2.0_wp * ktheta * dt2)
          end if
          if (abs_dt3 < phi) then
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
          e_tmp_pwmcos = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2 * ene_coef_t3
          epwmcos      = epwmcos + e_tmp_pwmcos

          ! ======
          ! Forces
          ! ======
          !
          ! --
          ! r0
          ! --
          !
          ! f_r0_coef_tmp = e_tmp_pwmcos / ene_coef_r0
          f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2 * ene_coef_t3
          f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_r0_tmp(1:3)
          !
          ! --
          ! t1
          ! --
          !
          ! f_t1_coef_tmp = e_tmp_pwmcos / ene_coef_t1 * grad_coef_t1
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
          ! f_t2_coef_tmp = e_tmp_pwmcos / ene_coef_t2 * grad_coef_t2
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
          ! f_t3_coef_tmp = e_tmp_pwmcos / ene_coef_t3 * grad_coef_t3
          f_t3_coef_tmp =  etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2 * grad_coef_t3
          f_t3_tmp(1:6) = - f_t3_coef_tmp * grad_t3(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t3_tmp(1:3)
          f_DB_0(1:3)   = f_DB_0(1:3) - f_t3_tmp(1:3)
          f_CA_C(1:3)   = f_CA_C(1:3) - f_t3_tmp(4:6)
          f_CA_N(1:3)   = f_CA_N(1:3) + f_t3_tmp(4:6)


          ! ------------------------
          ! update all the forces...
          ! ------------------------
          !
          force(1:3, i_CA_0,id) = force(1:3, i_CA_0,id) + f_CA_0(1:3)
          force(1:3, i_CA_N,id) = force(1:3, i_CA_N,id) + f_CA_N(1:3)
          force(1:3, i_CA_C,id) = force(1:3, i_CA_C,id) + f_CA_C(1:3)
          force(1:3, i_DB_0,id) = force(1:3, i_DB_0,id) + f_DB_0(1:3)
          force(1:3, i_DB_5,id) = force(1:3, i_DB_5,id) + f_DB_5(1:3)
          force(1:3, i_DB_3,id) = force(1:3, i_DB_3,id) + f_DB_3(1:3)
          force(1:3, i_DS_0,id) = force(1:3, i_DS_0,id) + f_DS_0(1:3)

          ! ------
          ! Virial
          ! ------
          do l = 1, 3
            do m = l+1, 3
              var_tmp &
                  = f_CA_0(l) * coord(m, i_CA_0) &
                  + f_CA_N(l) * coord(m, i_CA_N) &
                  + f_CA_C(l) * coord(m, i_CA_C) &
                  + f_DB_0(l) * coord(m, i_DB_0) &
                  + f_DB_5(l) * coord(m, i_DB_5) &
                  + f_DB_3(l) * coord(m, i_DB_3) &
                  + f_DS_0(l) * coord(m, i_DS_0)
              virial(m, l) = virial(m, l) - var_tmp
              virial(l, m) = virial(l, m) - var_tmp
            end do
            var_tmp &
                = f_CA_0(l) * coord(l, i_CA_0) &
                + f_CA_N(l) * coord(l, i_CA_N) &
                + f_CA_C(l) * coord(l, i_CA_C) &
                + f_DB_0(l) * coord(l, i_DB_0) &
                + f_DB_5(l) * coord(l, i_DB_5) &
                + f_DB_3(l) * coord(l, i_DB_3) &
                + f_DS_0(l) * coord(l, i_DS_0)
            virial(l,l) = virial(l,l) - var_tmp
          end do

        end do                  ! k

      end if                    ! fin_pwmcos > ini_pwmcos

    end do                      ! i_pwm
    !$omp end parallel

    call timer(TimerCGPWMcos, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_pwmcos_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_pwmcosns
  !> @brief        calculate Pwmcosns interactions between protein and DNA
  !! @authors      CT
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    pairlist   : pairlist information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] epwmcosns  : Pwmcosns energy
  !! @note         3SPN.2C + protein
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_pwmcosns(enefunc, pairlist, coord, force, &
                                        virial, epwmcosns)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epwmcosns

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: natom, npwmcosns
    integer           :: ini_pwmcosns, fin_pwmcosns
    integer           :: k, l, m
    integer           :: ipwm, ipair
    ! indicex + params
    integer           :: i_CA_0, i_CA_N, i_CA_C
    integer           :: i_DP_0
    integer           :: i_DS_0
    real(wp)          :: cutoff, cutoff_sqr
    ! native values
    real(wp)          :: r00, t10, t20
    ! vectors and distances
    real(wp)          :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)          :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)          :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)          :: t1, t2
    real(wp)          :: cos_t1, cos_t2
    real(wp)          :: sin_t1, sin_t2
    ! modulating functions
    real(wp)          :: dt1, dt2
    real(wp)          :: abs_dt1, abs_dt2
    real(wp)          :: cos_dt1, cos_dt2
    real(wp)          :: sin_dt1, sin_dt2
    real(wp)          :: sigma, phi, phi2
    real(wp)          :: sig_sqr
    real(wp)          :: ktheta_2, ktheta
    ! energy calc
    real(wp)          :: pwm_score
    real(wp)          :: ene_coef_r0
    real(wp)          :: grad_coef_t1, ene_coef_t1
    real(wp)          :: grad_coef_t2, ene_coef_t2
    real(wp)          :: etmp0
    real(wp)          :: e_tmp_pwmcosns
    ! force coef
    real(wp)          :: var_tmp
    real(wp)          :: grad_r0(3)
    real(wp)          :: grad_t1(6)
    real(wp)          :: grad_t2(6)
    ! force calc
    real(wp)          :: f_r0_coef_tmp
    real(wp)          :: f_r0_tmp(6)
    real(wp)          :: f_t1_coef_tmp
    real(wp)          :: f_t1_tmp(6)
    real(wp)          :: f_t2_coef_tmp
    real(wp)          :: f_t2_tmp(6)
    real(wp)          :: f_CA_0(3)
    real(wp)          :: f_CA_N(3)
    real(wp)          :: f_CA_C(3)
    real(wp)          :: f_DP_0(3)
    real(wp)          :: f_DS_0(3)

    integer,  pointer :: pwmcosns_list(:,:)
    integer,  pointer :: n_pwmcosns_pair(:)
    integer,  pointer :: pwmcosns_id(:)
    integer,  pointer :: pwmcosns_id_N(:)
    integer,  pointer :: pwmcosns_id_C(:)
    real(wp), pointer :: pwmcosns_r0(:)
    real(wp), pointer :: pwmcosns_theta1(:)
    real(wp), pointer :: pwmcosns_theta2(:)
    real(wp), pointer :: pwmcosns_ene(:)
    integer,  pointer :: pwmcosns_spec(:)
    integer,  pointer :: pwmcosns_to_pair(:)
    integer,  pointer :: base_type(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPwmcosns, TimerOn)

    natom   = size(coord(1,:))
    npwmcosns = enefunc%num_pwmcosns_terms

    pwmcosns_list    => pairlist%cg_pwmcosns_list
    n_pwmcosns_pair  => pairlist%num_cg_pwmcosns

    pwmcosns_id      => enefunc%pwmcosns_protein_id
    pwmcosns_id_N    => enefunc%pwmcosns_protein_id_N
    pwmcosns_id_C    => enefunc%pwmcosns_protein_id_C
    pwmcosns_r0      => enefunc%pwmcosns_r0
    pwmcosns_theta1  => enefunc%pwmcosns_theta1
    pwmcosns_theta2  => enefunc%pwmcosns_theta2
    pwmcosns_ene     => enefunc%pwmcosns_ene
    pwmcosns_spec    => enefunc%pwmcosns_specificity
    pwmcosns_to_pair => enefunc%pwmcosns_to_pairlist_id
    base_type        => enefunc%NA_base_type

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

    !$omp parallel default(none)                  &
    !$omp private(my_id, id,                      &
    !$omp         ini_pwmcosns, fin_pwmcosns,     &
    !$omp         k, l, m,                        &
    !$omp         ipwm, ipair,                    &
    !$omp         i_CA_0, i_CA_N, i_CA_C,         &
    !$omp         i_DP_0, i_DS_0,                 &
    !$omp         cutoff, cutoff_sqr,             &
    !$omp         r00, t10, t20,                  &
    !$omp         vbc, vbs, vcn,                  &
    !$omp         rbc, rbs, rcn,                  &
    !$omp         rbc_sqr, rbs_sqr,               &
    !$omp         rcn_sqr,                        &
    !$omp         rbc_inv, rbs_inv,               &
    !$omp         rcn_inv,                        &
    !$omp         ebc, ebs, ecn,                  &
    !$omp         t1, t2,                         &
    !$omp         cos_t1, cos_t2,                 &
    !$omp         sin_t1, sin_t2,                 &
    !$omp         dt1, dt2,                       &
    !$omp         abs_dt1, abs_dt2,               &
    !$omp         cos_dt1, cos_dt2,               &
    !$omp         sin_dt1, sin_dt2,               &
    !$omp         pwm_score,                      &
    !$omp         ene_coef_r0,                    &
    !$omp         grad_coef_t1, ene_coef_t1,      &
    !$omp         grad_coef_t2, ene_coef_t2,      &
    !$omp         etmp0, e_tmp_pwmcosns,          &
    !$omp         var_tmp, grad_r0,               &
    !$omp         grad_t1, grad_t2,               &
    !$omp         f_r0_coef_tmp, f_r0_tmp,        &
    !$omp         f_t1_coef_tmp, f_t1_tmp,        &
    !$omp         f_t2_coef_tmp, f_t2_tmp,        &
    !$omp         f_CA_0, f_CA_N, f_CA_C,         &
    !$omp         f_DP_0,                         &
    !$omp         f_DS_0                          &
    !$omp         )                               &
    !$omp shared(coord, my_city_rank, nproc_city, &
    !$omp        nthread, natom, npwmcosns,       &
    !$omp        sig_sqr,                         &
    !$omp        sigma, phi, phi2,                &
    !$omp        ktheta_2, ktheta,                &
    !$omp        pwmcosns_list, n_pwmcosns_pair,  &
    !$omp        pwmcosns_id, pwmcosns_id_N,      &
    !$omp        pwmcosns_id_C, pwmcosns_to_pair, &
    !$omp        pwmcosns_r0, pwmcosns_theta1,    &
    !$omp        pwmcosns_theta2,                 &
    !$omp        pwmcosns_ene,                    &
    !$omp        pwmcosns_spec, base_type,        &
    !$omp        force)                           &
    !$omp reduction(+:virial) reduction(+:epwmcosns)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do ipwm = 1, npwmcosns

      ! if (mod(ipwm - 1, nproc_city*nthread) .ne. my_id) cycle

      i_CA_0     = pwmcosns_id      (ipwm)
      ipair      = pwmcosns_to_pair (ipwm)

      if (mod(ipair - 1, nproc_city*nthread) .ne. my_id) cycle

      i_CA_N     = pwmcosns_id_N    (ipwm)
      i_CA_C     = pwmcosns_id_C    (ipwm)

      r00        = pwmcosns_r0      (ipwm)
      t10        = pwmcosns_theta1  (ipwm)
      t20        = pwmcosns_theta2  (ipwm)

      cutoff     = pwmcosns_r0(ipwm) + 5.0
      cutoff_sqr = cutoff * cutoff

      ! -------------
      ! Ca_C ==> Ca_N
      ! -------------
      !
      vcn(1:3)   = coord(1:3, i_CA_N) - coord(1:3, i_CA_C)
      rcn_sqr    = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
      rcn        = sqrt(rcn_sqr)
      rcn_inv    = 1.0_wp / rcn
      ecn(1:3)   = vcn(1:3) * rcn_inv

      ini_pwmcosns = 1
      fin_pwmcosns = n_pwmcosns_pair(ipair)

      do k = ini_pwmcosns, fin_pwmcosns

        i_DP_0   = pwmcosns_list(k, ipair)

        vbc(1:3) = coord(1:3, i_CA_0) - coord(1:3, i_DP_0)
        rbc_sqr  = vbc(1) * vbc(1) + vbc(2) * vbc(2) + vbc(3) * vbc(3)

        if (rbc_sqr >= cutoff_sqr) cycle

        f_CA_0(1:3) = 0.0_wp
        f_CA_N(1:3) = 0.0_wp
        f_CA_C(1:3) = 0.0_wp
        f_DP_0(1:3) = 0.0_wp
        f_DS_0(1:3) = 0.0_wp

        etmp0        = 0.0_wp
        e_tmp_pwmcosns = 0.0_wp

        ! ==================================
        ! Distance and angle calculations...
        ! ==================================

        i_DS_0 = i_DP_0 + 1

        ! --------------
        ! phos ==> sugar
        ! --------------
        !
        vbs(1:3) = coord(1:3, i_DS_0) - coord(1:3, i_DP_0)
        rbs_sqr  = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
        rbs      = sqrt(rbs_sqr)
        rbs_inv  = 1.0_wp / rbs
        ebs(1:3) = vbs(1:3) * rbs_inv

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
        t1     = acos(cos_t1)

        ! ---------------------------------------
        ! Angle t2: base0 -- Ca <==> Ca_N -- Ca_C
        ! ---------------------------------------
        !
        cos_t2 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
        if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
        if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
        sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
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
        !
        ! -----------------
        ! Simple angle test
        ! -----------------
        !
        if (abs_dt1 >= phi2 .or. &
            abs_dt2 >= phi2) then
          cycle
        end if

        ! ---------------------------------
        ! PWM energy: sequence specificity!
        ! ---------------------------------
        !
        etmp0 = pwmcosns_ene(ipwm)

        ! -------------------
        ! Distance modulating
        ! -------------------
        !
        call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, grad_r0)

        ! ----------------
        ! Angle modulating
        ! ----------------
        if (abs_dt1 < phi) then
          ene_coef_t1  = 1.0_wp
          grad_coef_t1 = 0.0_wp
          grad_t1(1:6) = 0.0_wp
        else
          sin_dt1      = sin(ktheta * dt1)
          ene_coef_t1  = sin_dt1 * sin_dt1
          call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
          grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
        end if

        if (abs_dt2 < phi) then
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
        e_tmp_pwmcosns = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2
        epwmcosns      = epwmcosns + e_tmp_pwmcosns

        ! ======
        ! Forces
        ! ======
        !
        ! --
        ! r0
        ! --
        !
        ! f_r0_coef_tmp = e_tmp_pwmcosns / ene_coef_r0
        f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2
        f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
        f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
        f_DP_0(1:3)   = f_DP_0(1:3) - f_r0_tmp(1:3)
        !
        ! --
        ! t1
        ! --
        !
        ! f_t1_coef_tmp = e_tmp_pwmcosns / ene_coef_t1 * grad_coef_t1
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
        ! f_t2_coef_tmp = e_tmp_pwmcosns / ene_coef_t2 * grad_coef_t2
        f_t2_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t1 * grad_coef_t2
        f_t2_tmp(1:6) = - f_t2_coef_tmp * grad_t2(1:6)
        f_CA_0(1:3)   = f_CA_0(1:3) + f_t2_tmp(1:3)
        f_DP_0(1:3)   = f_DP_0(1:3) - f_t2_tmp(1:3)
        f_CA_C(1:3)   = f_CA_C(1:3) - f_t2_tmp(4:6)
        f_CA_N(1:3)   = f_CA_N(1:3) + f_t2_tmp(4:6)


        ! ------------------------
        ! update all the forces...
        ! ------------------------
        !
        force(1:3, i_CA_0,id) = force(1:3, i_CA_0,id) + f_CA_0(1:3)
        force(1:3, i_CA_N,id) = force(1:3, i_CA_N,id) + f_CA_N(1:3)
        force(1:3, i_CA_C,id) = force(1:3, i_CA_C,id) + f_CA_C(1:3)
        force(1:3, i_DP_0,id) = force(1:3, i_DP_0,id) + f_DP_0(1:3)
        force(1:3, i_DS_0,id) = force(1:3, i_DS_0,id) + f_DS_0(1:3)

        ! ------
        ! Virial
        ! ------
        do l = 1, 3
          do m = l+1, 3
            var_tmp &
                = f_CA_0(l) * coord(m, i_CA_0) &
                + f_CA_N(l) * coord(m, i_CA_N) &
                + f_CA_C(l) * coord(m, i_CA_C) &
                + f_DP_0(l) * coord(m, i_DP_0) &
                + f_DS_0(l) * coord(m, i_DS_0)
            virial(m, l) = virial(m, l) - var_tmp
            virial(l, m) = virial(l, m) - var_tmp
          end do
          var_tmp &
              = f_CA_0(l) * coord(l, i_CA_0) &
              + f_CA_N(l) * coord(l, i_CA_N) &
              + f_CA_C(l) * coord(l, i_CA_C) &
              + f_DP_0(l) * coord(l, i_DP_0) &
              + f_DS_0(l) * coord(l, i_DS_0)
          virial(l,l) = virial(l,l) - var_tmp
        end do

      end do

    end do
    !$omp end parallel

    call timer(TimerCGPwmcosns, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_pwmcosns

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_CG_pwmcosns_pbc
  !> @brief        calculate Pwmcosns interactions between protein and DNA (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : PBC coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] epwmcosns  : Pwmcosns energy
  !! @note         3SPN.2C + protein
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_CG_pwmcosns_pbc(enefunc, boundary, pairlist, &
                                            coord, force, virial, epwmcosns)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: epwmcosns

    ! local variables
    integer           :: my_id, id
    integer           :: omp_get_num_threads, omp_get_thread_num
    integer           :: npwmcosns
    integer           :: num_pwmcosns, ini_pwmcosns, fin_pwmcosns
    integer           :: i, j, k, l, m
    integer           :: i1, i2, i3, k_pbc ! for PBC
    integer           :: i_atom, j_atom
    integer           :: ipwm
    ! indicex + params
    integer           :: i_CA_0, i_CA_N, i_CA_C
    integer           :: i_DP_0
    integer           :: i_DS_0
    real(wp)          :: cutoff, cutoff_sqr
    ! native values
    real(wp)          :: r00, t10, t20
    ! vectors and distances
    real(wp)          :: vbc(3), rbc, rbc_sqr, rbc_inv, ebc(3)
    real(wp)          :: vbs(3), rbs, rbs_sqr, rbs_inv, ebs(3)
    real(wp)          :: vcn(3), rcn, rcn_sqr, rcn_inv, ecn(3)
    ! angles
    real(wp)          :: t1, t2
    real(wp)          :: cos_t1, cos_t2
    real(wp)          :: sin_t1, sin_t2
    ! modulating functions
    real(wp)          :: dt1, dt2
    real(wp)          :: abs_dt1, abs_dt2
    real(wp)          :: cos_dt1, cos_dt2
    real(wp)          :: sin_dt1, sin_dt2
    real(wp)          :: sigma, phi, phi2
    real(wp)          :: sig_sqr
    real(wp)          :: ktheta_2, ktheta
    ! energy calc
    real(wp)          :: ene_coef_r0
    real(wp)          :: grad_coef_t1, ene_coef_t1
    real(wp)          :: grad_coef_t2, ene_coef_t2
    real(wp)          :: etmp0
    real(wp)          :: e_tmp_pwmcosns
    ! force coef
    real(wp)          :: var_tmp
    real(wp)          :: grad_r0(3)
    real(wp)          :: grad_t1(6)
    real(wp)          :: grad_t2(6)
    ! force calc
    real(wp)          :: f_r0_coef_tmp
    real(wp)          :: f_r0_tmp(6)
    real(wp)          :: f_t1_coef_tmp
    real(wp)          :: f_t1_tmp(6)
    real(wp)          :: f_t2_coef_tmp
    real(wp)          :: f_t2_tmp(6)
    real(wp)          :: f_CA_0(3)
    real(wp)          :: f_CA_N(3)
    real(wp)          :: f_CA_C(3)
    real(wp)          :: f_DP_0(3)
    real(wp)          :: f_DS_0(3)

    integer,  pointer :: cg_list_base(:)
    integer,  pointer :: pwmcosns_list(:,:)
    integer,  pointer :: n_pwmcosns_pair(:, :)
    integer,  pointer :: pwmcosns_id(:)
    integer,  pointer :: pwmcosns_id_N(:)
    integer,  pointer :: pwmcosns_id_C(:)
    real(wp), pointer :: pwmcosns_r0(:)
    real(wp), pointer :: pwmcosns_theta1(:)
    real(wp), pointer :: pwmcosns_theta2(:)
    real(wp), pointer :: pwmcosns_ene(:)
    integer,  pointer :: pwmcosns_spec(:)
    integer,  pointer :: base_type(:)

    real(wp)          :: bsize(3), half_bsize(3)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGPwmcosns, TimerOn)

    npwmcosns        = enefunc%num_pwmcosns_terms

    pwmcosns_list    => pairlist%cg_pwmcosns_list
    n_pwmcosns_pair  => pairlist%num_cg_pwmcosns_calc
    cg_list_base     => enefunc%cg_particle_DNA_base

    pwmcosns_id      => enefunc%pwmcosns_protein_id
    pwmcosns_id_N    => enefunc%pwmcosns_protein_id_N
    pwmcosns_id_C    => enefunc%pwmcosns_protein_id_C
    pwmcosns_r0      => enefunc%pwmcosns_r0
    pwmcosns_theta1  => enefunc%pwmcosns_theta1
    pwmcosns_theta2  => enefunc%pwmcosns_theta2
    pwmcosns_ene     => enefunc%pwmcosns_ene
    pwmcosns_spec    => enefunc%pwmcosns_specificity
    base_type        => enefunc%NA_base_type

    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

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

    num_pwmcosns = 0

    !$omp parallel default(none)                 &
    !$omp firstprivate(num_pwmcosns)             &
    !$omp private(my_id, id,                     &
    !$omp         ini_pwmcosns, fin_pwmcosns,    &
    !$omp         i, j, k, l, m,                 &
    !$omp         i1, i2, i3, k_pbc,             &
    !$omp         i_atom, j_atom,                &
    !$omp         ipwm,                          &
    !$omp         i_CA_0, i_CA_N, i_CA_C,        &
    !$omp         i_DP_0, i_DS_0,                &
    !$omp         cutoff, cutoff_sqr,            &
    !$omp         r00, t10, t20,                 &
    !$omp         vbc, vbs, vcn,                 &
    !$omp         rbc, rbs, rcn,                 &
    !$omp         rbc_sqr, rbs_sqr,              &
    !$omp         rcn_sqr,                       &
    !$omp         rbc_inv, rbs_inv,              &
    !$omp         rcn_inv,                       &
    !$omp         ebc, ebs, ecn,                 &
    !$omp         t1, t2,                        &
    !$omp         cos_t1, cos_t2,                &
    !$omp         sin_t1, sin_t2,                &
    !$omp         dt1, dt2,                      &
    !$omp         abs_dt1, abs_dt2,              &
    !$omp         cos_dt1, cos_dt2,              &
    !$omp         sin_dt1, sin_dt2,              &
    !$omp         ene_coef_r0,                   &
    !$omp         grad_coef_t1, ene_coef_t1,     &
    !$omp         grad_coef_t2, ene_coef_t2,     &
    !$omp         etmp0, e_tmp_pwmcosns,         &
    !$omp         var_tmp, grad_r0,              &
    !$omp         grad_t1, grad_t2,              &
    !$omp         f_r0_coef_tmp, f_r0_tmp,       &
    !$omp         f_t1_coef_tmp, f_t1_tmp,       &
    !$omp         f_t2_coef_tmp, f_t2_tmp,       &
    !$omp         f_CA_0, f_CA_N, f_CA_C,        &
    !$omp         f_DP_0,                        &
    !$omp         f_DS_0                         &
    !$omp         )                              &
    !$omp shared(coord,                          &
    !$omp        my_city_rank, nproc_city,       &
    !$omp        nthread, npwmcosns,             &
    !$omp        sig_sqr,                        &
    !$omp        sigma, phi, phi2,               &
    !$omp        ktheta_2, ktheta,               &
    !$omp        cg_list_base,                   &
    !$omp        pwmcosns_list, n_pwmcosns_pair, &
    !$omp        pwmcosns_id, pwmcosns_id_N,     &
    !$omp        pwmcosns_id_C,                  &
    !$omp        pwmcosns_r0, pwmcosns_theta1,   &
    !$omp        pwmcosns_theta2,                &
    !$omp        pwmcosns_ene,                   &
    !$omp        pwmcosns_spec, base_type,       &
    !$omp        bsize, half_bsize,              &
    !$omp        force)                          &
    !$omp reduction(+:virial) reduction(+:epwmcosns)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do ipwm = 1, npwmcosns

      ini_pwmcosns = num_pwmcosns + 1
      fin_pwmcosns = num_pwmcosns + n_pwmcosns_pair(ipwm, id)
      num_pwmcosns = fin_pwmcosns

      if (fin_pwmcosns >= ini_pwmcosns) then

        i_CA_0     = pwmcosns_id      (ipwm)
        i_CA_N     = pwmcosns_id_N    (ipwm)
        i_CA_C     = pwmcosns_id_C    (ipwm)

        r00        = pwmcosns_r0      (ipwm)
        t10        = pwmcosns_theta1  (ipwm)
        t20        = pwmcosns_theta2  (ipwm)

        cutoff     = pwmcosns_r0(ipwm) + 5.0_wp
        cutoff_sqr = cutoff * cutoff

        etmp0      = pwmcosns_ene(ipwm)

        ! -------------
        ! Ca_C ==> Ca_N
        ! -------------
        !
        vcn(1:3)   = coord(1:3, i_CA_N) - coord(1:3, i_CA_C)
        if (vcn(1) > half_bsize(1)) then
          vcn(1) = vcn(1) - bsize(1)
        else if (vcn(1) < -half_bsize(1)) then
          vcn(1) = vcn(1) + bsize(1)
        end if
        if (vcn(2) > half_bsize(2)) then
          vcn(2) = vcn(2) - bsize(2)
        else if (vcn(2) < -half_bsize(2)) then
          vcn(2) = vcn(2) + bsize(2)
        end if
        if (vcn(3) > half_bsize(3)) then
          vcn(3) = vcn(3) - bsize(3)
        else if (vcn(3) < -half_bsize(3)) then
          vcn(3) = vcn(3) + bsize(3)
        end if
        !
        rcn_sqr    = vcn(1) * vcn(1) + vcn(2) * vcn(2) + vcn(3) * vcn(3)
        rcn        = sqrt(rcn_sqr)
        rcn_inv    = 1.0_wp / rcn
        ecn(1:3)   = vcn(1:3) * rcn_inv

        do k = ini_pwmcosns, fin_pwmcosns

          ! i_DP_0   = pwmcosns_list(k, id)
          k_pbc = pwmcosns_list(k,id)
          i_DP_0 = k_pbc / 27
          k_pbc = k_pbc - i_DP_0*27
          i3 = k_pbc / 9
          k_pbc = k_pbc - i3*9
          i2 = k_pbc / 3
          i1 = k_pbc - i2*3
          i1 = i1 - 1
          i2 = i2 - 1
          i3 = i3 - 1

          vbc(1)  = coord(1, i_CA_0) - coord(1, i_DP_0) - bsize(1) * real(i1, wp)
          vbc(2)  = coord(2, i_CA_0) - coord(2, i_DP_0) - bsize(2) * real(i2, wp)
          vbc(3)  = coord(3, i_CA_0) - coord(3, i_DP_0) - bsize(3) * real(i3, wp)
          rbc_sqr = vbc(1) * vbc(1) + vbc(2) * vbc(2) + vbc(3) * vbc(3)

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

          i_DS_0 = i_DP_0 + 1

          ! --------------
          ! sugar <== base
          ! --------------
          !
          vbs(1:3) = coord(1:3, i_DS_0) - coord(1:3, i_DP_0)
          if (vbs(1) > half_bsize(1)) then
            vbs(1) = vbs(1) - bsize(1)
          else if (vbs(1) < -half_bsize(1)) then
            vbs(1) = vbs(1) + bsize(1)
          end if
          if (vbs(2) > half_bsize(2)) then
            vbs(2) = vbs(2) - bsize(2)
          else if (vbs(2) < -half_bsize(2)) then
            vbs(2) = vbs(2) + bsize(2)
          end if
          if (vbs(3) > half_bsize(3)) then
            vbs(3) = vbs(3) - bsize(3)
          else if (vbs(3) < -half_bsize(3)) then
            vbs(3) = vbs(3) + bsize(3)
          end if
          !
          rbs_sqr  = vbs(1) * vbs(1) + vbs(2) * vbs(2) + vbs(3) * vbs(3)
          rbs      = sqrt(rbs_sqr)
          rbs_inv  = 1.0_wp / rbs
          ebs(1:3) = vbs(1:3) * rbs_inv

          ! -----------
          ! base ==> Ca
          ! -----------
          !
          rbc      = sqrt(rbc_sqr)
          rbc_inv  = 1.0_wp / rbc
          ebc(1:3) = vbc(1:3) * rbc_inv

          ! -----------------------------
          ! Angle t1: sugar -- base -- Ca
          ! -----------------------------
          !
          cos_t1 = ebc(1) * ebs(1) + ebc(2) * ebs(2) + ebc(3) * ebs(3)
          if (cos_t1 >  1.0_wp) cos_t1 =  1.0_wp
          if (cos_t1 < -1.0_wp) cos_t1 = -1.0_wp
          sin_t1 = sqrt(1.0_wp - cos_t1 * cos_t1)
          t1     = acos(cos_t1)

          ! ---------------------------------------
          ! Angle t2: base0 -- Ca <==> Ca_N -- Ca_C
          ! ---------------------------------------
          !
          cos_t2 = ebc(1) * ecn(1) + ebc(2) * ecn(2) + ebc(3) * ecn(3)
          if (cos_t2 >  1.0_wp) cos_t2 =  1.0_wp
          if (cos_t2 < -1.0_wp) cos_t2 = -1.0_wp
          sin_t2 = sqrt(1.0_wp - cos_t2 * cos_t2)
          t2     = acos(cos_t2)

          ! ---------------
          ! Delta angles...
          ! ---------------
          !
          dt1     = t1 - t10
          dt2     = t2 - t20
          abs_dt1 = abs(dt1)
          abs_dt2 = abs(dt2)


          ! =======================================================================
          ! Energy/Force Calculation: gamma * (pwm_score + pwm_shift) * f * g1 * g2
          ! =======================================================================
          !
          ! -----------------
          ! Simple angle test
          ! -----------------
          !
          if (abs_dt1 >= phi2 .or. &
              abs_dt2 >= phi2) then
            cycle
          end if

          ! -------------------
          ! Distance modulating
          ! -------------------
          !
          call calculate_gaussian(vbc, rbc, r00, sig_sqr, ene_coef_r0, grad_r0)

          ! ----------------
          ! Angle modulating
          ! ----------------
          if (abs_dt1 < phi) then
            ene_coef_t1  = 1.0_wp
            grad_coef_t1 = 0.0_wp
            grad_t1(1:6) = 0.0_wp
          else
            sin_dt1      = sin(ktheta * dt1)
            ene_coef_t1  = sin_dt1 * sin_dt1
            call calculate_angle(vbc, vbs, rbc, rbs, cos_t1, sin_t1, grad_t1)
            grad_coef_t1 = ktheta * sin(2.0_wp * ktheta * dt1)
          end if
          if (abs_dt2 < phi) then
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
          e_tmp_pwmcosns = etmp0 * ene_coef_r0 * ene_coef_t1 * ene_coef_t2
          epwmcosns      = epwmcosns + e_tmp_pwmcosns

          ! ======
          ! Forces
          ! ======
          !
          ! --
          ! r0
          ! --
          !
          ! f_r0_coef_tmp = e_tmp_pwmcosns / ene_coef_r0
          f_r0_coef_tmp = etmp0 * ene_coef_t1 * ene_coef_t2
          f_r0_tmp(1:3) = - f_r0_coef_tmp * grad_r0(1:3)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_r0_tmp(1:3)
          f_DP_0(1:3)   = f_DP_0(1:3) - f_r0_tmp(1:3)
          !
          ! --
          ! t1
          ! --
          !
          ! f_t1_coef_tmp = e_tmp_pwmcosns / ene_coef_t1 * grad_coef_t1
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
          ! f_t2_coef_tmp = e_tmp_pwmcosns / ene_coef_t2 * grad_coef_t2
          f_t2_coef_tmp = etmp0 * ene_coef_r0 * ene_coef_t1 * grad_coef_t2
          f_t2_tmp(1:6) = - f_t2_coef_tmp * grad_t2(1:6)
          f_CA_0(1:3)   = f_CA_0(1:3) + f_t2_tmp(1:3)
          f_DP_0(1:3)   = f_DP_0(1:3) - f_t2_tmp(1:3)
          f_CA_C(1:3)   = f_CA_C(1:3) - f_t2_tmp(4:6)
          f_CA_N(1:3)   = f_CA_N(1:3) + f_t2_tmp(4:6)


          ! ------------------------
          ! update all the forces...
          ! ------------------------
          !
          force(1:3, i_CA_0,id) = force(1:3, i_CA_0,id) + f_CA_0(1:3)
          force(1:3, i_CA_N,id) = force(1:3, i_CA_N,id) + f_CA_N(1:3)
          force(1:3, i_CA_C,id) = force(1:3, i_CA_C,id) + f_CA_C(1:3)
          force(1:3, i_DP_0,id) = force(1:3, i_DP_0,id) + f_DP_0(1:3)
          force(1:3, i_DS_0,id) = force(1:3, i_DS_0,id) + f_DS_0(1:3)

          ! ------
          ! Virial
          ! ------
          do l = 1, 3
            do m = l+1, 3
              var_tmp &
                  = f_CA_0(l) * coord(m, i_CA_0) &
                  + f_CA_N(l) * coord(m, i_CA_N) &
                  + f_CA_C(l) * coord(m, i_CA_C) &
                  + f_DP_0(l) * coord(m, i_DP_0) &
                  + f_DS_0(l) * coord(m, i_DS_0)
              virial(m, l) = virial(m, l) - var_tmp
              virial(l, m) = virial(l, m) - var_tmp
            end do
            var_tmp &
                = f_CA_0(l) * coord(l, i_CA_0) &
                + f_CA_N(l) * coord(l, i_CA_N) &
                + f_CA_C(l) * coord(l, i_CA_C) &
                + f_DP_0(l) * coord(l, i_DP_0) &
                + f_DS_0(l) * coord(l, i_DS_0)
            virial(l,l) = virial(l,l) - var_tmp
          end do

        end do

      end if

    end do
    !$omp end parallel

    call timer(TimerCGPwmcosns, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_CG_pwmcosns_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_attractive
  !> @brief        calculate 3SPN.2C morse ( potential ) gradients
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
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: dist
    real(wp),                intent(in)    :: sigma
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: epsilon
    real(wp),                intent(inout) :: ene
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    real(wp)                 :: expon, exptmp, etmp
    real(wp)                 :: f_coef


    if (dist > sigma) then
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

    return

  end subroutine calculate_attractive

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_nonlocal_dihedral
  !> @brief        calculate dihedral angle
  !! @authors      CK; CT
  !! @param[in]    aindex  : atom indices
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] cos_dih : cosin of dihedral angles
  !! @param[inout] sin_dih : sin of dihedral angles
  !! @param[inout] grad    : gradient of dihedral angles
  !
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !! @note         Copied from src/lib/dihedral_libs.fpp
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_nonlocal_dihedral(aindex, coord, cos_dih, sin_dih, grad)

    ! formal arguments
    integer,   intent(in)    :: aindex(:)
    real(wp),  intent(in)    :: coord(:,:)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(1:9)

    ! local variable
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl


    dij(1:3) = coord(1:3,aindex(1)) - coord(1:3,aindex(2))
    djk(1:3) = coord(1:3,aindex(2)) - coord(1:3,aindex(3))
    dlk(1:3) = coord(1:3,aindex(4)) - coord(1:3,aindex(3))

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
  !  Subroutine    calculate_nonlocal_dihedral_pbc
  !> @brief        calculate dihedral angle considering PBC
  !! @authors      CK; CT
  !! @param[in]    aindex  : atom indices
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    bsize   : PBC box size
  !! @param[inout] cos_dih : cosin of dihedral angles
  !! @param[inout] sin_dih : sin of dihedral angles
  !! @param[inout] grad    : gradient of dihedral angles
  !
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !! @note         Copied from src/lib/dihedral_libs.fpp
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_nonlocal_dihedral_pbc(aindex, coord, bsize, &
                                             cos_dih, sin_dih, grad)

    ! formal arguments
    integer,   intent(in)    :: aindex(:)
    real(wp),  intent(in)    :: coord(:,:)
    real(wp),  intent(in)    :: bsize(:)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(1:9)

    ! local variable
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl
    real(wp)                 :: half_bsize(3)


    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    dij(1:3) = coord(1:3,aindex(1)) - coord(1:3,aindex(2))
    if (dij(1) > half_bsize(1)) then
      dij(1) = dij(1) - bsize(1)
    else if (dij(1) < -half_bsize(1)) then
      dij(1) = dij(1) + bsize(1)
    end if
    if (dij(2) > half_bsize(2)) then
      dij(2) = dij(2) - bsize(2)
    else if (dij(2) < -half_bsize(2)) then
      dij(2) = dij(2) + bsize(2)
    end if
    if (dij(3) > half_bsize(3)) then
      dij(3) = dij(3) - bsize(3)
    else if (dij(3) < -half_bsize(3)) then
      dij(3) = dij(3) + bsize(3)
    end if

    djk(1:3) = coord(1:3,aindex(2)) - coord(1:3,aindex(3))
    if (djk(1) > half_bsize(1)) then
      djk(1) = djk(1) - bsize(1)
    else if (djk(1) < -half_bsize(1)) then
      djk(1) = djk(1) + bsize(1)
    end if
    if (djk(2) > half_bsize(2)) then
      djk(2) = djk(2) - bsize(2)
    else if (djk(2) < -half_bsize(2)) then
      djk(2) = djk(2) + bsize(2)
    end if
    if (djk(3) > half_bsize(3)) then
      djk(3) = djk(3) - bsize(3)
    else if (djk(3) < -half_bsize(3)) then
      djk(3) = djk(3) + bsize(3)
    end if

    dlk(1:3) = coord(1:3,aindex(4)) - coord(1:3,aindex(3))
    if (dlk(1) > half_bsize(1)) then
      dlk(1) = dlk(1) - bsize(1)
    else if (dlk(1) < -half_bsize(1)) then
      dlk(1) = dlk(1) + bsize(1)
    end if
    if (dlk(2) > half_bsize(2)) then
      dlk(2) = dlk(2) - bsize(2)
    else if (dlk(2) < -half_bsize(2)) then
      dlk(2) = dlk(2) + bsize(2)
    end if
    if (dlk(3) > half_bsize(3)) then
      dlk(3) = dlk(3) - bsize(3)
    else if (dlk(3) < -half_bsize(3)) then
      dlk(3) = dlk(3) + bsize(3)
    end if

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

  end subroutine calculate_nonlocal_dihedral_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_angle
  !> @brief        calculate angle gradients
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
    real(wp) :: sin_ang_safe

    if (sin_ang < EPS) then
      sin_ang_safe = EPS
    else
      sin_ang_safe = sin_ang
    end if

    tmp1 = 1.0_wp  / (sin_ang_safe * r1 * r2)
    tmp2 = cos_ang / (sin_ang_safe * r1 * r1)
    tmp3 = cos_ang / (sin_ang_safe * r2 * r2)

    grad(1:3) = - tmp1 * dr2(1:3) + tmp2 * dr1(1:3)
    grad(4:6) = - tmp1 * dr1(1:3) + tmp3 * dr2(1:3)

    return

  end subroutine calculate_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_gaussian
  !> @brief        calculate gaussian ( potential ) gradients
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
    real(wp),                intent(in)    :: vec(:)
    real(wp),                intent(in)    :: r
    real(wp),                intent(in)    :: r0
    real(wp),                intent(in)    :: sigma_sqr
    real(wp),                intent(inout) :: enecoef
    real(wp),                intent(inout) :: grad(1:3)

    ! local variables
    real(wp)                 :: delta_r
    real(wp)                 :: expon
    real(wp)                 :: f_coef


    delta_r = r - r0
    expon  = - 0.5_wp * delta_r * delta_r / sigma_sqr
    enecoef   = exp(expon)
    f_coef = - delta_r * enecoef / sigma_sqr / r
    grad(1:3) = f_coef * vec(1:3)

    return

  end subroutine calculate_gaussian

end module at_energy_cg_nonlocal_mod
