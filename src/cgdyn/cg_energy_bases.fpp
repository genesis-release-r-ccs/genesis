!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_bases_mod
!> @brief   Calculates base-base interaction energy.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_bases_mod

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
  public :: compute_energy_DNA_base_stacking

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_DNA_base_stacking
  !> @brief        Calculates base stacking energy.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] ebase   : base energy of target systems
  !! @note         3SPN.2C
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_DNA_base_stacking(domain, enefunc, coord, force, &
                                              ebase, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: ebase(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    integer                  :: i, j, k, ix
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
    real(wp)                 :: grad_bell_coef_0, grad_bell_coef_1
    real(wp)                 :: grad_bell_coef_3
    real(wp)                 :: work(1:9), ebase_temp, viri(1:3)
    real(wp)                 :: g_attr_1(1:3), g_attr_3(1:3)
    real(wp)                 :: theta_bs_threshold_1, theta_bs_threshold_2

    integer                  :: icel1, icel2, icel3, i1, i2, i3
    integer                  :: start_icel1, start_icel2, start_icel3
    integer                  :: id

    real(wp),        pointer :: ene_coef(:), sigma(:), theta_bs(:)
    integer,         pointer :: list(:,:), func(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)

    call timer(TimerBaseStack, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    list        => enefunc%base_stack_list
    func        => enefunc%base_stack_func
    ene_coef    => enefunc%base_stack_epsilon
    sigma       => enefunc%base_stack_sigma
    theta_bs    => enefunc%base_stack_theta_bs

    K_theta     =  enefunc%base_stack_k
    alpha       =  enefunc%base_stack_alpha

    theta_bs_threshold_2 = PI / K_theta
    theta_bs_threshold_1 = theta_bs_threshold_2 / 2.0_wp

    ! calculate base stack energy
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, ix, icel1, icel2, icel3, i1, i2, i3, start_icel1,   &
    !$omp         start_icel2, start_icel3, j, k, i_sugar, i_base5, i_base3, &
    !$omp         v21, r21_sqr, r21, v23, r23_sqr, r23, dot_v21_v23,         &
    !$omp         cos_t_bs, sin_t_bs, delta_t_bs, abs_delta_t_bs, delta_r23, &
    !$omp         expnt, exp_expnt, m_exp_expnt, ene_morse, grad_morse,      &
    !$omp         ene_repl, grad_repl, ene_attr, grad_attr, cos_dt,          &
    !$omp         coef_bell, grad_bell_coef_0, grad_bell_coef_1, viri,       &
    !$omp         grad_bell_coef_3, g_attr_1, g_attr_3, work, ebase_temp)      
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    ebase_temp = 0.0_wp

    do i = id+1, enefunc%num_stack_domain, nthread

      ! initialization of work
      !
      work(1:9) = 0.0_wp

      ! indices
      !
      i1    = id_g2l(list(1,i))
      i2    = id_g2l(list(2,i))
      i3    = id_g2l(list(3,i))

      ! r21: distance between sugar and base5
      !
      v21(1:3) = coord(i1,1:3) - coord(i2,1:3)
      r21_sqr  = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3) 
      r21      = sqrt(r21_sqr)

      ! distance between base5 and base3
      v23(1:3) = coord(i3,1:3) - coord(i2,1:3)
      r23_sqr  = v23(1)*v23(1) + v23(2)*v23(2) + v23(3)*v23(3)
      r23      = sqrt(r23_sqr)

      ! cos_t_bs: cos(theta_bs)
      ! sin_t_bs: sin(theta_bs)
      !
      dot_v21_v23 = v21(1) * v23(1) + v21(2) * v23(2) + v21(3) * v23(3)
      cos_t_bs    = dot_v21_v23 / (r21 * r23)
      if ( cos_t_bs > 1.0_wp )  cos_t_bs =  1.0_wp
      if ( cos_t_bs < -1.0_wp ) cos_t_bs = -1.0_wp
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
      grad_morse  = 2.0_wp * alpha * ene_coef(i) * exp_expnt * m_exp_expnt / r23

      ! ---------
      ! Repulsion
      ! ---------
      if ( delta_r23 < 0.0_wp ) then
        ene_repl  = ene_morse
        grad_repl = grad_morse
      else
        ene_repl  = 0.0_wp
        grad_repl = 0.0_wp
      end if
      work(4:6) =   grad_repl * v23(1:3)
      work(7:9) = - grad_repl * v23(1:3)

      ! ----------
      ! Attraction
      ! ----------
      if ( delta_r23 < 0.0_wp ) then
        ene_attr  = - ene_coef(i)
        grad_attr = 0.0_wp
      else
        ene_attr  = ene_morse - ene_coef(i)
        grad_attr = grad_morse
      end if

      if (abs_delta_t_bs < theta_bs_threshold_1) then
        coef_bell = 1.0_wp
        work(4:6) = work(4:6) + grad_attr * v23(1:3)
        work(7:9) = work(7:9) - grad_attr * v23(1:3)
      else if (abs_delta_t_bs < theta_bs_threshold_2) then
        cos_dt    = cos(K_theta*delta_t_bs)
        coef_bell = 1.0_wp - cos_dt * cos_dt
        grad_bell_coef_0 = - ene_attr * K_theta    &
            * sin(2.0_wp * K_theta * delta_t_bs) &
            / (sin_t_bs * r21 * r23 )
        grad_bell_coef_1 = - grad_bell_coef_0 * dot_v21_v23 / r21_sqr
        grad_bell_coef_3 = - grad_bell_coef_0 * dot_v21_v23 / r23_sqr
        g_attr_1(1:3) = grad_bell_coef_0 * v23(1:3) + grad_bell_coef_1 * v21(1:3)
        g_attr_3(1:3) = grad_bell_coef_0 * v21(1:3) + grad_bell_coef_3 * v23(1:3)
        work(1:3) = work(1:3) - g_attr_1(1:3)
        work(4:6) = work(4:6) + g_attr_1(1:3) + g_attr_3(1:3) + coef_bell * grad_attr * v23(1:3)
        work(7:9) = work(7:9) - g_attr_3(1:3) - coef_bell * grad_attr * v23(1:3)
      else
        coef_bell = 0.0_wp
      end if

      ebase_temp = ebase_temp + ene_repl + coef_bell * ene_attr

      ! virial
      !
      viri(1:3) = work(1:3)*coord(i1,1:3) &
                + work(4:6)*coord(i2,1:3) &
                + work(7:9)*coord(i3,1:3)
      virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

      ! store force
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) + work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(4:6)
      force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(7:9)

    end do

    ebase(id+1) = ebase(id+1) + ebase_temp

    !$omp end parallel 
    
    call timer(TimerBaseStack, TimerOff)

    return

  end subroutine compute_energy_DNA_base_stacking

end module cg_energy_bases_mod
