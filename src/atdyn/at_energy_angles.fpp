!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_angles_mod
!> @brief   calculate angle energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_angles_mod

  use at_boundary_str_mod
  use at_enefunc_str_mod
  use table_libs_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_angle
  public :: compute_energy_angle_g96
  public :: compute_energy_angle_multi
  public :: compute_energy_flexible_angle
  public :: compute_energy_local_angle
  public :: compute_energy_urey
  public :: compute_energy_angle_pbc
  public :: compute_energy_local_angle_pbc
  public :: compute_energy_flexible_angle_pbc

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle
  !> @brief        calculate angle energy
  !! @authors      YS, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif, vtmp
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: func(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    theta0 => enefunc%angl_theta_min
    work   => enefunc%work
    ! ~CG~ to distinguish normal angle term and AICG2+ local term (functype = 21)
    func   => enefunc%angl_func


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                    &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32, &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123, t_dif,     &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp)                      &
    !$omp shared(istart, iend, fc, theta0, work, list, func, coord)    &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

       ! only calculate for type=1 angle terms??
       ! not sure... only exclude AICG2+ local potential for now....
       if (func(i) == 21) then
          ! work(1:6, i) = 0.0_wp
          cycle
       end if

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt(r12_2*r32_2)
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3) ) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t)
      cos_t  = max(-1.0_wp, cos_t)
      t123   = acos(cos_t)
      t_dif  = t123 - theta0(i)
      eangle = eangle + fc(i) * t_dif * t_dif

      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max(EPS, sin_t)
      cc_frc = - (2.0_wp * fc(i) * t_dif) / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
       
       ! ~CG~ 3SPN.2C DNA: not calculate type 21 angles...
       if (func(i) == 21) then
          cycle
       end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle_multi
  !> @brief        calculate multi angle energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_multi(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif, t_dif1, vtmp
    real(wp)                 :: vtmpa, vtmpb,vtmpmax
    real(wp)                 :: vea, veb, ctmp
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: fc1(:), theta01(:)
    real(wp),        pointer :: epsa(:), cgamma(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    if (enefunc%num_angles <= 0) return

    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list    => enefunc%angl_list
    fc      => enefunc%angl_force_const
    theta0  => enefunc%angl_theta_min
    fc1     => enefunc%angl_force_const1
    theta01 => enefunc%angl_theta_min1
    cgamma  => enefunc%angl_gamma
    epsa    => enefunc%angl_epsa
    work    => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                    &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32, &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123, t_dif,     &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp, t_dif1, vtmpa, vtmpb,&
    !$omp         vtmpmax, vea, veb, ctmp)                             &
    !$omp shared(istart, iend, fc, theta0, work, list, coord,          &
    !$omp        fc1, theta01, epsa, cgamma)                           &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt(r12_2*r32_2)
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t)
      cos_t  = max(-1.0_wp, cos_t)
      t123   = acos(cos_t)
      t_dif  = t123 - theta0(i)
      t_dif1 = t123 - theta01(i)
      vtmpa  = fc(i)  * t_dif  * t_dif + epsa(i)
      vtmpb  = fc1(i) * t_dif1 * t_dif1

      vea   = exp(-cgamma(i)*vtmpa)
      veb   = exp(-cgamma(i)*vtmpb)
      vtmp = vea + veb
      eangle = eangle-(1.0_wp/cgamma(i))*log(vtmp)

      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max(EPS, sin_t)
      ctmp = 2.0_wp * (fc(i) * t_dif * vea/vtmp +  fc1(i) * t_dif1 * veb/vtmp)
      cc_frc = - ctmp / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_angle_multi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle_g96
  !> @brief        calculate angle energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_g96(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, dif, cc_frc, cc_frc2, cc_frc3, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    theta0 => enefunc%angl_theta_min
    work   => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                         &
    !$omp private(i, j, k, d12, d32, r12_2, inv_r12_2, r32_2, inv_r32_2,    &
    !$omp         r12r32, inv_r12r32, cos_t, dif, cc_frc, cc_frc2, cc_frc3, &
    !$omp         vtmp)                                                     &
    !$omp shared(istart, iend, fc, theta0, work, list, coord)               &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))

      r12_2   = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2   = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32  = sqrt(r12_2*r32_2)

      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t   = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32

      dif     = cos_t - cos(theta0(i))
      cc_frc  = fc(i) * dif
      eangle  = eangle + (cc_frc * dif)

      ! gradient: dE/dX
      !
      cc_frc  = 2.0_wp * cc_frc
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle_g96

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_urey
  !> @brief        calculate Urey-Bradley energy
  !! @authors      YS, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_urey(enefunc, coord, force, virial, eurey)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eurey

    ! local variables
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc_ub(:), r0_ub(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_urey
    iend   = enefunc%iend_urey

    list   => enefunc%urey_list
    fc_ub  => enefunc%urey_force_const
    r0_ub  => enefunc%urey_rmin
    work   => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                             &
    !$omp private(i, j, k, d13, r13, ub_dif, cc_frc, vtmp)      &
    !$omp shared(istart, iend, fc_ub, r0_ub, work, coord, list) &
    !$omp reduction(+:eurey) reduction(+:virial)
    !
    do i = istart, iend

      ! urey-bradley energy: E=K[ub-ub0]^2
      !
      d13(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      r13    = sqrt(d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3))
      ub_dif = r13 - r0_ub(i)
      eurey  = eurey + fc_ub(i) * ub_dif * ub_dif

      ! gradient: dE/dx
      !
      cc_frc = (2.0_wp * fc_ub(i) * ub_dif) / r13
      work(1:3,i) = cc_frc * d13(1:3)

      ! virial from urey-bradley
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d13(k)*work(j,i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d13(j)*work(j,i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do


    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do


    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_urey

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flexible_angle
  !> @brief        calculate flexible angle energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_angle(enefunc, coord, force, virial, &
                                           eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif, vtmp, etmp, gradient
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    integer                  :: i, j, k, istart, iend, atype

!    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: work(:,:)
    real(wp),        pointer :: theta(:,:)
    real(wp),        pointer :: efunc(:,:)
    real(wp),        pointer :: d2func(:,:)
    real(wp),        pointer :: min_th(:,:)
    real(wp),        pointer :: max_th(:,:)
    real(wp),        pointer :: ener_corr(:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: at(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angflex
    iend   = enefunc%iend_angflex

    list   => enefunc%anglflex_list
    at     => enefunc%anglflex_type
    theta  => enefunc%anglflex_theta
    efunc  => enefunc%anglflex_efunc
    d2func => enefunc%anglflex_d2func
    min_th => enefunc%anglflex_min_th
    max_th => enefunc%anglflex_max_th
    ener_corr => enefunc%anglflex_ener_corr

    work   => enefunc%work

    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                           &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32,        &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123,                   &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp, etmp, gradient, atype )     &
    !$omp shared(istart, iend,  work, list, coord, at, theta, efunc, d2func,  &
    !$omp        min_th, max_th, ener_corr)                                   &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle tabled energy
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt(r12_2*r32_2)
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t)
      cos_t  = max(-1.0_wp, cos_t)
      t123   = acos(cos_t)

      atype = at(i)
      if (t123 < min_th(1,atype)) then

        etmp =  AICG2P_FBA_MIN_ANG_FORCE*(t123-min_th(1,atype))+min_th(2,atype)
       
        gradient = AICG2P_FBA_MIN_ANG_FORCE
     
      else if (t123 > max_th(1,atype)) then
        etmp =  AICG2P_FBA_MAX_ANG_FORCE*(t123-max_th(1,atype))+max_th(2,atype)
       
        gradient = AICG2P_FBA_MAX_ANG_FORCE
     
      else

        call table_flexibleangle(atype, t123, theta, efunc, d2func,     &
                               etmp, gradient)
      end if

      eangle   = eangle + AICG2P_K_ANG*(etmp-ener_corr(atype))
      gradient = AICG2P_K_ANG*gradient

      ! gradient: dE/dX
      !
      sin_t   = sin(t123)
      sin_t   = max(EPS, sin_t)
      ! d(theta)/dr = -1/sin(theta)*(dcos/dr) 
      cc_frc  = -gradient / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_flexible_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_angle
  !> @brief        calculate angle energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_angle(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d13(1:3), r13, r_dif, cc_frc, vtmp, etmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), r0(:), w(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: func(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    r0     => enefunc%angl_theta_min
    w      => enefunc%angl_w
    work   => enefunc%work
    ! ~CG~ to distinguish normal angle term and AICG2+ local term (functype = 21)
    func   => enefunc%angl_func

    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                 &
    !$omp private(i, j, k, d13, r13, r_dif, cc_frc, vtmp, etmp)     &
    !$omp shared(istart, iend, fc, r0, w, work, coord, func, list)  &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

       ! ~CG~ 3SPN.2C DNA: 
       ! only compute energy and force for functype 21: AICG2+ local angles
       if (func(i) /= 21) then
          ! work(1:3, i) = 0.0_wp
          cycle
       end if

      ! gaussian angle energy
      !
      d13(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(3,i))
      r13      = sqrt(d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3))
      r_dif    = (r13 - r0(i))/w(i)
      etmp     = fc(i)*exp(-0.5_wp*r_dif*r_dif)
      eangle   = eangle + etmp

      ! gradient: dE/dx
      !
      cc_frc   = -r_dif*etmp/w(i)/r13
      work(1:3,i) = cc_frc * d13(1:3)

      ! virial from urey-bradley
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d13(k)*work(j,i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d13(j)*work(j,i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do


    ! store force
    !
    do i = istart, iend
       
       if (func(i) /= 21) then
          cycle
       end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(1:3,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_local_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle_pbc
  !> @brief        calculate angle energy
  !! @authors      YS, TM, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eangle   : ene_angle angle energy of target systems
  !! @note         Here PBC coordinates should be used
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_pbc(enefunc, boundary, coord, force, virial, &
                                      eangle)

    ! formal arguments
    type(s_enefunc), target,  intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eangle

    ! local variables
    real(wp)          :: d12(1:3), r12_2, inv_r12_2
    real(wp)          :: d32(1:3), r32_2, inv_r32_2
    real(wp)          :: r12r32, inv_r12r32
    real(wp)          :: cos_t, sin_t, t123, t_dif, vtmp
    real(wp)          :: cc_frc, cc_frc2, cc_frc3
    integer           :: i, j, k, istart, iend
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: fc(:), theta0(:)
    real(wp), pointer :: work(:,:)
    integer,  pointer :: list(:,:)
    integer,  pointer :: func(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    theta0 => enefunc%angl_theta_min
    work   => enefunc%work
    func   => enefunc%angl_func

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                    &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32, &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123, t_dif,     &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp)                      &
    !$omp shared(istart, iend, fc, theta0, work, list, func, coord,    &
    !$omp        bsize, half_bsize)                                    &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

       if (func(i) == 21) then
          cycle
       end if

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      if (d12(1) > half_bsize(1)) then
        d12(1) = d12(1) - bsize(1)
      else if (d12(1) < -half_bsize(1)) then
        d12(1) = d12(1) + bsize(1)
      end if
      if (d12(2) > half_bsize(2)) then
        d12(2) = d12(2) - bsize(2)
      else if (d12(2) < -half_bsize(2)) then
        d12(2) = d12(2) + bsize(2)
      end if
      if (d12(3) > half_bsize(3)) then
        d12(3) = d12(3) - bsize(3)
      else if (d12(3) < -half_bsize(3)) then
        d12(3) = d12(3) + bsize(3)
      end if

      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      if (d32(1) > half_bsize(1)) then
        d32(1) = d32(1) - bsize(1)
      else if (d32(1) < -half_bsize(1)) then
        d32(1) = d32(1) + bsize(1)
      end if
      if (d32(2) > half_bsize(2)) then
        d32(2) = d32(2) - bsize(2)
      else if (d32(2) < -half_bsize(2)) then
        d32(2) = d32(2) + bsize(2)
      end if
      if (d32(3) > half_bsize(3)) then
        d32(3) = d32(3) - bsize(3)
      else if (d32(3) < -half_bsize(3)) then
        d32(3) = d32(3) + bsize(3)
      end if

      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt(r12_2*r32_2)
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t)
      cos_t  = max(-1.0_wp, cos_t)
      t123   = acos(cos_t)
      t_dif  = t123 - theta0(i)
      eangle = eangle + fc(i) * t_dif * t_dif

      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max(EPS, sin_t)
      cc_frc = - (2.0_wp * fc(i) * t_dif) / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
       
       if (func(i) == 21 ) then
          cycle
       end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_angle_pbc

  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flexible_angle_pbc
  !> @brief        calculate flexible angle energy
  !! @authors      CK, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eangle   : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_angle_pbc(enefunc, boundary, coord, &
                                               force, virial, eangle)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eangle

    ! local variables
    real(wp)          :: d12(1:3), r12_2, inv_r12_2
    real(wp)          :: d32(1:3), r32_2, inv_r32_2
    real(wp)          :: r12r32, inv_r12r32
    real(wp)          :: cos_t, sin_t, t123, t_dif, vtmp, etmp, gradient
    real(wp)          :: cc_frc, cc_frc2, cc_frc3
    integer           :: i, j, k, istart, iend, atype
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: work(:,:)
    real(wp), pointer :: theta(:,:)
    real(wp), pointer :: efunc(:,:)
    real(wp), pointer :: d2func(:,:)
    real(wp), pointer :: min_th(:,:)
    real(wp), pointer :: max_th(:,:)
    real(wp), pointer :: ener_corr(:)
    integer,  pointer :: list(:,:)
    integer,  pointer :: at(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angflex
    iend   = enefunc%iend_angflex

    list   => enefunc%anglflex_list
    at     => enefunc%anglflex_type
    theta  => enefunc%anglflex_theta
    efunc  => enefunc%anglflex_efunc
    d2func => enefunc%anglflex_d2func
    min_th => enefunc%anglflex_min_th
    max_th => enefunc%anglflex_max_th
    ener_corr => enefunc%anglflex_ener_corr

    work   => enefunc%work

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                           &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32,        &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123,                   &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp, etmp, gradient, atype )     &
    !$omp shared(istart, iend,  work, list, coord, at, theta, efunc, d2func,  &
    !$omp        min_th, max_th, ener_corr,                                   &
    !$omp        bsize, half_bsize)                                           &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle tabled energy
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      if (d12(1) > half_bsize(1)) then
        d12(1) = d12(1) - bsize(1)
      else if (d12(1) < -half_bsize(1)) then
        d12(1) = d12(1) + bsize(1)
      end if
      if (d12(2) > half_bsize(2)) then
        d12(2) = d12(2) - bsize(2)
      else if (d12(2) < -half_bsize(2)) then
        d12(2) = d12(2) + bsize(2)
      end if
      if (d12(3) > half_bsize(3)) then
        d12(3) = d12(3) - bsize(3)
      else if (d12(3) < -half_bsize(3)) then
        d12(3) = d12(3) + bsize(3)
      end if

      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      if (d32(1) > half_bsize(1)) then
        d32(1) = d32(1) - bsize(1)
      else if (d32(1) < -half_bsize(1)) then
        d32(1) = d32(1) + bsize(1)
      end if
      if (d32(2) > half_bsize(2)) then
        d32(2) = d32(2) - bsize(2)
      else if (d32(2) < -half_bsize(2)) then
        d32(2) = d32(2) + bsize(2)
      end if
      if (d32(3) > half_bsize(3)) then
        d32(3) = d32(3) - bsize(3)
      else if (d32(3) < -half_bsize(3)) then
        d32(3) = d32(3) + bsize(3)
      end if

      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt(r12_2*r32_2)
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t)
      cos_t  = max(-1.0_wp, cos_t)
      t123   = acos(cos_t)

      atype = at(i)
      if (t123 < min_th(1,atype)) then

        etmp =  AICG2P_FBA_MIN_ANG_FORCE*(t123-min_th(1,atype))+min_th(2,atype)
       
        gradient = AICG2P_FBA_MIN_ANG_FORCE
     
      else if (t123 > max_th(1,atype)) then
        etmp =  AICG2P_FBA_MAX_ANG_FORCE*(t123-max_th(1,atype))+max_th(2,atype)
       
        gradient = AICG2P_FBA_MAX_ANG_FORCE
     
      else

        call table_flexibleangle(atype, t123, theta, efunc, d2func,     &
                               etmp, gradient)
      end if

      eangle   = eangle + AICG2P_K_ANG*(etmp-ener_corr(atype))
      gradient = AICG2P_K_ANG*gradient

      ! gradient: dE/dX
      !
      sin_t   = sin(t123)
      sin_t   = max(EPS, sin_t)
      ! d(theta)/dr = -1/sin(theta)*(dcos/dr) 
      cc_frc  = -gradient / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6,i) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_flexible_angle_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_angle_pbc
  !> @brief        calculate angle energy
  !! @authors      CK, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eangle   : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_angle_pbc(enefunc, boundary, coord, force, &
                                            virial, eangle)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eangle

    ! local variables
    integer           :: i, j, k, istart, iend
    real(wp)          :: d13(1:3), r13, r_dif, cc_frc, vtmp, etmp
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: fc(:), r0(:), w(:)
    real(wp), pointer :: work(:,:)
    integer,  pointer :: list(:,:)
    integer,  pointer :: func(:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    r0     => enefunc%angl_theta_min
    w      => enefunc%angl_w
    work   => enefunc%work
    func   => enefunc%angl_func

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                 &
    !$omp private(i, j, k, d13, r13, r_dif, cc_frc, vtmp, etmp)     &
    !$omp shared(istart, iend, fc, r0, w, work, coord, func, list,  &
    !$omp        bsize, half_bsize)                                 &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

       ! ~CG~ 3SPN.2C DNA: 
       ! only compute energy and force for functype 21: AICG2+ local angles
       if (func(i) /= 21 ) then
          ! work(1:3, i) = 0.0_wp
          cycle
       end if

      ! gaussian angle energy
      !
      d13(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(3,i))
      if (d13(1) > half_bsize(1)) then
        d13(1) = d13(1) - bsize(1)
      else if (d13(1) < -half_bsize(1)) then
        d13(1) = d13(1) + bsize(1)
      end if
      if (d13(2) > half_bsize(2)) then
        d13(2) = d13(2) - bsize(2)
      else if (d13(2) < -half_bsize(2)) then
        d13(2) = d13(2) + bsize(2)
      end if
      if (d13(3) > half_bsize(3)) then
        d13(3) = d13(3) - bsize(3)
      else if (d13(3) < -half_bsize(3)) then
        d13(3) = d13(3) + bsize(3)
      end if

      r13      = sqrt(d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3))
      r_dif    = (r13 - r0(i))/w(i)
      etmp     = fc(i)*exp(-0.5_wp*r_dif*r_dif)
      eangle   = eangle + etmp

      ! gradient: dE/dx
      !
      cc_frc   = -r_dif*etmp/w(i)/r13
      work(1:3,i) = cc_frc * d13(1:3)

      ! virial from urey-bradley
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d13(k)*work(j,i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d13(j)*work(j,i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do


    ! store force
    !
    do i = istart, iend
       
       if (func(i) /= 21 ) then
          cycle
       end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(1:3,i)
    end do


    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_local_angle_pbc

end module at_energy_angles_mod
