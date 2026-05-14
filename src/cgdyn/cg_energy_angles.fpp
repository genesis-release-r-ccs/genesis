!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_angles_mod
!> @brief   Calculates angle energy.
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_angles_mod

  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use table_libs_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_energy_angle
  public  :: compute_energy_angle_g96
  public  :: compute_energy_flexible_angle
  public  :: compute_energy_local_angle

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_angle
  !> @brief        Calculates angle energy.
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle(domain, enefunc, coord, force, &
                                  eangle, eurey, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)
    real(dp),                intent(inout) :: eurey(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc_ub
    real(wp)                 :: eangle_temp, eurey_temp, work(9)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: id

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: fc_ub(:), r0_ub(:)
    integer,         pointer :: anglelist(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerAngle, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    anglelist   => enefunc%angl_list
    fc          => enefunc%angl_force_const
    theta0      => enefunc%angl_theta_min
    fc_ub       => enefunc%urey_force_const
    r0_ub       => enefunc%urey_rmin

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, j, ix, icel1, icel2, icel3, i1, i2, i3, d12, r12_2, &
    !$omp         inv_r12_2, d32, r32_2, inv_r32_2, r12r32, inv_r12r32,      &
    !$omp         cos_t, sin_t, t123, t_dif, cc_frc, cc_frc2, cc_frc3,       &
    !$omp         d13, r13, ub_dif, cc_frc_ub, eangle_temp, eurey_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    eangle_temp = 0.0_wp
    eurey_temp = 0.0_wp

    do i = id+1, enefunc%num_angle_domain, nthread

      ix = enefunc%num_angle_flexible_domain+enefunc%num_angle_local_domain+i

      i1    = id_g2l(anglelist(1,ix))
      i2    = id_g2l(anglelist(2,ix))
      i3    = id_g2l(anglelist(3,ix))

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(i1,1:3) - coord(i2,1:3)
      d32(1:3) = coord(i3,1:3) - coord(i2,1:3)
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt( r12_2*r32_2 )

      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = ( d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3) ) * inv_r12r32
      cos_t  = min(  1.0_wp, cos_t ) 
      cos_t  = max( -1.0_wp, cos_t ) 
      t123   = acos(cos_t)
      t_dif  = t123 - theta0(ix)
      eangle_temp = eangle_temp + fc(ix) * t_dif * t_dif
    
      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max ( EPS, sin_t )
      cc_frc = - ( 2.0_wp * fc(ix) * t_dif ) / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3) = cc_frc * ( d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2 )
      work(4:6) = cc_frc * ( d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3 )

      ! virial
      !
      virial(1,1,id+1) = virial(1,1,id+1) - d12(1)*work(1) - d32(1)*work(4)
      virial(2,2,id+1) = virial(2,2,id+1) - d12(2)*work(2) - d32(2)*work(5)
      virial(3,3,id+1) = virial(3,3,id+1) - d12(3)*work(3) - d32(3)*work(6)

      ! urey-bradley energy: E=K[r-r0]^2
      !
      if (abs(fc_ub(ix)) > EPS) then

        d13(1:3) = coord(i1,1:3) - coord(i3,1:3)
        r13    = sqrt( d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3) )
        ub_dif = r13 - r0_ub(ix)
        eurey_temp = eurey_temp + fc_ub(ix) * ub_dif * ub_dif

        ! gradient: dE/dx
        !
        cc_frc_ub = ( 2.0_wp * fc_ub(ix) * ub_dif ) / r13
        work(7:9) = cc_frc_ub * d13(1:3)

        ! virial
        !
        virial(1,1,id+1) = virial(1,1,id+1) - d13(1)*work(7)
        virial(2,2,id+1) = virial(2,2,id+1) - d13(2)*work(8)
        virial(3,3,id+1) = virial(3,3,id+1) - d13(3)*work(9)

      else

        work(7:9) = 0.0_wp

      end if

      ! store force
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1)-work(1:3)-work(7:9)
      force(i2,1:3,id+1) = force(i2,1:3,id+1)+work(1:3)+work(4:6)
      force(i3,1:3,id+1) = force(i3,1:3,id+1)-work(4:6)+work(7:9)

    end do

    eangle(id+1) = eangle(id+1) + eangle_temp
    eurey(id+1)  = eurey(id+1)  + eurey_temp

    !$omp end parallel 

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_flexible_angle
  !> @brief        Calculates flexible angle energy.
  !! @authors      CK< JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_angle(domain, enefunc, coord, &
                                           force, eangle, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123
    real(wp)                 :: etmp, gradient
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    real(wp)                 :: eangle_temp, work(9)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: atype
    integer                  :: id

    integer,         pointer :: anglelist(:,:), angle_type(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)

    real(wp),        pointer :: theta(:,:)
    real(wp),        pointer :: efunc(:,:)
    real(wp),        pointer :: d2func(:,:)
    real(wp),        pointer :: min_th(:,:)
    real(wp),        pointer :: max_th(:,:)
    real(wp),        pointer :: ener_corr(:)

    call timer(TimerAngle, TimerOn)

    ncell_local     => domain%num_cell_local
    id_g2l          => domain%id_g2l

    anglelist       => enefunc%angl_list
    angle_type      => enefunc%angl_kind

    theta           => enefunc%anglflex_theta
    efunc           => enefunc%anglflex_efunc
    d2func          => enefunc%anglflex_d2func
    min_th          => enefunc%anglflex_min_th
    max_th          => enefunc%anglflex_max_th
    ener_corr       => enefunc%anglflex_ener_corr

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, j, ix, icel1, icel2, icel3, i1, i2, i3, d12, d32,   &
    !$omp         r12_2, r32_2, r12r32, inv_r12r32, inv_r12_2, inv_r32_2,    &
    !$omp         cos_t, sin_t, t123, cc_frc, cc_frc2, cc_frc3, atype,       &
    !$omp         etmp, gradient, eangle_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    eangle_temp = 0.0_wp

    do i = id+1, enefunc%num_angle_flexible_domain, nthread

      i1    = id_g2l(anglelist(1,i))
      i2    = id_g2l(anglelist(2,i))
      i3    = id_g2l(anglelist(3,i))

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(i1,1:3) - coord(i2,1:3)
      d32(1:3) = coord(i3,1:3) - coord(i2,1:3)
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt( r12_2*r32_2 )

      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = (d12(1)*d32(1)+d12(2)*d32(2)+d12(3)*d32(3)) * inv_r12r32
      cos_t  = min( 1.0_wp, cos_t) 
      cos_t  = max(-1.0_wp, cos_t) 
      t123   = acos(cos_t)

      atype  = angle_type(i)

      if (t123 < min_th(1,atype)) then

        etmp = AICG2P_FBA_MIN_ANG_FORCE*(t123-min_th(1,atype)) &
             + min_th(2,atype)
        gradient = AICG2P_FBA_MIN_ANG_FORCE

      else if (t123 > max_th(1,atype)) then

        etmp = AICG2P_FBA_MAX_ANG_FORCE*(t123-max_th(1,atype)) &
             + max_th(2,atype)
        gradient = AICG2P_FBA_MAX_ANG_FORCE

      else

        call table_flexibleangle(atype, t123, theta, efunc, d2func, &
                                 etmp, gradient)

      endif
 
      eangle_temp = eangle_temp + AICG2P_K_ANG*(etmp-ener_corr(atype))
      gradient = AICG2P_K_ANG*gradient
    
      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max(EPS, sin_t)
      cc_frc = -gradient / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3) = cc_frc * (d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2)
      work(4:6) = cc_frc * (d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3)

      ! virial
      !
      virial(1,1,id+1) = virial(1,1,id+1) - d12(1)*work(1) - d32(1)*work(4)
      virial(2,2,id+1) = virial(2,2,id+1) - d12(2)*work(2) - d32(2)*work(5)
      virial(3,3,id+1) = virial(3,3,id+1) - d12(3)*work(3) - d32(3)*work(6)

      ! store force
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) + work(4:6)
      force(i3,1:3,id+1) = force(i3,1:3,id+1) - work(4:6)

    end do

    eangle(id+1) = eangle(id+1) + eangle_temp

    !$omp end parallel 

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_flexible_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_local_angle
  !> @brief        Calculates local angle energy (Gaussian type).
  !! @authors      CK, JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_angle(domain, enefunc, coord, force, &
                                        eangle, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: d13(1:3), r13, r_dif, r_dif2, cc_frc
    real(wp)                 :: eangle_temp, work(9)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: id

    real(wp),        pointer :: fc(:), r0(:), width(:)
    integer,         pointer :: anglelist(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerAngle, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    anglelist       => enefunc%angl_list
    fc              => enefunc%angl_force_const
    r0              => enefunc%angl_theta_min
    width           => enefunc%angl_width

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, j, ix, icel1, icel2, icel3, i1, i2, i3, d13, r13,   &
    !$omp         r_dif, r_dif2, cc_frc, eangle_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    eangle_temp = 0.0_wp

    do i = id+1, enefunc%num_angle_local_domain, nthread

      ix = enefunc%num_angle_flexible_domain + i

      i1    = id_g2l(anglelist(1,ix))
      i2    = id_g2l(anglelist(2,ix))
      i3    = id_g2l(anglelist(3,ix))

      ! angle energy: E=K[t-t0]^2
      !
      d13(1:3) = coord(i1,1:3) - coord(i3,1:3)
      r13      = sqrt(d13(1)*d13(1)+d13(2)*d13(2)+d13(3)*d13(3))
      r_dif    = (r13-r0(ix))/width(ix)
      r_dif2   = -0.5_wp*r_dif*r_dif
      r_dif2   = fc(ix)*exp(r_dif2)
      eangle_temp = eangle_temp + r_dif2
    
      ! gradient: dE/dX
      !
      cc_frc = -r_dif * r_dif2 / (width(ix)*r13)
      work(1:3) = cc_frc * d13(1:3)

      ! virial
      !
      virial(1,1,id+1) = virial(1,1,id+1) - d13(1)*work(1)
      virial(2,2,id+1) = virial(2,2,id+1) - d13(2)*work(2)
      virial(3,3,id+1) = virial(3,3,id+1) - d13(3)*work(3)

      ! store force
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
      force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(1:3)

    end do

    eangle(id+1) = eangle(id+1) + eangle_temp

    !$omp end parallel 

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_local_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_angle_g96
  !> @brief        Calculates angle energy for Martini.
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_g96(domain, enefunc, coord, force, &
                                      eangle, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif, coef
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc_ub
    real(wp)                 :: eangle_temp, work(9)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: id

    real(wp),        pointer :: fc(:), theta0(:)
    integer,         pointer :: anglelist(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerAngle, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    anglelist   => enefunc%angl_list
    fc          => enefunc%angl_force_const
    theta0      => enefunc%angl_theta_min

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, j, ix, icel1, icel2, icel3, i1, i2, i3, d12, r12_2, &
    !$omp         inv_r12_2, d32, r32_2, inv_r32_2, r12r32, inv_r12r32,      &
    !$omp         cos_t, sin_t, t123, t_dif, cc_frc, cc_frc2, cc_frc3,       &
    !$omp         d13, r13, ub_dif, cc_frc_ub, eangle_temp, work, coef)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    eangle_temp = 0.0_wp

    do ix = id+1, enefunc%num_angle_domain, nthread

      i1    = id_g2l(anglelist(1,ix))
      i2    = id_g2l(anglelist(2,ix))
      i3    = id_g2l(anglelist(3,ix))

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(i1,1:3) - coord(i2,1:3)
      d32(1:3) = coord(i3,1:3) - coord(i2,1:3)
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = max(EPS,r12_2*r32_2)
      r12r32   = sqrt(r12r32)

      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = ( d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3) ) * inv_r12r32
      t_dif  = cos_t - cos(theta0(ix))
      coef   = fc(ix)*t_dif
      eangle_temp = eangle_temp + coef*t_dif
    
      ! gradient: dE/dX
      !
      coef = 2.0_wp * coef
      work(1:3) = coef * (d32(1:3)*inv_r12r32 - d12(1:3)*cos_t*inv_r12_2)
      work(4:6) = coef * (d12(1:3)*inv_r12r32 - d32(1:3)*cos_t*inv_r32_2)

      ! virial
      !
      virial(1,1,id+1) = virial(1,1,id+1) - d12(1)*work(1) - d32(1)*work(4)
      virial(2,2,id+1) = virial(2,2,id+1) - d12(2)*work(2) - d32(2)*work(5)
      virial(3,3,id+1) = virial(3,3,id+1) - d12(3)*work(3) - d32(3)*work(6)

      ! store force
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1)-work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1)+work(1:3)+work(4:6)
      force(i3,1:3,id+1) = force(i3,1:3,id+1)-work(4:6)

    end do

    eangle(id+1) = eangle(id+1) + eangle_temp

    !$omp end parallel 

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle_g96

end module cg_energy_angles_mod
