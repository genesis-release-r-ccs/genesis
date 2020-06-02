!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_angles_mod
!> @brief   calculate angle energy
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_angles_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_angle
  public  :: compute_energy_angle_g96

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_angle
  !> @brief        calculate angle energy
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle(domain, enefunc, coord, force, eangle, eurey)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)
    real(dp),                intent(inout) :: eurey(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc_ub
    real(wp)                 :: eangle_temp, eurey_temp, work(9)
    integer                  :: i, j, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: list(3)
    integer                  :: id, omp_get_thread_num

    real(wp),        pointer :: fc(:,:), theta0(:,:)
    real(wp),        pointer :: fc_ub(:,:), r0_ub(:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)

    call timer(TimerAngle, TimerOn)

#ifdef PKTIMER
    call timer_sta(212)
#ifdef FJ_PROF_FAPP
    call fapp_start("compute_energy_angle",212,0)
#endif
#endif

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    nangle      => enefunc%num_angle
    anglelist   => enefunc%angle_list
    fc          => enefunc%angle_force_const
    theta0      => enefunc%angle_theta_min
    fc_ub       => enefunc%urey_force_const
    r0_ub       => enefunc%urey_rmin

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                           &
    !$omp private(id, i, j, ix, icel1, icel2, icel3, i1, i2, i3, d12, r12_2, &
    !$omp         inv_r12_2, d32, r32_2, inv_r32_2, r12r32, inv_r12r32,      &
    !$omp         cos_t, sin_t, t123, t_dif, cc_frc, cc_frc2, cc_frc3,       &
    !$omp         d13, r13, ub_dif, cc_frc_ub, eangle_temp, eurey_temp,      &
    !$omp         work, list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      eangle_temp = 0.0_wp
      eurey_temp = 0.0_wp

      do ix = 1, nangle(i)

        list(1) = anglelist(1,ix,i)
        list(2) = anglelist(2,ix,i)
        list(3) = anglelist(3,ix,i)
        icel1   = id_g2l(1,list(1))
        i1      = id_g2l(2,list(1))
        icel2   = id_g2l(1,list(2))
        i2      = id_g2l(2,list(2))
        icel3   = id_g2l(1,list(3))
        i3      = id_g2l(2,list(3))

        ! angle energy: E=K[t-t0]^2
        !
        d12(1) = coord(1,i1,icel1) - coord(1,i2,icel2)
        d12(2) = coord(2,i1,icel1) - coord(2,i2,icel2)
        d12(3) = coord(3,i1,icel1) - coord(3,i2,icel2)
        d32(1) = coord(1,i3,icel3) - coord(1,i2,icel2)
        d32(2) = coord(2,i3,icel3) - coord(2,i2,icel2)
        d32(3) = coord(3,i3,icel3) - coord(3,i2,icel2)
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
        t_dif  = t123 - theta0(ix,i)
        eangle_temp = eangle_temp + fc(ix,i) * t_dif * t_dif
    
        ! gradient: dE/dX
        !
!       sin_t  = sin(t123)
        sin_t  = sqrt(1.0_wp-cos_t*cos_t)
        sin_t  = max ( EPS, sin_t )
        cc_frc = - ( 2.0_wp * fc(ix,i) * t_dif ) / sin_t
        cc_frc2 = cos_t * inv_r12_2
        cc_frc3 = cos_t * inv_r32_2
        work(1) = cc_frc * ( d32(1) * inv_r12r32 - d12(1) * cc_frc2 )
        work(2) = cc_frc * ( d32(2) * inv_r12r32 - d12(2) * cc_frc2 )
        work(3) = cc_frc * ( d32(3) * inv_r12r32 - d12(3) * cc_frc2 )
        work(4) = cc_frc * ( d12(1) * inv_r12r32 - d32(1) * cc_frc3 )
        work(5) = cc_frc * ( d12(2) * inv_r12r32 - d32(2) * cc_frc3 )
        work(6) = cc_frc * ( d12(3) * inv_r12r32 - d32(3) * cc_frc3 )

        ! urey-bradley energy: E=K[r-r0]^2
        !
        if (abs(fc_ub(ix,i)) > EPS) then

          d13(1) = coord(1,i1,icel1) - coord(1,i3,icel3)
          d13(2) = coord(2,i1,icel1) - coord(2,i3,icel3)
          d13(3) = coord(3,i1,icel1) - coord(3,i3,icel3)
          r13    = sqrt( d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3) )
          ub_dif = r13 - r0_ub(ix,i)
          eurey_temp = eurey_temp + fc_ub(ix,i) * ub_dif * ub_dif

          ! gradient: dE/dx
          !
          cc_frc_ub = ( 2.0_wp * fc_ub(ix,i) * ub_dif ) / r13
          work(7) = cc_frc_ub * d13(1)
          work(8) = cc_frc_ub * d13(2)
          work(9) = cc_frc_ub * d13(3)

        else

          work(7) = 0.0_wp
          work(8) = 0.0_wp
          work(9) = 0.0_wp

        end if

        ! store force
        !
        force(1,i1,icel1,id+1) = force(1,i1,icel1,id+1)-work(1)-work(7)
        force(2,i1,icel1,id+1) = force(2,i1,icel1,id+1)-work(2)-work(8)
        force(3,i1,icel1,id+1) = force(3,i1,icel1,id+1)-work(3)-work(9)
        force(1,i2,icel2,id+1) = force(1,i2,icel2,id+1)+work(1)+work(4)
        force(2,i2,icel2,id+1) = force(2,i2,icel2,id+1)+work(2)+work(5)
        force(3,i2,icel2,id+1) = force(3,i2,icel2,id+1)+work(3)+work(6)
        force(1,i3,icel3,id+1) = force(1,i3,icel3,id+1)-work(4)+work(7)
        force(2,i3,icel3,id+1) = force(2,i3,icel3,id+1)-work(5)+work(8)
        force(3,i3,icel3,id+1) = force(3,i3,icel3,id+1)-work(6)+work(9)

      end do

      eangle(id+1) = eangle(id+1) + eangle_temp
      eurey(id+1)  = eurey(id+1)  + eurey_temp

    end do

    !angle energy for water (for the case of no constraint)
    !
    if (enefunc%table%water_bond_calc) then

      do i = id+1, ncell_local, nthread

        eangle_temp = 0.0_wp

        do ix = 1, domain%num_water(i)

          i2 = domain%water_list(1,ix,i)
          i1 = domain%water_list(2,ix,i)
          i3 = domain%water_list(3,ix,i)

          ! angle energy: E=K[t-t0]^2
          !
          d12(1) = coord(1,i1,i) - coord(1,i2,i)
          d12(2) = coord(2,i1,i) - coord(2,i2,i)
          d12(3) = coord(3,i1,i) - coord(3,i2,i)
          d32(1) = coord(1,i3,i) - coord(1,i2,i)
          d32(2) = coord(2,i3,i) - coord(2,i2,i)
          d32(3) = coord(3,i3,i) - coord(3,i2,i)
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
          t_dif  = t123 - enefunc%table%HOH_angle
          eangle_temp = eangle_temp + enefunc%table%HOH_force *t_dif*t_dif

          ! gradient: dE/dX
          !
!         sin_t  = sin(t123)
          sin_t  = sqrt(1.0_wp-cos_t*cos_t)
          sin_t  = max ( EPS, sin_t )
          cc_frc = - (2.0_wp*enefunc%table%HOH_force*t_dif) / sin_t
          cc_frc2 = cos_t * inv_r12_2
          cc_frc3 = cos_t * inv_r32_2
          work(1) = cc_frc * ( d32(1) * inv_r12r32 - d12(1) * cc_frc2 )
          work(2) = cc_frc * ( d32(2) * inv_r12r32 - d12(2) * cc_frc2 )
          work(3) = cc_frc * ( d32(3) * inv_r12r32 - d12(3) * cc_frc2 )
          work(4) = cc_frc * ( d12(1) * inv_r12r32 - d32(1) * cc_frc3 )
          work(5) = cc_frc * ( d12(2) * inv_r12r32 - d32(2) * cc_frc3 )
          work(6) = cc_frc * ( d12(3) * inv_r12r32 - d32(3) * cc_frc3 )
          work(7) = 0.0_wp
          work(8) = 0.0_wp
          work(9) = 0.0_wp

          ! store force
          !
          force(1,i1,i,1) = force(1,i1,i,1)-work(1)
          force(2,i1,i,1) = force(2,i1,i,1)-work(2)
          force(3,i1,i,1) = force(3,i1,i,1)-work(3)
          force(1,i2,i,1) = force(1,i2,i,1)+work(1)+work(4)
          force(2,i2,i,1) = force(2,i2,i,1)+work(2)+work(5)
          force(3,i2,i,1) = force(3,i2,i,1)+work(3)+work(6)
          force(1,i3,i,1) = force(1,i3,i,1)-work(4)
          force(2,i3,i,1) = force(2,i3,i,1)-work(5)
          force(3,i3,i,1) = force(3,i3,i,1)-work(6)

        end do

        eangle(id+1) = eangle(id+1) + eangle_temp

      end do
    end if

    !$omp end parallel 

#ifdef PKTIMER
    call timer_end(212)
#ifdef FJ_PROF_FAPP
    call fapp_stop("compute_energy_angle",212,0)
#endif
#endif

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle


  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    compute_energy_angle_g96
  !> @brief        calculate angle energy
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eangle  : angle energy of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_g96(domain, enefunc, coord, force, eangle)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eangle(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, dif, coef
    real(wp)                 :: eangle_temp, work(6)
    integer                  :: i, ix, icel1, icel2, icel3, i1, i2, i3
    integer                  :: id, omp_get_thread_num

    real(wp),        pointer :: fc(:,:), theta0(:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    call timer(TimerAngle, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    nangle      => enefunc%num_angle
    anglelist   => enefunc%angle_list
    fc          => enefunc%angle_force_const
    theta0      => enefunc%angle_theta_min

    ! calculation of angle energy and gradient
    !
    !$omp parallel default(shared)                                        &
    !$omp private(id, i, ix, icel1, icel2, icel3, i1, i2, i3, d12, d32,   &
    !$omp         r12_2, inv_r12_2, r32_2, inv_r32_2, r12r32, inv_r12r32, &
    !$omp         cos_t, dif, coef, eangle_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread
    
      eangle_temp = 0.0_wp

      do ix = 1, nangle(i)

        icel1 = id_g2l(1,anglelist(1,ix,i))
        i1    = id_g2l(2,anglelist(1,ix,i))
        icel2 = id_g2l(1,anglelist(2,ix,i))
        i2    = id_g2l(2,anglelist(2,ix,i))
        icel3 = id_g2l(1,anglelist(3,ix,i))
        i3    = id_g2l(2,anglelist(3,ix,i))

        ! angle energy: E=K[t-t0]^2
        !
        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i2,icel2)
        d32(1:3) = coord(1:3,i3,icel3) - coord(1:3,i2,icel2)

        r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
        r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
        r12r32   = sqrt( r12_2*r32_2 )

        inv_r12r32 = 1.0_wp / r12r32
        inv_r12_2  = 1.0_wp / r12_2
        inv_r32_2  = 1.0_wp / r32_2

        cos_t  = d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)
        cos_t  = cos_t * inv_r12r32

        dif    = cos_t - cos(theta0(ix, i))
        coef   = fc(ix,i) * dif

        eangle_temp = eangle_temp + ( coef * dif )

        ! gradient: dE/dX
        !
        coef = 2.0_wp * coef
        work(1:3) = coef*(d32(1:3) * inv_r12r32 - d12(1:3) * cos_t * inv_r12_2)
        work(4:6) = coef*(d12(1:3) * inv_r12r32 - d32(1:3) * cos_t * inv_r32_2)

        ! store force
        !
        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1)- work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)+ work(1:3)+work(4:6)
        force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)- work(4:6)

      end do

      eangle(id+1) = eangle(id+1) + eangle_temp

    end do

    !$omp end parallel

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle_g96

end module sp_energy_angles_mod
