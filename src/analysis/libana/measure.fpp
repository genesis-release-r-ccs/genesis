!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   measure_mod
!> @brief   measure distance, angle and torsion
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module measure_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_com
  public :: compute_dis
  public :: compute_ang
  public :: compute_dih

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      compute_com
  !> @brief        compute center of mass coordinate
  !! @authors      NT
  !! @return       center of mass
  !! @param[in]    coord   : coordinate of points
  !! @param[in]    mass    : mass of points
  !! @param[in]    group   : indices of points group
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function compute_com(coord, mass, group)

    ! return value
    real(wp)                 :: compute_com(3)

    ! formal arguments
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(in)    :: mass (:)
    integer,                 intent(in)    :: group(:)

    ! local variables
    real(wp)                 :: total_mass, w
    integer                  :: i


    if (size(group) == 0) &
      call error_msg('Compute_Com> atom count is ZERO.')

    total_mass = 0.0_wp
    do i = 1, size(group)
      total_mass = total_mass + mass(group(i))
    end do

    compute_com(1:3) = 0.0_wp
    do i = 1, size(group)
      w = mass(group(i)) / total_mass
      compute_com(1:3) = compute_com(1:3) + coord(1:3,group(i)) * w
    end do

    return

  end function compute_com

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      compute_dis
  !> @brief        compute distance of two coordinates
  !! @authors      NT
  !! @return       distance
  !! @param[in]    c1 : coordinate 1
  !! @param[in]    c2 : coordinate 2
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function compute_dis(c1, c2)

    ! return value
    real(wp)                 :: compute_dis

    ! formal arguments
    real(wp)                 :: c1(3)
    real(wp)                 :: c2(3)

    ! local variables
    real(wp)                 :: dis2


    dis2 = (c1(1) - c2(1))**2 + (c1(2) - c2(2))**2 + (c1(3) - c2(3))**2
    compute_dis = sqrt(dis2)

    return

  end function compute_dis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      compute_ang
  !> @brief        compute angle between three coordinates
  !! @authors      NT
  !! @return       angle
  !! @param[in]    c1 : coordinate 1
  !! @param[in]    c2 : coordinate 2
  !! @param[in]    c3 : coordinate 3
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function compute_ang(c1, c2, c3)

    ! return value
    real(wp)                 :: compute_ang

    ! formal arguments
    real(wp)                 :: c1(3)
    real(wp)                 :: c2(3)
    real(wp)                 :: c3(3)

    ! local variables
    real(wp)                 :: v12(3), v32(3)
    real(wp)                 :: d12, d32, cosp, ang


    v12(1:3) = c1(1:3) - c2(1:3)
    v32(1:3) = c3(1:3) - c2(1:3)

    d12 = sqrt(v12(1)*v12(1) + v12(2)*v12(2) + v12(3)*v12(3))
    d32 = sqrt(v32(1)*v32(1) + v32(2)*v32(2) + v32(3)*v32(3))

    if (d12 > 0.0_wp .and. d32 > 0.0_wp) then

      v12(1) = v12(1) / d12
      v12(2) = v12(2) / d12
      v12(3) = v12(3) / d12

      v32(1) = v32(1) / d32
      v32(2) = v32(2) / d32
      v32(3) = v32(3) / d32

      cosp = v12(1)*v32(1) + v12(2)*v32(2) + v12(3)*v32(3)
      if (cosp  >  1.0_wp) then
        cosp = 1.0_wp
      else if (cosp < -1.0_wp) then
        cosp = -1.0_wp
      end if

      ang = acos(cosp)
    end if

    ang = ang * (180.0_wp / PI)

    compute_ang = ang

    return

  end function compute_ang

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      compute_dih
  !> @brief        compute dihedral angle between four coordinates
  !! @authors      NT
  !! @return       dihedral angles
  !! @param[in]    c1 : coordinate 1
  !! @param[in]    c2 : coordinate 2
  !! @param[in]    c3 : coordinate 3
  !! @param[in]    c4 : coordinate 4
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function compute_dih(c1, c2, c3, c4)

    ! return value
    real(wp)                 :: compute_dih

    ! formal arguments
    real(wp)                 :: c1(3)
    real(wp)                 :: c2(3)
    real(wp)                 :: c3(3)
    real(wp)                 :: c4(3)

    ! local variables
    real(wp)                 :: d21(3), d32(3), d43(3), p12(3), p23(3), p123(3)
    real(wp)                 :: r12, r23, r12r23, s1223, s32123, cosp, tor

    
    d21(1:3) = c2(1:3) - c1(1:3)
    d32(1:3) = c3(1:3) - c2(1:3)
    d43(1:3) = c4(1:3) - c3(1:3)

    p12(1) = d21(2)*d32(3) - d32(2)*d21(3)
    p12(2) = d21(3)*d32(1) - d32(3)*d21(1)
    p12(3) = d21(1)*d32(2) - d32(1)*d21(2)
    p23(1) = d32(2)*d43(3) - d43(2)*d32(3)
    p23(2) = d32(3)*d43(1) - d43(3)*d32(1)
    p23(3) = d32(1)*d43(2) - d43(1)*d32(2)

    r12  = p12(1)*p12(1) + p12(2)*p12(2) + p12(3)*p12(3)
    r23  = p23(1)*p23(1) + p23(2)*p23(2) + p23(3)*p23(3)
    r12r23 = sqrt(r12 * r23)

    if (r12r23 > 0.0_wp) then

      ! in the case of bond angle = 0 or pi can not calculate
      ! torsional angle   then do not calculate torsional energy

      r12r23 = 1.0_wp / r12r23
      s1223  = p12(1) * p23(1) + p12(2) * p23(2) + p12(3) * p23(3)
      cosp   = s1223 * r12r23

      if (cosp >  1.0_wp) cosp =  1.0_wp
      if (cosp < -1.0_wp) cosp = -1.0_wp

      p123(1)  = p12(2) * p23(3) - p23(2) * p12(3)
      p123(2)  = p12(3) * p23(1) - p23(3) * p12(1)
      p123(3)  = p12(1) * p23(2) - p23(1) * p12(2) 
      s32123 = d32(1) * p123(1) + d32(2) * p123(2) + d32(3) * p123(3)

      tor = acos(cosp)                                
      if (s32123 < 0.0_wp) then                        
        tor = (2.0_wp * pi) - tor                        
      end if

    end if

    tor = tor * (180.0_wp / PI)
    if (tor .gt. 180.0_wp) tor = tor - 360.0_wp

    compute_dih = tor

    return

  end function compute_dih

end module measure_mod
