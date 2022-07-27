!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_bonds_mod
!> @brief   calculate bond energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ),
!           Cheng Tan (CT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_bonds_mod

  use at_enefunc_str_mod
  use at_boundary_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_bond
  public :: compute_energy_bond_quartic
  public :: compute_energy_bond_pbc
  public :: compute_energy_bond_quartic_pbc

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond
  !> @brief        calculate bond energy
  !! @authors      YS, TM, JJ
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ebond   : bond energy of target systems
  !!
  !! @note         Fujitsu suggested that force reduction in OpenMP is too slow.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond(enefunc, coord, force, virial, ebond)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ebond

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), r0(:), work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerBond, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_bond
    iend   = enefunc%iend_bond
    list   => enefunc%bond_list
    fc     => enefunc%bond_force_const
    r0     => enefunc%bond_dist_min
    work   => enefunc%work


    ! calculation of bond energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, d12, r12, r_dif, cc_frc, vtmp)    &
    !$omp shared(istart, iend, fc, r0, work, coord, list)    &
    !$omp reduction(+:ebond) reduction(+:virial) 
    !
    do i = istart, iend
    
      ! bond energy: E=K[b-b0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
      r_dif = r12 - r0(i)
      ebond = ebond + fc(i) * r_dif * r_dif

      ! gradient: dE/dX
      !
      cc_frc    = (2.0_wp * fc(i) * r_dif) / r12
      work(1:3,i) = cc_frc * d12(1:3)

      ! virial coefficient
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j, i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j, i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force: F=-dE/dX
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do

    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond_quartic
  !> @brief        calculate quartic bond energy for CG DNA 3SPN.2C model
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ebond   : quartic bond energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond_quartic(enefunc, coord, force, virial, ebond)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ebond

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc, vtmp
    real(wp)                 :: r_dif_2
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), r0(:), work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerBond, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_bond_quartic
    iend   = enefunc%iend_bond_quartic
    list   => enefunc%bond_quartic_list
    fc     => enefunc%bond_quartic_force_const
    r0     => enefunc%bond_quartic_dist_min
    work   => enefunc%work

    ! calculation of quartic bond energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, d12, r12, r_dif, r_dif_2, cc_frc, &
    !$omp         vtmp)                                      &
    !$omp shared(istart, iend, fc, r0, work, coord, list)    &
    !$omp reduction(+:ebond) reduction(+:virial) 
    !
    do i = istart, iend
    
      ! bond_quartic energy: E=K[b-b0]^2 + 100 K[b-b0]^4
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
      r_dif = r12 - r0(i)
      r_dif_2 = r_dif * r_dif
      ebond = ebond + fc(i) * r_dif_2 * (1.0_wp + r_dif_2 * 100.0_wp)

      ! gradient: dE/dX
      !
      cc_frc    = (2.0_wp * fc(i) * r_dif * (1.0_wp + 200.0_wp * r_dif_2)) / r12
      work(1:3,i) = cc_frc * d12(1:3)

      ! virial coefficient
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j, i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j, i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force: F=-dE/dX
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do

    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond_quartic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond_pbc
  !> @brief        calculate bond energy with pbc coordinates
  !! @authors      YS, TM, JJ, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : pbc coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ebond    : bond energy of target systems
  !! @note         Here `coord` should be the PBC coordinates (`coord_pbc`)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond_pbc(enefunc, boundary, coord, force, &
                                     virial, ebond)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ebond

    ! local variables
    real(wp)          :: d12(1:3), r12, r_dif, cc_frc, vtmp
    integer           :: i, j, k, istart, iend
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: fc(:), r0(:), work(:,:)
    integer,  pointer :: list(:,:)


    call timer(TimerBond, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_bond
    iend   = enefunc%iend_bond
    list   => enefunc%bond_list
    fc     => enefunc%bond_force_const
    r0     => enefunc%bond_dist_min
    work   => enefunc%work

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of bond energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, d12, r12, r_dif, cc_frc, vtmp)    &
    !$omp shared(istart, iend, fc, r0, work, coord, list,    &
    !$omp        bsize, half_bsize)                          &
    !$omp reduction(+:ebond) reduction(+:virial) 
    !
    do i = istart, iend
    
      ! bond energy: E=K[b-b0]^2
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

      r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
      r_dif = r12 - r0(i)
      ebond = ebond + fc(i) * r_dif * r_dif

      ! gradient: dE/dX
      !
      cc_frc    = (2.0_wp * fc(i) * r_dif) / r12
      work(1:3,i) = cc_frc * d12(1:3)

      ! virial coefficient
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j, i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j, i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force: F=-dE/dX
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do

    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond_quartic_pbc
  !> @brief        calculate quartic bond energy for CG DNA 3SPN.2C model
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : PBC coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] ebond    : quartic bond energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond_quartic_pbc(enefunc, boundary, coord, &
                                             force, virial, ebond)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: ebond

    ! local variables
    real(wp)          :: d12(1:3), r12, r_dif, cc_frc, vtmp
    real(wp)          :: r_dif_2
    integer           :: i, j, k, istart, iend
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: fc(:), r0(:), work(:,:)
    integer,  pointer :: list(:,:)


    call timer(TimerBond, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_bond_quartic
    iend   = enefunc%iend_bond_quartic
    list   => enefunc%bond_quartic_list
    fc     => enefunc%bond_quartic_force_const
    r0     => enefunc%bond_quartic_dist_min
    work   => enefunc%work

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of quartic bond energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, d12, r12, r_dif, r_dif_2, cc_frc, &
    !$omp         vtmp)                                      &
    !$omp shared(istart, iend, fc, r0, work, coord, list,    &
    !$omp        bsize, half_bsize)                          &
    !$omp reduction(+:ebond) reduction(+:virial) 
    !
    do i = istart, iend
    
      ! bond_quartic energy: E=K[b-b0]^2 + 100 K[b-b0]^4
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

      r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
      r_dif = r12 - r0(i)
      r_dif_2 = r_dif * r_dif
      ebond = ebond + fc(i) * r_dif_2 * (1.0_wp + r_dif_2 * 100.0_wp)

      ! gradient: dE/dX
      !
      cc_frc    = (2.0_wp * fc(i) * r_dif * (1.0_wp + 200.0_wp * r_dif_2)) / r12
      work(1:3,i) = cc_frc * d12(1:3)

      ! virial coefficient
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j, i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j, i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force: F=-dE/dX
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do

    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond_quartic_pbc

end module at_energy_bonds_mod
