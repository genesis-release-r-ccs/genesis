!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_dihedrals_mod
!> @brief   calculate dihedral energy
!! @authors Chigusa Kobayashi(CK), Takao Yoda (TY), Cheng Tan(CT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_dihedrals_mod

  use at_boundary_str_mod
  use at_enefunc_str_mod
  use dihedral_libs_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutine
  public :: compute_energy_dihed
  public :: compute_energy_rb_dihed
  public :: compute_energy_improp
  public :: compute_energy_improp_cos
  public :: compute_energy_cmap
  public :: compute_energy_flexible_dihed
  public :: compute_energy_local_dihed
  public :: compute_energy_cg_dihed_periodic_cos2mod
  public :: compute_energy_cg_dihed_periodic_sin3mod
  public :: compute_energy_cg_dihed_gaussian_cos2mod
  public :: compute_energy_cg_dihed_gaussian_sin3mod
  public :: compute_energy_cg_dihed_flexible_cos2mod
  ! PBC
  public :: compute_energy_dihed_pbc
  public :: compute_energy_local_dihed_pbc
  public :: compute_energy_flexible_dihed_pbc
  public :: compute_energy_cg_dihed_periodic_cos2mod_pbc
  public :: compute_energy_cg_dihed_gaussian_cos2mod_pbc
  public :: compute_energy_cg_dihed_flexible_cos2mod_pbc

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed
  !> @brief        calculate dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, id, krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  = enefunc%istart_dihedral
    iend    = enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    nperiod => enefunc%dihe_periodicity
    phase   => enefunc%dihe_phase
    work    => enefunc%work
!   edihe = 0.0_wp
    ! ~CG~ 3SPN.2C DNA: dihedral type 1
    func    => enefunc%dihe_func

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, func, coord)    &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      ! ~CG~ 3SPN.2C DNA: only calculate functype /= 1 dihedrals...
      if (func(i) == 21 .or. func(i) == 32 .or. func(i) == 33 .or. &
          func(i) == 41 .or. func(i) == 43) then
        ! work(1:9, i) = 0.0_wp
        cycle
      end if

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot = 0
      do while (krot < nperiod(i))
        tmp   = cosnt*cos_dih - sinnt*sin_dih
        sinnt = sinnt*cos_dih + cosnt*sin_dih
        cosnt = tmp
        krot = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      edihe = edihe + fc(i) * (1.0_wp + cospha*cosnt + sinnt*sinpha)

      grad_coef = fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      ! ~CG~ 3SPN.2C DNA: only calculate functype /= 1 dihedrals...
      if (func(i) == 21 .or. func(i) == 32 .or. func(i) == 33 .or. &
          func(i) == 41 .or. func(i) == 43) then
        cycle
      end if
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_rb_dihed
  !> @brief        calculate Ryckaert-Bellemans dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_rb_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, istart, iend, icn
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, coef
    real(wp)                 :: cos_dih, sin_dih, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:,:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerDihedral, TimerOn)

    list    => enefunc%rb_dihe_list
    fc      => enefunc%rb_dihe_c
    work    => enefunc%work

    istart  = enefunc%istart_rb_dihed
    iend    = enefunc%iend_rb_dihed

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k,  aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         icn,  vtmp, coef)               &
    !$omp shared(istart, iend, list, fc, work, coord)          &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

!     psi = phi - pi
      cos_dih = -cos_dih
      sin_dih = -sin_dih

      coef = 0.0_wp
      do icn = 1, 6
        coef = coef + real(icn-1,wp) * fc(icn,i) * cos_dih**(icn-2)
        edihe = edihe + fc(icn,i) * cos_dih**(icn-1)
      end do

      grad_coef = sin_dih * coef
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial 
      !
      do j = 1, 3
        do k = j + 1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_rb_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp
  !> @brief        calculate improper dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eimpr   : improper dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp(enefunc, coord, force, virial, eimpr)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eimpr

    ! local variables
    integer                  :: i, j, k
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi, vtmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosdif, sindif
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)


    call timer(TimerDihedral, TimerOn)

    istart = enefunc%istart_improper
    iend   = enefunc%iend_improper

    list   => enefunc%impr_list
    fc     => enefunc%impr_force_const
    phase  => enefunc%impr_phase
    work   => enefunc%work

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k,  aindex, cos_dih, sin_dih, grad, v, grad_coef,      &
    !$omp         cospha, sinpha, diffphi, cosdif, sindif, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, coord)          &
    !$omp reduction(+:eimpr) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      cosdif = cos_dih*cospha + sin_dih*sinpha
      ! In this case, sindif = -sindif (2019/07/03) 
      sindif = cos_dih*sinpha - sin_dih*cospha

      if (cosdif > 1.0E-1_wp) then
        diffphi = asin(sindif)
      else
        diffphi = sign(1.0_wp,sindif)*acos(cosdif)
      end if

      eimpr     = eimpr + fc(i)*diffphi*diffphi

      grad_coef = 2.0_wp*fc(i)*diffphi
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      i1 = list(1,i)
      i2 = list(2,i)
      i3 = list(3,i)
      i4 = list(4,i)
      force(1,i1) = force(1,i1) - work(1,i)
      force(2,i1) = force(2,i1) - work(2,i)
      force(3,i1) = force(3,i1) - work(3,i)
      force(1,i2) = force(1,i2) + work(1,i) - work(4,i)
      force(2,i2) = force(2,i2) + work(2,i) - work(5,i)
      force(3,i2) = force(3,i2) + work(3,i) - work(6,i)
      force(1,i3) = force(1,i3) + work(4,i) + work(7,i)
      force(2,i3) = force(2,i3) + work(5,i) + work(8,i)
      force(3,i3) = force(3,i3) + work(6,i) + work(9,i)
      force(1,i4) = force(1,i4) - work(7,i)
      force(2,i4) = force(2,i4) - work(8,i)
      force(3,i4) = force(3,i4) - work(9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_improp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp_cos
  !> @brief        calculate improper dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eimpr   : improper dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp_cos(enefunc, coord, force, virial, eimpr)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eimpr

    ! local variables
    integer                  :: i, j, k,  krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi, tmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosdif, sindif
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)


    call timer(TimerDihedral, TimerOn)

    istart = enefunc%istart_improper
    iend   = enefunc%iend_improper

    list    => enefunc%impr_list
    fc      => enefunc%impr_force_const
    nperiod => enefunc%impr_periodicity
    phase   => enefunc%impr_phase
    work    => enefunc%work

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, coord)          &
    !$omp reduction(+:eimpr) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot = 0
      do while (krot < nperiod(i))
        tmp   = cosnt*cos_dih - sinnt*sin_dih
        sinnt = sinnt*cos_dih + cosnt*sin_dih
        cosnt = tmp
        krot = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      eimpr = eimpr + fc(i) * (1.0_wp + cospha*cosnt + sinnt*sinpha)

      grad_coef = fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      i1 = list(1,i)
      i2 = list(2,i)
      i3 = list(3,i)
      i4 = list(4,i)
      force(1,i1) = force(1,i1) - work(1,i)
      force(2,i1) = force(2,i1) - work(2,i)
      force(3,i1) = force(3,i1) - work(3,i)
      force(1,i2) = force(1,i2) + work(1,i) - work(4,i)
      force(2,i2) = force(2,i2) + work(2,i) - work(5,i)
      force(3,i2) = force(3,i2) + work(3,i) - work(6,i)
      force(1,i3) = force(1,i3) + work(4,i) + work(7,i)
      force(2,i3) = force(2,i3) + work(5,i) + work(8,i)
      force(3,i3) = force(3,i3) + work(6,i) + work(9,i)
      force(1,i4) = force(1,i4) - work(7,i)
      force(2,i4) = force(2,i4) - work(8,i)
      force(3,i4) = force(3,i4) - work(9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_improp_cos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cmap
  !> @brief        calculate cmap energy
  !! @authors      CK, TY
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ecmap   : cmap energy of target systems
  !! @note         A.D.MacKerell et al., J.Comput.Chem., 25, 1400-1415 (2004).
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_cmap(enefunc, coord, force, virial, ecmap)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ecmap

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, dihed
    real(wp)                 :: delta, inv_delta, ctmp, vtmp
    integer                  :: aindex(1:8)
    integer                  :: ncmap, istart, iend
    integer                  :: icmap
    integer                  :: igrid(1:2), ngrid0, i, j, m, igr, ip, ictype
    real(wp)                 :: work(1:9,1:2), grad_coef(1:2), gradtot(1:2)
    real(wp)                 :: dgrid(1:2)
    real(wp)                 :: gradphi(1:9), vphi(1:3,1:3)
    real(wp)                 :: gradpsi(1:9), vpsi(1:3,1:3)
    real(wp)                 :: grid_power(1:2,1:4), dgrid_power(1:2,1:4)

    real(wp),        pointer :: coef(:,:,:,:,:)
    real(wp),        pointer :: f(:,:,:)
    integer,         pointer :: resol(:), list(:,:), ctype(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_cmap
    iend   = enefunc%iend_cmap
    ncmap  =  enefunc%num_cmaps
    list   => enefunc%cmap_list
    coef   => enefunc%cmap_coef
    resol  => enefunc%cmap_resolution
    ctype  => enefunc%cmap_type
    f      => enefunc%cmap_force

    ! calculate cmap energy and gradient
    !
    !$omp parallel do default(none)                                       &
    !$omp private(icmap, i, j, aindex, ngrid0, delta, inv_delta, ip,      &
    !$omp         cos_dih, sin_dih, dihed, igr, dgrid, igrid, vtmp,       &
    !$omp         grid_power, dgrid_power, work, m,  ctmp, grad_coef,     &
    !$omp         gradpsi, gradphi, vpsi, vphi, gradtot, ictype)          &
    !$omp shared(istart, iend, list, coord, coef, resol,  f, ctype)       &
    !$omp reduction(+:ecmap) reduction(+:virial)
    do icmap = istart, iend

      ictype    = ctype(icmap)
      ngrid0    = resol(ictype)
      delta     = 360.0_wp/real(ngrid0,wp)
      inv_delta = 1.0_wp/delta

      ! calc dih1 and gradient
      !
      aindex(1:8) = list(1:8,icmap)

      ! dih1 is backbone phi
      call calculate_dihedral(aindex(1:4), coord, cos_dih, sin_dih,          &
                              gradphi, vphi)

      gradphi(1:9) = gradphi(1:9)/RAD

      if (abs(cos_dih) > 1.0E-1_wp) then
        dihed = asin(sin_dih)/RAD
        if (cos_dih < 0.0_wp) then
          if (dihed > 0.0_wp) then
            dihed=180.0_wp-dihed
          else
            dihed=-180.0_wp-dihed
          end if
        end if
      else
        dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
      end if

      if (dihed < -180.0_wp) then
        dihed = dihed + 360.0_wp
      else if (dihed > 180.0_wp) then
        dihed = dihed - 360.0_wp
      end if

      igr      = int((dihed+180.0_wp)*inv_delta)
      dgrid(1) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
      igrid(1) = igr + 1

      !dbg
      if (igrid(1) < 1 .or. igrid(1) > 24) then
        write(6,'("debug-charmm-6-2")')
        write(6,'(2e20.4)') dihed, inv_delta
        write(6,'(i8)') igrid(1)
      end if
      !dbg

      ! dih2 is backbone psi
      call calculate_dihedral(aindex(5:8), coord, cos_dih, sin_dih,          &
                              gradpsi, vpsi)

      gradpsi(1:9) = gradpsi(1:9)/RAD

      if (abs(cos_dih) > 1.0E-1_wp) then
        dihed = asin(sin_dih)/RAD
        if (cos_dih < 0.0_wp) then
          if (dihed > 0.0_wp) then
            dihed=180.0_wp-dihed
          else
            dihed=-180.0_wp-dihed
          end if
        end if
      else
        dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
      end if

      if (dihed < -180.0_wp) then
        dihed = dihed + 360.0_wp
      else if (dihed > 180.0_wp) then
        dihed = dihed - 360.0_wp
      end if

      igr      = int((dihed+180.0_wp)*inv_delta)
      dgrid(2) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
      igrid(2) = igr + 1

      !dbg
      if (igrid(2) < 1 .or. igrid(2) > 24) then
        write(6,'("debug-charmm-6-3")')
        write(6,'(2e20.4)') dihed, inv_delta
        write(6,'(i8)') igrid(2)
      end if
      !dbg

      grid_power(1:2,1) = 1.0_wp

      ip = 1
      do while(ip < 4)
        ip = ip + 1
        grid_power(1:2,ip) = grid_power(1:2,ip-1)*dgrid(1:2)
      end do

      dgrid_power(1:2,1) = 0.0_wp
      dgrid_power(1:2,2) = inv_delta

      ip = 2
      do while(ip < 4)
        ip = ip + 1
        dgrid_power(1:2,ip) = grid_power(1:2,ip-1)*real(ip-1,wp)*inv_delta
      end do

      ! calculate Ecmap and gradient
      !
      work(1:9,1:2) = 0.0_wp
      gradtot(1:2) = 0.0_wp

      do j = 1, 4
        do i = 1, 4
          ctmp = coef(i,j,igrid(2),igrid(1),ictype)

          ! cmap energy
          !
          ecmap = ecmap + grid_power(2,i)*grid_power(1,j)*ctmp

          ! for gradient
          !
          grad_coef(1) = -dgrid_power(1,j)*grid_power(2,i)*ctmp
          grad_coef(2) = -dgrid_power(2,i)*grid_power(1,j)*ctmp

          gradtot(1:2) = gradtot(1:2) + grad_coef(1:2)
        end do
      end do
      work(1:9,1) = gradtot(1)*gradphi(1:9)
      work(1:9,2) = gradtot(2)*gradpsi(1:9)

      f(1:3,1,icmap) = -work(1:3,1)
      f(1:3,2,icmap) =  work(1:3,1) - work(4:6,1)
      f(1:3,3,icmap) =  work(4:6,1) + work(7:9,1)
      f(1:3,4,icmap) = -work(7:9,1)

      f(1:3,5,icmap) = -work(1:3,2)
      f(1:3,6,icmap) =  work(1:3,2) - work(4:6,2)
      f(1:3,7,icmap) =  work(4:6,2) + work(7:9,2)
      f(1:3,8,icmap) = -work(7:9,2)

      do i  = 1, 3
        do j  = i+1, 3
          vtmp =  (gradtot(1)*vphi(j,i) + gradtot(2)*vpsi(j,i))/RAD
          virial(j,i) =  virial(j,i)  + vtmp
          virial(i,j) =  virial(i,j)  + vtmp
        end do
        vtmp =  (gradtot(1)*vphi(i,i) + gradtot(2)*vpsi(i,i))/RAD
        virial(i,i) =  virial(i,i)  + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do icmap = istart, iend
      aindex(1:8) = list(1:8, icmap)

      force(1:3, aindex(1)) = force(1:3, aindex(1)) + f(1:3,1,icmap)
      force(1:3, aindex(2)) = force(1:3, aindex(2)) + f(1:3,2,icmap)
      force(1:3, aindex(3)) = force(1:3, aindex(3)) + f(1:3,3,icmap)
      force(1:3, aindex(4)) = force(1:3, aindex(4)) + f(1:3,4,icmap)
      force(1:3, aindex(5)) = force(1:3, aindex(5)) + f(1:3,5,icmap)
      force(1:3, aindex(6)) = force(1:3, aindex(6)) + f(1:3,6,icmap)
      force(1:3, aindex(7)) = force(1:3, aindex(7)) + f(1:3,7,icmap)
      force(1:3, aindex(8)) = force(1:3, aindex(8)) + f(1:3,8,icmap)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flexible_dihed
  !> @brief        calculate flexible dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, id, krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4), dtype
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, etmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cos_2dih, sin_2dih
    real(wp)                 :: cos_3dih, sin_3dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), c(7)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: coef(:,:)
    real(wp),        pointer :: ener_corr(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:), types(:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  = enefunc%istart_dihedflex
    iend    = enefunc%iend_dihedflex
    list    => enefunc%diheflex_list
    coef    => enefunc%diheflex_coef
    types   => enefunc%diheflex_type
    func    => enefunc%diheflex_func
!   c0      => enefunc%diheflex_c0
!   c1      => enefunc%diheflex_c1
!   c2      => enefunc%diheflex_c2
!   c3      => enefunc%diheflex_c3
!   c4      => enefunc%diheflex_c4
!   c5      => enefunc%diheflex_c5
!   c6      => enefunc%diheflex_c6
    work    => enefunc%work
    ener_corr => enefunc%diheflex_ener_corr

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, etmp,         &
    !$omp         cos_2dih, cos_3dih, sin_2dih, sin_3dih, dtype, c)            &
    !$omp shared(istart, iend, list, coef, types, func, work, coord, ener_corr)  &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      if (func(i) /= 22) then
        cycle
      end if

      aindex(1:4) = list(1:4,i)
      dtype = types(i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      !
      ! cos(2t)=2*cos(t)^2-1
      ! sin(2t)=2*cos(t)*sin(t)
      ! cos(3t)= 4*cos(t)^3-3*cos(t)
      ! sin(3t)=-4*sin(t)^3+3*sin(t)
      !
      cos_2dih =  2.0_wp*cos_dih*cos_dih-1.0_wp
      sin_2dih =  2.0_wp*cos_dih*sin_dih
      cos_3dih =  4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih
      sin_3dih = -4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih
      c(1:7) = coef(1:7,dtype)
      etmp =  c(1)                             &
            + c(2)*cos_dih  + c(3)*sin_dih     &
            + c(4)*cos_2dih + c(5)*sin_2dih    &
            + c(6)*cos_3dih + c(7)*sin_3dih    
     edihe = edihe + AICG2P_K_DIHE*(etmp-ener_corr(dtype))

     grad_coef = - c(2)*sin_dih                         & 
                 + c(3)*cos_dih                         &
                 - 2.0_wp*c(4)*sin_2dih &
                 + 2.0_wp*c(5)*cos_2dih &
                 - 3.0_wp*c(6)*sin_3dih &
                 + 3.0_wp*c(7)*cos_3dih 
      grad_coef = -AICG2P_K_DIHE*grad_coef

      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend

      if (func(i) /= 22) then
        cycle
      end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)

    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_flexible_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_dihed
  !> @brief        calculate local dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, id, krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: grad_coef, tmp
    real(wp)                 :: cos_dih, sin_dih, vtmp, dihe, etmp
    real(wp)                 :: cospha, sinpha
    real(wp)                 :: cosdif, sindif, diffphi
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), theta_min(:), w(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    ! ~CG~ 3SPN.2C DNA: for AICG2+ local and 3PSN DNA dih SPSP/PSPS
    integer,         pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart    = enefunc%istart_dihedral
    iend      = enefunc%iend_dihedral
    list      => enefunc%dihe_list
    fc        => enefunc%dihe_force_const
    theta_min => enefunc%dihe_theta_min
    w         => enefunc%dihe_w
    work      => enefunc%work
    ! ~CG~ 3SPN.2C DNA: dihedral functype
    func      => enefunc%dihe_func

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         vtmp,  dihe, etmp, cospha, sinpha, cosdif, sindif, diffphi)  &
    !$omp shared(istart, iend, list, fc, theta_min, w, work, func, coord)      &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

       if (func(i) /= 21) then
          cycle
       end if

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cospha = cos(theta_min(i))
      sinpha = sin(theta_min(i))

      cosdif = cos_dih*cospha + sin_dih*sinpha
      ! correct definintion (2019/03/27)
      sindif = sin_dih*cospha - cos_dih*sinpha

      if (cosdif > 1.0E-1_wp) then
        diffphi = asin(sindif)
      else
        diffphi = sign(1.0_wp,sindif)*acos(cosdif)
      end if

      dihe = diffphi/w(i)

      etmp = fc(i)*exp(-0.5_wp*dihe*dihe)
!      write(6,*) i, etmp, fc(i), dihe*dihe, exp(-0.5_wp*dihe*dihe)
      edihe = edihe + etmp
! from equation
!      grad_coef = -1*-1*(dihe*etmp)/w(i)
      grad_coef = (dihe*etmp)/w(i)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
       ! ~CG~ 3SPN.2C DNA: only calculate functype == 21 dihedrals...
       if (func(i) /= 21) then
          cycle
       end if
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_local_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_periodic_cos2mod
  !> @brief        calculate dihedral energy in a safe way
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_periodic_cos2mod(enefunc, coord, force, &
                                                      virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, l, m, n, id, omp_get_thread_num
    integer                  :: krot
    real(wp)                 :: tmp, etmp, var_tmp
    real(wp)                 :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)                 :: rji2, rkj2, rkl2
    real(wp)                 :: rji, rkj, rkl
    real(wp)                 :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)                 :: inv_rji, inv_rkj, inv_rkl
    real(wp)                 :: dotpro_ijk, dotpro_jkl
    real(wp)                 :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2
    real(wp)                 :: inv_raijkl
    real(wp)                 :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)                 :: v1_sqr, v2_sqr
    real(wp)                 :: u_t_1
    real(wp)                 :: u_t_2
    real(wp)                 :: A1, A2, A3, A4
    real(wp)                 :: B1
    real(wp)                 :: C1
    real(wp)                 :: cos_dih, sin_dih, u_dih
    real(wp)                 :: cospha, sinpha, cosnt, sinnt
    real(wp)                 :: coef_dih_shift
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)
    real(wp)                 :: ftmp(1:3,1:4)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: fc(:)
    real(wp),        pointer :: phase(:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    phase   => enefunc%dihe_phase
    nperiod => enefunc%dihe_periodicity
    func    => enefunc%dihe_func
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, krot, id,               &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih,            &
    !$omp         cospha, sinpha, cosnt, sinnt,       &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        fc, phase, nperiod, func, coord,     &
    !$omp        force, nthread)                      &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread
      if (func(i) /= 32) then
        cycle
      end if


      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      !
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      !
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      !
      raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

      ! -----------------------------------------
      ! calculate cos^2 theta 1 and cos^2 theta 2
      ! -----------------------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      !
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      !
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      !
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot  = 0

      do while (krot < nperiod(i))
        tmp   = cosnt * cos_dih - sinnt * sin_dih
        sinnt = sinnt * cos_dih + cosnt * sin_dih
        cosnt = tmp
        krot  = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      u_dih  = fc(i) * (coef_dih_shift + 1.0_wp + cospha*cosnt + sinnt*sinpha)

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_periodic_cos2mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_periodic_sin3mod
  !> @brief        calculate dihedral energy in a safe way (sin^3)
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         M Bulacu et al., JCTC, 2013, 9, 3282-3292.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_periodic_sin3mod(enefunc, coord, force, &
                                                      virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, l, m, n, id, omp_get_thread_num
    integer                  :: krot
    real(wp)                 :: tmp, etmp, var_tmp
    real(wp)                 :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)                 :: rji2, rkj2, rkl2
    real(wp)                 :: rji, rkj, rkl
    real(wp)                 :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)                 :: inv_rji, inv_rkj, inv_rkl
    real(wp)                 :: dotpro_ijk, dotpro_jkl
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2
    real(wp)                 :: inv_raijkl
    real(wp)                 :: cos_t_1, sin_t_1, t_1, u_t_1
    real(wp)                 :: cos_t_2, sin_t_2, t_2, u_t_2
    real(wp)                 :: cos_2_t, cos_2_t_m_1
    real(wp)                 :: cos_2_t_p_1, cos_2_t_p_1_sqr
    real(wp)                 :: A1, A2, A3, A4
    real(wp)                 :: B1, B2, B3, B4
    real(wp)                 :: C1, C2, C3, C4
    real(wp)                 :: cos_dih, sin_dih, u_dih
    real(wp)                 :: cospha, sinpha, cosnt, sinnt
    real(wp)                 :: coef_dih_shift
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)
    real(wp)                 :: ftmp(1:3,1:4)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: fc(:)
    real(wp),        pointer :: phase(:)
    integer,         pointer :: nperiod(:)
    integer,         pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = sin^3(t)
    ! 
    ! M'(t) = 3 sin^2(t) cos(t)
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    phase   => enefunc%dihe_phase
    nperiod => enefunc%dihe_periodicity
    func    => enefunc%dihe_func
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, krot, id,               &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         cos_t_1, sin_t_1, t_1, u_t_1,       &
    !$omp         cos_t_2, sin_t_2, t_2, u_t_2,       &
    !$omp         cos_2_t, cos_2_t_m_1,               &
    !$omp         cos_2_t_p_1, cos_2_t_p_1_sqr,       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, B2, B3, B4,                     &
    !$omp         C1, C2, C3, C4,                     &
    !$omp         cos_dih, sin_dih, u_dih,            &
    !$omp         cospha, sinpha, cosnt, sinnt,       &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        fc, phase, nperiod, func, coord,     &
    !$omp        force, nthread)                      &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread
      if (func(i) /= 33) then
        cycle
      end if

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      ! 
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      ! 
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      ! 
      raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      cos_t_1 = dotpro_ijk * inv_rji * inv_rkj
      cos_t_1 = min( 1.0_wp, cos_t_1)
      cos_t_1 = max(-1.0_wp, cos_t_1)
      sin_t_1 = sqrt(1.0_wp - cos_t_1 * cos_t_1)
      t_1     = acos(cos_t_1)
      !
      cos_t_2 = dotpro_jkl * inv_rkj * inv_rkl
      cos_t_2 = min( 1.0_wp, cos_t_2)
      cos_t_2 = max(-1.0_wp, cos_t_2)
      sin_t_2 = sqrt(1.0_wp - cos_t_2 * cos_t_2)
      t_2     = acos(cos_t_2)

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      u_t_1 = sin_t_1 * sin_t_1 * sin_t_1
      A1    = sin_t_1 * inv_rkj * inv_rji2
      A2    = sin_t_1 * inv_rji * inv_rkj2 * (-cos_t_1)
      B1    = 3.0_wp * cos_t_1 * sin_t_1

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      u_t_2 = sin_t_2 * sin_t_2 * sin_t_2
      A3    = sin_t_2 * inv_rkl * inv_rkj2 * cos_t_2
      A4    = sin_t_2 * inv_rkj * inv_rkl2
      C1    = 3.0_wp * cos_t_2 * sin_t_2

      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot  = 0
        
      do while (krot < nperiod(i))
        tmp   = cosnt * cos_dih - sinnt * sin_dih
        sinnt = sinnt * cos_dih + cosnt * sin_dih
        cosnt = tmp
        krot  = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      u_dih  = fc(i) * (coef_dih_shift + cospha*cosnt + sinnt*sinpha)

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      tmp = u_dih * u_t_2 * B1
      B2  = inv_rji * inv_rkj
      B3  = inv_rji2
      B4  = inv_rkj2
      Q(1:3, 1) = tmp * (- B2 * dkj(1:3) - B3 * cos_t_1 * dji(1:3))
      Q(1:3, 3) = tmp * (+ B2 * dji(1:3) + B4 * cos_t_1 * dkj(1:3))
      Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
      Q(1:3, 4) = 0.0_wp

      ! ------
      ! part 3
      ! ------
      tmp = u_dih * u_t_1 * C1
      C2  = inv_rkj * inv_rkl
      C3  = inv_rkj2
      C4  = inv_rkl2
      R(1:3, 2) = tmp * (C2 * dkl(1:3) - C3 * cos_t_2 * dkj(1:3))
      R(1:3, 4) = tmp * (C2 * dkj(1:3) - C4 * cos_t_2 * dkl(1:3))
      R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
      R(1:3, 1) = 0.0_wp

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_periodic_sin3mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_gaussian_cos2mod
  !> @brief        calculate Gaussian-type dihedral energy in a safe way
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_gaussian_cos2mod(enefunc, coord, force, &
                                                      virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, l, m, n, id, omp_get_thread_num
    real(wp)                 :: tmp, etmp, var_tmp
    real(wp)                 :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)                 :: rji2, rkj2, rkl2
    real(wp)                 :: rji, rkj, rkl
    real(wp)                 :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)                 :: inv_rji, inv_rkj, inv_rkl
    real(wp)                 :: dotpro_ijk, dotpro_jkl
    real(wp)                 :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2
    real(wp)                 :: inv_raijkl
    real(wp)                 :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)                 :: v1_sqr, v2_sqr
    real(wp)                 :: u_t_1
    real(wp)                 :: u_t_2
    real(wp)                 :: A1, A2, A3, A4
    real(wp)                 :: B1
    real(wp)                 :: C1
    real(wp)                 :: cos_dih, sin_dih, u_dih, dih
    real(wp)                 :: cos_tmin, sin_tmin
    real(wp)                 :: cos_d_dih, sin_d_dih, d_dih
    real(wp)                 :: coef_dih_shift
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)
    real(wp)                 :: ftmp(1:3,1:4)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: fc(:)
    real(wp),        pointer :: t_min(:)
    real(wp),        pointer :: w(:)
    integer,         pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    t_min   => enefunc%dihe_theta_min
    w       => enefunc%dihe_w
    func    => enefunc%dihe_func
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, id,                     &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih, dih,       &
    !$omp         cos_tmin, sin_tmin,                 &
    !$omp         cos_d_dih, sin_d_dih, d_dih,        &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        fc, t_min, w, func, coord, force,    &
    !$omp        nthread)                             &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread

      if (func(i) /= 41) then
        cycle
      end if

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      !
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      !
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      !
      raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      ! 
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cos_tmin = cos(t_min(i))
      sin_tmin = sin(t_min(i))
      cos_d_dih = cos_dih * cos_tmin + sin_dih * sin_tmin
      sin_d_dih = sin_dih * cos_tmin - cos_dih * sin_tmin
      cos_d_dih = min( 1.0_wp, cos_d_dih)
      cos_d_dih = max(-1.0_wp, cos_d_dih)
      sin_d_dih = min( 1.0_wp, sin_d_dih)
      sin_d_dih = max(-1.0_wp, sin_d_dih)

      if (cos_d_dih > 1.0E-1_wp) then
        d_dih = asin(sin_d_dih)
      else
        d_dih = sign(1.0_wp, sin_d_dih)*acos(cos_d_dih)
      end if

      tmp   = d_dih / w(i)
      u_dih = fc(i) * (coef_dih_shift + exp(-0.5_wp * tmp * tmp))

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - (u_dih - fc(i) * coef_dih_shift) * tmp / w(i)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_gaussian_cos2mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_gaussian_sin3mod
  !> @brief        calculate dihedral energy in a safe way (sin^3) for Gaussian
  !! @authors      CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         M Bulacu et al., JCTC, 2013, 9, 3282-3292.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_gaussian_sin3mod(enefunc, coord, force, &
                                                      virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, l, m, n, id, omp_get_thread_num
    real(wp)                 :: tmp, etmp, var_tmp
    real(wp)                 :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)                 :: rji2, rkj2, rkl2
    real(wp)                 :: rji, rkj, rkl
    real(wp)                 :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)                 :: inv_rji, inv_rkj, inv_rkl
    real(wp)                 :: dotpro_ijk, dotpro_jkl
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2
    real(wp)                 :: inv_raijkl
    real(wp)                 :: cos_t_1, sin_t_1, t_1, u_t_1
    real(wp)                 :: cos_t_2, sin_t_2, t_2, u_t_2
    real(wp)                 :: cos_2_t, cos_2_t_m_1
    real(wp)                 :: cos_2_t_p_1, cos_2_t_p_1_sqr
    real(wp)                 :: A1, A2, A3, A4
    real(wp)                 :: B1, B2, B3, B4
    real(wp)                 :: C1, C2, C3, C4
    real(wp)                 :: cos_dih, sin_dih, u_dih, dih
    real(wp)                 :: cos_tmin, sin_tmin
    real(wp)                 :: cos_d_dih, sin_d_dih, d_dih
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)
    real(wp)                 :: ftmp(1:3,1:4)

    integer,         pointer :: list(:,:)
    real(wp),        pointer :: fc(:)
    real(wp),        pointer :: t_min(:)
    real(wp),        pointer :: w(:)
    integer,         pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = sin^3(t)
    ! 
    ! M'(t) = 3 sin^2(t) cos(t)
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    t_min   => enefunc%dihe_theta_min
    w       => enefunc%dihe_w
    func    => enefunc%dihe_func

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, id,                     &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         cos_t_1, sin_t_1, t_1, u_t_1,       &
    !$omp         cos_t_2, sin_t_2, t_2, u_t_2,       &
    !$omp         cos_2_t, cos_2_t_m_1,               &
    !$omp         cos_2_t_p_1, cos_2_t_p_1_sqr,       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, B2, B3, B4,                     &
    !$omp         C1, C2, C3, C4,                     &
    !$omp         cos_dih, sin_dih, u_dih, dih,       &
    !$omp         cos_tmin, sin_tmin,                 &
    !$omp         cos_d_dih, sin_d_dih, d_dih,        &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list,                  &
    !$omp        fc, t_min, w, func, coord, force,    &
    !$omp        nthread)                             &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread

      if (func(i) /= 43) then
        cycle
      end if

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      ! 
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      ! 
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      ! 
      raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      cos_t_1 = dotpro_ijk * inv_rji * inv_rkj
      cos_t_1 = min( 1.0_wp, cos_t_1)
      cos_t_1 = max(-1.0_wp, cos_t_1)
      sin_t_1 = sqrt(1.0_wp - cos_t_1 * cos_t_1)
      t_1     = acos(cos_t_1)
      !
      cos_t_2 = dotpro_jkl * inv_rkj * inv_rkl
      cos_t_2 = min( 1.0_wp, cos_t_2)
      cos_t_2 = max(-1.0_wp, cos_t_2)
      sin_t_2 = sqrt(1.0_wp - cos_t_2 * cos_t_2)
      t_2     = acos(cos_t_2)

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      u_t_1 = sin_t_1 * sin_t_1 * sin_t_1
      A1    = sin_t_1 * inv_rkj * inv_rji2
      A2    = sin_t_1 * inv_rji * inv_rkj2 * (-cos_t_1)
      B1    = 3.0_wp * cos_t_1 * sin_t_1

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      u_t_2 = sin_t_2 * sin_t_2 * sin_t_2
      A3    = sin_t_2 * inv_rkl * inv_rkj2 * cos_t_2
      A4    = sin_t_2 * inv_rkj * inv_rkl2
      C1    = 3.0_wp * cos_t_2 * sin_t_2

      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      
      cos_tmin = cos(t_min(i))
      sin_tmin = sin(t_min(i))
      cos_d_dih = cos_dih * cos_tmin + sin_dih * sin_tmin
      sin_d_dih = sin_dih * cos_tmin - cos_dih * sin_tmin
      cos_d_dih = min( 1.0_wp, cos_d_dih)
      cos_d_dih = max(-1.0_wp, cos_d_dih)
      sin_d_dih = min( 1.0_wp, sin_d_dih)
      sin_d_dih = max(-1.0_wp, sin_d_dih)

      if (cos_d_dih > 1.0E-1_wp) then
        d_dih = asin(sin_d_dih)
      else
        d_dih = sign(1.0_wp, sin_d_dih)*acos(cos_d_dih)
      end if

      tmp = d_dih / w(i)
      u_dih  = fc(i) * exp(-0.5_wp * tmp * tmp)

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - u_dih * tmp / w(i)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      tmp = u_dih * u_t_2 * B1
      B2  = inv_rji * inv_rkj
      B3  = inv_rji2
      B4  = inv_rkj2
      Q(1:3, 1) = tmp * (- B2 * dkj(1:3) - B3 * cos_t_1 * dji(1:3))
      Q(1:3, 3) = tmp * (+ B2 * dji(1:3) + B4 * cos_t_1 * dkj(1:3))
      Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
      Q(1:3, 4) = 0.0_wp

      ! ------
      ! part 3
      ! ------
      tmp = u_dih * u_t_1 * C1
      C2  = inv_rkj * inv_rkl
      C3  = inv_rkj2
      C4  = inv_rkl2
      R(1:3, 2) = tmp * (C2 * dkl(1:3) - C3 * cos_t_2 * dkj(1:3))
      R(1:3, 4) = tmp * (C2 * dkj(1:3) - C4 * cos_t_2 * dkl(1:3))
      R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
      R(1:3, 1) = 0.0_wp

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel 

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_gaussian_sin3mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_flexible_cos2mod
  !> @brief        calculate flexible dihedral energy in a safe way
  !! @authors      CK, CT
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_flexible_cos2mod(enefunc, coord, force, &
                                                      virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: istart, iend
    integer                  :: i, l, m, n, id, omp_get_thread_num
    integer                  :: dtype
    real(wp)                 :: tmp, etmp, var_tmp
    real(wp)                 :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)                 :: rji2, rkj2, rkl2
    real(wp)                 :: rji, rkj, rkl
    real(wp)                 :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)                 :: inv_rji, inv_rkj, inv_rkl
    real(wp)                 :: dotpro_ijk, dotpro_jkl
    real(wp)                 :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: raijk2, rajkl2
    real(wp)                 :: inv_raijk2, inv_rajkl2
    real(wp)                 :: inv_raijkl
    real(wp)                 :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)                 :: v1_sqr, v2_sqr
    real(wp)                 :: u_t_1
    real(wp)                 :: u_t_2
    real(wp)                 :: A1, A2, A3, A4
    real(wp)                 :: B1
    real(wp)                 :: C1
    real(wp)                 :: cos_dih, sin_dih, u_dih
    real(wp)                 :: cos_2dih, sin_2dih
    real(wp)                 :: cos_3dih, sin_3dih
    real(wp)                 :: coef_dih_shift
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)
    real(wp)                 :: ftmp(1:3,1:4)
    real(wp)                 :: c(7)

    integer,         pointer :: list(:,:)
    integer,         pointer :: types(:)
    real(wp),        pointer :: coef(:,:)
    real(wp),        pointer :: ener_corr(:)
    integer,         pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart         = enefunc%istart_dihedflex
    iend           = enefunc%iend_dihedflex
    list           => enefunc%diheflex_list
    types          => enefunc%diheflex_type
    func           => enefunc%diheflex_func
    coef           => enefunc%diheflex_coef
    ener_corr      => enefunc%diheflex_ener_corr
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, dtype, id,              &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih,            &
    !$omp         cos_2dih, sin_2dih, c,              &
    !$omp         cos_3dih, sin_3dih,                 &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        coef, types, func, ener_corr, coord, &
    !$omp        force, nthread)                      &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread

      if (func(i) /= 52) then
        cycle
      end if

      dtype = types(i)

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      !
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      !
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      !
      raijk2     = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2     = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)

      ! -----------------------------------------
      ! calculate cos^2 theta 1 and cos^2 theta 2
      ! -----------------------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      !
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      !
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      !
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cos_2dih = 2.0_wp * cos_dih * cos_dih - 1.0_wp
      sin_2dih = 2.0_wp * cos_dih * sin_dih
      cos_3dih = cos_2dih * cos_dih - sin_2dih * sin_dih
      sin_3dih = sin_2dih * cos_dih + cos_2dih * sin_dih
      c(1:7) = coef(1:7,dtype)
      u_dih = coef_dih_shift + AICG2P_K_DIHE * (c(1) &
          + c(2) * cos_dih  + c(3) * sin_dih         &
          + c(4) * cos_2dih + c(5) * sin_2dih        &
          + c(6) * cos_3dih + c(7) * sin_3dih        &
          - ener_corr(dtype))

      ! ==============
      ! Compute energy
      ! ==============
      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================
      grad_dih_coef =  AICG2P_K_DIHE * ( &
          - c(2)*sin_dih                 &
          + c(3)*cos_dih                 &
          - 2.0_wp*c(4)*sin_2dih         &
          + 2.0_wp*c(5)*cos_2dih         &
          - 3.0_wp*c(6)*sin_3dih         &
          + 3.0_wp*c(7)*cos_3dih)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_flexible_cos2mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed_pbc
  !> @brief        calculate dihedral energy
  !! @authors      CK, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed_pbc(enefunc, boundary, coord, force, virial, &
                                      edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: i, j, k, id, krot
    integer           :: i1, i2, i3, i4
    integer           :: istart, iend
    integer           :: aindex(1:4)
    real(wp)          :: cospha, sinpha, grad_coef, tmp
    real(wp)          :: cos_dih, sin_dih
    real(wp)          :: cosnt, sinnt, vtmp
    real(wp)          :: grad(1:9), v(1:3,1:3)
    real(wp)          :: bsize(3)

    real(wp), pointer :: fc(:), phase(:)
    real(wp), pointer :: work(:,:)
    integer,  pointer :: list(:,:)
    integer,  pointer :: nperiod(:)
    integer,  pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  = enefunc%istart_dihedral
    iend    = enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    nperiod => enefunc%dihe_periodicity
    phase   => enefunc%dihe_phase
    work    => enefunc%work
    func    => enefunc%dihe_func

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef, &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp)             &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, func, coord,  &
    !$omp        bsize)                                                      &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      if (func(i) == 21 .or. func(i) == 32 .or. func(i) == 33 .or. &
          func(i) == 41 .or. func(i) == 43) then
        cycle
      end if

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral_pbc(aindex, coord, bsize, cos_dih, sin_dih, grad, v)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot = 0
      do while (krot < nperiod(i))
        tmp   = cosnt*cos_dih - sinnt*sin_dih
        sinnt = sinnt*cos_dih + cosnt*sin_dih
        cosnt = tmp
        krot = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      edihe = edihe + fc(i) * (1.0_wp + cospha*cosnt + sinnt*sinpha)

      grad_coef = fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      ! ~CG~ 3SPN.2C DNA: only calculate functype /= 1 dihedrals...
      if (func(i) == 21 .or. func(i) == 32 .or. func(i) == 33 .or. &
          func(i) == 41 .or. func(i) == 43) then
        cycle
      end if
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_dihed_pbc
  !> @brief        calculate local dihedral energy based on PBC coords
  !! @authors      CK, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_dihed_pbc(enefunc, boundary, coord, force, &
                                            virial, edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: i, j, k, id, krot
    integer           :: i1, i2, i3, i4
    integer           :: istart, iend
    integer           :: aindex(1:4)
    real(wp)          :: grad_coef, tmp
    real(wp)          :: cos_dih, sin_dih, vtmp, dihe, etmp
    real(wp)          :: cospha, sinpha
    real(wp)          :: cosdif, sindif, diffphi
    real(wp)          :: grad(1:9), v(1:3,1:3)
    real(wp)          :: bsize(3)

    real(wp), pointer :: fc(:), theta_min(:), w(:)
    real(wp), pointer :: work(:,:)
    integer,  pointer :: list(:,:)
    integer,  pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart    = enefunc%istart_dihedral
    iend      = enefunc%iend_dihedral
    list      => enefunc%dihe_list
    fc        => enefunc%dihe_force_const
    theta_min => enefunc%dihe_theta_min
    w         => enefunc%dihe_w
    work      => enefunc%work
    func      => enefunc%dihe_func

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                           &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,  &
    !$omp         vtmp,  dihe, etmp, cospha, sinpha, cosdif, sindif, diffphi) &
    !$omp shared(istart, iend, list, fc, theta_min, w, work, func, coord,     &
    !$omp        bsize)                                                       &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

       if (func(i) /= 21) then
          cycle
       end if

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral_pbc(aindex, coord, bsize, cos_dih, sin_dih, grad,v)

      cospha = cos(theta_min(i))
      sinpha = sin(theta_min(i))

      cosdif = cos_dih*cospha + sin_dih*sinpha
      sindif = sin_dih*cospha - cos_dih*sinpha

      if (cosdif > 1.0E-1_wp) then
        diffphi = asin(sindif)
      else
        diffphi = sign(1.0_wp,sindif)*acos(cosdif)
      end if

      dihe = diffphi/w(i)

      etmp = fc(i)*exp(-0.5_wp*dihe*dihe)
      edihe = edihe + etmp
      grad_coef = (dihe*etmp)/w(i)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
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
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_local_dihed_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flexible_dihed_pbc
  !> @brief        calculate flexible dihedral energy from PBC coords
  !! @authors      CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_dihed_pbc(enefunc, boundary, coord, &
                                               force, virial, edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: i, j, k, id, krot
    integer           :: i1, i2, i3, i4
    integer           :: istart, iend
    integer           :: aindex(1:4), dtype
    real(wp)          :: cospha, sinpha, grad_coef, tmp, etmp
    real(wp)          :: cos_dih, sin_dih
    real(wp)          :: cos_2dih, sin_2dih
    real(wp)          :: cos_3dih, sin_3dih
    real(wp)          :: cosnt, sinnt, vtmp
    real(wp)          :: grad(1:9), v(1:3,1:3), c(7)
    real(wp)          :: bsize(3)

    real(wp), pointer :: fc(:), phase(:)
    real(wp), pointer :: coef(:,:)
    real(wp), pointer :: ener_corr(:)
    real(wp), pointer :: work(:,:)
    integer,  pointer :: list(:,:), types(:)
    integer,  pointer :: nperiod(:)
    integer,  pointer :: func(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  = enefunc%istart_dihedflex
    iend    = enefunc%iend_dihedflex
    list    => enefunc%diheflex_list
    coef    => enefunc%diheflex_coef
    types   => enefunc%diheflex_type
    func    => enefunc%diheflex_func
    work    => enefunc%work
    ener_corr => enefunc%diheflex_ener_corr

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef, &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, etmp,       &
    !$omp         cos_2dih, cos_3dih, sin_2dih, sin_3dih, dtype, c)          &
    !$omp shared(istart, iend, list, coef, types, func, work, coord,         &
    !$omp        ener_corr, bsize)                                           &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      if (func(i) /= 22) then
        cycle
      end if

      aindex(1:4) = list(1:4,i)
      dtype = types(i)
      call calculate_dihedral_pbc(aindex, coord, bsize, cos_dih, sin_dih, grad, v)

      cos_2dih =  2.0_wp*cos_dih*cos_dih-1.0_wp
      sin_2dih =  2.0_wp*cos_dih*sin_dih
      cos_3dih =  4.0_wp*cos_dih*cos_dih*cos_dih-3.0_wp*cos_dih
      sin_3dih = -4.0_wp*sin_dih*sin_dih*sin_dih+3.0_wp*sin_dih
      c(1:7) = coef(1:7,dtype)
      etmp =  c(1)                             &
            + c(2)*cos_dih  + c(3)*sin_dih     &
            + c(4)*cos_2dih + c(5)*sin_2dih    &
            + c(6)*cos_3dih + c(7)*sin_3dih    
     edihe = edihe + AICG2P_K_DIHE*(etmp-ener_corr(dtype))

     grad_coef = - c(2)*sin_dih                         & 
                 + c(3)*cos_dih                         &
                 - 2.0_wp*c(4)*sin_2dih &
                 + 2.0_wp*c(5)*cos_2dih &
                 - 3.0_wp*c(6)*sin_3dih &
                 + 3.0_wp*c(7)*cos_3dih 
      grad_coef = -AICG2P_K_DIHE*grad_coef

      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend

      if (func(i) /= 22) then
        cycle
      end if

      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_flexible_dihed_pbc


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_periodic_cos2mod_pbc
  !> @brief        calculate dihedral energy in a safe way
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_periodic_cos2mod_pbc(enefunc, boundary, &
                                                    coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: istart, iend
    integer           :: i, l, m, n, id, omp_get_thread_num
    integer           :: krot
    real(wp)          :: tmp, etmp, var_tmp
    real(wp)          :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)          :: rji2, rkj2, rkl2
    real(wp)          :: rji, rkj, rkl
    real(wp)          :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)          :: inv_rji, inv_rkj, inv_rkl
    real(wp)          :: dotpro_ijk, dotpro_jkl
    real(wp)          :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)          :: aijk(1:3), ajkl(1:3)
    real(wp)          :: raijk2, rajkl2
    real(wp)          :: inv_raijk2, inv_rajkl2
    real(wp)          :: inv_raijkl
    real(wp)          :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)          :: v1_sqr, v2_sqr
    real(wp)          :: u_t_1
    real(wp)          :: u_t_2
    real(wp)          :: A1, A2, A3, A4
    real(wp)          :: B1
    real(wp)          :: C1
    real(wp)          :: cos_dih, sin_dih, u_dih
    real(wp)          :: cospha, sinpha, cosnt, sinnt
    real(wp)          :: coef_dih_shift
    real(wp)          :: grad_dih_coef
    real(wp)          :: P(1:3,1:4)
    real(wp)          :: Q(1:3,1:4)
    real(wp)          :: R(1:3,1:4)
    real(wp)          :: ftmp(1:3,1:4)
    real(wp)          :: bsize(3), half_bsize(3)

    integer,  pointer :: list(:,:)
    real(wp), pointer :: fc(:)
    real(wp), pointer :: phase(:)
    integer,  pointer :: nperiod(:)
    integer,  pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    phase   => enefunc%dihe_phase
    nperiod => enefunc%dihe_periodicity
    func    => enefunc%dihe_func
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    bsize(1)        = boundary%box_size_x
    bsize(2)        = boundary%box_size_y
    bsize(3)        = boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, krot, id,               &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih,            &
    !$omp         cospha, sinpha, cosnt, sinnt,       &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        fc, phase, nperiod, func, coord,     &
    !$omp        bsize, half_bsize,                   &
    !$omp        force, nthread)                      &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread
      if (func(i) /= 32) then
        cycle
      end if

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      if (dji(1) > half_bsize(1)) then
        dji(1) = dji(1) - bsize(1)
      else if (dji(1) < -half_bsize(1)) then
        dji(1) = dji(1) + bsize(1)
      end if
      if (dji(2) > half_bsize(2)) then
        dji(2) = dji(2) - bsize(2)
      else if (dji(2) < -half_bsize(2)) then
        dji(2) = dji(2) + bsize(2)
      end if
      if (dji(3) > half_bsize(3)) then
        dji(3) = dji(3) - bsize(3)
      else if (dji(3) < -half_bsize(3)) then
        dji(3) = dji(3) + bsize(3)
      end if

      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      if (dkj(1) > half_bsize(1)) then
        dkj(1) = dkj(1) - bsize(1)
      else if (dkj(1) < -half_bsize(1)) then
        dkj(1) = dkj(1) + bsize(1)
      end if
      if (dkj(2) > half_bsize(2)) then
        dkj(2) = dkj(2) - bsize(2)
      else if (dkj(2) < -half_bsize(2)) then
        dkj(2) = dkj(2) + bsize(2)
      end if
      if (dkj(3) > half_bsize(3)) then
        dkj(3) = dkj(3) - bsize(3)
      else if (dkj(3) < -half_bsize(3)) then
        dkj(3) = dkj(3) + bsize(3)
      end if

      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))
      if (dkl(1) > half_bsize(1)) then
        dkl(1) = dkl(1) - bsize(1)
      else if (dkl(1) < -half_bsize(1)) then
        dkl(1) = dkl(1) + bsize(1)
      end if
      if (dkl(2) > half_bsize(2)) then
        dkl(2) = dkl(2) - bsize(2)
      else if (dkl(2) < -half_bsize(2)) then
        dkl(2) = dkl(2) + bsize(2)
      end if
      if (dkl(3) > half_bsize(3)) then
        dkl(3) = dkl(3) - bsize(3)
      else if (dkl(3) < -half_bsize(3)) then
        dkl(3) = dkl(3) + bsize(3)
      end if

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      ! 
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      ! 
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      ! 
      raijk2  = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2  = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      ! 
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot  = 0
        
      do while (krot < nperiod(i))
        tmp   = cosnt * cos_dih - sinnt * sin_dih
        sinnt = sinnt * cos_dih + cosnt * sin_dih
        cosnt = tmp
        krot  = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      u_dih  = fc(i) * (coef_dih_shift + 1.0_wp + cospha*cosnt + sinnt*sinpha)

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_periodic_cos2mod_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_gaussian_cos2mod_pbc
  !> @brief        calculate Gaussian-type dihedral energy in a safe way
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_gaussian_cos2mod_pbc(enefunc, boundary, &
                                                    coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: istart, iend
    integer           :: i, l, m, n, id, omp_get_thread_num
    real(wp)          :: tmp, etmp, var_tmp
    real(wp)          :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)          :: rji2, rkj2, rkl2
    real(wp)          :: rji, rkj, rkl
    real(wp)          :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)          :: inv_rji, inv_rkj, inv_rkl
    real(wp)          :: dotpro_ijk, dotpro_jkl
    real(wp)          :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)          :: aijk(1:3), ajkl(1:3)
    real(wp)          :: raijk2, rajkl2
    real(wp)          :: inv_raijk2, inv_rajkl2
    real(wp)          :: inv_raijkl
    real(wp)          :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)          :: v1_sqr, v2_sqr
    real(wp)          :: u_t_1
    real(wp)          :: u_t_2
    real(wp)          :: A1, A2, A3, A4
    real(wp)          :: B1
    real(wp)          :: C1
    real(wp)          :: cos_dih, sin_dih, u_dih, dih
    real(wp)          :: cos_tmin, sin_tmin
    real(wp)          :: cos_d_dih, sin_d_dih, d_dih
    real(wp)          :: coef_dih_shift
    real(wp)          :: grad_dih_coef
    real(wp)          :: P(1:3,1:4)
    real(wp)          :: Q(1:3,1:4)
    real(wp)          :: R(1:3,1:4)
    real(wp)          :: ftmp(1:3,1:4)
    real(wp)          :: bsize(3), half_bsize(3)

    integer,  pointer :: list(:,:)
    real(wp), pointer :: fc(:)
    real(wp), pointer :: t_min(:)
    real(wp), pointer :: w(:)
    integer,  pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  =  enefunc%istart_dihedral
    iend    =  enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    t_min   => enefunc%dihe_theta_min
    w       => enefunc%dihe_w
    func    => enefunc%dihe_func
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, id,                     &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih, dih,       &
    !$omp         cos_tmin, sin_tmin,                 &
    !$omp         cos_d_dih, sin_d_dih, d_dih,        &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        fc, t_min, w, func, coord, force,    &
    !$omp        bsize, half_bsize,                   &
    !$omp        nthread)                             &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread

      if (func(i) /= 41) then
        cycle
      end if

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      if (dji(1) > half_bsize(1)) then
        dji(1) = dji(1) - bsize(1)
      else if (dji(1) < -half_bsize(1)) then
        dji(1) = dji(1) + bsize(1)
      end if
      if (dji(2) > half_bsize(2)) then
        dji(2) = dji(2) - bsize(2)
      else if (dji(2) < -half_bsize(2)) then
        dji(2) = dji(2) + bsize(2)
      end if
      if (dji(3) > half_bsize(3)) then
        dji(3) = dji(3) - bsize(3)
      else if (dji(3) < -half_bsize(3)) then
        dji(3) = dji(3) + bsize(3)
      end if

      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      if (dkj(1) > half_bsize(1)) then
        dkj(1) = dkj(1) - bsize(1)
      else if (dkj(1) < -half_bsize(1)) then
        dkj(1) = dkj(1) + bsize(1)
      end if
      if (dkj(2) > half_bsize(2)) then
        dkj(2) = dkj(2) - bsize(2)
      else if (dkj(2) < -half_bsize(2)) then
        dkj(2) = dkj(2) + bsize(2)
      end if
      if (dkj(3) > half_bsize(3)) then
        dkj(3) = dkj(3) - bsize(3)
      else if (dkj(3) < -half_bsize(3)) then
        dkj(3) = dkj(3) + bsize(3)
      end if

      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))
      if (dkl(1) > half_bsize(1)) then
        dkl(1) = dkl(1) - bsize(1)
      else if (dkl(1) < -half_bsize(1)) then
        dkl(1) = dkl(1) + bsize(1)
      end if
      if (dkl(2) > half_bsize(2)) then
        dkl(2) = dkl(2) - bsize(2)
      else if (dkl(2) < -half_bsize(2)) then
        dkl(2) = dkl(2) + bsize(2)
      end if
      if (dkl(3) > half_bsize(3)) then
        dkl(3) = dkl(3) - bsize(3)
      else if (dkl(3) < -half_bsize(3)) then
        dkl(3) = dkl(3) + bsize(3)
      end if

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      !
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      !
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      !
      raijk2  = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2  = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      ! 
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cos_tmin = cos(t_min(i))
      sin_tmin = sin(t_min(i))
      cos_d_dih = cos_dih * cos_tmin + sin_dih * sin_tmin
      sin_d_dih = sin_dih * cos_tmin - cos_dih * sin_tmin
      cos_d_dih = min( 1.0_wp, cos_d_dih)
      cos_d_dih = max(-1.0_wp, cos_d_dih)
      sin_d_dih = min( 1.0_wp, sin_d_dih)
      sin_d_dih = max(-1.0_wp, sin_d_dih)

      if (cos_d_dih > 1.0E-1_wp) then
        d_dih = asin(sin_d_dih)
      else
        d_dih = sign(1.0_wp, sin_d_dih)*acos(cos_d_dih)
      end if

      tmp   = d_dih / w(i)
      u_dih = fc(i) * (coef_dih_shift + exp(-0.5_wp * tmp * tmp))

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef = - (u_dih - fc(i) * coef_dih_shift) * tmp / w(i)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_gaussian_cos2mod_pbc


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cg_dihed_flexible_cos2mod_pbc
  !> @brief        calculate Gaussian-type dihedral energy in a safe way
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine compute_energy_cg_dihed_flexible_cos2mod_pbc(enefunc, boundary, &
                                                    coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: edihe

    ! local variables
    integer           :: istart, iend
    integer           :: i, l, m, n, id, omp_get_thread_num
    integer           :: dtype
    real(wp)          :: tmp, etmp, var_tmp
    real(wp)          :: dji(1:3), dkj(1:3), dkl(1:3)
    real(wp)          :: rji2, rkj2, rkl2
    real(wp)          :: rji, rkj, rkl
    real(wp)          :: inv_rji2, inv_rkj2, inv_rkl2
    real(wp)          :: inv_rji, inv_rkj, inv_rkl
    real(wp)          :: dotpro_ijk, dotpro_jkl
    real(wp)          :: inv_dotpro_ijk, inv_dotpro_jkl
    real(wp)          :: aijk(1:3), ajkl(1:3)
    real(wp)          :: raijk2, rajkl2
    real(wp)          :: inv_raijk2, inv_rajkl2
    real(wp)          :: inv_raijkl
    real(wp)          :: u1, u2, v1, v2, w1, w2 ! see C.Tan et al. 2020
    real(wp)          :: v1_sqr, v2_sqr
    real(wp)          :: u_t_1
    real(wp)          :: u_t_2
    real(wp)          :: A1, A2, A3, A4
    real(wp)          :: B1
    real(wp)          :: C1
    real(wp)          :: cos_dih, sin_dih, u_dih, dih
    real(wp)          :: cos_2dih, sin_2dih
    real(wp)          :: cos_3dih, sin_3dih
    real(wp)          :: coef_dih_shift
    real(wp)          :: grad_dih_coef
    real(wp)          :: P(1:3,1:4)
    real(wp)          :: Q(1:3,1:4)
    real(wp)          :: R(1:3,1:4)
    real(wp)          :: ftmp(1:3,1:4)
    real(wp)          :: c(7)
    real(wp)          :: bsize(3), half_bsize(3)

    integer,  pointer :: list(:,:)
    integer,  pointer :: types(:)
    real(wp), pointer :: coef(:,:)
    real(wp), pointer :: ener_corr(:)
    integer,  pointer :: func(:)


    ! =======================================================
    ! New dihedral potential with angle modulating functions:
    ! 
    ! U = f(phi) * M(t1) * M(t2)
    ! 
    ! M(t) = cos^2[3*(t-5pi/6)]
    !      = sin^2(t) * [2 cos(2t) + 1]^2
    ! 
    ! M'(t) = 6 sin(t) cos(t) [2 cos(2t) - 1] [2 cos(2t) + 1]
    ! =======================================================
    
    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart         = enefunc%istart_dihedflex
    iend           = enefunc%iend_dihedflex
    list           => enefunc%diheflex_list
    types          => enefunc%diheflex_type
    func           => enefunc%diheflex_func
    coef           => enefunc%diheflex_coef
    ener_corr      => enefunc%diheflex_ener_corr
    coef_dih_shift = enefunc%cg_safe_dih_ene_shift

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel default(none)                      &
    !$omp private(i, l, m, n, dtype, id,              &
    !$omp         tmp, etmp, var_tmp,                 &
    !$omp         dji, dkj, dkl,                      &
    !$omp         rji2, rkj2, rkl2, rji, rkj, rkl,    &
    !$omp         inv_rji2, inv_rkj2, inv_rkl2,       &
    !$omp         inv_rji, inv_rkj, inv_rkl,          &
    !$omp         dotpro_ijk, dotpro_jkl,             &
    !$omp         inv_dotpro_ijk, inv_dotpro_jkl,     &
    !$omp         aijk, ajkl, raijk2, rajkl2,         &
    !$omp         inv_raijk2, inv_rajkl2, inv_raijkl, &
    !$omp         u1, u2, v1, v2, w1, w2,             &
    !$omp         v1_sqr, v2_sqr,                     &
    !$omp         u_t_1, u_t_2,                       &
    !$omp         A1, A2, A3, A4,                     &
    !$omp         B1, C1,                             &
    !$omp         cos_dih, sin_dih, u_dih,            &
    !$omp         cos_2dih, sin_2dih, c,              &
    !$omp         cos_3dih, sin_3dih,                 &
    !$omp         grad_dih_coef, P, Q, R, ftmp)       &
    !$omp shared(istart, iend, list, coef_dih_shift,  &
    !$omp        coef, types, func, ener_corr,        &
    !$omp        bsize, half_bsize,                   &
    !$omp        coord, force,                        &
    !$omp        nthread)                             &
    !$omp reduction(+:virial) reduction(+:edihe)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    do i = istart+id, iend, nthread

      if (func(i) /= 52) then
        cycle
      end if

      dtype = types(i)

      dji(1:3) = coord(1:3, list(1, i)) - coord(1:3, list(2, i))
      if (dji(1) > half_bsize(1)) then
        dji(1) = dji(1) - bsize(1)
      else if (dji(1) < -half_bsize(1)) then
        dji(1) = dji(1) + bsize(1)
      end if
      if (dji(2) > half_bsize(2)) then
        dji(2) = dji(2) - bsize(2)
      else if (dji(2) < -half_bsize(2)) then
        dji(2) = dji(2) + bsize(2)
      end if
      if (dji(3) > half_bsize(3)) then
        dji(3) = dji(3) - bsize(3)
      else if (dji(3) < -half_bsize(3)) then
        dji(3) = dji(3) + bsize(3)
      end if

      dkj(1:3) = coord(1:3, list(2, i)) - coord(1:3, list(3, i))
      if (dkj(1) > half_bsize(1)) then
        dkj(1) = dkj(1) - bsize(1)
      else if (dkj(1) < -half_bsize(1)) then
        dkj(1) = dkj(1) + bsize(1)
      end if
      if (dkj(2) > half_bsize(2)) then
        dkj(2) = dkj(2) - bsize(2)
      else if (dkj(2) < -half_bsize(2)) then
        dkj(2) = dkj(2) + bsize(2)
      end if
      if (dkj(3) > half_bsize(3)) then
        dkj(3) = dkj(3) - bsize(3)
      else if (dkj(3) < -half_bsize(3)) then
        dkj(3) = dkj(3) + bsize(3)
      end if

      dkl(1:3) = coord(1:3, list(4, i)) - coord(1:3, list(3, i))
      if (dkl(1) > half_bsize(1)) then
        dkl(1) = dkl(1) - bsize(1)
      else if (dkl(1) < -half_bsize(1)) then
        dkl(1) = dkl(1) + bsize(1)
      end if
      if (dkl(2) > half_bsize(2)) then
        dkl(2) = dkl(2) - bsize(2)
      else if (dkl(2) < -half_bsize(2)) then
        dkl(2) = dkl(2) + bsize(2)
      end if
      if (dkl(3) > half_bsize(3)) then
        dkl(3) = dkl(3) - bsize(3)
      else if (dkl(3) < -half_bsize(3)) then
        dkl(3) = dkl(3) + bsize(3)
      end if

      rji2     = dji(1)*dji(1) + dji(2)*dji(2) + dji(3)*dji(3)
      rkj2     = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
      rkl2     = dkl(1)*dkl(1) + dkl(2)*dkl(2) + dkl(3)*dkl(3)
      inv_rji2 = 1.0_wp / rji2
      inv_rkj2 = 1.0_wp / rkj2
      inv_rkl2 = 1.0_wp / rkl2
      rji      = sqrt(rji2)
      rkj      = sqrt(rkj2)
      rkl      = sqrt(rkl2)
      inv_rji  = 1.0_wp / rji
      inv_rkj  = 1.0_wp / rkj
      inv_rkl  = 1.0_wp / rkl

      ! -----------------------------
      ! dot products for ji*kj, kj*kl
      ! -----------------------------
      !
      dotpro_ijk = - dji(1)*dkj(1) - dji(2)*dkj(2) - dji(3)*dkj(3)
      dotpro_jkl =   dkj(1)*dkl(1) + dkj(2)*dkl(2) + dkj(3)*dkl(3)
      inv_dotpro_ijk = 1.0_wp / dotpro_ijk
      inv_dotpro_jkl = 1.0_wp / dotpro_jkl

      ! ---------------------------
      ! cross products for ijk, jkl
      ! ---------------------------
      !
      aijk(1) = dji(2)*dkj(3) - dji(3)*dkj(2)
      aijk(2) = dji(3)*dkj(1) - dji(1)*dkj(3)
      aijk(3) = dji(1)*dkj(2) - dji(2)*dkj(1)
      !
      ajkl(1) = dkl(2)*dkj(3) - dkl(3)*dkj(2)
      ajkl(2) = dkl(3)*dkj(1) - dkl(1)*dkj(3)
      ajkl(3) = dkl(1)*dkj(2) - dkl(2)*dkj(1)
      !
      raijk2  = aijk(1)*aijk(1) + aijk(2)*aijk(2) + aijk(3)*aijk(3)
      rajkl2  = ajkl(1)*ajkl(1) + ajkl(2)*ajkl(2) + ajkl(3)*ajkl(3)
      if (raijk2 < EPS) then
        raijk2 = EPS
      end if
      if (rajkl2 < EPS) then
        rajkl2 = EPS
      end if
      inv_raijk2 = 1.0_wp / raijk2
      inv_rajkl2 = 1.0_wp / rajkl2
      inv_raijkl = sqrt(inv_raijk2*inv_rajkl2)
 
      ! -----------------------------
      ! calculate theta 1 and theta 2
      ! -----------------------------
      !
      u1 = dotpro_ijk * dotpro_ijk * inv_rji2 * inv_rkj2
      !
      u2 = dotpro_jkl * dotpro_jkl * inv_rkl2 * inv_rkj2

      ! -----------------------------
      ! calculate U(t_1) and grad t_1
      ! -----------------------------
      !
      if (u1 < 0.75_wp) then
        u_t_1 = 1
        !
        A1 = rkj * inv_raijk2
        A2 = - dotpro_ijk * inv_raijk2 * inv_rkj
        B1 = 0.0_wp
      else
        v1     = 4.0_wp * u1 - 1.0_wp
        v1_sqr = v1 * v1
        w1     = v1 - 2.0_wp
        u_t_1  = (1.0_wp - u1) * v1_sqr
        !
        A1 = v1_sqr * inv_rkj * inv_rji2
        A2 = - v1_sqr * u1 * inv_rkj * inv_dotpro_ijk
        B1 = 6.0_wp * u1 * v1 * w1
      end if

      ! -----------------------------
      ! calculate U(t_2) and grad t_2
      ! -----------------------------
      !
      if (u2 < 0.75_wp) then
        u_t_2 = 1
        !
        A3 = dotpro_jkl * inv_rajkl2 * inv_rkj
        A4 = rkj * inv_rajkl2
        C1 = 0.0_wp
      else
        v2     = 4.0_wp * u2 - 1.0_wp
        v2_sqr = v2 * v2
        w2     = v2 - 2.0_wp
        u_t_2  = (1.0_wp - u2) * v2_sqr
        !
        A3 = v2_sqr * u2 * inv_rkj * inv_dotpro_jkl
        A4 = v2_sqr * inv_rkj * inv_rkl2
        C1 = 6.0_wp * u2 * v2 * w2
      end if

      ! 
      ! ---------------------------
      ! compute cos_dih and sin_dih
      ! ---------------------------
      ! 
      cos_dih = (aijk(1)*ajkl(1) + aijk(2)*ajkl(2) + aijk(3)*ajkl(3))*inv_raijkl
      cos_dih = min( 1.0_wp, cos_dih)
      cos_dih = max(-1.0_wp, cos_dih)
      ! 
      tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
      sin_dih = tmp * rkj * inv_raijkl
      sin_dih = min( 1.0_wp, sin_dih)
      sin_dih = max(-1.0_wp, sin_dih)

      cos_2dih = 2.0_wp * cos_dih * cos_dih - 1.0_wp
      sin_2dih = 2.0_wp * cos_dih * sin_dih
      cos_3dih = cos_2dih * cos_dih - sin_2dih * sin_dih
      sin_3dih = sin_2dih * cos_dih + cos_2dih * sin_dih
      c(1:7) = coef(1:7,dtype)
      u_dih = coef_dih_shift + AICG2P_K_DIHE * (c(1) &
          + c(2) * cos_dih  + c(3) * sin_dih         &
          + c(4) * cos_2dih + c(5) * sin_2dih        &
          + c(6) * cos_3dih + c(7) * sin_3dih        &
          - ener_corr(dtype))

      ! ==============
      ! Compute energy
      ! ==============

      etmp  = u_dih * u_t_1 * u_t_2
      edihe = edihe + etmp

      ! ========================
      ! Calculate all the forces
      ! ========================

      grad_dih_coef =  AICG2P_K_DIHE * ( &
          - c(2)*sin_dih                 &
          + c(3)*cos_dih                 &
          - 2.0_wp*c(4)*sin_2dih         &
          + 2.0_wp*c(5)*cos_2dih         &
          - 3.0_wp*c(6)*sin_3dih         &
          + 3.0_wp*c(7)*cos_3dih)

      ! ------
      ! part 1
      ! ------
      P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
      P(1:3, 2) = (u_t_2 * grad_dih_coef * (- A1 - A2)) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
      P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
      P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)
 
      ! ------
      ! part 2
      ! ------
      if (B1 < EPS) then
        Q(1:3, 1:4) = 0.0_wp
      else
        tmp = u_dih * u_t_2 * B1
        Q(1:3, 1) = tmp * (- inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
        Q(1:3, 3) = tmp * (+ inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
        Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
        Q(1:3, 4) = 0.0_wp
      end if

      ! ------
      ! part 3
      ! ------
      if (C1 < EPS) then
        R(1:3, 1:4) = 0.0_wp
      else
        tmp =  u_dih * u_t_1 * C1
        R(1:3, 2) = tmp * (inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
        R(1:3, 4) = tmp * (inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
        R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
        R(1:3, 1) = 0.0_wp
      end if

      ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

      force(1:3, list(1,i),id+1) = force(1:3, list(1,i),id+1) + ftmp(1:3, 1)
      force(1:3, list(2,i),id+1) = force(1:3, list(2,i),id+1) + ftmp(1:3, 2)
      force(1:3, list(3,i),id+1) = force(1:3, list(3,i),id+1) + ftmp(1:3, 3)
      force(1:3, list(4,i),id+1) = force(1:3, list(4,i),id+1) + ftmp(1:3, 4)

      ! virial
      do l = 1, 3
        do m = l+1, 3
          var_tmp = 0.0_wp
          do n = 1, 4
            var_tmp = var_tmp + ftmp(l, n) * coord(m, n) 
          end do
          virial(m, l) = virial(m, l) - var_tmp
          virial(l, m) = virial(l, m) - var_tmp
        end do
        var_tmp = 0.0_wp
        do n = 1, 4
          var_tmp = var_tmp + ftmp(l, n) * coord(l, n) 
        end do
        virial(l,l) = virial(l,l) - var_tmp
      end do
      !
    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cg_dihed_flexible_cos2mod_pbc

end module at_energy_dihedrals_mod
