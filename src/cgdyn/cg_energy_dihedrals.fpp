!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_dihedrals_mod
!> @brief   Calculates dihedral energy.
!! @authors Chigusa Kobayashi (CK), Jaewoon Jung (JJ) , Takao Yoda (TY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_dihedrals_mod

  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use dihedral_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public :: compute_energy_dihed
  public :: compute_energy_local_dihed
  public :: compute_energy_flexible_dihed

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed
  !> @brief        Calculates dihedral energy.
  !! @authors      CK, JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed(domain, enefunc, coord, force, edihe, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    integer                  :: i, j, k, id
    integer                  :: ix
    integer                  :: i1, i2, i3, i4
    integer                  :: nperiod_temp
    real(wp)                 :: edihe_temp, etmp, fc_temp, phase_temp
    real(wp)                 :: cwork(1:3,1:4), ftmp(1:3,1:4)
    real(wp)                 :: work(1:9), viri(1:3)

    real(wp),        pointer :: fc(:), phase(:)
    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: nperiod(:), func(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    dihelist       => enefunc%dihe_list
    fc             => enefunc%dihe_force_const
    nperiod        => enefunc%dihe_periodicity
    phase          => enefunc%dihe_phase
    func           => enefunc%dihe_kind

    !$omp parallel default(shared)                                &
    !$omp private(i, j, k, id, cwork, ix, work, edihe_temp, ftmp, &
    !$omp         fc_temp, phase_temp, nperiod_temp, etmp, viri,  &
    !$omp         i1, i2, i3, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    edihe_temp = 0.0_wp

    do i = id+1, enefunc%num_dihe_domain, nthread

      ix = enefunc%num_dihe_flexible_domain+enefunc%num_dihe_local_domain+i

      i1    = id_g2l(dihelist(1,ix))
      i2    = id_g2l(dihelist(2,ix))
      i3    = id_g2l(dihelist(3,ix))
      i4    = id_g2l(dihelist(4,ix))

      cwork(1:3,1) = coord(i1,1:3)
      cwork(1:3,2) = coord(i2,1:3)
      cwork(1:3,3) = coord(i3,1:3)
      cwork(1:3,4) = coord(i4,1:3)

      fc_temp = fc(ix)
      phase_temp = phase(ix)
      nperiod_temp = nperiod(ix)

      if (func(ix) == 1 .or. func(ix) == 31) then

        call dihed(cwork, fc_temp, phase_temp, nperiod_temp, etmp, work, viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst1',my_city_rank,etmp
        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      else if (func(ix) == 32) then

        call dihed_cos2mod(cwork, fc_temp, phase_temp, nperiod_temp,  &
                           enefunc%cg_safe_dih_ene_shift, etmp, ftmp, viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst2',my_city_rank,etmp
        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) + ftmp(1:3,1)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + ftmp(1:3,2)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + ftmp(1:3,3)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) + ftmp(1:3,4)

      end if

      virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

    end do

    edihe(id+1) = edihe(id+1) + edihe_temp

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_dihed
  !> @brief        Calculates local dihedral energy.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_dihed(domain, enefunc, coord, force, edihe, &
                                        virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    integer                  :: i, j, k, id
    integer                  :: ix
    integer                  :: i1, i2, i3, i4
    real(wp)                 :: edihe_temp, etmp, viri(1:3)
    real(wp)                 :: fc_temp, phase_temp, width_temp
    real(wp)                 :: work(1:9), cwork(1:3,1:4), ftmp(1:3,1:4)

    real(wp),        pointer :: fc(:), phase(:), width(:)
    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: func(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    dihelist       => enefunc%dihe_list
    fc             => enefunc%dihe_force_const
    phase          => enefunc%dihe_phase
    width          => enefunc%dihe_width
    func           => enefunc%dihe_kind

    !$omp parallel default(shared)                            &
    !$omp private(i, j, k, id, work, etmp, cwork, edihe_temp, &
    !$omp         fc_temp, phase_temp, width_temp, ftmp,      &
    !$omp         viri, i1, i2, i3, i4, ix)
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    edihe_temp = 0.0_wp

    do i = id+1, enefunc%num_dihe_local_domain, nthread

      ix = enefunc%num_dihe_flexible_domain + i

      i1    = id_g2l(dihelist(1,ix))
      i2    = id_g2l(dihelist(2,ix))
      i3    = id_g2l(dihelist(3,ix))
      i4    = id_g2l(dihelist(4,ix))

      cwork(1:3,1) = coord(i1,1:3)
      cwork(1:3,2) = coord(i2,1:3)
      cwork(1:3,3) = coord(i3,1:3)
      cwork(1:3,4) = coord(i4,1:3)
      fc_temp = fc(ix)
      phase_temp = phase(ix)
      width_temp = width(ix)

      if (func(ix) == 21) then

        call local_dihed(cwork, fc_temp, phase_temp, width_temp, &
                         etmp, work, viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst3',my_city_rank,etmp
        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      else if (func(ix) == 41) then

        call local_dihed_cos2mod(cwork, fc_temp, phase_temp, width_temp,    &
                                 enefunc%cg_safe_dih_ene_shift, etmp, ftmp, &
                                 viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst4',my_city_rank,etmp

        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) + ftmp(1:3,1)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + ftmp(1:3,2)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + ftmp(1:3,3)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) + ftmp(1:3,4)

      end if        

      virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

    end do

    edihe(id+1) = edihe(id+1) + edihe_temp

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_local_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flexible_dihed
  !> @brief        Calculates flexible dihedral angle energy.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flexible_dihed(domain, enefunc, coord, force, &
                                           edihe, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    integer                  :: i, ix, k, id
    integer                  :: i1, i2, i3, i4
    integer                  :: dtype
    real(wp)                 :: edihe_temp, vtmp, etmp, viri(1:3)
    real(wp)                 :: cwork(1:3,1:4), ftmp(1:3,1:4)
    real(wp)                 :: work(1:9), c(1:7)

    integer,         pointer :: dihelist(:,:)
    integer,         pointer :: dihe_type(:), dihe_func(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)
    real(wp),        pointer :: coef(:,:), ener_corr(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    dihelist       => enefunc%dihe_list
    coef           => enefunc%diheflex_coef
    ener_corr      => enefunc%diheflex_ener_corr
    dihe_type      => enefunc%dihe_kind
    dihe_func      => enefunc%dihe_periodicity

    !$omp parallel default(shared)                                       &
    !$omp private(i, ix, k, id, cwork, vtmp, c, etmp, dtype, edihe_temp, &
    !$omp         ftmp, work, viri, i1, i2, i3, i4)
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    edihe_temp = 0.0_wp

    do i = id+1, enefunc%num_dihe_flexible_domain, nthread

      i1    = id_g2l(dihelist(1,i))
      i2    = id_g2l(dihelist(2,i))
      i3    = id_g2l(dihelist(3,i))
      i4    = id_g2l(dihelist(4,i))

      cwork(1:3,1) = coord(i1,1:3)
      cwork(1:3,2) = coord(i2,1:3)
      cwork(1:3,3) = coord(i3,1:3)
      cwork(1:3,4) = coord(i4,1:3)

      dtype  = dihe_type(i)
      c(1:7) = coef(1:7,dtype)
      vtmp   = ener_corr(dtype)

      if (dihe_func(i) == 22) then

        call dihed_flexible(cwork, vtmp, c, etmp, work, viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst5',my_city_rank,etmp
        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      else if (dihe_func(i) == 52) then

        call dihed_flexible_cos2mod(cwork, vtmp, c, &
                   enefunc%cg_safe_dih_ene_shift, etmp, ftmp, viri)
        if (.not.(etmp < 100.0_wp)) write(*,*) 'tst6',my_city_rank,etmp
        edihe_temp = edihe_temp + etmp
        force(i1,1:3,id+1) = force(i1,1:3,id+1) + ftmp(1:3,1)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + ftmp(1:3,2)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + ftmp(1:3,3)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) + ftmp(1:3,4)

      end if

      virial(1,1,id+1) = virial(1,1,id+1) + viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) + viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) + viri(3)

    end do

    edihe(id+1) = edihe(id+1) + edihe_temp

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_flexible_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dihed
  !> @brief        Calculates dihedral energy called from energy_local_dihed.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dihed(cwork, fc, phase, nperiod, etmp, grad, viri)

    ! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: fc, phase
    integer,                 intent(in)    :: nperiod
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: grad(:)
    real(wp),                intent(inout) :: viri(:)

    ! local variables
    integer                  :: aindex(1:4), krot
    real(wp)                 :: cospha, sinpha, cos_dih, sin_dih
    real(wp)                 :: grad_coef, tmp
    real(wp)                 :: cosnt, sinnt
    real(wp)                 :: v(1:3,1:3)

    aindex(1) = 1
    aindex(2) = 2
    aindex(3) = 3
    aindex(4) = 4
    call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

    cosnt = 1.0_wp
    sinnt = 0.0_wp
    krot  = 0
    do while (krot < nperiod)
      tmp   = cosnt*cos_dih - sinnt*sin_dih
      sinnt = sinnt*cos_dih + cosnt*sin_dih
      cosnt = tmp 
      krot  = krot + 1
    end do

    cospha = cos(phase)
    sinpha = sin(phase)

    etmp  = fc*(1.0_wp + cospha*cosnt + sinnt*sinpha)

    grad_coef = fc*real(nperiod,wp)*(cospha*sinnt-cosnt*sinpha)
    grad(1:9) = grad_coef*grad(1:9)
    viri(1)   = grad_coef*v(1,1)
    viri(2)   = grad_coef*v(2,2)
    viri(3)   = grad_coef*v(3,3)

    return

  end subroutine dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dihed
  !> @brief        Calculates dihedral energy called from energy_local_dihed.  
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dihed_flexible(cwork, vtmp, c, etmp, grad, viri)

! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: vtmp
    real(wp),                intent(in)    :: c(:)
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: grad(:)
    real(wp),                intent(inout) :: viri(:)

! local variables
    integer                  :: aindex(1:4)
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cos_2dih, sin_2dih
    real(wp)                 :: cos_3dih, sin_3dih
    real(wp)                 :: grad_coef, tmp
    real(wp)                 :: v(1:3,1:3)

    aindex(1) = 1
    aindex(2) = 2
    aindex(3) = 3
    aindex(4) = 4
    call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

    cos_2dih =  2.0_wp*cos_dih*cos_dih - 1.0_wp
    sin_2dih =  2.0_wp*cos_dih*sin_dih
    cos_3dih =  4.0_wp*cos_dih*cos_dih*cos_dih - 3.0_wp*cos_dih
    sin_3dih = -4.0_wp*sin_dih*sin_dih*sin_dih + 3.0_wp*sin_dih

    tmp = c(1) + c(2)*cos_dih  + c(3)*sin_dih     &
               + c(4)*cos_2dih + c(5)*sin_2dih    &
               + c(6)*cos_3dih + c(7)*sin_3dih

    etmp = AICG2P_K_DIHE * (tmp - vtmp)

    grad_coef = -        c(2)*sin_dih  +        c(3)*cos_dih  &
                - 2.0_wp*c(4)*sin_2dih + 2.0_wp*c(5)*cos_2dih &
                - 3.0_wp*c(6)*sin_3dih + 3.0_wp*c(7)*cos_3dih
    grad_coef = -AICG2P_K_DIHE * grad_coef
    grad(1:9) = grad_coef*grad(1:9)

    viri(1) = grad_coef*v(1,1)
    viri(2) = grad_coef*v(2,2)
    viri(3) = grad_coef*v(3,3)

    return

  end subroutine dihed_flexible

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    local_dihed
  !> @brief        local dihedral energy called from energy_local_dihed
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine local_dihed(cwork, fc, phase, width, etmp, grad, viri)

    ! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: fc, phase, width
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: grad(:)
    real(wp),                intent(inout) :: viri(:)

    ! local variables
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, cos_dih, sin_dih
    real(wp)                 :: grad_coef, dihe
    real(wp)                 :: cosdif, sindif, diffphi, theta
    real(wp)                 :: v(1:3,1:3)

    aindex(1) = 1
    aindex(2) = 2
    aindex(3) = 3
    aindex(4) = 4
    call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

    theta = phase
    cospha = cos(theta)
    sinpha = sin(theta)

    cosdif = cos_dih*cospha + sin_dih*sinpha
    sindif = sin_dih*cospha - cos_dih*sinpha
    sindif = min(sindif, 1.0_wp-EPS)
    sindif = max(sindif, -1.0_wp+EPS)
    cosdif = min(cosdif, 1.0_wp-EPS)
    cosdif = max(cosdif, -1.0_wp+EPS)

    if (cosdif > 0.1_wp) then
      diffphi = asin(sindif)
    else
      diffphi = sign(1.0_wp,sindif)*acos(cosdif)
    end if

    dihe = diffphi / width
    etmp  = fc*exp(-0.5_wp*dihe*dihe)

    grad_coef = (dihe*etmp)/width
    grad(1:9) = grad_coef*grad(1:9)

    viri(1) = grad_coef * v(1,1)
    viri(2) = grad_coef * v(2,2)
    viri(3) = grad_coef * v(3,3)

    return

  end subroutine local_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dihed_cos2mod
  !> @brief        local dihedral energy called from energy_local_dihed
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dihed_cos2mod(cwork, fc, phase, nperiod, coef_dih_shift, &
                           etmp, ftmp, viri)

    ! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: fc, phase
    integer,                 intent(in)    :: nperiod
    real(wp),                intent(in)    :: coef_dih_shift
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: ftmp(:,:)
    real(wp),                intent(inout) :: viri(:)

    ! local variables
    integer                  :: krot, i
    real(wp)                 :: tmp
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
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)

    dji(1:3) = cwork(1:3,1) - cwork(1:3,2)
    dkj(1:3) = cwork(1:3,2) - cwork(1:3,3)
    dkl(1:3) = cwork(1:3,4) - cwork(1:3,3)

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
    if ( raijk2 < EPS ) then
      raijk2 = EPS
    end if
    if ( rajkl2 < EPS ) then
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
    if ( u1 < 0.75_wp ) then
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
    if ( u2 < 0.75_wp ) then
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
    cos_dih = min(  1.0_wp-EPS, cos_dih )
    cos_dih = max( -1.0_wp+EPS, cos_dih )
      !
    tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
    sin_dih = tmp * rkj * inv_raijkl
    sin_dih = min(  1.0_wp-EPS, sin_dih )
    sin_dih = max( -1.0_wp+EPS, sin_dih )

    cosnt = 1.0_wp
    sinnt = 0.0_wp
    krot  = 0

    do while (krot < nperiod)
      tmp   = cosnt * cos_dih - sinnt * sin_dih
      sinnt = sinnt * cos_dih + cosnt * sin_dih
      cosnt = tmp
      krot  = krot+1
    end do

    cospha = cos(phase)
    sinpha = sin(phase)

    u_dih  = fc * ( coef_dih_shift + 1.0_wp + cospha*cosnt + sinnt*sinpha)

    ! ==============
    ! Compute energy
    ! ==============

    etmp  = u_dih * u_t_1 * u_t_2

    ! ========================
    ! Calculate all the forces
    ! ========================

    grad_dih_coef = - fc * real(nperiod,wp) * (cospha*sinnt - cosnt*sinpha)

    ! ------
    ! part 1
    ! ------
    P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
    P(1:3, 2) = ( u_t_2 * grad_dih_coef * (- A1 - A2) ) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
    P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
    P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)

    ! ------
    ! part 2
    ! ------
    if ( B1 < EPS ) then
      Q(1:3, 1:4) = 0.0_wp
    else
      tmp = u_dih * u_t_2 * B1
      Q(1:3, 1) = tmp * ( - inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
      Q(1:3, 3) = tmp * ( + inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
      Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
      Q(1:3, 4) = 0.0_wp
    end if

    ! ------
    ! part 3
    ! ------
    if ( C1 < EPS ) then
      R(1:3, 1:4) = 0.0_wp
    else
      tmp =  u_dih * u_t_1 * C1
      R(1:3, 2) = tmp * ( inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
      R(1:3, 4) = tmp * ( inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
      R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
      R(1:3, 1) = 0.0_wp
    end if

    ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

    viri(1:3) = 0.0_wp
    do i = 1, 4
      viri(1:3) = viri(1:3) + ftmp(1:3,i)*cwork(1:3,i)
    end do

    return

  end subroutine dihed_cos2mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dihed_flexible_cos2mod
  !> @brief        dihedral energy called from energy_local_dihed
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dihed_flexible_cos2mod(cwork, vtmp, c, coef_dih_shift, &
                                    etmp, ftmp, viri)

  ! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: vtmp
    real(wp),                intent(in)    :: c(:)
    real(wp),                intent(in)    :: coef_dih_shift
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: ftmp(:,:)
    real(wp),                intent(inout) :: viri(:)

  ! local variables
    integer                  :: i
    real(wp)                 :: tmp
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
    real(wp)                 :: grad_dih_coef
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)

    dji(1:3) = cwork(1:3,1) - cwork(1:3,2)
    dkj(1:3) = cwork(1:3,2) - cwork(1:3,3)
    dkl(1:3) = cwork(1:3,4) - cwork(1:3,3)

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
    if ( raijk2 < EPS ) then
      raijk2 = EPS
    end if
    if ( rajkl2 < EPS ) then
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
    if ( u1 < 0.75_wp ) then
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
    if ( u2 < 0.75_wp ) then
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
    cos_dih = min(  1.0_wp-EPS, cos_dih )
    cos_dih = max( -1.0_wp+EPS, cos_dih )
!
    tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
    sin_dih = tmp * rkj * inv_raijkl
    sin_dih = min(  1.0_wp-EPS, sin_dih )
    sin_dih = max( -1.0_wp+EPS, sin_dih )

    cos_2dih = 2.0_wp * cos_dih * cos_dih - 1.0_wp
    sin_2dih = 2.0_wp * cos_dih * sin_dih
    cos_3dih = cos_2dih * cos_dih - sin_2dih * sin_dih
    sin_3dih = sin_2dih * cos_dih + cos_2dih * sin_dih
    u_dih = coef_dih_shift + AICG2P_K_DIHE * (c(1) &
          + c(2) * cos_dih  + c(3) * sin_dih         &
          + c(4) * cos_2dih + c(5) * sin_2dih        &
          + c(6) * cos_3dih + c(7) * sin_3dih        &
          - vtmp)

  ! ==============
  ! Compute energy
  ! ==============
    etmp  = u_dih * u_t_1 * u_t_2

  ! ========================
  ! Calculate all the forces
  ! ========================
    grad_dih_coef = -        c(2)*sin_dih  +        c(3)*cos_dih    &
                    - 2.0_wp*c(4)*sin_2dih + 2.0_wp*c(5)*cos_2dih   &
                    - 3.0_wp*c(6)*sin_3dih + 3.0_wp*c(7)*cos_3dih
    grad_dih_coef = AICG2P_K_DIHE * grad_dih_coef

  ! ------
  ! part 1
  ! ------
    P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
    P(1:3, 2) = ( u_t_2 * grad_dih_coef * (- A1 - A2) ) * aijk(1:3) + u_t_1 * grad_dih_coef * A3 * ajkl(1:3)
    P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
    P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)

  ! ------
  ! part 2
  ! ------
    if ( B1 < EPS ) then
      Q(1:3, 1:4) = 0.0_wp
    else
      tmp = u_dih * u_t_2 * B1
      Q(1:3, 1) = tmp * ( - inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
      Q(1:3, 3) = tmp * ( + inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
      Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
      Q(1:3, 4) = 0.0_wp
    end if

  ! ------
  ! part 3
  ! ------
    if ( C1 < EPS ) then
      R(1:3, 1:4) = 0.0_wp
    else
      tmp =  u_dih * u_t_1 * C1
      R(1:3, 2) = tmp * ( inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
      R(1:3, 4) = tmp * ( inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
      R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
      R(1:3, 1) = 0.0_wp
    end if

    ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)
    viri(1:3) = 0.0_wp
    do i = 1, 4
      viri(1:3) = viri(1:3) + ftmp(1:3,i)*cwork(1:3,i)
    end do

    return

  end subroutine dihed_flexible_cos2mod

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    local_dihed_cos2mod
  !> @brief        local dihedral energy called from energy_local_dihed
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine local_dihed_cos2mod(cwork, fc, phase, width, coef_dih_shift, &
                                 etmp, ftmp, viri)

! formal arguments
    real(wp),                intent(in)    :: cwork(:,:)
    real(wp),                intent(in)    :: fc, phase, width
    real(wp),                intent(in)    :: coef_dih_shift
    real(wp),                intent(inout) :: etmp
    real(wp),                intent(inout) :: ftmp(:,:)
    real(wp),                intent(inout) :: viri(:)

    ! local variables
    integer                  :: i
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
    real(wp)                 :: cos_tmin, sin_tmin
    real(wp)                 :: cos_d_dih, sin_d_dih, d_dih
    real(wp)                 :: grad_dih_coef, tmp
    real(wp)                 :: P(1:3,1:4)
    real(wp)                 :: Q(1:3,1:4)
    real(wp)                 :: R(1:3,1:4)

    dji(1:3) = cwork(1:3,1) - cwork(1:3,2)
    dkj(1:3) = cwork(1:3,2) - cwork(1:3,3)
    dkl(1:3) = cwork(1:3,4) - cwork(1:3,3)

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
    if ( u1 < 0.75_wp ) then
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
    if ( u2 < 0.75_wp ) then
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
    cos_dih = min(  1.0_wp-EPS, cos_dih )
    cos_dih = max( -1.0_wp+EPS, cos_dih )
!
    tmp     = aijk(1)*dkl(1) + aijk(2)*dkl(2) + aijk(3)*dkl(3)
    sin_dih = tmp * rkj * inv_raijkl
    sin_dih = min(  1.0_wp-EPS, sin_dih )
    sin_dih = max( -1.0_wp+EPS, sin_dih )

    cos_tmin = cos(phase)
    sin_tmin = sin(phase)
    cos_d_dih = cos_dih * cos_tmin + sin_dih * sin_tmin
    sin_d_dih = sin_dih * cos_tmin - cos_dih * sin_tmin

    cos_d_dih = min( 1.0_wp-EPS, cos_d_dih)
    cos_d_dih = max(-1.0_wp+EPS, cos_d_dih)
    sin_d_dih = min( 1.0_wp-EPS, sin_d_dih)
    sin_d_dih = max(-1.0_wp+EPS, sin_d_dih)

    if (cos_d_dih > 1.0E-1_wp) then
      d_dih = asin(sin_d_dih)
    else
      d_dih = sign(1.0_wp, sin_d_dih)*acos(cos_d_dih)
    endif

    tmp = d_dih / width
    u_dih  = fc * (coef_dih_shift + exp(-0.5_wp * tmp * tmp))

    ! ==============
    ! Compute energy
    ! ==============
    etmp  = u_dih * u_t_1 * u_t_2

    ! ========================
    ! Calculate all the forces
    ! ========================

    grad_dih_coef = - u_dih * tmp / width

    ! ------
    ! part 1
    ! ------
    P(1:3, 1) = u_t_2 * grad_dih_coef * A1 * aijk(1:3)
    P(1:3, 2) = (u_t_2*grad_dih_coef*(-A1-A2))*aijk(1:3) &
              + u_t_1*grad_dih_coef*A3*ajkl(1:3)
    P(1:3, 4) = - u_t_1 * grad_dih_coef * A4 * ajkl(1:3)
    P(1:3, 3) = - P(1:3, 1) - P(1:3, 2) - P(1:3, 4)

    ! ------
    ! part 2
    ! ------
    if ( B1 < EPS ) then
      Q(1:3, 1:4) = 0.0_wp
    else
      tmp = u_dih * u_t_2 * B1
      Q(1:3, 1) = tmp * ( - inv_dotpro_ijk * dkj(1:3) - inv_rji2 * dji(1:3))
      Q(1:3, 3) = tmp * ( + inv_dotpro_ijk * dji(1:3) + inv_rkj2 * dkj(1:3))
      Q(1:3, 2) = - Q(1:3, 1) - Q(1:3, 3)
      Q(1:3, 4) = 0.0_wp
    end if

    ! ------
    ! part 3
    ! ------
    if ( C1 < EPS ) then
      R(1:3, 1:4) = 0.0_wp
    else
      tmp =  u_dih * u_t_1 * C1
      R(1:3, 2) = tmp * ( inv_dotpro_jkl * dkl(1:3) - inv_rkj2 * dkj(1:3))
      R(1:3, 4) = tmp * ( inv_dotpro_jkl * dkj(1:3) - inv_rkl2 * dkl(1:3))
      R(1:3, 3) = - R(1:3, 2) - R(1:3, 4)
      R(1:3, 1) = 0.0_wp
    end if

    ftmp(1:3, 1:4) = P(1:3, 1:4) + Q(1:3, 1:4) + R(1:3, 1:4)

    viri(1:3) = 0.0_wp
    do i = 1, 4
      viri(1:3) = viri(1:3) + ftmp(1:3,i)*cwork(1:3,i)
    end do

    return

  end subroutine local_dihed_cos2mod

end module cg_energy_dihedrals_mod
