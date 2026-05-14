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
  public :: compute_energy_dihed_localres
  public :: compute_energy_rb_dihed
  public :: compute_energy_improp
  public :: compute_energy_improp_cos

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

  subroutine compute_energy_dihed(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, k, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, edihe_temp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: ndihe(:), ndihe_local(:), ndihe_flexible(:)
    integer,         pointer :: dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    ndihe          => enefunc%num_dihedral
    ndihe_local    => enefunc%num_dihe_local
    ndihe_flexible => enefunc%num_dihe_flexible
    dihelist       => enefunc%dihe_list
    fc             => enefunc%dihe_force_const
    nperiod        => enefunc%dihe_periodicity
    phase          => enefunc%dihe_phase

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, cwork,        &
    !$omp         ix, work, edihe_temp, nrot,                                  &
    !$omp         icel1, i1, icel2, i2, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = ndihe_local(i)+ndihe_flexible(i)+1, ndihe(i)

        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)
        cosnt = 1.0_wp
        sinnt = 0.0_wp
        krot = 0
        nrot = nperiod(ix,i)
        if (enefunc%notation_14types > 0) &
        nrot = mod(nperiod(ix,i), enefunc%notation_14types)
        do while (krot < nrot)
          tmp   = cosnt*cos_dih - sinnt*sin_dih
          sinnt = sinnt*cos_dih + cosnt*sin_dih
          cosnt = tmp
          krot = krot+1
        end do

        cospha = cos(phase(ix, i))
        sinpha = sin(phase(ix, i))
        edihe_temp = edihe_temp + fc(ix, i) *  &
                                  (1.0_wp + cospha*cosnt + sinnt*sinpha)

        grad_coef = fc(ix, i) * real(nrot,wp) *  &
                                    (cospha*sinnt - cosnt*sinpha)
        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_dihed
  !> @brief        calculate local dihedral energy
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_dihed(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, k, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, cos_dih, sin_dih
    real(wp)                 :: grad_coef, tmp, edihe_temp, dihe
    real(wp)                 :: cosdif, sindif, diffphi, theta
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:), phase(:,:), width(:,:)
    integer,         pointer :: ndihe(:), ndihe_local(:), ndihe_flexible(:)
    integer,         pointer :: dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    ndihe          => enefunc%num_dihedral
    ndihe_local    => enefunc%num_dihe_local
    ndihe_flexible => enefunc%num_dihe_flexible
    dihelist       => enefunc%dihe_list
    fc             => enefunc%dihe_force_const
    phase          => enefunc%dihe_phase
    width          => enefunc%dihe_width

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         theta, cospha, sinpha, cosdif, sindif, diffphi, tmp, cwork,  &
    !$omp         ix, dihe, work, edihe_temp, icel1, i1, icel2, i2, icel3, i3, &
    !$omp         icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = ndihe_flexible(i)+1, ndihe_local(i)+ndihe_flexible(i)

        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)
        fc_temp = fc(ix,i)
        phase_temp = phase(ix,i)
        width_temp = width(ix,i)

        if (func(ix,i) == 21) then

          call local_dihed(cwork, fc_temp, phase_temp, width_temp, &
                           etmp, work)
          edihe_temp = edihe_temp + etmp

        else if (func(ix,i) == 41) then

          call local_dihed_cos2mod(aindex, cwork, phase(ix,i), width(ix,i), &
                                   etmp, work)

        end if        


        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

        theta = phase(ix,i)
        cospha = cos(theta)
        sinpha = sin(theta)

        cosdif = cos_dih*cospha + sin_dih*sinpha
        sindif = sin_dih*cospha - cos_dih*sinpha

        if (cosdif > 0.1_wp) then
          diffphi = asin(sindif)
        else
          diffphi = sign(1.0_wp,sindif)*acos(cosdif)
        end if

        dihe = diffphi / width(ix,i)
        tmp  = fc(ix,i)*exp(-0.5_wp*dihe*dihe)

        edihe_temp = edihe_temp + tmp

        grad_coef = (dihe*tmp)/width(ix,i)
        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

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

  subroutine compute_energy_flexible_dihed(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, k, id, krot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4), dtype
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, edihe_temp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cos_2dih, sin_2dih
    real(wp)                 :: cos_3dih, sin_3dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9), c(1:7)

    integer,         pointer :: ndihe_flexible(:)
    integer,         pointer :: dihelist(:,:,:)
    integer,         pointer :: dihe_type(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)
    real(wp),        pointer :: coef(:,:), ener_corr(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    ndihe_flexible => enefunc%num_dihe_flexible
    dihelist       => enefunc%dihe_list
    coef           => enefunc%diheflex_coef
    ener_corr      => enefunc%diheflex_ener_corr
    dihe_type      => enefunc%dihe_kind

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, cwork,        &
    !$omp         ix, work, edihe_temp, cos_2dih, sin_2dih, cos_3dih,          &
    !$omp         sin_3dih, c, dtype, icel1, i1, icel2, i2, icel3, i3, icel4,  &
    !$omp         i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = 1, ndihe_flexible(i)

        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

        cos_2dih =  2.0_wp*cos_dih*cos_dih - 1.0_wp
        sin_2dih =  2.0_wp*cos_dih*sin_dih
        cos_3dih =  4.0_wp*cos_dih*cos_dih*cos_dih - 3.0_wp*cos_dih
        sin_3dih = -4.0_wp*sin_dih*sin_dih*sin_dih + 3.0_wp*sin_dih

        dtype = dihe_type(ix,i)
        c(1:7)   = coef(1:7,dtype)

        tmp = c(1)                             &
            + c(2)*cos_dih  + c(3)*sin_dih     &
            + c(4)*cos_2dih + c(5)*sin_2dih    &
            + c(6)*cos_3dih + c(7)*sin_3dih
        edihe_temp = edihe_temp + AICG2P_K_DIHE*(tmp-ener_corr(dtype))
       
        grad_coef = -       c(2)*sin_dih  &
                  +         c(3)*cos_dih  &
                  -  2.0_wp*c(4)*sin_2dih &
                  +  2.0_wp*c(5)*cos_2dih &
                  -  3.0_wp*c(6)*sin_3dih &
                  +  3.0_wp*c(7)*cos_3dih
 
        grad_coef = -AICG2P_K_DIHE * grad_coef
        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_flexible_dihed


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed_localres
  !> @brief        Calculates dihedral energy.
  !! @authors      CK, JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Chem. Comput., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed_localres(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, edihe_temp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosdif, sindif, diffphi
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)
    integer,         pointer :: dkind(:,:)


    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    ndihe       => enefunc%num_dihedral
    dihelist    => enefunc%dihe_list
    fc          => enefunc%dihe_force_const
    nperiod     => enefunc%dihe_periodicity
    phase       => enefunc%dihe_phase
    dkind       => enefunc%dihe_kind

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,      &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, cwork,        &
    !$omp         sindif, cosdif, diffphi, ix, work, edihe_temp, nrot,         &
    !$omp         icel1, i1, icel2, i2, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = 1, ndihe(i)

        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4

        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)
        cospha = cos(phase(ix, i))
        sinpha = sin(phase(ix, i))

        if (dkind(ix,i) == 0) then
          cosnt = 1.0_wp
          sinnt = 0.0_wp
          krot = 0
          nrot = nperiod(ix,i)
          if (enefunc%notation_14types > 0) &
          nrot = mod(nperiod(ix,i), enefunc%notation_14types)
          do while (krot < nrot)
            tmp   = cosnt*cos_dih - sinnt*sin_dih
            sinnt = sinnt*cos_dih + cosnt*sin_dih
            cosnt = tmp
            krot = krot+1
          end do

          edihe_temp = edihe_temp + fc(ix, i) *  &
                                    (1.0_wp + cospha*cosnt + sinnt*sinpha)

          grad_coef = fc(ix, i) * real(nrot,wp) *  &
                                      (cospha*sinnt - cosnt*sinpha)
        else if (dkind(ix,i) == 1) then

          cosdif = cos_dih*cospha + sin_dih*sinpha
          sindif = cos_dih*sinpha - sin_dih*cospha

          if (cosdif > 1.0E-1_wp) then
            diffphi = asin(sindif)
          else
            diffphi = sign(1.0_wp,sindif)*acos(cosdif)
          endif
          edihe_temp = edihe_temp + fc(ix, i)*diffphi*diffphi
          grad_coef = 2.0_wp * fc(ix, i)*diffphi

        endif

        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed_localres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_rb_dihed
  !> @brief        Calculates dihedral energy.
  !! @authors      CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_rb_dihed(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, k, id, icn
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: grad_coef, edihe_temp, coef
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)

    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    ndihe       => enefunc%num_rb_dihedral
    dihelist    => enefunc%rb_dihe_list
    fc          => enefunc%rb_dihe_c

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         icn, coef, vtmp, cwork,  ix, work, edihe_temp, icel1, i1,    &
    !$omp         icel2, i2, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif


    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp
      do ix = 1, ndihe(i)
        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

!       psi = phi - pi
        cos_dih = -cos_dih
        sin_dih = -sin_dih
        coef = 0.0_wp
        do icn = 1, 6
          coef = coef + real(icn-1,wp) * fc(icn, ix, i) * cos_dih**(icn-2)
          edihe_temp = edihe_temp + fc(icn, ix, i) * cos_dih**(icn-1)
        end do
       
        grad_coef = sin_dih * coef
        work(1:9) = grad_coef*grad(1:9)

        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do
    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_rb_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp
  !> @brief        Calculates improper energy.
  !! @authors      CK, JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eimprop : improper dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp(domain, enefunc, coord, force, eimprop)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eimprop(nthread)

    ! local variables
    integer                  :: i, ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, eimp_temp
    real(wp)                 :: cosdif, sindif, diffphi
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)
    integer                  :: id

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: nimp(:), imprlist(:,:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    nimp        => enefunc%num_improper
    imprlist    => enefunc%impr_list
    fc          => enefunc%impr_force_const
    phase       => enefunc%impr_phase

    !$omp parallel default(shared)                                             &
    !$omp private(i, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,         &
    !$omp         cospha, sinpha, cosnt, sinnt, tmp, vtmp, cwork,              &
    !$omp         sindif, cosdif, diffphi, ix, work, eimp_temp,                &
    !$omp         icel1, i1, icel2, i2, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      eimp_temp = 0.0_wp

      do ix = 1, nimp(i)

        i1    = id_g2l(imprlist(1,ix,i))
        i2    = id_g2l(imprlist(2,ix,i))
        i3    = id_g2l(imprlist(3,ix,i))
        i4    = id_g2l(imprlist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4

        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)
        cospha = cos(phase(ix, i))
        sinpha = sin(phase(ix, i))

        cosdif = cos_dih*cospha + sin_dih*sinpha
        sindif = cos_dih*sinpha - sin_dih*cospha

        if (cosdif > 1.0E-1_wp) then
          diffphi = asin(sindif)
        else
          diffphi = sign(1.0_wp,sindif)*acos(cosdif)
        endif
        eimp_temp = eimp_temp + fc(ix, i)*diffphi*diffphi
        grad_coef = 2.0_wp*fc(ix, i)*diffphi

        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      eimprop(id+1) = eimprop(id+1) + eimp_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_improp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp_cos
  !> @brief        Calculates improper energy.
  !! @authors      CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eimprop : improper dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp_cos(domain, enefunc, coord, force, eimprop)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eimprop(nthread)

    ! local variables
    integer                  :: i, j, k, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp, eimpr_temp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: nimp(:), imprlist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    nimp        => enefunc%num_improper
    imprlist    => enefunc%impr_list
    fc          => enefunc%impr_force_const
    nperiod     => enefunc%impr_periodicity
    phase       => enefunc%impr_phase

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, cwork,        &
    !$omp         ix, work, eimpr_temp, nrot, icel1, i1, icel2, i2, icel3, i3, &
    !$omp         icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      eimpr_temp = 0.0_wp

      do ix = 1, nimp(i)

        i1    = id_g2l(imprlist(1,ix,i))
        i2    = id_g2l(imprlist(2,ix,i))
        i3    = id_g2l(imprlist(3,ix,i))
        i4    = id_g2l(imprlist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4

        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)
        cosnt = 1.0_wp
        sinnt = 0.0_wp
        krot = 0
        nrot = nperiod(ix,i)
        if (enefunc%notation_14types > 0) &
        nrot = mod(nperiod(ix,i), enefunc%notation_14types)

        do while (krot < nrot)
          tmp   = cosnt*cos_dih - sinnt*sin_dih
          sinnt = sinnt*cos_dih + cosnt*sin_dih
          cosnt = tmp
          krot = krot+1
        end do

        cospha = cos(phase(ix, i))
        sinpha = sin(phase(ix, i))
        eimpr_temp = eimpr_temp + fc(ix, i) *  &
                                  (1.0_wp + cospha*cosnt + sinnt*sinpha)

        grad_coef = fc(ix, i) * real(nrot,wp) *  &
                                    (cospha*sinnt - cosnt*sinpha)
        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      eimprop(id+1) = eimprop(id+1) + eimpr_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return


    return

  end subroutine compute_energy_improp_cos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_local_dihed
  !> @brief        calculate local dihedral energy
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_local_dihed(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: edihe(nthread)

    ! local variables
    integer                  :: i, j, k, id, krot, nrot
    integer                  :: ix, icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, cos_dih, sin_dih
    real(wp)                 :: grad_coef, tmp, edihe_temp, dihe
    real(wp)                 :: cosdif, sindif, diffphi, theta
    real(wp)                 :: grad(1:9), v(1:3,1:3), cwork(1:3,1:4)
    real(wp)                 :: work(1:9)

    real(wp),        pointer :: fc(:,:), phase(:,:), width(:,:)
    integer,         pointer :: ndihe(:), ndihe_local(:), ndihe_flexible(:)
    integer,         pointer :: dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerDihedral, TimerOn)

    ncell_local    => domain%num_cell_local
    id_g2l         => domain%id_g2l

    ndihe          => enefunc%num_dihedral
    ndihe_local    => enefunc%num_dihe_local
    ndihe_flexible => enefunc%num_dihe_flexible
    dihelist       => enefunc%dihe_list
    fc             => enefunc%dihe_force_const
    phase          => enefunc%dihe_phase
    width          => enefunc%dihe_width

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         theta, cospha, sinpha, cosdif, sindif, diffphi, tmp, cwork,  &
    !$omp         ix, dihe, work, edihe_temp, icel1, i1, icel2, i2, icel3, i3, &
    !$omp         icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = ndihe_flexible(i)+1, ndihe_local(i)+ndihe_flexible(i)

        i1    = id_g2l(dihelist(1,ix,i))
        i2    = id_g2l(dihelist(2,ix,i))
        i3    = id_g2l(dihelist(3,ix,i))
        i4    = id_g2l(dihelist(4,ix,i))

        cwork(1:3,1) = coord(i1,1:3)
        cwork(1:3,2) = coord(i2,1:3)
        cwork(1:3,3) = coord(i3,1:3)
        cwork(1:3,4) = coord(i4,1:3)
        fc_temp = fc(ix,i)
        phase_temp = phase(ix,i)
        width_temp = width(ix,i)

        if (func(ix,i) == 21) then

          call local_dihed(cwork, fc_temp, phase_temp, width_temp, &
                           etmp, work)
          edihe_temp = edihe_temp + etmp

        else if (func(ix,i) == 41) then

          call local_dihed_cos2mod(aindex, cwork, phase(ix,i), width(ix,i), &
                                   etmp, work)

        end if


        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
        call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)

        theta = phase(ix,i)
        cospha = cos(theta)
        sinpha = sin(theta)

        cosdif = cos_dih*cospha + sin_dih*sinpha
        sindif = sin_dih*cospha - cos_dih*sinpha

        if (cosdif > 0.1_wp) then
          diffphi = asin(sindif)
        else
          diffphi = sign(1.0_wp,sindif)*acos(cosdif)
        end if

        dihe = diffphi / width(ix,i)
        tmp  = fc(ix,i)*exp(-0.5_wp*dihe*dihe)

        edihe_temp = edihe_temp + tmp

        grad_coef = (dihe*tmp)/width(ix,i)
        work(1:9) = grad_coef*grad(1:9)


        force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
        force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3) - work(4:6)
        force(i3,1:3,id+1) = force(i3,1:3,id+1) + work(4:6) + work(7:9)
        force(i4,1:3,id+1) = force(i4,1:3,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return


end module cg_energy_dihedrals_mod
