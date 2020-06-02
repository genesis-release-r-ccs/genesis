!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_dihedrals_mod
!> @brief   calculate dihedral energy
!! @authors Chigusa Kobayashi (CK), Jaewoon Jung (JJ) , Takao Yoda (TY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_dihedrals_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use dihedral_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_dihed
  public :: compute_energy_dihed_localres
  public :: compute_energy_rb_dihed
  public :: compute_energy_improp
  public :: compute_energy_improp_cos
  public :: compute_energy_cmap

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed
  !> @brief        calculate dihedral energy
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
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
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
    integer                  :: omp_get_thread_num, list(4)

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    ndihe       => enefunc%num_dihedral
    dihelist    => enefunc%dihe_list
    fc          => enefunc%dihe_force_const
    nperiod     => enefunc%dihe_periodicity
    phase       => enefunc%dihe_phase

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp, cwork,        &
    !$omp         ix, work, edihe_temp, nrot,                                  &
    !$omp         icel1, i1, icel2, i2, icel3, i3, icel4, i4, list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      edihe_temp = 0.0_wp

      do ix = 1, ndihe(i)


        list(1) = dihelist(1,ix,i)
        list(2) = dihelist(2,ix,i)
        list(3) = dihelist(3,ix,i)
        list(4) = dihelist(4,ix,i)

        icel1 = id_g2l(1,list(1))
        i1    = id_g2l(2,list(1))
        icel2 = id_g2l(1,list(2))
        i2    = id_g2l(2,list(2))
        icel3 = id_g2l(1,list(3))
        i3    = id_g2l(2,list(3))
        icel4 = id_g2l(1,list(4))
        i4    = id_g2l(2,list(4))

        cwork(1,1) = coord(1,i1,icel1)
        cwork(2,1) = coord(2,i1,icel1)
        cwork(3,1) = coord(3,i1,icel1)
        cwork(1,2) = coord(1,i2,icel2)
        cwork(2,2) = coord(2,i2,icel2)
        cwork(3,2) = coord(3,i2,icel2)
        cwork(1,3) = coord(1,i3,icel3)
        cwork(2,3) = coord(2,i3,icel3)
        cwork(3,3) = coord(3,i3,icel3)
        cwork(1,4) = coord(1,i4,icel4)
        cwork(2,4) = coord(2,i4,icel4)
        cwork(3,4) = coord(3,i4,icel4)

        aindex(1) = 1
        aindex(2) = 2
        aindex(3) = 3
        aindex(4) = 4
!       call calculate_dihedral(aindex, cwork, cos_dih, sin_dih, grad, v)
        call calculate_dihedral_2(aindex, cwork, cos_dih, sin_dih, grad)
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
!       sinpha = sin(phase(ix, i))
        sinpha = sqrt(1.0_wp-cospha*cospha)
        edihe_temp = edihe_temp + fc(ix, i) *  &
                                  (1.0_wp + cospha*cosnt + sinnt*sinpha)

        grad_coef = fc(ix, i) * real(nrot,wp) *  &
                                    (cospha*sinnt - cosnt*sinpha)
        work(1) = grad_coef*grad(1)
        work(2) = grad_coef*grad(2)
        work(3) = grad_coef*grad(3)
        work(4) = grad_coef*grad(4)
        work(5) = grad_coef*grad(5)
        work(6) = grad_coef*grad(6)
        work(7) = grad_coef*grad(7)
        work(8) = grad_coef*grad(8)
        work(9) = grad_coef*grad(9)


        force(1,i1,icel1,id+1) = force(1,i1,icel1,id+1) - work(1)
        force(2,i1,icel1,id+1) = force(2,i1,icel1,id+1) - work(2)
        force(3,i1,icel1,id+1) = force(3,i1,icel1,id+1) - work(3)
        force(1,i2,icel2,id+1) = force(1,i2,icel2,id+1)  &
                               + work(1) - work(4)
        force(2,i2,icel2,id+1) = force(2,i2,icel2,id+1)  &
                               + work(2) - work(5)
        force(3,i2,icel2,id+1) = force(3,i2,icel2,id+1)  &
                               + work(3) - work(6)
        force(1,i3,icel3,id+1) = force(1,i3,icel3,id+1)  &
                               + work(4) + work(7)
        force(2,i3,icel3,id+1) = force(2,i3,icel3,id+1)  &
                               + work(5) + work(8)
        force(3,i3,icel3,id+1) = force(3,i3,icel3,id+1)  &
                               + work(6) + work(9)
        force(1,i4,icel4,id+1) = force(1,i4,icel4,id+1) - work(7)
        force(2,i4,icel4,id+1) = force(2,i4,icel4,id+1) - work(8)
        force(3,i4,icel4,id+1) = force(3,i4,icel4,id+1) - work(9)

!       work(1:9) = grad_coef*grad(1:9)


!       force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
!       force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)  &
!                                 + work(1:3) - work(4:6)
!       force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)  &
!                                 + work(4:6) + work(7:9)
!       force(1:3,i4,icel4,id+1) = force(1:3,i4,icel4,id+1) - work(7:9)

      end do

      edihe(id+1) = edihe(id+1) + edihe_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed_localres
  !> @brief        calculate dihedral energy
  !! @authors      CK, JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed_localres(domain, enefunc, coord, force, edihe)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
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
    integer                  :: omp_get_thread_num

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    integer(1),      pointer :: dkind(:,:)


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

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        cwork(1:3,1) = coord(1:3,i1,icel1)
        cwork(1:3,2) = coord(1:3,i2,icel2)
        cwork(1:3,3) = coord(1:3,i3,icel3)
        cwork(1:3,4) = coord(1:3,i4,icel4)

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


        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)  &
                                  + work(1:3) - work(4:6)
        force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)  &
                                  + work(4:6) + work(7:9)
        force(1:3,i4,icel4,id+1) = force(1:3,i4,icel4,id+1) - work(7:9)

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
  !> @brief        calculate dihedral energy
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
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
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
    integer                  :: omp_get_thread_num

    real(wp),        pointer :: fc(:,:,:)
    integer,         pointer :: ndihe(:), dihelist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)

    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    ndihe       => enefunc%num_rb_dihedral
    dihelist    => enefunc%rb_dihe_list
    fc          => enefunc%rb_dihe_c

    !$omp parallel default(shared)                                             &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         icn, coef, vtmp, cwork,  ix, work, edihe_temp,        &
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
        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel2 = id_g2l(1,dihelist(2,ix,i))
        i2    = id_g2l(2,dihelist(2,ix,i))
        icel3 = id_g2l(1,dihelist(3,ix,i))
        i3    = id_g2l(2,dihelist(3,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        cwork(1:3,1) = coord(1:3,i1,icel1)
        cwork(1:3,2) = coord(1:3,i2,icel2)
        cwork(1:3,3) = coord(1:3,i3,icel3)
        cwork(1:3,4) = coord(1:3,i4,icel4)

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

        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)  &
                                  + work(1:3) - work(4:6)
        force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)  &
                                  + work(4:6) + work(7:9)
        force(1:3,i4,icel4,id+1) = force(1:3,i4,icel4,id+1) - work(7:9)

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
  !> @brief        calculate improper energy
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
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
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
    integer                  :: id, omp_get_thread_num

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: nimp(:), imprlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


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

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel2 = id_g2l(1,imprlist(2,ix,i))
        i2    = id_g2l(2,imprlist(2,ix,i))
        icel3 = id_g2l(1,imprlist(3,ix,i))
        i3    = id_g2l(2,imprlist(3,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        cwork(1:3,1) = coord(1:3,i1,icel1)
        cwork(1:3,2) = coord(1:3,i2,icel2)
        cwork(1:3,3) = coord(1:3,i3,icel3)
        cwork(1:3,4) = coord(1:3,i4,icel4)

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


        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)  &
                                  + work(1:3) - work(4:6)
        force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)  &
                                  + work(4:6) + work(7:9)
        force(1:3,i4,icel4,id+1) = force(1:3,i4,icel4,id+1) - work(7:9)

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
  !> @brief        calculate improper energy
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
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
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
    integer                  :: omp_get_thread_num

    real(wp),        pointer :: fc(:,:), phase(:,:)
    integer,         pointer :: nimp(:), imprlist(:,:,:)
    integer,         pointer :: nperiod(:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


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
    !$omp         ix, work, eimpr_temp, nrot,                                  &
    !$omp         icel1, i1, icel2, i2, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      eimpr_temp = 0.0_wp

      do ix = 1, nimp(i)

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel2 = id_g2l(1,imprlist(2,ix,i))
        i2    = id_g2l(2,imprlist(2,ix,i))
        icel3 = id_g2l(1,imprlist(3,ix,i))
        i3    = id_g2l(2,imprlist(3,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        cwork(1:3,1) = coord(1:3,i1,icel1)
        cwork(1:3,2) = coord(1:3,i2,icel2)
        cwork(1:3,3) = coord(1:3,i3,icel3)
        cwork(1:3,4) = coord(1:3,i4,icel4)

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


        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1)  &
                                  + work(1:3) - work(4:6)
        force(1:3,i3,icel3,id+1) = force(1:3,i3,icel3,id+1)  &
                                  + work(4:6) + work(7:9)
        force(1:3,i4,icel4,id+1) = force(1:3,i4,icel4,id+1) - work(7:9)

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
  !  Subroutine    compute_energy_cmap
  !> @brief        calculate cmap energy
  !! @authors      CK,  JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eimprop : improper dihedral energy of target systems
  !! @note         A.D.MacKerell et al., J.Comput.Chem., 25, 1400-1415 (2004).
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_cmap(domain, enefunc, coord, force, ecmap)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: ecmap(nthread)

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, dihed
    real(wp)                 :: delta, inv_delta, ctmp, vtmp
    integer                  :: aindex(1:8)
    integer                  :: k, ix, icel(1:8), iatom(1:8)
    integer                  :: igrid(1:2), ngrid0, i, j, m, igr, ip, itype
    real(wp)                 :: work(1:9,1:2), grad_coef(1:2), gradtot(1:2)
    real(wp)                 :: dgrid(1:2), cwork(1:3,1:8)
    real(wp)                 :: gradphi(1:9), vphi(1:3,1:3),f(3,8)
    real(wp)                 :: gradpsi(1:9), vpsi(1:3,1:3)
    real(wp)                 :: grid_power(1:2,1:4), dgrid_power(1:2,1:4)
    real(wp)                 :: ecmap_temp
    integer                  :: id, omp_get_thread_num
    integer                  :: ic1, ic2, ia1, ia2

    real(wp),        pointer :: coef(:,:,:,:,:)
    integer,         pointer :: ncmap(:), cmaplist(:,:,:), cmaptype(:,:)
    integer,         pointer :: resol(:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    call timer(TimerDihedral, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    ncmap       => enefunc%num_cmap
    cmaplist    => enefunc%cmap_list
    cmaptype    => enefunc%cmap_type
    resol       => enefunc%cmap_resolution
    coef        => enefunc%cmap_coef

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, icel, iatom, aindex, ngrid0, delta, inv_delta,    &
    !$omp         ip, itype, k, j, f,    &
    !$omp         cos_dih, sin_dih, dihed, igr, dgrid, igrid, vtmp,       &
    !$omp         grid_power, dgrid_power, work, m,  ctmp, grad_coef,     &
    !$omp         gradpsi, gradphi, vpsi, vphi, gradtot, cwork,           &
    !$omp         ic1, ic2, ia1, ia2, ecmap_temp)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      ecmap_temp = 0.0_wp

      do ix = 1, ncmap(i)

        do k = 1, 8
          icel(k)  = id_g2l(1,cmaplist(k,ix,i))
          iatom(k) = id_g2l(2,cmaplist(k,ix,i))

          aindex(k) = k

          cwork(1:3,k) = coord(1:3,iatom(k),icel(k))
        end do

        itype  = cmaptype(ix,i)
        ngrid0  = resol(itype)
        delta     = 360.0_wp/real(ngrid0,wp)
        inv_delta = 1.0_wp/delta

      ! dih1 is backbone phi
        call calculate_dihedral(aindex(1:4), cwork, cos_dih, sin_dih,       &
                               gradphi, vphi)

        gradphi(1:9) = gradphi(1:9)/RAD
        if (abs(cos_dih) > 1.0E-1_wp) then
          dihed = asin(sin_dih)/RAD
          if (cos_dih < 0.0_wp) then
            if (dihed > 0.0_wp) then
              dihed=180.0_wp-dihed
            else
              dihed=-180.0_wp-dihed
            endif
          endif
        else
          dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
        endif

        if (dihed < -180.0_wp) then
          dihed = dihed + 360.0_wp
        else if (dihed > 180.0_wp) then
          dihed = dihed - 360.0_wp
        endif

        igr      = int((dihed+180.0_wp)*inv_delta)
        dgrid(1) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
        igrid(1) = igr + 1

        ! dih2 is backbone psi
        call calculate_dihedral(aindex(5:8), cwork, cos_dih, sin_dih,          &
                                gradpsi, vpsi)

        gradpsi(1:9) = gradpsi(1:9)/RAD
       
        if (abs(cos_dih) > 1.0E-1_wp) then
          dihed = asin(sin_dih)/RAD
          if (cos_dih < 0.0_wp) then
            if (dihed > 0.0_wp) then
              dihed=180.0_wp-dihed
            else
              dihed=-180.0_wp-dihed
            endif
          endif
        else
          dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
        endif

        if (dihed < -180.0_wp) then
          dihed = dihed + 360.0_wp
        else if (dihed > 180.0_wp) then
          dihed = dihed - 360.0_wp
        endif

        igr      = int((dihed+180.0_wp)*inv_delta)
        dgrid(2) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
        igrid(2) = igr + 1

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
          do k = 1, 4
            ctmp = coef(k,j,igrid(2),igrid(1),itype)

            ! cmap energy
            !
            ecmap_temp = ecmap_temp + grid_power(2,k)*grid_power(1,j)*ctmp

            ! for gradient
            !
            grad_coef(1) = -dgrid_power(1,j)*grid_power(2,k)*ctmp
            grad_coef(2) = -dgrid_power(2,k)*grid_power(1,j)*ctmp

            gradtot(1:2) = gradtot(1:2) + grad_coef(1:2)
          end do
        end do
        work(1:9,1) = gradtot(1)*gradphi(1:9)
        work(1:9,2) = gradtot(2)*gradpsi(1:9)

        f(1:3,1) = - work(1:3,1)
        f(1:3,2) =   work(1:3,1) - work(4:6,1)
        f(1:3,3) =   work(4:6,1) + work(7:9,1)
        f(1:3,4) = - work(7:9,1)
        f(1:3,5) = - work(1:3,2)
        f(1:3,6) =   work(1:3,2) - work(4:6,2)
        f(1:3,7) =   work(4:6,2) + work(7:9,2)
        f(1:3,8) = - work(7:9,2)

        do k = 1, 8
          ic1 = icel(k)
          ia1 = iatom(k)
          force(1:3,ia1,ic1,id+1) = force(1:3,ia1,ic1,id+1) + f(1:3,k)
        end do

      end do

      ecmap(id+1) = ecmap(id+1) + ecmap_temp

    end do

    !$omp end parallel

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cmap

  subroutine calculate_dihedral_2(aindex, coord, cos_dih, sin_dih, grad)

    ! formal arguments
    integer,   intent(in)    :: aindex(:)
    real(wp),  intent(in)    :: coord(:,:)
    real(wp),  intent(inout) :: cos_dih
    real(wp),  intent(inout) :: sin_dih
    real(wp),  intent(inout) :: grad(:)

    ! local variable
    integer                  :: i, j
    real(wp)                 :: dij(1:3), djk(1:3), dlk(1:3)
    real(wp)                 :: aijk(1:3), ajkl(1:3)
    real(wp)                 :: tmp(1:4)
    real(wp)                 :: raijk2, rajkl2, vtmp
    real(wp)                 :: inv_raijk2, inv_rajkl2, inv_raijkl
    real(wp)                 :: rjk, inv_rjk, dotpro_ijk, dotpro_jkl

    dij(1) = coord(1,aindex(1)) - coord(1,aindex(2))
    dij(2) = coord(2,aindex(1)) - coord(2,aindex(2))
    dij(3) = coord(3,aindex(1)) - coord(3,aindex(2))
    djk(1) = coord(1,aindex(2)) - coord(1,aindex(3))
    djk(2) = coord(2,aindex(2)) - coord(2,aindex(3))
    djk(3) = coord(3,aindex(2)) - coord(3,aindex(3))
    dlk(1) = coord(1,aindex(4)) - coord(1,aindex(3))
    dlk(2) = coord(2,aindex(4)) - coord(2,aindex(3))
    dlk(3) = coord(3,aindex(4)) - coord(3,aindex(3))

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

    grad(1) =  tmp(1)*aijk(1)
    grad(2) =  tmp(1)*aijk(2)
    grad(3) =  tmp(1)*aijk(3)
    grad(4) = -tmp(3)*aijk(1) + tmp(4)*ajkl(1)
    grad(5) = -tmp(3)*aijk(2) + tmp(4)*ajkl(2)
    grad(6) = -tmp(3)*aijk(3) + tmp(4)*ajkl(3)
    grad(7) =                 - tmp(2)*ajkl(1)
    grad(8) =                 - tmp(2)*ajkl(2)
    grad(9) =                 - tmp(2)*ajkl(3)

    return

  end subroutine calculate_dihedral_2

end module sp_energy_dihedrals_mod
