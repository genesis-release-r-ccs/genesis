!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_bonds_mod
!> @brief   calculate bond energy
!! @authors Jaewoon Jung (JJ), Yuji Sugia (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_bonds_mod

  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_bond
  public :: compute_energy_enm 
  ! FEP
  public :: compute_energy_bond_fep

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond
  !> @brief        calculate bond energy
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] ebond   : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond(domain, enefunc, coord, force, virial, ebond)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: ebond(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc
    real(wp)                 :: ebond_temp, work(3), viri(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id, omp_get_thread_num
    integer                  :: list(3), pbc_int, k1, k2, k3

    real(wp),        pointer :: fc(:,:), r0(:,:)
    real(wp),        pointer :: system_size(:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bond_pbc(:,:)


    call timer(TimerBond, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    system_size => domain%system_size

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    fc          => enefunc%bond_force_const
    r0          => enefunc%bond_dist_min
    bond_pbc    => enefunc%bond_pbc

    ! calculate bond energy
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         cc_frc, ebond_temp, work, list, viri, pbc_int, k1, k2, k3)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      ebond_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))
        pbc_int = bond_pbc(ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1

        ! bond energy: E=K[b-b0]^2
        !
        d12(1) = coord(1,i1,icel1) - coord(1,i2,icel2) &
               + system_size(1)*real(k1,wp)
        d12(2) = coord(2,i1,icel1) - coord(2,i2,icel2) &
               + system_size(2)*real(k2,wp)
        d12(3) = coord(3,i1,icel1) - coord(3,i2,icel2) &
               + system_size(3)*real(k3,wp)
        r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
        r_dif = r12 - r0(ix,i)
        ebond_temp = ebond_temp + fc(ix,i) * r_dif * r_dif

        ! gradient: dE/dX
        !
        cc_frc  = (2.0_wp * fc(ix,i) * r_dif) / r12
        work(1) = cc_frc * d12(1)
        work(2) = cc_frc * d12(2)
        work(3) = cc_frc * d12(3)

        ! virial
        !
        viri(1) = viri(1) + d12(1)*work(1)  
        viri(2) = viri(2) + d12(2)*work(2)  
        viri(3) = viri(3) + d12(3)*work(3)  
  
        ! store force: F=-dE/dX
        !
        force(1,i1,icel1,id+1) = force(1,i1,icel1,id+1) - work(1)
        force(2,i1,icel1,id+1) = force(2,i1,icel1,id+1) - work(2)
        force(3,i1,icel1,id+1) = force(3,i1,icel1,id+1) - work(3)
        force(1,i2,icel2,id+1) = force(1,i2,icel2,id+1) + work(1)
        force(2,i2,icel2,id+1) = force(2,i2,icel2,id+1) + work(2)
        force(3,i2,icel2,id+1) = force(3,i2,icel2,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      ebond(id+1) = ebond(id+1) + ebond_temp

    end do

    ! Bond energy for TIP2P case
    !
    if (enefunc%table%TIP2_bond_calc) then

      do i = id+1, ncell_local, nthread

        ebond_temp = 0.0_wp
        viri(1:3) = 0.0_wp

        do ix = 1, domain%num_water(i)

          list(1) = domain%water_list(1,ix,i)
          list(2) = domain%water_list(2,ix,i)

          d12(1) = coord(1,list(1),i) - coord(1,list(2),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(2),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(2),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%TIP2_bond
          ebond_temp = ebond_temp + enefunc%table%TIP2_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%TIP2_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)
          viri(2) = viri(2) + d12(2)*work(2)
          viri(3) = viri(3) + d12(3)*work(3)
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(2),i,1) = force(1,list(2),i,1) + work(1)
          force(2,list(2),i,1) = force(2,list(2),i,1) + work(2)
          force(3,list(2),i,1) = force(3,list(2),i,1) + work(3)

        end do

        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
        ebond(id+1) = ebond(id+1) + ebond_temp

      end do
    end if

    !O-H bond energy in water (for the case of no constraint)
    !
    if (enefunc%table%water_bond_calc) then

      do i = id+1, ncell_local, nthread

        ebond_temp = 0.0_wp
        viri(1:3) = 0.0_wp

        do ix = 1, domain%num_water(i)

          list(1) = domain%water_list(1,ix,i)
          list(2) = domain%water_list(2,ix,i)
          list(3) = domain%water_list(3,ix,i)

          ! first OH
          !
          d12(1) = coord(1,list(1),i) - coord(1,list(2),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(2),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(2),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)  
          viri(2) = viri(2) + d12(2)*work(2)  
          viri(3) = viri(3) + d12(3)*work(3)  
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(2),i,1) = force(1,list(2),i,1) + work(1)
          force(2,list(2),i,1) = force(2,list(2),i,1) + work(2)
          force(3,list(2),i,1) = force(3,list(2),i,1) + work(3)

          ! seoncd OH
          !
          d12(1) = coord(1,list(1),i) - coord(1,list(3),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(3),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(3),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)  
          viri(2) = viri(2) + d12(2)*work(2)  
          viri(3) = viri(3) + d12(3)*work(3)  
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(3),i,1) = force(1,list(3),i,1) + work(1)
          force(2,list(3),i,1) = force(2,list(3),i,1) + work(2)
          force(3,list(3),i,1) = force(3,list(3),i,1) + work(3)

          ! HH
          !
          if (enefunc%table%water_bond_calc_HH) then
            d12(1) = coord(1,list(2),i) - coord(1,list(3),i)
            d12(2) = coord(2,list(2),i) - coord(2,list(3),i)
            d12(3) = coord(3,list(2),i) - coord(3,list(3),i)
            r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
            r_dif = r12 - enefunc%table%HH_bond
            ebond_temp = ebond_temp + enefunc%table%HH_force*r_dif*r_dif
            cc_frc  = (2.0_wp*enefunc%table%HH_force*r_dif) / r12
            work(1) = cc_frc * d12(1)
            work(2) = cc_frc * d12(2)
            work(3) = cc_frc * d12(3)
            viri(1) = viri(1) + d12(1)*work(1)  
            viri(2) = viri(2) + d12(2)*work(2)  
            viri(3) = viri(3) + d12(3)*work(3)  
            force(1,list(2),i,1) = force(1,list(2),i,1) - work(1)
            force(2,list(2),i,1) = force(2,list(2),i,1) - work(2)
            force(3,list(2),i,1) = force(3,list(2),i,1) - work(3)
            force(1,list(3),i,1) = force(1,list(3),i,1) + work(1)
            force(2,list(3),i,1) = force(2,list(3),i,1) + work(2)
            force(3,list(3),i,1) = force(3,list(3),i,1) + work(3)
          end if

        end do

        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
        ebond(id+1) = ebond(id+1) + ebond_temp

      end do
    end if

    !$omp end parallel 
    
    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_enm 
  !> @brief        calculate enm energy
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] ebond   : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_enm(domain, enefunc, coord, force, virial, eenm)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: eenm(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc
    real(wp)                 :: eenm_temp, work(3), viri(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id, omp_get_thread_num
    integer                  :: list(3), pbc_int, k1, k2, k3

    real(wp),        pointer :: fc(:,:), r0(:,:)
    real(wp),        pointer :: system_size(:)
    integer,         pointer :: enm(:), enmlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    call timer(TimerBond, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    enm         => enefunc%num_enm
    enmlist     => enefunc%enm_list
    fc          => enefunc%enm_force_const
    r0          => enefunc%enm_dist_min

    ! calculate bond energy
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         cc_frc, eenm_temp, work, list, viri, pbc_int, k1, k2, k3)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      eenm_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do ix = 1, enm(i)

        icel1 = id_g2l(1,enmlist(1,ix,i))
        i1    = id_g2l(2,enmlist(1,ix,i))
        icel2 = id_g2l(1,enmlist(2,ix,i))
        i2    = id_g2l(2,enmlist(2,ix,i))

        ! bond energy: E=K[b-b0]^2
        !
        d12(1) = coord(1,i1,icel1) - coord(1,i2,icel2) 
        d12(2) = coord(2,i1,icel1) - coord(2,i2,icel2) 
        d12(3) = coord(3,i1,icel1) - coord(3,i2,icel2) 
        r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
        r_dif = r12 - r0(ix,i)
        eenm_temp = eenm_temp + fc(ix,i) * r_dif * r_dif

        ! gradient: dE/dX
        !
        cc_frc  = (2.0_wp * fc(ix,i) * r_dif) / r12
        work(1) = cc_frc * d12(1)
        work(2) = cc_frc * d12(2)
        work(3) = cc_frc * d12(3)

        ! virial
        !
        viri(1) = viri(1) + d12(1)*work(1)  
        viri(2) = viri(2) + d12(2)*work(2)  
        viri(3) = viri(3) + d12(3)*work(3)  
  
        ! store force: F=-dE/dX
        !
        force(1,i1,icel1,id+1) = force(1,i1,icel1,id+1) - work(1)
        force(2,i1,icel1,id+1) = force(2,i1,icel1,id+1) - work(2)
        force(3,i1,icel1,id+1) = force(3,i1,icel1,id+1) - work(3)
        force(1,i2,icel2,id+1) = force(1,i2,icel2,id+1) + work(1)
        force(2,i2,icel2,id+1) = force(2,i2,icel2,id+1) + work(2)
        force(3,i2,icel2,id+1) = force(3,i2,icel2,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      eenm(id+1) = eenm(id+1) + eenm_temp

    end do

    !$omp end parallel 
    
    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_enm
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond_fep
  !> @brief        calculate bond energy for FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial term of target systems
  !! @param[inout] ebond   : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond_fep(domain, enefunc, coord, force, virial, &
                                     ebond)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: ebond(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc
    real(wp)                 :: ebond_temp, work(3), viri(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id, omp_get_thread_num
    integer                  :: list(3), pbc_int, k1, k2, k3

    real(wp),        pointer :: fc(:,:), r0(:,:)
    real(wp),        pointer :: system_size(:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: bond_pbc(:,:)

    ! FEP
    integer                  :: fg1, fg2
    real(wp)                 :: lambbond

    call timer(TimerBond, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    system_size => domain%system_size

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    fc          => enefunc%bond_force_const
    r0          => enefunc%bond_dist_min
    bond_pbc    => enefunc%bond_pbc

    ! calculate bond energy
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         cc_frc, ebond_temp, work, list, viri, pbc_int, k1,   &
    !$omp         k2, k3, fg1, fg2, lambbond)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      ebond_temp = 0.0_wp
      viri(1:3) = 0.0_wp

      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))
        pbc_int = bond_pbc(ix,i)
        k3 = pbc_int / 9
        pbc_int = pbc_int - k3*9
        k2 = pbc_int / 3
        k1 = pbc_int - k2*3
        k1 = k1 - 1
        k2 = k2 - 1
        k3 = k3 - 1

        ! FEP: Determine lambbond
        fg1 = domain%fepgrp(i1,icel1)
        fg2 = domain%fepgrp(i2,icel2)
        lambbond = enefunc%table_bond_lambda(fg1,fg2)
        lambbond = lambbond &
          + real(enefunc%bond_singleB(ix,i), wp) &
          * (enefunc%lambbondB - enefunc%lambbondA)

        ! bond energy: E=K[b-b0]^2
        !
        d12(1) = coord(1,i1,icel1) - coord(1,i2,icel2) &
               + system_size(1)*real(k1,wp)
        d12(2) = coord(2,i1,icel1) - coord(2,i2,icel2) &
               + system_size(2)*real(k2,wp)
        d12(3) = coord(3,i1,icel1) - coord(3,i2,icel2) &
               + system_size(3)*real(k3,wp)
        r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
        r_dif = r12 - r0(ix,i)
        ebond_temp = ebond_temp + fc(ix,i) * r_dif * r_dif * lambbond

        ! gradient: dE/dX
        !
        cc_frc  = (2.0_wp * fc(ix,i) * r_dif) / r12 * lambbond
        work(1) = cc_frc * d12(1)
        work(2) = cc_frc * d12(2)
        work(3) = cc_frc * d12(3)

        ! virial
        !
        viri(1) = viri(1) + d12(1)*work(1)  
        viri(2) = viri(2) + d12(2)*work(2)  
        viri(3) = viri(3) + d12(3)*work(3)  
  
        ! store force: F=-dE/dX
        !
        force(1,i1,icel1,id+1) = force(1,i1,icel1,id+1) - work(1)
        force(2,i1,icel1,id+1) = force(2,i1,icel1,id+1) - work(2)
        force(3,i1,icel1,id+1) = force(3,i1,icel1,id+1) - work(3)
        force(1,i2,icel2,id+1) = force(1,i2,icel2,id+1) + work(1)
        force(2,i2,icel2,id+1) = force(2,i2,icel2,id+1) + work(2)
        force(3,i2,icel2,id+1) = force(3,i2,icel2,id+1) + work(3)

      end do

      virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
      virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
      virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
      ebond(id+1) = ebond(id+1) + ebond_temp

    end do

    ! Bond energy for TIP2P case
    !
    if (enefunc%table%TIP2_bond_calc) then

      do i = id+1, ncell_local, nthread

        ebond_temp = 0.0_wp
        viri(1:3) = 0.0_wp

        do ix = 1, domain%num_water(i)

          list(1) = domain%water_list(1,ix,i)
          list(2) = domain%water_list(2,ix,i)

          d12(1) = coord(1,list(1),i) - coord(1,list(2),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(2),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(2),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%TIP2_bond
          ebond_temp = ebond_temp + enefunc%table%TIP2_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%TIP2_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)
          viri(2) = viri(2) + d12(2)*work(2)
          viri(3) = viri(3) + d12(3)*work(3)
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(2),i,1) = force(1,list(2),i,1) + work(1)
          force(2,list(2),i,1) = force(2,list(2),i,1) + work(2)
          force(3,list(2),i,1) = force(3,list(2),i,1) + work(3)

        end do

        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
        ebond(id+1) = ebond(id+1) + ebond_temp

      end do
    end if

    !O-H bond energy in water (for the case of no constraint)
    !
    if (enefunc%table%water_bond_calc) then

      do i = id+1, ncell_local, nthread

        ebond_temp = 0.0_wp
        viri(1:3) = 0.0_wp

        do ix = 1, domain%num_water(i)

          list(1) = domain%water_list(1,ix,i)
          list(2) = domain%water_list(2,ix,i)
          list(3) = domain%water_list(3,ix,i)

          ! first OH
          !
          d12(1) = coord(1,list(1),i) - coord(1,list(2),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(2),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(2),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)  
          viri(2) = viri(2) + d12(2)*work(2)  
          viri(3) = viri(3) + d12(3)*work(3)  
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(2),i,1) = force(1,list(2),i,1) + work(1)
          force(2,list(2),i,1) = force(2,list(2),i,1) + work(2)
          force(3,list(2),i,1) = force(3,list(2),i,1) + work(3)

          ! seoncd OH
          !
          d12(1) = coord(1,list(1),i) - coord(1,list(3),i)
          d12(2) = coord(2,list(1),i) - coord(2,list(3),i)
          d12(3) = coord(3,list(1),i) - coord(3,list(3),i)
          r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1) = cc_frc * d12(1)
          work(2) = cc_frc * d12(2)
          work(3) = cc_frc * d12(3)
          viri(1) = viri(1) + d12(1)*work(1)  
          viri(2) = viri(2) + d12(2)*work(2)  
          viri(3) = viri(3) + d12(3)*work(3)  
          force(1,list(1),i,1) = force(1,list(1),i,1) - work(1)
          force(2,list(1),i,1) = force(2,list(1),i,1) - work(2)
          force(3,list(1),i,1) = force(3,list(1),i,1) - work(3)
          force(1,list(3),i,1) = force(1,list(3),i,1) + work(1)
          force(2,list(3),i,1) = force(2,list(3),i,1) + work(2)
          force(3,list(3),i,1) = force(3,list(3),i,1) + work(3)

          ! HH
          !
          if (enefunc%table%water_bond_calc_HH) then
            d12(1) = coord(1,list(2),i) - coord(1,list(3),i)
            d12(2) = coord(2,list(2),i) - coord(2,list(3),i)
            d12(3) = coord(3,list(2),i) - coord(3,list(3),i)
            r12   = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
            r_dif = r12 - enefunc%table%HH_bond
            ebond_temp = ebond_temp + enefunc%table%HH_force*r_dif*r_dif
            cc_frc  = (2.0_wp*enefunc%table%HH_force*r_dif) / r12
            work(1) = cc_frc * d12(1)
            work(2) = cc_frc * d12(2)
            work(3) = cc_frc * d12(3)
            viri(1) = viri(1) + d12(1)*work(1)  
            viri(2) = viri(2) + d12(2)*work(2)  
            viri(3) = viri(3) + d12(3)*work(3)  
            force(1,list(2),i,1) = force(1,list(2),i,1) - work(1)
            force(2,list(2),i,1) = force(2,list(2),i,1) - work(2)
            force(3,list(2),i,1) = force(3,list(2),i,1) - work(3)
            force(1,list(3),i,1) = force(1,list(3),i,1) + work(1)
            force(2,list(3),i,1) = force(2,list(3),i,1) + work(2)
            force(3,list(3),i,1) = force(3,list(3),i,1) + work(3)
          end if

        end do

        virial(1,1,id+1) = virial(1,1,id+1) - viri(1)
        virial(2,2,id+1) = virial(2,2,id+1) - viri(2)
        virial(3,3,id+1) = virial(3,3,id+1) - viri(3)
        ebond(id+1) = ebond(id+1) + ebond_temp

      end do
    end if

    !$omp end parallel 
    
    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond_fep

end module sp_energy_bonds_mod
