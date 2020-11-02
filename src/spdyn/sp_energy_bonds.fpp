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
  !! @param[inout] ebond   : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond(domain, enefunc, coord, force, ebond)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: ebond(nthread)

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc
    real(wp)                 :: ebond_temp, work(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id, omp_get_thread_num
    integer                  :: list(3)

    real(wp),        pointer :: fc(:,:), r0(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    call timer(TimerBond, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    fc          => enefunc%bond_force_const
    r0          => enefunc%bond_dist_min

    ! calculate bond energy
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         cc_frc, ebond_temp, work, list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      ebond_temp = 0.0_wp

      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))

        ! bond energy: E=K[b-b0]^2
        !
        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i2,icel2)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        r_dif = r12 - r0(ix,i)
        ebond_temp = ebond_temp + fc(ix,i) * r_dif * r_dif

        ! gradient: dE/dX
        !
        cc_frc  = (2.0_wp * fc(ix,i) * r_dif) / r12
        work(1:3) = cc_frc * d12(1:3)
  
        ! store force: F=-dE/dX
        !
        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1) + work(1:3)

      end do

      ebond(id+1) = ebond(id+1) + ebond_temp

    end do

    !O-H bond energy in water (for the case of no constraint)
    !
    if (enefunc%table%water_bond_calc) then

      do i = id+1, ncell_local, nthread

        ebond_temp = 0.0_wp

        do ix = 1, domain%num_water(i)

          list(1:3) = domain%water_list(1:3,ix,i)

          ! first OH
          !
          d12(1:3) = coord(1:3,list(1),i) - coord(1:3,list(2),i)
          r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1:3) = cc_frc * d12(1:3)
          force(1:3,list(1),i,1) = force(1:3,list(1),i,1) - work(1:3)
          force(1:3,list(2),i,1) = force(1:3,list(2),i,1) + work(1:3)

          ! seoncd OH
          !
          d12(1:3) = coord(1:3,list(1),i) - coord(1:3,list(3),i)
          r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
          r_dif = r12 - enefunc%table%OH_bond
          ebond_temp = ebond_temp + enefunc%table%OH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%OH_force*r_dif) / r12
          work(1:3) = cc_frc * d12(1:3)
          force(1:3,list(1),i,1) = force(1:3,list(1),i,1) - work(1:3)
          force(1:3,list(3),i,1) = force(1:3,list(3),i,1) + work(1:3)

          ! HH
          !
          d12(1:3) = coord(1:3,list(2),i) - coord(1:3,list(3),i)
          r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
          r_dif = r12 - enefunc%table%HH_bond
          ebond_temp = ebond_temp + enefunc%table%HH_force*r_dif*r_dif
          cc_frc  = (2.0_wp*enefunc%table%HH_force*r_dif) / r12
          work(1:3) = cc_frc * d12(1:3)
          force(1:3,list(2),i,1) = force(1:3,list(2),i,1) - work(1:3)
          force(1:3,list(3),i,1) = force(1:3,list(3),i,1) + work(1:3)

        end do

        ebond(id+1) = ebond(id+1) + ebond_temp

      end do
    end if

    !$omp end parallel 
    
    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond

end module sp_energy_bonds_mod
