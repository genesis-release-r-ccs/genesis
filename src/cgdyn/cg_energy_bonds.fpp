!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_bonds_mod
!> @brief   Calculates bond energy.
!! @authors Jaewoon Jung (JJ), Yuji Sugia (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_bonds_mod

  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public :: compute_energy_bond

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond
  !> @brief        Calculates bond energy.
  !! @authors      JJ, YS
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] ebond   : bond energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond(domain, enefunc, coord, force, ebond, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: ebond(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, r_dif_2, cc_frc
    real(wp)                 :: ebond_temp, work(3)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: id

    real(wp),        pointer :: fc(:), r0(:)
    integer,         pointer :: bondlist(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerBond, TimerOn)

    ncell_local   => domain%num_cell_local
    id_g2l        => domain%id_g2l

    bondlist      => enefunc%bond_list
    fc            => enefunc%bond_force_const
    r0            => enefunc%bond_dist_min

    ! calculate bond energy
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         r_dif_2, cc_frc, ebond_temp, work)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    ebond_temp = 0.0_wp
    do i = id+1, enefunc%num_bondsq_domain, nthread

      i1    = id_g2l(bondlist(1,i))
      i2    = id_g2l(bondlist(2,i))

      ! bond energy: E=K[b-b0]^2
      !
      d12(1:3) = coord(i1,1:3) - coord(i2,1:3)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      r_dif = r12 - r0(i)
      ebond_temp = ebond_temp + fc(i) * r_dif * r_dif

      ! gradient: dE/dX
      !
      cc_frc  = (2.0_wp * fc(i) * r_dif) / r12
      work(1:3) = cc_frc * d12(1:3)

      ! store force: F=-dE/dX
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3)

      virial(1,1,id+1) = virial(1,1,id+1) - d12(1)*work(1)
      virial(2,2,id+1) = virial(2,2,id+1) - d12(2)*work(2)
      virial(3,3,id+1) = virial(3,3,id+1) - d12(3)*work(3)
    end do

    do i = id+1, enefunc%num_bond_domain, nthread

      ix = i + enefunc%num_bondsq_domain
      i1 = id_g2l(bondlist(1,ix))
      i2 = id_g2l(bondlist(2,ix))

      ! bond energy: E=K[b-b0]^2+100 K[b-b0]^4
      !
      d12(1:3) = coord(i1,1:3) - coord(i2,1:3)
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      r_dif = r12 - r0(ix)
      r_dif_2 = r_dif * r_dif
      ebond_temp = ebond_temp + fc(ix) * r_dif_2 * (1.0_wp+r_dif_2*100.0_wp)

      ! gradient: dE/dX
      !
      cc_frc = (2.0_wp * fc(ix) * r_dif * (1.0_wp+200.0_wp*r_dif_2)) / r12 
      work(1:3) = cc_frc * d12(1:3)

      ! store force: F=-dE/dX
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3)

      virial(1,1,id+1) = virial(1,1,id+1) - d12(1)*work(1)
      virial(2,2,id+1) = virial(2,2,id+1) - d12(2)*work(2)
      virial(3,3,id+1) = virial(3,3,id+1) - d12(3)*work(3)
    end do

    ebond(id+1) = ebond(id+1) + ebond_temp

    !$omp end parallel 
   
    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond

end module cg_energy_bonds_mod
