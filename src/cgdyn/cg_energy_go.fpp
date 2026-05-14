!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_contacts_mod
!> @brief   Calculates native contact energy.
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_go_mod

  use cg_pairlist_str_mod
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
  public :: compute_energy_contact_126
  public :: compute_energy_contact_1210

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126
  !> @brief        Calculates contact energy.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] econtact: contact energy of target systems
  !! @param[inout] enoncontact: non-contact energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126(domain, enefunc, coord, force, &
                                        econtact, enoncontact)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: econtact(nthread)
    real(dp),                intent(inout) :: enoncontact(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, term_lj12, term_lj6
    real(wp)                 :: econt, encont, work(3), coef, cutoff2
    integer                  :: list(2)
    integer                  :: i, ix, icel1, icel2, i1, i2
    integer                  :: iatmcls, jatmcls
    integer                  :: id

    real(wp),        pointer :: contact_lj12(:), contact_lj6(:)
    real(wp),        pointer :: nonb_lj12(:,:)
    integer,         pointer :: contactlist(:,:)
    integer,         pointer :: atmcls(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerContact, TimerOn)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    atmcls      => domain%atom_cls_no

    contactlist => enefunc%contact_list
    contact_lj12=> enefunc%contact_lj12
    contact_lj6 => enefunc%contact_lj6
    nonb_lj12   => enefunc%nonb_lj12

    cutoff2     = enefunc%cutoffdist * enefunc%cutoffdist

    ! calculate bond energy
    !
    !$omp parallel default(shared) &
    !$omp private(id, i, ix, icel1, i1, icel2, i2, lj6, lj12, dij, rij2, &
    !$omp         inv_rij2, inv_rij6, inv_rij12, term_lj12, term_lj6,    &
    !$omp         econt, encont, iatmcls, jatmcls, coef, work, list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    econt  = 0.0_wp
    encont = 0.0_wp

    do i = id+1, enefunc%num_contact_domain, nthread

      list(1:2) = contactlist(1:2,i)
      i1    = id_g2l(list(1))
      i2    = id_g2l(list(2))
      dij(1:3)  = coord(i1,1:3) - coord(i2,1:3)

      ! contact energy
      !
      lj6       = contact_lj6(i)
      lj12      = contact_lj12(i)
      rij2      = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6

      term_lj12 = lj12 * inv_rij12
      term_lj6  = lj6 * inv_rij6
      econt     = econt + term_lj12 - term_lj6
      coef      = 12.0_wp*term_lj12 - 6.0_wp*term_lj6

      ! noncontact energy
      !
      if (rij2 < cutoff2) then
        iatmcls   = atmcls(i1)
        jatmcls   = atmcls(i2)
        lj12      = nonb_lj12(iatmcls, jatmcls)
        term_lj12 = lj12 * inv_rij12
        encont    = encont - term_lj12
        coef      = coef - 12.0_wp*term_lj12
      end if

      ! gradient: dE/dX
      !
      coef      = - inv_rij2 * coef
      work(1:3) = coef * dij(1:3)

      ! store force: F=-dE/dX
      !
      force(i1,1:3,id+1) = force(i1,1:3,id+1) - work(1:3)
      force(i2,1:3,id+1) = force(i2,1:3,id+1) + work(1:3)

    end do

    econtact(id+1) = econtact(id+1) + econt
    enoncontact(id+1) = enoncontact(id+1) + encont

    !$omp end parallel

    call timer(TimerContact, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_126

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_1210
  !> @brief        Calculates contact energy.
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] econtact: contact energy of target systems
  !! @param[inout] enoncontact: non-contact energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_1210(domain, enefunc, coord, force, &
                                         econtact, enoncontact, virial)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: econtact(nthread)
    real(dp),                intent(inout) :: enoncontact(nthread)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(wp)                 :: dij(3), rij2 
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj10, lj12, term_lj12, term_lj10
    real(wp)                 :: econt, encont, work(3), coef, cutoff2
    integer                  :: list(2)
    integer                  :: i, ix, icel1, icel2, i1, i2
    integer                  :: iatmcls, jatmcls
    integer                  :: id

    real(wp),        pointer :: contact_lj12(:), contact_lj10(:)
    real(wp),        pointer :: nonb_lj12(:,:)
    integer,         pointer :: contactlist(:,:)
    integer,         pointer :: atmcls(:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:)

    call timer(TimerNonBond, TimerOn)
    call timer(TimerContact, TimerOn)

    ncell_local  => domain%num_cell_local
    id_g2l       => domain%id_g2l
    atmcls       => domain%atom_cls_no

    contactlist  => enefunc%contact_list
    contact_lj12 => enefunc%contact_lj12
    contact_lj10 => enefunc%contact_lj10
    nonb_lj12    => enefunc%nonb_lj12

    cutoff2      = enefunc%cutoffdist * enefunc%cutoffdist

    ! calculate bond energy
    !
    !$omp parallel default(shared) &
    !$omp private(id, i, ix, icel1, i1, icel2, i2, lj10, lj12, dij, rij2, &
    !$omp         inv_rij2, inv_rij6, inv_rij10, inv_rij12, term_lj12,    &
    !$omp         term_lj10, econt, encont, iatmcls, jatmcls, coef, work, &
    !$omp         list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    econt  = 0.0_wp
    encont = 0.0_wp

    do i = id+1, enefunc%num_contact_domain, nthread

      list(1:2) = contactlist(1:2,i)
      i1    = id_g2l(list(1))
      i2    = id_g2l(list(2))
      dij(1)  = coord(i1,1) - coord(i2,1) 
      dij(2)  = coord(i1,2) - coord(i2,2) 
      dij(3)  = coord(i1,3) - coord(i2,3) 

      ! contact energy
      !
      lj10      = contact_lj10(i)
      lj12      = contact_lj12(i)
      rij2      = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6

      term_lj12 = lj12 * inv_rij12
      term_lj10 = lj10 * inv_rij10
      econt     = econt + term_lj12 - term_lj10
      coef      = 12.0_wp*term_lj12 - 10.0_wp*term_lj10

!       ! noncontact energy
!       !
!       if (rij2 < cutoff2) then
!         iatmcls   = atmcls(i1,icel1)
!         jatmcls   = atmcls(i2,icel2)
!         lj12      = nonb_lj12(iatmcls, jatmcls)
!         term_lj12 = lj12 * inv_rij12
!         encont    = encont - term_lj12
!         coef      = coef - 12.0_wp*term_lj12
!       end if

      ! gradient: dE/dX
      !
      coef      = - inv_rij2 * coef
      work(1) = coef * dij(1)
      work(2) = coef * dij(2)
      work(3) = coef * dij(3)

      ! store force: F=-dE/dX
      !
      force(i1,1,id+1) = force(i1,1,id+1) - work(1)
      force(i1,2,id+1) = force(i1,2,id+1) - work(2)
      force(i1,3,id+1) = force(i1,3,id+1) - work(3)
      force(i2,1,id+1) = force(i2,1,id+1) + work(1)
      force(i2,2,id+1) = force(i2,2,id+1) + work(2)
      force(i2,3,id+1) = force(i2,3,id+1) + work(3)

      ! virial
      !
      virial(1,1,id+1) = virial(1,1,id+1) - dij(1)*work(1)
      virial(2,2,id+1) = virial(2,2,id+1) - dij(2)*work(2)
      virial(3,3,id+1) = virial(3,3,id+1) - dij(3)*work(3)

    end do

    econtact(id+1) = econtact(id+1) + econt
    enoncontact(id+1) = enoncontact(id+1) + encont

    !$omp end parallel
    call timer(TimerNonBond, TimerOff)
    call timer(TimerContact, TimerOff)

    return

  end subroutine compute_energy_contact_1210


end module cg_energy_go_mod

