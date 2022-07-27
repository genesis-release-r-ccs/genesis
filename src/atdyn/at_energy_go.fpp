!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_go_mod
!> @brief   calculate contact and non-contact energy (Go potential)
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_go_mod

  use at_pairlist_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_contact_126
  public :: compute_energy_contact_126_multi
  public :: compute_energy_contact_1210
  public :: compute_energy_contact_12106
  public :: compute_energy_contact_126_pbc
  public :: compute_energy_contact_1210_pbc
  public :: compute_energy_contact_AICG2P_pbc
  public :: compute_energy_contact_12106_pbc
  public :: compute_energy_contact_12106_multi
  public :: compute_energy_noncontact_nobc_KBGO
  public :: compute_energy_general_exv_AICG2P
  public :: compute_energy_general_exv_AICG2P_pbc
  public :: compute_energy_noncontact_nobc
  public :: compute_energy_noncontact_pbc
  public :: compute_energy_noncontact_pbc_KBGO

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126
  !> @brief        calculate contact energy for all-atom Go model
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         P.C.Whitford et al., Proteins, 75, 430-441 (2009)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126(enefunc, coord, force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, coef
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj6(:), contact_lj12(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj12 => enefunc%contact_lj12

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij12, coef,   &
    !$omp         work, lj6, lj12,a1, a2, id)                             &
    !$omp shared(istart, iend, list, coord, contact_lj6, contact_lj12,    &
    !$omp        force, nthread)                                          &
    !$omp reduction(+:virial) reduction(+:econt)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-6.0_wp*lj6*inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel  
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_126

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126_multi
  !> @brief        calculate contact energy for AAGO
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126_multi(enefunc, coord, &
                                              force,  virial, energy_mb)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(:,:,:)
    real(wp),                intent(inout) :: energy_mb(:,:)

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, coef
    real(wp)                 :: work(3)
    real(wp)                 :: cutoff2, cutoff
    real(wp)                 :: veps, rmin
    integer                  :: i, j, l, istart, iend, imb
    integer                  :: id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj12(:), contact_lj6(:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: model(:)
    real(wp),        pointer :: nonb_eps(:), nonb_rmin(:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_multi_contact
    iend         =  enefunc%iend_multi_contact
    model        => enefunc%multi_contact_model
    list         => enefunc%multi_contact_list
    contact_lj6  => enefunc%multi_contact_lj6
    contact_lj12 => enefunc%multi_contact_lj12

    cutoff       =  enefunc%cutoffdist
    cutoff2      =  cutoff * cutoff
    nonb_eps     => enefunc%nonb_eps
    nonb_rmin    => enefunc%nonb_rmin

    ! calculate energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6,  inv_rij12,           &
    !$omp         coef, work, lj6, lj12, a1, a2, imb, veps, rmin, j)         &
    !$omp shared(istart, iend, list, coord, contact_lj6, model, cutoff2,     &
    !$omp       contact_lj12, nonb_eps, nonb_rmin)                           &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:energy_mb)
    !
    do i = istart, iend

      ! compute distance
      !
      a1   = list(1,i)
      a2   = list(2,i)
      imb  = model(i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      energy_mb(MBeneCntc,imb) = energy_mb(MBeneCntc,imb) &
         + (lj12 * inv_rij12 - lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12- 6.0_wp * lj6 * inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,imb) = force(1:3,a1,imb) - work(1:3)
      force(1:3,a2,imb) = force(1:3,a2,imb) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l,imb) = virial(1:3,l,imb) - dij(1:3)*work(l)
      end do

      if (rij2 < cutoff2) then
        veps = sqrt(nonb_eps(a1)*nonb_eps(a2))
        rmin = (nonb_rmin(a1)+nonb_rmin(a2))
        rmin = rmin * rmin * rmin
        rmin = rmin * rmin
        lj12 = veps * rmin * rmin
        energy_mb(MBeneNonb,imb) = energy_mb(MBeneNonb,imb) - (lj12 * inv_rij12)

        coef = - inv_rij2 * 12.0_wp*lj12*inv_rij12
        work(1:3) = -coef*dij(1:3)
        force(1:3,a1,imb) = force(1:3,a1,imb) - work(1:3)
        force(1:3,a2,imb) = force(1:3,a2,imb) + work(1:3)
        do l = 1, 3
          virial(1:3,l,imb) = virial(1:3,l,imb) + dij(1:3)*work(l)
        end do
      end if

    end do
    !$omp end parallel do
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_126_multi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_1210
  !> @brief        calculate contact energy for C alpha Go model
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_1210(enefunc, molecule, &
                                         coord, force, virial, econt, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: econt
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                 :: dij(3), rij2, rij
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj10, lj12, coef
    real(wp)                 :: term
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2
    integer                  :: goelec

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: charge(:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: molno(:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerContact, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12
    charge       => molecule%charge
    molno        => molecule%molecule_no

    el_fact      =  ELECOEF/enefunc%dielec_const
    debye        =  enefunc%debye
    goelec       =  enefunc%go_electrostatic

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                             &
    !$omp private(id, i, l, dij, rij2, rij, inv_rij2, inv_rij6, inv_rij10,   &
    !$omp         inv_rij12, coef, work, lj10, lj12, a1, a2, term)           &
    !$omp shared(istart, iend, list, coord, contact_lj10, contact_lj12,      &
    !$omp        charge, molno, el_fact, debye, goelec, nthread, force)      &
    !$omp reduction(+:virial) reduction(+:econt) reduction(+:eelec)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10)
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10)

      ! debye-huckel
      ! 
      if (goelec == GoElectrostaticALL) then
        rij = sqrt(rij2)
        term  = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
        eelec = eelec + term
        coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
      else if (goelec == GoElectrostaticINTER) then
        if (molno(a1) /= molno(a2)) then
          rij = sqrt(rij2)
          term  = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
          eelec = eelec + term
          coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
        end if
      end if

      ! gradient
      !
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel
    call timer(TimerContact, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_1210

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_12106
  !> @brief        calculate contact energy for KBGO
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         Karanicolas & Brooks, Protein Sci., 11, 2351-2361 (2002)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_12106(enefunc, coord, force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj6, lj10, lj12, coef
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: contact_lj6(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                             &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij10, inv_rij12, &
    !$omp         coef, work, lj6, lj10, lj12, a1, a2, id)                   &
    !$omp shared(istart, iend, list, coord, contact_lj6,                     &
    !$omp       contact_lj10, contact_lj12, force, nthread)                  &
    !$omp reduction(+:virial) reduction(+:econt)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10 + lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10 &
                          + 6.0_wp * lj6 * inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_12106

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126_pbc
  !> @brief        calculate contact energy for all-atom Go model
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         P.C.Whitford et al., Proteins, 75, 430-441 (2009)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126_pbc(enefunc, boundary, coord, &
                                            force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, coef
    real(wp)                 :: work(3)
    real(wp)                 :: bsize(3), inv_bsize(3)
    real(wp)                 :: cutoff, cutoff2
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj6(:), contact_lj12(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)

    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj12 => enefunc%contact_lj12

    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    cutoff       =  enefunc%cutoffdist

    bsize(1)     =  boundary%box_size_x
    bsize(2)     =  boundary%box_size_y
    bsize(3)     =  boundary%box_size_z

    inv_bsize(1:3) = 1.0_wp/bsize(1:3)
    cutoff2      =  cutoff * cutoff


    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij12, coef,   &
    !$omp         work, lj6, lj12,a1, a2, id)                             &
    !$omp shared(istart, iend, list, coord, contact_lj6, contact_lj12,    &
    !$omp         bsize, inv_bsize, cutoff2, force, nthread)              &
    !$omp reduction(+:virial) reduction(+:econt)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      dij(1:3) = dij(1:3) - bsize(1:3)*anint(dij(1:3)*inv_bsize(1:3))
      rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

!      if ( rij2 >= cutoff2 ) cycle

      lj6  = contact_lj6(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj6 * inv_rij6)

      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-6.0_wp*lj6*inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)

      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do

    end do

    !$omp end parallel

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_126_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_1210_pbc
  !> @brief        calculate contact energy for C alpha Go model
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_1210_pbc(enefunc, molecule, boundary, &
                                             coord, force, virial, econt, eelec)

    ! formal arguments
    type(s_enefunc), target, intent(in)   :: enefunc
    type(s_molecule),target, intent(in)   :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt
    real(wp),                intent(inout) :: eelec

    ! local variables
    real(wp)                 :: dij(3), rij2, rij
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj10, lj12, coef, term
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3)
    real(wp)                 :: bsize(3), inv_bsize(3)
    real(wp)                 :: cutoff, cutoff2
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2
    integer                  :: goelec

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: charge(:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: molno(:)


    call timer(TimerNonBond, TimerOn)

    list         => enefunc%contact_list
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12
    charge       => molecule%charge
    molno        => molecule%molecule_no

    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    cutoff       =  enefunc%cutoffdist

    bsize(1)     =  boundary%box_size_x
    bsize(2)     =  boundary%box_size_y
    bsize(3)     =  boundary%box_size_z

    inv_bsize(1:3) = 1.0_wp/bsize(1:3)
    cutoff2      =  cutoff * cutoff

    el_fact      =  ELECOEF/enefunc%dielec_const
    debye        =  enefunc%debye
    goelec       =  enefunc%go_electrostatic

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                             &
    !$omp private(id, i, l, dij, rij2, rij, inv_rij2, inv_rij6, inv_rij10,   &
    !$omp         inv_rij12, coef, work, lj10, lj12, a1, a2, term)           &
    !$omp shared(istart, iend, list, coord, contact_lj10, contact_lj12,      &
    !$omp        bsize, inv_bsize, cutoff2, charge, molno, el_fact,          &
    !$omp        debye, goelec, force, nthread)                              &
    !$omp reduction(+:virial) reduction(+:econt) reduction(+:eelec)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      dij(1:3) = dij(1:3) - bsize(1:3)*anint(dij(1:3)*inv_bsize(1:3))
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

!      if ( rij2 >= cutoff2 ) cycle

      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10)
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10)

      ! debye-huckel
      ! 
      if (goelec == GoElectrostaticALL) then
        rij = sqrt(rij2)
        term  = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
        eelec = eelec + term
        coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
      else if (goelec == GoElectrostaticINTER) then
        if (molno(a1) /= molno(a2)) then
          rij = sqrt(rij2)
          term  = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
          eelec = eelec + term
          coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
        end if
      end if

      ! gradient
      !
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)

      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do

    end do
    !$omp end parallel
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_1210_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_12106_pbc
  !> @brief        calculate contact energy for KBGO
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         Karanicolas & Brooks, Protein Sci., 11, 2351-2361 (2002)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_12106_pbc(enefunc, boundary, coord, &
                                              force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj6, lj10, lj12, coef
    real(wp)                 :: work(3)
    real(wp)                 :: bsize(3), inv_bsize(3)
    real(wp)                 :: cutoff, cutoff2
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: contact_lj6(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)

    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12

    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    cutoff       =  enefunc%cutoffdist

    bsize(1)     =  boundary%box_size_x
    bsize(2)     =  boundary%box_size_y
    bsize(3)     =  boundary%box_size_z

    inv_bsize(1:3) = 1.0_wp/bsize(1:3)
    cutoff2      =  cutoff * cutoff


    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                              &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij10, inv_rij12, &
    !$omp         coef, work, lj6, lj10, lj12, a1, a2, id)                   &
    !$omp shared(istart, iend, list, coord, contact_lj6,                     &
    !$omp       contact_lj10, contact_lj12,                                  &
    !$omp       bsize, inv_bsize, cutoff2, force, nthread)                   &
    !$omp reduction(+:virial) reduction(+:econt)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      dij(1:3) = dij(1:3) - bsize(1:3)*anint(dij(1:3)*inv_bsize(1:3))
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

!      if ( rij2 >= cutoff2 ) cycle

      lj6  = contact_lj6(i)
      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10 + lj6 * inv_rij6)

      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10 &
                          + 6.0_wp * lj6 * inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)

      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do

    end do

    !$omp end parallel

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_12106_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_12106_multi
  !> @brief        calculate contact energy for KBGO
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         Karanicolas & Brooks, Protein Sci., 11, 2351-2361 (2002)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_12106_multi(enefunc, coord, &
                                                force,  virial, energy_mb)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(:,:,:)
    real(wp),                intent(inout) :: energy_mb(:,:)

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj6, lj10, lj12, coef
    real(wp)                 :: work(3)
    real(wp)                 :: cutoff2, cutoff
    real(wp)                 :: veps, rmin
    integer                  :: i, j, l, istart, iend, imb
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: contact_lj6(:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: model(:)
    real(wp),        pointer :: nonb_eps(:), nonb_rmin(:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_multi_contact
    iend         =  enefunc%iend_multi_contact
    model        => enefunc%multi_contact_model
    list         => enefunc%multi_contact_list
    contact_lj6  => enefunc%multi_contact_lj6
    contact_lj10 => enefunc%multi_contact_lj10
    contact_lj12 => enefunc%multi_contact_lj12

    cutoff       =  enefunc%cutoffdist
    cutoff2      =  cutoff * cutoff
    nonb_eps     => enefunc%nonb_eps
    nonb_rmin    => enefunc%nonb_rmin

    ! calculate energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij10, inv_rij12, &
    !$omp         coef, work, lj6, lj10, lj12, a1, a2, imb, veps, rmin, j)   &
    !$omp shared(istart, iend, list, coord, contact_lj6, model, cutoff2,     &
    !$omp       contact_lj10, contact_lj12, nonb_eps, nonb_rmin)             &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:energy_mb)
    !
    do i = istart, iend

      ! compute distance
      !
      a1   = list(1,i)
      a2   = list(2,i)
      imb  = model(i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      energy_mb(MBeneCntc,imb) = energy_mb(MBeneCntc,imb) &
         + (lj12 * inv_rij12 - lj10 * inv_rij10 + lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10 &
                          + 6.0_wp * lj6 * inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,imb) = force(1:3,a1,imb) - work(1:3)
      force(1:3,a2,imb) = force(1:3,a2,imb) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l,imb) = virial(1:3,l,imb) - dij(1:3)*work(l)
      end do

      if (rij2 < cutoff2) then
        veps = sqrt(nonb_eps(a1)*nonb_eps(a2))
        rmin = (nonb_rmin(a1)+nonb_rmin(a2))
        rmin = rmin * rmin * rmin
        rmin = rmin * rmin
        lj12 = veps * rmin * rmin
        energy_mb(MBeneNonb,imb) = energy_mb(MBeneNonb,imb) &
           - (lj12 * inv_rij12)

        coef = - inv_rij2 * 12.0_wp*lj12*inv_rij12
        work(1:3) = -coef*dij(1:3)
        force(1:3,a1,imb) = force(1:3,a1,imb) - work(1:3)
        force(1:3,a2,imb) = force(1:3,a2,imb) + work(1:3)
        do l = 1, 3
          virial(1:3,l,imb) = virial(1:3,l,imb) + dij(1:3)*work(l)
        end do
      end if

    end do
    !$omp end parallel do
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_12106_multi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_nobc_KBGO
  !> @brief        calculate non-contact energy with pairlist (NOBC)
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_nobc_KBGO(enefunc, molecule, pairlist, &
                                                 coord, force, virial, encont)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12, term_lj12
    real(wp)                  :: lj12, coef, noncontact_lj12
    real(wp)                  :: work(3)
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: veps, rmin
    integer                   :: i, j, k, l, natom, id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   ::  nthread, my_id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: forcefield
         
    integer, pointer          :: num_nb15_calc(:,:), nb15_calc_list(:,:)
    integer, pointer          :: atmcls(:)
    real(wp), pointer         :: nonb_lj12(:,:)
    real(wp), pointer         :: nonb_eps(:), nonb_rmin(:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    natom          =  molecule%num_atoms
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list
    num_nb15       = 0
    atmcls         => molecule%atom_cls_no

    forcefield     =  enefunc%forcefield
    cutoff         =  enefunc%cutoffdist
    cutoff2        =  cutoff * cutoff
    nonb_eps       => enefunc%nonb_eps
    nonb_rmin      => enefunc%nonb_rmin

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_nb15)                                    & 
    !$omp private(ini_nb15, fin_nb15, i, k, j, l, dij, coef, work,  &
    !$omp         rij2, inv_rij2, inv_rij6, inv_rij12, term_lj12,   &
    !$omp         my_id, lj12, id, veps, rmin)           &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list, cutoff2,     &
    !$omp        coord, nthread,  my_city_rank, nproc_city, atmcls, &
    !$omp        forcefield, nonb_eps, nonb_rmin, force)            &
    !$omp reduction(+:virial) reduction(+:encont) 
    !
#ifdef OMP
    id      = omp_get_thread_num()
    nthread = omp_get_num_threads()
#else
    id      = 0
    nthread = 1
#endif
    my_id   = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      if (mod(i-1,nproc_city*nthread) .ne. my_id) ini_nb15 = fin_nb15 + 1

      do k = ini_nb15, fin_nb15
        
        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        veps = sqrt(nonb_eps(i)*nonb_eps(j))
        rmin = (nonb_rmin(i)+nonb_rmin(j))
        rmin = rmin * rmin * rmin
        rmin = rmin * rmin
        lj12 = veps * rmin * rmin
        ! cutoff
        !
        if (rij2 >= cutoff2) cycle

        ! non-contact energy
        !
        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij12 = inv_rij6 * inv_rij6
        term_lj12 = lj12 * inv_rij12
        
        encont = encont + term_lj12
        
        ! gradient
        !
        coef = - inv_rij2 * (12.0_wp*term_lj12)
        work(1:3) = coef*dij(1:3)
        
        ! store force
        !
        force(1:3,i,id+1) = force(1:3,i,id+1) - work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) + work(1:3)
        
        ! virial
        !
        do l = 1, 3
          virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        end do

      end do

    end do
    !$omp end parallel
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_noncontact_nobc_KBGO

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_AICG2P_pbc
  !> @brief        calculate contact energy for C alpha AICG2P Go model
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_AICG2P_pbc(enefunc, boundary, coord, &
                                               force, virial, econt)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: econt

    ! local variables
    integer           :: i, l, istart, iend, id, omp_get_thread_num
    integer           :: a1, a2
    real(wp)          :: dij(3), rij2, rij
    real(wp)          :: inv_rij2, inv_rij4, inv_rij10, inv_rij12
    real(wp)          :: lj10, lj12, coef, term
    real(wp)          :: work(3)
    real(wp)          :: bsize(3), half_bsize(3)

    real(wp), pointer :: contact_lj10(:), contact_lj12(:)
    integer,  pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerContact, TimerOn)

    list           => enefunc%contact_list
    contact_lj10   => enefunc%contact_lj10
    contact_lj12   => enefunc%contact_lj12

    istart         =  enefunc%istart_contact
    iend           =  enefunc%iend_contact

    bsize(1)        =  boundary%box_size_x
    bsize(2)        =  boundary%box_size_y
    bsize(3)        =  boundary%box_size_z
    half_bsize(1:3) = 0.5_wp * bsize(1:3)

    ! calculate energy and gradient
    !
    !$omp parallel default(none)             &
    !$omp private(i, l, a1, a2, id,          &
    !$omp         dij, rij2, rij,            &
    !$omp         inv_rij2, inv_rij4,        &
    !$omp         inv_rij10, inv_rij12,      &
    !$omp         lj10, lj12, term,          &
    !$omp         coef, work)                &
    !$omp shared(istart, iend, list, coord,  &
    !$omp        contact_lj10, contact_lj12, &
    !$omp        bsize, half_bsize, force,   &
    !$omp        nthread)                    &
    !$omp reduction(+:virial)                &
    !$omp reduction(+:econt) 
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = istart+id, iend, nthread

      ! compute distance
      !
      a1 = list(1, i)
      a2 = list(2, i)

      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      ! dij(1:3) = dij(1:3) - bsize(1:3)*anint(dij(1:3)*inv_bsize(1:3))
      if ( dij(1) > half_bsize(1) ) then
        dij(1) = dij(1) - bsize(1)
      else if ( dij(1) < -half_bsize(1) ) then
        dij(1) = dij(1) + bsize(1)
      end if
      if ( dij(2) > half_bsize(2) ) then
        dij(2) = dij(2) - bsize(2)
      else if ( dij(2) < -half_bsize(2) ) then
        dij(2) = dij(2) + bsize(2)
      end if
      if ( dij(3) > half_bsize(3) ) then
        dij(3) = dij(3) - bsize(3)
      else if ( dij(3) < -half_bsize(3) ) then
        dij(3) = dij(3) + bsize(3)
      end if

      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij4  = inv_rij2 * inv_rij2
      inv_rij10 = inv_rij4 * inv_rij4 * inv_rij2
      inv_rij12 = inv_rij10 * inv_rij2
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10)
      coef = - inv_rij2 * (12.0_wp * lj12 * inv_rij12 - 10.0_wp * lj10 * inv_rij10)

      ! gradient
      !
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,a1,id+1) = force(1:3,a1,id+1) - work(1:3)
      force(1:3,a2,id+1) = force(1:3,a2,id+1) + work(1:3)

      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do

    end do
    !$omp end parallel 

    call timer(TimerContact, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_AICG2P_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_general_exv_AICG2P
  !> @brief        calculate non-contact energy with pairlist (NOBC)
  !! @authors      TM, CK, CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @note         AICG
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_general_exv_AICG2P(enefunc, molecule, pairlist, &
                                               coord, force, virial, encont)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12, term_lj12
    real(wp)                  :: lj12, coef, noncontact_lj12
    real(wp)                  :: work(3)
    real(wp)                  :: cutoff2, cutoff
    real(wp)                  :: ei_sqrt, ej_sqrt
    real(wp)                  :: si_half, sj_half
    real(wp)                  :: veps, rmin
    integer                   :: i, j, k, l, natom, id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: nthread, my_id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: forcefield
         
    integer, pointer          :: num_nb15_calc(:,:), nb15_calc_list(:,:)
    real(wp), pointer         :: param_sigma(:), param_epsilon(:)
    real(wp), parameter       :: exv_cutoff = 2.0
    real(wp), parameter       :: exv_cutoff12 = exv_cutoff**(-12)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGexv, TimerOn)
    ! use pointers
    !
    natom          =  molecule%num_atoms
    num_nb15_calc  => pairlist%num_cg_exv_calc
    nb15_calc_list => pairlist%cg_exv_list
    num_nb15       = 0
    param_epsilon  => enefunc%param_epsilon
    param_sigma    => enefunc%param_sigma

    forcefield     =  enefunc%forcefield

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_nb15)                                    & 
    !$omp private(ini_nb15, fin_nb15, i, k, j, l, dij, coef, work,  &
    !$omp         si_half, sj_half,                                 &
    !$omp         ei_sqrt, ej_sqrt,                                 &
    !$omp         rij2, inv_rij2, inv_rij6, inv_rij12, term_lj12,   &
    !$omp         my_id, lj12, id, veps, rmin, cutoff2)             &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list,              &
    !$omp        coord, nthread,  my_city_rank, nproc_city,         &
    !$omp        param_sigma, param_epsilon,                        &
    !$omp        forcefield, force)                                 &
    !$omp reduction(+:virial) reduction(+:encont) 
    !
#ifdef OMP
    id      = omp_get_thread_num()
    nthread = omp_get_num_threads()
#else
    id      = 0
    nthread = 1
#endif
    my_id   = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      ! if (mod(i-1,nproc_city*nthread) .ne. my_id) ini_nb15 = fin_nb15 + 1

      ei_sqrt             = param_epsilon(i)
      si_half             = param_sigma(i)

      do k = ini_nb15, fin_nb15
        
        j = nb15_calc_list(k,id+1)

        ej_sqrt     = param_epsilon(j)
        sj_half     = param_sigma(j)

        veps   = ei_sqrt * ej_sqrt
        rmin   = si_half + sj_half
        cutoff2 = exv_cutoff*exv_cutoff*rmin*rmin

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! cutoff
        !
        if (rij2 >= cutoff2) cycle

        rmin = rmin * rmin * rmin
        rmin = rmin * rmin
        lj12 = rmin * rmin

        ! non-contact energy
        !
        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij12 = inv_rij6 * inv_rij6
        term_lj12 = veps *  (lj12*inv_rij12)
        
        encont = encont + term_lj12 - exv_cutoff12*veps
        
        ! gradient
        !
        coef = - inv_rij2 * (12.0_wp*term_lj12)
        work(1:3) = coef*dij(1:3)
        
        ! store force
        !
        force(1:3,i,id+1) = force(1:3,i,id+1) - work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) + work(1:3)
        
        ! virial
        !
        do l = 1, 3
          virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        end do

      end do

    end do
    !$omp end parallel
    call timer(TimerCGexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_general_exv_AICG2P

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_general_exv_AICG2P_pbc
  !> @brief        calculate excluded volume energy with pairlist (PBC)
  !! @authors      CT
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eexv     : nonnative contact energy of target systems
  !! @note         AICG2P
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_general_exv_AICG2P_pbc(enefunc, boundary, pairlist,&
                                                   coord, force, virial, eexv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eexv

    ! local variables
    integer             :: my_id, id
    integer             :: omp_get_num_threads, omp_get_thread_num
    !
    integer             :: natom
    integer             :: num_exv, ini_exv, fin_exv
    integer             :: i, j, k, k1, l, i1, i2, i3
    !
    logical             :: proceed
    real(wp)            :: ei_sqrt, ej_sqrt
    real(wp)            :: si_half, sj_half
    real(wp)            :: epsilon
    real(wp)            :: sigma, sigma_sqr
    real(wp)            :: dij(3), rij_sqr
    real(wp)            :: inv_rij_sqr
    real(wp)            :: sig_over_rij_sqr
    real(wp)            :: sig_over_rij_6th, sig_over_rij_12th
    real(wp)            :: grad_coef_exv, grad(3)
    real(wp)            :: eexv_tmp
    !
    real(wp)            :: coord_tmp(3)
    integer             :: atomclass_i, atomclass_j
    real(wp)            :: force_tmp(3), factor
    real(wp)            :: virial_tmp(3,3)
    real(wp)            :: virial_omp(3,3,nthread)
    real(wp)            :: eexv_omp(nthread), eexv_omp_tmp

    integer,  pointer   :: exv_list(:,:)
    integer,  pointer   :: num_exv_calc(:,:)
    real(wp), pointer   :: param_sigma(:), param_epsilon(:)
    real(wp), parameter :: exv_cutoff_factor = 4.0_wp
    real(wp), parameter :: exv_energy_shift = 1.0_wp / 4096.0_wp
    real(wp)            :: bsize(3), inv_bsize(3)


    call timer(TimerNonBond, TimerOn)
    call timer(TimerCGexv, TimerOn)

    natom          = size(coord(1,:))
    num_exv_calc   => pairlist%num_cg_exv_calc
    exv_list       => pairlist%cg_exv_list
    param_epsilon  => enefunc%param_epsilon
    param_sigma    => enefunc%param_sigma
    bsize(1)       =  boundary%box_size_x
    bsize(2)       =  boundary%box_size_y
    bsize(3)       =  boundary%box_size_z
    inv_bsize(1:3) = 1.0_wp/bsize(1:3)

    num_exv      = 0

    eexv_omp(1:nthread) = 0.0_wp
    virial_omp(1:3,1:3,1:nthread) = 0.0_wp

    !$omp parallel default(none)                       &
    !$omp firstprivate(num_exv)                        &
    !$omp private(my_id, id,                           &
    !$omp         ini_exv, fin_exv,                    &
    !$omp         i, k, j, l, epsilon,                 &
    !$omp         proceed,                             &
    !$omp         si_half, sj_half,                    &
    !$omp         ei_sqrt, ej_sqrt,                    &
    !$omp         sigma, sigma_sqr,                    &
    !$omp         dij, rij_sqr,                        &
    !$omp         inv_rij_sqr, sig_over_rij_sqr,       &
    !$omp         sig_over_rij_6th, sig_over_rij_12th, &
    !$omp         grad_coef_exv, grad, eexv_tmp,       &
    !$omp         coord_tmp, force_tmp, factor,        &
    !$omp         atomclass_i, atomclass_j,            &
    !$omp         virial_tmp, eexv_omp_tmp,            &
    !$omp         k1, i1, i2, i3)                      &
    !$omp shared(coord, my_city_rank, nproc_city,      &
    !$omp        nthread, natom,                       &
    !$omp        exv_list, num_exv_calc,               &
    !$omp        param_sigma, param_epsilon,           &
    !$omp        bsize, inv_bsize,                     &
    !$omp        virial_omp, eexv_omp,                 &
    !$omp        force)                                 
    !
#ifdef OMP
    id      = omp_get_thread_num()
    ! nthread = omp_get_num_threads()
#else
    id      = 0
    ! nthread = 1
#endif
    my_id   = my_city_rank * nthread + id
    id      = id + 1

    do i = 1, natom - 1

      proceed = .false.

      ini_exv = num_exv + 1
      fin_exv = num_exv + num_exv_calc(i, id)
      num_exv = fin_exv

      if (fin_exv >= ini_exv) proceed = .true.

      if (proceed) then

        coord_tmp(1:3)      = coord(1:3, i)
        ei_sqrt             = param_epsilon(i)
        si_half             = param_sigma(i)
        force_tmp(1:3)      = 0.0_wp
        virial_tmp(1:3,1:3) = 0.0_wp
        eexv_omp_tmp        = 0.0_wp

        do k = ini_exv, fin_exv

          k1          = exv_list(k, id)
          j           = k1 / 27
          k1          = k1 - j*27
          i3          = k1 / 9
          k1          = k1 - i3*9
          i2          = k1 / 3
          i1          = k1 - i2*3
          i1          = i1 - 1
          i2          = i2 - 1
          i3          = i3 - 1
          ej_sqrt     = param_epsilon(j)
          sj_half     = param_sigma(j)

          ! compute distance
          !
          dij(1)  = coord(1,j) - coord_tmp(1) + bsize(1)*real(i1,wp)
          dij(2)  = coord(2,j) - coord_tmp(2) + bsize(2)*real(i2,wp)
          dij(3)  = coord(3,j) - coord_tmp(3) + bsize(3)*real(i3,wp)
!         dij(1)  = dij(1) - bsize(1) * anint(dij(1) * inv_bsize(1))
!         dij(2)  = dij(2) - bsize(2) * anint(dij(2) * inv_bsize(2))
!         dij(3)  = dij(3) - bsize(3) * anint(dij(3) * inv_bsize(3))
          rij_sqr   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          epsilon   = ei_sqrt * ej_sqrt
          sigma     = si_half + sj_half
          sigma_sqr = sigma * sigma

!         factor = merge(0.0_wp,1.0_wp,rij_sqr>=sigma_sqr * exv_cutoff_factor)
          if (rij_sqr>=sigma_sqr * exv_cutoff_factor) cycle

          inv_rij_sqr       = 1.0_wp / rij_sqr
          sig_over_rij_sqr  = sigma_sqr * inv_rij_sqr

          sig_over_rij_6th  = sig_over_rij_sqr * sig_over_rij_sqr * sig_over_rij_sqr
          sig_over_rij_12th = sig_over_rij_6th * sig_over_rij_6th

          ! exv energy
          !
!         eexv_tmp     = epsilon * ( sig_over_rij_12th - exv_energy_shift ) * factor
          eexv_tmp     = epsilon * ( sig_over_rij_12th - exv_energy_shift ) 
          eexv_omp_tmp = eexv_omp_tmp + eexv_tmp
          
          ! gradient
          !
!         grad_coef_exv = - inv_rij_sqr * 12.0_wp * epsilon * (- sig_over_rij_12th) * factor
          grad_coef_exv = - inv_rij_sqr * 12.0_wp * epsilon * (- sig_over_rij_12th) 
          grad(1)       = grad_coef_exv * dij(1)
          grad(2)       = grad_coef_exv * dij(2)
          grad(3)       = grad_coef_exv * dij(3)

          ! store force
          !
          force_tmp(1) = force_tmp(1) - grad(1)
          force_tmp(2) = force_tmp(2) - grad(2)
          force_tmp(3) = force_tmp(3) - grad(3)
          force(1,j,id) = force(1,j,id) + grad(1)
          force(2,j,id) = force(2,j,id) + grad(2)
          force(3,j,id) = force(3,j,id) + grad(3)

          ! virial
          !
          virial_tmp(1,1) = virial_tmp(1,1) - dij(1)*grad(1)
          virial_tmp(2,1) = virial_tmp(2,1) - dij(2)*grad(1)
          virial_tmp(3,1) = virial_tmp(3,1) - dij(3)*grad(1)
          virial_tmp(1,2) = virial_tmp(1,2) - dij(1)*grad(2)
          virial_tmp(2,2) = virial_tmp(2,2) - dij(2)*grad(2)
          virial_tmp(3,2) = virial_tmp(3,2) - dij(3)*grad(2)
          virial_tmp(1,3) = virial_tmp(1,3) - dij(1)*grad(3)
          virial_tmp(2,3) = virial_tmp(2,3) - dij(2)*grad(3)
          virial_tmp(3,3) = virial_tmp(3,3) - dij(3)*grad(3)

        end do

        force(1:3,i,id) = force(1:3,i,id) + force_tmp(1:3)
        eexv_omp(id) = eexv_omp(id) + eexv_omp_tmp
        virial_omp(1:3,1:3,id) = virial_omp(1:3,1:3,id) + virial_tmp(1:3,1:3)

      end if
    end do
    !$omp end parallel

    do i = 1, nthread
      eexv = eexv + eexv_omp(i)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,i)
    end do

    call timer(TimerCGexv, TimerOff)
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_general_exv_AICG2P_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_nobc
  !> @brief        calculate non-contact energy with pairlist (NOBC)
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_nobc(enefunc, molecule, pairlist, &
                                            coord, force, virial, encont, eelec)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                  :: dij(3), rij2, rij
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12, term_lj12
    real(wp)                  :: lj12, coef, noncontact_lj12
    real(wp)                  :: term
    real(wp)                  :: work(3)
    real(wp)                  :: el_fact, debye
    real(wp)                  :: cutoff2, cutoff
    integer                   :: i, j, k, l, natom, id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: nthread, my_id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: forcefield
    integer                   :: goelec

    integer, pointer          :: num_nb15_calc(:,:), nb15_calc_list(:,:)
    integer, pointer          :: atmcls(:)
    integer, pointer          :: molno(:)
    real(wp), pointer         :: charge(:)
    real(wp), pointer         :: nonb_lj12(:,:)


    call timer(TimerNonBond, TimerOn)

    natom          =  molecule%num_atoms
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list
    num_nb15       = 0
    atmcls         => molecule%atom_cls_no
    charge         => molecule%charge
    molno          => molecule%molecule_no

    noncontact_lj12=  enefunc%noncontact_lj12
    forcefield     =  enefunc%forcefield
    cutoff         =  enefunc%cutoffdist
    cutoff2        =  cutoff * cutoff
    nonb_lj12      => enefunc%nonb_lj12

    el_fact      =  ELECOEF/enefunc%dielec_const
    debye        =  enefunc%debye
    goelec       =  enefunc%go_electrostatic

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_nb15)                                    & 
    !$omp private(ini_nb15, fin_nb15, i, k, j, l, dij, coef, work,  &
    !$omp         rij2, rij, inv_rij2, inv_rij6, inv_rij12,         &
    !$omp         term_lj12, my_id, lj12, id, term)                 &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list, cutoff2,     &
    !$omp        coord, nthread,  my_city_rank, nproc_city, atmcls, &
    !$omp        forcefield, noncontact_lj12, nonb_lj12,            &
    !$omp        charge, molno, el_fact, debye, goelec, force)      &
    !$omp reduction(+:virial) reduction(+:encont) reduction(+:eelec)
    !
#ifdef OMP
    id      = omp_get_thread_num()
    nthread = omp_get_num_threads()
#else
    id      = 0
    nthread = 1
#endif
    my_id   = my_city_rank * nthread + id
!    if (forcefield /= ForcefieldKBGO) lj12 = noncontact_lj12

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      if (mod(i-1,nproc_city*nthread) /= my_id) ini_nb15 = fin_nb15 + 1

      do k = ini_nb15, fin_nb15
        
        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        !        if (forcefield == ForcefieldKBGO)  &
           lj12 = nonb_lj12(atmcls(i),atmcls(j))
        ! cutoff
        !
        if (rij2 >= cutoff2) cycle

        ! non-contact energy
        !
        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij12 = inv_rij6 * inv_rij6
        term_lj12 = lj12 * inv_rij12

        encont = encont + term_lj12
        coef = - inv_rij2 * (12.0_wp*term_lj12)

        ! debye-huckel
        ! 
        if (goelec == GoElectrostaticALL) then
          rij = sqrt(rij2)
          term  = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
          eelec = eelec + term
          coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
        else if (goelec == GoElectrostaticINTER) then
          if (molno(i) /= molno(j)) then
            rij = sqrt(rij2)
            term  = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
            eelec = eelec + term
            coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
          end if
        end if

        ! gradient
        !
        work(1:3) = coef*dij(1:3)

        ! store force
        !
        force(1:3,i,id+1) = force(1:3,i,id+1) - work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) + work(1:3)

        ! virial
        !
        do l = 1, 3
          virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        end do

      end do

    end do
    !$omp end parallel
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_noncontact_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_pbc
  !> @brief        calculate non-contact energy with pairlist (PBC)
  !! @authors      JJ, TM, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : non-contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_pbc(enefunc, molecule,                &
                                           boundary, pairlist, coord, force, &
                                           virial, encont, eelec)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2, rij
    real(wp)                  :: cutoff2
    real(wp)                  :: coef
    real(wp)                  :: term
    real(wp)                  :: work(3)
    real(wp)                  :: el_fact, debye
    real(wp)                  :: lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: box_size(1:3)
    integer                   :: ii, i, j, k, m
    integer                   :: id, my_id
    integer                   :: goelec

    real(wp),         pointer :: cutoff
    real(wp),         pointer :: bsize_x, bsize_y, bsize_z
    real(wp),         pointer :: charge(:)
    integer,          pointer :: molno(:)
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    call timer(TimerNonBond, TimerOn)

    atmcls          => molecule%atom_cls_no
    nsolute         => enefunc%table%num_solute
    solute_list     => enefunc%table%solute_list
    cutoff          => enefunc%cutoffdist
    num_nb15_calc   => pairlist%table%num_nb15_calc
    nb15_calc_list  => pairlist%table%nb15_calc_list
    charge          => molecule%charge
    molno           => molecule%molecule_no

    box_size(1)     =  boundary%box_size_x
    box_size(2)     =  boundary%box_size_y
    box_size(3)     =  boundary%box_size_z

    el_fact      =  ELECOEF/enefunc%dielec_const
    debye        =  enefunc%debye
    goelec       =  enefunc%go_electrostatic

    cutoff2  = cutoff * cutoff

    ! calculate energy and gradient
    !
    !$omp parallel                                                    &
    !$omp private(id, my_id, i, ii, k, j, m, dij,                     &
    !$omp         rij2, rij, work, coef, rtmp,                        &
    !$omp         inv_rij2, inv_rij6, inv_rij12,                      &
    !$omp         term_lj12, term_lj6, lj12, lj6, term,               &
    !$omp         force_local)                                        &
    !$omp shared(charge, molno, el_fact, debye, goelec, force)        &
    !$omp reduction(+:virial) reduction(+:encont) reduction(+:eelec)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    ! solute-solute
    !
    do ii = 1, nsolute-1

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(ii)

        j = nb15_calc_list(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3)-box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))

          inv_rij2  = 1.0_wp / rij2
          inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
          inv_rij12 = inv_rij6 * inv_rij6
          term_lj12 = lj12 * inv_rij12

          encont = encont + term_lj12
          coef = - inv_rij2 * (12.0_wp*term_lj12)

          ! debye-huckel
          ! 
          if (goelec == GoElectrostaticALL) then
            rij = sqrt(rij2)
            term  = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
            eelec = eelec + term
            coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
          else if (goelec == GoElectrostaticINTER) then
            if (molno(i) /= molno(j)) then
              rij = sqrt(rij2)
              term  = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
              eelec = eelec + term
              coef  = coef - term/rij * (1.0_wp/debye + 1.0_wp/rij)
            end if
          end if

          ! gradient
          !
          work(1:3) = coef*dij(1:3)

          ! store force
          !
          force_local(1:3) = force_local(1:3) - work(1:3)
          force(1:3,j,id+1) = force(1:3,j,id+1) + work(1:3)

          ! virial
          !
          do m = 1, 3
            virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
          end do

        end if

      end do

      force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

    end do

    !$omp end parallel

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_noncontact_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_pbc_KBGO
  !> @brief        calculate non-contact energy with pairlist (PBC)
  !! @authors      JJ, TM, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : non-contact energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_pbc_KBGO(enefunc, molecule,           &
                                           boundary, pairlist, coord, force, &
                                           virial, encont)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont

    ! local variables
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: coef
    real(wp)                  :: work(3)
    real(wp)                  :: lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: box_size(1:3)
    integer                   :: ii, i, j, k, m
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff
    real(wp),         pointer :: bsize_x, bsize_y, bsize_z
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)

    ! KASA add the variables 
    real(wp),         pointer :: nonb_eps(:), nonb_rmin(:)
    real(wp),         pointer :: veps, rmin 
     
#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    call timer(TimerNonBond, TimerOn)

    atmcls          => molecule%atom_cls_no
    nsolute         => enefunc%table%num_solute
    solute_list     => enefunc%table%solute_list
    cutoff          => enefunc%cutoffdist
    num_nb15_calc   => pairlist%table%num_nb15_calc
    nb15_calc_list  => pairlist%table%nb15_calc_list

    box_size(1)     =  boundary%box_size_x
    box_size(2)     =  boundary%box_size_y
    box_size(3)     =  boundary%box_size_z

    cutoff2  = cutoff * cutoff

    ! KASA add the links 
    nonb_eps        => enefunc%nonb_eps
    nonb_rmin       => enefunc%nonb_rmin 


    ! calculate energy and gradient
    !
    !$omp parallel                                                    &
    !$omp private(id, my_id, i, ii, k, j, m, dij,                     &
    !$omp         rij2, work, coef, rtmp,                             &
    !$omp         inv_rij2, inv_rij6, inv_rij12,                      &
    !$omp         term_lj12, term_lj6, lj12, lj6, veps, rmin,         &
    !$omp         force_local)                                        &
    !$omp reduction(+:virial) reduction(+:encont)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    ! solute-solute
    !
    do ii = 1, nsolute-1

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(ii)

        j = nb15_calc_list(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3)-box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! KASA add the following codes
        veps  = sqrt(nonb_eps(i)*nonb_eps(j))
        rmin  = (nonb_rmin(i)+nonb_rmin(j))
        rmin  = rmin * rmin * rmin
        rmin  = rmin * rmin
        lj12  = veps * rmin * rmin 

        if (rij2 < cutoff2) then

          !lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))

          inv_rij2  = 1.0_wp / rij2
          inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
          inv_rij12 = inv_rij6 * inv_rij6
          term_lj12 = lj12 * inv_rij12

          encont = encont + term_lj12

          ! gradient
          !
          coef = - inv_rij2 * (12.0_wp*term_lj12)
          work(1:3) = coef*dij(1:3)

          ! store force
          !
          force_local(1:3) = force_local(1:3) - work(1:3)
          force(1:3,j,id+1)= force(1:3,j,id+1) + work(1:3)

          ! virial
          !
          do m = 1, 3
            virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
          end do

        end if

      end do

      force(1:3,i,id+1) = force(1:3,i,id+1) + force_local(1:3)

    end do

    !$omp end parallel

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_noncontact_pbc_KBGO

end module at_energy_go_mod
