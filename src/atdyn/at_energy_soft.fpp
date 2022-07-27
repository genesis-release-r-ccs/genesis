!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_soft_mod
!> @brief   calculate contact and non-contact energy (Soft potential)
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK), Tadashi Ando (TA)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_soft_mod

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
  public :: compute_energy_contact_soft_nobc
  public :: compute_energy_contact_soft_pbc
  public :: compute_energy_noncontact_soft_nobc
  public :: compute_energy_noncontact_soft_pbc

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_soft_nobc
  !> @brief        calculate contact energy for SOFT model (NOBC)
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_soft_nobc(enefunc, molecule, &
                                             coord, force, virial, econt, eelec)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_molecule),target, intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt
    real(wp),                intent(inout) :: eelec

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, rij
    real(wp)                 :: eps, sig, khh, coef
    real(wp)                 :: term0, term1, term2
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: eij(:), sij(:), kij(:), charge(:)
    integer,         pointer :: list(:,:), func(:)


    call timer(TimerNonBond, TimerOn)

    list    => enefunc%contact_list
    func    => enefunc%contact_func
    eij     => enefunc%contact_lj12
    sij     => enefunc%contact_lj10
    kij     => enefunc%contact_lj6
    charge  => molecule%charge

    istart  =  enefunc%istart_contact
    iend    =  enefunc%iend_contact

    cutoff  =  enefunc%cutoffdist
    cutoff2 =  cutoff * cutoff

    el_fact =  ELECOEF/enefunc%dielec_const
    debye   =  enefunc%debye

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, rij, eps, sig, khh, coef,    &
    !$omp         term0, term1, term2, work, a1, a2, id)                  &
    !$omp shared(istart, iend, list, func, coord, charge, cutoff, cutoff2,&
    !$omp        el_fact, eij, sij, kij, debye, force, nthread)           &
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
      rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      if (rij2 >= cutoff2) cycle

      inv_rij2  = 1.0_wp / rij2

      eps = eij(i)
      sig = sij(i)
      khh = kij(i)

      if      (func(i) == 10) then

        ! 12-6
        term0 = sig**2*inv_rij2
        term2 = term0**3
        term1 = term2*term2
        econt = econt + 4.0_wp*eps*(term1 - term2)
        coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)

      else if (func(i) == 11) then

        ! 12-10
        eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
        term0 = sig**2*inv_rij2
        term2 = term0**5
        term1 = term2*term0
        econt = econt + eps*(term1 - term2)
        coef  = - inv_rij2 * eps*(12.0_wp*term1 - 10.0_wp*term2)

      else if (func(i) == 14) then

        ! 12-6-hh
        if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
          rij = sqrt(rij2)
          econt = econt + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2 - eps
          coef  = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
        else
          term0 = sig**2*inv_rij2
          term2 = term0**3
          term1 = term2*term2
          econt = econt + 4.0_wp*eps*(term1 - term2)
          coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
        end if

      else if (func(i) == 15) then

        ! 12-10-hh
        if (rij2 < (sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2) then
          rij = sqrt(rij2)
          econt = econt + khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2 - eps
          coef  = 2.0_wp*khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp)) / rij
        else
          eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
          term0 = sig**2*inv_rij2
          term2 = term0**5
          term1 = term2*term0
          econt = econt + eps*(term1 - term2)
          coef  = - inv_rij2 * eps*(12.0_wp*term1 - 10.0_wp*term2)
        end if

      end if

      if (debye /= 0.0_wp) then
        rij = sqrt(rij2)
        term0 = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
        eelec = eelec + term0
        coef  = coef - term0/rij * (1.0_wp/debye + 1.0_wp/rij)
      end if

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

  end subroutine compute_energy_contact_soft_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_soft_pbc
  !> @brief        calculate contact energy for SOFT model (PBC)
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_soft_pbc(enefunc, molecule, boundary, &
                                             coord, force, virial, econt, eelec)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_molecule),target, intent(in)    :: molecule
    type(s_boundary),target, intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt
    real(wp),                intent(inout) :: eelec

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, rij
    real(wp)                 :: eps, sig, khh, coef
    real(wp)                 :: term0, term1, term2
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3)
    real(wp)                 :: box_size(3)
    integer                  :: i, l, istart, iend, id, omp_get_thread_num
    integer                  :: a1, a2

    real(wp),        pointer :: eij(:), sij(:), kij(:), charge(:)
    integer,         pointer :: list(:,:), func(:)


    call timer(TimerNonBond, TimerOn)

    list        => enefunc%contact_list
    func        => enefunc%contact_func
    eij         => enefunc%contact_lj12
    sij         => enefunc%contact_lj10
    kij         => enefunc%contact_lj6
    charge      => molecule%charge

    istart      =  enefunc%istart_contact
    iend        =  enefunc%iend_contact

    cutoff      =  enefunc%cutoffdist
    cutoff2     =  cutoff * cutoff

    el_fact     =  ELECOEF/enefunc%dielec_const
    debye       =  enefunc%debye

    box_size(1) =  boundary%box_size_x
    box_size(2) =  boundary%box_size_y
    box_size(3) =  boundary%box_size_z

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                       &
    !$omp private(i, l, dij, rij2, inv_rij2, rij, eps, sig, khh, coef,    &
    !$omp         term0, term1, term2, work, a1, a2, id)                  &
    !$omp shared(istart, iend, list, func, coord, charge, cutoff, cutoff2,&
    !$omp        el_fact, eij, sij, kij, debye, box_size, force, nthread) &
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
      dij(1:3) = dij(1:3)-box_size(1:3)*anint(dij(1:3)/box_size(1:3))
      rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      if (rij2 >= cutoff2) cycle

      inv_rij2  = 1.0_wp / rij2

      eps = eij(i)
      sig = sij(i)
      khh = kij(i)

      if      (func(i) == 10) then

        ! 12-6
        term0 = sig**2*inv_rij2
        term2 = term0**3
        term1 = term2*term2
        econt = econt + 4.0_wp*eps*(term1 - term2)
        coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)

      else if (func(i) == 11) then

        ! 12-10
        eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
        term0 = sig**2*inv_rij2
        term2 = term0**5
        term1 = term2*term0
        econt = econt + eps*(term1 - term2)
        coef  = - inv_rij2 * eps* (12.0_wp*term1 - 10.0_wp*term2)

      else if (func(i) == 14) then

        ! 12-6-hh
        if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
          rij = sqrt(rij2)
          econt = econt + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2 - eps
          coef  = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
        else
          term0 = sig**2*inv_rij2
          term2 = term0**3
          term1 = term2*term2
          econt = econt + 4.0_wp*eps*(term1 - term2)
          coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
        end if

      else if (func(i) == 15) then

        ! 12-10-hh
        if (rij2 < (sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2) then
          rij = sqrt(rij2)
          econt = econt + khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2 - eps
          coef  = 2.0_wp*khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp)) / rij
        else
          eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
          term0 = sig**2*inv_rij2
          term2 = term0**5
          term1 = term2*term0
          econt = econt + eps*(term1 - term2)
          coef  = - inv_rij2 * eps*(12.0_wp*term1 - 10.0_wp*term2)
        end if

      end if

      if (debye /= 0.0_wp) then
        rij = sqrt(rij2)
        term0 = el_fact*charge(a1)*charge(a2)*exp(-rij/debye)/rij
        eelec = eelec + term0
        coef  = coef - term0/rij * (1.0_wp/debye + 1.0_wp/rij)
      end if

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

  end subroutine compute_energy_contact_soft_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_soft_nobc
  !> @brief        calculate non-contact energy with pairlist (NOBC)
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_soft_nobc(enefunc, molecule, pairlist, &
                                            coord, force, virial, encont, eelec)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_molecule),target, intent(in)    :: molecule
    type(s_pairlist),target, intent(in)    :: pairlist
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: encont
    real(wp),                intent(inout) :: eelec


    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, rij
    real(wp)                 :: eps, sig, khh, coef
    real(wp)                 :: term0, term1, term2
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3)
    integer                  :: i, j, k, l, natom, id
    integer                  :: num_nb15, ini_nb15, fin_nb15
    integer                  :: nthread, my_id, nonb_func
    integer                  :: omp_get_num_threads, omp_get_thread_num

    integer,  pointer        :: num_nb15_calc(:,:), nb15_calc_list(:,:)
    integer,  pointer        :: atmcls(:)
    real(wp), pointer        :: eij(:,:), sij(:,:), kij(:,:), charge(:)


    call timer(TimerNonBond, TimerOn)

    atmcls         => molecule%atom_cls_no
    charge         => molecule%charge
    eij            => enefunc%nonb_lj12
    sij            => enefunc%nonb_lj10
    kij            => enefunc%nonb_lj6
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list

    natom          =  molecule%num_atoms
    nonb_func      =  enefunc%nonb_func
    cutoff         =  enefunc%cutoffdist
    cutoff2        =  cutoff * cutoff
    el_fact        =  ELECOEF/enefunc%dielec_const
    debye          =  enefunc%debye

    num_nb15       =  0


    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                     &
    !$omp firstprivate(num_nb15)                                     &
    !$omp private(ini_nb15, fin_nb15, i, k, j, l, dij, coef,         &
    !$omp         rij2, inv_rij2, rij, eps, sig, khh,                &
    !$omp         term0, term1, term2, my_id, id, work)              &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list, cutoff2,      &
    !$omp        coord, nthread, my_city_rank, nproc_city, atmcls,   &
    !$omp        el_fact, debye, eij, sij, kij, nonb_func, charge,   &
    !$omp        force),                                             &
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
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 >= cutoff2) cycle

        inv_rij2  = 1.0_wp / rij2

        eps = eij(atmcls(i),atmcls(j))
        sig = sij(atmcls(i),atmcls(j))
        khh = kij(atmcls(i),atmcls(j))

        if      (nonb_func == 10) then

          ! 12-6
          term0 = sig**2*inv_rij2
          term2 = term0**3
          term1 = term2*term2
          encont = encont + 4.0_wp*eps*(term1 - term2)
          coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)

        else if (nonb_func == 11) then

          ! 12-10
          eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
          term0 = sig**2*inv_rij2
          term2 = term0**5
          term1 = term2*term0
          encont = encont + eps*(term1 - term2)
          coef  = - inv_rij2 * eps* (12.0_wp*term1 - 10.0_wp*term2)

        else if (nonb_func == 12) then

          ! WCA
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            term0 = sig**2*inv_rij2
            term2 = term0**3
            term1 = term2*term2
            encont = encont + 4.0_wp*eps*(term1 - term2) + eps
            coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
          else
            coef = 0.0_wp
          end if

        else if (nonb_func == 13) then

          ! hh
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2
            coef  = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
          else
            coef  = 0.0_wp
          end if

        else if (nonb_func == 14) then

          ! 12-6-hh
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2 - eps
            coef  = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
          else
            term0 = sig**2*inv_rij2
            term2 = term0**3
            term1 = term2*term2
            encont = encont + 4.0_wp*eps*(term1 - term2)
            coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
          end if

        else if (nonb_func == 15) then

          ! 12-10-hh
          if (rij2 < (sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2 - eps
            coef  = 2.0_wp*khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp)) / rij
          else
            eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
            term0 = sig**2*inv_rij2
            term2 = term0**5
            term1 = term2*term0
            encont = encont + eps*(term1 - term2)
            coef  = - inv_rij2 * eps*(12.0_wp*term1 - 10.0_wp*term2)
          end if

        end if

        if (debye /= 0.0_wp) then
          rij = sqrt(rij2)
          term0 = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
          eelec = eelec + term0
          coef  = coef - term0/rij * (1.0_wp/debye + 1.0_wp/rij)
        end if

        ! store force
        !
        work(1:3) = coef*dij(1:3)
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

  end subroutine compute_energy_noncontact_soft_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_soft_pbc
  !> @brief        calculate non-contact energy with pairlist (PBC)
  !! @authors      TM, CK, TA
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_soft_pbc(enefunc, molecule, boundary, &
                                  pairlist, coord, force, virial, encont, eelec)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_molecule),target, intent(in)    :: molecule
    type(s_boundary),target, intent(in)    :: boundary
    type(s_pairlist),target, intent(in)    :: pairlist
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: encont
    real(wp),                intent(inout) :: eelec

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, rij
    real(wp)                 :: eps, sig, khh, coef
    real(wp)                 :: term0, term1, term2
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: el_fact, debye
    real(wp)                 :: work(3), box_size(3)
    integer                  :: i, j, k, l, natom, id
    integer                  :: nthread, my_id, nonb_func
    integer                  :: omp_get_num_threads, omp_get_thread_num

    integer,  pointer        :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,  pointer        :: atmcls(:)
    real(wp), pointer        :: eij(:,:), sij(:,:), kij(:,:), charge(:)


    call timer(TimerNonBond, TimerOn)

    atmcls         => molecule%atom_cls_no
    charge         => molecule%charge
    eij            => enefunc%nonb_lj12
    sij            => enefunc%nonb_lj10
    kij            => enefunc%nonb_lj6
    num_nb15_calc  => pairlist%table%num_nb15_calc
    nb15_calc_list => pairlist%table%nb15_calc_list

    natom          =  molecule%num_atoms
    nonb_func      =  enefunc%nonb_func
    cutoff         =  enefunc%cutoffdist
    cutoff2        =  cutoff * cutoff
    el_fact        =  ELECOEF/enefunc%dielec_const
    debye          =  enefunc%debye
    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z


    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp private(i, k, j, l, dij, coef, work, rij2, inv_rij2, rij, &
    !$omp         eps, sig, khh, term0, term1, term2, my_id, id)    &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list, cutoff2,     &
    !$omp        coord, nthread, my_city_rank, nproc_city, atmcls,  &
    !$omp        el_fact, debye, eij, sij, kij, nonb_func, charge,  &
    !$omp        box_size, force)                                   &
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

    do i = 1, natom - 1

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(i)

        j = nb15_calc_list(k,i)

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        dij(1:3) = dij(1:3)-box_size(1:3)*anint(dij(1:3)/box_size(1:3))
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 >= cutoff2) cycle

        inv_rij2  = 1.0_wp / rij2

        eps = eij(atmcls(i),atmcls(j))
        sig = sij(atmcls(i),atmcls(j))
        khh = kij(atmcls(i),atmcls(j))

        if      (nonb_func == 10) then

          ! 12-6
          term0 = sig**2*inv_rij2
          term2 = term0**3
          term1 = term2*term2
          encont = encont + 4.0_wp*eps*(term1 - term2)
          coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)

        else if (nonb_func == 11) then

          ! 12-10
          eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
          term0 = sig**2*inv_rij2
          term2 = term0**5
          term1 = term2*term0
          encont = encont + eps*(term1 - term2)
          coef  = - inv_rij2 * eps* (12.0_wp*term1 - 10.0_wp*term2)

        else if (nonb_func == 12) then

          ! WCA
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            term0 = sig**2*inv_rij2
            term2 = term0**3
            term1 = term2*term2
            encont = encont + 4.0_wp*eps*(term1 - term2) + eps
            coef = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
          else
            coef = 0.0_wp
          end if

        else if (nonb_func == 13) then

          ! hh
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2
            coef = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
          else
            coef = 0.0_wp
          end if

        else if (nonb_func == 14) then

          ! 12-6-hh
          if (rij2 < (sig*2.0_wp**(1.0_wp/6.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp))**2 - eps
            coef  = 2.0_wp*khh*(rij - sig*2.0_wp**(1.0_wp/6.0_wp)) / rij
          else
            term0 = sig**2*inv_rij2
            term2 = term0**3
            term1 = term2*term2
            encont = encont + 4.0_wp*eps*(term1 - term2)
            coef  = - inv_rij2 * 4.0_wp*eps*(12.0_wp*term1 - 6.0_wp*term2)
          end if

        else if (nonb_func == 15) then

          ! 12-10-hh
          if (rij2 < (sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2) then
            rij = sqrt(rij2)
            encont = encont + khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp))**2 - eps
            coef  = 2.0_wp*khh*(rij - sig*(12.0_wp/10.0_wp)**(1.0_wp/2.0_wp)) / rij
          else
            eps   = (6.0_wp*(6.0_wp/5.0_wp)**5)*eps
            term0 = sig**2*inv_rij2
            term2 = term0**5
            term1 = term2*term0
            encont = encont + eps*(term1 - term2)
            coef  = - inv_rij2 * eps*(12.0_wp*term1 - 10.0_wp*term2)
          end if

        end if

        if (debye /= 0.0_wp) then
          rij = sqrt(rij2)
          term0 = el_fact*charge(i)*charge(j)*exp(-rij/debye)/rij
          eelec = eelec + term0
          coef  = coef - term0/rij * (1.0_wp/debye + 1.0_wp/rij)
        end if

        ! store force
        !
        work(1:3) = coef*dij(1:3)
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

  end subroutine compute_energy_noncontact_soft_pbc

end module at_energy_soft_mod
