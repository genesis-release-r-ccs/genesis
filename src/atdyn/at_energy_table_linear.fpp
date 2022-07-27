!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_table_linear_mod
!> @brief   calculate nonbonded energy with table
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!!          Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_table_linear_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_pairlist_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! variables
  real(wp), save, allocatable :: force_i(:,:)
  real(wp), save, allocatable :: force_j(:,:,:)

  ! subroutines
  public  :: compute_energy_nonbond14_table_linear
  public  :: compute_energy_nonbond_table_linear
  public  :: compute_energy_nonbond_table_linear_check
  public  :: compute_force_nonbond_table_linear

  private :: compute_energy_nonbond14_table_linear_charmm
  private :: compute_energy_nonbond14_table_linear_gro_amber
  private :: compute_energy_nonbond14_table_linear_charmm_check
  private :: compute_energy_nonbond14_table_linear_gro_amber_check
  private :: compute_energy_nonbond_table_solute_linear
  private :: compute_energy_nonbond_table_water_linear
  private :: compute_energy_nonbond_table_solute_linear_check
  private :: compute_energy_nonbond_table_water_linear_check
  private :: compute_force_nonbond_table_solute_linear
  private :: compute_force_nonbond_table_water_linear

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear
  !> @brief        Calculate PME nonbond 1-4 term real part + Lennard-Jones
  !!               with linear 1/R**2 table
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear( &
                                          enefunc, molecule, coord, &
                                          force, virial, eelec, evdw)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw


    if (enefunc%nonb_limiter) then
      ! ==> Type 5
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldCHARMM19) then

        call compute_energy_nonbond14_table_linear_charmm_check( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)

      ! ==> Type 11
      else  if (enefunc%forcefield == ForcefieldAMBER  .or. &
                enefunc%forcefield == ForcefieldGROAMBER .or. &
                enefunc%forcefield == ForcefieldGROMARTINI) then


        call compute_energy_nonbond14_table_linear_gro_amber_check( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)
      end if 
    else 
      ! ==> Type 5
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldCHARMM19) then

        call compute_energy_nonbond14_table_linear_charmm( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)

      ! ==> Type 11
      else  if (enefunc%forcefield == ForcefieldAMBER  .or. &
                enefunc%forcefield == ForcefieldGROAMBER .or. &
                enefunc%forcefield == ForcefieldGROMARTINI) then


        call compute_energy_nonbond14_table_linear_gro_amber( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)
      end if
    end if

    return

  end subroutine compute_energy_nonbond14_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear
  !> @brief        Calculate nonbonded energy with lookup table
  !! @authors      NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    boundary : boundary information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear( &
                                          enefunc, molecule, pairlist,    &
                                          boundary, coord, force, virial, &
                                          eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    ! Compute solute-solute and solute-water interactions
    !
    call compute_energy_nonbond_table_solute_linear(enefunc, molecule, &
                                boundary, pairlist, coord, force,      &
                                virial, eelec, evdw)

    ! Compute water-water interactions
    !
    call compute_energy_nonbond_table_water_linear(enefunc, molecule,  &
                                boundary, pairlist, coord, force,      &
                                virial, eelec, evdw)

    return

  end subroutine compute_energy_nonbond_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_check
  !> @brief        Calculate nonbonded energy with lookup table
  !! @authors      NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    boundary : boundary information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_check( &
                                          enefunc, molecule, pairlist,    &
                                          boundary, coord, force, virial, &
                                          eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    ! Compute solute-solute and solute-water interactions
    !
    call compute_energy_nonbond_table_solute_linear_check(enefunc, molecule, &
                                boundary, pairlist, coord, force,      &
                                virial, eelec, evdw)

    ! Compute water-water interactions
    !
    call compute_energy_nonbond_table_water_linear_check(enefunc, molecule,  &
                                boundary, pairlist, coord, force,      &
                                virial, eelec, evdw)

    return

  end subroutine compute_energy_nonbond_table_linear_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_linear
  !> @brief        Calculate nonbonded force with lookup table
  !! @authors      NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    boundary : boundary information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_linear( &
                                         enefunc, molecule, pairlist, &
                                         boundary, coord, force, virial)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)


    ! Compute solute-solute and solute-water interactions
    !
    call compute_force_nonbond_table_solute_linear(enefunc, molecule, &
                                boundary, pairlist, coord, force, virial)

    ! Compute water-water interactions
    !
    call compute_force_nonbond_table_water_linear(enefunc, molecule,  &
                                boundary, pairlist, coord, force, virial)

    return

  end subroutine compute_force_nonbond_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm
  !> @brief        Calculate PME nonbond 1-4 term real part + Lennard-Jones
  !!               with linear 1/R**2 table
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm( &
                                         enefunc, molecule, coord, &
                                         force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, table(6)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: natom
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom          => enefunc%table%num_solute
    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    cutoff         => enefunc%cutoffdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    table_ene      => enefunc%table%table_ene
    table_grad     => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,enefunc%table%num_solute))
               !force_j(3,size(enefunc%nb14_calc_list)))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp


    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef, lj6, lj12, term_lj12, term_lj6,       &
    !$omp         term_elec, table, force_local)                               &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, natom-1

      force_local(1:3) = 0.0_wp

      i = enefunc%table%solute_list(ii)

      if (mod(ii-1,nproc_city*nthread) == my_id) then

        do k = 1, num_nb14_calc(ii)

          j = nb14_calc_list(k,ii)

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2 = cutoff2 * enefunc%table%density / rij2
            L = int(rij2)
            R = rij2 - L
            lj6  = enefunc%nb14_lj6(atmcls(i),atmcls(j))
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j))

            table(1:6) = table_ene(3*L-2:3*L+3)
            term_lj12 = table(1) + R*(table(4)-table(1))
            term_lj6  = table(2) + R*(table(5)-table(2))
            term_elec = table(3) + R*(table(6)-table(3))
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec

            table(1:6) = table_grad(3*L-2:3*L+3)
            term_lj12 = table(1) + R*(table(4)-table(1))
            term_lj6  = table(2) + R*(table(5)-table(2))
            term_elec = table(3) + R*(table(6)-table(3))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + charge(i)*charge(j)*term_elec

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k,ii) = work(1:3)

            ! virial
            !
            do m = 1, 3
              virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
            end do

          end if

        end do

        force_i(1:3,ii) = force_local(1:3)

      end if

    end do

    num_nb14_all  = num_nb14
    !$omp end parallel

    do ii = 1, natom-1
       i = enefunc%table%solute_list(ii)
       force(1,i,1) = force(1,i,1) + force_i(1,ii)
       force(2,i,1) = force(2,i,1) + force_i(2,ii)
       force(3,i,1) = force(3,i,1) + force_i(3,ii)
    end do

    do ii = 1, natom-1
       do k = 1, num_nb14_calc(ii)
          j = nb14_calc_list(k,ii)
          force(1,j,1) = force(1,j,1) + force_j(1,k,ii)
          force(2,j,1) = force(2,j,1) + force_j(2,k,ii)
          force(3,j,1) = force(3,j,1) + force_j(3,k,ii)
       end do
    end do

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber
  !> @brief        Calculate PME nonbond 1-4 term real part + Lennard-Jones
  !!               with linear 1/R**2 table
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber( &
                                            enefunc, molecule, coord, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, table(6)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: qq_scale, lj_scale
    real(wp)                  :: inv_r12, inv_r121, inv_r123, inv_r126,inv_r1212
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    real(wp),         pointer :: nb14_qq_scale(:,:), nb14_lj_scale(:,:)
    integer,          pointer :: natom
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom          => enefunc%table%num_solute
    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    cutoff         => enefunc%cutoffdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    nb14_qq_scale  => enefunc%nb14_qq_scale
    nb14_lj_scale  => enefunc%nb14_lj_scale
    table_ene      => enefunc%table%table_ene
    table_grad     => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,natom))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp


    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef, lj6, lj12, term_lj12, term_lj6,       &
    !$omp         term_elec, table, force_local, qq_scale, lj_scale,           &
    !$omp         inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212)            &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
   my_id = my_city_rank * nthread + id

    do ii = 1, natom - 1

      force_local(1:3) = 0.0_wp

      i = enefunc%table%solute_list(ii)

      if (mod(ii-1,nproc_city*nthread) == my_id) then

        do k = 1, num_nb14_calc(ii)

          j        = nb14_calc_list(k,ii)
          qq_scale = nb14_qq_scale(k,ii)
          lj_scale = nb14_lj_scale(k,ii)

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r12  = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2 = cutoff2 * enefunc%table%density / rij2
            L = int(rij2)
            R = rij2 - L

            lj6  = enefunc%nb14_lj6 (atmcls(i),atmcls(j)) * lj_scale
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j)) * lj_scale

            table(1:6) = table_ene(3*L-2:3*L+3)
            term_lj12 = inv_r1212 
            term_lj6  = inv_r126   
            term_elec = table(3) + R*(table(6)-table(3))
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec * qq_scale

            table(1:6) = table_grad(3*L-2:3*L+3)
            term_lj12 = -12.0_wp * inv_r1212 * inv_r12
            term_lj6  = -6.0_wp * inv_r126 * inv_r12
            term_elec = table(3) + R*(table(6)-table(3))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + charge(i)*charge(j)*term_elec * qq_scale

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k,ii) = work(1:3)

            ! virial
            !
            do m = 1, 3
              virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
            end do

          end if

        end do

        force_i(1:3,ii) = force_local(1:3)

      end if

    end do

    num_nb14_all  = num_nb14
    !$omp end parallel

    do ii = 1, natom-1
      i = enefunc%table%solute_list(ii)
      force(1,i,1) = force(1,i,1) + force_i(1,ii)
      force(2,i,1) = force(2,i,1) + force_i(2,ii)
      force(3,i,1) = force(3,i,1) + force_i(3,ii)
    end do

    do ii = 1, natom-1
      do k = 1, num_nb14_calc(ii)
        j = nb14_calc_list(k,ii)
        force(1,j,1) = force(1,j,1) + force_j(1,k,ii)
        force(2,j,1) = force(2,j,1) + force_j(2,k,ii)
        force(3,j,1) = force(3,j,1) + force_j(3,k,ii)
      end do
    end do

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm_check
  !> @brief        Calculate PME nonbond 1-4 term real part + Lennard-Jones
  !!               with linear 1/R**2 table
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm_check( &
                                                     enefunc, molecule, coord, &
                                                     force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, table(6)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: minimum_contact
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: natom
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom          => enefunc%table%num_solute
    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    cutoff         => enefunc%cutoffdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    table_ene      => enefunc%table%table_ene
    table_grad     => enefunc%table%table_grad

    minimum_contact =  enefunc%minimum_contact

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,enefunc%table%num_solute))
               !force_j(3,size(enefunc%nb14_calc_list)))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp


    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef, lj6, lj12, term_lj12, term_lj6,       &
    !$omp         term_elec, table, force_local)                               &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
   my_id = my_city_rank * nthread + id

    do ii = 1, natom - 1

      force_local(1:3) = 0.0_wp

      i = enefunc%table%solute_list(ii)

      if (mod(ii-1,nproc_city*nthread) == my_id) then

        do k = 1, num_nb14_calc(ii)

          j = nb14_calc_list(k,ii)

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2     = max(rij2, minimum_contact)
            rij2 = cutoff2 * enefunc%table%density / rij2
            L = int(rij2)
            R = rij2 - L
            lj6  = enefunc%nb14_lj6(atmcls(i),atmcls(j))
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j))

            table(1:6) = table_ene(3*L-2:3*L+3)
            term_lj12 = table(1) + R*(table(4)-table(1))
            term_lj6  = table(2) + R*(table(5)-table(2))
            term_elec = table(3) + R*(table(6)-table(3))
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec

            table(1:6) = table_grad(3*L-2:3*L+3)
            term_lj12 = table(1) + R*(table(4)-table(1))
            term_lj6  = table(2) + R*(table(5)-table(2))
            term_elec = table(3) + R*(table(6)-table(3))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + charge(i)*charge(j)*term_elec

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k,ii) = work(1:3)

            ! virial
            !
            do m = 1, 3
              virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
            end do

          end if

        end do

        force_i(1:3,ii) = force_local(1:3)

      end if

    end do

    num_nb14_all  = num_nb14
    !$omp end parallel

    do ii = 1, natom-1
      i = enefunc%table%solute_list(ii)
      force(1,i,1) = force(1,i,1) + force_i(1,ii)
      force(2,i,1) = force(2,i,1) + force_i(2,ii)
      force(3,i,1) = force(3,i,1) + force_i(3,ii)
    end do

    do ii = 1, natom-1
      do k = 1, num_nb14_calc(ii)
        j = nb14_calc_list(k,ii)
        force(1,j,1) = force(1,j,1) + force_j(1,k,ii)
        force(2,j,1) = force(2,j,1) + force_j(2,k,ii)
        force(3,j,1) = force(3,j,1) + force_j(3,k,ii)
      end do
    end do

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber_check
  !> @brief        Calculate PME nonbond 1-4 term real part + Lennard-Jones
  !!               with linear 1/R**2 table
  !! @authors      !TODO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber_check( &
                                                     enefunc, molecule, coord, &
                                                     force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, table(6)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: qq_scale, lj_scale
    real(wp)                  :: inv_r12, inv_r121, inv_r123, inv_r126,inv_r1212
    real(wp)                  :: minimum_contact
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    real(wp),         pointer :: nb14_qq_scale(:,:), nb14_lj_scale(:,:)
    integer,          pointer :: natom
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom          => enefunc%table%num_solute
    charge         => molecule%charge
    atmcls         => molecule%atom_cls_no
    cutoff         => enefunc%cutoffdist
    num_nb14_calc  => enefunc%num_nb14_calc
    nb14_calc_list => enefunc%nb14_calc_list
    nb14_qq_scale  => enefunc%nb14_qq_scale
    nb14_lj_scale  => enefunc%nb14_lj_scale
    table_ene      => enefunc%table%table_ene
    table_grad     => enefunc%table%table_grad

    minimum_contact =  enefunc%minimum_contact

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,natom))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp


    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef, lj6, lj12, term_lj12, term_lj6,       &
    !$omp         term_elec, table, force_local, qq_scale, lj_scale,           &
    !$omp         inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212)            &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
   my_id = my_city_rank * nthread + id

    do ii = 1, natom - 1

      force_local(1:3) = 0.0_wp

      i = enefunc%table%solute_list(ii)

      if (mod(ii-1,nproc_city*nthread) == my_id) then

        do k = 1, num_nb14_calc(ii)

          j        = nb14_calc_list(k,ii)
          qq_scale = nb14_qq_scale(k,ii)
          lj_scale = nb14_lj_scale(k,ii)

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = max(rij2, minimum_contact)
          inv_r12  = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2 = cutoff2 * enefunc%table%density / rij2
            L = int(rij2)
            R = rij2 - L

            lj6  = enefunc%nb14_lj6 (atmcls(i),atmcls(j)) * lj_scale
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j)) * lj_scale

            table(1:6) = table_ene(3*L-2:3*L+3)
            term_lj12 = inv_r1212 
            term_lj6  = inv_r126   
            term_elec = table(3) + R*(table(6)-table(3))
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec * qq_scale

            table(1:6) = table_grad(3*L-2:3*L+3)
            term_lj12 = -12.0_wp * inv_r1212 * inv_r12
            term_lj6  = -6.0_wp * inv_r126 * inv_r12
            term_elec = table(3) + R*(table(6)-table(3))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                      + charge(i)*charge(j)*term_elec * qq_scale

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k,ii) = work(1:3)

            ! virial
            !
            do m = 1, 3
              virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
            end do

          end if

        end do

        force_i(1:3,ii) = force_local(1:3)

      end if

    end do

    num_nb14_all  = num_nb14
    !$omp end parallel

    do ii = 1, natom-1
      i = enefunc%table%solute_list(ii)
      force(1,i,1) = force(1,i,1) + force_i(1,ii)
      force(2,i,1) = force(2,i,1) + force_i(2,ii)
      force(3,i,1) = force(3,i,1) + force_i(3,ii)
    end do

    do ii = 1, natom-1
      do k = 1, num_nb14_calc(ii)
        j = nb14_calc_list(k,ii)
        force(1,j,1) = force(1,j,1) + force_j(1,k,ii)
        force(2,j,1) = force(2,j,1) + force_j(2,k,ii)
        force(3,j,1) = force(3,j,1) + force_j(3,k,ii)
      end do
    end do

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_solute_linear
  !> @brief        Calculate PME real part + Lennard-Jones with lookup table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_solute_linear(enefunc, molecule,     &
                                             boundary, pairlist, coord, force, &
                                             virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, qtmp, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: box_size(1:3)
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calcw(:), nb15_calc_listw(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nsolute         => enefunc%table%num_solute
    charge          => molecule%charge
    atmcls          => molecule%atom_cls_no
    num_nb15_calc   => pairlist%table%num_nb15_calc
    num_nb15_calcw  => pairlist%table%num_nb15_calcw
    nb15_calc_list  => pairlist%table%nb15_calc_list
    nb15_calc_listw => pairlist%table%nb15_calc_listw
    solute_list     => enefunc%table%solute_list
    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    box_size(1)     =  boundary%box_size_x
    box_size(2)     =  boundary%box_size_y
    box_size(3)     =  boundary%box_size_z

    cutoff2  = cutoff * cutoff

    ! calculate energy and gradient
    !
    !$omp parallel                                                    &
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, ii, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         force_local)                                        &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do ii = 1, nsolute-1

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(ii)

        j = nb15_calc_list(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3) - box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! comute energy
          !
          L1 = 3*L - 2
          term_lj12 = table_ene(L1  ) + R*(table_ene(L1+3)-table_ene(L1  ))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw      = evdw + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    do ii = 1, nsolute

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calcw(ii)

        j = nb15_calc_listw(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3)-box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6 (atmcls(i),atmcls(j))


          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! compute energy
          !
          L1 = 3*L - 2
          term_lj12 = table_ene(L1  ) + R*(table_ene(L1+3)-table_ene(L1  ))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw      = evdw  + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_energy_nonbond_table_solute_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_water_linear
  !> @brief        Calculate Nonbonded water-water interaction with table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_water_linear(enefunc, molecule,     &
                                                   boundary, pairlist, coord, &
                                                   force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: rtmp(1:3,1:3)
    real(wp)                  :: force_local(1:3,1:3)
    real(wp)                  :: dij(1:3,1:3,1:3)
    real(wp)                  :: work(1:3,1:3,1:3)
    real(wp)                  :: rij2(1:3,1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: R, cutoff2
    real(wp)                  :: Delta_pbc(1:3), box_size(1:3)
    real(wp)                  :: grad_coef, term_lj, term_elec
    integer                   :: type(1:3,1:3), water_type
    integer                   :: i, j, k, L, L1
    integer                   :: i1(1:3), j1(1:3)
    integer                   :: ii, jj, ia, ja, m
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene_WW(:,:), table_de_WW(:,:)
    integer,          pointer :: nwater
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nwater         => enefunc%table%num_water
    num_nb15_calc  => pairlist%table%num_nb15_calc_water
    nb15_calc_list => pairlist%table%nb15_calc_list_water
    cutoff         => enefunc%cutoffdist
    table_ene_WW   => enefunc%table%table_ene_WW
    table_de_WW    => enefunc%table%table_de_WW
    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z

    cutoff2  = cutoff * cutoff

    type(1,1)     = 1
    type(1,2:3)   = 2
    type(2:3,1)   = 2
    type(2:3,2:3) = 3

    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb15)                                               &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, i1, j1, ia, ja, m,   &
    !$omp         ii, jj, dij, rij2, R, L, L1, term_lj, grad_coef, term_elec,  &
    !$omp         rtmp,  force_local, viri_local, work, Delta_pbc, water_type) &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do i = 1, nwater-1

      i1(1:3) = enefunc%table%water_list(1:3,i)

      do ii = 1, 3
        rtmp(1:3,ii) = coord(1:3,i1(ii))
      end do

      force_local(1:3,1:3) = 0.0_wp

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(i)

        j = nb15_calc_list(k,i)

        j1(1:3) = enefunc%table%water_list(1:3,j)

        ! compute distance
        !
        dij(1:3,1,1)   = rtmp(1:3,1) - coord(1:3,j1(1))
        Delta_pbc(1:3) = box_size(1:3)*anint(dij(1:3,1,1)/box_size(1:3))
        do ii = 1, 3
          do jj = 1, 3
            ja = j1(jj)
            dij(1:3,ii,jj) = rtmp(1:3,ii) - coord(1:3,ja) - Delta_pbc(1:3)
            rij2(ii,jj) = dij(1,ii,jj)*dij(1,ii,jj) +dij(2,ii,jj)*dij(2,ii,jj) &
                        + dij(3,ii,jj)*dij(3,ii,jj)
          end do
        end do

        viri_local(1:3,1:3) = 0.0_wp

        do ii = 1, 3
          do jj = 1, 3

            if (rij2(ii,jj) < cutoff2) then

              rij2(ii,jj) = cutoff2 * enefunc%table%density / rij2(ii,jj)

              L    = int(rij2(ii,jj))
              R    = rij2(ii,jj) - L
              water_type = type(ii,jj)

              ! compute energy and force
              !
              L1 = 6*L-5
              term_lj   =   table_ene_WW(L1,  water_type)
              term_elec =   table_ene_WW(L1+1,water_type)
              grad_coef =   table_ene_WW(L1+2,water_type)
              term_lj   = term_lj                         &
                        + R*table_ene_WW(L1+3,water_type)
              term_elec = term_elec                       &
                        + R*table_ene_WW(L1+4,water_type)
              grad_coef = grad_coef                       &
                        + R*table_ene_WW(L1+5,water_type)

              evdw  = evdw  + term_lj
              eelec = eelec + term_elec

              work(1:3,ii,jj) = grad_coef * dij(1:3,ii,jj)

              force_local(1:3,ii) = force_local(1:3,ii) - work(1:3,ii,jj)
              ja = j1(jj)
              force(1:3,ja,id+1) = force(1:3,ja,id+1) + work(1:3,ii,jj)

              do m = 1, 3
                viri_local(1:3,m) = viri_local(1:3,m)  &
                                  + dij(1:3,ii,jj)*work(m,ii,jj)
              end do

            end if

          end do
        end do

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) - viri_local(1:3,m)
        end do

      end do

      do ii = 1, 3
        ia = i1(ii)
        force(1:3,ia,id+1) = force(1:3,ia,id+1) + force_local(1:3,ii)
      end do

    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_water_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_solute_linear_check
  !> @brief        Calculate PME real part + Lennard-Jones with lookup table
  !! @authors      JJ, TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_solute_linear_check( &
                                            enefunc, molecule, boundary, &
                                            pairlist, coord, force,      &
                                            virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, qtmp, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: box_size(1:3)
    real(wp)                  :: minimum_contact
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calcw(:), nb15_calc_listw(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nsolute         => enefunc%table%num_solute
    charge          => molecule%charge
    atmcls          => molecule%atom_cls_no
    num_nb15_calc   => pairlist%table%num_nb15_calc
    num_nb15_calcw  => pairlist%table%num_nb15_calcw
    nb15_calc_list  => pairlist%table%nb15_calc_list
    nb15_calc_listw => pairlist%table%nb15_calc_listw
    solute_list     => enefunc%table%solute_list
    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    box_size(1)     =  boundary%box_size_x
    box_size(2)     =  boundary%box_size_y
    box_size(3)     =  boundary%box_size_z

    cutoff2  = cutoff * cutoff
    minimum_contact =  enefunc%minimum_contact

    ! calculate energy and gradient
    !
    !$omp parallel                                                    &
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, ii, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         force_local)                                        &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do ii = 1, nsolute - 1

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(ii)

        j = nb15_calc_list(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3) - box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2 = max(rij2, minimum_contact)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! comute energy
          !
          L1 = 3*L - 2
          term_lj12 = table_ene(L1  ) + R*(table_ene(L1+3)-table_ene(L1  ))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw      = evdw + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    do ii = 1, nsolute

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calcw(ii)

        j = nb15_calc_listw(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3) - box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rij2 = max(rij2, minimum_contact)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6 (atmcls(i),atmcls(j))


          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! compute energy
          !
          L1 = 3*L - 2
          term_lj12 = table_ene(L1  ) + R*(table_ene(L1+3)-table_ene(L1  ))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw      = evdw  + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_energy_nonbond_table_solute_linear_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_water_linear_check
  !> @brief        Calculate Nonbonded water-water interaction with table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_water_linear_check( &
                                            enefunc, molecule,         &
                                            boundary, pairlist, coord, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw

    ! local variables
    real(wp)                  :: rtmp(1:3,1:3)
    real(wp)                  :: force_local(1:3,1:3)
    real(wp)                  :: dij(1:3,1:3,1:3)
    real(wp)                  :: work(1:3,1:3,1:3)
    real(wp)                  :: rij2(1:3,1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: R, cutoff2
    real(wp)                  :: Delta_pbc(1:3), box_size(1:3)
    real(wp)                  :: grad_coef, term_lj, term_elec
    real(wp)                  :: minimum_contact
    integer                   :: type(1:3,1:3), water_type
    integer                   :: i, j, k, L, L1
    integer                   :: i1(1:3), j1(1:3)
    integer                   :: ii, jj, ia, ja, m
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene_WW(:,:), table_de_WW(:,:)
    integer,          pointer :: nwater
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nwater         => enefunc%table%num_water
    num_nb15_calc  => pairlist%table%num_nb15_calc_water
    nb15_calc_list => pairlist%table%nb15_calc_list_water
    cutoff         => enefunc%cutoffdist
    table_ene_WW   => enefunc%table%table_ene_WW
    table_de_WW    => enefunc%table%table_de_WW
    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z

    cutoff2  = cutoff * cutoff
    minimum_contact = enefunc%minimum_contact

    type(1,1)     = 1
    type(1,2:3)   = 2
    type(2:3,1)   = 2
    type(2:3,2:3) = 3

    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb15)                                               &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, i1, j1, ia, ja, m,   &
    !$omp         ii, jj, dij, rij2, R, L, L1, term_lj, grad_coef, term_elec,  &
    !$omp         rtmp,  force_local, viri_local, work, Delta_pbc, water_type) &
    !$omp reduction(+:virial) reduction(+:eelec)  reduction(+:evdw)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do i = 1, nwater-1

      i1(1:3) = enefunc%table%water_list(1:3,i)

      do ii = 1, 3
        rtmp(1:3,ii) = coord(1:3,i1(ii))
      end do

      force_local(1:3,1:3) = 0.0_wp

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(i)

        j = nb15_calc_list(k,i)

        j1(1:3) = enefunc%table%water_list(1:3,j)

        ! compute distance
        !
        dij(1:3,1,1)   = rtmp(1:3,1) - coord(1:3,j1(1))
        Delta_pbc(1:3) = box_size(1:3)*anint(dij(1:3,1,1)/box_size(1:3))
        do ii = 1, 3
          do jj = 1, 3
            ja = j1(jj)
            dij(1:3,ii,jj) = rtmp(1:3,ii) - coord(1:3,ja) - Delta_pbc(1:3)
            rij2(ii,jj) = dij(1,ii,jj)*dij(1,ii,jj) +dij(2,ii,jj)*dij(2,ii,jj) &
                        + dij(3,ii,jj)*dij(3,ii,jj)
            rij2(ii,jj) = max(rij2(ii,jj), minimum_contact)
          end do
        end do

        viri_local(1:3,1:3) = 0.0_wp

        do ii = 1, 3
          do jj = 1, 3

            if (rij2(ii,jj) < cutoff2) then

              rij2(ii,jj) = cutoff2 * enefunc%table%density / rij2(ii,jj)

              L    = int(rij2(ii,jj))
              R    = rij2(ii,jj) - L
              water_type = type(ii,jj)

              ! compute energy and force
              !
              L1 = 6*L-5
              term_lj   =   table_ene_WW(L1,  water_type)
              term_elec =   table_ene_WW(L1+1,water_type)
              grad_coef =   table_ene_WW(L1+2,water_type)
              term_lj   = term_lj                         &
                        + R*table_ene_WW(L1+3,water_type)
              term_elec = term_elec                       &
                        + R*table_ene_WW(L1+4,water_type)
              grad_coef = grad_coef                       &
                        + R*table_ene_WW(L1+5,water_type)

              evdw  = evdw  + term_lj
              eelec = eelec + term_elec

              work(1:3,ii,jj) = grad_coef * dij(1:3,ii,jj)

              force_local(1:3,ii) = force_local(1:3,ii) - work(1:3,ii,jj)
              ja = j1(jj)
              force(1:3,ja,id+1) = force(1:3,ja,id+1) + work(1:3,ii,jj)

              do m = 1, 3
                viri_local(1:3,m) = viri_local(1:3,m)  &
                                  + dij(1:3,ii,jj)*work(m,ii,jj)
              end do

            end if

          end do
        end do

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) - viri_local(1:3,m)
        end do

      end do

      do ii = 1, 3
        ia = i1(ii)
        force(1:3,ia,id+1) = force(1:3,ia,id+1) + force_local(1:3,ii)
      end do

    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_water_linear_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_solute_linear
  !> @brief        Calculate PME real part + Lennard-Jones with lookup table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_solute_linear(enefunc, molecule, &
                                       boundary, pairlist, coord, force, virial)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)

    ! local variables
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: R, qtmp, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: box_size(1:3)
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calcw(:), nb15_calc_listw(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nsolute         => enefunc%table%num_solute
    charge          => molecule%charge
    atmcls          => molecule%atom_cls_no
    num_nb15_calc   => pairlist%table%num_nb15_calc
    num_nb15_calcw  => pairlist%table%num_nb15_calcw
    nb15_calc_list  => pairlist%table%nb15_calc_list
    nb15_calc_listw => pairlist%table%nb15_calc_listw
    solute_list     => enefunc%table%solute_list
    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    box_size(1)     =  boundary%box_size_x
    box_size(2)     =  boundary%box_size_y
    box_size(3)     =  boundary%box_size_z

    cutoff2  = cutoff * cutoff

    ! calculate energy and gradient
    !
    !$omp parallel                                                    &
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, ii, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         force_local)                                        &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do ii = 1, nsolute-1

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(ii)

        j = nb15_calc_list(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3) - box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    do ii = 1, nsolute

      i = solute_list(ii)

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calcw(ii)

        j = nb15_calc_listw(k,ii)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)
        dij(1:3) = dij(1:3) - box_size(1:3)*anint(dij(1:3)/box_size(1:3))

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2  = cutoff2 * enefunc%table%density / rij2

          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6 (atmcls(i),atmcls(j))


          ! compute force
          !
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_force_nonbond_table_solute_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_water_linear
  !> @brief        Calculate Nonbonded water-water interaction with table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_water_linear(enefunc, molecule,      &
                                                   boundary, pairlist, coord, &
                                                   force, virial)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_boundary), target, intent(in)    :: boundary
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)

    ! local variables
    real(wp)                  :: rtmp(1:3,1:3)
    real(wp)                  :: force_local(1:3,1:3)
    real(wp)                  :: dij(1:3,1:3,1:3)
    real(wp)                  :: work(1:3,1:3,1:3)
    real(wp)                  :: rij2(1:3,1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: R, cutoff2
    real(wp)                  :: Delta_pbc(1:3), box_size(1:3)
    real(wp)                  :: grad_coef, term_lj, term_elec
    integer                   :: type(1:3,1:3), water_type
    integer                   :: i, j, k, L, L1
    integer                   :: i1(1:3), j1(1:3)
    integer                   :: ii, jj, ia, ja, m
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene_WW(:,:), table_de_WW(:,:)
    integer,          pointer :: nwater
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    ! use pointers
    nwater         => enefunc%table%num_water
    num_nb15_calc  => pairlist%table%num_nb15_calc_water
    nb15_calc_list => pairlist%table%nb15_calc_list_water
    cutoff         => enefunc%cutoffdist
    table_ene_WW   => enefunc%table%table_ene_WW
    table_de_WW    => enefunc%table%table_de_WW
    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z

    cutoff2  = cutoff * cutoff

    type(1,1)     = 1
    type(1,2:3)   = 2
    type(2:3,1)   = 2
    type(2:3,2:3) = 3

    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb15)                                               &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, i1, j1, ia, ja, m,   &
    !$omp         ii, jj, dij, rij2, R, L, L1, term_lj, grad_coef, term_elec,  &
    !$omp         rtmp,  force_local, viri_local, work, Delta_pbc, water_type) &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id   = my_city_rank * nthread + id

    do i = 1, nwater - 1

      i1(1:3) = enefunc%table%water_list(1:3,i)

      do ii = 1, 3
        rtmp(1:3,ii) = coord(1:3,i1(ii))
      end do

      force_local(1:3,1:3) = 0.0_wp

      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      do k = 1, num_nb15_calc(i)

        j = nb15_calc_list(k,i)

        j1(1:3) = enefunc%table%water_list(1:3,j)

        ! compute distance
        !
        dij(1:3,1,1)   = rtmp(1:3,1) - coord(1:3,j1(1))
        Delta_pbc(1:3) = box_size(1:3)*anint(dij(1:3,1,1)/box_size(1:3))
        do ii = 1, 3
          do jj = 1, 3
            ja = j1(jj)
            dij(1:3,ii,jj) = rtmp(1:3,ii) - coord(1:3,ja) - Delta_pbc(1:3)
            rij2(ii,jj) = dij(1,ii,jj)*dij(1,ii,jj) +dij(2,ii,jj)*dij(2,ii,jj) &
                        + dij(3,ii,jj)*dij(3,ii,jj)
          end do
        end do

        viri_local(1:3,1:3) = 0.0_wp

        do ii = 1, 3
          do jj = 1, 3

            if (rij2(ii,jj) < cutoff2) then

              rij2(ii,jj) = cutoff2 * enefunc%table%density / rij2(ii,jj)

              L    = int(rij2(ii,jj))
              R    = rij2(ii,jj) - L
              water_type = type(ii,jj)


              ! compute force
              !
              grad_coef = table_de_WW(L,water_type)       &
                        + R*(table_de_WW(L+1,water_type)  &
                           - table_de_WW(L,  water_type))

              work(1:3,ii,jj) = grad_coef * dij(1:3,ii,jj)

              force_local(1:3,ii) = force_local(1:3,ii) - work(1:3,ii,jj)
              ja = j1(jj)
              force(1:3,ja,id+1) = force(1:3,ja,id+1) + work(1:3,ii,jj)

              do m = 1, 3
                viri_local(1:3,m) = viri_local(1:3,m)  &
                                  + dij(1:3,ii,jj)*work(m,ii,jj)
              end do

            end if

          end do
        end do

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) - viri_local(1:3,m)
        end do

      end do

      do ii = 1, 3
        ia = i1(ii)
        force(1:3,ia,id+1) = force(1:3,ia,id+1) + force_local(1:3,ii)
      end do

    end do
    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_water_linear

end module at_energy_table_linear_mod
