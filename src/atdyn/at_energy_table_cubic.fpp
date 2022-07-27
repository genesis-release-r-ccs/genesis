!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_table_cubic_mod
!> @brief   calculate nonbonded energy with table
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_table_cubic_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_pairlist_mod
  use molecules_str_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! variables
  real(wp), save, allocatable :: force_i(:,:)
  real(wp), save, allocatable :: force_j(:,:,:)

  ! subroutines
  public  :: compute_energy_nonbond14_table
  public  :: compute_energy_nonbond_table
  public  :: compute_force_nonbond_table
  public  :: compute_energy_nonbond_nobc_table
  public  :: compute_force_nonbond_nobc_table

  private :: compute_energy_nonbond14_table_charmm
  private :: compute_energy_nonbond14_table_ecqm
  private :: compute_energy_nonbond14_table_gro_amber
  private :: compute_energy_nonbond_table_solute_cubic
  private :: compute_energy_nonbond_table_water_cubic
  private :: compute_force_nonbond_table_solute_cubic
  private :: compute_force_nonbond_table_water_cubic
  private :: compute_energy_nonbond_nobc_table_cubic
  private :: compute_energy_nonbond_nobc_table_cubic_ecqm
  private :: compute_force_nonbond_nobc_table_cubic
  private :: compute_force_nonbond_nobc_table_cubic_ecqm

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table
  !> @brief        Calculate nonbond 1-4 term real part + Lennard-Jones 
  !!               with table
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

  subroutine compute_energy_nonbond14_table(enefunc, molecule, coord, &
                                            force, virial, eelec, evdw)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw


    ! ==> Type 2/3
    if (enefunc%forcefield == ForcefieldCHARMM .or. &
        enefunc%forcefield == ForcefieldCHARMM19) then

      call compute_energy_nonbond14_table_charmm( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)
      
      if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%qm_classical    &
                               .and. enefunc%qmmm%num_qmmmbonds > 0) then
        call compute_energy_nonbond14_table_ecqm( enefunc, molecule, &
                                          coord, force, virial, eelec)
      end if

    ! ==> Type 8/9
    else  if (enefunc%forcefield == ForcefieldAMBER    .or. &
              enefunc%forcefield == ForcefieldGROAMBER .or. &
              enefunc%forcefield == ForcefieldGROMARTINI) then

      call compute_energy_nonbond14_table_gro_amber( &
                                    enefunc, molecule, &
                                    coord, force, virial, eelec, evdw)
    end if

    return
  end subroutine compute_energy_nonbond14_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table
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

  subroutine compute_energy_nonbond_table(enefunc, molecule, pairlist,    &
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
    call compute_energy_nonbond_table_solute_cubic(enefunc, molecule, &
                                boundary, pairlist, coord, force,     &
                                virial, eelec, evdw)

    ! Compute water-water interactions
    !
    call compute_energy_nonbond_table_water_cubic(enefunc, molecule,  &
                                boundary, pairlist, coord, force,     &
                                virial, eelec, evdw)

    return
 
  end subroutine compute_energy_nonbond_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table
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

  subroutine compute_force_nonbond_table(enefunc, molecule, pairlist,    &
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
    call compute_force_nonbond_table_solute_cubic(enefunc, molecule, &
                                boundary, pairlist, coord, force, virial)

    ! Compute water-water interactions
    !
    call compute_force_nonbond_table_water_cubic(enefunc, molecule,  &
                                boundary, pairlist, coord, force, virial)

    return
 
  end subroutine compute_force_nonbond_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_nobc_table
  !> @brief        Calculate nonbonded force with lookup table (NOBC)
  !! @authors      NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_nobc_table(enefunc, molecule, pairlist, &
                                               coord, force, virial,        &
                                               eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    ! Compute interactions
    !
    call compute_energy_nonbond_nobc_table_cubic(enefunc, molecule, &
                                pairlist, coord, force,     &
                                virial, eelec, evdw)
    if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%qm_classical    &
                             .and. enefunc%qmmm%num_qmmmbonds > 0) then
      call compute_energy_nonbond_nobc_table_cubic_ecqm( &
             enefunc, molecule, pairlist, coord, force, virial, eelec)
    end if

    return
 
  end subroutine compute_energy_nonbond_nobc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_nobc_table
  !> @brief        Calculate nonbonded energy with lookup table (NOBC)
  !! @authors      NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_nobc_table(enefunc, molecule, pairlist, &
                                              coord, force, virial)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)


    ! Compute interactions
    !
    call compute_force_nonbond_nobc_table_cubic(enefunc, molecule, &
                                pairlist, coord, force, virial)
    if (enefunc%qmmm%do_qmmm .and. enefunc%qmmm%qm_classical    &
                             .and. enefunc%qmmm%num_qmmmbonds > 0) then
      call compute_force_nonbond_nobc_table_cubic_ecqm( &
             enefunc, molecule, pairlist, coord, force, virial)
    end if

    return
 
  end subroutine compute_force_nonbond_nobc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_charmm
  !> @brief        Calculate nonbond 1-4 term real part + Lennard-Jones 
  !!               with table
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

  subroutine compute_energy_nonbond14_table_charmm(enefunc, molecule, coord, &
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
    real(wp)                  :: work(1:3), force_local(1:3)
    real(wp)                  :: R, h00, h01, h10, h11, table(12)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    real(wp),         pointer :: nb14_qq_scale_c19
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

    nb14_qq_scale_c19  => enefunc%nb14_qq_scale_c19

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(enefunc%num_nb14_calc(1:natom)))

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
    !$omp         term_elec, table, h00, h10, h01, h11, force_local)           &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw) 
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, natom - 1

      i = enefunc%table%solute_list(ii)

      force_local(1:3) = 0.0_wp

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

            rij2 = rij2 * enefunc%table%density 
            L    = int(rij2)
            R    = rij2 - L
            lj6  = enefunc%nb14_lj6(atmcls(i),atmcls(j))
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j))

            h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
            h10  = R*(1.0_wp-R)*(1.0_wp-R)
            h01  = R*R*(3.0_wp-2.0_wp*R)
            h11  = R*R*(R-1.0_wp)

            table(1:12) = table_ene(6*L-5:6*L+6)
            term_lj12 = table(1)*h00 + table(2)*h10 + table(7 )*h01 + table(8 )*h11
            term_lj6  = table(3)*h00 + table(4)*h10 + table(9 )*h01 + table(10)*h11
            term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 + table(12)*h11
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec*nb14_qq_scale_c19

            table(1:12) = table_grad(6*L-5:6*L+6)
            term_lj12 = table(1)*h00 + table(2)*h10 + table(7 )*h01 + table(8 )*h11
            term_lj6  = table(3)*h00 + table(4)*h10 + table(9 )*h01 + table(10)*h11
            term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 + table(12)*h11
            grad_coef = term_lj12*lj12 - term_lj6*lj6 +             &
                        charge(i)*charge(j)*term_elec*nb14_qq_scale_c19

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3)  = force_local(1:3) - work(1:3)
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

  end subroutine compute_energy_nonbond14_table_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_ecqm
  !> @brief        Remove Coulomb between QM-EC 1-4 terms
  !! @authors      KY
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_ecqm(enefunc, molecule, coord, &
                                                 force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    ! local variables
    integer,         pointer :: ecqm_num_nb14
    integer,         pointer :: ecqm_nb14_list(:,:)

    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: term_elec
    real(wp)                  :: R, h00, h01, h10, h11, table(12)
    integer                   :: ii, i, j, k, m, L
    integer                   :: id, my_id

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)

    integer                   :: qm_natoms, ec_natoms
    integer,          pointer :: qmatom_id(:), ecatom_id(:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    ecqm_num_nb14   => enefunc%table%ecqm_num_nb14
    ecqm_nb14_list  => enefunc%table%ecqm_nb14_list

    charge         => molecule%charge
    cutoff         => enefunc%cutoffdist
    table_ene      => enefunc%table%table_ene
    table_grad     => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff

    !$omp parallel                                                        &
    !$omp private(id, my_id, ii, i, j, m, dij, rij2, R, L, work,          &
    !$omp         grad_coef, term_elec, table, h00, h10, h01, h11)        &
    !$omp         reduction(+:virial) reduction(+:eelec) 
     
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, ecqm_num_nb14

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      i = ecqm_nb14_list(1,ii)
      j = ecqm_nb14_list(2,ii)

      ! compute distance
      !
      dij(1:3) = coord(1:3,i) - coord(1:3,j)
      rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      ! cutoff
      !
      if (rij2 < cutoff2) then

        rij2 = rij2 * enefunc%table%density 
        L    = int(rij2)
        R    = rij2 - L

        h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10  = R*(1.0_wp-R)*(1.0_wp-R)
        h01  = R*R*(3.0_wp-2.0_wp*R)
        h11  = R*R*(R-1.0_wp)

        table(1:12) = table_ene(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 + table(12)*h11
        eelec = eelec - charge(i)*charge(j)*term_elec

        table(1:12) = table_grad(6*L-5:6*L+6)
        term_elec = table(5)*h00 + table(6)*h10 + table(11)*h01 + table(12)*h11
        grad_coef = charge(i)*charge(j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        force(1:3,i,id+1) = force(1:3,i,id+1) + work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) - work(1:3)

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) + dij(1:3)*work(m)
        end do

      end if
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_ecqm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_gro_amber
  !> @brief        Calculate nonbond 1-4 term real part + Lennard-Jones 
  !!               with table
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

  subroutine compute_energy_nonbond14_table_gro_amber(enefunc, molecule, coord,&
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
    real(wp)                  :: work(1:3), force_local(1:3)
    real(wp)                  :: R, h00, h01, h10, h11, table(12)
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: qq_scale, lj_scale
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
    !$omp         term_elec, table, h00, h10, h01, h11, force_local, qq_scale, &
    !$omp         lj_scale)                                                    &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw) 
    !
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, natom-1

      i = enefunc%table%solute_list(ii)

      force_local(1:3) = 0.0_wp

      if (mod(ii-1,nproc_city*nthread) == my_id) then

        do k = 1, num_nb14_calc(ii)

          j        = nb14_calc_list(k,ii)
          qq_scale = nb14_qq_scale(k,ii)
          lj_scale = nb14_lj_scale(k,ii)

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2 = rij2 * enefunc%table%density 
            L    = int(rij2)
            R    = rij2 - L

            lj6  = enefunc%nb14_lj6 (atmcls(i),atmcls(j)) * lj_scale
            lj12 = enefunc%nb14_lj12(atmcls(i),atmcls(j)) * lj_scale

            h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
            h10  = R*(1.0_wp-R)*(1.0_wp-R)
            h01  = R*R*(3.0_wp-2.0_wp*R)
            h11  = R*R*(R-1.0_wp)

            table(1:12) = table_ene(6*L-5:6*L+6)
            term_lj12 = table(1)*h00+ table(2)*h10+ table(7 )*h01+ table(8 )*h11
            term_lj6  = table(3)*h00+ table(4)*h10+ table(9 )*h01+ table(10)*h11
            term_elec = table(5)*h00+ table(6)*h10+ table(11)*h01+ table(12)*h11
            evdw = evdw + term_lj12*lj12 - term_lj6*lj6
            eelec = eelec + charge(i)*charge(j)*term_elec * qq_scale

            table(1:12) = table_grad(6*L-5:6*L+6)
            term_lj12 = table(1)*h00+ table(2)*h10+ table(7 )*h01+ table(8 )*h11
            term_lj6  = table(3)*h00+ table(4)*h10+ table(9 )*h01+ table(10)*h11
            term_elec = table(5)*h00+ table(6)*h10+ table(11)*h01+ table(12)*h11
            grad_coef = term_lj12*lj12 - term_lj6*lj6 + &
                            charge(i)*charge(j)*term_elec * qq_scale

            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k,ii)   = work(1:3)

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

  end subroutine compute_energy_nonbond14_table_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_solute_cubic
  !> @brief        Calculate real part + Lennard-Jones with lookup table
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

  subroutine compute_energy_nonbond_table_solute_cubic(enefunc, molecule,     &
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
    real(wp)                  :: h00, h01, h10, h11
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    real(wp),         pointer :: bsize_x, bsize_y, bsize_z
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calcw(:), nb15_calc_listw(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


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
    num_nb15 = 0

    ! calculate energy and gradient
    !
    !$omp parallel                                                    & 
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, ii, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         h00, h01, h10, h11, force_local)                    &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
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
      qtmp = charge(i)
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

          rij2 = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute force
          !
          L1 = 6*L - 5
          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! compute energy
          !
          L1 = 6*L - 5
          term_lj12 = table_ene(L1  )*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_lj12 = term_lj12 + table_ene(L1+6 )*h01 + table_ene(L1+7 )*h11
          term_lj6  = term_lj6  + table_ene(L1+8 )*h01 + table_ene(L1+9 )*h11
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          evdw      = evdw  + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    ! solute-water
    !
    num_nb15 = 0
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

          rij2  = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute energy and force
          !
          L1 = 6*L - 5
          term_lj12 = table_ene(L1  )*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_lj12 = term_lj12 + table_ene(L1+6 )*h01 + table_ene(L1+7 )*h11
          term_lj6  = term_lj6  + table_ene(L1+8 )*h01 + table_ene(L1+9 )*h11
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          evdw      = evdw  + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_energy_nonbond_table_solute_cubic
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_water_cubic
  !> @brief        Calculate nonbonded water-water interaction with table
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

  subroutine compute_energy_nonbond_table_water_cubic(enefunc, molecule,     &
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
    real(wp)                  :: h00, h01, h10, h11
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


    nwater         => enefunc%table%num_water
    cutoff         => enefunc%cutoffdist
    table_ene_WW   => enefunc%table%table_ene_WW
    table_de_WW    => enefunc%table%table_de_WW
    num_nb15_calc  => pairlist%table%num_nb15_calc_water
    nb15_calc_list => pairlist%table%nb15_calc_list_water

    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z

    cutoff2  = cutoff * cutoff
    num_nb15 = 0

    type(1,1)     = 1
    type(1,2:3)   = 2
    type(2:3,1)   = 2
    type(2:3,2:3) = 3

    ! calculate energy and gradient  water-water
    !
    !$omp parallel                                                            & 
    !$omp firstprivate(num_nb15)                                              &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, i1, j1, ia, ja,     &
    !$omp         ii, jj, dij, rij2, R, L, L1, term_lj, grad_coef, term_elec, &
    !$omp         rtmp,  force_local, viri_local, work, Delta_pbc,            &
    !$omp         water_type, m, h00, h01, h10, h11)                          &
    !$omp reduction(+:virial) reduction(+:eelec)  reduction(+:evdw) 
    !
#ifdef OMP
    id    = omp_get_thread_num()
#else
    id    = 0
#endif
    my_id = my_city_rank * nthread + id

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
            rij2(ii,jj) = dij(1,ii,jj)*dij(1,ii,jj) + dij(2,ii,jj)*dij(2,ii,jj) &
                        + dij(3,ii,jj)*dij(3,ii,jj)
          end do
        end do

        viri_local(1:3,1:3) = 0.0_wp

        do ii = 1, 3
          do jj = 1, 3

            if (rij2(ii,jj) < cutoff2) then

              rij2(ii,jj) = rij2(ii,jj) * enefunc%table%density 

              L   = int(rij2(ii,jj))
              R   = rij2(ii,jj) - L
              water_type = type(ii,jj)
              h00 = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
              h10 = R*(1.0_wp-R)*(1.0_wp-R)
              h01 = R*R*(3.0_wp-2.0_wp*R)
              h11 = R*R*(R-1.0_wp)

              ! compute energy and force
              !
              L1 = 6*L - 5
              term_lj   = table_ene_WW(L1,   water_type)*h00 &
                        + table_ene_WW(L1+1, water_type)*h10
              term_elec = table_ene_WW(L1+2, water_type)*h00 &
                        + table_ene_WW(L1+3, water_type)*h10
              grad_coef = table_ene_WW(L1+4, water_type)*h00 &
                        + table_ene_WW(L1+5, water_type)*h10
              term_lj   = term_lj                            &
                        + table_ene_WW(L1+6, water_type)*h01 &
                        + table_ene_WW(L1+7, water_type)*h11
              term_elec = term_elec                          &
                        + table_ene_WW(L1+8, water_type)*h01 &
                        + table_ene_WW(L1+9, water_type)*h11
              grad_coef = grad_coef                          &
                        + table_ene_WW(L1+10,water_type)*h01 &
                        + table_ene_WW(L1+11,water_type)*h11

              evdw  = evdw  + term_lj
              eelec = eelec + term_elec

              work(1:3,ii,jj) = grad_coef * dij(1:3,ii,jj)

              ! store force
              !
              force_local(1:3,ii) = force_local(1:3,ii) - work(1:3,ii,jj)
              ja = j1(jj)
              force(1:3,ja,id+1) = force(1:3,ja,id+1) + work(1:3,ii,jj)

              ! virial
              !
              do m = 1, 3
                viri_local(1:3,m) = viri_local(1:3,m)  &
                                  + dij(1:3,ii,jj)*work(m,ii,jj)
              end do

            end if

          end do
        end do

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

  end subroutine compute_energy_nonbond_table_water_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_solute_cubic
  !> @brief        Calculate real part + Lennard-Jones with lookup table
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

  subroutine compute_force_nonbond_table_solute_cubic(enefunc, molecule, &
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
    real(wp)                  :: h00, h01, h10, h11
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    real(wp),         pointer :: bsize_x, bsize_y, bsize_z
    integer,          pointer :: nsolute, atmcls(:), solute_list(:)
    integer,          pointer :: num_nb15_calc(:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calcw(:), nb15_calc_listw(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


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
    num_nb15 = 0

    ! calculate energy and gradient
    !
    !$omp parallel                                                    & 
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, ii, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         h00, h01, h10, h11, force_local)                    &
    !$omp reduction(+:virial)
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
      qtmp = charge(i)
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

          rij2 = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute force
          !
          L1 = 6*L - 5
          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    ! solute-water
    !
    num_nb15 = 0
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

          rij2  = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute force only
          !
          L1 = 6*L-5
          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
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

  end subroutine compute_force_nonbond_table_solute_cubic
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_water_cubic
  !> @brief        Calculate nonbonded water-water interaction with table
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

  subroutine compute_force_nonbond_table_water_cubic(enefunc, molecule, &
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
    real(wp)                  :: rtmp(1:3,1:3)
    real(wp)                  :: force_local(1:3,1:3)
    real(wp)                  :: dij(1:3,1:3,1:3)
    real(wp)                  :: work(1:3,1:3,1:3)
    real(wp)                  :: rij2(1:3,1:3)
    real(wp)                  :: viri_local(1:3,1:3)
    real(wp)                  :: R, cutoff2
    real(wp)                  :: Delta_pbc(1:3), box_size(1:3)
    real(wp)                  :: grad_coef, term_lj, term_elec
    real(wp)                  :: h00, h01, h10, h11
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


    nwater         => enefunc%table%num_water
    cutoff         => enefunc%cutoffdist
    table_ene_WW   => enefunc%table%table_ene_WW
    table_de_WW    => enefunc%table%table_de_WW
    num_nb15_calc  => pairlist%table%num_nb15_calc_water
    nb15_calc_list => pairlist%table%nb15_calc_list_water

    box_size(1)    =  boundary%box_size_x
    box_size(2)    =  boundary%box_size_y
    box_size(3)    =  boundary%box_size_z

    cutoff2  = cutoff * cutoff
    num_nb15 = 0

    type(1,1)     = 1
    type(1,2:3)   = 2
    type(2:3,1)   = 2
    type(2:3,2:3) = 3

    ! calculate energy and gradient  water-water
    !
    !$omp parallel                                                            & 
    !$omp firstprivate(num_nb15)                                              &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, i1, j1, ia, ja,     &
    !$omp         ii, jj, dij, rij2, R, L, L1, term_lj, grad_coef, term_elec, &
    !$omp         rtmp,  force_local, viri_local, work, Delta_pbc,            &
    !$omp         water_type, m, h00, h01, h10, h11)                          &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id    = omp_get_thread_num()
#else
    id    = 0
#endif
    my_id = my_city_rank * nthread + id

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
            rij2(ii,jj) = dij(1,ii,jj)*dij(1,ii,jj) + dij(2,ii,jj)*dij(2,ii,jj) &
                        + dij(3,ii,jj)*dij(3,ii,jj)
          end do
        end do

        viri_local(1:3,1:3) = 0.0_wp

        do ii = 1, 3
          do jj = 1, 3

            if (rij2(ii,jj) < cutoff2) then

              rij2(ii,jj) = rij2(ii,jj) * enefunc%table%density 

              L   = int(rij2(ii,jj))
              R   = rij2(ii,jj) - L
              water_type = type(ii,jj)
              h00 = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
              h10 = R*(1.0_wp-R)*(1.0_wp-R)
              h01 = R*R*(3.0_wp-2.0_wp*R)
              h11 = R*R*(R-1.0_wp)

              L  = 2*L - 1
              grad_coef = table_de_WW(L,   water_type)*h00   &
                        + table_de_WW(L+1, water_type)*h10   &
                        + table_de_WW(L+2, water_type)*h01   &
                        + table_de_WW(L+3, water_type)*h11

              work(1:3,ii,jj) = grad_coef * dij(1:3,ii,jj)

              ! store force
              !
              force_local(1:3,ii) = force_local(1:3,ii) - work(1:3,ii,jj)
              ja = j1(jj)
              force(1:3,ja,id+1) = force(1:3,ja,id+1) + work(1:3,ii,jj)

              ! virial
              !
              do m = 1, 3
                viri_local(1:3,m) = viri_local(1:3,m)  &
                                  + dij(1:3,ii,jj)*work(m,ii,jj)
              end do

            end if

          end do
        end do

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

  end subroutine compute_force_nonbond_table_water_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_nobc_table_cubic
  !> @brief        Calculate real part + Lennard-Jones with lookup table
  !! @authors      JJ, NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information 
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_nobc_table_cubic(enefunc, molecule,&
                                            pairlist, coord, force, &
                                            virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
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
    real(wp)                  :: h00, h01, h10, h11
    integer                   :: i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: natom, atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    natom           => molecule%num_atoms
    charge          => molecule%charge
    atmcls          => molecule%atom_cls_no
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list

    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff

    num_nb15 = 0

    ! calculate energy and gradient
    !
    !$omp parallel                                                    & 
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         h00, h01, h10, h11, force_local)                    &
    !$omp reduction(+:virial) reduction(+:eelec) reduction(+:evdw)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(i-1,nproc_city*nthread) /= my_id) ini_nb15 = fin_nb15 + 1

      do k = ini_nb15, fin_nb15

        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2 = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute force
           
          L1 = 6*L - 5
          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          ! compute energy
           
          L1 = 6*L - 5
          term_lj12 = table_ene(L1  )*h00 + table_ene(L1+1)*h10
          term_lj6  = table_ene(L1+2)*h00 + table_ene(L1+3)*h10
          term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
          term_lj12 = term_lj12 + table_ene(L1+6 )*h01 + table_ene(L1+7 )*h11
          term_lj6  = term_lj6  + table_ene(L1+8 )*h01 + table_ene(L1+9 )*h11
          term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
          evdw      = evdw  + term_lj12*lj12 - term_lj6*lj6
          eelec     = eelec + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_energy_nonbond_nobc_table_cubic
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_nobc_table_cubic_ecqm
  !> @brief        Remove Coulomb between QM-EC 1-5 terms
  !! @authors      KY
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_nobc_table_cubic_ecqm(enefunc, molecule,&
                                         pairlist, coord, force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    ! local variables
    integer,         pointer :: ecqm_num_nb15
    integer,         pointer :: ecqm_nb15_list(:,:)

    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: term_elec
    real(wp)                  :: R, h00, h01, h10, h11
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: id, my_id
    integer                   :: qm_natoms, ec_natoms

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: qmatom_id(:), ecatom_id(:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    ecqm_num_nb15   => pairlist%ecqm_num_nb15
    ecqm_nb15_list  => pairlist%ecqm_nb15_list

    charge          => molecule%charge
    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff

    !$omp parallel                                                        &
    !$omp private(id, my_id, ii, i, j, m, dij, rij2, R, L, L1, work,      &
    !$omp         grad_coef, term_elec, h00, h10, h01, h11)               &
    !$omp reduction(+:virial) reduction(+:eelec) 
     
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, ecqm_num_nb15

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      i = ecqm_nb15_list(1,ii)
      j = ecqm_nb15_list(2,ii)

      ! compute distance
      !
      dij(1:3) = coord(1:3,i) - coord(1:3,j)

      rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      if (rij2 < cutoff2) then

        rij2 = rij2 * enefunc%table%density
        L    = int(rij2)
        R    = rij2 - L

        h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10  = R*(1.0_wp-R)*(1.0_wp-R)
        h01  = R*R*(3.0_wp-2.0_wp*R)
        h11  = R*R*(R-1.0_wp)

        ! compute force

        L1 = 6*L - 5
        term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
        term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
        grad_coef = charge(i)*charge(j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        ! compute energy

        L1 = 6*L - 5
        term_elec = table_ene(L1+4)*h00 + table_ene(L1+5)*h10
        term_elec = term_elec + table_ene(L1+10)*h01 + table_ene(L1+11)*h11
        eelec     = eelec - charge(i)*charge(j)*term_elec

        ! store force
        !
        force(1:3,i,id+1) = force(1:3,i,id+1) + work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) - work(1:3)

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
        end do

      end if
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_nobc_table_cubic_ecqm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_nobc_table_cubic
  !> @brief        Calculate real part + Lennard-Jones with lookup table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information 
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_nobc_table_cubic(enefunc, molecule, &
                                       pairlist, coord, force, virial)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
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
    real(wp)                  :: h00, h01, h10, h11
    integer                   ::  i, j, k, m, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: natom, atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    natom           => molecule%num_atoms
    charge          => molecule%charge
    atmcls          => molecule%atom_cls_no
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list  => pairlist%nb15_calc_list

    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff
    num_nb15 = 0

    ! calculate energy and gradient
    !
    !$omp parallel                                                    & 
    !$omp firstprivate(num_nb15)                                      &
    !$omp private(id, my_id, ini_nb15, fin_nb15, i, k, j, m, dij, &
    !$omp         rij2, R, L, L1, work, grad_coef, rtmp, qtmp,        &
    !$omp         term_lj12, term_lj6, term_elec, lj12, lj6,          &
    !$omp         h00, h01, h10, h11, force_local)                    &
    !$omp reduction(+:virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    ! solute-solute
    !
    do i = 1, natom - 1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      rtmp(1:3) = coord(1:3,i)
      qtmp = charge(i)
      force_local(1:3) = 0.0_wp

      if (mod(i-1,nproc_city*nthread) /= my_id) ini_nb15 = fin_nb15 + 1

      do k = ini_nb15, fin_nb15

        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = rtmp(1:3) - coord(1:3,j)

        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        if (rij2 < cutoff2) then

          rij2 = rij2 * enefunc%table%density 
          L    = int(rij2)
          R    = rij2 - L
          lj12 = enefunc%nonb_lj12(atmcls(i),atmcls(j))
          lj6  = enefunc%nonb_lj6(atmcls(i),atmcls(j))

          h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
          h10  = R*(1.0_wp-R)*(1.0_wp-R)
          h01  = R*R*(3.0_wp-2.0_wp*R)
          h11  = R*R*(R-1.0_wp)

          ! compute force
          !
          L1 = 6*L - 5
          term_lj12 = table_grad(L1  )*h00 + table_grad(L1+1)*h10
          term_lj6  = table_grad(L1+2)*h00 + table_grad(L1+3)*h10
          term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
          term_lj12 = term_lj12 + table_grad(L1+6 )*h01 + table_grad(L1+7 )*h11
          term_lj6  = term_lj6  + table_grad(L1+8 )*h01 + table_grad(L1+9 )*h11
          term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

    return

  end subroutine compute_force_nonbond_nobc_table_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_nobc_table_cubic_ecqm
  !> @brief        Remove Coulomb between QM-EC 1-5 terms
  !! @authors      KY
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_nobc_table_cubic_ecqm(enefunc, molecule,&
                                         pairlist, coord, force, virial)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)

    ! local variables
    integer,         pointer :: ecqm_num_nb15
    integer,         pointer :: ecqm_nb15_list(:,:)

    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: grad_coef
    real(wp)                  :: work(3)
    real(wp)                  :: term_elec
    real(wp)                  :: R, h00, h01, h10, h11
    integer                   :: ii, i, j, k, m, L, L1
    integer                   :: id, my_id
    integer                   :: qm_natoms, ec_natoms

    real(wp),         pointer :: cutoff, charge(:)
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: qmatom_id(:), ecatom_id(:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif


    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    ecqm_num_nb15   => pairlist%ecqm_num_nb15
    ecqm_nb15_list  => pairlist%ecqm_nb15_list

    charge          => molecule%charge
    cutoff          => enefunc%cutoffdist
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad

    cutoff2  = cutoff * cutoff

    !$omp parallel                                                        &
    !$omp private(id, my_id, ii, i, j, m, dij, rij2, R, L, L1, work,      &
    !$omp         grad_coef, term_elec, h00, h10, h01, h11)               &
    !$omp reduction(+:virial) 
     
#ifdef OMP
    id      = omp_get_thread_num()
#else
    id      = 0
#endif
    my_id = my_city_rank * nthread + id

    do ii = 1, ecqm_num_nb15

      if (mod(ii-1,nproc_city*nthread) /= my_id) cycle

      i = ecqm_nb15_list(1,ii)
      j = ecqm_nb15_list(2,ii)

      ! compute distance
      !
      dij(1:3) = coord(1:3,i) - coord(1:3,j)

      rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      if (rij2 < cutoff2) then

        rij2 = rij2 * enefunc%table%density
        L    = int(rij2)
        R    = rij2 - L

        h00  = (1.0_wp + 2.0_wp*R)*(1.0_wp-R)*(1.0_wp-R)
        h10  = R*(1.0_wp-R)*(1.0_wp-R)
        h01  = R*R*(3.0_wp-2.0_wp*R)
        h11  = R*R*(R-1.0_wp)

        ! compute force

        L1 = 6*L - 5
        term_elec = table_grad(L1+4)*h00 + table_grad(L1+5)*h10
        term_elec = term_elec + table_grad(L1+10)*h01 + table_grad(L1+11)*h11
        grad_coef = charge(i)*charge(j)*term_elec

        work(1:3) = grad_coef*dij(1:3)

        ! store force
        !
        force(1:3,i,id+1) = force(1:3,i,id+1) + work(1:3)
        force(1:3,j,id+1) = force(1:3,j,id+1) - work(1:3)

        ! virial
        !
        do m = 1, 3
          virial(1:3,m) = virial(1:3,m) - dij(1:3)*work(m)
        end do

      end if
    end do
    !$omp end parallel

    return

  end subroutine compute_force_nonbond_nobc_table_cubic_ecqm

end module at_energy_table_cubic_mod
