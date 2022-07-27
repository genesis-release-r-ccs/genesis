!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_table_linear_bondcorr_mod
!> @brief   Smooth particle mesh ewald method
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_table_linear_bondcorr_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use at_boundary_mod
  use math_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! variables
  integer,  public, parameter :: PmeMaxNspline = 10
  real(wp), save, allocatable :: force_i(:,:)
  real(wp), save, allocatable :: force_j(:,:,:)

  ! subroutines
  public   :: pme_bond_corr_linear
  private  :: pme_bond_corr_linear_general
  private  :: pme_bond_corr_linear_gro_amber
  private  :: pme_bond_corr_linear_general_check
  private  :: pme_bond_corr_linear_gro_amber_check

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear
  !> @brief        Calculate PME bonded energy correction with lookup table
  !! @authors      JJ, TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear(enefunc, molecule, coord, force, virial,eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec


    if (enefunc%nonb_limiter) then
      call pme_bond_corr_linear_general_check(enefunc, molecule, coord,       &
                                          force, virial, eelec)

      if (enefunc%forcefield /= ForcefieldCHARMM) then
   
        call pme_bond_corr_linear_gro_amber_check(enefunc, molecule, coord,   &
                                            force,  virial, eelec)
   
      end if
    else
      call pme_bond_corr_linear_general(enefunc, molecule, coord, force,      &
                                        virial, eelec)
     
      if (enefunc%forcefield /= ForcefieldCHARMM) then
     
        call pme_bond_corr_linear_gro_amber(enefunc, molecule, coord, force,  &
                                            virial, eelec)
     
      end if
    end if

    return

  end subroutine pme_bond_corr_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general
  !> @brief        Calculate PME bonded energy correction with lookup table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general(enefunc, molecule, coord, force,  &
                                          virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                  :: u, el_fact, alpha, alpha2m, alpha2sp
    real(wp)                  :: rij(3), rij2, qij, uij, fij(3)
    real(wp)                  :: vxx,vyx,vzx,vyy,vzy,vzz
    real(wp)                  :: Res, cutoff2
    integer                   :: i, j, k, jatm
    integer                   :: num_excl,ini_excl,fin_excl
    integer                   :: iatm, L

    real(wp),         pointer :: q(:), force_wk(:,:)
    real(wp),         pointer :: table_elec(:), table_delec(:)
    integer,          pointer :: natm, num_nonb_excl(:), nonb_excl_list(:)
    integer,          pointer :: i_list(:)


    ! use pointers
    !
    natm           => molecule%num_atoms
    q              => molecule%charge
    num_nonb_excl  => enefunc%num_nonb_excl
    nonb_excl_list => enefunc%nonb_excl_list
    table_elec     => enefunc%table%table_ecor
    table_delec    => enefunc%table%table_decor
    force_wk       => enefunc%pme%force_wk
    i_list         => enefunc%pme%i_list
    el_fact        =  ELECOEF/enefunc%dielec_const
    alpha          =  enefunc%pme%alpha
    alpha2m        =  enefunc%pme%alpha2m
    alpha2sp       =  enefunc%pme%alpha2sp
    cutoff2        =  enefunc%cutoffdist * enefunc%cutoffdist

    ! Initializing the energy and force 
    !
    u        = 0.0_wp
    num_excl = 0
    i_list(:)     = 0
    force_wk(:,:) = 0.0_wp

    ! Running only the 1-2 and 1-3 connected atoms
    !
    !$omp parallel default (shared)                                       &
    !$omp firstprivate(num_excl)                                          &
    !$omp private(i, ini_excl, fin_excl, j, jatm, rij2, k, rij, qij,      &
    !$omp         uij, fij, vxx, vyx, vzx, vyy, vzy, vzz, Res, L)         &
    !$omp shared(natm, alpha, alpha2sp, alpha2m, coord, q, el_fact,       &
    !$omp        num_nonb_excl, nonb_excl_list, nproc_city, my_city_rank, &
    !$omp        i_list, force_wk)                                        &
    !$omp reduction(+:eelec) reduction(+:virial)
    !
    do i = 1, natm - 1

      ini_excl = num_excl + 1
      fin_excl = num_excl + num_nonb_excl(i)
      num_excl = fin_excl

      if (mod(i-1,nproc_city) == my_city_rank) then

        !$omp do
        do j = ini_excl, fin_excl

          jatm      = nonb_excl_list(j)
          i_list(j) = i

          ! Calculating R_ij
          rij2 = 0.0_wp

          do k = 1, 3
            rij(k) = coord(k,i) - coord(k,jatm)
            rij2 = rij2 + rij(k)**2
          end do

          qij  = q(i) * q(jatm)

          rij2 = cutoff2 * enefunc%table%density /rij2
          L    = int(rij2)
          Res  = rij2 - L

          uij = table_elec(L) + Res*(table_elec(L+1)-table_elec(L))
          eelec = eelec + uij*qij

          ! Calculating F_ij and adding it to F_i (symmetry considered)
          !
          uij = table_delec(L) + Res*(table_delec(L+1)-table_delec(L))

          do k = 1, 3
            fij(k) = uij * rij(k) * qij
            force_wk(k,j) = fij(k)
          end do

          ! virial
          !
          vxx = rij(1) * fij(1)
          vyx = rij(2) * fij(1)
          vzx = rij(3) * fij(1)
          vyy = rij(2) * fij(2)
          vzy = rij(3) * fij(2)
          vzz = rij(3) * fij(3)
          virial(1,1) = virial(1,1) - vxx
          virial(2,1) = virial(2,1) - vyx
          virial(3,1) = virial(3,1) - vzx
          virial(1,2) = virial(1,2) - vyx
          virial(2,2) = virial(2,2) - vyy
          virial(3,2) = virial(3,2) - vzy
          virial(1,3) = virial(1,3) - vzx
          virial(2,3) = virial(2,3) - vzy
          virial(3,3) = virial(3,3) - vzz

        end do
        !$omp end do

      end if

    end do
    !$omp end parallel

    do j = 1, size(enefunc%nonb_excl_list(:))

      iatm = i_list(j)
      jatm = nonb_excl_list(j)

      if (iatm > 0) then
        do k = 1, 3
          force(k,iatm,1) = force(k,iatm,1) - force_wk(k,j)
          force(k,jatm,1) = force(k,jatm,1) + force_wk(k,j)
        end do
      end if

    end do

    return

  end subroutine pme_bond_corr_linear_general

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber
  !> @brief        Calculate PME bonded energy correction with lookup table
  !! @authors      CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber(enefunc, molecule, coord, &
                                            force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: work(3)
    real(wp)                  :: R
    real(wp)                  :: term_elec, grad_coef
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: qq_scale, cc
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_elec(:), table_delec(:)
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
    table_elec     => enefunc%table%table_ecor
    table_delec    => enefunc%table%table_decor

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,natom))
               !force_j(3,size(enefunc%nb14_calc_list)))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    !force_j(1:3,1:size(enefunc%nb14_calc_list)) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef,                                       &
    !$omp         term_elec, force_local, qq_scale, cc)                 &
    !$omp reduction(+:virial) reduction(+:eelec) 
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
!          qq_scale = nb14_qq_scale(k,ii)-1.0_wp
          qq_scale = -nb14_qq_scale(k,ii)+1.0_wp

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

            term_elec   = table_elec(L) + R*(table_elec(L+1)-table_elec(L))
            cc    = charge(i)*charge(j)*qq_scale
            eelec = eelec + term_elec*cc

            term_elec   = table_delec(L) + R*(table_delec(L+1)-table_delec(L))
            grad_coef = cc*term_elec

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

  end subroutine pme_bond_corr_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_check
  !> @brief        Calculate PME bonded energy correction with lookup table
  !! @authors      JJ, TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_check(enefunc, molecule, coord, &
                                                force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    ! local variables
    real(wp)                  :: u, el_fact, alpha, alpha2m, alpha2sp
    real(wp)                  :: rij(3), rij2, qij, uij, fij(3)
    real(wp)                  :: vxx,vyx,vzx,vyy,vzy,vzz
    real(wp)                  :: Res, cutoff2
    real(wp)                  :: minimum_contact
    integer                   :: i, j, k, jatm
    integer                   :: num_excl,ini_excl,fin_excl
    integer                   :: iatm, L

    real(wp),         pointer :: q(:), force_wk(:,:)
    real(wp),         pointer :: table_elec(:), table_delec(:)
    integer,          pointer :: natm, num_nonb_excl(:), nonb_excl_list(:)
    integer,          pointer :: i_list(:)


    ! use pointers
    !
    natm           => molecule%num_atoms
    q              => molecule%charge
    num_nonb_excl  => enefunc%num_nonb_excl
    nonb_excl_list => enefunc%nonb_excl_list
    table_elec     => enefunc%table%table_ecor
    table_delec    => enefunc%table%table_decor
    force_wk       => enefunc%pme%force_wk
    i_list         => enefunc%pme%i_list
    el_fact        =  ELECOEF/enefunc%dielec_const
    alpha          =  enefunc%pme%alpha
    alpha2m        =  enefunc%pme%alpha2m
    alpha2sp       =  enefunc%pme%alpha2sp
    cutoff2        =  enefunc%cutoffdist * enefunc%cutoffdist
    minimum_contact=  enefunc%minimum_contact

    ! Initializing the energy and force 
    !
    u        = 0.0_wp
    num_excl = 0
    i_list(:)     = 0
    force_wk(:,:) = 0.0_wp

    ! Running only the 1-2 and 1-3 connected atoms
    !
    !$omp parallel default (shared)                                       &
    !$omp firstprivate(num_excl)                                          &
    !$omp private(i, ini_excl, fin_excl, j, jatm, rij2, k, rij, qij,      &
    !$omp         uij, fij, vxx, vyx, vzx, vyy, vzy, vzz, Res, L)         &
    !$omp shared(natm, alpha, alpha2sp, alpha2m, coord, q, el_fact,       &
    !$omp        num_nonb_excl, nonb_excl_list, nproc_city, my_city_rank, &
    !$omp        i_list, force_wk)                                        &
    !$omp reduction(+:eelec) reduction(+:virial)
    !
    do i = 1, natm - 1

      ini_excl = num_excl + 1
      fin_excl = num_excl + num_nonb_excl(i)
      num_excl = fin_excl

      if (mod(i-1,nproc_city) == my_city_rank) then

        !$omp do
        do j = ini_excl, fin_excl

          jatm      = nonb_excl_list(j)
          i_list(j) = i

          ! Calculating R_ij
          rij2 = 0.0_wp

          do k = 1, 3
            rij(k) = coord(k,i) - coord(k,jatm)
            rij2 = rij2 + rij(k)**2
          end do
          rij2     = max(rij2, minimum_contact)

          qij  = q(i) * q(jatm)

          rij2 = cutoff2 * enefunc%table%density /rij2
          L    = int(rij2)
          Res  = rij2 - L

          uij = table_elec(L) + Res*(table_elec(L+1)-table_elec(L))
          eelec = eelec + uij*qij

          ! Calculating F_ij and adding it to F_i (symmetry considered)
          !
          uij = table_delec(L) + Res*(table_delec(L+1)-table_delec(L))

          do k = 1, 3
            fij(k) = uij * rij(k) * qij
            force_wk(k,j) = fij(k)
          end do

          ! virial
          !
          vxx = rij(1) * fij(1)
          vyx = rij(2) * fij(1)
          vzx = rij(3) * fij(1)
          vyy = rij(2) * fij(2)
          vzy = rij(3) * fij(2)
          vzz = rij(3) * fij(3)
          virial(1,1) = virial(1,1) - vxx
          virial(2,1) = virial(2,1) - vyx
          virial(3,1) = virial(3,1) - vzx
          virial(1,2) = virial(1,2) - vyx
          virial(2,2) = virial(2,2) - vyy
          virial(3,2) = virial(3,2) - vzy
          virial(1,3) = virial(1,3) - vzx
          virial(2,3) = virial(2,3) - vzy
          virial(3,3) = virial(3,3) - vzz

        end do
        !$omp end do

      end if

    end do
    !$omp end parallel

    do j = 1, size(enefunc%nonb_excl_list(:))

      iatm = i_list(j)
      jatm = nonb_excl_list(j)

      if (iatm > 0) then
        do k = 1, 3
          force(k,iatm,1) = force(k,iatm,1) - force_wk(k,j)
          force(k,jatm,1) = force(k,jatm,1) + force_wk(k,j)
        end do
      end if

    end do

    return

  end subroutine pme_bond_corr_linear_general_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber_check
  !> @brief        Calculate PME bonded energy correction with lookup table
  !! @authors      CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates 
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber_check(enefunc, molecule, coord, &
                                                  force, virial, eelec)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec

    real(wp)                  :: dij(3)
    real(wp)                  :: rij2
    real(wp)                  :: cutoff2
    real(wp)                  :: work(3)
    real(wp)                  :: R
    real(wp)                  :: term_elec, grad_coef
    real(wp)                  :: force_local(1:3)
    real(wp)                  :: qq_scale, cc
    real(wp)                  :: minimum_contact
    integer                   :: ii, i, j, k, m, L
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: id, my_id
    integer                   :: num_nb14_all, max_nb14

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_elec(:), table_delec(:)
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
    table_elec     => enefunc%table%table_ecor
    table_delec    => enefunc%table%table_decor
    minimum_contact=  enefunc%minimum_contact

    cutoff2  = cutoff * cutoff
    num_nb14 = 0
    max_nb14 = max(1,maxval(num_nb14_calc(1:natom)))

    if (.not. allocated(force_i)) then
      allocate(force_i(3,enefunc%table%num_solute), &
               force_j(3,max_nb14,natom))
               !force_j(3,size(enefunc%nb14_calc_list)))
    end if

    force_i(1:3,1:natom) = 0.0_wp
    !force_j(1:3,1:size(enefunc%nb14_calc_list)) = 0.0_wp
    force_j(1:3,1:max_nb14,1:natom) = 0.0_wp

    ! calculate energy and gradient
    !
    !$omp parallel                                                             &
    !$omp firstprivate(num_nb14)                                               &
    !$omp private(id, my_id, ini_nb14, fin_nb14, ii, i, k, j, m, dij, rij2,    &
    !$omp         R, L, work, grad_coef,                                       &
    !$omp         term_elec, force_local, qq_scale, cc)                 &
    !$omp reduction(+:virial) reduction(+:eelec) 
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
!          qq_scale = nb14_qq_scale(k,ii)-1.0_wp
          qq_scale = -nb14_qq_scale(k,ii)+1.0_wp

          ! compute distance
          !
          dij(1:3) = coord(1:3,i) - coord(1:3,j)
          rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2 = max(rij2, minimum_contact)

          ! cutoff
          !
          if (rij2 < cutoff2) then

            rij2 = cutoff2 * enefunc%table%density / rij2
            L = int(rij2)
            R = rij2 - L

            term_elec   = table_elec(L) + R*(table_elec(L+1)-table_elec(L))
            cc    = charge(i)*charge(j)*qq_scale
            eelec = eelec + term_elec*cc

            term_elec   = table_delec(L) + R*(table_delec(L+1)-table_delec(L))
            grad_coef = cc*term_elec

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

  end subroutine pme_bond_corr_linear_gro_amber_check

end module at_energy_table_linear_bondcorr_mod
