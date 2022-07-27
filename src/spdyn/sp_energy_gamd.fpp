!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_gamd_mod
!> @brief   compute energy for GaMD
!           Refs: Miao Y. et al., J. Chem. Theory Comput. 11, 3584 (2015)
!                 Pang Y. et al., J. Chem. Theory Comput. 13, 9 (2017)
!! @authors Hiraku Oshima (HO)
!  
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_gamd_mod

  use sp_energy_restraints_mod
  use sp_energy_dihedrals_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_restraints_str_mod
  use sp_domain_str_mod
  use sp_boundary_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
  use string_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public :: compute_boost_potential_gamd
  public :: compute_stat_gamd
  public :: boost_gamd
  public :: compute_energy_restraints_gamd
  public :: compute_energy_dihed_gamd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_boost_potential_gamd
  !> @brief        compute boost potential in GaMD
  !! @authors      HO
  !! @param[in]    ene     : input energy
  !! @param[in]    k       : boost parameter for GaMD
  !! @param[in]    ene_th  : energy threshold for GaMD
  !! @param[out]   e_gamd  : boost potential
  !! @param[put]   w_gamd  : boost coefficient for force
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_boost_potential_gamd(ene, k, ene_th, e_gamd, w_gamd)

    ! formal arguments
    real(dp),                intent(in)  :: ene
    real(dp),                intent(in)  :: k
    real(dp),                intent(in)  :: ene_th
    real(dp),                intent(out) :: e_gamd
    real(wp),                intent(out) :: w_gamd

    ! local variables
    real(dp) :: ev, factor


    e_gamd = 0.0_dp
    if (ene < ene_th) then
      ev = ene - ene_th
      factor = k * ev
      e_gamd = factor * 0.5_dp * ev
      w_gamd = 1.0_wp + real(factor, wp)
    else
      w_gamd = 1.0_wp
    end if

    return

  end subroutine compute_boost_potential_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_stat_gamd
  !> @brief        compute statistics for GaMD
  !! @authors      HO
  !! @param[in]    E       : input energy
  !! @param[inout] E_max   : maximum energy
  !! @param[inout] E_min   : minimum energy
  !! @param[inout] E_ave   : average of energy
  !! @param[inout] E_ave2  : average of square energy
  !! @param[inout] E_dev   : deviation of energy
  !! @param[inout] counter : number of samples
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_stat_gamd(E, E_max, E_min, E_ave, E_ave2, E_dev, counter)

    ! formal arguments
    real(dp),                intent(in)    :: E
    real(dp),                intent(inout) :: E_max
    real(dp),                intent(inout) :: E_min
    real(dp),                intent(inout) :: E_ave
    real(dp),                intent(inout) :: E_ave2
    real(dp),                intent(inout) :: E_dev
    integer,                 intent(inout) :: counter

    ! local variables
    real(dp) :: E_diff


    if (E > E_max) E_max = E
    if (E < E_min) E_min = E

    E_diff = E - E_ave
    E_ave  = E_ave  + E_diff/real(counter, dp)
    E_ave2 = E_ave2 + E_diff * (E - E_ave)
    E_dev = dsqrt(E_ave2/real(counter, dp))
    counter = counter + 1

    return

  end subroutine compute_stat_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    boost_gamd
  !> @brief        boost the potential energy using GaMD method
  !! @authors      HO
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[inout] energy     : energy information
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine boost_gamd(domain, enefunc, energy, force, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),target,  intent(inout) :: enefunc
    type(s_energy),          intent(inout) :: energy
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)

    ! local variable
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    real(wp)                 :: f_tmp(1:3)
    real(dp)                 :: v_tmp(3,3)
    real(wp)                 :: w_pot_gamd, w_dih_gamd
    real(dp)                 :: ene_pot, ene_dih
    real(dp)                 :: before_reduce(2), after_reduce(2)

    type(s_enefunc_gamd), pointer :: gamd


    gamd => enefunc%gamd

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! Reduce potential energy or dihedral energy.
    !
#ifdef HAVE_MPI_GENESIS
    if (gamd%boost_dih) then
      before_reduce(1) = energy%dihedral
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        before_reduce(1) = before_reduce(1) + energy%cmap
      end if
      call mpi_allreduce(before_reduce, after_reduce, 1, mpi_real8,  &
        mpi_sum, mpi_comm_country, ierror)
      ene_dih = after_reduce(1)
    else if (gamd%boost_dual) then
      before_reduce(1) = energy%bond          &
                       + energy%angle         &
                       + energy%urey_bradley  &
                       + energy%improper      &
                       + energy%electrostatic &
                       + energy%van_der_waals
      before_reduce(2) = energy%dihedral
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        before_reduce(2) = before_reduce(2) + energy%cmap
      end if
      call mpi_allreduce(before_reduce, after_reduce, 2, mpi_real8,  &
        mpi_sum, mpi_comm_country, ierror)
      ene_pot = after_reduce(1)
      ene_dih = after_reduce(2)
    else if (gamd%boost_pot) then
      before_reduce(1) = energy%bond          &
                       + energy%angle         &
                       + energy%urey_bradley  &
                       + energy%dihedral      &
                       + energy%improper      &
                       + energy%electrostatic &
                       + energy%van_der_waals
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        before_reduce(1) = before_reduce(1) + energy%cmap
      end if
      call mpi_allreduce(before_reduce, after_reduce, 1, mpi_real8,  &
        mpi_sum, mpi_comm_country, ierror)
      ene_pot = after_reduce(1)
    end if
#else
    if (gamd%boost_dih) then
      ene_dih = energy%dihedral
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        ene_dih = ene_dih + energy%cmap
      end if
    else if (gamd%boost_dual) then
      ene_pot = energy%bond          &
              + energy%angle         &
              + energy%urey_bradley  &
              + energy%improper      &
              + energy%electrostatic &
              + energy%van_der_waals
      ene_dih = energy%dihedral
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        ene_dih = ene_dih + energy%cmap
      end if
    else if (gamd%boost_pot) then
      ene_pot = energy%bond          &
              + energy%angle         &
              + energy%urey_bradley  &
              + energy%dihedral      &
              + energy%improper      &
              + energy%electrostatic &
              + energy%van_der_waals
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        ene_pot = ene_pot + energy%cmap
      end if
    end if
#endif
    ! Compute statistics
    !
    if (gamd%gamd_stat) then

      if (gamd%boost_dih) then

        ! Compute statistics of dihedral energy (i.e., max, min, ave, sd)
        ! for updating the gamd parameters.
        !
        call compute_stat_gamd(ene_dih, gamd%ene_dih_max, gamd%ene_dih_min, &
          gamd%ene_dih_ave, gamd%ene_dih_ave2, gamd%ene_dih_dev,       &
          gamd%count_dih)

      else if (gamd%boost_pot) then

        ! Compute statistics of total potential energy (i.e., max, min, ave, sd)
        ! for updating the gamd parameters.
        !
        call compute_stat_gamd(ene_pot, gamd%ene_pot_max, gamd%ene_pot_min, &
          gamd%ene_pot_ave, gamd%ene_pot_ave2, gamd%ene_pot_dev,       &
          gamd%count_pot)

      else if (gamd%boost_dual) then

        ! Compute statistics of dihedral energy (i.e., max, min, ave, sd)
        ! for updating the gamd parameters.
        !
        call compute_stat_gamd(ene_dih, gamd%ene_dih_max, gamd%ene_dih_min, &
          gamd%ene_dih_ave, gamd%ene_dih_ave2, gamd%ene_dih_dev,       &
          gamd%count_dih)

        ! Compute statistics of total potential energy (i.e., max, min, ave, sd)
        ! for updating the gamd parameters.
        !
        call compute_stat_gamd(ene_pot, gamd%ene_pot_max, gamd%ene_pot_min, &
          gamd%ene_pot_ave, gamd%ene_pot_ave2, gamd%ene_pot_dev,       &
          gamd%count_pot)

      end if

    end if

    ! Skip boosting or not
    !
    if (.not. gamd%gamd_boost) then
      return
    end if

    ! Boost potential and force
    !
    if (gamd%boost_dih) then

      ! Compute boost potential (= dihedral_gamd) and 
      ! weight of force (= w_dih_gamd) for dihedral.
      !
      call compute_boost_potential_gamd(ene_dih, gamd%k_dih, gamd%ene_dih_th, &
        energy%dihedral_gamd, w_dih_gamd)

      ! Reweight the force of dihedral (+cmap if charmm) (= f_dihe_omp).
      !
      !$omp parallel do default(shared) private(id, i, ix, f_tmp) &
      !$omp schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_tmp(1:3) = (w_dih_gamd-1.0_wp) * gamd%f_dihe_omp(1:3,ix,i,1)
          do id = 2, nthread
            f_tmp(1:3) = f_tmp(1:3) + (w_dih_gamd-1.0_wp) &
              * gamd%f_dihe_omp(1:3,ix,i,id)
          end do
          force(1:3,ix,i) = force(1:3,ix,i) + f_tmp(1:3)
        end do
      end do
      !$omp end parallel do

      ! Reweight virial.
      do id = 1, nthread
        virial(1:3,1:3) = virial(1:3,1:3) + (w_dih_gamd-1.0_wp) &
          * gamd%v_dihe_omp(1:3,1:3,id)
      end do

    else if (gamd%boost_pot) then

      ! Compute boost potential (= total_gamd) and 
      ! weight of force (= w_pot_gamd) for total potential.
      !
      call compute_boost_potential_gamd(ene_pot, gamd%k_pot, gamd%ene_pot_th, &
        energy%total_gamd, w_pot_gamd)

      ! Reweight force and virial of total potential.
      !
      !$omp parallel do default(shared) private(id, i, ix, f_tmp) &
      !$omp schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_tmp(1:3) = (1.0_wp - w_pot_gamd) * gamd%f_rest_omp(1:3,ix,i,1)
          do id = 2, nthread
            f_tmp(1:3) = f_tmp(1:3) &
              + (1.0_wp - w_pot_gamd) * gamd%f_rest_omp(1:3,ix,i,id)
          end do
          force(1:3,ix,i) = w_pot_gamd * force(1:3,ix,i) + f_tmp(1:3)
        end do
      end do
      !$omp end parallel do
      v_tmp(1:3,1:3) = (1.0_wp - w_pot_gamd) * gamd%v_rest_omp(1:3,1:3,1)
      do id = 2, nthread
        v_tmp(1:3,1:3) = v_tmp(1:3,1:3) &
          + (1.0_wp - w_pot_gamd)*gamd%v_rest_omp(1:3,1:3,id)
      end do
      virial(1:3,1:3) = w_pot_gamd*virial(1:3,1:3) + v_tmp(1:3,1:3)

    else if (gamd%boost_dual) then

      ! Compute boost potential (= dihedral_gamd) and 
      ! weight of force (= w_dih_gamd) for dihedral.
      !
      call compute_boost_potential_gamd(ene_dih, gamd%k_dih, gamd%ene_dih_th, &
        energy%dihedral_gamd, w_dih_gamd)

      ! Compute boost potential (= total_gamd) and 
      ! weight of force (= w_pot_gamd) for total potential.
      !
      call compute_boost_potential_gamd(ene_pot, gamd%k_pot, gamd%ene_pot_th, &
        energy%total_gamd, w_pot_gamd)

      ! Reweight force and virial of total potential.
      !
      !$omp parallel do default(shared) private(id, i, ix, f_tmp) &
      !$omp schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_tmp(1:3) = (w_dih_gamd - w_pot_gamd) * gamd%f_dihe_omp(1:3,ix,i,1)&
            + (1.0_wp - w_pot_gamd) * gamd%f_rest_omp(1:3,ix,i,1)
          do id = 2, nthread
            f_tmp(1:3) = f_tmp(1:3) &
              + (w_dih_gamd - w_pot_gamd) * gamd%f_dihe_omp(1:3,ix,i,id) &
              + (1.0_wp - w_pot_gamd) * gamd%f_rest_omp(1:3,ix,i,id)
          end do
          force(1:3,ix,i) = w_pot_gamd * force(1:3,ix,i) + f_tmp(1:3)
        end do
      end do
      !$omp end parallel do

      v_tmp(1:3,1:3) = (w_dih_gamd - w_pot_gamd) * gamd%v_dihe_omp(1:3,1:3,1)&
        + (1.0_wp - w_pot_gamd) * gamd%v_rest_omp(1:3,1:3,1)
      do id = 2, nthread
        v_tmp(1:3,1:3) = v_tmp(1:3,1:3) &
          + (w_dih_gamd - w_pot_gamd) * gamd%v_dihe_omp(1:3,1:3,id) &
          + (1.0_wp - w_pot_gamd)*gamd%v_rest_omp(1:3,1:3,id)
      end do
      virial(1:3,1:3) = w_pot_gamd*virial(1:3,1:3) + v_tmp(1:3,1:3)

    end if

  end subroutine boost_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_gamd
  !> @brief        calculate restraint energy for GaMD
  !! @authors      HO
  !! @param[in]    get_coord  : flag for use of coord
  !! @param[in]    calc_force : flag for calculation of force
  !! @param[in]    domain     : domain information
  !! @param[in]    boundary   : boundary information
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] virial_ext : extern virial term of target systems
  !! @param[inout] eposi      : point restraint energy of target systems
  !! @param[inout] energy     : energy information
  !! @detail       force and virial are saved to scale them after energy
  !!               evaluation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_gamd(get_coord, calc_force, domain, &
                      boundary, enefunc, coord, force, virial, virial_ext, &
                      eposi, energy)

    ! formal arguments
    logical,                 intent(in)    :: get_coord
    logical,                 intent(in)    :: calc_force
    type(s_domain),  target, intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: virial_ext(3,3,nthread)
    real(dp),                intent(inout) :: eposi(nthread)
    type(s_energy),          intent(inout) :: energy

    ! local variables
    integer           :: ncell, natom, id, i, j, k, ix, ic, jc
    integer           :: omp_get_thread_num
    real(wp), pointer :: f_rest_omp(:,:,:,:)
    real(dp), pointer :: v_rest_omp(:,:,:)


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

      f_rest_omp => enefunc%gamd%f_rest_omp
      v_rest_omp => enefunc%gamd%v_rest_omp

      ! initialize
      !$omp parallel do
      do id = 1, nthread
        f_rest_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      v_rest_omp(1:3,1:3,1:nthread) = 0.0_dp

      !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_rest_omp(1:3,ix,i,id+1) = force(1:3,ix,i,id+1)
        end do
      end do
      !$omp end parallel
      do id = 1, nthread
        v_rest_omp(1:3,1:3,id) = virial(1:3,1:3,id)
      end do

    end if

    call compute_energy_restraints(get_coord, calc_force, domain, boundary, &
                        enefunc, coord, &
                        force, virial, virial_ext, &
                        eposi, energy%restraint_rmsd, &
                        energy%rmsd, energy%restraint_distance, &
                        energy%restraint_emfit, energy%emcorr)

    if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

      !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_rest_omp(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) &
            - f_rest_omp(1:3,ix,i,id+1)
        end do
      end do
      !$omp end parallel
      do id = 1, nthread
        v_rest_omp(1:3,1:3,id) = virial(1:3,1:3,id) - &
          v_rest_omp(1:3,1:3,id)
      end do

    end if

    return

  end subroutine compute_energy_restraints_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed_gamd
  !> @brief        calculate dihedral, cmap, and improper energies for GaMD
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    npt      : flag for NPT or not
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] edihe    : dihedral energy of target systems
  !! @param[inout] ecmap    : cmap energy of target systems
  !! @param[inout] eimprop  : improper energy of target systems
  !! @detail       force and virial are saved to scale them after energy
  !!               evaluation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed_gamd(domain, enefunc, npt, nonb_ene, &
                                  coord, force, virial, edihed, ecmap, eimprop)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: edihed(nthread)
    real(dp),                intent(inout) :: ecmap(nthread)
    real(dp),                intent(inout) :: eimprop(nthread)

    ! local variables
    integer           :: ncell, natom, id, i, j, k, ix, ic, jc
    integer           :: omp_get_thread_num
    real(wp), pointer :: f_dihe_omp(:,:,:,:)
    real(dp), pointer :: v_dihe_omp(:,:,:)


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

      f_dihe_omp => enefunc%gamd%f_dihe_omp
      v_dihe_omp => enefunc%gamd%v_dihe_omp

      ! initialize
      !$omp parallel do
      do id = 1, nthread
        f_dihe_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
        v_dihe_omp(1,1,id) = virial(1,1,id)
        v_dihe_omp(2,2,id) = virial(2,2,id)
        v_dihe_omp(3,3,id) = virial(3,3,id)
      end do
      !$omp end parallel do

      !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_dihe_omp(1:3,ix,i,id+1) = force(1:3,ix,i,id+1)
        end do
      end do
      !$omp end parallel

    end if

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres(domain, enefunc, coord, &
        force, virial, edihed)
    else
      call compute_energy_dihed(domain, enefunc, coord, &
        force, virial, edihed)
    end if

    if (enefunc%forcefield == ForcefieldCHARMM .or. &
        enefunc%forcefield == ForcefieldAMBER) then
      ! cmap energy
      ! 
      call compute_energy_cmap(domain, enefunc, coord, force, virial, ecmap)
    end if

    if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

      !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          f_dihe_omp(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) &
            - f_dihe_omp(1:3,ix,i,id+1)
        end do
      end do
      v_dihe_omp(1,1,id+1) = virial(1,1,id+1) - v_dihe_omp(1,1,id+1)
      v_dihe_omp(2,2,id+1) = virial(2,2,id+1) - v_dihe_omp(2,2,id+1)
      v_dihe_omp(3,3,id+1) = virial(3,3,id+1) - v_dihe_omp(3,3,id+1)
      !$omp end parallel

    end if

    return

  end subroutine compute_energy_dihed_gamd

end module sp_energy_gamd_mod
