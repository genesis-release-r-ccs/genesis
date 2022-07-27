!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_gamd_mod
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

module at_energy_gamd_mod

  use at_energy_restraints_mod
  use at_energy_dihedrals_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_boundary_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
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
    real(wp),                intent(in)  :: ene
    real(wp),                intent(in)  :: k
    real(wp),                intent(in)  :: ene_th
    real(wp),                intent(out) :: e_gamd
    real(wp),                intent(out) :: w_gamd

    ! local variables
    real(wp) :: ev, factor


    e_gamd = 0.0_wp
    if (ene < ene_th) then
      ev = ene - ene_th
      factor = k * ev
      e_gamd = factor * 0.5_wp * ev
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
    real(wp),                intent(in)     :: E
    real(wp),                intent(inout)  :: E_max
    real(wp),                intent(inout)  :: E_min
    real(wp),                intent(inout)  :: E_ave
    real(wp),                intent(inout)  :: E_ave2
    real(wp),                intent(inout)  :: E_dev
    integer,                 intent(inout)  :: counter

    ! local variables
    real(wp) :: E_diff


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

  subroutine boost_gamd(enefunc, natom, energy, temporary, force, virial)

    ! formal arguments
    type(s_enefunc),target,  intent(inout) :: enefunc
    integer,                 intent(in)    :: natom
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)

    ! local variable
    real(wp)                     :: after_allreduce(999), before_allreduce(999)
    integer                      :: i, j, k, n
    real(wp)                     :: w_pot_gamd, w_dih_gamd
    real(wp)                     :: ene_pot, ene_dih
    integer                      :: ncycle, icycle, nlen, ixx

    type(s_enefunc_gamd),pointer :: gamd


    gamd => enefunc%gamd

    ! Calculate potential energy or dihedral energy.
    !
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
              + energy%van_der_waals &
              + energy%contact       &
              + energy%noncontact    &
              + energy%solvation
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
              + energy%van_der_waals &
              + energy%contact       &
              + energy%noncontact    &
              + energy%solvation
      if (enefunc%forcefield == ForcefieldCHARMM .or. &
          enefunc%forcefield == ForcefieldAMBER) then
        ene_pot = ene_pot + energy%cmap
      end if
    end if

#ifdef HAVE_MPI_GENESIS
    if (.not. gamd%boost_pot) then
      ! Allreduce force
      !
      ncycle = (natom - 1) / mpi_drain + 1
      nlen   = mpi_drain
      ixx    = 1

      temporary(1:3,1:natom) = 0.0_wp
      do i = 1, natom
        temporary(1,i) = gamd%f_dihe(1,i)
        temporary(2,i) = gamd%f_dihe(2,i)
        temporary(3,i) = gamd%f_dihe(3,i)
      end do
      do icycle = 1, ncycle
        if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
        call mpi_allreduce(temporary(1,ixx), gamd%f_dihe(1,ixx), 3*nlen, &
                        mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
        ixx = ixx + nlen
      end do

      ! Allreduce virial
      !
      n = 0
      do i = 1, 3
        do j = 1, 3
          n = n + 1
          before_allreduce(n) = gamd%v_dihe(i,j)
        end do
      end do

      call mpi_allreduce(before_allreduce, after_allreduce, 9, &
                         mpi_wp_real,  mpi_sum,                 &
                         mpi_comm_country, ierror)

      n = 0
      do i = 1, 3
        do j = 1, 3
          n = n + 1
          gamd%v_dihe(i,j) = after_allreduce(n)
        end do
      end do
    end if


    if (.not. gamd%boost_dih) then
      ! Allreduce force
      !
      ncycle = (natom - 1) / mpi_drain + 1
      nlen   = mpi_drain
      ixx    = 1

      temporary(1:3,1:natom) = 0.0_wp
      do i = 1, natom
        temporary(1,i) = gamd%f_rest(1,i)
        temporary(2,i) = gamd%f_rest(2,i)
        temporary(3,i) = gamd%f_rest(3,i)
      end do
      do icycle = 1, ncycle
        if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
        call mpi_allreduce(temporary(1,ixx), gamd%f_rest(1,ixx), 3*nlen, &
                        mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
        ixx = ixx + nlen
      end do

      ! Allreduce virial
      !
      n = 0
      do i = 1, 3
        do j = 1, 3
          n = n + 1
          before_allreduce(n) = gamd%v_rest(i,j)
        end do
      end do

      call mpi_allreduce(before_allreduce, after_allreduce, 9, &
                         mpi_wp_real,  mpi_sum,                 &
                         mpi_comm_country, ierror)

      n = 0
      do i = 1, 3
        do j = 1, 3
          n = n + 1
          gamd%v_rest(i,j) = after_allreduce(n)
        end do
      end do
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

      ! Reweight the force of dihedral (+cmap if charmm) (= f_dihe)
      ! and virial.
      !
      do i = 1, natom
        force(1:3,i) = force(1:3,i) + (w_dih_gamd-1.0_wp) * gamd%f_dihe(1:3,i)
      end do
      virial(1:3,1:3) = virial(1:3,1:3) + (w_dih_gamd-1.0_wp)*&
        gamd%v_dihe(1:3,1:3)

    else if (gamd%boost_pot) then

      ! Compute boost potential (= total_gamd) and 
      ! weight of force (= w_pot_gamd) for total potential.
      !
      call compute_boost_potential_gamd(ene_pot, gamd%k_pot, &
        gamd%ene_pot_th, energy%total_gamd, w_pot_gamd)

      ! Reweight force and virial of total potential
      !
      do i = 1, natom
        force(1:3,i) = w_pot_gamd * force(1:3,i) &
                       + (1.0_wp - w_pot_gamd) * gamd%f_rest(1:3,i)
      end do
      virial(1:3,1:3) = w_pot_gamd*virial(1:3,1:3) &
                        + (1.0_wp - w_pot_gamd) * gamd%v_rest(1:3,1:3)

    else if (gamd%boost_dual) then

      ! Compute boost potential (= dihedral_gamd) and 
      ! weight of force (= w_dih_gamd) for dihedral.
      !
      call compute_boost_potential_gamd(ene_dih, gamd%k_dih, &
        gamd%ene_dih_th, energy%dihedral_gamd, w_dih_gamd)

      ! Compute boost potential (= total_gamd) and 
      ! weight of force (= w_pot_gamd) for total potential.
      !
      call compute_boost_potential_gamd(ene_pot, gamd%k_pot, &
        gamd%ene_pot_th, energy%total_gamd, w_pot_gamd)

      ! Reweight force and virial of total potential.
      !
      do i = 1, natom
        force(1:3,i) = w_pot_gamd * force(1:3,i) &
                       + (w_dih_gamd - w_pot_gamd) * gamd%f_dihe(1:3,i) &
                       + (1.0_wp - w_pot_gamd) * gamd%f_rest(1:3,i)
      end do
      virial(1:3,1:3) = w_pot_gamd*virial(1:3,1:3) &
                        + (w_dih_gamd - w_pot_gamd) * gamd%v_dihe(1:3,1:3) &
                        + (1.0_wp - w_pot_gamd) * gamd%v_rest(1:3,1:3)

    end if

    return

  end subroutine boost_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_gamd
  !> @brief        calculate restraint energy for GaMD
  !! @authors      HO
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    boundary   : boundary information
  !! @param[in]    natom      : number of atoms
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] virial_ext : extern virial term of target systems
  !! @param[inout] energy     : energy information
  !! @detail       force and virial are saved to scale them after energy
  !!               evaluation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_gamd(enefunc, boundary, natom, coord, &
                                            force, virial, virial_ext, energy)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)
    type(s_energy),          intent(inout) :: energy

    ! local variables
    integer                  :: i, j, k, ix, ic, jc
    real(wp), pointer :: f_rest(:,:)
    real(wp), pointer :: v_rest(:,:)


    if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

      f_rest => enefunc%gamd%f_rest
      v_rest => enefunc%gamd%v_rest

      do i = 1, natom
        f_rest(1:3,i) = force(1:3,i)
      end do

      do j = 1, 3
        do k = 1, 3
          v_rest(j,k) = virial(j,k)
        end do
      end do

    end if

    call compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                        virial_ext, energy)

    if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

      do i = 1, natom
        f_rest(1:3,i) = force(1:3,i) - f_rest(1:3,i)
      end do

      do j = 1, 3
        do k = 1, 3
          v_rest(j,k) = virial(j,k) - v_rest(j,k)
        end do
      end do

    end if

    return

  end subroutine compute_energy_restraints_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed_gamd
  !> @brief        calculate dihedral, cmap, and improper energies for GaMD
  !! @authors      HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    natom    : number of atoms
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] energy   : energy information
  !! @detail       force and virial are saved to scale them after energy
  !!               evaluation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed_gamd(enefunc, natom, coord, force, &
                                       virial, energy)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    type(s_energy),          intent(inout) :: energy

    ! local variables
    integer                  :: id, i, j, k, ix, ic, jc
    real(wp), pointer :: f_dihe(:,:)
    real(wp), pointer :: v_dihe(:,:)


    if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

      f_dihe => enefunc%gamd%f_dihe
      v_dihe => enefunc%gamd%v_dihe

      do i = 1, natom
        f_dihe(1:3,i) = force(1:3,i)
      end do

      do j = 1, 3
        do k = 1, 3
          v_dihe(j,k) = virial(j,k)
        end do
      end do

    end if

    ! dihedral energy
    !
    call compute_energy_dihed(enefunc, coord, force, virial,  &
                          energy%dihedral)

    if (enefunc%forcefield == ForcefieldCHARMM .or. &
        enefunc%forcefield == ForcefieldAMBER) then
      ! cmap energy
      !
      call compute_energy_cmap(enefunc, coord, force, virial,   &
                            energy%cmap)
    end if

    if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

      do i = 1, natom
        f_dihe(1:3,i) = force(1:3,i) - f_dihe(1:3,i)
      end do

      do j = 1, 3
        do k = 1, 3
          v_dihe(j,k) = virial(j,k) - v_dihe(j,k)
        end do
      end do

    end if

    return

  end subroutine compute_energy_dihed_gamd

end module at_energy_gamd_mod
