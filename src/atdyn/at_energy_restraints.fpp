!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_restraints_mod
!> @brief   calculate restraints energy
!! @authors Chigusa Kobayashi (CK), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_restraints_mod

  use at_boundary_str_mod
  use at_rpath_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_experiments_mod
  use fitting_mod
  use fitting_str_mod
  use dihedral_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_energy_restraint
  public  :: compute_energy_restraint_each
  private :: compute_energy_positional_restraint
  private :: compute_energy_distance_restraint
  private :: compute_energy_angle_restraint
  private :: compute_energy_dihedral_restraint
  private :: compute_energy_rmsd_restraint
  private :: compute_energy_pc_restraint
  private :: compute_energy_repulsive_restraint
  private :: compute_energy_flatbottom_restraint
  private :: compute_energy_rg_restraint
 
contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraint
  !> @brief        calculate restraint energy
  !! @authors      TM, CK
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    boundary   : information of boundary
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : virial of target systems
  !! @param[inout] energy     : information of energy
  !! @param[inout] cv         : cv of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                                      virial_ext, energy)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)
    type(s_energy),          intent(inout) :: energy

    ! local variables
    real(wp)                 :: eneval
    real(wp)                 :: cv
    integer                  :: inum
    logical                  :: calc_force


    do inum = 1, enefunc%num_restraintfuncs

      ! only main_rank computes both energy and force
      ! and other ranks compute only energy
      ! After energy calculation, force is allreduced, while energy is not.
      ! However, only EM restraint energy is allreduced
      !
      if (main_rank .or. replica_main_rank) then
        calc_force = .true.
      else
        if (enefunc%restraint_kind(inum) == RestraintsFuncEM) then
          calc_force = .true.
        else
          calc_force = .false.
        end if
      end if


      call compute_energy_restraint_each(inum, calc_force, enefunc, boundary, &
                                  coord, force, virial, virial_ext, eneval, cv)
 
      energy%restraint(inum) = eneval
      energy%restraint_cv(inum) = cv

    end do

    return

  end subroutine compute_energy_restraint

!!!develop

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraint_each
  !> @brief        calculate restraint energy
  !! @authors      CK, TM
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    boundary   : information of boundary
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : virial of target systems
  !! @param[inout] eneval     : restraint energy of target systems
  !! @param[inout] cv         : cv of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraint_each(inum, calc_force, enefunc, boundary,&
                                   coord, force, virial, virial_ext, eneval, cv)
 
    ! formal arguments
    integer,                  intent(in)    :: inum
    logical,                  intent(in)    :: calc_force
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_boundary),         intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: virial_ext(3,3)
    real(wp),                 intent(inout) :: eneval
    real(wp),                 intent(inout) :: cv


    call timer(TimerRestraint, TimerOn)

    eneval = 0.0_wp
    cv     = 0.0_wp
 
    select case(enefunc%restraint_kind(inum)) 

    case (RestraintsFuncPOSI) 

      ! Positional restraint
      !
      call compute_energy_positional_restraint  &
        (enefunc, coord, inum, calc_force, force, virial, virial_ext, eneval)

    case (RestraintsFuncDIST:RestraintsFuncDISTCOM)

      ! Distance restraint
      !
      call compute_energy_distance_restraint    &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncANGLE:RestraintsFuncANGLECOM)

      ! Angle restraint
      !
      call compute_energy_angle_restraint       &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncDIHED:RestraintsFuncDIHEDCOM)

      ! Dihedral angle restraint
      !
      call compute_energy_dihedral_restraint    &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

      call compute_energy_rmsd_restraint        &
        (enefunc, coord, inum, calc_force, .false., force,  &
         virial, virial_ext, eneval, cv)

    case (RestraintsFuncPC:RestraintsFuncPCCOM)

      ! PC restraint
      !
      call compute_energy_pc_restraint        &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncEM)

      ! Experimental data restraint (only this goes to at_experiments.fpp)
      !
      call compute_energy_experimental_restraint &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncRG:RestraintsFuncRGWOMASS)

      call compute_energy_rg_restraint &
        (enefunc, coord, inum, calc_force, force, virial, eneval, cv)

    case (RestraintsFuncREPUL:RestraintsFuncREPULCOM)

      ! Repulsive restraint
      !
      call compute_energy_repulsive_restraint    &
        (enefunc, boundary, coord, inum, calc_force, force, virial, eneval)

    case (RestraintsFuncFB:RestraintsFuncFBCOM)

      ! flat-bottom restraint
      !
      call compute_energy_flatbottom_restraint    &
        (enefunc, boundary, coord, inum, calc_force, force, virial, eneval)

    end select

    call timer(TimerRestraint, TimerOff)

    return

  end subroutine compute_energy_restraint_each

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_positional_restraint
  !> @brief        calculate position restraint energy
  !! @authors      CK
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    inum        : pointer for restraint function
  !! @param[in]    const       : force constants
  !! @param[inout] force       : forces of target systems
  !! @param[inout] virial      : virial of target systems
  !! @param[out]   virial_ext  : external virial term of target systems
  !! @param[out]   eposi       : position constraint energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_positional_restraint(enefunc, coord, inum, &
                                   calc_force, force, virial, virial_ext, eposi)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)
    real(wp),                intent(out)   :: eposi

    ! local variables
    real(wp)                 :: d(1:3), d2(1:3)
    real(wp)                 :: rot_matrix(1:3,1:3)
    real(wp)                 :: const(0:3), etmp, rtmp, grad_coef, vtmp
    real(wp)                 :: work(1:3)
    real(wp)                 :: com_mov(1:3), com_ref(1:3)
    real(wp)                 :: rmsd
    real(wp)                 :: virial_local(1:3,1:3)
    integer                  :: iexpo
    integer                  :: i, j, k, iatm, ierr
    integer                  :: ikind

    real(wp),        pointer :: refcoord(:,:)
    real(wp),        pointer :: fit_refcoord(:,:), fit_work(:,:)
    real(wp),        pointer :: fit_mass(:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atomlist(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: expo(:)


    refcoord => enefunc%restraint_refcoord
    expo     => enefunc%restraint_exponent_func
    list     => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atomlist => enefunc%restraint_atomlist

    fit_refcoord => enefunc%fit_refcoord
    fit_work     => enefunc%fit_work
    fit_mass     => enefunc%fit_mass
    iexpo    =  expo(inum)

    ! from spdyn (2017/05/05)
    ierr  = 0

    if (enefunc%rpath_flag .and. enefunc%rpath_pos_func > 0) then
      if (allocated(enefunc%fit_work)) fit_work(:,:)=0.0_wp
      rot_matrix(:,:)=0.0_wp
      rot_matrix(1,1)=1.0_wp
      rot_matrix(2,2)=1.0_wp
      rot_matrix(3,3)=1.0_wp

      select case (enefunc%fitting_method)

      case (FittingMethodTR_ROT)

        if (enefunc%fitting_move == FittingMoveSYS) then
          fit_work(1:3,1:enefunc%num_atoms_ref) = &
                     coord(1:3,1:enefunc%num_atoms_ref) 
          call fit_trrot(enefunc%num_fitting, enefunc%fitting_atom,  &
                         fit_refcoord, fit_mass, coord, rot_matrix,  &
                         com_ref, com_mov, rmsd, ierr)

        else if (enefunc%fitting_move == FittingMoveREF) then
          fit_work(1:3,1:enefunc%num_atoms_ref) = &
                     refcoord(1:3,1:enefunc%num_atoms_ref) 
          call fit_trrot(enefunc%num_fitting, enefunc%fitting_atom,  &
                         coord, fit_mass, refcoord, rot_matrix,  &
                         com_ref, com_mov, rmsd, ierr)
        end if
         
      case (FittingMethodXYTR_ZROT)

        if (enefunc%fitting_move == FittingMoveSYS) then
          fit_work(1:3,1:enefunc%num_atoms_ref) = &
                     coord(1:3,1:enefunc%num_atoms_ref) 
          call fit_xytr_zrot2(enefunc%num_fitting, enefunc%fitting_atom,  &
                         fit_refcoord, coord, fit_mass, rot_matrix,  &
                         com_ref, com_mov, rmsd, ierr)

        else if (enefunc%fitting_move == FittingMoveREF) then
          fit_work(1:3,1:enefunc%num_atoms_ref) = &
                     refcoord(1:3,1:enefunc%num_atoms_ref) 
          call fit_xytr_zrot2(enefunc%num_fitting, enefunc%fitting_atom,  &
                         coord, refcoord, fit_mass, rot_matrix,  &
                         com_ref, com_mov, rmsd, ierr)
        end if

      case (FittingMethodNO)
         ! Do nothing

      end select

      if (ierr /= 0) then
        call error_msg("Compute_Energy_Positional_Restraint> fitting failed")
      end if
      if (enefunc%fitting_method /= FittingMethodNO)  then
        call transform(enefunc%num_atoms_ref, rot_matrix, com_ref, com_mov, &
                     fit_work)
      end if

    end if


    do j = 1, 4
      const(j-1) = enefunc%restraint_const(j,inum)
    end do

    eposi = 0.0_wp

    virial_local(1:3,1:3) = 0.0_wp

    do i = 1, numatoms(list(1,inum))

      iatm = atomlist(i,list(1,inum))

      ! Positional restraint energy: 
      !   E=(Kx(x-x0)^2+Ky(y-y0)^2+Kz(z-z0)^2)^(iexpo/2)
      !
      d(1:3) = coord(1:3,iatm) - refcoord(1:3,iatm)
      rtmp = const(1)*d(1)*d(1) + const(2)*d(2)*d(2) + const(3)*d(3)*d(3)

      if (iexpo == 2) then
        etmp      = const(0)*rtmp
        grad_coef = const(0)*2.0_wp
      else
        etmp      = const(0)*sqrt(rtmp)**iexpo
        grad_coef = const(0)*real(iexpo,wp)*(sqrt(rtmp)**(iexpo-2))
      end if

      eposi = eposi + etmp

      d2(1:3) = d(1:3)
      if (enefunc%rpath_sum_mf_flag .and. enefunc%rpath_pos_func > 0) then
        if ((enefunc%fitting_method == FittingMethodTR_ROT) .or. &
            (enefunc%fitting_method == FittingMethodXYTR_ZROT)) then
          d2(1:3) = fit_work(1:3,iatm) - refcoord(1:3,iatm)
        end if

        enefunc%stats_delta(3*i-2) = d2(1)
        enefunc%stats_delta(3*i-1) = d2(2)
        enefunc%stats_delta(3*i  ) = d2(3)

      end if

      if (calc_force) then
        ! gradient: dE/dX
        !
        work(1:3)       = grad_coef * const(1:3) * d(1:3)

        force(1:3,iatm) = force(1:3,iatm) - work(1:3)

        ! external virial
        !
        do j = 1, 3
          do k = j+1, 3
            vtmp = - coord(k,iatm) * work(j)
            virial_local(k,j) = virial_local(k,j) + vtmp
            virial_local(j,k) = virial_local(j,k) + vtmp
          end do
          vtmp = - coord(j,iatm) * work(j)
          virial_local(j,j) = virial_local(j,j) + vtmp
        end do

      end if

    end do

    if (enefunc%rpath_sum_mf_flag .and. enefunc%rpath_pos_func > 0) then
      enefunc%stats_count = enefunc%stats_count + 1
      do i = 1, enefunc%stats_dimension/3
        enefunc%stats_grad (1,i,3*i-2) = 1.0_wp
        enefunc%stats_grad (2,i,3*i-1) = 1.0_wp
        enefunc%stats_grad (3,i,3*i  ) = 1.0_wp
      end do
    end if

    if (calc_force) then
      if (.not. enefunc%pressure_position) then
        virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_local(1:3,1:3)
      else
        virial(1:3,1:3) = virial(1:3,1:3) + virial_local(1:3,1:3)
      end if
    end if

    return

  end subroutine compute_energy_positional_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_distance_restraint
  !> @brief        calculate distance restraint energy
  !! @authors      CK
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    inum      : pointer for restraint function
  !! @param[in]    const     : force constants
  !! @param[inout] force     : forces of target systems
  !! @param[inout] virial    : virial of target systems
  !! @param[out]   erestdist : distance restraint energy of target systems
  !! @param[out]   cv        : distance of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_distance_restraint(enefunc, coord, inum,      &
                                               calc_force, force, virial, &
                                               erestdist, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: erestdist
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                 :: r_sum, r_dif, coef, coefi, const(2), ref(2)
    real(wp)                 :: grad_coef, vtmp
    real(wp)                 :: work(1:3)
    integer                  :: numgrp, iexpo
    integer                  :: j, k, iatm, ind
    integer                  :: ndist, idist, id1, id2, id
    integer                  :: i_rpath_dim

    real(wp),        pointer :: wtind(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)
    integer,         pointer :: expoind(:,:)


    funcgrp  => enefunc%restraint_funcgrp 
    expo     => enefunc%restraint_exponent_func
    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    expoind  => enefunc%restraint_exponent_dist
    wtind    => enefunc%restraint_weight_dist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    coefgrp  => enefunc%restraint_wtmp
    com      => enefunc%restraint_wcom1
    dij      => enefunc%restraint_wcom2
    rij      => enefunc%restraint_wdrt

    numgrp   = funcgrp(inum)
    iexpo    = expo(inum)
    ndist    = int(numgrp/2)

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2)   = enefunc%restraint_ref(1:2,inum)

    !
    ! for steered MD
    if (inum == enefunc%steered_function) then
      ref(1:2)  = enefunc%target_value
    end if

    ! calculation of com
    !
    coefgrp(:,:) = 0.0_wp
    if (kind(inum) == RestraintsFuncDIST) then

      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = 1.0_wp/real(numatoms(grplist(j,inum)),wp)
        end do
      end do

    else

      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = masscoef(k,grplist(j,inum))
        end do
      end do

    end if

    do j = 1, numgrp
      com(1:3,j) = 0.0_wp

      do k = 1, numatoms(grplist(j,inum))
        iatm = atmlist(k,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + coord(1:3,iatm) * coefgrp(k,j)
      end do

    end do

    r_sum = 0.0_wp

    do idist = 1, ndist
      id1 = idist*2-1
      id2 = idist*2

      dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)

      rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) +  &
                            dij(2,idist)*dij(2,idist) +  &
                            dij(3,idist)*dij(3,idist))

      if (expoind(idist,inum) == 1) then
        r_sum = r_sum + wtind(idist,inum) * rij(idist)
      else
        r_sum = r_sum + wtind(idist,inum) * rij(idist)**(expoind(idist,inum))
      end if
    end do

    if (r_sum < ref(1)) then
      ind = 1
    else if (r_sum > ref(2)) then
      ind = 2
    else 
      ! ref(1) < r_sum < ref(2) : force and energy is zero 
      return
    end if

    r_dif = r_sum - ref(ind)
    cv    = r_sum
    if (iexpo == 2) then
      coef      = const(ind)
      erestdist = coef*r_dif*r_dif
      grad_coef = 2.0_wp * const(ind)*r_dif

    else
      coef      = const(ind) 
      erestdist = coef*r_dif**(iexpo)
      grad_coef = real(iexpo,wp) * const(ind)*(r_dif**(iexpo-1))

    end if
   

    ! Calculation of force
    !   coef  = const * (R - Ref)^(k-1)
    !   coefi = coef * u_i * e_n * r_n^(e_n -2)
    !   force(n,i) = coefi * dij(n) * sig(i) * coefgrp(i)   
    !
    if (calc_force) then
      do idist = 1, ndist

        if (rij(idist) < 1e-10_wp) then
          coefi = 0.0_wp
        else
          if (expoind(idist,inum) == 1) then
            coefi = grad_coef * wtind(idist,inum) / rij(idist)
          else
            coefi = expoind(idist,inum) * grad_coef  &
                   * rij(idist) ** (expoind(idist,inum)-2)
          end if
        end if

        id = 2*idist-1
        do j = 1, numatoms(grplist(id,inum))
          iatm = atmlist(j,grplist(id,inum))

          work(1:3) = - coefi * coefgrp(j,id) * dij(1:3,idist)

          force(1:3,iatm) = force(1:3,iatm) + work(1:3)

        end do

        id = 2*idist
        do j = 1, numatoms(grplist(id,inum))
          iatm = atmlist(j,grplist(id,inum))
          work(1:3) =   coefi * coefgrp(j,id) * dij(1:3,idist)
          force(1:3,iatm) = force(1:3,iatm) + work(1:3)
        end do

        ! virial
        ! 
        work(1:3) = coefi * dij(1:3,idist)

        do j = 1,3

          do k = j+1, 3
            vtmp = -work(k) * dij(j,idist)
            virial(k,j) = virial(k,j) + vtmp
            virial(j,k) = virial(j,k) + vtmp
          end do

          vtmp = -work(j) * dij(j,idist)
          virial(j,j) = virial(j,j) + vtmp
        end do

      end do
    end if


    ! rpath
    !
    if (enefunc%rpath_sum_mf_flag) then
      i_rpath_dim=enefunc%restraint_rpath_func(inum)
      if (i_rpath_dim >= 1) then
      if (i_rpath_dim == 1) enefunc%stats_count = enefunc%stats_count + 1
        enefunc%stats_delta(i_rpath_dim) = r_dif
        enefunc%stats_grad(1,1,i_rpath_dim) = dij(1, 1) / rij(1)
        enefunc%stats_grad(2,1,i_rpath_dim) = dij(2, 1) / rij(1)
        enefunc%stats_grad(3,1,i_rpath_dim) = dij(3, 1) / rij(1)
        enefunc%stats_grad(1,2,i_rpath_dim) = - dij(1, 1) / rij(1)
        enefunc%stats_grad(2,2,i_rpath_dim) = - dij(2, 1) / rij(1)
        enefunc%stats_grad(3,2,i_rpath_dim) = - dij(3, 1) / rij(1)
      end if
    end if

    return

  end subroutine compute_energy_distance_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle_restraint
  !> @brief        calculate angle restraint energy
  !! @authors      CK
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    inum      : pointer for restraint function
  !! @param[in]    const     : force constants
  !! @param[inout] force     : forces of target systems
  !! @param[inout] virial    : virial of target systems
  !! @param[out]   erestangle: angle restraint energy of target systems
  !! @param[out]   cv        : angle of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_restraint(enefunc, coord, inum,       &
                                            calc_force, force, virial,  &
                                            erestangle, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: erestangle
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                 :: dij(1:3), dkj(1:3)
    real(wp)                 :: rij, rkj, inv_rij2, inv_rkj2
    real(wp)                 :: rij2, rkj2, rijrkj,inv_rijrkj, vtmp
    real(wp)                 :: cos_t, t_dif, coef, theta, sin_t, grad_coef
    real(wp)                 :: v(1:3,1:3)
    real(wp)                 :: const(2), ref(2), work(1:6)
    real(wp)                 :: radref(1:2)
    integer                  :: numgrp, iexpo
    integer                  :: i, j, k, iatm, ind

    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)
    logical                  :: anglecheck


    funcgrp  => enefunc%restraint_funcgrp 
    expo     => enefunc%restraint_exponent_func
    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    coefgrp  => enefunc%restraint_wtmp
    com      => enefunc%restraint_wcom1

    numgrp   = funcgrp(inum)
    iexpo    = expo   (inum)

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2)   = enefunc%restraint_ref(1:2,inum)
    !
    ! for steered MD
    if (inum == enefunc%steered_function) then
      ref(1:2)   = enefunc%target_value
    end if

    coefgrp = 0.0_wp

    if (kind(inum) == RestraintsFuncANGLE) then
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = 1.0_wp/real(numatoms(grplist(j,inum)),wp)
        end do
      end do
    else
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = masscoef(k,grplist(j,inum))
        end do
      end do
    end if


    ! calculation of com
    !    
    com(1:3,1:3) = 0.0_wp

    do i = 1,3 
      do j = 1,numatoms(grplist(i,inum))
        iatm = atmlist(j,grplist(i,inum))
        com(1:3,i) = com(1:3,i) + coord(1:3,iatm) * coefgrp(j,i)
      end do
    end do
    
    radref(1:2) = ref(1:2)*RAD

    dij(1:3) = com(1:3,1) - com(1:3,2)
    dkj(1:3) = com(1:3,3) - com(1:3,2)

    rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
    rkj2  = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)

    rijrkj = sqrt(rij2*rkj2)

    inv_rijrkj = 1.0_wp /rijrkj
    inv_rij2 = 1.0_wp /rij2
    inv_rkj2 = 1.0_wp /rkj2

    cos_t = (dij(1)*dkj(1) + dij(2)*dkj(2) + dij(3)*dkj(3) )*inv_rijrkj
    cos_t  = min(  1.0_wp, cos_t )
    cos_t  = max( -1.0_wp, cos_t )
    theta = acos(cos_t)

    if (theta < radref(1)) then
      ind = 1
    else if (theta > radref(2)) then
      ind = 2
    else 
      return
    end if

    anglecheck = .true.
    do while(anglecheck)
      if (theta > pi)  theta = theta - 2.0_wp * pi
      if (theta < -pi) theta = theta + 2.0_wp * pi
      if (theta >= -pi  .and. theta <= pi)  anglecheck = .false.
    end do

    t_dif = theta - radref(ind)
    cv    = theta

    anglecheck = .true.
    do while(anglecheck)
      if (t_dif > pi)  t_dif = t_dif - 2.0_wp * pi
      if (t_dif < -pi) t_dif = t_dif + 2.0_wp * pi
      if (t_dif >= -pi  .and. t_dif <= pi)  anglecheck = .false.
    end do

    if (iexpo == 2) then
      coef       = const(ind) 
      erestangle = coef * t_dif * t_dif
      grad_coef  = 2.0_wp * const(ind) * t_dif
    else
      coef       = const(ind)
      erestangle = coef * t_dif**(iexpo)
      grad_coef  = real(iexpo,wp) * const(ind)*(t_dif**(iexpo-1))
    end if


    ! compute force
    !
    if (calc_force) then    
      sin_t = (1.0_wp-cos_t*cos_t)
      sin_t = sqrt(sin_t)
      sin_t  = max( EPS, sin_t )
      grad_coef = -grad_coef / sin_t

      work(1:3) = grad_coef * (dkj(1:3) * inv_rijrkj -  &
                               dij(1:3) * cos_t * inv_rij2)

      work(4:6) = grad_coef * (dij(1:3) * inv_rijrkj -  &
                               dkj(1:3) * cos_t * inv_rkj2)

      ! virial from angle
      !
      do j = 1, 3

        do k = j+1, 3
          vtmp = -(dij(k)*work(j) + dkj(k)*work(j+3))
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do

        vtmp = -(dij(j)*work(j) + dkj(j)*work(j+3))
        virial(j,j) = virial(j,j) + vtmp
      end do

      do j = 1, numatoms(grplist(1,inum))
        iatm = atmlist(j,grplist(1,inum))
        force(1:3,iatm) = force(1:3,iatm) - work(1:3)*coefgrp(j,1)
      end do

      do j = 1, numatoms(grplist(2,inum))
        iatm = atmlist(j,grplist(2,inum))
        force(1:3,iatm) = force(1:3,iatm)  &
                          + (work(1:3) + work(4:6)) *coefgrp(j,2)
      end do

      do j = 1, numatoms(grplist(3,inum))
        iatm = atmlist(j,grplist(3,inum))
        force(1:3,iatm) = force(1:3,iatm) - work(4:6)*coefgrp(j,3)
      end do
    end if

    return

  end subroutine compute_energy_angle_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihedral_restraint
  !> @brief        calculate dihedral restraint energy
  !! @authors      CK
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    inum      : pointer for restraint function
  !! @param[in]    const     : force constants
  !! @param[inout] force     : forces of target systems
  !! @param[inout] virial    : virial of target systems
  !! @param[out]   erestdihed: dihed restraint energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihedral_restraint(enefunc, coord, inum,      &
                                               calc_force, force, virial, &
                                               erestdihed, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: erestdihed
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, coef
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi
    real(wp)                 :: cosdif, sindif, theta, vtmp
    real(wp)                 :: const(2), ref(2), work(1:9)
    real(wp)                 :: radref(1:2)
    real(wp)                 :: grad(1:9), v(1:3,1:3)
    integer                  :: aindex(1:4)
    integer                  :: numgrp, iexpo
    integer                  :: i, j, k, iatm, ind
    integer                  :: i_rpath_dim

    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    funcgrp  => enefunc%restraint_funcgrp 
    expo     => enefunc%restraint_exponent_func
    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    coefgrp  => enefunc%restraint_wtmp
    com      => enefunc%restraint_wcom1

    numgrp   = funcgrp(inum)
    iexpo    = expo   (inum)

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2)   = enefunc%restraint_ref(1:2,inum)
    !
    ! for steered MD
    if (inum == enefunc%steered_function) then
      ref(1:2)  = enefunc%target_value
    end if

    coefgrp = 0.0_wp

    if (kind(inum) == RestraintsFuncDIHED) then
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = 1.0_wp/real(numatoms(grplist(j,inum)),wp)
        end do
      end do
    else
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = masscoef(k,grplist(j,inum))
        end do
      end do
    end if

    ! calculation of com
    !
    com(1:3,1:4)=0.0_wp

    do i = 1,4 
      do j = 1,numatoms(grplist(i,inum))
        iatm = atmlist(j,grplist(i,inum))
        com(1:3,i) = com(1:3,i) + coord(1:3,iatm) * coefgrp(j,i)
      end do
      aindex(i) = i
    end do
    radref(1:2) = ref(1:2)*RAD

    call calculate_dihedral(aindex, com, cos_dih, sin_dih, grad, v)

    if (cos_dih > 1.0E-1_wp) then
      theta = asin(sin_dih)
    else
      theta = sign(1.0_wp,sin_dih)*acos(cos_dih)
    end if
    cv = theta
   
    if (theta < radref(1)) then
      ind = 1
    else if (theta > radref(2)) then
      ind = 2
    else 
      return
    end if
    cospha = cos(radref(ind))
    sinpha = sin(radref(ind))

    cosdif = cos_dih*cospha + sin_dih*sinpha
    sindif = cos_dih*sinpha - sin_dih*cospha

    if (cosdif > 1.0E-1_wp) then
      diffphi = asin(sindif)
    else
      diffphi = sign(1.0_wp,sindif)*acos(cosdif)
    end if

    if (iexpo == 2) then
      coef       = const(ind) 
      erestdihed = coef * diffphi * diffphi
      grad_coef  = 2.0_wp * const(ind) * diffphi
    else
      coef       = const(ind)
      erestdihed = coef * diffphi **(iexpo)
      grad_coef  = real(iexpo,wp) * const(ind)*(diffphi**(iexpo-1))
    end if
   

    ! compute force
    !
    if (calc_force) then

      if (abs(theta) < 1e-10_wp) then
        grad_coef = 0.0_wp
      end if
      work(1:9) = grad_coef*grad(1:9)

      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)

          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do

        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

      do j = 1, numatoms(grplist(1,inum))
        iatm = atmlist(j,grplist(1,inum))
        force(1:3,iatm) = force(1:3,iatm) - work(1:3)*coefgrp(j,1)
      end do

      do j = 1, numatoms(grplist(2,inum))
        iatm = atmlist(j,grplist(2,inum))
        force(1:3,iatm) = force(1:3,iatm) + (work(1:3) - work(4:6))*coefgrp(j,2)
      end do

      do j = 1, numatoms(grplist(3,inum))
        iatm = atmlist(j,grplist(3,inum))
        force(1:3,iatm) = force(1:3,iatm) + (work(4:6) + work(7:9))*coefgrp(j,3)
      end do

      do j = 1, numatoms(grplist(4,inum))
        iatm = atmlist(j,grplist(4,inum))
        force(1:3,iatm) = force(1:3,iatm) - work(7:9)*coefgrp(j,4)
      end do
    end if

    if (enefunc%rpath_sum_mf_flag) then
      i_rpath_dim=enefunc%restraint_rpath_func(inum)
      if (i_rpath_dim >= 1) then
        if (i_rpath_dim == 1) enefunc%stats_count = enefunc%stats_count + 1
        theta = theta/RAD - ref(1)
        if (theta > 180.0_wp) then
          theta = theta - 360.0_wp
        else if (theta < -180.0_wp) then
          theta = theta + 360.0_wp
        end if
        enefunc%stats_delta(i_rpath_dim) = theta * RAD
        enefunc%stats_grad(1,1,i_rpath_dim) = grad(1) 
        enefunc%stats_grad(2,1,i_rpath_dim) = grad(2) 
        enefunc%stats_grad(3,1,i_rpath_dim) = grad(3) 
        enefunc%stats_grad(1,2,i_rpath_dim) = - grad(1) + grad(4) 
        enefunc%stats_grad(2,2,i_rpath_dim) = - grad(2) + grad(5) 
        enefunc%stats_grad(3,2,i_rpath_dim) = - grad(3) + grad(6) 
        enefunc%stats_grad(1,3,i_rpath_dim) = - grad(4) - grad(7) 
        enefunc%stats_grad(2,3,i_rpath_dim) = - grad(5) - grad(8) 
        enefunc%stats_grad(3,3,i_rpath_dim) = - grad(6) - grad(9) 
        enefunc%stats_grad(1,4,i_rpath_dim) = grad(7) 
        enefunc%stats_grad(2,4,i_rpath_dim) = grad(8) 
        enefunc%stats_grad(3,4,i_rpath_dim) = grad(9)       
      end if
    end if

    return

  end subroutine compute_energy_dihedral_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_rmsd_restraint
  !> @brief        calculate position restraint energy with RMSD
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    rep_md  : flag for replica MD or not
  !! @param[in]    const   : force constants
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   virial_ext  : external virial term of target systems
  !! @param[out]   eposi   : position constraint energy of target systems
  !! @param[out]   cv      : CV of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_rmsd_restraint(enefunc, coord, inum,              &
                                           calc_force, rep_md, force, virial, &
                                           virial_ext, eposi, cv)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    logical,                 intent(in)    :: rep_md
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)
    real(wp),                intent(out)   :: eposi
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                 :: r_dif, coef, grad_coef, vtmp
    real(wp)                 :: rmsd, Krmsd, factor, totmass
    real(wp)                 :: const(2), ref(2), work(1:3)
    real(wp)                 :: dcom(1:3), dref(1:3)
    real(wp)                 :: dadd(1:3), dsub(1:3)
    real(wp)                 :: crd_com(3), crd_com_ref(3)
    real(wp)                 :: rot_matrix(1:3,1:3)
    real(wp)                 :: virial_local(1:3,1:3)
    integer                  :: i, j, k, igrp, iatm, ind, ierr
    integer                  :: natmgrp, iexpo

    real(wp),        pointer :: refcoord(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: fit_mass(:)
    real(wp),        pointer :: totmass_grp(:)
    real(wp),        pointer :: masstmp(:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: expo(:)


    expo           => enefunc%restraint_exponent_func
    grplist        => enefunc%restraint_grouplist
    numatoms       => enefunc%restraint_numatoms
    atmlist        => enefunc%restraint_atomlist
    kind           => enefunc%restraint_kind
    masscoef       => enefunc%restraint_masscoef
    masstmp        => enefunc%restraint_masstmp
    totmass_grp    => enefunc%restraint_totmass_group
    refcoord       => enefunc%restraint_refcoord
    coefgrp        => enefunc%restraint_wtmp
    fit_mass       => enefunc%fit_mass

    Krmsd    =  enefunc%rmsd_force

    igrp     = grplist(1,inum)
    natmgrp  = numatoms(igrp)
    iexpo    = expo(inum)

    virial_local(1:3,1:3) = 0.0_wp

    rot_matrix(:,:) = 0.0_wp
    rot_matrix(1,1) = 1.0_wp
    rot_matrix(2,2) = 1.0_wp
    rot_matrix(3,3) = 1.0_wp

    masstmp(:)   = 0.0_wp
    coefgrp(:,:) = 0.0_wp
    do i = 1, natmgrp
      iatm = atmlist(i,igrp)
      if (kind(inum) == RestraintsFuncRMSD) then
        masstmp(iatm) = 1.0_wp/real(natmgrp,wp)
        coefgrp(i,1)  = 1.0_wp/real(natmgrp,wp)
      else
        masstmp(iatm) = masscoef(i,igrp)
        coefgrp(i,1)  = masscoef(i,igrp)
      end if
    end do
    if (kind(inum) == RestraintsFuncRMSD) then
       totmass = real(natmgrp,wp)
    else
       totmass = totmass_grp(igrp)
    end if

    select case (enefunc%fitting_method)
   
    case (FittingMethodTR_ROT)
   
        call fit_trrot(enefunc%num_fitting, enefunc%fitting_atom,  &
                       coord, fit_mass, refcoord, rot_matrix,  &
                       crd_com, crd_com_ref, rmsd, ierr)
       
    case (FittingMethodXYTR_ZROT)
   
        call fit_xytr_zrot2(enefunc%num_fitting, enefunc%fitting_atom,  &
                       coord, refcoord, fit_mass, rot_matrix,  &
                       crd_com, crd_com_ref, rmsd, ierr)
   
    case (FittingMethodNO)
       ! Do nothing
   
    end select

    rmsd = 0.0_wp
    do i = 1, natmgrp
      iatm = atmlist(i,igrp)
      dcom(1:3) = coord(1:3,iatm)    - crd_com(1:3) 
      dref(1:3) = refcoord(1:3,iatm) - crd_com_ref(1:3)
      dadd(1:3) = rot_matrix(1:3,1)*dref(1) + &
                  rot_matrix(1:3,2)*dref(2) + &
                  rot_matrix(1:3,3)*dref(3)
      dsub(1:3) = dcom(1:3)-dadd(1:3)
      rmsd      = rmsd + coefgrp(i,1)* &
                  (dsub(1)*dsub(1) +dsub(2)*dsub(2) + dsub(3)*dsub(3))
    end do

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2) = enefunc%restraint_ref(1:2,inum)

    rmsd = sqrt(rmsd)
    !
    ! for steered & targeted MD
    !
    if (inum == enefunc%target_function .or.  &
        inum == enefunc%steered_function) then
      ref(1:2)   = enefunc%target_value
    end if
    cv   = rmsd

    if (rmsd < ref(1)) then
      ind = 1
    else if (rmsd > ref(2)) then
      ind = 2
    else 
      return
    end if

    ! for steered MD or normal restraint
    if (inum /= enefunc%target_function) then
      r_dif = rmsd - ref(ind)

      if (iexpo == 2) then
        coef  = const(ind) 
        eposi = coef * r_dif  * r_dif
        grad_coef  = 2.0_wp * const(ind) * r_dif
      else
        coef  = const(ind) 
        eposi = coef * r_dif ** (iexpo)
        grad_coef  = real(iexpo,wp) * const(ind) * r_dif ** (iexpo-1) 
      end if
     
      ! Compute gradient 
      !
      if (calc_force) then
        if (rmsd < 1e-10_wp) then
          grad_coef = 0.0_wp
        else
          grad_coef = grad_coef/rmsd
        end if
     
        do i = 1, natmgrp
          iatm = atmlist(i,igrp)
     
          ! gradient: dE/dX
          !
          dcom(1:3) = coord(1:3,iatm) - crd_com(1:3)
          dref(1:3) = refcoord(1:3,iatm) - crd_com_ref(1:3)
          dadd(1:3) = matmul(rot_matrix(1:3,1:3), dref(1:3))
          dsub(1:3) = dcom(1:3)-dadd(1:3)
     
          ! gradient: dE/dX
          !
          work(1:3) = grad_coef * coefgrp(i,1) * dsub(1:3)
     
          force(1:3,iatm) = force(1:3,iatm) - work(1:3)
     
          ! external virial
          !
          do j = 1, 3
            do k = j+1, 3
     
              vtmp = -coord(k,iatm) * work(j)
              virial_local(k,j) = virial_local(k,j) + vtmp
              virial_local(j,k) = virial_local(j,k) + vtmp
            end do
     
            vtmp = -coord(j,iatm) * work(j)
            virial_local(j,j) = virial_local(j,j) + vtmp
          end do
        end do
      end if

    else
    ! for targeted MD
      eposi = 0.0_wp
      do i = 1, natmgrp
        iatm = atmlist(i,igrp)
        dcom(1:3) = coord(1:3,iatm) - crd_com(1:3)
        dref(1:3) = refcoord(1:3,iatm) - crd_com_ref(1:3)
        dadd(1:3) = matmul(rot_matrix(1:3,1:3), dref(1:3))
        dsub(1:3) = dcom(1:3) - dadd(1:3)
        factor    = Krmsd*totmass*coefgrp(i,1)*(1.0_wp-ref(1)/rmsd)
        eposi  = eposi + factor*dsub(1)*dsub(1) &
                       + factor*dsub(2)*dsub(2) &
                       + factor*dsub(3)*dsub(3)
        factor = 2.0_wp*factor
      
        ! gradient: dE/dX
        !
        if (calc_force) then
          work(1:3) = factor * dsub(1:3)
      
          force(1:3,iatm)  = force(1:3,iatm) - work(1:3)
      
          ! external virial
          !
          do j = 1, 3
            do k = j+1, 3
      
              vtmp = -coord(k,iatm) * work(j)
              virial_local(k,j) = virial_local(k,j) + vtmp
              virial_local(j,k) = virial_local(j,k) + vtmp
            end do
      
            vtmp = -coord(j,iatm) * work(j)
            virial_local(j,j) = virial_local(j,j) + vtmp
          end do
        end if
      end do

    end if

    if (calc_force) then
      if (.not. enefunc%pressure_rmsd) then
        virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_local(1:3,1:3)
      else
        virial(1:3,1:3) = virial(1:3,1:3) + virial_local(1:3,1:3)
      end if
    end if

    return

  end subroutine compute_energy_rmsd_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_pc_restraint
  !> @brief        calculate position restraint energy with PC
  !! @authors      YK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   erestpc : principal component restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_pc_restraint(enefunc, coord, inum, &
                                         calc_force, force, virial, erestpc, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: erestpc
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                 :: const, ref, tot_mass, work(1:3)
    real(wp)                 :: crd_com(1:3), crd_ref_com(1:3)
    real(wp)                 :: pc, grad_erestpc
    real(wp)                 :: sym_matrix(1:4,1:4)
    real(wp)                 :: rot_matrix(1:3,1:3), mass
    real(wp)                 :: eval(1:4), evec(1:4), lwork(1:11)
    real(wp)                 :: dref(1:3), dmov(1:3), dsub(1:3), dadd(1:3)
    real(wp)                 :: diff(1:3), diff_fit(1:3), rmsd_pc, rmsd_pc_fit
    real(wp)                 :: rmsd
    integer                  :: i, j, k, iatm, natmgrp, igrp
    integer                  :: pc_mode_sta, pc_mode_end
    integer                  :: ierr

    real(wp),        pointer :: refcoord(:,:), enefunc_pc_mode(:)
    real(wp),        pointer :: fit_mass(:)
    real(wp),        pointer :: totmass_grp(:)
    real(wp),        pointer :: pc_mode(:), pc_mode_fit(:), grad_pc(:,:)
    real(wp),        pointer :: crd(:), crd_ref(:), crd_ref_fit(:)
    integer,         pointer :: numatoms(:), atmlist(:,:), kind(:)
    integer,         pointer :: grplist(:,:)


    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind

    refcoord       => enefunc%restraint_refcoord
    fit_mass       => enefunc%fit_mass
    totmass_grp    => enefunc%restraint_totmass_group
    pc_mode        => enefunc%pc_mode_temp
    pc_mode_fit    => enefunc%pc_mode_fit
    grad_pc        => enefunc%grad_pc
    crd            => enefunc%pc_crd
    crd_ref        => enefunc%pc_crd_ref
    crd_ref_fit    => enefunc%pc_crd_ref_fit

    enefunc_pc_mode => enefunc%pc_mode

    igrp     = grplist(1,inum)
    natmgrp  = numatoms(igrp)
    const    = enefunc%restraint_const(1,inum)
    ref      = enefunc%restraint_ref  (1,inum)

    if (kind(inum) == RestraintsFuncPCCOM) then
      tot_mass = totmass_grp(igrp)
    else
      tot_mass = real(natmgrp,wp)
    end if

    !mass     = 1.0_wp
    rot_matrix(:,:) = 0.0_wp
    rot_matrix(1,1) = 1.0_wp
    rot_matrix(2,2) = 1.0_wp
    rot_matrix(3,3) = 1.0_wp

    ! allocate
    !allocate(crd(1:natmgrp*3), crd_ref(1:natmgrp*3), crd_ref_fit(1:natmgrp*3),&
    !         pc_mode(1:natmgrp*3), pc_mode_fit(1:natmgrp*3),&
    !         grad_pc(1:3,1:natmgrp))

    ! center of mass
    crd_com(:)     = 0.0_wp
    crd_ref_com(:) = 0.0_wp
    !tot_mass       = 0.0_wp

    select case (enefunc%fitting_method)
   
    case (FittingMethodTR_ROT)
   
      call fit_trrot(enefunc%num_fitting, enefunc%fitting_atom,  &
                     coord, fit_mass, refcoord, rot_matrix,  &
                     crd_com, crd_ref_com, rmsd, ierr)
       
    case (FittingMethodXYTR_ZROT)
   
      call fit_xytr_zrot2(enefunc%num_fitting, enefunc%fitting_atom,  &
                     coord, refcoord, fit_mass, rot_matrix,  &
                     crd_com, crd_ref_com, rmsd, ierr)
   
    case (FittingMethodNO)
       ! Do nothing
   
    end select

    ! pc
    pc_mode_sta = (enefunc%restraint_mode(inum)-1)*3*natmgrp + 1
    pc_mode_end =  enefunc%restraint_mode(inum)   *3*natmgrp
    pc_mode(1:3*natmgrp) = enefunc_pc_mode(pc_mode_sta:pc_mode_end)

    rmsd_pc     = 0.0_wp
    rmsd_pc_fit = 0.0_wp
    do i = 1, natmgrp
      iatm = atmlist(i,grplist(1,inum))
      crd((i-1)*3+1:i*3)     = coord(1:3,iatm)    - crd_com(1:3)
      crd_ref((i-1)*3+1:i*3) = refcoord(1:3,iatm) - crd_ref_com(1:3)
      crd_ref_fit((i-1)*3+1:i*3) =                                  &
              matmul(rot_matrix(1:3,1:3),crd_ref((i-1)*3+1:i*3))
      pc_mode_fit((i-1)*3+1:i*3) =                                  &
              matmul(rot_matrix(1:3,1:3),pc_mode((i-1)*3+1:i*3))
    end do
   
    rmsd_pc     = sqrt(rmsd_pc/tot_mass)
    rmsd_pc_fit = sqrt(rmsd_pc_fit/tot_mass)

    pc = dot_product(pc_mode_fit(1:3*natmgrp),&
                    (crd(1:3*natmgrp)-crd_ref_fit(1:3*natmgrp)))
    cv = pc

    ! energy of pc restraint
    erestpc = const * (pc - ref)**2

    ! gradient (dE/dx = dE/dpc * dpc/dx)
    ! dE/dpc
    grad_erestpc = 2.0_wp * const * (pc - ref)

    ! dpc/dx and dE/dx
    if (calc_force) then
      do i = 1, natmgrp
        iatm = atmlist(i,grplist(1,inum))
        grad_pc(1:3,i)  = pc_mode_fit((i-1)*3+1:3*i)
        work(1:3)       = grad_erestpc * grad_pc(1:3,i)
        force(1:3,iatm) = force(1:3,iatm) - work(1:3)
      end do
    
      ! virial (temporary)
      virial(1:3,1:3) = 0.0_wp
    end if

    if (enefunc%rpath_sum_mf_flag) then
      if (inum == 1) enefunc%stats_count = enefunc%stats_count + 1
      enefunc%stats_delta(inum) = pc - ref
      do i = 1, natmgrp
        enefunc%stats_grad(1:3,i,inum) = pc_mode_fit((i-1)*3+1:3*i)
      end do
    end if

    return

  end subroutine compute_energy_pc_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_energy_rg_restraint
  !> @brief      calculate rg restraint energy and forces
  !! @authors    MK
  !! @param[inout]  enefunc : potential energy functions [str]
  !! @param[in]  coord      : coordinates of target systems [dble]
  !! @param[in]  inum       : coordinates of target systems [dble]
  !! @param[in]  calc_force : coordinates of target systems [dble]
  !! @param[out] force      : forces of target systems [dble]
  !! @param[out] virial     : of target systems [dble]
  !! @param[out] eneval     : energy of target systems [dble]
  !! @param[out] cv         : Rg [dble]
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine compute_energy_rg_restraint(enefunc, coord, inum, calc_force, &
                                         force, virial, eneval, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: eneval
    real(wp),                intent(out)   :: cv

    ! local variables
    real(wp)                               :: com(3)
    real(wp)                               :: const, ref
    real(wp)                               :: msd, rg, rg_dif, vtmp
    real(wp)                               :: coef_force, local_coef
    real(wp)                               :: dcom(3), work(3)
    real(wp)                               :: flat_rad
    integer                                :: igrp
    integer                                :: natmgrp
    integer                                :: iatm, iexpo
    integer                                :: i, j, k
    logical                                :: caging

    ! pointers
    ! mass / total mass (of group)
    real(wp),                      pointer :: masscoef(:,:)
    real(wp),                      pointer :: masstmp(:)
    integer,                       pointer :: numatoms(:)
    integer,                       pointer :: grplist(:,:)
    integer,                       pointer :: atmlist(:,:)
    integer,                       pointer :: kind(:)
    integer,                       pointer :: expo(:)


    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    masstmp  => enefunc%restraint_masstmp
    expo     => enefunc%restraint_exponent_func

    igrp    = grplist(1,inum)
    natmgrp = numatoms(igrp)
    iexpo   = expo(inum)

    const    = enefunc%restraint_const(1,inum)
    ref      = enefunc%restraint_ref(1,inum)
    caging   = enefunc%restraint_caging(inum)
    flat_rad = enefunc%restraint_flat_radius(inum)

    do i = 1, natmgrp
      iatm = atmlist(i,igrp)
      if (kind(inum) == RestraintsFuncRGWOMASS) then
        masstmp(iatm) = 1.0_wp / real(natmgrp,wp)
      else
        masstmp(iatm) = masscoef(i,igrp)
      end if
    end do

    ! calculate center of mass (or geometric center)
    com(1:3) = 0.0_wp
    do i = 1, natmgrp
      iatm = atmlist(i,grplist(1,inum))
      com(1:3) = com(1:3) + masstmp(iatm) * coord(1:3,iatm)
    end do

    ! get Rg
    rg = 0.0_wp
    do i = 1, natmgrp
      iatm = atmlist(i,grplist(1,inum))
      msd = (coord(1,iatm) - com(1)) ** 2 + &
            (coord(2,iatm) - com(2)) ** 2 + &
            (coord(3,iatm) - com(3)) ** 2
      rg = rg + masstmp(iatm) * msd
    end do
    rg = sqrt( rg )
    cv = rg ! for output

    ! set constraint index
    if (rg < ref .and. caging) then
      eneval = 0.0_wp
      return
    end if

    ! energy
    rg_dif = rg - ref
    if (abs(rg_dif) < flat_rad) then
      eneval = 0.0_wp
      return
    end if
    if (rg_dif < 0.0_wp) then
      rg_dif = rg_dif + flat_rad
    else
      rg_dif = rg_dif - flat_rad
    end if
    eneval = const * (rg_dif ** iexpo)

    ! force
    if (rg < 1.0e-10_wp) then
      coef_force = 0.0_wp
    else
      coef_force = const * real(iexpo,wp) * ( rg_dif ** ( iexpo - 1 ) ) / rg
    end if

    if (calc_force) then
      do i = 1, natmgrp
        iatm = atmlist(i,grplist(1,inum))
        local_coef = coef_force * ( 1.0_wp - masstmp(iatm) ) * masstmp(iatm)
        dcom(1:3) = coord(1:3,iatm) - com(1:3)
        work(1:3) = local_coef * dcom(1:3)
        force(1:3,iatm) = force(1:3,iatm) - work(1:3)
        ! virial
        do j = 1, 3
          do k = j+1, 3
            vtmp = -dcom(k) * work(j)
            virial(k,j) = virial(k,j) + vtmp
            virial(j,k) = virial(j,k) + vtmp
          end do
          vtmp = -dcom(j) * work(j)
          virial(j,j) = virial(j,j) + vtmp
        end do
      end do
    end if

    return

  end subroutine compute_energy_rg_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_repulsive_restraint
  !> @brief        calculate repulsive restraint energy
  !! @authors      HO
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[in]    boundary  : information of boundary
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    inum      : pointer for restraint function
  !! @param[in]    const     : force constants
  !! @param[inout] force     : forces of target systems
  !! @param[inout] virial    : virial of target systems
  !! @param[out]   erepul    : repulsive restraint energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_repulsive_restraint(enefunc, boundary, coord, inum, &
                                                calc_force, force, virial, &
                                                erepul)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: erepul

    ! local variables
    real(wp)                 :: r_dif, const(2), ref(2)
    real(wp)                 :: vtmp
    real(wp)                 :: work(1:3)
    real(wp)                 :: coefi
    real(wp)                 :: box_size(3)
    integer                  :: numgrp, iexpo
    integer                  :: j, k, iatm, ind
    integer                  :: ndist, idist, id1, id2, id
    integer                  :: irpath_num

    real(wp), allocatable    :: grad_coef(:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    funcgrp  => enefunc%restraint_funcgrp
    expo     => enefunc%restraint_exponent_func
    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    coefgrp  => enefunc%restraint_wtmp
    com      => enefunc%restraint_rcom1
    dij      => enefunc%restraint_rcom2
    rij      => enefunc%restraint_rdrt
    box_size(1) = boundary%box_size_x
    box_size(2) = boundary%box_size_y
    box_size(3) = boundary%box_size_z

    numgrp   = funcgrp(inum)
    iexpo    = expo(inum)
    ndist    = int(numgrp*(numgrp-1)/2)

    !
    allocate(grad_coef(ndist))

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2)   = enefunc%restraint_ref(1:2,inum)

    ! 
    ! check rpath index for each replica function
    irpath_num = enefunc%restraint_replica_index(inum)

    !
    ! for steered MD
    if (inum == enefunc%steered_function) then
      ref(1:2)   = enefunc%target_value
    end if

    ! calculation of com
    !
    coefgrp(:,:) = 0.0_wp
    if (kind(inum) == RestraintsFuncREPUL) then
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = 1.0_wp/real(numatoms(grplist(j,inum)),wp)
        end do
      end do
    else
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = masscoef(k,grplist(j,inum))
        end do
      end do
    end if
    do j = 1, numgrp
      com(1:3,j) = 0.0_wp
      do k = 1, numatoms(grplist(j,inum))
        iatm = atmlist(k,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + coord(1:3,iatm) * coefgrp(k,j)
      end do
    end do


    ! Calculation of distance between two group
    !     U(r) = const(1) * (r - Ref(1))^k : r < Ref(1)
    !     U(r) = 0                               : Ref(1) <= r
    !
    !           ndist
    !   erepul = sum U(r_n)
    !             n
    !
    !     r_n  = |rcom_i - rcom_j|
    !
    !            N_i
    !   rcom_i = sum w_l * r_l ; if w_l = mass_l/totmass: mass weighted,
    !             l in group-i ;    w_l = 1     / N_i   : average
    !
    !            N_j
    !   rcom_j = sum w_m * r_m ; if w_m = mass_m/totmass: mass weighted,
    !             m in group-j ;    w_m = 1     / N_j   : average
    !
    ! Calculation of force
    ! grad_coef(n) = k*const(1)*(r_n - Ref(1))^(k-1)
    !   force(n,i) = grad_coef(n) * r_n^(-1) * dij(n) * sig(i) * coefgrp(i)   
    !

    idist = 0
    do id1 = 1, numgrp
      do id2 = id1+1, numgrp
        idist = idist + 1
        dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
        if (boundary%type == BoundaryTypePBC) then
            dij(1:3,idist) = dij(1:3,idist) - box_size(1:3)* &
                             anint(dij(1:3,idist)/box_size(1:3))
        end if
        rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) +  &
                              dij(2,idist)*dij(2,idist) +  &
                              dij(3,idist)*dij(3,idist))
      end do
    end do

    erepul = 0.0_wp
    do idist = 1, ndist
      if (rij(idist) < ref(1)) then
        r_dif = rij(idist) - ref(1)
        if (iexpo == 2) then
          erepul           = erepul + const(1)*r_dif*r_dif
          grad_coef(idist) = 2.0_wp*const(1)*r_dif
        else
          erepul           = erepul + const(1)*r_dif**(iexpo)
          grad_coef(idist) = real(iexpo,wp)*const(1)*r_dif**(iexpo-1)
        end if
      else
        grad_coef(idist) = 0.0_wp
      end if
    end do
   
    if (calc_force) then
      idist = 0
      do id1 = 1, numgrp
        do id2 = id1+1, numgrp
          idist = idist + 1

          if (rij(idist) < 1e-10_wp) then
            coefi = 0.0_wp
          else
            coefi = grad_coef(idist) * 1.0_wp / rij(idist)
          end if

          do j = 1, numatoms(grplist(id1,inum))
            iatm = atmlist(j,grplist(id1,inum))
            work(1:3) = - coefi * coefgrp(j,id1) * dij(1:3,idist)
            force(1:3,iatm) = force(1:3,iatm) + work(1:3)
          end do

          do j = 1, numatoms(grplist(id2,inum))
            iatm = atmlist(j,grplist(id2,inum))
            work(1:3) =   coefi * coefgrp(j,id2) * dij(1:3,idist)
            force(1:3,iatm) = force(1:3,iatm) + work(1:3)
          end do

          ! virial
          ! 
          work(1:3) = coefi * dij(1:3,idist)
          do j = 1,3
            do k = j+1, 3
              vtmp = -work(k) * dij(j,idist)
              virial(k,j) = virial(k,j) + vtmp
              virial(j,k) = virial(j,k) + vtmp
            end do
            vtmp = -work(j) * dij(j,idist)
            virial(j,j) = virial(j,j) + vtmp
          end do

        end do
      end do
    end if

    deallocate(grad_coef)

    return

  end subroutine compute_energy_repulsive_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_flatbottom_restraint
  !> @brief        calculate flat-bottom restraint energy
  !! @authors      HO
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[in]    boundary  : information of boundary
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    inum      : pointer for restraint function
  !! @param[in]    const     : force constants
  !! @param[inout] force     : forces of target systems
  !! @param[inout] virial    : virial of target systems
  !! @param[out]   efb       : flat-bottom restraint energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_flatbottom_restraint(enefunc, boundary, coord, &
                                                 inum, calc_force, force,  &
                                                 virial, efb)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(out)   :: efb

    ! local variables
    real(wp)                 :: r_dif, const(2), ref(2)
    real(wp)                 :: vtmp
    real(wp)                 :: work(1:3)
    real(wp)                 :: coefi
    real(wp)                 :: box_size(3)
    integer                  :: numgrp, iexpo
    integer                  :: j, k, iatm, ind
    integer                  :: ndist, idist, id1, id2, id
    integer                  :: irpath_num

    real(wp), allocatable    :: grad_coef(:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: numatoms(:)
    integer,         pointer :: atmlist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    funcgrp  => enefunc%restraint_funcgrp 
    expo     => enefunc%restraint_exponent_func
    grplist  => enefunc%restraint_grouplist
    numatoms => enefunc%restraint_numatoms
    atmlist  => enefunc%restraint_atomlist
    kind     => enefunc%restraint_kind
    masscoef => enefunc%restraint_masscoef
    coefgrp  => enefunc%restraint_wtmp
    com      => enefunc%restraint_rcom1
    dij      => enefunc%restraint_rcom2
    rij      => enefunc%restraint_rdrt
    box_size(1) = boundary%box_size_x
    box_size(2) = boundary%box_size_y
    box_size(3) = boundary%box_size_z

    numgrp   = funcgrp(inum)
    iexpo    = expo(inum)
    ndist    = int(numgrp*(numgrp-1)/2)

    !
    allocate(grad_coef(ndist))

    const(1:2) = enefunc%restraint_const(1:2,inum)
    ref(1:2)   = enefunc%restraint_ref(1:2,inum)

    ! 
    ! check rpath index for each replica function
    irpath_num = enefunc%restraint_replica_index(inum)

    !
    ! for steered MD
    if (inum == enefunc%steered_function) then
      ref(1:2)   = enefunc%target_value
    end if

    ! calculation of com
    !
    coefgrp(:,:) = 0.0_wp
    if (kind(inum) == RestraintsFuncFB) then
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = 1.0_wp/real(numatoms(grplist(j,inum)),wp)
        end do
      end do
    else
      do j = 1, numgrp
        do k = 1, numatoms(grplist(j,inum))
          coefgrp(k,j) = masscoef(k,grplist(j,inum))
        end do
      end do
    end if
    do j = 1, numgrp
      com(1:3,j) = 0.0_wp
      do k = 1, numatoms(grplist(j,inum))
        iatm = atmlist(k,grplist(j,inum))
        com(1:3,j) = com(1:3,j) + coord(1:3,iatm) * coefgrp(k,j)
      end do
    end do


    ! Calculation of distance between two group
    !     U(r) = const(1) * (r - Ref(1))^k : r > Ref(1)
    !     U(r) = 0                               : Ref(1) >= r
    !
    !           ndist
    !   efb    = sum U(r_n)
    !             n
    !
    !     r_n  = |rcom_i - rcom_j|
    !
    !            N_i
    !   rcom_i = sum w_l * r_l ; if w_l = mass_l/totmass: mass weighted,
    !             l in group-i ;    w_l = 1     / N_i   : average
    !
    !            N_j
    !   rcom_j = sum w_m * r_m ; if w_m = mass_m/totmass: mass weighted,
    !             m in group-j ;    w_m = 1     / N_j   : average
    !
    ! Calculation of force
    ! grad_coef(n) = k*const(1)*(r_n - Ref(1))^(k-1)
    !   force(n,i) = grad_coef(n) * r_n^(-1) * dij(n) * sig(i) * coefgrp(i)   
    !

    idist = 0
    do id1 = 1, numgrp
      do id2 = id1+1, numgrp
        idist = idist + 1
        dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
        if (boundary%type == BoundaryTypePBC) then
            dij(1:3,idist) = dij(1:3,idist) - box_size(1:3)* &
                             anint(dij(1:3,idist)/box_size(1:3))
        end if
        rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) +  &
                              dij(2,idist)*dij(2,idist) +  &
                              dij(3,idist)*dij(3,idist))
      end do
    end do

    efb = 0.0_wp
    do idist = 1, ndist
      if (rij(idist) > ref(1)) then
        r_dif = rij(idist) - ref(1)
        if (iexpo == 2) then
          efb              = efb + const(1)*r_dif*r_dif
          grad_coef(idist) = 2.0_wp*const(1)*r_dif
        else
          efb              = efb + const(1)*r_dif**(iexpo)
          grad_coef(idist) = real(iexpo,wp)*const(1)*r_dif**(iexpo-1)
        end if
      else
        grad_coef(idist) = 0.0_wp
      end if
    end do

    if (calc_force) then
      idist = 0
      do id1 = 1, numgrp
        do id2 = id1+1, numgrp
          idist = idist + 1

          if (rij(idist) < 1e-10_wp) then
            coefi = 0.0_wp
          else
            coefi = grad_coef(idist) * 1.0_wp / rij(idist)
          end if

          do j = 1, numatoms(grplist(id1,inum))
            iatm = atmlist(j,grplist(id1,inum))
            work(1:3) = - coefi * coefgrp(j,id1) * dij(1:3,idist)
            force(1:3,iatm) = force(1:3,iatm) + work(1:3)
          end do

          do j = 1, numatoms(grplist(id2,inum))
            iatm = atmlist(j,grplist(id2,inum))
            work(1:3) =   coefi * coefgrp(j,id2) * dij(1:3,idist)
            force(1:3,iatm) = force(1:3,iatm) + work(1:3)
          end do

          ! virial
          ! 
          work(1:3) = coefi * dij(1:3,idist)
          do j = 1,3
            do k = j+1, 3
              vtmp = -work(k) * dij(j,idist)
              virial(k,j) = virial(k,j) + vtmp
              virial(j,k) = virial(j,k) + vtmp
            end do
            vtmp = -work(j) * dij(j,idist)
            virial(j,j) = virial(j,j) + vtmp
          end do

        end do
      end do
    end if

    deallocate(grad_coef)

    return

  end subroutine compute_energy_flatbottom_restraint

end module at_energy_restraints_mod
