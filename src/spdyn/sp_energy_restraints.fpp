!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_restraints_mod
!> @brief   calculate restraints energy
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_restraints_mod

  use sp_rpath_str_mod
  use sp_enefunc_fit_mod
  use sp_restraints_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_experiments_mod
  use sp_boundary_str_mod
  use fitting_mod
  use fitting_str_mod
  use dihedral_libs_mod
  use mpi_parallel_mod
  use timers_mod
  use messages_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! subroutines
  public  :: compute_energy_restraints
  private :: compute_energy_restraints_pos
  private :: compute_energy_restraints_rmsd
  private :: compute_energy_restraints_dist
  private :: compute_energy_restraints_angle
  private :: compute_energy_restraints_dihed
  private :: compute_energy_restraints_repul
  private :: compute_energy_restraints_fb
  private :: get_restraint_coordinates
  private :: get_restraint_forces

  ! FEP
  public  :: compute_energy_restraints_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints
  !> @brief        calculate restraint energy
  !! @authors      JJ, TM
  !! @param[in]    get_coord  : flag for whether to get coordinates
  !! @param[in]    calc_force : flag for whether to calculate forces
  !! @param[in]    domain     : domain information
  !! @param[in]    boundary   : boundary information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] virial_ext : virial term of target systems
  !! @param[inout] eposi      : point restraint energy of target systems
  !! @param[inout] ermsd      : rmsd restraint energy of target systems
  !! @param[inout] rmsd       : rmsd
  !! @param[inout] ebonds     : 
  !! @param[inout] eemfit     :
  !! @param[inout] emcorr     :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints(get_coord, calc_force, domain, boundary,&
                                       enefunc, coord, force, virial,          &
                                       virial_ext, eposi, ermsd, rmsd, ebonds, &
                                       eemfit, emcorr)

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
    real(dp),                intent(inout) :: ermsd
    real(dp),                intent(inout) :: rmsd
    real(dp),                intent(inout) :: ebonds
    real(dp),                intent(inout) :: eemfit
    real(dp),                intent(inout) :: emcorr

    ! local variables
    real(dp)                 :: edist, eangle, edihed, etilt, eexp, erepul, efb
    real(wp),        pointer :: bonds_coord(:,:), bonds_force(:,:)
    integer,         pointer :: bonds_to_atom(:)


    call timer(TimerRestraint, TimerOn)

    bonds_to_atom => enefunc%restraint_bondslist_to_atomlist
    bonds_coord   => enefunc%restraint_bonds_coord
    bonds_force   => enefunc%restraint_bonds_force

    edist  = 0.0_dp
    eangle = 0.0_dp
    edihed = 0.0_dp
    ermsd  = 0.0_dp
    eemfit = 0.0_dp
    erepul = 0.0_dp
    efb    = 0.0_dp
    bonds_force(1:3,1:enefunc%num_atoms_bonds_restraint) = 0.0_wp

    if (get_coord) then
      call get_restraint_coordinates(domain, enefunc, coord, bonds_to_atom, &
                                     bonds_coord)
    end if

    ! point restraint energy
    !
    if (enefunc%restraint_posi) &
    call compute_energy_restraints_pos(calc_force, domain, enefunc, coord,  &
                                       force, virial, virial_ext, eposi)

    if (enefunc%restraint_rmsd) &
    call compute_energy_restraints_rmsd(domain, enefunc, coord,             &
                                       force, virial, virial_ext, ermsd, rmsd)

    call compute_energy_restraints_dist(calc_force, enefunc, bonds_coord,   &
                                       bonds_force, virial, edist)

    call compute_energy_restraints_angle(calc_force, enefunc, bonds_coord,  &
                                       bonds_force, virial, eangle)

    call compute_energy_restraints_dihed(calc_force, enefunc, bonds_coord,  &
                                       bonds_force, virial, edihed)

    call compute_energy_experimental_restraint(calc_force, domain, boundary,&
                                       enefunc, coord, force, virial_ext,   &
                                       eemfit, emcorr)

    call compute_energy_restraints_repul(calc_force, enefunc, boundary, &
                                    bonds_coord, bonds_force, virial, erepul)

    call compute_energy_restraints_fb(calc_force, enefunc, boundary, &
                                    bonds_coord, bonds_force, virial, efb)

    ebonds = ebonds + edist + eangle + edihed + erepul + efb 
 
    if (calc_force) then
      call get_restraint_forces(domain, enefunc, bonds_to_atom,             &
                                bonds_force, force)
    end if

    call timer(TimerRestraint, TimerOff)

    return

  end subroutine compute_energy_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_pos
  !> @brief        calculate position restraint energy
  !! @authors      JJ
  !! @param[in]    calc_force : flag for whether to calculate forces
  !! @param[in]    domain     : domain information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] virial_ext : virial term of target systems
  !! @param[out]   erest      : restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_pos(calc_force, domain, enefunc, coord, &
                                           force, virial, virial_ext, erest)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: virial_ext(3,3,nthread)
    real(dp),                intent(inout) :: erest(nthread)

    ! local variables
    real(wp)                 :: d(1:3), v(1:3, 1:3), work(1:3), f(1:3) 
    real(wp)                 :: enemax, erest_temp, grad_coef
    real(wp)                 :: rtmp, etmp, wtmp
    real(wp)                 :: tmp_wp(1:3)
    real(wip)                :: d2(1:3)
    integer                  :: i, ix, i1, j, k
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(wp),        pointer :: restraint_force(:,:,:), restraint_coord(:,:,:)
    real(wp),        pointer :: rotated_coord(:,:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)
    integer,         pointer :: id_atm2cv(:)


    id_g2l          => domain%id_g2l

    nrestraint      => enefunc%num_restraint
    restraint_atom  => enefunc%restraint_atom
    restraint_force => enefunc%restraint_force
    restraint_coord => enefunc%restraint_coord
    rotated_coord   => enefunc%rotated_coord
    id_atm2cv       => enefunc%stats_id_atm2cv

    ncell_local = domain%num_cell_local

    enemax = 1.0e+4_wp

    select case (enefunc%fitting_method)

    case (FittingMethodTR_ROT)

      rotated_coord   (1:3, 1:MaxRest, 1:ncell_local) = 0.0_wp
      call fitting_sel(domain, enefunc, coord, rotated_coord)
      
    case (FittingMethodXYTR_ZROT)

      rotated_coord   (1:3, 1:MaxRest, 1:ncell_local) = 0.0_wp
      call fitting_sel_2d(domain, enefunc, coord, rotated_coord)

    case (FittingMethodNO)

      ! nothing to do
    end select

    ! calculate restraint energy
    !
    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, i1, d, d2,  grad_coef, rtmp, v, etmp, wtmp, work, &
    !$omp         erest_temp, j, k, f, tmp_wp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell_local, nthread

      erest_temp = 0.0_wp
      v(1:3,1:3) = 0.0_wp

      do ix = 1, nrestraint(i)

        i1 = id_g2l(2,restraint_atom(ix,i))

        ! restraint energy: E=K[b-b0]^2
        !
        d(1:3)    = coord(1:3,i1,i) - restraint_coord(1:3,ix,i)
        f(1:3)    = restraint_force(2:4,ix,i)

        rtmp      = f(1)*d(1)*d(1) + f(2)*d(2)*d(2) + f(3)*d(3)*d(3)
        etmp      = restraint_force(1,ix,i)*rtmp
        grad_coef = restraint_force(1,ix,i)*2.0_wp
        work(1:3) = grad_coef*f(1:3)*d(1:3)

        d2(1:3) = real(d(1:3),wip)
        if ((enefunc%fitting_method == FittingMethodTR_ROT) .or. &
            (enefunc%fitting_method == FittingMethodXYTR_ZROT)) then
          tmp_wp(1:3)  = rotated_coord(1:3,ix,i) - restraint_coord(1:3,ix,i)
          d2(1:3)      = real(tmp_wp(1:3),wip)
        end if

        ! check too big energy
        !
        if (etmp > enemax) call error_msg('Posicon_ene> energy is too big')
        erest_temp = erest_temp + etmp

        if (calc_force) then
          ! force
          !
          force(1:3,i1,i,id+1) = force(1:3,i1,i,id+1) - work(1:3)

          ! virial
          !
          do j = 1, 3
            do k = j+1, 3
              wtmp = -coord(k,i1,i)*work(j)
              v(k,j) = v(k,j) + wtmp
              v(j,k) = v(j,k) + wtmp
            end do
            wtmp = -coord(j,i1,i)*work(j)
            v(j,j) = v(j,j) + wtmp
          end do
        end if

        ! rpath
        !
        if (enefunc%rpath_sum_mf_flag .and. enefunc%rpath_pos_func > 0) then
          enefunc%stats_delta(3*id_atm2cv(restraint_atom(ix,i)) - 2) = real(d2(1),wp)
          enefunc%stats_delta(3*id_atm2cv(restraint_atom(ix,i)) - 1) = real(d2(2),wp)
          enefunc%stats_delta(3*id_atm2cv(restraint_atom(ix,i))    ) = real(d2(3),wp)
        end if

      end do
      if (calc_force .and. .not. enefunc%pressure_position) &
        virial_ext(1:3,1:3,id+1) = virial_ext(1:3,1:3,id+1) + v(1:3,1:3)
      if (calc_force .and. enefunc%pressure_position) &
        virial(1:3,1:3,id+1) = virial(1:3,1:3,id+1) + v(1:3,1:3)
      erest(id+1) = erest(id+1) + erest_temp

    end do

    !$omp end parallel

    ! rpath
    !
    if (enefunc%rpath_sum_mf_flag .and. enefunc%rpath_pos_func > 0) then
      enefunc%stats_count = enefunc%stats_count + 1
      do i = 1, enefunc%stats_dimension/3
        enefunc%stats_grad (1,i,3*i-2) = 1.0_dp
        enefunc%stats_grad (2,i,3*i-1) = 1.0_dp
        enefunc%stats_grad (3,i,3*i  ) = 1.0_dp
      end do
    end if

    return

  end subroutine compute_energy_restraints_pos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_rmsd
  !> @brief        calculate position restraint energy
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[out]   force       : forces of target systems
  !! @param[out]   virial      : virial term of target systems
  !! @param[out]   virial_ext  : external virial term of target systems
  !! @param[out]   ermsd       : rmsd restraint energy
  !! @param[out]   rmsd        : rmsd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_rmsd(domain, enefunc, coord, force, &
                                            virial, virial_ext, ermsd, rmsd)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: virial_ext(3,3,nthread)
    real(dp),                intent(inout) :: ermsd
    real(dp),                intent(inout) :: rmsd

    ! local variables
    real(dp)                 :: com_ref(1:3), com_mov(1:3), tot_mass
    real(dp)                 :: sym_matrix(1:4,1:4), rot_matrix(1:3,1:3)
    real(dp)                 :: eval(1:4), evec(1:4), work(1:11)
    real(dp)                 :: dref(1:3), dmov(1:3), dsub(1:3), dadd(1:3)
    real(dp)                 :: values(1:2)
    real(dp)                 :: Krmsd, rmsd0, factor

    real(wp)                 :: viri_local(1:3)
    integer                  :: i, ix, ig, i1
    integer                  :: Natom
    integer                  :: ncell_local
    integer                  :: ierr
    logical                  :: withmass

    real(wp),        pointer :: ref_coord(:,:,:)
    real(wp),        pointer :: rotated_coord(:,:,:)
    real(wip),       pointer :: mass(:,:)
    integer(int2),   pointer :: id_g2l(:,:)
    integer,         pointer :: nrestraint(:), restraint_atom(:,:)


    id_g2l          => domain%id_g2l
    mass            => domain%mass

    nrestraint      => enefunc%num_restraint
    restraint_atom  => enefunc%restraint_atom
    ref_coord       => enefunc%restraint_coord
    rotated_coord   => enefunc%rotated_coord

    Krmsd           =  enefunc%rmsd_force
    rmsd0           =  enefunc%rmsd_target
    withmass        =  enefunc%rmsd_withmass

    ncell_local     = domain%num_cell_local
    Natom           = dble(enefunc%nrmsd)

    select case (enefunc%fitting_method)

    case (FittingMethodTR_ROT)

      rotated_coord   (1:3, 1:MaxRest, 1:ncell_local) = 0.0_wp
      call fitting_sel(domain, enefunc, coord, rotated_coord)
      
    case (FittingMethodXYTR_ZROT)

      rotated_coord   (1:3, 1:MaxRest, 1:ncell_local) = 0.0_wp
      call fitting_sel_2d(domain, enefunc, coord, rotated_coord)

    case (FittingMethodNO)

      rotated_coord => enefunc%restraint_coord

      ! nothing to do
    end select

    ! rmsd
    !
    rmsd = 0.0_dp 
    tot_mass = 0.0_dp 

    if (withmass) then
      do i = 1, ncell_local
        do ix = 1, nrestraint(i)
          ig = restraint_atom(ix,i)
          i1 = id_g2l(2,ig)
          dsub(1:3) = coord(1:3,i1,i) - rotated_coord(1:3,ix,i)
          if (enefunc%fitting_method == FittingMethodXYTR_ZROT) &
            dsub(3)=0.0_dp 

          rmsd = rmsd + mass(i1,i) &
                       *(dsub(1)*dsub(1)+dsub(2)*dsub(2)+dsub(3)*dsub(3))
          tot_mass = tot_mass + mass(i1,i)
        end do
      end do
      values(1) = rmsd
      values(2) = tot_mass
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, values, 2, mpi_real8, mpi_sum, &
                         mpi_comm_city, ierror)
#endif
      rmsd     = values(1)
      tot_mass = values(2)
      rmsd = sqrt(rmsd/tot_mass)
    else
      do i = 1, ncell_local
        do ix = 1, nrestraint(i)
          ig = restraint_atom(ix,i)
          i1 = id_g2l(2,ig)
          dsub(1:3) = coord(1:3,i1,i) - rotated_coord(1:3,ix,i)
          if (enefunc%fitting_method == FittingMethodXYTR_ZROT) &
            dsub(3)=0.0_dp 
          rmsd = rmsd + dsub(1)*dsub(1)+dsub(2)*dsub(2)+dsub(3)*dsub(3)
        end do
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, rmsd, 1, mpi_real8, mpi_sum, &
                         mpi_comm_city, ierror)
#endif
      rmsd = sqrt(rmsd/real(Natom,dp))
    end if

    if (enefunc%restraint_rmsd_target) then

      ! targeted md case
      !
      ermsd = 0.0_dp
      viri_local(1:3) = 0.0_dp

      if (withmass) then

        do i = 1, ncell_local
          do ix = 1, nrestraint(i)
            ig = restraint_atom(ix,i)
            i1 = id_g2l(2,ig)
            dsub(1:3) = coord(1:3,i1,i)-rotated_coord(1:3,ix,i)
            if (enefunc%fitting_method == FittingMethodXYTR_ZROT) &
              dsub(3)=0.0_dp
            factor = Krmsd*mass(i1,i)*(1.0_wp-rmsd0/rmsd)
            ermsd  = ermsd + factor*dsub(1)*dsub(1) &
                           + factor*dsub(2)*dsub(2) &
                           + factor*dsub(3)*dsub(3)
            factor = -2.0_wp*factor
            force(1:3,i1,i,1) = force(1:3,i1,i,1) + factor*dsub(1:3)
            viri_local(1:3) = viri_local(1:3) + coord(1:3,ix,i)*factor*dsub(1:3)
          end do
        end do

      else

        do i = 1, ncell_local
          do ix = 1, nrestraint(i)
            ig = restraint_atom(ix,i)
            i1 = id_g2l(2,ig)
            dsub(1:3) = coord(1:3,i1,i)-rotated_coord(1:3,ix,i)
            if (enefunc%fitting_method == FittingMethodXYTR_ZROT) &
              dsub(3)=0.0_dp
            factor = Krmsd*(1.0_wp-rmsd0/rmsd)
            ermsd  = ermsd + factor*dsub(1)*dsub(1) &
                           + factor*dsub(2)*dsub(2) &
                           + factor*dsub(3)*dsub(3)
            factor = -2.0_wp*factor
            force(1:3,i1,i,1) = force(1:3,i1,i,1) + factor*dsub(1:3)
            viri_local(1:3) = viri_local(1:3) + coord(1:3,ix,i)*factor*dsub(1:3)
          end do
        end do

      end if
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, ermsd, 1, mpi_real8, mpi_sum, &
                         mpi_comm_city, ierror)
#endif

    else

      ! steered MD or just rmsd restraint
      !
      ermsd = Krmsd * (rmsd - rmsd0)**2
      viri_local(1:3) = 0.0_dp 
      factor = -2.0_dp*Krmsd*(rmsd-rmsd0) / (Natom*rmsd) 
      do i = 1, ncell_local
        do ix = 1, nrestraint(i)

          ig = restraint_atom(ix,i)
          i1 = id_g2l(2,ig)
          dsub(1:3) = coord(1:3,i1,i)-rotated_coord(1:3,ix,i)
          if (enefunc%fitting_method == FittingMethodXYTR_ZROT) &
            dsub(3)=0.0_dp 
          force(1:3,i1,i,1) = force(1:3,i1,i,1) + factor*dsub(1:3)
          viri_local(1:3) = viri_local(1:3) + coord(1:3,ix,i)*factor*dsub(1:3)
        end do
      end do

    end if

    if (.not. enefunc%pressure_rmsd) then
      virial_ext(1,1,1) = virial_ext(1,1,1) + viri_local(1)
      virial_ext(2,2,1) = virial_ext(2,2,1) + viri_local(2)
      virial_ext(3,3,1) = virial_ext(3,3,1) + viri_local(3)
    else
      virial(1,1,1) = virial(1,1,1) + viri_local(1)
      virial(2,2,1) = virial(2,2,1) + viri_local(2)
      virial(3,3,1) = virial(3,3,1) + viri_local(3)
    end if

    return

  end subroutine compute_energy_restraints_rmsd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_dist
  !> @brief        calculate distance restraint energy
  !! @authors      JJ, CK
  !! @param[in]    calc_force : flag for whether to calculate forces
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    dist_coord : restraint bonds coordinates
  !! @param[inout] dist_force : restraint bonds forces
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] edist      : distance restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_dist(calc_force, enefunc, dist_coord, &
                                            dist_force, virial, edist)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: dist_coord(:,:)
    real(wp),                intent(inout) :: dist_force(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: edist

    ! local variable
    real(wp)                 :: sig(1:2)
    real(wp)                 :: work(1:3)
    real(wp)                 :: r_sum, dtmp, rtmp, r_dif, coef, coefi, grad_coef
    real(wp)                 :: tmp_wp
    integer                  :: numgrp, iexpo, idist, ndist
    integer                  :: inum, i, j, jx, index, id, id1, id2, ind
    integer                  :: iatm
    integer                  :: i_rpath_dim

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: wtind(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)
    integer,         pointer :: expoind(:,:)


    num_funcs => enefunc%num_restraintfuncs
    funcgrp   => enefunc%restraint_funcgrp
    const     => enefunc%restraint_const
    ref       => enefunc%restraint_ref
    expo      => enefunc%restraint_exponent_func
    grplist   => enefunc%restraint_grouplist
    numatoms  => enefunc%restraint_numatoms
    bondslist => enefunc%restraint_bondslist
    expoind   => enefunc%restraint_exponent_dist
    wtind     => enefunc%restraint_weight_dist
    kind      => enefunc%restraint_kind
    masscoef  => enefunc%restraint_masscoef
    coefgrp   => enefunc%restraint_wtmp
    com       => enefunc%restraint_wcom1
    dij       => enefunc%restraint_wcom2
    rij       => enefunc%restraint_wdrt

    do inum = 1, num_funcs

      if (kind(inum) == RestraintsFuncDIST .or. &
          kind(inum) == RestraintsFuncDISTCOM) then

        ! initialization
        !
        numgrp = funcgrp(inum)
        iexpo  = expo(inum)
        ndist  = int(numgrp/2)

        if (kind(inum) == RestraintsFuncDIST) then
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                    1.0_wp/real(numatoms(grplist(j,inum)),wp)
          end do
        else
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                    masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
          end do
        end if

        ! check atom group range and calculate center of mass of each group
        !
        do j = 1, numgrp
          com(1:3,j) = 0.0_wp
          do jx = 1, numatoms(grplist(j,inum))
            index = bondslist(jx,grplist(j,inum))
            com(1:3,j) = com(1:3,j) + dist_coord(1:3,index) * coefgrp(jx,j)
          end do
        end do

        ! calculation of distance energy
        !
        r_sum = 0.0_wp

        do idist = 1, ndist

          id1 = idist*2 - 1
          id2 = idist*2
          dtmp = 0.0_wp

          dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
          rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) + &
                                dij(2,idist)*dij(2,idist) + &
                                dij(3,idist)*dij(3,idist) )
          if (expoind(idist,inum) == 1) then
            r_sum = r_sum + wtind(idist,inum)*rij(idist)
          else
            r_sum = r_sum + wtind(idist,inum)*rij(idist)**(expoind(idist,inum))
          end if

        end do

        if (r_sum < ref(1,inum)) then
          ind = 1
        else if (r_sum > ref(2,inum)) then
          ind = 2
        else
          cycle
        end if

        r_dif = r_sum -ref(ind,inum)

        if (iexpo == 2) then
          coef       = const(ind,inum)
          edist      = edist + coef*r_dif*r_dif
          grad_coef  = 2.0_wp * const(ind,inum)*r_dif
        else
          coef       = const(ind,inum)
          edist      = edist + coef*r_dif**(iexpo)
          grad_coef  = real(iexpo,wp) * const(ind,inum)*r_dif**(iexpo-1)
        end if

        sig(1) = -1.0_wp
        sig(2) = 1.0_wp

        ! calculation of force
        !
        if (calc_force) then
          do idist = 1,ndist

            if (rij(idist) < 1e-10_wp) then
              coefi = 0.0_wp
            else
              if (expoind(idist,inum) == 1) then
                coefi = grad_coef * wtind(idist,inum) / rij(idist)
              else
                coefi = expoind(idist,inum) * grad_coef &
                       * rij(idist)**(expoind(idist,inum)-2)
              end if
            end if

            do i = 1, 2
              id = idist*2-2+i
              do j = 1,numatoms(grplist(id,inum))
                iatm = bondslist(j,grplist(id,inum))
                work(1:3) = sig(i) * coefi  * coefgrp(j,id) * dij(1:3,idist)
                dist_force(1:3,iatm) = dist_force(1:3,iatm) + work(1:3)
              end do
            end do

            ! calculation of virial
            !
            if (main_rank) then
              do i = 1, 3
                work(i)= coefi * dij(i,idist)
                do j = 1, 3
                  virial(j,i,1) = virial(j,i,1) - work(i) * dij(j,idist)
                end do
              end do
            end if
          end do
        end if

        ! rpath
        !
        if (enefunc%rpath_sum_mf_flag) then
          i_rpath_dim=enefunc%restraint_rpath_func(inum)
          if (i_rpath_dim >= 1) then
          if (i_rpath_dim == 1) enefunc%stats_count = enefunc%stats_count + 1
            enefunc%stats_delta(i_rpath_dim) = r_dif
            tmp_wp = dij(1, 1) / rij(1)
            enefunc%stats_grad(1,1,i_rpath_dim) = tmp_wp
            tmp_wp = dij(2, 1) / rij(1)
            enefunc%stats_grad(2,1,i_rpath_dim) = tmp_wp
            tmp_wp = dij(3, 1) / rij(1)
            enefunc%stats_grad(3,1,i_rpath_dim) = tmp_wp
            tmp_wp = -dij(1, 1) / rij(1)
            enefunc%stats_grad(1,2,i_rpath_dim) = tmp_wp
            tmp_wp = -dij(2, 1) / rij(1)
            enefunc%stats_grad(2,2,i_rpath_dim) = tmp_wp
            tmp_wp = -dij(3, 1) / rij(1)
            enefunc%stats_grad(3,2,i_rpath_dim) = tmp_wp
          end if
        end if

      end if

    end do

    return

  end subroutine compute_energy_restraints_dist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_angle
  !> @brief        calculate angle restraint energy
  !! @authors      CK
  !! @param[in]    calc_force  : flag for whether to calculate forces
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    angle_coord : restraint bonds coordinates
  !! @param[inout] angle_force : restraint bonds forces
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] erestangle  : angle restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_angle(calc_force, enefunc, angle_coord, &
                                             angle_force, virial, erestangle)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: angle_coord(:,:)
    real(wp),                intent(inout) :: angle_force(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: erestangle

    ! local variable
    real(wp)                 :: dij(1:3), dkj(1:3)
    real(wp)                 :: inv_rij2, inv_rkj2
    real(wp)                 :: rij2, rkj2, rijrkj,inv_rijrkj, vtmp
    real(wp)                 :: cos_t, t_dif, coef, theta, sin_t, grad_coef
    real(wp)                 :: work(1:6)
    real(wp)                 :: radref(1:2)
    integer                  :: numgrp, iexpo
    integer                  :: inum, j, k, jx, index, ind
    integer                  :: iatm
    logical                  :: anglecheck

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    num_funcs  => enefunc%num_restraintfuncs
    funcgrp    => enefunc%restraint_funcgrp
    const      => enefunc%restraint_const
    ref        => enefunc%restraint_ref
    expo       => enefunc%restraint_exponent_func
    grplist    => enefunc%restraint_grouplist
    numatoms   => enefunc%restraint_numatoms
    bondslist  => enefunc%restraint_bondslist
    kind       => enefunc%restraint_kind
    masscoef   => enefunc%restraint_masscoef
    coefgrp    => enefunc%restraint_wtmp
    com        => enefunc%restraint_wcom1


    do inum = 1, num_funcs

      if (kind(inum) == RestraintsFuncANGLE .or. &
          kind(inum) == RestraintsFuncANGLECOM) then

        ! initialization
        !
        numgrp = funcgrp(inum)
        iexpo  = expo(inum)

        if (kind(inum) == RestraintsFuncANGLE) then
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                    1.0_wp/real(numatoms(grplist(j,inum)),wp)
          end do
        else
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                    masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
          end do
        end if

        ! check atom group range and calculate center of mass of each group
        !
        do j = 1, numgrp
          com(1:3,j) = 0.0_wp
          do jx = 1, numatoms(grplist(j,inum))
            index = bondslist(jx,grplist(j,inum))
            com(1:3,j) = com(1:3,j) + angle_coord(1:3,index) * coefgrp(jx,j)
          end do
        end do

        ! calculation of angle energy
        !
        radref(1:2) = ref(1:2,inum)*RAD
     
        dij(1:3) = com(1:3,1) - com(1:3,2)
        dkj(1:3) = com(1:3,3) - com(1:3,2)
        rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        rkj2  = dkj(1)*dkj(1) + dkj(2)*dkj(2) + dkj(3)*dkj(3)
        rijrkj = sqrt(rij2*rkj2)

        inv_rijrkj = 1.0_wp /rijrkj
        inv_rij2 = 1.0_wp /rij2
        inv_rkj2 = 1.0_wp /rkj2

        cos_t = (dij(1)*dkj(1) + dij(2)*dkj(2) + dij(3)*dkj(3) )*inv_rijrkj
        cos_t  = min( 1.0_wp, cos_t)
        cos_t  = max(-1.0_wp, cos_t)
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

        anglecheck = .true.
        do while(anglecheck)
          if (t_dif > pi)  t_dif = t_dif - 2.0_wp * pi
          if (t_dif < -pi) t_dif = t_dif + 2.0_wp * pi
          if (t_dif >= -pi  .and. t_dif <= pi)  anglecheck = .false.
        end do

        if (iexpo == 2) then
          coef       = const(ind,inum) 
          erestangle = erestangle + coef * t_dif * t_dif
          grad_coef  = 2.0_wp * const(ind,inum) * t_dif
        else
          coef       = const(ind,inum) 
          erestangle = erestangle + coef * t_dif**(iexpo)
          grad_coef  = real(iexpo,wp) * const(ind,inum) * t_dif**(iexpo-1)
        end if
        
        sin_t = (1.0_wp-cos_t*cos_t)
        sin_t = sqrt(sin_t)
        sin_t  = max(EPS, sin_t)
        grad_coef = -grad_coef / sin_t
  
        work(1:3) = grad_coef * (dkj(1:3) * inv_rijrkj -  &
                                 dij(1:3) * cos_t * inv_rij2)

        work(4:6) = grad_coef * (dij(1:3) * inv_rijrkj -  &
                                 dkj(1:3) * cos_t * inv_rkj2)
     
        ! compute force
        !
        if (calc_force)  then    
          do j = 1,numatoms(grplist(1,inum))
            iatm = bondslist(j,grplist(1,inum))
            angle_force(1:3,iatm) = angle_force(1:3,iatm)                     &
                                   - work(1:3)*coefgrp(j,1)
          end do
     
          do j = 1,numatoms(grplist(2,inum))
            iatm = bondslist(j,grplist(2,inum))
            angle_force(1:3,iatm) = angle_force(1:3,iatm) +                   &
                                  (work(1:3) + work(4:6))*coefgrp(j,2)
          end do
     
          do j = 1,numatoms(grplist(3,inum))
            iatm = bondslist(j,grplist(3,inum))
            angle_force(1:3,iatm) = angle_force(1:3,iatm)                     &
                                   - work(4:6)*coefgrp(j,3)
          end do

          ! calculation of virial
          !
          if (main_rank) then
            do j = 1, 3
         
              do k = j+1, 3
                vtmp = -(dij(k)*work(j) + dkj(k)*work(j+3))
                virial(k,j,1) = virial(k,j,1) + vtmp
                virial(j,k,1) = virial(j,k,1) + vtmp
              end do
         
              vtmp = -(dij(j)*work(j) + dkj(j)*work(j+3))
              virial(j,j,1) = virial(j,j,1) + vtmp
            end do
          end if
        end if
      end if
    end do

    return

  end subroutine compute_energy_restraints_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_dihed
  !> @brief        calculate dihedral restraint energy
  !! @authors      CK
  !! @param[in]    calc_force  : flag for whether to calculate forces
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    dihed_coord : restraint bonds coordinates
  !! @param[inout] dihed_force : restraint bonds forces
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] erestdihed  : dihedral restraint energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_dihed(calc_force, enefunc, dihed_coord, &
                                             dihed_force, virial, erestdihed)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: dihed_coord(:,:)
    real(wp),                intent(inout) :: dihed_force(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: erestdihed

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, coef
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi
    real(wp)                 :: cosdif, sindif, theta, vtmp
    real(wp)                 :: work(1:9)
    real(wp)                 :: radref(1:2)
    real(wp)                 :: tmp_wp
    integer                  :: aindex(1:4)
    integer                  :: numgrp, iexpo, inum
    integer                  :: i, j, k, iatm, ind
    integer                  :: i_rpath_dim
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:)
    real(wp),        pointer :: coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    num_funcs  => enefunc%num_restraintfuncs
    funcgrp    => enefunc%restraint_funcgrp 
    expo       => enefunc%restraint_exponent_func
    grplist    => enefunc%restraint_grouplist
    numatoms   => enefunc%restraint_numatoms
    bondslist  => enefunc%restraint_bondslist
    kind       => enefunc%restraint_kind
    masscoef   => enefunc%restraint_masscoef
    coefgrp    => enefunc%restraint_wtmp
    com        => enefunc%restraint_wcom1
    const      => enefunc%restraint_const
    ref        => enefunc%restraint_ref

    do inum = 1, num_funcs

      if (kind(inum) == RestraintsFuncDIHED .or. &
          kind(inum) == RestraintsFuncDIHEDCOM) then
        coefgrp(1:enefunc%max_restraint_numatoms,1:enefunc%num_restraintgroups) &
          = 0.0_wp
        numgrp  = funcgrp(inum)
        iexpo   = expo   (inum)


        if (kind(inum) == RestraintsFuncDIHED) then
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                1.0_wp/real(numatoms(grplist(j,inum)),wp)
          end do
        else
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
                    masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
          end do
        end if
     
        ! calculation of com
        !
        
        com(1:3,1:4)=0.0_wp
     
        do i = 1,4 
          do j = 1,numatoms(grplist(i,inum))
            iatm = bondslist(j,grplist(i,inum))
            com(1:3,i) = com(1:3,i) + dihed_coord(1:3,iatm) * coefgrp(j,i)
          end do
          aindex(i) = i
        end do

        radref(1:2) = ref(1:2,inum)*RAD
        call calculate_dihedral(aindex, com, cos_dih, sin_dih, grad, v)

        if (cos_dih > 1.0E-1_wp) then
          theta = asin(sin_dih)
        else
          theta = sign(1.0_wp,sin_dih)*acos(cos_dih)
        end if
       
        if (theta < radref(1)) then
          ind = 1
        else if (theta > radref(2)) then
          ind = 2
        else 
          cycle
!          return
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
          coef       = const(ind,inum) 
          erestdihed = erestdihed + coef * diffphi * diffphi
          grad_coef  = 2.0_wp * const(ind,inum) * diffphi
        else
          coef       = const(ind,inum)
          erestdihed = erestdihed + coef * diffphi**(iexpo)
          grad_coef  = real(iexpo,wp) * const(ind,inum) * diffphi**(iexpo-1)
        end if


        if (calc_force) then       
          ! compute virial
          !
          if (main_rank) then
            do j = 1, 3
              do k = j+1, 3
                vtmp = grad_coef*v(k,j)
         
                virial(k,j,1) = virial(k,j,1) + vtmp
                virial(j,k,1) = virial(j,k,1) + vtmp
              end do
         
              vtmp = grad_coef*v(j,j)
              virial(j,j,1) = virial(j,j,1) + vtmp
            end do
          end if
          if (abs(theta) < 1e-10_wp) then
            grad_coef = 0.0_wp
          end if
          work(1:9) = grad_coef*grad(1:9)

          ! compute force
          !
          do j = 1,numatoms(grplist(1,inum))
            iatm = bondslist(j,grplist(1,inum))
            dihed_force(1:3,iatm) = dihed_force(1:3,iatm)             &
                                   - work(1:3)*coefgrp(j,1)
          end do
       
          do j = 1,numatoms(grplist(2,inum))
            iatm = bondslist(j,grplist(2,inum))
            dihed_force(1:3,iatm) = dihed_force(1:3,iatm) +           &
                                  (work(1:3) - work(4:6))*coefgrp(j,2)
          end do
       
          do j = 1,numatoms(grplist(3,inum))
            iatm = bondslist(j,grplist(3,inum))
            dihed_force(1:3,iatm) = dihed_force(1:3,iatm) +           &
                                  (work(4:6) + work(7:9))*coefgrp(j,3)
          end do
       
          do j = 1,numatoms(grplist(4,inum))
            iatm = bondslist(j,grplist(4,inum))
            dihed_force(1:3,iatm) = dihed_force(1:3,iatm)             &
                                   - work(7:9)*coefgrp(j,4)
          end do
        end if

        ! rpath
        !
        if (enefunc%rpath_sum_mf_flag) then
          i_rpath_dim=enefunc%restraint_rpath_func(inum)
          if (i_rpath_dim >= 1) then
            if (i_rpath_dim == 1) enefunc%stats_count = enefunc%stats_count + 1
            theta = theta/RAD - ref(1,inum)
            if (theta > 180.0_wp) then
              theta = theta - 360.0_wp
            else if (theta < -180.0_wp) then
              theta = theta + 360.0_wp
            end if
            enefunc%stats_delta(i_rpath_dim)    = theta * RAD
            enefunc%stats_grad(1,1,i_rpath_dim) = grad(1)
            enefunc%stats_grad(2,1,i_rpath_dim) = grad(2)
            enefunc%stats_grad(3,1,i_rpath_dim) = grad(3)
            tmp_wp = - grad(1) + grad(4) 
            enefunc%stats_grad(1,2,i_rpath_dim) = tmp_wp
            tmp_wp = - grad(2) + grad(5) 
            enefunc%stats_grad(2,2,i_rpath_dim) = tmp_wp
            tmp_wp = - grad(3) + grad(6) 
            enefunc%stats_grad(3,2,i_rpath_dim) = tmp_wp
            tmp_wp = - grad(4) - grad(7) 
            enefunc%stats_grad(1,3,i_rpath_dim) = tmp_wp
            tmp_wp = - grad(5) - grad(8) 
            enefunc%stats_grad(2,3,i_rpath_dim) = tmp_wp
            tmp_wp = - grad(6) - grad(9) 
            enefunc%stats_grad(3,3,i_rpath_dim) = tmp_wp
            enefunc%stats_grad(1,4,i_rpath_dim) = grad(7) 
            enefunc%stats_grad(2,4,i_rpath_dim) = grad(8) 
            enefunc%stats_grad(3,4,i_rpath_dim) = grad(9)    
          end if
        end if

      end if
    end do

    return

  end subroutine compute_energy_restraints_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_repul
  !> @brief        calculate repulsive restraint energy
  !! @authors      HO
  !! @param[in]    calc_force  : flag for whether to calculate forces
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    repul_coord : restraint bonds coordinates
  !! @param[inout] repul_force : restraint bonds forces
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] erepul      : repulsive restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_repul(calc_force, enefunc, boundary, &
                                       repul_coord, repul_force, virial, erepul)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: repul_coord(:,:)
    real(wp),                intent(inout) :: repul_force(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: erepul

    ! local variable
    real(wp)                 :: sig(1:2)
    real(wp)                 :: work(1:3)
    real(wp)                 :: r_sum, dtmp, rtmp, r_dif, coef, coefi
    integer                  :: numgrp, iexpo, idist, ndist
    integer                  :: inum, i, j, jx, index, id, id1, id2, ind
    integer                  :: iatm
    real(wp),    allocatable :: grad_coef(:)
    real(wp)                 :: box_size(3)

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    num_funcs => enefunc%num_restraintfuncs
    funcgrp   => enefunc%restraint_funcgrp
    const     => enefunc%restraint_const
    ref       => enefunc%restraint_ref
    expo      => enefunc%restraint_exponent_func
    grplist   => enefunc%restraint_grouplist
    numatoms  => enefunc%restraint_numatoms
    bondslist => enefunc%restraint_bondslist
    kind      => enefunc%restraint_kind
    masscoef  => enefunc%restraint_masscoef
    coefgrp   => enefunc%restraint_wtmp
    com       => enefunc%restraint_rcom1
    dij       => enefunc%restraint_rcom2
    rij       => enefunc%restraint_rdrt
    box_size(1) = boundary%box_size_x
    box_size(2) = boundary%box_size_y
    box_size(3) = boundary%box_size_z

    do inum = 1, num_funcs
      if (kind(inum) == RestraintsFuncREPUL .or. &
          kind(inum) == RestraintsFuncREPULCOM) then

        ! initialization
        !
        numgrp = funcgrp(inum)
        iexpo  = expo(inum)
        ndist  = int(numgrp*(numgrp-1)/2)
        allocate(grad_coef(ndist))
        if (kind(inum) == RestraintsFuncREPUL) then
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
              1.0_wp/real(numatoms(grplist(j,inum)),wp)
          end do
        else
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
              masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
          end do
        end if

        ! check atom group range and calculate center of mass of each group
        !
        do j = 1, numgrp
          com(1:3,j) = 0.0_wp
          do jx = 1, numatoms(grplist(j,inum))
            index = bondslist(jx,grplist(j,inum))
            com(1:3,j) = com(1:3,j) + repul_coord(1:3,index) * coefgrp(jx,j)
          end do
        end do

        ! calculation of repulsive distance energy
        !
        idist = 0
        do id1 = 1, numgrp
          do id2 = id1+1, numgrp
            idist = idist + 1
            dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
            if (boundary%type == BoundaryTypePBC) then
              dij(1:3,idist) = dij(1:3,idist) - box_size(1:3) * &
                anint(dij(1:3,idist)/box_size(1:3))
            end if
            rij(idist)     = sqrt(dij(1,idist)*dij(1,idist) +  &
              dij(2,idist)*dij(2,idist) +  &
              dij(3,idist)*dij(3,idist))
          end do
        end do

        do idist = 1, ndist
          if (rij(idist) < ref(1,inum)) then
            r_dif = rij(idist) - ref(1,inum)
            if (iexpo == 2) then
              erepul           = erepul + const(1,inum)*r_dif*r_dif
              grad_coef(idist) = 2.0_wp * const(1,inum)*r_dif
            else
              erepul           = erepul + const(1,inum)*r_dif**(iexpo)
              grad_coef(idist) = real(iexpo,wp) * const(1,inum)*r_dif**(iexpo-1)
            end if
          else
            grad_coef(idist) = 0.0_wp
          end if
        end do

        ! calculation of force
        !
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
                iatm = bondslist(j,grplist(id1,inum))
                work(1:3) = - coefi * coefgrp(j,id1) * dij(1:3,idist)
                repul_force(1:3,iatm) = repul_force(1:3,iatm) + work(1:3)
              end do
              do j = 1, numatoms(grplist(id2,inum))
                iatm = bondslist(j,grplist(id2,inum))
                work(1:3) =   coefi * coefgrp(j,id2) * dij(1:3,idist)
                repul_force(1:3,iatm) = repul_force(1:3,iatm) + work(1:3)
              end do

              ! virial
              ! 
              if (main_rank) then
                do i = 1, 3
                  work(i)= coefi * dij(i,idist)
                  do j = 1, 3
                    virial(j,i,1) = virial(j,i,1) - work(i) * dij(j,idist)
                  end do
                end do
              end if
            end do
          end do
        end if
        deallocate(grad_coef)
      end if
    end do

    return

  end subroutine compute_energy_restraints_repul

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_fb
  !> @brief        calculate flat bottom restraint energy
  !! @authors      HO
  !! @param[in]    calc_force : flag for whether to calculate forces
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    boundary   : boundary condition information
  !! @param[in]    fb_coord   : restraint bonds coordinates
  !! @param[inout] fb_force   : restraint bonds forces
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] efb        : flat bottom restraint energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_fb(calc_force, enefunc, boundary, &
                                          fb_coord, fb_force, virial, efb)

    ! formal arguments
    logical,                 intent(in)    :: calc_force
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: fb_coord(:,:)
    real(wp),                intent(inout) :: fb_force(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)
    real(dp),                intent(inout) :: efb

    ! local variable
    real(wp)                 :: sig(1:2)
    real(wp)                 :: work(1:3)
    real(wp)                 :: r_sum, dtmp, rtmp, r_dif, coef, coefi
    integer                  :: numgrp, iexpo, idist, ndist
    integer                  :: inum, i, j, jx, index, id, id1, id2, ind
    integer                  :: iatm
    real(wp),    allocatable :: grad_coef(:)
    real(wp)                 :: box_size(3)

    real(wp),        pointer :: const(:,:)
    real(wp),        pointer :: ref(:,:)
    real(wp),        pointer :: masscoef(:,:)
    real(wp),        pointer :: com(:,:), dij(:,:), rij(:), coefgrp(:,:)
    integer,         pointer :: num_funcs, numatoms(:)
    integer,         pointer :: bondslist(:,:)
    integer,         pointer :: kind(:)
    integer,         pointer :: grplist(:,:)
    integer,         pointer :: funcgrp(:)
    integer,         pointer :: expo(:)


    num_funcs => enefunc%num_restraintfuncs
    funcgrp   => enefunc%restraint_funcgrp
    const     => enefunc%restraint_const
    ref       => enefunc%restraint_ref
    expo      => enefunc%restraint_exponent_func
    grplist   => enefunc%restraint_grouplist
    numatoms  => enefunc%restraint_numatoms
    bondslist => enefunc%restraint_bondslist
    kind      => enefunc%restraint_kind
    masscoef  => enefunc%restraint_masscoef
    coefgrp   => enefunc%restraint_wtmp
    com       => enefunc%restraint_rcom1
    dij       => enefunc%restraint_rcom2
    rij       => enefunc%restraint_rdrt
    box_size(1) = boundary%box_size_x
    box_size(2) = boundary%box_size_y
    box_size(3) = boundary%box_size_z

    do inum = 1, num_funcs
      if (kind(inum) == RestraintsFuncFB .or. &
          kind(inum) == RestraintsFuncFBCOM) then

        ! initialization
        !
        numgrp = funcgrp(inum)
        iexpo  = expo(inum)
        ndist  = int(numgrp*(numgrp-1)/2)
        allocate(grad_coef(ndist))
        if (kind(inum) == RestraintsFuncFB) then
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
              1.0_wp/real(numatoms(grplist(j,inum)),wp)
          end do
        else
          do j = 1, numgrp
            coefgrp(1:numatoms(grplist(j,inum)),j) = &
              masscoef(1:numatoms(grplist(j,inum)),grplist(j,inum))
          end do
        end if

        ! check atom group range and calculate center of mass of each group
        !
        do j = 1, numgrp
          com(1:3,j) = 0.0_wp
          do jx = 1, numatoms(grplist(j,inum))
            index = bondslist(jx,grplist(j,inum))
            com(1:3,j) = com(1:3,j) + fb_coord(1:3,index) * coefgrp(jx,j)
          end do
        end do

        ! calculation of distance energy
        !
        idist = 0
        do id1 = 1, numgrp
          do id2 = id1+1, numgrp
            idist = idist + 1
            dij(1:3,idist) = com(1:3,id1) - com(1:3,id2)
            if (boundary%type == BoundaryTypePBC) then
              dij(1:3,idist) = dij(1:3,idist) - box_size(1:3) * &
                anint(dij(1:3,idist)/box_size(1:3))
            end if
            rij(idist)       = sqrt(dij(1,idist)*dij(1,idist) +  &
              dij(2,idist)*dij(2,idist) +  &
              dij(3,idist)*dij(3,idist))
          end do
        end do

        do idist = 1, ndist
          if (rij(idist) > ref(1,inum)) then
            r_dif = rij(idist) - ref(1,inum)
            if (iexpo == 2) then
              efb              = efb + const(1,inum)*r_dif*r_dif
              grad_coef(idist) = 2.0_wp * const(1,inum)*r_dif
            else
              efb              = efb + const(1,inum)*r_dif**(iexpo)
              grad_coef(idist) = real(iexpo,wp)*const(1,inum)*r_dif**(iexpo-1)
            end if
          else
            grad_coef(idist) = 0.0_wp
          end if
        end do

        ! calculation of force
        !
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
                iatm = bondslist(j,grplist(id1,inum))
                work(1:3) = - coefi * coefgrp(j,id1) * dij(1:3,idist)
                fb_force(1:3,iatm) = fb_force(1:3,iatm) + work(1:3)
              end do
              do j = 1, numatoms(grplist(id2,inum))
                iatm = bondslist(j,grplist(id2,inum))
                work(1:3) =   coefi * coefgrp(j,id2) * dij(1:3,idist)
                fb_force(1:3,iatm) = fb_force(1:3,iatm) + work(1:3)
              end do

              ! virial
              ! 
              if (main_rank) then
                do i = 1, 3
                  work(i)= coefi * dij(i,idist)
                  do j = 1, 3
                    virial(j,i,1) = virial(j,i,1) - work(i) * dij(j,idist)
                  end do
                end do
              end if

            end do
          end do
        end if
        deallocate(grad_coef)
      end if
    end do

    return

  end subroutine compute_energy_restraints_fb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_restraint_coordinates
  !> @brief        get restraint coordinates
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    bonds_to_atom : bonds list to atom list for restrains
  !! @param[inout] bonds_coord   : restraint bonds coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_restraint_coordinates(domain, enefunc, coord, bonds_to_atom, &
                                       bonds_coord)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    real(wip),               intent(in)    :: coord(:,:,:)
    integer,                 intent(in)    :: bonds_to_atom(:)
    real(wp),                intent(inout) :: bonds_coord(:,:)

    ! local variable
    integer                  :: i, ix, icel, ncell_local

    integer(int2),   pointer :: id_g2l(:,:)


    id_g2l => domain%id_g2l

    ncell_local = domain%num_cell_local

    bonds_coord(1:3,1:enefunc%num_atoms_bonds_restraint) = 0.0_wp

    do i = 1, enefunc%num_atoms_bonds_restraint

      icel = id_g2l(1,bonds_to_atom(i))

      if (icel > 0 .and. icel <= ncell_local) then

        ix = id_g2l(2,bonds_to_atom(i))
        bonds_coord(1:3,i) = coord(1:3,ix,icel)

      end if
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, bonds_coord, &
                       3*enefunc%num_atoms_bonds_restraint, &
                       mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#endif

    return

  end subroutine get_restraint_coordinates

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_restraint_forces
  !> @brief        get restraint forces
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_restraint_forces(domain, enefunc, dist_to_atom, &
                                  dist_force, force)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    integer,                 intent(in)    :: dist_to_atom(:)
    real(wp),                intent(inout) :: dist_force(:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)

    ! local variable
    integer                  :: i, ix, icel, ncell_local
    integer(int2),   pointer :: id_g2l(:,:)


    id_g2l => domain%id_g2l

    ncell_local = domain%num_cell_local

    do i = 1, enefunc%num_atoms_bonds_restraint

      icel = id_g2l(1,dist_to_atom(i))

      if (icel > 0 .and. icel <= ncell_local) then

        ix = id_g2l(2,dist_to_atom(i))
        force(1:3,ix,icel,1) = force(1:3,ix,icel,1) + dist_force(1:3,i)

      end if
    end do

    return

  end subroutine get_restraint_forces

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2, val3)

    ! formal arguments
    real(dp),                intent(inout) :: val1(:), val2(:), val3

    ! local variables
    real(dp)                      :: before_reduce(7), after_reduce(7)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4:6) = val2(1:3)
    before_reduce(7)   = val3

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_city, ierror)
#else
    after_reduce(1:7) = before_reduce(1:7)
#endif

    val1(1:3)    = after_reduce(1:3)
    val2(1:3)    = after_reduce(4:6)
    val3         = after_reduce(7)

    return

  end subroutine reduce_com

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_restraints_fep
  !> @brief        calculate restraint energy for FEP
  !! @authors      HO
  !! @param[in]    get_coord  : flag for whether to get coordinates
  !! @param[in]    calc_force : flag for whether to calculate forces
  !! @param[in]    domain     : domain information
  !! @param[in]    boundary   : boundary information
  !! @param[inout] enefunc    : potential energy functions information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !! @param[inout] virial_ext : virial term of target systems
  !! @param[inout] eposi      : point restraint energy of target systems
  !! @param[inout] ermsd      : rmsd restraint energy of target systems
  !! @param[inout] rmsd       : rmsd
  !! @param[inout] ebonds     : 
  !! @param[inout] eemfit     :
  !! @param[inout] emcorr     :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_restraints_fep(get_coord, calc_force, domain, boundary,&
                                       enefunc, coord, force, virial,          &
                                       virial_ext, eposi, ermsd, rmsd, &
                                       ebonds, eemfit, emcorr)

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
    real(dp),                intent(inout) :: ermsd
    real(dp),                intent(inout) :: rmsd
    real(dp),                intent(inout) :: ebonds
    real(dp),                intent(inout) :: eemfit
    real(dp),                intent(inout) :: emcorr

    ! local variables
    real(dp)                 :: edist, eangle, edihed, etilt, eexp, erepul, efb
    real(wp),        pointer :: bonds_coord(:,:), bonds_force(:,:)
    integer,         pointer :: bonds_to_atom(:)

    ! FEP
    integer                  :: ncell, natom, id, i, j, k, ix, ic, jc
    integer                  :: omp_get_thread_num
    real(wp),        pointer :: f_fep(:,:,:,:)
    real(dp)                 :: v_fep(3,3,nthread)
    real(dp)                 :: v_ext_fep(3,3,nthread)
    real(dp)                 :: eposi_fep(nthread)

    call timer(TimerRestraint, TimerOn)

    bonds_to_atom => enefunc%restraint_bondslist_to_atomlist
    bonds_coord   => enefunc%restraint_bonds_coord
    bonds_force   => enefunc%restraint_bonds_force

    edist  = 0.0_dp
    eangle = 0.0_dp
    edihed = 0.0_dp
    ermsd  = 0.0_dp
    eemfit = 0.0_dp
    erepul = 0.0_dp
    efb    = 0.0_dp
    bonds_force(1:3,1:enefunc%num_atoms_bonds_restraint) = 0.0_wp

    if (get_coord) then
      call get_restraint_coordinates(domain, enefunc, coord, bonds_to_atom, &
                                     bonds_coord)
    end if

    ! FEP: initialize temporary arrays
    f_fep => domain%f_fep_omp
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom
    !$omp parallel do
    do id = 1, nthread
      f_fep(1:3,1:natom,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do
    v_fep(1:3,1:3,1:nthread) = 0.0_dp
    v_ext_fep(1:3,1:3,1:nthread) = 0.0_dp
    eposi_fep(1:nthread) = 0.0_dp
    !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        f_fep(1:3,ix,i,id+1) = force(1:3,ix,i,id+1)
      end do
    end do
    !$omp end parallel
    do id = 1, nthread
      v_fep(1:3,1:3,id)     = virial(1:3,1:3,id)
      v_ext_fep(1:3,1:3,id) = virial_ext(1:3,1:3,id)
      eposi_fep(id) = eposi(id)
    end do

    ! point restraint energy
    !
    if (enefunc%restraint_posi) &
    call compute_energy_restraints_pos(calc_force, domain, enefunc, coord,  &
                                       force, virial, virial_ext, eposi)

    if (enefunc%restraint_rmsd) &
    call compute_energy_restraints_rmsd(domain, enefunc, coord,             &
                                       force, virial, virial_ext, ermsd, rmsd)

    call compute_energy_restraints_dist(calc_force, enefunc, bonds_coord,   &
                                       bonds_force, virial, edist)

    call compute_energy_restraints_angle(calc_force, enefunc, bonds_coord,  &
                                       bonds_force, virial, eangle)

    call compute_energy_restraints_dihed(calc_force, enefunc, bonds_coord,  &
                                       bonds_force, virial, edihed)

    call compute_energy_experimental_restraint(calc_force, domain, boundary,&
                                       enefunc, coord, force, virial_ext,   &
                                       eemfit, emcorr)

    call compute_energy_restraints_repul(calc_force, enefunc, boundary, &
                                    bonds_coord, bonds_force, virial, erepul)

    call compute_energy_restraints_fb(calc_force, enefunc, boundary, &
                                    bonds_coord, bonds_force, virial, efb)

    ! FEP: multiply lambrest by energies (excluding eposi)
    ermsd  = enefunc%lambrest*ermsd
    eemfit = enefunc%lambrest*eemfit
    emcorr = enefunc%lambrest*emcorr
    ebonds = ebonds + enefunc%lambrest*&
             (edist + eangle + edihed + erepul + efb)

    if (calc_force) then
      call get_restraint_forces(domain, enefunc, bonds_to_atom,             &
                                bonds_force, force)
    end if

    ! FEP: multiply lambrest by forces and virials (and eposi) using OMP
    !$omp parallel default(shared) private(id, i, ix) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        f_fep(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - f_fep(1:3,ix,i,id+1)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + (enefunc%lambrest-1.0_wp)*f_fep(1:3,ix,i,id+1)
      end do
    end do
    !$omp end parallel
    do id = 1, nthread
      v_fep(1:3,1:3,id)  = virial(1:3,1:3,id) - v_fep(1:3,1:3,id)
      virial(1:3,1:3,id) = virial(1:3,1:3,id) + (real(enefunc%lambrest,dp)-1.0_dp)*v_fep(1:3,1:3,id)
      v_ext_fep(1:3,1:3,id)  = virial_ext(1:3,1:3,id) - v_ext_fep(1:3,1:3,id)
      virial_ext(1:3,1:3,id) = virial_ext(1:3,1:3,id) + (real(enefunc%lambrest,dp)-1.0_dp)*v_ext_fep(1:3,1:3,id)
      eposi_fep(id) = eposi(id) - eposi_fep(id)
      eposi(id)     = eposi(id) + (real(enefunc%lambrest,dp)-1.0_dp)*eposi_fep(id)
    end do

    call timer(TimerRestraint, TimerOff)

    return

  end subroutine compute_energy_restraints_fep

end module sp_energy_restraints_mod
