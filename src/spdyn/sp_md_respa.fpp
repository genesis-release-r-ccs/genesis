!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_md_respa_mod
!> @brief   perform molecular dynamics simulation with velocity verlet
!!          and respa algorithm
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_md_respa_mod

  use sp_output_mod
  use sp_update_domain_mod
  use sp_assign_velocity_mod
  use sp_dynvars_mod
  use sp_constraints_mod
  use sp_communicate_mod
  use sp_energy_mod
  use sp_remd_str_mod
  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use random_mod
  use math_libs_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_str_mod
  use sp_fep_energy_mod
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
  public  :: vverlet_respa_dynamics
  private :: initial_vverlet
  private :: integrate_vv1
  private :: integrate_vv2
  private :: nve_vv1
  private :: nve_vv2
  private :: mtk_barostat_vv1
  private :: langevin_thermostat_vv1
  private :: langevin_thermostat_vv2
  private :: langevin_barostat_vv1
  private :: langevin_barostat_vv2
  private :: update_barostat
  private :: update_barostat_mtk

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vverlet_respa_dynamics
  !> @brief        velocity verlet integrator using respa
  !! @authors      JJ
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraint information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vverlet_respa_dynamics(output, domain, enefunc, dynvars, &
                                    dynamics, pairlist, boundary,     &
                                    constraints, ensemble, comm, remd,&
                                    alchemy)

    use Ctim

    ! formal arguments
    type(s_output),            intent(inout) :: output
    type(s_domain),  target,   intent(inout) :: domain
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars), target,   intent(inout) :: dynvars
    type(s_dynamics),target,   intent(inout) :: dynamics
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary
    type(s_constraints),       intent(inout) :: constraints
    type(s_ensemble),          intent(inout) :: ensemble
    type(s_comm),              intent(inout) :: comm
    type(s_remd),              intent(inout) :: remd
    type(s_alchemy), optional, intent(inout) :: alchemy

    ! local variables
    real(wip)                :: simtim, temperature
    real(wip)                :: dt_short, dt_long
    real(wip)                :: min_time
    integer                  :: i, j, k, jx, nsteps, ii
    integer                  :: istep, multistep
    integer                  :: iseed, istart, iend
    integer                  :: min_k, k_min, k_max
    logical                  :: npt, npt1

    integer,         pointer :: atmcls_pbc(:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wip),       pointer :: coord_old(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:), force_omp(:,:,:,:)
    real(wip),       pointer :: force(:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial_extern(:,:)
    real(dp),        pointer :: virial(:,:), virial_long(:,:)
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: ncell, natom(:)
    logical,         pointer :: XI_RESPA, XO_RESPA


    atmcls_pbc    => domain%atmcls_pbc
    ncell         => domain%num_cell_local
    natom         => domain%num_atom
    mass          => domain%mass
    inv_mass      => domain%inv_mass
    coord         => domain%coord
    coord_ref     => domain%coord_ref
    coord_old     => domain%coord_old
    coord_pbc     => domain%translated
    force         => domain%force    
    force_omp     => domain%force_omp
    force_short   => domain%force_short
    force_long    => domain%force_long
    vel           => domain%velocity
    vel_ref       => domain%velocity_ref
    force_pbc     => domain%force_pbc
    virial_cell   => domain%virial_cellpair
    virial        => dynvars%virial
    virial_long   => dynvars%virial_long
    virial_extern => dynvars%virial_extern

    temperature   =  ensemble%temperature
    nsteps        =  dynamics%nsteps
    istart        =  dynamics%istart_step
    iend          =  dynamics%iend_step
    multistep     =  dynamics%elec_long_period
    dt_short      =  dynamics%timestep/AKMA_PS
    dt_long       =  real(multistep,wip)*dt_short
    simtim        =  dynamics%initial_time
    iseed         =  dynamics%iseed_init_velocity
    npt           =  ensemble%use_barostat

    if (domain%fep_use) then
      ! FEP: Copy lambda of enefunc to domain, because they are needed for 
      ! calculation of virial in SHAKE.
      domain%lambljA   = enefunc%lambljA
      domain%lambljB   = enefunc%lambljB
      domain%lambelA   = enefunc%lambelA
      domain%lambelB   = enefunc%lambelB
      domain%lambbondA = enefunc%lambbondA
      domain%lambbondB = enefunc%lambbondB
      domain%lambrest  = enefunc%lambrest
    end if

    if (abs(dynamics%initial_rmsd) < 0.001_wip)  &
      dynamics%initial_rmsd = real(dynvars%energy%rmsd,wip)
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt_long*dt_long)

    ! Check restart
    !
    if (.not. dynamics%restart) then

      call initial_velocity(temperature,            &
                            domain%num_atom,        &
                            domain%num_deg_freedom, &
                            domain%num_cell_local,  &
                            domain%id_g2l,          &
                            domain%inv_mass,        &
                            iseed,                  &
                            domain%velocity)

      call stop_trans_rotation(domain%num_cell_local,         &
                               domain%num_atom,               &
                               dynamics%stop_com_translation, &
                               dynamics%stop_com_rotation,    &
                               domain%mass,                   &
                               domain%coord,                  &
                               domain%velocity)

      call initial_vverlet(npt, output, enefunc, dynamics,       &
                           pairlist, boundary, ensemble, constraints, &
                           domain, dynvars, comm)

    else

      ! After 2nd cycle of REMD simulation (istart /= 1), this is skipped
      !
      if (istart == 1) then

        call communicate_coor(domain, comm)

        call compute_energy(domain, enefunc, pairlist, boundary, coord,  &
                            npt, .false., .true., .false.,               &
                            enefunc%nonb_limiter,                        &
                            dynvars%energy, atmcls_pbc,                  &
                            coord_pbc, force_short,                      &
                            force_long, force_omp, force_pbc,            &
                            virial_cell, virial, virial_long, virial_extern)
  
        call communicate_force(domain, comm, force_short)
        call communicate_force(domain, comm, force_long)

        if (constraints%water_type == TIP4) then
          call water_force_redistribution(constraints, domain, &
                                          force_short, virial)
          call water_force_redistribution(constraints, domain, &
                                          force_long, virial_long)
        end if

      end if
    
    end if

    call mpi_barrier(mpi_comm_country, ierror)
    Timt=0.0d0
    Timb=0.0d0
    mpi_tot_tran=0.d0

    ! Main loop
    !
    do i = istart, iend-multistep+1, multistep

      mpi_tran = 0.0_dp
      mpi_bari = 0.0_dp
      dynvars%time = dynamics%timestep * real(i-1,dp)
      dynvars%step = i - 1
      enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

      ! output trajectory and restart data
      !
      if (i > istart) then
        call output_md(output, dynamics, boundary, pairlist, ensemble, &
                       constraints, dynvars, domain, enefunc, remd)
        ! FEP: Compute and output energy differnces between adjacent states
        if ((domain%fep_use).and.(present(alchemy))) then
          if (dynamics%fepout_period > 0) then
            if ((mod(i-1,dynamics%fepout_period) == 0) .and. &
              (mod(i-1,dynamics%eneout_period) == 0)) then
              if (i-1 > dynamics%equilsteps) then
                call compute_fep_energy(domain, enefunc, dynvars, pairlist, &
                  ensemble, boundary, alchemy, remd)
                call output_fep_energy(output, enefunc, dynvars)
              end if
            end if
          end if
        end if
      end if

      ! Inner integrator
      !
      do j = 1, multistep

        istep = i + j - 1

        if (dynamics%target_md .or. dynamics%steered_md) &
          enefunc%rmsd_target = dynamics%initial_rmsd &
                              + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                               *real(istep,wip)/real(nsteps,wip)

        ! save coordiantes for constraint
        !
        call timer(TimerIntegrator, TimerOn)

        ! VV1
        !
        call integrate_vv1(dynamics, istep, istart, j, dt_long, dt_short,  &
                           ensemble, domain, constraints, boundary, dynvars)

        call timer(TimerIntegrator, TimerOff)

        if (mod(istep-1, dynamics%eneout_period) == 0 .and. istep > istart) then
          call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                               dynvars)
          call output_dynvars (output, enefunc, dynvars, ensemble)
        end if

        ! Remove translational and rotational motion about COM(t + dt)
        !
        if (dynamics%stoptr_period > 0) then
          if (mod(istep-1,dynamics%stoptr_period) == 0) then
              call stop_trans_rotation(domain%num_cell_local,         &
                                       domain%num_atom,               &
                                       dynamics%stop_com_translation, &
                                       dynamics%stop_com_rotation,    &
                                       domain%mass,                   &
                                       domain%coord,                  &
                                       domain%velocity)
          end if
        end if

        ! cell migration and update cell pairlist
        !
        if (dynamics%nbupdate_period > 0 .and. &
            j == multistep .and. i > istart) then

          if (domain%fep_use) then
            ! FEP
            call domain_interaction_update_fep(istep, dynamics%nbupdate_period, &
                                         domain, enefunc, pairlist, boundary, &
                                         constraints, comm)
          else
            call domain_interaction_update(istep, dynamics%nbupdate_period,     &
                                         domain, enefunc, pairlist, boundary, &
                                         constraints, comm)
          end if

        end if

        ! Simulated annealing
        !
        if (dynamics%anneal_period > 0 .and. i > istart) then
          if (mod(istep-1,dynamics%anneal_period) == 0) &
            call simulated_annealing_vverlet(dynamics, ensemble)
        end if

        ! short range forces
        !
        call timer(TimerIntegrator, TimerOn)
        call timer(TimerComm1, TimerOn)

        call communicate_coor(domain, comm)

        call timer(TimerComm1, TimerOff)
        call timer(TimerIntegrator, TimerOff)

        if (j < multistep) then

          pairlist%univ_update = 0
          enefunc%rpath_sum_mf_flag = enefunc%rpath_flag
          npt1 = npt .and. mod(istep-1,dynamics%baro_period) == 0

          call compute_energy_short(domain, enefunc, pairlist, boundary, coord,&
                                    npt1,                                      &
                                    dynvars%energy,                            &
                                    atmcls_pbc,                                &
                                    coord_pbc,                                 &
                                    force_short,                               &
                                    force_omp,                                 &
                                    force_pbc,                                 &
                                    virial_cell,                               &
                                    virial, virial_extern)


          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)

          call communicate_force(domain, comm, force_short)

          if (constraints%water_type == TIP4) &
            call water_force_redistribution(constraints, domain, force_short, &
                                            virial)
          call timer(TimerComm2, TimerOff)
          call timer(TimerIntegrator, TimerOff)

        end if
        
        ! full step velocities with foce_short
        ! long range force for last inner step
        !
        if (j == multistep) then

          enefunc%rpath_sum_mf_flag = enefunc%rpath_flag
          npt1 = npt .and. mod(istep,dynamics%baro_period) == 0

          call compute_energy(domain, enefunc, pairlist, boundary, coord,      &
                              npt1,.false.,                                    &
                              mod(istep,dynamics%eneout_period) == 0,          &
                              .false.,                                         &
                              enefunc%nonb_limiter,                            &
                              dynvars%energy,                                  &
                              atmcls_pbc,                                      &
                              coord_pbc,                                       &
                              force_short,                                     &
                              force_long,                                      &
                              force_omp,                                       &
                              force_pbc,                                       &
                              virial_cell,                                     &
                              virial, virial_long, virial_extern)

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)

          call communicate_force(domain, comm, force_long)
          call communicate_force(domain, comm, force_short)

          if (constraints%water_type == TIP4) then
            call water_force_redistribution(constraints, domain, &
                                            force_short, virial)
            call water_force_redistribution(constraints, domain, &
                                            force_long, virial_long)
          end if
          call timer(TimerComm2, TimerOff)
          call timer(TimerIntegrator, TimerOff)

        end if

        ! VV2
        !
        call integrate_vv2(dynamics, istep, j, dt_long, dt_short, ensemble,    &
                           domain, constraints, boundary, dynvars)

      end do

    end do

    ! for final output
    !
    dynvars%time = dynamics%timestep * real(iend,dp)
    dynvars%step = iend

    ! output trajectory(t + dt) and restart data
    !
    call output_md(output, dynamics, boundary, pairlist, ensemble, &
                   constraints, dynvars, domain, enefunc, remd)

!   ! output parallel I/O restart
!   !
!   call output_prst_md(output, enefunc, dynamics, boundary, &
!                       dynvars, domain, constraints)

    ! FEP: Compute and output energy differnces between adjacent states
    if ((domain%fep_use).and.(present(alchemy))) then
      if (dynamics%fepout_period > 0) then
        call compute_fep_energy(domain, enefunc, dynvars, pairlist, &
          ensemble, boundary, alchemy, remd)
        call output_fep_energy(output, enefunc, dynvars)
      end if
    end if

    call integrate_vv1(dynamics, iend+1, istart, 1, dt_long, dt_short, ensemble, &
                       domain, constraints, boundary, dynvars)

    call coord_vel_ref(domain, dynvars)
    if (ensemble%tpcontrol == TpcontrolNHC) &
      dynvars%nh_velocity(1:5) = dynvars%nh_velocity_ref(1:5) 

    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                         dynvars)
    call output_dynvars (output, enefunc, dynvars, ensemble)

    return

  end subroutine vverlet_respa_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_vverlet
  !> @brief        compute the first step (0+dt)
  !! @authors      JJ
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynamics    : dynamics information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

   subroutine initial_vverlet(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, domain,  &
                              dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: temperature
    real(wip)                :: simtim, dt, half_dt, friction
    integer                  :: i, k, l, ix, ncell

    integer,         pointer :: atmcls_pbc(:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:)
    real(dp),        pointer :: virial_const(:,:), virial_extern(:,:)
    real(dp),        pointer :: virial(:,:), virial_long(:,:)
    real(dp),        pointer :: kin_half(:), ekin_half
    real(dp),        pointer :: kin_ref(:), ekin_ref
    real(dp),        pointer :: kin_full(:), ekin_full
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: natom(:)


    atmcls_pbc    => domain%atmcls_pbc
    natom         => domain%num_atom
    coord         => domain%coord
    coord_ref     => domain%coord_ref
    coord_pbc     => domain%translated
    force         => domain%force
    force_short   => domain%force_short
    force_long    => domain%force_long 
    force_omp     => domain%force_omp
    vel           => domain%velocity
    vel_ref       => domain%velocity_ref
    mass          => domain%mass
    inv_mass      => domain%inv_mass
    force_pbc     => domain%force_pbc
    virial_cell   => domain%virial_cellpair
    virial        => dynvars%virial
    virial_long   => dynvars%virial_long
    virial_extern => dynvars%virial_extern
    virial_const  => dynvars%virial_const 
    kin_half      => dynvars%kin_half
    ekin_half     => dynvars%ekin_half
    kin_full      => dynvars%kin_full
    ekin_full     => dynvars%ekin_full
    kin_ref       => dynvars%kin_ref
    ekin_ref      => dynvars%ekin_ref

    dt            =  dynamics%timestep/AKMA_PS
    half_dt       =  0.5_wip * dt
    simtim        =  dynamics%initial_time
    temperature   =  ensemble%temperature
    friction      =  ensemble%gamma_t * AKMA_PS
    ncell         =  domain%num_cell_local

    dynvars%time  = simtim
    dynvars%step  = 0

    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        vel_ref(1:3,ix,i)   = vel(1:3,ix,i)
      end do
    end do

    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref, &
                               domain, constraints, coord, vel,            &
                               virial_const)

      do i = 1, ncell
        do ix = 1, natom(i)
          coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        end do
      end do

    else if (constraints%water_type == TIP4) then

      call decide_dummy(domain, constraints, coord)

    end if

    ! calculate energy(0) and forces(0)
    !
    call communicate_coor(domain, comm)

    call compute_energy(domain, enefunc, pairlist, boundary, coord,     &
                        npt, .false., .true., .false.,                  &
                        enefunc%nonb_limiter,                           &
                        dynvars%energy, atmcls_pbc,                     &
                        coord_pbc, force_short,                         &
                        force_long, force_omp, force_pbc,               &
                        virial_cell, virial, virial_long, virial_extern)

    call communicate_force(domain, comm, force_short)
    call communicate_force(domain, comm, force_long)
    if (constraints%water_type == TIP4) then
      call water_force_redistribution(constraints, domain, &
                                      force_short, virial)
      call water_force_redistribution(constraints, domain, &
                                      force_long, virial_long)
    end if

    ! force
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        force(1:3,ix,i) = force_short(1:3,ix,i) + force_long(1:3,ix,i)
      end do
    end do

    ! velocities at 0 + dt/2
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = vel_ref(1:3,ix,i) &
                      + half_dt*force(1:3,ix,i)*inv_mass(ix,i)
        coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
      end do
    end do

    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref,    &
                               domain, constraints, coord, vel, virial_const)
      virial_const(1:3,1:3) = 2.0_dp * virial_const(1:3,1:3)
    end if

    call calc_kinetic(ncell, natom, mass, vel, kin_half, ekin_half)

    ! velocities at 0 + dt/2
    ! 
    do i = 1, ncell
      do ix = 1, natom(i)
        vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                                &
                          - half_dt*force(1:3,ix,i)*inv_mass(ix,i)
        coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
      end do
    end do

    call compute_constraints(ConstraintModeLEAP, .false., -dt, coord_ref,    &
                             domain, constraints, coord, vel_ref,            &
                             virial_const)
    virial_const(1:3,1:3) = 2.0_dp * virial_const(1:3,1:3)

    call calc_kinetic(ncell, natom, mass, vel_ref, kin_ref, ekin_ref)

    ! calculate velocity(0)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        vel_ref(1:3,ix,i) = 0.5_wip*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
      end do
    end do

    call calc_kinetic(ncell, natom, mass, vel_ref, kin_full, ekin_full)

    ! vel <= updated velocities (0) and coord <= constrained coordinates(0)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = vel_ref(1:3,ix,i)
        coord(1:3,ix,i) = coord_ref(1:3,ix,i)
      end do
    end do

    ! output dynvars(0)
    !
    dynvars%time = 0.0_wip
    dynvars%step = 0

    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble)

    dynamics%restart = .true.

    return

  end subroutine initial_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv1
  !> @brief        VV1 with thermostat/barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step        
  !! @param[in]    istart      : dynamics start step
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv1(dynamics, istep, istart, inner_step, dt_long, &
                           dt_short, ensemble, domain, constraints, boundary, &
                           dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: istart
    integer,                 intent(in)    :: inner_step
    real(wip),               intent(in)    :: dt_long
    real(wip),               intent(in)    :: dt_short
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars

    integer  :: alloc_stat, ncell
    real(wip):: dt_therm, dt_baro


    ! select ensemble, thermostat and barostat
    !
    if (.not. ensemble%group_tp .and. constraints%rigid_bond) &
      call error_msg('In VRES, ensemble%group_tp=NO does not work')

    select case (ensemble%ensemble)

    case (EnsembleNVE)

      call nve_vv1(dynamics, istep, dt_long, dt_short, domain, constraints, &
                   dynvars)

    case (EnsembleNVT)

      select case (ensemble%tpcontrol)

      case (TpcontrolBerendsen, TpcontrolBussi, TpcontrolNHC)

        call vel_rescaling_thermostat_vv1(inner_step, dt_long, dt_short,     &
                                          dynamics, istep, istart, ensemble, &
                                          domain, constraints, dynvars)

      case default

        call error_msg('Thermostat is not properly defined')

      end select

    case (EnsembleNPT:EnsembleNPgT)

      select case (ensemble%tpcontrol)

      case (TpcontrolBussi, TpcontrolBerendsen, TpcontrolNHC)

        call mtk_barostat_vv1(dynamics, istep, istart, dt_long, dt_short,    &
                              ensemble, domain, constraints,         &
                              boundary, dynvars)

      case default

        call error_msg('Available thermostat/barostat for NPT are ' //  &
                       'only Bussi/Berendsen/NHC')

      end select

    case default

      call error_msg('Ensemble is not defined properly')

    end select

    return

  end subroutine integrate_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv2
  !> @brief        VV2 with thermostat/barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step        
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv2(dynamics, istep, inner_step, dt_long, dt_short,  &
                           ensemble, domain, constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: inner_step
    real(wip),               intent(in)    :: dt_long 
    real(wip),               intent(in)    :: dt_short
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars

    real(wip):: dt_therm


    dt_therm = dt_short * real(dynamics%thermo_period,wip)

    select case (ensemble%ensemble)

    case (EnsembleNVE)

      call nve_vv2(dynamics, inner_step, dt_long, dt_short, domain, &
                   constraints, dynvars)

    case (EnsembleNVT)

      select case (ensemble%tpcontrol)

      case (TpcontrolLangevin)

        call langevin_thermostat_vv2(dynamics, istep, inner_step,    &
                                     dt_long, dt_short, ensemble,    &
                                     domain, constraints, dynvars)

      case (TpcontrolBussi, TpcontrolBerendsen, TpcontrolNHC)

        call nve_vv2(dynamics, inner_step, dt_long, dt_short, domain,  &
                     constraints, dynvars)

      case default

        call error_msg('Thermostat is not properly defined')

      end select

    case (EnsembleNPT:EnsembleNPgT)

      select case (ensemble%tpcontrol)

        case (TpcontrolBussi, TpcontrolBerendsen, TpcontrolNHC)

          call mtk_barostat_vv2(dynamics, istep, dt_long, dt_short, &
                                ensemble, domain,                   &
                                constraints, boundary, dynvars)
  
        case default

        call error_msg('Available thermostat/barostat for NPT are' // &
                       'only Bussi/Berendsen/NHC')

      end select

    case default

      call error_msg('Ensemble is not defined properly')

    end select

    return

  end subroutine integrate_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    coord_vel_ref
  !> @brief        output in the last step
  !! @authors      JJ
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine coord_vel_ref(domain, dynvars)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix

    integer,         pointer :: ncell, natom(:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), vel_half(:,:,:)


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_half    => domain%velocity_half

    ! kinetic energy and virial output
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        vel_half(1,ix,i) = vel(1,ix,i)
        vel_half(2,ix,i) = vel(2,ix,i)
        vel_half(3,ix,i) = vel(3,ix,i)
        vel(1,ix,i) = vel_ref(1,ix,i)
        vel(2,ix,i) = vel_ref(2,ix,i)
        vel(3,ix,i) = vel_ref(3,ix,i)
        coord(1,ix,i) = coord_ref(1,ix,i)
        coord(2,ix,i) = coord_ref(2,ix,i)
        coord(3,ix,i) = coord_ref(3,ix,i)
      end do
    end do

    return

  end subroutine coord_vel_ref

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv1
  !> @brief        VV1 with NVE
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @@aram[in]    istep       : dynamics step
  !! @param[in]    dt_long     : long time step        
  !! @param[in]    dt_short    : short time step        
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv1(dynamics, istep, dt_long, dt_short, domain, constraints, &
                     dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    real(wip),               intent(in)    :: dt_long
    real(wip),               intent(in)    :: dt_short
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix, k, l
    real(wip)                :: half_dt_long, half_dt_short
    real(wip)                :: factor

    integer,         pointer :: ncell, natom(:), nwater(:), water_list(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), vel_half(:,:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:), viri_group(:,:)
    real(dp),        pointer :: virial_long(:,:)
    real(dp),        pointer :: kin(:), kin_full(:), kin_half(:), kin_ref(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    nwater      => domain%num_water
    water_list  => domain%water_list
    mass        => domain%mass
    inv_mass    => domain%inv_mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_half    => domain%velocity_half
    force       => domain%force     
    force_long  => domain%force_long
    force_short => domain%force_short
    virial      => dynvars%virial
    virial_long => dynvars%virial_long
    viri_const  => dynvars%virial_const
    viri_group  => dynvars%virial_group
    ekin_full   => dynvars%ekin_full
    ekin_half   => dynvars%ekin_half
    ekin_ref    => dynvars%ekin_ref
    ekin        => dynvars%ekin
    kin         => dynvars%kin
    kin_full    => dynvars%kin_full
    kin_half    => dynvars%kin_half
    kin_ref     => dynvars%kin_ref

    half_dt_long   = 0.5_wip * dt_long  
    half_dt_short  = 0.5_wip * dt_short

    call timer(TimerUpdate, TimerOn)

    ! reference coordinates and velocities
    !
    call copy_coord_vel(ncell, natom, coord, vel, coord_ref, vel_ref)

    ! kinetic energy and virial output
    !
    if (mod(istep-1, dynamics%eneout_period) == 0) then

      !$omp parallel do private(i, ix, factor)
      do i = 1, ncell
        do ix = 1, natom(i)
          factor = half_dt_short * inv_mass(ix,i)
          force(1,ix,i) = force_long(1,ix,i)+force_short(1,ix,i)
          force(2,ix,i) = force_long(2,ix,i)+force_short(2,ix,i)
          force(3,ix,i) = force_long(3,ix,i)+force_short(3,ix,i)
          vel_half(1,ix,i) = factor*force(1,ix,i)
          vel_half(2,ix,i) = factor*force(2,ix,i)
          vel_half(3,ix,i) = factor*force(3,ix,i)
        end do
      end do
      !$omp end parallel do

      call compute_kin_group(constraints, ncell, nwater, water_list, &
                             mass, vel_half, kin_half, ekin_half)
      call compute_kin_group(constraints, ncell, nwater, water_list, &
                             mass, vel, kin_full, ekin_full)
      ekin = ekin_full + 2.0_dp*ekin_half/3.0_dp
      kin(1:3) = kin_full(1:3) + kin_half(1:3)

      viri_group(1:3,1:3) = 0.0_dp
      call compute_virial_group(constraints, ncell, nwater, water_list, &
                                mass, coord_ref, force, viri_group)
      virial(1:3,1:3) = virial(1:3,1:3) + viri_group(1:3,1:3)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_long(1:3,1:3)

    end if

    ! VV1
    !
    if (mod(istep-1, dynamics%elec_long_period) == 0) then

      !$omp parallel do private(i, ix, factor)
      do i = 1, ncell
        do ix = 1, natom(i)
          factor = half_dt_long * inv_mass(ix,i)
          vel(1,ix,i)   = vel(1,ix,i)   + factor*force_long(1,ix,i)
          vel(2,ix,i)   = vel(2,ix,i)   + factor*force_long(2,ix,i)
          vel(3,ix,i)   = vel(3,ix,i)   + factor*force_long(3,ix,i)
        end do
      end do
      !$omp end parallel do

    end if

    !$omp parallel do private(i, ix, factor)
    do i = 1, ncell
      do ix = 1, natom(i)
        factor = half_dt_short * inv_mass(ix,i)
        vel(1,ix,i)   = vel(1,ix,i)   + factor*force_short(1,ix,i)
        vel(2,ix,i)   = vel(2,ix,i)   + factor*force_short(2,ix,i)
        vel(3,ix,i)   = vel(3,ix,i)   + factor*force_short(3,ix,i)
        coord(1,ix,i) = coord(1,ix,i) + dt_short*vel(1,ix,i)
        coord(2,ix,i) = coord(2,ix,i) + dt_short*vel(2,ix,i)
        coord(3,ix,i) = coord(3,ix,i) + dt_short*vel(3,ix,i)
      end do
    end do
    !$omp end parallel do

    call timer(TimerUpdate, TimerOff)

    ! RATTLE VV1
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt_short,  &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
    end if

    return

  end subroutine nve_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv2
  !> @brief        VV2 with NVE
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : long time step        
  !! @param[in]    dt_short    : short time step        
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv2(dynamics, inner_step, dt_long, dt_short, domain, &
                     constraints, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: inner_step
    real(wip),               intent(in)    :: dt_long
    real(wip),               intent(in)    :: dt_short
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix, k, l, ncell
    real(wip)                :: half_dt_long, half_dt_short
    real(wip)                :: factor, vel_change(3)
    real(dp)                 :: virial_constraint(3)

    integer,         pointer :: natom(:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)


    ! use pointers
    !
    ncell       =  domain%num_cell_local
    natom       => domain%num_atom
    mass        => domain%mass
    inv_mass    => domain%inv_mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_full
    force_long  => domain%force_long
    force_short => domain%force_short
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const

    half_dt_long   =  0.5_wip * dt_long 
    half_dt_short  =  0.5_wip * dt_short

    call timer(TimerUpdate, TimerOn)

    ! VV2
    !
    if (inner_step == dynamics%elec_long_period) then

      !$omp parallel do private(i, ix, factor)
      do i = 1, ncell
        do ix = 1, natom(i)
          factor = half_dt_long * inv_mass(ix,i)
          vel(1,ix,i)   = vel(1,ix,i)   + factor*force_long(1,ix,i)
          vel(2,ix,i)   = vel(2,ix,i)   + factor*force_long(2,ix,i)
          vel(3,ix,i)   = vel(3,ix,i)   + factor*force_long(3,ix,i)
        end do
      end do
      !$omp end parallel do

    end if

    !$omp parallel do private(i, ix, factor)
    do i = 1, ncell
      do ix = 1, natom(i)
        factor = half_dt_short * inv_mass(ix,i)
        vel(1,ix,i)   = vel(1,ix,i)   + factor*force_short(1,ix,i)
        vel(2,ix,i)   = vel(2,ix,i)   + factor*force_short(2,ix,i)
        vel(3,ix,i)   = vel(3,ix,i)   + factor*force_short(3,ix,i)
      end do
    end do
    !$omp end parallel do

    call timer(TimerUpdate, TimerOff)

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               vel, viri_const)
    end if

    return

  end subroutine nve_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_rescaling_thermostat_vv1
  !> @brief        VV1 with Bussi's thermostat (group temp or no constraint)
  !! @authors      JJ
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : long time step
  !! @param[in]    dt_short    : short time step
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step        
  !! @param[in]    istart      : dynamics start step
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vel_rescaling_thermostat_vv1(inner_step, dt_long, dt_short,      &
                                          dynamics, istep, istart, ensemble,  &
                                          domain, constraints, dynvars)

    ! formal arguments
    integer,                 intent(in)    :: inner_step
    real(wip),               intent(in)    :: dt_long
    real(wip),               intent(in)    :: dt_short
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: istart
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    logical                  :: calc_thermostat
    integer                  :: i, ix, k, l
    integer                  :: num_degree
    real(wip)                :: temp0, tau_t
    real(wip)                :: half_dt_long, half_dt_short, dt_therm
    real(wip)                :: factor, rr
    real(wip)                :: scale_vel, scale_vel2
    real(dp)                 :: kin_full_ref(3), ekin_full_ref

    integer,         pointer :: ncell, natom(:), nwater(:), water_list(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), vel_half(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:), viri_group(:,:)
    real(dp),        pointer :: virial_long(:,:)
    real(dp),        pointer :: kin(:), kin_full(:), kin_half(:), kin_ref(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    nwater      => domain%num_water
    water_list  => domain%water_list
    mass        => domain%mass
    inv_mass    => domain%inv_mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_half    => domain%velocity_half
    force_long  => domain%force_long
    force_short => domain%force_short
    force       => domain%force      
    virial      => dynvars%virial
    virial_long => dynvars%virial_long
    viri_const  => dynvars%virial_const
    viri_group  => dynvars%virial_group
    ekin_full   => dynvars%ekin_full
    ekin_half   => dynvars%ekin_half
    ekin_ref    => dynvars%ekin_ref 
    ekin        => dynvars%ekin
    kin         => dynvars%kin
    kin_full    => dynvars%kin_full
    kin_half    => dynvars%kin_half
    kin_ref     => dynvars%kin_ref

    half_dt_long   = 0.5_wip * dt_long 
    half_dt_short  = 0.5_wip * dt_short 
    dt_therm = dt_short * real(dynamics%thermo_period,wip)
    temp0    = ensemble%temperature
    tau_t    = ensemble%tau_t / AKMA_PS

    calc_thermostat = mod(istep-1,dynamics%thermo_period) == 0 .and. istep > istart

    if (ensemble%group_tp) then
      num_degree = domain%num_group_freedom
    else
      num_degree = domain%num_deg_freedom
    end if

    scale_vel  = 1.0_wip
    scale_vel2 = 1.0_wip

    call timer(TimerUpdate, TimerOn)

    ! reference coordinates and velocities
    !
    call copy_coord_vel(ncell, natom, coord, vel, coord_ref, vel_ref)

    ! thermostat
    !
    if (calc_thermostat) then

      ! half-time step velocity
      do i = 1, ncell
        do ix = 1, natom(i)
          factor = half_dt_short * inv_mass(ix,i)
          force(1,ix,i) = force_long(1,ix,i) + force_short(1,ix,i)
          force(2,ix,i) = force_long(2,ix,i) + force_short(2,ix,i)
          force(3,ix,i) = force_long(3,ix,i) + force_short(3,ix,i)
          vel_half(1,ix,i) = factor*force(1,ix,i)
          vel_half(2,ix,i) = factor*force(2,ix,i)
          vel_half(3,ix,i) = factor*force(3,ix,i)
        end do
      end do
      if (ensemble%group_tp) then
        call compute_kin_group(constraints, ncell, nwater, water_list, mass, &
                               vel_half, kin_half, ekin_half)
        call compute_kin_group(constraints, ncell, nwater, water_list, mass, &
                               vel_ref, kin_full, ekin_full)
      else
        call calc_kinetic(ncell, natom, mass, vel_half, kin_half, ekin_half)
        call calc_kinetic(ncell, natom, mass, vel_ref, kin_full, ekin_full)
      end if

      ekin = ekin_full + 2.0_dp*ekin_half/3.0_dp
      kin(1:3) = kin_full(1:3) + kin_half(1:3)

      ! calculate scaling factor
      !
      if (ensemble%tpcontrol == TpcontrolBussi) then
        rr = real(random_get_gauss(),wip)
        call vel_scale_bussi(num_degree, dt_therm, tau_t, temp0, ekin, rr, &
                             scale_vel)
      else if (ensemble%tpcontrol == TpcontrolBerendsen) then
        call vel_scale_berendsen(num_degree, dt_therm, tau_t, temp0, ekin, &
                                 scale_vel)
      else if (ensemble%tpcontrol == TpcontrolNHC) then
        call vel_scale_nhc(num_degree, dt_therm, tau_t, temp0, ekin, &
                           ensemble, dynvars, scale_vel)
      end if

      scale_vel2 = scale_vel * scale_vel

    end if

    ekin = ((1.0_dp+2.0_dp*scale_vel2)*ekin_full + 2.0_dp*ekin_half)/3.0_dp
    kin(1:3) = 0.5_dp*((1.0_dp+scale_vel2)*kin_full(1:3)) + kin_half(1:3)

    if (ensemble%group_tp) then
      call update_vel_group(constraints, ncell, nwater, water_list, &
                            scale_vel, mass, vel)
    else
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1,ix,i) = vel(1,ix,i)*scale_vel
          vel(2,ix,i) = vel(2,ix,i)*scale_vel
          vel(3,ix,i) = vel(3,ix,i)*scale_vel
        end do
      end do
    end if
 
    ! VV1
    if (inner_step == 1) then

      !$omp parallel do private(i, ix, k, l, factor)
      do i = 1, ncell
        do ix = 1, natom(i)
          factor = half_dt_long * inv_mass(ix,i)
          vel(1,ix,i)   = vel(1,ix,i)   + factor*force_long(1,ix,i)
          vel(2,ix,i)   = vel(2,ix,i)   + factor*force_long(2,ix,i)
          vel(3,ix,i)   = vel(3,ix,i)   + factor*force_long(3,ix,i)
        end do
      end do

    end if

    !$omp parallel do private(i, l, k, ix, factor)
    do i = 1, ncell
      do ix = 1, natom(i)
        factor = half_dt_short * inv_mass(ix,i)
        vel(1,ix,i)   = vel(1,ix,i)   + factor*force_short(1,ix,i)
        vel(2,ix,i)   = vel(2,ix,i)   + factor*force_short(2,ix,i)
        vel(3,ix,i)   = vel(3,ix,i)   + factor*force_short(3,ix,i)
        coord(1,ix,i) = coord_ref(1,ix,i) + dt_short*vel(1,ix,i)
        coord(2,ix,i) = coord_ref(2,ix,i) + dt_short*vel(2,ix,i)
        coord(3,ix,i) = coord_ref(3,ix,i) + dt_short*vel(3,ix,i)
      end do
    end do
    !$omp end parallel do
    call timer(TimerUpdate, TimerOff)

    ! RATTLE VV1
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt_short,  &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
    end if

    if (mod(istep-1, dynamics%eneout_period) == 0 .and. ensemble%group_tp) then
      viri_group(1:3,1:3) = 0.0_dp
      call compute_virial_group(constraints, ncell, nwater, water_list, &
                                mass, coord_ref, force, viri_group)
      virial(1:3,1:3) = virial(1:3,1:3) + viri_group(1:3,1:3)
      virial(1:3,1:3) = virial(1:3,1:3) + virial_long(1:3,1:3)
    end if

    return

  end subroutine vel_rescaling_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv1
  !> @brief        control temperature using Langevin thermostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv1(dynamics, istep, inner_step, &
                                     dt_long, dt_short, ensemble, &
                                     domain, constraints, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(wip),                intent(in)    :: dt_long, dt_short
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: inv_dt
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: temp0, gamma_t, scale_v
    real(wip)                :: factor, sigma
    real(wip)                :: rsq, v1, v2, grandom(1:3)
    real(wip)                :: kBT, vel_tmp(1:3)
    integer                  :: j, jx, k, l, ncell

    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: random_f(:,:,:), mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: temporary(:,:,:)
    real(dp) ,       pointer :: virial(:,:), viri_const(:,:)
    integer,         pointer :: natom(:)


    inv_dt        =  1.0_wip / dt_short
    dt_therm      =  dt_short * real(dynamics%thermo_period, wip)
    half_dt_therm =  dt_therm / 2.0_wip
    temp0         =  ensemble%temperature
    gamma_t       =  ensemble%gamma_t *AKMA_PS
    random_f      => ensemble%random_force
    ncell         =  domain%num_cell_local

    natom         => domain%num_atom
    mass          => domain%mass
    inv_mass      => domain%inv_mass
    coord         => domain%coord
    coord_ref     => domain%coord_ref
    vel           => domain%velocity
    vel_ref       => domain%velocity_ref
    force_long    => domain%force_long
    force_short   => domain%force_short
    temporary     => domain%coord_old
    virial        => dynvars%virial
    viri_const    => dynvars%virial_const

    ! setup variables
    !
    kBT      = KBOLTZ * temp0

    ! scale factor for velocities
    !
    scale_v = exp(-gamma_t*half_dt_therm)

    ! random force
    !
    if (istep == 1) then

      factor  = 1.0_wip - scale_v*scale_v
      factor  = factor*KBOLTZ*temp0

      do j = 1, ncell
        do jx = 1, natom(j)

          sigma = sqrt(factor*inv_mass(jx,j))
          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

    end if

    ! Thermostat before VV1
    !
    if (mod(istep-1, dynamics%thermo_period) == 0) then

      ! Thermostat
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)*scale_v
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do

    end if

    ! VV1 with long range force
    !
    if (inner_step == 1) then
      do j = 1, ncell
        do jx = 1, natom(j)
          factor = 0.5_wip * dt_long * inv_mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
      end do
    end if

    ! VV1
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        factor = 0.5_wip * dt_short * inv_mass(jx,j)
        vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        coord(1:3,jx,j) = coord_ref(1:3,jx,j)+vel(1:3,jx,j)*dt_short
      end do
    end do

    ! Coordinate constraint (RATTLE VV1)
    !
    if (constraints%rigid_bond) then
      do j = 1, ncell
        do jx = 1, natom(j)
          temporary(1,jx,j) = coord(1,jx,j)
          temporary(2,jx,j) = coord(2,jx,j)
          temporary(3,jx,j) = coord(3,jx,j)
        end do
      end do

      call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
      viri_const(1:3,1:3) = 2.0_dp * viri_const(1:3,1:3)

      do j = 1, ncell
        do jx = 1, natom(j)
          vel_tmp(1:3)  = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
          vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
        end do
      end do
    end if

    return

  end subroutine langevin_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv2
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv2(dynamics, istep, inner_step, dt_long,    &
                                     dt_short, ensemble, domain, constraints, &
                                     dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(wip),                intent(in)    :: dt_long
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: inv_dt, temp0
    real(wip)                :: scale_v, factor
    real(wip)                :: gamma_t
    real(wip)                :: sigma
    real(wip)                :: v1, v2, rsq, grandom(1:3)
    real(wip)                :: dt_therm, half_dt_therm
    integer                  :: j, k, l, jx, ncell

    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(dp),        pointer :: viri_const(:,:)
    real(wip),       pointer :: random_f(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    integer,         pointer :: natom(:)


    inv_dt     =  1.0_wip/dt_short
    temp0      =  ensemble%temperature
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    random_f   => ensemble%random_force
    ncell      =  domain%num_cell_local

    mass       => domain%mass
    inv_mass   => domain%inv_mass
    natom      => domain%num_atom
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force_long => domain%force_long
    force_short=> domain%force_short
    viri_const => dynvars%virial_const

    ! time step
    !
    dt_therm       = dt_short * real(dynamics%thermo_period,wip)
    half_dt_therm  = dt_therm / 2.0_wip

    ! VV2
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        factor = 0.5_wip * dt_short*inv_mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
      end do
    end do

    if (inner_step == dynamics%elec_long_period) then

      ! VV2 (long range force)
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          factor = 0.5_wip * dt_long * inv_mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
      end do

      if (mod(istep,  dynamics%thermo_period) == 0) then

        ! random force
        !
        scale_v = exp(-gamma_t*half_dt_therm)
        factor   = 1.0_wip - scale_v*scale_v
        factor   = factor*KBOLTZ*temp0/2.0_wip

        do j = 1, ncell
          do jx = 1, natom(j)

            sigma = sqrt(factor*inv_mass(jx,j))

            rsq = 2.0_wip

            do while (rsq >= 1.0_wip)
              v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
            grandom(1) = rsq * v1
            rsq = 2.0_wip
  
            do while (rsq >= 1.0_wip)
              v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
            grandom(2) = rsq * v1
            rsq = 2.0_wip
  
            do while (rsq >= 1.0_wip)
              v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
            grandom(3) = rsq * v1
            random_f(1:3,jx,j) = sigma*grandom(1:3)
  
          end do
        end do

        ! Thermostat
        !
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * scale_v
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do

      end if

    end if

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) &
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               vel, viri_const)

    return

  end subroutine langevin_thermostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv1
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv1(dynamics, istep, dt_long, dt_short,  &
                                   ensemble, domain, constraints,       &
                                   boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    real(wip),                intent(in)    :: dt_long
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: inv_dt, temp0, press0, d_ndegf
    real(dp)                 :: kin(1:3), kin_temp(1:3), ekin
    real(wip)                :: delta_vel(1:3)
    real(wip)                :: volume, press(1:3)
    real(wip)                :: pressxy, pressxyz
    real(wip)                :: factor
    real(wip)                :: bmoment_ref(3), scale_b(1:3)
    real(wip)                :: gamma_t, gamma_p
    real(wip)                :: sigma
    real(wip)                :: v1, v2, rsq, grandom(1:3)
    real(wip)                :: half_dt_short
    real(wip)                :: dt_baro, half_dt_baro, quart_dt_baro
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: size_scale(1:3), vel_scale
    real(dp)                 :: virial_constraint(3,3), virial_sum(3)
    integer                  :: i, j, k, l, jx, i_ndegf, maxiter
    integer                  :: ncell, nboundary

    real(wip),       pointer :: pmass, pforce(:), random_f(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: temporary(:,:,:), temporary1(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp) ,       pointer :: virial_short(:,:), virial_long(:,:)
    real(dp) ,       pointer :: viri_const(:,:), kin_ref(:)
    real(wip),       pointer :: bmoment(:)
    integer,         pointer :: natom(:)


    inv_dt       =  1.0_wip/dt_short
    temp0        =  ensemble%temperature
    press0       =  ensemble%pressure * ATMOS_P
    gamma_t      =  ensemble%gamma_t * AKMA_PS
    gamma_p      =  ensemble%gamma_p * AKMA_PS
    pmass        => ensemble%pmass
    pforce       => ensemble%pforce
    random_f     => ensemble%random_force
    i_ndegf      =  domain%num_deg_freedom
    d_ndegf      =  real(i_ndegf,wip)
    ncell        =  domain%num_cell_local
    nboundary    =  domain%num_cell_boundary

    natom        => domain%num_atom
    mass         => domain%mass
    inv_mass     => domain%inv_mass
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    force_long   => domain%force_long
    force_short  => domain%force_short
    temporary    => domain%coord_old
    temporary1   => domain%velocity_full
    virial_short => dynvars%virial
    virial_long  => dynvars%virial_long
    viri_const   => dynvars%virial_const
    bmoment      => dynvars%barostat_momentum
    kin_ref      => dynvars%kin_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt_short = dt_short / 2.0_wip
    dt_therm      = dt_short * real(dynamics%thermo_period,wip)
    half_dt_therm = dt_therm / 2.0_wip
    dt_baro       = dt_short * real(dynamics%baro_period,wip)
    half_dt_baro  = dt_baro / 2.0_wip
    quart_dt_baro = half_dt_baro / 2.0_wip

    ! scale factor for veloctiy rescaling
    !
    vel_scale = exp(-gamma_t*half_dt_therm)

    ! maximum iteration
    !
    if (mod(istep-1,dynamics%baro_period) == 0) then
      if (constraints%rigid_bond) then
        maxiter = 4
      else
        maxiter = 1
      end if
    else
      maxiter = 1
    end if

    ! constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wip

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! pmass and stochastic force (pforce)
    !
    if (istep == 1) then

      if (ensemble%isotropy == IsotropyISO) then

        pmass = (d_ndegf+3.0_wip)*KBOLTZ*temp0 / (2.0_wip*PI*gamma_p)**2
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = pforce(1)
        pforce(3) = pforce(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        pmass = (d_ndegf+3.0_wip)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = pforce(1)
        pforce(3) = sigma * real(random_get_gauss(),wip)

      else if (ensemble%isotropy == IsotropyANISO) then

        pmass = (d_ndegf+3.0_wip)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = sigma * real(random_get_gauss(),wip)
        pforce(3) = sigma * real(random_get_gauss(),wip)

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        pmass = (d_ndegf+1.0_wip)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = 0.0_wip
        pforce(2) = 0.0_wip
        pforce(3) = sigma * real(random_get_gauss(),wip)

      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(pforce, 3, mpi_wip_real, 0, mpi_comm_country, ierror)
#endif

      ! random force
      !
      factor = 1.0_wip - vel_scale*vel_scale
      factor = factor*KBOLTZ*temp0/2.0_wip

      do j = 1, ncell
        do jx = 1, natom(j)

          sigma = sqrt(factor*inv_mass(jx,j))

          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_wip
          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

      kin_ref(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)
          kin_ref(1:3) = kin_ref(1:3) &
                       + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_ref, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)

    end if

    do j = 1, ncell
      do jx = 1, natom(j)
        vel_ref(1:3,jx,j) = vel(1:3,jx,j)
        temporary1(1:3,jx,j) = 0.0_wip
      end do
    end do

    if (mod(istep-1,dynamics%baro_period) == 0) then

      do i = 1, maxiter

        ! Barostat 1
        !
        kin_temp(1:3) = 0.0_dp
        do j = 1, ncell
          do jx = 1, natom(j)
            kin_temp(1:3) = kin_temp(1:3)  &
                          + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
          end do
        end do

        ! virial + virial_constraint
        !
        virial_sum(1) = virial_short(1,1) + virial_long(1,1)
        virial_sum(2) = virial_short(2,2) + virial_long(2,2)
        virial_sum(3) = virial_short(3,3) + virial_long(3,3)
        virial_sum(1) = virial_sum(1) &
                      + virial_constraint(1,1)
        virial_sum(2) = virial_sum(2) &
                      + virial_constraint(2,2)
        virial_sum(3) = virial_sum(3) &
                      + virial_constraint(3,3)

        call reduce_pres(kin_temp, ekin, virial_sum)
        kin(1:3) = 0.5_dp*(kin_ref(1:3)+kin_temp(1:3))
        ekin = 0.5_dp*(kin(1)+kin(2)+kin(3))

        press(1:3) = real(kin(1:3) + virial_sum(1:3),wip)/volume
        pressxyz = (press(1) + press(2) + press(3))/3.0_wip
        pressxy  = (press(1) + press(2))/2.0_wip

        ! update barostat
        !
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel_ref(1:3,jx,j)
          end do
        end do
        call update_barostat(ensemble, boundary, bmoment_ref, pforce,   &
                             press, pressxyz, pressxy, press0, volume,  &
                             d_ndegf, pmass, gamma_p, ekin,             &
                             dt_baro, half_dt_baro, natom, ncell, vel,  &
                             bmoment)
        size_scale(1:3) = exp(bmoment(1:3)*dt_baro)

        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do

        do j = 1, ncell
          do jx = 1, natom(j)
            factor = half_dt_short * inv_mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*temporary1(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel(1:3,jx,j)*dt_short
          end do
        end do

        ! RATTLE VV1
        !
        if (constraints%rigid_bond) then
          do j = 1, ncell
            do jx = 1, natom(j)
              temporary(1:3,jx,j) = coord(1:3,jx,j)
            end do
          end do
          call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                   coord_ref, domain, constraints, coord, &
                                   vel, viri_const)

          virial_constraint(1:3,1:3) = virial_constraint(1:3,1:3) &
                                     + 2.0_dp * viri_const(1:3,1:3)

          do j = 1, ncell
            do jx = 1, natom(j)
              delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
              vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
              temporary1(1:3,jx,j) = temporary1(1:3,jx,j) &
                                   + mass(jx,j)*delta_vel(1:3)/half_dt_short
            end do
          end do
        end if

      end do

      ! thermostat
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j) * vel_scale
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do

      ! VV1 (long range force)
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          factor = 0.5_wip * dt_long * inv_mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
      end do

      do j = 1, ncell
        do jx = 1, natom(j)
          factor = half_dt_short * inv_mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

    else

      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      ! Thermostat
      !
      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do
      end if

      do j = 1, ncell
        do jx = 1, natom(j)
          factor = half_dt_short * inv_mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

      kin_temp(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)
          kin_temp(1:3) = kin_temp(1:3) + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_temp, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)

      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      ! Thermostat
      !
      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do
      end if

      ! VV1 (long range force)
      !
      if (mod(istep-1,dynamics%elec_long_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            factor = 0.5_wip * dt_long * inv_mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do

      end if

      do j = 1, ncell
        do jx = 1, natom(j)
          factor = half_dt_short * inv_mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

    end if

    kin_ref(1:3) = kin_temp(1:3)

    if (mod(istep-1,dynamics%baro_period) == 0) then

      viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

      ! compute box size(t+2dt)
      !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
      !
      scale_b(1:3) = exp(bmoment(1:3)*dt_baro)
      boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
      boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
      boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

      call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                         boundary%box_size_z)

      boundary%cell_size_x = boundary%box_size_x /real(boundary%num_cells_x,wip)
      boundary%cell_size_y = boundary%box_size_y /real(boundary%num_cells_y,wip)
      boundary%cell_size_z = boundary%box_size_z /real(boundary%num_cells_z,wip)

      ! update boudary conditions
      !
      dynvars%barostat_momentum(1:3) = bmoment(1:3)
      do j = 1, ncell+nboundary
        do jx = 1, natom(j)
          domain%trans_vec(1:3,jx,j) = domain%trans_vec(1:3,jx,j) * scale_b(1:3)
        end do
      end do

      domain%system_size(1) = boundary%box_size_x
      domain%system_size(2) = boundary%box_size_y
      domain%system_size(3) = boundary%box_size_z
    
    end if

    return

  end subroutine langevin_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv2
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv2(dynamics, istep, inner_step, dt_long,    &
                                   dt_short, ensemble, domain, constraints, &
                                   boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(wip),                intent(in)    :: dt_long
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: inv_dt, temp0, press0, d_ndegf
    real(wip)                :: volume
    real(wip)                :: factor
    real(wip)                :: bmoment_ref(3)
    real(wip)                :: gamma_t, gamma_p
    real(wip)                :: sigma, virial_constraint(1:3,1:3)
    real(wip)                :: v1, v2, rsq, grandom(1:3)
    real(wip)                :: half_dt_baro, quart_dt_baro
    real(wip)                :: half_dt_short
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: vel_scale
    real(wip)                :: vel_change(1:3)
    integer                  :: j, k, l, jx, i_ndegf, ncell

    real(wip),       pointer :: pmass, pforce(:), random_f(:,:,:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), force_add(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: coord_deri(:,:,:)
    real(dp) ,       pointer :: virial_short(:,:), virial_long(:,:)
    real(dp) ,       pointer :: viri_const(:,:)
    real(wip),       pointer :: bmoment(:)
    integer,         pointer :: natom(:)


    inv_dt       =  1.0_wip/dt_short
    temp0        =  ensemble%temperature
    press0       =  ensemble%pressure * ATMOS_P
    gamma_t      =  ensemble%gamma_t * AKMA_PS
    gamma_p      =  ensemble%gamma_p * AKMA_PS
    pmass        => ensemble%pmass
    pforce       => ensemble%pforce
    random_f     => ensemble%random_force
    i_ndegf      =  domain%num_deg_freedom
    d_ndegf      =  real(i_ndegf,wip)
    ncell        =  domain%num_cell_local

    mass         => domain%mass
    inv_mass     => domain%inv_mass
    natom        => domain%num_atom
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    coord_deri   => domain%velocity_full
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    force_long   => domain%force_long
    force_short  => domain%force_short
    force_add    => domain%coord_old
    virial_short => dynvars%virial
    virial_long  => dynvars%virial_long
    viri_const   => dynvars%virial_const
    bmoment      => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt_short = dt_short / 2.0_wip
    dt_therm      = dt_short * real(dynamics%thermo_period,wip)
    half_dt_therm = dt_therm / 2.0_wip
    half_dt_baro  = dt_short * real(dynamics%baro_period,wip) / 2.0_wip
    quart_dt_baro = half_dt_baro / 2.0_wip

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wip

    ! initial constraint force
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        force_add(1:3,jx,j) = 0.0_wip
      end do
    end do

    ! Langevin piston force
    !
    if (mod(istep, dynamics%baro_period) == 0) then

      ! barostate coefficient
      !
      bmoment_ref(1:3) = bmoment(1:3)

      ! pmass and stochastic force (pforce)
      !
      if (ensemble%isotropy == IsotropyISO) then

        sigma = sqrt(2.0_wip*gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = pforce(1)
        pforce(3) = pforce(1)
  
      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = pforce(1)
        pforce(3) = sigma * real(random_get_gauss(),wip)

      else if (ensemble%isotropy == IsotropyANISO) then

        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * real(random_get_gauss(),wip)
        pforce(2) = sigma * real(random_get_gauss(),wip)
        pforce(3) = sigma * real(random_get_gauss(),wip)
  
      else if (ensemble%isotropy == IsotropyXY_Fixed) then
  
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = 0.0_wip
        pforce(2) = 0.0_wip
        pforce(3) = sigma * real(random_get_gauss(),wip)

      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(pforce, 3, mpi_wip_real, 0, mpi_comm_country, ierror)
#endif
    end if

    ! Langevin random force
    !
    if (mod(istep, dynamics%thermo_period) == 0) then

      vel_scale = exp(-gamma_t*half_dt_therm)
      factor    = 1.0_wip - vel_scale*vel_scale
      factor    = factor*KBOLTZ*temp0/2.0_wip
      do j = 1, ncell
        do jx = 1, natom(j)

          sigma = sqrt(factor*inv_mass(jx,j))

          rsq = 2.0_wip

          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_wip
  
          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do
  
          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_wip
  
          do while (rsq >= 1.0_wip)
            v1  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            v2  = 2.0_wip*real(random_get(),wip) - 1.0_wip
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

    end if

    ! VV2 (short range force)
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        factor = half_dt_short * inv_mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_add(1:3,jx,j)
      end do
    end do

    ! VV2 (long range force)
    !
    if (inner_step == dynamics%elec_long_period) then
      do j = 1, ncell
        do jx = 1, natom(j)
          factor = 0.5_wip * dt_long * inv_mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
      end do
    end if

    ! Thermostat
    !
    if (mod(istep, dynamics%thermo_period) == 0) then
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do
    end if

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      do j = 1, ncell
        do jx = 1, natom(j)
          coord_deri(1:3,jx,j) = vel(1:3,jx,j) &
                               + bmoment(1:3)*coord(1:3,jx,j)
          coord_ref(1:3,jx,j) = coord_deri(1:3,jx,j)
        end do
      end do
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               coord_deri, viri_const)

      if (mod(istep,dynamics%baro_period) == 0) then

        do j = 1, ncell
          do jx = 1, natom(j)
            vel_change(1:3) = coord_deri(1:3,jx,j) - coord_ref(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + vel_change(1:3)
            force_add(1:3,jx,j) = force_add(1:3,jx,j)                &
                                + mass(jx,j)*vel_change(1:3)/half_dt_short
          end do
        end do

        ! constraint virial
        ! 
        do j = 1, ncell
          do jx = 1, natom(j)
            virial_constraint(1,1) = virial_constraint(1,1)             + &
                                     force_add(1,jx,j)*coord(1,jx,j)
            virial_constraint(2,2) = virial_constraint(2,2)             + &
                                     force_add(2,jx,j)*coord(2,jx,j)
            virial_constraint(3,3) = virial_constraint(3,3)             + &
                                     force_add(3,jx,j)*coord(3,jx,j)
          end do
        end do
      end if

    end if

    if (mod(istep, dynamics%baro_period) == 0) &
      viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    return

  end subroutine langevin_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv1
  !> @brief        Bussi thermostat and barostat
  !! @authors      TA, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    istart      : dynamics start step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv1(dynamics, istep, istart, dt_long, dt_short, &
                              ensemble, domain, constraints,      &
                              boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: istart
    real(wip),                intent(in)    :: dt_long
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    logical                  :: calc_thermostat, calc_barostat
    logical                  :: calc_elec_long
    real(wip)                :: dt, inv_dt, press0, d_ndegf
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: dt_baro, half_dt_baro, quart_dt
    real(wip)                :: volume, press(3)
    real(wip)                :: pressxy, pressxyz
    real(wip)                :: factor
    real(wip)                :: bmoment_ref(3)
    real(wip)                :: tau_t, tau_p
    real(wip)                :: gr, random_gr, scale_vel
    real(wip)                :: half_dt, half_dt_long
    real(wip)                :: size_scale(3), scale_b(3)
    real(wip)                :: vel_scale(3)
    real(dp)                 :: virial_constraint(3,3), virial_sum(3)
    real(dp)                 :: kin_full_ref(3), ekin_full_ref
    integer                  :: i, j, k, l, jx, i_ndegf
    integer                  :: ncell, nboundary, num_degree 
    integer                  :: omp_get_thread_num, id

    real(wip),       pointer :: pmass, ekin0, temp0, degree
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: temporary(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), vel_half(:,:,:)
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wp),        pointer :: trans_vec(:,:,:)
    real(dp) ,       pointer :: virial(:,:), viri_const(:,:), viri_group(:,:)
    real(dp) ,       pointer :: virial_long(:,:)
    real(dp),        pointer :: kin_full(:), kin_half(:), kin_ref(:), kin(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin
    real(wip),       pointer :: bmoment(:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)


    dt          =  dt_short
    inv_dt      =  1.0_wip/dt
    temp0       => ensemble%temperature
    ekin0       => ensemble%kinetic
    degree      => ensemble%degree
    pmass       => ensemble%pmass
    press0      =  ensemble%pressure * ATMOS_P
    tau_t       =  ensemble%tau_t / AKMA_PS
    tau_p       =  ensemble%tau_p / AKMA_PS
    i_ndegf     =  domain%num_deg_freedom
    if (ensemble%group_tp) i_ndegf = domain%num_group_freedom
    d_ndegf     =  real(i_ndegf,wip)
    ncell       =  domain%num_cell_local
    nboundary   =  domain%num_cell_boundary
    natom       => domain%num_atom
    nwater      => domain%num_water
    water_list  => domain%water_list
    mass        => domain%mass
    inv_mass    => domain%inv_mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_half    => domain%velocity_half
    force       => domain%force     
    force_long  => domain%force_long
    force_short => domain%force_short
    temporary   => domain%velocity_full
    trans_vec   => domain%trans_vec
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    viri_group  => dynvars%virial_group
    virial_long => dynvars%virial_long
    kin_full   => dynvars%kin_full
    kin_half   => dynvars%kin_half
    kin_ref    => dynvars%kin_ref
    kin        => dynvars%kin
    ekin_full  => dynvars%ekin_full
    ekin_half  => dynvars%ekin_half
    ekin_ref   => dynvars%ekin_ref
    ekin       => dynvars%ekin
    bmoment     => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt       = 0.5_wip * dt
    half_dt_long  = 0.5_wip * dt_long 
    dt_therm      = dt_short * real(dynamics%thermo_period,wip)
    half_dt_therm = 0.5_wip * dt_therm 
    dt_baro       = dt_short * real(dynamics%baro_period,wip)
    half_dt_baro  = 0.5_wip * dt_baro

    calc_elec_long  = mod(istep-1, dynamics%elec_long_period) == 0
    calc_thermostat = mod(istep-1, dynamics%thermo_period) == 0 .and. &
                      istep > istart
    calc_barostat   = mod(istep-1, dynamics%baro_period) == 0 .and. &
                      istep > istart

    ! pmass and ekin0
    !
    if (istep == 1) then
      degree = d_ndegf + 3.0_wip
      if (ensemble%isotropy == IsotropyXY_Fixed) degree = d_ndegf + 1.0_wip
      ekin0 = 0.5_wip*KBOLTZ*temp0 * degree
      pmass  = degree*KBOLTZ*temp0 * tau_p*tau_p
    end if

    ! refereence
    !
    bmoment_ref(1:3) = bmoment(1:3)
    call copy_coord_vel(ncell, natom, coord, vel, coord_ref, vel_ref)

    ! From here, last part of VV2
    !

    ! scale factor from barostat momentum
    !
    if (calc_elec_long .and. istep > 1) then

      gr = bmoment(1)+bmoment(2)+bmoment(3)
      gr = gr/degree
      scale_b(1:3) = bmoment(1:3) + gr
      vel_scale(1:3) = exp(-scale_b(1:3)*half_dt_long)

      if (ensemble%group_tp) then
        call update_vel_group_3d(constraints, ncell, nwater, water_list, &
                                 vel_scale, mass, vel)
      else
        !$omp parallel do private(j,jx)
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1,jx,j) = vel(1,jx,j) * vel_scale(1)
            vel(2,jx,j) = vel(2,jx,j) * vel_scale(2)
            vel(3,jx,j) = vel(3,jx,j) * vel_scale(3)
          end do
        end do
        !$omp end parallel do
      end if

    end if

    if (calc_thermostat) then

      !$omp parallel do private(j, jx, factor)
      do j = 1, ncell
        do jx = 1, natom(j)
          factor = half_dt * inv_mass(jx,j)
          force(1,jx,j) = force_long(1,jx,j)+force_short(1,jx,j)
          force(2,jx,j) = force_long(2,jx,j)+force_short(2,jx,j)
          force(3,jx,j) = force_long(3,jx,j)+force_short(3,jx,j)
          vel_half(1,jx,j) = factor*force(1,jx,j)
          vel_half(2,jx,j) = factor*force(2,jx,j)
          vel_half(3,jx,j) = factor*force(3,jx,j)
        end do
      end do
      !$omp end parallel do

      if (ensemble%group_tp) then
        call compute_kin_group(constraints, ncell, nwater, water_list, mass, &
                               vel_half, kin_half, ekin_half)
        call compute_kin_group(constraints, ncell, nwater, water_list, mass, &
                               vel, kin_full, ekin_full)
      else
        call calc_kinetic(ncell, natom, mass, vel_half, kin_half, ekin_half)
        call calc_kinetic(ncell, natom, mass, vel     , kin_full, ekin_full)
      end if

      ekin = ekin_full + 2.0_dp*ekin_half/3.0_dp
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))
   
      ! calculate scaling factor
      !
      if (ensemble%tpcontrol == TpcontrolBussi) then
        random_gr = real(random_get_gauss(),wip)
        call vel_scale_bussi(int(degree), half_dt_therm, tau_t, temp0, ekin, &
                             random_gr, scale_vel) 
      else if (ensemble%tpcontrol == TpcontrolBerendsen) then
        call vel_scale_berendsen(int(degree), half_dt_therm, tau_t, temp0,   &
                                 ekin, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolNHC) then
        call vel_scale_nhc(int(degree), half_dt_therm, tau_t, temp0, ekin,   &
                           ensemble, dynvars, scale_vel)
      end if

      if (ensemble%group_tp) then
        call update_vel_group(constraints, ncell, nwater, water_list, &
                              scale_vel, mass, vel)
      else
        !$omp parallel do private(j,jx)
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1,jx,j) = vel(1,jx,j)*scale_vel
            vel(2,jx,j) = vel(2,jx,j)*scale_vel
            vel(3,jx,j) = vel(3,jx,j)*scale_vel
          end do
        end do
        !$omp end parallel do
      end if

      bmoment(1:3) = bmoment(1:3) * scale_vel

      kin_full(1:3) = kin_full(1:3)*scale_vel*scale_vel
      kin(1:3) = kin_full(1:3) + kin_half(1:3)
      ekin = 0.5_dp*(kin_full(1)+kin_full(2)+kin_full(3)) + 2.0_dp*ekin_half/3.0_dp

    end if

    ! update barostat of last VV2 and first VV1 (combination)
    !
    if (calc_barostat) then

      virial(1:3,1:3) = virial(1:3,1:3) + virial_long(1:3,1:3)
      if (ensemble%group_tp) then
        viri_group(1:3,1:3) = 0.0_dp
        call compute_virial_group(constraints, ncell, nwater, water_list, &
                                  mass, coord_ref, force, viri_group)
        virial_sum(1) = virial(1,1) + viri_group(1,1)
        virial_sum(2) = virial(2,2) + viri_group(2,2)
        virial_sum(3) = virial(3,3) + viri_group(3,3)
        virial(1,1) = virial_sum(1)
        virial(2,2) = virial_sum(2)
        virial(3,3) = virial_sum(3)
      else
        virial_sum(1) = virial(1,1)
        virial_sum(2) = virial(2,2)
        virial_sum(3) = virial(3,3)
      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, virial_sum, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)

#endif

      press(1:3) = (kin(1:3) + virial_sum(1:3))/volume
      pressxyz = (press(1)+press(2)+press(3))/3.0_dp
      pressxy  = (press(1)+press(2))/2.0_dp

      ! update barostat
      !
      call update_barostat_mtk(ensemble, press(1), press(2), press(3), &
                               pressxyz, pressxy, press0, volume,      &
                               d_ndegf, pmass, dt_baro, ekin, bmoment)

    end if

    ! from here, only VV1
    !
    if (calc_thermostat) then

      ! calculate scaling factor
      !
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))
      if (ensemble%tpcontrol == TpcontrolBussi) then
        random_gr = real(random_get_gauss(),wip)
        call vel_scale_bussi(int(degree), half_dt_therm, tau_t, temp0, ekin, &
                             random_gr, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolBerendsen) then
        call vel_scale_berendsen(int(degree), half_dt_therm, tau_t, temp0,   &
                                 ekin, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolNHC) then
        call vel_scale_nhc(int(degree), half_dt_therm, tau_t, temp0, ekin,   &
                           ensemble, dynvars, scale_vel)
      end if

      if (ensemble%group_tp) then
        call update_vel_group(constraints, ncell, nwater, water_list, &
                              scale_vel, mass, vel)
      else
        !$omp parallel do private(j, jx)
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1,jx,j) = vel(1,jx,j)*scale_vel
            vel(2,jx,j) = vel(2,jx,j)*scale_vel
            vel(3,jx,j) = vel(3,jx,j)*scale_vel
          end do
        end do
        !$omp end parallel do
      end if

      bmoment(1:3) = bmoment(1:3) * scale_vel

      kin_full(1:3) = kin_full(1:3)*scale_vel*scale_vel
      ekin = 0.5_dp * (kin_full(1)+kin_full(2)+kin_full(3)) + 2.0_dp*ekin_half/3.0_dp
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))

    end if
 
    ! VV1 (long range force)
    !
    if (calc_elec_long) then

      ! scale factor from barostat momentum
      !
      gr = bmoment(1)+bmoment(2)+bmoment(3)
      scale_b(1:3) = bmoment(1:3) + gr/degree
      vel_scale(1:3) = exp(-scale_b(1:3)*half_dt_long)

      ! scale velocity
      !
      if (ensemble%group_tp) then
        call update_vel_group_3d(constraints, ncell, nwater, water_list, &
                              vel_scale, mass, vel)
      else
        !$omp parallel do private(j, jx)
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1,jx,j) = vel(1,jx,j) * vel_scale(1)
            vel(2,jx,j) = vel(2,jx,j) * vel_scale(2)
            vel(3,jx,j) = vel(3,jx,j) * vel_scale(3)
          end do
        end do
        !$omp end parallel do
      end if

      !$omp parallel private(j, jx, factor, id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          factor = half_dt_long * inv_mass(jx,j)
          vel(1,jx,j) = vel(1,jx,j) + factor*force_long(1,jx,j)
          vel(2,jx,j) = vel(2,jx,j) + factor*force_long(2,jx,j)
          vel(3,jx,j) = vel(3,jx,j) + factor*force_long(3,jx,j)
        end do
      end do
      !$omp end parallel 

    end if

    ! VV1
    !
    !$omp parallel private(j, jx, factor, id) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
    do j = id+1, ncell, nthread
      do jx = 1, natom(j)
        factor = half_dt * inv_mass(jx,j)
        vel(1,jx,j) = vel(1,jx,j) + factor*force_short(1,jx,j)
        vel(2,jx,j) = vel(2,jx,j) + factor*force_short(2,jx,j)
        vel(3,jx,j) = vel(3,jx,j) + factor*force_short(3,jx,j)
      end do
    end do
    !$omp end parallel
    size_scale(1:3)  = exp(bmoment(1:3)*half_dt)
    if (ensemble%group_tp) then
      call compute_vv1_coord_group(constraints, ncell, natom, nwater,       &
                                   water_list, mass, inv_mass, force_short, &
                                   coord_ref, vel, size_scale, dt, coord)
    else
      !$omp parallel do private(j, jx)
      do j = 1, ncell
        do jx = 1, natom(j)
          coord(1,jx,j) = size_scale(1)*coord_ref(1,jx,j) + vel(1,jx,j)*dt
          coord(2,jx,j) = size_scale(2)*coord_ref(2,jx,j) + vel(2,jx,j)*dt
          coord(3,jx,j) = size_scale(3)*coord_ref(3,jx,j) + vel(3,jx,j)*dt
          coord(1,jx,j) = size_scale(1)*coord(1,jx,j)
          coord(2,jx,j) = size_scale(2)*coord(2,jx,j)
          coord(3,jx,j) = size_scale(3)*coord(3,jx,j)
        end do
      end do
      !$omp end parallel do
    end if

    ! constraints
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt_short, &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
    end if

    ! update barostat momentum
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)

    ! compute box size
    !
    scale_b(1:3) = exp(bmoment(1:3)*dt)
    boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
    boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
    boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

    call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                       boundary%box_size_z)

    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,wip)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,wip)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,wip)

    ! update boudary conditions
    !
    !$omp parallel private(j,jx,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
    do j = id+1, ncell+nboundary, nthread
      do jx = 1, natom(j)
        trans_vec(1,jx,j) = trans_vec(1,jx,j) * scale_b(1)
        trans_vec(2,jx,j) = trans_vec(2,jx,j) * scale_b(2)
        trans_vec(3,jx,j) = trans_vec(3,jx,j) * scale_b(3)
      end do
    end do
    !$omp end parallel

    domain%system_size(1) = boundary%box_size_x
    domain%system_size(2) = boundary%box_size_y
    domain%system_size(3) = boundary%box_size_z

    return

  end subroutine mtk_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv2
  !> @brief        VV2 with MTK barostat
  !! @authors      TA, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv2(dynamics, istep, dt_long, dt_short,    &
                              ensemble, domain, constraints,         &
                              boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    real(wip),                intent(in)    :: dt_long
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, inv_dt, temp0, press0, d_ndegf, degree
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: dt_baro, half_dt_baro
    real(dp)                 :: cm(1:8)
    real(wip)                :: volume
    real(wip)                :: factor
    real(wip)                :: tau_t, tau_p
    real(wip)                :: ekin0
    real(dp)                 :: virial_constraint(1:3,1:3)
    real(wip)                :: half_dt, half_dt_long
    real(wip)                :: vel_change(1:3)
    integer                  :: j, k, l, jx, i_ndegf, ncell
    integer                  :: omp_get_thread_num, id

    real(wip),       pointer :: pmass, pforce(:)
    real(wip),       pointer :: mass(:,:), inv_mass(:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(wip),       pointer :: force_add(:,:,:), coord_deri(:,:,:)
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp) ,       pointer :: virial(:,:), viri_const(:,:), virial_long(:,:)
    real(wip),       pointer :: bmoment(:)
    integer,         pointer :: natom(:)


    dt          =  dt_short
    inv_dt      =  1.0_wip/dt
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    tau_t       =  ensemble%tau_t / AKMA_PS
    tau_p       =  ensemble%tau_p / AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  domain%num_deg_freedom
    d_ndegf     =  real(i_ndegf,wip)
    ncell       =  domain%num_cell_local
    mass        => domain%mass
    inv_mass    => domain%inv_mass
    natom       => domain%num_atom
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_short => domain%force_short
    force_long  => domain%force_long
    force_add   => domain%coord_old
    coord_deri  => domain%velocity_full
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    virial_long => dynvars%virial_long
    bmoment     => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt = 0.5_wip * dt
    half_dt_long = 0.5_wip * dt_long
    dt_therm      = dt_short * real(dynamics%thermo_period,wip)
    half_dt_therm = 0.5_wip * dt_therm

    ! shift velocity, force
    !
    if (mod(istep, dynamics%elec_long_period) == 0) then

      cm(1:8) = 0.0_dp
      !$omp parallel private(j,jx,id) reduction(+:cm)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          cm(1)  = cm(1) + mass(jx,j)*vel(1,jx,j)
          cm(2)  = cm(2) + mass(jx,j)*vel(2,jx,j)
          cm(3)  = cm(3) + mass(jx,j)*vel(3,jx,j)
          cm(4)  = cm(4) + mass(jx,j)
          cm(5)  = cm(5) + force_long(1,jx,j)
          cm(6)  = cm(6) + force_long(2,jx,j)
          cm(7)  = cm(7) + force_long(3,jx,j)
          cm(8)    = cm(8)   + 1.0_dp
        end do
      end do
      !$omp end parallel 

      call mpi_allreduce(mpi_in_place, cm, 8, mpi_real8, mpi_sum, &
                         mpi_comm_city, ierror)


      !$omp parallel private(j,jx,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          vel(1,jx,j) = vel(1,jx,j) - real(cm(1)/cm(4),wip)
          vel(2,jx,j) = vel(2,jx,j) - real(cm(2)/cm(4),wip)
          vel(3,jx,j) = vel(3,jx,j) - real(cm(3)/cm(4),wip)
          force_long(1,jx,j) = force_long(1,jx,j) - real(cm(5)/cm(8),wip)
          force_long(2,jx,j) = force_long(2,jx,j) - real(cm(6)/cm(8),wip)
          force_long(3,jx,j) = force_long(3,jx,j) - real(cm(7)/cm(8),wip)
        end do
      end do
      !$omp end parallel 

    end if

    ! VV2 (long range force)
    !
    if (mod(istep, dynamics%elec_long_period) == 0) then
      !$omp parallel private(j,jx,factor,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          factor = half_dt_long * inv_mass(jx,j)
          vel(1,jx,j) = vel(1,jx,j) + factor*force_long(1,jx,j)
          vel(2,jx,j) = vel(2,jx,j) + factor*force_long(2,jx,j)
          vel(3,jx,j) = vel(3,jx,j) + factor*force_long(3,jx,j)
        end do
      end do
      !$omp end parallel 
    end if

    ! VV2 (short range force)
    !
    !$omp parallel private(j,jx,factor,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
    do j = id+1, ncell, nthread
      do jx = 1, natom(j)
        factor = half_dt * inv_mass(jx,j)
        vel(1,jx,j) = vel(1,jx,j) + factor*force_short(1,jx,j)
        vel(2,jx,j) = vel(2,jx,j) + factor*force_short(2,jx,j)
        vel(3,jx,j) = vel(3,jx,j) + factor*force_short(3,jx,j)
      end do
    end do
    !$omp end parallel

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then

      !$omp parallel private(j,jx,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          coord_deri(1,jx,j) = vel(1,jx,j) + bmoment(1)*coord(1,jx,j)
          coord_deri(2,jx,j) = vel(2,jx,j) + bmoment(2)*coord(2,jx,j)
          coord_deri(3,jx,j) = vel(3,jx,j) + bmoment(3)*coord(3,jx,j)
          coord_ref(1,jx,j)  = coord_deri(1,jx,j)
          coord_ref(2,jx,j)  = coord_deri(2,jx,j)
          coord_ref(3,jx,j)  = coord_deri(3,jx,j)
        end do
      end do
      !$omp end parallel

      call compute_constraints(ConstraintModeVVER2, .false., dt, coord_ref, &
                               domain, constraints, coord, coord_deri,      &
                               viri_const)

      !$omp parallel private(j,jx,vel_change,id)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do j = id+1, ncell, nthread
        do jx = 1, natom(j)
          vel_change(1) = coord_deri(1,jx,j) - coord_ref(1,jx,j)
          vel_change(2) = coord_deri(2,jx,j) - coord_ref(2,jx,j)
          vel_change(3) = coord_deri(3,jx,j) - coord_ref(3,jx,j)
          vel(1,jx,j) = vel(1,jx,j) + vel_change(1)
          vel(2,jx,j) = vel(2,jx,j) + vel_change(2)
          vel(3,jx,j) = vel(3,jx,j) + vel_change(3)
        end do
      end do
      !$omp end parallel

    end if

    return

  end subroutine mtk_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_barostat
  !> @brief        update barostat parameter bmoment (eta)
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat(ensemble, boundary, bmoment_ref, pforce, &
                             press, pressxyz, pressxy, press0, &
                             volume, d_ndegf, pmass, gamma_p,  &
                             ekin, half_dt, quart_dt, natom,   &
                             ncell, vel, bmoment)

    ! formal arguments
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: bmoment_ref(:)
    real(wip),               intent(in)    :: pforce(:)
    real(dp),                intent(in)    :: ekin
    real(wip),               intent(in)    :: press(:)
    real(wip),               intent(in)    :: pressxyz
    real(wip),               intent(in)    :: pressxy
    real(wip),               intent(in)    :: press0
    real(wip),               intent(in)    :: volume
    real(wip),               intent(in)    :: d_ndegf
    real(wip),               intent(in)    :: pmass
    real(wip),               intent(in)    :: gamma_p
    real(wip),               intent(in)    :: quart_dt
    real(wip),               intent(in)    :: half_dt
    integer,                 intent(in)    :: natom(:)
    integer,                 intent(in)    :: ncell
    real(wip),               intent(inout) :: vel(:,:,:)
    real(wip),               intent(inout) :: bmoment(:)

    ! local variable
    real(wip)        :: gamma0, pressxy0, trace
    real(wip)        :: vel_scale(1:3), ekin_real
    integer          :: i, ix


    gamma0 = ensemble%gamma*ATMOS_P*100.0_wip/1.01325_wip
    ekin_real = real(ekin, wip)

    ! eta(t+1/4dt)
    !
    bmoment(1:3) = exp(-gamma_p*quart_dt/2.0_wip)*bmoment_ref(1:3)

    ! eta(t+1/4dt) is scaled according to pressure
    !
    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + quart_dt*(3.0_wip*volume*(pressxyz - press0) &
                              + 6.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_wip
      bmoment(2) = 0.0_wip
      bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0) &
                              + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass

    end if

    bmoment(1:3)  = exp(-gamma_p*quart_dt/2.0_wip)*bmoment(1:3)

    ! velocities are scaled according to scaled eta(t+1/4dt)
    !
    trace = bmoment(1)+bmoment(2)+bmoment(3)
    vel_scale(1:3) = exp(-half_dt*(bmoment(1:3)+trace/d_ndegf))
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = vel(1:3,ix,i) * vel_scale(1:3)
      end do
    end do

    ! eta(t+1/4dt) is scaled
    !
    bmoment(1:3) = exp(-gamma_p*quart_dt/2.0_wip)*bmoment(1:3)

    ! eta(t+1/2dt)
    !
    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + quart_dt*(3.0_wip*volume*(pressxyz - press0) &
                              + 6.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_wip
      bmoment(2) = 0.0_wip
      bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0) &
                              + 2.0_wip*ekin_real/d_ndegf + pforce(3))/pmass

    end if

    bmoment(1:3)  = exp(-gamma_p*quart_dt/2.0_wip)*bmoment(1:3)

    return

  end subroutine update_barostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_barostat_mtk
  !> @brief        update barostat parameter bmoment for MTK
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat_mtk(ensemble, pressx, pressy, pressz, pressxyz, &
                                 pressxy, press0, volume, d_ndegf, pmass,    &
                                 half_dt, ekin, bmoment)

    ! formal arguments
    type(s_ensemble),        intent(in)    :: ensemble
    real(wip),               intent(in)    :: pressx
    real(wip),               intent(in)    :: pressy
    real(wip),               intent(in)    :: pressz
    real(wip),               intent(in)    :: pressxyz
    real(wip),               intent(in)    :: pressxy
    real(wip),               intent(in)    :: press0
    real(wip),               intent(in)    :: volume
    real(wip),               intent(in)    :: d_ndegf
    real(wip),               intent(in)    :: pmass
    real(wip),               intent(in)    :: half_dt
    real(dp),                intent(in)    :: ekin
    real(wip),               intent(inout) :: bmoment(:)

    real(wip)                :: ekin_real


    ekin_real = real(ekin,wip)

    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + half_dt*(3.0_wip*volume*(pressxyz - press0) &
                              + 6.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      bmoment(1) = bmoment(1) + half_dt*(volume*(pressxy - press0)   &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(2) = bmoment(2) + half_dt*(volume*(pressxy - press0)   &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz  - press0)   &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass

    else if (ensemble%isotropy == IsotropyANISO) then

      bmoment(1) = bmoment(1) + half_dt*(volume*(pressx - press0)    &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(2) = bmoment(2) + half_dt*(volume*(pressy - press0)    &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz - press0)    &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_wip
      bmoment(2) = 0.0_wip
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz - press0) &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass

    end if

    return

  end subroutine update_barostat_mtk

#ifdef HAVE_MPI_GENESIS
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine reduce_pres
  !> @brief
  !! @authors   JJ
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_pres(val1, val2, val3)

    ! formal arguments
    real(dp), intent(inout) :: val1(:)
    real(dp), intent(inout) :: val2
    real(dp), intent(inout) :: val3(:)

    ! local variables
    double precision              :: before_reduce(7), after_reduce(7)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4)   = val2
    before_reduce(5:7) = val3(1:3)

    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_double_precision,   &
                       mpi_sum, mpi_comm_country, ierror)

    val1(1:3) = after_reduce(1:3)
    val2      = after_reduce(4)
    val3(1:3) = after_reduce(5:7)

    return

  end subroutine reduce_pres
#endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_kinetic
  !> @brief        kinetic energy calculation
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_kinetic(ncell, natom, mass, vel, kin, ekin)

    ! formal arguments
    integer,                 intent(inout) :: ncell
    integer,                 intent(inout) :: natom(:)
    real(wip),               intent(inout) :: mass(:,:)
    real(wip),               intent(inout) :: vel (:,:,:)
    real(dp),                intent(inout) :: kin(:)
    real(dp),                intent(inout) :: ekin

    integer                  :: i, ix


    kin(1:3) = 0.0_dp
    !$omp parallel do private(i, ix) reduction(+:kin)
    do i = 1, ncell
      do ix = 1, natom(i)
        kin(1) = kin(1) + mass(ix,i)*vel(1,ix,i)*vel(1,ix,i)
        kin(2) = kin(2) + mass(ix,i)*vel(2,ix,i)*vel(2,ix,i)
        kin(3) = kin(3) + mass(ix,i)*vel(3,ix,i)*vel(3,ix,i)
      end do
    end do
    !$omp end parallel do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, kin, 3, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif
    ekin = kin(1) + kin(2) + kin(3)
    ekin = 0.5_dp * ekin

    return

  end subroutine calc_kinetic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !
  !  Subroutine    vel_scale_bussi
  !> @brief        scaling factor with Bussi's thermostat
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vel_scale_bussi(degree, dt, tau_t, temp0, ekin, rr, scale_vel)

    ! formal arguments
    integer,                 intent(in)    :: degree
    real(wip),               intent(in)    :: dt, tau_t
    real(wip),               intent(in)    :: temp0
    real(dp),                intent(in)    :: ekin
    real(wip),               intent(in)    :: rr
    real(wip),               intent(inout) :: scale_vel

    real(wip)                :: tempf, tempt, factor


    factor = exp(-dt/tau_t)
    tempf = 2.0_wip * real(ekin,wip)/(real(degree,wip)*KBOLTZ)
    tempt = tempf*factor    &
          + temp0/real(degree,wip)*(1.0_wip-factor)   &
            *(sum_gauss(degree-1)+rr*rr)              &
          + 2.0_wip*sqrt(tempf*temp0/real(degree,wip) &
            *(1.0_wip-factor)*factor)*rr
    scale_vel = sqrt(tempt/tempf)

#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(scale_vel, 1, mpi_wip_real, 0, mpi_comm_country, ierror)
#endif

    return

  end subroutine vel_scale_bussi

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_scale_berendsen
  !> @brief        scaling factor with Berendsen's thermostat
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vel_scale_berendsen(degree, dt, tau_t, temp0, ekin, scale_vel)

    ! formal arguments
    integer,                 intent(in)    :: degree
    real(wip),               intent(in)    :: dt, tau_t
    real(wip),               intent(in)    :: temp0
    real(dp),                intent(in)    :: ekin
    real(wip),               intent(inout) :: scale_vel

    real(wip)                :: tempf, tempt, factor


    factor = exp(-dt/tau_t)
    tempf = 2.0_wip * real(ekin,wip)/(real(degree,wip)*KBOLTZ)
    scale_vel = sqrt(1.0_wip + (dt/tau_t)*(temp0/tempf-1.0_wip))

    return

  end subroutine vel_scale_berendsen

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_scale_nhc
  !> @brief        scaling factor with NHC thermostat
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vel_scale_nhc(degree, dt, tau_t, temp0, ekin, ensemble, &
                           dynvars, scale_vel)

    ! formal arguments
    integer,                 intent(in)    :: degree
    real(wip),               intent(in)    :: dt, tau_t
    real(wip),               intent(in)    :: temp0
    real(dp),                intent(in)    :: ekin
    type(s_ensemble), target,intent(inout) :: ensemble
    type(s_dynvars),  target,intent(inout) :: dynvars
    real(wip),               intent(inout) :: scale_vel

    integer                  :: nh_length, nh_step
    integer                  :: i, j, k
    real(dp)                 :: ekf
    real(wip)                :: w(3)
    real(wip)                :: KbT
    real(wip)                :: dt_small, dt_1, dt_2, dt_4, dt_8
    real(wip)                :: tempf, tempt, scale_kin
    real(wip),       pointer :: nh_mass(:), nh_vel(:)
    real(wip),       pointer :: nh_force(:), nh_coef(:)


    nh_length   = ensemble%nhchain
    nh_step     = ensemble%nhmultistep
    KbT         = KBOLTZ * temp0

    nh_mass     => dynvars%nh_mass
    nh_vel      => dynvars%nh_velocity
    nh_force    => dynvars%nh_force
    nh_coef     => dynvars%nh_coef

    ! NH mass
    !
    nh_mass(2:nh_length) = KbT * tau_t*tau_t
    nh_mass(1)           = real(degree,wip) * nh_mass(2)

    ! Yoshida coefficient
    !
    w(1) = 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
    w(3) = w(1)
    w(2) = 1.0_dp - w(1) - w(3)

    dt_small  = dt / real(nh_step, wip)
    scale_vel = 1.0_wip
    ekf = 2.0_dp*ekin

    do i = 1, nh_step
      do j = 1, 3

        dt_1 = w(j) * dt_small
        dt_2 = dt_1 * 0.5_dp
        dt_4 = dt_2 * 0.5_dp
        dt_8 = dt_4 * 0.5_dp

        nh_force(nh_length) = nh_mass(nh_length-1) &
                             *nh_vel(nh_length-1)*nh_vel(nh_length-1)-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

        do k = nh_length-1, 2, -1
          nh_force(k) = nh_mass(k-1)*nh_vel(k-1)*nh_vel(k-1)-KbT
          nh_force(k) = nh_force(k) / nh_mass(k)
          nh_coef(k)  = exp(-nh_vel(k+1)*dt_8)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(1) = real(ekf,wip) - real(degree,wip)*KbT
        nh_force(1) = nh_force(1) / nh_mass(1)
        nh_coef(1)  = exp(-nh_vel(2)*dt_8)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        scale_kin = exp(-nh_vel(1)*dt_1)
        scale_vel = scale_vel * exp(-nh_vel(1)*dt_2)
        ekf = ekf * scale_kin

        nh_force(1) = (ekf - real(degree,wip)*KbT) / nh_mass(1)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        do k = 2, nh_length-1
          nh_force(k) = nh_mass(k-1)*nh_vel(k-1)*nh_vel(k-1)-KbT
          nh_force(k) = nh_force(k) / nh_mass(k)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(nh_length) = nh_mass(nh_length-1) &
                             *nh_vel(nh_length-1)*nh_vel(nh_length-1)-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

      end do
    end do

    return

  end subroutine vel_scale_nhc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_coord_vel
  !> @brief        copy coordinate and velocity
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_coord_vel(ncell, natom, coord_in, vel_in, coord_out, vel_out)

    ! formal arguments
    integer,                 intent(in)    :: ncell
    integer,                 intent(in)    :: natom(:)
    real(wip),               intent(in)    :: coord_in(:,:,:)
    real(wip),               intent(in)    :: vel_in(:,:,:)
    real(wip),               intent(inout) :: coord_out(:,:,:)
    real(wip),               intent(inout) :: vel_out(:,:,:)

    integer                  :: i, ix


    !$omp parallel do
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_out(1,ix,i) = coord_in(1,ix,i)
        coord_out(2,ix,i) = coord_in(2,ix,i)
        coord_out(3,ix,i) = coord_in(3,ix,i)
        vel_out(1,ix,i)   = vel_in(1,ix,i)
        vel_out(2,ix,i)   = vel_in(2,ix,i)
        vel_out(3,ix,i)   = vel_in(3,ix,i)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine copy_coord_vel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine bcast_boxsize
  !> @brief
  !! @authors   JJ
  !! @param[in]
  !! @date    2012/09/10 (JJ)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bcast_boxsize(val1, val2, val3)

    ! formal arguments
    real(wip),               intent(inout) :: val1, val2, val3

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(wip)                :: list(3)

    list(1) = val1
    list(2) = val2
    list(3) = val3

    call mpi_bcast(list, 3, mpi_wip_real, 0, mpi_comm_country, ierror)

    val1 = list(1)
    val2 = list(2)
    val3 = list(3)

#endif

    return

  end subroutine bcast_boxsize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_vverlet 
  !> @brief        change target temperature linearly
  !! @authors      TM, JJ
  !! @param[in]    dynamics: dynamics information
  !! @param[out]   ensemble: ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine simulated_annealing_vverlet(dynamics, ensemble)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble

    ! local variable
    real(wip)                :: old_temperature


    if (.not. dynamics%annealing) return

    old_temperature      = ensemble%temperature
    ensemble%temperature = ensemble%temperature + dynamics%dtemperature

    if (main_rank) then
      write(MsgOut,'(A,F10.3,A,F10.3)')                              &
            'Simulated_Annealing_Leapfrog> Anneal temperature from', &
            old_temperature, ' to ', ensemble%temperature

      if (ensemble%temperature < 0.0_wip) &
        call error_msg( &
        'Simulated_Annealing_Leapfrog> Error: Temperature is less than 0 K')

    end if

    return

  end subroutine simulated_annealing_vverlet

end module sp_md_respa_mod
