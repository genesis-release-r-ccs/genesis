!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_md_leapfrog_mod
!> @brief   Performs molecular dynamics simulations.
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_md_leapfrog_mod

  use cg_output_mod
  use cg_update_domain_mod
  use cg_assign_velocity_mod
  use cg_dynvars_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
  use cg_domain_mod
  use cg_energy_mod
  use cg_output_str_mod
  use cg_dynamics_str_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use random_mod
  use messages_mod
  use fileio_grotop_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: leapfrog_dynamics
  private :: initial_leapfrog
  private :: control_temp_pres_leap
  private :: berendsen_leapfrog
  private :: nose_hoover_leapfrog
  private :: gaussian_leapfrog
  private :: langevin_leapfrog_nvt
  private :: langevin_leapfrog_npt
  private :: simulated_annealing_leapfrog
  private :: reduce_pres
  private :: bcast_boxsize

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    leapfrog_dynamics
  !> @brief        Implements a leapfrog integrator with domain decomposition.
  !! @authors      JJ, TM
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine leapfrog_dynamics(grotop, molecule, output, domain, enefunc, &
                               dynvars, dynamics, pairlist, boundary,     &
                               ensemble, comm)

    ! formal arguments
    type(s_grotop),          intent(in   ) :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm   

    ! local variables
    real(wip)                :: simtim, dt, temperature
    integer                  :: i, nsteps, iseed
    integer                  :: istart, iend
    logical                  :: npt

    temperature    =  ensemble%temperature
    npt            =  ensemble%use_barostat
    nsteps         =  dynamics%nsteps
    istart         =  dynamics%istart_step
    iend           =  dynamics%iend_step
    dt             =  dynamics%timestep/AKMA_PS
    simtim         =  dynamics%initial_time
    iseed          =  dynamics%iseed_init_velocity


    if (abs(dynamics%initial_rmsd) < 0.001_wp)  &
      dynamics%initial_rmsd = real(dynvars%energy%rmsd,wp)
    if (dynamics%target_md) enefunc%rmsd_force = 1.0_wp / real(dt*dt,wp)

    ! Check restart
    !
    if (.not. dynamics%restart) then

      call initial_velocity(temperature, iseed, domain)

      call stop_trans_rotation(dynamics, domain)

      call initial_leapfrog(npt, output, enefunc, dynamics,   &
                            pairlist, boundary, ensemble,     &
                            domain, dynvars, comm)
    end if

#ifdef KCOMP
    ! Start performance check on K computer
    !
!   call fipp_start()
#endif

    ! Main loop 
    !   coord is at 0 +  dt and vel is at 0 + 1/2dt, if restart off
    !   coord is at t + 2dt and vel is at t + 3/2dt, if restart on
    !
    do i = istart, iend

      simtim = simtim + dynamics%timestep
      dynvars%time = simtim
      dynvars%step = i
      if (dynamics%target_md .or. dynamics%steered_md) &
        enefunc%rmsd_target = dynamics%initial_rmsd &
                            + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                             *real(dynvars%step,wp)/real(nsteps,wp)

      call timer(TimerIntegrator, TimerOn)

      ! Save coordinates(t + dt) and velocities(t + 1/2dt)
      !
      call coord_ref_save(domain)

      call timer(TimerIntegrator, TimerOff)

      ! send the coordinate data
      !
      call timer(TimerIntegrator, TimerOn)
      call timer(TimerComm1, TimerOn)

      call communicate_coor(domain, comm)

      call timer(TimerComm1, TimerOff)
      call timer(TimerIntegrator, TimerOff)

      ! Compute energy(t + dt), force(t + dt), and internal virial(t + dt)
      !
      call compute_energy(pairlist, boundary,                &
                          npt, .false.,                      &
                          mod(i,dynamics%eneout_period)==0,  &
                          enefunc, domain, dynvars)

      call timer(TimerIntegrator, TimerOn)
      call timer(TimerComm2, TimerOn)
                          
      call communicate_force(domain, comm)

      call timer(TimerComm2, TimerOff)

      if (ensemble%tpcontrol /= TpcontrolLangevin) then

        ! Newtonian dynamics
        !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
        !   r(t+2dt) = r(t+dt) + dt*v(t+3/2dt)
        !
        call update_coord_vel(dt, domain)

        ! Control temperature and pressure
        !   coord     is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref is at t +  dt, vel_ref is at t + 1/2dt
        !   scale velocities(t + 3/2dt) and coordinates(t + 2dt)
        !

        if (ensemble%ensemble /= EnsembleNVE) then

          call control_temp_pres_leap(dynamics, ensemble, domain, &
                                      boundary, dynvars)

        end if

      else

        ! Langevin dynamics
        !
        if (ensemble%ensemble == EnsembleNVT) then

          call langevin_leapfrog_nvt(dynamics, ensemble, domain, &
                                     dynvars)

        else if (ensemble%ensemble == EnsembleNPT  .or. &
                 ensemble%ensemble == EnsembleNPAT .or. &
                 ensemble%ensemble == EnsembleNPgT) then

          call langevin_leapfrog_npt(dynamics, ensemble, domain, &
                                     boundary, dynvars)

        end if

      end if

      call timer(TimerIntegrator, TimerOff)


      if (dynamics%stoptr_period > 0) then

        if (mod(i,dynamics%stoptr_period) == 0) then

          call coord_mid_save(domain)
          call stop_trans_rotation(dynamics, domain)
        end if

      end if

      ! ------------------------------------------------------------------------
      ! OUTPUT energy, trajectory, and restart file
      !   coord     is at t + 2dt, vel      is at t + 3/2dt
      !   coord_ref is at t +  dt, vel_ref  is at t +    dt
      !   volume    is at t +  dt, box_size is at t +   2dt

      call random_push_stock

      ! output dynvars(t + dt)
      !
      call output_md(output, enefunc, dynamics, boundary, ensemble,  &
                                dynvars, domain)

      if (mod(i,dynamics%eneout_period) == 0) then
        call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                             dynvars)
        call output_dynvars (output, enefunc, dynvars, ensemble, boundary)
      end if

      ! update interaction
      !
      if (dynamics%nbupdate_period > 0) &
        call domain_interaction_update_md(i, grotop, dynamics, molecule,  &
                                          domain, enefunc, pairlist,      &
                                          boundary, comm)

      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then
        if (mod(i,dynamics%anneal_period) == 0) then

          call simulated_annealing_leapfrog(dynamics, ensemble)

        end if
      end if

    end do

#ifdef KCOMP
    ! Start performance check on K computer
    !
!   call fipp_stop()
#endif

    return

  end subroutine leapfrog_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_leapfrog
  !> @brief        Computes the first step (0+dt).
  !! @authors      JJ
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    output      : output information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_leapfrog(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, domain, dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(wip)                :: factor, temperature
    real(wip)                :: imass, simtim, dt, friction
    integer                  :: i, ncell

    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wp),        pointer :: coord_pbc(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:)
    real(wip),       pointer :: force(:,:), force_long(:,:)
    real(wip),       pointer :: mass(:)
    real(wp),        pointer :: force_omp(:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    integer,         pointer :: natom(:)

    natom           => domain%num_atom
    coord           => domain%coord
    coord_ref       => domain%coord_ref
    coord_pbc       => domain%translated
    force           => domain%force
    force_long      => domain%force_long
    force_omp       => domain%force_omp
    force_pbc       => domain%force_pbc
    virial_cell     => domain%virial_cellpair
    vel             => domain%velocity
    vel_ref         => domain%velocity_ref
    mass            => domain%mass
    virial          => dynvars%virial

    ncell           =  domain%num_cell_local
    dt              =  dynamics%timestep/AKMA_PS
    simtim          =  dynamics%initial_time
    temperature     =  ensemble%temperature
    friction        =  ensemble%gamma_t * AKMA_PS

    dynvars%time    = simtim
    dynvars%step    = 0


    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do i = 1, domain%num_atom_domain
      coord_ref(i,1) = coord(i,1)
      coord_ref(i,2) = coord(i,2)
      coord_ref(i,3) = coord(i,3)
      vel_ref(i,1)   = vel(i,1)
      vel_ref(i,2)   = vel(i,2)
      vel_ref(i,3)   = vel(i,3)
    end do

    ! calculate energy(0) and forces(0)
    !
    call communicate_coor(domain, comm)

    call compute_energy(pairlist, boundary,        &
                        npt, .false., .true.,      &
                        enefunc, domain, dynvars)

    ! get forces by communication
    !
    call communicate_force(domain, comm)

    ! Calculate velocities(0 - dt/2)
    !
    do i = 1, domain%num_atom_domain
      imass = 1.0_wp/mass(i)
      vel(i,1) = vel(i,1) - 0.5_wp*dt*force(i,1)*imass
      vel(i,2) = vel(i,2) - 0.5_wp*dt*force(i,2)*imass
      vel(i,3) = vel(i,3) - 0.5_wp*dt*force(i,3)*imass
    end do

    ! first step dynamics
    !   Calculate coordinates(0 + dt) and velocities(0 + 1/2dt)
    !   Note: this routine corresponds to i = 0 in the main loop
    !

    ! coord_ref <= coord(0) and vel_ref <= vel(0 - 1/2dt)
    !
    do i = 1, domain%num_atom_domain
      coord_ref(i,1) = coord(i,1)
      coord_ref(i,2) = coord(i,2)
      coord_ref(i,3) = coord(i,3)
      vel_ref(i,1)   = vel(i,1)
      vel_ref(i,2)   = vel(i,2)
      vel_ref(i,3)   = vel(i,3)
    end do

    ! calculate energy(0) and forces(0)
    !
    call communicate_coor(domain, comm)

    call compute_energy(pairlist, boundary,        &
                        npt, .false., .true.,      &
                        enefunc, domain, dynvars)

    call communicate_force(domain, comm)

    ! Newtonian dynamics
    !   v(0+1/2dt) = v(0-1/2dt) + dt*F(0)/m
    !   r(0+dt) = r(0) + dt*v(0+1/2dt)
    !
    do i = 1, domain%num_atom_domain
      factor = dt/mass(i)
      vel(i,1) = vel(i,1) + factor*force(i,1)
      vel(i,2) = vel(i,2) + factor*force(i,2)
      vel(i,3) = vel(i,3) + factor*force(i,3)
      coord(i,1) = coord(i,1) + dt*vel(i,1)
      coord(i,2) = coord(i,2) + dt*vel(i,2)
      coord(i,3) = coord(i,3) + dt*vel(i,3)
    end do

    ! output dynvars(0)
    !
    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)

    dynamics%restart = .true.

    ! at this point
    !   coord is at 0 + dt, and vel is at 0 + 1/2dt

    return

  end subroutine initial_leapfrog
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_temp_pres_leap
  !> @brief        Controls temperature and pressure.
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_temp_pres_leap(dynamics, ensemble, domain, &
                                    boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars 


    select case (ensemble%tpcontrol)

    case (TpcontrolNoseHoover)

      call nose_hoover_leapfrog(dynamics, ensemble, domain, &
                                dynvars)

    case (TpcontrolBerendsen)

      call berendsen_leapfrog  (dynamics, ensemble, domain, &
                                boundary, dynvars)

    case (TpcontrolGauss)

      call gaussian_leapfrog   (dynamics, ensemble, domain)

    end select

    return

  end subroutine control_temp_pres_leap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nose_hoover_leapfrog
  !> @brief        Controls temperature and pressure using the Nose-Hoover 
  !                method.
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nose_hoover_leapfrog(dynamics, ensemble, domain, &
                                  dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, inv_dt, temp0, press0, tau_t
    real(dp)                 :: kin(1:3)
    real(wip)                :: vel_tmp(1:3), ekin
    real(wip)                :: factor
    real(wip)                :: sigma, qmass
    integer                  :: i, num_degree

    real(wip),       pointer :: fric_ref, fric
    real(wip),       pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:)
    real(dp),        pointer :: virial(:,:)
    real(wip),       pointer :: mass(:)
    integer,         pointer :: ncell, natom(:)


    natom      => domain%num_atom
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    mass       => domain%mass
    ncell      => domain%num_cell_local

    virial     => dynvars%virial
    fric_ref   => dynvars%thermostat_momentum_ref
    fric       => dynvars%thermostat_momentum

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    num_degree =  domain%num_deg_freedom
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure
    tau_t      =  ensemble%tau_t/AKMA_PS


    ! save friction coeff(t+1/2dt)
    !
    dynvars%thermostat_momentum_ref = dynvars%thermostat_momentum

    ! compute kinetic energy(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    kin(1:3) = 0.0_dp

    do i = 1, domain%num_atom_domain
      vel_tmp(1) = 0.5_wip*(vel(i,1) + vel_ref(i,1))
      vel_tmp(2) = 0.5_wip*(vel(i,2) + vel_ref(i,2))
      vel_tmp(3) = 0.5_wip*(vel(i,3) + vel_ref(i,3))
      kin(1) = kin(1) + mass(i)*vel_tmp(1)*vel_tmp(1)
      kin(2) = kin(2) + mass(i)*vel_tmp(2)*vel_tmp(2)
      kin(3) = kin(3) + mass(i)*vel_tmp(3)*vel_tmp(3)
    end do
    ekin   = 0.5_dp * (kin(1)+kin(2)+kin(3))

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! thermostat
    !   fric(t+3/2dt) = fric(t+1/2dt) + dt*[2Ekin(t+dt)-2s]/qmass
    !   fric(t+dt)    = [fric(t+3/2dt) + fric(t+1/2dt)]/2
    !   v(t+dt)    = [v(t+3/2dt) + v(t+1/2dt)]/2
    !   v(t+3/2dt) = v(t+1/2dt) + dt*[F(t+dt)/m - fric(t+dt)*v(t+dt)]
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    sigma = 0.5_wip * num_degree * KBOLTZ * temp0
    qmass = 2.0_wip * sigma * tau_t**2
    fric  = fric_ref + dt*(2.0_wip*ekin - 2.0_wip*sigma)/qmass

    factor = 0.5_wip*(fric_ref + fric)
    do i = 1, domain%num_atom_domain
      vel_tmp(1)   = 0.5_wip*(vel(i,1) + vel_ref(i,1))
      vel_tmp(2)   = 0.5_wip*(vel(i,2) + vel_ref(i,2))
      vel_tmp(3)   = 0.5_wip*(vel(i,3) + vel_ref(i,3))
      vel(i,1)  = vel_ref(i,1)   &
                + (force(i,1)/mass(i) - factor*vel_tmp(1))*dt
      vel(i,2)  = vel_ref(i,2)   &
                + (force(i,2)/mass(i) - factor*vel_tmp(2))*dt
      vel(i,3)  = vel_ref(i,3)   &
                + (force(i,3)/mass(i) - factor*vel_tmp(3))*dt
      coord(i,1) = coord_ref(i,1) + vel(i,1)*dt
      coord(i,2) = coord_ref(i,2) + vel(i,2)*dt
      coord(i,3) = coord_ref(i,3) + vel(i,3)*dt
    end do

    return

  end subroutine nose_hoover_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    berendsen_leapfrog
  !> @brief        Controls temperature and pressure using the Berendsen method.
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine berendsen_leapfrog(dynamics, ensemble, domain, &
                                boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, inv_dt, temp0, press0
    real(dp)                 :: kin(1:3)
    real(wip)                :: tau_t, tau_p, compress
    real(dp)                 :: volume, press(1:3)
    real(wip)                :: vel_tmp(1:3), tpress, tempt
    real(wip)                :: scale_v, scale_c(3), factor
    real(dp)                 :: ekin, virial_sum(1:3)
    integer                  :: i, num_degree

    real(wip),       pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),        pointer :: trans(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:)
    real(wip),       pointer :: mass(:)
    real(dp),        pointer :: virial(:,:)
    integer,         pointer :: ncell, nboundary
    integer,         pointer :: natom(:)


    ncell      => domain%num_cell_local
    nboundary  => domain%num_cell_boundary
    mass       => domain%mass
    natom      => domain%num_atom
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    virial     => dynvars%virial
    trans      => domain%trans_vec

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    num_degree =  domain%num_deg_freedom
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure 
    tau_t      =  ensemble%tau_t/AKMA_PS
    tau_p      =  ensemble%tau_p/AKMA_PS
    compress   =  ensemble%compressibility

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    kin(1:3) = 0.0_dp

    do i = 1, domain%num_atom_domain
      vel_tmp(1) = 0.5_wip*(vel(i,1) + vel_ref(i,1))
      vel_tmp(2) = 0.5_wip*(vel(i,2) + vel_ref(i,2))
      vel_tmp(3) = 0.5_wip*(vel(i,3) + vel_ref(i,3))
      kin(1) = kin(1) + mass(i)*vel_tmp(1)*vel_tmp(1)
      kin(2) = kin(2) + mass(i)*vel_tmp(2)*vel_tmp(2)
      kin(3) = kin(3) + mass(i)*vel_tmp(3)*vel_tmp(3)
    end do

    ekin   = 0.5_dp * (kin(1)+kin(2)+kin(3))

    call reduce_pres(kin, ekin, virial(1,1), virial(2,2), &
                     virial(3,3), virial_sum)
    tempt  = 2.0_wip * ekin / (num_degree * KBOLTZ)


    ! thermostat
    !   scale_v    = sqrt(1+dt/tau_t(temp_ref/temp(t+dt) -1))
    !   v(t+3/2dt) = (v(t+1/2dt) + dtF(t+dt)/m) * scale_v
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    scale_v = sqrt(1.0_wip + dt * (temp0/tempt - 1.0_wip)/tau_t)

    do i = 1, domain%num_atom_domain
      factor   = dt/mass(i)
      vel(i,1) = (vel_ref(i,1)+factor*force(i,1))*scale_v
      vel(i,2) = (vel_ref(i,2)+factor*force(i,2))*scale_v
      vel(i,3) = (vel_ref(i,3)+factor*force(i,3))*scale_v
      coord(i,1) = coord_ref(i,1) + vel(i,1)*dt
      coord(i,2) = coord_ref(i,2) + vel(i,2)*dt
      coord(i,3) = coord_ref(i,3) + vel(i,3)*dt
    end do

    if (ensemble%use_barostat) then

      ! compute pressure(t+dt)
      !
      press(1:3) = P_ATMOS * (kin(1:3) + virial_sum(1:3)) / volume
      tpress = (press(1) + press(2) + press(3))/3.0_wip

      ! barostat
      !   scale coordinates(t+2dt) and box_size(t+2dt)
      !   scale_c  = (1- compress*dt*(Pext - P(t+dt)/tau_p))**1/3
      !   r(t+2dt)    = scale_c * r(t+2dt)
      !   size(t+2dt) = scale_c * size(t+dt)
      !
      if (ensemble%isotropy == IsotropyISO) then

        scale_c(1) = (1.0_wip - compress*dt*(press0 - tpress)/tau_p)**ONE_THIRD
        scale_c(2) = scale_c(1)
        scale_c(3) = scale_c(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        scale_c(1) = (1.0_wip -    &
            compress*dt*(press0 - 0.5_wip*(press(1)+press(2)))/tau_p)**ONE_THIRD
        scale_c(2) = (1.0_wip -    &
            compress*dt*(press0 - 0.5_wip*(press(1)+press(2)))/tau_p)**ONE_THIRD
        scale_c(3) = (1.0_wip -    &
            compress*dt*(press0 - press(3))/tau_p)**ONE_THIRD

      else if (ensemble%isotropy == IsotropyANISO) then

        scale_c(1:3) = (1.0_wip -   &
           compress*dt*(press0 - press(1:3))/tau_p)**ONE_THIRD

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        scale_c(1) = 1.0_wip
        scale_c(2) = 1.0_wip
        scale_c(3) = (1.0_wip -    &
             compress*dt*(press0 - press(3))/tau_p)**ONE_THIRD

      end if

      do i = 1, domain%num_atom_domain
        coord(i,1) = scale_c(1) * coord(i,1)
        coord(i,2) = scale_c(2) * coord(i,2)
        coord(i,3) = scale_c(3) * coord(i,3)
        trans(i,1) = real(scale_c(1),wp) * trans(i,1)
        trans(i,2) = real(scale_c(2),wp) * trans(i,2)
        trans(i,3) = real(scale_c(3),wp) * trans(i,3)
      end do

      do i = domain%num_atom_domain+1, &
             domain%num_atom_domain+domain%num_atom_boundary
        trans(i,1) = real(scale_c(1),wp) * trans(i,1)
        trans(i,2) = real(scale_c(2),wp) * trans(i,2)
        trans(i,3) = real(scale_c(3),wp) * trans(i,3)
      end do

      boundary%box_size_x = scale_c(1) * boundary%box_size_x_ref
      boundary%box_size_y = scale_c(2) * boundary%box_size_y_ref
      boundary%box_size_z = scale_c(3) * boundary%box_size_z_ref

      call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                         boundary%box_size_z)

      boundary%cell_size_x = boundary%box_size_x/real(boundary%num_cells_x,wip)
      boundary%cell_size_y = boundary%box_size_y/real(boundary%num_cells_y,wip)
      boundary%cell_size_z = boundary%box_size_z/real(boundary%num_cells_z,wip)

      domain%system_size(1) = real(boundary%box_size_x,wp)
      domain%system_size(2) = real(boundary%box_size_y,wp)
      domain%system_size(3) = real(boundary%box_size_z,wp)

    end if

    return

  end subroutine berendsen_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gaussian_leapfrog
  !> @brief        Controls temperature using a Gaussian thermostat.
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gaussian_leapfrog(dynamics, ensemble, domain)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    real(wip)                :: dt, inv_dt, tempt, beta
    real(dp)                 :: ekin
    real(wip)                :: factor1, factor2, factor
    real(wip)                :: vel_tmp(1:3), temp0
    integer                  :: i, num_degree

    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    integer,         pointer :: ncell, natom(:)


    ncell      => domain%num_cell_local
    natom      => domain%num_atom
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    mass       => domain%mass

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    temp0      =  ensemble%temperature
    num_degree =  domain%num_deg_freedom


    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    ekin = 0.0_dp
    do i = 1, domain%num_atom_domain
      vel_tmp(1) = 0.5_wip*(vel(i,1)+vel_ref(i,1))
      vel_tmp(2) = 0.5_wip*(vel(i,2)+vel_ref(i,2))
      vel_tmp(3) = 0.5_wip*(vel(i,3)+vel_ref(i,3))
      ekin = ekin + mass(i) &
                  *(vel_tmp(1)*vel_tmp(1)+vel_tmp(2)*vel_tmp(2)  &
                   +           vel_tmp(3)*vel_tmp(3))
    end do
    ekin   = 0.5_dp * ekin

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif
    tempt  = 2.0_wip * ekin / (num_degree * KBOLTZ)

    ! calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
    !   v(t+3/2dt) = v(t+1/2dt)*(2*beta-1) + beta*dt*F(t+dt)/m
    !   r(t+2dt) = r(t+dt) + dt*vel(t+3/2dt)
    !
    beta    = sqrt(temp0/tempt)
    factor  = beta * dt
    factor1 = 2.0_wip*beta - 1.0_wip

    do i = 1, domain%num_atom_domain
      factor2  = factor/mass(i)
      vel(i,1) = factor1*vel_ref(i,1) + factor2*force(i,1)
      vel(i,2) = factor1*vel_ref(i,2) + factor2*force(i,2)
      vel(i,3) = factor1*vel_ref(i,3) + factor2*force(i,3)
      coord(i,1) = coord_ref(i,1) + dt*vel(i,1)
      coord(i,2) = coord_ref(i,2) + dt*vel(i,2)
      coord(i,3) = coord_ref(i,3) + dt*vel(i,3)
    end do

    return

  end subroutine gaussian_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_nvt
  !> @brief        Controls temperature using Langevin dynamics.
  !! @authors      JJ, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_leapfrog_nvt(dynamics, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, temp0, gamma_t, inv_dt
    real(wip)                :: scale_v, scale_f, factor
    real(wip)                :: sigma
    real(wip)                :: v1, v2, rsq, grandom(3)
    integer                  :: j, jx, ncell, start_j

    real(wip),       pointer :: coord(:,:), force(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:)
    real(dp),        pointer :: viri(:,:), viri_const(:,:)
    real(wip),       pointer :: mass(:)
    integer,         pointer :: natom(:)


    natom      => domain%num_atom
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    force      => domain%force
    viri       => dynvars%virial
    viri_const => dynvars%virial_const

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    temp0      =  ensemble%temperature
    ncell      =  domain%num_cell_local
    gamma_t    =  ensemble%gamma_t*AKMA_PS


    ! add random force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
    !
    factor = 2.0_wip*gamma_t*KBOLTZ*temp0/dt

    do j = 1, domain%num_atom_domain

      sigma = sqrt(mass(j) * factor)
      rsq = 2.0_wip

      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(1) = rsq * v1

      rsq = 2.0_wip

      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(2) = rsq * v1

      rsq = 2.0_wip

      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(3) = rsq * v1

      force(j,1) = force(j,1) + sigma * grandom(1)
      force(j,2) = force(j,2) + sigma * grandom(2)
      force(j,3) = force(j,3) + sigma * grandom(3)

    end do

    ! Langevin dynamics (Langevin thermostat)
    !   calculate velocities(t + 3/2dt)
    !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
    !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
    !
    factor  = 1.0_wip + 0.5_wip * dt * gamma_t
    scale_v = (2.0_wip/factor) - 1.0_wip
    scale_f = dt/factor

    !$omp parallel do default(shared)     &
    !$omp private(j, jx, start_j)                                               
    !
    do j = 1, domain%num_atom_domain
      vel(j,1) = scale_v*vel(j,1) + scale_f*force(j,1)/mass(j)
      vel(j,2) = scale_v*vel(j,2) + scale_f*force(j,2)/mass(j)
      vel(j,3) = scale_v*vel(j,3) + scale_f*force(j,3)/mass(j)
      coord(j,1) = coord(j,1) + dt*vel(j,1)
      coord(j,2) = coord(j,2) + dt*vel(j,2)
      coord(j,3) = coord(j,3) + dt*vel(j,3)
    end do
    !$omp end parallel do

    return

  end subroutine langevin_leapfrog_nvt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_npt
  !> @brief        Implements Langevin thermostat and barostat.
  !! @authors      JJ, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_leapfrog_npt(dynamics, ensemble, domain, &
                                   boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, inv_dt, temp0, press0, gamma0, pressxy0
    real(dp)                 :: d_ndegf, kin(1:3), ekin
    real(dp)                 :: virial_sum(1:3)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(wip)                :: coord_tmp(1:3)
    real(wip)                :: scale_f(3), scale_v(3), fact(3), factor
    real(wip)                :: bmoment_ref(3), bmoment2(3), scale_b(3)
    real(wip)                :: gamma_t, gamma_p, pmass, pforce(3)
    real(wip)                :: sigma
    real(wip)                :: v1, v2, rsq, grandom(3)
    integer                  :: i, j, jx, i_ndegf, ncell, nboundary
    integer                  :: start_j

    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: coord(:,:)
    real(wip),       pointer :: coord_ref(:,:), coord_old(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)
    real(wip),       pointer :: bmoment(:)
    integer,         pointer :: natom(:)


    mass       => domain%mass
    natom      => domain%num_atom
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    coord_old  => domain%coord_old
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    bmoment    => dynvars%barostat_momentum

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P
    gamma0     =  ensemble%gamma    * ATMOS_P*100.0_wip/1.01325_wip !TODO
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    gamma_p    =  ensemble%gamma_p * AKMA_PS
    i_ndegf    =  domain%num_deg_freedom
    d_ndegf    =  real(i_ndegf,dp)
    ncell      =  domain%num_cell_local
    nboundary  =  domain%num_cell_boundary


    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    bmoment_ref(1:3) = bmoment(1:3)

    ! Initial guess to calculate eta(t+dt)
    !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
    !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
    !
    !$omp parallel do default(shared) private(j, jx, start_j)                 
    !
    do j = 1, domain%num_atom_domain
      vel(j,1) = vel(j,1) + dt*force(j,1)/mass(j)
      vel(j,2) = vel(j,2) + dt*force(j,2)/mass(j)
      vel(j,3) = vel(j,3) + dt*force(j,3)/mass(j)
      coord(j,1) = coord(j,1) + dt*vel(j,1)
      coord(j,2) = coord(j,2) + dt*vel(j,2)
      coord(j,3) = coord(j,3) + dt*vel(j,3)
    end do
    !$omp end parallel do

    ! add random_force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
    !
    factor   = 2.0_wip*gamma_t*KBOLTZ*temp0/dt

    do j = 1, domain%num_atom_domain
      sigma = sqrt(mass(j) * factor)
      rsq = 2.0_wip

      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(1) = rsq * v1

      rsq = 2.0_wip

      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(2) = rsq * v1

      rsq = 2.0_wip
      do while (rsq >= 1.0_wip)
        v1  = 2.0_wip*random_get() - 1.0_wip
        v2  = 2.0_wip*random_get() - 1.0_wip
        rsq = v1*v1 + v2*v2
      end do

      rsq   = sqrt(-2.0_wip * log(rsq) / rsq)
      grandom(3) = rsq * v1

      force(j,1:3) = force(j,1:3) + sigma * grandom(1:3)

    end do

    ! calculate piston mass (pmass) and stochastic force (pforce)
    ! acting on the barostat
    !
    if (ensemble%isotropy == IsotropyISO) then

      pmass = real(i_ndegf+3,dp)*KBOLTZ*temp0 / (2.0_wip*PI*gamma_p)**2
      sigma = sqrt(2.0_wip*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = pforce(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      pmass = real(i_ndegf+3,dp)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
      sigma = sqrt(2.0_wip*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyANISO) then

      pmass = real(i_ndegf+3,wip)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
      sigma = sqrt(2.0_wip*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = sigma * random_get_gauss()
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyXY_FIXED) then

      pmass = real(i_ndegf+1,dp)*KBOLTZ*temp0/(3.0_wip*(2.0_wip*PI*gamma_p)**2)
      sigma = sqrt(2.0_wip*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = 0.0_wip
      pforce(2) = 0.0_wip
      pforce(3) = sigma * random_get_gauss()

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(pforce, 3, mpi_real8, 0, mpi_comm_country, ierror)
#endif

    ! iteration
    !   Since leapfrog algorithm cannot directly compute vel(t+dt),
    !   iteration scheme is required to obtain T(t+dt) and P(t+dt).
    !
    do i = 1, 8

      ! compute kinetic energy(t+dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      kin(1:3) = 0.0_dp

      do j = 1, domain%num_atom_domain
        kin(1:3) = kin(1:3) + mass(j) * vel(1:3,j)*vel(1:3,j)
        kin(1:3) = kin(1:3) + mass(j) * vel_ref(1:3,j)*vel_ref(1:3,j)
      end do
      kin(1:3) = kin(1:3) / 2.0_dp

!     ! scale factor for harmonic oscillators (R.Pastor et al., Mol.Phys. 1988)
!     ! it is not used in the barostat (Jaewoon Jung)
!     if (i /= 1) then
!       kin(1:3) = kin(1:3) * (1.0_dp + 0.5_dp*gamma_t*dt)
!     end if

      ekin   = 0.5_dp * (kin(1) + kin(2) + kin(3))
      call reduce_pres(kin, ekin, virial(1,1), virial(2,2), virial(3,3),    &
                       virial_sum)

      ! compute pressure(t+dt) in the unit of kcal/mol*A3
      !
      press(1:3) = (kin(1:3) + virial_sum(1:3)) / volume
      pressxyz   = (press(1) + press(2) + press(3))/3.0_dp
      pressxy    = (press(1) + press(2))/2.0_dp

      ! compute thermostat and barostat parameters
      !   for isotropic systems
      !     eta(t+3/2dt) = eta(t+1/2dt) + dt[dim*V(t+dt)*(P(t+dt)-Pext)
      !                                     + dim*2Ekin(t+dt)/Nf +Rp]/pmass
      !     eta(t+3/2dt) = eps(-gamma_p*dt)*eta(t+3/2dt)
      !     eta(t+dt)    = 0.5*(eta(t+3/2dt) + eta(t+1/2dt))
      !     factor       = 1 + [gamma_t+(1+dim/Nf)*eta(t+dt)]dt/2
      !
      if (ensemble%isotropy == IsotropyISO) then

        bmoment(1) = bmoment_ref(1) + dt*(3.0_wip*volume*(pressxyz - press0)   &
                                    + 6.0_wip*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(1)

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wip*(bmoment(1:3) + bmoment_ref(1:3))

        fact(1:3) = 1.0_wip+                                                   &
           (gamma_t+bmoment2(1:3)*(1.0_wip+3.0_wip/d_ndegf))*0.5_wip*dt

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(pressxy - press0)      &
                       + 2.0_wip*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3)   = bmoment_ref(3) + dt*(volume*(press(3)  - press0)      &
                       + 2.0_wip*ekin/d_ndegf + pforce(3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(pressxy - pressxy0)    &
                       + 2.0_wip*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3)  - press0)        &
                     + 2.0_wip*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wip*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wip+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wip*dt

      else if (ensemble%isotropy == IsotropyANISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1:3) = bmoment_ref(1:3) + dt*(volume*(press(1:3) - press0)   &
                       + 2.0_wip*ekin/d_ndegf + pforce(1:3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(press(1:2) - pressxy0) &
                       + 2.0_wip*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3) - press0)         &
                     + 2.0_wip*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wip*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wip+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wip*dt

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        bmoment(1) = 0.0_wip
        bmoment(2) = 0.0_wip
        bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3) - press0)           &
                                    + 2.0_wip*ekin/d_ndegf + pforce(3))/pmass

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wip*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wip+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wip*dt

      end if

      scale_v(1:3) = 2.0_wip/fact(1:3) - 1.0_wip
      scale_f(1:3) = dt/fact(1:3)


      ! Langevin dynamics
      !   calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
      !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
      !   r(t+3/2dt) = 0.5*(r(t+2dt) + r(t+dt))
      !   r(t+2dt)   = r(t+dt) + dt*(v(t+3/2dt) + eta(t+3/2dt)*r(t+3/2dt)))
      !
      do j = 1, domain%num_atom_domain+domain%num_atom_boundary
        vel(1:3,j) = scale_v(1:3)*vel_ref(1:3,j) &
                   + scale_f(1:3)*force(1:3,j)/mass(j)
        coord_tmp(1:3) = 0.5_wip*(coord_ref(1:3,j) + coord(1:3,j))
        coord(1:3,j) = coord_ref(1:3,j) &
                     + dt*(vel(1:3,j) + bmoment(1:3)*coord_tmp(1:3))
      end do

    end do


    ! compute box size(t+2dt)
    !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
    !

    scale_b(1:3) = exp(bmoment(1:3)*dt)
    boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
    boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
    boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

    call bcast_boxsize(boundary%box_size_x, &
                       boundary%box_size_y, &
                       boundary%box_size_z)

    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,wip)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,wip)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,wip)

    ! update boudary conditions
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)
    do j = 1, domain%num_atom_domain
      domain%trans_vec(1:3,j) = domain%trans_vec(1:3,j) * real(scale_b(1:3),wp)
    end do

    domain%system_size(1) = real(boundary%box_size_x,wp)
    domain%system_size(2) = real(boundary%box_size_y,wp)
    domain%system_size(3) = real(boundary%box_size_z,wp)

    return

  end subroutine langevin_leapfrog_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_leapfrog
  !> @brief        Changes target temperature linearly.
  !! @authors      TM
  !! @param[in]    dynamics: dynamics information
  !! @param[out]   ensemble: ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine simulated_annealing_leapfrog(dynamics, ensemble)
  
    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble
    
    ! local variable
    real(dp)                 :: old_temperature
    

    if (.not. dynamics%annealing) return

    old_temperature      = ensemble%temperature
    ensemble%temperature = ensemble%temperature + dynamics%dtemperature
    
    if (main_rank) then
      write(MsgOut,'(A,F10.3,A,F10.3)')                                     &
            'Simulated_Annealing_Leapfrog: Anneals temperature from', &
            old_temperature, ' to ', ensemble%temperature

      if (ensemble%temperature < 0.0_dp) &
        call error_msg( &
        'Simulated_Annealing_Leapfrog: Error: Temperature is less than 0 K.')

    end if

    return

  end subroutine simulated_annealing_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_pres
  !> @brief        Reduces pressure.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_pres(val1, val2, val3, val4, val5, val6)

    ! formal arguments
    real(dp),                intent(inout) :: val1(1:3), val2
    real(dp),                intent(inout) :: val3, val4, val5
    real(dp),                intent(inout) :: val6(1:3)

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: before_reduce(7), after_reduce(7)

    before_reduce(1:3)  = val1(1:3)
    before_reduce(4)    = val2
    before_reduce(5)    = val3
    before_reduce(6)    = val4
    before_reduce(7)    = val5

    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)

    val1(1:3) = after_reduce(1:3)
    val2      = after_reduce(4)
    val6(1:3) = after_reduce(5:7)

#else

    val6(1) = val3
    val6(2) = val4
    val6(3) = val5

#endif

    return

  end subroutine reduce_pres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bcast_boxsize
  !> @brief        Broadcasts box size.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bcast_boxsize(val1, val2, val3)

    ! formal arguments
    real(wip),                intent(inout) :: val1, val2, val3

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: list(3)


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
  !  Subroutine    coord_ref_save
  !> @brief        Saves reference coordinates and velocities.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine coord_ref_save(domain)

    ! formal arguments
    type(s_domain),  target,  intent(inout) :: domain

    real(wip),       pointer  :: coord(:,:), coord_ref(:,:) 
    real(wip),       pointer  :: vel(:,:), vel_ref(:,:) 
    integer                   :: i

    coord     => domain%coord
    coord_ref => domain%coord_ref
    vel       => domain%velocity
    vel_ref   => domain%velocity_ref

    !$omp parallel do default(shared) private(i)
    do i = 1, domain%num_atom_domain
      coord_ref(i,1) = coord(i,1)
      coord_ref(i,2) = coord(i,2)
      coord_ref(i,3) = coord(i,3)
      vel_ref  (i,1) = vel  (i,1)
      vel_ref  (i,2) = vel  (i,2)
      vel_ref  (i,3) = vel  (i,3)
    end do
    !$omp end parallel do

    return

  end subroutine coord_ref_save

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_coord_vel
  !> @brief        Updates coordinates and velocities.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_coord_vel(dt, domain)

    ! formal arguments
    real(wip),                intent(in   ) :: dt
    type(s_domain),  target,  intent(inout) :: domain

    real(wip),       pointer  :: coord(:,:), vel(:,:)
    real(wip),       pointer  :: mass(:), force(:,:)
    integer                   :: i
    real(wip)                 :: factor

    coord     => domain%coord
    vel       => domain%velocity
    mass      => domain%mass
    force     => domain%force

    !$omp parallel do default(shared) private(i, factor)
    do i = 1, domain%num_atom_domain
      factor = dt/mass(i)
      vel(i,1) = vel(i,1)+factor*force(i,1)
      vel(i,2) = vel(i,2)+factor*force(i,2)
      vel(i,3) = vel(i,3)+factor*force(i,3)
      coord(i,1) = coord(i,1) + dt*vel(i,1)
      coord(i,2) = coord(i,2) + dt*vel(i,2)
      coord(i,3) = coord(i,3) + dt*vel(i,3)
    end do
    !$omp end parallel do

    return

  end subroutine update_coord_vel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    coord_mid_save
  !> @brief        Saves midpoint coordinates.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine coord_mid_save(domain)

    ! formal arguments
    type(s_domain),  target,  intent(inout) :: domain

    real(wip),       pointer  :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer  :: coord_old(:,:)
    integer                   :: i

    coord     => domain%coord
    coord_ref => domain%coord_ref
    coord_old => domain%coord_old

    !$omp parallel do default(shared) private(i)
    do i = 1, domain%num_atom_domain
      coord_old(i,1) = 0.5_wip*(coord(i,1)+coord_ref(i,1))
      coord_old(i,2) = 0.5_wip*(coord(i,2)+coord_ref(i,2))
      coord_old(i,3) = 0.5_wip*(coord(i,3)+coord_ref(i,3))
    end do
    !$omp end parallel do

    return

  end subroutine coord_mid_save

end module cg_md_leapfrog_mod
