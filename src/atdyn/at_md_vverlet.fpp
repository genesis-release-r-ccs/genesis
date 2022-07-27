!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_md_vverlet_mod
!> @brief   perform molecular dynamics simulation with velocity verlet algorithm
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!           Tadashi Ando (TA), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_md_vverlet_mod

  use at_output_mod
  use at_dynvars_mod
  use at_ensemble_mod
  use at_constraints_mod
  use at_boundary_mod
  use at_pairlist_mod
  use at_energy_str_mod
  use at_energy_mod
  use at_output_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_remd_str_mod
  use at_rpath_str_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_enefunc_gbsa_mod
  use at_enefunc_gamd_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_pdb_mod
  use timers_mod
  use random_mod
  use math_libs_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: vverlet_dynamics
  private :: initial_vverlet
  private :: control_temp_pres_vver1
  private :: control_temp_pres_vver2
  private :: evans_thermostat_vverlet 
  private :: berendsen_thermostat_vverlet
  private :: andersen_thermostat_vverlet
  private :: nosehoover_thermostat_vverlet
  private :: langevin_thermostat_vv1
  private :: langevin_thermostat_vv2
  private :: langevin_barostat_vv1
  private :: langevin_barostat_vv2
  private :: update_barostat
  private :: bussi_thermostat_vverlet
  private :: bussi_barostat_vv1
  private :: bussi_barostat_vv2
  private :: simulated_annealing_vverlet
  private :: update_qm_charge_vverlet


contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vverlet_dynamics
  !> @brief        velocity verlet integrator
  !! @authors      JJ, TM, CK, TA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] ensemble    : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vverlet_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                              pairlist, boundary, constraints, ensemble)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics), target, intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble

    ! local variables
    real(wp)                  :: simtim, dt
    real(wp)                  :: shrink_ratio_x, shrink_ratio_y, shrink_ratio_z 
    integer                   :: i, j, natom
    integer                   :: nsteps, istart, iend
    integer                   :: count_qm
    integer                   :: icount
    character(256)            :: folder,basename
    character                 :: num*9
    logical                   :: savefile

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: force(:,:), mass(:), inv_mass(:)
    real(wp),         pointer :: force_omp(:,:,:)


    mass      => molecule%mass
    inv_mass  => molecule%inv_mass
    coord     => dynvars%coord
    coord_pbc => dynvars%coord_pbc
    coord_ref => dynvars%coord_ref
    force     => dynvars%force
    force_omp => dynvars%force_omp
    vel       => dynvars%velocity
    vel_ref   => dynvars%velocity_ref

    natom     =  molecule%num_atoms
    istart    =  dynamics%istart_step
    iend      =  dynamics%iend_step
    nsteps    =  dynamics%nsteps
    dt        =  dynamics%timestep/AKMA_PS
    simtim    =  dynamics%initial_time

    if (dynamics%target_md) enefunc%rmsd_force = 1.0_wp / (dt*dt)
    if (abs(dynamics%initial_value) < 0.001_wp)  then
      if (dynamics%target_md) &
        dynamics%initial_value =  &
          dynvars%energy%restraint_cv(enefunc%target_function)
      if (dynamics%steered_md) &
        dynamics%initial_value =  &
          dynvars%energy%restraint_cv(enefunc%steered_function)
    end  if
    enefunc%target_value = dynamics%initial_value 

    ! first-step MD
    ! 
    if (.not. dynamics%restart) then
      call initial_vverlet(output, molecule, enefunc, dynamics, pairlist, &
                           boundary, ensemble, constraints, dynvars)
    else
      ! Remove velocity of fixed atoms
      ! 
      if (constraints%num_fixatm > 0) &
        call clear_fixatm_component(constraints, natom, dynvars%velocity)
    end if

    !
    ! stop NPAT & NPgT
    ! 
    if (ensemble%ensemble == EnsembleNPAT .or.  &
        ensemble%ensemble == EnsembleNPT  .or.  &
        ensemble%ensemble == EnsembleNPgT)      & 
      call error_msg('Vverlet_dynamics> Barostats are not allowed in ATDYN')


    ! Reset qm_count
    ! 
    if(enefunc%qmmm%do_qmmm) then 
      count_qm = 0
      enefunc%qmmm%qm_count = istart
    end if

    ! Main MD loop
    ! at this point
    !   coord, vel, and force are at t = 0   , if restart off
    !   coord, vel, and force are at t = t+dt, if restart on
    !

    do i = istart, iend
 
      call timer(TimerIntegrator, TimerOn)

      simtim = simtim + dt*AKMA_PS
      dynvars%time = simtim
      dynvars%step = i

      if (dynamics%target_md .or. dynamics%steered_md) &
        enefunc%target_value = dynamics%initial_value &
                            + (dynamics%final_value-dynamics%initial_value) &
                             *real(dynvars%step,wp)/real(nsteps,wp)

      enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

      ! Save coordinates(t) and velocities(t)
      !
      do j = 1, natom
        coord_ref(1,j) = coord(1,j)
        coord_ref(2,j) = coord(2,j)
        coord_ref(3,j) = coord(3,j)
        vel_ref(1,j) = vel(1,j)
        vel_ref(2,j) = vel(2,j)
        vel_ref(3,j) = vel(3,j)
      end do


      ! Compute velocities(t + 1/2 dt) and coordinates(t + dt)
      !
      if (ensemble%tpcontrol /= TpcontrolNoseHoover .and. &
          ensemble%tpcontrol /= TpcontrolLangevin   .and. &
          ensemble%tpcontrol /= TpcontrolBussi      ) then

        ! Newtonian dynamics
        !   v(t+1/2dt) = v(t) + 0.5dt*F(t)/m
        !   r(t+dt)    = r(t) + dt*v(t+1/2dt)
        !
        do j = 1, natom
          vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
          vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
          vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)

          coord(1,j) = coord(1,j) + dt*vel(1,j)
          coord(2,j) = coord(2,j) + dt*vel(2,j)
          coord(3,j) = coord(3,j) + dt*vel(3,j)
        end do

        ! Coordinate constraint (RATTLE VV1)
        !
        if (constraints%rigid_bond) then
          call compute_constraints(ConstraintModeVVER1, dt, molecule, &
                                   dynvars, constraints)
        end if

      else

        call control_temp_pres_vver1(molecule, dynamics, ensemble,  &
                                     constraints, boundary, dynvars)

      end if

      call timer(TimerIntegrator, TimerOff)

      ! calculate potential energy(t + dt), force(t + dt), and virial(t + dt)
      !
      call compute_energy(molecule, enefunc, pairlist, boundary, &
                          mod(i,dynamics%eneout_period) == 0,    &
                          enefunc%nonb_limiter,  &
                          dynvars%coord,         &
                          dynvars%trans,         &
                          dynvars%coord_pbc,     &
                          dynvars%energy,        &
                          dynvars%temporary,     &
                          dynvars%force,         &
                          dynvars%force_omp,     &
                          dynvars%virial,        &
                          dynvars%virial_extern, &
                          constraints)


      call timer(TimerIntegrator, TimerOn)

      if (ensemble%tpcontrol == TpcontrolLangevin) then

        if (ensemble%ensemble == EnsembleNVT) then

          call langevin_thermostat_vv2(molecule, dynamics, ensemble, &
                                       boundary, constraints, dynvars)

        else

          call langevin_barostat_vv2  (molecule, dynamics, ensemble, &
                                       constraints, boundary, dynvars)

        end if

      else

        ! Compute velocities(t + dt)
        !   v(t+dt) = v(t+1/2dt) + 0.5dt*F(t+dt)/m
        !
        do j = 1, natom
          vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
          vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
          vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)
        end do

        ! Velocty constraint (RATTLE VV2)
        !   coord     is at t + dt, vel (unconstrained) is at t + dt
        !   coord_ref is at t,      vel_ref             is at t
        !   compute constrained velocities(t+dt)
        !   add constraint virial(t+dt)
        !   Note that constraint virial(t) is added to virial(t+dt) in leapfrog
        !             which is correct?
        !
        if (constraints%rigid_bond) then
  
          call compute_constraints(ConstraintModeVVER2, &
                                   dt, molecule, dynvars, constraints) 
  
        end if

        ! Control temperature VV2
        !
        if (ensemble%ensemble /= EnsembleNVE) then
  
          call control_temp_pres_vver2(molecule, dynamics, ensemble, dynvars)
  
        end if
      end if

      call timer(TimerIntegrator, TimerOff)

      ! Remove translational and rotational motion about COM(t + dt)
      !
      if (dynamics%stoptr_period > 0) then

        if (mod(i,dynamics%stoptr_period) == 0) then
          call stop_trans_rotation(molecule%num_atoms, molecule%mass, &
                                   dynamics%stop_com_translation,     &
                                   dynamics%stop_com_rotation,        &
                                   dynvars%coord,                     &
                                   constraints%fixatm,                &
                                   dynvars%velocity)
        end if

      end if


      ! Update nonbond pairlist for coordinates(t + 2dt)
      !
      if (dynamics%nbupdate_period > 0) then
        if (mod(i,dynamics%nbupdate_period) == 0 .and. real_calc) then

          ! Update boundary if pressure is controlled
          !
          if (ensemble%use_barostat) then
            if (enefunc%forcefield == ForcefieldRESIDCG) then
              call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
                  enefunc%cg_pairlistdist_126,                     &
                  enefunc%cg_pairlistdist_PWMcos,                  &
                  enefunc%cg_pairlistdist_DNAbp,                   &
                  enefunc%cg_pairlistdist_exv,                     &
                  boundary)
            else
              call update_boundary(enefunc%table%table,               &
                  enefunc%pairlistdist,                               &
                  boundary)
            end if
          end if

          ! Shrink box
          !
          if (dynamics%shrink_box) then
            if (.not. ensemble%use_barostat .and. &
                mod(i, dynamics%shrink_period) == 0) then
              boundary%box_size_x     = boundary%box_size_x - dynamics%dbox_x
              boundary%box_size_y     = boundary%box_size_y - dynamics%dbox_y
              boundary%box_size_z     = boundary%box_size_z - dynamics%dbox_z
              boundary%box_size_x_ref = boundary%box_size_x
              boundary%box_size_y_ref = boundary%box_size_y
              boundary%box_size_z_ref = boundary%box_size_z

              shrink_ratio_x = boundary%box_size_x &
                               / (boundary%box_size_x + dynamics%dbox_x)
              shrink_ratio_y = boundary%box_size_y &
                               / (boundary%box_size_y + dynamics%dbox_y)
              shrink_ratio_z = boundary%box_size_z &
                               / (boundary%box_size_z + dynamics%dbox_z)

              do j = 1, natom
                coord(1,j) = coord(1,j) * shrink_ratio_x
                coord(2,j) = coord(2,j) * shrink_ratio_y
                coord(3,j) = coord(3,j) * shrink_ratio_z
                vel(1,j)   = vel(1,j) * shrink_ratio_x
                vel(2,j)   = vel(2,j) * shrink_ratio_y
                vel(3,j)   = vel(3,j) * shrink_ratio_z
              end do

              if (enefunc%forcefield == ForcefieldRESIDCG) then
                call update_boundary_cg(enefunc%cg_pairlistdist_ele, &
                    enefunc%cg_pairlistdist_126,                     &
                    enefunc%cg_pairlistdist_PWMcos,                  &
                    enefunc%cg_pairlistdist_DNAbp,                   &
                    enefunc%cg_pairlistdist_exv,                     &
                    boundary)
              else
                call update_boundary(enefunc%table%table,               &
                    enefunc%pairlistdist,                               &
                    boundary)
              end if
            end if
          end if

          call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                               coord_pbc, pairlist)

          !if(enefunc%spot_use) &
          !  call update_spot_atomlist(molecule%num_atoms, coord, boundary)

        end if
      end if


      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then

        if (mod(i,dynamics%anneal_period) == 0) then
          call simulated_annealing_vverlet(dynamics, enefunc, ensemble)
        end if

      end if

      ! Update GAMD
      !
      if (enefunc%gamd%update_period > 0) then
        if (mod(i,enefunc%gamd%update_period) == 0) then
          call update_gamd_vverlet(output, enefunc, dynvars)
        end if
      end if

      ! Update QM charge for ESP-MM
      !
      if (dynamics%esp_mm .and. dynamics%calc_qm_period > 0) then
        if (mod(i,dynamics%calc_qm_period) == 0) then
          count_qm = count_qm + 1
          call update_qm_charge_vverlet(enefunc, molecule, pairlist, dynamics, i, &
              count_qm, dynvars, enefunc%qmmm)
        end if
      end if

      ! Output energy(t + dt) and dynamical variables(t + dt)
      !
      call output_md(output, molecule, enefunc, dynamics, boundary, &
                     ensemble, dynvars)

    end do
    
    return

  end subroutine vverlet_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_vverlet
  !> @brief        compute the first step (0+dt)
  !! @authors      JJ, CK, TM, KY
  !! @param[in]    output      : output information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    pairlist    : non-bond pair list information
  !! @param[in]    boundary    : boundary conditions information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_vverlet(output, molecule, enefunc, dynamics, pairlist, &
                             boundary, ensemble, constraints, dynvars)

    ! formal arguments
    type(s_output),           intent(in)    :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynamics), target, intent(inout) :: dynamics
    type(s_pairlist),         intent(in)    :: pairlist
    type(s_boundary),         intent(in)    :: boundary
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, simtim
    integer                   :: j, natom
    character(256)            :: folder,basename
    character                 :: num*9
    logical                   :: savefile

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: force(:,:), mass(:), inv_mass(:)
    real(wp),         pointer :: force_omp(:,:,:)
    integer,          pointer :: iseed


    ! use pointers
    natom     =  molecule%num_atoms
    coord     => dynvars%coord
    coord_ref => dynvars%coord_ref
    force     => dynvars%force
    force_omp => dynvars%force_omp
    vel       => dynvars%velocity
    vel_ref   => dynvars%velocity_ref
    mass      => molecule%mass
    inv_mass  => molecule%inv_mass
    iseed     => dynamics%iseed
    dt        =  dynamics%timestep/AKMA_PS
    simtim    =  dynamics%initial_time

    dynvars%time = simtim
    dynvars%step = 0


    ! generate initial velocities
    !
    call generate_velocity(ensemble%temperature, molecule%num_atoms, &
                           molecule%mass, iseed, dynvars%velocity)

    ! Remove velocity of fixed atoms
    ! 
    if (constraints%num_fixatm > 0) &
      call clear_fixatm_component(constraints, natom, dynvars%velocity)

    call stop_trans_rotation(molecule%num_atoms, molecule%mass,      &
                             dynamics%stop_com_translation,          &
                             dynamics%stop_com_rotation,             &
                             dynvars%coord,                          &
                             constraints%fixatm,                     &
                             dynvars%velocity)


    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do j = 1, natom
      coord_ref(1,j) = coord(1,j)
      coord_ref(2,j) = coord(2,j)
      coord_ref(3,j) = coord(3,j)
      vel_ref(1,j)   = vel(1,j)
      vel_ref(2,j)   = vel(2,j)
      vel_ref(3,j)   = vel(3,j)
    end do

    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)
      do j = 1, natom
        coord_ref(1,j) = coord(1,j)
        coord_ref(2,j) = coord(2,j)
        coord_ref(3,j) = coord(3,j)
      end do

    end if


    ! calculate energy(0) and forces(0)
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter,  &
                        dynvars%coord,         &
                        dynvars%trans,         &
                        dynvars%coord_pbc,     &
                        dynvars%energy,        &
                        dynvars%temporary,     &
                        dynvars%force,         &
                        dynvars%force_omp,     &
                        dynvars%virial,        &
                        dynvars%virial_extern, &
                        constraints)


    ! if rigid-body on, update velocity(0 + dt/2 and 0 - dt/2)
    !
    if (constraints%rigid_bond) then

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do j = 1, natom
        vel(1,j) = vel_ref(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
        vel(2,j) = vel_ref(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
        vel(3,j) = vel_ref(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)

        coord(1,j) = coord_ref(1,j) + dt*vel(1,j)
        coord(2,j) = coord_ref(2,j) + dt*vel(2,j)
        coord(3,j) = coord_ref(3,j) + dt*vel(3,j)
      end do

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)
      dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)

      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))/dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))/dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))/dt
      end do


      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do j = 1, natom
        vel_ref(1,j) = vel_ref(1,j) - 0.5_wp*dt*force(1,j)*inv_mass(j)
        vel_ref(2,j) = vel_ref(2,j) - 0.5_wp*dt*force(2,j)*inv_mass(j)
        vel_ref(3,j) = vel_ref(3,j) - 0.5_wp*dt*force(3,j)*inv_mass(j)

        coord(1,j) = coord_ref(1,j) - dt*vel_ref(1,j)
        coord(2,j) = coord_ref(2,j) - dt*vel_ref(2,j)
        coord(3,j) = coord_ref(3,j) - dt*vel_ref(3,j)
      end do

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)
      dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)

      do j = 1, natom
        vel_ref(1,j) = (coord_ref(1,j) - coord(1,j))/dt
        vel_ref(2,j) = (coord_ref(2,j) - coord(2,j))/dt
        vel_ref(3,j) = (coord_ref(3,j) - coord(3,j))/dt
      end do

      ! calculate velocity(0)
      !
      do j = 1, natom
        vel_ref(1,j) = 0.5_wp*(vel(1,j) + vel_ref(1,j))
        vel_ref(2,j) = 0.5_wp*(vel(2,j) + vel_ref(2,j))
        vel_ref(3,j) = 0.5_wp*(vel(3,j) + vel_ref(3,j))
      end do

      ! vel <= updated velocities(0) and coord <= constrained coordinates(0)
      !
      do j = 1, natom
        vel(1,j)   = vel_ref(1,j)
        vel(2,j)   = vel_ref(2,j)
        vel(3,j)   = vel_ref(3,j)
        coord(1,j) = coord_ref(1,j)
        coord(2,j) = coord_ref(2,j)
        coord(3,j) = coord_ref(3,j)
      end do

    end if


    ! output dynvars(0)
    !
    call compute_dynvars(molecule, enefunc, dynamics, boundary, ensemble, &
                         dynvars)
    call output_dynvars (output, enefunc, dynvars, ensemble)

    dynamics%restart = .true.

    ! at this point
    !   coord, velocity, and force are at 0 


    return

  end subroutine initial_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_temp_pres_vver1
  !> @brief        driver to control temperature and pressure
  !! @authors      JJ, TA, TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         Nose-Hoover thermostating may have to be skipped at mdstep=1
  !                and Langevin dynamics may have to be Newtonian dynamics 
  !                at mdstep = 1 to get the identical trajectory with LEAP
  !                please ask TM before you edit.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_temp_pres_vver1(molecule, dynamics, ensemble,   &
                                     constraints, boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variable
    real(wp)                  :: dt
    integer                   :: natom, j, alloc_stat

    real(wp),         pointer :: coord(:,:), vel(:,:), force(:,:), inv_mass(:)


    ! use pointers
    !
    inv_mass => molecule%inv_mass
    coord    => dynvars%coord
    force    => dynvars%force
    vel      => dynvars%velocity
    natom    =  molecule%num_atoms
    dt       =  dynamics%timestep/AKMA_PS


    if (ensemble%tpcontrol == TpcontrolLangevin) then

      if (ensemble%ensemble == EnsembleNVT) then

        call langevin_thermostat_vv1(molecule, dynamics, ensemble, &
                                     boundary, constraints, dynvars)

      else

        call langevin_barostat_vv1  (molecule, dynamics, ensemble, &
                                     constraints, boundary, dynvars)

      end if

    else if (ensemble%tpcontrol == TpcontrolNoseHoover) then

      call nosehoover_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

      do j = 1, natom
        vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
        vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
        vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)
        coord(1,j) = coord(1,j) + dt*vel(1,j)
        coord(2,j) = coord(2,j) + dt*vel(2,j)
        coord(3,j) = coord(3,j) + dt*vel(3,j)
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        call compute_constraints(ConstraintModeVVER1, dt, molecule, &
                                 dynvars, constraints)
        dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)
      end if

    else if (ensemble%tpcontrol == TpcontrolBussi) then

      if (ensemble%ensemble == EnsembleNVT) then

        do j = 1, natom
          vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
          vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
          vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)
          coord(1,j) = coord(1,j) + dt*vel(1,j)
          coord(2,j) = coord(2,j) + dt*vel(2,j)
          coord(3,j) = coord(3,j) + dt*vel(3,j)
        end do

        ! RATTLE VV1
        !
        if (constraints%rigid_bond) then
          call compute_constraints(ConstraintModeVVER1, dt, molecule, &
                                   dynvars, constraints)
          dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)
        end if

      else 

        call bussi_barostat_vv1(molecule, dynamics, ensemble, &
                                constraints, boundary, dynvars)

      end if

    end if

    return

  end subroutine control_temp_pres_vver1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_temp_pres_vver2
  !> @brief        driver to control temperature and pressure
  !! @authors      JJ, TM, TA
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         This subroutine may have to be skipped at mdstep = 1
  !                to get the identical trajectory with LEA
  !                please ask TM before you edit this subroutine.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_temp_pres_vver2(molecule, dynamics, ensemble, dynvars)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_dynvars),         intent(inout) :: dynvars


    ! Berendsen thermostat
    !
    if (ensemble%tpcontrol == TpcontrolBerendsen) then 
      call berendsen_thermostat_vverlet (molecule, dynamics, ensemble, dynvars)
      
    ! Andersen thermostat
    !
    else if (ensemble%tpcontrol == TpcontrolAndersen) then
      call andersen_thermostat_vverlet  (molecule, dynamics, ensemble, dynvars)
      
    ! Evans thermostat
    !
    else if (ensemble%tpcontrol == TpcontrolEvans) then
      call evans_thermostat_vverlet     (molecule, dynamics, ensemble, dynvars)
        
    ! Nose-Hoover thermostat
    ! 
    else if (ensemble%tpcontrol == TpcontrolNoseHoover) then
      call nosehoover_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    ! Bussi thermostat
    !
    else if (ensemble%tpcontrol == TpcontrolBussi ) then
 
      call bussi_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    end if

    return
        
  end subroutine control_temp_pres_vver2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    evans_thermostat_vverlet
  !> @brief        control temperature using Evans thermostat
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine evans_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: ekin, dot_vf, factor, temp0, dt
    integer                   :: j, natom

    real(wp),         pointer :: force(:,:)
    real(wp),         pointer :: vel(:,:)
    real(wp),         pointer :: mass(:)


    ! use pointers
    !
    dt    =  dynamics%timestep/AKMA_PS
    temp0 =  ensemble%temperature
    mass  => molecule%mass
    natom =  molecule%num_atoms
    force => dynvars%force
    vel   => dynvars%velocity


    ! calculate kinetic energy (t+dt)
    !
    ekin = 0.0_wp
    do j = 1, natom
      ekin = ekin + mass(j)*(vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)
    end do

    
    ! dot product between velocities and forces
    !
    dot_vf = 0.0_wp
    do j = 1, natom
      dot_vf = dot_vf + vel(1,j)*force(1,j)  
      dot_vf = dot_vf + vel(2,j)*force(2,j)  
      dot_vf = dot_vf + vel(3,j)*force(3,j)  
    end do


    ! calculate the factor 
    ! it is the ratio of the scalar product to the kinetic energy 
    !
    factor = dot_vf/ekin

    ! get the modified velocities
    !
    do j = 1, natom
      vel(1,j) = vel(1,j)*exp(-factor*dt)
      vel(2,j) = vel(2,j)*exp(-factor*dt)
      vel(3,j) = vel(3,j)*exp(-factor*dt)
    end do

    return

  end subroutine evans_thermostat_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    berendsen_thermostat_vverlet
  !> @brief        control temperature using Berendsen thermostat
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine berendsen_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: temp0, tau_t, dt
    real(wp)                  :: ekf, tempf, lambda
    integer                   :: j
    integer                   :: num_degree, natom

    real(wp),         pointer :: vel(:,:)
    real(wp),         pointer :: mass(:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    temp0      =  ensemble%temperature
    tau_t      =  ensemble%tau_t/AKMA_PS
    natom      =  molecule%num_atoms
    num_degree =  molecule%num_deg_freedom
    vel        => dynvars%velocity
    mass       => molecule%mass


    ! calculate temperature(t + dt) 
    !
    ekf = 0.0_wp
    do j = 1, natom
      ekf = ekf + mass(j)*(vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)
    end do
    tempf = ekf/(num_degree*KBOLTZ)


    ! calculate modified velocities(t + dt)
    !
    lambda = sqrt(1.0_wp + (dt/tau_t)*(temp0/tempf - 1.0_wp))
    do j = 1, natom
      vel(1,j) = lambda*vel(1,j)
      vel(2,j) = lambda*vel(2,j)
      vel(3,j) = lambda*vel(3,j)
    end do

    return

  end subroutine berendsen_thermostat_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    andersen_thermostat_vverlet
  !> @brief        control temperature using Andersen thermostat
  !! @authors      JJ
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine andersen_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics), target, intent(inout) :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: kBT, TWO_PI, dt
    real(wp)                  :: temp0, softness, tau_t, cutoff
    real(wp)                  :: rand, rand1, rand2 
    real(wp)                  :: factor, factor1, factor2, alpha, beta
    real(wp)                  :: vel1, vel2, vel3
    integer                   :: j, switch
    integer                   :: num_degree, natom

    real(wp),         pointer :: vel(:,:)
    real(wp),         pointer :: mass(:)
    integer,          pointer :: iseed


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    temp0      =  ensemble%temperature
    softness   =  ensemble%softness
    tau_t      =  ensemble%tau_t
    natom      =  molecule%num_atoms
    num_degree =  molecule%num_deg_freedom
    vel        => dynvars%velocity
    mass       => molecule%mass
    iseed      => dynamics%iseed

    ! initial setup
    !
    kBT     = KBOLTZ * temp0
    TWO_PI  = 2.0_wp * PI
    cutoff  = 1.0_wp - exp(-dt / tau_t)
    factor1 = softness
    factor2 = sqrt(1.0_wp - softness**2)

    ! update velocities
    !
    switch = 0
    do j = 1, natom

      rand = random_get_legacy(iseed)

      ! criteria for velocity updates
      ! 
      if (rand <= cutoff) then
        
        factor = sqrt(kBT / (2.0_wp*mass(j)))
       
        if (switch == 0) then
          
          ! generate random velocities
          !
          rand1 = random_get_legacy(iseed)
          rand2 = random_get_legacy(iseed)
          alpha = factor * sqrt(-2.0_wp * log(rand1)) 
          beta  = TWO_PI * rand2
          vel1  = alpha * cos(beta) 
          vel2  = alpha * sin(beta) 

          rand1 = random_get_legacy(iseed)
          rand2 = random_get_legacy(iseed)
          alpha = factor * sqrt(-2.0_wp * log(rand1)) 
          beta  = TWO_PI * rand2
          vel3  = alpha * cos(beta)

          ! update velocities
          !
          vel(1,j) = factor1*vel(1,j) + factor2*vel1      
          vel(2,j) = factor1*vel(2,j) + factor2*vel2      
          vel(3,j) = factor1*vel(3,j) + factor2*vel3

          switch = 1
        
        else       
       
          vel1  = alpha * sin(beta)

          rand1 = random_get_legacy(iseed)      
          rand2 = random_get_legacy(iseed)     
          alpha = factor * sqrt(-2.0_wp * log(rand1)) 
          beta  = TWO_PI * rand2
           
          vel2  = alpha * cos(beta)          
          vel3  = alpha * sin(beta)          

          ! update velocities
          !
          vel(1,j) = factor1*vel(1,j) + factor2*vel1      
          vel(2,j) = factor1*vel(2,j) + factor2*vel2      
          vel(3,j) = factor1*vel(3,j) + factor2*vel3
          
          switch = 0
        end if

      end if
    
    end do

    return
  
  end subroutine andersen_thermostat_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nosehoover_thermostat_vverlet
  !> @brief        control temperature using Nose-Hoover thermostat
  !! @authors      JJ, TM
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !! @note         S.Nose, J.Chem.Phys., 81, 511-519 (1984)
  !!               S.Nose, Mol.Phys., 52, 255-268 (1984)
  !!               W.G.Hoover, Phys.Rev.A, 31, 1695-1697 (1985)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nosehoover_thermostat_vverlet(molecule, dynamics, ensemble,       &
                                           dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: temp0, tau_t, dt
    real(wp)                  :: ekf, tempf, factor1, scale_v
    integer                   :: j
    integer                   :: num_degree, natom

    real(wp),         pointer :: vel(:,:)
    real(wp),         pointer :: mass(:)
    real(wp),         pointer :: moment


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    temp0      =  ensemble%temperature
    tau_t      =  ensemble%tau_t/AKMA_PS
    natom      =  molecule%num_atoms
    num_degree =  molecule%num_deg_freedom
    vel        => dynvars%velocity
    mass       => molecule%mass
    moment     => dynvars%thermostat_momentum

    ! calculate temperature(t or t + dt)
    !
    ekf = 0.0_wp
    do j = 1, natom
      ekf = ekf + mass(j)*(vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)
    end do
    tempf = ekf/(real(num_degree,wp)*KBOLTZ)

    ! calculate modified velocities(t or t + dt)
    !
    factor1 = tempf/temp0 - 1
    factor1 = factor1/(tau_t**2)

    ! update the Nose-Hoover factor
    !
    moment = moment + factor1*dt/4.0_wp

    ! scale velocities according to the updated Nose-Hoover factor
    ! 
    scale_v = exp(-moment*dt/2.0_wp)
    do j = 1, natom
      vel(1:3,j) = scale_v*vel(1:3,j)
    end do

    ! calculate temperature with updated velocities
    !
    ekf = 0.0_wp
    do j = 1, natom
      ekf = ekf + mass(j)*(vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)
    end do
    tempf = ekf/(real(num_degree,wp)*KBOLTZ)

    ! calculate modified velocities(t or t + dt)
    !
    factor1 = tempf/temp0  - 1
    factor1 = factor1/(tau_t**2)
    factor1 = factor1

    ! update the Nose-Hoover factor
    !
    moment = moment + factor1*dt/4.0_wp

    return

  end subroutine nosehoover_thermostat_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv1
  !> @brief        control temperature using Langevin thermostat
  !! @authors      JJ, KY
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv1(molecule, dynamics, ensemble, &
                                     boundary, constraints, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_molecule), target, intent(in)    :: molecule
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_boundary),         intent(in)    :: boundary
    type(s_constraints), target, intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, h_dt
    real(wp)                  :: temperature, gamma_t, scale_v
    real(wp)                  :: factor, sigma, temp0
    real(wp)                  :: v1, v2, rsq, grandom(1:3)
    real(wp)                  :: kBT, TWO_PI
    real(wp)                  :: tvel(3)
    integer                   :: j, natom

    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),         pointer :: mass(:), inv_mass(:), temporary(:,:)
    real(wp),         pointer :: viri(:,:), viri_const(:,:)
    real(wp),         pointer :: random_force(:,:)
    logical,          pointer :: fixatm(:)


    ! use pointers
    !
    dt            =  dynamics%timestep/AKMA_PS
    temp0         =  ensemble%temperature
    gamma_t       =  ensemble%gamma_t*AKMA_PS
    natom         =  molecule%num_atoms
    mass          => molecule%mass
    inv_mass      => molecule%inv_mass
    coord         => dynvars%coord
    coord_ref     => dynvars%coord_ref
    vel           => dynvars%velocity
    vel_ref       => dynvars%velocity_ref
    force         => dynvars%force
    random_force  => dynvars%random_force
    temporary     => dynvars%temporary
    viri          => dynvars%virial
    viri_const    => dynvars%virial_const
    fixatm        => constraints%fixatm

    ! setup variables
    !
    kBT      = KBOLTZ * temp0
    TWO_PI   = 2.0_wp * PI
    h_dt     = dt / 2.0_wp

    ! random force
    !
    if (dynvars%step == 1) then

      factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/h_dt
      do j = 1, natom
        if (fixatm(j)) cycle
  
        sigma = sqrt(mass(j) * factor)
  
        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(1) = rsq * v1
  
        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(2) = rsq * v1
  
        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(3) = rsq * v1
  
        random_force(1:3,j) = sigma*grandom(1:3)

      end do
 
    end if

    ! scale factor for velocities
    !
    scale_v = exp(-gamma_t * h_dt)
     
    ! thermostat
    !
    do j = 1, natom
      vel(1:3,j) = vel_ref(1:3,j)*scale_v
    end do

    ! VV1
    !
    do j = 1, natom
      vel(1:3,j) = vel(1:3,j) + h_dt*force(1:3,j)*inv_mass(j)
      vel(1:3,j) = vel(1:3,j) + h_dt*random_force(1:3,j)*inv_mass(j)
      coord(1:3,j) = coord_ref(1:3,j) + vel(1:3,j)*dt
    end do

    ! RATTLE VV1
    !
    if (constraints%rigid_bond) then
      do j = 1, natom
        temporary(1,j) = coord(1,j)
        temporary(2,j) = coord(2,j)
        temporary(3,j) = coord(3,j)
      end do

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)
      dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)

      do j = 1, natom
        tvel(1) = (coord(1,j) - temporary(1,j))/dt
        tvel(2) = (coord(2,j) - temporary(2,j))/dt
        tvel(3) = (coord(3,j) - temporary(3,j))/dt
        vel(1,j) = vel(1,j) + tvel(1)
        vel(2,j) = vel(2,j) + tvel(2)
        vel(3,j) = vel(3,j) + tvel(3)
      end do
    end if

    return

  end subroutine langevin_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv2
  !> @brief        control temperature using Langevin thermostat and barostat
  !! @authors      JJ, KY
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    boundary    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv2(molecule, dynamics, ensemble, &
                                     boundary, constraints, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_boundary),         intent(in)    :: boundary
    type(s_constraints), target, intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0
    real(wp)                  :: factor, sigma
    real(wp)                  :: gamma_t
    real(wp)                  :: v1, v2, rsq, grandom(1:3), fcm(0:4)
    real(wp)                  :: h_dt, vel_scale
    integer                   :: j, natom

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: temporary(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: random_f(:,:)
    logical,          pointer :: fixatm(:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    temp0      =  ensemble%temperature
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    mass       => molecule%mass
    inv_mass   => molecule%inv_mass
    natom      =  molecule%num_atoms
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    force      => dynvars%force
    random_f   => dynvars%random_force
    temporary  => dynvars%temporary
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    fixatm     => constraints%fixatm

    ! time step
    !
    h_dt = dt / 2.0_wp

    ! random force
    !
    fcm(0:3) = 0.0_wp
    factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/dt

    do j = 1, natom
      if (fixatm(j)) cycle

      sigma = sqrt(mass(j) * factor)

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(1) = rsq * v1

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(2) = rsq * v1

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(3) = rsq * v1

      random_f(1:3,j) = sigma*grandom(1:3)

    end do

    ! VV2
    !
    do j = 1, natom
      vel(1:3,j) = vel(1:3,j) + h_dt*(force(1:3,j)+random_f(1:3,j))*inv_mass(j)
    end do

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeVVER2, &
                               dt, molecule, dynvars, constraints)
    end if

    ! termostat
    !
    vel_scale = exp(-gamma_t*h_dt)
    do j = 1, natom
      vel(1:3,j) = vel(1:3,j)*vel_scale
    end do

    return

  end subroutine langevin_thermostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv1
  !> @brief        control temperature using Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv1(molecule, dynamics, ensemble,  &
                                   constraints, boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0, press0, d_ndegf
    real(wp)                  :: ke(3), ekin, tvel(3)
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: factor
    real(wp)                  :: gamma_t, gamma_p
    real(wp)                  :: sigma
    real(wp)                  :: v1, v2, rsq, grandom(1:3), fcm(0:4)
    real(wp)                  :: h_dt, q_dt, size_scale(1:3), vel_scale
    real(wp)                  :: virial_constraint(3,3)
    integer                   :: i, j, i_ndegf, natom, maxiter

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), temporary(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: bmoment(:), bmoment_ref(:)
    real(wp),         pointer :: pmass, pforce(:)
    real(wp),         pointer :: random_f(:,:)


    ! use pointers
    !
    dt          =  dynamics%timestep/AKMA_PS
    inv_dt      =  1.0_wp/dt
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    gamma_t     =  ensemble%gamma_t * AKMA_PS
    gamma_p     =  ensemble%gamma_p * AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  molecule%num_deg_freedom
    d_ndegf     =  real(i_ndegf,wp)
    mass        => molecule%mass
    inv_mass    => molecule%inv_mass
    natom       =  molecule%num_atoms
    coord       => dynvars%coord
    coord_ref   => dynvars%coord_ref
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    force       => dynvars%force
    random_f    => dynvars%random_force
    temporary   => dynvars%temporary
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    bmoment     => dynvars%barostat_momentum
    bmoment_ref => dynvars%barostat_momentum_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wp
 
    ! time step
    !
    h_dt = dt / 2.0_wp
    q_dt = dt / 4.0_wp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    if (dynvars%step == 1) then

      ! pmass and stochastic force (pforce)
      !
      if (ensemble%isotropy == IsotropyISO) then
        pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0/(2.0_wp*PI*gamma_p)**2
        sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = pforce(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then
        pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0/(3.0_wp*(2.0_wp*PI*gamma_p)**2)
        sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = sigma * random_get_gauss()

      else if (ensemble%isotropy == IsotropyANISO) then
        pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0/(3.0_wp*(2.0_wp*PI*gamma_p)**2)
        sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = sigma * random_get_gauss()
        pforce(3) = sigma * random_get_gauss()

      else if (ensemble%isotropy == IsotropyXY_Fixed) then
        pmass = real(i_ndegf+1,wp)*KBOLTZ*temp0/(3.0_wp*(2.0_wp*PI*gamma_p)**2)
        sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
        pforce(1) = 0.0_wp
        pforce(2) = 0.0_wp
        pforce(3) = sigma * random_get_gauss()

      end if

      ! random force
      !
      factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/h_dt

      do j = 1, natom
        sigma = sqrt(mass(j) * factor)
  
        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(1) = rsq * v1

        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(2) = rsq * v1
  
        rsq = 2.0_wp
        do while (rsq >= 1.0_wp)
          v1  = 2.0_wp*random_get() - 1.0_wp
          v2  = 2.0_wp*random_get() - 1.0_wp
          rsq = v1*v1 + v2*v2
        end do
        rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
        grandom(3) = rsq * v1

        random_f(1:3,j) = sigma*grandom(1:3)

      end do
 
    end if
    
    ! iteration of barostats
    !
    do i = 1, maxiter

      ! thermostat1 (scale velocities)
      ! 
      vel_scale = exp(-gamma_t*q_dt)
      do j = 1, natom
        vel(1:3,j) = vel_ref(1:3,j)*vel_scale
      end do

      ! compute kinetic energy(t+dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      ke(1:3) = 0.0_wp
      do j = 1, natom
        ke(1)  = ke(1) + mass(j)*vel(1,j)*vel(1,j)
        ke(2)  = ke(2) + mass(j)*vel(2,j)*vel(2,j)
        ke(3)  = ke(3) + mass(j)*vel(3,j)*vel(3,j)
      end do
      ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))

      ! compute pressure(t+dt) in the unit of kcal/mol*A3
      !
      pressx = (ke(1) + virial(1,1) + virial_constraint(1,1))/volume
      pressy = (ke(2) + virial(2,2) + virial_constraint(2,2))/volume
      pressz = (ke(3) + virial(3,3) + virial_constraint(3,3))/volume
      pressxyz = (pressx + pressy + pressz)/3.0_wp
      pressxy  = (pressx + pressy)/2.0_wp

      ! update barostat
      !
      call update_barostat(molecule, ensemble, dynamics, boundary,    &
                           pressx, pressy, pressz, pressxyz, pressxy, &
                           volume, ekin, dynvars)

      ! thermostat2 (scale velocities)
      !
      vel_scale = exp(-gamma_t*q_dt)
      do j = 1, natom
        vel(1:3,j) = vel(1:3,j)*vel_scale
      end do

      ! VV1
      !
      size_scale(1:3) = exp(bmoment(1:3)*dt)
      do j = 1, natom
        vel(1,j) = vel(1,j) + h_dt*force(1,j)*inv_mass(j)
        vel(2,j) = vel(2,j) + h_dt*force(2,j)*inv_mass(j)
        vel(3,j) = vel(3,j) + h_dt*force(3,j)*inv_mass(j)

        vel(1,j) = vel(1,j) + h_dt*random_f(1,j)*inv_mass(j)
        vel(2,j) = vel(2,j) + h_dt*random_f(2,j)*inv_mass(j)
        vel(3,j) = vel(3,j) + h_dt*random_f(3,j)*inv_mass(j)

        coord(1,j) = size_scale(1)*coord_ref(1,j) + vel(1,j)*dt
        coord(2,j) = size_scale(2)*coord_ref(2,j) + vel(2,j)*dt
        coord(3,j) = size_scale(3)*coord_ref(3,j) + vel(3,j)*dt
      end do 

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, natom
          temporary(1,j) = coord(1,j)
          temporary(2,j) = coord(2,j)
          temporary(3,j) = coord(3,j)
        end do

        call compute_constraints(ConstraintModeLEAP, &
                                 dt, molecule, dynvars, constraints)

        virial_constraint(1:3,1:3) = &
             virial_constraint(1:3,1:3) + 2.0_wp * viri_const(1:3,1:3)

        do j = 1, natom
          tvel(1) = (coord(1,j) - temporary(1,j))/dt
          tvel(2) = (coord(2,j) - temporary(2,j))/dt
          tvel(3) = (coord(3,j) - temporary(3,j))/dt

          vel(1,j) = vel(1,j) + tvel(1)
          vel(2,j) = vel(2,j) + tvel(2)
          vel(3,j) = vel(3,j) + tvel(3)

          force(1,j) = force(1,j) + mass(j)*tvel(1)/h_dt
          force(2,j) = force(2,j) + mass(j)*tvel(2)/h_dt
          force(3,j) = force(3,j) + mass(j)*tvel(3)/h_dt
        end do
      end if

    end do

    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    ! compute box size(t+2dt)
    !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
    !
    boundary%box_size_x = exp(bmoment(1)*dt) * boundary%box_size_x_ref
    boundary%box_size_y = exp(bmoment(2)*dt) * boundary%box_size_y_ref
    boundary%box_size_z = exp(bmoment(3)*dt) * boundary%box_size_z_ref

    ! update barostat momentum
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)

    return

  end subroutine langevin_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv2
  !> @brief        control temperature using Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv2(molecule, dynamics, ensemble,  &
                                   constraints, boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0, press0, gamma0, d_ndegf
    real(wp)                  :: ke(3), ekin
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: factor
    real(wp)                  :: gamma_t, gamma_p
    real(wp)                  :: sigma
    real(wp)                  :: v1, v2, rsq, grandom(1:3)
    real(wp)                  :: virial_constraint(1:3,1:3)
    real(wp)                  :: h_dt, q_dt, vel_scale
    integer                   :: i, j, i_ndegf, natom, maxiter

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: force_add(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), vel_old(:,:)
    real(wp),         pointer :: coord(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: bmoment(:), bmoment_ref(:)
    real(wp),         pointer :: pmass, pforce(:)
    real(wp),         pointer :: random_f(:,:)


    ! use pointers
    !
    dt          =  dynamics%timestep/AKMA_PS
    inv_dt      =  1.0_wp/dt
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    gamma0      =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    gamma_t     =  ensemble%gamma_t * AKMA_PS
    gamma_p     =  ensemble%gamma_p * AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  molecule%num_deg_freedom
    d_ndegf     =  real(i_ndegf,wp)
    mass        => molecule%mass
    inv_mass    => molecule%inv_mass
    natom       =  molecule%num_atoms
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    vel_old     => dynvars%coord_ref
    coord       => dynvars%coord
    force       => dynvars%force
    random_f    => dynvars%random_force
    force_add   => dynvars%temporary
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    bmoment     => dynvars%barostat_momentum
    bmoment_ref => dynvars%barostat_momentum_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    h_dt = dt / 2.0_wp
    q_dt = dt / 4.0_wp

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wp

    ! initial constraint force
    !
    force_add(1:3,1:natom) = 0.0_wp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! pmass and stochastic force (pforce)
    !
    if (ensemble%isotropy == IsotropyISO) then
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = pforce(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyANISO) then
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
      pforce(1) = sigma * random_get_gauss()
      pforce(2) = sigma * random_get_gauss()
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)
      pforce(1) = 0.0_wp
      pforce(2) = 0.0_wp
      pforce(3) = sigma * random_get_gauss()

    end if

    ! random force
    !
    factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/dt

    do j = 1, natom
      sigma = sqrt(mass(j) * factor)

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(1) = rsq * v1

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(2) = rsq * v1

      rsq = 2.0_wp
      do while (rsq >= 1.0_wp)
        v1  = 2.0_wp*random_get() - 1.0_wp
        v2  = 2.0_wp*random_get() - 1.0_wp
        rsq = v1*v1 + v2*v2
      end do
      rsq   = sqrt(-2.0_wp * log(rsq) / rsq)
      grandom(3) = rsq * v1

      random_f(1:3,j) = sigma*grandom(1:3)
    end do
   
    vel_ref(1:3,1:natom) = vel(1:3,1:natom)

    ! iteraction of barostats
    !
    do i = 1, maxiter

      vel(1:3,1:natom) = vel_ref(1:3,1:natom)

      ! VV2
      !
      do j = 1, natom
        vel(1:3,j) = vel(1:3,j) + h_dt*force(1:3,j)*inv_mass(j)
        vel(1:3,j) = vel(1:3,j) + h_dt*random_f(1:3,j)*inv_mass(j)
        vel(1:3,j) = vel(1:3,j) + h_dt*force_add(1:3,j)*inv_mass(j)
      end do

      ! RATTLE VV2
      !
      if (constraints%rigid_bond) then
 
        vel_old(1:3,1:natom) = vel(1:3,1:natom) 
        call compute_constraints(ConstraintModeVVER2, &
                                 dt, molecule, dynvars, constraints)
        do j = 1, natom
          force_add(1:3,j) = force_add(1:3,j)                          &
                           + mass(j)*(vel(1:3,j)-vel_old(1:3,j))/h_dt
        end do

        ! constraint virial
        !
        if (i == maxiter) then
          do j = 1, natom
            virial_constraint(1,1) = virial_constraint(1,1)  &
                                   + force_add(1,j)*coord(1,j)
            virial_constraint(2,2) = virial_constraint(2,2)  &
                                   + force_add(2,j)*coord(2,j)
            virial_constraint(3,3) = virial_constraint(3,3)  &
                                   + force_add(3,j)*coord(3,j)
          end do
        end if

      end if

      ! termostat
      ! 
      vel_scale = exp(-gamma_t*q_dt)
      do j = 1, natom
        vel(1:3,j) = vel(1:3,j)*vel_scale
      end do

      ! compute kinetic energy(t+dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      ke(1:3) = 0.0_wp
      do j = 1, natom
        ke(1) = ke(1) + mass(j)*vel(1,j)*vel(1,j)
        ke(2) = ke(2) + mass(j)*vel(2,j)*vel(2,j)
        ke(3) = ke(3) + mass(j)*vel(3,j)*vel(3,j)
      end do
      ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))

      ! compute pressure(t+dt) in the unit of kcal/mol*A3
      !
      pressx   = (ke(1) + virial(1,1) + virial_constraint(1,1))/volume
      pressy   = (ke(2) + virial(2,2) + virial_constraint(2,2))/volume
      pressz   = (ke(3) + virial(3,3) + virial_constraint(3,3))/volume
      pressxyz = (pressx + pressy + pressz)/3.0_wp
      pressxy  = (pressx + pressy)/2.0_wp

      ! update barostat
      !
      call update_barostat(molecule, ensemble, dynamics, boundary,    &
                           pressx, pressy, pressz, pressxyz, pressxy, &
                           volume, ekin, dynvars)

      ! termostat
      !
      do j = 1, natom
        vel(1:3,j) = vel(1:3,j)*vel_scale
      end do

    end do

    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)
 
    return

  end subroutine langevin_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  update_barostat
  !> @brief      update barostat parameter bmoment (eta)
  !! @authors    JJ, TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat(molecule, ensemble, dynamics, boundary,    &
                             pressx, pressy, pressz, pressxyz, pressxy, &
                             volume, ekin, dynvars)

    ! formal arguments
    type(s_molecule),         intent(in)    :: molecule
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_boundary),         intent(in)    :: boundary
    real(wp),                 intent(in)    :: ekin
    real(wp),                 intent(in)    :: pressx, pressy, pressz
    real(wp),                 intent(in)    :: pressxyz, pressxy
    real(wp),                 intent(in)    :: volume
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variable
    real(wp)                  :: vel_scale(1:3)
    real(wp)                  :: d_ndegf
    real(wp)                  :: press0, gamma0, pressxy0
    real(wp)                  :: gamma_p, pmass
    real(wp)                  :: q_dt, h_dt
    integer                   :: i, natom

    real(wp),        pointer  :: pforce(:), vel(:,:)
    real(wp),        pointer  :: bmoment(:), bmoment_ref(:)


    press0      =  ensemble%pressure * ATMOS_P
    gamma0      =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    gamma_p     =  ensemble%gamma_p * AKMA_PS
    pmass       =  ensemble%pmass
    pforce      => ensemble%pforce
    natom       =  molecule%num_atoms
    d_ndegf     =  real(molecule%num_deg_freedom,wp)
    h_dt        =  (dynamics%timestep/AKMA_PS)/2.0_wp
    q_dt        =  (dynamics%timestep/AKMA_PS)/4.0_wp
    vel         => dynvars%velocity
    bmoment     => dynvars%barostat_momentum
    bmoment_ref => dynvars%barostat_momentum_ref

    ! eta(t+1/4dt)
    !
    bmoment(1:3) = exp(-gamma_p*q_dt/2.0_wp)*bmoment_ref(1:3)

    ! eta(t+1/4dt) is scaled according to pressure
    !
    if (ensemble%isotropy == IsotropyISO) then
      bmoment(1) = bmoment(1) + q_dt*(3.0_wp*volume*(pressxyz - press0)        &
                              + 6.0_wp*ekin/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)
      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then

        bmoment(1) = bmoment(1) + q_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

      else if (ensemble%ensemble == EnsembleNPgT) then

        pressxy0 = press0 - gamma0 / boundary%box_size_z

        bmoment(1) = bmoment(1) + q_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

      end if

      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then

        bmoment(1) = bmoment(1) + q_dt*(volume*(pressx - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressy - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

      else if (ensemble%ensemble == EnsembleNPgT) then

        pressxy0 = press0 - gamma0 / boundary%box_size_z

        bmoment(1) = bmoment(1) + q_dt*(volume*(pressx - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressy - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

      end if

      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      bmoment(1) = 0.0_wp
      bmoment(2) = 0.0_wp
      bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)                 &
                              + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    end if

    ! velocities are scaled according to scaled eta(t+1/4dt)
    !
    vel_scale(1:3) = exp(-h_dt*bmoment(1:3)*(1.0_wp+3.0_wp/d_ndegf))
    do i = 1, natom
      vel(1:3,i) = vel(1:3,i) * vel_scale(1:3)
    end do

    ! eta(t+1/4dt) is scaled 
    !
    bmoment(1:3) = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    ! eta(t+1/2dt) 
    !
    if (ensemble%isotropy == IsotropyISO) then
      bmoment(1) = bmoment(1) + q_dt*(3.0_wp*volume*(pressxyz - press0)        &
                               + 6.0_wp*ekin/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)
      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + q_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0   = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + q_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      end if

      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropyANISO) then
      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + q_dt*(volume*(pressx - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressy - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0   = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + q_dt*(volume*(pressx - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + q_dt*(volume*(pressy - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      end if

      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      bmoment(1) = 0.0_wp
      bmoment(2) = 0.0_wp
      bmoment(3) = bmoment(3) + q_dt*(volume*(pressz - press0)                 &
                              + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass
      bmoment(1:3)  = exp(-gamma_p*q_dt/2.0_wp)*bmoment(1:3)

    end if

    return

  end subroutine update_barostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_thermostat_vverlet
  !> @brief        control temperature using Bussi's stochastic re-scaling
  !! @authors      TA
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_thermostat_vverlet(molecule, dynamics, ensemble, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: temp0, tau_t, dt
    real(wp)                  :: ekf, tempf, lambda
    real(wp)                  :: s2, tempt, dtemp
    real(wp)                  :: gr, factor
    integer                   :: j
    integer                   :: num_degree, natom

    real(wp),         pointer :: vel(:,:)
    real(wp),         pointer :: mass(:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    temp0      =  ensemble%temperature
    tau_t      =  ensemble%tau_t/AKMA_PS
    natom      =  molecule%num_atoms
    num_degree =  molecule%num_deg_freedom
    vel        => dynvars%velocity
    mass       => molecule%mass


    ! calculate temperature(t + dt)
    !
    ekf = 0.0_wp
    do j = 1, natom
      ekf = ekf + mass(j)*(vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)
    end do
    tempf = ekf/(num_degree*KBOLTZ)
    ! calculate scaling factor
    !
    s2 = 2.0_wp*temp0**2/num_degree

    if ( tau_t > 0.0_wp ) then

! Ando
!      dtemp = 1.0_wp/(tau_t*s2)*(temp0-tempf)*dt + sqrt(2.0_wp*dt/tau_t)*random_get_gauss()
!      tempt = tempf + dtemp

! Bussi
      factor = exp(-dt/tau_t)
      gr = random_get_gauss()
      tempt = tempf*factor + temp0/num_degree*(1.0_wp-factor)                &
             *(sum_gauss(num_degree-1)+gr**2)                                &
             + 2.0_wp*sqrt(tempf*temp0/num_degree*(1.0_wp-factor)*factor)*gr
      dtemp = tempt - tempf

    else

      tempt = temp0 + sqrt(s2)*random_get_gauss()
      dtemp = tempt - tempf

    end if

    ! calculate modified velocities(t + dt)
    !
    lambda = sqrt(tempt/tempf)
    do j = 1, natom
      vel(1,j) = lambda*vel(1,j)
      vel(2,j) = lambda*vel(2,j)
      vel(3,j) = lambda*vel(3,j)
    end do

    return

  end subroutine bussi_thermostat_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_barostat_vv1
  !> @brief        control temperature using Bussi thermostat and barostat
  !! @authors      TA, TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_barostat_vv1(molecule, dynamics, ensemble,  &
                                constraints, boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, temp0, press0, pressxy0, gamma0, d_ndegf
    real(wp)                  :: ke(3), ekin, tvel(3)
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: factor
    real(wp)                  :: tau_t, tau_p
    real(wp)                  :: gr, ekin0
    real(wp)                  :: h_dt, size_scale(1:3), vel_scale
    real(wp)                  :: vel_scale_2(1:3)
    real(wp)                  :: virial_constraint(3,3)
    real(wp)                  :: ekin_ref, ke_ref(3)
    integer                   :: i, j, i_ndegf, natom, maxiter, ndof0

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), temporary(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: bmoment(:), bmoment_ref(:)
    real(wp),         pointer :: pmass, pforce(:)


    ! use pointers
    !
    dt          =  dynamics%timestep/AKMA_PS
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    gamma0      =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    tau_t       =  ensemble%tau_t / AKMA_PS
    tau_p       =  ensemble%tau_p / AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  molecule%num_deg_freedom
    d_ndegf     =  real(i_ndegf,wp)
    mass        => molecule%mass
    inv_mass    => molecule%inv_mass
    natom       =  molecule%num_atoms
    coord       => dynvars%coord
    coord_ref   => dynvars%coord_ref
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    force       => dynvars%force
    temporary   => dynvars%temporary
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    bmoment     => dynvars%barostat_momentum
    bmoment_ref => dynvars%barostat_momentum_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wp

    ! time step
    !
    h_dt = dt / 2.0_wp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! pmass and ekin0
    !
    if (ensemble%isotropy == IsotropyISO) then
      ndof0 = i_ndegf+1
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      ndof0 = i_ndegf+3
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2 / 3.0_wp
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    else if (ensemble%isotropy == IsotropyANISO) then
      ndof0 = i_ndegf+3
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2 / 3.0_wp
      ekin0 = 0.5_wp*KBOLTZ*temp0 * real(i_ndegf+3,wp)
    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      ndof0 = i_ndegf+1
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    end if

    ! compute kinetic energy
    !
    ke(1:3) = 0.0_wp
    do j = 1, natom
      ke(1)  = ke(1) + mass(j)*vel_ref(1,j)*vel_ref(1,j)
      ke(2)  = ke(2) + mass(j)*vel_ref(2,j)*vel_ref(2,j)
      ke(3)  = ke(3) + mass(j)*vel_ref(3,j)*vel_ref(3,j)
    end do
    ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))
    ke_ref(1:3) = ke(1:3)
    ekin_ref = ekin

    ! compute stochastic scaling factor
    !
    factor = exp(-h_dt/tau_t)
    ekin  = ekin + 0.5_wp*pmass*dot_product(bmoment_ref(1:3),bmoment_ref(1:3))
    gr = random_get_gauss()
    vel_scale = factor &
              + (1.0_wp-factor)*ekin0/(ndof0*ekin)*(sum_gauss(ndof0-1)+gr**2) &
              + 2.0_wp*gr*sqrt(ekin0/(ndof0*ekin)*factor*(1.0_wp-factor))
    vel_scale = sqrt(vel_scale)
    !NOTE: Most of the cases, sign(vel_scale) is +.
    !      For a harmonic osillator case, the value may have -.
    vel_scale = sign(vel_scale,                                               &
                     gr+sqrt(factor*ndof0*ekin/((1.0_wp-factor)*ekin0)))

    ! iteration of barostat
    !
    do i = 1, maxiter

      ! scale velocity
      !
      do j = 1, natom
        vel(1:3,j) = vel_ref(1:3,j)*vel_scale
      end do

      ! scale bmoment (eta) for barostat 
      !
      bmoment(1:3) = bmoment_ref(1:3)*vel_scale

      ! compute kinetic energy
      !
      ke(1:3) = ke_ref(1:3)*vel_scale*vel_scale
      ekin = ekin_ref*vel_scale*vel_scale

      ! compute pressure in the unit of kcal/mol*A3
      !
      pressx = (ke(1) + virial(1,1) + virial_constraint(1,1))/volume
      pressy = (ke(2) + virial(2,2) + virial_constraint(2,2))/volume
      pressz = (ke(3) + virial(3,3) + virial_constraint(3,3))/volume
      pressxyz = (pressx + pressy + pressz)/3.0_wp
      pressxy  = (pressx + pressy)/2.0_wp

      ! update barostat
      !
      if (ensemble%isotropy == IsotropyISO) then
        bmoment(1) = bmoment(1) + h_dt*(3.0_wp*volume*(pressxyz - press0)     &
                                + 6.0_wp*ekin/d_ndegf)/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1) = bmoment(1) + h_dt*(volume*(pressxy - press0)           &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(2) = bmoment(2) + h_dt*(volume*(pressxy - press0)           &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(3) = bmoment(3) + h_dt*(volume*(pressz  - press0)           &
                                  + 2.0_wp*ekin/d_ndegf)/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1) = bmoment(1) + h_dt*(volume*(pressxy - pressxy0)         &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(2) = bmoment(2) + h_dt*(volume*(pressxy - pressxy0)         &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(3) = bmoment(3) + h_dt*(volume*(pressz  - press0)           &
                                  + 2.0_wp*ekin/d_ndegf)/pmass

        end if

      else if (ensemble%isotropy == IsotropyANISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1) = bmoment(1) + h_dt*(volume*(pressx - press0)            &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(2) = bmoment(2) + h_dt*(volume*(pressy - press0)            &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)            &
                                  + 2.0_wp*ekin/d_ndegf)/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1) = bmoment(1) + h_dt*(volume*(pressx - pressxy0)          &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(2) = bmoment(2) + h_dt*(volume*(pressy - pressxy0)          &
                                  + 2.0_wp*ekin/d_ndegf)/pmass
          bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)            &
                                  + 2.0_wp*ekin/d_ndegf)/pmass

        end if

      else if (ensemble%isotropy == IsotropyXY_Fixed) then
        bmoment(1) = 0.0_wp
        bmoment(2) = 0.0_wp
        bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)              &
                                + 2.0_wp*ekin/d_ndegf)/pmass

      end if

      ! VV1
      !
      do j = 1, natom
        vel(1,j) = vel(1,j) + h_dt*force(1,j)*inv_mass(j)
        vel(2,j) = vel(2,j) + h_dt*force(2,j)*inv_mass(j)
        vel(3,j) = vel(3,j) + h_dt*force(3,j)*inv_mass(j)
      end do

      size_scale(1:3)  = exp(bmoment(1:3)*dt)
!      vel_scale_2(1:3) = sinh(bmoment(1:3)*dt)/bmoment(1:3)
      vel_scale_2(1:3) = powersinh(bmoment(1:3)*dt)*dt
      do j = 1, natom

        coord(1,j) = size_scale(1)*coord_ref(1,j) + vel_scale_2(1)*vel(1,j)
        coord(2,j) = size_scale(2)*coord_ref(2,j) + vel_scale_2(2)*vel(2,j)
        coord(3,j) = size_scale(3)*coord_ref(3,j) + vel_scale_2(3)*vel(3,j)

        vel(1,j) = vel(1,j)/size_scale(1)
        vel(2,j) = vel(2,j)/size_scale(2)
        vel(3,j) = vel(3,j)/size_scale(3)

      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then

        do j = 1, natom
          temporary(1,j) = coord(1,j)
          temporary(2,j) = coord(2,j)
          temporary(3,j) = coord(3,j)
        end do

        call compute_constraints(ConstraintModeLEAP, &
                                 dt, molecule, dynvars, constraints)

        virial_constraint(1:3,1:3) = &
             virial_constraint(1:3,1:3) + 2.0_wp * viri_const(1:3,1:3)

        do j = 1, natom
          tvel(1) = (coord(1,j) - temporary(1,j))/dt
          tvel(2) = (coord(2,j) - temporary(2,j))/dt
          tvel(3) = (coord(3,j) - temporary(3,j))/dt

          vel(1,j) = vel(1,j) + tvel(1)
          vel(2,j) = vel(2,j) + tvel(2)
          vel(3,j) = vel(3,j) + tvel(3)

          force(1,j) = force(1,j) + mass(j)*tvel(1)/h_dt
          force(2,j) = force(2,j) + mass(j)*tvel(2)/h_dt
          force(3,j) = force(3,j) + mass(j)*tvel(3)/h_dt
        end do

      end if

    end do

    ! update virial constraint
    !
    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    ! compute box size(t+dt)
    !
    boundary%box_size_x = exp(bmoment(1)*dt) * boundary%box_size_x_ref
    boundary%box_size_y = exp(bmoment(2)*dt) * boundary%box_size_y_ref
    boundary%box_size_z = exp(bmoment(3)*dt) * boundary%box_size_z_ref

    ! update barostat momentum
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)

    return

  end subroutine bussi_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_barostat_vv2
  !> @brief        control temperature using Bussi thermostat and barostat
  !! @authors      TA, TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_barostat_vv2(molecule, dynamics, ensemble,  &
                                constraints, boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, temp0, press0, pressxy0, gamma0, d_ndegf
    real(wp)                  :: ke(3), ekin
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: factor
    real(wp)                  :: tau_t, tau_p
    real(wp)                  :: gr, ekin0
    real(wp)                  :: virial_constraint(1:3,1:3)
    real(wp)                  :: h_dt, vel_scale
    integer                   :: i, j, i_ndegf, natom, maxiter, ndof0

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: force_add(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), vel_old(:,:)
    real(wp),         pointer :: coord(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: bmoment(:), bmoment_ref(:)
    real(wp),         pointer :: pmass, pforce(:)


    ! use pointers
    !
    dt          =  dynamics%timestep/AKMA_PS
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    gamma0      =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    tau_t       =  ensemble%tau_t   / AKMA_PS
    tau_p       =  ensemble%tau_p   / AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  molecule%num_deg_freedom
    d_ndegf     =  real(i_ndegf,wp)
    mass        => molecule%mass
    inv_mass    => molecule%inv_mass
    natom       =  molecule%num_atoms
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    vel_old     => dynvars%coord_ref
    coord       => dynvars%coord
    force       => dynvars%force
    force_add   => dynvars%temporary
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    bmoment     => dynvars%barostat_momentum
    bmoment_ref => dynvars%barostat_momentum_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    h_dt = dt / 2.0_wp

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_wp

    ! initial constraint force
    !
    force_add(1:3,1:natom) = 0.0_wp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! velocity reference
    !
    vel_ref(1:3,1:natom) = vel(1:3,1:natom)

    ! pmass and ekin0
    !
    if (ensemble%isotropy == IsotropyISO) then
      ndof0 = i_ndegf+1
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      ndof0 = i_ndegf+3
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2 / 3.0_wp
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    else if (ensemble%isotropy == IsotropyANISO) then
      ndof0 = i_ndegf+3
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2 / 3.0_wp
      ekin0 = 0.5_wp*KBOLTZ*temp0 * real(i_ndegf+3,wp)
    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      ndof0 = i_ndegf+1
      pmass = ndof0*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_wp*KBOLTZ*temp0 * ndof0
    end if

    ! iteration of barostat
    !
    do i = 1, maxiter

      vel(1:3,1:natom) = vel_ref(1:3,1:natom)

      ! VV2
      !
      do j = 1, natom
        vel(1:3,j) = vel(1:3,j) + h_dt*force(1:3,j)*inv_mass(j)
        vel(1:3,j) = vel(1:3,j) + h_dt*force_add(1:3,j)*inv_mass(j)
      end do

      ! RATTLE VV2
      !
      if (constraints%rigid_bond) then

        vel_old(1:3,1:natom) = vel(1:3,1:natom)
        call compute_constraints(ConstraintModeVVER2, &
                                 dt, molecule, dynvars, constraints)
        do j = 1, natom
          force_add(1:3,j) = force_add(1:3,j)                          &
                           + mass(j)*(vel(1:3,j)-vel_old(1:3,j))/h_dt
        end do

        ! constraint virial
        !
        if (i == maxiter) then
          do j = 1, natom
            virial_constraint(1,1) = virial_constraint(1,1)  &
                                   + force_add(1,j)*coord(1,j)
            virial_constraint(2,2) = virial_constraint(2,2)  &
                                   + force_add(2,j)*coord(2,j)
            virial_constraint(3,3) = virial_constraint(3,3)  &
                                   + force_add(3,j)*coord(3,j)
          end do
        end if

      end if

    end do

    ! compute kinetic energy
    !
    ke(1:3) = 0.0_wp
    do j = 1, natom
      ke(1) = ke(1) + mass(j)*vel(1,j)*vel(1,j)
      ke(2) = ke(2) + mass(j)*vel(2,j)*vel(2,j)
      ke(3) = ke(3) + mass(j)*vel(3,j)*vel(3,j)
    end do
    ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))

    ! compute pressure in the unit of kcal/mol*A3
    !
    pressx   = (ke(1) + virial(1,1) + virial_constraint(1,1))/volume
    pressy   = (ke(2) + virial(2,2) + virial_constraint(2,2))/volume
    pressz   = (ke(3) + virial(3,3) + virial_constraint(3,3))/volume
    pressxyz = (pressx + pressy + pressz)/3.0_wp
    pressxy  = (pressx + pressy)/2.0_wp

    ! update barostat
    !
    bmoment(1:3) = bmoment_ref(1:3)
    if (ensemble%isotropy == IsotropyISO) then
      bmoment(1) = bmoment(1) + h_dt*(3.0_wp*volume*(pressxyz - press0)        &
                              + 6.0_wp*ekin/d_ndegf)/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then

        bmoment(1) = bmoment(1) + h_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(2) = bmoment(2) + h_dt*(volume*(pressxy - press0)              &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(3) = bmoment(3) + h_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf)/pmass

      else if (ensemble%ensemble == EnsembleNPgT) then

        pressxy0 = press0 - gamma0 / boundary%box_size_z

        bmoment(1) = bmoment(1) + h_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(2) = bmoment(2) + h_dt*(volume*(pressxy - pressxy0)            &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(3) = bmoment(3) + h_dt*(volume*(pressz  - press0)              &
                                + 2.0_wp*ekin/d_ndegf)/pmass

      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then

        bmoment(1) = bmoment(1) + h_dt*(volume*(pressx - press0)               &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(2) = bmoment(2) + h_dt*(volume*(pressy - press0)               &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf)/pmass

      else if (ensemble%ensemble == EnsembleNPgT) then

        pressxy0 = press0 - gamma0 / boundary%box_size_z

        bmoment(1) = bmoment(1) + h_dt*(volume*(pressx - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(2) = bmoment(2) + h_dt*(volume*(pressy - pressxy0)             &
                                + 2.0_wp*ekin/d_ndegf)/pmass
        bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)               &
                                + 2.0_wp*ekin/d_ndegf)/pmass

      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      bmoment(1) = 0.0_wp
      bmoment(2) = 0.0_wp
      bmoment(3) = bmoment(3) + h_dt*(volume*(pressz - press0)                 &
                              + 2.0_wp*ekin/d_ndegf)/pmass

    end if

    ! scale velocity
    !
    factor = exp(-h_dt/tau_t)
    ekin  = ekin + 0.5_wp*pmass*dot_product(bmoment(1:3),bmoment(1:3))
    gr = random_get_gauss()
    vel_scale = factor &
              + (1.0_wp-factor)*ekin0/(ndof0*ekin)*(sum_gauss(ndof0-1)+gr**2) &
              + 2.0_wp*gr*sqrt(ekin0/(ndof0*ekin)*factor*(1.0_wp-factor))
    vel_scale = sqrt(vel_scale)
!NOTE: Most of the cases, sign(vel_scale) is +.
!      For a harmonic osillator case, the value may have -.
    vel_scale = sign(vel_scale,                                               &
                gr+sqrt(factor*ndof0*ekin/((1.0_wp-factor)*ekin0)))
    do j = 1, natom
      vel(1:3,j) = vel(1:3,j)*vel_scale
    end do

    ! scale bmoment (eta) for barostat 
    !
    bmoment(1:3) = bmoment(1:3)*vel_scale

    ! update virial constraint
    !
    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    ! update barostat momentum
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)

    return

  end subroutine bussi_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_vverlet
  !> @brief        change target temperature linearly
  !! @authors      TM
  !! @param[in]    dynamics : dynamics information
  !! @param[inout] enefunc  : potential energy function information
  !! @param[inout] ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine simulated_annealing_vverlet(dynamics, enefunc, ensemble)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_ensemble),        intent(inout) :: ensemble

    ! local variable
    real(wp)                 :: old_temperature


    if (.not. dynamics%annealing) &
      return

    old_temperature      = ensemble%temperature
    ensemble%temperature = ensemble%temperature + dynamics%dtemperature

    if (enefunc%eef1_use) then
      call setup_eef1_temperature(ensemble%temperature, enefunc)
    else if (enefunc%gbsa_use) then
      enefunc%gbsa%temperature = ensemble%temperature
    end if
    if (enefunc%cg_ele_calc) then
      enefunc%cg_ele_sol_T = ensemble%temperature
    end if

    if (ensemble%temperature < 0.0_wp) &
      call error_msg( &
        'Simulated_Annealing_Vverlet> Temperature is less than 0 K')

    if (main_rank) &
      write(MsgOut,'(A,F10.3,A,F10.3)')                             &
            'Simulated_Annealing_Vverlet> Anneal temperature from', &
            old_temperature, ' to ', ensemble%temperature

    return

  end subroutine simulated_annealing_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_gamd_vverlet
  !> @brief        update GaMD parameters
  !! @authors      HO
  !! @param[inout] output  : output information
  !! @param[inout] enefunc : potential energy functions information
  !! @param[inout] dynvars : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine update_gamd_vverlet(output, enefunc, dynvars)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    type(s_enefunc_gamd), pointer    :: gamd


    gamd => enefunc%gamd

    ! Compute and update gamd parameters
    !
    call setup_enefunc_gamd(gamd)

    ! output gamd parameters
    !
    call output_gamd(output, dynvars, enefunc)

    return

  end subroutine update_gamd_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_qm_charge_vverlet
  !> @brief        Perform QM/MM to obtain QM charges
  !! @authors      KY
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[in]    molecule : molecular information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    step     : step number
  !! @param[in]    ncount   : number of QM calc
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[inout] qmmm     : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine update_qm_charge_vverlet(enefunc, molecule, pairlist, dynamics, &
                                      step, ncount, dynvars, qmmm)

    ! formal arguments
    type(s_enefunc),    intent(inout) :: enefunc
    type(s_molecule),   intent(inout) :: molecule
    type(s_pairlist),   intent(in)    :: pairlist
    type(s_dynamics),   intent(inout) :: dynamics
    integer,            intent(in)    :: step
    integer,            intent(in)    :: ncount
    type(s_dynvars),    intent(inout) :: dynvars
    type(s_qmmm),       intent(inout) :: qmmm

    ! local variables
    type(s_energy) :: energy
    real(wp)       :: dummy(3,1)


    qmmm%qm_charge_save = qmmm%qm_charge

    qmmm%qm_count     = step
    qmmm%qm_classical = .false.
    call compute_energy_qmmm(enefunc, molecule, pairlist, dynvars%coord, qmmm, energy, dummy)
    qmmm%qm_classical = .true.

    if (dynamics%avg_qm_charge) then
      qmmm%qm_charge = (qmmm%qm_charge_save*real(ncount) + qmmm%qm_charge) &
                         /real(ncount+1)
    end if

    return

  end subroutine update_qm_charge_vverlet

end module at_md_vverlet_mod
