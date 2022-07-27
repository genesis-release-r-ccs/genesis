!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_md_leapfrog_mod
!> @brief   perform molecular dynamics simulation with leapfrog method algorithm
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!!          Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_md_leapfrog_mod
  
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
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

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

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    leapfrog_dynamics
  !> @brief        leapfrog integrator
  !! @authors      TM, CK, JJ, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraint information
  !! @param[inout] ensemble    : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine leapfrog_dynamics(output, molecule, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble

    ! local variables
    real(wp)                  :: simtim, dt, inv_dt, factor
    integer                   :: i, j, natom, nsteps
    integer                   :: istart, iend
    integer                   :: count_qm

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: force(:,:), mass(:), inv_mass(:)
    real(wp),         pointer :: force_omp(:,:,:)
    real(wp),         pointer :: viri(:,:), viri_const(:,:)
    real(wp),         pointer :: viri_ext(:,:)
    character(256)            :: folder,basename
    character                 :: num*9
    logical                   :: savefile

    integer :: icount


    mass       => molecule%mass
    inv_mass   => molecule%inv_mass
    coord      => dynvars%coord
    coord_pbc  => dynvars%coord_pbc
    coord_ref  => dynvars%coord_ref
    force      => dynvars%force
    force_omp  => dynvars%force_omp
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    viri       => dynvars%virial
    viri_const => dynvars%virial_const
    viri_ext   => dynvars%virial_extern

    natom      =  molecule%num_atoms
    nsteps     =  dynamics%nsteps
    istart     =  dynamics%istart_step
    iend       =  dynamics%iend_step
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    simtim     =  dynamics%initial_time

    if (abs(dynamics%initial_value) < 0.001_wp)  then
      if (dynamics%target_md) &
        dynamics%initial_value =  &
          dynvars%energy%restraint_cv(enefunc%target_function)
      if (dynamics%steered_md) &
        dynamics%initial_value =  &
          dynvars%energy%restraint_cv(enefunc%steered_function)
    end if
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt*dt)
    enefunc%target_value = dynamics%initial_value 

    ! first-step MD
    !
    if (.not. dynamics%restart) then
      call initial_leapfrog(output, molecule, enefunc, dynamics, pairlist, &
                            boundary, ensemble, constraints, dynvars)
    else
      ! Remove velocity of fixed atoms
      ! 
      if (constraints%num_fixatm > 0) &
        call clear_fixatm_component(constraints, natom, dynvars%velocity)

    end if

    ! Reset qm_count
    ! 
    if(enefunc%qmmm%do_qmmm) then
      count_qm = 0
      enefunc%qmmm%qm_count = istart
    end if

    ! Main MD loop 
    !   coord is at 0 +  dt and vel is at 0 + 1/2dt, if restart off
    !   coord is at t + 2dt and vel is at t + 3/2dt, if restart on
    !

    do i = istart, iend

      simtim = simtim + dynamics%timestep
      dynvars%time = simtim
      dynvars%step = i

      if (dynamics%target_md .or. dynamics%steered_md) &
        enefunc%target_value = dynamics%initial_value &
                            + (dynamics%final_value-dynamics%initial_value) &
                             *real(dynvars%step,wp)/real(nsteps,wp)

      enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

      ! Save coordinates(t + dt) and velocities(t + 1/2dt)
      !
      do j = 1, natom
        coord_ref(1,j) = coord(1,j)
        coord_ref(2,j) = coord(2,j)
        coord_ref(3,j) = coord(3,j)
        vel_ref(1,j)   = vel(1,j)
        vel_ref(2,j)   = vel(2,j)
        vel_ref(3,j)   = vel(3,j)
      end do


      ! Compute energy(t + dt), force(t + dt), and internal virial(t + dt)
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

      ! Compute velocities(t + 3/2dt) and coordinates(t + 2dt)
      !
      if (ensemble%tpcontrol /= TpcontrolLangevin) then

        ! Newtonian dynamics
        !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
        !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
        !
        do j = 1, natom
          vel(1,j) = vel(1,j) + dt*force(1,j)*inv_mass(j)
          vel(2,j) = vel(2,j) + dt*force(2,j)*inv_mass(j)
          vel(3,j) = vel(3,j) + dt*force(3,j)*inv_mass(j)

          coord(1,j) = coord(1,j) + dt*vel(1,j)
          coord(2,j) = coord(2,j) + dt*vel(2,j)
          coord(3,j) = coord(3,j) + dt*vel(3,j)
        end do


        ! Bond constraint
        !   coord (unconstrained) is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref             is at t +  dt, vel_ref is at t + 1/2dt
        !   compute constrained coordinates(t+2dt) and constraint virial(t+dt)
        !   update velocities(t + 3/2dt):
        !     v(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
        !   add constraint virial(t + dt)
        !
        if (constraints%rigid_bond) then

          call compute_constraints(ConstraintModeLEAP, &
                                   dt, molecule, dynvars, constraints)

          do j = 1, natom
            vel(1,j) = (coord(1,j) - coord_ref(1,j))*inv_dt
            vel(2,j) = (coord(2,j) - coord_ref(2,j))*inv_dt
            vel(3,j) = (coord(3,j) - coord_ref(3,j))*inv_dt
          end do
          viri(1:3,1:3) = viri(1:3,1:3) + viri_const(1:3,1:3)

        end if


        ! Control temperature and pressure
        !   coord     is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref is at t +  dt, vel_ref is at t + 1/2dt
        !   scale velocities(t + 3/2dt) and coordinates(t + 2dt)
        !

        if (ensemble%ensemble /= EnsembleNVE) then
          call control_temp_pres_leap(molecule, dynamics, ensemble, &
                                      constraints, boundary, dynvars)
        end if


      else if (ensemble%tpcontrol == TpcontrolLangevin) then

        ! Langevin dynamics
        !
        if (ensemble%ensemble == EnsembleNVT) then

          call langevin_leapfrog_nvt(molecule, dynamics, ensemble, &
                                     boundary, constraints, dynvars)

        else

          call langevin_leapfrog_npt(molecule, dynamics, ensemble, &
                                     constraints, boundary, dynvars)

        end if

      end if


      ! Remove translational and rotational motion about COM(t + 3/2dt)
      !
      if (dynamics%stoptr_period > 0) then
        if (mod(i,dynamics%stoptr_period) == 0) then

          do j = 1, natom
            dynvars%temporary(1,j) = 0.5_wp*(coord(1,j) + coord_ref(1,j))
            dynvars%temporary(2,j) = 0.5_wp*(coord(2,j) + coord_ref(2,j))
            dynvars%temporary(3,j) = 0.5_wp*(coord(3,j) + coord_ref(3,j))
          end do

          call stop_trans_rotation(molecule%num_atoms, molecule%mass,   &
                                   dynamics%stop_com_translation,       &
                                   dynamics%stop_com_rotation,          &
                                   dynvars%temporary,                   &
                                   constraints%fixatm,                  &
                                   dynvars%velocity)
        end if
      end if

      call timer(TimerIntegrator, TimerOff)


      ! Update boundary and nonbonded pairlist for coordinates(t + 2dt)
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

          if (dynamics%shrink_box) then
            if (.not. ensemble%use_barostat .and. &
                mod(i, dynamics%shrink_period) == 0) then
              boundary%box_size_x     = boundary%box_size_x - dynamics%dbox_x
              boundary%box_size_y     = boundary%box_size_y - dynamics%dbox_y
              boundary%box_size_z     = boundary%box_size_z - dynamics%dbox_z
              boundary%box_size_x_ref = boundary%box_size_x
              boundary%box_size_y_ref = boundary%box_size_y
              boundary%box_size_z_ref = boundary%box_size_z
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

        end if
      end if


      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then
        if (mod(i,dynamics%anneal_period) == 0) then

          call simulated_annealing_leapfrog(dynamics, enefunc, ensemble)

        end if
      end if

      ! Update GAMD
      !
      if (enefunc%gamd%update_period > 0) then
        if (mod(i,enefunc%gamd%update_period) == 0) then
          call update_gamd_leapfrog(output, enefunc, dynvars)
        end if
      end if

      ! Update QM charge for ESP-MM
      !
      if (dynamics%esp_mm .and. dynamics%calc_qm_period > 0) then
        if (mod(i,dynamics%calc_qm_period) == 0) then
          count_qm = count_qm + 1
          call update_qm_charge_leapfrog(enefunc, molecule, pairlist, dynamics, &
                  i, count_qm, dynvars, enefunc%qmmm)
        end if
      end if

      ! Output energy(t + dt) and dynamical variables(t + dt)
      !   coord     is at t +   2dt, coord_ref    is at t +    dt
      !   vel       is at t + 3/2dt, vel_ref      is at t + 1/2dt
      !   box_size  is at t +   2dt, box_size_ref is at t +    dt
      !
      call output_md(output, molecule, enefunc, dynamics, boundary, &
                     ensemble, dynvars)


    end do
    
    return

  end subroutine leapfrog_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_leapfrog
  !> @brief        compute the first step (0+dt)
  !! @authors      TM, CK, KY
  !! @param[in]    output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynamics    : dynamics information
  !! @param[in]    pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] constraints : bond constraint information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine initial_leapfrog(output, molecule, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, dynvars)

    ! formal arguments
    type(s_output),           intent(in)    :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynamics), target, intent(inout) :: dynamics
    type(s_pairlist),         intent(in)    :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_ensemble),         intent(inout) :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: friction, temperature, dt, simtim
    integer                   :: j, natom
    character(256)            :: folder,basename
    character                 :: num*9
    logical                   :: savefile

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_pbc(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: force(:,:), mass(:), inv_mass(:)
    real(wp),         pointer :: force_omp(:,:,:)
    real(wp),         pointer :: viri(:,:), viri_const(:,:)
    integer,          pointer :: iseed


    ! use pointers
    natom       =  molecule%num_atoms
    coord       => dynvars%coord
    coord_pbc   => dynvars%coord_pbc
    coord_ref   => dynvars%coord_ref
    force       => dynvars%force
    force_omp   => dynvars%force_omp
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    viri        => dynvars%virial
    viri_const  => dynvars%virial_const
    inv_mass    => molecule%inv_mass
    mass        => molecule%mass
    temperature =  ensemble%temperature
    iseed       => dynamics%iseed
    dt          =  dynamics%timestep/AKMA_PS
    friction    =  ensemble%gamma_t*AKMA_PS
    simtim      = dynamics%initial_time

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


    ! Calculate velocities(0 - dt/2)
    !
    if (.not. constraints%rigid_bond) then
      
      do j = 1, natom
        vel(1,j) = vel(1,j) - 0.5_wp*dt*force(1,j)*inv_mass(j)
        vel(2,j) = vel(2,j) - 0.5_wp*dt*force(2,j)*inv_mass(j)
        vel(3,j) = vel(3,j) - 0.5_wp*dt*force(3,j)*inv_mass(j)
      end do

    else

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

      do j = 1, natom
        vel_ref(1,j) = (coord_ref(1,j) - coord(1,j))/dt
        vel_ref(2,j) = (coord_ref(2,j) - coord(2,j))/dt
        vel_ref(3,j) = (coord_ref(3,j) - coord(3,j))/dt
      end do


      ! vel <= vel(0 - dt/2) and coord <= coord(0)
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


    ! first step dynamics
    !   Calculate coordinates(0 + dt) and velocities(0 + 1/2dt)
    !   Note: this routine corresponds to i = 0 in the main loop
    !

    ! coord_ref <= coord(0) and vel_ref <= vel(0 - 1/2dt)
    !
    do j = 1, natom
      coord_ref(1,j) = coord(1,j)
      coord_ref(2,j) = coord(2,j)
      coord_ref(3,j) = coord(3,j)
      vel_ref(1,j)   = vel(1,j)
      vel_ref(2,j)   = vel(2,j)
      vel_ref(3,j)   = vel(3,j)
    end do


    ! calculate potential energy(0) and force(0)
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

    ! calculate v(0+1/2dt) and r(0+dt)
    !   v(0+1/2dt) = v(0-1/2dt) + dt*F(0)/m
    !   r(0+dt) = r(0) + dt*v(0+1/2dt)
    !
    do j = 1, natom
      vel(1,j) = vel(1,j) + dt*force(1,j)*inv_mass(j)
      vel(2,j) = vel(2,j) + dt*force(2,j)*inv_mass(j)
      vel(3,j) = vel(3,j) + dt*force(3,j)*inv_mass(j)

      coord(1,j) = coord(1,j) + dt*vel(1,j)
      coord(2,j) = coord(2,j) + dt*vel(2,j)
      coord(3,j) = coord(3,j) + dt*vel(3,j)
    end do

    ! SHAKE and SETTLE
    !   calculate constrained coordinates(0 + dt) and constraint virial(0)
    !   update velocities(0 + 1/2dt):
    !     v(0+1/2dt) = (r(0+dt) - r(0))/dt
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)

      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))/dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))/dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))/dt
      end do

      viri(1:3,1:3) = viri(1:3,1:3) + viri_const(1:3,1:3)

    end if


    ! output dynvars(0)
    !
    call compute_dynvars(molecule, enefunc, dynamics, boundary, ensemble, dynvars)
    call output_dynvars (output, enefunc, dynvars, ensemble)

    dynamics%restart = .true.

    ! at this point
    !   coord is at 0 + dt, and vel is at 0 + 1/2dt


    return

  end subroutine initial_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_temp_pres_leap
  !> @brief        driver to control temperature and pressure
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_temp_pres_leap(molecule, dynamics, ensemble, constraints, &
                                    boundary, dynvars)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars
    

    if (ensemble%tpcontrol == TpcontrolBerendsen) then

      call berendsen_leapfrog  (molecule, dynamics, ensemble, constraints, &
                                boundary, dynvars)

    else if (ensemble%tpcontrol == TpcontrolNoseHoover) then

      call nose_hoover_leapfrog(molecule, dynamics, ensemble, constraints, &
                                dynvars)

    else if (ensemble%tpcontrol == TpcontrolGauss) then

      call gaussian_leapfrog   (molecule, dynamics, ensemble, constraints, &
                                dynvars)

    end if

    return

  end subroutine control_temp_pres_leap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    berendsen_leapfrog
  !> @brief        control temperature and pressure using Berendsen method
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         H.J.C.Berendsen et al., J.Chem.Phys., 81, 3684-3690 (1984)
  !!               Y.Zhang et al., J.Chem.Phys., 103, 10252-10266 (1995)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine berendsen_leapfrog(molecule, dynamics, ensemble, constraints, &
                                boundary, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0, press0, gamma0, pressxy0
    real(wp)                  :: ke(3), tau_t, tau_p, compress
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: tvel(3), ekin, tempt
    real(wp)                  :: scale_v, scale_c(3), factor
    integer                   :: j
    integer                   :: num_degree, natom

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), mass(:), virial(:,:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P
    gamma0     =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    tau_t      =  ensemble%tau_t/AKMA_PS
    tau_p      =  ensemble%tau_p/AKMA_PS
    compress   =  ensemble%compressibility/ATMOS_P
    num_degree =  molecule%num_deg_freedom
    mass       => molecule%mass
    natom      =  molecule%num_atoms
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    force      => dynvars%force
    virial     => dynvars%virial


    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z


    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    ke(1:3) = 0.0_wp
    do j = 1, natom
      tvel(1) = 0.5_wp*(vel(1,j) + vel_ref(1,j))
      tvel(2) = 0.5_wp*(vel(2,j) + vel_ref(2,j))
      tvel(3) = 0.5_wp*(vel(3,j) + vel_ref(3,j))
      ke(1) = ke(1) + mass(j)*tvel(1)*tvel(1)
      ke(2) = ke(2) + mass(j)*tvel(2)*tvel(2)
      ke(3) = ke(3) + mass(j)*tvel(3)*tvel(3)
    end do
    ekin   = 0.5_wp*(ke(1) + ke(2) + ke(3))
    tempt  = 2.0_wp*ekin/(num_degree*KBOLTZ)


    ! compute pressure(t+dt)
    !
    pressx = (ke(1) + virial(1,1))/volume
    pressy = (ke(2) + virial(2,2))/volume
    pressz = (ke(3) + virial(3,3))/volume
    pressxyz = (pressx + pressy + pressz)/3.0_wp
    pressxy  = (pressx + pressy)/2.0_wp


    ! thermostat
    !   scale_v    = sqrt(1+dt/tau_t(temp_ref/temp(t+dt) -1))
    !   v(t+3/2dt) = (v(t+1/2dt) + dtF(t+dt)/m) * scale_v
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    scale_v = sqrt(1.0_wp + dt * (temp0/tempt - 1.0_wp)/tau_t)

    do j = 1, natom
      factor   = dt/mass(j)
      vel(1,j) = (vel_ref(1,j) + factor*force(1,j))*scale_v
      vel(2,j) = (vel_ref(2,j) + factor*force(2,j))*scale_v
      vel(3,j) = (vel_ref(3,j) + factor*force(3,j))*scale_v

      coord(1,j) = coord_ref(1,j) + vel(1,j)*dt
      coord(2,j) = coord_ref(2,j) + vel(2,j)*dt
      coord(3,j) = coord_ref(3,j) + vel(3,j)*dt
    end do


    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)

      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))*inv_dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))*inv_dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))*inv_dt
      end do

    end if


    if (ensemble%use_barostat) then

      ! barostat
      !   scale coordinates(t+2dt) and box_size(t+2dt)
      !   scale_c  = (1- compress*dt*(Pext - P(t+dt)/tau_p))**1/3
      !   r(t+2dt)    = scale_c * r(t+2dt)
      !   size(t+2dt) = scale_c * size(t+dt)
      !
      if (ensemble%isotropy == IsotropyISO) then
        scale_c(1) = (1.0_wp - compress*dt*(press0 - pressxyz)/tau_p)**ONE_THIRD
        scale_c(2) = scale_c(1)
        scale_c(3) = scale_c(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        if (ensemble%ensemble == EnsembleNPT) then
          scale_c(1) = (1.0_wp - compress*dt*(press0 - pressxy)/tau_p)**ONE_THIRD
          scale_c(2) = scale_c(1)
          scale_c(3) = (1.0_wp - compress*dt*(press0 - pressz) /tau_p)**ONE_THIRD

        else if (ensemble%ensemble == EnsembleNPgT) then
          pressxy0   = press0 - gamma0 / boundary%box_size_z
          scale_c(1) = (1.0_wp - compress*dt*(pressxy0 - pressxy)/tau_p)**ONE_THIRD
          scale_c(2) = scale_c(1)
          scale_c(3) = (1.0_wp - compress*dt*(press0   - pressz )/tau_p)**ONE_THIRD

        end if

      else if (ensemble%isotropy == IsotropyANISO) then

        if (ensemble%ensemble == EnsembleNPT) then
          scale_c(1) = (1.0_wp - compress*dt*(press0 - pressx)/tau_p)**ONE_THIRD
          scale_c(2) = (1.0_wp - compress*dt*(press0 - pressy)/tau_p)**ONE_THIRD
          scale_c(3) = (1.0_wp - compress*dt*(press0 - pressz)/tau_p)**ONE_THIRD

        else if (ensemble%ensemble == EnsembleNPgT) then
          pressxy0   = press0 - gamma0 / boundary%box_size_z
          scale_c(1) = (1.0_wp - compress*dt*(pressxy0 - pressx)/tau_p)**ONE_THIRD
          scale_c(2) = (1.0_wp - compress*dt*(pressxy0 - pressy)/tau_p)**ONE_THIRD
          scale_c(3) = (1.0_wp - compress*dt*(press0   - pressz)/tau_p)**ONE_THIRD

        end if

      else if (ensemble%isotropy == IsotropyXY_Fixed) then
        scale_c(1) = 1.0_wp
        scale_c(2) = 1.0_wp
        scale_c(3) = (1.0_wp - compress*dt*(press0 - pressz)/tau_p)**ONE_THIRD

      end if

      do j = 1, natom
        coord(1,j) = scale_c(1) * coord(1,j)
        coord(2,j) = scale_c(2) * coord(2,j)
        coord(3,j) = scale_c(3) * coord(3,j)
      end do

      boundary%box_size_x = scale_c(1) * boundary%box_size_x_ref
      boundary%box_size_y = scale_c(2) * boundary%box_size_y_ref
      boundary%box_size_z = scale_c(3) * boundary%box_size_z_ref

    end if

    return

  end subroutine berendsen_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nose_hoover_leapfrog
  !> @brief        control temperature and pressure using Nose-Hoover method
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         S.Nose, J.Chem.Phys., 81, 511-519 (1984)
  !!               S.Nose, Mol.Phys., 52, 255-268 (1984)
  !!               W.G.Hoover, Phys.Rev.A, 31, 1695-1697 (1985)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine nose_hoover_leapfrog(molecule, dynamics, ensemble, constraints, &
                                  dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0, press0, tau_t
    real(wp)                  :: ke(3), tvel(3), ekin
    real(wp)                  :: factor
    real(wp)                  :: sigma, qmass
    integer                   :: j, num_degree, natom

    real(wp),         pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), virial(:,:)
    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: moment, moment_ref


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure
    tau_t      =  ensemble%tau_t/AKMA_PS
    natom      =  molecule%num_atoms
    num_degree =  molecule%num_deg_freedom
    mass       => molecule%mass
    inv_mass   => molecule%inv_mass
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    force      => dynvars%force
    virial     => dynvars%virial
    moment_ref => dynvars%thermostat_momentum_ref
    moment     => dynvars%thermostat_momentum


    ! save momentum(t+1/2dt)
    !
    dynvars%thermostat_momentum_ref = dynvars%thermostat_momentum 


    ! compute kinetic energy(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    ke(1:3) = 0.0_wp
    do j = 1, natom
      tvel(1) = 0.5_wp*(vel(1,j) + vel_ref(1,j))
      tvel(2) = 0.5_wp*(vel(2,j) + vel_ref(2,j))
      tvel(3) = 0.5_wp*(vel(3,j) + vel_ref(3,j))
      ke(1) = ke(1) + mass(j)*tvel(1)*tvel(1)
      ke(2) = ke(2) + mass(j)*tvel(2)*tvel(2) 
      ke(3) = ke(3) + mass(j)*tvel(3)*tvel(3)
    end do
    ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))


    ! thermostat
    !   moment(t+3/2dt) = moment(t+1/2dt) + dt*[2Ekin(t+dt)-2s]/qmass
    !   moment(t+dt)    = [moment(t+3/2dt) + moment(t+1/2dt)]/2
    !   v(t+dt)    = [v(t+3/2dt) + v(t+1/2dt)]/2
    !   v(t+3/2dt) = v(t+1/2dt) + dt*[F(t+dt)/m - moment(t+dt)*v(t+dt)]
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    sigma  = 0.5_wp*num_degree*KBOLTZ*temp0
    qmass  = 2.0_wp*sigma*tau_t*tau_t
    moment = moment_ref + dt*(2.0_wp*ekin - 2.0_wp*sigma)/qmass

    factor = 0.5_wp*(moment_ref + moment)
    do j = 1, natom
      tvel(1) = 0.5_wp*(vel(1,j) + vel_ref(1,j))
      tvel(2) = 0.5_wp*(vel(2,j) + vel_ref(2,j))
      tvel(3) = 0.5_wp*(vel(3,j) + vel_ref(3,j))

      vel(1,j) = vel_ref(1,j) + (force(1,j)*inv_mass(j) - factor*tvel(1))*dt
      vel(2,j) = vel_ref(2,j) + (force(2,j)*inv_mass(j) - factor*tvel(2))*dt
      vel(3,j) = vel_ref(3,j) + (force(3,j)*inv_mass(j) - factor*tvel(3))*dt

      coord(1,j) = coord_ref(1,j) + vel(1,j)*dt
      coord(2,j) = coord_ref(2,j) + vel(2,j)*dt
      coord(3,j) = coord_ref(3,j) + vel(3,j)*dt
    end do


    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)

      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))*inv_dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))*inv_dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))*inv_dt
      end do

    end if

    return

  end subroutine nose_hoover_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gaussian_leapfrog
  !> @brief        control temperature using Gaussian thermostat method
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         D.Brown & J.H.R.Clarke, Mol.Phys., 51, 1243-1252 (1984)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gaussian_leapfrog(molecule, dynamics, ensemble, constraints, &
                               dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, ekin, tempt, beta
    real(wp)                  :: factor1, factor2, factor
    real(wp)                  :: tvel(3), ke(3), temp0
    integer                   :: j, num_degree, natom

    real(wp),         pointer :: mass(:), coord(:,:), coord_ref(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    num_degree =  molecule%num_deg_freedom
    natom      =  molecule%num_atoms
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    force      => dynvars%force
    mass       => molecule%mass
    temp0      =  ensemble%temperature


    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    ke(1:3) = 0.0_wp
    do j = 1, natom
      tvel(1) = 0.5_wp*(vel(1,j) + vel_ref(1,j))
      tvel(2) = 0.5_wp*(vel(2,j) + vel_ref(2,j))
      tvel(3) = 0.5_wp*(vel(3,j) + vel_ref(3,j))
      ke(1) = ke(1) + mass(j)*tvel(1)*tvel(1)
      ke(2) = ke(2) + mass(j)*tvel(2)*tvel(2)
      ke(3) = ke(3) + mass(j)*tvel(3)*tvel(3)
    end do
    ekin  = 0.5_wp*(ke(1) + ke(2) + ke(3))
    tempt = 2.0_wp*ekin/(num_degree*KBOLTZ)


    ! calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
    !   v(t+3/2dt) = v(t+1/2dt)*(2*beta-1) + beta*dt*F(t+dt)/m
    !   r(t+2dt) = r(t+dt) + dt*vel(t+3/2dt)
    !
    beta    = sqrt(temp0/tempt)
    factor  = beta*dt
    factor1 = 2.0_wp*beta - 1.0_wp

    do j = 1, natom
      factor2  = factor/mass(j)
      vel(1,j) = factor1*vel_ref(1,j) + factor2*force(1,j)
      vel(2,j) = factor1*vel_ref(2,j) + factor2*force(2,j)
      vel(3,j) = factor1*vel_ref(3,j) + factor2*force(3,j)

      coord(1,j) = coord_ref(1,j) + dt*vel(1,j)
      coord(2,j) = coord_ref(2,j) + dt*vel(2,j)
      coord(3,j) = coord_ref(3,j) + dt*vel(3,j)
    end do


    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)

      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))*inv_dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))*inv_dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))*inv_dt
      end do

    end if


    return

  end subroutine gaussian_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_nvt
  !> @brief        control temperature using Langevin
  !! @authors      TM, KY
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         A.Brunger et al., Chem.Phys.Lett., 105, 495-500 (1984)
  !!               R.Pastor et al., Mol.Phys., 65, 1409-1419 (1988)
  !!               R.J.Loncharich et al., Biopolymers, 32, 523-535 (1992)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_leapfrog_nvt(molecule, dynamics, ensemble, boundary,  &
                                   constraints, dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_boundary),         intent(in)    :: boundary
    type(s_constraints), target, intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, temp0, gamma_t, inv_dt
    real(wp)                  :: scale_v, scale_f, factor
    real(wp)                  :: sigma
    real(wp)                  :: v1, v2, rsq, grandom(3)
    integer                   :: j, natom

    real(wp),         pointer :: coord(:,:), force(:,:), coord_ref(:,:)
    real(wp),         pointer :: vel(:,:), viri(:,:), viri_const(:,:)
    real(wp),         pointer :: mass(:), inv_mass(:)
    logical,          pointer :: fixatm(:)


    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    temp0      =  ensemble%temperature
    natom      =  molecule%num_atoms
    mass       => molecule%mass
    inv_mass   => molecule%inv_mass
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    vel        => dynvars%velocity
    force      => dynvars%force
    viri       => dynvars%virial
    viri_const => dynvars%virial_const
    fixatm     => constraints%fixatm

    gamma_t    = ensemble%gamma_t*AKMA_PS


    ! add random force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
    !   TM notes that the function gauss() should not be used in this loop
    !
    factor = 2.0_wp*gamma_t*KBOLTZ*temp0/dt

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

      force(1,j) = force(1,j) + sigma * grandom(1) 
      force(2,j) = force(2,j) + sigma * grandom(2) 
      force(3,j) = force(3,j) + sigma * grandom(3) 
    end do


    ! Langevin dynamics (Langevin thermostat)
    !   calculate velocities(t + 3/2dt)
    !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
    !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
    !
    factor  = 1.0_wp + 0.5_wp * dt * gamma_t
    scale_v = (2.0_wp/factor) - 1.0_wp
    scale_f = dt/factor

    !$omp parallel do default(none)           &
    !$omp private(j)                          &
    !$omp shared(natom, dt, scale_v, scale_f, &
    !$omp        vel, inv_mass, force, coord)
    !
    do j = 1, natom
      vel(1,j) = scale_v*vel(1,j) + scale_f*force(1,j)*inv_mass(j)
      vel(2,j) = scale_v*vel(2,j) + scale_f*force(2,j)*inv_mass(j)
      vel(3,j) = scale_v*vel(3,j) + scale_f*force(3,j)*inv_mass(j)

      coord(1,j) = coord(1,j) + dt*vel(1,j)
      coord(2,j) = coord(2,j) + dt*vel(2,j)
      coord(3,j) = coord(3,j) + dt*vel(3,j)
    end do
    !$omp end parallel do


    ! Bond constraint
    !   coord (unconstrained) is at t + 2dt, vel     is at t + 3/2dt
    !   coord_ref             is at t +  dt, vel_ref is at t + 1/2dt
    !   compute constrained coordinates(t+2dt) and constraint virial(t+dt)
    !   update velocities(t + 3/2dt):
    !     v(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !   add constraint virial(t + dt)
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, &
                               dt, molecule, dynvars, constraints)

      !$omp parallel do default(none)                   &
      !$omp private(j)                                  &
      !$omp shared(natom, vel, coord, coord_ref, inv_dt)
      !
      do j = 1, natom
        vel(1,j) = (coord(1,j) - coord_ref(1,j))*inv_dt
        vel(2,j) = (coord(2,j) - coord_ref(2,j))*inv_dt
        vel(3,j) = (coord(3,j) - coord_ref(3,j))*inv_dt
      end do
      !$omp end parallel do

      viri(1:3,1:3) = viri(1:3,1:3) + viri_const(1:3,1:3)

    end if

    return

  end subroutine langevin_leapfrog_nvt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_npt
  !> @brief        Langevin thermostat and barostat
  !! @authors      TM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] boundary    : boundary condition information
  !! @param[inout] dynvars     : dynamic variables information
  !! @note         A.Brunger et al., Chem.Phys.Lett., 105, 495-500 (1984)
  !!               R.Pastor et al., Mol.Phys., 65, 1409-1419 (1988)
  !!               R.J.Loncharich et al., Biopolymers, 32, 523-535 (1992)
  !!               G.J.Martyna et al., J.Chem.Phys., 101, 4177-4189 (1994)
  !!               S.E.Feller et al., J.Chem.Phys., 103, 4613-4621 (1995)
  !!               Y.Zhang et al., J.Chem.Phys., 103, 10252-10266 (1995)
  !!               D.Quigley & M.I.J.Probert, J.Chem.Phys., 120, 11432 (2004)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_leapfrog_npt(molecule, dynamics, ensemble, constraints, &
                                   boundary, dynvars)
    
    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars
    
    ! local variables
    real(wp)                  :: dt, inv_dt, temp0, press0, gamma0, pressxy0
    real(wp)                  :: d_ndegf, ke(3), ekin
    real(wp)                  :: volume, pressx, pressy, pressz
    real(wp)                  :: pressxy, pressxyz
    real(wp)                  :: tvel(3), tcrd(3)
    real(wp)                  :: scale_f(3), scale_v(3), fact(3), factor
    real(wp)                  :: bmoment_ref(3), bmoment2(3)
    real(wp)                  :: gamma_t, gamma_p, pmass, pforce(3)
    real(wp)                  :: sigma
    real(wp)                  :: v1, v2, rsq, grandom(3)
    integer                   :: i, j, i_ndegf, natom

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), coord_old(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: bmoment(:)


    ! use pointers 
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P 
    gamma0     =  ensemble%gamma    * ATMOS_P*100.0_wp/1.01325_wp
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    gamma_p    =  ensemble%gamma_p * AKMA_PS
    i_ndegf    =  molecule%num_deg_freedom
    d_ndegf    =  real(i_ndegf,wp)
    mass       => molecule%mass
    inv_mass   => molecule%inv_mass
    natom      =  molecule%num_atoms
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    vel        => dynvars%velocity
    vel_ref    => dynvars%velocity_ref
    force      => dynvars%force
    coord_old  => dynvars%temporary
    virial     => dynvars%virial 
    viri_const => dynvars%virial_const
    bmoment    => dynvars%barostat_momentum


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
    !$omp parallel do default(none)                      &
    !$omp shared(natom, vel, force, inv_mass, coord, dt) &
    !$omp private(j) 
    !
    do j = 1, natom
      vel(1,j) = vel(1,j) + dt*force(1,j)*inv_mass(j)
      vel(2,j) = vel(2,j) + dt*force(2,j)*inv_mass(j)
      vel(3,j) = vel(3,j) + dt*force(3,j)*inv_mass(j)

      coord(1,j) = coord(1,j) + dt*vel(1,j)
      coord(2,j) = coord(2,j) + dt*vel(2,j)
      coord(3,j) = coord(3,j) + dt*vel(3,j)
    end do
    !$omp end parallel do


    ! add random_force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
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

      force(1,j) = force(1,j) + sigma * grandom(1)
      force(2,j) = force(2,j) + sigma * grandom(2)
      force(3,j) = force(3,j) + sigma * grandom(3)
    end do


    ! calculate piston mass (pmass) and stochastic force (pforce)
    ! acting on the barostat
    !
    if (ensemble%isotropy == IsotropyISO) then
      pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0 / (2.0_wp*PI*gamma_p)**2
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = pforce(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0 / (3.0_wp*(2.0_wp*PI*gamma_p)**2)
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyANISO) then
      pmass = real(i_ndegf+3,wp)*KBOLTZ*temp0 / (3.0_wp*(2.0_wp*PI*gamma_p)**2)
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = sigma * random_get_gauss()
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      pmass = real(i_ndegf+1,wp)*KBOLTZ*temp0 / (3.0_wp*(2.0_wp*PI*gamma_p)**2)
      sigma = sqrt(2.0_wp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = 0.0_wp
      pforce(2) = 0.0_wp
      pforce(3) = sigma * random_get_gauss()

    end if


    ! iteration
    !   Since leapfrog algorithm cannot directly compute vel(t+dt),
    !   iteration scheme is required to obtain T(t+dt) and P(t+dt).
    !
    do i = 1, 8

      ! compute kinetic energy(t+dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      ke(1:3) = 0.0_wp
      do j = 1, natom
        ke(1) = ke(1) + mass(j)*vel(1,j)*vel(1,j)
        ke(1) = ke(1) + mass(j)*vel_ref(1,j)*vel_ref(1,j)
        ke(2) = ke(2) + mass(j)*vel(2,j)*vel(2,j)
        ke(2) = ke(2) + mass(j)*vel_ref(2,j)*vel_ref(2,j)
        ke(3) = ke(3) + mass(j)*vel(3,j)*vel(3,j)
        ke(3) = ke(3) + mass(j)*vel_ref(3,j)*vel_ref(3,j)
      end do
      ke(1:3) = ke(1:3) / 2.0_dp

      ! scale factor for harmonic oscillators (R.Pastor et al., Mol.Phys. 1988)
      if (i /= 1) then
        ke(1) = ke(1)*(1.0_wp + 0.5_wp*gamma_t*dt)
        ke(2) = ke(2)*(1.0_wp + 0.5_wp*gamma_t*dt)
        ke(3) = ke(3)*(1.0_wp + 0.5_wp*gamma_t*dt)
      end if
      ekin = 0.5_wp*(ke(1) + ke(2) + ke(3))


      ! compute pressure(t+dt) in the unit of kcal/mol*A3
      !
      pressx = (ke(1) + virial(1,1))/volume
      pressy = (ke(2) + virial(2,2))/volume
      pressz = (ke(3) + virial(3,3))/volume
      pressxyz = (pressx + pressy + pressz)/3.0_wp
      pressxy  = (pressx + pressy)/2.0_wp


      ! compute thermostat and barostat parameters
      !   for isotropic systems
      !     eta(t+3/2dt) = eta(t+1/2dt) + dt[dim*V(t+dt)*(P(t+dt)-Pext) 
      !                                     + dim*2Ekin(t+dt)/Nf +Rp]/pmass
      !     eta(t+3/2dt) = eps(-gamma_p*dt)*eta(t+3/2dt)
      !     eta(t+dt)    = 0.5*(eta(t+3/2dt) + eta(t+1/2dt))
      !     factor       = 1 + [gamma_t+(1+dim/Nf)*eta(t+dt)]dt/2
      !
      if (ensemble%isotropy == IsotropyISO) then

        bmoment(1) = bmoment_ref(1) + dt*(3.0_wp*volume*(pressxyz - press0)   &
                                    + 6.0_wp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(1)

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wp*(bmoment(1:3) + bmoment_ref(1:3))

        fact(1:3) = 1.0_wp+(gamma_t+bmoment2(1:3)*(1.0_wp+3.0_wp/d_ndegf))*0.5_wp*dt

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1) = bmoment_ref(1) + dt*(volume*(pressxy - press0)          &
                                      + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
          bmoment(2) = bmoment_ref(2) + dt*(volume*(pressxy - press0)          &
                                      + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(pressz  - press0)          &
                                      + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1) = bmoment_ref(1) + dt*(volume*(pressxy - pressxy0)        &
                                      + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
          bmoment(2) = bmoment_ref(2) + dt*(volume*(pressxy - pressxy0)        &
                                      + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(pressz  - press0)          &
                                      + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wp*dt

      else if (ensemble%isotropy == IsotropyANISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1) = bmoment_ref(1) + dt*(volume*(pressx - press0)           &
                                      + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
          bmoment(2) = bmoment_ref(2) + dt*(volume*(pressy - press0)           &
                                      + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(pressz - press0)           &
                                      + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1) = bmoment_ref(1) + dt*(volume*(pressx - pressxy0)         &
                                      + 2.0_wp*ekin/d_ndegf + pforce(1))/pmass
          bmoment(2) = bmoment_ref(2) + dt*(volume*(pressy - pressxy0)         &
                                      + 2.0_wp*ekin/d_ndegf + pforce(2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(pressz - press0)           &
                                      + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wp*dt

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        bmoment(1) = 0.0_wp
        bmoment(2) = 0.0_wp
        bmoment(3) = bmoment_ref(3) + dt*(volume*(pressz - press0)             &
                                    + 2.0_wp*ekin/d_ndegf + pforce(3))/pmass

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_wp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_wp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_wp*dt

      end if

      scale_v(1:3) = 2.0_wp/fact(1:3) - 1.0_wp
      scale_f(1:3) = dt/fact(1:3)


      ! Langevin dynamics
      !   calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
      !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
      !   r(t+3/2dt) = 0.5*(r(t+2dt) + r(t+dt))
      !   r(t+2dt)   = r(t+dt) + dt*(v(t+3/2dt) + eta(t+3/2dt)*r(t+3/2dt)))
      !
      !$omp parallel do default(none)                         &
      !$omp private(j, tcrd)                                  &
      !$omp shared(natom, dt, vel, vel_ref, scale_v, scale_f, &
      !$omp        inv_mass, force, coord, coord_ref, bmoment)
      !
      do j = 1, natom
        vel(1,j) = scale_v(1)*vel_ref(1,j) + scale_f(1)*force(1,j)*inv_mass(j)
        vel(2,j) = scale_v(2)*vel_ref(2,j) + scale_f(2)*force(2,j)*inv_mass(j)
        vel(3,j) = scale_v(3)*vel_ref(3,j) + scale_f(3)*force(3,j)*inv_mass(j)

        tcrd(1)  = 0.5_wp*(coord_ref(1,j) + coord(1,j))
        tcrd(2)  = 0.5_wp*(coord_ref(2,j) + coord(2,j))
        tcrd(3)  = 0.5_wp*(coord_ref(3,j) + coord(3,j))

        coord(1,j) = coord_ref(1,j) + dt*(vel(1,j) + bmoment(1)*tcrd(1))
        coord(2,j) = coord_ref(2,j) + dt*(vel(2,j) + bmoment(2)*tcrd(2))
        coord(3,j) = coord_ref(3,j) + dt*(vel(3,j) + bmoment(3)*tcrd(3))
      end do
      !$omp end parallel do


      ! Bond constraint
      !   store unconstrained coord(t + 2dt)
      !   compute constrained coordinates(t + 2dt)
      !   add contribution of constraints to virial
      !   correct velocities(t + 3/2dt) and forces(t + dt)
      !
      if (constraints%rigid_bond) then

        do j = 1, natom
          coord_old(1,j) = coord(1,j)
          coord_old(2,j) = coord(2,j)
          coord_old(3,j) = coord(3,j)
        end do

        call compute_constraints(ConstraintModeLEAP, &
                                 dt, molecule, dynvars, constraints)

        virial(1:3,1:3) = virial(1:3,1:3) + viri_const(1:3,1:3)

        !$omp parallel do default(none)                                 &
        !$omp private(j, tvel)                                          &
        !$omp shared(natom, coord, coord_old, vel, force, mass, inv_dt)
        !
        do j = 1, natom
          tvel(1)  = (coord(1,j) - coord_old(1,j))*inv_dt
          tvel(2)  = (coord(2,j) - coord_old(2,j))*inv_dt
          tvel(3)  = (coord(3,j) - coord_old(3,j))*inv_dt

          vel(1,j) = vel(1,j) + tvel(1)
          vel(2,j) = vel(2,j) + tvel(2)
          vel(3,j) = vel(3,j) + tvel(3)

          force(1,j) = force(1,j) + mass(j)*tvel(1)*inv_dt
          force(2,j) = force(2,j) + mass(j)*tvel(2)*inv_dt
          force(3,j) = force(3,j) + mass(j)*tvel(3)*inv_dt
        end do
        !$omp end parallel do

      end if

    end do


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

  end subroutine langevin_leapfrog_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_leapfrog
  !> @brief        change target temperature linearly
  !! @authors      TM
  !! @param[in]    dynamics : dynamics information
  !! @param[inout] ensemble : ensemble information
  !! @note         S.Kirkpatrick et al., Science, 220, 671-680 (1983)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine simulated_annealing_leapfrog(dynamics, enefunc, ensemble)

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
        'Simulated_Annealing_Leapfrog> Temperature is less than 0 K')

    if (main_rank) &
      write(MsgOut,'(A,F10.3,A,F10.3)')                              &
            'Simulated_Annealing_Leapfrog> Anneal temperature from', &
            old_temperature, ' to ', ensemble%temperature

    return

  end subroutine simulated_annealing_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_gamd_leapfrog
  !> @brief        update GaMD parameters
  !! @authors      HO
  !! @param[inout] output  : output information
  !! @param[inout] enefunc : potential energy functions information
  !! @param[inout] dynvars : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine update_gamd_leapfrog(output, enefunc, dynvars)

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

  end subroutine update_gamd_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_qm_charge_leapfrog
  !> @brief        Perform QM/MM to obtain QM charges
  !! @authors      KY
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[in]    molecule  : molecular information
  !! @param[in]    pairlist  : information of nonbonded pair list
  !! @param[in]    dynamics  : dynamics information
  !! @param[in]    step      : step number
  !! @param[in]    ncount    : number of QM calc
  !! @param[in]    dynvars   : dynamic variables information
  !! @param[inout] qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine update_qm_charge_leapfrog(enefunc, molecule, pairlist, dynamics, &
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

    if (dynamics%avg_qm_charge) then
      qmmm%qm_charge = (qmmm%qm_charge_save*real(ncount) + qmmm%qm_charge) &
                         /real(ncount+1)
    end if

    return

  end subroutine update_qm_charge_leapfrog

end module at_md_leapfrog_mod
