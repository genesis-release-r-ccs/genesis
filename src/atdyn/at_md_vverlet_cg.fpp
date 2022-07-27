!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_md_vverlet_cg_mod
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

module at_md_vverlet_cg_mod

  use at_output_mod
  use at_dynvars_mod
  use at_ensemble_mod
  use at_constraints_mod
  use at_boundary_mod
  use at_pairlist_mod
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
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: vverlet_dynamics_cg
  private :: initial_vverlet
  private :: control_temp_pres_vver1
  private :: langevin_thermostat_vv1
  private :: langevin_thermostat_vv2
  private :: langevin_barostat_vv1
  private :: langevin_barostat_vv2
  private :: update_barostat
  private :: simulated_annealing_vverlet

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

  subroutine vverlet_dynamics_cg(output, molecule, enefunc, dynvars, dynamics, &
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

    if (ensemble%tpcontrol /= TpcontrolLangevin) &
      call error_msg( &
        'Vverlet_dynamics_cg> VVER_CG is available only with Langevin')

    if (dynamics%target_md) enefunc%rmsd_force = 1.0_wp / (dt*dt)
    if (abs(dynamics%initial_value) .lt. 0.001_wp)  then
      if (dynamics%target_md) &
      dynamics%initial_value =  &
        dynvars%energy%restraint_cv(enefunc%target_function)
      if (dynamics%steered_md) &
      dynamics%initial_value =  &
        dynvars%energy%restraint_cv(enefunc%steered_function)
    end if
    enefunc%target_value = dynamics%initial_value 

    ! first-step MD
    !
    if (.not. dynamics%restart) then
      call initial_vverlet(output, molecule, enefunc, dynamics, pairlist, &
                           boundary, ensemble, constraints, dynvars)
    end if

    !
    ! stop NPAT & NPgT
    ! 
    if (ensemble%ensemble == EnsembleNPAT .or.  &
        ensemble%ensemble == EnsembleNPT  .or.  &
        ensemble%ensemble == EnsembleNPgT)      & 
      call error_msg('Vverlet_dynamics> Barostats are not allowed in ATDYN')


    ! Reset qm_count
    if(enefunc%qmmm%do_qmmm) enefunc%qmmm%qm_count = istart

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
                                   boundary%fixatm,                   &
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

      ! Output energy(t + dt) and dynamical variables(t + dt)
      !
      call output_md(output, molecule, enefunc, dynamics, boundary, &
                     ensemble, dynvars)

    end do
    
    return

  end subroutine vverlet_dynamics_cg

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
    if (boundary%num_fixatm > 0) &
      call clear_fixatm_component(constraints, natom, dynvars%velocity)

    call stop_trans_rotation(molecule%num_atoms, molecule%mass,      &
                             dynamics%stop_com_translation,          &
                             dynamics%stop_com_rotation,             &
                             dynvars%coord,                          &
                             boundary%fixatm,                        &
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
    type(s_dynamics),         intent(inout) :: dynamics
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

    end if

    return

  end subroutine control_temp_pres_vver1

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
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_molecule), target, intent(in)    :: molecule
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_boundary), target, intent(in)    :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, h_dt
    real(wp)                  :: temperature, gamma_t, scale_v
    real(wp)                  :: factor, sigma, temp0
    real(wp)                  :: v1, v2, rsq, grandom(1:3)
    real(wp)                  :: kBT, TWO_PI
    real(wp)                  :: tvel(3)
    integer                   :: j, natom, is, ie, quotient, remainder, count

    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),         pointer :: mass(:), inv_mass(:), temporary(:,:)
    real(wp),         pointer :: viri(:,:), viri_const(:,:)
    real(wp),         pointer :: random_force(:,:), random_force1(:,:)
    logical,          pointer :: fixatm(:)
    integer,          pointer :: istart, iend


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
    random_force1 => dynvars%random_force1
    temporary     => dynvars%temporary
    viri          => dynvars%virial
    viri_const    => dynvars%virial_const
    istart        => dynvars%istart_atom
    iend          => dynvars%iend_atom
    fixatm        => boundary%fixatm

    ! setup variables
    !
    kBT      = KBOLTZ * temp0
    TWO_PI   = 2.0_wp * PI
    h_dt     = dt / 2.0_wp

    ! random force
    !
    if (dynvars%step == 1) then

      factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/h_dt

      dynamics%iseed = dynamics%iseed + my_country_rank

      if (dynamics%multiple_random_seed) then
        call random_term
        call random_init(dynamics%iseed)
      end if

      random_force(1:3,1:natom) = 0.0_wp

#ifdef HAVE_MPI_GENESIS
      allocate(dynvars%displs(0:nproc_country-1),     &
               dynvars%recv_count(0:nproc_country-1), &
               stat = ierror)

      do j = 0, nproc_country - 1

        quotient  = natom / nproc_country
        remainder = mod(natom, nproc_country)

        if (j <= (remainder - 1)) then
          quotient = quotient + 1
          is       = quotient * j + 1
          ie       = is + quotient - 1
        else
          is       = (quotient + 1)*remainder + &
                     quotient*(j - remainder) + 1
          ie       = is + quotient - 1
        end if
        dynvars%recv_count(j) = 3*(ie-is+1)
        if (my_country_rank == j) then
          istart = is
          iend   = ie
        end if

      end do

      dynvars%displs(0) = 0
      do j = 0, nproc_country - 2
        dynvars%displs(j+1) = dynvars%displs(j) + dynvars%recv_count(j)
      end do

      do j = istart, iend

        if (fixatm(j)) cycle

        sigma = sqrt(mass(j) * factor)

        grandom(1) = random_get_gauss()
        grandom(2) = random_get_gauss()
        grandom(3) = random_get_gauss()
        random_force1(1:3,j-istart+1) = sigma*grandom(1:3)

      end do
      count = dynvars%recv_count(my_country_rank)
      call mpi_allgatherv(random_force1, count, mpi_wp_real,    &
                          random_force, dynvars%recv_count,     &
                          dynvars%displs, mpi_wp_real,          &
                          mpi_comm_country, ierror)

#else
      do j = 1, natom

        if (fixatm(j)) cycle

        sigma = sqrt(mass(j) * factor)

        grandom(1) = random_get_gauss()
        grandom(2) = random_get_gauss()
        grandom(3) = random_get_gauss()
        random_force(1:3,j) = sigma*grandom(1:3)

      end do
#endif

    end if

    ! scale factor for velocities
    !
    scale_v = exp(-gamma_t * h_dt)
     
    ! thermostat
    !
    !$omp parallel do default(shared) private(j)
    do j = 1, natom
      vel(1:3,j) = vel_ref(1:3,j)*scale_v
    end do
    !$omp end parallel do

    ! VV1
    !
    !$omp parallel do default(shared) private(j)
    do j = 1, natom
      vel(1:3,j) = vel(1:3,j) + h_dt*force(1:3,j)*inv_mass(j)
      vel(1:3,j) = vel(1:3,j) + h_dt*random_force(1:3,j)*inv_mass(j)
      coord(1:3,j) = coord_ref(1:3,j) + vel(1:3,j)*dt
    end do
    !$omp end parallel do

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
    type(s_boundary), target, intent(in)    :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: dt, inv_dt, temp0
    real(wp)                  :: factor, sigma
    real(wp)                  :: gamma_t
    real(wp)                  :: v1, v2, rsq, grandom(1:3), fcm(0:4)
    real(wp)                  :: h_dt, vel_scale
    integer                   :: j, natom
    integer                   :: count

    real(wp),         pointer :: mass(:), inv_mass(:)
    real(wp),         pointer :: temporary(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: virial(:,:), viri_const(:,:)
    real(wp),         pointer :: random_f(:,:), random_f1(:,:)
    logical,          pointer :: fixatm(:)
    integer,          pointer :: istart, iend


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
    random_f1  => dynvars%random_force1
    random_f   => dynvars%random_force
    temporary  => dynvars%temporary
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    istart     => dynvars%istart_atom
    iend       => dynvars%iend_atom
    fixatm     => boundary%fixatm

    ! time step
    !
    h_dt = dt / 2.0_wp

    ! random force
    !
    fcm(0:3) = 0.0_wp
    factor   = 2.0_wp*gamma_t*KBOLTZ*temp0/dt

#ifdef HAVE_MPI_GENESIS
    do j = istart, iend
      if (fixatm(j)) cycle

      sigma = sqrt(mass(j) * factor)

      grandom(1) = random_get_gauss()
      grandom(2) = random_get_gauss()
      grandom(3) = random_get_gauss()
      random_f1(1:3,j-istart+1) = sigma*grandom(1:3)

    end do
    count = dynvars%recv_count(my_country_rank)
    call mpi_allgatherv(random_f1, count, mpi_wp_real,    &
                        random_f, dynvars%recv_count,     &
                        dynvars%displs, mpi_wp_real,      &
                        mpi_comm_country, ierror)
#else
    do j = 1, natom
      if (fixatm(j)) cycle

      sigma = sqrt(mass(j) * factor)

      grandom(1) = random_get_gauss()
      grandom(2) = random_get_gauss()
      grandom(3) = random_get_gauss()
      random_f(1:3,j) = sigma*grandom(1:3)
    end do
#endif

!   !$omp end parallel do

    ! VV2
    !
    !$omp parallel do default(shared) private(j)
    do j = 1, natom
      vel(1,j) = vel(1,j) + h_dt*(force(1,j)+random_f(1,j))*inv_mass(j)
      vel(2,j) = vel(2,j) + h_dt*(force(2,j)+random_f(2,j))*inv_mass(j)
      vel(3,j) = vel(3,j) + h_dt*(force(3,j)+random_f(3,j))*inv_mass(j)
    end do
    !$omp end parallel do

    ! termostat
    !
    vel_scale = exp(-gamma_t*h_dt)
    !$omp parallel do default(shared) private(j)
    do j = 1, natom
      vel(1,j) = vel(1,j)*vel_scale
      vel(2,j) = vel(2,j)*vel_scale
      vel(3,j) = vel(3,j)*vel_scale
    end do
    !$omp end parallel do

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
    type(s_dynamics),         intent(inout) :: dynamics
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
    maxiter = 1

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    if (dynvars%step == 1) then

      dynamics%iseed = dynamics%iseed + my_country_rank

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

    maxiter = 1

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
  !  Subroutine  simulated_annealing_vverlet
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

end module at_md_vverlet_cg_mod
