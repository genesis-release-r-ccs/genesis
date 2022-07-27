!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_md_mts_mod
!> @brief   perform molecular dynamics simulation with velocity verlet
!!          and mts algorithm
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_md_mts_mod

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
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: pverlet_mts_dynamics
  private :: initial_vverlet

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pverlet_mts_dynamics
  !> @brief        velocity verlet integrator using mts
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

  subroutine pverlet_mts_dynamics(output, domain, enefunc, dynvars, &
                                  dynamics, pairlist, boundary,     &
                                  constraints, ensemble, comm, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_dynamics),target, intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd

    ! local variables
    real(wip)                :: simtim, temperature
    real(wip)                :: dt_short, dt_long
    integer                  :: i, j, k, jx, nsteps
    integer                  :: istep, multistep
    integer                  :: iseed, alloc_stat
    logical                  :: npt

    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wip),       pointer :: coord_old(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:)
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: ncell, natom(:)
    logical,         pointer :: XI_RESPA, XO_RESPA


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_old   => domain%coord_old
    coord_pbc   => domain%translated
    force       => domain%force    
    force_omp   => domain%force_omp
    force_short => domain%force_short
    force_long  => domain%force_long
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    XI_RESPA    => dynamics%xi_respa
    XO_RESPA    => dynamics%xo_respa

    temperature =  ensemble%temperature
    nsteps      =  dynamics%nsteps
    multistep   =  dynamics%elec_long_period
    dt_short    =  dynamics%timestep/AKMA_PS
    dt_long     =  real(multistep,wip)*dt_short
    simtim      =  dynamics%initial_time
    iseed       =  dynamics%iseed_init_velocity


    npt = (ensemble%ensemble == EnsembleNPT)

    ! Open output files
    !
    call open_output(output)

    ! Restart or not
    !   restart on : dynvars has been already replaced with restart data
    !                see subroutine "setup_restart"
    !   restart off: 1. generate initial velocities(t = 0)
    !                2. remove trans and rotat motions from
    !                velocities(0)
    !                3. update the number of degrees of freedom
    !                4. compute coordinates(0 + dt) and velocities(0 +
    !                dt)
    !
    if (dynamics%restart) then

      call communicate_coor(domain, comm)

      call compute_energy_short(domain, enefunc, pairlist, boundary, coord, &
                                mod(0,dynamics%eneout_period)== 0,     &
                                dynvars%energy, &
                                domain%atmcls_pbc, &
                                coord_pbc,      &
                                force_short,    &
                                force_omp,      &
                                force_pbc,      &
                                virial_cell,    &
                                dynvars%virial, &
                                dynvars%virial_extern)

      call compute_energy_long (domain, enefunc, boundary,       &
                                npt, dynvars%energy, force_long, &
                                dynvars%virial_long)

      call communicate_force(domain, comm, force_short)

      call communicate_force(domain, comm, force_long)

    else

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

      call initial_vverlet(npt, output, enefunc, dynamics,            &
                           pairlist, boundary, ensemble, constraints, &
                           domain, dynvars, comm)

    end if

    ! Outer integrator
    !
    do i = 1, nsteps / multistep

      simtim = simtim + dt_long * AKMA_PS
      dynvars%time = simtim
      dynvars%step = i * multistep

      ! Inner integrator
      !
      do j = 1, multistep

        call timer(TimerIntegrator, TimerOn)

        ! decide the current steps
        !
        istep = (i-1) * multistep + j

        ! save coordiantes for constraint
        !
        do k = 1, ncell
          do jx = 1, natom(k)
            coord_ref(1:3,jx,k) = coord(1:3,jx,k)
          end do
        end do

        ! PV1
        !
        call integrate_pv(dt_short, domain, constraints, dynvars)

        call timer(TimerIntegrator, TimerOff)

        ! cell migration and update cell pairlist
        !
        call timer(TimerComm1, TimerOn)

        if (dynamics%nbupdate_period > 0 .and. &
            j == multistep .and. i > 0) then

          call domain_interaction_update(istep, dynamics%nbupdate_period, &
                                         domain, enefunc, pairlist,       &
                                         boundary, constraints, comm)

        end if

        call timer(TimerComm1, TimerOff)

        ! short range forces
        !
        call communicate_coor(domain, comm)

        call compute_energy_short(domain, enefunc, pairlist, boundary, coord,  &
                                  mod(istep,dynamics%eneout_period) == 0,      &
                                  dynvars%energy,                              &
                                  domain%atmcls_pbc,                           &
                                  coord_pbc,                                   &
                                  force_short,                                 &
                                  force_omp,                                   &
                                  force_pbc,                                   &
                                  virial_cell,                                 &
                                  dynvars%virial,                              &
                                  dynvars%virial_extern)

        call communicate_force(domain, comm, force_short)
        
        ! full step velocities with foce_short
        !
        call timer(TimerIntegrator, TimerOn)

        ! VV
        if (istep == 1) then
          allocate(ensemble%random_force(3,MaxAtom,ncell), stat=alloc_stat)
          if (alloc_stat /= 0) call error_msg_alloc
        end if

        call langevin_thermostat_vv(dt_short, ensemble, domain, dynvars)

        do k = 1, ncell
          do jx = 1, natom(k)
            coord_ref(1:3,jx,k) = coord(1:3,jx,k)
          end do
        end do

        ! PV2
        call integrate_pv(dt_short, domain, constraints, dynvars)

        call timer(TimerIntegrator, TimerOff)

      end do

      ! energy output
      !
      if (dynamics%eneout_period > 0) then
        if (mod(istep,dynamics%eneout_period) == 0 ) then
          do k = 1, ncell
            do jx = 1, natom(k)
              force(1:3,jx,k) = force_short(1:3,jx,k) + force_long(1:3,jx,k)
            end do
          end do
        end if
      end if

      call output_md(output, dynamics, boundary, pairlist, ensemble, &
                     constraints, dynvars, domain, enefunc, remd)

      ! long range force 
      !
      do k = 1, ncell
        do jx = 1, natom(k)
          coord(1:3,jx,k) = coord(1:3,jx,k) + 0.5_wip*dt_long*vel(1:3,jx,k)
        end do
      end do
        
      call communicate_coor(domain, comm)
      call compute_energy_long(domain, enefunc, boundary,       &
                               npt, dynvars%energy, force_long, &
                               dynvars%virial_long)
      call communicate_force(domain, comm, force_long)

      do k = 1, ncell
        do jx = 1, natom(k)
          coord(1:3,jx,k) = coord(1:3,jx,k) - 0.5_wip*dt_long*vel(1:3,jx,k)
        end do
      end do

    end do

    ! Close output files
    !
    call close_output(output)

    return

  end subroutine pverlet_mts_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_vverlet
  !> @brief        compute the first step (0+dt)
  !! @authors      JJ
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
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
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(wip)                :: temperature
    real(wip)                :: imass, simtim, dt, friction
    integer                  :: i, ix, ncell

    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(wip),       pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(wip),       pointer :: force(:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:), force_omp(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:)
    real(wip),       pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: natom(:)


    natom       => domain%num_atom
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_pbc   => domain%translated
    force       => domain%force
    force_short => domain%force_short
    force_long  => domain%force_long 
    force_omp   => domain%force_omp
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    mass        => domain%mass
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair

    dt          =  dynamics%timestep/AKMA_PS
    simtim      =  dynamics%initial_time
    temperature =  ensemble%temperature
    friction    =  ensemble%gamma_t * AKMA_PS
    ncell       =  domain%num_cell_local

    dynvars%time = simtim
    dynvars%step = 0

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
                               dynvars%virial_const)

      do i = 1, ncell
        do ix = 1, natom(i)
          coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        end do
      end do
    end if

    ! calculate energy(0) and forces(0)
    !
    call communicate_coor(domain, comm)

    call compute_energy_short(domain, enefunc, pairlist, boundary, coord, &
                              .true.,         &
                              dynvars%energy, &
                              domain%atmcls_pbc, &
                              coord_pbc,      &
                              force_short,    &
                              force_omp,      &
                              force_pbc,      &
                              virial_cell,    &
                              dynvars%virial, &
                              dynvars%virial_extern)

    call compute_energy_long (domain, enefunc, boundary,        &
                              npt, dynvars%energy, force_long,  &
                              dynvars%virial_long)

    call communicate_force(domain, comm, force_short)

    call communicate_force(domain, comm, force_long)

    do i = 1, ncell
      do ix = 1, natom(i)
        force(1:3,ix,i) = force_short(1:3,ix,i) + force_long(1:3,ix,i)
      end do
    end do

    ! if rigid-body on, update velocity(0 + dt/2 and 0 - dt/2)
    !
    if (constraints%rigid_bond) then

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          imass = 1.0_wip/mass(ix,i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_wip*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., dt, coord_ref, &
                               domain, constraints, coord, vel,             &
                               dynvars%virial_const)

      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          imass = 1.0_wip/mass(ix,i)
          vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                             &
                             - 0.5_wip*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., -dt, coord_ref, &
                               domain, constraints, coord, vel_ref,          &
                               dynvars%virial_const)

      ! calculate velocity(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel_ref(1:3,ix,i) = 0.5_wip*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        end do
      end do

      ! vel <= updated velocities (0) and coord <= constrained
      ! coordinates(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! output dynvars(0)
    !
    dynvars%time = 0.0_wip
    dynvars%step = 0

    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain,        &
                         dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble)

    ! at this point
    !   coord, velocity, and force are at 0

    return

  end subroutine initial_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_pv
  !> @brief        PV1 
  !! @authors      JJ
  !! @param[in]    dt_short    : short time step        
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_pv(dt_short, domain, constraints, dynvars)

    ! formal arguments
    real(wip),               intent(in)    :: dt_short
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix, ncell
    real(wip)                :: h_dt_short
    real(wip)                :: factor

    integer,         pointer :: natom(:)
    real(wip),       pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wip),       pointer :: vel(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: mass(:,:)
    real(dp),        pointer :: viri_const(:,:)


    ! use pointers
    !
    ncell       =  domain%num_cell_local
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    force_long  => domain%force_long
    force_short => domain%force_short
    viri_const  => dynvars%virial_const

    h_dt_short  = dt_short / 2.0_wip

    do i = 1, ncell
      do ix = 1, natom(i)
        factor = h_dt_short / mass(ix,i)
        coord(1:3,ix,i) = coord(1:3,ix,i) + h_dt_short*vel(1:3,ix,i)
      end do
    end do

    ! RATTLE VV1
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., h_dt_short,  &
                               coord_ref, domain, constraints, coord,   &
                               vel, viri_const)
      dynvars%virial_const(1:3,1:3) = 2.0_wip * dynvars%virial_const(1:3,1:3)

    end if

    return

  end subroutine integrate_pv

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dt_short    : short time step        
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv(dt_short, ensemble, domain, dynvars)

    ! formal arguments
    real(wip),                intent(in)    :: dt_short
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: inv_dt, temp0
    real(wip)                :: scale_v, factor
    real(wip)                :: gamma_t
    real(wip)                :: sigma
    real(wip)                :: v1, v2, rsq, grandom(1:3)
    real(wip)                :: half_dt
    integer                  :: j, jx, ncell

    real(wip),       pointer :: mass(:,:)
    real(dp),        pointer :: viri_const(:,:)
    real(wip),       pointer :: random_f(:,:,:)
    real(wip),       pointer :: vel(:,:,:)
    real(wip),       pointer :: force_long(:,:,:), force_short(:,:,:)
    real(wip),       pointer :: coord(:,:,:)
    integer,         pointer :: natom(:)


    inv_dt     =  1.0_wip/dt_short
    temp0      =  ensemble%temperature
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    random_f   => ensemble%random_force
    ncell      =  domain%num_cell_local

    mass       => domain%mass
    natom      => domain%num_atom
    coord      => domain%coord
    vel        => domain%velocity
    force_long => domain%force_long
    force_short=> domain%force_short
    viri_const => dynvars%virial_const

    ! time step
    !
    half_dt  = dt_short / 2.0_wip

    ! random force
    !
    factor   = 2.0_wip*gamma_t*KBOLTZ*temp0/dt_short
    do j = 1, ncell
      do jx = 1, natom(j)

        sigma = sqrt(mass(jx,j) * factor)
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
        random_f(1:3,jx,j) = sigma*grandom(1:3)

      end do
    end do

    ! thermostat
    !
!   scale_v = 1.0_wp - 0.5_wp*gamma_t*dt_short
!   scale_v = 1.0_wp - 1.0_wp*gamma_t*dt_short
!   do j = 1, ncell
!     do jx = 1, natom(j)
!       vel(1:3,jx,j) = vel(1:3,jx,j)*scale_v
!     end do
!   end do

    ! VV
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        factor = dt_short/mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*random_f(1:3,jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
      end do
    end do

    scale_v = 1.0_wip / (1.0_wip + gamma_t*dt_short)
    do j = 1, ncell
      do jx = 1, natom(j)
        vel(1:3,jx,j) = vel(1:3,jx,j)*scale_v
      end do
    end do

  end subroutine langevin_thermostat_vv

end module sp_md_mts_mod
