!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_md_vverlet_mod
!> @brief   Performs molecular dynamics simulations with the velocity Verlet algorithm.
!! @authors Jaewoon Jung (JJ), Tadashi Ando (TA)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_md_vverlet_mod

  use cg_output_mod
  use cg_update_domain_mod
  use cg_assign_velocity_mod
  use cg_dynvars_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
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
  use math_libs_mod
  use messages_mod
  use fileio_grotop_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! subroutines
  public  :: vverlet_dynamics
  private :: initial_vverlet
  private :: integrate_vv1
  private :: integrate_vv2
  private :: nve_vv1
  private :: nve_vv2
  private :: langevin_thermostat_vv1


contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vverlet_dynamics
  !> @brief        Implements a velocity Verlet integrator.
  !! @authors      JJ
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

  subroutine vverlet_dynamics(grotop, molecule, output, domain, enefunc, &
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
    integer                  :: i, nsteps
    integer                  :: iseed
    integer                  :: istart, iend
    logical                  :: npt

    real(dp),        pointer :: virial_cell(:,:), virial(:,:)

    virial_cell    => domain%virial_cellpair
    virial         => dynvars%virial

    nsteps         =  dynamics%nsteps
    istart         =  dynamics%istart_step
    iend           =  dynamics%iend_step
    dt             =  dynamics%timestep/AKMA_PS
    simtim         =  dynamics%initial_time
    iseed          =  dynamics%iseed_init_velocity
    temperature    =  ensemble%temperature
    npt            =  ensemble%use_barostat


    if (abs(dynamics%initial_rmsd) < 0.001_wip)  &
      dynamics%initial_rmsd = real(dynvars%energy%rmsd,wp)
    if (dynamics%target_md) enefunc%rmsd_force = 1.0_wp / real(dt*dt,wp)

    ! Check restart
    !
    if (.not. dynamics%restart .and. dynamics%gen_velocity) then

      call initial_velocity(temperature, iseed, domain)
      call stop_trans_rotation(dynamics, domain)
      call initial_vverlet(npt, output, enefunc, dynamics,    &
                           pairlist, boundary, ensemble,      &
                           domain, dynvars, comm)

    else

      if (istart == 1) then

        call communicate_coor(domain, comm)
        call compute_energy(pairlist, boundary,               &
                            npt, .false.,                     &
                            mod(i,dynamics%eneout_period)==0, &
                            enefunc, domain, dynvars)
        call communicate_force(domain, comm)

      end if

    end if

    ! Main loop
    !
    do i = istart, iend

      dynvars%time = dynamics%timestep * real(i-1,wip)
      dynvars%step = i - 1

      if (dynamics%target_md .or. dynamics%steered_md)       &
        enefunc%rmsd_target = real(dynamics%initial_rmsd,wp) &
                            + real(dynamics%final_rmsd -     &
                                   dynamics%initial_rmsd,wp) &
                             *real(dynvars%step,wip)/real(nsteps,wip)

      if (i > istart) then

        call output_md(output, enefunc, dynamics, boundary, ensemble, &
                       dynvars, domain)

      end if

      ! VV1
      !
      call timer(TimerIntegrator, TimerOn)

      call integrate_vv1(dynamics, i, ensemble, domain, boundary, dynvars)

      call timer(TimerIntegrator, TimerOff)

      ! update cell and pairlist
      !
      if (mod(i-1,dynamics%eneout_period) == 0 .and. i > istart) then
        call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                             dynvars)
        call output_dynvars (output, enefunc, dynvars, ensemble, boundary)
      end if

      ! Remove translational and rotational motion about COM(t + dt)
      !
      if (dynamics%stoptr_period > 0) then
        if (mod(i-1,dynamics%stoptr_period) == 0) then

          call stop_trans_rotation(dynamics, domain)

        end if
      end if

      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then
        if (i > 1 .and. mod(i-1,dynamics%anneal_period) == 0) then
          call simulated_annealing_vverlet(dynamics, enefunc, ensemble)
        end if
      end if

      if (dynamics%nbupdate_period > 0 .and. i > 1) &
        call domain_interaction_update_md(i-1, grotop, dynamics, molecule,  &
                                          domain, enefunc, pairlist,        &
                                          boundary, comm)

      ! calculate potential energy(t + dt), force(t + dt), and virial(t + dt)
      !
      call timer(TimerIntegrator, TimerOn)

      call timer(TimerComm1, TimerOn)

      call communicate_coor(domain, comm)

      call timer(TimerComm1, TimerOff)
      call timer(TimerIntegrator, TimerOff)

      call compute_energy(pairlist, boundary,               &
                          npt, .false.,                     &
                          mod(i,dynamics%eneout_period)==0, &
                          enefunc, domain, dynvars)

      call timer(TimerIntegrator, TimerOn)
      call timer(TimerComm2, TimerOn)

      call communicate_force(domain, comm)

      call timer(TimerComm2, TimerOff)

!     call random_push_stock

      call integrate_vv2(dynamics, i, ensemble, domain, dynvars)

      call timer(TimerIntegrator, TimerOff)

    end do

    ! for final output
    !
    dynvars%time = dynamics%timestep * real(iend,dp)
    dynvars%step = iend

    ! output trajectory(t + dt) and restart data
    !
    call output_md(output, enefunc, dynamics, boundary, ensemble, &
                   dynvars, domain)

    call timer(TimerIntegrator, TimerOn)

    call integrate_vv1(dynamics, iend+1, ensemble, domain, boundary, dynvars)
    call coord_vel_ref(domain, dynvars)

    call timer(TimerIntegrator, TimerOff)
    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                         dynvars)
    call output_dynvars (output, enefunc, dynvars, ensemble, boundary)

    return

  end subroutine vverlet_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_vverlet
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

   subroutine initial_vverlet(npt, output, enefunc, dynamics, pairlist, &
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
    real(dp)                 :: temperature
    real(wip)                :: simtim, dt, friction
    integer                  :: j

    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wp),        pointer :: coord_pbc(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:)
    real(wip),       pointer :: force(:,:), force_long(:,:)
    real(wip),       pointer :: mass(:)
    real(wp),        pointer :: force_omp(:,:,:), force_pbc(:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    integer,         pointer :: natom(:)

    coord          => domain%coord
    coord_ref      => domain%coord_ref
    vel            => domain%velocity
    vel_ref        => domain%velocity_ref
    virial_cell    => domain%virial_cellpair
    virial         => dynvars%virial

    temperature    =  ensemble%temperature
    friction       =  ensemble%gamma_t * AKMA_PS
    dt             =  dynamics%timestep/AKMA_PS
    simtim         =  dynamics%initial_time

    dynvars%time    = simtim
    dynvars%step    = 0


    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do j = 1, domain%num_atom_domain
      coord_ref(j,1) = coord(j,1)
      coord_ref(j,2) = coord(j,2)
      coord_ref(j,3) = coord(j,3)
      vel_ref  (j,1) = vel  (j,1)
      vel_ref  (j,2) = vel  (j,2)
      vel_ref  (j,3) = vel  (j,3)
    end do

    ! calculate energy(0) and forces(0)
    !
    call communicate_coor(domain, comm)

    call compute_energy(pairlist, boundary,    &
                        npt, .false., .true.,  & 
                        enefunc, domain, dynvars)

    call communicate_force(domain, comm)

    ! output dynvars(0)
    !
    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)

    dynamics%restart = .true.

    return

  end subroutine initial_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv1
  !> @brief        Implements VV1 with thermostat/barostat.
  !! @authors      JJ, TA, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv1(dynamics, istep, ensemble, domain, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars

    integer  :: alloc_stat

    if (ensemble%tpcontrol == TpcontrolLangevin .and. &
        istep == 1 .and. .not. allocated(ensemble%random_force)) then
      allocate(ensemble%random_force(MaxAtom_domain,3), stat=alloc_stat)
      if (alloc_stat /= 0) call error_msg('Integrate_VV1: Failed to allocate random_force array.')
    end if

    ! select ensemble, thermostat and barostat
    !
    select case (ensemble%ensemble)

    case (EnsembleNVE)

      call nve_vv1(dynamics, istep, domain, dynvars)

    case (EnsembleNVT) 
    
      select case (ensemble%tpcontrol)
    
      case (TpcontrolBerendsen, TpcontrolBussi, TpcontrolNHC)
    
        call vel_rescaling_thermostat_vv1(dynamics, istep, ensemble, &
                                          domain, dynvars)
    
      case (TpcontrolLangevin)
    
        call langevin_thermostat_vv1(dynamics, istep, ensemble, domain, &
                                     dynvars)
    
      case default
    
        call error_msg('Integrate_VV1: Thermostat is not properly defined.')
    
    end select

    case (EnsembleNPT:EnsembleNPgT) 
    
      select case(ensemble%tpcontrol)
      
      case (TpcontrolBerendsen, TpcontrolBussi, TpcontrolNHC)
    
        call mtk_barostat_vv1(dynamics, istep, ensemble, domain, &
                              boundary, dynvars)
    
      case (TpcontrolLangevin)
    
    !   call langevin_barostat_vv1(dynamics, istep, ensemble, domain, &
    !                              boundary, dynvars)
    
      case default
    
        call error_msg('Integrate_VV1: Available thermostat/barostat '//&
                       'for NPT are only Langevin/Bussi/Berendsen/NHC.')
    
    end select

    case default

      call error_msg('Ensemble is not defined properly')

    end select

    return

  end subroutine integrate_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv2
  !> @brief        Implements VV2 with thermostat/barostat.
  !! @authors      JJ, TA
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv2(dynamics, istep, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_dynvars),         intent(inout) :: dynvars

    select case (ensemble%ensemble)

    case (EnsembleNVE)

      call nve_vv2(dynamics, domain, dynvars)

    case (EnsembleNVT)

      select case(ensemble%tpcontrol)

      case (TpcontrolBussi, TpcontrolBerendsen, TpcontrolNHC, TpcontrolLangevin)

        call nve_vv2(dynamics, domain, dynvars)

      case default

        call error_msg('Thermostat is not properly defined')

      end select

    case (EnsembleNPT:EnsembleNPgT) 

      select case(ensemble%tpcontrol)

      case (TpcontrolLangevin)

!       call langevin_barostat_vv2(dynamics, istep, ensemble, domain,  &
!                                  constraints, boundary, dynvars)

      case (TpcontrolBerendsen, TpcontrolBussi, TpcontrolNHC)

        call mtk_barostat_vv2(dynamics, ensemble, domain, dynvars)

      case default
          
        call error_msg('Integrate_VV1: Available thermostat/barostat for '//&
                       'NPT are only Langevin/Bussi/Berendsen/NHC.')
          
      end select

    case default

      call error_msg('Ensemble is not defined properly')

    end select

    return

  end subroutine integrate_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    coord_vel_ref
  !> @brief        Saves coordinate and velocity at reference values.
  !! @authors      JJ
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine coord_vel_ref(domain, dynvars)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars),         intent(inout) :: dynvars

    ! local variables
    integer                  :: i

    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)

    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    vel_half   => domain%velocity_half

    !$omp parallel do private(i)
    do i = 1, domain%num_atom_domain
      vel_half(i,1) = vel(i,1)
      vel_half(i,2) = vel(i,2)
      vel_half(i,3) = vel(i,3)
      vel(i,1) = vel_ref(i,1)
      vel(i,2) = vel_ref(i,2)
      vel(i,3) = vel_ref(i,3)
      coord(i,1) = coord_ref(i,1)
      coord(i,2) = coord_ref(i,2)
      coord(i,3) = coord_ref(i,3)
    end do
    !$omp end parallel do
    dynvars%nh_velocity(1:5) = dynvars%nh_velocity_ref(1:5)

    return

  end subroutine coord_vel_ref

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv1
  !> @brief        Implements VV1 with NVE.
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv1(dynamics, istep, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, id
    real(wip)                :: dt, half_dt
    real(wip)                :: factor
    real(dp)                 :: kin_temp(3)

    integer,         pointer :: natom(:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: mass(:)
    real(dp),        pointer :: kin(:), kin_full(:), kin_half(:), kin_ref(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin

    natom      => domain%num_atom
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    vel_half   => domain%velocity_half
    force      => domain%force
    ekin_full  => dynvars%ekin_full
    ekin_half  => dynvars%ekin_half
    ekin_ref   => dynvars%ekin_ref
    ekin       => dynvars%ekin
    kin        => dynvars%kin
    kin_full   => dynvars%kin_full
    kin_half   => dynvars%kin_half
    kin_ref    => dynvars%kin_ref

    dt         =  dynamics%timestep/AKMA_PS
    half_dt    =  dt/2.0_wip

    ! reference
    !
    !$omp parallel do private(i)
    do i = 1, domain%num_atom_domain
      vel_ref(i,1)   = vel(i,1)   
      vel_ref(i,2)   = vel(i,2)   
      vel_ref(i,3)   = vel(i,3)   
      coord_ref(i,1) = coord(i,1) 
      coord_ref(i,2) = coord(i,2) 
      coord_ref(i,3) = coord(i,3) 
    end do
    !$omp end parallel do

    if (mod(istep-1, dynamics%eneout_period) == 0) then

      !$omp parallel do private(i, factor)
      do i = 1, domain%num_atom_domain
        factor = half_dt / mass(i)
        vel_half(i,1) = factor*force(i,1)
        vel_half(i,2) = factor*force(i,2)
        vel_half(i,3) = factor*force(i,3)
      end do
      !$omp end parallel do
   
      ! kinetic energy at full time step
      !
      kin_temp(1:3) = 0.0_dp  
      !$omp parallel do private(i) reduction(+:kin_temp)
      do i = 1, domain%num_atom_domain
        kin_temp(1) = kin_temp(1) + mass(i)*vel(i,1)*vel(i,1) 
        kin_temp(2) = kin_temp(2) + mass(i)*vel(i,2)*vel(i,2) 
        kin_temp(3) = kin_temp(3) + mass(i)*vel(i,3)*vel(i,3)
      end do
       !$omp end parallel do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(kin_temp, kin_full, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)
#endif
      ekin_full = 0.5_dp * (kin_full(1)+kin_full(2)+kin_full(3))

      ! kinetic energy at half time step
      !
      kin_temp(1:3) = 0.0_dp  
      !$omp parallel private(i) reduction(+:kin_temp)
      do i = 1, domain%num_atom_domain
        kin_temp(1) = kin_temp(1) + mass(i)*vel_half(i,1)*vel_half(i,1)
        kin_temp(2) = kin_temp(2) + mass(i)*vel_half(i,2)*vel_half(i,2)
        kin_temp(3) = kin_temp(3) + mass(i)*vel_half(i,3)*vel_half(i,3)
      end do
      !$omp end parallel
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(kin_temp, kin_half, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)
#endif
      ekin_half = 0.5_dp * (kin_half(1)+kin_half(2)+kin_half(3))

      !optimal temperature
      !
      ekin = ekin_full + 2.0_dp*ekin_half/3.0_dp
      kin(1:3) = kin_full(1:3) + kin_half(1:3)

    end if

    ! VV1
    !$omp parallel do private(i,factor)
    do i = 1, domain%num_atom_domain
      factor = half_dt / mass(i)
      vel(i,1)   = vel(i,1)   + factor*force(i,1)
      vel(i,2)   = vel(i,2)   + factor*force(i,2)
      vel(i,3)   = vel(i,3)   + factor*force(i,3)
      coord(i,1) = coord(i,1) + dt*vel(i,1)
      coord(i,2) = coord(i,2) + dt*vel(i,2)
      coord(i,3) = coord(i,3) + dt*vel(i,3)
    end do
    !$omp end parallel do

    return

  end subroutine nve_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv2
  !> @brief        Implements VV2 with NVE.
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv2(dynamics, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i
    real(wip)                :: dt, half_dt
    real(dp)                 :: factor

    integer,         pointer :: natom(:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: mass(:)
    real(dp),        pointer :: viri_const(:,:)


    natom      => domain%num_atom
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_half   => domain%velocity_half
    vel_ref    => domain%velocity_full
    force      => domain%force
    viri_const => dynvars%virial_const

    dt         =  dynamics%timestep/AKMA_PS
    half_dt    =  0.5_wip * dt

    ! VV2
    !$omp parallel do private(i, factor)
    do i = 1, domain%num_atom_domain
      factor = half_dt / mass(i)
      vel_half(i,1) = vel(i,1)
      vel_half(i,2) = vel(i,2)
      vel_half(i,3) = vel(i,3)
      vel(i,1) = vel(i,1) + factor*force(i,1)
      vel(i,2) = vel(i,2) + factor*force(i,2)
      vel(i,3) = vel(i,3) + factor*force(i,3)
    end do
    !$omp end parallel do

    return

  end subroutine nve_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_rescaling_thermostat_vv1
  !> @brief        Implements VV1 with Bussi's thermostat (group temperature or no constraint).
  !! @authors      JJ
  !! @param[in]    dt          : time step
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vel_rescaling_thermostat_vv1(dynamics, istep, ensemble, domain,  &
                                          dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    logical                  :: calc_thermostat
    integer                  :: i, natom
    integer                  :: num_degree
    real(wip)                :: temp0, tau_t, dt, half_dt, dt_therm
    real(wip)                :: factor, rr
    real(wip)                :: scale_vel, scale_vel2
    real(dp)                 :: kin_temp(3)

    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: mass(:)
    real(dp),        pointer :: kin(:), kin_full(:), kin_half(:), kin_ref(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin

    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_half    => domain%velocity_half
    force       => domain%force
    ekin_full   => dynvars%ekin_full
    ekin_half   => dynvars%ekin_half
    ekin_ref    => dynvars%ekin_ref
    ekin        => dynvars%ekin
    kin         => dynvars%kin
    kin_full    => dynvars%kin_full
    kin_half    => dynvars%kin_half
    kin_ref     => dynvars%kin_ref

    dt          =  dynamics%timestep / AKMA_PS
    half_dt     =  dt / 2.0_wip
    dt_therm    =  dt * real(dynamics%thermo_period,wip)
    temp0       =  ensemble%temperature
    tau_t       =  ensemble%tau_t/AKMA_PS

    calc_thermostat = mod(istep-1,dynamics%thermo_period) == 0 .and. istep > 1
    num_degree      = domain%num_deg_freedom
    natom           = domain%num_atom_domain

    scale_vel  = 1.0_wip
    scale_vel2 = 1.0_wip

    ! reference coordinates and velocities
    !
    call copy_coord_vel(natom, coord, vel, coord_ref, vel_ref)

    ! thermostat
    !
    if (calc_thermostat) then

      !$omp parallel do private(i)
      do i = 1, natom
        vel_half(i,1) = vel_half(i,1) - vel(i,1)
        vel_half(i,2) = vel_half(i,2) - vel(i,2)
        vel_half(i,3) = vel_half(i,3) - vel(i,3)
      end do
 
      ! kinetic energy
      ! 
      call calc_kinetic(natom, mass, vel_half, kin_half, ekin_half)
      call calc_kinetic(natom, mass, vel_ref , kin_full, ekin_full)

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

    !$omp parallel do private(i)
    do i = 1, natom
      vel(i,1) = vel(i,1)*scale_vel
      vel(i,2) = vel(i,2)*scale_vel
      vel(i,3) = vel(i,3)*scale_vel
    end do
    !$omp end parallel do

    !$omp parallel do private(i, factor)
    do i = 1, natom
      factor = half_dt / mass(i)
      vel(i,1)   = vel(i,1)   + factor*force(i,1)
      vel(i,2)   = vel(i,2)   + factor*force(i,2)
      vel(i,3)   = vel(i,3)   + factor*force(i,3)
      coord(i,1) = coord_ref(i,1) + dt*vel(i,1)
      coord(i,2) = coord_ref(i,2) + dt*vel(i,2)
      coord(i,3) = coord_ref(i,3) + dt*vel(i,3)
    end do
    !$omp end parallel do

    return

  end subroutine vel_rescaling_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv1
  !> @brief        Bussi thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv1(dynamics, istep, ensemble, domain, boundary,  &
                              dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    logical                  :: calc_thermostat, calc_barostat
    real(wip)                :: dt, inv_dt, press0, d_ndegf
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: dt_baro, half_dt_baro, quart_dt
    real(wip)                :: vel_tmp(1:3)
    real(wip)                :: volume, press(1:3)
    real(wip)                :: pressxy, pressxyz
    real(wip)                :: factor
    real(wip)                :: bmoment_ref(3), scale_b(1:3), vel_scale_2(1:3)
    real(wip)                :: force_change(1:3)
    real(wip)                :: vel_scale(1:3), force_scale_2(1:3)
    real(wip)                :: tau_t, tau_p
    real(wip)                :: gr, random_gr
    real(wip)                :: half_dt, size_scale(1:3), scale_vel
    real(dp)                 :: virial_constraint(3,3), virial_sum(3)
    real(dp)                 :: kin_full_ref(3), ekin_full_ref
    integer                  :: i, j, k, l, jx, iter, maxiter
    integer                  :: i_ndegf
    integer                  :: natom, num_degree

    real(wip),       pointer :: pmass, ekin0, temp0, degree
    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: temporary(:,:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)
    real(wip),       pointer :: force(:,:)
    real(dp),        pointer :: virial(:,:)
    real(dp),        pointer :: kin_full(:), kin_half(:), kin_ref(:), kin(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin
    real(wip),       pointer :: bmoment(:)

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_wip/dt
    temp0      => ensemble%temperature
    ekin0      => ensemble%kinetic
    pmass      => ensemble%pmass
    press0     =  ensemble%pressure * ATMOS_P
    tau_t      =  ensemble%tau_t / AKMA_PS
    tau_p      =  ensemble%tau_p / AKMA_PS
    i_ndegf    =  domain%num_deg_freedom
    natom      =  domain%num_atom_domain
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    vel_half   => domain%velocity_half
    force      => domain%force
    temporary  => domain%velocity_full
    virial     => dynvars%virial
    kin_full   => dynvars%kin_full
    kin_half   => dynvars%kin_half
    kin_ref    => dynvars%kin_ref
    kin        => dynvars%kin
    ekin_full  => dynvars%ekin_full
    ekin_half  => dynvars%ekin_half
    ekin_ref   => dynvars%ekin_ref
    ekin       => dynvars%ekin
    bmoment    => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt  = dt / 2.0_wip
    quart_dt = half_dt / 2.0_wip
    dt_therm = dt * real(dynamics%thermo_period,wip)
    half_dt_therm = dt_therm / 2.0_wip
    dt_baro  = dt * real(dynamics%baro_period,wip)
    half_dt_baro = dt_baro / 2.0_wip

    calc_thermostat = mod(istep-1, dynamics%thermo_period) == 0 .and. istep > 1
    calc_barostat = mod(istep-1, dynamics%baro_period) == 0 .and. istep > 1

    ! pmass and ekin0
    !
    if (istep == 1) then
      num_degree = i_ndegf + 3
      if (ensemble%isotropy == IsotropyXY_Fixed) num_degree = i_ndegf + 1
      ekin0 = 0.5_wip*KBOLTZ*temp0 * real(num_degree,dp)
      pmass  = degree*KBOLTZ*temp0 * tau_p*tau_p
    end if

    ! initial values
    !
    random_gr = real(random_get_gauss(),wip)

    ! reference
    !
    bmoment_ref(1:3) = bmoment(1:3)
    call copy_coord_vel(natom, coord, vel, coord_ref, vel_ref)

    ! temperature at full time step
    !
    if (calc_thermostat) then

      !$omp parallel do private(j)
      do j = 1, natom
        vel_half(j,1) = vel_half(j,1) - vel(j,1)
        vel_half(j,2) = vel_half(j,2) - vel(j,2)
        vel_half(j,3) = vel_half(j,3) - vel(j,3)
      end do
      !$omp end parallel do

      call calc_kinetic(natom, mass, vel_half, kin_half, ekin_half)
      call calc_kinetic(natom, mass, vel_ref , kin_full, ekin_full)

      ekin = ekin_full + 2.0_dp*ekin_half/3.0_dp
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))

      ! calculate scaling factor
      !
      random_gr = real(random_get_gauss(),wip)
      if (ensemble%tpcontrol == TpcontrolBussi) then
        call vel_scale_bussi(num_degree, half_dt_therm, tau_t, temp0, ekin, &
                             random_gr, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolBerendsen) then
        call vel_scale_berendsen(num_degree, half_dt_therm, tau_t, temp0,   &
                                 ekin, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolNHC) then
        call vel_scale_nhc(num_degree, half_dt_therm, tau_t, temp0, ekin,   &
                           ensemble, dynvars, scale_vel)
      end if

      !$omp parallel do private(j)
      do j = 1, natom
        vel(j,1) = vel(j,1)*scale_vel
        vel(j,2) = vel(j,2)*scale_vel
        vel(j,3) = vel(j,3)*scale_vel
      end do
      bmoment(1:3) = bmoment(1:3) * scale_vel

      kin_full(1:3) = kin_full(1:3)*scale_vel*scale_vel
      kin(1:3) = kin_full(1:3) + kin_half(1:3)
      ekin = 0.5_dp*(kin_full(1)+kin_full(2)+kin_full(3)) &
           + 2.0_dp*ekin_half/3.0_dp

    end if

    if (calc_barostat) then

      virial_sum(1) = virial(1,1)
      virial_sum(2) = virial(2,2)
      virial_sum(3) = virial(3,3)
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, virial_sum, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)
#endif
      press(1:3) = (kin(1:3) + virial_sum(1:3))/volume
      pressxyz = (press(1)+press(2)+press(3))/3.0_dp
      pressxy  = (press(1)+press(2))/2.0_dp

      ! update barostat
      !
      call update_barostat(ensemble, boundary, press, pressxyz, pressxy,  &
                           press0, volume, d_ndegf, pmass, dt_baro, ekin, &
                           bmoment)
    end if

    if (calc_thermostat) then

      ! calculate scaling factor
      !
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))
      random_gr = real(random_get_gauss(),wip)
      if (ensemble%tpcontrol == TpcontrolBussi) then
        call vel_scale_bussi(num_degree, half_dt_therm, tau_t, temp0, ekin, &
                             random_gr, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolBerendsen) then
        call vel_scale_berendsen(num_degree, half_dt_therm, tau_t, temp0,   &
                                 ekin, scale_vel)
      else if (ensemble%tpcontrol == TpcontrolNHC) then
        call vel_scale_nhc(num_degree, half_dt_therm, tau_t, temp0, ekin,   &
                           ensemble, dynvars, scale_vel)
      end if

      !$omp parallel do private(j)
      do j = 1, natom
        vel(j,1) = vel(j,1)*scale_vel
        vel(j,2) = vel(j,2)*scale_vel
        vel(j,3) = vel(j,3)*scale_vel
      end do
      !$omp end parallel do

      bmoment(1:3) = bmoment(1:3) * scale_vel

      kin_full(1:3) = kin_full(1:3)*scale_vel*scale_vel
      ekin = 0.5_dp * (kin_full(1)+kin_full(2)+kin_full(3)) &
           + 2.0_dp*ekin_half/3.0_dp
      ekin = ekin + 0.5_dp*pmass*dot_product(bmoment(1:3),bmoment(1:3))

    end if

    ! VV1
    !
    size_scale(1:3)  = exp(bmoment(1:3)*dt)
    gr = bmoment(1)+bmoment(2)+bmoment(3)
    scale_b(1:3) = bmoment(1:3) + gr/degree
    vel_scale(1:3) = exp(-scale_b(1:3)*half_dt)
    vel_scale_2(1:3) = exp(bmoment(1:3)*half_dt)
    vel_scale_2(1:3) = vel_scale_2(1:3)*powersinh(bmoment(1:3)*half_dt)
    force_scale_2(1:3) = exp(-scale_b(1:3)*quart_dt)
    force_scale_2(1:3) = force_scale_2(1:3)*powersinh(scale_b(1:3)*quart_dt)

    !$omp parallel do private(j, factor)
    do j = 1, natom
      factor = half_dt / mass(j)
      vel(j,1:3)   = vel_scale(1:3)*vel(j,1:3) &
                   + factor*force_scale_2(1:3)*force(j,1:3)
      coord(j,1:3) = size_scale(1:3)*coord_ref(j,1:3) &
                   + dt*vel_scale_2(1:3)*vel(j,1:3)
    end do
    !$omp end parallel do

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
    natom = domain%num_atom_domain + domain%num_atom_boundary
    !$omp parallel do private(j, natom)
    do j = 1, natom
      domain%trans_vec(j,1:3) = domain%trans_vec(j,1:3) * scale_b(1:3)
    end do
    !$omp end parallel do

    domain%system_size(1) = boundary%box_size_x
    domain%system_size(2) = boundary%box_size_y
    domain%system_size(3) = boundary%box_size_z

    return

  end subroutine mtk_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv2
  !> @brief        Bussi thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv2(dynamics, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, half_dt, quart_dt
    real(dp)                 :: cm(8), viri_c(3)
    real(wip)                :: force_change(3)
    real(wip)                :: factor, alpha
    real(wip)                :: scale_b(3)
    real(wip)                :: vel_scale(3), force_scale_2(3)
    real(wip)                :: vel_change(3)
    integer                  :: j, jx, k, natom

    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: vel(:,:), vel_ref(:,:), vel_half(:,:)
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(dp),        pointer :: virial(:,:)
    real(wip),       pointer :: bmoment(:), degree

    dt         =  dynamics%timestep/AKMA_PS
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    vel_half   => domain%velocity_half
    force      => domain%force
    virial     => dynvars%virial
    bmoment    => dynvars%barostat_momentum

    natom      =  domain%num_atom_domain

    ! time step
    !
    half_dt = 0.5_wip * dt
    quart_dt = 0.25_wip * dt

    !$omp parallel do private(j)
    do j = 1, natom
      vel_half(j,1) = vel(j,1)
      vel_half(j,2) = vel(j,2)
      vel_half(j,3) = vel(j,3)
    end do

    alpha = bmoment(1)+bmoment(2)+bmoment(3)
    scale_b(1:3) = bmoment(1:3) + alpha/degree
    vel_scale(1:3) = exp(-scale_b(1:3)*half_dt)
    force_scale_2(1:3) = exp(-scale_b(1:3)*quart_dt)
    force_scale_2(1:3) = force_scale_2(1:3)*powersinh(scale_b(1:3)*quart_dt)

    ! VV2
    !
    !$omp parallel do private(j)
    do j = 1, natom
      factor = half_dt / mass(j)
      vel(j,1:3) = vel_scale(1:3)*vel(j,1:3) &
                 + factor*force_scale_2(1:3)*force(j,1:3)
    end do

    return

  end subroutine mtk_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_scale_bussi
  !> @brief        Calculates the scaling factor with Bussi's thermostat.
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
  !> @brief        Calculates the scaling factor with Berendsen's thermostat.
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

    real(wip)                :: tempf, factor

    factor = exp(-dt/tau_t)
    tempf = 2.0_wip * real(ekin,wip)/(real(degree,wip)*KBOLTZ)
    scale_vel = sqrt(1.0_wip + (dt/tau_t)*(temp0/tempf-1.0_wip))

    return

  end subroutine vel_scale_berendsen

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vel_scale_nhc
  !> @brief        Calculates the scaling factor with NHC thermostat.
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
    real(wip)                :: scale_kin
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
  !> @brief        Copies coordinate and velocity.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_coord_vel(natom, coord_in, vel_in, coord_out, vel_out)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wip),               intent(in)    :: coord_in(:,:)
    real(wip),               intent(in)    :: vel_in(:,:)
    real(wip),               intent(inout) :: coord_out(:,:)
    real(wip),               intent(inout) :: vel_out(:,:)

    integer                  :: i

    !$omp parallel do
    do i = 1, natom
      coord_out(i,1) = coord_in(i,1)
      coord_out(i,2) = coord_in(i,2)
      coord_out(i,3) = coord_in(i,3)
      vel_out  (i,1) = vel_in  (i,1)
      vel_out  (i,2) = vel_in  (i,2)
      vel_out  (i,3) = vel_in  (i,3)
    end do
    !$omp end parallel do

    return

  end subroutine copy_coord_vel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_barostat
  !> @brief        Updates barostat parameter bmoment.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat(ensemble, boundary, press, pressxyz, pressxy, &
                             press0, volume, d_ndegf, pmass, dt, ekin, bmoment)

    ! formal arguments
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: press(:)
    real(wip),               intent(in)    :: pressxyz
    real(wip),               intent(in)    :: pressxy
    real(wip),               intent(in)    :: press0
    real(wip),               intent(in)    :: volume
    real(wip),               intent(in)    :: d_ndegf
    real(wip),               intent(in)    :: pmass
    real(wip),               intent(in)    :: dt
    real(dp),                intent(in)    :: ekin
    real(wip),               intent(inout) :: bmoment(:)

    real(wip)                              :: ekin_real, gamma0, pressxy0

    gamma0    =  ensemble%gamma*ATMOS_P*100.0_wip/1.01325_wip
    ekin_real = real(ekin,wip)

    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + dt*(volume*(pressxyz - press0)  &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + dt*(volume*(pressxy - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf)/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(3) + dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf)/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + dt*(volume*(pressxy - pressxy0) &
                                + 2.0_wip*ekin_real/d_ndegf)/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(3) + dt*(volume*(press(3) - press0)   &
                                + 2.0_wip*ekin_real/d_ndegf)/pmass
      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1:3) = bmoment(1:3) + dt*(volume*(press(1:3) - press0)    &
                                    + 2.0_wip*ekin_real/d_ndegf)/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1:2) = bmoment(1:2) + dt*(volume*(press(1:2) - pressxy0)  &
                                    + 2.0_wip*ekin_real/d_ndegf)/pmass
        bmoment(3)   = bmoment(3) + dt*(volume*(press(3) - press0)   &
                                  + 2.0_wip*ekin_real/d_ndegf)/pmass
      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_wip
      bmoment(2) = 0.0_wip
      bmoment(3) = bmoment(3) + dt*(volume*(press(3) - press0) &
                              + 2.0_wip*ekin_real/d_ndegf)/pmass

    end if

    return

  end subroutine update_barostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bcast_boxsize
  !> @brief        Broadcasts box size.
  !! @authors      JJ
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
  !  Subroutine    langevin_thermostat_vv1
  !> @brief        Controls temperature using a Langevin thermostat.
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv1(dynamics, istep, ensemble, domain,  &
                                     dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wip)                :: dt, half_dt, inv_dt
    real(wip)                :: dt_therm, half_dt_therm
    real(wip)                :: temp0, gamma_t, scale_v
    real(wip)                :: factor, sigma
    real(wip)                :: grandom(1:3)
    real(wip)                :: kBT
    real(dp)                 :: kin_temp(3)
    integer                  :: i, natom, id

    real(wip),       pointer :: vel(:,:), vel_ref(:,:)
    real(wip),       pointer :: coord(:,:), coord_ref(:,:)
    real(wip),       pointer :: force(:,:)
    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: random_f(:,:)
    real(dp),        pointer :: kin(:), kin_full(:), kin_half(:), kin_ref(:)
    real(dp),        pointer :: ekin_full, ekin_half, ekin_ref, ekin

    dt            =  dynamics%timestep/AKMA_PS
    half_dt       =  dt / 2.0_wip
    dt_therm      =  dt * real(dynamics%thermo_period, dp)
    half_dt_therm =  dt_therm / 2.0_wip
    inv_dt        =  1.0_wip/dt
    temp0         =  ensemble%temperature
    gamma_t       =  ensemble%gamma_t *AKMA_PS
    random_f      => ensemble%random_force
    natom         =  domain%num_atom_domain

    mass          => domain%mass
    coord         => domain%coord
    coord_ref     => domain%coord_ref
    vel           => domain%velocity
    vel_ref       => domain%velocity_ref
    force         => domain%force
    ekin_full     => dynvars%ekin_full
    ekin_half     => dynvars%ekin_half
    ekin_ref      => dynvars%ekin_ref
    ekin          => dynvars%ekin
    kin           => dynvars%kin
    kin_full      => dynvars%kin_full
    kin_half      => dynvars%kin_half
    kin_ref       => dynvars%kin_ref

    ! setup variables
    !
    kBT      = KBOLTZ * temp0

    ! reference
    !
    !$omp parallel do private(i)
    do i = 1, natom
      vel_ref(i,1) = vel(i,1)
      vel_ref(i,2) = vel(i,2)
      vel_ref(i,3) = vel(i,3)
      coord_ref(i,1) = coord(i,1)
      coord_ref(i,2) = coord(i,2)
      coord_ref(i,3) = coord(i,3)
    end do
    !$omp end parallel do

    ! random force
    !
    if (mod(istep-1, dynamics%thermo_period) == 0) then
      scale_v = exp(- gamma_t*dt_therm)
      factor  = 1.0_wip - scale_v*scale_v
      factor  = factor*KBOLTZ*temp0

      !$omp parallel private(id, i, sigma, grandom)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
      do i = id+1, domain%num_atom_domain, nthread
        if (abs(mass(i)) < EPS) cycle
        sigma = sqrt(factor/mass(i))
        grandom(1) = random_get_gauss_omp(id+1)
        grandom(2) = random_get_gauss_omp(id+1)
        grandom(3) = random_get_gauss_omp(id+1)
        random_f(i,1) = sigma*grandom(1)
        random_f(i,2) = sigma*grandom(2)
        random_f(i,3) = sigma*grandom(3)
      end do
      !$omp end parallel
      if (mod(istep-1, dynamics%eneout_period) == 0 ) &
      call calc_kinetic(natom, mass, vel_ref, kin_full, ekin_full)
    end if

    ! VV1 (B, A)
    !
    !$omp parallel do private(i, factor)
    do i = 1, natom
      factor = half_dt / mass(i)
      vel(i,1) = vel_ref(i,1) + factor*force(i,1) 
      vel(i,2) = vel_ref(i,2) + factor*force(i,2) 
      vel(i,3) = vel_ref(i,3) + factor*force(i,3) 
      coord(i,1) = coord_ref(i,1) + vel(i,1)*half_dt
      coord(i,2) = coord_ref(i,2) + vel(i,2)*half_dt
      coord(i,3) = coord_ref(i,3) + vel(i,3)*half_dt
    end do
    !$omp end parallel do

    if (mod(istep-1, dynamics%thermo_period) == 0) then

      if (mod(istep-1, dynamics%eneout_period) == 0 ) &
      call calc_kinetic(natom, mass, vel, kin_half, ekin_half)

      ! Thermostat
      !
      !$omp parallel do private(i)
      do i = 1, domain%num_atom_domain
        vel(i,1) = vel(i,1) * scale_v
        vel(i,2) = vel(i,2) * scale_v
        vel(i,3) = vel(i,3) * scale_v
        vel(i,1) = vel(i,1) + random_f(i,1)
        vel(i,2) = vel(i,2) + random_f(i,2)
        vel(i,3) = vel(i,3) + random_f(i,3)
      end do
      !$omp end parallel do
      if (mod(istep-1, dynamics%eneout_period) == 0 ) &
      call calc_kinetic(natom, mass, vel, kin_half, ekin_half)
      ekin = ekin_half
      kin(1:3) = kin_half(1:3)

    end if

    ! VV1 (A)
    !
    !$omp parallel do private(i, factor)
    do i = 1, natom
      coord(i,1) = coord(i,1) + vel(i,1)*half_dt
      coord(i,2) = coord(i,2) + vel(i,2)*half_dt
      coord(i,3) = coord(i,3) + vel(i,3)*half_dt
    end do
    !$omp end parallel do

    return

  end subroutine langevin_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_vverlet 
  !> @brief        Changes target temperature linearly.
  !! @authors      TM
  !! @param[in]    dynamics: dynamics information
  !! @param[out]   ensemble: ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine simulated_annealing_vverlet(dynamics, enefunc, ensemble)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_enefunc ),        intent(inout) :: enefunc 
    type(s_ensemble),        intent(inout) :: ensemble

    ! local variable
    real(dp)                 :: old_temperature


    if (.not. dynamics%annealing) return

    old_temperature      = ensemble%temperature
    ensemble%temperature = ensemble%temperature + dynamics%dtemperature

    if (enefunc%cg_ele_calc) then
      enefunc%cg_ele_sol_T = ensemble%temperature
    end if

    if (main_rank) then
      write(MsgOut,'(A,F10.3,A,F10.3)')                                     &
            'Simulated_Annealing_Vverlet: Anneals temperature from', &
            old_temperature, ' to ', ensemble%temperature

      if (ensemble%temperature < 0.0_dp) &
        call error_msg( &
        'Simulated_Annealing_Leapfrog> Error: Temperature is less than 0 K')

    end if

    return

  end subroutine simulated_annealing_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_kinetic    
  !> @brief        Performs kinetic energy calculation.
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_kinetic(natom, mass, vel, kin, ekin)

    ! formal arguments
    integer,                 intent(in   ) :: natom
    real(wip),               intent(inout) :: mass(:)
    real(wip),               intent(inout) :: vel (:,:)
    real(dp),                intent(inout) :: kin(:)
    real(dp),                intent(inout) :: ekin

    integer                  :: i, ix

    kin(1:3) = 0.0_dp
    !$omp parallel do private(i) reduction(+:kin)
    do i = 1, natom
      kin(1) = kin(1) + mass(i)*vel(i,1)*vel(i,1)
      kin(2) = kin(2) + mass(i)*vel(i,2)*vel(i,2)
      kin(3) = kin(3) + mass(i)*vel(i,3)*vel(i,3)
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

end module cg_md_vverlet_mod
