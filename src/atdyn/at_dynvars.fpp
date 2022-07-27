!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_dynvars_mod
!> @brief   define dynamic variables
!! @authors Takaharu Mori (TM), Yuji Sugita (YS), Chigusa Kobayashi (CK),
!!          Kiyoshi Yagi (KY), Cheng Tan (CT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_dynvars_mod

  use at_output_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use at_restraints_str_mod
  use math_libs_mod
  use fileio_rst_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use random_mod

  implicit none
  private

  ! variables
  integer, public, save :: DynvarsOut = 6 ! this is necessary for REMD output.

  logical,         save :: etitle = .true.
  logical,         save :: vervose = .true.

  ! subroutines
  public  :: setup_dynvars
  public  :: output_dynvars
  public  :: compute_dynvars
  private :: output_dynvars_genesis
  private :: output_dynvars_charmm
  private :: output_dynvars_namd
  private :: output_dynvars_gromacs
  private :: compute_pressure
  public  :: generate_velocity
  public  :: stop_trans_rotation

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynvars
  !> @brief        copy molecule information into dynvars
  !! @authors      TM, YS
  !! @param[in]    molecule  : molecular information
  !! @param[in]    rst       : restart information
  !! @param[inout] dynvars   : dynamic variables
  !! @param[inout] dynamics  : dynamics information    (optional)
  !! @param[in]    tpcontrol : thermostat and barostat (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynvars(molecule, rst, dynvars, dynamics, tpcontrol)

    ! formal arguments
    type(s_molecule),                   intent(in)    :: molecule
    type(s_rst),                        intent(in)    :: rst
    type(s_dynvars),                    intent(inout) :: dynvars
    type(s_dynamics), target, optional, intent(inout) :: dynamics
    integer,                  optional, intent(in)    :: tpcontrol

    ! local variables
    integer                  :: natom


    natom = molecule%num_atoms

    !
    ! error check for tpcontrol & integrator
    !
    if (present(dynamics)) then
      if (dynamics%integrator == IntegratorLEAP .and. &
          tpcontrol           == TpcontrolBUSSI) then
        call error_msg('Setup_Dynvars> BUSSI is not allowed in Leap-frog')
      end if
    endif

    ! initialize dynvars
    !
    call init_dynvars(dynvars)

    ! allocate dynvars
    !
    call alloc_dynvars(dynvars, DynvarsGeneral, natom)

    ! setup coordinates
    !
    if (rst%rstfile_type == RstfileTypeUndef) then
      dynvars%coord    (1:3,1:natom) = molecule%atom_coord(1:3,1:natom)
      dynvars%coord_ref(1:3,1:natom) = dynvars%coord      (1:3,1:natom)
    else
      dynvars%coord    (1:3,1:natom) = rst%coord          (1:3,1:natom)
      dynvars%coord_ref(1:3,1:natom) = dynvars%coord      (1:3,1:natom)
    end if

    ! setup velocities
    !
    dynvars%velocity(1:3,1:natom)  = 0.0_wp
    dynvars%thermostat_momentum    = 0.0_wp
    dynvars%barostat_momentum(1:3) = 0.0_wp

    if (rst%rstfile_type == RstfileTypeUndef .or. &
        rst%rstfile_type == RstfileTypeMin ) then
      ! do nothing (or generate initial velocities at 0 step in MD)
    else
      if (present(dynamics)) then
        dynvars%velocity(1:3,1:natom)  = rst%velocity(1:3,1:natom)
        dynvars%thermostat_momentum    = rst%thermostat_momentum
        dynvars%barostat_momentum(1:3) = rst%barostat_momentum(1:3)
      end if
    end if
    dynvars%velocity_ref(1:3,1:natom)  = dynvars%velocity(1:3,1:natom)

    ! setup random system
    !
    if (present(dynamics)) then

      call random_init(dynamics%iseed)

      if (allocated(rst%random) .and. dynamics%random_restart) then
        call random_stock_mpi_frombyte(mpi_comm_country, rst%random)
        call random_pull_stock
      end if

    end if

    ! setup forces
    !
    dynvars%force(1:3,1:natom) = 0.0_wp

    ! setup random force
    !
    if (present(dynamics)) then
      if ((dynamics%integrator == IntegratorVVER .or.      &
           dynamics%integrator == IntegratorVVER_CG) .and. &
           tpcontrol           == TpcontrolLangevin) then
        call alloc_dynvars(dynvars, DynvarsLangevin,natom, nproc_country)
      end if
    end if

    ! change step number
    !
    if (present(dynamics)) then
      if (dynamics%restart .and. dynamics%integrator == IntegratorLEAP) then
        dynvars%step = 1
      else
        dynvars%step = 0
      end if
    end if

    return

  end subroutine setup_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars
  !> @brief        output dynamical variables
  !! @authors      YS
  !! @param[in]    output   : output information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamical variables information
  !! @param[in]    ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars(output, enefunc, dynvars, ensemble)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_ensemble), optional, intent(in) :: ensemble


    ! If REMD, output is done on replica_main_rank
    ! If MD, output is done on main_rank
    !
    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else if (output%vib) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if


    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_dynvars_genesis(output, enefunc, dynvars, ensemble)

    case (OutputStyleCHARMM)

      call output_dynvars_charmm (enefunc, dynvars)

    case (OutputStyleNAMD)

      call output_dynvars_namd   (enefunc, dynvars)

    case (OutputStyleGROMACS)

      call output_dynvars_gromacs(enefunc, dynvars)

    end select

    return

  end subroutine output_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_dynvars
  !> @brief        calculate other dynamical variables
  !! @authors      TM, CK
  !! @param[in]    molecule : molecular information
  !! @param[in]    enefunc  : potential energy function
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamics variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_dynvars(molecule, enefunc, dynamics, boundary, ensemble, &
                             dynvars)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_boundary),         intent(in)    :: boundary
    type(s_ensemble),         intent(in)    :: ensemble
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(wp)                  :: ekin, ekin_old, ekin_new
    real(wp)                  :: vel_oldx, vel_oldy, vel_oldz
    real(wp)                  :: rmsg, viri_ke
    real(wp)                  :: kex, key, kez
    real(wp)                  :: viri_ke_ext
    real(wp)                  :: gamma_t, dt, scale_v
    real(wp)                  :: totresene
    real(wp)                  :: virial_sum(1:3,1:3)
    integer                   :: j, natom

    real(wp),         pointer :: mass(:), coord_ref(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:), force(:,:)
    real(wp),         pointer :: vel_full(:,:)
    real(wp),         pointer :: virial(:,:)
    real(wp),         pointer :: virial_ext(:,:)
    real(wp),         pointer :: viri_const(:,:)


    mass        => molecule%mass
    coord_ref   => dynvars%coord_ref
    vel         => dynvars%velocity
    vel_ref     => dynvars%velocity_ref
    vel_full    => dynvars%temporary
    force       => dynvars%force
    virial      => dynvars%virial
    virial_ext  => dynvars%virial_extern
    natom       =  molecule%num_atoms
    gamma_t     =  ensemble%gamma_t*AKMA_PS
    dt          =  dynamics%timestep/AKMA_PS
    viri_const  => dynvars%virial_const


    if ((dynamics%integrator /= IntegratorBDEM) .and.  &
        (dynamics%integrator /= IntegratorBD2N) .and.  &
        (dynamics%integrator /= IntegratorSDMP)) then

      ! Compute velocities(t + dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      if (dynamics%integrator == IntegratorLEAP) then

        if (ensemble%tpcontrol == TpcontrolLangevin) then
          scale_v = sqrt(1.0_wp+0.5_wp*gamma_t*dt)
        else
          scale_v = 1.0_wp
        end if

        vel_full(1:3,1:natom) =                                          &
              0.5_wp*(vel(1:3,1:natom) + vel_ref(1:3,1:natom))*scale_v
      else
        vel_full(1:3,1:natom) = vel(1:3,1:natom)
      end if


      ! kinetic energy(t+dt) and temperature(t+dt)
      ! 
      kex = 0.0_wp
      key = 0.0_wp
      kez = 0.0_wp

      do j = 1, natom
        kex = kex + mass(j)*vel_full(1,j)*vel_full(1,j)
        key = key + mass(j)*vel_full(2,j)*vel_full(2,j)
        kez = kez + mass(j)*vel_full(3,j)*vel_full(3,j)
      end do

      ekin  = 0.5_wp * (kex + key + kez)
      dynvars%total_kene  = ekin
      dynvars%temperature = 2.0_wp*ekin/(molecule%num_deg_freedom*KBOLTZ)

    else

      kex = 0.0_wp
      key = 0.0_wp
      kez = 0.0_wp

      ekin  = 0.0_wp
      dynvars%total_kene  = ekin
      dynvars%temperature = ensemble%temperature

    end if


    totresene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      totresene = sum(dynvars%energy%restraint(1:enefunc%num_restraintfuncs))
    endif
    ! potential energy(t+dt) and total energy(t+dt)
    !
    dynvars%total_pene =  &
               dynvars%energy%bond               + &
               dynvars%energy%angle              + &
               dynvars%energy%urey_bradley       + &
               dynvars%energy%dihedral           + &
               dynvars%energy%improper           + &
               dynvars%energy%cmap               + &
               dynvars%energy%electrostatic      + &
               dynvars%energy%van_der_waals      + &
               dynvars%energy%contact            + &
               dynvars%energy%noncontact         + &
               dynvars%energy%solvation          + &
               dynvars%energy%base_stacking      + &
               dynvars%energy%base_pairing       + &
               dynvars%energy%cg_DNA_exv         + &
               dynvars%energy%cg_IDR_HPS         + &
               dynvars%energy%cg_IDR_KH          + &
               dynvars%energy%cg_KH_inter_pro    + &
               dynvars%energy%cg_exv             + &
               dynvars%energy%PWMcos             + &
               dynvars%energy%PWMcosns           + &
               totresene

    if (enefunc%dispersion_corr .ne. Disp_Corr_NONE) then

      dynvars%total_pene = dynvars%total_pene + &
                           dynvars%energy%disp_corr_energy

    end if

    if (enefunc%qmmm%do_qmmm) then

      dynvars%total_pene = dynvars%total_pene + &
                           dynvars%energy%qm_ene * CONV_UNIT_ENE

    end if

    dynvars%total_energy = &
               dynvars%total_pene + dynvars%total_kene

    ! GaMD
    !
    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_dih) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%dihedral_gamd
      else if (enefunc%gamd%boost_pot) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%total_gamd
      else if (enefunc%gamd%boost_dual) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%total_gamd + &
                             dynvars%energy%dihedral_gamd
      end if
    end if

    ! volume(t+dt) and pressure(t+dt)
    !   note that boundary%box_size     is at t + 2dt
    !             boundary%box_size_ref is at t +  dt
    !
    if (boundary%type == BoundaryTypePBC) then

      dynvars%box_size_x = boundary%box_size_x_ref
      dynvars%box_size_y = boundary%box_size_y_ref
      dynvars%box_size_z = boundary%box_size_z_ref

      dynvars%volume     = boundary%box_size_x_ref  &
                         * boundary%box_size_y_ref  &
                         * boundary%box_size_z_ref
    
    else if (boundary%type == BoundaryTypeNOBC) then

      dynvars%box_size_x = 0.0_wp
      dynvars%box_size_y = 0.0_wp
      dynvars%box_size_z = 0.0_wp

      dynvars%volume     = 0.0_wp
    
    end if


    if (boundary%type /= BoundaryTypeNOBC) then

      if (dynamics%integrator == IntegratorLEAP) then
        call compute_pressure(molecule%num_deg_freedom, dynvars%volume, &
                            dynvars%temperature, dynvars%virial,        &
                            dynvars%internal_pressure)
      else
        virial_sum(1:3,1:3) = dynvars%virial(1:3,1:3) + viri_const(1:3,1:3)
        call compute_pressure(molecule%num_deg_freedom, dynvars%volume, &
                            dynvars%temperature, virial_sum,            &
                            dynvars%internal_pressure)
      end if

      call compute_pressure(molecule%num_deg_freedom, dynvars%volume,   &
                            dynvars%temperature, dynvars%virial_extern, &
                            dynvars%external_pressure)
    else

      dynvars%internal_pressure = 0.0_wp
      dynvars%external_pressure = 0.0_wp

      virial_sum(1:3,1:3) = 0.0_wp
    end if

    if (dynvars%volume /= 0.0_wp) then
      if (dynamics%integrator == IntegratorLEAP) then
        dynvars%pressure_xx = P_ATMOS*(kex + virial(1,1))/dynvars%volume
        dynvars%pressure_yy = P_ATMOS*(key + virial(2,2))/dynvars%volume
        dynvars%pressure_zz = P_ATMOS*(kez + virial(3,3))/dynvars%volume

      else
        dynvars%pressure_xx = P_ATMOS*(kex + virial_sum(1,1))/dynvars%volume
        dynvars%pressure_yy = P_ATMOS*(key + virial_sum(2,2))/dynvars%volume
        dynvars%pressure_zz = P_ATMOS*(kez + virial_sum(3,3))/dynvars%volume

      end if
    else
      dynvars%pressure_xx = 0.0_wp
      dynvars%pressure_yy = 0.0_wp
      dynvars%pressure_zz = 0.0_wp
    end if

    ! other dynvars(t+dt)
    !
    ekin_old = 0.0_wp
    ekin_new = 0.0_wp
    rmsg     = 0.0_wp

    if (dynamics%integrator == IntegratorLEAP) then

      do j = 1, natom
        vel_oldx = 2.0_wp * vel_full(1,j) - vel(1,j)
        vel_oldy = 2.0_wp * vel_full(2,j) - vel(2,j)
        vel_oldz = 2.0_wp * vel_full(3,j) - vel(3,j)
        ekin_old = ekin_old + (vel_oldx**2 + vel_oldy**2 + vel_oldz**2)*mass(j)
        ekin_new = ekin_new + (vel(1,j)**2 + vel(2,j)**2 + vel(3,j)**2)*mass(j)
        rmsg     = rmsg     + force(1,j)**2 + force(2,j)**2 + force(3,j)**2
      end do

      ekin_old = 0.5_wp * ekin_old
      ekin_new = 0.5_wp * ekin_new

    else if (dynamics%integrator == IntegratorVVER .or.       &
             dynamics%integrator == IntegratorVVER_CG) then

      do j = 1, natom
        rmsg = rmsg + force(1,j)**2 + force(2,j)**2 + force(3,j)**2
      end do

    end if
    
    dynvars%hfc_kene     = 0.5_wp*(ekin_old + ekin_new)
    dynvars%rms_gradient = sqrt(rmsg/real(3*natom,wp))

    if (dynamics%integrator == IntegratorLEAP) then
      viri_ke = virial(1,1) + virial(2,2) + virial(3,3)
    else
      viri_ke = virial_sum(1,1) + virial_sum(2,2) + virial_sum(3,3)
    end if 

    dynvars%internal_virial = viri_ke/3.0_wp
    dynvars%virial_kene     = -0.5_wp*viri_ke

    viri_ke_ext = virial_ext(1,1) + virial_ext(2,2) + virial_ext(3,3)
    dynvars%external_virial = viri_ke_ext/3.0_wp

    return

  end subroutine compute_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_genesis
  !> @brief        output dynvars in GENESIS style
  !! @authors      YS, TM, YM
  !! @param[in]    output   : output information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    dynvars  : dynamical variables information
  !! @param[in]    ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_genesis(output, enefunc, dynvars, ensemble)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_ensemble), optional, intent(in) :: ensemble

    ! local variables
    integer,parameter        :: clength=16, flength=4
    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, frmt_res*10, rfrmt*7
    character                :: rfrmt_cont*9,frmt_cont*7
    real(wp)                 :: values(999)
    real(wp)                 :: ene_restraint


    ! common to all FF
    title = 'INFO:       STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    if (present(ensemble)) then
      write(category(ifm),frmt) 'TIME'
      values(ifm) = dynvars%time
      ifm = ifm+1

      write(category(ifm),frmt) 'TOTAL_ENE'
      values(ifm) = dynvars%total_energy
      ifm = ifm+1
    endif

    if (enefunc%morph_flag) then
      if (dynvars%iterations == 1) then
        etitle = .true.
      end if
      write(category(ifm),frmt) 'ITERATION'
      values(ifm) = dynvars%iterations
      ifm = ifm+1
    end if

    write(category(ifm),frmt) 'POTENTIAL_ENE'
    values(ifm) = dynvars%total_pene
    ifm = ifm+1

    if (present(ensemble)) then
      write(category(ifm),frmt) 'KINETIC_ENE'
      values(ifm) = dynvars%total_kene
      ifm = ifm+1
    end if

    write(category(ifm),frmt) 'RMSG'
    values(ifm) = dynvars%rms_gradient
    ifm = ifm+1

    if (.not. present(ensemble)) then
      write(category(ifm),frmt) 'MAXG'
      values(ifm) = dynvars%max_gradient
      ifm = ifm+1
    end if

    if (enefunc%num_bonds > 0 .or. enefunc%num_bonds_quartic > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = dynvars%energy%bond
      ifm = ifm+1
    end if

    if (enefunc%num_angles  > 0 .or. &
        enefunc%num_ureys   > 0 .or. &
        enefunc%num_angflex > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = dynvars%energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = dynvars%energy%urey_bradley
        ifm = ifm+1
      end if
    end if

    if (enefunc%num_dihedrals    > 0 .or. &
        enefunc%num_rb_dihedrals > 0 .or. &
        enefunc%num_dihedflex    > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = dynvars%energy%dihedral
      ifm = ifm+1
    end if

    if (enefunc%forcefield /= ForcefieldCAGO .and.  &
        enefunc%forcefield /= ForcefieldKBGO .and. &
        enefunc%forcefield /= ForcefieldRESIDCG) then

      if (enefunc%num_impropers > 0) then
        write(category(ifm),frmt) 'IMPROPER'
        values(ifm) = dynvars%energy%improper
        ifm = ifm+1
      end if
    end if

    if (enefunc%forcefield == ForcefieldCHARMM   .or. &
        enefunc%forcefield == ForcefieldCHARMM19 .or. &
        enefunc%forcefield == ForcefieldAMBER    .or. &
        enefunc%forcefield == ForcefieldGROAMBER) then

      if (enefunc%num_cmaps > 0) then
        write(category(ifm),frmt) 'CMAP'
        values(ifm) = dynvars%energy%cmap
        ifm = ifm+1
      end if
    end if

    ! ~CG~ 3SPN.2C DNA
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      if (enefunc%num_base_stack > 0) then
        write(category(ifm),frmt) 'BASE_STACK'
        values(ifm) = dynvars%energy%base_stacking
        ifm = ifm+1
      end if
      if (enefunc%cg_DNA_base_pair_calc) then
        write(category(ifm),frmt) 'BASE_PAIR'
        values(ifm) = dynvars%energy%base_pairing
        ifm = ifm+1
      end if
      if (enefunc%cg_DNA_exv_calc) then
        write(category(ifm),frmt) 'DNA_exv'
        values(ifm) = dynvars%energy%cg_DNA_exv
        ifm = ifm+1
      end if
    end if
    
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      if (enefunc%cg_IDR_HPS_calc) then
        write(category(ifm),frmt) 'IDR_HPS'
        values(ifm) = dynvars%energy%cg_IDR_HPS
        ifm = ifm+1
      end if
      if ( enefunc%cg_IDR_KH_calc ) then
        write(category(ifm),frmt) 'IDR_KH'
        values(ifm) = dynvars%energy%cg_IDR_KH
        ifm = ifm+1
      end if
      if (enefunc%cg_KH_calc) then
        write(category(ifm),frmt) 'pro_pro_KH'
        values(ifm) = dynvars%energy%cg_KH_inter_pro
        ifm = ifm+1
      end if
      if (enefunc%cg_pwmcos_calc) then
        write(category(ifm),frmt) 'PWMcos'
        values(ifm) = dynvars%energy%PWMcos
        ifm = ifm+1
      end if
      if (enefunc%cg_pwmcosns_calc) then
        write(category(ifm),frmt) 'PWMcosns'
        values(ifm) = dynvars%energy%PWMcosns
        ifm = ifm+1
      end if
    end if
    
    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldKBGO) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = dynvars%energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = dynvars%energy%noncontact
      ifm = ifm+1

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = dynvars%energy%electrostatic
      ifm = ifm+1

    else if (enefunc%forcefield == ForcefieldRESIDCG) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = dynvars%energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'CG_EXV'
      values(ifm) = dynvars%energy%cg_exv
      ifm = ifm+1

      if ( enefunc%cg_ele_calc) then
        write(category(ifm),frmt) 'ELECT'
        values(ifm) = dynvars%energy%electrostatic
        ifm = ifm+1
      end if

    else if (enefunc%forcefield == ForcefieldSOFT) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = dynvars%energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = dynvars%energy%noncontact
      ifm = ifm+1

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = dynvars%energy%electrostatic
      ifm = ifm+1

    else

      write(category(ifm),frmt) 'VDWAALS'
      values(ifm) = dynvars%energy%van_der_waals
      ifm = ifm+1

      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(category(ifm),frmt) 'DISP-CORR_ENE'
        values(ifm) = dynvars%energy%disp_corr_energy
        ifm = ifm+1
      end if

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = dynvars%energy%electrostatic
      ifm = ifm+1

      if (enefunc%eef1_use .or. enefunc%gbsa_use) then
        write(category(ifm),frmt) 'SOLVATION'
        values(ifm) = dynvars%energy%solvation
        ifm = ifm+1
      end if
    endif

    if (enefunc%morph_flag .and. .not. enefunc%morph_restraint_flag) then
      write(category(ifm),frmt) 'DRMS1'
      values(ifm) = dynvars%energy%drms(1)
      ifm = ifm+1

      write(category(ifm),frmt) 'DRMS2'
      values(ifm) = dynvars%energy%drms(2)
      ifm = ifm+1

      write(category(ifm),frmt) 'MORPH'
      values(ifm) = dynvars%energy%morph
      ifm = ifm+1
    end if

    if (enefunc%ez_membrane_flag) then
      write(category(ifm),frmt) 'EZ-MEMBRANE'
      values(ifm) = dynvars%energy%membrane
      ifm = ifm+1
    end if


    if (enefunc%num_restraintfuncs > 0) then
      ene_restraint =  &
         sum(dynvars%energy%restraint(1:enefunc%num_restraintfuncs))
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

      do i = 1, enefunc%num_restraintfuncs
        write(category(ifm),frmt_res) 'RESTRAINT', i
        values(ifm) = dynvars%energy%restraint(i)
        ifm = ifm+1
      end do
      do i = 1, enefunc%num_restraintfuncs
        write(category(ifm),frmt_res) 'RESTR_CVS', i
        values(ifm) = dynvars%energy%restraint_cv(i)
        ifm = ifm+1
      end do
    end if

    if (enefunc%steered_function > 0) then
      write(category(ifm),frmt) 'SMD_CV'
      values(ifm) = enefunc%target_value
      ifm = ifm+1
    end if
    if (enefunc%target_function > 0) then
      write(category(ifm),frmt) 'TMD_CV'
      values(ifm) = enefunc%target_value
      ifm = ifm+1
    end if

    if (enefunc%num_basins > 1) then
      do i = 1, enefunc%num_basins
        write(category(ifm),frmt_res) 'MUL_BASIN', i
        values(ifm) = dynvars%energy%basin_ratio(i)
        ifm = ifm+1
      end do
      do i = 1, enefunc%num_basins
        write(category(ifm),frmt_res) 'BASIN_ENE', i
        values(ifm) = dynvars%energy%basin_energy(i)
        ifm = ifm+1
      end do

    endif

    if (present(ensemble)) then

      write(category(ifm),frmt) 'TEMPERATURE'
      values(ifm) = dynvars%temperature
      ifm = ifm+1

      write(category(ifm),frmt) 'VOLUME'
      values(ifm) = dynvars%volume
      ifm = ifm+1

      if (output%verbose                    .or. &
          ensemble%ensemble == EnsembleNPT  .or. &
          ensemble%ensemble == EnsembleNPAT .or. &
          ensemble%ensemble == EnsembleNPgT ) then
     
        write(category(ifm),frmt) 'BOXX'
        values(ifm) = dynvars%box_size_x
        ifm = ifm+1
     
        write(category(ifm),frmt) 'BOXY'
        values(ifm) = dynvars%box_size_y
        ifm = ifm+1
     
        write(category(ifm),frmt) 'BOXZ'
        values(ifm) = dynvars%box_size_z
        ifm = ifm+1
     
        write(category(ifm),frmt) 'VIRIAL'
        values(ifm) = dynvars%internal_virial
        ifm = ifm+1
     
        if (enefunc%dispersion_corr == Disp_Corr_EPRESS) then
          write(category(ifm),frmt) 'DISP-CORR_VIR'
          values(ifm) = dynvars%energy%disp_corr_virial
          ifm = ifm+1
        end if
     
        write(category(ifm),frmt) 'PRESSURE'
        values(ifm) = dynvars%internal_pressure
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSXX'
        values(ifm) = dynvars%pressure_xx
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSYY'
        values(ifm) = dynvars%pressure_yy
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSZZ'
        values(ifm) = dynvars%pressure_zz
        ifm = ifm+1
     
      end if
    end if

    if (enefunc%qmmm%do_qmmm) then
      write(category(ifm),frmt) 'QM'
      values(ifm) = dynvars%energy%qm_ene * CONV_UNIT_ENE
      ifm = ifm + 1
    end if

    if (enefunc%spot_use) then
      write(category(ifm),frmt) 'SPOT'
      values(ifm) = dynvars%energy%spot
      ifm = ifm + 1
    end if

    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_pot) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = dynvars%energy%total_gamd
        ifm = ifm+1
      else if (enefunc%gamd%boost_dih) then
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = dynvars%energy%dihedral_gamd
        ifm = ifm+1
      else if (enefunc%gamd%boost_dual) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = dynvars%energy%total_gamd
        ifm = ifm+1
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = dynvars%energy%dihedral_gamd
        ifm = ifm+1
      end if
    end if

    if (etitle) then

      write(DynvarsOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(DynvarsOut,frmt) category(i)
        else
          write(DynvarsOut,frmt_cont) category(i)
        end if
      end do

      write(DynvarsOut,'(A80)') ' --------------- --------------- --------------- --------------- ---------------'
      if (nrep_per_proc == 1) vervose = .false.
      if (.not. vervose) etitle = .false.
      vervose = .false.
    end if

    write(DynvarsOut,'(A5,1x,I10,$)') 'INFO:', dynvars%step

    do i = 1, ifm-1
      if (i == ifm-1) then
        write(DynvarsOut,rfrmt) values(i)
      else
        write(DynvarsOut,rfrmt_cont) values(i)
      end if
    end do
    write(DynvarsOut,'(A)') ''

    return

  end subroutine output_dynvars_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_charmm
  !> @brief        output dynvars in CHARMM style
  !! @authors      YS, CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_charmm(enefunc, dynvars)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wp)                 :: time, totener, totke, energy, temperature
    real(wp)                 :: grms, hfctote, hfcke, ehfcor, virke
    real(wp)                 ::                hbonds, asp, user
    real(wp)                 :: vire, viri, presse, pressi, volume
    real(wp)                 :: boxx, boxy, boxz, pixx, piyy, pizz
    real(wp)                 :: cdihe, cintcr, noe
    real(wp)                 :: posicon, restdist
    integer                  :: i, step

    real(wp),        pointer :: bonds, angles, urey_b, dihedrals, impropers
    real(wp),        pointer :: cmaps, disp_corr
    real(wp),        pointer :: vdwaals, elec
    real(wp),        pointer :: contact, noncontact
    real(wp),        pointer :: base_stack
    real(wp),        pointer :: base_pair
    real(wp),        pointer :: cg_DNA_exv
    real(wp),        pointer :: cg_IDR_HPS
    real(wp),        pointer :: cg_IDR_KH
    real(wp),        pointer :: cg_KH
    real(wp),        pointer :: cg_exv
    real(wp),        pointer :: PWMcos
    real(wp),        pointer :: PWMcosns


    step        = 0
    time        = 0.0_wp
    hbonds      = 0.0_wp
    totener     = 0.0_wp
    totke       = 0.0_wp
    energy      = 0.0_wp
    temperature = 0.0_wp
    asp         = 0.0_wp
    user        = 0.0_wp
    vire        = 0.0_wp
    viri        = 0.0_wp
    presse      = 0.0_wp
    pressi      = 0.0_wp
    volume      = 0.0_wp
    boxx        = 0.0_wp
    boxy        = 0.0_wp
    boxz        = 0.0_wp
    pixx        = 0.0_wp
    piyy        = 0.0_wp
    pizz        = 0.0_wp
    hfctote     = 0.0_wp
    ehfcor      = 0.0_wp
    volume      = 0.0_wp


    step        = dynvars%step
    time        = dynvars%time
    energy      = dynvars%total_pene
    totke       = dynvars%total_kene
    totener     = dynvars%total_energy
    temperature = dynvars%temperature
    grms        = dynvars%rms_gradient
    hfcke       = dynvars%hfc_kene
    virke       = dynvars%virial_kene

    cdihe       = 0.0_wp
    cintcr      = 0.0_wp
    noe         = 0.0_wp

    viri        = dynvars%internal_virial
    vire        = dynvars%external_virial
    pressi      = dynvars%internal_pressure
    presse      = dynvars%external_pressure
    volume      = dynvars%volume
    boxx        = dynvars%box_size_x
    boxy        = dynvars%box_size_y
    boxz        = dynvars%box_size_z
    pixx        = dynvars%pressure_xx 
    piyy        = dynvars%pressure_yy
    pizz        = dynvars%pressure_zz

    bonds      => dynvars%energy%bond
    angles     => dynvars%energy%angle
    urey_b     => dynvars%energy%urey_bradley
    dihedrals  => dynvars%energy%dihedral
    impropers  => dynvars%energy%improper
    base_stack => dynvars%energy%base_stacking
    base_pair  => dynvars%energy%base_pairing
    cg_DNA_exv => dynvars%energy%cg_DNA_exv
    cg_IDR_HPS => dynvars%energy%cg_IDR_HPS
    cg_IDR_KH  => dynvars%energy%cg_IDR_KH
    cg_KH      => dynvars%energy%cg_KH_inter_pro
    cg_exv     => dynvars%energy%cg_exv
    PWMcos     => dynvars%energy%PWMcos
    PWMcosns   => dynvars%energy%PWMcosns
    elec       => dynvars%energy%electrostatic
    vdwaals    => dynvars%energy%van_der_waals
    cmaps      => dynvars%energy%cmap
    contact    => dynvars%energy%contact
    noncontact => dynvars%energy%noncontact
    disp_corr  => dynvars%energy%disp_corr_energy

    restdist = 0.0_wp
    posicon  = 0.0_wp
    do i = 1, enefunc%num_restraintfuncs
      select case (enefunc%restraint_kind(i))

      case (RestraintsFuncPOSI,                         &
            RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

        posicon  = posicon  + dynvars%energy%restraint(i)

      case (RestraintsFuncDIST:RestraintsFuncDISTCOM,   &
            RestraintsFuncANGLE:RestraintsFuncANGLECOM, &
            RestraintsFuncDIHED:RestraintsFuncDIHEDCOM, &
            RestraintsFuncREPUL:RestraintsFuncREPULCOM, &
            RestraintsFuncFB:RestraintsFuncFBCOM)

        restdist = restdist + dynvars%energy%restraint(i)

      end select
    end do


    ! write title if necessary
    !
    if (etitle) then
      write(DynvarsOut,'(A)') 'Output_Energy> CHARMM_Style is used'
      write(DynvarsOut,'(A)') ' '
      write(DynvarsOut,'(A79)') 'DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature'
      write(DynvarsOut,'(A79)') 'DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe'
      write(DynvarsOut,'(A79)') 'DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers'

      if (enefunc%num_cmaps > 0) then
        write(DynvarsOut,'(A79)') 'DYNA CROSS:           CMAPs                                                    '
      end if
      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(DynvarsOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      endif
      if (enefunc%num_base_stack > 0) then
        write(DynvarsOut,'(A79)') 'DYNA DNA:     BASE_STACKING                                                    '
      end if
      if (enefunc%cg_DNA_base_pair_calc) then
        write(DynvarsOut,'(A79)') 'DYNA DNA:     BASE_PAIRING                                                     '
      end if
      if (enefunc%cg_DNA_exv_calc) then
        write(DynvarsOut,'(A79)') 'DYNA DNA:               EXV                                                    '
      end if
      if (enefunc%cg_IDR_HPS_calc) then
        write(DynvarsOut,'(A79)') 'DYNA IDR:               HPS                                                    '
      end if
      if (enefunc%cg_IDR_KH_calc) then
        write(DynvarsOut,'(A79)') 'DYNA IDR:                KH                                                    '
      end if
      if (enefunc%cg_KH_calc) then
        write(DynvarsOut,'(A79)') 'DYNA pro-pro:            KH                                                    '
      end if
      if (enefunc%cg_pwmcos_calc) then
        write(DynvarsOut,'(A79)') 'DYNA protein-DNA:  PWMcos                                                      '
      end if
      if (enefunc%cg_pwmcosns_calc) then
        write(DynvarsOut,'(A79)') 'DYNA protein-DNA:  PWMcosns                                                    '
      end if

      if (enefunc%forcefield == ForcefieldAAGO .or. &
          enefunc%forcefield == ForcefieldCAGO .or. &
          enefunc%forcefield == ForcefieldRESIDCG .or. &
          enefunc%forcefield == ForcefieldKBGO) then
        write(DynvarsOut,'(A79)') 'DYNA GO:            CONTACT     NCONTACT                                       '
        write(DynvarsOut,'(A79)') 'DYNA CG:                EXV                                                    '
      else if (enefunc%forcefield == ForcefieldSOFT) then
        write(DynvarsOut,'(A79)') 'DYNA SOFT:          CONTACT     NCONTACT        ELEC                           '
      else
        write(DynvarsOut,'(A79)') 'DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER'
      end if

      if (enefunc%restraint_flag) then
        write(DynvarsOut,'(A79)') 'DYNA CONSTR:       HARMonic    CDIHedral        CIC     RESDistance       NOE'
      end if

      write(DynvarsOut,'(A79)') 'DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme'
      write(DynvarsOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
      if (nrep_per_proc == 1) vervose = .false.
      if (.not. vervose) etitle = .false.
      vervose = .false.
    end if


    ! write energy in CHARMM-style
    !
    write(DynvarsOut,'(A5,I9,5F13.5)') 'DYNA>', step, time, totener, totke, energy, temperature
    write(DynvarsOut,'(A14,5F13.5)')   'DYNA PROP>    ', grms, hfctote, hfcke, ehfcor, virke
    write(DynvarsOut,'(A14,5F13.5)')   'DYNA INTERN>  ', bonds, angles, urey_b, dihedrals, impropers

    if (enefunc%num_cmaps > 0) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA CROSS>   ', cmaps
    end if
    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA DISP>    ', disp_corr
    end if
    if (enefunc%num_base_stack > 0) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA DNA>     ', base_stack
    end if
    if (enefunc%cg_DNA_base_pair_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA DNA>     ', base_pair
    end if
    if (enefunc%cg_DNA_exv_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA DNA>     ', cg_DNA_exv
    end if
    if (enefunc%cg_IDR_HPS_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA IDR>     ', cg_IDR_HPS
    end if
    if (enefunc%cg_IDR_KH_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA IDR>     ', cg_IDR_KH
    end if
    if (enefunc%cg_KH_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA pro-pro> ', cg_KH
    end if
    if (enefunc%cg_pwmcos_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA PWMcos>  ', PWMcos
    end if
    if (enefunc%cg_pwmcosns_calc) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA PWMcosns>', PWMcosns
    end if

    if (enefunc%forcefield == ForcefieldAAGO .or.  &
        enefunc%forcefield == ForcefieldCAGO .or.  &
        enefunc%forcefield == ForcefieldRESIDCG .or.  &
        enefunc%forcefield == ForcefieldKBGO) then
      write(DynvarsOut,'(A14,2F13.5)')   'DYNA GO>      ', contact, noncontact
      write(DynvarsOut,'(A14,F13.5)')    'DYNA CG>      ', cg_exv
    else if (enefunc%forcefield == ForcefieldSOFT) then
      write(DynvarsOut,'(A14,3F13.5)')   'DYNA SOFT>    ', contact, noncontact, elec
    else
      write(DynvarsOut,'(A14,5F13.5)')   'DYNA EXTERN>  ', vdwaals, elec, hbonds, asp, user
    end if
!    write(DynvarsOut,'(A14,5F13.5)')   'DYNA IMAGES>  ', imnbvdw, imelec, imhbnd, rxnfield, extelec
!    write(DynvarsOut,'(A14,5F13.5)')   'DYNA EWALD>   ', ewksum, ewself, ewexcl, ewqcor, ewutil

    if (enefunc%restraint_flag) then
      write(DynvarsOut,'(A14,5F13.5)')   'DYNA CONSTR>  ', posicon, cdihe, cintcr, restdist, noe
    end if
    write(DynvarsOut,'(A14,5F13.5)')   'DYNA PRESS>   ', vire, viri, presse, pressi, volume
    write(DynvarsOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'

!    write(DynvarsOut,'(A10,3(A8,F11.5))') '      DYNA',' A     =',boxx,' B     =',boxy,' C     =',boxz
!    write(DynvarsOut,'(A10,3(A8,F11.2))') '      DYNA',' PIXX  =',pixx,' PIYY  =',piyy,' PIZZ  =',pizz

    return

  end subroutine output_dynvars_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_namd
  !> @brief        output dynvars in NAMD style
  !! @authors      YS, CK
  !! @param[in]    enefunc : potential energy function information
  !! @param[inout] dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_namd(enefunc, dynvars)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wp)                 :: ebond, eangle
    real(wp)                 ::               misc, kinetic
    real(wp)                 :: total, temp, total2, total3, tempavg
    real(wp)                 :: pressure, gpressure, volume
    real(wp)                 :: pressavg, gpressavg
    real(wp)                 :: posicon, restdist, eboundary
    integer                  :: i,nstep

    real(wp),        pointer :: edihed, eimprp
    real(wp),        pointer :: eelect, evdw
    real(wp),        pointer :: ebase_stack
    real(wp),        pointer :: ebase_pair
    real(wp),        pointer :: ecg_DNA_exv
    real(wp),        pointer :: ecg_IDR_HPS
    real(wp),        pointer :: ecg_IDR_KH
    real(wp),        pointer :: ecg_KH
    real(wp),        pointer :: ecg_exv
    real(wp),        pointer :: ePWMcos
    real(wp),        pointer :: ePWMcosns


    ! restraint energy
    !
    restdist = 0.0_wp
    posicon  = 0.0_wp
    do i = 1, enefunc%num_restraintfuncs
      select case (enefunc%restraint_kind(i))

      case (RestraintsFuncPOSI,                         &
            RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

        posicon = posicon + dynvars%energy%restraint(i)

      case (RestraintsFuncDIST:RestraintsFuncDISTCOM,   &
            RestraintsFuncANGLE:RestraintsFuncANGLECOM, &
            RestraintsFuncDIHED:RestraintsFuncDIHEDCOM, &
            RestraintsFuncREPUL:RestraintsFuncREPULCOM, &
            RestraintsFuncFB:RestraintsFuncFBCOM)

        restdist = restdist + dynvars%energy%restraint(i)

      end select
    end do

    nstep       =  dynvars%step
    ebond       =  dynvars%energy%bond  + restdist
    eangle      =  dynvars%energy%angle + dynvars%energy%urey_bradley
    edihed      => dynvars%energy%dihedral
    eimprp      => dynvars%energy%improper
    eelect      => dynvars%energy%electrostatic
    evdw        => dynvars%energy%van_der_waals
    ebase_stack => dynvars%energy%base_stacking
    ebase_pair  => dynvars%energy%base_pairing
    ecg_DNA_exv => dynvars%energy%cg_DNA_exv
    ecg_IDR_HPS => dynvars%energy%cg_IDR_HPS
    ecg_IDR_KH  => dynvars%energy%cg_IDR_KH
    ecg_KH      => dynvars%energy%cg_KH_inter_pro
    ecg_exv     => dynvars%energy%cg_exv
    ePWMcos     => dynvars%energy%PWMcos
    ePWMcosns   => dynvars%energy%PWMcosns
    eboundary   =  posicon
    
    kinetic   =  dynvars%total_kene
    total     =  dynvars%total_energy
    temp      =  dynvars%temperature

    ! write title if necessary
    !
    if (etitle) then
      write(DynvarsOut,'(A)') 'Output_Energy> NAMD_Style is used'
      write(DynvarsOut,'(A)') ' '
      write(DynvarsOut,'(A75)') 'ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP'
      write(DynvarsOut,'(A75)') '  BASE_STACKING   BASE_PAIRING        DNA_EXV         PWMCOS       PWMcosns'
      write(DynvarsOut,'(A75)') '        IDR_HPS         IDR_KH     PRO_PRO_KH            EXV               '
      write(DynvarsOut,'(A75)') '          ELECT            VDW       BOUNDARY           MISC        KINETIC'
      write(DynvarsOut,'(A75)') '          TOTAL           TEMP         TOTAL2         TOTAL3        TEMPAVG'
      write(DynvarsOut,'(A75)') '       PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG'
      write(DynvarsOut,'(A)') ' '
      etitle = .false.
    end if


    ! write energy in NAMD-style
    !
    write(DynvarsOut,'(A7,I8,4F15.4)')'ENERGY:', nstep, ebond, eangle, edihed, eimprp
    write(DynvarsOut,'(5F15.4)') ebase_stack, ebase_pair, ecg_DNA_exv, ePWMcos, ePWMcosns
    write(DynvarsOut,'(4F15.4)') ecg_IDR_HPS, ecg_IDR_KH, ecg_KH, ecg_exv
    write(DynvarsOut,'(5F15.4)') eelect, evdw, eboundary, misc, kinetic
    write(DynvarsOut,'(5F15.4)') total, temp, total2, total3, tempavg
    write(DynvarsOut,'(5F15.4)') pressure, gpressure, volume, pressavg, gpressavg
    write(DynvarsOut,'(A)') ' '


    return

  end subroutine output_dynvars_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_gromacs
  !> @brief        output dynvars in GROMACS style
  !! @authors      YS, CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_gromacs(enefunc, dynvars)

     ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wp)                 :: time, bonds, angles, urey_b
    real(wp)                 :: dihedrals, impropers, vdwaals, elec
    real(wp)                 :: restdist, posicon
    real(wp)                 :: energy, totke, totener
    real(wp)                 :: temperature, pressi, presse
    real(wp)                 :: disp_corr, pres_dc
    real(wp)                 :: noncontact, contact
    real(wp)                 :: base_stack
    real(wp)                 :: base_pair
    real(wp)                 :: cg_DNA_exv
    real(wp)                 :: cg_IDR_HPS
    real(wp)                 :: cg_IDR_KH
    real(wp)                 :: cg_KH
    real(wp)                 :: cg_exv
    real(wp)                 :: PWMcos
    real(wp)                 :: PWMcosns
    integer                  :: i, step


    step        = dynvars%step
    time        = dynvars%time
    bonds       = dynvars%energy%bond          * CAL2JOU
    angles      = dynvars%energy%angle         * CAL2JOU
    urey_b      = dynvars%energy%urey_bradley  * CAL2JOU
    dihedrals   = dynvars%energy%dihedral      * CAL2JOU
    impropers   = dynvars%energy%improper      * CAL2JOU
    vdwaals     = dynvars%energy%van_der_waals * CAL2JOU
    elec        = dynvars%energy%electrostatic * CAL2JOU
    noncontact  = dynvars%energy%noncontact    * CAL2JOU
    contact     = dynvars%energy%contact       * CAL2JOU
    base_stack  = dynvars%energy%base_stacking * CAL2JOU
    base_pair   = dynvars%energy%base_pairing  * CAL2JOU
    cg_DNA_exv  = dynvars%energy%cg_DNA_exv    * CAL2JOU
    cg_IDR_HPS  = dynvars%energy%cg_IDR_HPS    * CAL2JOU
    cg_IDR_KH   = dynvars%energy%cg_IDR_KH     * CAL2JOU
    cg_KH       = dynvars%energy%cg_KH_inter_pro * CAL2JOU
    cg_exv      = dynvars%energy%cg_exv        * CAL2JOU
    PWMcos      = dynvars%energy%PWMcos        * CAL2JOU
    PWMcosns    = dynvars%energy%PWMcosns      * CAL2JOU

    restdist = 0.0_wp
    posicon  = 0.0_wp
    do i = 1, enefunc%num_restraintfuncs
      select case (enefunc%restraint_kind(i))

      case (RestraintsFuncPOSI,                         &
            RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

        posicon  = posicon  + dynvars%energy%restraint(i)

      case (RestraintsFuncDIST:RestraintsFuncDISTCOM,   &
            RestraintsFuncANGLE:RestraintsFuncANGLECOM, &
            RestraintsFuncDIHED:RestraintsFuncDIHEDCOM, &
            RestraintsFuncREPUL:RestraintsFuncREPULCOM, &
            RestraintsFuncFB:RestraintsFuncFBCOM)

        restdist = restdist + dynvars%energy%restraint(i)

      end select
    end do
    posicon  = posicon  * CAL2JOU
    restdist = restdist * CAL2JOU

    energy      = dynvars%total_pene           * CAL2JOU
    totke       = dynvars%total_kene           * CAL2JOU
    totener     = dynvars%total_energy         * CAL2JOU
    temperature = dynvars%temperature
    pressi      = dynvars%internal_pressure    * ATM2BAR
    presse      = dynvars%external_pressure    * ATM2BAR
    disp_corr   = dynvars%energy%disp_corr_energy * CAL2JOU
    pres_dc     = dynvars%energy%disp_corr_virial * P_ATMOS/dynvars%volume &
                  * ATM2BAR


    write(DynvarsOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(DynvarsOut,'(I15,2F15.5)') &
         step, time, 0.0_wp
    write(DynvarsOut,'(A)') &
         ' '
    write(DynvarsOut,'(A)') &
         '   Energies (kJ/mol)'

    write(DynvarsOut,'(5A15)') &
         'Bond', 'Angle', 'Urey-bradley', 'Dihedral', 'Improper Dih.'
    write(DynvarsOut,'(5ES15.5E2)') &
         bonds,angles,urey_b,dihedrals,impropers

    if (enefunc%forcefield == ForcefieldRESIDCG) then
      if (enefunc%num_base_stack > 0) then
        write(DynvarsOut,'(A15)') 'Base Stacking'
        write(DynvarsOut,'(1ES15.5E2)') &
            base_stack
      end if
      if (enefunc%cg_DNA_base_pair_calc) then
        write(DynvarsOut,'(A15)') 'Base Pairing'
        write(DynvarsOut,'(1ES15.5E2)') &
            base_pair
      end if
      if (enefunc%cg_DNA_exv_calc) then
        write(DynvarsOut,'(A15)') 'DNA exv'
        write(DynvarsOut,'(1ES15.5E2)') cg_DNA_exv
      end if
      if (enefunc%cg_IDR_HPS_calc) then
        write(DynvarsOut,'(A15)') 'IDR HPS'
        write(DynvarsOut,'(1ES15.5E2)') cg_IDR_HPS
      end if
      if (enefunc%cg_IDR_KH_calc) then
        write(DynvarsOut,'(A15)') 'IDR KH'
        write(DynvarsOut,'(1ES15.5E2)') cg_IDR_KH
      end if
      if (enefunc%cg_KH_calc) then
        write(DynvarsOut,'(A15)') 'pro-pro KH'
        write(DynvarsOut,'(1ES15.5E2)') cg_KH
      end if
      if (enefunc%cg_pwmcos_calc) then
        write(DynvarsOut,'(A15)') 'PWMcos'
        write(DynvarsOut,'(1ES15.5E2)') PWMcos
      end if
      if (enefunc%cg_pwmcosns_calc) then
        write(DynvarsOut,'(A15)') 'PWMcosns'
        write(DynvarsOut,'(1ES15.5E2)') PWMcosns
      end if
      write(DynvarsOut,'(A15)') 'general exv'
      write(DynvarsOut,'(1ES15.5E2)') cg_exv
    end if

    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldRESIDCG .or. &
        enefunc%forcefield == ForcefieldKBGO) then
      write(DynvarsOut,'(5A15)') &
           ' Contact  ', ' Noncontact    ', 'Position Rest.', 'Potential',   &
           'Kinetic En.'
      write(DynvarsOut,'(5ES15.5E2)') &
           contact,noncontact,posicon,energy,totke

    else if (enefunc%forcefield == ForcefieldSOFT) then
      write(DynvarsOut,'(5A15)') &
           'Contact', 'Noncontact', 'Coulomb(SR)', 'Position Rest.', 'Potential'
      write(DynvarsOut,'(5ES15.5E2)') &
           contact,noncontact,elec,posicon,energy
      write(DynvarsOut,'(1A15)') &
           'Kinetic En.'
      write(DynvarsOut,'(1ES15.5E2)') &
           totke

    else if (enefunc%dispersion_corr == Disp_Corr_NONE) then
      write(DynvarsOut,'(5A15)') &
           'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Position Rest.', 'Potential',   &
           'Kinetic En.'
      write(DynvarsOut,'(5ES15.5E2)') &
           vdwaals,elec,posicon,energy,totke

    else
      write(DynvarsOut,'(5A15)') &
         'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Disper. corr.', 'Position Rest.', &
         'Pres. DC (bar)'
      write(DynvarsOut,'(5ES15.5E2)') &
           vdwaals,elec,disp_corr,posicon,pres_dc
      write(DynvarsOut,'(3A15)') &
         'Position Rest.', 'Potential', 'Kinetic En.'
      write(DynvarsOut,'(3ES15.5E2)') &
          posicon,energy,totke

    end if

    write(DynvarsOut,'(5A15)') &
         'Total Energy', 'Temperature', 'Pressure(int.)', 'Pressure(ext.)'
    write(DynvarsOut,'(5ES15.5E2)') &
         totener,temperature,pressi,presse

    write(DynvarsOut,'(A)')  ' '

    return

  end subroutine output_dynvars_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_pressure
  !> @brief        compute internal pressure
  !! @authors      TM
  !! @param[in]    num_degree  : number of degrees of freedom
  !! @param[in]    volume      : volume (A^3)
  !! @param[in]    temperature : temperature (K)
  !! @param[in]    virial      : internal virial
  !! @param[out]   pressure    : pressure (atm)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_pressure(num_degree, volume, temperature, virial, pressure)

    ! formal arguments
    integer,                 intent(in)    :: num_degree
    real(wp),                intent(in)    :: volume
    real(wp),                intent(in)    :: temperature
    real(wp),                intent(in)    :: virial(3,3)
    real(wp),                intent(out)   :: pressure

    ! local variables
    real(wp)                 :: viri_ke


    ! compute internal pressure
    !
    viri_ke  = virial(1,1) + virial(2,2) + virial(3,3)
    pressure = P_ATMOS / (3.0_wp * volume)                      &
              * (num_degree * KBOLTZ * temperature + viri_ke)

    return

  end subroutine compute_pressure

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate_velocity
  !> @brief        generate velocities
  !! @authors      TM
  !! @param[in]    temperature : initial temperature
  !! @param[in]    natom       : number of atoms
  !! @param[in]    mass        : atom mass
  !! @param[inout] iseed       : random number seed
  !! @param[inout] velocity    : atom velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine generate_velocity(temperature, natom, mass, iseed, velocity)

    ! formal arguments
    real(wp),                intent(in)    :: temperature
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: mass(:)
    integer,                 intent(inout) :: iseed
    real(wp),                intent(out)   :: velocity(:,:)

    ! local variables
    real(wp)                 :: alpha, beta, sigma, kBT
    integer                  :: i


    if (main_rank .or. replica_main_rank) then
      write(DynvarsOut,'(A)') 'Generate_Velocity> Generate initial velocities'
      write(DynvarsOut,'(A20,I10)')   '  iseed           = ', iseed
      write(DynvarsOut,'(A20,F10.3)') '  temperature     = ', temperature
      write(DynvarsOut,'(A)') ' '
    end if

    ! generate velocities
    !   Gaussian distribution was obtained by the Box-Muller method
    !
    kBT = KBOLTZ * temperature

    do i = 1, natom
      sigma = sqrt( kBT / mass(i) )

      alpha = sqrt( -2.0_wp * log(random_get_legacy(iseed)) )
      beta  = cos( 2.0_wp * PI * random_get_legacy(iseed) )
      velocity(1,i) = sigma * alpha * beta

      alpha = sqrt( -2.0_wp * log(random_get_legacy(iseed)) )
      beta  = cos( 2.0_wp * PI * random_get_legacy(iseed) )
      velocity(2,i) = sigma * alpha * beta

      alpha = sqrt( -2.0_wp * log(random_get_legacy(iseed)) )
      beta  = cos( 2.0_wp * PI * random_get_legacy(iseed) )
      velocity(3,i) = sigma * alpha * beta
    end do

    return 

  end subroutine generate_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    stop_trans_rotation
  !> @brief        remove trans and rotational motion about the center of mass
  !! @authors      TM, KY
  !! @param[in]    natom      : number of atoms
  !! @param[in]    mass       : atom mass
  !! @param[in]    stop_trans : flag for stop translational motion
  !! @param[in]    stop_rot   : flag for stop rotational motion
  !! @param[in]    coord      : atom coordinates
  !! @param[in]    fixatm     : fixatom if true
  !! @param[inout] velocity   : atom velocities
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine stop_trans_rotation(natom, mass, stop_trans, stop_rot, coord, &
                                 fixatm, velocity)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: mass(:)
    logical,                 intent(in)    :: stop_trans
    logical,                 intent(in)    :: stop_rot
    real(wp),                intent(in)    :: coord(:,:)
    logical,                 intent(in)    :: fixatm(:)
    real(wp),                intent(inout) :: velocity(:,:)

    ! local variables
    real(wp)                 :: inertia(3,3), inv_inertia(3,3)
    real(wp)                 :: total_mass, ccom(3), vcom(3)
    real(wp)                 :: lx, ly, lz
    real(wp)                 :: c(3), v(3), omega(3)
    real(wp)                 :: xx, xy, xz, yy, yz, zz, x, y, z
    integer                  :: i, j, k


    ! calculate coordinates and velocities about the center of mass
    !
    ccom(1:3)  = 0.0_wp
    vcom(1:3)  = 0.0_wp
    total_mass = 0.0_wp
    do i = 1, natom
      if (fixatm(i)) cycle
      ccom(1:3)  = ccom(1:3) + coord   (1:3,i)*mass(i)
      vcom(1:3)  = vcom(1:3) + velocity(1:3,i)*mass(i)
      total_mass = total_mass + mass(i)
    end do
    ccom(1:3) = ccom(1:3)/total_mass
    vcom(1:3) = vcom(1:3)/total_mass

    ! calculate the angular momentum (L) about the center of mass
    !
    lx = 0.0_wp
    ly = 0.0_wp
    lz = 0.0_wp
    do i = 1, natom
      if (fixatm(i)) cycle
      c(1:3) = coord   (1:3,i) - ccom(1:3)
      v(1:3) = velocity(1:3,i) - vcom(1:3)
      lx = lx + (c(2)*v(3) - c(3)*v(2))*mass(i)
      ly = ly + (c(3)*v(1) - c(1)*v(3))*mass(i)
      lz = lz + (c(1)*v(2) - c(2)*v(1))*mass(i)
    end do

    ! calculate the inertia tensor (I)
    !
    xx = 0.0_wp
    xy = 0.0_wp
    xz = 0.0_wp
    yy = 0.0_wp
    yz = 0.0_wp
    zz = 0.0_wp
    do i = 1, natom
      if (fixatm(i)) cycle
      c(1:3) = coord(1:3,i) - ccom(1:3)
      xx = xx + c(1)*c(1)*mass(i)
      xy = xy + c(1)*c(2)*mass(i)
      xz = xz + c(1)*c(3)*mass(i)
      yy = yy + c(2)*c(2)*mass(i)
      yz = yz + c(2)*c(3)*mass(i)
      zz = zz + c(3)*c(3)*mass(i)
    end do

    inertia(1,1) =   yy + zz
    inertia(1,2) = - xy
    inertia(1,3) = - xz
    inertia(2,1) = - xy
    inertia(2,2) =   xx + zz
    inertia(2,3) = - yz
    inertia(3,1) = - xz
    inertia(3,2) = - yz
    inertia(3,3) =   xx + yy

    ! calculate the inverse matrix of the inertia tensor
    !
    call compute_inverse_matrix(3, inertia, inv_inertia)

    ! calculate the angular velocity (OMG)
    !
    omega(1) = inv_inertia(1,1)*lx + inv_inertia(1,2)*ly + inv_inertia(1,3)*lz
    omega(2) = inv_inertia(2,1)*lx + inv_inertia(2,2)*ly + inv_inertia(2,3)*lz
    omega(3) = inv_inertia(3,1)*lx + inv_inertia(3,2)*ly + inv_inertia(3,3)*lz

    ! remove translational motion
    !
    if (stop_trans) then
      do i = 1, natom
        if (fixatm(i)) cycle
        velocity(1:3,i) = velocity(1:3,i) - vcom(1:3)
      end do
    end if

    ! remove rotational motion
    !
    if (stop_rot) then
      do i = 1, natom
        if (fixatm(i)) cycle
        c(1:3) = coord(1:3,i) - ccom(1:3)
        velocity(1,i) = velocity(1,i) - (omega(2)*c(3) - omega(3)*c(2))
        velocity(2,i) = velocity(2,i) - (omega(3)*c(1) - omega(1)*c(3))
        velocity(3,i) = velocity(3,i) - (omega(1)*c(2) - omega(2)*c(1))
      end do
    end if

    return

  end subroutine stop_trans_rotation

end module at_dynvars_mod
