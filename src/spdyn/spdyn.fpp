!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  SPDYN
!! @brief   Molecular dynamics Simulation of BioMolecules using
!!          Spacial Decomposition Scheme
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Yuji Sugita (YS), 
!!          Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

program spdyn

  use sp_setup_mpi_mod
  use sp_md_vverlet_mod
  use sp_md_leapfrog_mod
  use sp_dynamics_mod
  use sp_minimize_mod
  use sp_setup_spdyn_mod
  use sp_control_mod
  use sp_energy_mod
  use sp_energy_pme_mod
  use sp_communicate_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_output_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_remd_str_mod
  use sp_rpath_str_mod
  use sp_remd_mod
  use sp_rpath_mod
  use molecules_str_mod
  use fileio_control_mod
  use hardwareinfo_mod
  use timers_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

  integer             :: genesis_run_mode
  integer, parameter  :: GenesisMD    = 1
  integer, parameter  :: GenesisMIN   = 2
  integer, parameter  :: GenesisREMD  = 3
  integer, parameter  :: GenesisRPATH = 4
  integer, parameter  :: GenesisFEPMD   = 5
  integer, parameter  :: GenesisFEPMIN  = 6
  integer, parameter  :: GenesisFEPREMD = 7

  ! local variables
  character(MaxFilename)      :: ctrl_filename
  type(s_ctrl_data)           :: ctrl_data
  type(s_molecule)            :: molecule
  type(s_enefunc)             :: enefunc
  type(s_dynvars)             :: dynvars
  type(s_pairlist)            :: pairlist
  type(s_boundary)            :: boundary
  type(s_constraints)         :: constraints
  type(s_ensemble)            :: ensemble
  type(s_dynamics)            :: dynamics
  type(s_minimize)            :: minimize
  type(s_output)              :: output
  type(s_domain)              :: domain
  type(s_comm)                :: comm
  type(s_remd)                :: remd
  type(s_rpath)               :: rpath
  integer                     :: omp_get_max_threads, i
  real(dp)                    :: sas, eae
  type(s_alchemy)             :: alchemy

#ifdef HAVE_MPI_GENESIS
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, my_world_rank, ierror)
  call mpi_comm_size(mpi_comm_world, nproc_world,   ierror)
  main_rank = (my_world_rank == 0)
#else
  my_world_rank = 0
  nproc_world   = 1
  main_rank     = .true.
#endif

#ifdef OMP
  nthread = omp_get_max_threads()
#else
  nthread = 1
#endif
 
  ! show usage
  !
  call usage(ctrl_filename)


  ! get run mode from control file
  !
  call get_genesis_mode(ctrl_filename, genesis_run_mode)


  ! run genesis
  !
  call domain_decomposition_genesis(ctrl_filename, genesis_run_mode)

#ifdef HAVE_MPI_GENESIS
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_finalize(ierror)
#endif

  stop

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_genesis_mode
  !> @brief        get genesis run mode
  !! @authors      TM
  !! @param[in]    ctrl_filename    : control file name
  !! @param[out]   genesis_run_mode : run MD, MIN, REMD, RPATH
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_genesis_mode(ctrl_filename, genesis_run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename
    integer,                 intent(inout) :: genesis_run_mode


    if (find_ctrlfile_section(ctrl_filename, 'ALCHEMY')) then
      if (find_ctrlfile_section(ctrl_filename, 'MINIMIZE')) then
        genesis_run_mode = GenesisFEPMIN
      else if (find_ctrlfile_section(ctrl_filename, 'REMD')) then
        genesis_run_mode = GenesisFEPREMD
      else if (find_ctrlfile_section(ctrl_filename, 'DYNAMICS')) then
        genesis_run_mode = GenesisFEPMD
      else
        call error_msg('Get_Genesis_Mode> ERROR: Unknown control file format.')
      end if
    else  if (find_ctrlfile_section(ctrl_filename, 'REMD')) then
      genesis_run_mode = GenesisREMD

    else if (find_ctrlfile_section(ctrl_filename, 'RPATH')) then
      genesis_run_mode = GenesisRPATH

    else if (find_ctrlfile_section(ctrl_filename, 'DYNAMICS')) then
      genesis_run_mode = GenesisMD

    else if (find_ctrlfile_section(ctrl_filename, 'MINIMIZE')) then
      genesis_run_mode = GenesisMIN

    else
      call error_msg('Get_Genesis_Mode> ERROR: Unknown control file format.')

    end if

    return

  end subroutine get_genesis_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    domain_decomposition_genesis
  !> @brief        run genesis using domain decomposition scheme
  !! @authors      TM
  !! @param[in]    ctrl_filename   : control file name
  !! @param[in]    genesis_run_mod : run MD, MIN, REMD, RPATH
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine domain_decomposition_genesis(ctrl_filename, genesis_run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename
    integer,                 intent(in)    :: genesis_run_mode

#ifdef USE_GPU
    integer :: my_device_id
#endif


    ! set timer
    !
    call timer(TimerTotal, TimerOn)

    ! [Step0] Architecture & Compiler information
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP0] Architecture and Compiler Information'
      write(MsgOut,'(A)') ' '

      call hw_information

#ifdef USE_GPU
    else
      ! assign GPU
      call assign_gpu(my_device_id)
#endif /* USE_GPU */
    end if

    ! [Step1] Read control file
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP1] Read Control Parameters'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD,GenesisFEPMD)

      call control_md  (ctrl_filename, ctrl_data)

    case (GenesisMIN,GenesisFEPMIN)

      call control_min (ctrl_filename, ctrl_data)

    case (GenesisREMD,GenesisFEPREMD)

      call control_remd(ctrl_filename, ctrl_data)

    case (GenesisRPATH)

      call control_rpath(ctrl_filename, ctrl_data)

    end select


    ! [Step2] Setup MPI
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2] Setup MPI'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD,GenesisMIN,GenesisFEPMD,GenesisFEPMIN)

      call setup_mpi_md  (ctrl_data%ene_info)

    case (GenesisREMD,GenesisFEPREMD)

      call setup_mpi_remd(ctrl_data%ene_info, ctrl_data%rep_info, &
                          ctrl_data%bound_info)

    case (GenesisRPATH)

      call setup_mpi_rpath(ctrl_data%ene_info, ctrl_data%rpath_info)

    end select

    ! [Step3] Set relevant variables and structures 
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD)

      call setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                         dynvars, dynamics, constraints, ensemble, boundary, &
                         domain, comm)

    case (GenesisMIN)

      call setup_spdyn_min(ctrl_data, output, molecule, enefunc, pairlist,   &
                         dynvars, minimize, constraints, boundary, domain, comm)

    case (GenesisREMD)

      call setup_spdyn_remd(ctrl_data, output, molecule, enefunc, pairlist,  &
                          dynvars, dynamics, constraints, ensemble, boundary,&
                          domain, comm, remd)

    case (GenesisRPATH)

      call setup_spdyn_rpath(ctrl_data, output, molecule, enefunc, pairlist,  &
                          dynvars, dynamics, constraints, ensemble, boundary,&
                          domain, comm, rpath)

    case (GenesisFEPMD)

      call setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,  &
                         dynvars, dynamics, constraints, ensemble, boundary, &
                         domain, comm, alchemy)

    case (GenesisFEPMIN)

      call setup_spdyn_min(ctrl_data, output, molecule, enefunc, pairlist, &
                         dynvars, minimize, constraints, boundary, domain,     &
                         comm, alchemy)

    case (GenesisFEPREMD)

      call setup_spdyn_remd(ctrl_data, output, molecule, enefunc, pairlist,&
                          dynvars, dynamics, constraints, ensemble, boundary,&
                          domain, comm, remd, alchemy)

    end select


    ! [Step4] Compute single point energy for molecules
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
      write(MsgOut,'(A)') ' '
    end if

    call compute_energy(domain, enefunc, pairlist, boundary, domain%coord, &
                        .false., .true., .true., .true., & 
                        enefunc%nonb_limiter,            &
                        dynvars%energy,                  &
                        domain%atmcls_pbc,               &
                        domain%translated,               &
                        domain%force,                    &
                        domain%force_long,               &
                        domain%force_omp,                &
                        domain%force_pbc,                &
                        domain%virial_cellpair,          &
                        dynvars%virial,                  &
                        dynvars%virial_long,             &
                        dynvars%virial_extern)

    call output_energy(dynvars%step, enefunc, dynvars%energy)


    ! [Step5] Perform MD/REMD/RPATH simulation or Energy minimization
    !
#ifdef HAVE_MPI_GENESIS
    call mpi_barrier(mpi_comm_world, ierror)
#endif
    call timer(TimerDynamics, TimerOn)

    select case (genesis_run_mode)

    case (GenesisMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Molecular Dynamics Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_md(output, domain, enefunc, dynvars, dynamics, pairlist, &
                  boundary, constraints, ensemble, comm, remd)
    
    case (GenesisMIN)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Energy Minimization'
        write(MsgOut,'(A)') ' '
      end if

      call run_min(output, domain, enefunc, dynvars, minimize,       &
                     pairlist, boundary, constraints, comm)

    case (GenesisREMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Replica-Exchange MD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_remd(output, domain, enefunc, dynvars, dynamics, pairlist, &
                    boundary, constraints, ensemble, comm, remd)

    case (GenesisRPATH)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Replica Path MD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_rpath(output, domain, enefunc, dynvars, dynamics, pairlist, &
                    boundary, constraints, ensemble, comm, rpath, remd)

    case (GenesisFEPMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform FEP MD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_md_fep(output, domain, enefunc, dynvars, dynamics, pairlist, &
                  boundary, constraints, ensemble, comm, remd, alchemy)
    
    case (GenesisFEPMIN)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform FEP Energy Minimization'
        write(MsgOut,'(A)') ' '
      end if

      call run_min(output, domain, enefunc, dynvars, minimize,       &
                   pairlist, boundary, constraints, comm)

    case (GenesisFEPREMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform FEP REMD Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_remd_fep(output, domain, enefunc, dynvars, dynamics, pairlist, &
                    boundary, constraints, ensemble, comm, remd, alchemy)

    end select

    call timer(TimerDynamics, TimerOff)


    ! [Step6] Deallocate arrays
    !
    if (main_rank) then
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
      write(MsgOut,'(A)') ' '
    end if

    call dealloc_pme(enefunc)
    call dealloc_constraints_all(constraints)
    call dealloc_boundary_all   (boundary)
    call dealloc_pairlist_all   (pairlist)
    call dealloc_enefunc_all    (enefunc)
    call dealloc_domain_all     (domain)
    if ((genesis_run_mode == GenesisFEPMD) .or.    &
        (genesis_run_mode == GenesisFEPMIN) .or.   &
        (genesis_run_mode == GenesisFEPREMD)) then
        call dealloc_alchemy(alchemy)
    end if

    call timer(TimerTotal, TimerOff)


    ! output process time
    !
    call output_time

    return

  end subroutine domain_decomposition_genesis

end program spdyn
