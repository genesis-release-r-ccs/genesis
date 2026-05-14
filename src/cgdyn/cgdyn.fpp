!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  cgdyn
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

program cgdyn

  use cg_setup_mpi_mod
  use cg_md_vverlet_mod
  use cg_md_leapfrog_mod
  use cg_dynamics_mod
  use cg_minimize_mod
  use cg_setup_cgdyn_mod
  use cg_control_mod
  use cg_energy_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
  use cg_minimize_str_mod
  use cg_dynamics_str_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_output_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use fileio_control_mod
  use hardwareinfo_mod
  use timers_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_max_threads
#endif

  implicit none

  integer             :: genesis_run_mode
  integer, parameter  :: GenesisMD    = 1
  integer, parameter  :: GenesisMIN   = 2

  ! local variables
  character(MaxFilename)      :: ctrl_filename
  type(s_ctrl_data)           :: ctrl_data
  type(s_molecule)            :: molecule
  type(s_grotop)              :: grotop
  type(s_enefunc)             :: enefunc
  type(s_dynvars)             :: dynvars
  type(s_pairlist)            :: pairlist
  type(s_boundary)            :: boundary
  type(s_ensemble)            :: ensemble
  type(s_dynamics)            :: dynamics
  type(s_minimize)            :: minimize
  type(s_output)              :: output
  type(s_domain)              :: domain
  type(s_comm)                :: comm

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
  !! @param[out]   genesis_run_mode : run MD, MIN
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_genesis_mode(ctrl_filename, genesis_run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_filename
    integer,                 intent(inout) :: genesis_run_mode


    if (find_ctrlfile_section(ctrl_filename, 'DYNAMICS')) then
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
  !! @param[in]    genesis_run_mod : run MD, MIN
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

    case (GenesisMD)

      call control_md  (ctrl_filename, ctrl_data)

    case (GenesisMIN)

      call control_min (ctrl_filename, ctrl_data)

    end select


    ! [Step2] Setup MPI
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP2] Setup MPI'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD,GenesisMIN)

      call setup_mpi_md  (ctrl_data%ene_info)

    end select


    ! [Step3] Set relevant variables and structures 
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
      write(MsgOut,'(A)') ' '
    end if

    select case (genesis_run_mode)

    case (GenesisMD)

      call setup_cgdyn_md(ctrl_data, output, grotop, molecule, enefunc, &
                          pairlist, dynvars, dynamics, ensemble,        &
                          boundary, domain, comm)

    case (GenesisMIN)

      call setup_cgdyn_min(ctrl_data, output, grotop, molecule, enefunc, &
                          pairlist, dynvars, minimize, boundary, domain, &
                          comm)

    end select


    ! [Step4] Compute single point energy for molecules
    !
    if (main_rank) then
      write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
      write(MsgOut,'(A)') ' '
    end if

    call compute_energy(pairlist, boundary,      &
                        .false., .true., .true., &
                        enefunc, domain, dynvars)

    call output_energy(dynvars%step, enefunc, dynvars%energy)


    ! [Step5] Perform MD simulation or Energy minimization
    !
    call timer(TimerDynamics, TimerOn)

    select case (genesis_run_mode)

    case (GenesisMD)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Molecular Dynamics Simulation'
        write(MsgOut,'(A)') ' '
      end if

      call run_md(output, grotop, molecule, domain, enefunc, dynvars, &
                  dynamics, pairlist, boundary, ensemble, comm)
     
    case (GenesisMIN)

      if (main_rank) then
        write(MsgOut,'(A)') '[STEP5] Perform Energy Minimization'
        write(MsgOut,'(A)') ' '
      end if

      call run_min(output, domain, enefunc, dynvars, minimize,       &
                     pairlist, boundary, comm)

    end select

    call timer(TimerDynamics, TimerOff)


    ! [Step6] Deallocate arrays
    !
    if (main_rank) then
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
      write(MsgOut,'(A)') ' '
    end if

    call dealloc_boundary_all   (boundary)
    call dealloc_pairlist_all   (pairlist)
    call dealloc_enefunc_all    (enefunc)
    call dealloc_domain_all     (domain)

    call timer(TimerTotal, TimerOff)


    ! output process time
    !
    call output_time


    return

  end subroutine domain_decomposition_genesis

end program cgdyn
