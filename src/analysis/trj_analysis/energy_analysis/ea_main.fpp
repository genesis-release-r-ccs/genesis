!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ea_main
!! @brief   energy analysis (supporting REST)
!! @authors Motoshi Kamiya (MK)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program ea_main

  use ea_control_mod
  use ea_option_str_mod
  use ea_option_mod
  use ea_setup_mod
  use ea_analyze_mod

  use at_enefunc_str_mod
  use at_dynvars_str_mod
  use at_pairlist_str_mod
  use at_boundary_str_mod
  use at_constraints_str_mod
  use at_energy_str_mod
  use at_remd_str_mod
  use at_dynamics_str_mod
  use at_ensemble_str_mod

  use trajectory_str_mod
  use output_str_mod
  use fileio_control_mod
  use hardwareinfo_mod
  use timers_mod
  use molecules_str_mod
  use string_mod
  use mpi_parallel_mod
  use messages_mod

#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_molecule)       :: molecule
  type(s_trj_list)       :: trj_list
  type(s_trajectory)     :: trajectory
  type(s_output)         :: output
  type(s_option)         :: option

  type(s_enefunc)        :: enefunc
  type(s_dynvars)        :: dynvars
  type(s_pairlist)       :: pairlist
  type(s_boundary)       :: boundary
  type(s_constraints)    :: constraints
  type(s_remd)           :: remd
  type(s_ensemble)       :: ensemble
  type(s_dynamics)       :: dynamics
  integer                :: omp_get_max_threads

#ifdef HAVE_MPI_GENESIS
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, my_world_rank, ierror)
  call mpi_comm_size(mpi_comm_world, nproc_world,   ierror)
  main_rank = (my_world_rank .eq. 0)
#else
  my_world_rank = 0
  nproc_world   = 1
  main_rank     = .true.
#endif

  ! openmp
  nthread = 1
#ifdef OMP
  nthread = omp_get_max_threads()
#endif

  ! show usage
  !
  call usage(ctrl_filename)

  ! [Step1] Read control file
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '
  end if

  call control(ctrl_filename, ctrl_data)

  ! [Step2] Set relevant variables and structures
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '
  end if

  call setup(ctrl_data, output, molecule, enefunc, pairlist, &
             dynvars, boundary, remd, ensemble, dynamics, &
             trj_list, trajectory, option)

  ! [Step3] Analyze trajectory
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
  end if

  call analyze(output, molecule, enefunc, pairlist, dynvars, boundary, &
               ensemble, remd, trj_list, trajectory, option)

  ! [Step4] Deallocate memory
  !
  if (main_rank) then
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '
  end if

  call dealloc_constraints_all(constraints)
  call dealloc_boundary_all(boundary)
  call dealloc_pairlist_all(pairlist)
  call dealloc_energy_all(dynvars%energy)
  call dealloc_dynvars_all(dynvars)
  call dealloc_enefunc_all(enefunc)
  call dealloc_molecules_all(molecule)

#ifdef HAVE_MPI_GENESIS
  call mpi_finalize(ierror)
#endif


  stop

end program ea_main
