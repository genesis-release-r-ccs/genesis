!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  rg_main
!! @brief   RG analysis
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program rg_main

  use rg_analyze_mod
  use rg_setup_mod
  use rg_control_mod
  use rg_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_molecule)       :: molecule
  type(s_trj_list)       :: trj_list
  type(s_trajectory)     :: trajectory
  type(s_output)         :: output
  type(s_option)         :: option


  my_city_rank = 0
  nproc_city   = 1
  main_rank    = .true.

  ! show usage
  !
  call usage(ctrl_filename)

  ! [STEP1] Read control file
  !
  write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
  write(MsgOut,'(A)') ' '

  call control(ctrl_filename, ctrl_data)

  ! [STEP2] Set variables and structures
  !
  write(MsgOut,'(A)') '[STEP2] Set Variables and Structures'
  write(MsgOut,'(A)') ' '

  call setup(ctrl_data, molecule, trj_list, trajectory, output, option)

  ! [STEP3] Analyze trajectory
  !
  write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
  write(MsgOut,'(A)') ' '

  call analyze(molecule, trj_list, output, option, trajectory)

  ! [STEP4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP4] Deallocate memory'
  write(MsgOut,'(A)') ' '

  call dealloc_trajectory(trajectory)
  call dealloc_option(option)
  call dealloc_trj_list(trj_list)
  call dealloc_molecules_all(molecule)

  stop

end program rg_main
