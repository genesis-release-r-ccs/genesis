!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  pn_main
!! @brief   PATHCV analysis tool
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program pn_main

  use pn_analyze_mod
  use pn_setup_mod
  use pn_control_mod
  use pn_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use input_str_mod
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
  type(s_option)         :: option
  type(s_input)          :: input
  type(s_output)         :: output


  my_city_rank = 0
  nproc_city   = 1
  main_rank    = .true.


  ! show usage
  !
  call usage(ctrl_filename)


  ! [Step1] Read control file
  !
  write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
  write(MsgOut,'(A)') ' '

  call control(ctrl_filename, ctrl_data)


  ! [Step2] Set relevant variables and structures 
  !
  write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
  write(MsgOut,'(A)') ' '

  call setup(ctrl_data, option, input, output, molecule, trj_list, trajectory)
!  call setup(ctrl_data, option, input, output)


  ! [Step3] Analyze trajectory
  !
  write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
  write(MsgOut,'(A)') ' '

  if (option%trajectory) then
    call analyze(input, output, option,  molecule, trj_list, trajectory)
  else
    call analyze(input, output, option)
  endif


  ! [Step4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP4] Deallocate memory'
  write(MsgOut,'(A)') ' '

  stop

end program pn_main
