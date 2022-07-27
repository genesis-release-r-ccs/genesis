!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  pd_main
!! @brief   analysis trajectory files
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program pd_main

  use pd_draw_mod
  use pd_setup_mod
  use pd_control_mod
  use pd_option_str_mod
  use input_str_mod
  use output_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_input)          :: input
  type(s_output)         :: output
  type(s_option)         :: option


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

  call setup(ctrl_data,  &
             option,     &
             input,      &
             output)


  ! [Step3] Analyze trajectory
  !
  write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
  write(MsgOut,'(A)') ' '

  call draw(input, option, output)


  ! [Step4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP4] Deallocate memory'
  write(MsgOut,'(A)') ' '

  stop

end program pd_main
