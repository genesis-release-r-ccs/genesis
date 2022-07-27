!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ma_main
!! @brief   MSD analysis
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program da_main

  use da_analyze_mod
  use da_setup_mod
  use da_control_mod
  use da_option_str_mod
  use input_str_mod
  use output_str_mod
  use string_mod
  use messages_mod

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_input)          :: input
  type(s_output)         :: output
  type(s_option)         :: option


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

  call setup(ctrl_data, input, output, option)


  ! [Step3] Analysis of mean square dispacement
  !
  write(MsgOut,'(A)') '[STEP3] Analysis of mean square dispacement'
  call analyze(input, output, option)


  ! [Step4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP5] Deallocate memory'
  write(MsgOut,'(A)') ' '

  call dealloc_option(option)

  stop

end program da_main
