!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  rs_main
!! @brief   convert restart file
!! @authors Norio Takase (NT), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program rs_main

  use rs_convert_mod
  use rs_setup_mod
  use rs_control_mod
  use rs_option_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_input)          :: input
  type(s_molecule)       :: molecule
  type(s_option)         :: option
  type(s_output)         :: output


  my_city_rank = 0
  nproc_city   = 1
  main_rank    = .true.


  ! show usage
  !
  call usage(ctrl_filename)


  ! [Step1] Read control file
  !
  write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Convert'
  write(MsgOut,'(A)') ' '

  call control(ctrl_filename, ctrl_data)


  ! [Step2] Set relevant variables and structures 
  !
  write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
  write(MsgOut,'(A)') ' '

  call setup(ctrl_data, input, molecule, option, output)


  ! [Step3] Convert restart file
  !
  write(MsgOut,'(A)') '[STEP3] Convert restart file'
  write(MsgOut,'(A)') ' '

  call convert(input, molecule, option, output)


  ! [Step4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP4] Deallocate memory'
  write(MsgOut,'(A)') ' '


  stop

end program rs_main
