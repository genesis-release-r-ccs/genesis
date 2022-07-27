!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_setup_mod
!> @brief   setup variables and structures in eigmat analysis
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_setup_mod

  use ea_control_mod
  use output_mod
  use output_str_mod
  use input_mod
  use input_str_mod

  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data,  &
                   input,      &
                   output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module ea_setup_mod
