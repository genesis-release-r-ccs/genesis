!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_setup_mod
!> @brief   setup variables and structures
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_setup_mod

  use da_control_mod
  use da_option_str_mod
  use da_option_mod
  use output_mod
  use input_mod
  use input_str_mod
  use output_str_mod

  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, input, output, option)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! local variables

    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, option, output)

    return

  end subroutine setup

end module da_setup_mod
