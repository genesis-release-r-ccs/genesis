!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pm_setup_mod
!> @brief   setup variables and structures in PMF_ANALYSIS
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pm_setup_mod

  use pm_control_mod
  use pm_option_mod
  use pm_option_str_mod
  use output_mod
  use input_mod
  use output_str_mod
  use input_str_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in PMF_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, option, input, output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, option)


    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)


    return

  end subroutine setup

end module pm_setup_mod
