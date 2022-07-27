!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pd_setup_mod
!> @brief   setup variables and structures in pcavec drawer
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pd_setup_mod

  use pd_control_mod
  use pd_option_mod
  use pd_option_str_mod
  use input_mod
  use input_str_mod
  use output_mod
  use output_str_mod
  use trajectory_mod
  use trajectory_str_mod
  use select_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_pdb_mod

  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data,  &
                   option,     &
                   input,      &
                   output)

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

end module pd_setup_mod
