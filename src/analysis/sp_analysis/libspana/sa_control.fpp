!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_control_mod
!> @brief   read parameters and data for MD trajectory analysis
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_control_mod

  use sa_boundary_mod
  use sa_ensemble_mod
  use sa_option_mod
  use trajectory_mod
  use fitting_mod
  use output_mod
  use input_mod
  use fileio_control_mod
  use select_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_ctrl_data

    ! data for section input
    type(s_inp_info)      :: inp_info

    ! data for section output
    type(s_out_info)      :: out_info

    ! data for section trajectory
    type(s_trj_info)      :: trj_info

    ! data for section selection
    type(s_sel_info)      :: sel_info

    ! data for section selection
    type(s_fit_info)      :: fit_info

    ! data for section option
    type(s_opt_info)      :: opt_info

    ! data for section ensemble
    type(s_ens_info)      :: ens_info

    ! data for section boundary
    type(s_boundary_info) :: bound_info


  end type s_ctrl_data

  ! subroutines
  public  :: control

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control
  !> @brief        open/read/close control files
  !! @authors      TM
  !! @param[in]    filename  : file name of control file
  !! @param[inout] ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control(filename, ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    !
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control> File IO Error')


    ! read input section
    ! these functions are defined in /libana/input.fpp but not
    ! for MPI. All ranks output the values.
    call read_ctrl_input(handle, ctrl_data%inp_info)

    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)

    ! read trajectory section
    !
    call read_ctrl_trajectory(handle, ctrl_data%trj_info)

    ! read selection section
    !
    call read_ctrl_selection(handle, ctrl_data%sel_info)

    ! read fitting section
    !
    call read_ctrl_fitting(handle, ctrl_data%fit_info)

    ! read option section
    !
    call read_ctrl_option(handle, ctrl_data%opt_info)

    ! read ensemble section
    !
    call read_ctrl_ensemble(handle, ctrl_data%ens_info)

    ! read boundary section
    !
    call read_ctrl_boundary(handle, ctrl_data%bound_info)

    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine control

end module sa_control_mod
