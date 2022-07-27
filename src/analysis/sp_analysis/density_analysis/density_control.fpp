!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   density_control_mod
!> @brief   read parameters and data for MD trajectory analysis
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module density_control_mod

  use density_option_mod
  use sa_boundary_mod
  use sa_ensemble_mod
  use sa_option_mod, only : show_ctrl_option
  use fitting_mod
  use trajectory_mod
  use output_mod
  use input_mod
  use select_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_density_ctrl_data

    ! data for fitting section
    type(s_fit_info) :: fit_info

    ! data for section option
    type(s_density_opt_info) :: density_opt_info

  end type s_density_ctrl_data

  ! subroutines
  public  :: density_usage
  public  :: density_control

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    density_usage
  !> @brief        show density_usage of trj_analysis
  !! @authors      TM
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine density_usage(arg1)

    ! formal arguments
    character(*),            intent(out)   :: arg1

    ! local variables
    character(Maxline)       :: arg2, arg3
    integer                  :: iargc


    call getarg(1,arg1)

    if (iargc() < 1       .or. &
        arg1 .eq. '-help' .or. &
        arg1 .eq. '-h'    .or. &
        arg1 .eq. '-HELP' .or. &
        arg1 .eq. '-H') then

      call getarg(2, arg2)
      call tolower(arg2)

      if (iargc() < 2 .or. arg2 .ne. 'ctrl') then

        if (main_rank) then
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# normal usage'
          write(MsgOut,'(A)') '  % ./density_analysis INP'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# usage to see control parameters'
          write(MsgOut,'(A)') '  % ./density_analysis -h ctrl'
          write(MsgOut,'(A)') ' '
        end if
      else

        if (main_rank) then
           write(MsgOut,'(A)') '# control parameters in density_analysis'
           write(MsgOut,'(A)') ' '

          call show_ctrl_input ('psf,ref,pdb,top,armtop,ambref,ambcrd,grotop,groref,grocrd')
          call show_ctrl_output('map,pdb')
          call show_ctrl_trajectory
          call show_ctrl_ensemble
          call show_ctrl_boundary
          call show_ctrl_selection
          call show_ctrl_fitting
          call show_ctrl_option
          call show_density_ctrl_option
        end if
      end if

      stop

    else

      if (main_rank) then
       write(MsgOut,'(A)')'****************************************************'
       write(MsgOut,'(A)')'*            /==\_               _/==\             *'
       write(MsgOut,'(A)')'*                 ===============                  *'
       write(MsgOut,'(A)')'*                _===============_                 *'
       write(MsgOut,'(A)')'*            \==/      SPANA      \==/             *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'*            SPatial decomposition ANAlysis        *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'*                 DENSITY analysis                 *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'*              Developed by RIKEN TMS              *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'****************************************************'
       write(MsgOut,'(A)')' '
      end if
    end if

    return

  end subroutine density_usage

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    density_control
  !> @brief        open/read/close control files
  !! @authors      TM
  !! @param[in]    filename          : file name of control file
  !! @param[inout] density_ctrl_data : density information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine density_control(filename, density_ctrl_data)

    ! formal arguments
    character(*),              intent(in)    :: filename
    type(s_density_ctrl_data), intent(inout) :: density_ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    !
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control> File IO Error')

    ! read fitting section
    !
    call read_ctrl_fitting(handle, density_ctrl_data%fit_info)

    ! read option section
    !
    call read_ctrl_option (handle, density_ctrl_data%density_opt_info)

    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine density_control

end module density_control_mod
