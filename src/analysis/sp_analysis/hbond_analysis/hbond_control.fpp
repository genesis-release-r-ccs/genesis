!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hbond_control_mod
!> @brief   read parameters and data for MD trajectory analysis
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hbond_control_mod

  use hbond_option_mod
  use sa_boundary_mod
  use sa_ensemble_mod
  use sa_option_mod, only : show_ctrl_option
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
  type, public :: s_hbond_ctrl_data

    ! data for section option
    type(s_hbond_opt_info) :: hbond_opt_info

  end type s_hbond_ctrl_data

  ! subroutines
  public  :: hbond_usage
  public  :: hbond_control

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hbond_usage
  !> @brief        show hbond_usage of trj_analysis
  !! @authors      TM
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hbond_usage(arg1)

    ! formal arguments
    character(*),            intent(out)   :: arg1

    ! local variables
    character(Maxline)       :: arg2, arg3
    integer                  :: iargc


    call getarg(1,arg1)

    if (iargc() < 1        .or. &
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
          write(MsgOut,'(A)') '  % ./hbond_analysis INP'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# usage to see control parameters'
          write(MsgOut,'(A)') '  % ./hbond_analysis -h ctrl'
          write(MsgOut,'(A)') ' '
        end if
      else

        if (main_rank) then
          write(MsgOut,'(A)') '# control parameters in hbond_analysis'
          write(MsgOut,'(A)') ' '

          call show_ctrl_input ('psf,ref,armtop,ambcrd,grotop,grocrd')
          call show_ctrl_output('txt, HBLIST')
          call show_ctrl_trajectory
          call show_ctrl_ensemble
          call show_ctrl_boundary
          call show_ctrl_selection
          call show_ctrl_option
          call show_hbond_ctrl_option
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
       write(MsgOut,'(A)')'*                  H-bond analysis                 *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'*              Developed by RIKEN TMS              *'
       write(MsgOut,'(A)')'*                                                  *'
       write(MsgOut,'(A)')'****************************************************'
       write(MsgOut,'(A)')' '
      end if
    end if

    return

  end subroutine hbond_usage

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    hbond_control
  !> @brief        open/read/close control files
  !! @authors      TM
  !! @param[in]    filename  : file name of control file
  !! @param[inout] ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hbond_control(filename, hbond_ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_hbond_ctrl_data), intent(inout) :: hbond_ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    !
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control> File IO Error')

    ! read option section
    !
    call read_ctrl_option (handle, hbond_ctrl_data%hbond_opt_info)

    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine hbond_control

end module hbond_control_mod
