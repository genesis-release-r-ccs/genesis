!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_control_mod
!> @brief   read parameters and data for MD trajectory analysis
!! @authors Motoshi Kamiya (MK)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_control_mod

  use ea_option_mod
  use at_input_mod
  use output_mod

  use at_energy_mod
  use at_constraints_mod
  use at_boundary_mod
  use at_ensemble_mod
  use at_dynamics_mod
  use at_remd_mod
  use at_qmmm_mod

  use select_mod
  use at_restraints_mod
  use fitting_mod
  use trajectory_mod
  use mpi_parallel_mod
  use string_mod
  use messages_mod
  use fileio_control_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  character(MaxLine) :: commandline_options = "t:"

  ! structure
  type, public :: s_ctrl_data

    type(s_inp_info)      :: inp_info
    type(s_out_info)      :: out_info
    type(s_ene_info)      :: ene_info
    type(s_cons_info)     :: cons_info
    type(s_boundary_info) :: bound_info
    type(s_sel_info)      :: sel_info
    type(s_res_info)      :: res_info
    type(s_rep_info)      :: rep_info
    type(s_dyn_info)      :: dyn_info
    type(s_ens_info)      :: ens_info
    ! analysis specific
    type(s_trj_info)      :: trj_info
    type(s_opt_info)      :: opt_info
    ! for QM/MM
    type(s_qmmm_info)     :: qmmm_info

  end type s_ctrl_data

  public  :: usage
  public  :: control
  private :: show_ctrl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    usage
  !> @brief        show usage of energy_analysis
  !! @authors      MK
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine usage(arg1)

    ! formal arguments
    character(*), intent(out) :: arg1

    ! local variables
    character(MaxLine)        :: arg2
    integer                   :: iargc

    call getarg(1,arg1)

    if ( iargc() < 1     .or. &
         arg1 == '-help' .or. &
         arg1 == '-h'    .or. &
         arg1 == '-HELP' .or. &
         arg1 == '-H' ) then

      if (main_rank) then

        call getarg(2, arg2)
        call tolower(arg2)

        if (iargc() < 2 .or. arg2 /= 'ctrl' .and. arg2 /= 'ctrl_all') then

          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# normal usage'
          write(MsgOut,'(A)') '  % ./energy_analysis INP'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters'
          write(MsgOut,'(A)') '  % ./energy_analysis -h ctrl'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters'
          write(MsgOut,'(A)') '  % ./energy_analysis -h ctrl_all'
          write(MsgOut,'(A)') ' '

          stop

        else

          call show_ctrl(arg2)

        end if

      end if

    else

      if (main_rank) then
        write(MsgOut,'(A)')'****************************************************'
        write(MsgOut,'(A)')'*                                                  *'
        write(MsgOut,'(A)')'*                 ENERGY_ANALYSIS                  *'
        write(MsgOut,'(A)')'*                                                  *'
        write(MsgOut,'(A)')'*                 Energy analyzer                  *'
        write(MsgOut,'(A)')'*                                                  *'
        write(MsgOut,'(A)')'*               Developed by RIKEN                 *'
        write(MsgOut,'(A)')'*                                                  *'
        write(MsgOut,'(A)')'****************************************************'
        write(MsgOut,'(A)')' '
      end if

    end if

  end subroutine usage

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    show_ctrl
  !> @brief        show usage of energy_analysis control file
  !! @authors      MK
  !! @param[in]    ctrl_str : ctrl_string ("ctrl" or "ctrl_all")
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl(ctrl_str)

    ! formal arguments
    character(*), intent(in) :: ctrl_str

    ! local variables
    logical :: show_all

    show_all = (ctrl_str == 'ctrl_all')

    call show_ctrl_input(show_all,'md')
    call show_ctrl_output('enefile')
    call show_ctrl_remd(show_all,'remd')
    call show_ctrl_energy(show_all,'md')
    call show_ctrl_constraints(show_all,'md')
    call show_ctrl_boundary(show_all,'md')
    call show_ctrl_selection
    call show_ctrl_restraints(show_all,'md')
    call show_ctrl_trajectory
    call show_ctrl_option
    call show_ctrl_qmmm(show_all, 'md')

    stop

  end subroutine show_ctrl

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    control
  !> @brief        open/read/close control files for MD
  !! @authors      MK
  !! @param[in]    filename  : file name of control file
  !! @param[in]    ctrl_data : control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control(filename, ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data

    ! local variables
    integer :: handle, i
    character(1) :: c

    ! open control file
    !
    if (len_trim(filename) > 0) then
      call open_ctrlfile(filename, handle)

      if (handle == 0) &
        call error_msg('Control> File IO Error')

      ! read input section
      !
      call read_ctrl_input(handle, ctrl_data%inp_info)

      ! read output section
      !
      if (main_rank) then
        call read_ctrl_output(handle, ctrl_data%out_info)
      end if

      ! read replica section
      !
      if ( find_ctrlfile_section(filename, 'REMD') ) then
        call read_ctrl_remd(handle, ctrl_data%rep_info)
      end if

      ! read dynamics section
      !
      call read_ctrl_dynamics(handle, ctrl_data%dyn_info)

      ! read energy section
      !
      call read_ctrl_energy(handle, ctrl_data%ene_info)

      ! read constraints section
      !
      call read_ctrl_constraints(handle, ctrl_data%cons_info)

      ! read ensemble section
      !
      call read_ctrl_ensemble(handle, ctrl_data%ens_info)

      ! read boundary section
      !
      call read_ctrl_boundary(handle, ctrl_data%bound_info)

      ! read selection section
      !
      call read_ctrl_selection(handle, ctrl_data%sel_info)

      ! read restraints section
      !
      call read_ctrl_restraints(handle, ctrl_data%res_info)

      ! read trajectory section
      !
      call read_ctrl_trajectory(handle, ctrl_data%trj_info)

      ! read option section
      !
      call read_ctrl_option(handle, ctrl_data%opt_info)

      ! read QM/MM section
      !
      call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)

      ! close control file
      !
      call close_ctrlfile(handle)
      
    end if

    return

  end subroutine control

end module ea_control_mod
