!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_control_md
!> @brief   Reads parameters and data for MD simulation.
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!!          Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_control_mod

  use cg_output_mod
  use cg_input_mod
  use cg_minimize_mod
  use cg_dynamics_mod
  use cg_ensemble_mod
  use cg_restraints_mod
  use cg_boundary_mod
  use cg_energy_mod
  use select_mod
  use fitting_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef RICC
  use service_routines
#endif

  implicit none
  private

  ! structures
  type, public :: s_ctrl_data

    ! data for section input
    type(s_inp_info)      :: inp_info

    ! data for section output
    type(s_out_info)      :: out_info

    ! data for section energy
    type(s_ene_info)      :: ene_info

    ! data for section dynamics
    type(s_dyn_info)      :: dyn_info

    ! data for section minimize
    type(s_min_info)      :: min_info

    ! data for section ensemble
    type(s_ens_info)      :: ens_info

    ! data for section boundary
    type(s_boundary_info) :: bound_info

    ! data for section selection
    type(s_sel_info)      :: sel_info

    ! data for section restraints
    type(s_res_info)      :: res_info

    ! data for section fitting
    type(s_fit_info)      :: fit_info

    ! add more parameters here

  end type s_ctrl_data

  ! subroutines
  public  :: usage
  public  :: control_md
  public  :: control_min
  private :: show_ctrl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    usage
  !> @brief        Displays usage information for CGDYN.
  !! @authors      YS, TM
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine usage(arg1)

    ! formal arguments
    character(*),            intent(out) :: arg1

    ! local variables
    character(MaxLine)       :: arg2, arg3
    integer                  :: iargc


    call getarg(1,arg1)

    if (iargc() < 1      .or. &
         arg1 == '-help' .or. &
         arg1 == '-h'    .or. &
         arg1 == '-HELP' .or. &
         arg1 == '-H') then

      if (main_rank) then

        call getarg(2, arg2)
        call tolower(arg2)

        if (iargc() < 2 .or. arg2 /= 'ctrl' .and. arg2 /= 'ctrl_all') then
          write(MsgOut,'(A)') ' '
                  write(MsgOut,'(A)') '# Normal usage'
                  write(MsgOut,'(A)') '  % mpirun -np XX ./cgdyn INP'
                  write(MsgOut,'(A)') ' '
                  write(MsgOut,'(A)') '# Check control parameters of MD.'
                  write(MsgOut,'(A)') '  % ./cgdyn -h ctrl md'
                  write(MsgOut,'(A)') ' '
                  write(MsgOut,'(A)') '# Check control parameters of MIN.'
                  write(MsgOut,'(A)') '  % ./cgdyn -h ctrl min'
                  write(MsgOut,'(A)') ' '
                  write(MsgOut,'(A)') '# Check all control parameters of MD.'
                  write(MsgOut,'(A)') '  % ./cgdyn -h ctrl_all md'
                  write(MsgOut,'(A)') ' '
                  write(MsgOut,'(A)') '# Check all control parameters of MIN.'
                  write(MsgOut,'(A)') '  % ./cgdyn -h ctrl_all min'
                  write(MsgOut,'(A)') ' '
        else if (arg2 == 'ctrl' .or. arg2 == 'ctrl_all') then

          call getarg(3, arg3)
          call tolower(arg3)

          if (arg3 == '') &
            arg3 = 'md'

          call show_ctrl(arg2, arg3)

        end if

      end if

#ifdef HAVE_MPI_GENESIS
      call MPI_Finalize(ierror)
#endif

      stop

    else

      if (.not. main_rank) &
        return

      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*                  GENESIS CGDYN                   *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*        A Molecular Dynamics Simulator with       *'
      write(MsgOut,'(A)')'*           Spatial Decomposition Scheme           *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*               Developed by RIKEN                 *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')' '

    end if


    return

  end subroutine usage

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_md
  !> @brief        Opens, reads, and closes control files for MD.
  !! @authors      YS, TM, CK
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine control_md(filename, ctrl_data)
  
    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    
    ! local variables
    integer                  :: handle


    ! open control file
    ! 
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Error: Control_Md - File I/O Error.')



    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read dynamics section
    !
    call read_ctrl_dynamics(handle, ctrl_data%dyn_info)


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


    ! read fitting section
    !
    call read_ctrl_fitting_md(handle, ctrl_data%fit_info)


    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_min
  !> @brief        Opens, reads, and closes control files for minimization.
  !! @authors      TM, CK
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine control_min(filename, ctrl_data)
  
    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    
    ! local variables
    integer                  :: handle


    ! open control file
    ! 
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Error: Control_Min - File I/O Error.')


    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read minimize section
    !
    call read_ctrl_minimize(handle, ctrl_data%min_info)

    ! read boundary section
    !
    call read_ctrl_boundary(handle, ctrl_data%bound_info)


    ! read selection section
    !
    call read_ctrl_selection(handle, ctrl_data%sel_info)


    ! read restraints section
    !
    call read_ctrl_restraints(handle, ctrl_data%res_info)


    ! read fitting section
    !
    call read_ctrl_fitting_md(handle, ctrl_data%fit_info)


    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl
  !> @brief        Displays usage information for the GENESIS CGDYN control file.
  !! @authors      NT
  !! @param[in]    ctrl_str : ctrl string     : "ctrl", "ctrl_all"
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl(ctrl_str, run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_str
    character(*),            intent(in)    :: run_mode

    ! local variables
    logical                  :: show_all


    if (run_mode /= 'md'  .and. &
        run_mode /= 'min') then

      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# Check control parameters of MD.'
      write(MsgOut,'(A)') '  % ./cgdyn -h '//trim(ctrl_str)//' md'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# Check control parameters of MIN.'
      write(MsgOut,'(A)') '  % ./cgdyn -h '//trim(ctrl_str)//' min'
      write(MsgOut,'(A)') ' '
      return

    end if
      

    show_all = (ctrl_str == 'ctrl_all')


    ! show input section
    !
    call show_ctrl_input(show_all, run_mode)


    ! show output section
    !
    call show_ctrl_output(show_all, run_mode)


    ! show energy section
    !
    call show_ctrl_energy(show_all, run_mode)


    ! show dynamics section
    !
    call show_ctrl_dynamics(show_all, run_mode)


    ! show minimize section
    !
    call show_ctrl_minimize(show_all, run_mode)


    ! show ensemble section
    !
    call show_ctrl_ensemble(show_all, run_mode)


    ! show boundary section
    !
    call show_ctrl_boundary(show_all, run_mode)


    ! show restraints section
    !
    call show_ctrl_restraints(show_all, run_mode)


    return
  
  end subroutine show_ctrl

end module cg_control_mod
