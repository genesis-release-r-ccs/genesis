!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_control_md
!> @brief   read parameters and data for md simulation
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!!          Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_control_mod

  use sp_output_mod
  use sp_input_mod
  use sp_minimize_mod
  use sp_dynamics_mod
  use sp_ensemble_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_boundary_mod
  use sp_energy_mod
  use sp_remd_mod
  use sp_rpath_mod
  use sp_experiments_mod
  use sp_gamd_mod
  use select_mod
  use fitting_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_mod
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

    ! data for section replica
    type(s_rep_info)      :: rep_info

    ! data for section rpath
    type(s_rpath_info)    :: rpath_info

    ! data for section energy
    type(s_ene_info)      :: ene_info

    ! data for section dynamics
    type(s_dyn_info)      :: dyn_info

    ! data for section minimize
    type(s_min_info)      :: min_info

    ! data for section constraints
    type(s_cons_info)     :: cons_info

    ! data for section ensemble
    type(s_ens_info)      :: ens_info

    ! data for section boundary
    type(s_boundary_info) :: bound_info

    ! data for section selection
    type(s_sel_info)      :: sel_info

    ! data for section restraints
    type(s_res_info)      :: res_info

    ! data for section restraints
    type(s_exp_info)      :: exp_info

    ! data for section fitting
    type(s_fit_info)      :: fit_info

    ! data for section gamd
    type(s_gamd_info)     :: gamd_info

    ! data for section alchemy
    type(s_alch_info)     :: alch_info

    ! add more parameters here

  end type s_ctrl_data

  ! subroutines
  public  :: usage
  public  :: control_md
  public  :: control_min
  public  :: control_remd
  public  :: control_rpath
  private :: show_ctrl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    usage
  !> @brief        show usage of SPDYN
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
         arg1 .eq. '-help' .or. &
         arg1 .eq. '-h'    .or. &
         arg1 .eq. '-HELP' .or. &
         arg1 .eq. '-H') then

      if (main_rank) then

        call getarg(2, arg2)
        call tolower(arg2)

        if (iargc() < 2 .or. arg2 .ne. 'ctrl' .and. arg2 .ne. 'ctrl_all') then
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# normal usage'
          write(MsgOut,'(A)') '  % mpirun -np XX ./spdyn INP'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of md'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl md'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of min'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl min'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of remd'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl remd'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of rpath'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl rpath'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of md'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl_all md'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of min'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl_all min'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of remd'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl_all remd'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of rpath'
          write(MsgOut,'(A)') '  % ./spdyn -h ctrl_all rpath'
          write(MsgOut,'(A)') ' '

        else if (arg2 .eq. 'ctrl' .or. arg2 .eq. 'ctrl_all') then

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
      write(MsgOut,'(A)')'*                  GENESIS SPDYN                   *'
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
  !> @brief        open/read/close control files for MD
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
      call error_msg('Control_Md> File IO Error')


    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read gamd section
    !
    call read_ctrl_gamd(handle, ctrl_data%gamd_info)


    ! read alchemy section
    !
    if (find_ctrlfile_section(filename, 'ALCHEMY')) then
      call read_ctrl_alchemy(handle, ctrl_data%alch_info)
    end if


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read dynamics section
    !
    call read_ctrl_dynamics(handle, ctrl_data%dyn_info)


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


    ! read fitting section
    !
    call read_ctrl_fitting_md(handle, ctrl_data%fit_info)


    ! read experiments section
    !
    call read_ctrl_experiments(handle, ctrl_data%exp_info)

    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_min
  !> @brief        open/read/close control files for Minimization
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
      call error_msg('Control_Min> File IO Error')


    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read alchemy section
    !
    if (find_ctrlfile_section(filename, 'ALCHEMY')) then
      call read_ctrl_alchemy(handle, ctrl_data%alch_info)
    end if


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read minimize section
    !
    call read_ctrl_minimize(handle, ctrl_data%min_info)

    ! read constraints section
    !
    call read_ctrl_constraints(handle, ctrl_data%cons_info)

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


    ! read experiments section
    !
    call read_ctrl_experiments(handle, ctrl_data%exp_info)


    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_remd
  !> @brief        open/read/close control files for REMD
  !! @authors      TM, CK
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_remd(filename, ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    ! 
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control_Remd> File IO Error')


    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read replica section
    !
    call read_ctrl_remd(handle, ctrl_data%rep_info)


    ! read gamd section
    !
    call read_ctrl_gamd(handle, ctrl_data%gamd_info)


    ! read alchemy section
    !
    if (find_ctrlfile_section(filename, 'ALCHEMY')) then
      call read_ctrl_alchemy(handle, ctrl_data%alch_info)
    end if


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read dynamics section
    !
    call read_ctrl_dynamics(handle, ctrl_data%dyn_info)


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

    ! read fitting section
    !
    call read_ctrl_fitting_md(handle, ctrl_data%fit_info)


    ! read experiments section
    !
    call read_ctrl_experiments(handle, ctrl_data%exp_info)


    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine control_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_rpath
  !> @brief        open/read/close control files for RPATH
  !! @authors      TM, CK
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_rpath(filename, ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    ! 
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control_Rpath> File IO Error')


    ! read input section
    !
    call read_ctrl_input(handle, ctrl_data%inp_info)


    ! read output section
    !
    call read_ctrl_output(handle, ctrl_data%out_info)


    ! read rpath section
    !
    call read_ctrl_rpath(handle, ctrl_data%rpath_info)


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read dynamics section
    !
    call read_ctrl_dynamics(handle, ctrl_data%dyn_info)


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


    ! read fitting section
    !
    call read_ctrl_fitting_md(handle, ctrl_data%fit_info)


    ! read experiments section
    !
    call read_ctrl_experiments(handle, ctrl_data%exp_info)


    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine control_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl
  !> @brief        show usage of GENESIS SPDYN control file
  !! @authors      NT
  !! @param[in]    ctrl_str : ctrl string     : "ctrl", "ctrl_all"
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl(ctrl_str, run_mode)

    ! formal arguments
    character(*),            intent(in)    :: ctrl_str
    character(*),            intent(in)    :: run_mode

    ! local variables
    logical                  :: show_all


    if (run_mode .ne. 'md'  .and. &
        run_mode .ne. 'min' .and. &
        run_mode .ne. 'remd' .and. &
        run_mode .ne. 'rpath') then

      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of md'
      write(MsgOut,'(A)') '  % ./spdyn -h '//trim(ctrl_str)//' md'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of min'
      write(MsgOut,'(A)') '  % ./spdyn -h '//trim(ctrl_str)//' min'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of remd'
      write(MsgOut,'(A)') '  % ./spdyn -h '//trim(ctrl_str)//' remd'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of rpath'
      write(MsgOut,'(A)') '  % ./spdyn -h '//trim(ctrl_str)//' rpath'
      write(MsgOut,'(A)') ' '
      return

    end if
      

    show_all = (ctrl_str .eq. 'ctrl_all')


    ! show input section
    !
    call show_ctrl_input(show_all, run_mode)


    ! show output section
    !
    call show_ctrl_output(show_all, run_mode)


    ! show remd section
    !
    call show_ctrl_remd(show_all, run_mode)


    ! show rpath section
    !
    call show_ctrl_rpath(show_all, run_mode)


    ! show gamd section
    !
    call show_ctrl_gamd(show_all, run_mode)


    ! show alchemy section
    !
    call show_ctrl_alchemy(show_all, run_mode)


    ! show energy section
    !
    call show_ctrl_energy(show_all, run_mode)


    ! show dynamics section
    !
    call show_ctrl_dynamics(show_all, run_mode)


    ! show minimize section
    !
    call show_ctrl_minimize(show_all, run_mode)


    ! show constraints section
    !
    call show_ctrl_constraints(show_all, run_mode)


    ! show ensemble section
    !
    call show_ctrl_ensemble(show_all, run_mode)


    ! show boundary section
    !
    call show_ctrl_boundary(show_all, run_mode)


    ! show restraints section
    !
    call show_ctrl_restraints(show_all, run_mode)

    ! show experiments section
    !
    call show_ctrl_experiments(show_all, run_mode)

    return
  
  end subroutine show_ctrl

end module sp_control_mod
