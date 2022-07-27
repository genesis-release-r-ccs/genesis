!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_control_mod
!> @brief   read parameters and data for md simulation
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK), 
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_control_mod

  use at_output_mod
  use at_input_mod
  use at_morph_mod
  use at_remd_mod
  use at_rpath_mod
  use at_rpath_str_mod
  use at_minimize_mod
  use at_vibration_mod
  use at_dynamics_mod
  use at_ensemble_mod
  use at_restraints_mod
  use at_constraints_mod
  use at_boundary_mod
  use at_energy_mod
  use at_qmmm_mod
  use at_experiments_mod
  use at_gamd_mod
  use select_mod
  use fitting_mod
  use fileio_control_mod
  use string_mod
  use fitting_mod
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

    ! data for section QM/MM
    type(s_qmmm_info)      :: qmmm_info

    ! data for section restraints
    type(s_res_info)      :: res_info

    ! data for section restraints
    type(s_exp_info)      :: exp_info

    ! data for section fitting
    type(s_fit_info)      :: fit_info

    ! data for section vibration
    type(s_vib_info)      :: vib_info

    ! data for section morphing
    type(s_morph_info)    :: morph_info

    ! data for section gamd
    type(s_gamd_info)     :: gamd_info

    ! add more parameters here

  end type s_ctrl_data

  ! subroutines
  public  :: usage
  public  :: control_md
  public  :: control_min
  public  :: control_remd
  public  :: control_rpath
  public  :: control_vib
  public  :: control_morph
  private :: show_ctrl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    usage
  !> @brief        show usage of ATDYN
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

    if (iargc() < 1        .or. &
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
          write(MsgOut,'(A)') '  % ./atdyn INP'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of md'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl md'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of min'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl min'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of remd'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl remd'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of rpath'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl rpath'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check control parameters of vibration'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl vib'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of md'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl_all md'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of min'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl_all min'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of remd'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl_all remd'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of rpath'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl_all rpath'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# check all control parameters of vibration'
          write(MsgOut,'(A)') '  % ./atdyn -h ctrl_all vib'
          write(MsgOut,'(A)') ' '

        else if (arg2 .eq. 'ctrl' .or. arg2 .eq. 'ctrl_all') then

          call getarg(3, arg3)
          call tolower(arg3)

          if (arg3 .eq. '') &
            arg3 = 'md'

          call show_ctrl(arg2, arg3)

        end if

      end if

#ifdef HAVE_MPI_GENESIS
#ifdef QSIMULATE
      call BLACS_Exit(1)
#endif
      call MPI_Finalize(ierror)
#endif

      stop

    else

      if (.not. main_rank) &
        return

      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*                  GENESIS ATDYN                   *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*        A Molecular Dynamics Simulator with       *'
      write(MsgOut,'(A)')'*            Atomic Decomposition Scheme           *'
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


    ! read QM/MM section
    !
    call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)


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


    ! read QM/MM section
    !
    call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)


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


    ! read gamd section
    !
    call read_ctrl_gamd(handle, ctrl_data%gamd_info)


    ! read replica section
    !
    call read_ctrl_remd(handle, ctrl_data%rep_info)


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


    ! read QM/MM section
    !
    call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)


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
  !! @authors      TM, CK, KY
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


    if ((ctrl_data%rpath_info%rpathmode == RpathmodeMFEP)  .or. &
        (ctrl_data%rpath_info%rpathmode == RpathmodeMEPMD) .or. &
        (ctrl_data%rpath_info%rpathmode == RpathmodeFEP)) then

      ! read dynamics section
      !
      call read_ctrl_dynamics(handle, ctrl_data%dyn_info)

      ! read ensemble section
      !
      call read_ctrl_ensemble(handle, ctrl_data%ens_info)

    else if (ctrl_data%rpath_info%rpathmode == RpathmodeMEP) then

      ! read minimize section
      !
      call read_ctrl_minimize(handle, ctrl_data%min_info)

    end if


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


    ! read QM/MM section
    !
    call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)


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
  !  Subroutine    control_vib
  !> @brief        open/read/close control files for vibration
  !! @authors      KY
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine control_vib(filename, ctrl_data)
  
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


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read minimize section
    !
    call read_ctrl_vibration(handle, ctrl_data%vib_info)


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
    call read_ctrl_fitting(handle, ctrl_data%fit_info)

    ! read QM/MM section
    !
    call read_ctrl_qmmm(handle, filename, ctrl_data%qmmm_info)


    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_morph
  !> @brief        open/read/close control files for Minimization
  !! @authors      CK
  !! @param[in]    filename  : file name of control file
  !! @param[out]   ctrl_data : information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine control_morph(filename, ctrl_data)
  
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


    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)


    ! read morphing section
    !
    call read_ctrl_morph(handle, ctrl_data%morph_info)


    ! read boundary section
    !
    call read_ctrl_boundary(handle, ctrl_data%bound_info)


    ! read selection section
    !
    call read_ctrl_selection(handle, ctrl_data%sel_info)


    ! read restraints section
    !
    call read_ctrl_restraints(handle, ctrl_data%res_info)


    ! read restraints section
    !
    call read_ctrl_fitting(handle, ctrl_data%fit_info)

    ! close control file
    !
    call close_ctrlfile(handle)


    return

  end subroutine control_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl
  !> @brief        show usage of GENESIS ATDYN control file
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


    if (run_mode .ne. 'md'    .and. &
        run_mode .ne. 'min'   .and. &
        run_mode .ne. 'remd'  .and. &
        run_mode .ne. 'rpath' .and. &
        run_mode .ne. 'vib' ) then

      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of md'
      write(MsgOut,'(A)') '  % ./atdyn -h '//trim(ctrl_str)//' md'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of min'
      write(MsgOut,'(A)') '  % ./atdyn -h '//trim(ctrl_str)//' min'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of remd'
      write(MsgOut,'(A)') '  % ./atdyn -h '//trim(ctrl_str)//' remd'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of rpath'
      write(MsgOut,'(A)') '  % ./atdyn -h '//trim(ctrl_str)//' rpath'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') '# check control parameters of vib'
      write(MsgOut,'(A)') '  % ./atdyn -h '//trim(ctrl_str)//' vib'
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


    ! show gamd section
    !
    call show_ctrl_gamd(show_all, run_mode)


    ! show remd section
    !
    call show_ctrl_remd(show_all, run_mode)


    ! show rpath section
    !
    call show_ctrl_rpath(show_all, run_mode)


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


    ! show QM/MM section
    !
    call show_ctrl_qmmm(show_all, run_mode)


    ! show vibration section
    !
    call show_ctrl_vibration(show_all, run_mode)

    ! show experiments section
    !
    call show_ctrl_experiments(show_all, run_mode)


    return
  
  end subroutine show_ctrl

end module at_control_mod
