!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_output_mod
!> @brief   definition of output
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK), Jaewoon Jung (JJ),
!!          Ryuhei Harada (RH), Norio Takase (NT), Kiyoshi Yagi (KY),
!!          Yoshinobu Akinaga (YA)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_output_mod

  use at_dynvars_mod
  use at_boundary_mod
  use at_energy_mod
  use at_qmmm_mod
  use at_output_str_mod
  use at_minimize_str_mod
  use at_vibration_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_remd_str_mod
  use at_rpath_str_mod
  use at_morph_str_mod
  use at_boundary_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use fileio_rstmep_mod
  use fileio_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_out_info
    character(MaxFilename) :: logfile    = ''
    character(MaxFilename) :: dcdfile    = ''
    character(MaxFilename) :: dcdvelfile = ''
    character(MaxFilename) :: rstfile    = ''
    character(MaxFilename) :: pdbfile    = ''
    character(MaxFilename) :: remfile    = ''
    character(MaxFilename) :: rpathfile  = ''
    character(MaxFilename) :: rstmepfile = ''
    character(MaxFilename) :: minfofile  = ''
    character(MaxFilename) :: gamdfile   = ''
    character(MaxFilename) :: rpathlogfile  = ''
  end type s_out_info

  ! subroutines
  public  :: show_ctrl_output
  public  :: read_ctrl_output
  public  :: setup_output_md
  public  :: setup_output_min
  public  :: setup_output_remd
  public  :: setup_output_rpath
  public  :: setup_output_vib
  public  :: setup_output_bd
  public  :: setup_output_morph
  public  :: open_output
  public  :: close_output
  public  :: output_md
  public  :: output_min
  public  :: output_remd
  public  :: output_rpath
  public  :: output_vib
  public  :: output_bd
  public  :: output_morph
  public  :: output_gamd

  private :: output_restart_md
  private :: output_restart_min
  private :: output_restart_remd
  private :: output_restart_rpath
  private :: output_restart_bd
  private :: output_restart_morph
  private :: output_restart_pdb
  private :: write_restart_pdb
  private :: write_trajectory_dcd
  private :: write_trajectory_dcdvel
  private :: include_id_to_filename

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_output
  !> @brief        show OUTPUT section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_output(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') '# dcdfile    = sample.dcd   # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample.dcd   # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample.rst   # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample.pdb   # PDB file'
        write(MsgOut,'(A)') '# gamdfile   = sample.gamd  # gamd file'
        write(MsgOut,'(A)') ' '

      case ('min')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') '# dcdfile    = sample.dcd   # DCD trajectory file'
        write(MsgOut,'(A)') '# rstfile    = sample.rst   # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample.pdb   # PDB file'
        write(MsgOut,'(A)') ' '

      case ('vib')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'minfofile    = vib.minfo    # minfo file'
        write(MsgOut,'(A)') ' '

      case ('remd')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'logfile    = sample{}.log # log file of each replica'
        write(MsgOut,'(A)') '# dcdfile    = sample{}.dcd # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample{}.dcd # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample{}.rst # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample{}.pdb # PDB file'
        write(MsgOut,'(A)') '# remfile    = sample{}.rem # replica exchange ID file'
        write(MsgOut,'(A)') ' '

      case ('rpath')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'logfile    = sample{}.log # log file of each replica'
        write(MsgOut,'(A)') '# dcdfile    = sample{}.dcd # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample{}.dcd # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample{}.rst # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample{}.pdb # PDB file'
        write(MsgOut,'(A)') '# rpathfile  = sample{}.rpath # replica path ID file'
        write(MsgOut,'(A)') '# rpathlogfile = sample.rpathlog # rpathlog file'
        !write(MsgOut,'(A)') '# rstmepfile = sample{}.rstmep # restart file for MEP/FEP'
        write(MsgOut,'(A)') ' '

      case ('bd')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') '# dcdfile    = sample.dcd   # DCD trajectory file'
        write(MsgOut,'(A)') '# rstfile    = sample.rst   # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample.pdb   # PDB file'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('vib')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'minfofile    = vib.minfo    # minfo file'
        write(MsgOut,'(A)') ' '

      case ('remd', 'rpath')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'logfile    = sample{}.log # log file of each replica'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_output
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_output
  !> @brief        read OUTPUT section in the control file
  !! @authors      TM, JJ
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   out_info : OUTPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_output(handle, out_info)

    ! parameters
    character(*),            parameter     :: Section = 'Output'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_out_info),        intent(inout) :: out_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'logfile',   out_info%logfile)
    call read_ctrlfile_string(handle, Section, 'dcdfile',   out_info%dcdfile)
    call read_ctrlfile_string(handle, Section, 'dcdvelfile',out_info%dcdvelfile)
    call read_ctrlfile_string(handle, Section, 'rstfile',   out_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'pdbfile',   out_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'remfile',   out_info%remfile)
    call read_ctrlfile_string(handle, Section, 'rpathfile', out_info%rpathfile)
    call read_ctrlfile_string(handle, Section, 'rstmepfile',out_info%rstmepfile)
    call read_ctrlfile_string(handle, Section, 'minfofile' ,out_info%minfofile)
    call read_ctrlfile_string(handle, Section, 'gamdfile',  out_info%gamdfile)
    call read_ctrlfile_string(handle, Section, 'rpathlogfile',  out_info%rpathlogfile)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Output> Output Files'
      if (out_info%logfile /= '')    &
        write(MsgOut,*) ' logfile    = ', trim(out_info%logfile) 
      if (out_info%dcdfile /= '')    &
        write(MsgOut,*) ' dcdfile    = ', trim(out_info%dcdfile)
      if (out_info%dcdvelfile /= '') &
        write(MsgOut,*) ' dcdvelfile = ', trim(out_info%dcdvelfile)
      if (out_info%rstfile /= '')    &
        write(MsgOut,*) ' rstfile    = ', trim(out_info%rstfile)
      if (out_info%pdbfile /= '')    &
        write(MsgOut,*) ' pdbfile    = ', trim(out_info%pdbfile)
      if (out_info%remfile /= '')    &
        write(MsgOut,*) ' remfile    = ', trim(out_info%remfile)
      if (out_info%rpathfile /= '')  &
        write(MsgOut,*) ' rpathfile  = ', trim(out_info%rpathfile)
      if (out_info%rstmepfile /= '') &
        write(MsgOut,*) ' rstmepfile = ', trim(out_info%rstmepfile)
      if (out_info%minfofile /= '')  &
        write(MsgOut,*) ' minfofile  = ', trim(out_info%minfofile)
      if (out_info%gamdfile /= '')   &
        write(MsgOut,*) ' gamdfile   = ', trim(out_info%gamdfile)
      if (out_info%rpathfile /= '')  &
        write(MsgOut,*) ' rpathlogfile  = ', trim(out_info%rpathlogfile)
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_md
  !> @brief        setup output information for MD
  !! @authors      TM
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_md(out_info, dynamics, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_output),          intent(inout) :: output


    output%verbose = dynamics%verbose

    output%replica = .false.
    output%rpath   = .false.
    output%logout  = .false.

    if (dynamics%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') then
        call error_msg('Setup_Output_Md> Error: dcdfile name is not'//       &
                 ' specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      else
        output%dcdout  = .true.
        output%dcdfile = out_info%dcdfile
      end if
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') then
        call error_msg('Setup_Output_Md> Error: dcdvelfile name is not'//   &
                 ' specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      else
        output%dcdvelout  = .true.
        output%dcdvelfile = out_info%dcdvelfile
      end if
    end if

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') then
        call error_msg('Setup_Output_Md> Error: rstfile name is not'//&
                ' specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      else
        output%rstout  = .true.
        output%rstfile = out_info%rstfile
      end if

      if (out_info%pdbfile .eq. '') then
        output%pdbout  = .false.
      else
        output%pdbout  = .true.
        output%pdbfile = out_info%pdbfile
      end if
    end if

    if (out_info%gamdfile .ne. '') then
      output%gamdout  = .true.
      output%gamdfile = out_info%gamdfile
    end if

    return

  end subroutine setup_output_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_min
  !> @brief        setup output information for minimization
  !! @authors      TM
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    minimize : minimize information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_min(out_info, minimize, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_minimize),        intent(in)    :: minimize
    type(s_output),          intent(inout) :: output


    output%verbose    = minimize%verbose

    output%replica    = .false.
    output%rpath      = .false.
    output%logout     = .false.
    output%dcdvelout  = .false.

    if (minimize%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') then
        call error_msg('Setup_Output_Min> Error: dcdfile name is not'//     &
                 ' specified in [OUTPUT] (crdout_period > 0 in [MINIMIZE])')
      else
        output%dcdout  = .true.
        output%dcdfile = out_info%dcdfile
      end if
    end if

    if (minimize%rstout_period > 0) then
      if (out_info%rstfile .eq. '') then
        call error_msg('Setup_Output_Min> Error: rstfile name is not'//    &
                 ' specified in [OUTPUT] (rstout_period > 0 in [MINIMIZE])')
      else
        output%rstout  = .true.
        output%rstfile = out_info%rstfile
      end if

      if (out_info%pdbfile .eq. '') then
        output%pdbout  = .false.
      else
        output%pdbout  = .true.
        output%pdbfile = out_info%pdbfile
      end if
    end if

    return

  end subroutine setup_output_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_remd
  !> @brief        setup output information for REMD
  !! @authors      TM
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    remd     : REMD information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_remd(out_info, dynamics, remd, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_remd),            intent(in)    :: remd
    type(s_output),          intent(inout) :: output


    output%verbose = dynamics%verbose

    output%logfile = out_info%logfile
    call include_id_to_filename(output%logfile)
    output%logfile = output%logfile
    output%logout  = .true.

    if (dynamics%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: dcdfile name is not'//     &
                 ' specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      output%dcdfile = out_info%dcdfile
      call include_id_to_filename(output%dcdfile)
      output%dcdfile = output%dcdfile
      output%dcdout  = .true.

      if (out_info%remfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: remfile name is not'//     &
                 ' specified in [OUTPUT]')
      output%remfile = out_info%remfile
      call include_id_to_filename(output%remfile)
      output%remfile = output%remfile
      output%remout  = .true.
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: dcdvelfile name is not'//  &
                 ' specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      output%dcdvelfile = out_info%dcdvelfile
      call include_id_to_filename(output%dcdvelfile)
      output%dcdvelfile = output%dcdvelfile
      output%dcdvelout  = .true.
    end if

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: rstfile name is not'//     &
                 'specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      output%rstfile = out_info%rstfile
      call include_id_to_filename(output%rstfile)
      output%rstfile = output%rstfile
      output%rstout  = .true.

      if (remd%exchange_period > 0) then
        if (mod(dynamics%rstout_period,remd%exchange_period) /= 0) &
        call error_msg('Setup_Output_Remd> Error in rstout_period'//&
        '  mod(rstout_period, exchange_period) is not ZERO')
      end if

      if (out_info%pdbfile .ne. '') then
        output%pdbfile = out_info%pdbfile
        call include_id_to_filename(output%pdbfile)
        output%pdbfile = output%pdbfile
        output%pdbout  = .true.
      end if
    end if

    output%replica = .true.

    return

  end subroutine setup_output_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_rpath
  !> @brief        setup output information for RPATH
  !! @authors      TM, KY, YA
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[inout] rpath    : rpath information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_rpath(out_info, dynamics, rpath, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_rpath),           intent(inout) :: rpath
    type(s_output),          intent(inout) :: output


    output%verbose = dynamics%verbose

    output%logfile = out_info%logfile
    if (output%logfile .eq. '') then
      call error_msg('Setup_Output_Rpath> log filename is blank')
    end if
    call include_id_to_filename(output%logfile)
    output%logout  = .true.

    if (dynamics%crdout_period > 0 .or. rpath%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') &
        call error_msg('Setup_Output_Rpath> dcd filename is blank')
      output%dcdfile = out_info%dcdfile
      call include_id_to_filename(output%dcdfile)
      output%dcdout  = .true.
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') &
        call error_msg('Setup_Output_Rpath> dcdvel filename is blank')
      output%dcdvelfile = out_info%dcdvelfile
      call include_id_to_filename(output%dcdvelfile)
      output%dcdvelout  = .true.
    end if

    if (dynamics%rstout_period > 0 .or. rpath%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Rpath> rst filename is blank')
      output%rstfile = out_info%rstfile
      call include_id_to_filename(output%rstfile)
      output%rstout  = .true.

      if (rpath%rpath_period > 0) then
        if (mod(dynamics%rstout_period,rpath%rpath_period) /= 0) &
        call error_msg('Setup_Output_Rpath> Error in rstout_period'//&
        '  mod(rstout_period, rpath_period) is not ZERO')
      end if

      if (out_info%pdbfile .ne. '') then
        output%pdbfile = out_info%pdbfile
        call include_id_to_filename(output%pdbfile)
        output%pdbout  = .true.
      end if
    end if

    if (dynamics%crdout_period > 0 .or. rpath%crdout_period > 0) then
      if (out_info%rpathfile .ne. '') then
        output%rpathfile = out_info%rpathfile
        call include_id_to_filename(output%rpathfile)
        output%rpathout  = .true.
      end if
    end if

    if (rpath%rpathmode == RpathmodeFEP) then
      if (out_info%rstmepfile .eq. '') &
        call error_msg('Setup_Output_Rpath> rstmep filename is blank')
      output%rstmepfile = out_info%rstmepfile
      call include_id_to_filename(output%rstmepfile)
      output%rstmepout  = .true.
    end if
    if (dynamics%crdout_period > 0) then
      if (out_info%rpathlogfile .ne. '') then
        output%rpathlogfile = out_info%rpathlogfile
        output%rpathlogout  = .true.
      end if
    end if

    output%rpath = .true.

    return

  end subroutine setup_output_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_vib
  !> @brief        setup output information for minimization
  !! @authors      KY
  !! @param[in]    out_info  : OUTPUT section control parameters information
  !! @param[in]    vibration : vibration information
  !! @param[out]   output    : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_vib(out_info, vibration, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_vibration),       intent(inout) :: vibration
    type(s_output),          intent(inout) :: output


    ! local variables
    integer    :: file
    integer    :: ierr
    logical    :: ex


    output%replica    = .false.
    output%rpath      = .false.
    output%dcdvelout  = .false.
    output%vib        = .true.

    if (out_info%logfile .ne. '') then
      output%logfile = out_info%logfile
      call include_id_to_filename(output%logfile)
      output%logout  = .true.

    else
      output%logout  = .false.

    end if

    if (.not. vibration%gengrid) then
      if (out_info%minfofile .eq. '') then
        call error_msg('Setup_Output_Vib> Error: minfofile name is not'//     &
                 ' specified in [OUTPUT] ')
      else
        vibration%minfofile = out_info%minfofile
      end if

      ! Check whether minfofile already exists, and stop with error
      ! to avoid overwrite.
      !
      inquire(file=trim(vibration%minfofile), exist=ex)
      if(ex) then
        call error_msg('Setup_Output_Vib> Error: minfo file [ '// &
                 trim(vibration%minfofile)// &
                 ' ] already exists. ')
      end if

    else if(vibration%grid_ene_only) then

      ! Check whether datafile already exist, and stop with error
      ! to avoid overwrite.
      !
      !inquire(file=trim(vibration%datafile), exist=ex)
      !if(ex) then
      !  call error_msg('Setup_Output_Vib> Error: data file [ '// &
      !           trim(vibration%datafile)// &
      !           ' ] already exists. ')
      !end if

    end if

    return

  end subroutine setup_output_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_morph
  !> @brief        setup output information for minimization
  !! @authors      CK
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    morph    : morph information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_morph(out_info, morph, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_morph),           intent(in)    :: morph
    type(s_output),          intent(inout) :: output


    output%replica    = .false.
    output%rpath      = .false.
    output%logout     = .false.
    output%dcdvelout  = .false.

    if (morph%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') then
        call error_msg('Setup_Output_Morph> Error: dcdfile name is not specified in [OUTPUT] (crdout_period > 0 in [MORPH])')
      else
        output%dcdout  = .true.
        output%dcdfile = out_info%dcdfile
      end if
    end if

    if (morph%rstout_period > 0) then
      if (out_info%rstfile .eq. '') then
        call error_msg('Setup_Output_Morph> Error: rstfile name is not specified in [OUTPUT] (rstout_period > 0 in [MORPH])')
      else
        output%rstout  = .true.
        output%rstfile = out_info%rstfile
      end if

      if (out_info%pdbfile .eq. '') then
        output%pdbout  = .false.
      else
        output%pdbout  = .true.
        output%pdbfile = out_info%pdbfile
      end if
    end if

    return

  end subroutine setup_output_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_output
  !> @brief        open output file
  !! @authors      TM, SI
  !! @param[inout] outout : information of output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_output(output)

    ! formal arguments
    type(s_output),          intent(inout) :: output

    ! local variables
    integer                  :: file
    integer                  :: access


    ! open logfile 
    !
    if (output%logout) then
      if (replica_main_rank) then
        if (access(output%logfile, 'rw' ) == 0) then
          call open_file(file, output%logfile, IOFileOutputAppend)
        else
          call open_file(file, output%logfile, IOFileOutputNew)
        end if
        output%logunit = file
        DynvarsOut = output%logunit
      end if
    end if

    ! open dcdfile
    !
    if (output%dcdout) then
      if (main_rank .or. replica_main_rank) then
        if(access(output%dcdfile, 'rw' ) == 0) then
          call open_binary_file(file, output%dcdfile, IOFileOutputAppend)
        else
          call open_binary_file(file, output%dcdfile, IOFileOutputNew)
        end if
        output%dcdunit = file
      end if
    end if

    ! open dcdvelfile
    !
    if (output%dcdvelout) then
      if (main_rank .or. replica_main_rank) then            
        if(access(output%dcdvelfile, 'rw' ) == 0) then      
          call open_binary_file(file, output%dcdvelfile, IOFileOutputAppend)
        else
          call open_binary_file(file, output%dcdvelfile, IOFileOutputNew)
        end if
        output%dcdvelunit = file
      end if
    end if

    if (nrep_per_proc == 1) then
      ! open rstfile
      !   Just check whether rstfile is already existing.
      !
      if (output%rstout) then
        if (main_rank .or. replica_main_rank) then
          call open_file (file, output%rstfile, IOFileOutputNew)
          call close_file(file)
        end if
      end if

      ! open pdbfile
      !   Just check whether pdbfile is already existing.
      !
      if (output%pdbout) then
        if (main_rank .or. replica_main_rank) then
          call open_file (file, output%pdbfile, IOFileOutputNew)
          call close_file(file)
        end if
      end if
    end if

    ! open remfile
    !
    if (output%remout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%remfile, IOFileOutputNew)
        output%remunit = file
      end if
    end if

    ! open rpathfile
    !
    if (output%rpathout) then
      if (main_rank .or. replica_main_rank) then
        if(access(output%rpathfile, 'rw' ) == 0) then    
          call open_file(file, output%rpathfile, IOFileOutputAppend)
        else
          call open_file(file, output%rpathfile, IOFileOutputNew)
        end if 
      output%rpathunit = file
      end if
    end if

    ! open rstmepfile
    !   Just check whether rstmepfile is already existing.
    !
    if (output%rstmepout) then
      if (main_rank .or. replica_main_rank) then
        call open_file (file, output%rstmepfile, IOFileOutputNew)
        call close_file(file)
      end if
    end if

    ! open gamdfile
    !
    if (output%gamdout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%gamdfile, IOFileOutputNew)
        output%gamdunit = file
      end if
    end if

    ! open rpathlogfile (only main_rank)
    !
    if (output%rpathlogout) then
      if (main_rank) then
        call open_file(file, output%rpathlogfile, IOFileOutputNew)
        output%rpathlogunit = file
      end if
    end if

    return

  end subroutine open_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    close_output
  !> @brief        close output files
  !! @authors      TM
  !! @param[in]    output : information of output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine close_output(output)

    ! formal arguments
    type(s_output),          intent(in)    :: output


    ! close remfile
    !
    if (output%remout) &
      call close_file(output%remunit)

    ! close rpathfile
    !
    if (output%rpathout) &
      call close_file(output%rpathunit)

    ! close dcdfile
    !
    if (output%dcdout) &
      call close_file(output%dcdunit)

    ! close dcdvelfile
    !
    if (output%dcdvelout) &
      call close_file(output%dcdvelunit)

    ! close logfile
    !
    if (output%logout) &
      call close_file(output%logunit)

    ! close rpathlogfile
    !
    if (output%rpathlogout) &
      call close_file(output%rpathlogunit)

    return

  end subroutine close_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_md
  !> @brief        output trajectory and restart data for MD
  !! @authors      TM, NT
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_md(output, molecule, enefunc, dynamics, boundary, &
                       ensemble, dynvars)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(wp)                 :: dt, factor
    integer                  :: natom, istep, j

    real(wp),        pointer :: vel(:,:), vel_ref(:,:)


    natom = molecule%num_atoms
    istep = dynvars%step

    if (dynamics%eneout_period > 0) then
      if (mod(istep,dynamics%eneout_period) == 0) then

        call compute_dynvars(molecule, enefunc, dynamics, boundary, &
                             ensemble, dynvars)
        call output_dynvars (output, enefunc, dynvars, ensemble)

      end if
    end if

    if (dynamics%crdout_period > 0) then
      if (mod(istep,dynamics%crdout_period) == 0) then

        if (dynamics%integrator == IntegratorLEAP) then

          do j = 1, natom
            dynvars%temporary(1,j) = dynvars%coord_ref(1,j)
            dynvars%temporary(2,j) = dynvars%coord_ref(2,j)
            dynvars%temporary(3,j) = dynvars%coord_ref(3,j)
          end do

        else

          do j = 1, natom
            dynvars%temporary(1,j) = dynvars%coord(1,j)
            dynvars%temporary(2,j) = dynvars%coord(2,j)
            dynvars%temporary(3,j) = dynvars%coord(3,j)
          end do

        end if

        call wrap_molecules(molecule, boundary, dynvars%temporary)
        call write_trajectory_dcd(output, molecule, boundary, istep, &
                                  dynamics%nsteps, &
                                  dynamics%crdout_period, &
                                  dynamics%timestep, &
                                  dynvars%temporary)
      end if

    end if

    if (dynamics%velout_period > 0) then
      if (mod(istep,dynamics%velout_period) == 0) then

        if (dynamics%integrator == IntegratorLEAP) then

          dt =  dynamics%timestep/AKMA_PS
          if (ensemble%tpcontrol == TpcontrolLangevin) then
            factor = sqrt(1.0_wp+0.5_wp*ensemble%gamma_t*dt)
          else
            factor = 1.0_wp
          end if

          vel     => dynvars%velocity
          vel_ref => dynvars%velocity_ref

          do j = 1, natom
            dynvars%temporary(1,j) = 0.5_wp*(vel(1,j)+vel_ref(1,j))*factor
            dynvars%temporary(2,j) = 0.5_wp*(vel(2,j)+vel_ref(2,j))*factor
            dynvars%temporary(3,j) = 0.5_wp*(vel(3,j)+vel_ref(3,j))*factor
          end do

        else

          do j = 1, natom
            dynvars%temporary(1,j) = dynvars%velocity(1,j)
            dynvars%temporary(2,j) = dynvars%velocity(2,j)
            dynvars%temporary(3,j) = dynvars%velocity(3,j)
          end do

        end if

        call write_trajectory_dcdvel(output, molecule, boundary, istep, &
                                     dynamics%nsteps, &
                                     dynamics%velout_period, &
                                     dynamics%timestep, &
                                     dynvars%temporary)

      end if
    end if

    if (dynamics%rstout_period > 0) then
      if (mod(istep,dynamics%rstout_period) == 0) then

        call output_restart_md(output, molecule, dynamics, boundary, dynvars, &
                               enefunc%qmmm)

      end if
    end if

    ! output QM charge
    !
    if (dynamics%esp_mm) then
      if (istep == 1) then
          call output_qm_charge(molecule, enefunc%qmmm, output)
      else if (dynamics%calc_qm_period /= 0) then
        if (mod(istep,dynamics%calc_qm_period) == 0) &
          call output_qm_charge(molecule, enefunc%qmmm, output)
      end if
    end if

    return

  end subroutine output_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_min
  !> @brief        output trajectory and restart data for miniization
  !! @authors      TM, NT, KY
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    minimize : minimization information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    delta_r  : delta-r
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_min(output, molecule, enefunc, minimize, boundary, &
                        delta_r, maxg_id, dynvars)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_minimize),        intent(in)    :: minimize
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: delta_r
    integer,                 intent(in)    :: maxg_id
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: natom, istep, j, outunit


    natom = molecule%num_atoms
    istep = dynvars%step

    if (minimize%eneout_period > 0 .and. replica_main_rank) then
      if (mod(istep,minimize%eneout_period) == 0 .or. minimize%eneout) then

        if (output%logout) then
          outunit = output%logunit
        else
          outunit = MsgOut
        end if

        if (.not. minimize%eneout_short) then
          call output_dynvars(output, enefunc, dynvars)
          write(outunit, &
            '("  maxg = ",f10.4," at ",i10,x,a4,x,i6,x,a4,x,a6,x,a6,/)') &
            dynvars%max_gradient, maxg_id,  &
            molecule%segment_name(maxg_id), &
            molecule%residue_no(maxg_id),   &
            molecule%residue_name(maxg_id), &
            molecule%atom_name(maxg_id),    &
            molecule%atom_cls_name(maxg_id)

        else
          write(outunit, '(5X, "STEP ", i8, " POTENTIAL_ENE = ", F20.8, &
                           "  RMSG = ", F12.4, "  MAXG = ", F12.4)')    &
            istep, dynvars%energy%total, dynvars%rms_gradient,          &
            dynvars%max_gradient

        end if

      end if
    end if

    ! Do Not print the initial structure
    if (istep == 0) return

    if (minimize%crdout_period > 0) then
      if (mod(istep,minimize%crdout_period) == 0 .or. &
          minimize%crdout) then

        do j = 1, natom
          dynvars%temporary(1,j) = dynvars%coord_ref(1,j)
          dynvars%temporary(2,j) = dynvars%coord_ref(2,j)
          dynvars%temporary(3,j) = dynvars%coord_ref(3,j)
        end do

        call wrap_molecules(molecule, boundary, dynvars%temporary)
        call write_trajectory_dcd(output, molecule, boundary, istep, &
                                  minimize%nsteps, &
                                  minimize%crdout_period, &
                                  0.001_wp, &
                                  dynvars%temporary)
      end if
    end if

    if (minimize%rstout_period > 0) then
      if (mod(istep,minimize%rstout_period) == 0 .or. &
          minimize%rstout) then

        call output_restart_min(output, molecule, minimize, boundary, dynvars, &
                                enefunc%qmmm, delta_r)

      end if
    end if

    return

  end subroutine output_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_remd
  !> @brief        output restart data for REMD
  !! @authors      TM, NT
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[inout] remd     : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_remd(icycle, output, molecule, dynamics, dynvars, boundary,&
                         remd)

    ! formal arguments
    integer,                 intent(in)    :: icycle
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_remd),            intent(inout) :: remd


    ! output remfile
    !
    if (replica_main_rank .and. output%remout) then
      if (icycle /= remd%ncycles) then
        write(output%remunit,'(2I10)')          &
          dynvars%step,                         &
          remd%repid2parmsetid(my_country_no+1)
      end if
    end if

    ! output restart
    !
    if (dynamics%rstout_period > 0) then
      if (dynvars%step > 0 .and. mod(dynvars%step,dynamics%rstout_period) == 0) then
        call output_restart_remd(output, molecule, dynamics, dynvars, &
                                 boundary, remd)
      end if
    end if

    return

  end subroutine output_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_rpath
  !> @brief        output restart data for RPATH
  !! @authors      NT, YA, KY
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[in]    rpath     : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_rpath(output, molecule, enefunc, dynamics, dynvars, &
                          boundary, rpath)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_rpath),   target, intent(inout) :: rpath

    ! local variables
    type(s_rstmep)     :: rstmep
    integer            :: repid_i, i, ii, j, k


    ! wrire log (main rank only)
    !
    if (main_rank .and. output%rpathlogout) then
      write(output%rpathlogunit, '(I10,3F15.8)') &
           dynvars%step, rpath%sum_distance, &
           rpath%distance_init,rpath%distance_prev
    end if

    repid_i = my_country_no + 1
    if (nrep_per_proc > 1) repid_i = my_replica_no 

    if (rpath%rpathmode == RpathmodeMFEP) then

      ! write collective variables per replica
      !
      if (replica_main_rank .and. output%rpathout) then
        write(output%rpathunit, '(I10,$)') dynvars%step
        do k = 1, rpath%dimension
          write(output%rpathunit, '(E15.5,$)') rpath%rest_reference(1,k,repid_i)
        end do
        write(output%rpathunit, *)
      end if

      ! output restart
      !
      if (dynamics%rstout_period > 0) then
        if (mod(dynvars%step,dynamics%rstout_period) == 0) then
          call output_restart_rpath(output, molecule, dynamics, &
                                    dynvars, boundary, enefunc%qmmm, rpath)
        end if
      end if

    else if (rpath%rpathmode == RpathmodeMEPMD) then

      ! write collective variables per replica
      !
      if (replica_main_rank .and. output%rpathout) then
        write(output%rpathunit, '(I10,$)') dynvars%step
        do k = 1, rpath%dimension
          write(output%rpathunit, '(E15.5,$)') rpath%mep_coord(k,repid_i)
        end do
        write(output%rpathunit, *)
      end if

      ! output restart
      !
      if (dynamics%rstout_period > 0) then
        if (mod(dynvars%step,dynamics%rstout_period) == 0) then
          call output_restart_rpath(output, molecule, dynamics, &
                                    dynvars, boundary, enefunc%qmmm, rpath)
        end if
      end if

    else if (rpath%rpathmode == RpathmodeMEP) then

      ! write collective variables per replica
      !
      if (replica_main_rank .and. output%rpathout .and. &
          rpath%neb_output) then
        write(output%rpathunit, '(I10,$)') dynvars%step
        do k = 1, rpath%dimension
          write(output%rpathunit, '(E15.5,$)') rpath%mep_coord(k,repid_i)
        end do
        write(output%rpathunit, *)
      end if

      ! output energy
      !
      if (rpath%neb_output .and. rpath%eneout_period > 0) then
        if (mod(dynvars%step,rpath%eneout_period) == 0 .or. &
            rpath%eneout) then
          dynvars%total_pene = dynvars%energy%total
          call output_dynvars(output, enefunc, dynvars)
        end if
      end if

      ! output restart
      !
      if (rpath%rstout_period > 0) then
        if (mod(dynvars%step,rpath%rstout_period) == 0 .or. &
            rpath%rstout) then
          call output_restart_rpath(output, molecule, dynamics, &
                                    dynvars, boundary, enefunc%qmmm, rpath)
        end if
      end if

      ! output coordinates
      !
      if (rpath%neb_output .and. rpath%crdout_period > 0) then
        if ((mod(dynvars%step,rpath%crdout_period) == 0 .or. rpath%crdout) &
            .and. (.not. rpath%first_iter)) then
          do j = 1, molecule%num_atoms
            dynvars%temporary(1,j) = dynvars%coord(1,j)
            dynvars%temporary(2,j) = dynvars%coord(2,j)
            dynvars%temporary(3,j) = dynvars%coord(3,j)
          end do
          call wrap_molecules(molecule, boundary, dynvars%temporary)
          call write_trajectory_dcd(output, molecule, boundary, dynvars%step, &
                 rpath%ncycle, rpath%crdout_period, 0.0_wp, dynvars%temporary)
        end if
      end if

    else if (rpath%rpathmode == RpathmodeFEP) then

      ! output rstmep information
      !
      if (replica_main_rank .and. output%rstmepout) then
        rstmep%mep_natoms = rpath%mep_natoms
        allocate(rstmep%mep_coord(3*rpath%mep_natoms))
        rstmep%mep_coord  = rpath%mep_coord(:,repid_i)
        if (rpath%esp_energy .or. rpath%esp_md) then
          rstmep%qm_natoms = enefunc%qmmm%qm_natoms
          allocate(rstmep%qm_charge(rstmep%qm_natoms))
          rstmep%qm_charge = rpath%qm_charge(:,repid_i)
          rstmep%qm_energy = rpath%qm_energy(repid_i)
        end if
        allocate(rstmep%pfunc(2))
        rstmep%num_fep   = rpath%num_fep
        rstmep%pfunc(1)  = rpath%pfunc_f
        rstmep%pfunc(2)  = rpath%pfunc_b

        call output_rstmep(output%rstmepfile, rstmep)

        deallocate(rstmep%mep_coord, rstmep%pfunc)
        if (rpath%opt_micro) deallocate(rstmep%qm_charge)
        
      end if

      ! output restart
      !
      if (dynamics%rstout_period > 0 .and. &
          mod(dynvars%step,dynamics%rstout_period) == 0) then
        call output_restart_rpath(output, molecule, dynamics, &
                                  dynvars, boundary, enefunc%qmmm, rpath)
      end if

    else
      call error_msg('output_path> unrecognized rpathmode')

    end if

    return

  end subroutine output_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_vib
  !> @brief        output trajectory and restart data for miniization
  !! @authors      KY
  !! @param[inout] output    : output information
  !! @param[in]    molecule  : molecule information
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    vibration : vibration information
  !! @param[in]    boundary  : boundary condition information
  !! @param[inout] dynvars   : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_vib(output, molecule, enefunc, vibration, boundary, dynvars)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_vibration),       intent(in)    :: vibration
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(inout) :: dynvars


    if (.not. output%logout ) return

    dynvars%total_pene = dynvars%energy%total
    call output_dynvars(output, enefunc, dynvars)

    return

  end subroutine output_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_bd
  !> @brief        output trajectory and restart data for BD
  !! @authors      TA
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_bd(output, molecule, enefunc, dynamics, boundary, &
                       ensemble, dynvars)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_dynvars),         intent(inout) :: dynvars

    ! local variables
    real(wp)                 :: dt, factor
    integer                  :: natom, istep, j


    natom = molecule%num_atoms
    istep = dynvars%step

    if (dynamics%eneout_period > 0) then
      if (mod(istep,dynamics%eneout_period) == 0) then

        call compute_dynvars(molecule, enefunc, dynamics, boundary, &
                             ensemble, dynvars)
        call output_dynvars (output, enefunc, dynvars, ensemble)

      end if
    end if

    if (dynamics%crdout_period > 0) then
      if (mod(istep,dynamics%crdout_period) == 0) then

        do j = 1, natom
          dynvars%temporary(1,j) = dynvars%coord(1,j)
          dynvars%temporary(2,j) = dynvars%coord(2,j)
          dynvars%temporary(3,j) = dynvars%coord(3,j)
        end do

        call wrap_molecules(molecule, boundary, dynvars%temporary)
        call write_trajectory_dcd(output, molecule, boundary, istep, &
                                  dynamics%nsteps, &
                                  dynamics%crdout_period, &
                                  dynamics%timestep, &
                                  dynvars%temporary)
      end if

    end if

    if (dynamics%velout_period > 0) then
      if (mod(istep,dynamics%velout_period) == 0) then

        do j = 1, natom
          dynvars%temporary(1,j) = dynvars%velocity(1,j)
          dynvars%temporary(2,j) = dynvars%velocity(2,j)
          dynvars%temporary(3,j) = dynvars%velocity(3,j)
        end do

        call write_trajectory_dcdvel(output, molecule, boundary, istep, &
                                     dynamics%nsteps, &
                                     dynamics%velout_period, &
                                     dynamics%timestep, &
                                     dynvars%temporary)

      end if
    end if

    if (dynamics%rstout_period > 0) then
      if (mod(istep,dynamics%rstout_period) == 0) then

        call output_restart_bd(output, molecule, dynamics, boundary, dynvars)

      end if
    end if

    return

  end subroutine output_bd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_bd
  !> @brief        setup output information for BD
  !! @authors      TA
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_bd(out_info, dynamics, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_output),          intent(inout) :: output


    output%replica = .false.
    output%logout  = .false.

    if (dynamics%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') then
        call error_msg('Setup_Output_Bd> Error: dcdfile name is not specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      else
        output%dcdout  = .true.
        output%dcdfile = out_info%dcdfile
      end if
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') then
        call error_msg('Setup_Output_Bd> Error: dcdvelfile name is not specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      else
        output%dcdvelout  = .true.
        output%dcdvelfile = out_info%dcdvelfile
      end if
    end if

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') then
        call error_msg('Setup_Output_Bd> Error: rstfile name is not specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      else
        output%rstout  = .true.
        output%rstfile = out_info%rstfile
      end if

      if (out_info%pdbfile .eq. '') then
        output%pdbout  = .false.
      else
        output%pdbout  = .true.
        output%pdbfile = out_info%pdbfile
      end if
    end if

    return

  end subroutine setup_output_bd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_bd
  !> @brief        output restart data for BD
  !! @authors      TA
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_bd(output, molecule, dynamics, boundary, dynvars)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars

    ! local variables
    type(s_rst)              :: rst
    integer                  :: natom


    if (output%replica) &
      return

    if (output%rpath) &
      return

    call random_push_stock
    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    if (.not. main_rank) &
      return

    ! Output RST file
    !
    natom = molecule%num_atoms
    call alloc_rst(rst, RestartAtom, natom)

    rst%rstfile_type = RstfileTypeBd

    rst%iseed                  = dynamics%iseed
    rst%num_atoms              = molecule%num_atoms
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars%thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars%barostat_momentum(1:3)
    rst%coord(1:3,1:natom)     = dynvars%coord(1:3,1:natom)
    rst%velocity(1:3,1:natom)  = dynvars%velocity(1:3,1:natom)

    call output_rst(output%rstfile, rst)

    if (dynvars%step == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, dynvars%coord)

    return

  end subroutine output_restart_bd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_morph
  !> @brief        output restart data for MOPRH
  !! @authors      CK
  !! @param[inout] output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    morph    : RPATH information
  !! @param[in]    boundary : boundary information
  !! @param[in]    enefunc  : potential energy function information
  !! @param[in]    dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_morph(output, molecule, morph, boundary, enefunc, dynvars)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_morph),           intent(inout) :: morph
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars),         intent(inout) :: dynvars

    ! local variables
    integer     :: j, natom
    type(s_rst) :: rst


    if (.not. main_rank) return


    ! use pointers
    !
    natom = molecule%num_atoms

    ! output energy
    !

    if (dynvars%step > 0) then
      dynvars%total_pene = dynvars%energy%total
      call output_dynvars(output, enefunc, dynvars)
    end if


    ! output coordinates
    !
    if (morph%crdout_period > 0) then
      if (mod(morph%ncycles,morph%crdout_period) == 0) then

        do j = 1, natom
          dynvars%temporary(1,j) = dynvars%coord_ref(1,j)
          dynvars%temporary(2,j) = dynvars%coord_ref(2,j)
          dynvars%temporary(3,j) = dynvars%coord_ref(3,j)
        end do

        call wrap_molecules(molecule, boundary, dynvars%temporary)
        call write_trajectory_dcd(output, molecule, boundary, dynvars%step, &
                                  morph%ncycles, &
                                  morph%crdout_period, &
                                  morph%delta_r, &
                                  dynvars%temporary)
      end if
    end if


    ! output restart data
    !
    if (morph%rstout_period > 0) then
      if (mod(dynvars%step,morph%rstout_period) == 0) then

        call output_restart_morph(output, molecule, morph, dynvars, boundary)

      end if
    end if

    return

  end subroutine output_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_gamd
  !> @brief        output GaMD boost parameters
  !! @authors      HO
  !! @param[inout] output  : output information
  !! @param[in]    dynvars : dynamic variables information
  !! @param[in]    enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_gamd(output, dynvars, enefunc)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_enefunc),target,  intent(in)    :: enefunc

    ! local variables
    integer,parameter        :: clength=16, flength=4

    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, rfrmt*7
    character                :: rfrmt_cont*9, frmt_cont*7
    real(wp)                 :: values(999)
    logical, save            :: gamdtitle = .true.

    type(s_enefunc_gamd), pointer :: gamd


    gamd => enefunc%gamd

    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    if (output%replica) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if

    if (output%gamdout) then

      ifm = 1

      if (gamd%boost_pot .or. gamd%boost_dual) then

        write(category(ifm),frmt) 'POT_MAX'
        values(ifm) = gamd%ene_pot_max
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_MIN'
        values(ifm) = gamd%ene_pot_min
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_AVE'
        values(ifm) = gamd%ene_pot_ave
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_DEV'
        values(ifm) = gamd%ene_pot_dev
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_TH'
        values(ifm) = gamd%ene_pot_th
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_KC'
        values(ifm) = gamd%k_pot
        ifm = ifm+1

        write(category(ifm),frmt) 'POT_KC0'
        values(ifm) = gamd%k0_pot
        ifm = ifm+1
      end if

      if (gamd%boost_dih .or. gamd%boost_dual) then

        write(category(ifm),frmt) 'DIH_MAX'
        values(ifm) = gamd%ene_dih_max
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_MIN'
        values(ifm) = gamd%ene_dih_min
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_AVE'
        values(ifm) = gamd%ene_dih_ave
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_DEV'
        values(ifm) = gamd%ene_dih_dev
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_TH'
        values(ifm) = gamd%ene_dih_th
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_KC'
        values(ifm) = gamd%k_dih
        ifm = ifm+1

        write(category(ifm),frmt) 'DIH_KC0'
        values(ifm) = gamd%k0_dih
        ifm = ifm+1
      end if

      if (gamdtitle) then
        write(output%gamdunit,'(A,$)') title
        do i = 1, ifm-1
          if (i == ifm-1) then
            write(output%gamdunit,frmt) category(i)
          else
            write(output%gamdunit,frmt_cont) category(i)
          end if
        end do
        write(output%gamdunit,'(A80)') ' --------------- --------------- --------------- --------------- ---------------'
        gamdtitle = .false.
      end if

      write(output%gamdunit,'(6x,I10,$)') dynvars%step

      do i = 1, ifm-1
        if (i == ifm-1) then
          write(output%gamdunit,rfrmt) values(i)
        else
          write(output%gamdunit,rfrmt_cont) values(i)
        end if
      end do

      write(output%gamdunit,'(A)') ''
    end if

    return

  end subroutine output_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_md
  !> @brief        output restart data for MD
  !! @authors      TM
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    qmmm     : QMMM variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_md(output, molecule, dynamics, boundary, dynvars, &
                               qmmm)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    type(s_rst)              :: rst
    integer                  :: natom, nrep, ndim


    if (output%replica) &
      return

    if (output%rpath) &
      return

    call random_push_stock
    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    if (.not. main_rank) &
      return

    ! Output RST file
    !
    natom = molecule%num_atoms
    call alloc_rst(rst, RestartAtom, natom)


    rst%rstfile_type = RstfileTypeMd

    rst%iseed                  = dynamics%iseed
    rst%num_atoms              = molecule%num_atoms
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars%thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars%barostat_momentum(1:3)
    rst%coord(1:3,1:natom)     = dynvars%coord(1:3,1:natom)
    rst%velocity(1:3,1:natom)  = dynvars%velocity(1:3,1:natom)

    rst%sph_pot = boundary%sph_pot
    if (boundary%sph_pot) then
      call alloc_rst(rst, RestartSpot, boundary%nfunctions, natom)
      rst%nfunctions = boundary%nfunctions
      rst%radius     = boundary%radius
      rst%center     = boundary%center
      rst%const      = boundary%const
      rst%exponent   = boundary%exponent
      rst%fix_layer  = boundary%fix_layer
      rst%num_fixatm = boundary%num_fixatm
      rst%fixatm     = boundary%fixatm
    end if

    if (qmmm%do_qmmm .and. qmmm%is_qm_charge) then
      call alloc_rst(rst, RestartQMMM, qmmm%qm_natoms)
      rst%qm_charge = qmmm%qm_charge
    end if

    call output_rst(output%rstfile, rst)

    if (dynvars%step == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, dynvars%coord)

    return

  end subroutine output_restart_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_min
  !> @brief        output restart data for minimization
  !! @authors      TM
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    minimize : minimization information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    qmmm     : QMMM variables information
  !! @param[in]    delta_r  : delta-r
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_min(output, molecule, minimize, boundary, dynvars, &
                                qmmm, delta_r)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_minimize),        intent(in)    :: minimize
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_qmmm),            intent(in)    :: qmmm
    real(wp),                intent(in)    :: delta_r

    ! local variables
    type(s_rst)              :: rst
    integer                  :: natom


    if (output%replica) &
      return

    if (output%rpath) &
      return

    if (.not. main_rank) &
      return

    ! Output RST file
    !
    natom = molecule%num_atoms
    call alloc_rst(rst, RestartAtom, natom)

    rst%rstfile_type       = RstfileTypeMin
    rst%num_atoms          = molecule%num_atoms
    rst%box_size_x         = boundary%box_size_x
    rst%box_size_y         = boundary%box_size_y
    rst%box_size_z         = boundary%box_size_z
    rst%energy             = dynvars%energy%total
    rst%delta_r            = delta_r
    rst%coord(1:3,1:natom) = dynvars%coord_ref(1:3,1:natom)

    rst%sph_pot = boundary%sph_pot
    if (boundary%sph_pot) then
      call alloc_rst(rst, RestartSpot, boundary%nfunctions, natom)
      rst%nfunctions = boundary%nfunctions
      rst%radius     = boundary%radius
      rst%center     = boundary%center
      rst%const      = boundary%const
      rst%exponent   = boundary%exponent
      rst%fix_layer  = boundary%fix_layer
      rst%num_fixatm = boundary%num_fixatm
      rst%fixatm     = boundary%fixatm
    end if

    if (qmmm%do_qmmm .and. qmmm%is_qm_charge) then
      call alloc_rst(rst, RestartQMMM, qmmm%qm_natoms)
      rst%qm_charge = qmmm%qm_charge
    end if

    call output_rst(output%rstfile, rst)

    if (dynvars%step == minimize%nsteps) &
      call dealloc_rst_all(rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, dynvars%coord)

    return

  end subroutine output_restart_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_remd
  !> @brief        output restart data for REMD
  !! @authors      TM
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[inout] remd     : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_remd(output, molecule, dynamics, dynvars, &
                                 boundary, remd)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_remd),            intent(inout) :: remd

    ! local variables
    type(s_rst)              :: rst
    integer                  :: k, natom, nrep, ndim


    call random_push_stock
    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    if (.not. replica_main_rank) &
      return

    natom = molecule%num_atoms
    call alloc_rst(rst, RestartAtom, natom)

    nrep = remd%total_nreplicas
    ndim = remd%dimension
    call alloc_rst(rst, RestartReplica, nrep, ndim)


    ! Output RST file
    !
    rst%iseed                  = dynamics%iseed
    rst%num_atoms              = molecule%num_atoms
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars %thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)
    rst%coord(1:3,1:natom)     = dynvars %coord(1:3,1:natom)
    rst%velocity(1:3,1:natom)  = dynvars %velocity(1:3,1:natom)

    if (remd%equilibration_only) then
      rst%rstfile_type                     = RstfileTypeMd
    else
      rst%rstfile_type                     = RstfileTypeRemd
      rst%iseed_remd                       = remd%iseed
      rst%nreplicas                        = nrep
      rst%dimension                        = ndim
      rst%repid2parmsetid(1:nrep)          = remd%repid2parmsetid(1:nrep)
      rst%num_criteria (1:nrep,1:ndim,1:2) = &
                             remd%num_criteria (1:nrep,1:ndim,1:2)
      rst%num_exchanges(1:nrep,1:ndim,1:2) = &
                             remd%num_exchanges(1:nrep,1:ndim,1:2)
    end if

    rst%sph_pot = boundary%sph_pot
    if (boundary%sph_pot) then
      call alloc_rst(rst, RestartSpot, boundary%nfunctions, natom)
      rst%nfunctions = boundary%nfunctions
      rst%radius     = boundary%radius
      rst%center     = boundary%center
      rst%const      = boundary%const
      rst%exponent   = boundary%exponent
      rst%fix_layer  = boundary%fix_layer
      rst%num_fixatm = boundary%num_fixatm
      rst%fixatm     = boundary%fixatm
    end if

    call output_rst(output%rstfile, rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, &
                              dynvars%coord)


    if (dynvars%step == dynamics%nsteps) &
        call dealloc_rst_all(rst)

    return

  end subroutine output_restart_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_rpath
  !> @brief        output restart data for RPATH
  !! @authors      TM, YA, KY
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[in]    qmmm     : QMMM variables information
  !! @param[inout] rpath    : RPATH information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_rpath(output, molecule, dynamics, dynvars, &
                                  boundary, qmmm, rpath)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_qmmm),            intent(in)    :: qmmm
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    type(s_rst)              :: rst
    integer                  :: k, natom, nrep, ndim


    call random_push_stock
    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    if (.not. replica_main_rank) &
      return

    natom = molecule%num_atoms
    call alloc_rst(rst, RestartAtom, natom)

    nrep = rpath%nreplica
    ndim = rpath%dimension
    call alloc_rst(rst, RestartRpath, ndim, nrep)

    ! Output RST file
    !
    if (rpath%rpathmode == RpathmodeMFEP) then

      rst%iseed                  = dynamics%iseed
      rst%num_atoms              = molecule%num_atoms
      rst%box_size_x             = boundary%box_size_x
      rst%box_size_y             = boundary%box_size_y
      rst%box_size_z             = boundary%box_size_z
      rst%thermostat_momentum    = dynvars %thermostat_momentum
      rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)
      rst%coord(1:3,1:natom)     = dynvars %coord(1:3,1:natom)
      rst%velocity(1:3,1:natom)  = dynvars %velocity(1:3,1:natom)

      rst%rstfile_type           = RstfileTypeRpath
      rst%nreplicas              = nrep
      rst%dimension              = ndim
      rst%rest_reference(1:2,1:ndim,1:nrep) =  &
        rpath%rest_reference(1:2,1:ndim,1:nrep)

    else if (rpath%rpathmode == RpathmodeFEP) then

      rst%rstfile_type           = RstfileTypeMd
      rst%iseed                  = dynamics%iseed
      rst%num_atoms              = molecule%num_atoms
      rst%box_size_x             = boundary%box_size_x
      rst%box_size_y             = boundary%box_size_y
      rst%box_size_z             = boundary%box_size_z
      rst%thermostat_momentum    = dynvars %thermostat_momentum
      rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)
      rst%coord(1:3,1:natom)     = dynvars %coord(1:3,1:natom)
      rst%velocity(1:3,1:natom)  = dynvars %velocity(1:3,1:natom)

    else if (rpath%rpathmode == RpathmodeMEPMD) then

      rst%rstfile_type       = RstfileTypeMd
      rst%iseed              = dynamics%iseed
      rst%num_atoms          = molecule%num_atoms
      rst%box_size_x         = boundary%box_size_x
      rst%box_size_y         = boundary%box_size_y
      rst%box_size_z         = boundary%box_size_z
      rst%thermostat_momentum    = dynvars %thermostat_momentum
      rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)
      rst%coord(1:3,1:natom)     = dynvars %coord(1:3,1:natom)
      rst%velocity(1:3,1:natom)  = dynvars %velocity(1:3,1:natom)

    else if (rpath%rpathmode == RpathmodeMEP) then

      rst%rstfile_type       = RstfileTypeMIN
      rst%num_atoms          = molecule%num_atoms
      rst%box_size_x         = boundary%box_size_x
      rst%box_size_y         = boundary%box_size_y
      rst%box_size_z         = boundary%box_size_z
      rst%energy             = dynvars%energy%total
      rst%delta_r            = 0.0D+00
      rst%coord(1:3,1:natom) = dynvars%coord(1:3,1:natom)

    else

      call error_msg('output_restart_path> unrecognized rpathmode')

    end if

    rst%sph_pot = boundary%sph_pot
    if (boundary%sph_pot) then
      call alloc_rst(rst, RestartSpot, boundary%nfunctions, natom)
      rst%nfunctions = boundary%nfunctions
      rst%radius     = boundary%radius
      rst%center     = boundary%center
      rst%const      = boundary%const
      rst%exponent   = boundary%exponent
      rst%fix_layer  = boundary%fix_layer
      rst%num_fixatm = boundary%num_fixatm
      rst%fixatm     = boundary%fixatm
    end if

    if (qmmm%do_qmmm .and. qmmm%is_qm_charge) then
      call alloc_rst(rst, RestartQMMM, qmmm%qm_natoms)
      rst%qm_charge = qmmm%qm_charge
    end if

    call output_rst(output%rstfile, rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, dynvars%coord)

    if (dynvars%step == dynamics%nsteps) call dealloc_rst_all(rst)

    return

  end subroutine output_restart_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_morph
  !> @brief        output restart data for MOPRH
  !! @authors      CK
  !! @param[in]    output   : output information
  !! @param[in]    molecule : molecule information
  !! @param[in]    morph    : morph information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_morph(output, molecule, morph, dynvars, boundary)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_morph),           intent(in)    :: morph
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    type(s_rst)              :: rst
    integer                  :: k, natom


    natom = molecule%num_atoms

    call alloc_rst(rst, RestartAtom, natom)

    rst%rstfile_type       = RstfileTypeMin
    rst%num_atoms          = molecule%num_atoms
    rst%box_size_x         = boundary%box_size_x
    rst%box_size_y         = boundary%box_size_y
    rst%box_size_z         = boundary%box_size_z
    rst%energy             = dynvars%energy%total
    rst%delta_r            = morph%delta_r
    rst%coord(1:3,1:natom) = dynvars%coord_ref(1:3,1:natom)

    call output_rst(output%rstfile, rst)

    if (dynvars%step == morph%ncycles) &
      call dealloc_rst_all(rst)

    ! Output Restart PDB file
    !
    if (output%pdbout) &
      call output_restart_pdb(output%pdbfile, molecule, dynvars%coord)

    return

  end subroutine output_restart_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_pdb 
  !> @brief        a driver subroutine for writing PDB file
  !! @authors      CK
  !! @param[in]    pdb_filename : filename of PDB file
  !! @param[in]    molecule     : structure of molecule
  !! @param[in]    coord        : coordinate
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_pdb(pdb_filename, molecule, coord)

    ! formal arguments
    character(*),            intent(in)    :: pdb_filename
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)

    ! local variables
    integer                  :: file
    

    ! open PDB file
    !
    call open_file(file, pdb_filename, IOFileOutputReplace)

    ! write coordinate data from MD
    !
    call write_restart_pdb(file, molecule, coord)

    ! close PDB file
    !
    call close_file(file)

    return

  end subroutine output_restart_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_restart_pdb
  !> @brief        write restart PDB file
  !! @authors      CK
  !! @param[in]    file     : unit number of PDB file
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_restart_pdb(file, molecule, coord)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)

    ! local variables
    integer                  :: i, j, len
    character(80)            :: fmt_a, fmt_t
    character(6)             :: crec
    character(4)             :: catm, cres, cstr, cseg
    logical                  :: use_cid
    logical                  :: res_col6, res_col5
    logical                  :: atom_col7

    
    res_col6  = .false.
    res_col5  = .false.
    atom_col7 = .false.

    if (molecule%residue_no(molecule%num_atoms) >= 100000) then
      res_col6 = .true.
    else if (molecule%residue_no(molecule%num_atoms) >= 10000) then
      res_col5 = .true.
    end if

    if (molecule%atom_no(molecule%num_atoms) >= 1000000) then
      atom_col7 = .true.
    end if

    if (atom_col7) then
      if (res_col6) then
        use_cid = .false.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,i6,48x,a1)'
      else if (res_col5) then
        use_cid = .false.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a4,i7,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,i7,6x,a4,a1,i4,49x,a1)'
      end if
    else
      if (res_col6) then
        use_cid = .false.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,i6,48x,a1)'
      else if (res_col5) then
        use_cid = .false.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a6,i5,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,i5,6x,a4,a1,i4,49x,a1)'
      end if
    end if

    do i = 1, molecule%num_atoms
      crec = 'ATOM  '

      read(molecule%atom_name(i), *) cstr
      len = len_trim(cstr)
      if (len < 4) then
        write(catm, fmt='(1x,a3)') cstr
      else
        catm = molecule%atom_name(i)
      end if

      read(molecule%residue_name(i),*) cstr
      len = len_trim(cstr)
      if (len == 2) then
        write(cres, fmt='(1x,a2)') cstr
      else if (len == 1) then
        write(cres, fmt='(2x,a1)') cstr
      else
        cres = cstr
      end if

      cseg = molecule%segment_name(i)

      if (use_cid) then
        write(file, fmt=fmt_a)              &
              crec,                         &
              molecule%atom_no(i),          &
              catm,                         &
              cres,                         &
              molecule%chain_id(i),         &
              molecule%residue_no(i),       &
              (coord(j,i), j=1,3),          &
              molecule%atom_occupancy(i),   &
              molecule%atom_temp_factor(i), &
              cseg
      else
        write(file, fmt=fmt_a)              &
              crec,                         &
              molecule%atom_no(i),          &
              catm,                         &
              cres,                         &
              molecule%residue_no(i),       &
              (coord(j,i), j=1,3),          &
              molecule%atom_occupancy(i),   &
              molecule%atom_temp_factor(i), &
              cseg
      end if

    end do

    write(file, fmt='(a6,69x,a1)') 'END   ', ' '

    return 

  end subroutine write_restart_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcd
  !> @brief        write coordinates in DCD format      
  !! @authors      RH
  !! @param[inout] output    : output_information
  !! iparam[in]    molecule  : molecule information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    istep     : current step number
  !! @param[in]    nstep     : total # of steps
  !! @param[in]    outperiod : coordinate output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    coord     : coordinates at t + dt
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcd(output, molecule, boundary, istep, &
                                  nstep, outperiod, timestep, coord)

    ! parameters
    character(4),            parameter     :: HEADER  = 'CORD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wp),                parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    real(wp),                intent(in)    :: timestep
    real(wp),                intent(in)    :: coord(:,:)

    ! local variable
    real(dp)                 :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: i, file, natom
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date


    ! If REMD, output is done on replica_main_rank
    ! If MD, output is done on main_rank
    !
    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if


    file    = output%dcdunit

    ! write header
    !
    if (output%out_dcd_header) then

      ts_namd = real(timestep * 1000_wp / TIMEFAC, sp)
      natom   = molecule%num_atoms

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      if (boundary%type == BoundaryTypePBC) &
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER,icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom

      output%out_dcd_header = .false.

    end if


    ! write coordinates
    !
    if (boundary%type == BoundaryTypePBC) then
      dcd_cell    = 0.0_dp
      dcd_cell(1) = boundary%box_size_x_ref
      dcd_cell(3) = boundary%box_size_y_ref
      dcd_cell(6) = boundary%box_size_z_ref
      write(file) dcd_cell(:)
    end if

    write(file) real(coord(1,:),sp)
    write(file) real(coord(2,:),sp)
    write(file) real(coord(3,:),sp)

    return

  end subroutine write_trajectory_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcdvel
  !> @brief        write velocities in DCD format
  !! @authors      RH
  !! @param[inout] output    : output information
  !! @param[in]    molecule  : molecule information
  !! @param[in]    boundary  : boundary conditions information
  !! @param[in]    istep     : current step number
  !! @param[in]    nstep     : total # of steps
  !! @param[in]    outperiod : velocity output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    velocity  : velocities at t + dt
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcdvel(output, molecule, boundary, istep, &
                                     nstep, outperiod, timestep, velocity)

    ! parameters
    character(4),            parameter     :: HEADER  = 'VELD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wp),                parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    real(wp),                intent(in)    :: timestep
    real(wp),                intent(in)    :: velocity(:,:)

    ! local variable
    real(dp)                 :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: i, file, natom
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date


    ! If REMD, output is done on replica_main_rank
    ! If MD, output is done on main_rank
    !
    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if


    file    = output%dcdvelunit

    ! write header
    !
    if (output%out_dcdvel_header) then

      ts_namd = real(timestep * 1000_wp / TIMEFAC,sp)
      natom   = molecule%num_atoms

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      if (boundary%type == BoundaryTypePBC) &
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER,icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom

      output%out_dcdvel_header = .false.

    end if


    ! write velocities
    !
    if (boundary%type == BoundaryTypePBC) then
      dcd_cell    = 0.0_dp
      dcd_cell(1) = boundary%box_size_x_ref
      dcd_cell(3) = boundary%box_size_y_ref
      dcd_cell(6) = boundary%box_size_z_ref
      write(file) dcd_cell(:)
    end if

    write(file) real(velocity(1,:),sp)
    write(file) real(velocity(2,:),sp)
    write(file) real(velocity(3,:),sp)

    return

  end subroutine write_trajectory_dcdvel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    include_id_to_filename
  !> @brief        include id to filename
  !! @authors      TM
  !! @param[inout] filename : replicate filename
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine include_id_to_filename(filename)

    ! formal arguments
    character(MaxFilename),  intent(inout) :: filename

    ! local variables
    integer                  :: comp, ndigit, id
    integer                  :: i, j, ci1, ci2, cnumber
    character(MaxFilename)   :: filename_ori
    character(10)            :: frmt, cd, cid


    ! get replica id
    !
    id = my_country_no + 1
    if (nrep_per_proc > 1) id = my_replica_no
    do i = 1, 100
      comp = 10**i
      if (id < comp) then
        ndigit = i
        exit
      end if
    end do

    ! check filename
    !
    filename_ori = filename

    ci1 = scan(filename, '{')
    ci2 = scan(filename, '}')

    if (ci1 == 0 .or. ci2 ==0) then
      call error_msg('Include_Id_To_Filename> {} is not found in [OUTPUT] file')
    end if

    if (ci1 > 0 .and. ci2 > ci1) then

      write(cd,'(i10)') ndigit
      frmt = '(i' // trim(adjustl(cd)) // '.' // trim(adjustl(cd)) // ')'
      write(cid,frmt) id

      cnumber = len_trim(filename_ori)
      if (cnumber + ndigit > MaxFilename) &
         call error_msg('Error: too long filename'//filename_ori)

      j = 0
      do i = 1, ci1 - 1
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do
      do i = 1, ndigit
        j = j + 1
        filename(j:j) = cid(i:i)
      end do
      do i = ci2+1, MaxFilename
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do

    end if

    return

  end subroutine include_id_to_filename

end module at_output_mod
