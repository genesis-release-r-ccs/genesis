!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_output_md
!> @brief   definition of output
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Yasuhiro Matsunaga (YM),
!!          Ryuhei Harada (RH), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_output_mod

  use sp_dynvars_mod
  use sp_parallel_io_mod
  use sp_output_str_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_remd_str_mod
  use sp_rpath_str_mod
  use sp_constraints_mod
  use sp_alchemy_str_mod
  use sp_grest_energy_mod
  use molecules_str_mod
  use fileio_control_mod
  use fileio_rst_mod
  use fileio_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_out_info
    character(MaxFilename) :: logfile    = ''
    character(MaxFilename) :: dcdfile    = ''
    character(MaxFilename) :: enefile    = ''
    character(MaxFilename) :: dcdvelfile = ''
    character(MaxFilename) :: selfile    = ''
    character(MaxFilename) :: rstfile    = ''
    character(MaxFilename) :: pdbfile    = ''
    character(MaxFilename) :: remfile    = ''
    character(MaxFilename) :: rpathfile  = ''
    character(MaxFilename) :: mfrcfile   = ''
    character(MaxFilename) :: rpathlogfile  = ''
    character(MaxFilename) :: gamdfile   = ''
    character(MaxFilename) :: fepfile    = ''
  end type s_out_info

  ! valiables
  real(wip),        private, allocatable :: tmp_coord1(:,:)
  real(wip),        private, allocatable :: tmp_coord2(:,:)

  ! subroutines
  public  :: show_ctrl_output
  public  :: read_ctrl_output
  public  :: setup_output_md
  public  :: setup_output_min
  public  :: setup_output_remd
  public  :: setup_output_rpath
  public  :: open_output
  public  :: close_output
  public  :: output_md
  public  :: output_min
  public  :: output_remd
  public  :: output_rpath
  public  :: output_prst_md
  public  :: output_prst_min
  public  :: output_gamd

  private :: output_restart_md
  private :: output_restart_min
  private :: output_restart_remd
  private :: output_restart_rpath
  private :: write_trajectory_dcd
  private :: write_trajectory_dcdvel
  private :: reduce_coord
  private :: include_id_to_filename
  ! FEP
  public  :: output_fep_energy
  private :: match_fep_coord

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
        write(MsgOut,'(A)') '# dcdfile    = sample.dcd    # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample.dcd    # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample.rst    # restart file'
        write(MsgOut,'(A)') '# rstfile    = sample().rst  # parallel I/O restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample.pdb    # PDB file'
        write(MsgOut,'(A)') '# gamdfile   = sample.gamd   # gamd file'
        write(MsgOut,'(A)') '# fepfile    = sample.fepout # fep file'
        write(MsgOut,'(A)') ' '

      case ('min')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') '# dcdfile    = sample.dcd   # DCD trajectory file'
        write(MsgOut,'(A)') '# rstfile    = sample.rst   # restart file'
        write(MsgOut,'(A)') '# rstfile    = sample().rst # parallel I/O restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample.pdb   # PDB file'
        write(MsgOut,'(A)') ' '

      case ('remd')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'logfile    = sample{}.log # log file of each replica'
        write(MsgOut,'(A)') '# dcdfile    = sample{}.dcd  # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample{}.dcd  # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample{}.rst  # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample{}.pdb  # PDB file'
        write(MsgOut,'(A)') '# remfile    = sample{}.rem  # replica exchange ID file'
        write(MsgOut,'(A)') '# enefile    = sample{}.ene  # energy output for grest (only when analysis_grest = YES'
        write(MsgOut,'(A)') '# fepfile    = sample.fepout # fep file'
        write(MsgOut,'(A)') ' '

      case ('rpath')

        write(MsgOut,'(A)') '[OUTPUT]'
        write(MsgOut,'(A)') 'logfile    = sample{}.log # log file of each replica'
        write(MsgOut,'(A)') '# dcdfile    = sample{}.dcd # DCD trajectory file'
        write(MsgOut,'(A)') '# dcdvelfile = sample{}.dcd # DCD velocity file'
        write(MsgOut,'(A)') '# rstfile    = sample{}.rst # restart file'
        write(MsgOut,'(A)') '# pdbfile    = sample{}.pdb # PDB file'
        write(MsgOut,'(A)') '# rpathfile  = sample{}.rpath # replica path ID file'
        write(MsgOut,'(A)') '# mfrcfile   = sample{}.mfrc # mean force file'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

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
    call read_ctrlfile_string(handle, Section, 'enefile',   out_info%enefile)
    call read_ctrlfile_string(handle, Section, 'dcdvelfile',out_info%dcdvelfile)
    call read_ctrlfile_string(handle, Section, 'selfile',   out_info%selfile)
    call read_ctrlfile_string(handle, Section, 'rstfile',   out_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'pdbfile',   out_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'remfile',   out_info%remfile)
    call read_ctrlfile_string(handle, Section, 'rpathfile', out_info%rpathfile)
    call read_ctrlfile_string(handle, Section, 'rpathlogfile', out_info%rpathlogfile)
    call read_ctrlfile_string(handle, Section, 'mfrcfile', out_info%mfrcfile)
    call read_ctrlfile_string(handle, Section, 'gamdfile',  out_info%gamdfile)
    call read_ctrlfile_string(handle, Section, 'fepfile',   out_info%fepfile)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Output> Output Files'
      if (out_info%logfile .ne. '') &
        write(MsgOut,*) ' logfile    = ', trim(out_info%logfile) 
      if (out_info%dcdfile .ne. '') &
        write(MsgOut,*) ' dcdfile    = ', trim(out_info%dcdfile)
      if (out_info%enefile .ne. '') &
        write(MsgOut,*) ' enefile    = ', trim(out_info%enefile)
      if (out_info%dcdvelfile .ne. '') &
        write(MsgOut,*) ' dcdvelfile = ', trim(out_info%dcdvelfile)
      if (out_info%selfile .ne. '') &
        write(MsgOut,*) ' selfile    = ', trim(out_info%selfile)
      if (out_info%rstfile .ne. '') &
        write(MsgOut,*) ' rstfile    = ', trim(out_info%rstfile)
      if (out_info%pdbfile .ne. '') &
        write(MsgOut,*) ' pdbfile    = ', trim(out_info%pdbfile)
      if (out_info%remfile .ne. '') &
        write(MsgOut,*) ' remfile    = ', trim(out_info%remfile)
      if (out_info%rpathfile .ne. '') &
        write(MsgOut,*) ' rpathfile  = ', trim(out_info%rpathfile)
      if (out_info%rpathlogfile .ne. '') &
        write(MsgOut,*) ' rpathlogfile  = ', trim(out_info%rpathlogfile)
      if (out_info%mfrcfile .ne. '') &
        write(MsgOut,*) ' mfrcfile   = ', trim(out_info%mfrcfile)
      if (out_info%gamdfile .ne. '') &
        write(MsgOut,*) ' gamdfile   = ', trim(out_info%gamdfile)
      if (out_info%fepfile /= '') &
        write(MsgOut,*) ' fepfile    = ', trim(out_info%fepfile)
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
      if (out_info%dcdfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: dcdfile name is not'//&
                  ' specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      output%dcdout  = .true.
      output%dcdfile = out_info%dcdfile
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: dcdvelfile name is not'//&
                  ' specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      output%dcdvelout  = .true.
      output%dcdvelfile = out_info%dcdvelfile
    end if

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Md> Error: rstfile name is not'//&
                       'specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      output%rstout  = .true.
      output%rstfile = out_info%rstfile
    end if

    if (out_info%pdbfile .ne. '') then
      output%pdbout  = .true.
      output%pdbfile = out_info%pdbfile
    end if

    if (out_info%gamdfile .ne. '') then
      output%gamdout  = .true.
      output%gamdfile = out_info%gamdfile
    end if

    if (dynamics%fepout_period > 0) then
      if (out_info%fepfile == '') &
        call error_msg('Setup_Output_Md> Error: fepfile name is not'//&
                  ' specified in [OUTPUT] (fepout_period > 0 in [ALCHEMY])')
      output%fepout  = .true.
      output%fepfile = out_info%fepfile
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
    output%logout     = .false.
    output%dcdvelout  = .false.

    if (minimize%crdout_period > 0) then
      if (out_info%dcdfile .eq. '') &
        call error_msg('Setup_Output_Min> Error: dcdfile name is not'//&
                  ' specified in [OUTPUT] (crdout_period > 0 in [MINIMIZE])')
      output%dcdout  = .true.
      output%dcdfile = out_info%dcdfile
    end if

    if (minimize%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Min> Error: rstfile name is not'//&
                  ' specified in [OUTPUT] (rstout_period > 0 in [MINIMIZE])')
      output%rstout  = .true.
      output%rstfile = out_info%rstfile
    end if

    if (out_info%pdbfile .ne. '') then
      output%pdbout  = .true.
      output%pdbfile = out_info%pdbfile
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
  !! @param[in]    remd     : remd information
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
        call error_msg('Setup_Output_Remd> Error: dcdfile name is not'//&
                  ' specified in [OUTPUT] (crdout_period > 0 in [DYNAMICS])')
      output%dcdfile = out_info%dcdfile
      call include_id_to_filename(output%dcdfile)
      output%dcdfile = output%dcdfile
      output%dcdout  = .true.
      if (remd%analysis_grest) then
        if (out_info%enefile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: enefile name is not'//&
                       ' specified in [OUTPUT] (crdout_period > 0 in' //&
                       ' [DYNAMICS]) and analysis_grest = YES')
        output%enefile = out_info%enefile
        call include_id_to_filename(output%enefile)
        output%enefile = output%enefile
        output%eneout  = .true.
      end if
    end if

    if (dynamics%velout_period > 0) then
      if (out_info%dcdvelfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: dcdvelfile name is not'//&
                  ' specified in [OUTPUT] (velout_period > 0 in [DYNAMICS])')
      output%dcdvelfile = out_info%dcdvelfile
      call include_id_to_filename(output%dcdvelfile)
      output%dcdvelfile = output%dcdvelfile
      output%dcdvelout  = .true.
    end if


    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Remd> Error: rstfile name is not'//&
                  ' specified in [OUTPUT] (rstout_period > 0 in [DYNAMICS])')
      output%rstfile = out_info%rstfile
      call include_id_to_filename(output%rstfile)
      output%rstfile = output%rstfile
      output%rstout  = .true.

      if (.not. remd%equilibration_only .and. remd%exchange_period > 0) then
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

    if (dynamics%eneout_period > 0) then
      if (mod(remd%exchange_period,dynamics%eneout_period) /= 0) &
      call error_msg('Setup_Output_Remd> Error in eneout_period'//&
      '  mod(exchange_period, eneout_period) is not ZERO')
    end if

    if (out_info%remfile .ne. '') then
      output%remfile = out_info%remfile
      call include_id_to_filename(output%remfile)
      output%remfile = output%remfile
      output%remout  = .true.
    end if

    if (dynamics%fepout_period > 0) then
      if (out_info%fepfile == '') &
        call error_msg('Setup_Output_Remd> Error: fepfile name is not'//&
                  ' specified in [OUTPUT] (fepout_period > 0 in [ALCHEMY])')
      output%fepfile = out_info%fepfile
      call include_id_to_filename(output%fepfile)
      output%fepfile = output%fepfile
      output%fepout  = .true.

      ! output remfile to sort fepfile
      if (out_info%remfile == '') &
        call error_msg('Setup_Output_Remd> Error: remfile name is not'//&
                  ' specified in [OUTPUT]')
      output%remfile = out_info%remfile
      call include_id_to_filename(output%remfile)
      output%remfile = output%remfile
      output%remout  = .true.
    end if

    output%replica = .true.

    return

  end subroutine setup_output_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output_rpath
  !> @brief        setup output information for RPATH
  !! @authors      TM
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    rpath    : information of rpath
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output_rpath(out_info, dynamics, rpath, output)

    ! formal arguments
    type(s_out_info),        intent(in)    :: out_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_rpath),           intent(in)    :: rpath     
    type(s_output),          intent(inout) :: output


    output%verbose = dynamics%verbose

    output%logfile = out_info%logfile
    call include_id_to_filename(output%logfile)
    output%logout  = .true.

    if (dynamics%crdout_period > 0) then
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

    if (dynamics%rstout_period > 0) then
      if (out_info%rstfile .eq. '') &
        call error_msg('Setup_Output_Rpath> rst filename is blank')
      output%rstfile = out_info%rstfile
      call include_id_to_filename(output%rstfile)
      output%rstout  = .true.

      if (.not. rpath%equilibration_only .and. rpath%rpath_period > 0) then
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

    if (out_info%rpathfile .ne. '') then
      output%rpathfile = out_info%rpathfile
      call include_id_to_filename(output%rpathfile)
      output%rpathout  = .true.
    end if

    if (dynamics%crdout_period > 0) then
      if (out_info%mfrcfile .ne. '') then
        output%mfrcfile = out_info%mfrcfile
        call include_id_to_filename(output%mfrcfile)
        output%mfrcout  = .true.
      end if
    end if

    if (dynamics%eneout_period > 0) then
      if (mod(rpath%rpath_period,dynamics%eneout_period) /= 0) &
      call error_msg('Setup_Output_Rpath> Error in eneout_period'//&
      '  mod(rpath_period, eneout_period) is not ZERO')
    end if

    if (out_info%rpathlogfile .ne. '') then
      output%rpathlogfile = out_info%rpathlogfile
      output%rpathlogout  = .true.
    end if

    output%rpath = .true.

    return

  end subroutine setup_output_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    open_output
  !> @brief        open output file
  !! @authors      TM
  !! @param[inout] outout : information of output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine open_output(output)

    ! formal arguments
    type(s_output),          intent(inout) :: output

    ! local variables
    integer                  :: file


    ! open logfile 
    !
    if (output%logout) then
      if (replica_main_rank) then
        call open_file(file, output%logfile, IOFileOutputNew)
        output%logunit = file
        DynvarsOut = output%logunit
      end if
    end if


    ! open dcdfile
    !
    if (output%dcdout) then

      output%dcd_para = pio_check_ranked_file(output%dcdfile)

      if (output%dcd_para) then

        call pio_open_file(file, pio_get_ranked_filename(output%dcdfile,   &
                                                         my_country_rank), &
                           IOFileOutputNew)
        output%dcdunit = file

      else if (main_rank .or. replica_main_rank) then

        call open_binary_file(file, output%dcdfile, IOFileOutputNew)
        output%dcdunit = file

      end if

    end if

    ! open dcdvelfile
    !
    if (output%dcdvelout) then

      output%dcdvel_para = pio_check_ranked_file(output%dcdvelfile)

      if (output%dcdvel_para) then

        call pio_open_file(file, pio_get_ranked_filename(output%dcdvelfile, &
                                                         my_country_rank),  &
                           IOFileOutputNew)
        output%dcdvelunit = file

      else if (main_rank .or. replica_main_rank) then

        call open_binary_file(file, output%dcdvelfile, IOFileOutputNew)
        output%dcdvelunit = file

      end if

    end if

    ! open rstfile
    !   Just check whether rstfile is already existing.
    !
    if (output%rstout) then

      output%rst_para = pio_check_ranked_file(output%rstfile)
                           
      if (output%rst_para) then

        output%rstfile = pio_get_ranked_filename(output%rstfile, &
                                                 my_country_rank)
        call open_binary_file(file, output%rstfile, IOFileOutputNew)
        call close_file(file)

      else if (main_rank .or. replica_main_rank) then

        call open_file(file, output%rstfile, IOFileOutputNew)
        call close_file(file)

      end if

    end if

    ! open pdbfile
    !   Just check whether pdbfile is already existing.
    !
    if (output%pdbout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%pdbfile, IOFileOutputNew)
        call close_file(file)
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

    ! open enefile
    !
    if (output%eneout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%enefile, IOFileOutputNew)
        output%eneunit = file
      end if
    end if

    ! open rpathfile
    !
    if (output%rpathout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%rpathfile, IOFileOutputNew)
        output%rpathunit = file
      end if
    end if

    ! open rpathlogfile (only main rank)
    !
    if (output%rpathlogout) then
      if (main_rank) then
        call open_file(file, output%rpathlogfile, IOFileOutputNew)
        output%rpathlogunit = file
      end if
    end if

    ! open mfrcfile
    !
    if (output%mfrcout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%mfrcfile, IOFileOutputNew)
        output%mfrcunit = file
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

    ! open fepfile
    !
    if (output%fepout) then
      if (main_rank .or. replica_main_rank) then
        call open_file(file, output%fepfile, IOFileOutputNew)
        output%fepunit = file
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

    if (output%eneout) &
      call close_file(output%eneunit)

    ! close dcdfile
    !
    if (output%dcdout) &
      call close_file(output%dcdunit)

    ! close rpathlogfile
    !
    if (output%rpathout) &
      call close_file(output%rpathunit)

    ! close rpathlogfile
    !
    if (output%rpathlogout) &
      call close_file(output%rpathlogunit)

    ! close dcdvelfile
    !
    if (output%dcdvelout) &
      call close_file(output%dcdvelunit)

    ! close logfile
    !
    if (output%logout) &
      call close_file(output%logunit)


    return

  end subroutine close_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_md
  !> @brief        output trajectory and restart data for MD
  !! @authors      JJ, TM, NT
  !! @param[inout] output      : output information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] remd        : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_md(output, dynamics, boundary, pairlist, ensemble, &
                       constraints, dynvars, domain, enefunc, remd)

    ! formakl arguments
    type(s_output),          intent(inout) :: output
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(in)    :: constraints
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),            intent(inout) :: remd  

    ! local variables
    logical                   :: alloc_ref, dealloc_ref
    real(wip)                 :: scale_v, gamma_t, dt
    real(wip), allocatable    :: tmp_coord(:,:,:)
    integer                   :: i, ix, istep

    integer,          pointer :: ncell, natom(:)


    ncell        => domain%num_cell_local
    natom        => domain%num_atom

    gamma_t      =  ensemble%gamma_t*AKMA_PS
    dt           =  dynamics%timestep/AKMA_PS
    istep        =  dynvars%step

    ! Output coordinates (Wrap_all is not done)
    !
    if (dynamics%crdout_period > 0 .and. istep > 0) then

      if (mod(istep,dynamics%crdout_period) == 0) then

        if (.not. allocated(tmp_coord)) &
          allocate(tmp_coord(3, MaxAtom, ncell))

        tmp_coord(1:3,1:MaxAtom,1:ncell) = 0.0_wip

        if (dynamics%integrator == IntegratorLEAP) then
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%coord_ref(1:3,ix,i)
            end do
          end do
        else
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%coord(1:3,ix,i)
            end do
          end do
        end if

        if (output%dcd_para) then

          call pio_write_domain_trj(output%dcdunit,   &
                                    boundary,         &
                                    domain,           &
                                    tmp_coord,        &
                                int(dynamics%nsteps/dynamics%crdout_period),&
                                    .true.)

        else

          call write_trajectory_dcd(output,            &
                                    boundary,          &
                                    domain,            &
                                    istep,             &
                                    dynamics%nsteps,   &
                                    dynamics%crdout_period, &
                                    dynamics%timestep, &
                                    tmp_coord)

        end if

        if (remd%analysis_grest) then
          alloc_ref = (istep == dynamics%crdout_period)
          dealloc_ref = (istep == dynamics%nsteps)
          if (dynamics%integrator == IntegratorLEAP) then
            call compute_grest_energy_output(1, alloc_ref, dealloc_ref, &
                                             boundary, pairlist,     &
                                             ensemble, domain,       &
                                             enefunc, remd)
          else
            call compute_grest_energy_output(2, alloc_ref, dealloc_ref, &
                                             boundary, pairlist,     &
                                             ensemble, domain,       &
                                             enefunc, remd)
          end if
          call write_grest_energy_output(output, alloc_ref, istep,   &
                                         dynamics%nsteps, remd)
        end if

      end if
    end if


    ! Output velocities (Wrap_all is not done)
    !
    if (dynamics%velout_period > 0 .and. istep > 0) then
      if (mod(istep,dynamics%velout_period) == 0) then
        if (.not. allocated(tmp_coord)) &
          allocate(tmp_coord(3, MaxAtom, ncell))

        tmp_coord(1:3,1:MaxAtom,1:ncell) = 0.0_wip

        if (dynamics%integrator == IntegratorLEAP) then
          if (ensemble%tpcontrol == TpcontrolLangevin) then
            scale_v = sqrt(1.0_wip+0.5_wip*gamma_t*dt)
          else
            scale_v = 1.0_wip
          end if

          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = 0.5_wip * &
                                    (domain%velocity_ref(1:3,ix,i) + &
                                     domain%velocity    (1:3,ix,i))*scale_v
            end do
          end do
        else
          do i = 1, ncell
            do ix = 1, natom(i)
              tmp_coord(1:3,ix,i) = domain%velocity(1:3,ix,i)
            end do
          end do
        end if

        if (output%dcdvel_para) then

          call pio_write_domain_trj(output%dcdvelunit,   &
                                    boundary,            &
                                    domain,              &
                                    tmp_coord, &
                                int(dynamics%nsteps/dynamics%crdout_period),&
                                    .false.)

        else

          call write_trajectory_dcdvel(output,            &
                                       boundary,          &
                                       domain,            &
                                       istep,             &
                                       dynamics%nsteps,   &
                                       dynamics%velout_period, &
                                       dynamics%timestep, &
                                       tmp_coord)
        end if

      end if
    end if


    ! Output restart data
    !
    if (dynamics%rstout_period > 0 .and. istep > 0) then
      if (mod(istep,dynamics%rstout_period) == 0) then

        if (output%rst_para) then
          call output_prst_md(output, enefunc, dynamics, boundary, &
                              dynvars, domain, constraints)
        else
          call output_restart_md(output, domain, dynamics, boundary, dynvars)
        end if

      end if
    end if

    return

  end subroutine output_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_min
  !> @brief        output trajectory and restart data for minimization
  !! @authors      TM
  !! @param[inout] output      : output information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    minimize    : minimize information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    dynvars     : dynamic variables information
  !! @param[in]    delta_r     : delta r
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_min(output, enefunc, minimize, boundary, dynvars, delta_r, &
                        constraints, domain)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_minimize),        intent(in)    :: minimize
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars
    real(wip),               intent(in)    :: delta_r
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: istep


    istep        =  dynvars%step

    ! output energy
    !
    if (minimize%eneout_period > 0) then
      if (mod(istep,minimize%eneout_period) == 0) then
        call output_dynvars(output, enefunc, dynvars)
      end if
    end if

    ! Output coordinates
    !
    if (minimize%crdout_period > 0) then

      if (mod(istep,minimize%crdout_period) == 0) then

        if (output%dcd_para) then

          call pio_write_domain_trj(output%dcdunit,   &
                                    boundary,         &
                                    domain,           &
                                    domain%coord,     &
                                int(minimize%nsteps/minimize%crdout_period),&
                                    .true.)

        else

          call write_trajectory_dcd(output,          &
                                    boundary,        &
                                    domain,          &
                                    istep,           &
                                    minimize%nsteps, &
                                    minimize%crdout_period, &
                                    0.001_wip,        &
                                    domain%coord)
        end if
 
      end if
    end if


    ! output restart data
    !
    if (minimize%rstout_period > 0) then
      if (mod(istep,minimize%rstout_period) == 0) then
        if (output%rst_para) then
          call output_prst_min(output, enefunc, minimize, boundary, dynvars, &
                               domain, constraints)
        else
          call output_restart_min(output, domain, minimize, boundary, dynvars, &
                                  delta_r)
        end if
      end if
    end if

    return

  end subroutine output_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_remd
  !> @brief        output restart data for REMD
  !! @authors      TM, NT
  !! @param[in]    icycle      : number of remd cycle
  !! @param[inout] output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] domain      : domain information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    dynvars     : dynamic variables information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    constraints : constraints information
  !! @param[inout] remd        : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_remd(icycle, output, enefunc, domain, dynamics, dynvars, &
                         boundary, constraints, remd)

    ! formal arguments
    integer,                 intent(in)    :: icycle
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_domain),          intent(inout) :: domain
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_constraints),     intent(in)    :: constraints
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
      if (dynvars%step > 0 .and.                                    &
          mod(dynvars%step,dynamics%rstout_period) == 0) then
        if (output%rst_para) then
          call output_prst_remd(output, enefunc, dynamics, boundary,  &
                                dynvars, domain, constraints, remd)
        else
          call output_restart_remd(output, domain, dynamics, dynvars, &
                                   boundary, remd)
        end if
      end if
    end if

    return

  end subroutine output_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_rpath
  !> @brief        output restart data for RPATH
  !! @authors      NT
  !! @param[inout] output   : output information
  !! @param[inout] domain   : domain information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    rpath    : RPATH information
  !! @param[in]    enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_rpath(output, domain, dynamics, dynvars, boundary, rpath, &
                                                                       enefunc)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),          intent(in)    :: domain
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_boundary),        intent(in)    :: boundary
    type(s_rpath),   target, intent(inout) :: rpath
    type(s_enefunc),         intent(in)    :: enefunc

    ! local variables
    integer                  :: k, repid_i

    if (main_rank .and. output%rpathlogout) then
      write(output%rpathlogunit, '(I10,3F15.8)') dynvars%step, rpath%sum_distance, &
                                 rpath%distance_init,rpath%distance_prev
    endif

    ! wrire collective variables per replica
    !
    repid_i = my_country_no + 1

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
        call output_restart_rpath(output, domain, dynamics, dynvars,  &
                                  boundary, rpath)
      end if
    end if

    ! wrire mean forces per replica
    !
    if (replica_main_rank .and. output%mfrcout) then
      write(output%mfrcunit, '(I10,$)') dynvars%step
      do k = 1, rpath%dimension
        write(output%mfrcunit, '(E15.5,$)') enefunc%stats_force_save(k)
      end do
      write(output%mfrcunit, *)
    end if

    return

  end subroutine output_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_gamd
  !> @brief        output GaMD boost parameters
  !! @authors      HO
  !! @param[inout] output   : output information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_gamd(output, dynvars, enefunc)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_dynvars),         intent(in)    :: dynvars
    type(s_enefunc),target,  intent(in)    :: enefunc

    ! local variables
    integer,parameter            :: clength=16, flength=4
    integer                      :: i, ifm
    character(16)                :: title
    character(16)                :: category(999)
    character                    :: frmt*5, rfrmt*7
    character                    :: rfrmt_cont*9, frmt_cont*7
    real(wp)                     :: values(999)
    logical, save                :: gamdtitle = .true.
    type(s_enefunc_gamd),pointer :: gamd


    gamd => enefunc%gamd

    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'


    if (.not. main_rank) return

    if (output%gamdout) then

        ifm = 1

        if ((gamd%boost_pot) .or. &
            (gamd%boost_dual)) then
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

        if ((gamd%boost_dih) .or. &
            (gamd%boost_dual)) then
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
                if (i .eq. ifm-1) then
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
            if (i .eq. ifm-1) then
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
  !  Subroutine    output_prst_md
  !> @brief        output Parallel I/O restart data for MD
  !! @authors      NT
  !! @param[in]    output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    dynvars     : dynamic variables information
  !! @param[inout] domain      : domain information
  !! @param[in]    constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_prst_md(output, enefunc, dynamics, boundary, dynvars, &
                            domain, constraints)

    ! formal arguments
    type(s_output),          intent(in)   :: output
    type(s_enefunc),         intent(inout):: enefunc
    type(s_dynamics),        intent(in)   :: dynamics
    type(s_boundary),        intent(in)   :: boundary
    type(s_dynvars),         intent(in)   :: dynvars
    type(s_domain),          intent(inout):: domain
    type(s_constraints),     intent(in)   :: constraints

    ! local variables
    integer                  :: istep


    if (dynamics%rstout_period > 0 .and. output%rst_para) then

      if (mod(dynvars%step, dynamics%rstout_period) == 0) then

        pio_restart = .true.

        if (dynamics%hydrogen_mr) &
        call setup_mass_repartitioning_back(dynamics, constraints, domain, &
                                            enefunc)
 
        call pio_write_domain_rst(output%rstfile, &
                                  boundary,       &
                                  domain,         &
                                  enefunc,        &
                                  constraints,    &
                                  dynvars,        &
                                  dynamics)

        if (dynamics%hydrogen_mr) &
        call setup_mass_repartitioning(dynamics, constraints, domain, enefunc)

      end if

    end if

    return

  end subroutine output_prst_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_prst_min
  !> @brief        output Parallel I/O restart data for minimization
  !! @authors      NT
  !! @param[in]    output      : output information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    minimize    : minimize information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    dynvars     : dynamic variables information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_prst_min(output, enefunc, minimize, boundary, dynvars, &
                             domain, constraints)

    ! formal arguments
    type(s_output),          intent(in)   :: output
    type(s_enefunc),         intent(in)   :: enefunc
    type(s_minimize),        intent(in)   :: minimize
    type(s_boundary),        intent(in)   :: boundary
    type(s_dynvars),         intent(in)   :: dynvars
    type(s_domain),          intent(in)   :: domain
    type(s_constraints),     intent(in)   :: constraints

    ! local variables
    integer                  :: istep


    istep = dynvars%step

    if (minimize%rstout_period > 0 .and. output%rst_para) then

      if (mod(istep, minimize%rstout_period) == 0) then

        pio_restart = .false.

        call pio_write_domain_rst(output%rstfile, &
                                  boundary,       &
                                  domain,         &
                                  enefunc,        &
                                  constraints)

      end if

    end if

    return

  end subroutine output_prst_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_prst_remd
  !> @brief        output Parallel I/O restart data for MD
  !! @authors      NT
  !! @param[in]    output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    boundary    : boundary condition information
  !! @param[in]    dynvars     : dynamic variables information
  !! @param[inout] domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[in]    rem         : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_prst_remd(output, enefunc, dynamics, boundary, dynvars, &
                              domain, constraints, remd)

    ! formal arguments
    type(s_output),          intent(in)   :: output
    type(s_enefunc),         intent(inout):: enefunc
    type(s_dynamics),        intent(in)   :: dynamics
    type(s_boundary),        intent(in)   :: boundary
    type(s_dynvars),         intent(in)   :: dynvars
    type(s_domain),          intent(inout):: domain
    type(s_constraints),     intent(in)   :: constraints
    type(s_remd),            intent(in)   :: remd

    ! local variables
    integer                  :: istep


    if (dynamics%rstout_period > 0 .and. output%rst_para) then

      if (dynvars%step >0 .and. & 
          mod(dynvars%step,dynamics%rstout_period) == 0) then

        pio_restart = .true.
        if (dynamics%hydrogen_mr) &
        call setup_mass_repartitioning_back(dynamics, constraints, domain, &
                                            enefunc)

        call pio_write_domain_rst(output%rstfile, &
                                  boundary,       &
                                  domain,         &
                                  enefunc,        &
                                  constraints,    &
                                  dynvars,        &
                                  dynamics,       &
                                  remd)

        if (dynamics%hydrogen_mr) &
        call setup_mass_repartitioning(dynamics, constraints, domain, enefunc)

      end if

    end if

    return

  end subroutine output_prst_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_md
  !> @brief        output restart data for MD
  !! @authors      NT
  !! @param[in]    output   : output information
  !! @param[in]    domain   : domain information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_md(output, domain, dynamics, boundary, dynvars)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_domain),  target, intent(in)    :: domain
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars

    ! local variables
    type(s_rst)              :: rst
    integer                  :: i, ix, istep
    
    real(wip),       pointer :: coord(:,:,:), vel(:,:,:)
    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    coord     => domain%coord
    vel       => domain%velocity
    istep     =  dynvars%step

    if (output%replica) &
      return

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! setup restart information
    !
    call alloc_rst(rst, RestartAtom, int(natom_all))

    rst%rstfile_type = RstfileTypeMd
    rst%iseed                  = dynamics%iseed_init_velocity
    rst%num_atoms              = natom_all
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars%thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars%barostat_momentum(1:3)

    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (domain%fep_use) then
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    rst%coord(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! reduce velocities
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = vel(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))
    if (domain%fep_use) then
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    rst%velocity(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! output restart information
    !
    if (main_rank .or. replica_main_rank) &
      call output_rst(output%rstfile, rst)

    if (istep == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    return

  end subroutine output_restart_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_min
  !> @brief        output restart data for minimization
  !! @authors      NT
  !! @param[in]    output   : output information
  !! @param[in]    domain   : domain information
  !! @param[in]    minimize : minimization information
  !! @param[in]    boundary : boundary condition information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    delta_r  : delta-r
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_min(output, domain, minimize, boundary, dynvars, &
                                delta_r)

    ! formal arguments
    type(s_output),          intent(in)    :: output
    type(s_domain),  target, intent(in)    :: domain
    type(s_minimize),        intent(in)    :: minimize
    type(s_boundary),        intent(in)    :: boundary
    type(s_dynvars),         intent(in)    :: dynvars
    real(wip),               intent(in)    :: delta_r

    ! local variables
    type(s_rst)              :: rst
    integer                  :: i, ix, istep
    
    real(wip),       pointer :: coord(:,:,:)
    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    coord     => domain%coord
    istep     =  dynvars%step

    if (output%replica) &
      return

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! setup restart information
    !
    call alloc_rst(rst, RestartAtom, int(natom_all))

    rst%rstfile_type = RstfileTypeMin
    rst%num_atoms    = natom_all
    rst%box_size_x   = boundary%box_size_x
    rst%box_size_y   = boundary%box_size_y
    rst%box_size_z   = boundary%box_size_z
    rst%energy       = real(dynvars%energy%total,wp)
    rst%delta_r      = delta_r

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (domain%fep_use) then 
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    rst%coord(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! output restart information
    !
    if (main_rank) &
      call output_rst(output%rstfile, rst)

    if (istep == minimize%nsteps) &
      call dealloc_rst_all(rst)

    return

  end subroutine output_restart_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_remd
  !> @brief        output restart data for REMD
  !! @authors      TM
  !! @param[in]    output   : output information
  !! @param[in]    domain   : dmain information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[inout] remd     : REMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_remd(output, domain, dynamics, dynvars, &
                                 boundary, remd)

    ! formal arguments
    type(s_output),           intent(in)    :: output
    type(s_domain),   target, intent(in)    :: domain
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_dynvars),          intent(in)    :: dynvars
    type(s_boundary),         intent(in)    :: boundary
    type(s_remd),             intent(inout) :: remd

    ! local variables
    type(s_rst)              :: rst
    real(wip),       pointer :: coord(:,:,:), vel(:,:,:)
    integer                  :: i, ix, nrep, ndim
    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    coord     => domain%coord
    vel       => domain%velocity

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    call alloc_rst(rst, RestartAtom, int(natom_all))

    nrep = remd%total_nreplicas
    ndim = remd%dimension
    call alloc_rst(rst, RestartReplica, nrep, ndim)

    ! Output RST file
    !
    rst%iseed                  = dynamics%iseed_init_velocity
    rst%num_atoms              = natom_all
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars %thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)

    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))
    
    if (domain%fep_use) then
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    rst%coord(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! reduce velocities
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = vel(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (domain%fep_use) then
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    rst%velocity(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    if (remd%equilibration_only) then
      rst%rstfile_type                     = RstfileTypeMd
    else
      rst%rstfile_type                     = RstfileTypeRemd
      rst%iseed_remd                       = remd%iseed
      rst%nreplicas                        = nrep
      rst%dimension                        = ndim
      rst%repid2parmsetid(1:nrep)          = remd%repid2parmsetid(1:nrep)
      rst%num_criteria (1:nrep,1:ndim,1:2) = remd%num_criteria (1:nrep,1:ndim,1:2)
      rst%num_exchanges(1:nrep,1:ndim,1:2) = remd%num_exchanges(1:nrep,1:ndim,1:2)
    end if

    if (replica_main_rank) &
      call output_rst(output%rstfile, rst)

    if (dynvars%step == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    return

  end subroutine output_restart_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_restart_rpath
  !> @brief        output restart data for RPATH
  !! @authors      CK
  !! @param[in]    output   : output information
  !! @param[in]    domain   : dmain information
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dynvars  : dynamic variables information
  !! @param[in]    boundary : boundary information
  !! @param[inout] rpath    : rpath information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_restart_rpath(output, domain, dynamics, dynvars, &
                                 boundary, rpath)

    ! formal arguments
    type(s_output),           intent(in)    :: output
    type(s_domain),   target, intent(in)    :: domain
    type(s_dynamics),         intent(in)    :: dynamics
    type(s_dynvars),          intent(in)    :: dynvars
    type(s_boundary),         intent(in)    :: boundary
    type(s_rpath),            intent(inout) :: rpath

    ! local variables
    type(s_rst)              :: rst
    real(wip),       pointer :: coord(:,:,:), vel(:,:,:)
    integer                  :: i, ix, nrep, ndim
    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    coord     => domain%coord
    vel       => domain%velocity

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    call alloc_rst(rst, RestartAtom, int(natom_all))

    nrep = rpath%nreplica
    ndim = rpath%dimension
    call alloc_rst(rst, RestartRpath, ndim, nrep)

    ! Output RST file
    !
    rst%rstfile_type           = RstfileTypeRpath
    rst%iseed                  = dynamics%iseed_init_velocity
    rst%num_atoms              = natom_all
    rst%box_size_x             = boundary%box_size_x
    rst%box_size_y             = boundary%box_size_y
    rst%box_size_z             = boundary%box_size_z
    rst%thermostat_momentum    = dynvars %thermostat_momentum
    rst%barostat_momentum(1:3) = dynvars %barostat_momentum(1:3)

    call random_stock_mpi_tobyte(mpi_comm_country, rst%random)

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))
    rst%coord(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    ! reduce velocities
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = vel(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))
    rst%velocity(1:3,1:natom_all) = tmp_coord1(1:3,1:natom_all)

    rst%rstfile_type  = RstfileTypeRpath
    rst%nreplicas     = nrep
    rst%dimension     = ndim
    rst%rest_reference(1:2,1:ndim,1:nrep) =  &
      rpath%rest_reference(1:2,1:ndim,1:nrep)

    if (replica_main_rank) &
      call output_rst(output%rstfile, rst)

    if (dynvars%step == dynamics%nsteps) &
      call dealloc_rst_all(rst)

    return

  end subroutine output_restart_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcd
  !> @brief        write coordinates in DCD format with domain decomposition
  !! @authors      RH, TM, JJ
  !! @param[inout] output    : output information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    domain    : dmain information
  !! @param[in]    istep     : dynamics step
  !! @param[in]    nstep     : number of dynamics steps
  !! @param[in]    outperiod : coordinate output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    coord     : coordinates per domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcd(output, boundary, domain, istep, nstep, &
                                  outperiod, timestep, coord)

    ! parameters
    character(4),            parameter     :: HEADER  = 'CORD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wip),               parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),  target, intent(in)    :: domain
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    real(wip),               intent(in)    :: timestep
    real(wip),               intent(in)    :: coord(:,:,:)

    ! local variable
    real(dp)                 :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: file
    integer                  :: i, ix
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date

    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g
    file      =  output%dcdunit

    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = coord(1:3,ix,i)
      end do
    end do
    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

    if (domain%fep_use) then
      ! FEP: copy singleA atoms to singleB atoms
      call match_fep_coord(domain, tmp_coord1)
    end if

    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if
  

    ! write header
    !
    if (output%out_dcd_header) then

      ts_namd = real(timestep * 1000_wip / TIMEFAC, sp)

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER, icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom_all

      output%out_dcd_header = .false.

    end if

    ! write coordinates
    !
    dcd_cell(1:6) = 0.0_dp
    dcd_cell(1)   = boundary%box_size_x_ref
    dcd_cell(3)   = boundary%box_size_y_ref
    dcd_cell(6)   = boundary%box_size_z_ref
    write(file) dcd_cell(:)

    write(file) real(tmp_coord1(1,:),sp)
    write(file) real(tmp_coord1(2,:),sp)
    write(file) real(tmp_coord1(3,:),sp)

    return

  end subroutine write_trajectory_dcd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_grest_energy_output
  !> @brief        write grest energy output for mbar analysis
  !! @authors      JJ
  !! @param[inout] output    : unit number of crdfile
  !! @param[in]    initial   : initial step or not
  !! @param[in]    istep     : current step number
  !! @param[in]    nstep     : total # of steps
  !! @param[in]    remd      : remd information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_grest_energy_output(output, initial, istep, nstep, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    logical,                 intent(in)    :: initial
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    type(s_remd),            intent(in)    :: remd

    ! local variable
    real(dp)                 :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: file
    integer                  :: i, k
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date


    if (output%replica) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if

    if (initial) then
      do i = 1, remd%dimension
        if (remd%types(i) == RemdSoluteTempering) then
          write(output%eneunit,'("# Column   1: STEP")')
          do k = 1, remd%nreplicas(i)
            write(output%eneunit, '("# Column ",i3,": POTENTIAL ENERGY of", &
                  & i3,"-th parameter (kcal/mol)" )') k+1, k
          end do
          write(output%eneunit,'("# unit of output energy is kcal/mol")')
          write(output%eneunit,'("#")')
        end if
      end do
    end if

    write(output%eneunit, '(I16,2x)',advance="NO") istep
    do i = 1, remd%dimension
      if (remd%types(i) == RemdSoluteTempering) then
        do k = 1, remd%nreplicas(i)
          write(output%eneunit, '(F16.4)',advance="NO") &
                remd%grest_energy%potential_energy(k)
        end do
      end if
    end do
    write(output%eneunit,*)

    return

  end subroutine write_grest_energy_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_trajectory_dcdvel
  !> @brief        write velocities in DCD format with domain decomposition
  !! @authors      RH, TM, YM
  !! @param[inout] output    : output information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    domain    : dmain information
  !! @param[in]    istep     : dynamics step
  !! @param[in]    nstep     : number of dynamics steps
  !! @param[in]    outperiod : coordinate output period
  !! @param[in]    timestep  : time step in PS
  !! @param[in]    velocity  : velocities per domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_trajectory_dcdvel(output, boundary, domain, istep, nstep, &
                                     outperiod, timestep, velocity)

    ! parameters
    character(4),            parameter     :: HEADER  = 'VELD'
    integer(4),              parameter     :: NTITLE  = 4
    real(wip),               parameter     :: TIMEFAC = 48.88821

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: nstep
    integer,                 intent(in)    :: outperiod
    type(s_domain),  target, intent(in)    :: domain
    real(wip),               intent(in)    :: timestep
    real(wip),               intent(in)    :: velocity(:,:,:)

    ! local variable
    real(wip)                :: dcd_cell(6)
    real(sp)                 :: ts_namd
    integer                  :: file
    integer                  :: i, ix
    integer(4)               :: icntrl(20)
    character(80)            :: title(4)
    character(24)            :: name, date

    integer,         pointer :: ncell, natom(:), id_l2g(:,:)
    integer(iintegers), pointer :: natom_all


    ncell     => domain%num_cell_local
    natom_all => domain%num_atom_all
    natom     => domain%num_atom
    id_l2g    => domain%id_l2g


    if (.not. allocated(tmp_coord1)) &
      allocate(tmp_coord1(3, natom_all), tmp_coord2(3, natom_all))

    ! reduce coordinates
    !
    tmp_coord1(1:3,1:natom_all) = 0.0_wip
    do i = 1, ncell
      do ix = 1, natom(i)
        tmp_coord1(1:3,id_l2g(ix,i)) = velocity(1:3,ix,i)
      end do
    end do

    call reduce_coord(tmp_coord1, tmp_coord2, int(natom_all))

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

      ts_namd = real(timestep * 1000_wip / TIMEFAC,sp)

      icntrl(:)  = 0
      icntrl(1)  = nstep / outperiod            ! total # of frames   (NSET)
      icntrl(2)  = istep                        ! first time step     (NPRIV)
      icntrl(3)  = outperiod                    ! output period       (NSAVC)
      icntrl(4)  = nstep                        ! number of time step (NSTEP)
      icntrl(10) = transfer(ts_namd,icntrl(10)) ! length of time step (DELTA)
      icntrl(11) = 1                            ! flag for with unit-cell
      icntrl(20) = 24                           ! PRETEND TO BE CHARMM24 -JCP

      call fdate (date)
      call getlog(name)

      title(1) = 'REMARKS CREATED BY GENESIS                                                      '
      title(2) = 'REMARKS DATE: ' // date // ' CREATED BY USER: ' // name
      title(3) = '                                                                                '
      title(4) = '                                                                                '

      write(file) HEADER, icntrl(1:20)
      write(file) NTITLE,(title(i),i=1,NTITLE)
      write(file) natom_all

      output%out_dcdvel_header = .false.

    end if

    ! write coordinates
    !
    dcd_cell(1:6) = 0.0_wip
    dcd_cell(1)   = boundary%box_size_x_ref
    dcd_cell(3)   = boundary%box_size_y_ref
    dcd_cell(6)   = boundary%box_size_z_ref
    write(file) dcd_cell(:)

    write(file) real(tmp_coord1(1,:),sp)
    write(file) real(tmp_coord1(2,:),sp)
    write(file) real(tmp_coord1(3,:),sp)

    return

  end subroutine write_trajectory_dcdvel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_coord
  !> @brief        reduce coord
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_coord(coord, temporary, natom)

    ! formal arguments
    real(wip),               intent(inout) :: coord(:,:)
    real(wip),               intent(inout) :: temporary(:,:)
    integer,                 intent(in)    :: natom

#ifdef HAVE_MPI_GENESIS

    ! local variables
    integer                  :: j
    integer                  :: ncycle, icycle, nlen, ixx


    ! Reduce coordinate
    !
    temporary(1:3,1:natom) = 0.0_wip
    do j = 1, natom
      temporary(1,j) = coord(1,j)
      temporary(2,j) = coord(2,j)
      temporary(3,j) = coord(3,j)
    end do

    ncycle = (natom - 1) / mpi_drain + 1
    nlen = mpi_drain
    ixx  = 1

    do icycle = 1, ncycle
      if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
      call mpi_reduce(temporary(1,ixx), coord(1,ixx), 3*nlen,    &
                      mpi_wip_real, mpi_sum, 0, mpi_comm_country, &
                      ierror)
      ixx = ixx + nlen
    end do

#endif

    return

  end subroutine reduce_coord

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

    if (ci1 == 0 .or. ci2 == 0) then
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

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_fep_energy
  !> @brief        output energy difference between adjacent states in FEP
  !! @authors      HO
  !! @param[in]    output   : output information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    dynvars  : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_fep_energy(output, enefunc, dynvars)

    ! formal arguments
    type(s_output),             intent(in) :: output
    type(s_enefunc),            intent(in) :: enefunc
    type(s_dynvars),            intent(in) :: dynvars

    ! local variables
    integer,parameter        :: clength=16, flength=4
    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, frmt_res*10, rfrmt*7
    character                :: rfrmt_cont*9,frmt_cont*7
    real(dp)                 :: values(999)
    real(dp)                 :: ene_restraint
    logical, save            :: fep_title = .true.

    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if

    title = '#     STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    write(category(ifm),frmt) 'Total_E_ref'
    values(ifm) = dynvars%energy%deltU_fep(1)
    ifm = ifm+1

    if (enefunc%num_fep_neighbor == 2) then

      write(category(ifm),frmt) 'Delta_E_rev'
      values(ifm) = dynvars%energy%deltU_fep(2)
      ifm = ifm+1

      write(category(ifm),frmt) 'Delta_E_fwd'
      values(ifm) = dynvars%energy%deltU_fep(3)
      ifm = ifm+1

    else if (enefunc%num_fep_neighbor == 1) then

      if (enefunc%fep_direction == FEP_Forward) then

        write(category(ifm),frmt) 'Delta_E_fwd'
        values(ifm) = dynvars%energy%deltU_fep(2)
        ifm = ifm+1

      else if (enefunc%fep_direction == FEP_Reverse) then

        write(category(ifm),frmt) 'Delta_E_rev'
        values(ifm) = dynvars%energy%deltU_fep(2)
        ifm = ifm+1

      end if

    end if

    if (fep_title) then

      write(output%fepunit,'(A10,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(output%fepunit,frmt) category(i)
        else
          write(output%fepunit,frmt_cont) category(i)
        endif
      end do

      fep_title = .false.
    end if

    write(output%fepunit,'(I10,$)') dynvars%step

    do i = 1, ifm-1
      if (i == ifm-1) then
        write(output%fepunit,rfrmt) values(i)
      else
        write(output%fepunit,rfrmt_cont) values(i)
      endif
    end do

    return

  end subroutine output_fep_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    match_fep_single
  !> @brief        match up coordinates of single-topology parts
  !! @authors      HO
  !! @param[in]    domain : domain information
  !! @param[inout] coord  : coordinates to be matched
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine match_fep_coord(domain, coord)

    ! formal arguments
    type(s_domain), intent(in)    :: domain
    real(wip),      intent(inout) :: coord(:,:)

    ! local variables
    integer                       :: i, iatomA, iatomB 

    do i = 1, domain%num_atom_single_all
      iatomA = domain%id_singleA(i)
      iatomB = domain%id_singleB(i)
      coord(1:3,iatomB) = coord(1:3,iatomA)
    end do

    return

  end subroutine match_fep_coord

end module sp_output_mod
