!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_option_mod
!> @brief   module for analysis options
!! @authors Motoshi Kamiya (MK)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_option_mod

  use ea_option_str_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structure
  type, public :: s_opt_info
    logical                :: check_only     = .false.
    logical                :: allow_backup   = .false.
    logical                :: component      = .false.
    logical                :: rest_component = .false.
    character(MaxFilename) :: remfile        = '' ! hidden parameter
    integer                :: tgt_parmid     = 1  ! hidden parameter
    logical                :: mbar           = .false. ! for MBAR, hidden parameter
  end type s_opt_info

  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      MK
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_option

    if (main_rank) then
      write(MsgOut,'(A)') '[OPTION]'
      write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
      write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
      write(MsgOut,'(A)') 'component      = NO              # (YES/NO)'
!      write(MsgOut,'(A)') 'mbar           = NO              # for MBAR analysis'
    end if

    return

  end subroutine

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      MK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, opt_info)

    ! parameters
    character(*), parameter :: Section = 'Option'

    ! formal arguments
    integer,          intent(in)    :: handle
    type(s_opt_info), intent(inout) :: opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                               'check_only', opt_info%check_only)
    call read_ctrlfile_logical(handle, Section, &
                               'allow_backup', opt_info%allow_backup)
    call read_ctrlfile_logical(handle, Section, &
                               'component', opt_info%component)
    call read_ctrlfile_logical(handle, Section, &
                               'rest_component', opt_info%rest_component)
    call read_ctrlfile_logical(handle, Section, &
                               'mbar', opt_info%mbar)
    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options'
      if (opt_info%check_only) then
        write(MsgOut,'(A20,A3)')     '  check only      = ', 'yes'
      else
        write(MsgOut,'(A20,A3)')     '  check only      = ', 'no'
      end if
      if (opt_info%allow_backup) then
        write(MsgOut,'(A20,A3)')     '  allow backup    = ', 'yes'
      else
        write(MsgOut,'(A20,A2)')     '  allow backup    = ', 'no'
      end if
      if (opt_info%component) then
        write(MsgOut,'(A20,A3)')     '  component       = ', 'yes'
      else
        write(MsgOut,'(A20,A3)')     '  component       = ', 'no'
      end if
      if (opt_info%rest_component) then
        call error_msg('Read_Ctrl_Option> rest_component is not available. Please use analysis_grest in spdyn.')
      end if
      if (opt_info%mbar) then
        write(MsgOut,'(A20,A3)')     '  mbar            = ', 'yes'
      else
        write(MsgOut,'(A20,A3)')     '  mbar            = ', 'no'
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        read control parameters in OPTION section
  !! @authors      MK
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_option(opt_info, molecule, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option


    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! component (energy component)
    !
    option%component = opt_info%component

    ! rest component (rest component; Eu-u, Eu-v, Ev-v)
    ! it is not available (20220615)
    !
    option%rest_component = opt_info%rest_component

    ! remfile name
    !
    option%remfile = opt_info%remfile

    ! target id
    !
    option%tgt_parmid = opt_info%tgt_parmid

    ! mbar
    !
    option%mbar = opt_info%mbar

    return

  end subroutine setup_option

end module ea_option_mod
