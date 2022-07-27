!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pn_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pn_option_mod

  use pn_option_str_mod
  use fileio_control_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only      = .false.
    logical                         :: trajectory      = .false.
    integer                         :: nreplicas       = 1
    integer                         :: nfiles          = 1
  end type s_opt_info

  ! subroutines
  public  :: show_ctrl_option
  public  :: read_ctrl_option
  public  :: setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'nreplica       = 1'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'OPTION'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info


    ! read control parameters
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                               'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                               'trajectory', opt_info%trajectory)

    call read_ctrlfile_integer(handle, Section, &
                               'nreplica', opt_info%nreplicas)

    opt_info%nfiles = opt_info%nreplicas
    call read_ctrlfile_integer(handle, Section, &
                               'nfiles', opt_info%nfiles)

    call end_ctrlfile_section(handle)


    ! write parameter to MsgOut
    !

    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of OPTION'
    write(MsgOut,'(A)') ''

    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if

    if (opt_info%trajectory) then
      write(MsgOut,'(A20,A3)') '  trajectory      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  trajectory      = ', 'no'
    end if

    write(MsgOut,'(A20,I10)') &
         '     nreplica     = ', opt_info%nreplicas

    if (opt_info%nfiles .ne. opt_info%nreplicas) then
      write(MsgOut,'(A20,I10)') &
         '     nfiles       = ', opt_info%nfiles
    endif

    write(MsgOut,'(A)') ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_option),          intent(inout) :: option


    option%check_only = opt_info%check_only
    option%trajectory = opt_info%trajectory
    option%nreplicas  = opt_info%nreplicas
    option%nfiles     = opt_info%nfiles

    return

  end subroutine setup_option
  
end module pn_option_mod
