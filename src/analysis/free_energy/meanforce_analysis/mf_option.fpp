!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   mf_option_mod
!> @brief   module for analysis options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module mf_option_mod

  use mf_option_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only      = .false.
    logical                         :: allow_backup    = .false.
    integer                         :: cv_atom         = 2
    integer                         :: nimage          = 100
    real(wp)                        :: force_constant  = 100

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

    write(MsgOut,'(A)') '[MEANFORCE]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'cv_atom        = 1'
    write(MsgOut,'(A)') 'nimage         = 100'
    write(MsgOut,'(A)') 'force_constant = 100'
    write(MsgOut,'(A)') ''

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in MEANFORCE section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: SectionMeanforce = 'MEANFORCE'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info


    ! read control parameters
    !


    ! MEANFORCE section
    call begin_ctrlfile_section(handle, SectionMeanforce)

    call read_ctrlfile_logical(handle, SectionMeanforce, &
                               'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, SectionMeanforce, &
                               'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, SectionMeanforce, &
                               'cv_atom', opt_info%cv_atom)

    call read_ctrlfile_integer(handle, SectionMeanforce, &
                               'nimage', opt_info%nimage)

    call read_ctrlfile_real(handle, SectionMeanforce, &
                               'force_constant', opt_info%force_constant)

    call end_ctrlfile_section(handle)


    ! write parameter to MsgOut
    !

    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of MEANFORCE'
    write(MsgOut,'(A)') ''
    if (opt_info%check_only) then
      write(MsgOut,'(A20,A3)') '  check only      = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  check only      = ', 'no'
    end if
    if (opt_info%allow_backup) then
      write(MsgOut,'(A20,A3)') '  allow backup    = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  allow backup    = ', 'no'
    end if
    write(MsgOut,'(A20,I0)')   '  cv atoms        = ', opt_info%cv_atom
    write(MsgOut,'(A20,I0)')   '  # of images     = ', opt_info%nimage
    write(MsgOut,'(A20,F10.3)')   '  force constant  = ', opt_info%force_constant
    write(MsgOut,'(A)') ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, sel_info, molecule, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option


    option%check_only = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! cv atoms
    if (opt_info%cv_atom > size(sel_info%groups)) &
      call error_msg('Setup_Option> cv atom seleciton is out of range.')

    call select_atom(molecule, &
                     sel_info%groups(opt_info%cv_atom), &
                     option%cv_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> CV atom count: ', &
             size(option%cv_atom%idx)
    write(MsgOut,'(A)') ' '
 
    ! nimage
    option%nimage = opt_info%nimage

    ! force constant
    option%force_constant = opt_info%force_constant

    return

  end subroutine setup_option
  
end module mf_option_mod
