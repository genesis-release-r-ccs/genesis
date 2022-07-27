!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fa_option_mod
!> @brief   module for analysis options
!! @authors Daisuke Matsuoka (DM), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fa_option_mod

  use fa_option_str_mod
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
    logical    :: check_only      = .false.
    logical    :: allow_backup    = .false.
    integer    :: analysis_atom_1 = 1
    integer    :: analysis_atom_2 = 2
    real(wp)   :: Forster_radius  = 10.0_wp
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
  !! @authors      DM, NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option


    write(MsgOut,'(A)') '[SELECTION]'
    write(MsgOut,'(A)') 'group1           = resno:1         # atom group 1'
    write(MsgOut,'(A)') 'group2           = resno:3         # atom group 2'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only       = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup     = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'Forster_radius   = 10.0            # R_0'
    write(MsgOut,'(A)') 'analysis_atom_1  = 1               # group number'
    write(MsgOut,'(A)') 'analysis_atom_2  = 2               # group number'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DM, NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   opt_info : OPTION section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_option(handle, opt_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info

    ! local variables
    integer           :: i
    character(len=20) :: name


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                              'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                              'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, 'analysis_atom_1', opt_info%analysis_atom_1)
    call read_ctrlfile_integer(handle, Section, 'analysis_atom_2', opt_info%analysis_atom_2)
    call read_ctrlfile_real   (handle, Section, 'Forster_radius',  opt_info%Forster_radius)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options'
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

    write(MsgOut,'(A20, F7.3)')   &
                   '  Forster radius  = ', opt_info%Forster_radius

    write(MsgOut,'(A20,a5,i0)')   &
                   '  analysis_atom_1 = ', 'group', opt_info%analysis_atom_1

    write(MsgOut,'(A20,a5,i0)')   &
                   '  analysis_atom_2 = ', 'group', opt_info%analysis_atom_2

    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DM, NT
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
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option

    ! local variables
    integer           :: i
    character(len=20) :: name


    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! Forster radius
    !
    option%Forster_radius = opt_info%Forster_radius

    ! analysis atom
    !
    if (opt_info%analysis_atom_1 > size(sel_info%groups)) &
      call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp(1) = sel_info%groups(opt_info%analysis_atom_1)

    call select_atom(molecule,  &
                     option%analysis_atom_exp(1), &
                     option%analysis_atom(1))

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis_atom count: ', &
                            size(option%analysis_atom(1)%idx)

    if (opt_info%analysis_atom_2 > size(sel_info%groups)) &
      call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp(2) = sel_info%groups(opt_info%analysis_atom_2)

    call select_atom(molecule,  &
                     option%analysis_atom_exp(2), &
                     option%analysis_atom(2))

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis_atom count: ', &
                            size(option%analysis_atom(2)%idx)

    write(MsgOut,'(A)') ''

    return

  end subroutine setup_option

end module fa_option_mod
