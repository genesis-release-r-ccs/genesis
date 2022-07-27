!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_option_mod
!> @brief   module for tilt analysis options
!! @authors Daisuke Matsuoka (DM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_option_mod

  use ta_option_str_mod
  use trajectory_str_mod
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
    logical                         :: check_only     = .false.
    logical                         :: allow_backup   = .false.
    integer                         :: TM_helix_atom  = 1
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
  !! @authors      DM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'TM_helix_atom  = 1               # transmembrane helix atoms'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DM
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


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, &
                               'check_only', opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, &
                               'allow_backup', opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, &
                               'TM_helix_atom', opt_info%TM_helix_atom)
 
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
    write(MsgOut,'(A20,A5,I0)') &
                    '  TM_helix atom   = ', 'group', opt_info%TM_helix_atom
 
    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DM
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[out]   option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, &
                          sel_info, &
                          molecule, &
                          option)
  
    ! formal argments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option


    ! check only
    option%check_only    = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! TM helix atoms
    if (opt_info%TM_helix_atom > size(sel_info%groups)) &
      call error_msg('Setup_Option> TM_helix atom selection is out of range.')

    option%TM_helix_atom_exp = sel_info%groups(opt_info%TM_helix_atom)

    call select_atom(molecule, &
                     option%TM_helix_atom_exp, &
                     option%TM_helix_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> TM helix atom count: ', &
          size(option%TM_helix_atom%idx)

    write(MsgOut,'(A)') ' '

    return

  end subroutine setup_option

end module ta_option_mod
