!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ca_option_mod
!> @brief   module for analysis options
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ca_option_mod

  use ca_option_str_mod
  use pbc_correct_mod
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
    integer                         :: analysis_atom  = 0
    integer                         :: analysis_for   = AnalysisForAll
    integer                         :: output_coord   = OutputCoordXYZ
    logical                         :: output_atomno  = .true.
    integer                         :: pbc_correct    = PBCCModeNo
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
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_option

    write(MsgOut,'(A)') '[OPTION]'
    write(MsgOut,'(A)') 'check_only     = NO       # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO       # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'analysis_atom  = 2        # atom group'
    write(MsgOut,'(A)') 'analysis_for   = ALL      # analyze COM of (ALL/MOLECULE) in the selected group'
    write(MsgOut,'(A)') '                          # ALL     : Treat the selected group as a single molecule'
    write(MsgOut,'(A)') '                          # MOLECULE: Analyze each molecule in the selected group'
    write(MsgOut,'(A)') 'output_coord   = XYZ      # (XYZ/XY)'
    write(MsgOut,'(A)') 'pbc_correct    = NO       # (NO/MOLECULE) (MOLECULE requres psf/prmtop/grotop)'
    write(MsgOut,'(A)') ' '

    return

  end subroutine show_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      TM
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
                              'analysis_atom', opt_info%analysis_atom)

    call read_ctrlfile_type   (handle, Section, 'analysis_for', &
                               opt_info%analysis_for, AnalysisForTypes)

    call read_ctrlfile_type   (handle, Section, 'output_coord', &
                               opt_info%output_coord, OutputCoordTypes)

    call read_ctrlfile_logical(handle, Section, &
                              'output_atomno', opt_info%output_Atomno)

    call read_ctrlfile_type   (handle, Section, 'pbc_correct',  &
                               opt_info%pbc_correct, PBCCModeTypes)

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
                    '  analysis atom   = ', 'group', opt_info%analysis_atom

    write(MsgOut,'(A20,A10)')   &
                    '  analysis for    = ', AnalysisForTypes(opt_info%analysis_for)

    write(MsgOut,'(A20,A10)')   &
                    '  output coord    = ', OutputCoordTypes(opt_info%output_coord)

    if (opt_info%output_atomno) then
      write(MsgOut,'(A20,A3)') '  output atomno   = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  output atomno   = ', 'no'
    end if

    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      TM
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


    ! check only
    option%check_only = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    ! analysis atom
    if (opt_info%analysis_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp = sel_info%groups(opt_info%analysis_atom)

    call select_atom(molecule, &
                     option%analysis_atom_exp, &
                     option%analysis_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count : ', &
         size(option%analysis_atom%idx)

    option%analysis_for  = opt_info%analysis_for

    option%output_coord = opt_info%output_coord

    option%output_atomno = opt_info%output_atomno

    ! pbc correct
    option%pbcc_mode = opt_info%pbc_correct

    ! setup for PBC-correct mode "molecule"
    call setup_pbc_correct(option%pbcc_mode, molecule)

    write(MsgOut,'(A)') ''

    return

  end subroutine setup_option

end module ca_option_mod
