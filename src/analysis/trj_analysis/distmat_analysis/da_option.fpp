!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_option_mod
!> @brief   module for analysis options
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_option_mod

  use da_option_str_mod
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
    integer                         :: analysis_atom_1 = 1
    integer                         :: analysis_atom_2 = 2
    integer                         :: matrix_shape    = MatrixShapeHalf
    integer                         :: calc_mode       = intramolecule
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
    write(MsgOut,'(A)') 'check_only      = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup    = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'calc_mode       = intra           # (intra/inter)'
    write(MsgOut,'(A)') 'analysis_atom   = 1               # atom group'
    write(MsgOut,'(A)') '# analysis_atom_2 = 2             # for the inter-molecule distance'
    write(MsgOut,'(A)') 'matrix_shape    = HALF            # (HALF/FULL)'
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

    call read_ctrlfile_type   (handle, Section, 'calc_mode',       &
                               opt_info%calc_mode, calculation_mode)

    call read_ctrlfile_type   (handle, Section, 'matrix_shape',    &
                               opt_info%matrix_shape, MatrixShapeTypes)

    call read_ctrlfile_integer(handle, Section, &
                              'analysis_atom', opt_info%analysis_atom_1)

    call read_ctrlfile_integer(handle, Section, &
                              'analysis_atom_2', opt_info%analysis_atom_2)

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

    if (opt_info%calc_mode == intramolecule) then
      write(MsgOut,'(A20,A23)')   &
                   '  calc_mode       = ', 'intra-molecule distance'

      write(MsgOut,'(A20,A5,I0)') &
                   '  analysis atom   = ', 'group', opt_info%analysis_atom_1

      if (opt_info%matrix_shape == MatrixShapeHalf) then
        write(MsgOut,'(A20,A4)') '  matrix shape    = ', 'half'

      else if (opt_info%matrix_shape == MatrixShapeFull) then
        write(MsgOut,'(A20,A4)') '  matrix shape    = ', 'full'

      end if

    else if (opt_info%calc_mode == intermolecule) then
      write(MsgOut,'(A20,A23)')   &
                   '  calc_mode       = ', 'inter-molecule distance'

      write(MsgOut,'(A20,A5,I0)') &
                   '  analysis atom 1 = ', 'group', opt_info%analysis_atom_1

      write(MsgOut,'(A20,A5,I0)') &
                   '  analysis atom 2 = ', 'group', opt_info%analysis_atom_2

      write(MsgOut,'(A20,A4)') '  matrix shape    = ', 'full'

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
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(inout) :: option


    ! check only
    !
    option%check_only   = opt_info%check_only
    option%matrix_shape = opt_info%matrix_shape
    option%calc_mode    = opt_info%calc_mode

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! analysis atom
    !
    if (opt_info%analysis_atom_1 > size(sel_info%groups)) &
        call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp_1 = sel_info%groups(opt_info%analysis_atom_1)

    call select_atom(molecule, &
                     option%analysis_atom_exp_1, &
                     option%analysis_atom_1)

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count: ', &
         size(option%analysis_atom_1%idx)

    if (opt_info%calc_mode == intermolecule) then
      option%matrix_shape = MatrixShapeFull

      if (opt_info%analysis_atom_2 > size(sel_info%groups)) &
          call error_msg('Setup_Option> analysis atom selection is out of range.')

      option%analysis_atom_exp_2 = sel_info%groups(opt_info%analysis_atom_2)

      call select_atom(molecule, &
                       option%analysis_atom_exp_2, &
                       option%analysis_atom_2)

      write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count: ', &
           size(option%analysis_atom_2%idx)

    end if

    return

  end subroutine setup_option

end module da_option_mod
