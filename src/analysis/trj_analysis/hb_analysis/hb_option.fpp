!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hb_option_mod
!> @brief   module for hydrogen bond analysis options
!! @authors Daisuke Matsuoka (DM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hb_option_mod

  use hb_option_str_mod
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
    logical            :: check_only     = .false.
    logical            :: allow_backup   = .false.
    integer            :: output_type    = HBOutputModeCountSnap
    real(wp)           :: hb_distance    =   3.4_wp
    real(wp)           :: dha_angle      = 120.0_wp
    real(wp)           :: hda_angle      =  30.0_wp
    integer            :: analysis_atom  = 1
    integer            :: target_atom    = 1
    integer            :: boundary_Type  = BoundaryTypePBC
    character(MaxLine) :: solvent_list   = ''
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
    write(MsgOut,'(A)') 'check_only     = NO           # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO           # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'analysis_atom  = 2            # atom group for searching HB partners'
    write(MsgOut,'(A)') 'target_atom    = 1            # search the HB partners from this group'
    write(MsgOut,'(A)') 'hb_distance    =   3.4        # the upper limit of (D .. A) distance (default: 3.4 A)'
    write(MsgOut,'(A)') 'DHA_angle      = 120.0        # the lower limit of (D-H .. A) angle (default: 120 deg.)'
    write(MsgOut,'(A)') 'HDA_angle      =  30.0        # the upper limit of (H-D .. A) angle (default:  30 deg.)'
    write(MsgOut,'(A)') 'boundary_type  = PBC          # (PBC / NOBC)'
    write(MsgOut,'(A)') 'output_type    = count_snap   # (count_atom / count_snap)'
    write(MsgOut,'(A)') '                              #   count_atom: number of H-bonds is output for each pair'
    write(MsgOut,'(A)') '                              #   count_snap: number of H-bonds is output every snapshot'
    write(MsgOut,'(A)') '# solvent_list = TIP3 POPC    # molecule names treated as solvent (only for count_atom)'
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

    call read_ctrlfile_type   (handle, Section, &
                              'output_type', opt_info%output_type, HBOutputType)

    call read_ctrlfile_integer(handle, Section, &
                              'analysis_atom', opt_info%analysis_atom)

    call read_ctrlfile_integer(handle, Section, &
                              'target_atom', opt_info%target_atom)

    call read_ctrlfile_real   (handle, Section, &
                              'hb_distance', opt_info%hb_distance)

    call read_ctrlfile_real   (handle, Section, &
                              'dha_angle', opt_info%dha_angle)

    call read_ctrlfile_real   (handle, Section, &
                              'hda_angle', opt_info%hda_angle)

    call read_ctrlfile_type   (handle, Section, 'boundary_type',    &
                               opt_info%boundary_type, BoundaryTypeTypes)

    call read_ctrlfile_string (handle, Section, 'solvent_list',   &
                               opt_info%solvent_list)

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

    write(MsgOut,'(A20,A)')    '  solvent list    = ', trim(opt_info%solvent_list)

    write(MsgOut,'(A20,A8)')   '  output type     = ', HBOutputType(opt_info%output_type)

    write(MsgOut,'(A20,A5,I0)') &
                    '  analysis atom   = ', 'group', opt_info%analysis_atom

    write(MsgOut,'(A20,A5,I0)') &
                    '  target atom     = ', 'group', opt_info%target_atom

    write(MsgOut,'(A20,F7.3)') &
                    '  HB distance     = ', opt_info%hb_distance

    write(MsgOut,'(A20,F7.3)') &
                    '  D-H .. A angle  = ', opt_info%dha_angle

    write(MsgOut,'(A20,F7.3)') &
                    '  H-D .. A angle  = ', opt_info%hda_angle

    write(MsgOut,'(A20,A10)')     &
                    '  boundary type   = ', BoundaryTypeTypes(opt_info%boundary_type)

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
    integer :: i, ndata


    ! check only
    !
    option%check_only = opt_info%check_only

    ! backup output file
    !
    call setup_backup(opt_info%allow_backup)

    ! solvent molecule list
    !
    if (allocated(option%solvent_list))  &
      deallocate(option%solvent_list)

    if (opt_info%solvent_list == '' .or.  &
        index(opt_info%solvent_list, 'NONE') >= 1) then
      ndata = 0

    else
      ndata = split_num(opt_info%solvent_list)
      allocate(option%solvent_list(ndata))

      call split(ndata, ndata, opt_info%solvent_list, option%solvent_list)
    end if

    ! output type
    !
    option%output_type = opt_info%output_type

    ! pbc correct
    !
    option%boundary_type = opt_info%boundary_type

    ! analysis atom
    !
    if (opt_info%analysis_atom > size(sel_info%groups)) &
        call error_msg('Setup_Option> analysis atom selection is out of range.')

    option%analysis_atom_exp = sel_info%groups(opt_info%analysis_atom)

    call select_atom(molecule, &
                     option%analysis_atom_exp, &
                     option%analysis_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> analysis atom count: ', &
         size(option%analysis_atom%idx)

    ! HB target atom
    !
    option%target_atom_exp = sel_info%groups(opt_info%target_atom)

    call select_atom(molecule, &
                     option%target_atom_exp, &
                     option%target_atom)

    write(MsgOut,'(A,I8)') 'Setup_Option> target atom count: ', &
         size(option%target_atom%idx)

    ! criteria of hydrogen bonding
    !
    option%hb_distance = opt_info%hb_distance
    option%dha_angle   = opt_info%dha_angle
    option%hda_angle   = opt_info%hda_angle

    write(MsgOut,'(A)') ' '

    return

  end subroutine setup_option

end module hb_option_mod
