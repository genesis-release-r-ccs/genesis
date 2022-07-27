!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hbond_option_mod
!> @brief   module for analysis options
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hbond_option_mod

  use hbond_option_str_mod
  use select_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif


  implicit none
  private

  ! structures
  type, public :: s_hbond_opt_info

    integer            :: recenter      = 1
    integer            :: analysis_atom = 1
    integer            :: target_atom   = 1
    integer            :: output_type   = HBoutputModeCountSnap
    real(wp)           :: hb_distance   = 3.4_wp
    real(wp)           :: dha_angle     = 120.0_wp
    real(wp)           :: hda_angle     =  30.0_wp
    character(MaxLine) :: solvent_list  = ''

  end type s_hbond_opt_info

  ! subroutines
  public  :: show_hbond_ctrl_option
  public  :: read_ctrl_option
  public  :: hbond_setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_hbond_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      DM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_hbond_ctrl_option

    write(MsgOut,'(A)') '[HBOND_OPTION]'
    write(MsgOut,'(A)') 'recenter      = 1                # recenter'
    write(MsgOut,'(A)') 'output_type   = Count_snap       # (Count_snap / Count_Atom)'
    write(MsgOut,'(A)') 'analysis_atom = 1                # atom group for searching H-bond partners'
    write(MsgOut,'(A)') 'target_atom   = 2                # serach the HB partners from this atom group'
    write(MsgOut,'(A)') 'solvent_list  = TIP3 POPC        # molecule names treated as solvent (only for count_atom)'
    write(MsgOut,'(A)') '# HB_distance = 3.4              # the upper limit of (D .. A) distance (default: 3.4 A)'
    write(MsgOut,'(A)') '# DHA_angle   = 120.0            # the lower limit of (D-H .. A) angle (default: 120.0 deg)'
    write(MsgOut,'(A)') '# HDA_angle   =  20.0            # the upper limit of (H-D .. A) angle (default:  20.0 deg)'
  
    write(MsgOut,'(A)') ''

    return

  end subroutine show_hbond_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DM
  !! @param[in]    handle         : unit number of control file
  !! @param[out]   hbond_opt_info : hbond OPTION section control parameters 
  !!                                information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, hbond_opt_info)

    ! parameters
    character(*),           parameter     :: Section = 'HBOND_OPTION'

    ! formal argments
    integer,                intent(in)    :: handle
    type(s_hbond_opt_info), intent(inout) :: hbond_opt_info

    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    !---- options for group name (hbond analysis) -----
    call read_ctrlfile_integer(handle, Section, &
                              'analysis_atom', hbond_opt_info%analysis_atom)

    call read_ctrlfile_integer(handle, Section, &
                              'target_atom',   hbond_opt_info%target_atom)

    call read_ctrlfile_integer(handle, Section, &
                              'recenter',      hbond_opt_info%recenter)

    call read_ctrlfile_type   (handle, Section, &
                              'output_type', hbond_opt_info%output_type, &
                              HBOutputType)

    call read_ctrlfile_real   (handle, Section, &
                              'hb_distance', hbond_opt_info%hb_distance)

    call read_ctrlfile_real   (handle, Section, &
                              'dha_angle', hbond_opt_info%dha_angle)

    call read_ctrlfile_real   (handle, Section, &
                              'hda_angle', hbond_opt_info%hda_angle)

    call read_ctrlfile_string (handle, Section, 'solvent_list',   &
                               hbond_opt_info%solvent_list)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options (H-bond)'

      write(MsgOut,'(A20,A5,I0)') &
           '  recenter        = ', 'group', hbond_opt_info%recenter

      write(MsgOut,'(A20,A5,I0)') &
           '  analysis atom   = ', 'group', hbond_opt_info%analysis_atom

      write(MsgOut,'(A20,A5,I0)') &
           '  target atom     = ', 'group', hbond_opt_info%target_atom

      write(MsgOut,'(A20,A)') &
           '  output type     = ', HBOutputType(hbond_opt_info%output_type)

      write(MsgOut,'(A20,A)') &
           '  solvent list    = ', trim(hbond_opt_info%solvent_list)

      write(MsgOut,'(A20,F8.3)') &
           '  H-bond distance = ', hbond_opt_info%hb_distance

      write(MsgOut,'(A20,F8.3)') &
           '  D-H .. A angle  = ', hbond_opt_info%dha_angle

      write(MsgOut,'(A20,F8.3)') &
           '  H-D .. A angle  = ', hbond_opt_info%hda_angle

      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DM
  !! @param[in]    hbond_opt_info : hbond OPTION section control parameters 
  !!                                information
  !! @param[inout] hbond_option   : hbond option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine hbond_setup_option(hbond_opt_info, hbond_option)

    ! formal arguments
    type(s_hbond_opt_info), intent(in)    :: hbond_opt_info
    type(s_hbond_option),   intent(inout) :: hbond_option

    ! local variables
    integer :: ndata


    hbond_option%recenter      = hbond_opt_info%recenter
    hbond_option%analysis_atom = hbond_opt_info%analysis_atom
    hbond_option%target_atom   = hbond_opt_info%target_atom
    hbond_option%output_type   = hbond_opt_info%output_type
    hbond_option%hb_distance   = hbond_opt_info%hb_distance
    hbond_option%dha_angle     = hbond_opt_info%dha_angle
    hbond_option%hda_angle     = hbond_opt_info%hda_angle

    ! solvent molecule list
    !
    if (allocated(hbond_option%solvent_list))  &
      deallocate(hbond_option%solvent_list)

    if (hbond_opt_info%solvent_list .eq. '' .or.  &
        index(hbond_opt_info%solvent_list, 'NONE') >= 1) then
      ndata = 0

    else
      ndata = split_num(hbond_opt_info%solvent_list)
      allocate(hbond_option%solvent_list(ndata))

      call split(ndata, ndata, hbond_opt_info%solvent_list, &
                 hbond_option%solvent_list)
    end if

    return

  end subroutine hbond_setup_option

end module hbond_option_mod
