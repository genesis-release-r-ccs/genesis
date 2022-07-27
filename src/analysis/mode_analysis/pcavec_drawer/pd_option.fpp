!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pd_option_mod
!> @brief   module for analysis options
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pd_option_mod

  use pd_option_str_mod
  use select_mod
  use select_atoms_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only       = .false.
    logical                         :: allow_backup     = .false.
    integer                         :: mode_no          = 1
    real(wp)                        :: expand_vector    = 25.0
    real(wp)                        :: arrow_length     = 1.5
    logical                         :: arrow_reverse    = .false.
    character(50)                   :: vector_color_vmd = 'yellow'
    character(50)                   :: vector_color_pml = '0, 1, 1'
    real(wp)                        :: cylinder_radius  = 0.15
    real(wp)                        :: cone_radius      = 0.45
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
    write(MsgOut,'(A)') 'check_only       = NO            # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup     = NO            # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'mode_no          = 1'
    write(MsgOut,'(A)') 'expand_vector    = 25.0'
    write(MsgOut,'(A)') 'arrow_length     = 1.5'
    write(MsgOut,'(A)') 'cylinder_radius  = 0.15'
    write(MsgOut,'(A)') 'cone_radius      = 0.45'
    write(MsgOut,'(A)') 'arrow_reverse    = NO            # (YES/NO)'
    write(MsgOut,'(A)') 'vector_color_vmd = yellow        # color name for [VMD]'
    write(MsgOut,'(A)') 'vector_color_pml = 0, 1, 1       # RGB color (0.0-1.0) for [PyMOL]'
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

    call read_ctrlfile_logical(handle, Section, 'check_only',       &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup',     &
                               opt_info%allow_backup)

    call read_ctrlfile_integer(handle, Section, 'mode_no',          &
                               opt_info%mode_no)

    call read_ctrlfile_real   (handle, Section, 'expand_vector',    &
                               opt_info%expand_vector)

    call read_ctrlfile_real   (handle, Section, 'arrow_length',     &
                               opt_info%arrow_length)

    call read_ctrlfile_logical(handle, Section, 'arrow_reverse',    &
                               opt_info%arrow_reverse)

    call read_ctrlfile_string (handle, Section, 'vector_color_vmd', &
                               opt_info%vector_color_vmd)

    call read_ctrlfile_string (handle, Section, 'vector_color_pml', &
                               opt_info%vector_color_pml)

    call read_ctrlfile_real   (handle, Section, 'cylinder_radius',  &
                               opt_info%cylinder_radius)

    call read_ctrlfile_real   (handle, Section, 'cone_radius',      &
                               opt_info%cone_radius)

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

    write(MsgOut,'(A20,I0)') &
                    '  mode number     = ', opt_info%mode_no

    if (opt_info%arrow_reverse) then
      write(MsgOut,'(A20,A3)') '  arrow reverse   = ', 'yes'
    else
      write(MsgOut,'(A20,A2)') '  arrow_reverse   = ', 'no'
    end if

    write(MsgOut,'(A20,F10.3)') &
                    '  expand_vector   = ', opt_info%expand_vector

    write(MsgOut,'(A20,F10.3)') &
                    '  arrow_length    = ', opt_info%arrow_length

    write(MsgOut,'(A20,F10.5)') &
                    '  cylinder_radius = ', opt_info%cylinder_radius

    write(MsgOut,'(A20,F10.5)') &
                    '  cone_radius     = ', opt_info%cone_radius

    write(MsgOut,'(A20,A)') &
                    '  vector_color_vmd= ', opt_info%vector_color_vmd

    write(MsgOut,'(A20,A)') &
                    '  vector_color_pml= ', opt_info%vector_color_pml

    write(MsgOut,'(A)') ''

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      TM
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_option(opt_info, option)

    ! formal arguments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_option),          intent(inout) :: option


    ! check only
    option%check_only       = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    option%mode_no          = opt_info%mode_no
    option%arrow_reverse    = opt_info%arrow_reverse
    option%expand_vector    = opt_info%expand_vector
    option%arrow_length     = opt_info%arrow_length
    option%cylinder_radius  = opt_info%cylinder_radius
    option%cone_radius      = opt_info%cone_radius
    option%vector_color_vmd = opt_info%vector_color_vmd
    option%vector_color_pml = opt_info%vector_color_pml

    write(MsgOut,'(A)') ''

    return

  end subroutine setup_option

end module pd_option_mod
