!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   eg_option_mod
!> @brief   module for trajectory conversion options
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module eg_option_mod

  use eg_option_str_mod
  use fileio_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_opt_info
    logical                         :: check_only    = .false.
    logical                         :: allow_backup  = .false.
    integer                         :: map_format    = MapFormatSITUS
    real(wp)                        :: voxel_size    = 1.0_wp
    real(wp)                        :: x0            = -100.0_wp
    real(wp)                        :: y0            = -100.0_wp
    real(wp)                        :: z0            = -100.0_wp
    real(wp)                        :: box_size_x    =  200.0_wp
    real(wp)                        :: box_size_y    =  200.0_wp
    real(wp)                        :: box_size_z    =  200.0_wp
    real(wp)                        :: sigma         = 2.5_wp
    real(wp)                        :: tolerance     = 0.001_wp
    logical                         :: auto_margin   = .true.
    real(wp)                        :: margin_size_x =  20.0_wp
    real(wp)                        :: margin_size_y =  20.0_wp
    real(wp)                        :: margin_size_z =  20.0_wp
  end type s_opt_info

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
    write(MsgOut,'(A)') 'check_only     = NO              # only checking input files (YES/NO)'
    write(MsgOut,'(A)') 'allow_backup   = NO              # backup existing output files (YES/NO)'
    write(MsgOut,'(A)') 'map_format     = SITUS           # (SITUS)'
    write(MsgOut,'(A)') 'auto_margin    = YES             # (YES/NO)'
    write(MsgOut,'(A)') 'margin_size_x  = 20.0            # margin size'
    write(MsgOut,'(A)') 'margin_size_y  = 20.0            # margin size'
    write(MsgOut,'(A)') 'margin_size_z  = 20.0            # margin size'
    write(MsgOut,'(A)') 'voxel_size     =  1.0            # voxel size'
    write(MsgOut,'(A)') 'sigma          =  2.5            # resolution'
    write(MsgOut,'(A)') 'tolerance      = 0.001           # tolerance'
    write(MsgOut,'(A)') '# x0             = -100            # x-coordinates of first voxel'
    write(MsgOut,'(A)') '# y0             = -100            # y-coordinates of first voxel'
    write(MsgOut,'(A)') '# z0             = -100            # z-coordinates of first voxel'
    write(MsgOut,'(A)') '# box_size_x     =  200            # box size in x-axis'
    write(MsgOut,'(A)') '# box_size_y     =  200            # box size in y-axis'
    write(MsgOut,'(A)') '# box_size_z     =  200            # box size in z-axis'
    write(MsgOut,'(A)') ' '

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
    character(*),            parameter     :: Section = 'Option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_opt_info),        intent(inout) :: opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'check_only',   &
                               opt_info%check_only)

    call read_ctrlfile_logical(handle, Section, 'allow_backup', &
                               opt_info%allow_backup)

    call read_ctrlfile_type   (handle, Section, 'map_format',   &
                               opt_info%map_format, MapFormatTypes)

    call read_ctrlfile_real   (handle, Section, 'voxel_size',   &
                               opt_info%voxel_size)

    call read_ctrlfile_real   (handle, Section, 'x0',           &
                               opt_info%x0)

    call read_ctrlfile_real   (handle, Section, 'y0',           &
                               opt_info%y0)

    call read_ctrlfile_real   (handle, Section, 'z0',           &
                               opt_info%z0)

    call read_ctrlfile_real   (handle, Section, 'box_size_x',   &
                               opt_info%box_size_x)

    call read_ctrlfile_real   (handle, Section, 'box_size_y',   &
                               opt_info%box_size_y)

    call read_ctrlfile_real   (handle, Section, 'box_size_z',   &
                               opt_info%box_size_z)

    call read_ctrlfile_real   (handle, Section, 'sigma',        &
                               opt_info%sigma)

    call read_ctrlfile_real   (handle, Section, 'tolerance',    &
                               opt_info%tolerance)

    call read_ctrlfile_logical(handle, Section, 'auto_margin',  &
                               opt_info%auto_margin)

    call read_ctrlfile_real   (handle, Section, 'margin_size_x',&
                               opt_info%margin_size_x)

    call read_ctrlfile_real   (handle, Section, 'margin_size_y',&
                               opt_info%margin_size_y)

    call read_ctrlfile_real   (handle, Section, 'margin_size_z',&
                               opt_info%margin_size_z)

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
    write(MsgOut,'(A20,A10)') &
                    '  map format      = ', MapFormatTypes(opt_info%map_format)
    write(MsgOut,'(A)') ' '

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      NT
  !! @param[in]    opt_info : OPTION section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[inout] molecule : molecule information
  !! @param[out]   option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_option(opt_info, option)
  
    ! formal argments
    type(s_opt_info),        intent(in)    :: opt_info
    type(s_option),          intent(inout) :: option

    ! check only
    option%check_only    = opt_info%check_only

    ! backup output file
    call setup_backup(opt_info%allow_backup)

    option%map_format    = opt_info%map_format
    option%voxel_size    = opt_info%voxel_size
    option%x0            = opt_info%x0
    option%y0            = opt_info%y0
    option%z0            = opt_info%z0
    option%box_size_x    = opt_info%box_size_x
    option%box_size_y    = opt_info%box_size_y
    option%box_size_z    = opt_info%box_size_z
    option%sigma         = opt_info%sigma
    option%tolerance     = opt_info%tolerance
    option%auto_margin   = opt_info%auto_margin
    option%margin_size_x = opt_info%margin_size_x
    option%margin_size_y = opt_info%margin_size_y
    option%margin_size_z = opt_info%margin_size_z

    return

  end subroutine setup_option

end module eg_option_mod
