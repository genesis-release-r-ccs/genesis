!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sasa_option_mod
!> @brief   module for analysis options
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sasa_option_mod

  use sasa_option_str_mod
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
  type, public :: s_sasa_opt_info
    integer                 :: gid_solute    = 0
    real(wp)                :: probe_radius  = 1.4_wp
    real(wp)                :: delta_z       = 0.2_wp
    integer                 :: recenter      = 0
    character(MaxFilename)  :: radi_file     = ""
    integer                 :: out_style     = OutStyleHistory
  end type s_sasa_opt_info

  ! subroutines
  public  :: show_sasa_ctrl_option
  public  :: read_ctrl_option
  public  :: sasa_setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_sasa_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_sasa_ctrl_option

    write(MsgOut,'(A)') '[SASA_OPTION]'
    write(MsgOut,'(A)') 'solute          = 1          # target group number. If you do not specify,'
    write(MsgOut,'(A)') '                             SASA was calculated for all groups in the '
    write(MsgOut,'(A)') '                             section [SELECTION]'
    write(MsgOut,'(A)') 'recenter        = 1          # recenter group number'
    write(MsgOut,'(A)') 'output_style    = radial     # history, atomic, atomic+history'
    write(MsgOut,'(A)') 'radi_file       = radius.txt # atomic radius file'
    write(MsgOut,'(A)') 'probe_radius    = 1.4        # solvent probe radius'
    write(MsgOut,'(A)') 'delta_z         = 0.2        # delta-z'
     
    return

  end subroutine show_sasa_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      IY
  !! @param[in]    handle        : unit number of control file
  !! @param[out]   sasa_opt_info : sasa OPTION section control parameters 
  !!                               information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, sasa_opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Sasa_option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_sasa_opt_info),   intent(inout) :: sasa_opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, &
                              'solute', sasa_opt_info%gid_solute)

    call read_ctrlfile_integer(handle, Section, &
                              'recenter',sasa_opt_info%recenter)

    call read_ctrlfile_string(handle, Section, &
                              'radi_file', sasa_opt_info%radi_file)

    call read_ctrlfile_real(handle, Section, &
                              'probe_radius', sasa_opt_info%probe_radius)

    call read_ctrlfile_real(handle, Section, &
                              'delta_z', sasa_opt_info%delta_z)

    call read_ctrlfile_type(handle, Section, &
                              'output_style', sasa_opt_info%out_style, &
                              OutStyleTypes)

    call end_ctrlfile_section(handle)

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      IY
  !! @param[in]    sasa_opt_info : sasa OPTION section control parameters 
  !!                               information
  !! @param[inout] sasa_option   : sasa option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sasa_setup_option(sasa_opt_info, sasa_option)

    ! formal arguments
    type(s_sasa_opt_info),   intent(in)    :: sasa_opt_info
    type(s_sasa_option),     intent(inout) :: sasa_option


    sasa_option%recenter     = sasa_opt_info%recenter
    sasa_option%gid_solute   = sasa_opt_info%gid_solute
    sasa_option%radi_file    = sasa_opt_info%radi_file
    sasa_option%probe_radius = sasa_opt_info%probe_radius
    sasa_option%delta_z      = sasa_opt_info%delta_z
    sasa_option%out_style    = sasa_opt_info%out_style

    return

  end subroutine sasa_setup_option

end module sasa_option_mod
