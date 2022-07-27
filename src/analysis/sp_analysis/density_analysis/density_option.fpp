!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   density_option_mod
!> @brief   module for analysis options
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module density_option_mod

  use density_option_str_mod
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
  type, public :: s_density_opt_info
    integer    :: analysis_atom_solute
    integer    :: analysis_atom_solvent

    integer    :: density_type  = DensityTypeNumber
    integer    :: output_format = DensityFormatXPLOR
    real(wp)   :: ana_range     = 3.0_wp
    real(wp)   :: voxel_size    = 0.0_wp
    integer    :: recenter      = 0
    real(wp)   :: magnification = 1.0_wp
    logical    :: verbose       = .true.

  end type s_density_opt_info

  ! subroutines
  public  :: show_density_ctrl_option
  public  :: read_ctrl_option
  public  :: density_setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_density_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      DM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_density_ctrl_option

    write(MsgOut,'(A)') '[DENSITY_OPTION]'
    write(MsgOut,'(A)') 'density_type  = ELECTRON         # (NUMBER / ELECTRON)'
    write(MsgOut,'(A)') 'verbose       = yes              # (verbose log output)'
    write(MsgOut,'(A)') 'output_format = XPLOR            # (XPLOR / CCP4 / DX)'
    write(MsgOut,'(A)') 'solute        = 1                # solute group number'
    write(MsgOut,'(A)') 'solvent       = 2                # solvent group number'
    write(MsgOut,'(A)') 'recenter      = 1                # recenter group number'
    write(MsgOut,'(A)') 'range         = 5.0              # radius range'
    write(MsgOut,'(A)') 'voxel_size    = 0.0              # voxel size'
    write(MsgOut,'(A)') 'magnification = 5.0              # magnification in density output'
  
    write(MsgOut,'(A)') ''

    return

  end subroutine show_density_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      DM
  !! @param[in]    handle           : unit number of control file
  !! @param[out]   density_opt_info : density OPTION section control 
  !!                                  parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, density_opt_info)

    ! parameters
    character(*),             parameter     :: Section = 'Density_option'

    ! formal argments
    integer,                  intent(in)    :: handle
    type(s_density_opt_info), intent(inout) :: density_opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    !---- options for group name (density analysis) -----
    call read_ctrlfile_logical(handle, Section, &
                              'verbose', density_opt_info%verbose)

    call read_ctrlfile_integer(handle, Section, &
                              'solute',  density_opt_info%analysis_atom_solute)

    call read_ctrlfile_integer(handle, Section, &
                              'solvent', density_opt_info%analysis_atom_solvent)

    call read_ctrlfile_type   (handle, Section, &
                              'density_type', density_opt_info%density_type, &
                              DensityTypes)

    call read_ctrlfile_type   (handle, Section, &
                              'output_format', density_opt_info%output_format, &
                              DensityFormatTypes)

    call read_ctrlfile_integer(handle, Section, &
                              'recenter', density_opt_info%recenter)

    call read_ctrlfile_real   (handle, Section, &
                              'range', density_opt_info%ana_range)

    call read_ctrlfile_real   (handle, Section, &
                              'voxel_size', density_opt_info%voxel_size)

    call read_ctrlfile_real   (handle, Section, &
                              'magnification', density_opt_info%magnification)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Option> Parameters of Options (density)'

      if (density_opt_info%verbose) then
         write(MsgOut,'(A20,A3)') &
                         '  verbose         = ', 'yes'
      else
         write(MsgOut,'(A20,A2)') &
                         '  verbose         = ', 'no'
      end if

      write(MsgOut,'(A20,A5,I0)') &
           '  solvent atom    = ', 'group', &
           density_opt_info%analysis_atom_solvent

      write(MsgOut,'(A20,A5,I0)') &
           '  solute atom     = ', 'group', &
           density_opt_info%analysis_atom_solute

      write(MsgOut,'(A20,A5,I0)') &
           '  centering atom  = ', 'group', &
           density_opt_info%recenter

      write(MsgOut,'(A20,A)') &
           '  density type    = ', &
           DensityTypes(density_opt_info%density_type)

      write(MsgOut,'(A20,A)') &
           '  output format   = ', &
           DensityFormatTypes(density_opt_info%output_format)

      write(MsgOut,'(A20,F8.3)') &
           '  cutoff distance = ', &
           density_opt_info%ana_range

      write(MsgOut,'(A20,F8.3)') &
           '  voxel size      = ', &
           density_opt_info%voxel_size

      write(MsgOut,'(A20,F8.3)') &
           '  magnify         = ', &
           density_opt_info%magnification

      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      DM
  !! @param[in]    density_opt_info : density OPTION section control 
  !!                                  parameters information
  !! @param[inout] densiyt_option   : density option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine density_setup_option(density_opt_info, density_option)

    ! formal arguments
    type(s_density_opt_info), intent(in)    :: density_opt_info
    type(s_density_option),   intent(inout) :: density_option


    density_option%ana_range     = density_opt_info%ana_range
    density_option%voxel_size    = density_opt_info%voxel_size
    density_option%recenter      = density_opt_info%recenter
    density_option%gid_solvent   = density_opt_info%analysis_atom_solvent
    density_option%gid_solute    = density_opt_info%analysis_atom_solute
    density_option%density_type  = density_opt_info%density_type
    density_option%output_format = density_opt_info%output_format
    density_option%magnification = density_opt_info%magnification
    density_option%verbose       = density_opt_info%verbose

    return

  end subroutine density_setup_option

end module density_option_mod
