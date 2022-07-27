!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rdf_option_mod
!> @brief   module for analysis options
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rdf_option_mod

  use rdf_option_str_mod
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
  type, public :: s_rdf_opt_info
    integer    :: analysis_atom_solute  = 0
    integer    :: analysis_atom_solvent = 0
    integer    :: rmode                 = RModeRadial
    real(wp)   :: ana_range             = 0.0_wp
    real(wp)   :: binsize               = 0.0_wp
    real(wp)   :: bulk_value            = 0.0_wp 
    real(wp)   :: bulk_region           = 10.0_wp
    real(wp)   :: voxel_size            = 0.0_wp
    integer    :: recenter              = 0
  end type s_rdf_opt_info

  ! subroutines
  public  :: show_rdf_ctrl_option
  public  :: read_ctrl_option
  public  :: rdf_setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_rdf_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_rdf_ctrl_option

    write(MsgOut,'(A)') '[RDF_OPTION]'
    write(MsgOut,'(A)') 'solute        = 1              # solute group number'
    write(MsgOut,'(A)') 'solvent       = 2              # solvent group number'
    write(MsgOut,'(A)') 'rmode         = radial         # radial/proximal'
    write(MsgOut,'(A)') 'range         = 10             # range for output'
    write(MsgOut,'(A)') 'binsize       = 0.5            # bin size'
    write(MsgOut,'(A)') 'bulk_region   = 10             # fraction of the last profile used to calculate bulk density' 
    write(MsgOut,'(A)') 'voxel_size    = 0.0            # voxel size'
    write(MsgOut,'(A)') 'recenter      = 1              # recenter group number'

    return

  end subroutine show_rdf_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      IY
  !! @param[in]    handle       : unit number of control file
  !! @param[out]   rdf_opt_info : rdf OPTION section control parameters 
  !!                              information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, rdf_opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Rdf_option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_rdf_opt_info),    intent(inout) :: rdf_opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)

    !---- options for group name (Proximal distribution function) -----
    call read_ctrlfile_integer(handle, Section, &
                              'solute', rdf_opt_info%analysis_atom_solute)

    call read_ctrlfile_integer(handle, Section, &
                              'solvent', rdf_opt_info%analysis_atom_solvent)

    call read_ctrlfile_type(handle, Section, &
                              'rmode', rdf_opt_info%rmode, RModeTypes)

    call read_ctrlfile_real(handle, Section, &
                              'range', rdf_opt_info%ana_range)

    call read_ctrlfile_real(handle, Section, &
                              'binsize', rdf_opt_info%binsize)

    call read_ctrlfile_real(handle, Section, &
                              'bulk_region', rdf_opt_info%bulk_region)

    call read_ctrlfile_real(handle, Section, &
                              'bulk_value', rdf_opt_info%bulk_value)

    call read_ctrlfile_real(handle, Section, &
                              'voxel_size', rdf_opt_info%voxel_size)

    call read_ctrlfile_integer(handle, Section, &
                              'recenter', rdf_opt_info%recenter)

    call end_ctrlfile_section(handle)

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      IY
  !! @param[in]    rdf_opt_info : rdf OPTION section control parameters 
  !!                              information
  !! @param[inout] rdf_option   : rdf option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rdf_setup_option(rdf_opt_info, rdf_option)

    ! formal arguments
    type(s_rdf_opt_info),    intent(in)    :: rdf_opt_info
    type(s_rdf_option),      intent(inout) :: rdf_option


    rdf_option%gid_solvent = rdf_opt_info%analysis_atom_solvent
    rdf_option%gid_solute  = rdf_opt_info%analysis_atom_solute
    rdf_option%rmode       = rdf_opt_info%rmode
    rdf_option%ana_range   = rdf_opt_info%ana_range
    rdf_option%binsize     = rdf_opt_info%binsize
    rdf_option%bulk_region = rdf_opt_info%bulk_region
    rdf_option%bulk_value  = rdf_opt_info%bulk_value
    rdf_option%voxel_size  = rdf_opt_info%voxel_size
    rdf_option%recenter    = rdf_opt_info%recenter

    return

  end subroutine rdf_setup_option

end module rdf_option_mod
