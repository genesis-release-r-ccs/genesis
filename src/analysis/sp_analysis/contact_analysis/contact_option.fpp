!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   contact_option_mod
!> @brief   module for analysis options
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module contact_option_mod

  use contact_option_str_mod
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
  type, public :: s_contact_opt_info
    real(wp)   :: ana_range = 0.0_wp
    integer    :: recenter  = 0
    integer    :: mode      = ModeNumber
  end type s_contact_opt_info

  ! subroutines
  public  :: show_contact_ctrl_option
  public  :: read_ctrl_option
  public  :: contact_setup_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_contact_ctrl_option
  !> @brief        show control parameters in OPTION section
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_contact_ctrl_option

    write(MsgOut,'(A)') '[CONTACT_OPTION]'
    write(MsgOut,'(A)') 'range    = 10                  # range'
    write(MsgOut,'(A)') 'mode     = number              # mindist/number'

    return

  end subroutine show_contact_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_option
  !> @brief        read control parameters in OPTION section
  !! @authors      IY
  !! @param[in]    handle           : unit number of control file
  !! @param[out]   contact_opt_info : Contact OPTION section control parameters 
  !!               information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_option(handle, contact_opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Contact_option'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_contact_opt_info),intent(inout) :: contact_opt_info


    ! read parameters
    !
    call begin_ctrlfile_section(handle, Section)


    call read_ctrlfile_integer(handle, Section, &
                              'recenter', contact_opt_info%recenter)

    call read_ctrlfile_real(handle, Section, &
                              'range', contact_opt_info%ana_range)

    call read_ctrlfile_type(handle, Section, &
                              'mode', contact_opt_info%mode, ModeTypes)

    call end_ctrlfile_section(handle)

    return

  end subroutine read_ctrl_option

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_option
  !> @brief        setup option information
  !! @authors      IY
  !! @param[in]    contact_opt_info : Contact OPTION section control parameters 
  !!               information
  !! @param[inout] contact_option   : Contact option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine contact_setup_option(contact_opt_info, contact_option)

    ! formal arguments
    type(s_contact_opt_info),intent(in)    :: contact_opt_info
    type(s_contact_option),  intent(inout) :: contact_option


    contact_option%ana_range = contact_opt_info%ana_range
    contact_option%recenter  = contact_opt_info%recenter
    contact_option%mode      = contact_opt_info%mode

    return

  end subroutine contact_setup_option

end module contact_option_mod
