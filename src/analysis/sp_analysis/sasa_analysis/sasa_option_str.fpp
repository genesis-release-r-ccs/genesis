!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sasa_option_str_mod
!> @brief   structure of option information
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sasa_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_sasa_option

    integer                :: gid_solute    = 0
    real(wp)               :: probe_radius  = 0.0_wp
    real(wp)               :: delta_z       = 0.0_wp
    integer                :: recenter      = 0
    character(MaxFilename) :: radi_file     = ''
    integer                :: out_style     = 0

  end type s_sasa_option

  ! parameters
  integer,      public, parameter :: OutStyleHistory       = 1
  integer,      public, parameter :: OutStyleAtomic        = 2
  integer,      public, parameter :: OutStyleAtomicHistory = 3

  character(*), public, parameter :: OutStyleTypes(3) = (/'HISTORY       ', &
                                                          'ATOMIC        ', &
                                                          'ATOMIC+HISTORY'/)

  ! subroutines
  public :: dealloc_sasa_option

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_sasa_option
  !> @brief        deallocate option
  !! @authors      IY
  !! @param[inout] sasa_option : sasa option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_sasa_option(sasa_option)

    ! formal argments
    type(s_sasa_option),          intent(inout) :: sasa_option

    return

  end subroutine dealloc_sasa_option

end module sasa_option_str_mod
