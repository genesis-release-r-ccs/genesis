!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rdf_option_str_mod
!> @brief   structure of option information
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rdf_option_str_mod

  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_rdf_option

    integer    :: gid_solute
    integer    :: gid_solvent
    integer    :: rmode
    real(wp)   :: ana_range
    real(wp)   :: binsize
    real(wp)   :: bulk_value
    real(wp)   :: bulk_region
    real(wp)   :: voxel_size
    integer    :: recenter

  end type s_rdf_option

  ! parameters
  integer,      public, parameter :: RModeRadial   = 1
  integer,      public, parameter :: RModeProximal = 2

  character(*), public, parameter :: RModeTypes(2) = (/'RADIAL  ', &
                                                       'PROXIMAL'/)

  ! subroutines
  public :: dealloc_rdf_option

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rdf_option
  !> @brief        deallocate option
  !! @authors      IY
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rdf_option(rdf_option)

    ! formal argments
    type(s_rdf_option),          intent(inout) :: rdf_option

    return

  end subroutine dealloc_rdf_option

end module rdf_option_str_mod
