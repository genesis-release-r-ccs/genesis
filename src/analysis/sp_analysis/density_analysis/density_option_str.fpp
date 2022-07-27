!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   density_option_str_mod
!> @brief   structure of option information
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module density_option_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_density_option

    logical    :: verbose
    integer    :: output_format
    integer    :: density_type
    integer    :: gid_solute
    integer    :: gid_solvent

    real(wp)   :: ana_range
    real(wp)   :: voxel_size
    integer    :: recenter
    real(wp)   :: magnification

  end type s_density_option

  ! parameters

  ! density type
  integer,          public, parameter :: DensityTypeNumber   = 1
  integer,          public, parameter :: DensityTypeElectron = 2
  
  character(len=8), public, parameter :: DensityTypes(2) = (/'NUMBER  ', &
                                                             'ELECTRON'/)

  ! density file format
  integer,          public, parameter :: DensityFormatXPLOR = 1
  integer,          public, parameter :: DensityFormatCCP4  = 2
  integer,          public, parameter :: DensityFormatDX    = 3
  
  character(len=5), public, parameter :: DensityFormatTypes(3) = (/'XPLOR',  &
                                                                   'CCP4 ',  &
                                                                   'DX   '/)

end module density_option_str_mod
