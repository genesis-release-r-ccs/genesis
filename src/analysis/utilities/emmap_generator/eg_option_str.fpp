!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   eg_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module eg_option_str_mod

  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical             :: check_only
    integer             :: map_format
    real(wp)            :: voxel_size
    real(wp)            :: sigma
    real(wp)            :: tolerance
    logical             :: auto_margin
    real(wp)            :: margin_size_x
    real(wp)            :: margin_size_y
    real(wp)            :: margin_size_z
    real(wp)            :: x0
    real(wp)            :: y0
    real(wp)            :: z0
    real(wp)            :: box_size_x
    real(wp)            :: box_size_y
    real(wp)            :: box_size_z

  end type s_option

  integer,      public, parameter :: MapFormatSITUS       = 1

  character(*), public, parameter :: MapFormatTypes(1)  = (/'SITUS'/)

end module eg_option_str_mod
