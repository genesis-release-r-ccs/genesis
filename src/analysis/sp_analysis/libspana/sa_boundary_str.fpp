!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_boundary_str_mod
!> @brief   structure of boundary
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_boundary_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_boundary

    integer                       :: type           = 0
    real(wp)                      :: origin_x       = 0.0_wp
    real(wp)                      :: origin_y       = 0.0_wp
    real(wp)                      :: origin_z       = 0.0_wp

    integer                       :: num_cells_x    = 0
    integer                       :: num_cells_y    = 0
    integer                       :: num_cells_z    = 0
    real(wp)                      :: box_size_x     = 0.0_wp
    real(wp)                      :: box_size_y     = 0.0_wp
    real(wp)                      :: box_size_z     = 0.0_wp
    real(wp)                      :: cell_size_x    = 0.0_wp
    real(wp)                      :: cell_size_y    = 0.0_wp
    real(wp)                      :: cell_size_z    = 0.0_wp

    integer                       :: num_domain(3)

  end type s_boundary

  ! parameters
  integer,      public, parameter :: BoundaryTypeNOBC     = 1
  integer,      public, parameter :: BoundaryTypePBC      = 2

  character(*), public, parameter :: BoundaryTypeTypes(2) = (/'NOBC', &
                                                              'PBC '/)

end module sa_boundary_str_mod
