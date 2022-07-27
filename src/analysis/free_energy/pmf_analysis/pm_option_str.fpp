!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pm_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pm_option_str_mod

  use select_atoms_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    integer                       :: nreplicas
    integer                       :: dimension
    real(wp)                      :: temperature
    real(wp)                      :: cutoff
    real(wp), allocatable         :: grid_min(:)
    real(wp), allocatable         :: grid_max(:)
    integer,  allocatable         :: num_grids(:)
    real(wp), allocatable         :: center(:,:)
    real(wp), allocatable         :: band_width(:)
    logical,  allocatable         :: is_periodic(:)
    real(wp), allocatable         :: box_size(:)
    integer                       :: output_type
    real(wp), allocatable         :: delta_grid(:)

  end type s_option

  ! parameters for output type
  integer,      public, parameter :: OutputTypeGNUPLOT  = 1
  integer,      public, parameter :: OutputTypeMATLAB   = 2
  
  character(*), public, parameter :: OutputType(2) = (/'GNUPLOT', &
                                                      'MATLAB '/)

end module pm_option_str_mod
