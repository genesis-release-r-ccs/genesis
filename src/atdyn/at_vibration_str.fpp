!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_vibration_str_mod
!> @brief   structure of molecular vibration
!! @authors Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_vibration_str_mod

  use fileio_minfo_mod
  use messages_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_vibration
    integer                :: nreplica
    integer                :: vib_natoms
    integer, allocatable   :: vibatom_id(:)
    integer                :: minfo_natoms
    integer, allocatable   :: minfoatom_id(:)
    character(MaxFilename) :: minfofile
    character(MaxFilename) :: minfo_folder
    real(wp)               :: diff_stepsize
    real(wp)               :: cutoff
    type(s_minfo)          :: minfo_in

    logical                :: gengrid
    logical                :: grid_ene_only
    character(MaxFilename) :: gridfile
    character(MaxFilename) :: datafile
  end type s_vibration

  ! parameters
  integer,      public, parameter :: RunModeHarm = 1
  integer,      public, parameter :: RunModeQFF  = 2
  integer,      public, parameter :: RunModeGRID = 3
  
  character(*), public, parameter :: RunModeTypes(3) = (/'HARM',&
                                                         'QFF ',&
                                                         'GRID'/)
  ! subroutines

end module at_vibration_str_mod

