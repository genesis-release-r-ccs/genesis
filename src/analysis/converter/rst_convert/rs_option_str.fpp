!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rs_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rs_option_str_mod

  use select_atoms_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical  :: check_only
    integer  :: rstin_format
    integer  :: rstout_format
    integer  :: rstout_type
    integer  :: nreplicas
    integer  :: pbcc_mode
    real(wp) :: shift_coord(3)
    real(dp) :: min_energy
    real(dp) :: min_delta_r
    integer  :: md_iseed
    integer  :: step
    real(dp) :: md_thermo_moment
    real(dp) :: md_baro_moment(3)
    real(dp) :: origin_x
    real(dp) :: origin_y
    real(dp) :: origin_z

  end type s_option

  ! parameters for input restart file format
  integer,      public, parameter :: RstFormatNAMD    = 1
  integer,      public, parameter :: RstFormatGENESIS = 2
  integer,      public, parameter :: RstFormatPDB     = 3
  integer,      public, parameter :: RstFormatCHARMM  = 4
  integer,      public, parameter :: RstFormatAMBER   = 5

  character(*), public, parameter :: RstFormatTypes(5) = (/'NAMD   ', &
                                                           'GENESIS', &
                                                           'PDB    ', &
                                                           'CHARMM ', &
                                                           'AMBER  '/)

  ! parameters for output genesis restart file type
  integer,      public, parameter :: RstTypeAuto   = 1
  integer,      public, parameter :: RstTypeMin    = 2
  integer,      public, parameter :: RstTypeMD     = 3
  integer,      public, parameter :: RstTypeREMD   = 4

  character(*), public, parameter :: RstTypeTypes(4) = (/'AUTO', 'Min ', 'MD  ','REMD'/)

end module rs_option_str_mod
