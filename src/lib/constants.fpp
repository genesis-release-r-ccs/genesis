!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   constants_mod
!> @brief   constant variables (this module can be used like common in F77)
!! @authors Yuji Sugita (YS), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module constants_mod

#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! constants for real precision
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: dp = selected_real_kind(15, 307)

  ! constants for integer types
  integer,  public, parameter :: int1 = selected_int_kind(2)
  integer,  public, parameter :: int2 = selected_int_kind(4)
  integer,  public, parameter :: int4 = selected_int_kind(8)
  integer,  public, parameter :: int8 = selected_int_kind(16)

#if defined(_LARGE_INT)

  integer,  public, parameter :: iintegers = int8

#else

  integer,  public, parameter :: iintegers = int4

#endif

#if defined(_SINGLE)

  character(6), public, parameter :: precision_char = "single"

  integer,  public, parameter :: wp  = sp
  integer,  public, parameter :: wip = sp

#ifdef MPI

  integer,  public, parameter :: mpi_wp_real     = MPI_Real4
  integer,  public, parameter :: mpi_wp_complex  = MPI_Complex
  integer,  public, parameter :: mpi_wip_real    = MPI_Real4
  integer,  public, parameter :: mpi_wip_complex = MPI_Complex

#endif

#elif defined(_MIXED)

  character(6), public, parameter :: precision_char = "mixed"

  integer,  public, parameter :: wp  = sp
  integer,  public, parameter :: wip = dp

#ifdef MPI

  integer,  public, parameter :: mpi_wp_real     = MPI_Real4
  integer,  public, parameter :: mpi_wp_complex  = MPI_Complex
  integer,  public, parameter :: mpi_wip_real    = MPI_Real8
  integer,  public, parameter :: mpi_wip_complex = MPI_Double_complex

#endif

#else

  character(6), public, parameter :: precision_char = "double"

  integer,  public, parameter :: wp  = dp
  integer,  public, parameter :: wip = dp

#ifdef MPI

  integer,  public, parameter :: mpi_wp_real     = MPI_Real8
  integer,  public, parameter :: mpi_wp_complex  = MPI_Double_complex
  integer,  public, parameter :: mpi_wip_real    = MPI_Real8
  integer,  public, parameter :: mpi_wip_complex = MPI_Double_complex

#endif

#endif

  ! constants
  real(wp), public, parameter :: PI               = 3.14159265358979323846_wp
  real(wp), public, parameter :: SQRT_PI          = sqrt(PI)
  real(wp), public, parameter :: RAD              = PI/180.0_wp
  real(wp), public, parameter :: EPS              = 1.0e-10_wp
  real(wp), public, parameter :: HALF             = 0.5_wp
  real(wp), public, parameter :: ONE_THIRD        = 1.0_wp/3.0_wp
  real(wp), public, parameter :: AKMA_PS          = 4.88882129e-02_wp
  real(wp), public, parameter :: ELECOEF_CHARMM   = 332.0716_wp
  real(wp), public, parameter :: ELECOEF_AMBER    = 332.05221729_wp
  real(wp), public, parameter :: ELECOEF_GROMACS  = 332.0636930_wp
  real(wp), public, parameter :: KBOLTZ           = 1.987191e-03_wp
  real(wp), public, parameter :: ATMOS_P          = 1.4584007e-05_wp
  real(wp), public, parameter :: P_ATMOS          = 1.0_wp/ATMOS_P
  real(wp), public, parameter :: CAL2JOU          = 4.184_wp
  real(wp), public, parameter :: JOU2CAL          = 1.0_wp / CAL2JOU
  real(wp), public, parameter :: ATM2BAR          = 1.01325_wp
  real(wp), public, parameter :: BAR2ATM          = 1.0_wp / ATM2BAR
  real(wp), public, parameter :: CONV_UNIT_ENE    = 627.5095_wp
  real(wp), public, parameter :: CONV_UNIT_LEN    = 0.52917721092_wp
  real(wp), public, parameter :: CONV_UNIT_FORCE  = CONV_UNIT_ENE / CONV_UNIT_LEN
  real(wp), public, parameter :: ELMASS           = 1.660538921e-27_wp / 9.10938291e-31_wp
  real(wp), public, parameter :: HARTREE_WAVENUM  = 2.19474631e+05_wp
  real(wp), public, parameter :: VLIGHT_IN_AU     = 1.3703599918E+02_wp
  real(wp), public, parameter :: AVOGADRO         = 6.02214129e+23_wp

  real(wp), public, save      :: ELECOEF          = ELECOEF_CHARMM

  real(wp), public, save      :: LIGHT_ATOM_MASS_LIMIT = 2.1_wp

#ifdef _SINGLE
  real(wp), public, parameter :: EPS_CMAP         = 1.0e-3_wp
#else
  real(wp), public, parameter :: EPS_CMAP         = 1.0e-4_wp
#endif

end module constants_mod
