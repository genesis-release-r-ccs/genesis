!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_ensemble_str_mod
!> @brief   structure of ensemble
!! @authors Takaharu Mori (TM), Tadashi Ando (TA)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_ensemble_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_ensemble

    integer               :: ensemble
    real(wp)              :: temperature
    real(wp)              :: pressure
    real(wp)              :: gamma
    integer               :: tpcontrol

    ! parameters for Berendsen
    real(wp)              :: tau_t
    real(wp)              :: tau_p
    real(wp)              :: compressibility

    ! parameters for Langevin
    real(wp)              :: gamma_t
    real(wp)              :: gamma_p
    real(wp)              :: pforce(1:3)

    ! parameters for Andersen
    real(wp)              :: softness

    ! parameters for NPT and NPAT
    real(wp)              :: pmass
    integer               :: isotropy

    logical               :: use_barostat

  end type s_ensemble

  ! parameters
  integer,      public, parameter :: EnsembleNVE         = 1
  integer,      public, parameter :: EnsembleNVT         = 2
  integer,      public, parameter :: EnsembleNPT         = 3
  integer,      public, parameter :: EnsembleNPAT        = 4
  integer,      public, parameter :: EnsembleNPgT        = 5

  integer,      public, parameter :: TpcontrolNO         = 1
  integer,      public, parameter :: TpcontrolBerendsen  = 2
  integer,      public, parameter :: TpcontrolAndersen   = 3
  integer,      public, parameter :: TpcontrolEvans      = 4
  integer,      public, parameter :: TpcontrolNoseHoover = 5
  integer,      public, parameter :: TpcontrolLangevin   = 6
  integer,      public, parameter :: TpcontrolGauss      = 7
  integer,      public, parameter :: TpcontrolBussi      = 8
  
  integer,      public, parameter :: IsotropyISO         = 1
  integer,      public, parameter :: IsotropySEMI_ISO    = 2
  integer,      public, parameter :: IsotropyANISO       = 3
  integer,      public, parameter :: IsotropyXY_Fixed    = 4

  character(*), public, parameter :: EnsembleTypes(5)  = (/'NVE        ', &
                                                           'NVT        ', &
                                                           'NPT        ', &
                                                           'NPAT       ', &
                                                           'NPgT       '/)

  character(*), public, parameter :: TpcontrolTypes(8) = (/'NO         ', &
                                                           'BERENDSEN  ', &
                                                           'ANDERSEN   ', &
                                                           'EVANS      ', &
                                                           'NOSE-HOOVER', &
                                                           'LANGEVIN   ', &
                                                           'GAUSS      ', &
                                                           'BUSSI      '/)

  character(*), public, parameter :: IsotropyTypes(4)  = (/'ISO        ', &
                                                           'SEMI-ISO   ', &
                                                           'ANISO      ', &
                                                           'XY-FIXED   '/)

  ! subroutines
  public  :: init_ensemble

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_ensemble
  !> @brief        initialize ensemble information
  !! @authors      TM
  !! @param[out]   ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_ensemble(ensemble)

    ! formal arguments
    type(s_ensemble),        intent(inout) :: ensemble


    ensemble%ensemble        = 0
    ensemble%temperature     = 0.0_wp
    ensemble%pressure        = 0.0_wp
    ensemble%gamma           = 0.0_wp
    ensemble%tpcontrol       = 0
    ensemble%tau_t           = 0.0_wp
    ensemble%tau_p           = 0.0_wp
    ensemble%compressibility = 0.0_wp
    ensemble%gamma_t         = 0.0_wp
    ensemble%gamma_p         = 0.0_wp
    ensemble%softness        = 0.0_wp
    ensemble%isotropy        = 0
    ensemble%use_barostat    = .false.

    return

  end subroutine init_ensemble

end module at_ensemble_str_mod
