!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_ensemble_str_mod
!> @brief   structure of ensemble
!! @authors Takaharu Mori (TM), Tadashi Ando (TA)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_ensemble_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_ensemble

    integer             :: ensemble
    real(wip)           :: temperature
    real(wip)           :: kinetic    
    real(wip)           :: degree
    real(wip)           :: pressure
    real(wip)           :: gamma
    integer             :: tpcontrol

    ! parameters for Berendsen and Bussi
    real(wip)           :: tau_t
    real(wip)           :: tau_p
    real(wip)           :: compressibility

    ! parameters for Langevin
    real(wip)           :: gamma_t
    real(wip)           :: gamma_p
    real(wip)           :: pforce(1:3)
    real(wip), allocatable :: random_force(:,:,:)

    ! parameters for Nose-Hoover or MTK
    integer             :: nhmultistep
    integer             :: nhchain

    ! parameters for RESPA with MTK
    real(wP)            :: pressure_short
    real(wip)           :: pressure_long

    ! parameters for NPT and NPAT
    real(wip)           :: pmass
    integer             :: isotropy

    logical             :: use_barostat
    logical             :: group_tp

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
  integer,      public, parameter :: TpcontrolMTK        = 8
  integer,      public, parameter :: TpcontrolBussi      = 9
  integer,      public, parameter :: TpcontrolNHC        = 10
  
  integer,      public, parameter :: IsotropyISO         = 1
  integer,      public, parameter :: IsotropySEMI_ISO    = 2
  integer,      public, parameter :: IsotropyANISO       = 3
  integer,      public, parameter :: IsotropyXY_Fixed    = 4

  character(*), public, parameter :: EnsembleTypes(5)  = (/'NVE        ', &
                                                           'NVT        ', &
                                                           'NPT        ', &
                                                           'NPAT       ', &
                                                           'NPgT       '/)

  character(*), public, parameter :: TpcontrolTypes(10) = (/'NO         ', &
                                                            'BERENDSEN  ', &
                                                            'ANDERSEN   ', &
                                                            'EVANS      ', &
                                                            'NOSE-HOOVER', &
                                                            'LANGEVIN   ', &
                                                            'GAUSS      ', &
                                                            'MTK        ', &
                                                            'BUSSI      ', &
                                                            'NHC        '/)

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
    ensemble%temperature     = 0.0_wip
    ensemble%pressure        = 0.0_wip
    ensemble%gamma           = 0.0_wip
    ensemble%tpcontrol       = 0
    ensemble%tau_t           = 0.0_wip
    ensemble%tau_p           = 0.0_wip
    ensemble%compressibility = 0.0_wip
    ensemble%gamma_t         = 0.0_wip
    ensemble%gamma_p         = 0.0_wip
    ensemble%isotropy        = 0
    ensemble%use_barostat    = .false.

    return

  end subroutine init_ensemble

end module sp_ensemble_str_mod
