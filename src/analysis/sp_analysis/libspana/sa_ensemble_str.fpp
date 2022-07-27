!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_ensemble_str_mod
!> @brief   structure of ensemble
!! @authors Takaharu Mori (TM), Tadashi Ando (TA), Isseki Yu(IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_ensemble_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_ensemble

    integer             :: ensemble

  end type s_ensemble

  ! parameters
  integer,      public, parameter :: EnsembleNVE         = 1
  integer,      public, parameter :: EnsembleNVT         = 2
  integer,      public, parameter :: EnsembleNPT         = 3
  integer,      public, parameter :: EnsembleNPAT        = 4
  integer,      public, parameter :: EnsembleNPgT        = 5

  character(*), public, parameter :: EnsembleTypes(5)  = (/'NVE        ', &
                                                           'NVT        ', &
                                                           'NPT        ', &
                                                           'NPAT       ', &
                                                           'NPgT       '/)

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

    return

  end subroutine init_ensemble

end module sa_ensemble_str_mod
