!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pn_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pn_option_str_mod

  use select_atoms_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    logical                       :: trajectory
    integer                       :: nreplicas
    integer                       :: nfiles

  end type s_option

end module pn_option_str_mod
