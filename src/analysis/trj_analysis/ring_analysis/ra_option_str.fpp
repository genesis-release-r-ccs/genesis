!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_option_str_mod
!> @brief   structure of option information
!! @authors Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_option_str_mod

  use select_atoms_str_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    real(wp)                      :: min_contact
    integer                       :: num_ring
    integer                       :: num_all_ring
    integer                       :: num_ring_type
    integer                       :: num_other
    integer,      allocatable     :: atomflag(:) 
    integer,      allocatable     :: ring_atoms(:) 
    integer,      allocatable     :: ring_all_atoms(:) 
    integer,      allocatable     :: other_atoms(:) 
    character(6), allocatable     :: ring_atoms_in(:) 

  end type s_option

end module ra_option_str_mod
