!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pd_option_str_mod
!> @brief   structure of option information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pd_option_str_mod

  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option
    logical                       :: check_only
    integer                       :: mode_no
    real(wp)                      :: expand_vector
    real(wp)                      :: arrow_length
    logical                       :: arrow_reverse
    character(50)                 :: vector_color_vmd
    character(50)                 :: vector_color_pml
    real(wp)                      :: cylinder_radius
    real(wp)                      :: cone_radius
  end type s_option

end module pd_option_str_mod
