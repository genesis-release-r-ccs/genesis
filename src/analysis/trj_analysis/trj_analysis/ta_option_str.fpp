!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_option_str_mod

  use select_atoms_str_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only

    logical                       :: out_dis
    integer,          allocatable :: dist_list(:,:)
    real(wp),         allocatable :: distance(:)

    logical                       :: out_ang
    integer,          allocatable :: angl_list(:,:)
    real(wp),         allocatable :: angle(:)

    logical                       :: out_tor
    integer,          allocatable :: tors_list(:,:)
    real(wp),         allocatable :: torsion(:)

    logical                       :: out_cdis
    integer,          allocatable :: cdist_group(:,:)
    real(wp),         allocatable :: cdistance(:)

    logical                       :: out_cang
    integer,          allocatable :: cangl_group(:,:)
    real(wp),         allocatable :: cangle(:)

    logical                       :: out_ctor
    integer,          allocatable :: ctor_group(:,:)
    real(wp),         allocatable :: ctorsion(:)

    logical                       :: out_sum_dis
    integer,          allocatable :: dist_num(:)
    real(wp),         allocatable :: dist_weight(:,:)    

    type(s_selatoms), allocatable :: selatoms(:)

  end type s_option

end module ta_option_str_mod
