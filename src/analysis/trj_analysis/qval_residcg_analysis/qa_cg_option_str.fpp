!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qa_cg_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qa_cg_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    ! check only
    logical                       :: check_only

    ! parameters for Q-value cacluation
    logical                :: verbose
    real(wp)               :: lambda
    real(wp)               :: beta
    type(s_selatoms)       :: analysis_atom
    character(MaxLineLong) :: analysis_atom_exp

    ! atom selection
    type(s_selatoms),  allocatable :: selatoms(:)

  end type s_option

  ! structure
  type, public :: s_contact
    integer               :: n_pair
    real(wp), allocatable :: r0_ij(:)
    integer,  allocatable :: cnt_pair(:,:)
  end type s_contact


end module qa_cg_option_str_mod
