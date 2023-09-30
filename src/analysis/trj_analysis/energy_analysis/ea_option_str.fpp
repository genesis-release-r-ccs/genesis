!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_option_str_mod
!> @brief   structure of option information
!! @authors Motoshi Kamiya (MK)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_option_str_mod

  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                :: check_only
    logical                :: component
    logical                :: rest_component
    character(MaxFilename) :: remfile
    integer                :: tgt_parmid
    logical                :: mbar

  end type s_option

  public :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      MK
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal arguments
    type(s_option), intent(inout) :: option


    ! currently nothing to do

    return
  end subroutine dealloc_option

end module ea_option_str_mod
