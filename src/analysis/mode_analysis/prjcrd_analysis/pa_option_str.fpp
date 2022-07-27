!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pa_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pa_option_str_mod

  use select_atoms_str_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                       :: check_only
    integer                       :: vcv_matrix
    integer                       :: num_pca
    type(s_selatoms)              :: analysis_atom

  end type s_option

  ! parameters for Vcv Matrix
  integer, public, parameter      :: VcvMatrixLocal  = 1
  integer, public, parameter      :: VcvMatrixGlobal = 2

  character(*), public, parameter :: VcvMatrices(2)  = (/'LOCAL ',&
                                                         'GLOBAL'/)
  ! subroutines
  public  :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      NT
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option


    call dealloc_selatoms(option%analysis_atom)

    return

  end subroutine dealloc_option

end module pa_option_str_mod
