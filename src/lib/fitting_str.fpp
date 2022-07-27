!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fitting_str_mod
!> @brief   structure of fitting information
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fitting_str_mod

  use select_atoms_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_fitting

    integer          :: fitting_method
    type(s_selatoms) :: fitting_atom
    integer          :: num_grids
    real(wp)         :: grid_size
    logical          :: mass_weight = .false.

    logical          :: z_translation  = .false.
    real(wp)         :: rot_matrix(3,3)
    real(wp)         :: com_ref(3)
    real(wp)         :: com_mov(3)
    real(wp)         :: rmsd
    integer          :: ret_code

  end type s_fitting

  ! parameters for fitting method
  integer,      public, parameter :: FittingMethodNO        = 1
  integer,      public, parameter :: FittingMethodTR_ROT    = 2
  integer,      public, parameter :: FittingMethodTR        = 3
  integer,      public, parameter :: FittingMethodTR_ZROT   = 4
  integer,      public, parameter :: FittingMethodXYTR      = 5
  integer,      public, parameter :: FittingMethodXYTR_ZROT = 6

  character(*), public, parameter :: FittingMethodTypes(6)  = (/'NO       ', &
                                                                'TR+ROT   ', &
                                                                'TR       ', &
                                                                'TR+ZROT  ', &
                                                                'XYTR     ', &
                                                                'XYTR+ZROT'/)

  ! parameters for fitting method
  integer,      public, parameter :: FittingFileREF         = 1
  integer,      public, parameter :: FittingFileFIT         = 2

  character(*), public, parameter :: FittingFileTypes(2)  = (/'REF', &
                                                              'FIT'/)

  ! parameters for fitting method
  integer,      public, parameter :: FittingMoveREF         = 1
  integer,      public, parameter :: FittingMoveSYS         = 2

  character(*), public, parameter :: FittingMoveTypes(2)  = (/'REF', &
                                                              'SYS'/)

  ! subroutines
  public  :: dealloc_fitting

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_fitting
  !> @brief        deallocate fitting
  !! @authors      NT
  !! @param[inout] fitting : fitting information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_fitting(fitting)

    ! formal argments
    type(s_fitting),         intent(inout) :: fitting


    call dealloc_selatoms(fitting%fitting_atom)

    return

  end subroutine dealloc_fitting

end module fitting_str_mod
