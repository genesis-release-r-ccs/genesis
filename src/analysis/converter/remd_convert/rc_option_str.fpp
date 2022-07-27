!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rc_option_str_mod
!> @brief   structure of option information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rc_option_str_mod

  use constants_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical                :: check_only
    integer                :: convert_type
    integer, allocatable   :: convert_ids(:)
    integer                :: num_replicas
    integer                :: nsteps
    integer                :: exchange_period
    integer                :: crdout_period
    integer                :: eneout_period
    integer                :: logout_period
    integer                :: trjout_period
    integer                :: trjout_format
    integer                :: trjout_type
    character(MaxLineLong) :: trjout_atom_exp
    type(s_selatoms)       :: trjout_atom
    logical                :: centering
    character(MaxLineLong) :: centering_atom_exp
    type(s_selatoms)       :: centering_atom
    real(wp)               :: center_coord(3)
    integer                :: pbcc_mode

  end type s_option

  ! parameters
  integer,      public, parameter :: ConvertTypeReplica   = 1
  integer,      public, parameter :: ConvertTypeParameter = 2

  character(*), public, parameter :: ConvertTypeTypes(2)  = &
                                             (/'REPLICA  ','PARAMETER'/)

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

    ! local variables
    integer                  :: dealloc_stat


    if (allocated(option%convert_ids)) then

      dealloc_stat = 0
      deallocate(option%convert_ids, stat = dealloc_stat)

      if (dealloc_stat /= 0) &
        call error_msg_dealloc

    end if

    call dealloc_selatoms(option%trjout_atom)
    call dealloc_selatoms(option%centering_atom)

    return

  end subroutine dealloc_option

end module rc_option_str_mod
