!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_restraints_str_mod
!> @brief   structure of restraints information
!! @authors Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_restraints_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  type, public :: s_restraints

    logical                       :: restraint_flag = .false.
    logical                       :: verbose
    logical                       :: pressure_position
    logical                       :: pressure_rmsd
    integer                       :: nfunctions
    integer                       :: num_groups
    integer                       :: max_atoms

    ! restraints function (size = nfunctions)
    integer,              allocatable :: function(:)
    character(MaxLine),   allocatable :: constant(:)
    character(MaxLine),   allocatable :: reference(:)
    character(MaxLine),   allocatable :: select_index(:)
    integer,              allocatable :: direction(:)
    integer,              allocatable :: mode(:)
    integer,              allocatable :: exponent(:)
    character(MaxLine),   allocatable :: exponent_dist(:)
    character(MaxLine),   allocatable :: weight_dist(:)

    ! restraints group (size = num_group)
    character(MaxLine),   allocatable :: group(:)

    ! restraints list (size = num_group)
    integer,          allocatable :: atomlist(:,:)
    integer,          allocatable :: num_atoms(:)

  end type s_restraints

  ! parameters for allocatable variables
  integer,      public, parameter :: RestraintsFunc         = 1
  integer,      public, parameter :: RestraintsGroup        = 2
  integer,      public, parameter :: RestraintsList         = 3

  ! parameters for kind of functions
  integer,      public, parameter :: RestraintsFuncPOSI     = 1
  integer,      public, parameter :: RestraintsFuncDIST     = 2
  integer,      public, parameter :: RestraintsFuncDISTCOM  = 3
  integer,      public, parameter :: RestraintsFuncRMSD     = 4
  integer,      public, parameter :: RestraintsFuncRMSDCOM  = 5
  integer,      public, parameter :: RestraintsFuncANGLE    = 6
  integer,      public, parameter :: RestraintsFuncANGLECOM = 7
  integer,      public, parameter :: RestraintsFuncDIHED    = 8
  integer,      public, parameter :: RestraintsFuncDIHEDCOM = 9
  integer,      public, parameter :: RestraintsFuncPC       = 10
  integer,      public, parameter :: RestraintsFuncPCCOM    = 11
  integer,      public, parameter :: RestraintsFuncEM       = 12
  integer,      public, parameter :: RestraintsFuncREPUL    = 13
  integer,      public, parameter :: RestraintsFuncREPULCOM = 14
  integer,      public, parameter :: RestraintsFuncFB       = 15
  integer,      public, parameter :: RestraintsFuncFBCOM    = 16

  ! parameters for maximum number
  integer,      public, parameter :: RestraintsMaxConst     = 4
  integer,      public, parameter :: RestraintsMaxRef       = 2

  ! parameters for kind of direction
  integer,      public, parameter :: RestraintsDirALL       = 1
  integer,      public, parameter :: RestraintsDirX         = 2
  integer,      public, parameter :: RestraintsDirY         = 3
  integer,      public, parameter :: RestraintsDirZ         = 4

  ! restraints func type strings
  character(*), public, parameter :: RestraintsFuncTypes(16) = (/'POSI     ', &
                                                                 'DIST     ', &
                                                                 'DISTMASS ', &
                                                                 'RMSD     ', &
                                                                 'RMSDMASS ', &
                                                                 'ANGLE    ', &
                                                                 'ANGLEMASS', &
                                                                 'DIHED    ', &
                                                                 'DIHEDMASS', &
                                                                 'PC       ', &
                                                                 'PCMASS   ', &
                                                                 'EM       ', &
                                                                 'REPUL    ', &
                                                                 'REPULMASS', &
                                                                 'FB       ', &
                                                                 'FBMASS   '/)
 
  character(*), public, parameter :: RestraintsDirTypes(4) = (/'ALL ', &
                                                               'X   ', &
                                                               'Y   ', &
                                                               'Z   '/)

  ! subroutines
  public  :: init_restraints
  public  :: alloc_restraints
  public  :: dealloc_restraints
  public  :: dealloc_restraints_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_restraints
  !> @brief        initialize restraints information
  !! @authors      CK
  !! @param[out]   restraints : structure of restraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_restraints(restraints)

    ! formal arguments
    type(s_restraints),      intent(inout) :: restraints


    restraints%restraint_flag    = .false.
    restraints%verbose           = .false.
    restraints%pressure_position = .false.
    restraints%pressure_rmsd     = .false.
    restraints%nfunctions        = 0
    restraints%num_groups        = 0
    restraints%max_atoms         = 0

    return

  end subroutine init_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_restraints
  !> @brief        allocate restraints information
  !! @authors      CK
  !! @param[inout] restraints : structure of restraints information
  !! @param[in]    variable   : selected variable
  !! @param[in]    var_size   : size of the selected variable
  !! @param[in]    var_size2  : 2nd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_restraints(restraints, variable, var_size, var_size2)

    ! formal arguments
    type(s_restraints),      intent(inout) :: restraints
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer,      optional,  intent(in)    :: var_size2

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (RestraintsFunc)

      if (allocated(restraints%function)) then
        if (size(restraints%function(:)) /= var_size) &
          deallocate(restraints%function,      &
                     restraints%constant,      &
                     restraints%reference,     &
                     restraints%select_index,  &
                     restraints%direction,     &
                     restraints%mode,          &
                     restraints%exponent,      &
                     restraints%exponent_dist, &
                     restraints%weight_dist,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(restraints%function)) &
        allocate(restraints%function(var_size),      &
                 restraints%constant(var_size),      &
                 restraints%reference(var_size),     &
                 restraints%select_index(var_size),  &
                 restraints%direction(var_size),     &
                 restraints%mode(var_size),          &
                 restraints%exponent(var_size),      &
                 restraints%exponent_dist(var_size), &
                 restraints%weight_dist(var_size),   &
                 stat = alloc_stat)

      restraints%function     (1:var_size) = 0
      restraints%constant     (1:var_size) = ''
      restraints%reference    (1:var_size) = ''
      restraints%select_index (1:var_size) = ''
      restraints%mode         (1:var_size) = 0
      restraints%direction    (1:var_size) = RestraintsDirALL
      restraints%exponent     (1:var_size) = 0
      restraints%exponent_dist(1:var_size) = ''
      restraints%weight_dist  (1:var_size) = ''

    case (RestraintsGroup)

      if (allocated(restraints%group)) then
        if (size(restraints%group(:)) /= var_size) &
          deallocate(restraints%group, stat = dealloc_stat)
      end if

      if (.not. allocated(restraints%group)) &
        allocate(restraints%group(var_size), stat = alloc_stat)

      restraints%group(1:var_size) = ''

    case (RestraintsList)

      if (allocated(restraints%atomlist)) then
        if (size(restraints%num_atoms(:)) /= var_size) &
          deallocate(restraints%atomlist,  &
                     restraints%num_atoms, &
                     stat = dealloc_stat)
       end if

       if (.not. allocated(restraints%atomlist)) &
         allocate(restraints%atomlist(1:var_size2,1:var_size), &
                  restraints%num_atoms(1:var_size),            &
                  stat = alloc_stat)

       restraints%atomlist(1:var_size2,1:var_size) = 0
       restraints%num_atoms           (1:var_size) = 0

     case default

       call error_msg('Alloc_Restraints> bad variable')

     end select

     if (alloc_stat /= 0)   call error_msg_alloc
     if (dealloc_stat /= 0) call error_msg_dealloc

     return

  end subroutine alloc_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_restraints
  !> @brief        deallocate restraints information
  !! @authors      CK
  !! @param[inout] restraints : structure of restraints information
  !! @param[in]    variable   : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_restraints(restraints, variable)

    ! formatl arguments
    type(s_restraints),      intent(inout) :: restraints
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    select case (variable)

    case (RestraintsFunc)

      if (allocated(restraints%function)) then
        deallocate (restraints%function,      &
                    restraints%constant,      &
                    restraints%reference,     &
                    restraints%select_index,  &
                    restraints%mode,          &
                    restraints%direction,     &
                    restraints%exponent,      &
                    restraints%exponent_dist, &
                    restraints%weight_dist,   &
                    stat = dealloc_stat)

      end if

    case (RestraintsGroup)

      if (allocated(restraints%group)) then
        deallocate (restraints%group, &
                    stat = dealloc_stat)

      end if

    case (RestraintsList)

      if (allocated(restraints%atomlist)) then
        deallocate (restraints%atomlist,  &
                    restraints%num_atoms, &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Restraints> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_restraints_all
  !> @brief        deallocate all restraints information
  !! @authors      CK
  !! @param[inout] restraints : structure of restraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_restraints_all(restraints)

    ! formatl arguments
    type(s_restraints),      intent(inout) :: restraints


    call dealloc_restraints(restraints, RestraintsFunc)
    call dealloc_restraints(restraints, RestraintsGroup)
    call dealloc_restraints(restraints, RestraintsList)

    return

  end subroutine dealloc_restraints_all

end module sp_restraints_str_mod
