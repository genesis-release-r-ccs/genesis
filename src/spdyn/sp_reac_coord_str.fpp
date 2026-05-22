!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_reac_coord_str_mod
!> @brief   structure of reaction coordinate information
!! @authors Junichiro Koga (JK)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_reac_coord_str_mod

  use constants_mod
  use sp_restraints_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  type, public :: s_reac_coord

    logical                       :: reac_coord_flag = .false.
    logical                       :: verbose
    integer                       :: nreac_coord
    integer                       :: num_groups
    integer                       :: max_atoms

    integer,              allocatable :: function(:)
    integer,              allocatable :: direction(:)
    character(MaxLine),   allocatable :: select_index(:)
    character(MaxLine),   allocatable :: unitv(:)
    character(MaxLine),   allocatable :: exponent_dist(:)
    character(MaxLine),   allocatable :: weight_dist(:)

    ! reac coord group (size = num_group)
    character(MaxLine),   allocatable :: group(:)

    ! reac coord list (size = num_group)
    integer,          allocatable :: atomlist(:,:)
    integer,          allocatable :: num_atoms(:)

    integer,  allocatable           :: nxi(:)
    real(wp), allocatable           :: ximin(:), ximax(:), dxi(:)
    integer                         :: monitor_period
    real(wp), allocatable           :: uboundary(:)
    logical,  allocatable           :: periodicity(:)
    logical                         :: chk_periodicity  
    
    real(wp)                        :: time_tot=0.0

  end type s_reac_coord

  integer, public, parameter :: ReacCoord = 21

  ! parameters for allocatable variables
  integer,      public, parameter :: ReacCoordsFunc         = 1
  integer,      public, parameter :: ReacCoordsGroup        = 2
  integer,      public, parameter :: ReacCoordsList         = 3

  ! parameters for kind of functions
  integer,      public, parameter :: ReacCoordsFuncPOSI     = 1
  integer,      public, parameter :: ReacCoordsFuncPOSICOM  = 2
  integer,      public, parameter :: ReacCoordsFuncDIST     = 3
  integer,      public, parameter :: ReacCoordsFuncDISTCOM  = 4
  integer,      public, parameter :: ReacCoordsFuncANGLE    = 5
  integer,      public, parameter :: ReacCoordsFuncANGLECOM = 6
  integer,      public, parameter :: ReacCoordsFuncDIHED    = 7
  integer,      public, parameter :: ReacCoordsFuncDIHEDCOM = 8

    ! parameters for kind of direction
  integer,      public, parameter :: ReacCoordsDirALL       = 1
  integer,      public, parameter :: ReacCoordsDirX         = 2
  integer,      public, parameter :: ReacCoordsDirY         = 3
  integer,      public, parameter :: ReacCoordsDirZ         = 4
  integer,      public, parameter :: ReacCoordsDirX2        = 5 
  integer,      public, parameter :: ReacCoordsDirY2        = 6 
  integer,      public, parameter :: ReacCoordsDirZ2        = 7 

  ! reaction coordinate func type strings
  character(*), public, parameter :: ReacCoordsFuncTypes(8) = (/'POSI     ', &
                                                                'POSIMASS ', &
                                                                'DIST     ', &
                                                                'DISTMASS ', &
                                                                'ANGLE    ', &
                                                                'ANGLEMASS', &
                                                                'DIHED    ', &
                                                                'DIHEDMASS'/)

  character(*), public, parameter :: ReacCoordsDirTypes(7)  = (/'ALL ', &
                                                                'X   ', &
                                                                'Y   ', &
                                                                'Z   ', &
                                                                'X2  ', &
                                                                'Y2  ', &
                                                                'Z2  '/)                                                                 
                                                                
  public  :: init_reac_coord
  public  :: alloc_reac_coord
  public  :: dealloc_reac_coord
  public  :: dealloc_reac_coord_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_reac_coord
  !> @brief        initialize reaction coordinate information
  !! @authors      JK
  !! @param[out]   reac_coord : structure of reaction coordinate information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_reac_coord(reac_coord)

    ! formal arguments
    type(s_reac_coord),      intent(inout) :: reac_coord

    reac_coord%reac_coord_flag   = .false.
    reac_coord%verbose           = .false.
    reac_coord%nreac_coord       = 0
    reac_coord%num_groups        = 0
    reac_coord%max_atoms         = 0

    return

  end subroutine init_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_reac_coord
  !> @brief        allocate reaction coordinate information
  !! @authors      JK
  !! @param[inout] reac_coord : structure of reaction coordinate information
  !! @param[in]    variable   : selected variable
  !! @param[in]    var_size   : size of the selected variable
  !! @param[in]    var_size2  : 2nd size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_reac_coord(reac_coord, variable, var_size, var_size2)
    ! formal arguments
    type(s_reac_coord),     intent(inout) :: reac_coord
    integer,                intent(in)    :: variable
    integer,                intent(in)    :: var_size
    integer,      optional, intent(in)    :: var_size2
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case (RestraintsFunc)

      if (allocated(reac_coord%function)) then
        if (size(reac_coord%function(:)) /= var_size) &
          deallocate(reac_coord%function,      &
                     reac_coord%select_index,  &
                     reac_coord%direction,     &
                     reac_coord%unitv,         &
                     reac_coord%exponent_dist, &
                     reac_coord%weight_dist,   &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(reac_coord%function)) &
        allocate(reac_coord%function(var_size),      &
                 reac_coord%select_index(var_size),  &
                 reac_coord%direction(var_size),     &
                 reac_coord%unitv(var_size),         &
                 reac_coord%exponent_dist(var_size), &
                 reac_coord%weight_dist(var_size),   &
                 stat = alloc_stat)

      reac_coord%function     (1:var_size) = 0
      reac_coord%select_index (1:var_size) = ''
      reac_coord%direction    (1:var_size) = RestraintsDirALL
      reac_coord%unitv        (1:var_size) = ''
      reac_coord%exponent_dist(1:var_size) = ''
      reac_coord%weight_dist  (1:var_size) = ''

    case (RestraintsGroup)

      if (allocated(reac_coord%group)) then
        if (size(reac_coord%group(:)) /= var_size) &
          deallocate(reac_coord%group, stat = dealloc_stat)
      end if

      if (.not. allocated(reac_coord%group)) &
        allocate(reac_coord%group(var_size), stat = alloc_stat)

      reac_coord%group(1:var_size) = ''

    case (RestraintsList)

      if (allocated(reac_coord%atomlist)) then
        if (size(reac_coord%num_atoms(:)) /= var_size) &
          deallocate(reac_coord%atomlist,  &
                     reac_coord%num_atoms, &
                     stat = dealloc_stat)
       end if

       if (.not. allocated(reac_coord%atomlist)) &
         allocate(reac_coord%atomlist(1:var_size2,1:var_size), &
                  reac_coord%num_atoms(1:var_size),            &
                  stat = alloc_stat)

       reac_coord%atomlist(1:var_size2,1:var_size) = 0
       reac_coord%num_atoms           (1:var_size) = 0

     case (ReacCoord)
      if (allocated(reac_coord%nxi)) then
        if (size(reac_coord%nxi(:)) /= var_size) &
          deallocate(reac_coord%nxi,             &
                     reac_coord%ximin,           &
                     reac_coord%ximax,           &
                     reac_coord%dxi,             &
                     reac_coord%uboundary,       &
                     reac_coord%periodicity,      &
                     stat = dealloc_stat)
      end if

      if (.not. allocated(reac_coord%nxi))   &
        allocate(reac_coord%nxi      (var_size),   &
                 reac_coord%ximin    (var_size),   &
                 reac_coord%ximax    (var_size),   &
                 reac_coord%dxi      (var_size),   &
                 reac_coord%uboundary(var_size),   &
                 reac_coord%periodicity(var_size), &
                 stat = alloc_stat)

        reac_coord%nxi      (1:var_size) = 0
        reac_coord%ximin    (1:var_size) = 0_wp
        reac_coord%ximax    (1:var_size) = 0_wp
        reac_coord%dxi      (1:var_size) = 0_wp
        reac_coord%uboundary(1:var_size) = 0_wp
        reac_coord%periodicity (1:var_size) = .false.

     case default

       call error_msg('Alloc_ReacCoord> Invalid variable')

     end select

     if (alloc_stat /= 0)   call error_msg_alloc
     if (dealloc_stat /= 0) call error_msg_dealloc

     return
  end subroutine alloc_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_reac_coord
  !> @brief        deallocate reaction coordinate information
  !! @authors      JK
  !! @param[inout] reac_coord : structure of reaction coordinate information
  !! @param[in]    variable   : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_reac_coord(reac_coord, variable)

    ! formatl arguments
    type(s_reac_coord),      intent(inout) :: reac_coord
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    select case (variable)

    case (RestraintsFunc)

      if (allocated(reac_coord%function)) then
        deallocate (reac_coord%function,      &
                    reac_coord%select_index,  &
                    reac_coord%direction,     &
                    reac_coord%unitv,         &
                    reac_coord%exponent_dist, &
                    reac_coord%weight_dist,   &
                    stat = dealloc_stat)

      end if

    case (RestraintsGroup)

      if (allocated(reac_coord%group)) then
        deallocate (reac_coord%group, &
                    stat = dealloc_stat)

      end if

    case (RestraintsList)

      if (allocated(reac_coord%atomlist)) then
        deallocate (reac_coord%atomlist,  &
                    reac_coord%num_atoms, &
                    stat = dealloc_stat)
      end if

    case (ReacCoord)
      if (allocated(reac_coord%nxi)) then
        deallocate (reac_coord%nxi,          &
                    reac_coord%ximin,        &
                    reac_coord%ximax,        &
                    reac_coord%dxi,          &
                    reac_coord%uboundary,    &
                    reac_coord%periodicity,  &
                    stat = dealloc_stat)
      endif

    case default

      call error_msg('Dealloc_Restraints> Invalid variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_reac_coord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_reac_coord_all
  !> @brief        deallocate all reaction coordinate information
  !! @authors      JK
  !! @param[inout] reac_coord : structure of reaction coordinate information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_reac_coord_all(reac_coord)

    ! formatl arguments
    type(s_reac_coord),      intent(inout) :: reac_coord


    call dealloc_reac_coord(reac_coord, RestraintsFunc)
    call dealloc_reac_coord(reac_coord, RestraintsGroup)
    call dealloc_reac_coord(reac_coord, RestraintsList)
    call dealloc_reac_coord(reac_coord, ReacCoord)

    return

  end subroutine dealloc_reac_coord_all

end module sp_reac_coord_str_mod

