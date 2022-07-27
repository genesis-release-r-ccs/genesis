!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_boundary_str_mod
!> @brief   structure of boundary
!! @authors Takaharu Mori (TM), Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_boundary_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_boundary

    integer                       :: type
    real(wp)                      :: origin_x
    real(wp)                      :: origin_y
    real(wp)                      :: origin_z

    logical                       :: use_cell_linked_list
    logical                       :: wrap_all
    integer                       :: num_cells_x
    integer                       :: num_cells_y
    integer                       :: num_cells_z
    real(wp)                      :: box_size_x
    real(wp)                      :: box_size_y
    real(wp)                      :: box_size_z
    real(wp)                      :: box_size_x_ref
    real(wp)                      :: box_size_y_ref
    real(wp)                      :: box_size_z_ref
    real(wp)                      :: cell_size_x
    real(wp)                      :: cell_size_y
    real(wp)                      :: cell_size_z
    real(wp)                      :: pairlist_grid

    integer                       :: num_cells
    integer, allocatable          :: neighbor_cells(:,:)
    integer, allocatable          :: cell_head_atom(:)
    integer, allocatable          :: cell_head_atomw(:)

    integer                       :: num_neighbor_cells
    integer, allocatable          :: neighbor_cell_common_x(:)
    integer, allocatable          :: neighbor_cell_common_y(:)
    integer, allocatable          :: neighbor_cell_common_z(:)

    ! spherical boundary condition
    logical                       :: sph_pot
    integer                       :: nfunctions
    real(wp), allocatable         :: radius(:)
    real(wp), allocatable         :: center(:,:)
    real(wp), allocatable         :: const(:)
    integer, allocatable          :: exponent(:)
    real(wp)                      :: fix_layer
    integer                       :: num_fixatm
    logical, allocatable          :: fixatm(:)

    ! currently not in use
    integer                       :: num_nospot
    integer, allocatable          :: nospotlist(:)
    integer, allocatable          :: atomlist(:)

    ! CG: neighbor cell lists
    integer                       :: num_neighbor_cells_CG_ele
    integer                       :: num_neighbor_cells_CG_126
    integer                       :: num_neighbor_cells_CG_PWMcos
    integer                       :: num_neighbor_cells_CG_DNAbp
    integer                       :: num_neighbor_cells_CG_exv
    integer, allocatable          :: neighbor_cells_CG_ele(:,:)
    integer, allocatable          :: neighbor_cells_CG_126(:,:)
    integer, allocatable          :: neighbor_cells_CG_PWMcos(:,:)
    integer, allocatable          :: neighbor_cells_CG_DNAbp(:,:)
    integer, allocatable          :: neighbor_cells_CG_exv(:,:)
    integer, allocatable          :: neighbor_cell_common_x_ele(:)
    integer, allocatable          :: neighbor_cell_common_y_ele(:)
    integer, allocatable          :: neighbor_cell_common_z_ele(:)
    integer, allocatable          :: neighbor_cell_common_x_126(:)
    integer, allocatable          :: neighbor_cell_common_y_126(:)
    integer, allocatable          :: neighbor_cell_common_z_126(:)
    integer, allocatable          :: neighbor_cell_common_x_PWMcos(:)
    integer, allocatable          :: neighbor_cell_common_y_PWMcos(:)
    integer, allocatable          :: neighbor_cell_common_z_PWMcos(:)
    integer, allocatable          :: neighbor_cell_common_x_DNAbp(:)
    integer, allocatable          :: neighbor_cell_common_y_DNAbp(:)
    integer, allocatable          :: neighbor_cell_common_z_DNAbp(:)
    integer, allocatable          :: neighbor_cell_common_x_exv(:)
    integer, allocatable          :: neighbor_cell_common_y_exv(:)
    integer, allocatable          :: neighbor_cell_common_z_exv(:)
    logical                       :: calc_local_pbc

  end type s_boundary

  ! parameters for allocatable variables
  integer,      public, parameter :: BoundaryCells         = 1
  integer,      public, parameter :: BoundarySphericalPot  = 2
  integer,      public, parameter :: BoundaryCellsCGele    = 3
  integer,      public, parameter :: BoundaryCellsCG126    = 4
  integer,      public, parameter :: BoundaryCellsCGPWMcos = 5
  integer,      public, parameter :: BoundaryCellsCGDNAbp  = 6
  integer,      public, parameter :: BoundaryCellsCGexv    = 7

  ! parameters
  integer,      public, parameter :: BoundaryTypeNOBC     = 1
  integer,      public, parameter :: BoundaryTypeNOBC1    = 2
  integer,      public, parameter :: BoundaryTypePBC      = 3

  character(*), public, parameter :: BoundaryTypeTypes(3) = (/'NOBC  ', &
                                                              'NOBC1 ', &
                                                              'PBC   '/)

  ! subroutines
  public :: init_boundary
  public :: alloc_boundary
  public :: dealloc_boundary
  public :: dealloc_boundary_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_boundary
  !> @brief        initialize boundary conditions information
  !! @authors      NT
  !! @param[out]   boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_boundary(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    boundary%type                 = 0
    boundary%origin_x             = 0.0_wp
    boundary%origin_y             = 0.0_wp
    boundary%origin_z             = 0.0_wp
    boundary%use_cell_linked_list = .false.
    boundary%wrap_all             = .false.
    boundary%num_cells_x          = 0
    boundary%num_cells_y          = 0
    boundary%num_cells_z          = 0
    boundary%box_size_x           = 0.0_wp
    boundary%box_size_y           = 0.0_wp
    boundary%box_size_z           = 0.0_wp
    boundary%box_size_x_ref       = 0.0_wp
    boundary%box_size_y_ref       = 0.0_wp
    boundary%box_size_z_ref       = 0.0_wp
    boundary%cell_size_x          = 0.0_wp
    boundary%cell_size_y          = 0.0_wp
    boundary%cell_size_z          = 0.0_wp
    boundary%num_cells            = 0
    boundary%pairlist_grid        = 0.0_wp
    boundary%sph_pot              = .false.
    boundary%nfunctions           = 0
    boundary%num_nospot           = 0
    boundary%fix_layer            = 0.0_wp
    boundary%calc_local_pbc       = .false.

    return

  end subroutine init_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_boundary
  !> @brief        allocate boundary conditions information
  !! @authors      TM, KY
  !! @param[inout] boundary : boundary conditions information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_boundary(boundary, variable, var_size)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(BoundaryCells)

      if (allocated(boundary%neighbor_cell_common_x)) then

        if (size(boundary%neighbor_cell_common_x) <        &
            boundary%num_neighbor_cells ) then
          deallocate(boundary%neighbor_cells,              &
                     boundary%cell_head_atom,              &
                     boundary%cell_head_atomw,             &
                     boundary%neighbor_cell_common_x,      &
                     boundary%neighbor_cell_common_y,      &
                     boundary%neighbor_cell_common_z,      &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells)) then

        if (size(boundary%neighbor_cells(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells,  &
                     boundary%cell_head_atom,  &
                     boundary%cell_head_atomw, &
                     stat = dealloc_stat)
           if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells)) then

        allocate(boundary%neighbor_cells(boundary%num_neighbor_cells,var_size),&
                 boundary%cell_head_atom(var_size),                            &
                 boundary%cell_head_atomw(var_size),                           &
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x)) then

        allocate(boundary%neighbor_cell_common_x(boundary%num_neighbor_cells), &
                 boundary%neighbor_cell_common_y(boundary%num_neighbor_cells), &
                 boundary%neighbor_cell_common_z(boundary%num_neighbor_cells), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells(1:boundary%num_neighbor_cells,1:var_size) = 0
      boundary%neighbor_cell_common_x(1:boundary%num_neighbor_cells)    = 0
      boundary%neighbor_cell_common_y(1:boundary%num_neighbor_cells)    = 0
      boundary%neighbor_cell_common_z(1:boundary%num_neighbor_cells)    = 0
      boundary%cell_head_atom (1:var_size) = 0
      boundary%cell_head_atomw(1:var_size) = 0

    case(BoundaryCellsCGele)

      ! ======
      ! CG ELE
      ! ======
      ! 
      if (allocated(boundary%neighbor_cell_common_x_ele)) then

        if (size(boundary%neighbor_cell_common_x_ele) <   &
            boundary%num_neighbor_cells_CG_ele ) then
          deallocate(boundary%neighbor_cells_CG_ele,  &
                     boundary%neighbor_cell_common_x_ele, &
                     boundary%neighbor_cell_common_y_ele, &
                     boundary%neighbor_cell_common_z_ele, &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells_CG_ele)) then

        if (size(boundary%neighbor_cells_CG_ele(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells_CG_ele,  &
                     stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells_CG_ele)) then

        allocate(boundary%neighbor_cells_CG_ele(boundary%num_neighbor_cells_CG_ele,var_size),&
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x_ele)) then

        allocate(boundary%neighbor_cell_common_x_ele(boundary%num_neighbor_cells_CG_ele), &
                 boundary%neighbor_cell_common_y_ele(boundary%num_neighbor_cells_CG_ele), &
                 boundary%neighbor_cell_common_z_ele(boundary%num_neighbor_cells_CG_ele), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells_CG_ele(1: boundary%num_neighbor_cells_CG_ele,1:var_size) = 0
      boundary%neighbor_cell_common_x_ele(1:boundary%num_neighbor_cells_CG_ele)        = 0
      boundary%neighbor_cell_common_y_ele(1:boundary%num_neighbor_cells_CG_ele)        = 0
      boundary%neighbor_cell_common_z_ele(1:boundary%num_neighbor_cells_CG_ele)        = 0

    case(BoundaryCellsCG126)

      ! ======
      ! CG 126
      ! ======
      ! 
      if (allocated(boundary%neighbor_cell_common_x_126)) then

        if (size(boundary%neighbor_cell_common_x_126) <   &
            boundary%num_neighbor_cells_CG_126 ) then
          deallocate(boundary%neighbor_cells_CG_126,  &
                     boundary%neighbor_cell_common_x_126, &
                     boundary%neighbor_cell_common_y_126, &
                     boundary%neighbor_cell_common_z_126, &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells_CG_126)) then

        if (size(boundary%neighbor_cells_CG_126(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells_CG_126,  &
                     stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells_CG_126)) then

        allocate(boundary%neighbor_cells_CG_126(boundary%num_neighbor_cells_CG_126,var_size),&
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x_126)) then

        allocate(boundary%neighbor_cell_common_x_126(boundary%num_neighbor_cells_CG_126), &
                 boundary%neighbor_cell_common_y_126(boundary%num_neighbor_cells_CG_126), &
                 boundary%neighbor_cell_common_z_126(boundary%num_neighbor_cells_CG_126), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells_CG_126(1: boundary%num_neighbor_cells_CG_126,1:var_size) = 0
      boundary%neighbor_cell_common_x_126(1:boundary%num_neighbor_cells_CG_126)        = 0
      boundary%neighbor_cell_common_y_126(1:boundary%num_neighbor_cells_CG_126)        = 0
      boundary%neighbor_cell_common_z_126(1:boundary%num_neighbor_cells_CG_126)        = 0

    case(BoundaryCellsCGPWMcos)

      ! =========
      ! CG PWMCOS
      ! =========
      ! 
      if (allocated(boundary%neighbor_cell_common_x_PWMcos)) then

        if (size(boundary%neighbor_cell_common_x_PWMcos) <   &
            boundary%num_neighbor_cells_CG_PWMcos ) then
          deallocate(boundary%neighbor_cells_CG_PWMcos,  &
                     boundary%neighbor_cell_common_x_PWMcos, &
                     boundary%neighbor_cell_common_y_PWMcos, &
                     boundary%neighbor_cell_common_z_PWMcos, &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells_CG_PWMcos)) then

        if (size(boundary%neighbor_cells_CG_PWMcos(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells_CG_PWMcos,  &
                     stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells_CG_PWMcos)) then

        allocate(boundary%neighbor_cells_CG_PWMcos(boundary%num_neighbor_cells_CG_PWMcos,var_size),&
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x_PWMcos)) then

        allocate(boundary%neighbor_cell_common_x_PWMcos(boundary%num_neighbor_cells_CG_PWMcos), &
                 boundary%neighbor_cell_common_y_PWMcos(boundary%num_neighbor_cells_CG_PWMcos), &
                 boundary%neighbor_cell_common_z_PWMcos(boundary%num_neighbor_cells_CG_PWMcos), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells_CG_PWMcos(1: boundary%num_neighbor_cells_CG_PWMcos,1:var_size) = 0
      boundary%neighbor_cell_common_x_PWMcos(1:boundary%num_neighbor_cells_CG_PWMcos)        = 0
      boundary%neighbor_cell_common_y_PWMcos(1:boundary%num_neighbor_cells_CG_PWMcos)        = 0
      boundary%neighbor_cell_common_z_PWMcos(1:boundary%num_neighbor_cells_CG_PWMcos)        = 0

    case(BoundaryCellsCGDNAbp)

      ! ======
      ! CG DNABP
      ! ======
      ! 
      if (allocated(boundary%neighbor_cell_common_x_DNAbp)) then

        if (size(boundary%neighbor_cell_common_x_DNAbp) <   &
            boundary%num_neighbor_cells_CG_DNAbp ) then
          deallocate(boundary%neighbor_cells_CG_DNAbp,  &
                     boundary%neighbor_cell_common_x_DNAbp, &
                     boundary%neighbor_cell_common_y_DNAbp, &
                     boundary%neighbor_cell_common_z_DNAbp, &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells_CG_DNAbp)) then

        if (size(boundary%neighbor_cells_CG_DNAbp(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells_CG_DNAbp,  &
                     stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells_CG_DNAbp)) then

        allocate(boundary%neighbor_cells_CG_DNAbp(boundary%num_neighbor_cells_CG_DNAbp,var_size),&
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x_DNAbp)) then

        allocate(boundary%neighbor_cell_common_x_DNAbp(boundary%num_neighbor_cells_CG_DNAbp), &
                 boundary%neighbor_cell_common_y_DNAbp(boundary%num_neighbor_cells_CG_DNAbp), &
                 boundary%neighbor_cell_common_z_DNAbp(boundary%num_neighbor_cells_CG_DNAbp), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells_CG_DNAbp(1: boundary%num_neighbor_cells_CG_DNAbp,1:var_size) = 0
      boundary%neighbor_cell_common_x_DNAbp(1:boundary%num_neighbor_cells_CG_DNAbp)        = 0
      boundary%neighbor_cell_common_y_DNAbp(1:boundary%num_neighbor_cells_CG_DNAbp)        = 0
      boundary%neighbor_cell_common_z_DNAbp(1:boundary%num_neighbor_cells_CG_DNAbp)        = 0

    case(BoundaryCellsCGexv)

      ! ======
      ! CG EXV
      ! ======
      ! 
      if (allocated(boundary%neighbor_cell_common_x_exv)) then

        if (size(boundary%neighbor_cell_common_x_exv) <   &
            boundary%num_neighbor_cells_CG_exv ) then
          deallocate(boundary%neighbor_cells_CG_exv,  &
                     boundary%neighbor_cell_common_x_exv, &
                     boundary%neighbor_cell_common_y_exv, &
                     boundary%neighbor_cell_common_z_exv, &
                     stat = dealloc_stat)
        end if
      end if

      if (allocated(boundary%neighbor_cells_CG_exv)) then

        if (size(boundary%neighbor_cells_CG_exv(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells_CG_exv,  &
                     stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
      end if

      if (.not. allocated(boundary%neighbor_cells_CG_exv)) then

        allocate(boundary%neighbor_cells_CG_exv(boundary%num_neighbor_cells_CG_exv,var_size),&
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      end if

      if (.not. allocated(boundary%neighbor_cell_common_x_exv)) then

        allocate(boundary%neighbor_cell_common_x_exv(boundary%num_neighbor_cells_CG_exv), &
                 boundary%neighbor_cell_common_y_exv(boundary%num_neighbor_cells_CG_exv), &
                 boundary%neighbor_cell_common_z_exv(boundary%num_neighbor_cells_CG_exv), &
                 stat = alloc_stat)
      end if

      boundary%neighbor_cells_CG_exv(1: boundary%num_neighbor_cells_CG_exv,1:var_size) = 0
      boundary%neighbor_cell_common_x_exv(1:boundary%num_neighbor_cells_CG_exv)        = 0
      boundary%neighbor_cell_common_y_exv(1:boundary%num_neighbor_cells_CG_exv)        = 0
      boundary%neighbor_cell_common_z_exv(1:boundary%num_neighbor_cells_CG_exv)        = 0

    case (BoundarySphericalPot)

      if (allocated(boundary%radius)) then
        if (size(boundary%radius) == var_size) return
        
        deallocate(boundary%radius, &
                   boundary%center, &
                   boundary%const,  &
                   boundary%exponent, &
                   stat = dealloc_stat)
         if (dealloc_stat /= 0) call error_msg_dealloc
      end if

      allocate(boundary%radius(var_size),    &
               boundary%center(3, var_size), &
               boundary%const(var_size),     &
               boundary%exponent(var_size),  &
               stat = alloc_stat)

    case default

      call error_msg('Alloc_Boundary> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_boundary
  !> @brief        deallocate boundary conditions information
  !! @authors      TM, KY
  !! @param[inout] boundary : boundary conditions information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_boundary(boundary, variable)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(BoundaryCells)

      if (allocated(boundary%neighbor_cells)) then
        deallocate (boundary%neighbor_cells,  &
                    boundary%neighbor_cell_common_x, &
                    boundary%neighbor_cell_common_y, &
                    boundary%neighbor_cell_common_z, &
                    boundary%cell_head_atom,  &
                    boundary%cell_head_atomw, &
                    stat = dealloc_stat)
      end if

    case (BoundarySphericalPot)

      if (allocated(boundary%radius)) then
        deallocate(boundary%radius, &
                   boundary%center, &
                   boundary%const,  &
                   boundary%exponent, &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Boundary> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_boundary_all
  !> @brief        deallocate all boundary conditions information
  !! @authors      TM, KY
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_boundary_all(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    call dealloc_boundary(boundary, BoundaryCells)
    call dealloc_boundary(boundary, BoundarySphericalPot)

    return

  end subroutine dealloc_boundary_all

end module at_boundary_str_mod
