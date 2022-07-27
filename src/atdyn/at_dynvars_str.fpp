!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_dynvars_str_mod
!> @brief   structure of dynvars
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!! @note    dynvars include coord, force, energy, boxsize, pressure and so on.
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_dynvars_str_mod

  use at_energy_str_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_dynvars
    type(s_energy)                :: energy

    real(wp),         allocatable :: coord(:,:)
    real(wp),         allocatable :: trans(:,:)
    real(wp),         allocatable :: coord_pbc(:,:)
    real(wp),         allocatable :: coord_ref(:,:)
    real(wp),         allocatable :: velocity(:,:)
    real(wp),         allocatable :: velocity_ref(:,:)
    real(wp),         allocatable :: force(:,:)
    real(wp),         allocatable :: force_omp(:,:,:)
    real(wp),         allocatable :: random_force(:,:), random_force1(:,:)
    real(wp),         allocatable :: temporary(:,:)
    real(wp)                      :: virial(3,3)
    real(wp)                      :: virial_const(3,3)
    real(wp)                      :: virial_extern(3,3)

    integer                       :: step
    integer                       :: iterations
    integer                       :: istart_atom, iend_atom
    integer,          allocatable :: recv_count(:)
    integer,          allocatable :: displs(:)
    real(wp)                      :: time
    real(wp)                      :: total_pene
    real(wp)                      :: total_kene
    real(wp)                      :: total_energy
    real(wp)                      :: temperature
    real(wp)                      :: rms_gradient
    real(wp)                      :: max_gradient
    real(wp)                      :: hfc_kene
    real(wp)                      :: virial_kene

    real(wp)                      :: volume
    real(wp)                      :: internal_virial
    real(wp)                      :: internal_pressure
    real(wp)                      :: external_virial
    real(wp)                      :: external_pressure

    real(wp)                      :: thermostat_momentum
    real(wp)                      :: thermostat_momentum_ref
    real(wp)                      :: barostat_momentum(3)
    real(wp)                      :: barostat_momentum_ref(3)

    real(wp)                      :: box_size_x
    real(wp)                      :: box_size_y
    real(wp)                      :: box_size_z
    real(wp)                      :: pressure_xx
    real(wp)                      :: pressure_yy
    real(wp)                      :: pressure_zz
    logical                       :: verbose
  end type s_dynvars

  ! parameters
  integer,      public, parameter :: DynvarsGeneral  = 1
  integer,      public, parameter :: DynvarsLangevin = 2

  ! subroutines
  public :: init_dynvars
  public :: alloc_dynvars
  public :: dealloc_dynvars
  public :: dealloc_dynvars_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_dynvars
  !> @brief        initialize energy terms in dynvars
  !! @authors      YS, TM, CK
  !! @param[out]   dynvars : dynamic variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_dynvars(dynvars)

    ! formal arguments
    type(s_dynvars),         intent(inout) :: dynvars


    call init_energy (dynvars%energy)

    dynvars%virial(1:3,1:3)            = 0.0_wp
    dynvars%virial_const(1:3,1:3)      = 0.0_wp
    dynvars%virial_extern(1:3,1:3)     = 0.0_wp

    dynvars%step                       = 0
    dynvars%iterations                 = 0
    dynvars%time                       = 0.0_wp
    dynvars%total_pene                 = 0.0_wp
    dynvars%total_kene                 = 0.0_wp
    dynvars%total_energy               = 0.0_wp
    dynvars%temperature                = 0.0_wp
    dynvars%rms_gradient               = 0.0_wp
    dynvars%hfc_kene                   = 0.0_wp
    dynvars%virial_kene                = 0.0_wp

    dynvars%volume                     = 0.0_wp
    dynvars%internal_virial            = 0.0_wp
    dynvars%internal_pressure          = 0.0_wp
    dynvars%external_virial            = 0.0_wp
    dynvars%external_pressure          = 0.0_wp

    dynvars%thermostat_momentum        = 0.0_wp
    dynvars%thermostat_momentum_ref    = 0.0_wp
    dynvars%barostat_momentum(1:3)     = 0.0_wp
    dynvars%barostat_momentum_ref(1:3) = 0.0_wp

    dynvars%box_size_x                 = 0.0_wp
    dynvars%box_size_y                 = 0.0_wp
    dynvars%box_size_z                 = 0.0_wp
    dynvars%pressure_xx                = 0.0_wp
    dynvars%pressure_yy                = 0.0_wp
    dynvars%pressure_zz                = 0.0_wp
    
    dynvars%verbose                    = .false.

    return

  end subroutine init_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_dynvars
  !> @brief        allocate dynamic variables
  !! @authors      YS, TM
  !! @param[inout] dynvars  : dynamic variables
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_dynvars(dynvars, variable, var_size, var_size1)

    ! formal arguments
    type(s_dynvars),         intent(inout) :: dynvars
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size
    integer, optional,       intent(in)    :: var_size1

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(DynvarsGeneral)

      if (allocated(dynvars%coord)) then
        if (size(dynvars%coord(1,:)) == var_size) return
        deallocate(dynvars%coord,        &
                   dynvars%coord_pbc,    &
                   dynvars%force,        &
                   dynvars%force_omp,    &
                   dynvars%velocity,     &
                   dynvars%velocity_ref, &
                   dynvars%coord_ref,    &
                   dynvars%temporary,    &
                   stat = dealloc_stat)
      end if

      allocate(dynvars%coord(3,var_size),              &
               dynvars%trans(3,var_size),              &
               dynvars%coord_pbc(3,var_size),          &
               dynvars%force(3,var_size),              &
               dynvars%force_omp(3,var_size,nthread),  &
               dynvars%velocity(3,var_size),           &
               dynvars%velocity_ref(3,var_size),       &
               dynvars%coord_ref(3,var_size),          &
               dynvars%temporary(3,var_size),          &
               stat = alloc_stat)

      dynvars%coord       (1:3,1:var_size) = 0.0_wp
      dynvars%trans       (1:3,1:var_size) = 0.0_wp
      dynvars%coord_pbc   (1:3,1:var_size) = 0.0_wp
      dynvars%force       (1:3,1:var_size) = 0.0_wp
      dynvars%force_omp   (1:3,1:var_size,1:nthread) = 0.0_wp
      dynvars%velocity    (1:3,1:var_size) = 0.0_wp
      dynvars%velocity_ref(1:3,1:var_size) = 0.0_wp
      dynvars%coord_ref   (1:3,1:var_size) = 0.0_wp
      dynvars%temporary   (1:3,1:var_size) = 0.0_wp

    case(DynvarsLangevin)

      if (allocated(dynvars%random_force)) then
        if (size(dynvars%random_force(1,:)) == var_size) return
        deallocate(dynvars%random_force,         &
                   dynvars%random_force1,        &
                   stat = dealloc_stat)
      end if

      allocate(dynvars%random_force(3,var_size),              &
               dynvars%random_force1(3,var_size/var_size1+1), &
               stat = alloc_stat)

      dynvars%random_force (1:3,1:var_size) = 0.0_wp

    case default

      call error_msg('Alloc_Dynvars> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_dynvars
  !> @brief        deallocate dynvars information
  !! @authors      TM
  !! @param[inout] dynvars  : dynamic variables
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_dynvars(dynvars, variable)

    ! formal arguments
    type(s_dynvars),         intent(inout) :: dynvars
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(DynvarsGeneral)

      if (allocated(dynvars%coord)) then
        deallocate (dynvars%coord,        &
                    dynvars%coord_pbc,    &
                    dynvars%force,        &
                    dynvars%force_omp,    &
                    dynvars%velocity,     &
                    dynvars%velocity_ref, &
                    dynvars%coord_ref,    &
                    dynvars%temporary,    &
                    stat = dealloc_stat)
      end if

    case(DynvarsLangevin)

      if (allocated(dynvars%random_force)) then
        deallocate (dynvars%random_force, &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Dynvars> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_dynvars_all
  !> @brief        deallocate all dynvars information
  !! @authors      TM, YS
  !! @param[inout] dynvars : dynamic variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_dynvars_all(dynvars)

    ! formal arguments
    type(s_dynvars),         intent(inout) :: dynvars

    call dealloc_dynvars(dynvars, DynvarsGeneral )
    call dealloc_dynvars(dynvars, DynvarsLangevin)

    return

  end subroutine dealloc_dynvars_all

end module at_dynvars_str_mod
