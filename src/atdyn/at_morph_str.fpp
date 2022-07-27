!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_morph_str_mod
!> @brief   structure of replica path information
!! @authors Yasuaki Komuro (YK), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_morph_str_mod

  use at_output_str_mod
  use molecules_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_morph
    integer                           :: method 
    integer                           :: ncycles
    integer                           :: iterations
    integer                           :: iterations_2nd
    integer                           :: eneout_period
    integer                           :: crdout_period
    integer                           :: velout_period
    integer                           :: rstout_period
    integer                           :: stoptr_period
    integer                           :: nbupdate_period
    real(wp)                          :: morph_min_rmsd
    real(wp)                          :: morph_coef
    real(wp)                          :: morph_rmsd_prev
    real(wp)                          :: morph_rmsd
    real(wp)                          :: morph_drmsd
    real(wp)                          :: morph_step_ratio    = 0.5_wp
    real(wp)                          :: decrease = 0.9_wp
    real(wp)                          :: increase = 1.1_wp
    real(wp)                          :: minimize_cutoff = 1.0_wp-3
    real(wp)                          :: delta_r
    integer        ,allocatable       :: list_rmsd(:)
    integer        ,allocatable       :: morph_atom_list(:)
    real(wp)       ,allocatable       :: crd_mov(:,:)
    real(wp)       ,allocatable       :: crd_rmov(:,:)
    real(wp)       ,allocatable       :: crd_rmov_fit(:,:)
    real(wp)       ,allocatable       :: rmass(:)
    logical                           :: verbose 
    real(wp),       allocatable       :: vec(:)
    real(wp),       allocatable       :: upper(:)
    real(wp),       allocatable       :: lower(:)
    real(wp),       allocatable       :: gradient(:)
    real(wp),       allocatable       :: work_lbfgs(:)
    integer,        allocatable       :: list_bound(:)
    integer,        allocatable       :: iwork_lbfgs(:)
  end type s_morph

  ! parameters
  integer,      public, parameter :: MorphMinMethodSD     = 1
  integer,      public, parameter :: MorphMinMethodLBFGS  = 2
  integer,      public, parameter :: MorphMinMethodSTEPS  = 3

  integer,      public, parameter :: MorphBackbone   = 1
  integer,      public, parameter :: MorphSidechain  = 2
  integer,      public, parameter :: MorphOther      = 3
  
  character(*), public, parameter :: MorphMinMethodTypes(3) = (/'SD   ', &
                                                                'LBFGS', &
                                                                'STEPS'/)

  integer, public, parameter :: MorphListRMSD   = 1
  integer, public, parameter :: MorphLBFGS      = 2
  integer, public, parameter :: MorphListBBSC   = 3

  public  :: init_morph
  public  :: alloc_morph
  public  :: dealloc_morph

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_morph
  !> @brief        initialize morphing information
  !! @authors      CK
  !! @param[out]   morph : structure of restraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_morph(morph)

    ! formal arguments
    type(s_morph),      intent(inout) :: morph


    morph%morph_rmsd_prev     = 0.0_wp
    morph%morph_rmsd          = 0.0_wp
    morph%morph_drmsd         = 0.0_wp
    morph%delta_r             = 0.0_wp
    morph%verbose             = .false.
    morph%method              = 0

    return

  end subroutine init_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_morph
  !> @brief        allocate morph information
  !! @authors      CK
  !! @param[inout] morph : structure of restraints information
  !! @param[in]    variable   : selected variable
  !! @param[in]    var_size   : size of the selected variable
  !! @param[in]    var_size2  : size of the selected variable (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_morph(morph, variable, var_size, var_size2)

    ! formal arguments
    type(s_morph),           intent(inout) :: morph
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

    case (MorphListRMSD)

      if (allocated(morph%list_rmsd)) then
        if (size(morph%list_rmsd(:)) == var_size) return
          deallocate(morph%list_rmsd,            &
                     morph%crd_mov,              &
                     morph%crd_rmov,             &
                     morph%crd_rmov_fit,         &
                     morph%rmass,                &
                     stat = dealloc_stat)
      end if

      allocate(morph%list_rmsd(var_size),      &
               morph%crd_mov(1:3,var_size),         &
               morph%crd_rmov(1:3,var_size),        &
               morph%crd_rmov_fit(1:3,var_size),    &
               morph%rmass(var_size),               &
               stat = alloc_stat)

      morph%list_rmsd        (1:var_size) = 0
      morph%crd_mov      (1:3,1:var_size) = 0.0_wp
      morph%crd_rmov     (1:3,1:var_size) = 0.0_wp
      morph%crd_rmov_fit (1:3,1:var_size) = 0.0_wp
      morph%rmass            (1:var_size) = 0.0_wp

     case default

       call error_msg('Alloc_Morph> bad variable')

    case (MorphLBFGS)

      if (allocated(morph%vec)) then
        if (size(morph%vec(:)) == var_size) return
        deallocate(morph%vec,                  &
                   morph%upper,                &
                   morph%lower,                &
                   morph%gradient,             &
                   morph%work_lbfgs,           &
                   morph%list_bound,           &
                   morph%iwork_lbfgs,          &
                   stat = dealloc_stat)
      end if

      allocate(morph%vec(var_size),              &
               morph%upper(var_size),            &
               morph%lower(var_size),            &
               morph%gradient(var_size),         &
               morph%work_lbfgs(var_size2),      &
               morph%iwork_lbfgs(var_size*3),    &
               morph%list_bound(var_size),       &
               stat = alloc_stat)

      morph%vec        (1:var_size) = 0.0_wp
      morph%upper      (1:var_size) = 0.0_wp
      morph%lower      (1:var_size) = 0.0_wp
      morph%gradient   (1:var_size) = 0.0_wp
      morph%work_lbfgs (1:var_size2) = 0.0_wp
      morph%iwork_lbfgs(1:var_size*3) = 0
      morph%list_bound(1:var_size) = 0

    case (MorphListBBSC)

      if (allocated(morph%morph_atom_list)) then
        if (size(morph%morph_atom_list(:)) == var_size) return
          deallocate(morph%morph_atom_list,            &
                     stat = dealloc_stat)
      end if

      allocate(morph%morph_atom_list(var_size),      &
               stat = alloc_stat)

      morph%morph_atom_list   (1:var_size) = 0

     end select

     if (alloc_stat /= 0)   call error_msg_alloc
     if (dealloc_stat /= 0) call error_msg_dealloc

     return

  end subroutine alloc_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_morph
  !> @brief        deallocate morph information
  !! @authors      CK
  !! @param[inout] morph : structure of restraints information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_morph(morph)

    ! formal arguments
    type(s_morph),           intent(inout) :: morph

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    if (allocated(morph%list_rmsd)) then
      deallocate(morph%list_rmsd,            &
                 morph%crd_mov,              &
                 morph%crd_rmov,             &
                 morph%crd_rmov_fit,         &
                 morph%rmass,                &
                 stat = dealloc_stat)
    end if

    if (allocated(morph%vec)) then
      deallocate(morph%vec,                  &
                 morph%upper,                &
                 morph%lower,                &
                 morph%gradient,             &
                 morph%work_lbfgs,           &
                 morph%list_bound,           &
                 morph%iwork_lbfgs,          &
                 stat = dealloc_stat)
    end if

    if (allocated(morph%morph_atom_list)) then
      deallocate(morph%morph_atom_list,            &
                 stat = dealloc_stat)
    end if

     if (dealloc_stat /= 0) call error_msg_dealloc

     return

  end subroutine dealloc_morph

end module at_morph_str_mod
