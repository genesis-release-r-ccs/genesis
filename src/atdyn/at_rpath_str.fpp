!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_rpath_str_mod
!> @brief   structure of replica path information
!! @authors Yasuaki Komuro (YK), Takaharu Mori (TM)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_rpath_str_mod

  use at_output_str_mod
  use molecules_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rpath
  ! General variables
    ! Input parameters
    integer                           :: rpathmode
    integer                           :: ncycle
    integer                           :: nreplica
    real(wp)                          :: delta
    logical                           :: fix_terminal

    ! module private
    integer                           :: dimension
    real(dp),             allocatable :: force(:)
    real(dp),             allocatable :: before_gather(:)
    real(dp),             allocatable :: after_gather(:)
    logical                           :: equilibration_only

  ! MFEP variables
    ! Input
    logical                           :: avoid_shrinkage
    integer                           :: rpath_period
    real(wp)                          :: smooth
    logical                           :: use_restart
    real(wp)                          :: sum_distance
    real(wp)                          :: distance_prev
    real(wp)                          :: distance_init
    integer,              allocatable :: rest_function(:)

    ! module private
    real(wp),             allocatable :: rest_constants(:,:,:)
    real(wp),             allocatable :: rest_reference(:,:,:)
    real(dp),             allocatable :: metric(:,:)
    real(dp),             allocatable :: rest_reference_prev(:,:)
    real(dp),             allocatable :: rest_reference_init(:,:)
    integer                           :: fitting_method
    character(256)                    :: fitting_atom
    type(s_output)                    :: output

  ! MEP variables
    ! Input parameters
    logical                           :: mep_partial_opt
    integer                           :: eneout_period
    integer                           :: crdout_period
    integer                           :: rstout_period
    integer                           :: method

    ! String
    logical                           :: massWeightCoord
    real(wp)                          :: tol_energy
    real(wp)                          :: tol_path

    ! NEB
    real(wp)                          :: tol_rmsg
    real(wp)                          :: tol_maxg
    real(wp)                          :: k_spring
    logical                           :: climbing_image
    real(wp)                          :: tol_rmsg_cineb
    integer                           :: ncorrection
    logical                           :: lbfgs_bnd
    logical                           :: lbfgs_bnd_qmonly
    real(wp)                          :: lbfgs_bnd_maxmove
    logical                           :: verbose

    ! to be deprecated
    character(256)                    :: optimizer
    real(wp)                          :: force_scale_init
    real(wp)                          :: force_scale_max 
    ! to be deprecated

    ! module private
    real(wp)                          :: energy
    real(wp)                          :: energy_prev
    real(wp)                          :: pathlength
    real(wp)                          :: pathlength_prev
    integer                           :: mep_natoms
    integer                           :: neb_cycle = 0                  
    integer,              allocatable :: isave(:,:)             
    integer,              allocatable :: mepatom_id(:)
    integer,              allocatable :: iwa(:,:)               
    real(wp),             allocatable :: dsave(:,:)             
    real(wp),             allocatable :: mep_coord(:,:)
    real(wp),             allocatable :: mep_force(:,:)
    real(wp),             allocatable :: mep_energy(:)
    real(wp),             allocatable :: mep_energy_prev(:)
    real(wp),             allocatable :: mep_length(:)
    real(wp),             allocatable :: recv_buff(:,:)
    real(wp),             allocatable :: micro_force(:,:,:)     
    real(wp),             allocatable :: wa(:,:)                           
    logical                           :: opt_micro
    logical                           :: do_cineb
    logical                           :: eneout
    logical                           :: crdout
    logical                           :: rstout
    logical                           :: first_iter     = .true.  
    logical                           :: neb_output     = .true.  
    logical                           :: first_replica  = .true.  
    logical                           :: mep_exit       = .false.    
    logical,              allocatable :: lsave(:,:)               
    character(60)                     :: task                     
    character(60),        allocatable :: csave(:)                 

    real(wp),             allocatable :: mepmd_qmforce(:,:)

  ! FEP variables
    ! Input parameters
    integer                           :: fep_period
    logical                           :: esp_energy
    logical                           :: esp_md

    ! module private
    integer                           :: num_fep = 0
    real(wp),             allocatable :: qm_energy(:)
    real(wp),             allocatable :: qm_charge(:,:)
    real(wp)                          :: pfunc_f, pfunc_b
    real(wp)                          :: dfene_f, dfene_b

  end type s_rpath

  ! parameters for allocatable variables
  integer, public, parameter      :: RpathReplicas  = 1
  integer, public, parameter      :: RpathUmbrellas = 2

  integer, public, parameter      :: RpathFitNO            = 1
  integer, public, parameter      :: RpathFitTR_ROT        = 2
  integer, public, parameter      :: RpathFitXYTR_ZROT     = 3

  character(*), public, parameter :: RpathFitTypes(3)  = (/'NO         ', &
                                                           'TR+ROT     ', &
                                                           'XYTR+ZROT  '/)

  integer,      public, parameter :: RpathmodeMFEP  = 1
  integer,      public, parameter :: RpathmodeMEP   = 2
  integer,      public, parameter :: RpathmodeMEPMD = 3
  integer,      public, parameter :: RpathmodeFEP   = 4
  
  character(*), public, parameter :: RpathmodeTypes(4) = (/'MFEP  ', &
                                                           'MEP   ', &
                                                           'MEP/MD', &
                                                           'FEP   '/)

  integer,      public, parameter :: MEPmethod_String   = 1
  integer,      public, parameter :: MEPmethod_NEB      = 2

  character(*), public, parameter :: MEPmethodTypes(2) = (/'STRING', 'NEB   '/)

  ! subroutines
  public :: alloc_rpath
  public :: alloc_rpath_mep
  public :: alloc_rpath_mepmd
  public :: alloc_rpath_fep
  public :: dealloc_rpath
  public :: dealloc_rpath_mep
  public :: dealloc_rpath_mepmd
  public :: dealloc_rpath_fep
  public :: dealloc_rpath_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rpath
  !> @brief        allocate rpath information
  !! @authors      TM
  !! @param[inout] rpath     : information of rpath
  !! @param[in]    variable  : allocatable variable
  !! @param[in]    var_size1 : size of variables (REMD dimension)
  !! @param[in]    var_size2 : size of variables (Total number of replicas)
  !! @param[in]    var_size3 : size of variables (max number of replica)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rpath(rpath, variable, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size1
    integer,      optional,  intent(in)    :: var_size2
    integer,      optional,  intent(in)    :: var_size3

    ! local variables
    integer                  :: alloc_stat, dealloc_stat


    alloc_stat = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(RpathReplicas)
      if (allocated(rpath%rest_function)) then
        if (size(rpath%rest_function(:)) == var_size1) return
        deallocate(rpath%rest_function,       &
                   rpath%force,               &
                   rpath%metric,              &
                   rpath%before_gather,       &
                   rpath%after_gather,        &
                   stat = dealloc_stat)
      end if

      allocate(rpath%rest_function(var_size1),             &
               rpath%force(var_size1),                     &
               rpath%metric(var_size1,var_size1),          &
               rpath%before_gather(var_size1),             &
               rpath%after_gather(var_size1*var_size2),    &
               stat = alloc_stat)

      rpath%rest_function(1:var_size1)                 = 0
      rpath%force(1:var_size1)                         = 0.0_wp
      rpath%metric(1:var_size1,1:var_size1)            = 0.0_wp
      rpath%before_gather(1:var_size1)                 = 0.0_wp
      rpath%after_gather(1:var_size1*var_size2)        = 0.0_wp

    case(RpathUmbrellas)

      if (allocated(rpath%rest_constants)) then
        if (size(rpath%rest_constants(:,:,:)) == (4*var_size1*var_size2)) return
        deallocate(rpath%rest_constants,          &
                   rpath%rest_reference,          &
                   rpath%rest_reference_prev,     &
                   rpath%rest_reference_init,     &
                   stat = dealloc_stat)
      end if

      allocate(rpath%rest_constants(4,var_size1,var_size2), &
               rpath%rest_reference(2,var_size1,var_size2), &
               rpath%rest_reference_prev(var_size1,var_size2), &
               rpath%rest_reference_init(var_size1,var_size2), &
               stat = alloc_stat)

      rpath%rest_constants(1:4,1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference(1:2,1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference_prev(1:var_size1,1:var_size2) = 0.0_wp
      rpath%rest_reference_init(1:var_size1,1:var_size2) = 0.0_wp

    case default

      call error_msg('Alloc_Rpath> bad variable')

    end select

    if (alloc_stat /= 0) call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rpath_mep
  !> @brief        allocate rpath information
  !! @authors      YA, KY, SI
  !! @param[inout] rpath     : information of rpath
  !! @param[in]    var_size1 : MEP dimension
  !! @param[in]    var_size2 : Number of replicas
  !! @param[in]    var_size3 : Number of QM atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rpath_mep(rpath, var_size1, var_size2, var_size3, &
                             var_size4, var_size5)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: var_size1 ! rpath%mep_natoms * 3
    integer,                 intent(in)    :: var_size2 ! rpath%nreplica
    integer,                 intent(in)    :: var_size3 ! minimize%num_optatoms_micro
    integer,                 intent(in)    :: var_size4 ! nrep_per_proc   
    integer,      optional,  intent(in)    :: var_size5 ! qmmm%qm_natoms

    ! local variables
    integer                  :: alloc_stat
    integer                  :: tmp_var1, tmp_var2, tmp_var3, tmp_var4


    alloc_stat = 0

    ! allocate selected variables
    !
    ! Avoid error by re-allocation 
    if (.not. allocated(rpath%mep_coord)) then 
      allocate(rpath%force(var_size1),                        &
               rpath%mep_coord(var_size1, var_size2),         &
               rpath%mep_force(var_size1, var_size2),         &
               rpath%recv_buff(var_size1, var_size2),         &
               rpath%mep_energy(var_size2),                   &
               rpath%mep_energy_prev(var_size2),              &
               rpath%mep_length(var_size2),                   &
               rpath%before_gather(var_size1),                &
               rpath%after_gather(var_size1*var_size2),       &
               rpath%micro_force(3, var_size3, var_size4),    &
               stat = alloc_stat)
      if ((rpath%method == MEPmethod_NEB) .and. &
          (var_size4 > 1)) then
        allocate(rpath%isave(44, var_size2))
        allocate(rpath%dsave(29, var_size2))
        allocate(rpath%lsave(4,  var_size2))
        allocate(rpath%csave(var_size2))
        allocate(rpath%iwa(3*var_size1*var_size2, var_size2))
        tmp_var1 = rpath%ncorrection
        tmp_var2 = var_size1
        tmp_var3 = var_size1 * var_size2
        tmp_var4 = (2*tmp_var3 + 11*tmp_var1 +8) * tmp_var1 + 5*tmp_var2
        allocate(rpath%wa(tmp_var4, var_size2))
      end if

      rpath%force(1:var_size1)                  = 0.0_wp
      rpath%mep_coord(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%mep_force(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%recv_buff(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%mep_energy(1:var_size2)             = 0.0_wp
      rpath%mep_energy_prev(1:var_size2)        = 0.0_wp
      rpath%mep_length(1:var_size2)             = 0.0_wp
      rpath%before_gather(1:var_size1)          = 0.0_wp
      rpath%after_gather(1:var_size1*var_size2) = 0.0_wp
      rpath%micro_force(1:3, 1:var_size3, 1:var_size4) =  0.0_wp
    end if  

    if (.not. allocated(rpath%qm_energy)) then
      if (present(var_size5)) then
        allocate(rpath%qm_energy(var_size2), stat = alloc_stat)
        rpath%qm_energy = 0.0_wp
      end if
    end if

    if (alloc_stat /= 0) call error_msg_alloc

    return

  end subroutine alloc_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rpath_mepmd
  !> @brief        allocate rpath information
  !! @authors      KY
  !! @param[inout] rpath     : information of rpath
  !! @param[in]    var_size1 : MEP dimension
  !! @param[in]    var_size2 : Number of replicas
  !! @param[in]    var_size5 : Number of QM atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rpath_mepmd(rpath, var_size1, var_size2, var_size5)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: var_size1 ! rpath%mep_natoms * 3
    integer,                 intent(in)    :: var_size2 ! rpath%nreplica
    integer,      optional,  intent(in)    :: var_size5 ! qmmm%qm_natoms

    ! local variables
    integer                  :: alloc_stat
    integer                  :: tmp_var1, tmp_var2, tmp_var3, tmp_var4


    alloc_stat = 0

    ! allocate selected variables
    !
    ! Avoid error by re-allocation 
    if (.not. allocated(rpath%mep_coord)) then 
      allocate(rpath%force(var_size1),                        &
               rpath%mep_coord(var_size1, var_size2),         &
               rpath%mep_force(var_size1, var_size2),         &
               rpath%recv_buff(var_size1, var_size2),         &
               rpath%mep_energy(var_size2),                   &
               rpath%mep_energy_prev(var_size2),              &
               rpath%mep_length(var_size2),                   &
               rpath%before_gather(var_size1),                &
               rpath%after_gather(var_size1*var_size2),       &
               stat = alloc_stat)
      if (rpath%method == MEPmethod_NEB) then
        allocate(rpath%isave(44, var_size2))
        allocate(rpath%dsave(29, var_size2))
        allocate(rpath%lsave(4,  var_size2))
        allocate(rpath%csave(var_size2))
        allocate(rpath%iwa(3*var_size1*var_size2, var_size2))
        tmp_var1 = rpath%ncorrection
        tmp_var2 = var_size1
        tmp_var3 = var_size1 * var_size2
        tmp_var4 = (2*tmp_var3 + 11*tmp_var1 +8) * tmp_var1 + 5*tmp_var2
        allocate(rpath%wa(tmp_var4, var_size2))
      end if

      rpath%force(1:var_size1)                  = 0.0_wp
      rpath%mep_coord(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%mep_force(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%recv_buff(1:var_size1,1:var_size2)  = 0.0_wp
      rpath%mep_energy(1:var_size2)             = 0.0_wp
      rpath%mep_energy_prev(1:var_size2)        = 0.0_wp
      rpath%mep_length(1:var_size2)             = 0.0_wp
      rpath%before_gather(1:var_size1)          = 0.0_wp
      rpath%after_gather(1:var_size1*var_size2) = 0.0_wp
    end if  

    if(.not. allocated(rpath%qm_energy)) then
      if (present(var_size5)) then
        allocate(rpath%qm_energy(var_size2), stat = alloc_stat)
        rpath%qm_energy = 0.0_wp
      end if
    end if

    if (alloc_stat /= 0) call error_msg_alloc

    return

  end subroutine alloc_rpath_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rpath_fep
  !> @brief        allocate rpath information
  !! @authors      KY
  !! @param[inout] rpath     : information of rpath
  !! @param[in]    var_size1 : MEP dimension
  !! @param[in]    var_size2 : Number of replicas
  !! @param[in]    var_size3 : Number of QM atoms (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rpath_fep(rpath, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2
    integer,      optional,  intent(in)    :: var_size3

    ! local variables
    integer                  :: alloc_stat


    alloc_stat = 0

    ! allocate selected variables
    !
    allocate(rpath%mep_coord(var_size1, var_size2),         &
             stat = alloc_stat)

    rpath%mep_coord(1:var_size1,1:var_size2) = 0.0_wp

    if (present(var_size3)) then
      allocate(                                  &
        rpath%qm_charge(var_size3, var_size2),   &
        rpath%qm_energy(var_size2),              &
        stat = alloc_stat)
      rpath%qm_charge = 0.0_wp
      rpath%qm_energy = 0.0_wp
    end if

    if (alloc_stat /= 0) call error_msg_alloc

    return

  end subroutine alloc_rpath_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath
  !> @brief        deallocate rpath information
  !! @authors      TM
  !! @param[inout] rpath     : rpath information
  !! @param[in]    variable : allocatable variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath(rpath, variable)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case (RpathReplicas)

      if (allocated(rpath%rest_function)) then
        deallocate(rpath%rest_function,       &
                   rpath%force,               &
                   rpath%metric,              &
                   rpath%before_gather,       &
                   rpath%after_gather,        &
                   stat = dealloc_stat)
      end if

    case (RpathUmbrellas)

      if (allocated(rpath%rest_constants)) then
        deallocate(rpath%rest_constants,          &
                   rpath%rest_reference,          &
                   rpath%rest_reference_prev,     &
                   rpath%rest_reference_init,     &
                   stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Rpath> bad variable')

    end select 


    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath_mep
  !> @brief        deallocate rpath information (MEP)
  !! @authors      KY
  !! @param[inout] rpath     : rpath information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath_mep(rpath)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    if (allocated(rpath%mep_coord)) then
      deallocate(rpath%force,            &
                 rpath%mepatom_id,       &
                 rpath%mep_coord,        &
                 rpath%mep_force,        &
                 rpath%recv_buff,        &
                 rpath%mep_energy,       &
                 rpath%mep_energy_prev,  &
                 rpath%mep_length,       &
                 rpath%before_gather,    &
                 rpath%after_gather,     &
                 stat = dealloc_stat)
    end if

    ! QM/MM 
    if (allocated(rpath%qm_energy)) then
      deallocate(rpath%qm_energy,        &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine dealloc_rpath_mep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath_mepmd
  !> @brief        deallocate rpath information (MEP)
  !! @authors      KY
  !! @param[inout] rpath     : rpath information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath_mepmd(rpath)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    integer                  :: dealloc_stat

    dealloc_stat = 0

    if (allocated(rpath%mep_coord)) then
      deallocate(rpath%force,            &
                 rpath%mepatom_id,       &
                 rpath%mep_coord,        &
                 rpath%mep_force,        &
                 rpath%recv_buff,        &
                 rpath%mep_energy,       &
                 rpath%mep_energy_prev,  &
                 rpath%mep_length,       &
                 rpath%before_gather,    &
                 rpath%after_gather,     &
                 stat = dealloc_stat)
    end if

    ! QM/MM 
    if (allocated(rpath%qm_energy)) then
      deallocate(rpath%qm_energy,        &
                 stat = dealloc_stat)
    end if
    if (allocated(rpath%mepmd_qmforce)) then
      deallocate(rpath%mepmd_qmforce,    &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

  end subroutine dealloc_rpath_mepmd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath_fep
  !> @brief        deallocate rpath information (FEP)
  !! @authors      KY
  !! @param[inout] rpath     : rpath information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath_fep(rpath)

    ! formal arguments
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    if (allocated(rpath%mep_coord)) then
      deallocate(rpath%mepatom_id,       &
                 rpath%mep_coord,        &
                 stat = dealloc_stat)
    end if

    ! QM/MM 
    if (allocated(rpath%qm_charge)) then
      deallocate(rpath%qm_charge,        &
                 rpath%qm_energy,        &
                 stat = dealloc_stat)
    end if

  end subroutine dealloc_rpath_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rpath_all
  !> @brief        deallocate all rpath information
  !! @authors      TM
  !! @param[inout] rpath : information of rpath
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rpath_all(rpath)

    ! formal arguments
    type(s_rpath),         intent(inout) :: rpath


    if (rpath%rpathmode == RpathmodeMFEP) then
      call dealloc_rpath(rpath, RpathReplicas)
      call dealloc_rpath(rpath, RpathUmbrellas)

    else if (rpath%rpathmode == RpathmodeMEP) then
      call dealloc_rpath_mep(rpath)

    else if (rpath%rpathmode == RpathmodeMEPMD) then
      call dealloc_rpath_mepmd(rpath)

    else if (rpath%rpathmode == RpathmodeFEP) then
      call dealloc_rpath_fep(rpath)

    end if

    return

  end subroutine dealloc_rpath_all

end module at_rpath_str_mod
