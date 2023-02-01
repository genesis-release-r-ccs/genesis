!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_pme_notable_mod
!> @brief   calculate nonbonded energy without lookup table
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_pme_notable_mod

  use sp_energy_nonbond_gpu_mod
  use sp_energy_nonbond_fugaku_mod
  use sp_energy_nonbond_intel_mod
  use sp_energy_nonbond_generic_mod
  use sp_energy_nonbond_fep_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_pme_notbl
  public  :: compute_energy_nonbond_pme_notbl_check
  public  :: compute_force_nonbond_pme_notbl
  ! FEP
  public  :: compute_energy_nonbond_pme_notbl_fep
  public  :: compute_energy_nonbond_pme_notbl_check_fep
  public  :: compute_force_nonbond_pme_notbl_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_notbl  
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !! @authors      JJ, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_notbl( &
                                              domain, enefunc, pairlist, &
                                              npt, atmcls_pbc,           &
                                              coord_pbc, force_pbc,      &
                                              virial_cell, virial,       &
                                              eelec, evdw, ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)


    select case(domain%nonbond_kernel)

    case (NBK_Generic)
      call compute_energy_nonbond_notbl_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    case (NBK_Fugaku)      
      call compute_energy_nonbond_notbl_fugaku( &
                                     domain, enefunc, pairlist,          &
                                     atmcls_pbc,                         &
                                     coord_pbc, force_pbc, virial,       &
                                     eelec, evdw)

    case (NBK_Intel)
      call compute_energy_nonbond_notbl_intel( &
                                     domain, enefunc, pairlist,          &
                                     atmcls_pbc,                         &
                                     coord_pbc, force_pbc, virial,       &
                                     eelec, evdw)

    case (NBK_GPU)
      call compute_energy_nonbond_notbl_gpu( &
                                     domain, enefunc, pairlist,          &
                                     npt, coord_pbc, force_pbc,          &
                                     virial, eelec, evdw, ene_virial)

    end select

    return

  end subroutine compute_energy_nonbond_pme_notbl  

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_notbl_check
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !! @authors      JJ, CK, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_notbl_check( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force_pbc,      &
                                                 virial_cell, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)


    call compute_energy_nonbond_notbl_chk_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    return

  end subroutine compute_energy_nonbond_pme_notbl_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_pme_notbl
  !> @brief        calculate nonbonded force without lookup table of vdw
  !! @authors      JJ, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    cpu_calc    : flag for cpu calculation or not
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force       : forces for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_pme_notbl( &
                                             domain, enefunc, pairlist, &
                                             npt, cpu_calc, atmcls_pbc, &
                                             coord_pbc,                 &
                                             force, force_pbc,          &
                                             virial_cell, virial,       &
                                             ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)


    select case(domain%nonbond_kernel)

    case (NBK_Generic)
      call compute_force_nonbond_notbl_generic( &
                                     domain, enefunc, pairlist, &
                                     coord_pbc, force_pbc, virial_cell)

    case (NBK_Fugaku)
      call compute_force_nonbond_notbl_fugaku( &
                                     domain, enefunc, pairlist, &
                                     npt, atmcls_pbc,           &
                                     coord_pbc, force_pbc, virial)

    case (NBK_Intel)
      call compute_force_nonbond_notbl_intel( &
                                     domain, enefunc, pairlist, &
                                     npt, atmcls_pbc,           &
                                     coord_pbc, force_pbc, virial)

    case (NBK_GPU)
      call compute_force_nonbond_notbl_gpu( &
                                     domain, enefunc, pairlist, &
                                     npt, cpu_calc, coord_pbc, force, &
                                     force_pbc, virial, ene_virial)

    end select

    return

  end subroutine compute_force_nonbond_pme_notbl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_notbl_fep
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_notbl_fep( &
                                              domain, enefunc, pairlist, &
                                              npt, atmcls_pbc,           &
                                              coord_pbc, force_pbc,      &
                                              virial_cell, virial,       &
                                              eelec, evdw, ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)


    select case(domain%nonbond_kernel)

    case (NBK_Generic)
      call compute_energy_nonbond_notbl_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    case (NBK_Fugaku)      
      call compute_energy_nonbond_notbl_fugaku( &
                                     domain, enefunc, pairlist,          &
                                     atmcls_pbc,                         &
                                     coord_pbc, force_pbc, virial,       &
                                     eelec, evdw)

    case (NBK_Intel)
      call compute_energy_nonbond_notbl_intel( &
                                     domain, enefunc, pairlist,          &
                                     atmcls_pbc,                         &
                                     coord_pbc, force_pbc, virial,       &
                                     eelec, evdw)

    case (NBK_GPU)
      call compute_energy_nonbond_notbl_gpu_fep( &
                                     domain, enefunc, pairlist,          &
                                     npt, coord_pbc, force_pbc,          &
                                     virial, eelec, evdw, ene_virial)

    end select

    if (domain%nonbond_kernel /= NBK_GPU) then
      call compute_energy_nonbond_notbl_generic_fep( &
        domain, enefunc, pairlist,          &
        virial_cell,  &
        eelec, evdw)
    end if

    return

  end subroutine compute_energy_nonbond_pme_notbl_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_notbl_check_fep
  !> @brief        calculate nonbonded energy without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_notbl_check_fep( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force_pbc,      &
                                                 virial_cell, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)


    call compute_energy_nonbond_notbl_chk_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    call compute_energy_nonbond_notbl_chk_generic_fep( &
                                     domain, enefunc, pairlist,          &
                                     virial_cell,  &
                                     eelec, evdw)

    return

  end subroutine compute_energy_nonbond_pme_notbl_check_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_pme_notbl_fep
  !> @brief        calculate nonbonded force without lookup table of vdw
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    cpu_calc    : flag for cpu calculation or not
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force       : forces for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_pme_notbl_fep( &
                                             domain, enefunc, pairlist, &
                                             npt, cpu_calc, atmcls_pbc, &
                                             coord_pbc,                 &
                                             force, force_pbc,          &
                                             virial_cell, virial,       &
                                             ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    integer,                  intent(inout) :: atmcls_pbc(:)
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(wp),                 intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                 intent(inout) :: virial_cell(:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)


    select case(domain%nonbond_kernel)

    case (NBK_Generic)
      call compute_force_nonbond_notbl_generic( &
                                     domain, enefunc, pairlist, &
                                     coord_pbc, force_pbc, virial_cell)

    case (NBK_Fugaku)
      call compute_force_nonbond_notbl_fugaku( &
                                     domain, enefunc, pairlist, &
                                     npt, atmcls_pbc,           &
                                     coord_pbc, force_pbc, virial)

    case (NBK_Intel)
      call compute_force_nonbond_notbl_intel( &
                                     domain, enefunc, pairlist, &
                                     npt, atmcls_pbc,           &
                                     coord_pbc, force_pbc, virial)

    case (NBK_GPU)
      call compute_force_nonbond_notbl_gpu_fep( &
                                     domain, enefunc, pairlist, &
                                     npt, cpu_calc, coord_pbc, force, &
                                     force_pbc, virial, ene_virial)

    end select

    if (domain%nonbond_kernel /= NBK_GPU) then
      call compute_force_nonbond_notbl_generic_fep( &
        domain, enefunc, pairlist, &
        virial_cell)
    end if

    return

  end subroutine compute_force_nonbond_pme_notbl_fep

end module sp_energy_nonbond_pme_notable_mod
