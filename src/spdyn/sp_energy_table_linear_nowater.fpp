!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_nowater_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_nowater_mod

  use sp_energy_nonbond_gpu_mod
  use sp_energy_nonbond_intel_knl_mod
  use sp_energy_nonbond_k_generic_mod
  use sp_energy_nonbond_fugaku_mod
  use sp_energy_nonbond_generic_mod
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
  public  :: compute_energy_nonbond_table_linear
  public  :: compute_energy_nonbond_table_linear_check
  public  :: compute_force_nonbond_table_linear

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear( &
                                                 domain, enefunc, pairlist,   &
                                                 npt, atmcls_pbc,             &
                                                 coord_pbc, force_pbc,        &
                                                 virial_cell, virial,         &
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
      call compute_energy_nonbond_tbl_lnr_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    case (NBK_KGeneric)
      call compute_energy_nonbond_tbl_lnr_k_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    case (NBK_Fugaku)
      call compute_energy_nonbond_tbl_lnr_fugaku( &
                                     domain, enefunc, pairlist,          &
                                     atmcls_pbc,                         &
                                     coord_pbc, force_pbc, virial,       &
                                     eelec, evdw)

    case (NBK_IntelKnl)
      call compute_energy_nonbond_tbl_lnr_intel_knl( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    case (NBK_GPU)
      call compute_energy_nonbond_tbl_lnr_gpu( &
                                     domain, enefunc, pairlist,          &
                                     npt, coord_pbc, force_pbc,          &
                                     virial, eelec, evdw, ene_virial)

    end select

    return

  end subroutine compute_energy_nonbond_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_check
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ, CK, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_check( &
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


    call compute_energy_nonbond_tbl_lnr_chk_generic( &
                                     domain, enefunc, pairlist,          &
                                     coord_pbc, force_pbc, virial_cell,  &
                                     eelec, evdw)

    return

  end subroutine compute_energy_nonbond_table_linear_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_linear
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ, NT
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions
  !! @param[in]    pairlist    : interaction list in each domain
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    cpu_calc    : flag for cpu calculation or not
  !! @param[inout] coord_pbc   : coordinates for each cell
  !! @param[inout] force       : forces for each cell
  !! @param[inout] force_pbc   : forces for each cell
  !! @param[inout] virial_cell : virial term of target systems
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] ene_virial  : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_linear( &
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
      call compute_force_nonbond_tbl_lnr_generic( &
                                     domain, enefunc, pairlist, &
                                     coord_pbc, force_pbc, virial_cell)

    case (NBK_KGeneric)
      call compute_force_nonbond_tbl_lnr_k_generic( &
                                     domain, enefunc, pairlist, &
                                     coord_pbc, force_pbc, virial_cell)

    case (NBK_Fugaku)
      call compute_force_nonbond_tbl_lnr_fugaku( &
                                     domain, enefunc, pairlist, npt, &
                                     atmcls_pbc,                     &
                                     coord_pbc, force_pbc, virial)


    case (NBK_IntelKnl)
      call compute_force_nonbond_tbl_lnr_intel_knl( &
                                     domain, enefunc, pairlist, &
                                     coord_pbc, force_pbc, virial_cell)


    case (NBK_GPU)
      call compute_force_nonbond_tbl_lnr_gpu( &
                                     domain, enefunc, pairlist, &
                                     npt, cpu_calc, coord_pbc, force, &
                                     force_pbc, virial, ene_virial)

    end select

    return

  end subroutine compute_force_nonbond_table_linear

end module sp_energy_table_linear_nowater_mod
