!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbond_ljpme_mod
!> @brief   calculate nonbonded energy with ljpme
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbond_ljpme_mod

  use sp_energy_nonbond_gpu_mod
  use sp_energy_nonbond_fugaku_mod
  use sp_energy_nonbond_generic_mod
  use sp_energy_nonbond_intel_mod
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
  public  :: compute_energy_nonbond_table_ljpme
  public  :: compute_force_nonbond_table_ljpme

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_ljpme 
  !> @brief        calculate nonbonded energy with ljpme
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

  subroutine compute_energy_nonbond_table_ljpme(domain, enefunc, pairlist,  &
                                          npt, atmcls_pbc, coord_pbc,       &
                                          force_pbc, virial_cell, virial,   &
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
      call compute_energy_nonbond_tbl_ljpme_generic( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)

    case (NBK_Fugaku)
      call compute_energy_nonbond_tbl_ljpme_fugaku( &
                              domain, enefunc, pairlist, atmcls_pbc, &
                              coord_pbc, force_pbc, virial, eelec, evdw)

    case (NBK_Intel)
      call compute_energy_nonbond_tbl_ljpme_intel( &
                              domain, enefunc, pairlist, atmcls_pbc, &
                              coord_pbc, force_pbc, virial, eelec, evdw)

    case (NBK_GPU)
      call compute_energy_nonbond_tbl_ljpme_gpu( &
                              domain, enefunc, pairlist, npt,        &
                              coord_pbc, force_pbc, virial,          &
                              eelec, evdw, ene_virial)

    end select

    return

  end subroutine compute_energy_nonbond_table_ljpme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_ljpme
  !> @brief        calculate nonbonded force with lj pme
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

  subroutine compute_force_nonbond_table_ljpme(domain, enefunc, pairlist, &
                                          npt, cpu_calc, atmcls_pbc,      &
                                          coord_pbc, force, force_pbc,    &
                                          virial_cell, virial, ene_virial)

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
      call compute_force_nonbond_tbl_ljpme_generic( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell)

    case (NBK_Fugaku)
      call compute_force_nonbond_tbl_ljpme_fugaku( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc, virial)


    case (NBK_Intel)
      call compute_force_nonbond_tbl_ljpme_intel( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc, virial)


    case (NBK_GPU)
      call compute_force_nonbond_tbl_ljpme_gpu( &
                              domain, enefunc, pairlist, npt,        &
                              cpu_calc, coord_pbc, force,            &
                              force_pbc, virial, ene_virial)

    end select

    return

  end subroutine compute_force_nonbond_table_ljpme

end module sp_energy_nonbond_ljpme_mod
