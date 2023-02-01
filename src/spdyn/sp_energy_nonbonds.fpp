!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbonds_mod
!> @brief   calculate nonbond energy
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbonds_mod

  use sp_energy_pme_mod
  use sp_energy_table_cubic_mod
  use sp_energy_table_linear_bondcorr_mod
  use sp_energy_nonbond_pme_table_mod
  use sp_energy_table_linear_mod
  use sp_energy_nonbond_pme_notable_mod
  use sp_energy_nonbond_ljpme_mod
  use sp_energy_nonbond_generic_mod
  use sp_energy_nonbond_vacuum_mod
  use sp_pairlist_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_nonbond_cutoff
  public  :: compute_energy_nonbond_pme
  public  :: compute_energy_nonbond_pme_short
  public  :: compute_energy_nonbond_pme_long
  ! FEP
  public  :: compute_energy_nonbond_cutoff_fep
  public  :: compute_energy_nonbond_pme_fep
  public  :: compute_energy_nonbond_pme_short_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_cutoff
  !> @brief        compute nonbond energy with cutoff 
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pair-list information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                                           nonb_ene, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)


    call timer(TimerNonBond, TimerOn)

    if (enefunc%vacuum) then
      ! For vacuum without lookup table

      call compute_energy_nonbond14_vacuum( &
        domain, enefunc, &
        force, virial, eelec, evdw)

      if (nonb_ene) then

        call compute_energy_nonbond_vacuum( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)

      else

        call compute_force_nonbond_vacuum( &
          domain, enefunc, pairlist, &
          force, virial)

      end if

    else

      call compute_energy_nonbond14_table( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)
   
      if (nonb_ene) then

        call compute_energy_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial, eelec, evdw)
      
      else

        call compute_force_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial)
      
      end if

    end if

    call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_cutoff

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme
  !> @brief        Calculate nonbond energy by PME
  !! @authors      JJ, CK
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for contact checker
  !! @param[inout] coord_pbc    : pbc oriented coordinates
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : forces for each cell
  !! @param[inout] virial_cell  : virial term of target system in cell
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] eelec        : electrostatic energy of target systems
  !! @param[inout] evdw         : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme(domain, enefunc, pairlist, boundary, &
                                        npt, nonb_ene, nonb_limiter,         &
                                        atmcls_pbc,                          &
                                        coord_pbc, force_long, force,        &
                                        force_pbc, virial_cell, virial,      &
                                        eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    !local variables
    real(dp)                 :: ene_virial(5)


    call timer(TimerNonBond, TimerOn)

    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

      if (nonb_limiter) then

        if (enefunc%vdw == VDWCutoff) then
          call compute_energy_nonbond_pme_notbl_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else if (enefunc%vdw == VDWPME) then
          call compute_energy_nonbond_tbl_ljpme_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else
          call compute_energy_nonbond_table_linear_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        end if

      else

        if (nonb_ene) then

          if (enefunc%vdw == VDWCutoff) then
            call compute_energy_nonbond_pme_notbl( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc,      &
                              virial_cell, virial, eelec, evdw,      &
                              ene_virial)
          else if (enefunc%vdw == VDWPME) then
            call compute_energy_nonbond_table_ljpme( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc,      &
                              virial_cell, virial, eelec, evdw,      &
                              ene_virial)

          else 
            call compute_energy_nonbond_table_linear( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc,      &
                              virial_cell, virial, eelec, evdw,      &
                              ene_virial)
          end if

        else

          if (enefunc%vdw == VDWCutoff) then
            call compute_force_nonbond_pme_notbl( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
          else if (enefunc%vdw == VDWPME) then
            call compute_force_nonbond_table_ljpme( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
          else
            call compute_force_nonbond_table_linear( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
          end if

        end if
      end if

      if (domain%nonbond_kernel /= NBK_GPU) &
        call timer(TimerPmeReal, TimerOff)

    end if


    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then
 
      call timer(TimerPmeRecip, TimerOn)

      if (npt) then
        if (enefunc%vdw == VDWPME) then
          call pme_pre_lj(domain, boundary)
        else
          call pme_pre(domain, boundary)
        end if
      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_barrier(mpi_comm_city, ierror)
#endif

      if (enefunc%vdw == VDWPME) then
        call pme_recip_lj(domain, enefunc, force_long, virial, eelec, evdw)
      else
        call pme_recip(domain, force_long, virial, eelec)
      end if

      call timer(TimerPmeRecip, TimerOff)

    end if

    ! Add self energy
    !
    eelec(1) = eelec(1) + u_self

    if (domain%nonbond_kernel /= NBK_GPU) &
      call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_short
  !> @brief        Calculate nonbond energy (real part) by PME
  !! @authors      JJ
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[inout] coord_pbc    : pbc oriented coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : forces for each cell
  !! @param[inout] virial_cell  : virial term of target system in cell
  !! @param[inout] virial       : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_short(domain, enefunc, pairlist,       &
                                              npt,                             &
                                              atmcls_pbc,                      &
                                              coord_pbc, force, force_pbc,     &
                                              virial_cell, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: npt
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(dp)                 :: ene_virial(5), eelec(nthread), evdw(nthread)


    call timer(TimerNonBond, TimerOn)
 
    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

      if (enefunc%nonb_limiter) then

        if (enefunc%vdw == VDWCutoff) then
          call compute_energy_nonbond_pme_notbl_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else if (enefunc%vdw == VDWPME) then
          call compute_energy_nonbond_tbl_ljpme_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else
          call compute_energy_nonbond_table_linear_check( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        end if

      else

        if (enefunc%vdw == VDWCutoff) then
          call compute_force_nonbond_pme_notbl( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
        else if (enefunc%vdw == VDWPME) then
          call compute_force_nonbond_table_ljpme( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
        else
          call compute_force_nonbond_table_linear( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
        end if

        if (domain%nonbond_kernel /= NBK_GPU) &
          call timer(TimerPmeReal, TimerOff)

      end if

    end if

    if (domain%nonbond_kernel /= NBK_GPU) &
      call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_pme_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_long
  !> @brief        Calculate nonbond energy (reciprocal part) by PME
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary information
  !! @param[in]    npt      : flag for NPT or not
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_long(domain, enefunc, boundary, &
                                             npt, force, virial, eelec, &
                                             evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw (nthread)


    call timer(TimerNonBond, TimerOn)

    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then

      call timer(TimerPmeRecip, TimerOn)

      if (npt) then
        if (enefunc%vdw == VDWPME) then
          call pme_pre_lj(domain, boundary)
        else
          call pme_pre(domain, boundary)
        end if
      end if

      if (enefunc%vdw == VDWPME) then
        call pme_recip_lj(domain, enefunc, force, virial, eelec, evdw)
      else
        call pme_recip(domain, force, virial, eelec)
      end if

      call timer(TimerPmeRecip, TimerOff)

    end if

    ! Add self energy
    !
    eelec(1) = eelec(1) + u_self

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_pme_long

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_cutoff_fep
  !> @brief        compute nonbond energy with cutoff for FEP 
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pair-list information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_cutoff_fep(domain, enefunc, pairlist, &
                                           nonb_ene, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)


    call timer(TimerNonBond, TimerOn)

    if (enefunc%vacuum) then
      ! For vacuum without lookup table

      call compute_energy_nonbond14_vacuum_fep( &
        domain, enefunc, &
        force, virial, eelec, evdw)

      if (nonb_ene) then

        call compute_energy_nonbond_vacuum( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)

        call compute_energy_nonbond_vacuum_fep( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)

      else

        call compute_force_nonbond_vacuum( &
          domain, enefunc, pairlist, &
          force, virial)

        call compute_force_nonbond_vacuum_fep( &
          domain, enefunc, pairlist, &
          force, virial)

      end if

    else

      call compute_energy_nonbond14_table_fep( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)
   
      if (nonb_ene) then

        call compute_energy_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial, eelec, evdw)
   
        call compute_energy_nonbond_table_fep( &
                                    domain, enefunc, pairlist, &
                                    force, virial, eelec, evdw)
      
      else

        call compute_force_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial)
      
        call compute_force_nonbond_table_fep( &
                                  domain, enefunc, pairlist, &
                                  force, virial)
   
      end if

    end if

    call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_cutoff_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_fep
  !> @brief        Calculate nonbond energy by PME for FEP
  !! @authors      HO
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for contact checker
  !! @param[inout] coord_pbc    : pbc oriented coordinates
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : forces for each cell
  !! @param[inout] virial_cell  : virial term of target system in cell
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] eelec        : electrostatic energy of target systems
  !! @param[inout] evdw         : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_fep(domain, enefunc, pairlist, &
                                        boundary, &
                                        npt, nonb_ene, nonb_limiter,         &
                                        atmcls_pbc,                          &
                                        coord_pbc, force_long, force,        &
                                        force_pbc, virial_cell, virial,      &
                                        eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(inout)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    !local variables
    real(dp)                 :: ene_virial(5)


    call timer(TimerNonBond, TimerOn)

    !
    ! VDWPME is not allowed in FEP
    !

    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

      if (nonb_limiter) then

        if (enefunc%vdw == VDWCutoff) then
          call compute_energy_nonbond_pme_notbl_check_fep( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else if (enefunc%vdw /= VDWPME) then
          call compute_energy_nonbond_table_linear_check_fep( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        end if

      else

        if (nonb_ene) then

          if (enefunc%vdw == VDWCutoff) then
            call compute_energy_nonbond_pme_notbl_fep( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc,      &
                              virial_cell, virial, eelec, evdw,      &
                              ene_virial)
          else if (enefunc%vdw /= VDWPME) then
            call compute_energy_nonbond_table_linear_fep( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc, force_pbc,      &
                              virial_cell, virial, eelec, evdw,      &
                              ene_virial)
          end if

        else

          if (enefunc%vdw == VDWCutoff) then
            call compute_force_nonbond_pme_notbl_fep( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
          else if (enefunc%vdw /= VDWPME) then
            call compute_force_nonbond_table_linear_fep( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
          end if

        end if
      end if

      if (domain%nonbond_kernel /= NBK_GPU) &
        call timer(TimerPmeReal, TimerOff)

    end if


    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then
 
      call timer(TimerPmeRecip, TimerOn)

      if (npt) then
        call pme_pre(domain, boundary)
      end if

      call mpi_barrier(mpi_comm_city, ierror)

      call pme_recip(domain, force_long, virial, eelec)

      call timer(TimerPmeRecip, TimerOff)

    end if

    ! Add self energy
    !
    eelec(1) = eelec(1) + u_self

    if (domain%nonbond_kernel /= NBK_GPU) &
      call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_pme_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_short_fep
  !> @brief        Calculate nonbond energy (real part) by PME for FEP
  !! @authors      HO
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[inout] coord_pbc    : pbc oriented coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : forces for each cell
  !! @param[inout] virial_cell  : virial term of target system in cell
  !! @param[inout] virial       : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_short_fep(domain, enefunc, pairlist,   &
                                              npt,                             &
                                              atmcls_pbc,                      &
                                              coord_pbc, force, force_pbc,     &
                                              virial_cell, virial)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: npt
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(:,:,:)

    ! local variables
    real(dp)                 :: ene_virial(5), eelec(nthread), evdw(nthread)


    call timer(TimerNonBond, TimerOn)
 
    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

      if (enefunc%nonb_limiter) then

        if (enefunc%vdw == VDWCutoff) then
          call compute_energy_nonbond_pme_notbl_check_fep( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        else if (enefunc%vdw /= VDWPME) then
          call compute_energy_nonbond_table_linear_check_fep( &
                              domain, enefunc, pairlist, coord_pbc,  &
                              force_pbc, virial_cell, eelec, evdw)
        end if

      else

        if (enefunc%vdw == VDWCutoff) then
          call compute_force_nonbond_pme_notbl_fep( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
        else if (enefunc%vdw /= VDWPME) then
          call compute_force_nonbond_table_linear_fep( &
                              domain, enefunc, pairlist, npt,        &
                              .false., atmcls_pbc, coord_pbc,        &
                              force, force_pbc, virial_cell, virial, &
                              ene_virial)
        end if

        if (domain%nonbond_kernel /= NBK_GPU) &
          call timer(TimerPmeReal, TimerOff)

      end if

    end if

    if (domain%nonbond_kernel /= NBK_GPU) &
      call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_pme_short_fep

end module sp_energy_nonbonds_mod
