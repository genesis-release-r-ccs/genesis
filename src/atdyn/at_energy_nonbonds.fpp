!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_nonbonds_mod
!> @brief   calculate nonbonded energy
!! @authors Takaharu Mori (TM), Yuji Sugita (YS), Jaewoon Jung (JJ), 
!!          Takashi Imai (TI),  Chigusa Kobayashi (CK), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_nonbonds_mod

  use at_energy_pme_mod
  use at_energy_table_cubic_mod
  use at_energy_table_linear_bondcorr_mod
  use at_energy_table_linear_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_nonbond_nobc
  public  :: compute_energy_nonbond_cutoff 
  public  :: compute_energy_nonbond_pme

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_nobc
  !> @brief        Calculate nonbond energy in no boundary condition
  !! @authors      JJ, NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene,&
                                         coord, force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    call timer(TimerNonBond, TimerOn)

    call timer(TimerPmeReal, TimerOn)

    call compute_energy_nonbond14_table( &
                                  enefunc, molecule, &
                                  coord, force, virial, eelec, evdw)

    ! ==> Type 3 & 9
    if (nonb_ene) then
    !CK: use routines for charmm

      call compute_energy_nonbond_nobc_table( &
                                  enefunc, molecule, pairlist, &
                                  coord, force, virial, eelec, evdw)
    
    else

      call compute_force_nonbond_nobc_table( &
                                  enefunc, molecule, pairlist, &
                                  coord, force, virial)
    
    end if

    call timer(TimerPmeReal, TimerOff)

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_cutoff
  !> @brief        Calculate nonbond energy by cutoff
  !! @authors      JJ, NT, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy 
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_cutoff(enefunc, molecule, pairlist, &
                                           boundary, nonb_ene, coord,   &
                                           force, virial, eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    call timer(TimerNonBond, TimerOn)

    call timer(TimerPmeReal, TimerOn)

    call compute_energy_nonbond14_table( &
                                  enefunc, molecule, &
                                  coord, force, virial, eelec, evdw)

    ! ==> Type 2 and 8
    if (nonb_ene) then

      call compute_energy_nonbond_table( &
                                  enefunc, molecule, pairlist, &
                                  boundary, coord, force, virial, &
                                  eelec, evdw)

    else

      call compute_force_nonbond_table( &
                                  enefunc, molecule, pairlist, &
                                  boundary, coord, force, virial)

    end if


    call timer(TimerPmeReal, TimerOff)

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_cutoff

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme
  !> @brief        Calculate nonbond energy by PME
  !! @authors      TM, NT, CK
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    molecule     : molecule information
  !! @param[in]    pairlist     : pairlist information
  !! @param[in]    boundary     : boundary conditions information
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy 
  !! @param[in]    nonb_limiter : flag for contact checker
  !! @param[in]    coord        : coordinates of target systems
  !! @param[inout] force        : forces of target systems
  !! @param[inout] virial       : virial of target systems
  !! @param[inout] eelec        : electrostatic energy of target systems
  !! @param[inout] evdw         : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme(enefunc, molecule, pairlist, boundary, &
                                        nonb_ene, nonb_limiter,                &
                                        coord, force, virial,                  &
                                        eelec, evdw)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw


    call timer(TimerNonBond, TimerOn)


    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

      call compute_energy_nonbond14_table_linear( &
                                  enefunc, molecule, &
                                  coord, force, virial, eelec, evdw)

      ! ==> Type 5 & 11
      if (nonb_limiter) then
          call compute_energy_nonbond_table_linear_check( &
                                    enefunc, molecule, pairlist, &
                                    boundary, coord, force, virial, &
                                    eelec, evdw)

      else
        if (nonb_ene) then
      
          call compute_energy_nonbond_table_linear( &
                                    enefunc, molecule, pairlist, &
                                    boundary, coord, force, virial, &
                                    eelec, evdw)
        
        else
       
          call compute_force_nonbond_table_linear( &
                                    enefunc, molecule, pairlist, &
                                    boundary, coord, force, virial)
        
        end if
      end if

      call pme_bond_corr_linear(enefunc, molecule, coord, &
                                  force, virial, eelec)

      call timer(TimerPmeReal, TimerOff)

    end if

    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then

      call timer(TimerPmeRecip, TimerOn)

      call pme_pre(boundary, enefunc)
      call pme_recip(enefunc, boundary, molecule, &
                     coord, force, virial, eelec)

      call timer(TimerPmeRecip, TimerOff)

    end if

    ! Add self energy
    eelec = eelec + enefunc%pme%u_self

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_pme

end module at_energy_nonbonds_mod
