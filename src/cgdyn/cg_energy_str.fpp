!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_str_mod
!> @brief   Defines the structure of energy.
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Takaharu Mori (TM)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_energy
    real(dp)         :: total
    real(dp)         :: bond
    real(dp)         :: angle
    real(dp)         :: urey_bradley
    real(dp)         :: dihedral
    real(dp)         :: improper
    real(dp)         :: cmap

    ! ~CG~ 3SPN.2C DNA
    real(dp)         :: base_stacking
    real(dp)         :: base_pairing
    real(dp)         :: cg_DNA_exv

    real(dp)         :: electrostatic
    real(dp)         :: van_der_waals
    real(dp)         :: contact
    real(dp)         :: noncontact
    ! ~CG~ protein-DNA sequence-specific
    real(dp)         :: PWMcos
    real(dp)         :: PWMcosns
    ! ~CG~ general excluded volume
    real(dp)         :: cg_exv
    ! ~CG~ IDR models
    real(dp)         :: cg_IDR_HPS
    real(dp)         :: cg_IDR_KH
    real(dp)         :: cg_IDR_MPIPI
    real(dp)         :: cg_KH_inter_pro

    ! restraint
    real(dp)         :: restraint_position
    real(dp)         :: restraint_rmsd
    real(dp)         :: restraint_emfit
    real(dp)         :: restraint_distance
    real(dp)         :: rmsd
    real(dp)         :: emcorr
    ! dispersion correction
    real(dp)         :: disp_corr_energy
    real(dp)         :: disp_corr_virial
  end type s_energy

  ! parameters
  integer,      public, parameter :: ElectrostaticCutoff = 1
  integer,      public, parameter :: ElectrostaticPME    = 2
  integer,      public, parameter :: ElectrostaticRF     = 3

  character(*), public, parameter :: ElectrostaticTypes(3) = (/'CUTOFF', &
                                                               'PME   ', &
                                                               'RF    '/)

  integer,      public, parameter :: StructureCheckNone    = 1
  integer,      public, parameter :: StructureCheckFirst   = 2
  integer,      public, parameter :: StructureCheckDomain  = 3

  character(*), public, parameter :: StructureCheckTypes(3) = (/'NONE  ', &
                                                                'FIRST ', &
                                                                'DOMAIN'/)

  ! subroutines
  public  :: init_energy

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_energy
  !> @brief        Initializes potential energy.
  !! @authors      YS, TM
  !! @param[out]   energy : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_energy(energy)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy


    energy%total              = 0.0_dp
    energy%bond               = 0.0_dp
    energy%angle              = 0.0_dp
    energy%urey_bradley       = 0.0_dp
    energy%dihedral           = 0.0_dp
    energy%improper           = 0.0_dp
    energy%cmap               = 0.0_dp
    energy%base_stacking      = 0.0_dp
    energy%base_pairing       = 0.0_dp
    energy%cg_DNA_exv         = 0.0_dp
    energy%electrostatic      = 0.0_dp
    energy%van_der_waals      = 0.0_dp
    energy%contact            = 0.0_dp
    energy%noncontact         = 0.0_dp
    energy%PWMcos             = 0.0_dp
    energy%PWMcosns           = 0.0_dp
    energy%cg_exv             = 0.0_dp
    energy%cg_KH_inter_pro    = 0.0_dp
    energy%cg_IDR_KH          = 0.0_dp
    energy%cg_IDR_HPS         = 0.0_dp
    energy%cg_IDR_MPIPI       = 0.0_dp
    energy%restraint_position = 0.0_dp
    energy%restraint_distance = 0.0_dp
    energy%restraint_rmsd     = 0.0_dp
    energy%restraint_emfit    = 0.0_dp
    energy%disp_corr_energy   = 0.0_dp
    energy%disp_corr_virial   = 0.0_dp
    energy%rmsd               = 0.0_dp
    energy%emcorr             = 0.0_dp

    return

  end subroutine init_energy

end module cg_energy_str_mod
