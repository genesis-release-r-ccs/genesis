!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_str_mod
!> @brief   structure of energy
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Takaharu Mori (TM)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_energy
    real(dp)         :: total
    ! charmm
    real(dp)         :: bond
    real(dp)         :: enm 
    real(dp)         :: angle
    real(dp)         :: urey_bradley
    real(dp)         :: dihedral
    real(dp)         :: improper
    real(dp)         :: cmap
    real(dp)         :: electrostatic
    real(dp)         :: van_der_waals
    real(dp)         :: electric_field
    real(dp)         :: contact
    real(dp)         :: noncontact
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
    ! gamd
    real(dp)         :: total_gamd
    real(dp)         :: dihedral_gamd
    ! FEP
    real(dp)         :: deltU_fep(3)
  end type s_energy

  ! parameters
  integer,      public, parameter :: ElectrostaticCutoff = 1
  integer,      public, parameter :: ElectrostaticPME    = 2

  character(*), public, parameter :: ElectrostaticTypes(2) = (/'CUTOFF', &
                                                               'PME   '/)

  integer,      public, parameter :: VDWCutoff           = 1
  integer,      public, parameter :: VDWSwitch           = 2
  integer,      public, parameter :: VDWFSW              = 3
  integer,      public, parameter :: VDWShift            = 4
  integer,      public, parameter :: VDWPME              = 5

  character(*), public, parameter :: VDW_Types(5) = (/'CUTOFF', &
                                                      'SWITCH', &
                                                      'FSW   ', &
                                                      'SHIFT ', &
                                                      'PME   '/)
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
  !> @brief        initialize potential energy
  !! @authors      YS, TM
  !! @param[out]   energy : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_energy(energy)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy


    energy%total              = 0.0_dp
    energy%bond               = 0.0_dp
    energy%enm                = 0.0_dp
    energy%angle              = 0.0_dp
    energy%urey_bradley       = 0.0_dp
    energy%dihedral           = 0.0_dp
    energy%improper           = 0.0_dp
    energy%cmap               = 0.0_dp
    energy%electrostatic      = 0.0_dp
    energy%van_der_waals      = 0.0_dp
    energy%electric_field     = 0.0_dp
    energy%contact            = 0.0_dp
    energy%noncontact         = 0.0_dp
    energy%restraint_position = 0.0_dp
    energy%restraint_distance = 0.0_dp
    energy%restraint_rmsd     = 0.0_dp
    energy%restraint_emfit    = 0.0_dp
    energy%disp_corr_energy   = 0.0_dp
    energy%disp_corr_virial   = 0.0_dp
    energy%rmsd               = 0.0_dp
    energy%emcorr             = 0.0_dp
    energy%total_gamd         = 0.0_dp
    energy%dihedral_gamd      = 0.0_dp
    ! FEP
    energy%deltU_fep          = 0.0_dp

    return

  end subroutine init_energy

end module sp_energy_str_mod
