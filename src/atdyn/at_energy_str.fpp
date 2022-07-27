!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_str_mod
!> @brief   structure of energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Takashi Imai (TI), 
!!          Jaewoon Jung (JJ), Chigusa Kobayashi(CK), Kiyoshi Yagi (KY)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_str_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_energy

    real(wp)              :: total

    ! charmm            
    real(wp)              :: bond
    real(wp)              :: angle
    real(wp)              :: urey_bradley
    real(wp)              :: dihedral
    real(wp)              :: improper
    real(wp)              :: cmap
    real(wp)              :: electrostatic
    real(wp)              :: van_der_waals

    ! go model          
    real(wp)              :: contact
    real(wp)              :: noncontact
    real(wp)              :: morph
    real(wp)              :: membrane

    ! ~CG~ 3SPN.2C DNA
    real(wp)              :: base_stacking
    real(wp)              :: base_pairing
    real(wp)              :: cg_DNA_exv
    ! ~CG~ IDR models
    real(wp)              :: cg_IDR_HPS
    real(wp)              :: cg_IDR_KH
    ! ~CG~ protein-protein KH model
    real(wp)              :: cg_KH_inter_pro
    ! ~CG~ protein-DNA sequence-specific
    real(wp)              :: PWMcos
    ! ~CG~ protein-DNA sequence-nonspecific
    real(wp)              :: PWMcosns
    ! ~CG~ general excluded volume
    real(wp)              :: cg_exv

    ! restraint
    real(wp), allocatable :: restraint(:)
    real(wp), allocatable :: restraint_cv(:)

    ! drms for morph
    real(wp)              :: drms(1:2)

    ! dipersion correction
    real(wp)              :: disp_corr_energy
    real(wp)              :: disp_corr_virial
    ! multi-basin 
    real(wp), allocatable :: basin_ratio(:)
    real(wp), allocatable :: basin_energy(:)

    ! solvation free energy
    real(wp)              :: solvation

    ! QM energy in atomic unit
    real(wp)              :: qm_ene       = 0.0_wp

    ! sphere boundary condition
    real(wp)              :: spot

    ! GaMD
    real(wp)         :: total_gamd
    real(wp)         :: dihedral_gamd

  end type s_energy

  ! parameters
  integer,      public, parameter :: ElectrostaticCutoff = 1
  integer,      public, parameter :: ElectrostaticPME    = 2

  integer,      public, parameter :: ImplicitSolventNONE = 1
  integer,      public, parameter :: ImplicitSolventEEF1 = 2
  integer,      public, parameter :: ImplicitSolventIMM1 = 3
  integer,      public, parameter :: ImplicitSolventIMIC = 4
  integer,      public, parameter :: ImplicitSolventGBSA = 5

  character(*), public, parameter :: ElectrostaticTypes(2) = (/'CUTOFF', &
                                                               'PME   '/)

  character(*), public, parameter :: ImplicitSolventTypes(5) = (/'NONE  ', &
                                                                 'EEF1  ', &
                                                                 'IMM1  ', &
                                                                 'IMIC  ', &
                                                                 'GBSA  '/)

  ! parameters for allocatable variables
  integer,      public, parameter :: EneRestraints = 1
  integer,      public, parameter :: EneMultiBasin = 2

  ! subroutines
  public  :: init_energy
  public  :: alloc_energy
  public  :: dealloc_energy
  public  :: dealloc_energy_all

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


    energy%total              = 0.0_wp
    energy%bond               = 0.0_wp
    energy%angle              = 0.0_wp
    energy%urey_bradley       = 0.0_wp
    energy%dihedral           = 0.0_wp
    energy%improper           = 0.0_wp
    energy%cmap               = 0.0_wp
    energy%electrostatic      = 0.0_wp
    energy%van_der_waals      = 0.0_wp
    energy%contact            = 0.0_wp
    energy%noncontact         = 0.0_wp
    energy%morph              = 0.0_wp
    energy%membrane           = 0.0_wp
    energy%base_stacking      = 0.0_wp
    energy%base_pairing       = 0.0_wp
    energy%cg_DNA_exv         = 0.0_wp
    energy%cg_IDR_HPS         = 0.0_wp
    energy%cg_IDR_KH          = 0.0_wp
    energy%cg_KH_inter_pro    = 0.0_wp
    energy%PWMcos             = 0.0_wp
    energy%PWMcosns           = 0.0_wp
    energy%cg_exv             = 0.0_wp

    if (allocated(energy%restraint) .and. size(energy%restraint) > 0) then
      energy%restraint(1:size(energy%restraint)) = 0.0_wp
      energy%restraint_cv(1:size(energy%restraint)) = 0.0_wp
    endif

    energy%disp_corr_energy   = 0.0_wp
    energy%disp_corr_virial   = 0.0_wp

    energy%solvation          = 0.0_wp
    energy%spot               = 0.0_wp

    ! GaMD
    energy%total_gamd         = 0.0_wp
    energy%dihedral_gamd      = 0.0_wp

    return

  end subroutine init_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_energy
  !> @brief        allocate energy information
  !! @authors      CK
  !! @param[inout] energy : energy information
  !! @param[in]    variable : allocatable variables
  !! @param[in]    var_size : size of variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_energy(energy, variable, var_size)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(EneRestraints)

       if (var_size <= 0) return
       if (allocated(energy%restraint)) then
          if (size(energy%restraint(:)) == var_size) return
          deallocate(energy%restraint,           &
                     energy%restraint_cv,        &
                     stat = dealloc_stat)
       end if

       allocate(energy%restraint(var_size),        &
                energy%restraint_cv(var_size),     &
                stat = alloc_stat)

       energy%restraint    (1:var_size) = 0.0_wp
       energy%restraint_cv (1:var_size) = 0.0_wp

    case(EneMultiBasin)

       if (var_size <= 0) return
       if (allocated(energy%basin_ratio)) then
          if (size(energy%basin_ratio(:)) == var_size) return
          deallocate(energy%basin_ratio,        &
                     energy%basin_energy,       &
                     stat = dealloc_stat)
       end if

       allocate(energy%basin_ratio(var_size),        &
                energy%basin_energy(var_size),       &
                stat = alloc_stat)

       energy%basin_ratio    (1:var_size) = 0.0_wp
       energy%basin_energy   (1:var_size) = 0.0_wp


    case default

      call error_msg('Alloc_Energy> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine alloc_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_energy
  !> @brief        deallocate energy information
  !! @authors      CK
  !! @param[inout] energy : energy information
  !! @param[in]    variable : allocatable variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_energy(energy, variable)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat


    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(EneRestraints)

       if (allocated(energy%restraint)) then
          deallocate(energy%restraint,        &
                     energy%restraint_cv,     &
                     stat = dealloc_stat)
       end if

    case(EneMultiBasin)

       if (allocated(energy%basin_ratio)) then
          deallocate(energy%basin_ratio,        &
                     energy%basin_energy,       &
                     stat = dealloc_stat)
       end if

    case default

      call error_msg('Dealloc_Energy> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc    

    return

  end subroutine dealloc_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_energy_all
  !> @brief        deallocate all energy information
  !! @authors      CK
  !! @param[inout] energy : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_energy_all(energy)

    ! formal arguments
    type(s_energy),         intent(inout) :: energy


    call dealloc_energy(energy, EneRestraints)
    call dealloc_energy(energy, EneMultiBasin)

    return

  end subroutine dealloc_energy_all

end module at_energy_str_mod
