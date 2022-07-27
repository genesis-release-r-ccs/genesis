!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_gamd_mod
!> @brief   GaMD method
!           Refs: Miao Y. et al., J. Chem. Theory Comput. 11, 3584 (2015)
!                 Pang Y. et al., J. Chem. Theory Comput. 13, 9 (2017)
!! @authors Hiraku Oshima (HO)
! 
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_gamd_mod

  use at_dynamics_str_mod
  use at_enefunc_str_mod
  use at_enefunc_gamd_mod
  use at_remd_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters for GAMD type
  integer, public, parameter      :: GamdBoostPotential = 1
  integer, public, parameter      :: GamdBoostDihedral  = 2
  integer, public, parameter      :: GamdBoostDual      = 3
  character(*), public, parameter :: GamdBoostTypes(3)  = (/'POTENTIAL      ',&
                                                            'DIHEDRAL       ',&
                                                            'DUAL           '/)

  integer, public, parameter      :: GamdThreshLower = 1
  integer, public, parameter      :: GamdThreshUpper = 2
  character(*), public, parameter :: GamdThreshTypes(2)  = (/'LOWER           ', &
                                                             'UPPER           '/)

  ! structures
  type, public :: s_gamd_info
    logical  :: gamd          = .false.
    logical  :: boost         = .true.
    integer  :: boost_type    = GamdBoostDual
    integer  :: thresh_type   = GamdThreshLower
    real(wp) :: sigma0_pot    = 6.0_wp
    real(wp) :: sigma0_dih    = 6.0_wp
    integer  :: update_period = 0
    real(wp) :: pot_max       = -99999999.0_wp
    real(wp) :: pot_min       =  99999999.0_wp
    real(wp) :: pot_ave       =  0.0_wp
    real(wp) :: pot_dev       =  0.0_wp
    real(wp) :: dih_max       = -99999999.0_wp
    real(wp) :: dih_min       =  99999999.0_wp
    real(wp) :: dih_ave       =  0.0_wp
    real(wp) :: dih_dev       =  0.0_wp
  end type s_gamd_info

  ! public subroutines
  public  :: show_ctrl_gamd
  public  :: read_ctrl_gamd
  public  :: setup_gamd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_gamd
  !> @brief        show GAMD section usage
  !! @authors      HO
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "gamd"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine show_ctrl_gamd(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('gamd')

        write(MsgOut,'(A)') '[GAMD]'
        write(MsgOut,'(A)') 'gamd          = no    # Gaussian accelarated MD (GaMD) [yes, no]'
        write(MsgOut,'(A)') 'boost_type    = DUAL  # boost type of GaMD [DIHEDRAL, POTENTIAL, DUAL]'
        write(MsgOut,'(A)') 'thresh_type   = LOWER # threshold of GaMD [LOWER, UPPER]'
        write(MsgOut,'(A)') 'sigma0_pot    = 6.0   #'
        write(MsgOut,'(A)') 'sigma0_dih    = 6.0   #'
        write(MsgOut,'(A)') 'update_period = 0     # update period for GaMD parameter'
        write(MsgOut,'(A)') '# pot_max       = 0.0 # maximum of potential energy (kcal/mol)'
        write(MsgOut,'(A)') '# pot_min       = 0.0 # minimum of potential energy (kcal/mol)'
        write(MsgOut,'(A)') '# pot_ave       = 0.0 # average of potential energy (kcal/mol)'
        write(MsgOut,'(A)') '# pot_dev       = 0.0 # standard deviation of potential energy (kcal/mol)'
        write(MsgOut,'(A)') '# dih_max       = 0.0 # maximum of dihedral energy (kcal/mol)'
        write(MsgOut,'(A)') '# dih_min       = 0.0 # minimum of dihedral energy (kcal/mol)'
        write(MsgOut,'(A)') '# dih_ave       = 0.0 # average of dihedral energy (kcal/mol)'
        write(MsgOut,'(A)') '# dih_dev       = 0.0 # standard deviation of dihedral energy (kcal/mol)'
        write(MsgOut,'(A)') '# boost         = yes # Boosting [yes, no]'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('gamd')

        write(MsgOut,'(A)') '[GAMD]'
        write(MsgOut,'(A)') 'gamd          = no    # Gaussian accelarated MD (GaMD) [yes, no]'
        write(MsgOut,'(A)') 'boost_type    = DUAL  # boost type of GaMD [DIHEDRAL, POTENTIAL, DUAL]'
        write(MsgOut,'(A)') 'thresh_type   = LOWER # threshold of GaMD [LOWER, UPPER]'
        write(MsgOut,'(A)') 'sigma0_pot    = 6.0   #'
        write(MsgOut,'(A)') 'sigma0_dih    = 6.0   #'
        write(MsgOut,'(A)') 'update_period = 0     # update period for GaMD parameter'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_gamd
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_gamd
  !> @brief        read GAMD section in the control file
  !! @authors      HO
  !! @param[in]    handle    : unit number of control file
  !! @param[out]   gamd_info : GAMD section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_gamd(handle, gamd_info)

    ! parameters
    character(*),            parameter     :: Section = 'Gamd'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_gamd_info),       intent(inout) :: gamd_info

    ! local variables


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'gamd',  &
                               gamd_info%gamd)
    call read_ctrlfile_logical(handle, Section, 'boost',  &
                               gamd_info%boost)
    call read_ctrlfile_type   (handle, Section, 'boost_type',  &
                               gamd_info%boost_type, GamdBoostTypes)
    call read_ctrlfile_type   (handle, Section, 'thresh_type', &
                               gamd_info%thresh_type, GamdThreshTypes)
    call read_ctrlfile_real   (handle, Section, 'sigma0_pot',  &
                               gamd_info%sigma0_pot)
    call read_ctrlfile_real   (handle, Section, 'sigma0_dih',  &
                               gamd_info%sigma0_dih)
    call read_ctrlfile_integer(handle, Section, 'update_period',&
                               gamd_info%update_period)

    call read_ctrlfile_real(handle, Section, 'pot_max',       &
                               gamd_info%pot_max)
    call read_ctrlfile_real(handle, Section, 'pot_min',       &
                               gamd_info%pot_min)
    call read_ctrlfile_real(handle, Section, 'pot_ave',       &
                               gamd_info%pot_ave)
    call read_ctrlfile_real(handle, Section, 'pot_dev',       &
                               gamd_info%pot_dev)

    call read_ctrlfile_real(handle, Section, 'dih_max',       &
                               gamd_info%dih_max)
    call read_ctrlfile_real(handle, Section, 'dih_min',       &
                               gamd_info%dih_min)
    call read_ctrlfile_real(handle, Section, 'dih_ave',       &
                               gamd_info%dih_ave)
    call read_ctrlfile_real(handle, Section, 'dih_dev',       &
                               gamd_info%dih_dev)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      if (gamd_info%gamd) then

        write(MsgOut,'(A)') 'Read_Ctrl_Gamd> Gamd information'

        if (gamd_info%boost_type == GamdBoostPotential) then
          write(MsgOut,'(A)')   '  boost_type       = POTENTIAL'
          write(MsgOut,'(A20,F10.5)') &
            '  sigma0_pot       = ', gamd_info%sigma0_pot
          write(MsgOut,'(A20,F18.10)') &
            '  pot_max          = ', gamd_info%pot_max
          write(MsgOut,'(A20,F18.10)') &
            '  pot_min          = ', gamd_info%pot_min
          write(MsgOut,'(A20,F18.10)') &
            '  pot_ave          = ', gamd_info%pot_ave
          write(MsgOut,'(A20,F18.10)') &
            '  pot_dev          = ', gamd_info%pot_dev
        else if (gamd_info%boost_type == GamdBoostDihedral) then
          write(MsgOut,'(A)')   '  boost_type       = DIHEDRAL'
          write(MsgOut,'(A20,F10.5)') &
            '  sigma0_dih       = ', gamd_info%sigma0_dih
          write(MsgOut,'(A20,F18.10)') &
            '  dih_max          = ', gamd_info%dih_max
          write(MsgOut,'(A20,F18.10)') &
            '  dih_min          = ', gamd_info%dih_min
          write(MsgOut,'(A20,F18.10)') &
            '  dih_ave          = ', gamd_info%dih_ave
          write(MsgOut,'(A20,F18.10)') &
            '  dih_dev          = ', gamd_info%dih_dev
        else if (gamd_info%boost_type == GamdBoostDual) then
          write(MsgOut,'(A)')   '  boost_type       = DUAL'
          write(MsgOut,'(A20,F10.5)') &
            '  sigma0_pot       = ', gamd_info%sigma0_pot
          write(MsgOut,'(A20,F10.5)') &
            '  sigma0_dih       = ', gamd_info%sigma0_dih
          write(MsgOut,'(A20,E18.10)') &
            '  pot_max          = ', gamd_info%pot_max
          write(MsgOut,'(A20,E18.10)') &
            '  pot_min          = ', gamd_info%pot_min
          write(MsgOut,'(A20,E18.10)') &
            '  pot_ave          = ', gamd_info%pot_ave
          write(MsgOut,'(A20,E18.10)') &
            '  pot_dev          = ', gamd_info%pot_dev
          write(MsgOut,'(A20,E18.10)') &
            '  dih_max          = ', gamd_info%dih_max
          write(MsgOut,'(A20,E18.10)') &
            '  dih_min          = ', gamd_info%dih_min
          write(MsgOut,'(A20,E18.10)') &
            '  dih_ave          = ', gamd_info%dih_ave
          write(MsgOut,'(A20,E18.10)') &
            '  dih_dev          = ', gamd_info%dih_dev

        end if

        if (gamd_info%thresh_type == GamdThreshLower) then
          write(MsgOut,'(A)')   '  thresh_type      = LOWER'
        else if (gamd_info%thresh_type == GamdThreshUpper) then
          write(MsgOut,'(A)')   '  thresh_type      = UPPER'
        end if

        write(MsgOut,'(A20,I10)') &
          '  update_period    = ', gamd_info%update_period

        if (gamd_info%boost) then
          write(MsgOut,'(A)')   '  boost            = YES'
        else
          write(MsgOut,'(A)')   '  boost            = NO'
        end if

        write(MsgOut,'(A)') ''

      end if

    end if

    return

  end subroutine read_ctrl_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_gamd
  !> @brief        setup GAMD
  !! @authors      HO
  !! @param[in]    gamd_info : GAMD section control parameters information
  !! @param[inout] dynamics  : dynamics information
  !! @param[inout] molecule  : molecule information
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[inout] remd      : REMD information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_gamd(gamd_info, dynamics, molecule, enefunc, remd)

    ! formal arguments
    type(s_gamd_info),       intent(in)    :: gamd_info
    type(s_dynamics),target, intent(inout) :: dynamics
    type(s_molecule),target, intent(inout) :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_remd),  optional, intent(inout) :: remd

    ! local varialbles
    type(s_enefunc_gamd), pointer :: gamd


    gamd => enefunc%gamd

    if (.not. gamd_info%gamd) then
      return
    end if

    ! initialize variables in gamd
    !
    call init_enefunc_gamd(gamd)
    enefunc%gamd_use = gamd_info%gamd
    gamd%gamd_boost  = gamd_info%boost
    gamd%sigma0_pot  = gamd_info%sigma0_pot
    gamd%sigma0_dih  = gamd_info%sigma0_dih
    gamd%ene_pot_max = gamd_info%pot_max
    gamd%ene_pot_min = gamd_info%pot_min
    gamd%ene_pot_ave = gamd_info%pot_ave
    gamd%ene_pot_dev = gamd_info%pot_dev
    gamd%ene_dih_max = gamd_info%dih_max
    gamd%ene_dih_min = gamd_info%dih_min
    gamd%ene_dih_ave = gamd_info%dih_ave
    gamd%ene_dih_dev = gamd_info%dih_dev

    ! setup boost type
    !
    if (gamd_info%boost_type == GamdBoostDihedral) then
      gamd%boost_dih = .true.
    else  if (gamd_info%boost_type == GamdBoostPotential) then
      gamd%boost_pot = .true.
    else if (gamd_info%boost_type == GamdBoostDual) then
      gamd%boost_dual = .true.
    end if

    ! setup threshold type
    !
    if (gamd_info%thresh_type == GamdThreshLower) then
      gamd%thresh_lower = .true.
    else if (gamd_info%thresh_type == GamdThreshUpper) then
      gamd%thresh_upper = .true.
    end if

    ! setup update period
    !
    if (gamd_info%update_period > 0) then

      if (mod(dynamics%nsteps, gamd_info%update_period) /= 0) then
        call error_msg('Setup_Gamd> mod(nsteps,update_period) must be zero')
      else
        gamd%update_period = gamd_info%update_period
      end if

      if (present(remd)) then
        gamd%update_period = 0
        gamd%gamd_stat = .false.
        if (main_rank) then
          write(MsgOut,'(A)') 'Setup_Gamd> Parameter update is not allowed in &
                               &GaMD-REUS'
          write(MsgOut,'(A)') 'Setup_Gamd> GaMD parameters will not be updated'
        end if
      end if

    else if (gamd_info%update_period == 0) then

      gamd%update_period = 0
      gamd%gamd_stat = .false.
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Gamd> GaMD parameters will not be updated'
      end if

    else

      call error_msg('Setup_Gamd> update_period should be non-negative integer')

    end if

    ! write the summary of setup GaMD parameters
    ! (energy thresholds and force constants of boost potential)
    !
    if (main_rank) write(MsgOut,'(A)') 'Setup_Gamd> Gamd initial parameters'

    call alloc_enefunc_gamd(gamd, molecule, enefunc)

    if (gamd%boost_pot) then

      call setup_enefunc_gamd(gamd)

      if (main_rank) then
        write(MsgOut,'(A20,F14.5)') '  ene_pot_th  = ', gamd%ene_pot_th
        write(MsgOut,'(A20,F14.5)') '  k_pot       = ', gamd%k_pot
      end if

    else if (gamd%boost_dih) then

      call setup_enefunc_gamd(gamd)

      if (main_rank) then
        write(MsgOut,'(A20,F14.5)') '  ene_dih_th  = ', gamd%ene_dih_th
        write(MsgOut,'(A20,F14.5)') '  k_dih       = ', gamd%k_dih
      end if

    else if (gamd%boost_dual) then

      call setup_enefunc_gamd(gamd)

      if (main_rank) then
        write(MsgOut,'(A20,F14.5)') '  ene_pot_th  = ', gamd%ene_pot_th
        write(MsgOut,'(A20,F14.5)') '  k_pot       = ', gamd%k_pot
        write(MsgOut,'(A20,F14.5)') '  ene_dih_th  = ', gamd%ene_dih_th
        write(MsgOut,'(A20,F14.5)') '  k_dih       = ', gamd%k_dih
      end if

    end if

    if (main_rank) write(MsgOut,'(A)') ''

    if (.not. gamd%gamd_boost) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Gamd> Gamd boosting is not applied.'
        write(MsgOut,'(A)') ''
      end if
    end if

    if (gamd_info%update_period > 0) then
      gamd%ene_pot_ave = 0.0_wp
      gamd%ene_pot_dev = 0.0_wp
      gamd%ene_dih_ave = 0.0_wp
      gamd%ene_dih_dev = 0.0_wp
    end if

    return
  end subroutine setup_gamd

end module at_gamd_mod
