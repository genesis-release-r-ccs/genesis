!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_energy_mod
!> @brief   Computes energy.
!! @authors Jaewoon Jung(JJ), Yuji Sugita (YS), Takaharu Mori (TM),
!!          Chigusa Kobayashi (CK), Norio Takase (NT)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_energy_mod

  use cg_energy_bases_mod
  use cg_energy_dihedrals_mod
  use cg_energy_angles_mod
  use cg_energy_bonds_mod
  use cg_energy_go_mod
  use cg_energy_nonlocal_mod
  use cg_energy_martini_mod
  use cg_energy_restraints_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_energy_str_mod
  use cg_dynvars_str_mod
  use cg_restraints_str_mod
  use cg_domain_str_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
  use string_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  integer, parameter      :: FakeDefault      = 0

  ! structures
  type, public :: s_ene_info
    integer               :: forcefield       = ForcefieldRESIDCG
    integer               :: electrostatic    = ElectrostaticCUTOFF
    real(wp)              :: dielec_const     = 1.0_wp
    real(wp)              :: epsilon_rf       = 0.0_wp
    integer               :: output_style     = OutputStyleGENESIS
    logical               :: user_def_table   = .false.
    logical               :: assign_force_max = .false.
    real(wp)              :: upper_force_value= 100.0_wp

    real(wp)              :: cg_cellsize_min              = 12.0_wp
    real(wp)              :: cg_cutoffdist_vdw            = 12.0_wp
    real(wp)              :: cg_switchdist_vdw            = 10.0_wp
    real(wp)              :: cg_pairlistdist_vdw          = 13.5_wp
    real(wp)              :: cg_cutoffdist_ele            = 52.0_wp
    real(wp)              :: cg_cutoffdist_126            = 39.0_wp
    real(wp)              :: cg_cutoffdist_DNAbp          = 18.0_wp
    real(wp)              :: cg_pairlistdist_ele          = 57.0_wp
    real(wp)              :: cg_pairlistdist_126          = 44.0_wp
    real(wp)              :: cg_pairlistdist_WF           = 25.0_wp
    real(wp)              :: cg_pairlistdist_PWMcos       = 23.0_wp
    real(wp)              :: cg_pairlistdist_DNAbp        = 23.0_wp
    real(wp)              :: cg_pairlistdist_exv          = 15.0_wp
    real(wp)              :: cg_sol_temperature           = 300.0_wp 
    real(wp)              :: cg_sol_ionic_strength        = 0.15_wp
    real(wp)              :: cg_pro_DNA_ele_scale_Q       = -1.0_wp
    real(wp)              :: cg_PWMcos_sigma              = 1.0_wp
    real(wp)              :: cg_PWMcos_phi                = 10.0_wp
    real(wp)              :: cg_PWMcosns_sigma            = 1.0_wp
    real(wp)              :: cg_PWMcosns_phi              = 10.0_wp
    real(wp)              :: cg_IDR_HPS_epsilon           = 0.2_wp
    real(wp)              :: cg_exv_sigma_scaling         = 1.0_wp
    logical               :: cg_infinite_DNA              = .false.
    real(wp)              :: cg_IDR_MPIPI_cutoff_scale    = 3.0_wp

  end type s_ene_info

  ! varibles
  logical, save           :: etitle = .true.

  ! subroutines
  public  :: show_ctrl_energy
  public  :: read_ctrl_energy
  public  :: compute_energy
  public  :: output_energy
  private :: output_energy_genesis
  private :: output_energy_charmm
  private :: output_energy_namd
  private :: output_energy_gromacs

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_energy
  !> @brief        Displays ENERGY section usage.
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_energy(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = RESIDCG   # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = CUTOFF    # [CUTOFF,PME]'
        write(MsgOut,'(A)') '# output_style  = GENESIS   # format of energy output [GENESIS,CHARMM,NAMD,GROMACS]'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = RESIDCG   # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = CUTOFF    # [CUTOFF,PME]'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_energy
  !> @brief        Reads ENERGY section in the control file.
  !! @authors      YS, TM, JJ
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   ene_info : ENERGY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_energy(handle, ene_info)

    ! parameters
    character(*),            parameter     :: Section = 'Energy'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_ene_info),        intent(inout) :: ene_info

    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'forcefield',    &
                               ene_info%forcefield, ForceFieldTypes)
    call read_ctrlfile_type   (handle, Section, 'electrostatic', &
                               ene_info%electrostatic, ElectrostaticTypes)
    call read_ctrlfile_real   (handle, Section, 'dielec_const',  &
                               ene_info%dielec_const)
    call read_ctrlfile_real   (handle, Section, 'epsilon_rf',    &
                               ene_info%epsilon_rf)
    call read_ctrlfile_type   (handle, Section, 'output_style',  &
                               ene_info%output_style, OutputStyleTypes)
    call read_ctrlfile_logical(handle, Section, 'user_def_table',   &
                               ene_info%user_def_table)

    call read_ctrlfile_real   (handle, Section, 'cg_cellsize_min',        &
        ene_info%cg_cellsize_min)
    call read_ctrlfile_real   (handle, Section, 'cg_cutoffdist_ele',      &
        ene_info%cg_cutoffdist_ele)
    call read_ctrlfile_real   (handle, Section, 'cg_cutoffdist_vdw',      &
        ene_info%cg_cutoffdist_vdw)
    call read_ctrlfile_real   (handle, Section, 'cg_switchdist_vdw',      &
        ene_info%cg_switchdist_vdw)
    call read_ctrlfile_real   (handle, Section, 'cg_cutoffdist_126',      &
        ene_info%cg_cutoffdist_126)
    call read_ctrlfile_real   (handle, Section, 'cg_cutoffdist_DNAbp',    &
        ene_info%cg_cutoffdist_DNAbp)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_ele',    &
        ene_info%cg_pairlistdist_ele)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_vdw',    &
        ene_info%cg_pairlistdist_vdw)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_126',    &
        ene_info%cg_pairlistdist_126)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_WF',    &
        ene_info%cg_pairlistdist_WF)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_PWMcos', &
        ene_info%cg_pairlistdist_PWMcos)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_DNAbp',  &
        ene_info%cg_pairlistdist_DNAbp)
    call read_ctrlfile_real   (handle, Section, 'cg_pairlistdist_exv',    &
        ene_info%cg_pairlistdist_exv)

    call read_ctrlfile_real   (handle, Section, 'cg_sol_temperature',      &
        ene_info%cg_sol_temperature)
    call read_ctrlfile_real   (handle, Section, 'cg_sol_ionic_strength',   &
        ene_info%cg_sol_ionic_strength)
    call read_ctrlfile_real   (handle, Section, 'cg_pro_DNA_ele_scale_Q',  &
        ene_info%cg_pro_DNA_ele_scale_Q)
    call read_ctrlfile_real   (handle, Section, 'cg_PWMcos_sigma',         &
        ene_info%cg_PWMcos_sigma)
    call read_ctrlfile_real   (handle, Section, 'cg_PWMcos_phi',           &
        ene_info%cg_PWMcos_phi)
    call read_ctrlfile_real   (handle, Section, 'cg_PWMcosns_sigma',       &
        ene_info%cg_PWMcosns_sigma)
    call read_ctrlfile_real   (handle, Section, 'cg_PWMcosns_phi',         &
        ene_info%cg_PWMcosns_phi)
    call read_ctrlfile_real   (handle, Section, 'cg_IDR_HPS_epsilon',      &
        ene_info%cg_IDR_HPS_epsilon)
    call read_ctrlfile_real   (handle, Section, 'cg_exv_sigma_scaling',    &
        ene_info%cg_exv_sigma_scaling)
    call read_ctrlfile_logical(handle, Section, 'cg_infinite_DNA',         &
        ene_info%cg_infinite_DNA)
    call read_ctrlfile_logical(handle, Section, 'assign_force_max',        &
        ene_info%assign_force_max)
    call read_ctrlfile_real   (handle, Section, 'upper_force_value',       &
        ene_info%upper_force_value)
    call read_ctrlfile_real (handle, Section, 'cg_IDR_MPIPI_cutoff_scale', &
        ene_info%cg_IDR_MPIPI_cutoff_scale)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)')  'Read_Ctrl_Energy> Parameters of Energy Calculations'
      write(MsgOut,'(A20,A10)')                            &
            '  forcefield      = ',                        &
                trim(ForceFieldTypes(ene_info%forcefield))

      write(MsgOut,'(A20,A10)')                             &
            '  output_style    = ',                         &
            trim(OutputStyleTypes(ene_info%output_style))

      write(MsgOut,'(A)') ' '

    end if

    return
  
  end subroutine read_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy
  !> @brief        Computes potential energy.
  !! @authors      JJ
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[inout] enefunc       : potential energy functions information
  !! @param[inout] domain        : domain information
  !! @param[inout] dynvars       : dynvars information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy(pairlist, boundary,      &
                            npt, reduce, nonb_ene,   &
                            enefunc, domain, dynvars)

    ! formal arguments
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    real(dp)                 :: volume

    real(wip),       pointer :: coord(:,:), vel(:,:)
    real(wip),       pointer :: force(:,:), force_long(:,:)
    real(wp),        pointer :: coord_pbc(:,:)
    real(wp),        pointer :: force_omp(:,:,:)
    real(dp),        pointer :: virial(:,:)
    real(dp),        pointer :: virial_long(:,:), virial_ext(:,:)
    integer,         pointer :: natom(:)

    coord          => domain%coord
    coord_pbc      => domain%translated
    force          => domain%force
    force_long     => domain%force_long
    force_omp      => domain%force_omp
    vel            => domain%velocity
    virial         => dynvars%virial
    virial_long    => dynvars%virial_long
    virial_ext     => dynvars%virial_extern

    call timer(TimerEnergy, TimerOn)

    select case (enefunc%forcefield)

    case (ForcefieldAAGO, ForcefieldCAGO, ForcefieldRESIDCG)

      call compute_energy_go(domain, enefunc, pairlist, boundary, coord,       &
                             npt, reduce, nonb_ene, dynvars%energy, coord_pbc, &
                             force, force_omp, virial, virial_ext)

    case (ForceFieldGROMartini)

      call compute_energy_martini(domain, enefunc, pairlist, boundary, coord,  &
                             npt, reduce, nonb_ene, dynvars%energy, coord_pbc, &
                             force, force_omp, virial, virial_ext)

    end select

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy
  !> @brief        Outputs energy.
  !! @authors      YS
  !! @param[in]    step    : step 
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(in)    :: energy


    if (.not. main_rank) return

    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_energy_genesis(step, enefunc, energy)

    case (OutputStyleCHARMM)

      call output_energy_charmm(step, enefunc, energy)

    case (OutputStyleNAMD)

      call output_energy_namd(step, energy)

    case (OutputStyleGROMACS)

      call output_energy_gromacs(step, energy)

    end select

    return

  end subroutine output_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_go
  !> @brief        Computes potential energy with GROMACS-AMBER force field.
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[inout] energy        : energy information
  !! @param[inout] force         : forces of target systems
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_go(domain, enefunc, pairlist, boundary, coord, &
                               npt, reduce, nonb_ene, energy, coord_pbc,   &
                               force, force_omp, virial, virial_ext)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:)
    real(wip),               intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    integer,         pointer :: num_atom(:)
    real(wp),        pointer :: trans(:,:)
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eimprop_omp (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: econtact_omp(nthread)
    real(dp)                 :: enoncontact_omp (nthread)
    real(dp)                 :: ebase_omp   (nthread)
    real(dp)                 :: estack_omp  (nthread)
    real(dp)                 :: ednaexv_omp (nthread)
    real(dp)                 :: epwmcos_omp (nthread)
    real(dp)                 :: epwmcosns_omp (nthread)
    real(dp)                 :: eexv_omp    (nthread)
    real(dp)                 :: ekh_omp     (nthread)
    real(dp)                 :: eidrkh_omp  (nthread)
    real(dp)                 :: eidrhps_omp (nthread)
    real(dp)                 :: eidrmpipi_omp (nthread)
    real(wp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, start_i
    integer                  :: ilist

    ! pointer
    !
    num_atom  => domain%num_atom
    trans     => domain%trans_vec

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    virial     (1:3,1:3)               = 0.0_dp
    virial_ext (1:3,1:3)               = 0.0_dp

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    eurey_omp     (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    eimprop_omp   (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp
    econtact_omp  (1:nthread)         = 0.0_dp
    enoncontact_omp(1:nthread)        = 0.0_dp
    ebase_omp     (1:nthread)         = 0.0_dp
    estack_omp    (1:nthread)         = 0.0_dp
    ednaexv_omp   (1:nthread)         = 0.0_dp
    epwmcos_omp   (1:nthread)         = 0.0_dp
    epwmcosns_omp (1:nthread)         = 0.0_dp
    elec_omp      (1:nthread)         = 0.0_dp
    eexv_omp      (1:nthread)         = 0.0_dp
    ekh_omp       (1:nthread)         = 0.0_dp
    eidrkh_omp    (1:nthread)         = 0.0_dp
    eidrhps_omp   (1:nthread)         = 0.0_dp
    eidrmpipi_omp (1:nthread)         = 0.0_dp
 
    !$omp parallel private(id, i, start_i, ix, ilist)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      force_omp(i,1,id+1) = 0.0_wp
      force_omp(i,2,id+1) = 0.0_wp
      force_omp(i,3,id+1) = 0.0_wp
    end do

    do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
      force(i,1) = 0.0_wip
      force(i,2) = 0.0_wip
      force(i,3) = 0.0_wip
    end do

    ! pbc coordinates
    !
    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
        coord_pbc(i,1) = real(coord(i,1),wp) + trans(i,1)
        coord_pbc(i,2) = real(coord(i,2),wp) + trans(i,2)
        coord_pbc(i,3) = real(coord(i,3),wp) + trans(i,3)
      end do
    else if (boundary%type == boundaryTypeNOBC) then
      do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
        coord_pbc(i,1) = real(coord(i,1),wp)
        coord_pbc(i,2) = real(coord(i,2),wp)
        coord_pbc(i,3) = real(coord(i,3),wp)
      end do
    end if

    !$omp end parallel
    
    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord_pbc, &
                             force_omp, ebond_omp, virial_omp)

    ! angle energy
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then
      call compute_energy_flexible_angle(domain, enefunc, coord_pbc, &
                                 force_omp, eangle_omp, virial_omp)
      call compute_energy_local_angle(domain, enefunc, coord_pbc, &
                                  force_omp, eangle_omp, virial_omp)
    end if
    call compute_energy_angle(domain, enefunc, coord_pbc, &
                              force_omp, eangle_omp, eurey_omp, virial_omp)

    ! dihedral energy
    !
    if (enefunc%forcefield == ForcefieldRESIDCG) then

      call compute_energy_flexible_dihed(domain, enefunc, coord_pbc, &
                                 force_omp, edihed_omp, virial_omp)
      call compute_energy_local_dihed(domain, enefunc, coord_pbc, &
                                  force_omp, edihed_omp, virial_omp)
    end if

    call compute_energy_dihed(domain, enefunc, coord_pbc, &
                              force_omp, edihed_omp, virial_omp)

    ! base stacking energy
    !
    if (enefunc%forcefield == ForceFieldRESIDCG .and. &
        enefunc%num_base_stack_all > 0) then
      call compute_energy_DNA_base_stacking(domain, enefunc, coord_pbc, &
                                            force_omp, estack_omp,      &
                                            virial_omp)
    end if

    ! contact energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%forcefield == ForcefieldAAGO) then

        call error_msg('Compute_Energy_Go: PBC is not supported for AAGO.')

      else if (enefunc%forcefield == ForcefieldCAGO) then

        
      else if (enefunc%forcefield == ForcefieldRESIDCG) then

        call compute_energy_contact_1210(domain, enefunc, coord_pbc,  &
                                         force_omp, econtact_omp,     &
                                         enoncontact_omp, virial_omp)

        
      else if (enefunc%forcefield == ForcefieldKBGO) then

      
      end if

    case default

      if (enefunc%forcefield == ForceFieldAAGO) then

        call compute_energy_contact_126(domain, enefunc, coord_pbc,  &
                                        force_omp, econtact_omp,     &
                                        enoncontact_omp)

      else if (enefunc%forcefield == ForcefieldCAGO) then

      else if (enefunc%forcefield == ForcefieldRESIDCG) then

        call compute_energy_contact_1210(domain, enefunc, coord_pbc, &
                                         force_omp, econtact_omp,    &
                                         enoncontact_omp, virial_omp)

      else if (enefunc%forcefield == ForcefieldKBGO) then

      end if

    end select

    ! non-contact energy
    !
    if (enefunc%forcefield == ForceFieldKBGO .and. &
        boundary%type == BoundaryTypePBC) then
      call error_msg('Compute_Energy_Go: PBC is not supported for KBGO model.')

    else if (enefunc%forcefield == ForcefieldRESIDCG) then

      ! ~CG~ protein-IDR: HPS model
      if ( enefunc%cg_IDR_HPS_calc ) then
        if (nonb_ene) then
          call compute_energy_CG_IDR_HPS(domain, enefunc, pairlist, &
                               coord_pbc, force_omp, eidrhps_omp, virial_omp)
        else
          call compute_force_CG_IDR_HPS(domain, enefunc, pairlist, &
                               coord_pbc, force_omp, eidrhps_omp, virial_omp)
        end if
      end if

      ! ! ~CG~ protein-IDR: KH model
      if ( enefunc%cg_IDR_KH_calc ) then
        call compute_energy_CG_IDR_KH(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, eidrkh_omp, virial_omp)
      end if

      ! ! ~CG~ protein-protein: KH model
      if ( enefunc%cg_KH_calc ) then
        call compute_energy_CG_KH(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, ekh_omp, virial_omp)
      end if

      ! ~CG~ protein-IDR: MPIPI model
      if ( enefunc%cg_IDR_MPIPI_calc ) then
        if (nonb_ene) then
          call compute_energy_CG_IDR_MPIPI(domain, enefunc, pairlist, &
              coord_pbc, force_omp, eidrmpipi_omp, virial_omp)
        else
          call compute_force_CG_IDR_MPIPI(domain, enefunc, pairlist, &
              coord_pbc, force_omp, eidrmpipi_omp, virial_omp)
        end if
      end if

    end if


    ! -------------------------------------
    ! CG electrostatics: Debye-Huckel model
    ! -------------------------------------
    !
    if (enefunc%forcefield == ForceFieldRESIDCG) then
      if ( enefunc%cg_ele_calc ) &
        call compute_energy_CG_ele(domain, enefunc, pairlist, &
                                   coord_pbc, force_omp, elec_omp, virial_omp)
    end if

    ! ---------------------------------------------
    ! CG model: general excluded volume interaction
    ! ---------------------------------------------
    !
    if (enefunc%forcefield == ForceFieldRESIDCG) then

      call compute_energy_general_exv_AICG2P(domain, enefunc, pairlist, &
                               coord_pbc, force_omp, eexv_omp, virial_omp)

    end if

    ! ---------------------------
    ! CG protein-DNA interactions
    ! ---------------------------
    !
    if (enefunc%forcefield == ForceFieldRESIDCG) then
      if ( enefunc%cg_pwmcos_calc ) then
        call compute_energy_pwmcos(domain, enefunc, pairlist,     &
                               coord_pbc, force_omp, epwmcos_omp, &
                               virial_omp)
      end if
      if ( enefunc%cg_pwmcosns_calc ) then
        call compute_energy_pwmcosns(domain, enefunc, pairlist,    &
                               coord_pbc, force_omp, epwmcosns_omp,&
                               virial_omp)
      end if
    end if

    ! -----------------------------
    ! CG DNA non-local interactions
    ! -----------------------------
    !
    if (enefunc%forcefield == ForceFieldRESIDCG) then

      if ( enefunc%cg_DNA_exv_calc ) &
        call compute_energy_DNA_exv(domain, enefunc, pairlist, &
                                    coord_pbc, force_omp,      &
                                    ednaexv_omp, virial_omp)

      if ( enefunc%cg_DNA_base_pair_calc ) &
        call compute_energy_DNA_base_pairing(domain, enefunc, pairlist, &
                                    coord_pbc, force_omp, ebase_omp,    &
                                    virial_omp)

    end if

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain,                 &
                                     enefunc, coord, force_omp,              &
                                     virial_omp, virial_ext_omp,             &
                                     eposi_omp, energy%restraint_rmsd,       &
                                     energy%rmsd, energy%restraint_distance, &
                                     energy%restraint_emfit)

!   ! gather values
!   !
    !$omp parallel default(shared) private(id, i, ix, force_tmp, start_i) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ix = 1, nthread

      !$omp do private(i)
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        force(i,1) = force(i,1) + force_omp(i,1,ix)
        force(i,2) = force(i,2) + force_omp(i,2,ix)
        force(i,3) = force(i,3) + force_omp(i,3,ix)
      end do
      !$omp end do
      !$omp barrier
    end do
    if (enefunc%assign_force_max) then
      !$omp do 
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        if (force(i,1) > enefunc%upper_force_value) &
          force(i,1) = enefunc%upper_force_value
        if (force(i,1) < -enefunc%upper_force_value) &
          force(i,1) = -enefunc%upper_force_value
        if (force(i,2) > enefunc%upper_force_value) &
          force(i,2) = enefunc%upper_force_value
        if (force(i,2) < -enefunc%upper_force_value) &
          force(i,2) = -enefunc%upper_force_value
        if (force(i,3) > enefunc%upper_force_value) &
          force(i,3) = enefunc%upper_force_value
        if (force(i,3) < -enefunc%upper_force_value) &
          force(i,3) = -enefunc%upper_force_value
      end do
    end if

    !$omp end parallel

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%base_stacking      = energy%base_stacking      + estack_omp(id)
      energy%cg_DNA_exv         = energy%cg_DNA_exv         + ednaexv_omp(id)
      energy%cg_KH_inter_pro    = energy%cg_KH_inter_pro    + ekh_omp(id)
      energy%cg_IDR_KH          = energy%cg_IDR_KH          + eidrkh_omp(id)
      energy%cg_IDR_HPS         = energy%cg_IDR_HPS         + eidrhps_omp(id)
      energy%cg_IDR_MPIPI       = energy%cg_IDR_MPIPI       + eidrmpipi_omp(id)
      energy%contact            = energy%contact            + econtact_omp(id)
      energy%noncontact         = energy%noncontact         + enoncontact_omp(id)
      energy%PWMcos             = energy%PWMcos             + epwmcos_omp(id)
      energy%PWMcosns           = energy%PWMcosns           + epwmcosns_omp(id)
      energy%cg_exv             = energy%cg_exv             + eexv_omp(id)
      energy%base_pairing       = energy%base_pairing       + ebase_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond               &
                 + energy%angle              &
                 + energy%dihedral           &
                 + energy%improper           &
                 + energy%base_stacking      &
                 + energy%cg_DNA_exv         &
                 + energy%contact            &
                 + energy%noncontact         &
                 + energy%base_pairing       &
                 + energy%restraint_position &
                 + energy%PWMcos             &
                 + energy%PWMcosns           &
                 + energy%cg_exv             &
                 + energy%electrostatic      &
                 + energy%cg_KH_inter_pro    &
                 + energy%cg_IDR_KH          &
                 + energy%cg_IDR_HPS         &
                 + energy%cg_IDR_MPIPI

    if (reduce) &
      call reduce_ene_go(energy)
    call mpi_barrier(mpi_comm_country, ierror)

    return

  end subroutine compute_energy_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_martini
  !> @brief        Computes Martini potential energy.
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[inout] energy        : energy information
  !! @param[inout] force         : forces of target systems
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_martini(domain, enefunc, pairlist, boundary,   &
                               coord, npt, reduce, nonb_ene, energy,       &
                               coord_pbc, force, force_omp, virial, virial_ext)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:)
    real(wip),               intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    integer,         pointer :: num_atom(:)
    real(wp),        pointer :: trans(:,:)
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eimprop_omp (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(wp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, start_i
    integer                  :: ilist

    ! pointer
    !
    num_atom  => domain%num_atom
    trans     => domain%trans_vec

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    virial     (1:3,1:3)               = 0.0_dp
    virial_ext (1:3,1:3)               = 0.0_dp

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp
    elec_omp      (1:nthread)         = 0.0_dp
    evdw_omp      (1:nthread)         = 0.0_dp
 
    !$omp parallel private(id, i, start_i, ix, ilist)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = 1, domain%num_atom_domain+domain%num_atom_boundary
      force_omp(i,1,id+1) = 0.0_wp
      force_omp(i,2,id+1) = 0.0_wp
      force_omp(i,3,id+1) = 0.0_wp
    end do

    do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
      force(i,1) = 0.0_wip
      force(i,2) = 0.0_wip
      force(i,3) = 0.0_wip
    end do

    ! pbc coordinates
    !
    if (boundary%type == BoundaryTypePBC) then
      do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
        coord_pbc(i,1) = real(coord(i,1),wp) + trans(i,1)
        coord_pbc(i,2) = real(coord(i,2),wp) + trans(i,2)
        coord_pbc(i,3) = real(coord(i,3),wp) + trans(i,3)
      end do
    else if (boundary%type == boundaryTypeNOBC) then
      do i = id+1, domain%num_atom_domain+domain%num_atom_boundary, nthread
        coord_pbc(i,1) = real(coord(i,1),wp)
        coord_pbc(i,2) = real(coord(i,2),wp)
        coord_pbc(i,3) = real(coord(i,3),wp)
      end do
    end if

    !$omp end parallel
    !$omp barrier
    
    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord_pbc, &
                             force_omp, ebond_omp, virial_omp)

    ! angle energy
    !
    call compute_energy_angle_g96(domain, enefunc, coord_pbc, &
                              force_omp, eangle_omp, virial_omp)

    ! dihedral energy
    !
    call compute_energy_dihed(domain, enefunc, coord_pbc, &
                              force_omp, edihed_omp, virial_omp)

    ! nonlocal energy
    !
    if (nonb_ene) then
      if (enefunc%electrostatic == ElectrostaticRF) then
        call compute_energy_martini_rf(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, elec_omp, virial_omp)
      else
        call compute_energy_martini_ele(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, elec_omp, virial_omp)
      end if
      call compute_energy_martini_vdw(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, evdw_omp, virial_omp)
    else
      if (enefunc%electrostatic == ElectrostaticRF) then
!       call compute_energy_martini_rf(domain, enefunc, pairlist, &
!                            coord_pbc, force_omp, elec_omp, virial_omp)
        call compute_force_martini_rf(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, virial_omp)
      else
        call compute_force_martini_ele(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, virial_omp)
      end if
!     call compute_energy_martini_vdw(domain, enefunc, pairlist, &
!                            coord_pbc, force_omp, evdw_omp, virial_omp)
      call compute_force_martini_vdw(domain, enefunc, pairlist, &
                             coord_pbc, force_omp, virial_omp)
    end if

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain,                 &
                                     enefunc, coord, force_omp,              &
                                     virial_omp, virial_ext_omp,             &
                                     eposi_omp, energy%restraint_rmsd,       &
                                     energy%rmsd, energy%restraint_distance, &
                                     energy%restraint_emfit)

    ! gather values
    !
    !$omp parallel default(shared) private(id, i, ix, force_tmp, start_i) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do ix = 1, nthread

      !$omp do private(i)
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        force(i,1) = force(i,1) + force_omp(i,1,ix)
        force(i,2) = force(i,2) + force_omp(i,2,ix)
        force(i,3) = force(i,3) + force_omp(i,3,ix)
      end do
      !$omp end do
      !$omp barrier
    end do
    if (enefunc%assign_force_max) then
      !$omp do 
      do i = 1, domain%num_atom_domain+domain%num_atom_boundary
        if (force(i,1) > enefunc%upper_force_value) &
          force(i,1) = enefunc%upper_force_value
        if (force(i,1) < -enefunc%upper_force_value) &
          force(i,1) = -enefunc%upper_force_value
        if (force(i,2) > enefunc%upper_force_value) &
          force(i,2) = enefunc%upper_force_value
        if (force(i,2) < -enefunc%upper_force_value) &
          force(i,2) = -enefunc%upper_force_value
        if (force(i,3) > enefunc%upper_force_value) &
          force(i,3) = enefunc%upper_force_value
        if (force(i,3) < -enefunc%upper_force_value) &
          force(i,3) = -enefunc%upper_force_value
      end do
    end if

    !$omp end parallel

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    energy%electrostatic = energy%electrostatic + enefunc%rf_self

    ! total energy
    !
    energy%total = energy%bond               &
                 + energy%angle              &
                 + energy%dihedral           &
                 + energy%van_der_waals      &
                 + energy%electrostatic       

    if (reduce) &
      call reduce_ene(energy, virial)
    call mpi_barrier(mpi_comm_country, ierror)

    return

  end subroutine compute_energy_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_genesis
  !> @brief        Outputs energy in GENESIS style.
  !! @authors      TM, CK
  !! @param[in]    step    : step 
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_genesis(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(in)    :: energy

    ! local variables
    integer,parameter        :: clength=16, flength=4
    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, frmt_res*10, rfrmt*7
    character                :: rfrmt_cont*9,frmt_cont*7
    real(dp)                 :: values(999)
    real(dp)                 :: ene_restraint


    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    write(category(ifm),frmt) 'ENERGY'
    values(ifm) = energy%total
    ifm = ifm+1

    if (enefunc%num_bonds + enefunc%num_bonds_quartic > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = energy%bond
      ifm = ifm+1
    endif

    if (enefunc%num_angles + enefunc%num_angflex > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = energy%urey_bradley
        ifm = ifm+1
      endif
    endif

    if (enefunc%num_dihe_all > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = energy%dihedral
      ifm = ifm+1
    endif

    if (enefunc%num_base_stack_all > 0) then
      write(category(ifm),frmt) 'BASE STACKING'
      values(ifm) = energy%base_stacking
      ifm = ifm+1
    end if

    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldRESIDCG) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = energy%noncontact
      ifm = ifm+1

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldRESIDCG) then
        ! if ( enefunc%num_base_stack > 0 ) then
        !   write(category(ifm),frmt) 'BASE_STACK'
        !   values(ifm) = energy%base_stacking
        !   ifm = ifm+1
        ! end if
        if ( enefunc%cg_DNA_base_pair_calc ) then
          write(category(ifm),frmt) 'BASE_PAIRING'
          values(ifm) = energy%base_pairing
          ifm = ifm+1
        endif

        if ( enefunc%cg_DNA_exv_calc ) then
          write(category(ifm),frmt) 'DNA_exv'
          values(ifm) = energy%cg_DNA_exv
          ifm = ifm+1
        end if

        if ( enefunc%cg_IDR_HPS_calc ) then
          write(category(ifm),frmt) 'IDR_HPS'
          values(ifm) = energy%cg_IDR_HPS
          ifm = ifm+1
        end if

        if ( enefunc%cg_IDR_MPIPI_calc ) then
          write(category(ifm),frmt) 'IDR_MPIPI'
          values(ifm) = energy%cg_IDR_MPIPI
          ifm = ifm+1
        end if

        if ( enefunc%cg_IDR_KH_calc ) then
          write(category(ifm),frmt) 'IDR_KH'
          values(ifm) = energy%cg_IDR_KH
          ifm = ifm+1
        end if

        if ( enefunc%cg_KH_calc ) then
          write(category(ifm),frmt) 'pro_pro_KH'
          values(ifm) = energy%cg_KH_inter_pro
          ifm = ifm+1
        end if

        if ( enefunc%cg_pwmcos_calc ) then
          write(category(ifm),frmt) 'PWMcos'
          values(ifm) = energy%PWMcos
          ifm = ifm+1
        end if

        if ( enefunc%cg_pwmcosns_calc ) then
          write(category(ifm),frmt) 'PWMcosns'
          values(ifm) = energy%PWMcosns
          ifm = ifm+1
        end if

        write(category(ifm),frmt) 'CG_EXV'
        values(ifm) = energy%cg_exv
        ifm = ifm+1

      end if

    else if (enefunc%forcefield == ForcefieldSOFT) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = energy%noncontact
      ifm = ifm+1

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

    else

      write(category(ifm),frmt) 'VDWAALS'
      values(ifm) = energy%van_der_waals
      ifm = ifm+1

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

    end if

    if (enefunc%num_restraintfuncs > 0) then
      if (enefunc%restraint_rmsd) then
        write(category(ifm),frmt) 'RMSD'
        values(ifm) = energy%rmsd
        ifm = ifm+1
      endif

      if (enefunc%restraint_emfit) then
        write(category(ifm),frmt) 'EMCORR'
        values(ifm) = energy%emcorr
        ifm = ifm+1
      end if

      ene_restraint =   energy%restraint_distance &
                      + energy%restraint_position &
                      + energy%restraint_rmsd     &
                      + energy%restraint_emfit
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

    endif


    if (etitle) then

      write(MsgOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(MsgOut,frmt) category(i)
        else
          write(MsgOut,frmt_cont) category(i)
        endif
      end do

      write(MsgOut,'(A80)') ' --------------- --------------- --------------- --------------- ---------------'
      etitle = .false.
    end if

    write(MsgOut,'(6x,I10,$)') step

    do i = 1, ifm-1
      if (i == ifm-1) then
        write(MsgOut,rfrmt) values(i)
      else
        write(MsgOut,rfrmt_cont) values(i)
      endif
    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine output_energy_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_charmm
  !> @brief        Outputs energy in CHARMM style.
  !! @authors      YS, CK
  !! @param[in]    step    : step 
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_charmm(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(dp)                 :: time, totener, totke, energy_, temperature
    real(dp)                 :: grms, hfctote, hfcke, ehfcor, virke
    real(dp)                 :: hbonds, asp, user
    real(dp)                 :: imnbvdw, imelec, imhbnd, rxnfield, extelec
    real(dp)                 :: ewksum, ewself, ewexcl, ewqcor, ewutil
    real(dp)                 :: vire, viri, presse, pressi, volume
    real(dp)                 :: cdihe, cintcr, noe

    real(dp),        pointer :: bonds, angles, urey_b, dihedrals, impropers
    real(dp),        pointer :: vdwaals, elec, cmaps, disp_corr
    real(dp),        pointer :: posicon, restdist


    time        = 0.0_dp
    hbonds      = 0.0_dp
    totener     = 0.0_dp
    totke       = 0.0_dp
    energy_     = 0.0_dp
    temperature = 0.0_dp
    asp         = 0.0_dp
    user        = 0.0_dp
    imnbvdw     = 0.0_dp
    imelec      = 0.0_dp
    imhbnd      = 0.0_dp
    rxnfield    = 0.0_dp
    extelec     = 0.0_dp
    ewksum      = 0.0_dp
    ewself      = 0.0_dp
    ewexcl      = 0.0_dp
    ewqcor      = 0.0_dp
    ewutil      = 0.0_dp
    vire        = 0.0_dp
    viri        = 0.0_dp
    presse      = 0.0_dp
    pressi      = 0.0_dp
    volume      = 0.0_dp
    grms        = 0.0_dp
    hfctote     = 0.0_dp
    hfcke       = 0.0_dp
    ehfcor      = 0.0_dp
    virke       = 0.0_dp
    volume      = 0.0_dp
    noe         = 0.0_dp
    cdihe       = 0.0_dp
    cintcr      = 0.0_dp

    ! write title if necessary
    !
    if (etitle) then
      write(MsgOut,'(A)') 'Output_Energy> CHARMM_Style is used'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A79)') 'DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature'
      write(MsgOut,'(A79)') 'DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe'
      write(MsgOut,'(A79)') 'DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers'

      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(MsgOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      endif

      write(MsgOut,'(A79)') 'DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER'

      write(MsgOut,'(A79)') 'DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec'
      write(MsgOut,'(A79)') 'DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil'

      if (enefunc%restraint) then
        write(MsgOut,'(A79)') 'DYNA CONSTR:       HARMonic    CDIHedral          CIC     RESDistance       NOE'
      end if

      write(MsgOut,'(A79)') 'DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme'
      write(MsgOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
      etitle = .false.

    end if


    bonds      => energy%bond
    angles     => energy%angle
    urey_b     => energy%urey_bradley
    dihedrals  => energy%dihedral
    impropers  => energy%improper
    elec       => energy%electrostatic
    vdwaals    => energy%van_der_waals
    cmaps      => energy%cmap
    restdist   => energy%restraint_distance
    posicon    => energy%restraint_position 
    disp_corr  => energy%disp_corr_energy

    ! write energy in CHARMM-style
    !
    write(MsgOut,'(A5,I9,5F13.5)') 'DYNA>', step, time, totener, totke, energy_, temperature
    write(MsgOut,'(A14,5F13.5)')   'DYNA PROP>    ', grms, hfctote, hfcke, ehfcor, virke
    write(MsgOut,'(A14,5F13.5)')   'DYNA INTERN>  ', bonds, angles, urey_b, dihedrals, impropers

    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(MsgOut,'(A14,F13.5)')    'DYNA DISP>    ', disp_corr
    end if

    write(MsgOut,'(A14,5F13.5)')   'DYNA EXTERN>  ', vdwaals, elec, hbonds, asp, user

    write(MsgOut,'(A14,5F13.5)')   'DYNA IMAGES>  ', imnbvdw, imelec, imhbnd, rxnfield, extelec
    write(MsgOut,'(A14,5F13.5)')   'DYNA EWALD>   ', ewksum, ewself, ewexcl, ewqcor, ewutil

    if (enefunc%restraint) then
      write(MsgOut,'(A14,5F13.5)')   'DYNA CONSTR>  ', posicon, cdihe, cintcr, restdist, noe
    end if

    write(MsgOut,'(A14,5F13.5)')   'DYNA PRESS>   ', vire, viri, presse, pressi, volume
    write(MsgOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
    write(MsgOut,'(A)') ' '

    return

  end subroutine output_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_namd
  !> @brief        Outputs energy in NAMD style.
  !! @authors      YS, CK
  !! @param[in]    step    : step 
  !! @param[in]    energy : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_namd(step, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(dp)                 :: ebond, eangle
    real(dp)                 :: misc, kinetic
    real(dp)                 :: total, temp, total2, total3, tempavg
    real(dp)                 :: pressure, gpressure, volume
    real(dp)                 :: pressavg, gpressavg

    real(dp),        pointer :: edihed, eimprp
    real(dp),        pointer :: eelect, evdw
    real(dp),        pointer :: eboundary


    ! write title if necessary
    !
    if (etitle) then
      write(MsgOut,'(A)') 'Output_Energy> NAMD_Style is used'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A75)') 'ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP'
      write(MsgOut,'(A75)') '          ELECT            VDW       BOUNDARY           MISC        KINETIC'
      write(MsgOut,'(A75)') '          TOTAL           TEMP         TOTAL2         TOTAL3        TEMPAVG'
      write(MsgOut,'(A75)') '       PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG'
      write(MsgOut,'(A)') ' '
      etitle = .false.
    end if


    ebond     =  energy%bond  + energy%restraint_distance
    eangle    =  energy%angle + energy%urey_bradley
    edihed    => energy%dihedral
    eimprp    => energy%improper
    eelect    => energy%electrostatic
    evdw      => energy%van_der_waals
    eboundary => energy%restraint_position

    ! write energy in NAMD-style
    !
    write(MsgOut,'(A7,I8,4F15.4)')'ENERGY:', step, ebond, eangle, edihed,eimprp
    write(MsgOut,'(5F15.4)') eelect, evdw, eboundary, misc, kinetic
    write(MsgOut,'(5F15.4)') total, temp, total2, total3, tempavg
    write(MsgOut,'(5F15.4)') pressure, gpressure, volume, pressavg, gpressavg
    write(MsgOut,'(A)') ' '

    return

  end subroutine output_energy_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_gromacs
  !> @brief        Outputs energy in GROMACS style.
  !! @authors      NT
  !! @param[in]    step    : step
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_gromacs(step, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(dp)                 :: time, bonds, angles, urey_b
    real(dp)                 :: dihedrals, impropers, vdwaals, elec
    real(dp)                 :: restdist, posicon
    real(dp)                 :: energy_, totke, totener
    real(dp)                 :: temperature, pressi, presse
    real(dp)                 :: disp_corr


    time        = 0.0_dp
    bonds       = energy%bond          * CAL2JOU
    angles      = energy%angle         * CAL2JOU
    urey_b      = energy%urey_bradley  * CAL2JOU
    dihedrals   = energy%dihedral      * CAL2JOU
    impropers   = energy%improper      * CAL2JOU
    vdwaals     = energy%van_der_waals * CAL2JOU
    elec        = energy%electrostatic * CAL2JOU
    disp_corr   = energy%disp_corr_energy * CAL2JOU
    restdist    = energy%restraint_distance
    posicon     = energy%restraint_position+energy%restraint_rmsd
    energy_     = 0.0_dp
    totke       = 0.0_dp
    totener     = 0.0_dp
    temperature = 0.0_dp
    pressi      = 0.0_dp
    presse      = 0.0_dp


    write(MsgOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(MsgOut,'(I15,2F15.5)') &
         step, time, 0.0_dp
    write(MsgOut,'(A)') &
         ' '
    write(MsgOut,'(A)') &
         '   Energies (kJ/mol)'

    write(MsgOut,'(5A15)') &
         'Bond', 'Angle', 'Urey-bradley', 'Dihedral', 'Improper Dih.'
    write(MsgOut,'(5ES15.5E2)') &
         bonds,angles,urey_b,dihedrals,impropers

    write(MsgOut,'(4A15)') &
       'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Disper. corr.', 'Position Rest.'
    write(MsgOut,'(4ES15.5E2)') &
         vdwaals,elec,disp_corr,posicon
    write(MsgOut,'(2A15)') &
        'Potential', 'Kinetic En.'
    write(MsgOut,'(2ES15.5E2)') &
        energy_,totke


    write(MsgOut,'(5A15)') &
         'Total Energy', 'Temperature', 'Pressure(int.)', 'Pressure(ext.)'
    write(MsgOut,'(5ES15.5E2)') &
         totener,temperature,pressi,presse
    
    write(MsgOut,'(A)')  ' '

    return

  end subroutine output_energy_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    reduce_ene_go
  !> @brief        Reduces energy and virial.
  !! @authors      JJ
  !! @param[inout] energy : energy information
  !! @param[inout] virial : virial term of 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_ene_go(energy)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: before_reduce(19), after_reduce(19)

    ! Allreduce energy components
    !
    before_reduce(1)  = energy%bond
    before_reduce(2)  = energy%angle
    before_reduce(3)  = energy%dihedral
    before_reduce(4)  = energy%improper
    before_reduce(5)  = energy%base_stacking
    before_reduce(6)  = energy%contact
    before_reduce(7)  = energy%noncontact
    before_reduce(8)  = energy%restraint_position
    before_reduce(9)  = energy%electrostatic
    before_reduce(10) = energy%base_pairing
    before_reduce(11) = energy%PWMcos
    before_reduce(12) = energy%PWMcosns
    before_reduce(13) = energy%cg_DNA_exv
    before_reduce(14) = energy%cg_exv
    before_reduce(15) = energy%cg_KH_inter_pro
    before_reduce(16) = energy%cg_IDR_KH
    before_reduce(17) = energy%cg_IDR_HPS
    before_reduce(18) = energy%cg_IDR_MPIPI
    before_reduce(19) = energy%total

    call mpi_reduce(before_reduce, after_reduce, 19, mpi_real8,  &
                    mpi_sum, 0, mpi_comm_country, ierror)

    energy%bond               = after_reduce(1)
    energy%angle              = after_reduce(2)
    energy%dihedral           = after_reduce(3)
    energy%improper           = after_reduce(4)
    energy%base_stacking      = after_reduce(5)
    energy%contact            = after_reduce(6)
    energy%noncontact         = after_reduce(7)
    energy%restraint_position = after_reduce(8)
    energy%electrostatic      = after_reduce(9)
    energy%base_pairing       = after_reduce(10)
    energy%PWMcos             = after_reduce(11)
    energy%PWMcosns           = after_reduce(12)
    energy%cg_DNA_exv         = after_reduce(13)
    energy%cg_exv             = after_reduce(14)
    energy%cg_KH_inter_pro    = after_reduce(15)
    energy%cg_IDR_KH          = after_reduce(16)
    energy%cg_IDR_HPS         = after_reduce(17)
    energy%cg_IDR_MPIPI       = after_reduce(18)
    energy%total              = after_reduce(19)

#endif

    return

  end subroutine reduce_ene_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_ene
  !> @brief        Reduces energy and virial.
  !! @authors      JJ
  !! @param[inout] energy : energy information
  !! @param[inout] virial : virial term of
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_ene(energy, virial)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy
    real(dp),                intent(inout) :: virial(3,3)

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: before_reduce(23), after_reduce(23)
    integer                  :: i, j, n


    ! Allreduce virial and energy components
    !
    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        before_reduce(n) = virial(i,j)
      end do
    end do

    before_reduce(10) = energy%bond
    before_reduce(11) = energy%angle
    before_reduce(12) = energy%dihedral
    before_reduce(13) = energy%electrostatic
    before_reduce(14) = energy%van_der_waals
    before_reduce(16) = energy%restraint_position
    before_reduce(19) = energy%total

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 19, mpi_real8,  &
                       mpi_sum, mpi_comm_country, ierror)
#else
    after_reduce(1:19) = before_reduce(1:19)
#endif

    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        virial(i,j) = after_reduce(n)
      end do
    end do

    energy%bond               = after_reduce(10)
    energy%angle              = after_reduce(11)
    energy%dihedral           = after_reduce(12)
    energy%electrostatic      = after_reduce(13)
    energy%van_der_waals      = after_reduce(14)
    energy%restraint_position = after_reduce(16)
    energy%total              = after_reduce(19)

#endif

    return

  end subroutine reduce_ene

end module cg_energy_mod
