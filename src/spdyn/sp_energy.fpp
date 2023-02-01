!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_mod
!> @brief   compute energy
!! @authors Jaewoon Jung(JJ), Yuji Sugita (YS), Takaharu Mori (TM),
!!          Chigusa Kobayashi (CK), Norio Takase (NT), Hiraku Oshima (HO)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_mod

  use sp_energy_table_linear_mod
  use sp_energy_14_notable_mod
  use sp_energy_table_linear_bondcorr_mod
  use sp_energy_nonbonds_mod
  use sp_energy_restraints_mod
  use sp_energy_dihedrals_mod
  use sp_energy_angles_mod
  use sp_energy_bonds_mod
  use sp_energy_efield_mod
  use sp_energy_gamd_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_restraints_str_mod
  use sp_domain_str_mod
  use sp_experiments_mod
  use sp_enefunc_gamd_mod
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

  implicit none
  private

  integer, parameter      :: FakeDefault      = 0

  ! structures
  type, public :: s_ene_info
    integer               :: forcefield       = ForcefieldCHARMM
    character(20)         :: forcefield_char  = 'CHARMM'
    integer               :: electrostatic    = ElectrostaticPME
    integer               :: vdw              = VdwPME
    real(wp)              :: switchdist       = 10.0_wp
    real(wp)              :: cutoffdist       = 12.0_wp
    real(wp)              :: pairlistdist     = 13.5_wp
    real(wp)              :: dielec_const     = 1.0_wp
    real(wp)              :: dmin_size_cg     = 20.0_wp
    logical               :: dispersion_pme   = .false.
    logical               :: vdw_force_switch = .false.
    logical               :: vdw_shift        = .false.
    logical               :: cmap_pspline     = .false.
    logical               :: contact_check    = .false.
    logical               :: nonb_limiter     = .false.
    integer               :: structure_check  = StructureCheckNone
    logical               :: dsize_cg         = .false.
    real(wp)              :: pme_alpha ! do not initialize here
    real(wp)              :: pme_alpha_tol    = 1.0e-5
    integer               :: pme_ngrid_x      = FakeDefault
    integer               :: pme_ngrid_y      = FakeDefault
    integer               :: pme_ngrid_z      = FakeDefault
    integer               :: pme_nspline      = 4
    real(wp)              :: pme_max_spacing  = 1.2_wp
    integer               :: pme_scheme       = FFT_AutoSelect
    integer               :: nonbond_kernel   = NBK_AutoSelect
    logical               :: table            = .true.
    integer               :: table_order      = 1
    real(wp)              :: table_density    = 20.0_wp
    character(5)          :: water_model      = 'NONE'
    integer               :: output_style     = OutputStyleGENESIS
    integer               :: dispersion_corr  = Disp_corr_NONE
    real(wp)              :: minimum_contact  = 0.5_wp
    real(wp)              :: err_minimum_contact = 0.3_wp
    logical               :: use_knl_generic  = .false.
    real(wp)              :: scale_pairlist_fugaku  = 1.0_wp
    real(wp)              :: scale_pairlist_generic = 1.0_wp
    real(wp)              :: efield_x = 0.0_wp
    real(wp)              :: efield_y = 0.0_wp
    real(wp)              :: efield_z = 0.0_wp
    logical               :: efield_virial = .false.
    logical               :: efield_normal = .false.
    logical               :: vacuum           = .false.
  end type s_ene_info

  ! varibles
  logical, save           :: etitle = .true.

  ! subroutines
  public  :: show_ctrl_energy
  public  :: read_ctrl_energy
  public  :: compute_energy
  public  :: compute_energy_short
  public  :: compute_energy_long
  public  :: output_energy
  private :: compute_energy_charmm
  private :: compute_energy_amber
  private :: compute_energy_gro_amber
  private :: compute_energy_gro_martini
  private :: compute_energy_charmm_short
  private :: compute_energy_amber_short
  private :: compute_energy_gro_amber_short
  private :: compute_energy_general_long
  private :: output_energy_genesis
  private :: output_energy_charmm
  private :: output_energy_namd
  private :: output_energy_gromacs
  private :: reduce_ene
  private :: compute_stats
  ! FEP
  private :: compute_energy_charmm_fep
  private :: compute_energy_amber_fep
  private :: compute_energy_charmm_short_fep
  private :: compute_energy_amber_short_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_energy
  !> @brief        show ENERGY section usage
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

      case ('md', 'min', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = CHARMM    # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = PME       # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist    = 10.0      # switch distance'
        write(MsgOut,'(A)') 'cutoffdist    = 12.0      # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist  = 13.5      # pair-list distance'
        write(MsgOut,'(A)') '# nonbond_kernel = autoselect  # kernel type for real-space non-bonded interactions'
        write(MsgOut,'(A)') '# dielec_const  = 1.0       # dielectric constant'
        write(MsgOut,'(A)') '# vdw_force_switch = NO     # force switch option for van der Waals'
        write(MsgOut,'(A)') '# dispersion_pme   = NO     # PME for dispersion term in van der Waals'
        write(MsgOut,'(A)') '# vdw_shift        = NO     # shift option for van der Waals'
        write(MsgOut,'(A)') '# pme_scheme   = autoselect # reciprocal space calculation scheme in [PME]'
        write(MsgOut,'(A)') '# pme_alpha     = auto      # width of Gaussian distribution in [PME]'
        write(MsgOut,'(A)') '# pme_alpha_tol = 1.0e-5    # param for auto pme_alpha determination'
        write(MsgOut,'(A)') '# pme_nspline   = 4         # order of B-spline in [PME]'
        write(MsgOut,'(A)') '# pme_max_spacing  = 1.2    # Max grid spacing allowed '
        write(MsgOut,'(A)') '# table_density = 20.0      # number of bins used for lookup table'
        write(MsgOut,'(A)') '# output_style  = GENESIS   # format of energy output [GENESIS,CHARMM,NAMD,GROMACS]'
        write(MsgOut,'(A)') '# dispersion_corr = NONE    # dispersion correction [NONE,Energy,EPress]'
        if (run_mode == 'min') then
        write(MsgOut,'(A)') '# contact_check   = YES     # check atomic clash'
        write(MsgOut,'(A)') '# nonb_limiter    = YES     # avoid failure due to atomic clash'
        write(MsgOut,'(A)') '# minimum_contact = 0.5     # definition of atomic clash distance'
        end if
        write(MsgOut,'(A)') '# nonbond_kernel  = Generic # Nonbond kernel type &
                              &[autoselect, generic, kgeneric, intelknl]'
        write(MsgOut,'(A)') '# efield_x     # external electric field in x direction (V/anstrom)'
        write(MsgOut,'(A)') '# efield_y     # external electric field in y direction (V/anstrom)'
        write(MsgOut,'(A)') '# efield_z     # external electric field in z direction (V/anstrom)'
        write(MsgOut,'(A)') '# efield_virial# calculate virial from electric field in NPT'
        write(MsgOut,'(A)') '# efield_normal# normalized efield magnitue during NPT'
        write(MsgOut,'(A)') '# vacuum = NO               # vacuum option'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = CHARMM    # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = PME       # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist    = 10.0      # switch distance'
        write(MsgOut,'(A)') 'cutoffdist    = 12.0      # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist  = 13.5      # pair-list distance'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_energy
  !> @brief        read ENERGY section in the control file
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

    character(MaxLine)                     :: pme_alpha = "auto"
    real(wp)                               :: cutoff, cutoff2, mind
    real(wp)                               :: cutoff2_water
    integer                                :: cutoff_int2


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'forcefield',    &
                               ene_info%forcefield, ForceFieldTypes)
    call read_ctrlfile_string (handle, Section, 'forcefield',    &
                               ene_info%forcefield_char)
    call read_ctrlfile_type   (handle, Section, 'electrostatic', &
                               ene_info%electrostatic, ElectrostaticTypes)
    call read_ctrlfile_real   (handle, Section, 'switchdist',    &
                               ene_info%switchdist)
    call read_ctrlfile_real   (handle, Section, 'cutoffdist',    &
                               ene_info%cutoffdist)
    call read_ctrlfile_real   (handle, Section, 'pairlistdist',  &
                               ene_info%pairlistdist)
    call read_ctrlfile_real   (handle, Section, 'dmin_size_cg',  &
                               ene_info%dmin_size_cg)
    call read_ctrlfile_real   (handle, Section, 'dielec_const',  &
                               ene_info%dielec_const)
    call read_ctrlfile_logical(handle, Section, 'dispersion_pme', &
                               ene_info%dispersion_pme)
    call read_ctrlfile_logical(handle, Section, 'vdw_force_switch', &
                               ene_info%vdw_force_switch)
    call read_ctrlfile_logical(handle, Section, 'vdw_shift',     &
                               ene_info%vdw_shift)
    call read_ctrlfile_string (handle, Section, 'pme_alpha',     &
                               pme_alpha)
    call read_ctrlfile_real   (handle, Section, 'pme_alpha_tol', &
                               ene_info%pme_alpha_tol)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_x',   &
                               ene_info%pme_ngrid_x)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_y',   &
                               ene_info%pme_ngrid_y)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_z',   &
                               ene_info%pme_ngrid_z)
    call read_ctrlfile_integer(handle, Section, 'pme_nspline',   &
                               ene_info%pme_nspline)
    call read_ctrlfile_real   (handle, Section, 'pme_max_spacing', &
                               ene_info%pme_max_spacing)
    call read_ctrlfile_type   (handle, Section, 'PME_scheme',    &
                               ene_info%pme_scheme, FFT_Types)
    call read_ctrlfile_real   (handle, Section, 'table_density', &
                               ene_info%table_density)
    call read_ctrlfile_string (handle, Section, 'water_model',   &
                               ene_info%water_model)
    call read_ctrlfile_type   (handle, Section, 'output_style',  &
                               ene_info%output_style, OutputStyleTypes)
    call read_ctrlfile_type   (handle, Section, 'dispersion_corr',  &
                               ene_info%dispersion_corr, Disp_corr_Types)
    call read_ctrlfile_type   (handle, Section, 'structure_check',  &
                               ene_info%structure_check, StructureCheckTypes)
    call read_ctrlfile_logical(handle, Section, 'contact_check',  &
                               ene_info%contact_check)
    call read_ctrlfile_real   (handle, Section, 'scale_pairlist_fugaku', &
                               ene_info%scale_pairlist_fugaku)
    call read_ctrlfile_real   (handle, Section, 'scale_pairlist_generic', &
                               ene_info%scale_pairlist_generic)
    call read_ctrlfile_real   (handle, Section, 'efield_x', &
                               ene_info%efield_x)
    call read_ctrlfile_real   (handle, Section, 'efield_y', &
                               ene_info%efield_y)
    call read_ctrlfile_real   (handle, Section, 'efield_z', &
                               ene_info%efield_z)
    call read_ctrlfile_logical(handle, Section, 'efield_virial', &
                               ene_info%efield_virial)
    call read_ctrlfile_logical(handle, Section, 'efield_normal', &
                               ene_info%efield_normal)
    if (ene_info%contact_check) then
      ene_info%nonb_limiter=.true.
    end if
    call read_ctrlfile_logical(handle, Section, 'nonb_limiter',  &
                               ene_info%nonb_limiter)
    call read_ctrlfile_real   (handle, Section, 'minimum_contact', &
                               ene_info%minimum_contact)
!   call read_ctrlfile_logical(handle, Section, 'use_knl_generic', &
!                              ene_info%use_knl_generic)
    call read_ctrlfile_type   (handle, Section, 'Nonbond_kernel', &
                               ene_info%nonbond_kernel, NBK_Types)
    call read_ctrlfile_logical(handle, Section, 'vacuum',     &
                               ene_info%vacuum)

    call end_ctrlfile_section(handle)

    ! check vacuum
    !
    if (ene_info%vacuum) then
      ! Tables are not used.
      ene_info%table = .false.
      ! Switch functions are not used.
      ene_info%vdw_force_switch = .false.
    end if

    ! check table
    !
    if (ene_info%table) then
      if (ene_info%electrostatic == ElectrostaticCutoff) &
          ene_info%table_order = 3
      if (ene_info%electrostatic == ElectrostaticPME) &
          ene_info%table_order = 1
    end if

    ! error check for inputs
    !
    if (ene_info%switchdist > ene_info%cutoffdist) then
      call error_msg( &
         'Read_Ctrl_Energy> switchdist must be less than cutoffdist. '// &
         'Please set switchdist shorter than cutoffdist.')
    end if

    if (ene_info%cutoffdist >= ene_info%pairlistdist) then
      call error_msg( &
         'Read_Ctrl_Energy> cutoffdist must be less than pairlistdist. '// &
         'Please set cutoffdist shorter than pairlistdist.')
    end if

    if (ene_info%pme_nspline < 3 ) then
      call error_msg( &
         'Read_Ctrl_Energy> "pme_nspline" is too small. '// &
         'Please set pme_nspline > 3.')

    else if (mod(ene_info%pme_nspline,2) == 1 ) then
      call error_msg( &
         'Read_Ctrl_Energy> "pme_nspline" should be even (in the current version)')

    end if

    if (ene_info%water_model /= "NONE" .and. &
        ene_info%electrostatic == ElectrostaticCUTOFF) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: water_model is not available'

          ene_info%water_model = "NONE"

    end if

#ifdef USE_GPU
    if (ene_info%electrostatic == ElectrostaticCUTOFF) &
      call error_msg( &
         'Read_Ctrl_Energy> cutoff is not available with GPU. '// &
         'Please use GENESIS compiled without enabling GPU.')

    if (ene_info%water_model /= "NONE") ene_info%water_model = 'NONE'

    if (ene_info%nonb_limiter) &
      call error_msg( &
         'Read_Ctrl_Energy> nonb_limiter is not available with GPU. '// &
         'Please use GENESIS compiled without enabling GPU.')

    if (ene_info%structure_check /= StructureCheckNone) &
      call error_msg( &
         'Read_Ctrl_Energy> structure_check is not available with GPU. '// &
         'Please use GENESIS compiled without enabling GPU.')
#endif

    call tolower(pme_alpha)
    if (trim(pme_alpha) .eq. "auto") then
      ene_info%pme_alpha = get_ewald_alpha(ene_info%cutoffdist, &
                                           ene_info%pme_alpha_tol)
    else
      read(pme_alpha,*) ene_info%pme_alpha
    end if

    ! error check for each FFs
    !
    if (ene_info%forcefield == ForcefieldAMBER) then

      if (ene_info%electrostatic /= ElectrostaticPME) then
        if (.not. ene_info%vacuum) then
          call error_msg('Read_Ctrl_Energy> Electrostatic cutoff is not allowed in amber unless vacuum is used.')
        end if
      end if

      if (ene_info%switchdist /= ene_info%cutoffdist) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: switchdist should be equal to '//  &
          'cutoffdist if forcefield is AMBER'
          ene_info%switchdist = ene_info%cutoffdist
      end if

      if (ene_info%dispersion_corr /= Disp_corr_EPress) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dispersion correction should be '//&
          'set ene_info%cutoffdist if forcefield is AMBER'
          ene_info%dispersion_corr = Disp_corr_EPress
      end if

    end if

    if (ene_info%forcefield == ForcefieldCHARMM) then

      if (ene_info%dispersion_corr /= Disp_corr_NONE) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dispersion correction can not be'//&
          ' set ene_info%cutoffdist if forcefield is CHARMM'
          ene_info%dispersion_corr = Disp_corr_NONE
      end if

    end if

    if (ene_info%forcefield == ForcefieldGROMARTINI) then

      ene_info%dsize_cg  = .true.
      if (ene_info%dmin_size_cg < ene_info%pairlistdist) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dmin_size_cg is automatically'//     &
          ' set to pairlistdist if forcefield is GROMARTINI and '//         &
          'dmin_size_cg < pairlistdist'
        ene_info%dmin_size_cg = ene_info%pairlistdist

      end if

      if (ene_info%electrostatic == ElectrostaticCUTOFF .and. &
          .not. ene_info%vdw_shift) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW shift is automatically'//     &
          ' set if forcefield is GROMARTINI and CUTOFF'
          ene_info%vdw_shift = .true.
      end if

    end if

    ! In the case of LJ PME
    !
    if (ene_info%dispersion_pme) then
      if (ene_info%vdw_force_switch) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW_force_switch can not be'//     &
          ' set if dispersion_pme is YES'
        ene_info%vdw_force_switch = .false.
      end if
      if (ene_info%vdw_shift) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW_shift can not be set'//        &
          ' if dispersion_pme is YES'
        ene_info%vdw_shift = .false.
      end if
      if (ene_info%electrostatic == ElectrostaticCutoff) then
        call error_msg( &
          'Read_Ctrl_Energy>  WARNING: Electrostatic should be PME'//     &
          ' if dispersion_pme is YES')
      end if
    end if

    ! Case of force switch
    !
    if (ene_info%forcefield /= ForcefieldCHARMM) then
      if (ene_info%vdw_force_switch) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW force switch can not be'//     &
          ' set if forcefield is not CHARMM'
          ene_info%vdw_force_switch = .false.
      end if
    end if

    if (ene_info%vdw_force_switch) then
      if (ene_info%vdw_shift) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW_shift can not be set'//        &
          ' if vdW_force_switch is YES'
        ene_info%vdw_shift = .false.
      end if
    end if

    ! Case of shift
    !
    if (ene_info%forcefield /= ForcefieldGROAMBER .and. &
        ene_info%forcefield /= ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW shift can not be set'//        &
          ' if forcefield is not GROAMBER'
          ene_info%vdw_shift = .false.
      end if
    end if

    ! error check
    !

    if (.not. ene_info%table) then
      if (ene_info%nonb_limiter ) then
        if (main_rank) &
          write(MsgOut,'(A)') &
       'Read_Ctrl_Energy>  WARNING: nonb_limiter is only in table'
        ene_info%nonb_limiter = .false.
      end if
    end if

    if (ene_info%nonb_limiter .and. .not. ene_info%contact_check) then
      if (main_rank)       &
      write(MsgOut,'(A)') &
       'Read_Ctrl_Energy>  WARNING: contact_check is available in '//&
       'nonb_limiter=YES'
      ene_info%contact_check = .true.
    end if
    if (ene_info%contact_check .and.                            &
        ene_info%structure_check == StructureCheckNone) then
      ene_info%structure_check = StructureCheckFirst
    end if

    ! Decision of VDW types
    !
    if (ene_info%switchdist /= ene_info%cutoffdist) then
      ene_info%vdw = VDWSwitch
      if (ene_info%vdw_shift) ene_info%vdw = VDWShift
      if (ene_info%vdw_force_switch) ene_info%vdw = VDWFSW
    else
      ene_info%vdw = VDWCutoff
    end if
    if (ene_info%dispersion_pme) ene_info%vdw = VDWPME   

    ! error check for CUTOFF and vacuum
    !
    if (ene_info%vacuum) then
      if (ene_info%electrostatic /= ElectrostaticCutoff) then
        call error_msg( &
         'Read_Ctrl_Energy> vacuum is not allowed unless CUTOFF is used')
      end if
    end if
    
    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)')  'Read_Ctrl_Energy> Parameters of Energy Calculations'
      write(MsgOut,'(A20,A15)')                            &
            '  forcefield      = ',                        &
                trim(ForceFieldTypes(ene_info%forcefield))

      write(MsgOut,'(A20,F15.3)')                          &
            '  switchdist      = ', ene_info%switchdist     
      write(MsgOut,'(A20,F15.3)')                          &
            '  cutoffdist      = ', ene_info%cutoffdist
      write(MsgOut,'(A20,F15.3)')                          &
            '  pairlistdist    = ', ene_info%pairlistdist
      write(MsgOut,'(A20,F15.3)')                          &
            '  dielec_const    = ', ene_info%dielec_const

      if (ene_info%forcefield == ForcefieldGROMARTINI) then
        write(MsgOut,'(A20,F15.3)')                &
             '  dmin_size_cg    = ', ene_info%dmin_size_cg
      end if

      write(MsgOut,'(A20,A15)')                            &
            '  VDW type        = ', trim(VDW_types(ene_info%vdw))

      ! PME Case
      !
      if (ene_info%electrostatic == ElectrostaticPME) then

        write(MsgOut,'(A20,A15)')                            &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))
        if (ene_info%pme_ngrid_x /= FakeDefault .or.  &
            ene_info%pme_ngrid_y /= FakeDefault .or.  &
            ene_info%pme_ngrid_z /= FakeDefault) then
          write(MsgOut,'(A20,3I5)')                            &
                '  pme_ngrid(x,y,z)= ', ene_info%pme_ngrid_x   &
                                      , ene_info%pme_ngrid_y   &
                                      , ene_info%pme_ngrid_z
        end if
        write(MsgOut,'(A20,I15)')                            &
              '  pme_nspline     = ', ene_info%pme_nspline
        if (trim(pme_alpha) == "auto") then
          write(MsgOut,'(A20,F15.5)')                        &
              '  pme_alpha       = ', ene_info%pme_alpha
          write(MsgOut,'(A20,F15.5)')                        &
              '  pme_alpha_tol   = ', ene_info%pme_alpha_tol
        else
          write(MsgOut,'(A20,F15.5)')                &
              '  pme_alpha       = ', ene_info%pme_alpha
        end if
        write(MsgOut,'(A20,A15)')                            &
              '  pme_scheme      = ', trim(FFT_types(ene_info%pme_scheme))

        write(MsgOut,'(A20,A15)')                            &
              '  nonbond_kernel  = ', trim(NBK_types(ene_info%nonbond_kernel))

      ! CUTOFF Electrostatic
      !
      else if (ene_info%electrostatic == ElectrostaticCutoff) then

        write(MsgOut,'(A20,A15)')                            &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))

      end if

      if (ene_info%table) then
        write(MsgOut,'(A20,I15)')   '  table_order     = ', &
             ene_info%table_order
        write(MsgOut,'(A20,F15.3)') '  table_density   = ', &
             ene_info%table_density
      end if

      write(MsgOut,'(A20,A15)')                             &
            '  output_style    = ',                         &
            trim(OutputStyleTypes(ene_info%output_style))

      if (ene_info%dispersion_corr == Disp_corr_NONE) then
        write(MsgOut,'(A)') '  dispersion_corr =            none'
      else if (ene_info%dispersion_corr == Disp_corr_Energy) then
        write(MsgOut,'(A)') '  dispersion_corr =          energy'
      else if (ene_info%dispersion_corr == Disp_corr_EPress) then
        write(MsgOut,'(A)') '  dispersion_corr =          epress'
      end if

      if (ene_info%nonb_limiter) then
        write(MsgOut,'(A)') '  nonb_limiter    =             yes'
        write(MsgOut,'(A,F15.3)') ' minimum_contact  = ', &
             ene_info%minimum_contact
      else
        write(MsgOut,'(A)') '  nonb_limiter    =              no'
      end if
      if (ene_info%contact_check) then
        write(MsgOut,'(A)') '  contact_check   =             yes'
        write(MsgOut,'(A,F15.3)') '  minimum_contact = ', &
             ene_info%minimum_contact
      else
        write(MsgOut,'(A)') '  contact_check   =              no'
      end if
      if (ene_info%structure_check /= StructureCheckNone) then
        write(MsgOut,'(A20,A6)')                             &
              '  structure_check = ',                        &
              trim(StructureCheckTypes(ene_info%structure_check))
      end if

      write(MsgOut,'(A20,F15.3)') '  efield_x        = ', &
             ene_info%efield_x
      write(MsgOut,'(A20,F15.3)') '  efield_y        = ', &
             ene_info%efield_y
      write(MsgOut,'(A20,F15.3)') '  efield_z        = ', &
             ene_info%efield_z
      if (ene_info%efield_normal) then
        write(MsgOut,'(A)') '  efield_normal   =             yes'
      else
        write(MsgOut,'(A)') '  efield_normal   =              no'
      end if
      if (ene_info%efield_virial) then
        write(MsgOut,'(A)') '  efield_virial   =             yes'
      else
        write(MsgOut,'(A)') '  efield_virial   =              no'
      end if

      ! if vacuum
      if (ene_info%vacuum) then
        write(MsgOut,'(A)') '  vacuum          =             yes'
      else
        write(MsgOut,'(A)') '  vacuum          =              no'
      end if

      write(MsgOut,'(A)') ' '

    end if
    ene_info%minimum_contact=ene_info%minimum_contact*ene_info%minimum_contact

    if (ene_info%table_order == 1) then
      cutoff = ene_info%cutoffdist
      cutoff2 = cutoff*cutoff
      cutoff2_water  = (cutoff+4.5_wp) * (cutoff+4.5_wp)
      cutoff_int2    = int(cutoff2_water*ene_info%table_density)
      mind=cutoff2*ene_info%table_density/real(cutoff_int2,wp)+0.001_wp
      ene_info%err_minimum_contact=mind
    end if

    ene_info%err_minimum_contact=ene_info%err_minimum_contact* &
                                 ene_info%err_minimum_contact

    return

  end subroutine read_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy
  !> @brief        compute potential energy
  !! @authors      JJ, HO
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    merge_force   : flag for merge force
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[inout] atmcls_pbc    : atom class number
  !! @param[input] coord_pbc     : coordinates
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_long    : forces of target systems in long range
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : forces
  !! @param[inout] virial_cell   : virial term of target systems in cell
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_long   : virial term of target systems in long range
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy(domain, enefunc, pairlist, boundary, coord,  &
                            npt, reduce, nonb_ene_orig, merge_force,     &
                            nonb_limiter, energy, atmcls_pbc,            &
                            coord_pbc, force, force_long, force_omp,     &
                            force_pbc, virial_cell, virial, virial_long, &
                            virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene_orig
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter 
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    real(wip)                :: volume
    logical                  :: nonb_ene


    call timer(TimerEnergy, TimerOn)

    if (enefunc%gamd_use) then
      nonb_ene = .true.
    else
      nonb_ene = nonb_ene_orig
    end if

    if (domain%fep_use .and. enefunc%gamd_use) &
      call error_msg('Compute_Energy> GaMD is not available in FEP')

    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      if (domain%fep_use) then

        call compute_energy_charmm_fep( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, reduce, nonb_ene, merge_force,          &
                              nonb_limiter, energy, atmcls_pbc,            &
                              coord_pbc, force, force_long, force_omp,     &
                              force_pbc, virial_cell, virial, virial_long, &
                              virial_ext)

      else

        call compute_energy_charmm( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, reduce, nonb_ene, merge_force,          &
                              nonb_limiter, energy, atmcls_pbc,            &
                              coord_pbc, force, force_long, force_omp,     &
                              force_pbc, virial_cell, virial, virial_long, &
                              virial_ext)

      end if



    case (ForcefieldAMBER)

      if (domain%fep_use) then

        call compute_energy_amber_fep( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, reduce, nonb_ene, merge_force,          &
                              nonb_limiter, energy, atmcls_pbc,            &
                              coord_pbc, force, force_long, force_omp,     &
                              force_pbc, virial_cell, virial, virial_long, &
                              virial_ext)

      else

        call compute_energy_amber( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, reduce, nonb_ene, merge_force,          &
                              nonb_limiter, energy, atmcls_pbc,            &
                              coord_pbc, force, force_long, force_omp,     &
                              force_pbc, virial_cell, virial, virial_long, &
                              virial_ext)

      end if

    case (ForcefieldGROAMBER)

      if (domain%fep_use) &
        call error_msg('Compute_Energy> groamber is not available in FEP')

      call compute_energy_gro_amber( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, reduce, nonb_ene, merge_force,          &
                              nonb_limiter, energy, atmcls_pbc,            &
                              coord_pbc, force, force_long, force_omp,     &
                              force_pbc, virial_cell, virial, virial_long, &
                              virial_ext)


    case (ForcefieldGROMARTINI)

      if (domain%fep_use) &
        call error_msg('Compute_Energy> gro_martini is not available in FEP')

      call compute_energy_gro_martini( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, nonb_limiter,        &
                              energy, atmcls_pbc, coord_pbc,              &
                              force, force_long, force_omp, force_pbc,    &
                              virial_cell, virial, virial_ext)


    end select

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume =  boundary%box_size_x_ref * &
                boundary%box_size_y_ref * &
                boundary%box_size_z_ref
      if (domain%fep_use) then
        enefunc%dispersion_energy = &
          enefunc%dispersion_energy_CC + &
          enefunc%dispersion_energy_AC * enefunc%lambljA + &
          enefunc%dispersion_energy_BC * enefunc%lambljB + &
          enefunc%dispersion_energy_AA * enefunc%lambljA * enefunc%lambljA + &
          enefunc%dispersion_energy_BB * enefunc%lambljB * enefunc%lambljB + &
          enefunc%dispersion_energy_AB * enefunc%lambljA * enefunc%lambljB
      end if
      energy%disp_corr_energy = enefunc%dispersion_energy / volume
      if (enefunc%dispersion_corr == Disp_corr_EPress) then
        if (domain%fep_use) then
          enefunc%dispersion_virial = &
            enefunc%dispersion_virial_CC + &
            enefunc%dispersion_virial_AC * enefunc%lambljA + &
            enefunc%dispersion_virial_BC * enefunc%lambljB + &
            enefunc%dispersion_virial_AA * enefunc%lambljA * enefunc%lambljA + &
            enefunc%dispersion_virial_BB * enefunc%lambljB * enefunc%lambljB + &
            enefunc%dispersion_virial_AB * enefunc%lambljA * enefunc%lambljB
        end if
        energy%disp_corr_virial = enefunc%dispersion_virial / volume
        if (replica_main_rank .or. main_rank) then
          virial(1,1) = virial(1,1) + real(energy%disp_corr_virial,dp)
          virial(2,2) = virial(2,2) + real(energy%disp_corr_virial,dp)
          virial(3,3) = virial(3,3) + real(energy%disp_corr_virial,dp)
        end if
      end if

    end if

    if (enefunc%rpath_sum_mf_flag) then
      call compute_stats(enefunc)
    end if

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_short
  !> @brief        compute potential energy with short range interaction
  !! @authors      JJ, HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_short(domain, enefunc, pairlist, boundary, coord, &
                                  npt, energy, atmcls_pbc,                    &
                                  coord_pbc, force, force_omp, force_pbc,     &
                                  virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)


    call timer(TimerEnergy, TimerOn)


    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      if (domain%fep_use) then
        call compute_energy_charmm_short_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, energy, atmcls_pbc, coord_pbc,         &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)
      else
        call compute_energy_charmm_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, energy, atmcls_pbc, coord_pbc,         &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)
      end if

    case (ForcefieldAMBER)

      if (domain%fep_use) then
        call compute_energy_amber_short_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, energy, atmcls_pbc, coord_pbc,         &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)
      else
        call compute_energy_amber_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, energy, atmcls_pbc, coord_pbc,         &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)
      end if

    case (ForcefieldGROAMBER)

      if (domain%fep_use) &
        call error_msg('GROAMBER field is not available in FEP')

      call compute_energy_gro_amber_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, energy, atmcls_pbc, coord_pbc,         &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    end select


    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_long
  !> @brief        compute potential energy
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary information
  !! @param[in]    npt      : flag for NPT or not
  !! @param[inout] energy   : energy information
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_long(domain, enefunc, boundary, &
                                 npt, energy, force, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)


    call timer(TimerEnergy, TimerOn)

    call compute_energy_general_long( &
                              domain, enefunc, boundary, &
                              npt, energy, force, virial)

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_long

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy
  !> @brief        output energy
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
  !  Subroutine    compute_energy_charmm
  !> @brief        compute potential energy with charmm force field
  !! @authors      JJ
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    merge_force  : flag for merge force
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : forces
  !! @param[inout] virial_cell  : virial term of target systems in cell
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_long  : virial term of target systems in long range
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm(domain, enefunc, pairlist, boundary, coord, &
                                   npt, reduce, nonb_ene, merge_force,         &
                                   nonb_limiter, energy, atmcls_pbc,           &
                                   coord_pbc, force, force_long, force_omp,    &
                                   force_pbc, virial_cell, virial,             &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: trans(1:3)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, k, ix, ic, jc, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    !$omp parallel do private(i,ix,id)
    do i = 1, ncell
      do ix = 1, natom
        force     (1,ix,i) = 0.0_wip
        force     (2,ix,i) = 0.0_wip
        force     (3,ix,i) = 0.0_wip
        force_long(1,ix,i) = 0.0_wip
        force_long(2,ix,i) = 0.0_wip
        force_long(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end parallel do

    virial     (1:3,1:3) = 0.0_dp
    virial_long(1:3,1:3) = 0.0_dp
    virial_ext (1:3,1:3) = 0.0_dp

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do
  
      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      k = domain%num_atom_domain
      !$omp parallel do
      do i = 1, 3*k
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)


      else

        call compute_energy_nonbond_cutoff( &
                              domain, enefunc, pairlist,             &
                              nonb_ene, force_pbc, virial_omp,       &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Charmm> Unknown boundary condition')

    end select


    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

      ! dihedral and cmap energies
      !
      if (enefunc%gamd_use) then

        if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

          call compute_energy_dihed_gamd( &
                              domain, enefunc, npt, nonb_ene, coord, &
                              force_omp, virial_omp,                 &
                              edihed_omp, ecmap_omp, eimprop_omp)

        else

          if (enefunc%local_restraint) then
            call compute_energy_dihed_localres( &
                              domain, enefunc, coord,                &
                              force_omp, virial_omp, edihed_omp)
          else
            call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
          end if
      
          call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

        end if

      else

        if (enefunc%local_restraint) then
          call compute_energy_dihed_localres( &
                              domain, enefunc, coord,                &
                              force_omp, virial_omp, edihed_omp)
        else
          call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
        end if
    
        call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

      end if

      ! improper energy
      !
      call compute_energy_improp( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

      ! 1-4 interaction 
      !
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_table_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      ! external electric field
      !
      if (enefunc%use_efield) &
        call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

      ! restraint energy
      !
      if (enefunc%restraint) then

        if (enefunc%gamd_use) then

          if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

            call compute_energy_restraints_gamd( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp, eposi_omp, energy)

          else

            call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

          end if

        else

          call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

        end if

      end if

    end if

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, npt, nonb_ene, coord_pbc,    &
                              force_omp, force_pbc, virial_cell,     &
                              virial_omp, elec_omp, evdw_omp)
    end if

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt) .and. &
        domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (domain%nonbond_kernel /= NBK_GPU) then

        do i = 1, ncell
          do k = 1, 3
            do ix = 1, domain%num_atom(i)
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                   coord_pbc(ix,k,i)*force_pbc(ix,k,i,id+1)
            end do
          end do
        end do
        do i = id+1, maxcell, nthread
          ic = domain%cell_pairlist1(1,i)
          jc = domain%cell_pairlist1(2,i)
          if (domain%virial_check(jc,ic)==1) then
            trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                   - trans(k)*virial_cell(k,i)
            end do
          end if
        end do

      else

        k = domain%num_atom_domain
        do i = id+1, domain%num_atom_domain, nthread
          virial_omp(1,1,id+1) = virial_omp(1,1,id+1) &
                               + coord_pbc(i,1,1)*force_pbc(i,1,1,1)
          virial_omp(2,2,id+1) = virial_omp(2,2,id+1) &
                               + coord_pbc(k+i,1,1)*force_pbc(k+i,1,1,1)
          virial_omp(3,3,id+1) = virial_omp(3,3,id+1) &
                               + coord_pbc(2*k+i,1,1)*force_pbc(2*k+i,1,1,1)
        end do
      end if
      !$omp end parallel

    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)

          do ic = 2, nthread
            force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          if (merge_force) force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
          if (merge_force) force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
          if (merge_force) force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, ic, force_tmp) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel 

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, ic, force_tmp)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    do id = 1, nthread

      virial    (1,1) = virial    (1,1) + virial_omp    (1,1,id)
      virial    (2,2) = virial    (2,2) + virial_omp    (2,2,id)
      virial    (3,3) = virial    (3,3) + virial_omp    (3,3,id)
      virial_ext(1,1) = virial_ext(1,1) + virial_ext_omp(1,1,id)
      virial_ext(2,2) = virial_ext(2,2) + virial_ext_omp(2,2,id)
      virial_ext(3,3) = virial_ext(3,3) + virial_ext_omp(3,3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)
      energy%electric_field     = energy%electric_field     + efield_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond            &
                 + energy%angle           &
                 + energy%urey_bradley    &
                 + energy%dihedral        &
                 + energy%cmap            &
                 + energy%improper        &
                 + energy%electrostatic   &
                 + energy%van_der_waals   &
                 + energy%electric_field  &
                 + energy%restraint_position

    ! GaMD boost and statistics
    !
    if (enefunc%gamd_use) then
      call boost_gamd(domain, enefunc, energy, force, virial)
    end if

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber
  !> @brief        compute potential energy with AMBER99 force field
  !! @authors      JJ
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    merge_force  : flag for merge force
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : pbc forces
  !! @param[inout] virial_cell  : virial correction due to pbc
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_long  : virial term of target systems in long range
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber(domain, enefunc, pairlist, boundary, coord,  &
                                  npt, reduce, nonb_ene, merge_force,          &
                                  nonb_limiter, energy, atmcls_pbc,            &
                                  coord_pbc, force, force_long,                &
                                  force_omp, force_pbc, virial_cell, virial,   &
                                  virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: trans(1:3)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    !$omp parallel private(i,ix,id)
    !$omp do
    do i = 1, ncell
      do ix = 1, natom
        force     (1,ix,i) = 0.0_wip
        force     (2,ix,i) = 0.0_wip
        force     (3,ix,i) = 0.0_wip
        force_long(1,ix,i) = 0.0_wip
        force_long(2,ix,i) = 0.0_wip
        force_long(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end parallel

    virial     (1:3,1:3) = 0.0_dp
    virial_long(1:3,1:3) = 0.0_dp
    virial_ext (1:3,1:3) = 0.0_dp

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do 
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do

      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp 

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
      else

        call compute_energy_nonbond_cutoff( &
                              domain, enefunc, pairlist,             &
                              nonb_ene, force_pbc, virial_omp,       &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Amber> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      if (enefunc%gamd_use) then
        if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

          call compute_energy_dihed_gamd( &
                              domain, enefunc, npt, nonb_ene,        &
                              coord, force_omp, virial_omp,          &
                              edihed_omp, ecmap_omp, eimprop_omp)

        else

          if (enefunc%local_restraint) then
            call compute_energy_dihed_localres( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
          else
            call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
          end if

          call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

        end if

      else

        if (enefunc%local_restraint) then
          call compute_energy_dihed_localres( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
        else
          call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
        end if

        call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

      end if

      ! improper energy
      !
      call compute_energy_improp_cos( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

      ! 1-4 interaction with linear table
      !
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_notable( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      if (enefunc%use_efield) &
        call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

      ! restraint energy
      !
      if (enefunc%restraint) then

        if (enefunc%gamd_use) then

          if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

            call compute_energy_restraints_gamd( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy)

          else

            call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

          end if

        else

          call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

        end if
      end if
    end if

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, npt, nonb_ene, coord_pbc,    &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt) .and. &
        domain%nonbond_kernel /= NBK_Fugaku .and.       &
        domain%nonbond_kernel /= NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (domain%nonbond_kernel /= NBK_GPU) then

        do i = 1, ncell
          do k = 1, 3
            do ix = 1, domain%num_atom(i)
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                   coord_pbc(ix,k,i)*force_pbc(ix,k,i,id+1)
            end do
          end do
        end do
        do i = id+1, maxcell, nthread
          ic = domain%cell_pairlist1(1,i)
          jc = domain%cell_pairlist1(2,i)
          if (domain%virial_check(jc,ic)==1) then
            trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                   - trans(k)*virial_cell(k,i)
            end do
          end if
        end do

      else

        k = domain%num_atom_domain
        do i = id+1, domain%num_atom_domain, nthread
          virial_omp(1,1,id+1) = virial_omp(1,1,id+1) &
                               + coord_pbc(i,1,1)*force_pbc(i,1,1,1)
          virial_omp(2,2,id+1) = virial_omp(2,2,id+1) &
                               + coord_pbc(k+i,1,1)*force_pbc(k+i,1,1,1)
          virial_omp(3,3,id+1) = virial_omp(3,3,id+1) &
                               + coord_pbc(2*k+i,1,1)*force_pbc(2*k+i,1,1,1)
        end do

      end if
      !$omp end parallel

    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
          if (domain%nonbond_kernel /= NBK_GPU) then
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do
          else
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do
          end if
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          if (merge_force) force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
          if (merge_force) force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
          if (merge_force) force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if
!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    do id = 1, nthread

      virial    (1,1) = virial    (1,1) + virial_omp    (1,1,id)
      virial    (2,2) = virial    (2,2) + virial_omp    (2,2,id)
      virial    (3,3) = virial    (3,3) + virial_omp    (3,3,id)
      virial_ext(1,1) = virial_ext(1,1) + virial_ext_omp(1,1,id)
      virial_ext(2,2) = virial_ext(2,2) + virial_ext_omp(2,2,id)
      virial_ext(3,3) = virial_ext(3,3) + virial_ext_omp(3,3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do


    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%cmap          &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    ! GaMD boost and statistics
    !
    if (enefunc%gamd_use) then 
      call boost_gamd(domain, enefunc, energy, force, virial)
    end if

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !! @authors      JJ
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    merge_force  : flag for merge force
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : pbc forces
  !! @param[inout] virial_cell  : virial correction due to pbc
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_long  : virial term of target systems in long range
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber(domain, enefunc, pairlist, boundary,     &
                                      coord, npt, reduce, nonb_ene,            &
                                      merge_force, nonb_limiter,               &
                                      energy, atmcls_pbc, coord_pbc,           &
                                      force, force_long, force_omp, force_pbc, &
                                      virial_cell, virial, virial_long,        &
                                      virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: efield_omp  (nthread)
    real(wip)                :: trans(1:3)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    !$omp parallel private(i,ix,id)
    !$omp do
    do i = 1, ncell
      do ix = 1, natom
        force     (1,ix,i) = 0.0_wip
        force     (2,ix,i) = 0.0_wip
        force     (3,ix,i) = 0.0_wip
        force_long(1,ix,i) = 0.0_wip
        force_long(2,ix,i) = 0.0_wip
        force_long(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end parallel

    virial     (1:3,1:3) = 0.0_dp
    virial_long(1:3,1:3) = 0.0_dp
    virial_ext (1:3,1:3) = 0.0_dp

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do

      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    erbdihed_omp  (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp
    efield_omp    (1:nthread) = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)

      else

        call compute_energy_nonbond_cutoff( &
                              domain, enefunc, pairlist,             &
                              nonb_ene, force_pbc, virial_omp,       &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

      ! Ryckaert-Bellemans dihedral energy
      !
      call compute_energy_rb_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, erbdihed_omp)

      if (enefunc%pme_use) then

        if (enefunc%vdw_no_switch) then
          call compute_energy_nonbond14_notable( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)
        else
          call compute_energy_nonbond14_table_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)
        end if

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      ! external electric field
      !
      if (enefunc%use_efield) &
        call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)
    end if

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, npt, nonb_ene, coord_pbc,    &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt) .and. &
        domain%nonbond_kernel /= NBK_Fugaku .and.       &
        domain%nonbond_kernel /= NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (domain%nonbond_kernel /= NBK_GPU) then

        do i = 1, ncell
          do k = 1, 3
            do ix = 1, domain%num_atom(i)
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                   coord_pbc(ix,k,i)*force_pbc(ix,k,i,id+1)
            end do
          end do
        end do
        do i = id+1, maxcell, nthread
          ic = domain%cell_pairlist1(1,i)
          jc = domain%cell_pairlist1(2,i)
          if (domain%virial_check(jc,ic)==1) then
            trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                   - trans(k)*virial_cell(k,i)
            end do
          end if
        end do

      else

        k = domain%num_atom_domain
        do i = id+1, domain%num_atom_domain, nthread
          virial_omp(1,1,id+1) = virial_omp(1,1,id+1) &
                               + coord_pbc(i,1,1)*force_pbc(i,1,1,1)
          virial_omp(2,2,id+1) = virial_omp(2,2,id+1) &
                               + coord_pbc(k+i,1,1)*force_pbc(k+i,1,1,1)
          virial_omp(3,3,id+1) = virial_omp(3,3,id+1) &
                               + coord_pbc(2*k+i,1,1)*force_pbc(2*k+i,1,1,1)
        end do

      end if
      !$omp end parallel

    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)

          if (domain%nonbond_kernel /= NBK_GPU) then

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do

          else

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do

          end if

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          if (merge_force) force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
          if (merge_force) force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
          if (merge_force) force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    do id = 1, nthread

      virial    (1,1) = virial    (1,1) + virial_omp    (1,1,id)
      virial    (2,2) = virial    (2,2) + virial_omp    (2,2,id)
      virial    (3,3) = virial    (3,3) + virial_omp    (3,3,id)
      virial_ext(1,1) = virial_ext(1,1) + virial_ext_omp(1,1,id)
      virial_ext(2,2) = virial_ext(2,2) + virial_ext_omp(2,2,id)
      virial_ext(3,3) = virial_ext(3,3) + virial_ext_omp(3,3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%electric_field     = energy%electric_field     + efield_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond            &
                 + energy%angle           &
                 + energy%urey_bradley    &
                 + energy%dihedral        &
                 + energy%electrostatic   &
                 + energy%van_der_waals   &
                 + energy%electric_field  &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_martini
  !> @brief        compute potential energy with GROMACS-MARTINI force field
  !! @authors      JJ
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : pbc forces
  !! @param[inout] virial_cell  : virial correction due to pbc
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_martini(domain, enefunc, pairlist, boundary,   &
                                        coord, npt, reduce, nonb_ene,          &
                                        nonb_limiter, energy, atmcls_pbc,      &
                                        coord_pbc, force, force_long,          &
                                        force_omp, force_pbc,                  &
                                        virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    integer                  :: ncell, natom, id, i, ix, k


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell)   = 0.0_wip
    force_long(1:3,1:natom,1:ncell)   = 0.0_wip
    virial    (1:3,1:3)               = 0.0_dp 
    virial_ext(1:3,1:3)               = 0.0_dp 

    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:natom,1:3,1:ncell,1:nthread) = 0.0_wp
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    ebond_omp     (1:nthread)         = 0.0_dp 
    eangle_omp    (1:nthread)         = 0.0_dp 
    edihed_omp    (1:nthread)         = 0.0_dp 
    elec_omp      (1:nthread)         = 0.0_dp 
    evdw_omp      (1:nthread)         = 0.0_dp 
    eposi_omp     (1:nthread)         = 0.0_dp 

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_g96( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp)

      ! dihedral energy
      !
      call compute_energy_dihed(domain, enefunc, coord, force_omp,   &
                              virial_omp, edihed_omp)

      if (enefunc%pme_use) then

        call compute_energy_nonbond14_table_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)
    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
      else

        call compute_energy_nonbond_cutoff( &
                              domain, enefunc, pairlist, nonb_ene,   &
                              force_pbc, virial_omp, elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Gro_Martini> Unknown boundary condition')

    end select


    ! gather values
    !
    do id = 1, nthread

      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          force(1:3,ix,i) = force(1:3,ix,i) + force_omp(1:3,ix,i,id)
          force(1:3,ix,i) = force(1:3,ix,i) + force_pbc(ix,1:3,i,id)
        end do
      end do

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_gro_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_short
  !> @brief        compute potential energy with charmm force field
  !! @authors      JJ

  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_short(domain, enefunc, pairlist, boundary, &
                                         coord, npt, energy,                  &
                                         atmcls_pbc, coord_pbc,               &
                                         force, force_omp, force_pbc,         &
                                         virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, k, ki, ic, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    virial    (1,1)             = 0.0_dp
    virial    (2,2)             = 0.0_dp
    virial    (3,3)             = 0.0_dp
    virial_ext(1,1)             = 0.0_dp
    virial_ext(2,2)             = 0.0_dp
    virial_ext(3,3)             = 0.0_dp

    !$omp parallel private(i,ix,id)

    !$omp do 
    do i = 1, ncell
      do ix = 1, natom
        force(1,ix,i) = 0.0_wip
        force(2,ix,i) = 0.0_wip
        force(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end do

    if (domain%nonbond_kernel /= NBK_GPU) then
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        !$omp do
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, domain%num_atom(i)
              force_omp(1, ix, i, id) = 0.0_wp
              force_omp(2, ix, i, id) = 0.0_wp
              force_omp(3, ix, i, id) = 0.0_wp
              force_pbc(ix, 1, i, id) = 0.0_wp
              force_pbc(ix, 2, i, id) = 0.0_wp
              force_pbc(ix, 3, i, id) = 0.0_wp 
            end do
          end do
        end do
        !$omp end do
      else
        !$omp do
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
      end if
    else
      !$omp do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end do

      k = domain%num_atom_domain

      !$omp do
      do i = 1, 3*k
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end do
    end if

    !$omp do
    do id = 1, nthread
      virial_omp    (1,1,id) = 0.0_dp
      virial_omp    (2,2,id) = 0.0_dp
      virial_omp    (3,3,id) = 0.0_dp
      virial_ext_omp(1,1,id) = 0.0_dp
      virial_ext_omp(2,2,id) = 0.0_dp
      virial_ext_omp(3,3,id) = 0.0_dp
      ebond_omp     (id) = 0.0_dp
      eangle_omp    (id) = 0.0_dp
      eurey_omp     (id) = 0.0_dp
      edihed_omp    (id) = 0.0_dp
      eimprop_omp   (id) = 0.0_dp
      ecmap_omp     (id) = 0.0_dp
      elec_omp      (id) = 0.0_dp
      evdw_omp      (id) = 0.0_dp
      eposi_omp     (id) = 0.0_dp
      efield_omp    (id) = 0.0_dp
    end do

    !$omp end parallel

    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel)       &
      virial_cell   (1:3,1:maxcell) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc,                 &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp)

      else

        call error_msg( &
        'Compute_Energy_Charmm_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Charmm_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then

      call compute_energy_dihed_localres( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

    else

      call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

    end if

    ! improper energy
    !
    call compute_energy_improp( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

    ! cmap energy
    !
    call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

    if (enefunc%use_efield) &
      call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             & 
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, .false., .false.,            &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (domain%nonbond_kernel /= NBK_GPU) then
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do
          else
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do
          end if
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel
   
    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, k, ix, ki, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)

          ki = k + ix
          force_tmp(1) = 0.0_wip
          force_tmp(2) = 0.0_wip
          force_tmp(3) = 0.0_wip
!ocl nosimd
          do ic = 1, nthread
  
            force_tmp(1) = force_tmp(1) &
                         + force_omp(1,ix,i,ic)+force_pbc(1,ki,1,ic)
            force_tmp(2) = force_tmp(2) &
                         + force_omp(2,ix,i,ic)+force_pbc(2,ki,1,ic)
            force_tmp(3) = force_tmp(3) &
                         + force_omp(3,ix,i,ic)+force_pbc(3,ki,1,ic)
          end do
          force(1,ix,i) = force_tmp(1)
          force(2,ix,i) = force_tmp(2)
          force(3,ix,i) = force_tmp(3)
        end do
      end do
      !$omp end parallel 

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, k, ix, ki, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)

          ki = k + ix
          force_tmp(1) = 0.0_wip
          force_tmp(2) = 0.0_wip
          force_tmp(3) = 0.0_wip
!ocl nosimd
          do ic = 1, nthread

            force_tmp(1) = force_tmp(1) &
                         + force_omp(1,ix,i,ic)+force_pbc(ki,1,1,ic)
            force_tmp(2) = force_tmp(2) &
                         + force_omp(2,ix,i,ic)+force_pbc(ki,2,1,ic)
            force_tmp(3) = force_tmp(3) &
                         + force_omp(3,ix,i,ic)+force_pbc(ki,3,1,ic)
          end do
          force(1,ix,i) = force_tmp(1)
          force(2,ix,i) = force_tmp(2)
          force(3,ix,i) = force_tmp(3)
        end do
      end do
      !$omp end parallel

    end if

    return

  end subroutine compute_energy_charmm_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_short
  !> @brief        compute potential energy with AMBER99 force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber_short(domain, enefunc, pairlist, boundary,   &
                                        coord, npt, energy, atmcls_pbc,        &
                                        coord_pbc, force, force_omp, force_pbc,&
                                        virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_wip
    virial    (1:3,1:3)             = 0.0_dp 
    virial_ext(1:3,1:3)             = 0.0_dp 

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do
          
      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp 

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short( &
                              domain, enefunc, pairlist, npt,        & 
                              atmcls_pbc, coord_pbc,                 &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp)

      else

        call error_msg( &
        'Compute_Energy_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
    else
      call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp_cos( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

    ! cmap energy
    !
    call compute_energy_cmap( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_notable( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

    if (enefunc%use_efield) &
      call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, .false., .false.,            &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)

          if (domain%nonbond_kernel /= NBK_GPU) then

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do

          else

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do

          end if

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    return

  end subroutine compute_energy_amber_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber_short
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber_short(domain, enefunc, pairlist,         &
                                            boundary, coord, npt, energy,      &
                                            atmcls_pbc, coord_pbc,             &
                                            force, force_omp, force_pbc,       &
                                            virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: efield_omp  (nthread)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_wip
    virial    (1:3,1:3)             = 0.0_dp 
    virial_ext(1:3,1:3)             = 0.0_dp 

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do

      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    erbdihed_omp  (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp 

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc,                 &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp)

      else

        call error_msg( &
        'Compute_Energy_Gro_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle( & 
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    call compute_energy_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

    ! Ryckaert-Bellemans dihedral energy
    !
    call compute_energy_rb_dihed( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, erbdihed_omp)

    ! 1-4 interaction with linear table
    !
    if (enefunc%vdw_no_switch) then
      call compute_energy_nonbond14_notable( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)
    else
      call compute_energy_nonbond14_table_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)
    end if

    call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

    if (enefunc%use_efield) &
      call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints( &
                              .true., .true., domain, boundary,      &
                              enefunc, coord, force_omp,             &
                              virial_omp, virial_ext_omp,            &
                              eposi_omp, energy%restraint_rmsd,      &
                              energy%rmsd, energy%restraint_distance,&
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, .false., .false., coord_pbc, &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)

          if (domain%nonbond_kernel /= NBK_GPU) then

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do

          else

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do

          end if

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    return

  end subroutine compute_energy_gro_amber_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_general_long
  !> @brief        compute long range interaction
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary information
  !! @param[in]    npt      : flag for NPT or not
  !! @param[inout] energy   : energy information
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_general_long(domain, enefunc, boundary, &
                                         npt, energy, force, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    real(wip),               intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: elec_long
    real(dp)                 :: elec_omp(nthread)
    real(dp)                 :: evdw_omp(nthread)
    integer                  :: ncell, natom, id


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    force(1:3,1:natom,1:ncell) = 0.0_wip
    virial(1:3,1:3)            = 0.0_dp 
    elec_long                  = 0.0_dp 

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp 
    elec_omp  (1:nthread) = 0.0_dp 
    evdw_omp  (1:nthread) = 0.0_dp 

    if (enefunc%pme_use) &
      call compute_energy_nonbond_pme_long(domain, enefunc, boundary, &
                                           npt, force, virial_omp,    &
                                           elec_omp, evdw_omp)

    do id = 1, nthread

      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,id)

      elec_long = elec_long + elec_omp(id)

    end do

    energy%electrostatic = energy%electrostatic + elec_long

    ! total energy
    energy%total = energy%total + elec_long

    return

  end subroutine compute_energy_general_long

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_wait
  !> @brief        get the information from GPU to CPU
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_wait(pairlist,      &
                                             npt, nonb_ene, coord_pbc,       &
                                             force, force_pbc, virial_cell,  &
                                             virial, eelec, evdw)

    ! formal arguments
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)


#ifdef USE_GPU
    ! local
    real(dp) :: ene_virial(1:5)

    !
    ! GPU is used and overlapping is enabled
    !
    if (nonb_ene) then
      call gpu_wait_compute_energy_nonbond_table_linear_univ( &
                              coord_pbc, force_pbc(:,:,:,1),         &
                              ene_virial)
      eelec(1) = eelec(1) + ene_virial(1)
      evdw (1) = evdw(1)  + ene_virial(2)
      virial(1,1,1) = virial(1,1,1) + ene_virial(3)
      virial(2,2,1) = virial(2,2,1) + ene_virial(4)
      virial(3,3,1) = virial(3,3,1) + ene_virial(5)
    else
      call gpu_wait_compute_force_nonbond_table_linear_univ( &
                              coord_pbc, force_pbc(:,:,:,1),         &
                              ene_virial, npt)
      virial(1,1,1) = virial(1,1,1) + ene_virial(1)
      virial(2,2,1) = virial(2,2,1) + ene_virial(2)
      virial(3,3,1) = virial(3,3,1) + ene_virial(3)
    end if
    call timer(TimerPmeReal, TimerOff)
    call timer(TimerNonBond, TimerOff)
#endif

    return

  end subroutine compute_energy_nonbond_pme_wait

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_genesis
  !> @brief        output energy in GENESIS style
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

    if (enefunc%num_bond_all > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = energy%bond
      ifm = ifm+1
    end if

    if (enefunc%num_angl_all > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = energy%urey_bradley
        ifm = ifm+1
      end if

    end if

    if (enefunc%num_dihe_all > 0 .or. enefunc%num_rb_dihe_all > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = energy%dihedral
      ifm = ifm+1
    end if

    if (enefunc%num_impr_all > 0 ) then
      write(category(ifm),frmt) 'IMPROPER'
      values(ifm) = energy%improper
      ifm = ifm+1
    end if

    if (enefunc%forcefield == ForcefieldCHARMM .or. &
        enefunc%forcefield == ForcefieldAMBER) then
      if (enefunc%num_cmap_all > 0 ) then
        write(category(ifm),frmt) 'CMAP'
        values(ifm) = energy%cmap
        ifm = ifm+1
      end if
    end if

    write(category(ifm),frmt) 'VDWAALS'
    values(ifm) = energy%van_der_waals
    ifm = ifm+1

    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(category(ifm),frmt) 'DISP-CORR_ENE'
      values(ifm) = energy%disp_corr_energy
      ifm = ifm+1
    end if

    write(category(ifm),frmt) 'ELECT'
    values(ifm) = energy%electrostatic
    ifm = ifm+1

    if (enefunc%use_efield) then
      write(category(ifm),frmt) 'ELECTRIC FIELD'
      values(ifm) = energy%electric_field
      ifm = ifm+1
    end if

    if (enefunc%num_restraintfuncs > 0) then
      if (enefunc%restraint_rmsd) then
        write(category(ifm),frmt) 'RMSD'
        values(ifm) = energy%rmsd
        ifm = ifm+1
      end if

      ene_restraint =   energy%restraint_distance &
                      + energy%restraint_position &
                      + energy%restraint_rmsd
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

    end if

    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_pot) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = energy%total_gamd
        ifm = ifm+1
      else if (enefunc%gamd%boost_dih) then
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = energy%dihedral_gamd
        ifm = ifm+1
      else if (enefunc%gamd%boost_dual) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = energy%total_gamd
        ifm = ifm+1
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = energy%dihedral_gamd
        ifm = ifm+1
      end if
    end if

    if (etitle) then

      write(MsgOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(MsgOut,frmt) category(i)
        else
          write(MsgOut,frmt_cont) category(i)
        end if
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
      end if
    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine output_energy_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_charmm
  !> @brief        output energy in CHARMM style
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

      if (enefunc%num_cmap_all > 0) then
        write(MsgOut,'(A79)') 'DYNA CROSS:           CMAPs                                                    '
      end if
      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(MsgOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      end if

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

    if (enefunc%num_cmap_all > 0) then
      write(MsgOut,'(A14, F13.5)')   'DYNA CROSS>   ', cmaps
    end if
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
  !> @brief        output energy in NAMD style
  !! @authors      YS, CK
  !! @param[in]    step   : step
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
  !> @brief        output energy in GROMACS style
  !! @authors      NT
  !! @param[in]    step   : step
  !! @param[in]    energy : energy information
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


    time        = 0.0_wip
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
    energy_     = 0.0_wip
    totke       = 0.0_wip
    totener     = 0.0_wip
    temperature = 0.0_wip
    pressi      = 0.0_wip
    presse      = 0.0_wip


    write(MsgOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(MsgOut,'(I15,2F15.5)') &
         step, time, 0.0_wip
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
  !  Subroutine    reduce_ene
  !> @brief        reduce energy and virial
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
    real(dp)                 :: before_reduce(20), after_reduce(20)
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
    before_reduce(12) = energy%urey_bradley
    before_reduce(13) = energy%dihedral
    before_reduce(14) = energy%improper
    before_reduce(15) = energy%cmap
    before_reduce(16) = energy%electrostatic
    before_reduce(17) = energy%van_der_waals
    before_reduce(18) = energy%electric_field
    before_reduce(19) = energy%restraint_position
    before_reduce(20) = energy%total

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 20, mpi_real8,  &
                       mpi_sum, mpi_comm_country, ierror)
#else
    after_reduce(1:20) = before_reduce(1:20)
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
    energy%urey_bradley       = after_reduce(12)
    energy%dihedral           = after_reduce(13)
    energy%improper           = after_reduce(14)
    energy%cmap               = after_reduce(15)
    energy%electrostatic      = after_reduce(16)
    energy%van_der_waals      = after_reduce(17)
    energy%electric_field     = after_reduce(18)
    energy%restraint_position = after_reduce(19)
    energy%total              = after_reduce(20)

#endif

    return

  end subroutine reduce_ene

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_stats
  !> @brief        compute statistical quantities for RPATH
  !! @authors      YM
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_stats(enefunc)

    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, ifunc
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j
    real(dp)                 :: etmp, stmp, dtmp
    real(dp)                 :: d(1:3)
    real(wp),    allocatable :: collection(:)


    if (enefunc%rpath_pos_func > 0) then

      allocate(collection(enefunc%stats_dimension))

      collection(:) = 0.0_wp 

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(enefunc%stats_delta, collection, enefunc%stats_dimension,&
                         mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
#endif

    end if

    if (.not. replica_main_rank) then
      if (allocated(collection)) deallocate(collection)
      return
    end if

    do i = 1, enefunc%stats_dimension

      if (enefunc%rpath_pos_func > 0) then

        ifunc = enefunc%rpath_pos_func
        dtmp  = real(collection(i),dp)
        enefunc%stats_force(i) = enefunc%stats_force(i) + &
        2.0_dp * real(enefunc%restraint_const(1,ifunc),dp) * dtmp

      else

        ifunc = enefunc%rpath_rest_function(i)
        dtmp  = real(enefunc%stats_delta(i),dp)
        enefunc%stats_force(i) = enefunc%stats_force(i) + &
        2.0_dp * real(enefunc%restraint_const(1,ifunc),dp) * dtmp

      end if

    end do

    ifunc = enefunc%rpath_rest_function(1)

    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI .or. &
        enefunc%restraint_kind(ifunc) == RestraintsFuncPC .or. &
        enefunc%restraint_kind(ifunc) == RestraintsFuncPCCOM) then
      do dimno_i = 1, enefunc%stats_dimension
        etmp = 0.0_dp
        do i = 1, enefunc%stats_natom
          d(1:3) = real(enefunc%stats_grad(1:3,i,dimno_i),dp)
          do k = 1, 3
            etmp = etmp + (1.0_dp/real(enefunc%stats_mass(i,dimno_i),dp))*d(k)*d(k)
          end do
        end do
        enefunc%stats_metric(dimno_i, dimno_i) =  &
           enefunc%stats_metric(dimno_i, dimno_i) + etmp
      end do
    else
      do dimno_i = 1, enefunc%stats_dimension
        do dimno_j = 1, enefunc%stats_dimension
          etmp = 0.0_dp
          do i = 1, enefunc%stats_natom
            atom_i = enefunc%stats_atom(i,dimno_i)
            stmp = (1.0_dp / real(enefunc%stats_mass(i,dimno_i),dp))
            d(1:3) = real(enefunc%stats_grad(1:3,i,dimno_i),dp)
            do j = 1, enefunc%stats_natom
              atom_j = enefunc%stats_atom(j,dimno_j)
              if (atom_i == atom_j) then
                do k = 1, 3
!                  enefunc%stats_metric(dimno_i,dimno_j) = &
!                    enefunc%stats_metric(dimno_i,dimno_j) &
!                    + (1.0_wp / enefunc%stats_mass(i,dimno_i)) &
!                    * enefunc%stats_grad(k,i,dimno_i) * enefunc%stats_grad(k,j,dimno_j)
                  etmp = etmp + stmp * d(k) *  &
                     real(enefunc%stats_grad(k,j,dimno_j),dp)
                end do
              end if
            end do
          end do
          enefunc%stats_metric(dimno_i, dimno_j) =  &
             enefunc%stats_metric(dimno_i, dimno_j) + etmp
        end do
      end do
    end if

    if (enefunc%rpath_pos_func > 0) then
      deallocate(collection)
    end if

    return

  end subroutine compute_stats

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_fep
  !> @brief        compute potential energy with charmm force field for FEP
  !! @authors      HO
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    merge_force  : flag for merge force
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : forces
  !! @param[inout] virial_cell  : virial term of target systems in cell
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_long  : virial term of target systems in long range
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_fep(domain, enefunc, pairlist, boundary,    &
                                   coord,                                      &
                                   npt, reduce, nonb_ene, merge_force,         &
                                   nonb_limiter, energy, atmcls_pbc,           &
                                   coord_pbc, force, force_long, force_omp,    &
                                   force_pbc, virial_cell, virial,             &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: trans(1:3)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, k, ix, ic, jc, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    !$omp parallel do private(i,ix,id)
    do i = 1, ncell
      do ix = 1, natom
        force     (1,ix,i) = 0.0_wip
        force     (2,ix,i) = 0.0_wip
        force     (3,ix,i) = 0.0_wip
        force_long(1,ix,i) = 0.0_wip
        force_long(2,ix,i) = 0.0_wip
        force_long(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end parallel do

    virial     (1:3,1:3)  = 0.0_dp 
    virial_long(1:3,1:3)  = 0.0_dp 
    virial_ext (1:3,1:3)  = 0.0_dp 

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do
  
      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      k = domain%num_atom_domain
      !$omp parallel do
      do i = 1, 3*k
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    !$omp parallel do
    do id = 1, nthread
      domain%force_pbc_fep(1:natom,1:3,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp

    if (enefunc%gamd_use) &
      call error_msg('Compute_Energy> GaMD is not available in FEP')

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_fep( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)

      else

        call compute_energy_nonbond_cutoff_fep( &
                              domain, enefunc, pairlist,             &
                              nonb_ene, force_pbc, virial_omp,       &
                              elec_omp, evdw_omp)

      end if

    case default

      call error_msg('Compute_Energy_Charmm> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

      ! dihedral and cmap energies
      !

      if (enefunc%local_restraint) then
        call compute_energy_dihed_localres_fep( &
                            domain, enefunc, coord,                &
                            force_omp, virial_omp, edihed_omp)
      else
        call compute_energy_dihed_fep( &
                            domain, enefunc, coord, force_omp,     &
                            virial_omp, edihed_omp)
      end if
  
      call compute_energy_cmap_fep( &
                            domain, enefunc, coord, force_omp,     &
                            virial_omp, ecmap_omp)

      ! improper energy
      !
      call compute_energy_improp_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

      ! 1-4 interaction 
      !
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_table_linear_fep( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      ! external electric field
      !
      if (enefunc%use_efield) &
        call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

      ! restraint energy
      !
      if (enefunc%restraint) then

        call compute_energy_restraints_fep( &
                            .true., .true., domain, boundary,      &
                            enefunc, coord, force_omp,             &
                            virial_omp, virial_ext_omp,            &
                            eposi_omp, energy%restraint_rmsd,      &
                            energy%rmsd, energy%restraint_distance,&
                            energy%restraint_emfit, energy%emcorr)

      end if

    end if

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, npt, nonb_ene, coord_pbc,    &
                              force_omp, force_pbc, virial_cell,     &
                              virial_omp, elec_omp, evdw_omp)
    end if

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt) .and. &
        domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (domain%nonbond_kernel /= NBK_GPU) then

        do i = 1, ncell
          do k = 1, 3
            do ix = 1, domain%num_atom(i)
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                   coord_pbc(ix,k,i)*force_pbc(ix,k,i,id+1)
            end do
          end do
        end do
        do i = id+1, maxcell, nthread
          ic = domain%cell_pairlist1(1,i)
          jc = domain%cell_pairlist1(2,i)
          if (domain%virial_check(jc,ic)==1) then
            trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                   - trans(k)*virial_cell(k,i)
            end do
          end if
        end do

      else

        k = domain%num_atom_domain
        do i = id+1, domain%num_atom_domain, nthread
          virial_omp(1,1,id+1) = virial_omp(1,1,id+1) &
                               + coord_pbc(i,1,1)*force_pbc(i,1,1,1)
          virial_omp(2,2,id+1) = virial_omp(2,2,id+1) &
                               + coord_pbc(k+i,1,1)*force_pbc(k+i,1,1,1)
          virial_omp(3,3,id+1) = virial_omp(3,3,id+1) &
                               + coord_pbc(2*k+i,1,1)*force_pbc(2*k+i,1,1,1)
        end do
      end if
      !$omp end parallel

    end if

    ! FEP
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then
      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do k = 1, 3
          do ix = 1, domain%num_atom(i)
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
              domain%translated_fep(ix,k,i)*domain%force_pbc_fep(ix,k,i,id+1)
          end do
        end do
      end do
      !$omp end parallel
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)

          do ic = 2, nthread
            force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          if (merge_force) force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
          if (merge_force) force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
          if (merge_force) force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, ic, force_tmp) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel 

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, ic, force_tmp)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    ! FEP
    !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = domain%force_pbc_fep(ix,1:3,i,1)
        do ic = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + domain%force_pbc_fep(ix,1:3,i,ic)
        end do
        force(1:3,ix,i) = force(1:3,ix,i) + force_tmp(1:3)
      end do
    end do
    !$omp end parallel


    do id = 1, nthread

      virial    (1,1) = virial    (1,1) + virial_omp    (1,1,id)
      virial    (2,2) = virial    (2,2) + virial_omp    (2,2,id)
      virial    (3,3) = virial    (3,3) + virial_omp    (3,3,id)
      virial_ext(1,1) = virial_ext(1,1) + virial_ext_omp(1,1,id)
      virial_ext(2,2) = virial_ext(2,2) + virial_ext_omp(2,2,id)
      virial_ext(3,3) = virial_ext(3,3) + virial_ext_omp(3,3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)
      energy%electric_field     = energy%electric_field     + efield_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond            &
                 + energy%angle           &
                 + energy%urey_bradley    &
                 + energy%dihedral        &
                 + energy%cmap            &
                 + energy%improper        &
                 + energy%electrostatic   &
                 + energy%van_der_waals   &
                 + energy%electric_field  &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_charmm_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_fep
  !> @brief        compute potential energy with AMBER99 force field for FEP
  !! @authors      HO
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    coord        : coordinates of target systems
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    reduce       : flag for reduce energy and virial
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    merge_force  : flag for merge force
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[inout] energy       : energy information
  !! @param[inout] atmcls_pbc   : atom class number
  !! @param[input] coord_pbc    : coordinates
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_long   : forces of target systems in long range
  !! @param[inout] force_omp    : temprary forces of target systems
  !! @param[inout] force_pbc    : pbc forces
  !! @param[inout] virial_cell  : virial correction due to pbc
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] virial_long  : virial term of target systems in long range
  !! @param[inout] virial_ext   : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8


  subroutine compute_energy_amber_fep(domain, enefunc, pairlist, boundary, coord,  &
                                  npt, reduce, nonb_ene, merge_force,          &
                                  nonb_limiter, energy, atmcls_pbc,            &
                                  coord_pbc, force, force_long,                &
                                  force_omp, force_pbc, virial_cell, virial,   &
                                  virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wip),               intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: trans(1:3)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    !$omp parallel private(i,ix,id)
    !$omp do
    do i = 1, ncell
      do ix = 1, natom
        force     (1,ix,i) = 0.0_wip
        force     (2,ix,i) = 0.0_wip
        force     (3,ix,i) = 0.0_wip
        force_long(1,ix,i) = 0.0_wip
        force_long(2,ix,i) = 0.0_wip
        force_long(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end parallel

    virial     (1:3,1:3)  = 0.0_dp
    virial_long(1:3,1:3)  = 0.0_dp
    virial_ext (1:3,1:3)  = 0.0_dp

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do 
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do

      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    ! FEP
    !$omp parallel do
    do id = 1, nthread
      domain%force_pbc_fep(1:natom,1:3,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp 

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    if (enefunc%gamd_use) &
      call error_msg('Compute_Energy> GaMD is not available in FEP')

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_fep( &
                              domain, enefunc, pairlist, boundary,   &
                              npt, nonb_ene, nonb_limiter,           &
                              atmcls_pbc, coord_pbc,                 &
                              force_long, force_omp, force_pbc,      &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
      else

        call compute_energy_nonbond_cutoff_fep( &
                              domain, enefunc, pairlist,             &
                              nonb_ene, force_pbc, virial_omp,       &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Amber> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      if (enefunc%local_restraint) then
        call compute_energy_dihed_localres_fep( &
                            domain, enefunc, coord, force_omp,     &
                            virial_omp, edihed_omp)
      else
        call compute_energy_dihed_fep( &
                            domain, enefunc, coord, force_omp,     &
                            virial_omp, edihed_omp)
      end if

      call compute_energy_cmap_fep( &
                            domain, enefunc, coord, force_omp,     &
                            virial_omp, ecmap_omp)


      ! improper energy
      !
      call compute_energy_improp_cos_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

      ! 1-4 interaction with linear table
      !
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_notable_fep( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

      end if

      if (enefunc%use_efield) &
        call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

      ! restraint energy
      !
      if (enefunc%restraint) then

        call compute_energy_restraints_fep( &
                            .true., .true., domain, boundary,      &
                            enefunc, coord, force_omp,             &
                            virial_omp, virial_ext_omp,            &
                            eposi_omp, energy%restraint_rmsd,      &
                            energy%rmsd, energy%restraint_distance,&
                            energy%restraint_emfit, energy%emcorr)

      end if
    end if

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, npt, nonb_ene, coord_pbc,    &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt) .and. &
        domain%nonbond_kernel /= NBK_Fugaku .and.       &
        domain%nonbond_kernel /= NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      if (domain%nonbond_kernel /= NBK_GPU) then

        do i = 1, ncell
          do k = 1, 3
            do ix = 1, domain%num_atom(i)
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                   coord_pbc(ix,k,i)*force_pbc(ix,k,i,id+1)
            end do
          end do
        end do
        do i = id+1, maxcell, nthread
          ic = domain%cell_pairlist1(1,i)
          jc = domain%cell_pairlist1(2,i)
          if (domain%virial_check(jc,ic)==1) then
            trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                   - trans(k)*virial_cell(k,i)
            end do
          end if
        end do

      else

        k = domain%num_atom_domain
        do i = id+1, domain%num_atom_domain, nthread
          virial_omp(1,1,id+1) = virial_omp(1,1,id+1) &
                               + coord_pbc(i,1,1)*force_pbc(i,1,1,1)
          virial_omp(2,2,id+1) = virial_omp(2,2,id+1) &
                               + coord_pbc(k+i,1,1)*force_pbc(k+i,1,1,1)
          virial_omp(3,3,id+1) = virial_omp(3,3,id+1) &
                               + coord_pbc(2*k+i,1,1)*force_pbc(2*k+i,1,1,1)
        end do

      end if
      !$omp end parallel

    end if

    ! FEP
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then
      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = 1, ncell
        do k = 1, 3
          do ix = 1, domain%num_atom(i)
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
              domain%translated_fep(ix,k,i)*domain%force_pbc_fep(ix,k,i,id+1)
          end do
        end do
      end do
      !$omp end parallel
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
          if (domain%nonbond_kernel /= NBK_GPU) then
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do
          else
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do
          end if
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          if (merge_force) force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
          if (merge_force) force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
          if (merge_force) force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic) 
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if
!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)
          if (merge_force) then
            force_tmp(1) = force_tmp(1) + force_long(1,ix,i)
            force_tmp(2) = force_tmp(2) + force_long(2,ix,i)
            force_tmp(3) = force_tmp(3) + force_long(3,ix,i)
          end if
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    ! FEP
    !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = domain%force_pbc_fep(ix,1:3,i,1)
        do ic = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + domain%force_pbc_fep(ix,1:3,i,ic)
        end do
        force(1:3,ix,i) = force(1:3,ix,i) + force_tmp(1:3)
      end do
    end do
    !$omp end parallel

    do id = 1, nthread

      virial    (1,1) = virial    (1,1) + virial_omp    (1,1,id)
      virial    (2,2) = virial    (2,2) + virial_omp    (2,2,id)
      virial    (3,3) = virial    (3,3) + virial_omp    (3,3,id)
      virial_ext(1,1) = virial_ext(1,1) + virial_ext_omp(1,1,id)
      virial_ext(2,2) = virial_ext(2,2) + virial_ext_omp(2,2,id)
      virial_ext(3,3) = virial_ext(3,3) + virial_ext_omp(3,3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do


    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%cmap          &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_amber_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_short_fep
  !> @brief        compute potential energy with charmm force field for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_short_fep(domain, enefunc, pairlist, boundary, &
                                         coord, npt, energy,                  &
                                         atmcls_pbc, coord_pbc,               &
                                         force, force_omp, force_pbc,         &
                                         virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, k, ki, ic, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

    !$omp parallel private(i,ix,id)

    !$omp do 
    do i = 1, ncell
      do ix = 1, natom
        force(1,ix,i) = 0.0_wip
        force(2,ix,i) = 0.0_wip
        force(3,ix,i) = 0.0_wip
      end do
    end do
    !$omp end do

    if (domain%nonbond_kernel /= NBK_GPU) then
      if (domain%nonbond_kernel /= NBK_Fugaku .and. &
          domain%nonbond_kernel /= NBK_Intel) then
        !$omp do
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, domain%num_atom(i)
              force_omp(1, ix, i, id) = 0.0_wp
              force_omp(2, ix, i, id) = 0.0_wp
              force_omp(3, ix, i, id) = 0.0_wp
              force_pbc(ix, 1, i, id) = 0.0_wp
              force_pbc(ix, 2, i, id) = 0.0_wp
              force_pbc(ix, 3, i, id) = 0.0_wp 
            end do
          end do
        end do
        !$omp end do
      else
        !$omp do
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
      end if
    else
      !$omp do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end do

      k = domain%num_atom_domain

      !$omp do
      do i = 1, 3*k
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end do
    end if

    ! FEP
    !$omp parallel do
    do id = 1, nthread
      domain%force_pbc_fep(1:natom,1:3,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do

    !$omp do
    do id = 1, nthread
      virial_omp    (1,1,id) = 0.0_dp
      virial_omp    (2,2,id) = 0.0_dp
      virial_omp    (3,3,id) = 0.0_dp
      virial_ext_omp(1,1,id) = 0.0_dp
      virial_ext_omp(2,2,id) = 0.0_dp
      virial_ext_omp(3,3,id) = 0.0_dp
      ebond_omp     (id) = 0.0_dp
      eangle_omp    (id) = 0.0_dp
      eurey_omp     (id) = 0.0_dp
      edihed_omp    (id) = 0.0_dp
      eimprop_omp   (id) = 0.0_dp
      ecmap_omp     (id) = 0.0_dp
      elec_omp      (id) = 0.0_dp
      evdw_omp      (id) = 0.0_dp
      eposi_omp     (id) = 0.0_dp
      efield_omp    (id) = 0.0_dp
    end do

    !$omp end parallel

    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel)       &
      virial_cell   (1:3,1:maxcell) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short_fep( &
                              domain, enefunc, pairlist, npt,        &
                              atmcls_pbc, coord_pbc,                 &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp)

      else

        call error_msg( &
        'Compute_Energy_Charmm_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Charmm_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then

      call compute_energy_dihed_localres_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

    else

      call compute_energy_dihed_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)

    end if

    ! improper energy
    !
    call compute_energy_improp_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

    ! cmap energy
    !
    call compute_energy_cmap_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear_fep( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

    if (enefunc%use_efield) &
      call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints_fep( &
                          .true., .true., domain, boundary,      &
                          enefunc, coord, force_omp,             &
                          virial_omp, virial_ext_omp,            &
                          eposi_omp, energy%restraint_rmsd,      &
                          energy%rmsd, energy%restraint_distance,&
                          energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, .false., .false.,            &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)
          if (domain%nonbond_kernel /= NBK_GPU) then
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do
          else
            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do
          end if
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel
   
    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, k, ix, ki, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)

          ki = k + ix
          force_tmp(1) = 0.0_wip
          force_tmp(2) = 0.0_wip
          force_tmp(3) = 0.0_wip
!ocl nosimd
          do ic = 1, nthread
  
            force_tmp(1) = force_tmp(1) &
                         + force_omp(1,ix,i,ic)+force_pbc(1,ki,1,ic)
            force_tmp(2) = force_tmp(2) &
                         + force_omp(2,ix,i,ic)+force_pbc(2,ki,1,ic)
            force_tmp(3) = force_tmp(3) &
                         + force_omp(3,ix,i,ic)+force_pbc(3,ki,1,ic)
          end do
          force(1,ix,i) = force_tmp(1)
          force(2,ix,i) = force_tmp(2)
          force(3,ix,i) = force_tmp(3)
        end do
      end do
      !$omp end parallel 

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, k, ix, ki, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)

          ki = k + ix
          force_tmp(1) = 0.0_wip
          force_tmp(2) = 0.0_wip
          force_tmp(3) = 0.0_wip
!ocl nosimd
          do ic = 1, nthread

            force_tmp(1) = force_tmp(1) &
                         + force_omp(1,ix,i,ic)+force_pbc(ki,1,1,ic)
            force_tmp(2) = force_tmp(2) &
                         + force_omp(2,ix,i,ic)+force_pbc(ki,2,1,ic)
            force_tmp(3) = force_tmp(3) &
                         + force_omp(3,ix,i,ic)+force_pbc(ki,3,1,ic)
          end do
          force(1,ix,i) = force_tmp(1)
          force(2,ix,i) = force_tmp(2)
          force(3,ix,i) = force_tmp(3)
        end do
      end do
      !$omp end parallel

    end if

    ! FEP
    !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = domain%force_pbc_fep(ix,1:3,i,1)
        do ic = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + domain%force_pbc_fep(ix,1:3,i,ic)
        end do
        force(1:3,ix,i) = force(1:3,ix,i) + force_tmp(1:3)
      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_charmm_short_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_short_fep
  !> @brief        compute potential energy with AMBER99 force field for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[inout] energy      : energy information
  !! @param[inout] atmcls_pbc  : atom class number
  !! @param[input] coord_pbc   : coordinates
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : forces
  !! @param[inout] virial_cell : virial term of target systems in cell
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber_short_fep(domain, enefunc, pairlist, boundary,   &
                                        coord, npt, energy, atmcls_pbc,        &
                                        coord_pbc, force, force_omp, force_pbc,&
                                        virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wip),               intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    integer,                 intent(inout) :: atmcls_pbc(:)
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wip),               intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: efield_omp (nthread)
    real(wip)                :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, k, start_i
    integer                  :: omp_get_thread_num


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_wip
    virial    (1:3,1:3)             = 0.0_dp 
    virial_ext(1:3,1:3)             = 0.0_dp 

    if (domain%nonbond_kernel /= NBK_GPU) then

      if (domain%nonbond_kernel == NBK_Fugaku .or. &
          domain%nonbond_kernel == NBK_Intel) then

        !$omp parallel do private(k, id, i, ix)
        do id = 1, nthread
          do i = 1, ncell
            do ix = 1, natom
              force_omp(1,ix,i,id) = 0.0_wp
              force_omp(2,ix,i,id) = 0.0_wp
              force_omp(3,ix,i,id) = 0.0_wp
            end do
          end do
          if (domain%nonbond_kernel == NBK_Fugaku) then
            do i = 1, domain%num_atom_domain
              force_pbc(1,i,1,id) = 0.0_wp
              force_pbc(2,i,1,id) = 0.0_wp
              force_pbc(3,i,1,id) = 0.0_wp
            end do
          else if (domain%nonbond_kernel == NBK_Intel) then
            do i = 1, domain%num_atom_domain
              force_pbc(i,1,1,id) = 0.0_wp
              force_pbc(i,2,1,id) = 0.0_wp
              force_pbc(i,3,1,id) = 0.0_wp
            end do
          end if
        end do
        !$omp end parallel do
          
      else

        !$omp parallel do
        do id = 1, nthread
          force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
          force_pbc(1:natom,1:3,1:ncell,id) = 0.0_wp
        end do
        !$omp end parallel do

      end if

    else
      !$omp parallel do
      do id = 1, nthread
        force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      end do
      !$omp end parallel do
      !$omp parallel do
      do i = 1, 3*domain%num_atom_domain
        force_pbc(i,1,1,1) = 0.0_wp
      end do
      !$omp end parallel do
    end if

    ! FEP
    !$omp parallel do
    do id = 1, nthread
      domain%force_pbc_fep(1:natom,1:3,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp 
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp 
    virial_cell   (1:3,1:maxcell) = 0.0_dp 
    ebond_omp     (1:nthread) = 0.0_dp 
    eangle_omp    (1:nthread) = 0.0_dp 
    eurey_omp     (1:nthread) = 0.0_dp 
    edihed_omp    (1:nthread) = 0.0_dp 
    eimprop_omp   (1:nthread) = 0.0_dp 
    ecmap_omp     (1:nthread) = 0.0_dp 
    elec_omp      (1:nthread) = 0.0_dp 
    evdw_omp      (1:nthread) = 0.0_dp 
    eposi_omp     (1:nthread) = 0.0_dp 
    efield_omp    (1:nthread) = 0.0_dp 

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short_fep( &
                              domain, enefunc, pairlist, npt,        & 
                              atmcls_pbc, coord_pbc,                 &
                              force_omp, force_pbc,                  &
                              virial_cell, virial_omp)

      else

        call error_msg( &
        'Compute_Energy_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
    else
      call compute_energy_dihed_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp_cos_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, eimprop_omp)

    ! cmap energy
    !
    call compute_energy_cmap_fep( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_notable_fep( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear( &
                              domain, enefunc, atmcls_pbc,           &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_omp, evdw_omp, elec_omp)

    if (enefunc%use_efield) &
      call compute_energy_efield( &
                              domain, enefunc, coord, force_omp,     &
                              virial_omp, efield_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints_fep( &
                          .true., .true., domain, boundary,      &
                          enefunc, coord, force_omp,             &
                          virial_omp, virial_ext_omp,            &
                          eposi_omp, energy%restraint_rmsd,      &
                          energy%rmsd, energy%restraint_distance,&
                          energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
    if (domain%nonbond_kernel == NBK_GPU .and. enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
                              pairlist, .false., .false.,            &
                              coord_pbc, force_omp, force_pbc,       &
                              virial_cell, virial_omp,               &
                              elec_omp, evdw_omp)
    end if

    ! gather values
    !
    if (domain%nonbond_kernel /= NBK_Fugaku .and. &
        domain%nonbond_kernel /= NBK_Intel  .and. &
        domain%nonbond_kernel /= NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        do ix = 1, domain%num_atom(i)
          force_tmp(1:3) = force_omp(1:3,ix,i,1)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,1)

          if (domain%nonbond_kernel /= NBK_GPU) then

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
              force_tmp(1:3) = force_tmp(1:3) + force_pbc(ix,1:3,i,ic)
            end do

          else

            do ic = 2, nthread
              force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,ic)
            end do

          end if

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_GPU) then

      !$omp parallel default(shared) private(id, i, ix, force_tmp, ic, k, start_i)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      k = domain%num_atom_domain
      do i = id+1, ncell, nthread
        start_i = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(    start_i+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(  k+start_i+ix,1,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(2*k+start_i+ix,1,1,1)
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
          end do
          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Fugaku) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(1,k+ix,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(2,k+ix,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(3,k+ix,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    else if (domain%nonbond_kernel == NBK_Intel) then

      !$omp parallel default(shared) private(id, i, ix, k, force_tmp, ic)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do i = id+1, ncell, nthread
        k = domain%start_atom(i)
        do ix = 1, domain%num_atom(i)
          force_tmp(1) = force_omp(1,ix,i,1)
          force_tmp(2) = force_omp(2,ix,i,1)
          force_tmp(3) = force_omp(3,ix,i,1)
          force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,1)
          force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,1)
          force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,1)

!ocl nosimd
          do ic = 2, nthread
            force_tmp(1) = force_tmp(1) + force_omp(1,ix,i,ic)
            force_tmp(2) = force_tmp(2) + force_omp(2,ix,i,ic)
            force_tmp(3) = force_tmp(3) + force_omp(3,ix,i,ic)
            force_tmp(1) = force_tmp(1) + force_pbc(k+ix,1,1,ic)
            force_tmp(2) = force_tmp(2) + force_pbc(k+ix,2,1,ic)
            force_tmp(3) = force_tmp(3) + force_pbc(k+ix,3,1,ic)
          end do

          force(1:3,ix,i) = force_tmp(1:3)
        end do
      end do
      !$omp end parallel

    end if

    ! FEP
    !$omp parallel default(shared) private(id, i, ix, force_tmp, ic) 
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = domain%force_pbc_fep(ix,1:3,i,1)
        do ic = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + domain%force_pbc_fep(ix,1:3,i,ic)
        end do
        force(1:3,ix,i) = force(1:3,ix,i) + force_tmp(1:3)
      end do
    end do
    !$omp end parallel

    return

  end subroutine compute_energy_amber_short_fep

end module sp_energy_mod
