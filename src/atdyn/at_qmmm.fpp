!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_qmmm_mod
!> @brief   perform QM/MM calculation
!! @authors Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY), Kenta Yamada (KYMD)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_qmmm_mod

  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_boundary_str_mod
  use at_output_str_mod
  use molecules_str_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use fileio_mod
  use fileio_control_mod
  use fileio_rst_mod
  use fileio_psf_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! atomic symbols
  character(3)            :: qm_atomic_symbol(118) = (/ &
      'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ', 'O  ', 'F  ', 'Ne ', &
      'Na ', 'Mg ', 'Al ', 'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', &
      'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ', &
      'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', 'Rr ', 'Sr ', 'Y  ', 'Zr ', &
      'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ', &
      'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ', &
      'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ', &
      'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ', &
      'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', &
      'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ', &
      'Md ', 'No ', 'Lr ', 'Rf ', 'Db ', 'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ', &
      'Rg ', 'Cn ', 'Uut', 'Fl ', 'Uup', 'Lv ', 'Uus', 'Uuo' /)

  ! atomic mass (Taken from Standard Atomic Weights 2013)
  real(wp)               :: qm_atomic_mass(118) = (/ &
      1.0080_wp,   4.0026_wp,   6.9675_wp,   9.0122_wp,  10.8135_wp, &
     12.0106_wp,  14.0069_wp,  15.9994_wp,  18.9984_wp,  20.1797_wp, &
     22.9898_wp,  24.3050_wp,  26.9815_wp,  28.0850_wp,  30.9738_wp, &
     32.0675_wp,  35.4515_wp,  39.948_wp,   39.0983_wp,  40.078_wp,  &
     44.9559_wp,  47.867_wp,   50.9415_wp,  51.9961_wp,  54.9380_wp, &
     55.845_wp,   58.9332_wp,  58.6934_wp,  63.546_wp,   65.38_wp,   &
     69.723_wp,   72.63_wp,    74.9216_wp,  78.96_wp,    79.904_wp,  &
     83.798_wp,   85.4678_wp,  87.62_wp,    88.9059_wp,  91.224_wp,  &
     92.9064_wp,  95.96_wp,     0.0000_wp, 101.07_wp,   102.9055_wp, &
    106.42_wp,   107.8682_wp, 112.411_wp,  114.818_wp,  118.710_wp,  &
    121.760_wp,  127.60_wp,   126.9044_wp, 131.293_wp,  132.9055_wp, &
    137.327_wp,  138.9055_wp, 140.116_wp,  140.9077_wp, 144.242_wp,  &
      0.0000_wp, 150.36_wp,   151.964_wp,  157.25_wp,   158.9254_wp, &
    162.500_wp,  164.9303_wp, 167.259_wp,  168.9342_wp, 173.054_wp,  &
    174.9668_wp, 178.49_wp,   180.9479_wp, 183.84_wp,   186.207_wp,  &
    190.23_wp,   192.217_wp,  195.084_wp,  196.9666_wp, 200.59_wp,   &
    204.384_wp,  207.2_wp,    208.9804_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, 232.0381_wp, &
    231.0359_wp, 238.0289_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp,   0.0000_wp, &
      0.0000_wp,   0.0000_wp,   0.0000_wp /)

  real(wp), parameter     :: deuterium_mass = 2.0141_wp
  !
  real(wp), parameter     :: mass_error = 0.05_wp
  !
  !   1 Debye = 1/299792458 x 1e-21 = 3.335641e-30 Cm = 3.934303e-1 au
  !     using 1 Bohr  = 5.291772109e-11 m,   1 e     = 1.602176634e-10 C
  real(wp), parameter     :: conv_debye_au = 3.934303E-01

  integer, parameter      :: DEFAULT_INT = 10000

  !! For GBSA read
  real(wp) :: salt_cons, temperature, eps_solute, eps_solvent

  ! structures
  type, public :: s_qmmm_info
    logical               :: do_qmmm = .false.
    integer               :: qmtyp = 0
    integer               :: qmsave_period       = 10
    integer               :: exclude_charge      = ExcludeChargeGroup
    integer               :: qmmaxtrial          = 1
    logical               :: qm_debug            = .false.
    logical               :: qminfo              = .false.
    real(wp)              :: linkbond            = -1.0_wp
    integer               :: qm_total_charge     = DEFAULT_INT

    character(MaxFilename) :: qmcnt = ''
    character(MaxFilename) :: qmexe = ''
    character(MaxFilename) :: qmatm_select_index  = ''
    character(MaxFilename) :: workdir  = ''
    character(MaxFilename) :: savedir  = ''
    character(MaxFilename) :: basename = ''

    ! parameters not used now
    real(wp)              :: mm_cutoffdist       = -1.0_wp
    logical               :: mm_cutoff_byresidue = .true.
    character(256)        :: qm_exemode = 'system'
    integer               :: qm_nprocs  = 1
  end type s_qmmm_info

  ! subroutines
  public  :: show_ctrl_qmmm
  public  :: read_ctrl_qmmm
  public  :: setup_qmmm
  public  :: setup_qmmm_directory
  private :: qm_check_input
  public  :: qm_atomic_number
  public  :: qm_generate_input
  public  :: qm_read_output
  public  :: qm_post_process
  public  :: write_qminfo
  public  :: output_qm_charge
  public  :: qm_monitor
  public  :: qmmm_finalize
  private :: setup_mmlist
  public  :: update_mmlist
#ifdef QSIMULATE
  public :: qm_fill_input_qsimulate
  private :: qm_fill_input_qsimulate_gbsa
  private :: qm_fill_input_qsimulate_retry
#endif
  private :: qm_check_input_qsimulate
  private :: qm_check_input_qchem
  private :: qm_check_input_gaussian
  private :: qm_check_input_terachem
  private :: qm_check_input_dftbplus
  private :: qm_check_input_molpro
  private :: qm_check_input_gaussian_frwf
  private :: qm_check_input_orca
  private :: qm_generate_input_qchem
  private :: qm_generate_input_gaussian
  private :: qm_generate_input_terachem
  private :: qm_generate_input_dftbplus
  private :: qm_generate_input_molpro
  private :: qm_generate_input_gaussian_frwf
  private :: qm_generate_input_orca
  private :: qm_read_output_qchem
  private :: qm_read_output_gaussian
  private :: qm_read_output_terachem
  private :: qm_read_output_dftbplus
  private :: qm_read_output_molpro
  private :: qm_read_output_gaussian_frwf
  private :: qm_read_output_orca
  private :: count_qmmmbonds
  private :: generate_linkatoms
  private :: check_linkmm
  private :: rearrange_term_list
  private :: rearrange_impr_list
  private :: merge_qm_and_linkh
  private :: separate_qm_and_linkh

contains


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_qmmm
  !> @brief        show QMMM section usage
  !! @authors      YA
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "remd", "min", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_qmmm(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd', 'min', 'vib', 'rpath')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') '# group1        = an:CA              # QM group 1'
        write(MsgOut,'(A)') '# group2        = an:CA and resno:1  # QM group 2'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[QMMM]'
        write(MsgOut,'(A)') '## Example of QM/MM'
        write(MsgOut,'(A)') '# qmtyp              = QChem      # QM solver'
        write(MsgOut,'(A)') '# qmcnt              = qchem.inp  # template file for qmtyp input'
        write(MsgOut,'(A)') '# qmexe              = runqc.sh   # QM run script'
        write(MsgOut,'(A)') '# qmatm_select_index = 1 2        # QM groups'
        write(MsgOut,'(A)') '# exclude_charge     = atom       # atom or group'
        write(MsgOut,'(A)') '# qmmaxtrial         = 1          # maximum number of QM restarts'
        write(MsgOut,'(A)') '# workdir            = qmmm       # a working directory for QM jobs'
        write(MsgOut,'(A)') '# savedir            = qmmm_save  # a directory to save QM jobs'
        write(MsgOut,'(A)') '# basename           = job        # the name of QM files' 
        write(MsgOut,'(A)') '# qmsave_period      = 1          # period for saving QM files'
        write(MsgOut,'(A)') '# qm_total_charge    = 0          # QM total charge (must be integer)'
        write(MsgOut,'(A)') '# linkbond           = 1.0        # QM-MM link bond distance'
!ky        write(MsgOut,'(A)') '# qm_debug           = yes   # enable debug mode'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_qmmm
  !> @brief        read QMMM section in the control file
  !! @authors      YA, KY
  !! @param[in]    handle   : unit number of control fileha
  !! @param[in]    min_info : QMMM section in control parameters 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_qmmm(handle, ctrlfile, qmmm_info)

    ! parameters
    character(*),            parameter     :: Section = 'QMMM'

    ! formal arguments
    integer,                 intent(in)    :: handle
    character(*),            intent(in)    :: ctrlfile
    type(s_qmmm_info),       intent(inout) :: qmmm_info

    ! local variables
    character(256)           :: line
    character(6)             :: word
    logical                  :: found
    integer                  :: unit_no, nsta, nend


    ! search QMMM section in control file
    !
    qmmm_info%do_qmmm = find_ctrlfile_section(ctrlfile, Section)
    if (.not. qmmm_info%do_qmmm) then
      return
    end if

    call begin_ctrlfile_section(handle, Section)
    ! read parameters from control file
    !
    call read_ctrlfile_type   (handle, Section, 'qmtyp',              &
                               qmmm_info%qmtyp, QMtypTypes)
    call read_ctrlfile_string (handle, Section, 'qmcnt',              &
                               qmmm_info%qmcnt)
    call read_ctrlfile_string (handle, Section, 'qmexe',              &
                               qmmm_info%qmexe)
    call read_ctrlfile_string (handle, Section, 'qmatm_select_index', &
                               qmmm_info%qmatm_select_index)
    call read_ctrlfile_string (handle, Section, 'workdir',            &
                               qmmm_info%workdir)
    call read_ctrlfile_string (handle, Section, 'savedir',            &
                               qmmm_info%savedir)
    call read_ctrlfile_string (handle, Section, 'basename',           &
                               qmmm_info%basename)
    call read_ctrlfile_integer(handle, Section, 'qmsave_period',      &
                               qmmm_info%qmsave_period)
    call read_ctrlfile_type   (handle, Section, 'exclude_charge',     &
                               qmmm_info%exclude_charge, ExcludeChargeTypes)
    call read_ctrlfile_string (handle, Section, 'qm_exemode',         &
                               qmmm_info%qm_exemode)
    call read_ctrlfile_integer(handle, Section, 'qm_nprocs',          &
                               qmmm_info%qm_nprocs)
    call read_ctrlfile_integer(handle, Section, 'qmmaxtrial',         &
                               qmmm_info%qmmaxtrial)
    call read_ctrlfile_logical(handle, Section, 'qm_debug',           &
                               qmmm_info%qm_debug)
    call read_ctrlfile_string (handle, Section, 'qmexe',              &
                               qmmm_info%qmexe)
    call read_ctrlfile_logical(handle, Section, 'qminfo',             &
                               qmmm_info%qminfo)
    call read_ctrlfile_real   (handle, Section, 'mm_cutoffdist',      &
                               qmmm_info%mm_cutoffdist)
    call read_ctrlfile_logical(handle, Section, 'mm_cutoff_byresidue', &
                               qmmm_info%mm_cutoff_byresidue)
    call read_ctrlfile_real   (handle, Section, 'linkbond',           &
                               qmmm_info%linkbond)
    call read_ctrlfile_integer(handle, Section, 'qm_total_charge',    &
                               qmmm_info%qm_total_charge)
    call end_ctrlfile_section(handle)

    if (qmmm_info%qmtyp /= QMtypQSimulate .and. &
        qmmm_info%qmtyp /= QMtypQCHEM     .and. &
        qmmm_info%qmtyp /= QMtypG09       .and. &
        qmmm_info%qmtyp /= QMtypTERACHEM  .and. &
        qmmm_info%qmtyp /= QMtypDFTBPLUS  .and. &
        qmmm_info%qmtyp /= QMtypG09_FRWF)     &
      call error_msg('Read_Ctrl_QMMM> QM program is not defined')

    if (trim(qmmm_info%qmcnt) .eq. '') &
      call error_msg('Read_Ctrl_QMMM> QM template file is not defined')

    if (qmmm_info%qmtyp /= QMtypQSimulate .and. &
        trim(qmmm_info%qmexe) .eq. '') &
      call error_msg('Read_Ctrl_QMMM> QM run script is not defined')

#ifndef QSIMULATE
    if (qmmm_info%qmtyp == QMtypQSimulate) then
      call error_msg('Read_Ctrl_QMMM> qmtyp = qsimulate, but GENESIS is not '//&
              'configured with --enable-qsimulate. Please recompile the program')
    end if
#endif

    if (trim(qmmm_info%qmatm_select_index) .eq. '') &
      call error_msg('Read_Ctrl_QMMM> No QM atoms defined')

    if (qmmm_info%qmmaxtrial < 0) &
      call error_msg('Read_Ctrl_QMMM> Illegal value of qmmaxtrial')

    if (qmmm_info%exclude_charge /= ExcludeChargeAtom .and. &
        qmmm_info%exclude_charge /= ExcludeChargeGroup .and. &
        qmmm_info%exclude_charge /= ExcludeChargeAtomAndDistribute)    &
      call error_msg('Read_Ctrl_QMMM> Illegal value of exclude_charge')

    if (qmmm_info%workdir == '') &
      qmmm_info%workdir = 'qmmm'

! -- 
    select case(qmmm_info%qmtyp)
      case (QMtypQSimulate, QMtypQCHEM, QMtypG09, QMtypTERACHEM, QMtypDFTBPLUS, &
            QMtypG09_FRWF)
        if (trim(qmmm_info%qm_exemode) .eq. 'spawn') then
          write(MsgOut, '(a,a)') &
          'Read_Ctrl_QMMM> Warning: qm_exemode = spawn is not available for ', &
          trim(qmmm_info%qm_exemode)
          qmmm_info%qm_exemode = 'system'
        end if
    end select

    if (trim(qmmm_info%qm_exemode) .eq. 'system') then
      if (qmmm_info%qm_nprocs > 1) then
        write(MsgOut, '(a,i0,a)') 'Read_Ctrl_QMMM> Warning: qm_nprocs = ', &
                                  qmmm_info%qm_nprocs, ' is ignored'
        qmmm_info%qm_nprocs = 1
      end if
#ifndef HAVE_MPI_GENESIS
    else if (trim(qmmm_info%qm_exemode) .eq. 'spawn') then
      call error_msg('Read_Ctrl_QMMM> qm_exemode = spawn requires MPI')
#endif
    else
      call error_msg('Read_Ctrl_QMMM> Illegal value of qm_exemode')
    end if
! -- 

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(a)') 'Read_Ctrl_QMMM> Parameters of QM/MM'

      write(MsgOut,'(a,a)')       &
            '  qmtyp              = ', trim(QMtypTypes(qmmm_info%qmtyp))
      write(MsgOut,'(a,a)')       &
            '  qmcnt              = ', trim(qmmm_info%qmcnt)
      if (len(trim(qmmm_info%qmexe)) > 0) then
        write(MsgOut,'(a,a)')       &
              '  qmexe              = ', trim(qmmm_info%qmexe)
      end if
      write(MsgOut,'(a,a)')       &
            '  workdir            = ', trim(qmmm_info%workdir)
      if (qmmm_info%savedir /= '') then
        write(MsgOut,'(a,a)')     &
            '  savedir            = ', trim(qmmm_info%savedir)
      else
        write(MsgOut,'(a)')       &
            '  savedir            = none'
      end if
      if (qmmm_info%basename /= '') then
        write(MsgOut,'(a,a)')     &
            '  basename           = ', trim(qmmm_info%basename)
      else
        write(MsgOut,'(a)')       &
            '  basename           = blank'
      end if
      write(MsgOut,'(A,I0)')      &
            '  qmsave_period      = ', qmmm_info%qmsave_period
      write(MsgOut,'(a,I0)')      &
            '  qmmaxtrial         = ', qmmm_info%qmmaxtrial
      write(MsgOut,'(a,a)')       &
            '  qmatm_select_index = ', trim(qmmm_info%qmatm_select_index)
      write(MsgOut,'(a,a)')       &
            '  exclude_charge     = ', &
              trim(ExcludeChargeTypes(qmmm_info%exclude_charge))
      if (qmmm_info%qm_total_charge /= DEFAULT_INT) then
        write(MsgOut,'(a,I0)')   &
              '  qm_total_charge    = ', qmmm_info%qm_total_charge
      end if
!ky      if (qmmm_info%mm_cutoffdist > 0.0_wp) then
!ky        write(MsgOut,'(a,f10.3)') &
!ky            '  mm_cutoffdist       = ', qmmm_info%mm_cutoffdist
!ky        if (qmmm_info%mm_cutoff_byresidue) then
!ky          write(MsgOut,'(a)') '  mm_cutoff_byresidue = yes'
!ky        else
!ky          write(MsgOut,'(a)') '  mm_cutoff_byresidue = no'
!ky        end if
!ky      end if
!ky      write(MsgOut,'(a,a)')       &
!ky            '  qm_exemode         = ', trim(qmmm_info%qm_exemode)
!ky      write(MsgOut,'(a,i0)')      &
!ky            '  qm_nprocs          = ', qmmm_info%qm_nprocs
!ky      if (qmmm_info%qm_debug) then
!ky        write(MsgOut,'(a)')       &
!ky            '  qm_debug           = yes'
!ky      end if

      write(MsgOut,'(a)') ' '
    end if

    return

  end subroutine read_ctrl_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_qmmm
  !> @brief        setup QM/MM information
  !> @authors      YA, KY
  !! @param[in]    qmmm_info  : QMMM section control parameters information
  !! @param[in]    sel_info   : SELECTION section control parameters information
  !! @param[in]    boundary   : boundary information
  !! @param[in]    psf        : PSF information
  !! @param[in]    rst        : rst information
  !! @param[in]    forcefield : forcefield
  !! @param[in]    molecule   : molecular information
  !! @param[out]   qmmm       : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_qmmm(qmmm_info, sel_info, boundary, psf, rst, forcefield, &
                        molecule, qmmm)

    ! formal arguments
    type(s_qmmm_info),       intent(in)    :: qmmm_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_psf),             intent(in)    :: psf
    type(s_rst),             intent(in)    :: rst
    type(s_boundary),        intent(in)    :: boundary
    integer,                 intent(in)    :: forcefield
    type(s_molecule),        intent(inout) :: molecule
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    real(wp)                 :: mm_mass
    integer                  :: atomic_no
    integer                  :: i, j, offset
    integer                  :: ngroup, igroup, natom, iatom
    integer                  :: istart, iend, mm_start, mm_end, mm_block
    integer                  :: num_partner, list_partner(10)
    integer                  :: temp
    integer                  :: kalloc_stat, kdealloc_stat
    character(4)             :: mm_name
    character                :: num*5
    character(MaxFilename)   :: path
    logical                  :: found_atom, setup_qmmm_error

    type(s_selatoms), allocatable :: selatoms(:)
    integer,          allocatable :: group_list(:), atom_no(:)


    ! flag of QM/MM
    !
    qmmm%do_qmmm = qmmm_info%do_qmmm

    if (.not. qmmm%do_qmmm) return

    if ((forcefield /= ForcefieldCHARMM) .and. \
        (forcefield /= ForcefieldAMBER)) then
        if (main_rank) then
          call error_msg("Setup_QMMM> Check force field")
        end if
    end if

    setup_qmmm_error = .false.

    ! setup qmmm
    !
    ! read from the input
    qmmm%qmtyp               = qmmm_info%qmtyp
    qmmm%qmcnt               = qmmm_info%qmcnt
    qmmm%qmmaxtrial          = qmmm_info%qmmaxtrial
    qmmm%exclude_charge      = qmmm_info%exclude_charge
    qmmm%save_qminfo         = qmmm_info%qminfo
    qmmm%qmsave_period       = qmmm_info%qmsave_period
    qmmm%mm_cutoffdist       = qmmm_info%mm_cutoffdist
    qmmm%mm_cutoff_byresidue = qmmm_info%mm_cutoff_byresidue
    qmmm%qm_exemode          = qmmm_info%qm_exemode
    qmmm%qm_nprocs           = qmmm_info%qm_nprocs
    qmmm%linkbond            = qmmm_info%linkbond
    qmmm%qm_total_charge     = dble(qmmm_info%qm_total_charge)

    ! Set forcefield-specifc parameters
    !
    if (forcefield .eq. ForcefieldAMBER) then
      if (qmmm%exclude_charge /= ExcludeChargeAtomAndDistribute) then
        if (main_rank) write(MsgOut, '(a)') "Setup_QMMM> Warning: Force to exclude_charge = amber for AMBER forcefifeld"
        qmmm%exclude_charge = ExcludeChargeAtomAndDistribute
      end if
    end if
    if (forcefield .eq. ForcefieldCHARMM) then
      if (qmmm%exclude_charge == ExcludeChargeAtomAndDistribute) then
        if (main_rank) write(MsgOut, '(a)') "Setup_QMMM> Warning: exclude_charge = amber is not valid for CHARMM forcefield"
        if (main_rank) write(MsgOut, '(a)') "                     exclude_charge = group will be used."
        qmmm%exclude_charge = ExcludeChargeGroup
      end if
    end if
    if (qmmm%linkbond < 0.0_wp) then
      if (forcefield .eq. ForcefieldAMBER) then
        if (main_rank) write(MsgOut, '(a,f5.3,a)') "Setup_QMMM> Link bond length is set to ", LinkBondAMBER, " angstrom"
        qmmm%linkbond = LinkBondAMBER
      else if (forcefield .eq. ForcefieldCHARMM) then
        if (main_rank) write(MsgOut, '(a,f5.3,a)') "Setup_QMMM> Link bond length is set to ", LinkBondCHARMM, " angstrom"
        qmmm%linkbond = LinkBondCHARMM
      end if
    end if

    if (forcefield .eq. ForcefieldAMBER .and.  &
        qmmm_info%qm_total_charge == DEFAULT_INT) then
      call error_msg('Fatal ERROR: AMBER ff requires qm_total_charge as an input.')
    end if

    ! directory and file names of QM jobs
    write(num,'(i0)') (my_country_no+1)
    if (nrep_per_proc > 1) write(num,'(i0)') (my_replica_no+1)
    qmmm%workdir = trim(qmmm_info%workdir)//'.'//num
    if (replica_main_rank) &
      call system('mkdir -p '//trim(qmmm%workdir)//' > /dev/null 2>&1')

    if (qmmm_info%savedir .ne. '') then
      qmmm%savedir  = trim(qmmm_info%savedir)//'.'//num
      if (replica_main_rank) &
        call system('mkdir -p '//trim(qmmm%savedir)//' > /dev/null 2>&1')
    else
      qmmm%savedir  = qmmm%workdir
    end if
    qmmm%qmbase0    = qmmm_info%basename

    ! qmexe with an absolute path
    if (qmmm_info%qmexe(1:1) .eq. "/") then 
      qmmm%qmexe = trim(qmmm_info%qmexe)

    else if (qmmm_info%qmexe(1:2) .eq. "~/") then
      call getenv('HOME', path)
      qmmm%qmexe = trim(path)//trim(qmmm_info%qmexe(2:))

    else
      call getcwd(path)
      qmmm%qmexe = trim(path)//'/'//trim(qmmm_info%qmexe)
    end if
    
    ! initialize
    qmmm%dryrun       = .false.
    qmmm%ene_only     = .false.
    qmmm%qm_count     = 0
    qmmm%savefile     = .false.
    qmmm%qm_classical = .false.
    qmmm%qm_get_esp   = .false.

    ! debug option
    qmmm%qm_debug    = qmmm_info%qm_debug


    if (main_rank) write(MsgOut,'(a)') "Setup_QMMM> Setup QM region"

    ! Boundary for QM/MM run is currently limited to NOBC
    if (boundary%type /= BoundaryTypeNOBC .and. main_rank) then
      write(MsgOut,*)
      write(MsgOut,'(a)') '  Error in [BOUNDARY]!'
      write(MsgOut,'(a)') '  type = '// &
        trim(BoundaryTypeTypes(boundary%type))//' is not supported for QMMM.'
      write(MsgOut,'(a)') '  type = NOBC must be used.'
      write(MsgOut,*)
      call error_msg('Setup_QMMM> Error while setting up QMMM.')

    end if

    ! number of QM and MM atoms
    ngroup = split_num(trim(qmmm_info%qmatm_select_index))

    allocate(group_list(ngroup), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    call split(ngroup, ngroup, qmmm_info%qmatm_select_index, group_list)

    allocate(selatoms(ngroup), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    ! Avoid segmentation falut due to the lack of [SELECTION]
    !
    if (size(sel_info%groups) == 0) then
      deallocate(selatoms, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      deallocate(group_list, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      call error_msg('Setup_QMMM> Selection for QM region was NOT found!')
    end if

    natom = 0
    do i = 1, ngroup
      igroup = group_list(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom  = natom + size(selatoms(i)%idx)
    end do

    qmmm%qm_natoms = natom
    qmmm%mm_natoms = molecule%num_atoms - qmmm%qm_natoms

    if (qmmm%qm_natoms == 0) then
      deallocate(group_list, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      deallocate(selatoms, stat = kdealloc_stat)
      if (kdealloc_stat /= 0) call error_msg_dealloc
      call error_msg('Setup_QMMM> Number of QM atoms is zero')
    end if

    ! allocate array to save forces on QM and MM atoms
    allocate(qmmm%qm_force(3, qmmm%qm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    allocate(qmmm%mm_force(3, qmmm%mm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    ! allocate array to save QM charges
    allocate(qmmm%qm_charge(qmmm%qm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    allocate(qmmm%qm_charge_save(qmmm%qm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    qmmm%qm_charge      = 0.0_wp
    qmmm%qm_charge_save = 0.0_wp
    qmmm%is_qm_charge   = .false.

    if (allocated(rst%qm_charge)) then
      ! retrieve QM charges from rstfile
      if (size(rst%qm_charge) == qmmm%qm_natoms) then
        qmmm%qm_charge = rst%qm_charge
        qmmm%is_qm_charge = .true.
      end if
    end if

    ! list of QM atom indices
    !
    allocate(qmmm%qmatom_id(qmmm%qm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    allocate(qmmm%qm_atomic_no(qmmm%qm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    allocate(qmmm%mmatom_id(qmmm%mm_natoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    offset = 0
    do i = 1, ngroup
      igroup = group_list(i)
      natom  = size(selatoms(i)%idx)
      qmmm%qmatom_id(offset+1:offset+natom)     = selatoms(i)%idx(1:natom)
      offset = offset + natom
    end do

    deallocate(selatoms, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(group_list, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc

    ! sort QM atom indices in ascending order
    !
    do i = qmmm%qm_natoms, 2, -1
      do j = 1, i - 1
        if (qmmm%qmatom_id(j) > qmmm%qmatom_id(j+1)) then
          temp = qmmm%qmatom_id(j)
          qmmm%qmatom_id(j)   = qmmm%qmatom_id(j+1)
          qmmm%qmatom_id(j+1) = temp
        end if
      end do
    end do

    ! QM atomic numbers
    !
    found_atom = .true.
    do i = 1, qmmm%qm_natoms
      iatom = qmmm%qmatom_id(i)
      mm_mass = molecule%mass(iatom)
      call qm_atomic_number(mm_mass, atomic_no)
      qmmm%qm_atomic_no(i) = atomic_no
      if (atomic_no < 0) found_atom = .false.
    end do

    if (.not. found_atom .and. main_rank) then
      write(MsgOut,'(a)') '  Error! Unknown atom is detected.'
      write(MsgOut,'(a)') &
        '  The following MM atoms cannnot be assigned to QM atoms.'
      do i = 1, qmmm%qm_natoms
        if (qmmm%qm_atomic_no(i) < 0) then
          iatom = qmmm%qmatom_id(i)
          write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6, "mass = ",f12.4)') &
            iatom,                         &
            molecule%segment_name(iatom),  &
            molecule%residue_no(iatom),    &
            molecule%residue_name(iatom),  &
            molecule%atom_name(iatom),     &
            molecule%atom_cls_name(iatom), &
            molecule%mass(iatom)
        end if
      end do
      write(MsgOut,'(a)') &
        '  QM atoms are assigned based on the mass of MM atoms.'

      setup_qmmm_error = .true.

    end if

    ! Punch out QM atom info.
    !
    if (main_rank) then
      write(MsgOut,'(a)') "  QM assignment info"
      do i = 1, qmmm%qm_natoms
        iatom   = qmmm%qmatom_id(i)
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6," assigned to QM atom ",i4, &
        & " of element: ",a3,x,i3)')     &
          iatom,                         &
          molecule%segment_name(iatom),  &
          molecule%residue_no(iatom),    &
          molecule%residue_name(iatom),  &
          molecule%atom_name(iatom),     &
          molecule%atom_cls_name(iatom), &
          i,                             &
          qm_atomic_symbol(qmmm%qm_atomic_no(i)),     &
          qmmm%qm_atomic_no(i)
      end do
      write(MsgOut,'(a,i0)') "  number of QM atoms = ", qmmm%qm_natoms
      write(MsgOut, '(a)') ' '
      if (qmmm%is_qm_charge) &
         write(MsgOut,'(a)') "  QM charges retreived from rstfile."
      write(MsgOut, '(a)') ' '
    end if

    ! MM atom indices
    !
    allocate(atom_no(molecule%num_atoms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    do i = 1, molecule%num_atoms
       atom_no(i) = i
    end do

    mm_end = qmmm%qmatom_id(1) - 1
    iend   = mm_end
    qmmm%mmatom_id(1:mm_end) = atom_no(1:iend)

    do i = 1, qmmm%qm_natoms - 1
      mm_block = qmmm%qmatom_id(i+1) - qmmm%qmatom_id(i) - 1
      if (mm_block > 0) then
        mm_start = mm_end + 1
        mm_end = mm_start + mm_block - 1
        istart = qmmm%qmatom_id(i) + 1
        iend   = qmmm%qmatom_id(i+1) - 1
        qmmm%mmatom_id(mm_start:mm_end) = atom_no(istart:iend)
      end if
    end do

    mm_block = molecule%num_atoms - qmmm%qmatom_id(qmmm%qm_natoms) - 1
    if (mm_block > 0) then
      mm_start = mm_end + 1
      mm_end = mm_start + mm_block
      istart = qmmm%qmatom_id(qmmm%qm_natoms) + 1
      iend   = molecule%num_atoms
      qmmm%mmatom_id(mm_start:mm_end) = atom_no(istart:iend)    
    end if

    deallocate(atom_no, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc

    ! number of QM-MM bonds
    !
    call count_qmmmbonds(molecule, qmmm)

    if (qmmm%num_qmmmbonds /= 0) then

      ! generate link atom coordinates
      ! remove MM charges
      allocate(qmmm%linkatom_coord(3, qmmm%num_qmmmbonds), stat = kalloc_stat)  
      if (kalloc_stat /= 0) call error_msg_alloc
      allocate(qmmm%linkatom_force(3, qmmm%num_qmmmbonds), stat = kalloc_stat)  
      if (kalloc_stat /= 0) call error_msg_alloc
      allocate(qmmm%qmmmbond_list(2, qmmm%num_qmmmbonds), stat = kalloc_stat)
      if (kalloc_stat /= 0) call error_msg_alloc
      allocate(qmmm%linkatom_charge(qmmm%num_qmmmbonds), stat = kalloc_stat)  
      qmmm%linkatom_charge(1:qmmm%num_qmmmbonds) = 0.0_wp
      if (kalloc_stat /= 0) call error_msg_alloc
      allocate(qmmm%linkatom_global_address(qmmm%num_qmmmbonds), stat = kalloc_stat)  
      if (kalloc_stat /= 0) call error_msg_alloc

      call generate_linkatoms(molecule, psf, qmmm)

      call check_linkmm(molecule, qmmm, setup_qmmm_error)

    else
      qmmm%ec_natoms = 0

    end if

    ! setup list of MM atoms for QM-MM interaction
    if (qmmm%mm_cutoffdist > 0.0_wp) &
      call setup_mmlist(molecule, qmmm)

    ! rearrange the force field terms

    ! bonds
    call rearrange_term_list(2,                              &
             qmmm%qm_natoms, qmmm%qmatom_id, qmmm%qm_nbonds, &
             molecule%num_bonds, forcefield, molecule%bond_list)

    ! angles
    if (molecule%num_angles > 0) then
      call rearrange_term_list(3,                                &
             qmmm%qm_natoms, qmmm%qmatom_id, qmmm%qm_nangles,    &
             molecule%num_angles, forcefield, molecule%angl_list)
    end if

    ! dihedrals
    if (molecule%num_dihedrals > 0) then
      call rearrange_term_list(4,                                &
             qmmm%qm_natoms, qmmm%qmatom_id, qmmm%qm_ndihedrals, &
             molecule%num_dihedrals, forcefield, molecule%dihe_list)

    end if

    ! impropers
    if (molecule%num_impropers > 0) then
      call rearrange_impr_list(                                  &
             qmmm%qm_natoms, qmmm%qmatom_id, qmmm%qm_nimpropers, &
             molecule%num_impropers, forcefield, molecule%impr_list)
    end if

    ! cmaps
    if (molecule%num_cmaps > 0) then
      call rearrange_term_list(8,                                &
             qmmm%qm_natoms, qmmm%qmatom_id, qmmm%qm_ncmaps,     &
             molecule%num_cmaps, forcefield, molecule%cmap_list)
    end if

    ! set the charges of QM region to zero
    do i = 1, qmmm%qm_natoms
      molecule%charge(qmmm%qmatom_id(i)) = 0.0D+00
    end do

    if (setup_qmmm_error) &
      call error_msg('Stop with error while setting up QMMM.')

    ! Check template for QM input
    call qm_check_input(qmmm, setup_qmmm_error)
    if (setup_qmmm_error) &
      call error_msg('Setup_QMMM> Error while checking '//trim(qmmm%qmcnt)//'.')

#ifdef HAVE_MPI_GENESIS
    ! Synchronize MPI processes here, so as to avoid non-main rank processes to 
    ! start QM jobs when error.
    call mpi_barrier(mpi_comm_world, i)
#endif

    return

  end subroutine setup_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_qmmm_directory
  !> @brief        setup directory for QM/MM information
  !> @authors      SI
  !! @param[in]    qmmm_info : QMMM section control parameters information
  !! @param[in]    sel_info  : SELECTION section control parameters information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    psf       : PSF information
  !! @param[in]    rst       : rst information
  !! @param[in]    molecule  : molecular information
  !! @param[out]   qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_qmmm_directory(qmmm_info, qmmm)

    ! formal arguments
    type(s_qmmm_info),       intent(in)    :: qmmm_info
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    character                :: num*5
    character(MaxFilename)   :: path


    ! directory and file names of QM jobs
    write(num,'(i0)') my_replica_no
    qmmm%workdir = trim(qmmm_info%workdir)//'.'//num

    if (replica_main_rank) &
      call system('mkdir -p '//trim(qmmm%workdir)//' > /dev/null 2>&1')

    if (qmmm_info%savedir .ne. '') then
      qmmm%savedir  = trim(qmmm_info%savedir)//'.'//num
      if (replica_main_rank) &
        call system('mkdir -p '//trim(qmmm%savedir)//' > /dev/null 2>&1')
    else
      qmmm%savedir  = qmmm%workdir
    end if
    qmmm%qmbase0    = qmmm_info%basename

    ! qmexe with an absolute path
    if (qmmm_info%qmexe(1:1) .eq. "/") then 
      qmmm%qmexe = trim(qmmm_info%qmexe)

    else if (qmmm_info%qmexe(1:2) .eq. "~/") then
      call getenv('HOME', path)
      qmmm%qmexe = trim(path)//trim(qmmm_info%qmexe(2:))

    else
      call getcwd(path)
      qmmm%qmexe = trim(path)//'/'//trim(qmmm_info%qmexe)
    end if

    return

  end subroutine setup_qmmm_directory

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input
  !> @brief        Check a template file for QM input
  !> @authors      KY
  !! @param[out]   qmmm             : QM/MM information
  !! @param[out]   setup_qmmm_error : returns true when error
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input(qmmm, setup_qmmm_error)

#ifdef QSIMULATE
    use json_module
#endif

    type(s_qmmm),  intent(inout) :: qmmm
    logical,       intent(inout) :: setup_qmmm_error
#ifdef QSIMULATE
    character(kind=json_CK,len=:), allocatable :: input
    integer                      :: csize
#endif

    ! Check template for QM input
    if (main_rank) then
      select case (qmmm%qmtyp)
        case (QMtypQCHEM)
          call qm_check_input_qchem(qmmm, setup_qmmm_error)
        case (QMtypG09)
          call qm_check_input_gaussian(qmmm, setup_qmmm_error)
        case (QMtypTERACHEM)
          call qm_check_input_terachem(qmmm, setup_qmmm_error)
        case (QMtypDFTBPLUS)
          call qm_check_input_dftbplus(qmmm, setup_qmmm_error)
        case (QMtypMOLPRO)
          call qm_check_input_molpro(qmmm, setup_qmmm_error)
        case (QMtypG09_FRWF)
          call qm_check_input_gaussian_frwf(qmmm, setup_qmmm_error)
        case (QMtypORCA)
          call qm_check_input_orca(qmmm, setup_qmmm_error)
      end select
    end if

#ifdef HAVE_MPI_GENESIS
    ! broadcast QM information to other processes
    if (qmmm%qmtyp == QMtypTERACHEM) then
      call mpi_bcast(qmmm%tcscr0, MaxLine, mpi_character, 0, mpi_comm_world, ierror)
    end if

#ifdef QSIMULATE
    if (qmmm%qmtyp == QMtypQSimulate) then
      call qm_check_input_qsimulate(qmmm, setup_qmmm_error)
    end if
#endif
#endif
    
    return

  end subroutine qm_check_input

#ifdef QSIMULATE
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_fill_input_qsimulate
  !> @brief        from the qmmm%qs_input previously read in, it adds other 
  !                relevant parameters for the qmmm input such as parameters for
  !                  - Changes optimization parameters for QM when it is not the first attempt
  !                  - GBSA
  !! @authors      Klaas Gunst
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] inputjson   : copy of the qs_input with any other relevant parameters
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_fill_input_qsimulate(qmmm, enefunc, json, inputjson)
    use json_module

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    type(s_enefunc), intent(in) :: enefunc
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: inputjson
    logical                                  :: found
    ! logical, intent(inout)      :: setup_error

    call json%add_by_path(inputjson, 'silent', .true.)
    call json%add_by_path(inputjson, 'bagel(1).angstrom', .true.)
    if (qmmm%savefile) then
      call json%add_by_path(inputjson, 'file', trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
    end if

    if (qmmm%qmtrial /= 0) call qm_fill_input_qsimulate_retry(qmmm, json, inputjson)
    if (enefunc%gbsa_use) call qm_fill_input_qsimulate_gbsa(qmmm, enefunc, json, inputjson)

    return

  end subroutine qm_fill_input_qsimulate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_fill_input_qsimulate_retry
  !> @brief        Changes optimization parameters for QM when it is not the first attempt
  !! @authors      Klaas Gunst
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] inputjson   : copy of the qs_input with any other relevant parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_fill_input_qsimulate_retry(qmmm, json, inputjson)
    use json_module

    ! formal arguments
    type(s_qmmm), intent(inout)    :: qmmm
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: inputjson

    call json%traverse(inputjson, set_retry_callback)
    return

    contains

  end subroutine qm_fill_input_qsimulate_retry

  subroutine set_retry_callback(json, p, finished)
    use json_module

    class(json_core), intent(inout)       :: json
    type(json_value), pointer, intent(in) :: p
    logical, intent(out)                  :: finished

    character(len=:), allocatable         :: title_name
    logical                               :: found
    finished = .false.

    call json%get(p, 'title', title_name, found)
    if (.not. found) return

    if (title_name .eq. 'xtb' .or. title_name .eq. 'dftb') then
      finished = .true.
      call json%update(p, 'broyden_alpha', 0.2_dp, found)
    endif
    if (&
      title_name .eq. 'hf' .or. title_name .eq. 'ks' .or. &
      title_name .eq. 'rhf' .or. title_name .eq. 'rks' .or. &
      title_name .eq. 'uhf' .or. title_name .eq. 'uks' .or. &
      title_name .eq. 'cis' .or. title_name .eq. 'mp2' &
      ) then
      finished = .true.
      call json%update(p, 'qc', .true., found)
    endif

    deallocate(title_name)
  end subroutine set_retry_callback


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_fill_input_qsimulate_gbsa
  !> @brief        Changes optimization parameters for gbsa
  !! @authors      Klaas Gunst
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] inputjson   : copy of the qs_input with any other relevant parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_fill_input_qsimulate_gbsa(qmmm, enefunc, json, inputjson)
    use json_module

    ! formal arguments
    type(s_qmmm), intent(inout)    :: qmmm
    type(s_enefunc), intent(in)    :: enefunc
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: inputjson


    character(len=:), allocatable  :: basisname
    logical                        :: found

    call json%get(inputjson, 'bagel(1).basis', basisname, found)
    if (.not. found) then
      call error_msg( &
        'INPUT_QSIMULATE> "basis" should be present in first block of BAGEL json input')
    endif
    if (basisname .ne. 'xtb' .and. basisname .ne. 'dftb') then
      call error_msg( &
        'INPUT_QSIMULATE> QMMM-GBSA works only for xtb or dftb calculations.')
    endif
    deallocate(basisname)

    ! Absurd workaround but it seems that for gfortran 10.2.1 with some optimization settings
    ! including set_gbsa_callback as `contains` gives rise to segmentation faults when called.
    ! So need to do this through private variables of the module
    salt_cons = enefunc%gbsa%salt_cons
    temperature = enefunc%gbsa%temperature
    eps_solvent = enefunc%gbsa%eps_solvent
    eps_solute = enefunc%gbsa%eps_solute
    call json%traverse(inputjson, set_gbsa_callback)

  end subroutine qm_fill_input_qsimulate_gbsa

  subroutine set_gbsa_callback(json, p, finished)
    use json_module

    class(json_core), intent(inout)       :: json
    type(json_value), pointer, intent(in) :: p
    logical, intent(out)                  :: finished

    character(len=:), allocatable         :: title_name
    logical                               :: found
    finished = .false.

    call json%get(p, 'title', title_name, found)
    if (.not. found) return

    if (title_name .eq. 'xtb' .or. title_name .eq. 'dftb') then
      finished = .true.
      call json%update(p, 'solvation', .true., found)
      call json%update(p, 'gbsa_for_qmmm', .true., found)
      call json%update(p, 'salt_cons', salt_cons, found)
      call json%update(p, 'gbsa_temperature', temperature, found)
      call json%update(p, 'epsilon_solvent', eps_solvent, found)
      call json%update(p, 'epsilon_solute', eps_solute, found)
    endif
    deallocate(title_name)
  end subroutine set_gbsa_callback
#endif


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_mmlist
  !> @brief        setup for MM list
  !> @authors      KY
  !! @param[in]    molecule  : molecular information
  !! @param[out]   qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_mmlist(molecule, qmmm)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    integer :: i, id, nres, natm, max_natm, res_no1, res_no2


    qmmm%mm_all_natoms = qmmm%mm_natoms
    allocate(qmmm%mm_natoms_res(molecule%num_residues))

    max_natm = 0
    nres     = 1
    natm     = 1
    id       = qmmm%mmatom_id(1)
    res_no1  = molecule%residue_no(id)

    do i = 2, qmmm%mm_natoms
      id      = qmmm%mmatom_id(i)
      res_no2 = molecule%residue_no(id)
      if (res_no2 == res_no1) then
        natm = natm + 1
      else
        qmmm%mm_natoms_res(nres) = natm
        nres = nres + 1
        if (natm > max_natm) max_natm = natm

        natm = 1
        res_no1 = res_no2
      end if

    end do
    qmmm%mm_natoms_res(nres) = natm
    qmmm%mm_nres = nres

    !dbg write(MsgOut,'(i5)') max_natm
    !dbg write(MsgOut,'(i5)') qmmm%mm_nres
    !dbg do i = 1, qmmm%mm_nres
    !dbg   write(MsgOut,'("resid=",i5," natom=",i5)') i,qmmm%mm_natoms_res(i)
    !dbg end do

    allocate(qmmm%mmatom_id_res(max_natm, qmmm%mm_nres))

    i = 1
    do nres = 1, qmmm%mm_nres
      do natm = 1, qmmm%mm_natoms_res(nres)
        qmmm%mmatom_id_res(natm,nres) = qmmm%mmatom_id(i)
        i = i + 1
      end do
    end do

    !dbg do nres = 1, qmmm%mm_nres
    !dbg   write(MsgOut,'("resid=",i5)') nres
    !dbg   natm = qmmm%mm_natoms_res(nres)
    !dbg   write(MsgOut,'(5x,10i5)') qmmm%mmatom_id_res(1:natm,nres)
    !dbg end do

    return

  end subroutine setup_mmlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_qsimulate
  !> @brief        it actually read the template and store in qmmm%qm_template
  !! @authors      Toru Shiozaki, Kiyoshi Yagi
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_qsimulate(qmmm, setup_error)

#ifdef QSIMULATE
    use iso_c_binding
    use json_module
#endif
    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error

#ifdef QSIMULATE
    character(kind=json_CK,len=:), allocatable :: input
    integer                      :: csize

    if (main_rank) then
      write(MsgOut,*)
      write(MsgOut,'(a)') '  Check the input file for QSimulate [ ' &
                           //trim(qmmm%qmcnt)//' ]'

      call qmmm%qs_input%initialize()
      call qmmm%qs_input%load(filename = trim(qmmm%qmcnt))
      call qmmm%qs_input%serialize(input)
      csize = len(input)
    endif
    call mpi_bcast(csize, 1, mpi_integer, 0, mpi_comm_world, ierror)

    if (.not. main_rank) allocate(character(csize)::input)
    call mpi_bcast(input, csize, mpi_character, 0, mpi_comm_world, ierror)

    if (.not. main_rank) call qmmm%qs_input%deserialize(input)

    ! TODO : perform basic check of the input parameters
#endif

    return

  end subroutine qm_check_input_qsimulate

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_qchem
  !> @brief        check template for qchem input
  !! @authors      KYMD
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_qchem(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error


    return

  end subroutine qm_check_input_qchem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_gaussian
  !> @brief        check template for gaussian input
  !! @authors      KY
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_gaussian(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical, intent(inout)      :: setup_error

    ! local variables
    logical                     :: defined_coord, defined_charge
    integer                     :: file1, ierror
    character(256)              :: string


    write(MsgOut,*) 
    write(MsgOut,'(a)') '  Check the control file for Gaussian [ ' &
                         //trim(qmmm%qmcnt)//' ]'

    defined_coord   = .false.
    defined_charge  = .false.

    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    do
      read(file1, '(a)', iostat = ierror) string
      if (ierror /= 0) exit
      call tolower(string)

      if (index(string, '#coord') > 0) defined_coord   = .true.
      if (index(string, '#charg') > 0) defined_charge  = .true.
    end do

    if (.not. defined_coord .or. .not. defined_charge) then
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      if (.not. defined_coord) &
        write(MsgOut,'("ERROR>>> #coordinate# is not found!")')
      if (.not. defined_charge) &
        write(MsgOut,'("ERROR>>> #charge# is not found!")')

      setup_error = .true.

    else
      write(MsgOut,'(a)') '  Passed the check!'

    end if

    call close_file(file1)

    return

  end subroutine qm_check_input_gaussian

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_terachem
  !> @brief        check template for terachem input
  !! @authors      KY
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_terachem(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error

    ! local variables
    integer                  :: file1
    character(1024)          :: string1, string2
    character(256)           :: key, value, guess
    logical                  :: scratch, amber, atomname
    character(10)            :: linkHname
    integer                  :: i, num_qmatom


    write(MsgOut,*) 
    write(MsgOut,'(a)') '  Check the control file for TeraChem [ ' &
                         //trim(qmmm%qmcnt)//' ]'
    amber = .true.
    scratch = .false.
    linkHname = ''
    atomname  = .false.
    num_qmatom = 0

    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    do while (.true.)

      read(file1, '(a)', end=100) string1
      i = index(string1,'#')
      if (i > 0) string1 = string1(1:i)
      string2 = string1
      call tolower(string1)

      !if (index(string1,'end') > 0) exit

      if (index(string1,'scrdir') > 0) then
        scratch = .true.
        read(string2,*) key, value
        qmmm%tcscr0 = trim(value)

      else if (index(string1,'amber') > 0) then
        ! check amber = yes
        read(string2,*) key, value
        if (trim(value) == 'no') amber = .false.

      else if (index(string1,'guess') > 0) then
        i = index(string1, 'guess')
        value = adjustl(string1(i+5:))
        i = index(value,'/', back=.true.)
        guess = value(1:i-1)

      else if (index(string1,'linkhname') >0) then
        read(string2,*) key, linkHname

      else if (index(string1,'atomname') >0) then
        atomname = .true.
        read(string2,*) key, num_qmatom

      end if

    end do
100 continue

    call close_file(file1)

    if (.not. scratch) then
      ! default scratch dir of TeraChem
      qmmm%tcscr0 = 'scr'
    end if
    write(MsgOut,'(a)') '    scratch directory  = '//trim(qmmm%tcscr0)

!    if (trim(guess) /= trim(qmmm%tcscr)) then
!      write(MsgOut,'(a)') '  !! WARNING !!'
!      write(MsgOut,'(a)') '  Guess MO directory is not the same as the scratch directory.'
!      write(MsgOut,'(a)') '    o guess dir   = '//trim(guess)
!      write(MsgOut,'(a)') '    o scratch dir = '//trim(qmmm%tcscr)
!    end if

    if (len(trim(linkHname)) /= 0) then
      if (atomname) then
        write(MsgOut,*) 
        write(MsgOut,'(a)') 'WARNING>>> linkHname and atomname are both present in ' &
                            //trim(qmmm%qmcnt)
        write(MsgOut,'(a)') 'WARNING>>> linkHname is ignored.'
        write(MsgOut,*) 

      else if (linkHname(1:1) .ne. 'H') then
        write(MsgOut,*) 
        write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
        write(MsgOut,'(a)') 'ERROR>>> linkHname must start with "H"!'
        write(MsgOut,*) 
        setup_error = .true.

      else
        write(MsgOut,'(a)') '    link hydrogen name = '//trim(linkHname)

      end if

    else
      write(MsgOut,'(a)') '    link hydrogen name = H'

    end if

    if (.not. amber) then
      write(MsgOut,*) 
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      write(MsgOut,'(a)') 'ERROR>>> amber must be yes in TeraChem input!'
      setup_error = .true.
    end if

    if (atomname .and. &
        num_qmatom /= qmmm%qm_natoms + qmmm%num_qmmmbonds) then
      write(MsgOut,*) 
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      write(MsgOut,'(a)') 'ERROR>>> Number of atomname does not match the number',&
                          ' of QM atoms + link hydrogen atoms!'
      write(MsgOut,'(a,i0)') 'ERROR>>> atomname      ',num_qmatom
      write(MsgOut,'(a,i0)') 'ERROR>>> QM atom       ',qmmm%qm_natoms
      write(MsgOut,'(a,i0)') 'ERROR>>> link hydrogen ',qmmm%num_qmmmbonds
      setup_error = .true.
    end if

    if (.not. setup_error) then
      write(MsgOut,'(a)') '  Passed the check!'
    end if

    write(MsgOut,*) 
    
    return

  end subroutine qm_check_input_terachem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_dftb+
  !> @brief        check template for dftb+ input
  !! @authors      KY
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_dftbplus(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical, intent(inout)      :: setup_error

    ! local variables
    integer                  :: file1, i
    character(1024)          :: string1, string2
    logical                  :: geom, efield


    write(MsgOut,*) 
    write(MsgOut,'(a)') '  Check the control file for DFTB+ [ ' &
                         //trim(qmmm%qmcnt)//' ]'

    geom   = .false.
    efield = .false.
    qmmm%charges_bin = .true.

    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    do while (.true.)

      read(file1, '(a)', end=100) string1
      i = index(string1,'#')
      if (i > 0) string1 = string1(1:i)
      string2 = string1
      call tolower(string1)

      if (index(string1,'geometry') > 0) then
        geom = .true.

      else if (index(string1,'electricfield') > 0) then
        efield = .true.

      else if (index(string1,'writechargesastext') > 0 &
               .and. index(string1,'yes') > 0) then
        qmmm%charges_bin = .false.

      end if

    end do
100 continue

    call close_file(file1)

    if (geom) then
      write(MsgOut,*) 
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      write(MsgOut,'(a)') 'ERROR>>> Entry for Geometry exists.'
      setup_error = .true.
    end if

    if (efield) then
      write(MsgOut,*) 
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      write(MsgOut,'(a)') 'ERROR>>> Entry for ElectricField exists.'
      setup_error = .true.
    end if

    if (.not. geom .and. .not. efield) then
      write(MsgOut,'(a)') '  Passed the check!'
    end if

    write(MsgOut,*) 
    
    return

  end subroutine qm_check_input_dftbplus

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_gaussian_frwf
  !> @brief        check template for gaussian input
  !! @authors      KYMD, KY
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_gaussian_frwf(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error

    ! local variables
    integer                     :: file1, ipos, ierror, is_exist
    logical                     :: is_route_section
    logical                     :: defined_gen, defined_iop6, defined_iop4
    logical                     :: defined_coord, defined_charge
    logical                     :: defined_elfield
    character(256)              :: string


    is_route_section  = .false.
    defined_gen       = .false.
    defined_iop6      = .false.
    defined_iop4      = .false.
    defined_coord     = .false.
    defined_charge    = .false.
    defined_elfield   = .false.

    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    do
      read(file1, '(a)', iostat = ierror) string
      if (ierror /= 0) exit
      call tolower(string)

      if (.not.is_route_section .and. (index(string, '#') /= 0)) &
        is_route_section = .true.
      if (is_route_section .and. len_trim(string) == 0) exit
      if (.not.is_route_section) cycle

      ! whether in qmmm%qmcnt, "gen" and its analogous keywords are specified or not 
      !
      ipos = index(string, '!') - 1
      if (ipos == -1) ipos = len_trim(string)
      is_exist = index(string(1: ipos), 'gen') + index(string(1: ipos), 'genecp') + index(string(1: ipos), 'pseudo')
      if (is_exist /= 0) defined_gen = .true.
      is_exist = index(string(1: ipos), '6/13=1')
      if (is_exist /= 0) defined_iop6 = .not.defined_iop6
      is_exist = index(string(1: ipos), '4/5=100')
      if (is_exist /= 0) defined_iop4 = .not.defined_iop4
    end do
    do
      read(file1, '(a)', iostat = ierror) string
      if (ierror /= 0) exit
      call tolower(string)

      ipos = index(string, '!') - 1
      if (ipos == -1) ipos = len_trim(string)

      is_exist = index(string(1: ipos), '#coordinate#')
      if (is_exist /= 0) defined_coord   = .not.defined_coord
      is_exist = index(string(1: ipos), '#charge#')
      if (is_exist /= 0) defined_charge  = .not.defined_charge
      is_exist = index(string(1: ipos), '#elec_field#')
      if (is_exist /= 0) defined_elfield = .not.defined_elfield
      if (defined_coord .and. (.not.defined_gen .or. defined_elfield)) exit
    end do

    call close_file(file1)

    ! Flag to save electric field on MM atoms, arising from QM density, in read-write file
    ! ... In Single-Point calculation, qmmm%ene_only has to be YES ...
    !if (.not.defined_iop6 .and. .not.qmmm%ene_only)                           &
    ! call error_msg('Setup_QMMM> iop(6/13=1) must be given in your '//trim(qmmm%qmcnt))

    ! Flag to read initial guess from previous results if possible
    ! 
    if (.not.defined_iop4)                                                    &
     write(MsgOut, '(/a/)') "NOTICE: It is a good idea to add iop(4/5=100) to your "//trim(qmmm%qmcnt)

    ! Sign to place QM and MM charges coordinates
    !
    !if (.not.defined_coord)                                                   &
    ! call error_msg('Setup_QMMM> #coordinate# must be given in your '//trim(qmmm%qmcnt))

    ! Sign to place MM charges coordinates for electric field
    !
    !if (defined_gen .and. .not.defined_elfield)                               &
    ! call error_msg('Setup_QMMM> #elec_field# must be given in your '//trim(qmmm%qmcnt))

    if (.not. defined_coord .or. .not. defined_charge .or. &
        .not. defined_coord .or. .not. (defined_iop6 .and. &
        .not.qmmm%ene_only)) then
      write(MsgOut,'(a)') 'ERROR>>> Fatal error in '//trim(qmmm%qmcnt)
      if (.not.defined_iop6 .and. .not.qmmm%ene_only) &
        write(MsgOut,'("ERROR>>> iop(6/13=1) is not found!")')
      if (.not. defined_coord) &
        write(MsgOut,'("ERROR>>> #coordinate# is not found!")')
      if (.not. defined_charge) &
        write(MsgOut,'("ERROR>>> #charge# is not found!")')
      if (.not. defined_coord) &
        write(MsgOut,'("ERROR>>> #elec_field# is not found!")')

      setup_error = .true.

    end if
 
    qmmm%defined_gen = defined_gen

    return

  end subroutine qm_check_input_gaussian_frwf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_molpro
  !> @brief        check template for molpro input
  !! @authors      YA
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_molpro(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error


    return

  end subroutine qm_check_input_molpro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_check_input_orca
  !> @brief        check template for orca input
  !! @authors      
  !! @param[in]    qmmm        : information on QM/MM setting
  !! @param[inout] setup_error : return true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_check_input_orca(qmmm, setup_error)

    ! formal arguments
    type(s_qmmm), intent(inout) :: qmmm
    logical,      intent(inout) :: setup_error


    return

  end subroutine qm_check_input_orca

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_atomic_number
  !> @brief        get atomic number
  !! @authors      YA, KY
  !! @param[in]    mm_mass   : atom mass
  !! @param[in]    atomic_no : atomic number
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_atomic_number(mm_mass, atomic_no)

    ! formal arguments
    real(wp),                intent(in)  :: mm_mass
    integer,                 intent(out) :: atomic_no

    ! local variables
    integer                  :: ielem


    do ielem = 1, 118
      if (mm_mass >= qm_atomic_mass(ielem) - mass_error .and. &
        mm_mass <= qm_atomic_mass(ielem) + mass_error) then
        atomic_no = ielem
        return
      end if
    end do

    if (mm_mass >= deuterium_mass - mass_error .and. &
      mm_mass <= deuterium_mass + mass_error) then
      atomic_no = 1
      return
    end if

    atomic_no = -1

    return

  end subroutine qm_atomic_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input
  !> @brief        generate QM input file
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of whole molecule
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    character(4)  :: string
    integer       :: i, qm_atom, mm_atom
    real(wp)      :: rr01, r0(3), r1(3), r01(3)


    ! generate link atom coordinates
    if (qmmm%num_qmmmbonds > 0) then
      do i = 1, qmmm%num_qmmmbonds

        qm_atom = qmmm%qmmmbond_list(1,i) 
        mm_atom = qmmm%qmmmbond_list(2,i)

        r0(1:3) = coord(1:3,qm_atom)
        r01(1:3) = coord(1:3,mm_atom) - r0(1:3)

        rr01 = sqrt(r01(1) * r01(1) + r01(2) * r01(2) + r01(3) * r01(3))
        qmmm%linkatom_coord(1:3,i) = &
          r0(1:3) + r01(1:3) * (qmmm%linkbond / rr01)
        qmmm%linkatom_global_address(i) = mm_atom

      end do
    end if

    ! update list of MM atoms
    call update_mmlist(coord, qmmm)

    if (replica_main_rank) then
      select case (qmmm%qmtyp)
        case (QMtypQCHEM)
          call qm_generate_input_qchem(molecule, coord, qmmm)
        case (QMtypG09)
          call qm_generate_input_gaussian(molecule, coord, qmmm)
        case (QMtypTERACHEM)
          call qm_generate_input_terachem(molecule, coord, qmmm)
        case (QMtypDFTBPLUS)
          call qm_generate_input_dftbplus(molecule, coord, qmmm)
        case (QMtypMOLPRO)
          call qm_generate_input_molpro(molecule, coord, qmmm)
        case (QMtypG09_FRWF)
          call qm_generate_input_gaussian_frwf(molecule, coord, qmmm)
        case (QMtypORCA)
          call qm_generate_input_orca(molecule, coord, qmmm)
      end select
    end if

    return

  end subroutine qm_generate_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output
  !> @brief        Retrieve energy and gradients from QM output files
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] energy    : energy information
  !! @param[inout] qm_error  : error flag
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output(molecule, qmmm, energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    type(s_energy),          intent(inout) :: energy
    logical,                 intent(inout) :: qm_error


    ! local variables
    real(wp)            :: diff
    real(wp)            :: rab, vab(3), rab1, rab3, ff, lafn(3), la_chg
    integer             :: i, j, n, address, qm_atom, mm_atom
    integer             :: errorcode
    integer(4)          :: istat, unlink
    character(len_exec) :: exec

    real(wp),   pointer :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp),   pointer :: qm_charge(:), la_charge(:)
    real(wp),   pointer :: qm_dipole(:)


    ! use pointers
    !
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge
    qm_dipole => qmmm%qm_dipole
 
    ! Read QM energy (hartree) and gradient (hartree/bohr) 
    if (replica_main_rank) then
      select case (qmmm%qmtyp)
        case (QMtypQCHEM)
          call qm_read_output_qchem(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypG09)
          call qm_read_output_gaussian(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypTERACHEM)
          call qm_read_output_terachem(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypDFTBPLUS)
          call qm_read_output_dftbplus(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypMOLPRO)
          call qm_read_output_molpro(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypG09_FRWF)
          call qm_read_output_gaussian_frwf(molecule, qmmm, energy%qm_ene, qm_error)
        case (QMtypORCA)
          call qm_read_output_orca(molecule, qmmm, energy%qm_ene, qm_error)
      end select

      if (.not. qmmm%savefile .and. .not. qm_error) then
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') &
                       trim(qmmm%workdir)//'/'//trim(qmmm%qminp)
        istat = unlink(trim(qmmm%workdir)//'/'//trim(qmmm%qminp))
        if (istat /= 0) then
          call error_msg ('Compute_Energy_QMMM> Error while removing '// &
                           trim(qmmm%workdir)//'/'//trim(qmmm%qminp))
        end if

        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') &
                       trim(qmmm%workdir)//'/'//trim(qmmm%qmout)
        istat = unlink(trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
        if (istat /= 0) then
          call error_msg ('Compute_Energy_QMMM> Error while removing '// &
                           trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
        end if

      else if (qmmm%savedir /= qmmm%workdir) then
        exec = 'mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qminp)// &
               '   '//trim(qmmm%workdir)//'/'//trim(qmmm%qmout)//' '//trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

      end if

      if (qmmm%savefile .and. qmmm%save_qminfo) call write_qminfo(qmmm, energy%qm_ene)

    end if

#ifdef HAVE_MPI_GENESIS
    ! broadcast QM information to other processes
    call mpi_bcast(qm_error , 1, mpi_logical, 0, mpi_comm_country, ierror)
    call mpi_bcast(energy%qm_ene, 1, mpi_wp_real, 0, mpi_comm_country, ierror)
    call mpi_bcast(qm_force(1,1), qmmm%qm_natoms*3, mpi_wp_real, 0,     &
                   mpi_comm_country, ierror)
    call mpi_bcast(mm_force(1,1), qmmm%mm_natoms*3, mpi_wp_real, 0,     &
                   mpi_comm_country, ierror)
    call mpi_bcast(qm_charge(1) , qmmm%qm_natoms, mpi_wp_real, 0,       &
                   mpi_comm_country, ierror)
    call mpi_bcast(qmmm%is_qm_charge, 1, mpi_logical, 0,                &
                   mpi_comm_country, ierror)
    call mpi_bcast(qm_dipole(1), 3, mpi_wp_real, 0, mpi_comm_country, ierror)

    if (qmmm%num_qmmmbonds > 0) then
      call mpi_bcast(la_force(1,1), qmmm%num_qmmmbonds*3, mpi_wp_real, 0, &
                     mpi_comm_country, ierror)
      call mpi_bcast(la_charge(1) , qmmm%num_qmmmbonds, mpi_wp_real, 0,   &
                     mpi_comm_country, ierror)
    end if
#endif

    return

  end subroutine qm_read_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_post_process
  !> @brief        Post process after QM info is obtained
  !! @authors      YA, KY
  !! @param[in]    coord     : coordinates
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] energy    : energy information
  !! @param[inout] force     : atomic forces
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qm_post_process(coord, qmmm, energy, force)

    ! formal arguments
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm), target,    intent(inout) :: qmmm
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: force(:,:)

    ! local variables
    real(wp)          :: rab, vab(3), rab1, rab3, ff, lafn(3), la_chg
    real(wp)          :: diff
    integer           :: i, j, n, address, qm_atom, mm_atom
    integer           :: errorcode
    character(len_exec) :: exec

    real(wp), pointer :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer :: qm_charge(:), la_charge(:)


    ! use pointers
    !
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    ! ------------------------------------
    ! QM energy / gradient obtained safely

    energy%total  = energy%total + energy%qm_ene * CONV_UNIT_ENE

    if (.not. qmmm%ene_only) then
      do i = 1, qmmm%mm_natoms
        address = qmmm%mmatom_id(i)
        force(1:3,address) = force(1:3,address) &
                           + mm_force(1:3,i) * CONV_UNIT_FORCE
      end do

      do i = 1, qmmm%qm_natoms
        address = qmmm%qmatom_id(i)
        force(1:3,address) = force(1:3,address) &
                           + qm_force(1:3,i) * CONV_UNIT_FORCE
      end do

      if (qmmm%num_qmmmbonds /= 0) then
        !
        ! gradient correction due to linkatoms
        do n = 1, qmmm%num_qmmmbonds
          lafn    = la_force(1:3,n) * CONV_UNIT_FORCE
          qm_atom = qmmm%qmmmbond_list(1,n)
          mm_atom = qmmm%qmmmbond_list(2,n)

          vab(1:3) = coord(1:3,mm_atom) - coord(1:3,qm_atom)

          rab = 0.0_wp
          do i = 1, 3
            rab = rab +  vab(i)*vab(i)
          end do
          rab = sqrt(rab)
          rab1 = 1.0_wp / rab
          rab3 = 1.0_wp / rab / rab / rab

          do i = 1, 3

            ff = 0.0_wp
            do j = 1, 3
               ff = ff + vab(i)*vab(j)*lafn(j)
            end do
            ff = ff * qmmm%linkbond * rab3

            ! gradient correction to QM atom
            force(i,qm_atom) = force(i,qm_atom)               &
              + (1.0_wp - qmmm%linkbond * rab1) * lafn(i)     &
              + ff

            ! gradient correction to MM atom
            force(i,mm_atom) = force(i,mm_atom)               &
              + qmmm%linkbond * rab1 * lafn(i)                &
              - ff

          end do

        end do

        ! link-atom charge is distributed among rest QM atoms
        do n = 1, qmmm%num_qmmmbonds
          la_chg = la_charge(n) / real(qmmm%qm_natoms, wp)
          do i = 1, qmmm%qm_natoms
            qm_charge(i) = qm_charge(i) + la_chg
          end do
        end do

      end if

    end if

    return

  end subroutine qm_post_process

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_qminfo
  !> @brief        write information of QM calc. to file
  !! @authors      KY
  !! @param[in]    qmmm
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine write_qminfo(qmmm, qm_energy)

    ! formal arguments
    type(s_qmmm),            intent(in) :: qmmm
    real(wp),                intent(in) :: qm_energy

    ! local variables
    integer                  :: file1, i
    character(256)           :: qminfo


    if (.not. replica_main_rank) return

    qminfo = trim(qmmm%savedir)//'/'//trim(qmmm%qminfo)
    call open_file(file1, qminfo, IOFileOutputReplace)

    write(file1, '("QM energy")')
    write(file1, '(e30.20)') qm_energy
    write(file1, '("QM dipole")')
    write(file1, '(3e30.20)') qmmm%qm_dipole

    if (qmmm%is_qm_charge) then
      write(file1, '("QM charge")')
      write(file1, '(i8)') qmmm%qm_natoms
      do i = 1, qmmm%qm_natoms
        write(file1, '(i8,e30.20)') qmmm%qmatom_id(i), qmmm%qm_charge(i)
      end do
      if (qmmm%num_qmmmbonds /= 0) then
        write(file1, '("Link atom charge")')
        write(file1, '(i8)') qmmm%num_qmmmbonds
        do i = 1, qmmm%num_qmmmbonds
          write(file1, '(2i8,e30.20)') &
            qmmm%qmmmbond_list(:,i), qmmm%linkatom_charge(i)
        end do
      end if
    end if

    if (.not. qmmm%ene_only) then
      write(file1, '("QM gradient")')
      write(file1, '(i8)') qmmm%qm_natoms
      do i = 1, qmmm%qm_natoms
        write(file1, '(i8,3e30.20)') qmmm%qmatom_id(i), qmmm%qm_force(:,i)
      end do

      if (qmmm%num_qmmmbonds /= 0) then
        write(file1, '("Link atom gradient")')
        write(file1, '(i8)') qmmm%num_qmmmbonds
        do i = 1, qmmm%num_qmmmbonds
          write(file1, '(2i8,3e30.20)') &
            qmmm%qmmmbond_list(:,i), qmmm%linkatom_force(:,i)
        end do
      end if

      write(file1, '("MM gradient")')
      write(file1, '(i8)') qmmm%mm_natoms
      do i = 1, qmmm%mm_natoms
        write(file1, '(i8,3e30.20)') qmmm%mmatom_id(i), qmmm%mm_force(:,i)
      end do
    end if

    call close_file(file1)

    return

  end subroutine write_qminfo

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_qm_charge
  !> @brief        print QM charge to output
  !! @authors      KY
  !! @param[in]    molecule  : molecular information
  !! @param[in]    qmmm      : QM/MM information
  !! @param[in]    output   : output information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine output_qm_charge(molecule, qmmm, output)

    ! formal arguments
    type(s_molecule),         intent(in) :: molecule
    type(s_qmmm),             intent(in) :: qmmm
    type(s_output), optional, intent(in) :: output

    ! local variables
    integer :: i, iatom, iout


    if (present(output)) then
      ! output
      if (.not. replica_main_rank) return
      if (output%logout) then
        iout = output%logunit
      else
        iout = MsgOut
      end if

    else
      ! setup
      if (.not. main_rank) return
      iout = MsgOut

    end if

    write(iout,'(4x," QM charge:",17x,"atom info",16x,&
                    "(new)",7x,"(old)",7x,"delta")')
    do i = 1, qmmm%qm_natoms
      iatom   = qmmm%qmatom_id(i)
      write(iout,'(4x," QM charge: ",2i6,x,a4,x,i6,x,a4,x,a6,3f12.4)')  &
        i,                             &
        iatom,                         &
        molecule%segment_name(iatom),  &
        molecule%residue_no(iatom),    &
        molecule%residue_name(iatom),  &
        molecule%atom_name(iatom),     &
        qmmm%qm_charge(i),             &
        qmmm%qm_charge_save(i),        &
       (qmmm%qm_charge(i) - qmmm%qm_charge_save(i))
    end do
    write(iout,*)

    return

  end subroutine output_qm_charge

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_monitor
  !> @brief        Monitor QM execution
  !! @authors      YA
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_monitor(qmmm)

    ! formal arguments
    type(s_qmmm),            intent(in) :: qmmm

    ! local variables
    integer,    parameter  :: monitor_interval = 1
    integer,    parameter  :: monitor_max = 3600

    integer                :: i
    character(len_exec)    :: exec
    character(MaxFilename) :: filename
    logical                :: done


    filename = trim(qmmm%workdir)//'/qm_done'
    do i = 1, monitor_max
      inquire(file=trim(filename), exist=done)
      if (done) exit
      call sleep(monitor_interval)
    end do

    if (done) then
      exec = 'rm '//trim(filename)
      call system(exec)
    end if

    return

  end subroutine qm_monitor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qmmm_finalize
  !> @brief        Finalize QM/MM calculation
  !! @authors      KY
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine qmmm_finalize(qmmm)

    ! formal arguments
    type(s_qmmm),  intent(inout) :: qmmm

    ! local variables
    integer            :: ilen
    character(len_exec):: exec


    if (qmmm%workdir /= qmmm%savedir) then
      if (nrep_per_proc == 1) then
        exec = 'rm -rf '//trim(qmmm%workdir)
        call system(exec)
      else
        if (main_rank) then
          ilen = len_trim(qmmm%workdir)
          if (nrep_per_proc < 10) then
            qmmm%workdir = qmmm%workdir(1:ilen-1)//'*'
            write(*,*) "SITO: qmdir=", qmmm%workdir
          else
            qmmm%workdir = qmmm%workdir(1:ilen-2)//'*'
          end if
          exec = 'rm -rf '//trim(qmmm%workdir)
          call system(exec)
        end if
      end if 
    end if

    return

  end subroutine qmmm_finalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_mmlist
  !> @brief        Update list of MM atoms that are used for QM calc.
  !! @authors      KY
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_mmlist(coord, qmmm)

    ! formal arguments
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm), target,    intent(inout) :: qmmm

    ! local variables
    real(wp)       :: cutoffdist2, dij(3), rij2
    integer        :: n, i, j, imm, jqm, ii
    integer        :: id, my_id, nthread, nsize
#ifdef OMP
    integer        :: omp_get_thread_num, omp_get_max_threads
#endif
    logical        :: add_res
    logical        :: byresidue

    integer, allocatable :: mm_natoms(:), mmatom_id_n(:,:)
    integer, pointer     :: mmatom_id(:), qmatom_id(:)
    integer, pointer     :: mm_natoms_res(:), mmatom_id_res(:,:)


    if (qmmm%mm_cutoffdist < 0.0_wp) return

    mmatom_id => qmmm%mmatom_id
    qmatom_id => qmmm%qmatom_id
    mm_natoms_res => qmmm%mm_natoms_res
    mmatom_id_res => qmmm%mmatom_id_res
    cutoffdist2 = qmmm%mm_cutoffdist * qmmm%mm_cutoffdist
    byresidue = qmmm%mm_cutoff_byresidue

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    nsize = int(qmmm%mm_all_natoms/nthread) + 1
    allocate(mm_natoms(nthread), mmatom_id_n(nsize,nthread))

    !$omp parallel &
    !$omp private(id, my_id, n, add_res, i, j, imm, jqm, dij, rij2, ii)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    !for hybrid OMP/MPI parallel
    !my_id = my_city_rank * nthread + id
    my_id = id
    id = id + 1

    mm_natoms(id) = 0
    do n = 1, qmmm%mm_nres
      add_res = .false.

      !for hybrid OMP/MPI parallel
      !if (mod(n-1,nproc_city*nthread) /= my_id) cycle
      if (mod(n-1,nthread) /= my_id) cycle

      do i = 1, mm_natoms_res(n)
        imm = mmatom_id_res(i,n)

        do j = 1, qmmm%qm_natoms
          jqm = qmatom_id(j)

          ! compute distance
          !
          dij(1:3) = coord(1:3,imm) - coord(1:3,jqm)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoffdist2) then
            if (byresidue) then
              do ii = 1, mm_natoms_res(n)
                mm_natoms(id) = mm_natoms(id) + 1
                mmatom_id_n(mm_natoms(id),id) = mmatom_id_res(ii,n)
              end do
              add_res = .true.

            else
              mm_natoms(id) = mm_natoms(id) + 1
              mmatom_id_n(mm_natoms(id),id) = imm

            end if

            ! exit the loop over QM atoms
            exit

          end if

        end do
        if (add_res) exit
        
      end do
    end do
    !$omp end parallel

    qmmm%mm_natoms = 0
    do n = 1, nthread
      mmatom_id(qmmm%mm_natoms+1:qmmm%mm_natoms + mm_natoms(n)) = &
        mmatom_id_n(1:mm_natoms(n),n)
      qmmm%mm_natoms = qmmm%mm_natoms + mm_natoms(n)
    end do

    deallocate(mm_natoms, mmatom_id_n)

    !dbg write(MsgOut,'("mm_natoms = ",i10)')  qmmm%mm_natoms

    return

  end subroutine update_mmlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_qchem
  !> @brief        generate QChem input file
  !! @authors      YA
  !! @param[in]    molecule  : information of whole molecule
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_qchem(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    character(256)           :: qminp
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    character(1024)          :: string, string2
    character(1)             :: conv
    integer                  :: file1, file2
    integer                  :: i, cnt
    integer                  :: address, atomic_no, charge, spin


    ! open files
    !
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    qminp = trim(qmmm%workdir)//'/'//trim(qmmm%qminp)
    call open_file(file2, qminp, IOFileOutputReplace)

    do while (.true.)

      read(file1, '(a)', end=100) string
      string2 = string
      call tolower(string)
      if (string(1:9) .eq. '$molecule') then
        write(file2, '(a)') '$molecule'
        read(file1, '(a)', end=100) string
        write(file2, '(a)') trim(string)

        ! QM atom coordinates
        !
        cnt = 0
        do i = 1, qmmm%qm_natoms
          atomic_no = qmmm%qm_atomic_no(i)
          address = qmmm%qmatom_id(i)
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                        qm_atomic_symbol(atomic_no), coord(1:3,address)
        end do

        ! link atom coordinates
        !
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                       "H  ", qmmm%linkatom_coord(1:3,i)
        end do

        !
        ! Sort link atoms
        call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

        !
        ! Write QM coordinates
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          write(file2, *) trim(coord_txt(i))
        end do

        ! MM point charges
        !
        if (qmmm%mm_natoms > 0) then

          write(file2, '(a4,/)') '$end'
          write(file2, '(a17)') '$external_charges'
          do i = 1, qmmm%mm_natoms
            address = qmmm%mmatom_id(i)
            write(file2, '(3(f18.10,x),f18.10)') &
              coord(1:3,address), molecule%charge(address)
          end do

        end if

      else if (index(string,"scf_convergence") > 0) then

        write(conv,'(i1)') 8-qmmm%qmtrial*2
        string2='scf_convergence    '//conv
        write(file2,'(a)') trim(string2)

      else

        write(file2, '(a)') trim(string2)

      end if
    end do

100 continue

    call close_file(file1)
    call close_file(file2)

    return

  end subroutine qm_generate_input_qchem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_gaussian
  !> @brief        generate Gaussian input file
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of whole molecule
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_gaussian(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    integer                  :: file1, file2
    integer                  :: i, cnt
    integer                  :: address, atomic_no, charge, spin
    character(256)           :: qminp
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    character(1024)          :: string1, string2
    character(1)             :: conv


    ! open files
    !
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    qminp=trim(qmmm%workdir)//'/'//trim(qmmm%qminp)
    call open_file(file2, qminp, IOFileOutputReplace)

    do while (.true.)

      read(file1, '(a)', end=100) string1
      string2 = string1
      call tolower(string1)
      if (index(string1,'#coord') > 0) then

        ! QM atom coordinates
        !
        cnt = 0
        do i = 1, qmmm%qm_natoms
          atomic_no = qmmm%qm_atomic_no(i)
          address = qmmm%qmatom_id(i)
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                        qm_atomic_symbol(atomic_no), coord(1:3,address)
        end do

        ! link atom coordinates
        !
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                       "H  ", qmmm%linkatom_coord(1:3,i)
        end do
        
        !
        ! Sort link atoms
        call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

        !
        ! Write QM coordinates
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          write(file2, *) trim(coord_txt(i))
        end do

        ! end of QM coordinates
        !
        write(file2, *) ''

      else if (index(string1,'#charg') > 0) then

        ! MM point charges
        !
        if (qmmm%mm_natoms > 0) then

          do i = 1, qmmm%mm_natoms
            address = qmmm%mmatom_id(i)
            write(file2, '(3(f18.10,x),f18.10)') &
              coord(1:3,address), molecule%charge(address)
          end do

          write(file2, *) ''

        end if

      else if (index(string1,'#elec') > 0) then

        ! Electric fields on MM atoms
        !
        if (qmmm%mm_natoms > 0) then

          do i = 1, qmmm%mm_natoms
            address = qmmm%mmatom_id(i)
            write(file2, '(3(f18.10,x))') &
              coord(1:3,address)
          end do

          write(file2, *) ''

        end if

      else if (index(string1,"conver=") > 0) then

        write(conv,'(i1)') 8-qmmm%qmtrial*2
        i=index(string1,"conver=")
        string2=string1(1:i+6)//conv//string1(i+8:)
        write(file2,'(a)') trim(string2)

      else

        write(file2, '(a)') trim(string2)

      end if
    end do

100 continue

    call close_file(file1)
    call close_file(file2)

    return

  end subroutine qm_generate_input_gaussian

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_terachem
  !> @brief        generate TeraChem input file
  !! @authors      KY
  !! @param[in]    molecule  : information of whole molecule
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_terachem(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    real(wp)                 :: conv
    integer                  :: file1, file2, file3, file4
    integer                  :: i, cnt
    integer                  :: address, atomic_no, charge
    character(256)           :: qmxyz
    character(256)           :: pcxyz
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    character(1024)          :: string1, string2
    character(80)            :: key
    character(10)            :: linkHname

    character(10), allocatable :: atomname(:)


    ! initialize
    linkHname = ''

    ! open files
    !
    qmxyz=trim(qmmm%qmbasename)//'_qm.xyz'
    pcxyz=trim(qmmm%qmbasename)//'_pc.xyz'
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    call open_file(file2, trim(qmmm%workdir)//'/'//trim(qmmm%qminp), IOFileOutputReplace)

    write(file2, '(a)') 'jobname      '//trim(qmmm%qmbasename)
    qmmm%tcscr = trim(qmmm%tcscr0)//'_'//trim(qmmm%qmbasename)
    write(file2, '(a)') 'scrdir       '//trim(qmmm%tcscr)

    write(file2, '(a)') 'coordinates  '//trim(qmxyz)
    if (qmmm%mm_natoms > 0) write(file2, '(a)') 'pointcharges '//trim(pcxyz)
    write(file2, '(a)') 'amber        yes'

    do while (.true.)

      read(file1, '(a)', end=100) string1

      i = index(trim(adjustl(string1)),'#')
      if (i > 0) then
        if (i == 1) then
          cycle
        else
          string1 = string1(1:i-1)
        end if
      end if

      string2 = string1
      call tolower(string1)

      if (index(string1,'coordinates' ) > 0 .or. &
          index(string1,'pointcharges') > 0 .or. &
          index(string1,'scrdir'      ) > 0 .or. &
          index(string1,'amber'       ) > 0) then
        cycle

      else if (index(string1,"convthre") > 0) then
        read(string2,*) key, conv
        !dbg write(MsgOut,'("debug: key=",a)') trim(key)
        do i = 1, qmmm%qmtrial
           conv = conv * 1.D+01
        end do
        !dbg write(MsgOut,'("debug: conv=",e12.4)') conv
        write(file2,'(a,3x,e12.4)') trim(key), conv

      else if (index(string1,"linkhname") > 0) then
        read(string2,*) key, linkHname

      else if (index(string1,"atomname") > 0) then
        allocate(atomname(qmmm%qm_natoms + qmmm%num_qmmmbonds))
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, '(a)') atomname(i)
        end do

      else

        write(file2, '(a)') trim(string2)

      end if
    end do

100 continue

    call close_file(file1)
    call close_file(file2)

    ! QM atom coordinates
    !
    call open_file(file3, trim(qmmm%workdir)//'/'//trim(qmxyz), IOFileOutputReplace)

    write(file3, '(i0)') qmmm%qm_natoms + qmmm%num_qmmmbonds
    write(file3, *) ''

    if (.not. allocated(atomname)) then
      cnt = 0
      do i = 1, qmmm%qm_natoms
        atomic_no = qmmm%qm_atomic_no(i)
        address = qmmm%qmatom_id(i)
        cnt = cnt + 1
        write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                      qm_atomic_symbol(atomic_no), coord(1:3,address)
      end do

      ! link atom coordinates
      !
      if (len(trim(linkHname)) /= 0) then
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a,3(x,f18.10))') &
                       trim(linkHname), qmmm%linkatom_coord(1:3,i)
        end do

      else
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                       "H  ", qmmm%linkatom_coord(1:3,i)
        end do
      end if

    else
      cnt = 0
      do i = 1, qmmm%qm_natoms
        address = qmmm%qmatom_id(i)
        cnt = cnt + 1
        write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                      atomname(i), coord(1:3,address)
      end do

      ! link atom coordinates
      !
      do i = 1, qmmm%num_qmmmbonds
        cnt = cnt + 1
        write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
              atomname(i + qmmm%qm_natoms), qmmm%linkatom_coord(1:3,i)
      end do

    end if

    !
    ! Sort link atoms
    call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
      qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

    !
    ! Write QM coordinates
    do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
      write(file3, *) trim(coord_txt(i))
    end do


    write(file3, *) ''
    call close_file(file3)

    ! end of QM coordinates
    !


    ! MM point charges
    !
    if (qmmm%mm_natoms > 0) then
      call open_file(file4, trim(qmmm%workdir)//'/'//trim(pcxyz), IOFileOutputReplace)

      write(file4, '(i0)') qmmm%mm_natoms
      write(file4, *) ''

      do i = 1, qmmm%mm_natoms
        address = qmmm%mmatom_id(i)
        write(file4, '(4(f21.16))') &
          molecule%charge(address), coord(1:3,address)
      end do

      write(file4, *) ''

      call close_file(file4)
    end if

    return

  end subroutine qm_generate_input_terachem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_dftbplus
  !> @brief        generate DFTB+ input file
  !! @authors      KY
  !! @param[in]    molecule : information of whole molecule
  !! @param[in]    coord    : atomic coordinates
  !! @param[in]    qmmm     : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_dftbplus(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    integer                  :: file1, file2, mm_charges
    integer                  :: i, j, cnt
    integer                  :: ntyp, ityp, address, atomic_no
    character(1024)          :: string1, string2
    character(256)           :: fname_mm
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    logical                  :: found

    integer, allocatable     :: qm_atomtypes(:)


    ! open files
    !
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    call open_file(file2, trim(qmmm%workdir)//'/'//trim(qmmm%qminp), IOFileOutputReplace)

    ! file name for MM charges
    fname_mm='mm_charges.dat'

    ! QM atom coordinates
    !
    allocate(qm_atomtypes(qmmm%qm_natoms + 1))

    if (qmmm%num_qmmmbonds > 0) then
      ntyp = 1
      qm_atomtypes(ntyp) = 1
    else
      ntyp = 0
    end if

    do i = 1, qmmm%qm_natoms
      atomic_no = qmmm%qm_atomic_no(i)
      found = .false.
      do j = 1, ntyp
        if (atomic_no == qm_atomtypes(j)) then
          found = .true. 
          exit
        end if
      end do
      if (.not. found) then
        ntyp = ntyp + 1
        qm_atomtypes(ntyp) = atomic_no
      end if
    end do

    write(file2, '("Geometry = GenFormat {")')
    write(file2, '(i0," C")') qmmm%qm_natoms + qmmm%num_qmmmbonds
    do i = 1, ntyp
      write(file2, '(a3,$)') qm_atomic_symbol(qm_atomtypes(i))
    end do
    write(file2, *) ''
    write(file2, *) ''

    cnt = 0
    do i = 1, qmmm%qm_natoms
      atomic_no = qmmm%qm_atomic_no(i)
      address = qmmm%qmatom_id(i)
      do j = 1, ntyp
        if (atomic_no == qm_atomtypes(j)) then
          ityp = j
          exit
        end if
      end do
      cnt = cnt + 1
!      write(coord_txt(cnt), '(x,i6,i4,3(x,f18.10))') &
!                    i, ityp, coord(1:3,address)
      write(coord_txt(cnt), '(i4,3(x,f18.10))') &
                    ityp, coord(1:3,address)
    end do

    ! link atom coordinates
    !
    do i = 1, qmmm%num_qmmmbonds
      ityp = 1
      cnt = cnt + 1
!      write(coord_txt(cnt), '(x,i6,i4,3(x,f18.10))') &
!                   i + qmmm%qm_natoms, ityp, qmmm%linkatom_coord(1:3,i)
      write(coord_txt(cnt), '(i4,3(x,f18.10))') &
                   ityp, qmmm%linkatom_coord(1:3,i)
    end do

    !
    ! Sort link atoms
    call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
      qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

    !
    ! Write QM coordinates
    do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
!      write(file2, *) trim(coord_txt(i))
      write(file2, '(x,i6,a)') i, trim(coord_txt(i))
    end do

    write(file2, '("}")')

    ! end of QM coordinates
    !

    do while (.true.)

      read(file1, '(a)', end=100) string1

      i = index(trim(adjustl(string1)),'#')
      if (i > 0) then
        if (i == 1) then
          cycle
        else
          string1 = string1(1:i-1)
        end if
      end if

      string2 = string1
      call tolower(string1)

      if (index(string1,'hamiltonian' ) > 0) then
        write(file2, '(a)') trim(string2)

        ! MM point charges
        !
        if (qmmm%mm_natoms > 0) then
          write(file2,'(2x,"ElectricField = {")')
          write(file2,'(4x,"PointCharges = {")')
          write(file2,'(6x,"CoordsAndCharges [Angstrom] = DirectRead {")')
          write(file2,'(8x,"Records = ",i0)') qmmm%mm_natoms
          write(file2,'(8x,"File    = ",a)') trim(fname_mm)
          call open_file(mm_charges, trim(qmmm%workdir)//'/'//trim(fname_mm), IOFileOutputReplace)
          do i = 1, qmmm%mm_natoms
            address = qmmm%mmatom_id(i)
            write(mm_charges, '(f21.16, 1x, f21.16, 1x, f21.16, 1x, f21.16)') &
              coord(1:3,address), molecule%charge(address)
          end do
          call close_file(mm_charges)


          write(file2, '(6x,"}")')
          write(file2, '(4x,"}")')
          write(file2, '(2x,"}")')
        end if

      else

        write(file2, '(a)') trim(string2)

      end if

    end do

100 continue

    call close_file(file1)
    call close_file(file2)

    return

  end subroutine qm_generate_input_dftbplus

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_gaussian_frwf
  !> @brief        generate Gaussian input file
  !! @authors      KYMD, KY
  !! @param[in]    molecule  : information of whole molecule
  !! @param[in]    coord     : atomic coordinates
  !! @param[in]    qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_gaussian_frwf(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    integer                  :: file1, file2
    integer                  :: i,cnt ! , ii
    integer                  :: address, atomic_no, charge, spin
    integer                  :: ierror
    character(256)           :: qminp
    character(256)           :: string, string_org
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    logical                  :: defined_gen


    defined_gen = qmmm%defined_gen

    ! open files
    !
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    rewind(file1)

    qminp = trim(qmmm%workdir)//'/'//trim(qmmm%qminp)
    call open_file(file2, qminp, IOFileOutputReplace)

    do 
      read(file1, '(a)', iostat = ierror) string
      if (ierror /= 0) exit
      write(string_org, '(a)') string
      call tolower(string)
      if (index(string,'#coord') > 0) then

        ! QM atom coordinates
        !
        cnt = 0
        do i = 1, qmmm%qm_natoms
          atomic_no = qmmm%qm_atomic_no(i)
          address   = qmmm%qmatom_id(i)
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') qm_atomic_symbol(atomic_no), coord(1:3,address)
        end do

        ! link atom coordinates
        !
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') "H  ", qmmm%linkatom_coord(1:3,i)
        end do

        !
        ! Sort link atoms
        call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

        !
        ! Write QM coordinates
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          write(file2, *) trim(coord_txt(i))
        end do

        ! end of QM coordinates
        !
        write(file2, *) ''

      else if (index(string,'#charg') > 0) then

        ! MM point charges
        !
        do i = 1, qmmm%mm_natoms
          address = qmmm%mmatom_id(i)
          write(file2, '(3(f18.10,x),f18.10)') coord(1:3,address), molecule%charge(address)

        end do
        write(file2, *) ''

      else if (index(string,'#elec') > 0) then

        ! Electric fields on MM atoms
        !
        do i = 1, qmmm%mm_natoms
          address = qmmm%mmatom_id(i)
          write(file2, '(3(f18.10,x))') coord(1:3,address)

        end do
        write(file2, *) ''

      else

        write(file2, '(a)') trim(string_org)

      end if
    end do

    call close_file(file1)
    call close_file(file2)

    return

  end subroutine qm_generate_input_gaussian_frwf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_molpro
  !> @brief        generate Molpro input file
  !! @authors      YA
  !! @param[in]    molecule : information of whole molecule
  !! @param[in]    coord    : atomic coordinates
  !! @param[in]    qmmm     : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_molpro(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(in)    :: qmmm

    ! local variables
    integer                  :: file1, file2, file3
    integer                  :: i, cnt
    integer                  :: address, atomic_no, charge, spin
    character(256)           :: qminp, qminppc
    character(256)           :: string, string2
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)


    ! open files
    !
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    qminp = trim(qmmm%workdir)//'/'//trim(qmmm%qminp)
    call open_file(file2, qminp, IOFileOutputReplace)
    qminppc = trim(qmmm%qminp)//'.pc'
    call open_file(file3, trim(qmmm%workdir)//'/'//qminppc, IOFileOutputReplace)

    do while (.true.)

      read(file1, '(a)', end=100) string
      call tolower(string)
      if (string(1:12) .eq. '#coordinate#') then

        ! QM atom coordinates
        !
        cnt = 0
        do i = 1, qmmm%qm_natoms
          atomic_no = qmmm%qm_atomic_no(i)
          address = qmmm%qmatom_id(i)
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                        qm_atomic_symbol(atomic_no), coord(1:3,address)
        end do

        ! link atom coordinates
        !
        do i = 1, qmmm%num_qmmmbonds
          cnt = cnt + 1
          write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
                       "H  ", qmmm%linkatom_coord(1:3,i)
        end do

        !
        ! Sort link atoms
        call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

        !
        ! Write QM coordinates
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          write(file2, *) trim(coord_txt(i))
        end do

      else if (index(string,"#pcfile#") > 0) then

        i = index(string,"#pcfile#")
        string2 = string(1:i-1)//trim(qminppc)//string(i+8:)
        write(file2,'(a)') trim(string2)

      else

        write(file2, '(a)') trim(string)

      end if
    end do

100 continue

    call close_file(file1)
    call close_file(file2)

    ! MM point charges
    !

    if (qmmm%mm_natoms > 0) then

      write(file3, '(a)') "Generated by GENESIS ATDYN"
      write(file3, '(i0)') qmmm%mm_natoms
      do i = 1, qmmm%mm_natoms
        address = qmmm%mmatom_id(i)
        write(file3, '(3(f18.10,x),f18.10," 1")') &
          coord(1:3,address), molecule%charge(address)
      end do
      write(file3, *) ''

    end if

    call close_file(file3)

    return

  end subroutine qm_generate_input_molpro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_generate_input_orca
  !> @brief        generate ORCA input file
  !! @authors      YA, KY
  !! @param[in]    molecule : information of whole molecule
  !! @param[in]    coord    : atomic coordinates
  !! @param[in]    qmmm     : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_generate_input_orca(molecule, coord, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    real(wp)                 :: conv
    integer                  :: file1, file2, file3, file4
    integer                  :: i, i1, cnt
    integer                  :: address, atomic_no, charge
    character(256)           :: qmxyz
    character(256)           :: pcxyz
    character(1024)          :: string1, string2
    character(128)           :: coord_txt(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    character(32)            :: coordfile_kwd = "xyzfile"
    character(10)            :: linkHname

    character(10), allocatable :: atomname(:)


    ! open files
    !
    qmxyz = trim(qmmm%qmbasename)//'_qm.xyz'
    pcxyz = trim(qmmm%qmbasename)//'_pc.xyz'
    call open_file(file1, trim(qmmm%qmcnt), IOFileInput)
    call open_file(file2, trim(qmmm%workdir)//'/'//trim(qmmm%qminp), IOFileOutputReplace)

    write(file2, '(a)') '# Input for ORCA'
    write(file2, '(a)') '# generated by GENESIS'
    write(file2, '(a)') '#'

    do while (.true.)

      read(file1, '(a)', end=100) string1

      string2 = string1
      call tolower(string2)
      i1 = index(string2, trim(coordfile_kwd))
      if (i1 > 0) then
        ! Coordinates of QM atoms
        write(file2, '(a)') trim(string1)//' '//trim(qmxyz)
      else
        write(file2, '(a)') trim(string1)
      end if

    end do

100 continue

    ! Coordinates and charges of MM atoms
    write(file2, '(a)') '%pointcharges "'//trim(pcxyz)//'"'


    call close_file(file1)
    call close_file(file2)

    ! QM atom coordinates
    !
    call open_file(file3, trim(qmmm%workdir)//'/'//trim(qmxyz), IOFileOutputReplace)

    write(file3, '(i0)') qmmm%qm_natoms + qmmm%num_qmmmbonds
    write(file3, *) ''

    cnt = 0
    do i = 1, qmmm%qm_natoms
      atomic_no = qmmm%qm_atomic_no(i)
      address = qmmm%qmatom_id(i)
      cnt = cnt + 1
      write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
        qm_atomic_symbol(atomic_no), coord(1:3,address)
    end do

    ! link atom coordinates
    !
    do i = 1, qmmm%num_qmmmbonds
      cnt = cnt + 1
      write(coord_txt(cnt), '(x,a3,3(x,f18.10))') &
        "H  ", qmmm%linkatom_coord(1:3,i)
    end do

    !
    ! Sort link atoms
    call merge_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
      qmmm%qmatom_id, qmmm%linkatom_global_address, coord_txt)

    !
    ! Write QM coordinates
    do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
      write(file3, *) trim(coord_txt(i))
    end do

    write(file3, *) ''
    call close_file(file3)

    ! MM point charges
    !
    if (qmmm%mm_natoms > 0) then
      call open_file(file4, trim(qmmm%workdir)//'/'//trim(pcxyz), IOFileOutputReplace)

      write(file4, '(i0)') qmmm%mm_natoms
      write(file4, *) ''

      do i = 1, qmmm%mm_natoms
        address = qmmm%mmatom_id(i)
        write(file4, '(4(f21.16))') &
          molecule%charge(address), coord(1:3,address)
      end do

      write(file4, *) ''

      call close_file(file4)
    end if

    return

  end subroutine qm_generate_input_orca

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_qchem
  !> @brief        retrieve energy and gradients from QChem output files
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_qchem(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)                 :: total_energy, mm_energy, force_tmp(3)
    real(wp)                 :: force_tmp_all(3,qmmm%qm_natoms+qmmm%num_qmmmbonds)
    integer                  :: file1, fchk
    integer                  :: i
    integer                  :: nsta, nend
    integer                  :: address, atomic_no, charge, spin
    integer(4)               :: istat, unlink
    character(len_exec)      :: exec
    character(MaxLine)       :: line
    character(MaxFilename)   :: efield, fname
    logical                  :: exist_energy, exist_fchk, exist_force, ex
    logical                  :: exist_fchk42
    logical                  :: found_energy, found_dipole, found_force

    real(wp), pointer :: qm_force(:,:), mm_force(:,:), la_force(:,:)


    ! use pointers
    qm_force => qmmm%qm_force
    mm_force => qmmm%mm_force
    la_force => qmmm%linkatom_force

    qm_error = .true.

    ! inquire files
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_energy)
    if (.not. exist_energy) then
      write(MsgOut,'(a,i0)') &
       'Fatal error in Q-Chem calculation. '//trim(qmmm%qmout)// &
       ' not found! Rank = ', my_country_no
      qm_error = .true.
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check
    !
    qm_error = .true.
    do while (.true.)
      read(file1, '(a)', end=100) line
      i = index(line, 'Have a nice day')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .false.
        exit
      end if

      i = index(line, 'Q-Chem fatal error occurred')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .true.
        exit
      end if
    end do

100 continue

    if (qm_error) then
      write(MsgOut,'(a,i0)') &
       'Fatal error in Q-Chem calculation while reading '//trim(qmmm%qmout)// &
       '! Rank = ', my_country_no
      return
    end if

    ! retrieve QM energy
    !
    total_energy = 0.0_wp
    mm_energy    = 0.0_wp
    found_energy = .false.
    rewind(file1)
    do while (.true.)
      read(file1, '(a)', end=200) line

      ! QMMM energy
      if (line(1:29) .eq. ' The QM part of the Energy is') then
        nsta = 1
        nend = len(line) - 30
        call read_real(line(30:256), nsta, nend, total_energy)
        found_energy = .true.
 
      ! SCF energy
      !else if (line(1:38) .eq. ' Total energy in the final basis set =') then
      !  nsta = 1
      !  nend = len(line) - 39
      !  call read_real(line(39:256), nsta, nend, total_energy)
      !  found_energy = .true.

      ! MP2 energy
      !else if (line(1:34) .eq. '        MP2         total energy =') then
      !  nsta = 1
      !  nend = len(line) - 35
      !  call read_real(line(35:256), nsta, nend, total_energy)
      !  found_energy = .true.
 
      ! RI-MP2 energy
      !else if (line(1:33) .eq. 'RI-MP2 TOTAL ENERGY             =') then
      !  nsta = 1
      !  nend = len(line) - 34
      !  call read_real(line(34:256), nsta, nend, total_energy)
      !  found_energy = .true.
 
      ! Charge-Charge energy
      else if (line(1:27) .eq. ' Charge-charge energy     =') then
        nsta = 1
        nend = len(line) - 28
        call read_real(line(28:256), nsta, nend, mm_energy)
      end if

    end do

200 continue

    qm_error = (.not. found_energy)
    if (qm_error) then
      write(MsgOut,'(a,i0)') &
       'Fatal error in Q-Chem calculation while reading '//trim(qmmm%qmout)// &
       '. Energy is not found! Rank = ', my_country_no
      return
    end if

    if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                 & ", MM electrostatic energy = ", &
                 & f18.10)') total_energy, mm_energy

    qm_energy = total_energy - mm_energy

    ! retrieve dipole moment
    !
    qmmm%qm_dipole = 0.0_wp
    found_dipole   = .false.

    ! Q-Chem 4.3
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.FChk', &
            exist=exist_fchk)
    ! Q-Chem 4.2
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.inp.0.fchk', &
            exist=exist_fchk42)

    if (exist_fchk) then
      call open_file(fchk, &
           trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.FChk', IOFileInput)

      do while (.true.)
        read(fchk, '(a)', end=210) line
        ! QM dipole
        if (line(1:11) .eq. 'Dipole_Data') then
          read(fchk, '(a)') line
          nsta = 1
          nend = 16
          call read_real(line, nsta, nend, qmmm%qm_dipole(1))
          nsta = 17
          nend = 32
          call read_real(line, nsta, nend, qmmm%qm_dipole(2))
          nsta = 33
          nend = 48
          call read_real(line, nsta, nend, qmmm%qm_dipole(3))
          found_dipole = .true. 
          exit
        end if
      end do
      call close_file(fchk)

      if (.not. qmmm%savefile) then
        fname = trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.FChk'
        if (qmmm%qm_debug) &
          write(MsgOut,'(''QMMM_debug> remove '',a)') trim(fname)
        istat = unlink(trim(fname))
        if (istat /= 0) &
          call error_msg ('qm_read_output> Error while removing '//trim(fname))

        fname = trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.0.FChk'
        inquire(file=trim(fname), exist=ex)
        if (ex) then
          if (qmmm%qm_debug) &
            write(MsgOut,'(''QMMM_debug> remove '',a)') trim(fname)
          istat = unlink(trim(fname))
          if (istat /= 0) &
            call error_msg ('qm_read_output> Error while removing '//trim(fname))
        end if

      else if (qmmm%workdir /= qmmm%savedir) then
        fname = trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)
        inquire(file=trim(fname)//'.0.FChk', exist=ex)
        if (ex) then
          exec='mv '//trim(fname)//'.FChk '//trim(fname)//'.0.FChk '//trim(qmmm%savedir)
        else
          exec='mv '//trim(fname)//'.FChk '//trim(qmmm%savedir)
        end if
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

      end if

    else if (exist_fchk42) then
      call open_file(fchk, &
           trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.inp.0.fchk', IOFileInput)

      do while (.true.)
        read(fchk, '(a)', end=210) line
        ! QM dipole
        if (line(1:11) .eq. 'Dipole_Data') then
          read(fchk, '(a)') line
          nsta = 1
          nend = 16
          call read_real(line, nsta, nend, qmmm%qm_dipole(1))
          nsta = 17
          nend = 32
          call read_real(line, nsta, nend, qmmm%qm_dipole(2))
          nsta = 33
          nend = 48
          call read_real(line, nsta, nend, qmmm%qm_dipole(3))
          found_dipole = .true. 
          exit
        end if
      end do
      call close_file(fchk)

      if (.not. qmmm%savefile) then
        fname = trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.inp.0.fchk'
        if (qmmm%qm_debug) &
          write(MsgOut,'(''QMMM_debug> remove '',a)') trim(fname)
        istat = unlink(trim(fname))
        if (istat /= 0) &
          call error_msg ('qm_read_output> Error while removing '//trim(fname))

      else if (qmmm%workdir /= qmmm%savedir) then
        fname = trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)
        exec='mv '//trim(fname)//'.inp.0.fchk '//trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

      end if

    else
      ! read from log file - the digit is only 4.
      rewind(file1)
      do while (.true.)
        read(file1, '(a)', end=210) line
        if (line(21:47) .eq. 'Cartesian Multipole Moments') then
          read(file1, '(a)') line
          read(file1, '(a)') line
          read(file1, '(a)') line
          read(file1, '(a)') line
          read(file1, '(a)') line
          nsta = 11
          nend = 23
          call read_real(line, nsta, nend, qmmm%qm_dipole(1))
          nsta = 31
          nend = 43
          call read_real(line, nsta, nend, qmmm%qm_dipole(2))
          nsta = 51
          nend = 63
          call read_real(line, nsta, nend, qmmm%qm_dipole(3))
          found_dipole = .true. 
          exit
        end if
      end do

    end if

    qmmm%qm_dipole = qmmm%qm_dipole * conv_debye_au

    if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole

210 continue

    call close_file(file1)

    if (qmmm%ene_only) return

    ! retrieve atomic forces
    !
    efield=trim(qmmm%workdir)//'/efield.dat'
    inquire(file=efield, exist=exist_force)

    if (.not. exist_force) then
      write(MsgOut,'(a,i0)') &
        'Fatal error in Q-Chem. efield.dat for '//trim(qmmm%qmbasename)// &
        ' does not exist! Rank =',my_country_no
      qm_error = .true.
      return
    end if

    call open_file(file1, efield, IOFileInput)

    found_force = .false.
    do i = 1, qmmm%mm_natoms
      read(file1, *, err=300) force_tmp(1:3)
      address = qmmm%mmatom_id(i)
      mm_force(1:3,i) = force_tmp(1:3) * molecule%charge(address)
    end do

    do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
      read(file1, *, err=300) force_tmp_all(1:3,i)
    end do

    !
    ! Separate to QM and link atom groups
    call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
      qmmm%qmatom_id, qmmm%linkatom_global_address, &
      force_tmp_all, 3)

    do i = 1, qmmm%qm_natoms
      qm_force(1:3,i) = - force_tmp_all(1:3,i)
    end do
    do i = 1, qmmm%num_qmmmbonds
      la_force(1:3,i) = - force_tmp_all(1:3,qmmm%qm_natoms+i)
    end do

    found_force = .true.
    if (qmmm%savefile) then
      exec='mv '//trim(efield)//' ' &
                //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'_efield.dat'
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

    else
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') trim(efield)
      istat = unlink(trim(efield))
      if (istat /= 0) &
        call error_msg ('qm_read_output> Error while removing '//trim(efield))
    end if

300 continue
    call close_file(file1)

    if (.not. found_force) then
      write(MsgOut,'(a,i0)') &
        'Unexpected EOF while reading efield.dat for '//trim(qmmm%qmbasename)// & 
        '. Rank =', my_country_no
      qm_error = .true.
    end if

    ! ESP charge is not provided by Q-Chem
    qmmm%is_qm_charge = .false.
    if (qmmm%qm_get_esp) then
      qmmm%qmtrial = qmmm%qmmaxtrial
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error while reading data  for '//trim(qmmm%qmbasename)// &
        '! Rank = ', my_country_no
      write(MsgOut,'(a)') 'ESP charge is not supported by Q-Chem.'
      return
    end if

    return

  end subroutine qm_read_output_qchem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_gaussian
  !> @brief        retrieve energy and gradients from Gaussian output files
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_gaussian(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)                 :: rdum
    real(wp)                 :: total_energy, mm_energy
    real(wp)                 :: force_tmp((qmmm%qm_natoms+qmmm%num_qmmmbonds)*3)
    real(wp)                 :: charge_tmp(qmmm%qm_natoms+qmmm%num_qmmmbonds) 
    integer                  :: cnt
    integer                  :: file1, file2
    integer                  :: i, j
    integer                  :: nsta, nend
    integer                  :: address
    integer                  :: ndata, nlines
    integer                  :: idum
    integer(4)               :: istat, unlink
    character(len_exec)      :: exec
    character(MaxLine)       :: line
    character(MaxFilename)   :: fchk
    logical                  :: exist_log, exist_fchk
    logical                  :: found_energy_total, found_energy_charge
    logical                  :: found_force, found_field, found_charge

    real(wp), pointer        :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer        :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    qm_error            = .true.
    found_energy_total  = .false.
    found_energy_charge = .false.
    found_force         = .false.
    found_charge        = .false.
    found_field         = .false.

    ! inquire log file
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
       'Fatal error in Gaussian. '//trim(qmmm%qmout)// &
       ' not found! Rank = ',my_country_no
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check
    !
    qm_error = .true.
    do while (.true.)
      read(file1, '(a)', end=100) line
      i = index(line, 'Normal termination of Gaussian')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .false.
        exit
      end if
    end do

100 continue

    if (qm_error) then
      write(MsgOut,'(a,i0)') &
        'Fatal error in Gaussian calculation while reading ' &
        //trim(qmmm%qmout)//'. Rank = ', my_country_no
      return
    end if

    ! read Gaussian log file
    !
    mm_energy    = 0.0_wp
    rewind(file1)
    do while (.true.)
      read(file1, '(a)', end=200) line
      ! retrieve charge-charge energy
      !
      if (line(2:28) .eq. 'Self energy of the charges') then
        nsta = 30
        nend = 51
        call read_real(line, nsta, nend, mm_energy)
        found_energy_charge = .true.
      end if

      ! retrieve electric field
      !
      if (line(33:64) .eq. '-------- Electric Field --------') then
        do j = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds + 2
          read(file1, '(a)') line
        end do
        do j = 1, qmmm%mm_natoms
          !read(file1, *, err=200) idum, rdum, force_tmp(1:3)
          read(file1, '(a)') line
          read(line(25:66), *, err=200) force_tmp(1:3)
          address = qmmm%mmatom_id(j)
          mm_force(1:3,j) = molecule%charge(address) * force_tmp(1:3)
        end do
        found_field = .true.
      end if

    end do

200 continue

    call close_file(file1)

    ! read Fchk file
    !
    fchk = trim(qmmm%workdir)//'/gaussian.Fchk'

    inquire(file=trim(fchk), exist=exist_fchk)
    if (.not. exist_fchk) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error. gaussian.Fchk for '//trim(qmmm%qmbasename)// &
        ' is not found! Rank = ', my_country_no
      return
    end if

    call open_file(file2, trim(fchk), IOFileInput)

    total_energy = 0.0_wp
    do while (.true.)
      read(file2, '(a)', end=300) line

      ! QMMM total energy
      if (line(1:12) .eq. 'Total Energy') then
        nsta = 50
        nend = 71
        call read_real(line, nsta, nend, total_energy)
        found_energy_total = .true.
      end if
 
      ! retrieve atomic forces
      !
      if (line(1:18) .eq. 'Cartesian Gradient') then
        nlines = (qmmm%qm_natoms + qmmm%num_qmmmbonds) * 3 / 5
        ndata = 5
        nsta = 0
        do i = 1, nlines
          read(file2, *, end=300) force_tmp(nsta+1:nsta+ndata)
          nsta = nsta + 5
        end do
        ndata = mod((qmmm%qm_natoms+qmmm%num_qmmmbonds) * 3, 5)
        if (ndata > 0) then
          read(file2, *, end=300) force_tmp(nsta+1:nsta+ndata)
        end if
        found_force = .true.
      end if

      ! retrieve dipole moment
      !
      if (line(1:13) .eq. 'Dipole Moment') then
        read(file2, *, end=300) qmmm%qm_dipole(1:3)
      end if

      ! retrieve QM ESP charges
      !
      if (line(1:10) .eq. 'ESP Charge') then
        nlines = (qmmm%qm_natoms + qmmm%num_qmmmbonds) / 5
        ndata = 5
        nsta = 0
        do i = 1, nlines
          read(file2, *, end=300) charge_tmp(nsta+1:nsta+ndata)
          nsta = nsta + 5
        end do
        ndata = mod((qmmm%qm_natoms+qmmm%num_qmmmbonds), 5)
        if (ndata > 0) then
          read(file2, *, end=300) charge_tmp(nsta+1:nsta+ndata)
        end if
        found_charge = .true.
      end if

    end do

300 continue

    call close_file(file2)

    qm_error = .not. (found_energy_total )
    if (qm_error) then
      write(MsgOut,'(a,i0)') &
        'Fatal error while reading gaussian.Fchk for ' &
        //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
      if (.not. found_energy_total)  &
        write(MsgOut,'(a)') 'QM energy is not found!'
      if (.not. found_energy_charge) &
        write(MsgOut,'(a)') 'MM energy is not found!'
    end if
    qm_energy = total_energy - mm_energy

    if (.not. qmmm%ene_only) then
      qm_error = .not. (found_force) 
      if (qm_error) then
        write(MsgOut,'(a,i0)') &
          'Fatal error while reading gaussian.Fchk for ' &
          //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
        if (.not. found_force) write(MsgOut,'(a)') 'Force is not found!'
        if (.not. found_field) write(MsgOut,'(a)') 'Field is not found!'
      end if

      !
      ! Separate to QM and link atom groups
      call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
        qmmm%qmatom_id, qmmm%linkatom_global_address, &
        force_tmp, 3)
      call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
        qmmm%qmatom_id, qmmm%linkatom_global_address, &
        charge_tmp, 1)

      nsta = 0
      do i = 1, qmmm%qm_natoms
        qm_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do
      do i = 1, qmmm%num_qmmmbonds
        la_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do
    end if

    if (found_charge) then
      qmmm%is_qm_charge = .true.
      nsta = 1
      do i = 1, qmmm%qm_natoms
        qm_charge(i) = charge_tmp(nsta)
        nsta = nsta + 1
      end do
      do i = 1, qmmm%num_qmmmbonds
        la_charge(i) = charge_tmp(nsta)
        nsta = nsta + 1
      end do

    else
      qmmm%is_qm_charge = .false.
      if (qmmm%qm_get_esp) then
        qm_error = .true.
        write(MsgOut,'(a,i0)') &
          'Fatal error while reading gaussian.Fchk for '//trim(qmmm%qmbasename)// &
          '! Rank = ', my_country_no
        write(MsgOut,'(a)') 'ESP charge is not found!'
      end if

    end if

    if (qm_error) return

    if (qmmm%qm_debug) then
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                     ", MM electrostatic energy = ", f18.10)') &
                      total_energy, mm_energy
      write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole

    end if

    if (qmmm%savefile) then
      exec='mv '//trim(fchk)//' ' &
                //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'.Fchk'
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

    else
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') trim(fchk)
      istat = unlink(trim(fchk))
      if (istat /= 0) &
        call error_msg ('qm_read_output> Error while removing '//trim(fchk))

    end if

    return

  end subroutine qm_read_output_gaussian

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_terachem
  !> @brief        retrieve energy and gradients from TeraChem output files
  !! @authors      KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_terachem(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)               :: force_tmp(3)
    real(wp)               :: force_tmp_all(3,qmmm%qm_natoms+qmmm%num_qmmmbonds)
    real(wp)               :: charge_tmp(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    integer                :: file1, file2
    integer                :: i, j, idx1, idx2
    character(len_exec)    :: exec
    character(MaxLine)     :: line
    character(MaxFilename) :: dat, grad
    logical                :: exist_log, exist_dat, exist_grad, exist_wf
    logical                :: found_energy_total
    logical                :: found_force, found_mm_force, found_charge

    real(wp), pointer      :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer      :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    qm_error            = .true.
    found_energy_total  = .false.
    found_force         = .false.
    found_charge        = .false.
    found_mm_force      = .false.

    ! inquire log file
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
       'Fatal error in TeraChem '//trim(qmmm%qmout)// &
       ' not found! Rank = ',my_country_no
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check
    !
    qm_error = .true.
    do while (.true.)
      read(file1, '(a)', end=100) line
      i = index(line, 'Job finished')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .false.
        exit
      end if
    end do

100 continue

    if (qm_error) then
      write(MsgOut,'(a,i0)') &
        'Fatal error in TeraChem calculation while reading ' &
        //trim(qmmm%qmout)//'. Rank = ', my_country_no
      return
    end if

    ! read TeraChem log file
    !
    rewind(file1)
    do while (.true.)
      read(file1, '(a)', end=200) line

      ! retrieve dipole moment
      !
      if (line(1:17) .eq. 'QM  DIPOLE MOMENT') then
        idx1 = index(line,'{')+1
        idx2 = index(line,'}')-1
        read(line(idx1:idx2), *) qmmm%qm_dipole(1:3)
        qmmm%qm_dipole(1:3) = qmmm%qm_dipole(1:3) * conv_debye_au
      end if

      ! retrieve MM gradient
      !
      if (index(line, 'Point charge part') > 0) then
        do j = 1, qmmm%mm_natoms
          read(file1, *, err=200) force_tmp(1:3)
          mm_force(1:3,j) = -force_tmp(1:3)
        end do
        found_mm_force = .true.
      end if

      ! retrieve QM ESP charges
      !
      if (line(1:23) .eq. 'ESP unrestraint charges' .or. &
          line(1:21) .eq. 'ESP restraint charges') then
        read(file1, '(a)') line
        read(file1, '(a)') line
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, '(a)') line
          read(line(38:48), *) charge_tmp(i)
        end do
        !
        ! Separate to QM and link atom groups
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          charge_tmp, 1)

        qm_charge(:) = charge_tmp(1:qmmm%qm_natoms)
        la_charge(:) = charge_tmp(qmmm%qm_natoms+1:)
        
        found_charge = .true.
      end if

    end do

200 continue

    call close_file(file1)

    ! read scr/results.dat
    !
    dat = trim(qmmm%workdir)//'/'//trim(qmmm%tcscr)//'/results.dat'

    inquire(file=trim(dat), exist=exist_dat)
    if (.not. exist_dat) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error. '//trim(qmmm%tcscr)//'/results.dat for ' &
         //trim(qmmm%qmbasename)//' is not found! Rank = ', my_country_no
      return
    end if

    call open_file(file2, trim(dat), IOFileInput)

    qm_energy = 0.0_wp
    do while (.true.)
      read(file2, '(a)', end=300) line

      ! QMMM total energy
      if (line(1:19) .eq. 'Ground state energy') then
        read(file2, *) qm_energy
        found_energy_total = .true.
      end if
    end do

300 continue

    call close_file(file2)

    ! retrieve atomic forces from scr/grad.xyz
    !
    grad = trim(qmmm%workdir)//'/'//trim(qmmm%tcscr)//'/grad.xyz'

    inquire(file=trim(grad), exist=exist_grad)
    if (exist_grad) then
      found_force = .true.

      call open_file(file2, trim(grad), IOFileInput)

      read(file2, *)
      read(file2, *)
      do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
        read(file2, '(a)') line
        read(line(4:), *) force_tmp_all(1:3,i)
      end do
      !
      ! Separate to QM and link atom groups
      call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
        qmmm%qmatom_id, qmmm%linkatom_global_address, &
        force_tmp_all, 3)

      qm_force(:,:) = -force_tmp_all(:,1:qmmm%qm_natoms)
      la_force(:,:) = -force_tmp_all(:,qmmm%qm_natoms+1:)

      call close_file(file2)

    end if

    ! error check, and copy data
    !
    if (.not. found_energy_total) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error. '//trim(qmmm%tcscr)//'/results.dat for ' &
        //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
      write(MsgOut,'(a)') 'QM energy is not found!'
    end if

    if (.not. qmmm%ene_only .and. .not. found_force) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error while reading '//trim(qmmm%tcscr)//'/grad.xyz for ' &
        //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
      if (.not. found_force) write(MsgOut,'(a)') 'Force is not found!'
    end if

    if (.not. qmmm%ene_only .and. .not. found_mm_force) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error while reading '//trim(qmmm%qmout)//'. Rank = ', my_country_no
      write(MsgOut,'(a)') 'MM Force is not found!'
    end if

    if (found_charge) then
      qmmm%is_qm_charge = .true.

    else
      qmmm%is_qm_charge = .false.
      if (qmmm%qm_get_esp) then
        qm_error = .true.
        write(MsgOut,'(a,i0)') &
          'Fatal error in TeraChem calculation while reading ' &
          //trim(qmmm%qmout)//'. Rank = ', my_country_no
        write(MsgOut,'(a)') 'ESP charge is not found!'
      end if

    end if
    if (qm_error) return

    if (qmmm%qm_debug) then
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                     ", MM electrostatic energy = ", f18.10)') &
                      qm_energy, 0.0_wp
      write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole

    end if

    if (qmmm%savefile) then

      if (qmmm%workdir /= qmmm%savedir) then

        ! save scrdir
        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%tcscr)//' '//trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

        ! save xyz files
        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_qm.xyz '& 
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_pc.xyz ' &
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

      end if

    else

      ! remove scrdir (including result.dat and grad.xyz)
      exec='rm -rf '//trim(qmmm%workdir)//'/'//trim(qmmm%tcscr)
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

      ! xyz files
      exec='rm '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_qm.xyz ' & 
                //trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_pc.xyz'
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

    end if

    return

  end subroutine qm_read_output_terachem

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_dftbplus
  !> @brief        retrieve energy and gradients from DFTB+ output files
  !! @authors      KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_dftbplus(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)               :: force_tmp(4)
    real(wp)               :: force_tmp_all(3,qmmm%qm_natoms+qmmm%num_qmmmbonds)
    real(wp)               :: charge_tmp(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    integer                :: file1
    integer                :: i, j
    integer(4)             :: istat, unlink
    character(len_exec)    :: exec
    character(MaxLine)     :: line
    character(MaxFilename) :: fname
    logical                :: exist_log, ex
    logical                :: found_energy_total
    logical                :: found_force, found_mm_force, found_charge
    logical                :: found_version 

    real(wp), pointer      :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer      :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    qm_error            = .false.
    found_energy_total  = .false.
    found_force         = .false.
    found_charge        = .false.
    found_mm_force      = .false.
    found_version       = .false.
    qm_energy = 0.0_wp
    force_tmp = 0.0_wp 

    ! inquire log file
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
       'Fatal error in DFTB+. '//trim(qmmm%qmout)// &
       ' not found! Rank = ',my_country_no
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! convergence check
    !
    do while (.true.)
      read(file1, '(a)', end=100) line
      if (index(line, 'SCC is NOT converged') > 0) then
        qm_error = .true.
        write(MsgOut,'(a,i0)') &
          'Fatal error in DFTB+. SCC is not converged. ' &
          //' Rank = ', my_country_no
        return
      end if
    end do

100 continue

    call close_file(file1)

    ! inquire detailed.out file
    !
    inquire(file=trim(qmmm%workdir)//'/detailed.out', exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
       'Fatal error in DFTB+. detailed.out not found! Rank = ',my_country_no
      return
    end if

    ! open detailed.out file
    !
    call open_file(file1, trim(qmmm%workdir)//'/detailed.out', IOFileInput)

    do while (.true.)
      read(file1, '(a)', end=200) line

      ! retrieve self-consistent charges
      !
      if (index(line,'Atomic gross charges') > 0) then
        read(file1, *)
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, *) j, charge_tmp(i)
        end do
        !
        ! Separate to QM and link atom groups
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          charge_tmp, 1)

        qm_charge(:) = charge_tmp(1:qmmm%qm_natoms)
        la_charge(:) = charge_tmp(qmmm%qm_natoms+1:)

        found_charge = .true.
      end if

      ! QMMM total energy
      !
      if (index(line,'Total Mermin free energy') > 0) then
        read(line(26:49), *) qm_energy
        found_energy_total = .true.
      end if

      ! Check version of DFTB+
      !
      if (index(line,'Atom Sh.   l   m       Population  Label') > 0) then
        ! This program is more than version 19.1
        found_version = .true.
      end if    

      ! retrieve QM gradient
      !
      if (index(line,'Total Forces') > 0) then
        found_force = .true.

        if (found_version) then
          ! For dftb+ ver.19.1 & 20.1 
          do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
            read(file1, *) force_tmp(1:4)
            force_tmp_all(1:3,i) = force_tmp(2:4)
          end do          
        else
          ! For dftb+ ver.18.1 
          do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
            read(file1, *) force_tmp(1:3)
            force_tmp_all(1:3,i) = force_tmp(1:3)
          end do
        end if

        !
        ! Separate to QM and link atom groups
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          force_tmp_all, 3)

        qm_force(:,:) = force_tmp_all(:,1:qmmm%qm_natoms)
        la_force(:,:) = force_tmp_all(:,qmmm%qm_natoms+1:)

      end if

      ! retrieve MM gradient
      !
      if (index(line,'Forces on external charges') > 0) then
        found_mm_force = .true.

        do j = 1, qmmm%mm_natoms
          read(file1, *, err=200) force_tmp(1:3)
          mm_force(1:3,j) = force_tmp(1:3)
        end do

      end if

      ! retrieve dipole moment
      !
      if (index(line,'Dipole moment') > 0 .and. index(line,'au') > 0) then
        read(line(15:56), *) qmmm%qm_dipole(1:3)
      end if

    end do

200 continue

    call close_file(file1)

    ! error check, and copy data
    !
    if (.not. found_energy_total) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error in detailed.out for ' &
        //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
      write(MsgOut,'(a)') 'QM energy is not found!'
    end if

    if (.not. qmmm%ene_only .and. .not. found_force) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
        'Fatal error while reading detailed.out for ' &
        //trim(qmmm%qmbasename)//'! Rank = ', my_country_no
      if (.not. found_force) write(MsgOut,'(a)') 'Force is not found!'
    end if

    if (found_charge) then
      qmmm%is_qm_charge = .true.

    else
      qmmm%is_qm_charge = .false.
      if (qmmm%qm_get_esp) then
        qm_error = .true.
        write(MsgOut,'(a,i0)') &
          'Fatal error while reading detailed.out. ' &
          //' Rank = ', my_country_no
        write(MsgOut,'(a)') 'SCC charge is not found!'
      end if

    end if
    if (qm_error) return

    if (qmmm%qm_debug) then
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                     ", MM electrostatic energy = ", f18.10)') &
                      qm_energy, 0.0_wp
      write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole

    end if

    ! detailed.out
    fname = trim(qmmm%workdir)//'/detailed.out'
    if (qmmm%savefile) then
      exec='mv '//trim(fname)//' ' &
                //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'_detailed.out'
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

    else
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') trim(fname)
      istat = unlink(trim(fname))
      if (istat /= 0) &
        call error_msg ('qm_read_output> Error while removing '//trim(fname))

    end if

    ! charges.bin or charges.dat
    if (qmmm%savefile) then
      if (qmmm%charges_bin) then
        ! save charges.bin
        inquire(file=trim(qmmm%workdir)//'/charges.bin', exist=ex)
        if (ex) then
          exec='cp '//trim(qmmm%workdir)//'/charges.bin ' &
               //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'_charges.bin'
          if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
          call system(exec)
        end if

      else
        ! save charges.dat
        inquire(file=trim(qmmm%workdir)//'/charges.dat', exist=ex)
        if (ex) then
          exec='cp '//trim(qmmm%workdir)//'/charges.dat ' &
               //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'_charges.dat'
          if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
          call system(exec)

        end if

      end if
    end if

    return

  end subroutine qm_read_output_dftbplus

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_gaussian_frwf
  !> @brief        retrieve energy and gradients from Gaussian FORMATTED rwf,
  !                which is generated with system-called formrwf.
  !                Please contact KYMD to have formrwf program. 
  !! @authors      YA, KY, KYMD
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_gaussian_frwf(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),     intent(in)    :: molecule
    type(s_qmmm), target, intent(inout) :: qmmm
    real(wp),             intent(inout) :: qm_energy
    logical,              intent(inout) :: qm_error

    ! local variables
    real(wp)               :: rdum
    real(wp)               :: total_energy, mm_energy
    real(wp)               :: force_tmp((qmmm%qm_natoms+qmmm%num_qmmmbonds)*3)
    real(wp)               :: charge_tmp(qmmm%qm_natoms+qmmm%num_qmmmbonds) 
    real(wp)            :: field_tmp(1: 5)
    integer                :: file1
    integer                :: i, j
    integer                :: nsta, nend
    integer                :: address
    integer                :: ndata, nlines
    integer                :: idum
    integer             :: ierror, natom_qmlink, jatom_mm
    character(len_exec)    :: exec
    character(MaxLine)     :: line
    character(MaxFilename) :: rwf, frwf
    logical                :: exist_log, exist_rwf, exist_frwf
    logical                :: found_energy, found_force, found_charge
    logical                :: found_dipole, found_selfe, found_efield

    real(wp), pointer      :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer      :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force => qmmm%qm_force
    mm_force => qmmm%mm_force
    la_force => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge


    ! inquire log file
    !
    qm_error = .false.
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist = exist_log)
    if (.not.exist_log) then
      if (main_rank .or. replica_main_rank) &
        write(MsgOut,'(a)') 'Fatal error in Gaussian: NO Output file was found!'
      qm_error = .true.
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check & read Self energy with Output file
    !
    qm_error     = .true.
    found_selfe  = .false.
    do 
      read(file1, '(a)', iostat = ierror) line
      if (ierror /= 0) exit
      if (line(2:27) == 'Self energy of the charges') then
        nsta = 30
        nend = 50
        call read_real(line, nsta, nend, mm_energy)
        found_selfe = .true.
        cycle
      end if
      i = index(line, 'Normal termination of Gaussian')
      if (i /= 0) then
        qm_error = .false.
        exit
      end if
    end do

    ! close log file
    !
    call close_file(file1)
    if (qm_error) then 
      if (main_rank .or. replica_main_rank) &
        write(MsgOut,'(a)') 'Fatal error in Gaussian execution!'
      return
    end if

    ! inquire rwf file
    !> When frwf had been generated normally, this rwf was removed.
    !
    rwf  = trim(qmmm%workdir)//'/'//'gaussian.rwf'
    frwf = trim(qmmm%workdir)//'/'//'gaussian.Frwf'
    inquire(file = trim(rwf), exist = exist_rwf)
    if (exist_rwf) then
      if (main_rank .or. replica_main_rank) then
        write(MsgOut,'(a)') "Fatal error in Gaussian: NO Gaussian 'Formatted' RWF was found!"
        write(MsgOut,'(a)') '.....   This routine requires formrwf convert program   .....'
      end if
      qm_error = .true.
      return
    end if

    ! inquire Frwf file
    !> This file was generated in 'rungaussian.sh'
    !
    inquire(file = trim(frwf), exist = exist_frwf)
    if (.not.exist_frwf) then
      if (main_rank .or. replica_main_rank) &
        write(MsgOut,'(a)') 'Fatal error in Gaussian: NO Gaussian Formatted RWF was found!'
      qm_error = .true.
      return
    end if

    ! open frwf with read only
    !
    call open_file(file1, trim(frwf), IOFileInput)

    ! retrieve QM energy
    !
    found_energy = .false.
    found_force  = .false.
    found_dipole = .false.
    !found_selfe  = .false.
    found_efield = .false.
    found_charge = .false.
    do 
      read(file1, '(a)', iostat = ierror) line
      if (ierror /= 0) exit

      ! total energy of system consisting of QM and point charges
      !
      if (line(1:12) == 'Total Energy') then
        nsta = 50
        nend = 71
        call read_real(line, nsta, nend, total_energy)
        found_energy = .true.
      end if
 
      ! retrieve QM gradients
      !
      natom_qmlink = (qmmm%qm_natoms+qmmm%num_qmmmbonds)
      if (line(1:18) == 'Cartesian Gradient') then
        ndata = 5
        nend  = 0
        do 
          nsta = nend + 1
          nend = nend + ndata
          if (nend > natom_qmlink * 3) nend = natom_qmlink * 3
          read(file1, '(5e16.8)') force_tmp(nsta: nend)
          if (nend == natom_qmlink * 3) exit
        end do
        ! Separate to QM and link atom groups
        !
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          force_tmp, 3)

        found_force = .true.
      end if

      ! retrieve dipole moment
      !
      if (line(1:13) == 'Dipole Moment') then
        read(file1, *) qmmm%qm_dipole(1:3)
        found_dipole = .true.
      end if

      ! MM (point charges) self energy of system consisting of QM and point charges
      !
      ! ===== Obtained from Output file
      !if (line(1:14) == 'MM Self Energy') then
      !  nsta = 50
      !  nend = 71
      !  call read_real(line, nsta, nend, mm_energy)
      !  found_selfe = .true.
      !end if

      ! retrieve QM ESP charges
      !
      if (line(1:10) .eq. 'ESP Charge') then
        nlines = (qmmm%qm_natoms + qmmm%num_qmmmbonds) / 5
        ndata = 5
        nsta = 0
        do i = 1, nlines
          read(file1, *) charge_tmp(nsta+1:nsta+ndata)
          nsta = nsta + 5
        end do
        ndata = mod((qmmm%qm_natoms+qmmm%num_qmmmbonds), 5)
        if (ndata > 0) then
          read(file1, *) charge_tmp(nsta+1:nsta+ndata)
        end if
        !
        ! Separate to QM and link atom groups
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          charge_tmp, 1)

        found_charge = .true.
      end if

      ! retrieve MM electric fields
      !
      if (line(1:18) == 'MM Electric Fields') then
        ndata    = 5
        nend     = 0
        jatom_mm = 0
        do 
          nsta = nend + 1
          nend = nend + ndata
          if (nend > qmmm%mm_natoms * 3) nend = qmmm%mm_natoms * 3
          read(file1, '(5e16.8)') field_tmp(1: 5)
          !call pass_value(qmmm%mm_natoms, nsta, nend, jatom_mm, field_tmp, mm_force)
          call pass_value2(qmmm%mm_natoms, nsta, nend, field_tmp, mm_force)
          if (nend == qmmm%mm_natoms * 3) exit
        end do
        found_efield = .true.
      end if
    end do
    call close_file(file1)

    if (.not. qmmm%ene_only) then
      qm_error = .not.(found_energy .and. found_force .and. found_dipole .and. &
                       found_selfe .and. found_efield)
      if (qm_error) then
        if (main_rank .or. replica_main_rank) then
          if (.not.found_energy) write(MsgOut,'(a)') 'Fatal error in Gaussian. Energy not found!'
          if (.not.found_force)  write(MsgOut,'(a)') 'Fatal error in Gaussian. Grad not found!'
          if (.not.found_dipole) write(MsgOut,'(a)') 'Fatal error in Gaussian. Dipole not found!'
          if (.not.found_selfe)  write(MsgOut,'(a)') 'Fatal error in Gaussian. Self Energy not found!'
          if (.not.found_efield) write(MsgOut,'(a)') 'Fatal error in Gaussian. Electric field not found!'
        end if
        return
      end if

      nsta = 0
      do i = 1, qmmm%qm_natoms
        qm_force(1:3, i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do
      do i = 1, qmmm%num_qmmmbonds
        la_force(1:3, i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do

      do j = 1, qmmm%mm_natoms
        address = qmmm%mmatom_id(j)
        mm_force(1: 3, j) = molecule%charge(address) * mm_force(1: 3, j)
      end do

      if (found_charge) then
        qmmm%is_qm_charge = .true.
        nsta = 1
        do i = 1, qmmm%qm_natoms
          qm_charge(i) = charge_tmp(nsta)
          nsta = nsta + 1
        end do
        do i = 1, qmmm%num_qmmmbonds
          la_charge(i) = charge_tmp(nsta)
          nsta = nsta + 1
        end do

      else
        qmmm%is_qm_charge = .false.
        if (qmmm%qm_get_esp) then
          qm_error = .true.
          write(MsgOut,'(a,i0)') 'Fatal error in Gaussian. Rank = ', my_country_no
          write(MsgOut,'(a)') 'ESP charge is not found!'
          return
        end if

      end if

    else
      qm_error = .not.(found_energy .and. found_dipole .and. found_selfe )
      if (qm_error) then
        if (main_rank .or. replica_main_rank) then
          if (.not.found_energy) write(MsgOut,'(a)') 'Fatal error in Gaussian. Energy not found!'
          if (.not.found_dipole) write(MsgOut,'(a)') 'Fatal error in Gaussian. Dipole not found!'
          if (.not.found_selfe) write(MsgOut,'(a)') 'Fatal error in Gaussian. Self Energy not found!'
        end if
        return
      end if

    end if

    qm_energy      = total_energy - mm_energy
    !qmmm%qm_dipole = qmmm%qm_dipole * conv_debye_au


    if (qmmm%savefile) then
      exec='mv '//trim(frwf)//' ' &
                //trim(qmmm%savedir)//'/'//trim(qmmm%qmbasename)//'.Frwf'
    else
      exec='rm -f '//trim(frwf)
    end if
    if (replica_main_rank) call system(exec)

    if (qmmm%qm_debug) then
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, ", MM electrostatic energy = ", f18.10)') &
       total_energy, mm_energy
      write(MsgOut,'("QMMM_debug> QM dipole in a.u. = ", 3f12.4)') qmmm%qm_dipole
      write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
    end if

    return

  contains

    !====1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    pass_value
    !> @brief        pass value in 1D vector to 2D matrix
    !! @authors      KYMD
    !! @param[in]    natom_mm   : # of MM atoms
    !! @param[in]    istart     : # of first item in read line
    !! @param[in]    iend       : # of last item in read line
    !! @param[in]    jpos_save  : Position of column in 2D matrix
    !! @param[in]    field_tmp  : 1D vector
    !! @param[out]   mm_force   : 2D matrix
    !
    !====1=========2=========3=========4=========5=========6=========7=========8

    subroutine pass_value(natom_mm, istart, iend, jpos_save, field_tmp, &
                          mm_force)

      ! formal arguments
      integer,               intent(in)    :: natom_mm
      integer,               intent(in)    :: istart
      integer,               intent(in)    :: iend
      integer,               intent(inout) :: jpos_save
      real(wp),              intent(in)    :: field_tmp(5)
      real(wp),              intent(out)   :: mm_force(3, natom_mm)

      ! local variables
      integer :: ipass_pattern, iend_pattern
      integer :: i, j, jend


      ipass_pattern = mod(istart, 3)
      iend_pattern  = mod(iend, 5)
      if (ipass_pattern == 1) then
        do i = 1, 3
            mm_force(i, jpos_save + 1) = field_tmp(i)
        end do
        if (iend_pattern == 3) return
        do i = 1, 2
            mm_force(i, jpos_save + 2) = field_tmp(i + 3)
        end do
        jpos_save = jpos_save + 1
      else if (ipass_pattern == 0) then
        do i = 3, 3
            mm_force(i, jpos_save + 1) = field_tmp(i - 2)
        end do
        if (iend_pattern == 1) return
        do i = 1, 3
            mm_force(i, jpos_save + 2) = field_tmp(i + 1)
        end do
        if (iend_pattern == 4) return
        do i = 1, 1
            mm_force(i, jpos_save + 3) = field_tmp(i + 4)
        end do
        jpos_save = jpos_save + 2
      else if (ipass_pattern == 2) then
        do i = 2, 3
            mm_force(i, jpos_save + 1) = field_tmp(i - 1)
        end do
        if (iend_pattern == 2) return
        do i = 1, 3
            mm_force(i, jpos_save + 2) = field_tmp(i + 2)
        end do
        jpos_save = jpos_save + 1
      end if

      return

    end subroutine pass_value

    !====1=========2=========3=========4=========5=========6=========7=========8
    !
    !  Subroutine    pass_value2
    !> @brief        pass value in 1D vector to 2D matrix
    !! @authors      KYMD
    !! @param[in]    natom_mm   : # of MM atoms
    !! @param[in]    istart     : # of first item in read line
    !! @param[in]    iend       : # of last item in read line
    !! @param[in]    field_tmp  : 1D vector
    !! @param[out]   mm_force   : 2D matrix
    !
    !====1=========2=========3=========4=========5=========6=========7=========8

    subroutine pass_value2(natom_mm, istart, iend, field_tmp, mm_force)

      ! formal arguments
      integer,               intent(in)    :: natom_mm
      integer,               intent(in)    :: istart
      integer,               intent(in)    :: iend
      real(wp),              intent(in)    :: field_tmp(5)
      real(wp),              intent(out)   :: mm_force(1: 3*natom_mm, 1: 1)

      ! local variables
      integer :: jend


      jend = max(5 - mod(iend, 5), mod(iend, 5))
      mm_force(istart: iend, 1) = field_tmp(1: jend)

      return

    end subroutine pass_value2

  end subroutine qm_read_output_gaussian_frwf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_molpro
  !> @brief        retrieve energy and gradients from Molpro output files
  !! @authors      YA, KY
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_molpro(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)             :: total_energy, mm_energy
    real(wp)             :: force_tmp((qmmm%qm_natoms+qmmm%num_qmmmbonds)*3)
    integer              :: file1
    integer              :: i, j
    integer              :: nsta, nend
    integer              :: idum
    character(256)       :: line
    logical              :: exist_log
    logical              :: found_energy, found_force, found_mmforce
    logical              :: found_dipole

    real(wp), pointer    :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer    :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    qm_error = .true.

    ! inquire log file
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut, *) trim(qmmm%qmout)//" not found !"
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check
    !
    qm_error = .true.
    do while (.true.)
      read(file1, '(a)', end=100) line
      i = index(line, 'Variable memory released')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .false.
        exit
      end if
    end do

100 continue

    if (qm_error) then
      write(MsgOut,'(a)') 'Fatal error in Molpro calculation for ' &
        //trim(qmmm%qmout)//'!'
      return
    end if

    ! read molpro log file
    !
    found_energy  = .false.
    found_force   = .false.
    found_mmforce = .false.
    found_dipole  = .false.
    rewind(file1)
    do while (.true.)
      read(file1, '(a)', end=200) line

      ! retrieve total energy
      !
      if (line(2:19) .eq. 'ENERGY           =') then
        read(line(20:), *, err=200) total_energy
        found_energy = .true.
      end if

      ! retrieve dipole moment
      !
      if (line(2:29) .eq. '[DMX,DMY,DMZ](1:3)       = [') then
        i = index(line, '] AU')
        read(line(30:i-1), *, err=200) qmmm%qm_dipole(1:3)
        found_dipole = .true.
      end if

      ! retrieve charge-charge energy
      !
      ! total_energy in Molpro doesn't include charge-charge interaction
      mm_energy     = 0.0_wp

      ! retrieve QM forces
      !
      if (line(2:20) .eq. 'Atom          dE/dx') then
        read(file1, '(a)') line
        nsta = 0
        do j = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, *, err=200) idum, force_tmp(nsta+1:nsta+3)
          nsta = nsta + 3
        end do
        ! Separate to QM and link atom groups
        !
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          force_tmp, 3)

        found_force = .true.
      end if

      ! retrieve point charge forces
      !
      if (line(2:26) .eq. 'Point Charge        dE/dx') then
        read(file1, '(a)') line
        do j = 1, qmmm%mm_natoms
           read(file1, *, err=200) idum, mm_force(1:3,j)
           mm_force(1:3,j) = - mm_force(1:3,j)
        end do
        found_mmforce = .true.
      end if

    end do

200 continue

    call close_file(file1)

    qm_error = .not. found_energy
    if (qm_error) then
      write(MsgOut,'(a)') &
        'Fatal error in Molpro while reading' &
        //trim(qmmm%qmout)//'! Rank = ', my_country_no
      write(MsgOut,'(a)') 'QM energy is not found!'
      return
    end if
    qm_energy = total_energy

    if (.not. qmmm%ene_only) then
      qm_error = .not. (found_force .and. found_mmforce)
      if (qm_error) then
        write(MsgOut,'(a)') &
          'Fatal error in Molpro while reading' &
          //trim(qmmm%qmout)//'! Rank = ', my_country_no
        if (.not. found_force) write(MsgOut,'(a)') 'Force  is not found!'
        if (.not. found_mmforce) write(MsgOut,'(a)') 'MM-Force  is not found!'
        return
      end if

      nsta = 0
      do i = 1, qmmm%qm_natoms
        qm_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do
      do i = 1, qmmm%num_qmmmbonds
        la_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        nsta = nsta + 3
      end do

    end if

    if (qmmm%qm_debug) then 
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                 & ", MM electrostatic energy = ", &
                 & f18.10)') total_energy, mm_energy
      write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole 
    end if

    return

  end subroutine qm_read_output_molpro

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    qm_read_output_orca
  !> @brief        retrieve energy and gradients from ORCA output files
  !! @authors      YA
  !! @param[in]    molecule  : information of molecules
  !! @param[in]    qmmm      : QM/MM information
  !! @param[inout] qm_energy : QM energy
  !! @param[inout] qm_error  : error flag
  !! @note  : This routine is called only from replica_main_rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine qm_read_output_orca(molecule, qmmm, qm_energy, qm_error)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm), target,    intent(inout) :: qmmm
    real(wp),                intent(inout) :: qm_energy
    logical,                 intent(inout) :: qm_error

    ! local variables
    real(wp)             :: force_tmp((qmmm%qm_natoms+qmmm%num_qmmmbonds)*3)
    real(wp)             :: charge_tmp(qmmm%qm_natoms+qmmm%num_qmmmbonds)
    integer              :: file1
    integer              :: i, j, nline, nsta, nend
    character(len_exec)  :: exec
    character(MaxLine)   :: line
    logical              :: exist_log
    logical              :: found_energy, found_dipole
    logical              :: found_force, found_mm_force, found_charge

    real(wp), pointer    :: qm_force(:,:), mm_force(:,:), la_force(:,:)
    real(wp), pointer    :: qm_charge(:), la_charge(:)


    ! use pointers
    qm_force  => qmmm%qm_force
    mm_force  => qmmm%mm_force
    la_force  => qmmm%linkatom_force
    qm_charge => qmmm%qm_charge
    la_charge => qmmm%linkatom_charge

    qm_error       = .true.
    found_energy   = .false.
    found_force    = .false.
    found_charge   = .false.
    found_mm_force = .false.

    ! inquire log file
    !
    inquire(file=trim(qmmm%workdir)//'/'//trim(qmmm%qmout), exist=exist_log)
    if (.not. exist_log) then
      qm_error = .true.
      write(MsgOut,'(a,i0)') &
       'Fatal error in ORCA '//trim(qmmm%qmout)// &
       ' not found! Rank = ',my_country_no
      return
    end if

    ! open log file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmout), IOFileInput)

    ! error check
    !
    qm_error = .true.
    do while (.true.)
      read(file1, '(a)', end=100) line
      i = index(line, 'ORCA TERMINATED NORMALLY')
      if (i > 0 .and. i <= len(line)) then
        qm_error = .false.
        exit
      end if
    end do

100 continue

    if (qm_error) then
      write(MsgOut,'(a,i0)') &
        'Fatal error in ORCA calculation while reading ' &
        //trim(qmmm%qmout)//'. Rank = ', my_country_no
      return
    end if

    ! read ORCA log file
    !
    found_energy  = .false.
    found_force   = .false.
    found_mm_force = .false.
    found_dipole  = .false.
    found_charge  = .false.
    rewind(file1)
    do while (.true.)
      read(file1, '(a)', end=200) line

      ! retrieve dipole moment
      !
      if (index(line, "Total Dipole Moment    :") > 0) then
        nsta = 25
        nend = 37
        call read_real(line, nsta, nend, qmmm%qm_dipole(1))
        nsta = 39
        nend = 51
        call read_real(line, nsta, nend, qmmm%qm_dipole(2))
        nsta = 53
        nend = 65
        call read_real(line, nsta, nend, qmmm%qm_dipole(3))
        found_dipole = .true.

      ! retrieve QM ESP charges
      !
      else if (index(line, "CHELPG Charges") > 0) then
        read(file1, *) line
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, '(a)') line
          nsta = 12
          nend = 26
          call read_real(line, nsta, nend, charge_tmp(i))          
        end do
        !
        ! Separate to QM and link atom groups
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          charge_tmp, 1)

        found_charge = .true.

      end if
      if (found_dipole .and. found_charge) exit

    end do
200 continue
    call close_file(file1)
    
    ! open engrad file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//".engrad", IOFileInput)
    do while (.true.)
      read(file1, '(a)', end=210) line

      ! retrieve QM energy
      !
      if (index(line, "The current total energy in Eh") > 0) then
        read(file1, *) line
        read(file1, *) qm_energy
        found_energy = .true.

      ! retrieve QM gradient
      !
      else if (index(line, "The current gradient in Eh/bohr") > 0) then
        read(file1, *) line
        nsta = 0
        do i = 1, qmmm%qm_natoms + qmmm%num_qmmmbonds
          read(file1, *) force_tmp(nsta+1)
          read(file1, *) force_tmp(nsta+2)
          read(file1, *) force_tmp(nsta+3)
          nsta = nsta + 3
        end do
        ! Separate to QM and link atom groups
        !
        call separate_qm_and_linkh(qmmm%qm_natoms, qmmm%num_qmmmbonds, &
          qmmm%qmatom_id, qmmm%linkatom_global_address, &
          force_tmp, 3)

        found_force = .true.

      end if
      if (found_energy .and. found_force) exit

    end do
210 continue
    call close_file(file1)

    ! open pcgrad file
    !
    call open_file(file1, trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//".pcgrad", IOFileInput)
    read(file1, *) nline

    ! retrieve MM gradient
    !
    do i = 1, nline
      read(file1, *) mm_force(1:3,i)
    end do
    mm_force(:,:) = -mm_force(:,:)
    found_mm_force = .true.
    call close_file(file1)

    ! Error check
    !
    qm_error = .not. found_energy
    if (qm_error) then
      write(MsgOut,'(a)') &
        'Fatal error in ORCA while reading' &
        //trim(qmmm%qmout)//'! Rank = ', my_country_no
      write(MsgOut,'(a)') 'QM energy is not found!'
      return
    end if

    if (.not. qmmm%ene_only) then
      qm_error = .not. (found_force .and. found_mm_force)
      if (qm_error) then
        write(MsgOut,'(a)') &
          'Fatal error in ORCA while reading' &
          //trim(qmmm%qmout)//'! Rank = ', my_country_no
        if (.not. found_force) write(MsgOut,'(a)') 'Force  is not found!'
        if (.not. found_mm_force) write(MsgOut,'(a)') 'MM-Force  is not found!'
        return
      end if

      nsta = 0
      do i = 1, qmmm%qm_natoms
        qm_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        qm_charge(i) = charge_tmp((nsta+3)/3)
        nsta = nsta + 3
      end do
      do i = 1, qmmm%num_qmmmbonds
        la_force(1:3,i) = -force_tmp(nsta+1:nsta+3)
        la_charge(i) = charge_tmp((nsta+3)/3)
        nsta = nsta + 3
      end do

    end if

    if (qmmm%qm_debug) then
      write(MsgOut,'("QMMM_debug> QM energy = ", f18.10, &
                     ", MM electrostatic energy = ", f18.10)') &
                      qm_energy, 0.0_wp
      write(MsgOut,'("QMMM_debug> QM dipole = ", 3f12.4)') &
                      qmmm%qm_dipole

    end if

    ! Save QM reults
    !
    if (qmmm%savefile) then

      if (qmmm%workdir /= qmmm%savedir) then

        ! save engrad file
        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.engrad ' & 
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

        ! save pcgrad file
        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.pcgrad ' &
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

        ! save xyz files
        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_qm.xyz ' & 
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

        exec='mv '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_pc.xyz ' &
                  //trim(qmmm%savedir)
        if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
        call system(exec)

      end if

    else

      exec='rm '//trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.engrad ' &
                //trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'.pcgrad ' &
                //trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_qm.xyz ' &
                //trim(qmmm%workdir)//'/'//trim(qmmm%qmbasename)//'_pc.xyz '
      if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> exec '',a)') trim(exec)
      call system(exec)

    end if

    return

  end subroutine qm_read_output_orca

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_qmmmbonds
  !> @brief        get number of QM-MM bonds
  !! @authors      YA
  !! @param[in]    molecule  : molecular information
  !! @param[inout] qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine count_qmmmbonds(molecule, qmmm)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    integer                  :: list_bond(10), list_bond_0(10)
    integer                  :: num_bond
    integer                  :: i, j, k, l, iatom, ibond


    ! initialization
    !
    qmmm%num_qmmmbonds = 0

    do i = 1, qmmm%qm_natoms

      iatom = qmmm%qmatom_id(i)

      ! find bonds involving QM atoms
      !
      num_bond = 0
      do ibond = 1, molecule%num_bonds
        if (iatom == molecule%bond_list(1,ibond)) then
          num_bond = num_bond + 1
          list_bond(num_bond) = molecule%bond_list(2,ibond)
        else if (iatom == molecule%bond_list(2,ibond)) then
          num_bond = num_bond + 1
          list_bond(num_bond) = molecule%bond_list(1,ibond)
        end if
      end do

      ! eliminate QM-QM bonds
      !
      list_bond_0(1:num_bond) = list_bond(1:num_bond)
      do j = 1, num_bond
        do k = 1, qmmm%qm_natoms
          if (list_bond(j) == qmmm%qmatom_id(k)) then
            list_bond_0(j) = -1
          end if
        end do
      end do

      l=0
      do j = 1, num_bond
         if (list_bond_0(j) /= -1) then
           l = l + 1
           list_bond(l) = list_bond_0(j) 
         end if
      end do
      num_bond = l

      qmmm%num_qmmmbonds = qmmm%num_qmmmbonds + num_bond
    
    end do

    return

  end subroutine count_qmmmbonds

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    generate_linkatoms
  !> @brief        generate link atom bond_list and exclude MM atoms
  !! @authors      YA, KY
  !! @param[in]    molecule  : molecular information
  !! @param[in]    psf       : PSF information
  !! @param[inout] qmmm      : QM/MM information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine generate_linkatoms(molecule, psf, qmmm)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_psf),             intent(in)    :: psf
    type(s_qmmm),            intent(inout) :: qmmm

    ! local variables
    real(wp)                 :: rr01, r0(3), r1(3), r01(3)
    integer                  :: list_bond(10), list_bond_0(10)
    integer                  :: num_bond, num_linkatoms
    integer                  :: i, j, k, l, ista, iend, iatom, ibond
    integer                  :: igrp, natom_grp, iatom_rm
    integer                  :: lkatm
    integer                  :: max_ecatoms, num_ecatoms, nec0
    integer                  :: kalloc_stat
    integer, parameter       :: num_ecatm_per_link = 20
    real(wp)                 :: total_charge, shift_charge
!    integer                  :: locsta, locend
!    real(wp)                 :: rr01, r0(3), r1(3), r01(3)


    if (main_rank) then
      write(MsgOut,'(/a)') "Generate_LinkAtoms> QM-MM interface info. Link hydrogen is set between:"
      write(MsgOut,'(14x,''[ QM atom ]'',13x,''-'',13x,''[ MM atom ]'')')
    end if

    if (qmmm%exclude_charge == ExcludeChargeAtom) then
      max_ecatoms = qmmm%num_qmmmbonds
    else if (qmmm%exclude_charge == ExcludeChargeGroup) then
      max_ecatoms = qmmm%num_qmmmbonds*num_ecatm_per_link
    else if (qmmm%exclude_charge == ExcludeChargeAtomAndDistribute) then
      max_ecatoms = qmmm%num_qmmmbonds
    end if
    allocate(qmmm%ecatom_id(max_ecatoms), stat = kalloc_stat)  
    if (kalloc_stat /= 0) call error_msg_alloc

    num_linkatoms = 0
    num_ecatoms = 1
    do i = 1, qmmm%qm_natoms

      iatom = qmmm%qmatom_id(i)

      ! find bonds involving QM atoms
      !
      num_bond = 0
      do ibond = 1, molecule%num_bonds
        if (iatom == molecule%bond_list(1,ibond)) then
          num_bond = num_bond + 1
          list_bond(num_bond) = molecule%bond_list(2,ibond)
        else if (iatom == molecule%bond_list(2,ibond)) then
          num_bond = num_bond + 1
          list_bond(num_bond) = molecule%bond_list(1,ibond)
        end if
      end do

      ! eliminate QM-QM bonds
      !
      list_bond_0(1:num_bond) = list_bond(1:num_bond)
      do j = 1, num_bond
        do k = 1, qmmm%qm_natoms
          if (list_bond(j) == qmmm%qmatom_id(k)) then
            list_bond_0(j) = -1
          end if
        end do
      end do

      l=0
      do j = 1, num_bond
         if (list_bond_0(j) /= -1) then
           l = l + 1
           list_bond(l) = list_bond_0(j) 
         end if
      end do
      num_bond = l

      ! Save QM-MM bond pair and punch out link atom info.
      !
      do j = 1, num_bond
        qmmm%qmmmbond_list(1,num_linkatoms+j) = iatom
        qmmm%qmmmbond_list(2,num_linkatoms+j) = list_bond(j)

        lkatm = iatom
        if (main_rank) &
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6'' - '',$)') &
          lkatm,                           &
          molecule%segment_name(lkatm),    &
          molecule%residue_no(lkatm),      &
          molecule%residue_name(lkatm),    &
          molecule%atom_name(lkatm),       &
          molecule%atom_cls_name(lkatm) 

        lkatm = list_bond(j)
        if (main_rank) &
        write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
          lkatm,                           &
          molecule%segment_name(lkatm),    &
          molecule%residue_no(lkatm),      &
          molecule%residue_name(lkatm),    &
          molecule%atom_name(lkatm),       &
          molecule%atom_cls_name(lkatm) 
      end do

      num_linkatoms = num_linkatoms + num_bond
    
      ! remove excluded MM charges
      !
      !if (qmmm%exclude_charge == ExcludeChargeAtom) then
      if (qmmm%exclude_charge == ExcludeChargeAtom .or. &
          qmmm%exclude_charge == ExcludeChargeAtomAndDistribute) then

        do j = 1, num_bond
          qmmm%ecatom_id(num_ecatoms) = list_bond(j)
          num_ecatoms = num_ecatoms + 1

          do k = 1, qmmm%mm_natoms
            if (list_bond(j) == qmmm%mmatom_id(k)) then
              if (k < qmmm%mm_natoms) then
                qmmm%mmatom_id(k:qmmm%mm_natoms-1) = &
                  qmmm%mmatom_id(k+1:qmmm%mm_natoms)
              end if
              qmmm%mm_natoms = qmmm%mm_natoms - 1

              exit
            end if
          end do
        end do

      else if (qmmm%exclude_charge == ExcludeChargeGroup) then

        do j = 1, num_bond

          do igrp = 1, psf%num_groups
            ista = psf%grp_list(1,igrp) + 1
            if (igrp < psf%num_groups) then
              iend = psf%grp_list(1,igrp+1) 
            else
              iend = molecule%num_atoms 
            end if
            natom_grp = iend - ista + 1 

            if (list_bond(j) >= ista .and. list_bond(j) <= iend) then 

              nec0 = num_ecatoms
              do iatom_rm = ista, iend
                do k = 1, qmmm%mm_natoms
                  if (iatom_rm == qmmm%mmatom_id(k)) then
                    if (k < qmmm%mm_natoms) then
                      qmmm%mmatom_id(k:qmmm%mm_natoms-1) = &
                        qmmm%mmatom_id(k+1:qmmm%mm_natoms)
                    end if
                    qmmm%mm_natoms = qmmm%mm_natoms - 1

                    qmmm%ecatom_id(num_ecatoms) = iatom_rm
                    num_ecatoms = num_ecatoms + 1
                    if (num_ecatoms > max_ecatoms) then
                      call error_msg('Generate_LinkAtoms> Too many excluded charge atoms.')
                    end if

                    exit  ! loop over mmatoms

                  end if
                end do
              end do

              ! Punch out excluded MM atom info.
              !
              do l = nec0, num_ecatoms - 1
                lkatm = qmmm%ecatom_id(l)
                if (lkatm == list_bond(j)) cycle
                if (main_rank) &
                  write(MsgOut,'(40x,i6,x,a4,x,i6,x,a4,x,a6,x,a6," excluded")')&
                               lkatm, molecule%segment_name(lkatm),            &
                               molecule%residue_no(lkatm),                     &
                               molecule%residue_name(lkatm),                   &
                               molecule%atom_name(lkatm),                      &
                               molecule%atom_cls_name(lkatm)
              end do

              exit  ! loop over group

            end if

          end do
        end do

      end if

    end do

    qmmm%ec_natoms = num_ecatoms - 1

    if (qmmm%exclude_charge == ExcludeChargeAtomAndDistribute) then

      ! total charge of the QM region
      !
      total_charge = 0.0_wp
      do i = 1, qmmm%qm_natoms
        total_charge = total_charge + molecule%charge(qmmm%qmatom_id(i))
      end do

      ! total charge of the excluded atoms
      !
      do i = 1, qmmm%ec_natoms
        total_charge = total_charge + molecule%charge(qmmm%ecatom_id(i))
      end do

      if (main_rank) then 
        write(MsgOut,*)
        write(MsgOut, '(a,f8.4, a, i6, a)') " Distributing charge ", &
        total_charge - qmmm%qm_total_charge, " over ", qmmm%mm_natoms, " MM atoms"
      end if

      ! externaly modify MM charges
      !
      shift_charge = (total_charge - qmmm%qm_total_charge) / dble(qmmm%mm_natoms)
      do i = 1, qmmm%mm_natoms
        j = qmmm%mmatom_id(i)
        molecule%charge(j) = molecule%charge(j) + shift_charge
      end do

      ! eliminate link MM atom charge
      !
      do i = 1, qmmm%ec_natoms
        molecule%charge(qmmm%ecatom_id(i)) = 0.0_wp
      end do
    end if

    if (main_rank) then
      write(MsgOut,'(/,a,i0)') "  number of link atoms           = ", &
        qmmm%num_qmmmbonds
      write(MsgOut,'(a,i0/)') "  number of external MM charges  = ", &
        qmmm%mm_natoms
    end if

    return

  end subroutine generate_linkatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_linkmm
  !> @brief        check that link MM atom is NOT hydrogen 
  !! @authors      KY
  !! @param[in]    molecule    : molecular information
  !! @param[inout] qmmm        : QM/MM information
  !! @param[inout] setup_error : true, when error is detected
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine check_linkmm(molecule, qmmm, setup_error)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    type(s_qmmm),     intent(inout) :: qmmm
    logical,          intent(inout) :: setup_error

    ! local variables
    real(wp) :: mm_mass
    integer  :: i, lkmm, ano
    logical  :: is_hydrogen


    is_hydrogen = .false.
    do i = 1, qmmm%num_qmmmbonds
      lkmm = qmmm%qmmmbond_list(2,i)
      mm_mass = molecule%mass(lkmm)
      call qm_atomic_number(mm_mass, ano)
      if (ano == 1) then
        is_hydrogen = .true.
        exit
      end if
    end do

    if (is_hydrogen .and. main_rank) then
      write(MsgOut,'(a)') &
        'ERROR>>> Hydrogen atom(s) is detected as link MM atom!'
      write(MsgOut,*)
      do i = 1, qmmm%num_qmmmbonds
        lkmm = qmmm%qmmmbond_list(2,i)
        mm_mass = molecule%mass(lkmm)
        call qm_atomic_number(mm_mass, ano)
        if (ano == 1) then
          write(MsgOut,'(i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
            lkmm,                         &
            molecule%segment_name(lkmm),  &
            molecule%residue_no(lkmm),    &
            molecule%residue_name(lkmm),  &
            molecule%atom_name(lkmm),     &
            molecule%atom_cls_name(lkmm)
        end if
      end do
      write(MsgOut,*)
      write(MsgOut,'(a)') &
       'ERROR>>> It is most likely that you forgot to include the &
       &hydrogen in a selection group for QM atoms.' 
      write(MsgOut,*)

      setup_error = .true.

    end if

    return

  end subroutine check_linkmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rearrange_term_list
  !> @brief        rearrange the list of terms in such a way that 
  !!                     list(:, 1:mm_terms)                   ... MM terms
  !!                     list(:, mm_terms:mm_terms + qm_terms) ... QM terms
  !! @authors      KY
  !! @param[in]    ityp      : type of force field terms
  !!                           2 = bond, 3 = angle, 4 = dihed, CMAP
  !!                           5 = imprpr
  !! @param[in]    qm_natoms : number of QM atoms
  !! @param[in]    qmatom_id : ID of QM atoms
  !! @param[out]   qm_terms  : number of terms in QM region
  !! @param[in]    num_terms : number of all terms 
  !! @param[in]    forcefield : forcefield identifier
  !! @param[out]   list      : list of terms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine rearrange_term_list(ityp, qm_natoms, qmatom_id, &
                                 qm_terms, num_terms, forcefield, list)

    ! formal arguments
    integer, intent(in)    :: ityp
    integer, intent(in)    :: qm_natoms
    integer, intent(in)    :: qmatom_id(qm_natoms)
    integer, intent(out)   :: qm_terms
    integer, intent(in)    :: num_terms
    integer, intent(in)    :: forcefield
    integer, intent(inout) :: list(ityp,num_terms)

    ! local variables
    integer                :: i, n, mm_terms
    integer                :: j, jj, jmax
    integer                :: kalloc_stat, kdealloc_stat

    integer,   allocatable :: qm_index(:)
    integer,   allocatable :: mm_index(:)
    integer,   allocatable :: qm_list (:,:)
    integer,   allocatable :: jlist(:)


    if (ityp == 2) then
       jmax = 2
       allocate(jlist(2), stat = kalloc_stat)
       if (kalloc_stat /= 0) call error_msg_alloc
       jlist(1) = 1
       jlist(2) = 2

    else if (ityp == 3) then
       jmax = 3
       allocate(jlist(3), stat = kalloc_stat)
       if (kalloc_stat /= 0) call error_msg_alloc
       jlist(1) = 1
       jlist(2) = 2
       jlist(3) = 3

    else if (mod(ityp, 4) == 0) then
      if (forcefield == ForcefieldCHARMM) then
        jmax = 2
        allocate(jlist(2), stat = kalloc_stat)
        if (kalloc_stat /= 0) call error_msg_alloc
        jlist(1) = 2
        jlist(2) = 3
      else if (forcefield == ForcefieldAMBER) then
        jmax = 4
        allocate(jlist(4), stat = kalloc_stat)
        if (kalloc_stat /= 0) call error_msg_alloc
        jlist(1) = 1
        jlist(2) = 2
        jlist(3) = 3
        jlist(4) = 4
      end if

    end if

    allocate(qm_index(num_terms), mm_index(num_terms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    qm_terms = 0
    mm_terms = 0

    do n = 1, num_terms

      jj = 0
      do i = 1, qm_natoms
        do j = 1, jmax
          if (list(jlist(j),n) == qmatom_id(i)) then
            jj = jj + 1
          end if
        end do
        if (jj == jmax) exit
      end do

      if (jj == jmax) then
        qm_terms = qm_terms + 1
        qm_index(qm_terms) = n

      else
        mm_terms = mm_terms + 1
        mm_index(mm_terms) = n

      end if

    end do

    allocate(qm_list(ityp, qm_terms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    do n = 1, qm_terms
      qm_list(:,n) = list(:,qm_index(n))
    end do

    do n = 1, mm_terms
      list(:,n) = list(:,mm_index(n))
    end do

    i = 1
    do n = mm_terms+1, num_terms
      list(:,n) = qm_list(:,i)
      i = i + 1
    end do

    deallocate(jlist, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(qm_index, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(mm_index, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(qm_list, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine rearrange_term_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rearrange_impr_list
  !> @brief        rearrange the list of terms in such a way that 
  !!                     list(:, 1:mm_terms)                   ... MM terms
  !!                     list(:, mm_terms:mm_terms + qm_terms) ... QM terms
  !! @authors      KY
  !! @param[in]    qm_natoms : number of QM atoms
  !! @param[in]    qmatom_id : ID of QM atoms
  !! @param[out]   qm_terms  : number of terms in QM region
  !! @param[in]    num_terms : number of all terms 
  !! @param[in]    forcefield : forcefield identifier
  !! @param[out]   list      : list of terms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine rearrange_impr_list(qm_natoms, qmatom_id, &
                                 qm_terms, num_terms, forcefield, list)

    ! formal arguments
    integer, intent(in)    :: qm_natoms
    integer, intent(in)    :: qmatom_id(qm_natoms)
    integer, intent(out)   :: qm_terms
    integer, intent(in)    :: num_terms
    integer, intent(in)    :: forcefield
    integer, intent(inout) :: list(4,num_terms)

    ! local variables
    integer                :: i, n, mm_terms
    integer                :: j, jj, jmax
    integer                :: kalloc_stat, kdealloc_stat

    integer, allocatable   :: qm_index(:)
    integer, allocatable   :: mm_index(:)
    integer, allocatable   :: qm_list (:,:)
    !integer, allocatable   :: jlist(:)


    allocate(qm_index(num_terms), mm_index(num_terms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc

    if (forcefield == ForcefieldCHARMM) then
      jmax = 3
    else if (forcefield == ForcefieldAMBER) then
      jmax = 4
    end if
    qm_terms = 0
    mm_terms = 0

    do n = 1, num_terms

      jj = 0
      do i = 1, qm_natoms
        do j = 1, 4
          if (list(j,n) == qmatom_id(i)) then
            jj = jj + 1
          end if
        end do
        if (jj >= jmax) exit
      end do

      if (jj >= jmax) then
        qm_terms = qm_terms + 1
        qm_index(qm_terms) = n

      else
        mm_terms = mm_terms + 1
        mm_index(mm_terms) = n

      end if

    end do

    allocate(qm_list(4, qm_terms), stat = kalloc_stat)
    if (kalloc_stat /= 0) call error_msg_alloc
    do n = 1, qm_terms
      qm_list(:,n) = list(:,qm_index(n))
    end do

    do n = 1, mm_terms
      list(:,n) = list(:,mm_index(n))
    end do

    i = 1
    do n = mm_terms+1, num_terms
      list(:,n) = qm_list(:,i)
      i = i + 1
    end do

    deallocate(qm_index, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(mm_index, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc
    deallocate(qm_list, stat = kdealloc_stat)
    if (kdealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine rearrange_impr_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    merge_qm_and_linkh
  !> @brief        Merge QM atom list and link-H list
  !! @authors      YA
  !! @param[in]    qm_natoms : number of QM atoms
  !! @param[in]    num_linkh : number of link hydrogens
  !! @param[in]    qm_address : QM atom global addresses
  !! @param[in]    lh_address : link atom global addresses
  !! @param[inout] array      : array to sort
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine merge_qm_and_linkh(qm_natoms, num_linkh, qm_address, &
                                lh_address, array)

    ! formal arguments
    integer,                 intent(in)    :: qm_natoms
    integer,                 intent(in)    :: num_linkh
    integer,                 intent(in)    :: qm_address(qm_natoms)
    integer,                 intent(in)    :: lh_address(num_linkh)
    character(128),          intent(inout) :: array(*)

    ! Local variables
    integer        :: address(qm_natoms+num_linkh)
    integer        :: idx(1,qm_natoms+num_linkh)
    integer        :: i, minv, cnt, idum(1)
    character(128) :: tmp(qm_natoms+num_linkh)


    cnt = 0
    do i = 1, qm_natoms
      cnt = cnt + 1
      address(cnt) = qm_address(i)
    end do
    do i = 1, num_linkh
      cnt = cnt + 1
      address(cnt) = lh_address(i)
    end do
    
    !
    minv = 0
    do i = 1, qm_natoms + num_linkh
      idx(:,i) = MINLOC(address, MASK=address > minv)
      minv = MINVAL(address, MASK=address > minv)
    end do
    
    !
    do i = 1, qm_natoms + num_linkh
      tmp(i) = array(idx(1,i))
    end do

    array(1:qm_natoms+num_linkh) = tmp(1:qm_natoms+num_linkh)
    
    return

  end subroutine merge_qm_and_linkh

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    separate_qm_and_linkh
  !> @brief        Separate atom list into QM atoms and link hydrogens
  !! @authors      YA
  !! @param[in]    qm_natoms : number of QM atoms
  !! @param[in]    num_linkh : number of link hydrogens
  !! @param[in]    qm_address : QM atom global addresses
  !! @param[inout] lh_address : link atom global addresses
  !! @param[in]    array      : array to separate
  !! @param[in]    ld         : leading dimension of array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine separate_qm_and_linkh(qm_natoms, num_linkh, qm_address, &
                                   lh_address, array, ld)

    ! formal arguments
    integer,                 intent(in)    :: qm_natoms
    integer,                 intent(in)    :: num_linkh
    integer,                 intent(in)    :: qm_address(qm_natoms)
    integer,                 intent(in)    :: lh_address(num_linkh)
    integer,                 intent(in)    :: ld
    real(wp),                intent(inout) :: array(ld,*)

    ! Local variables
    real(wp) :: tmp(ld,qm_natoms+num_linkh)
    integer  :: address(qm_natoms+num_linkh)
    integer  :: idx(1,qm_natoms+num_linkh)
    integer  :: i, minv, cnt


    cnt = 0
    do i = 1, qm_natoms
      cnt = cnt + 1
      address(cnt) = qm_address(i)
    end do
    do i = 1, num_linkh
      cnt = cnt + 1
      address(cnt) = lh_address(i)
    end do
    
    !
    minv = 0
    do i = 1, qm_natoms + num_linkh
      idx(:,i) = MINLOC(address, MASK=address > minv)
      minv = MINVAL(address, MASK=address > minv)
    end do
    
    !
    do i = 1, qm_natoms + num_linkh
      tmp(1:ld,idx(1,i)) = array(1:ld,i)
    end do

    array(1:ld,1:qm_natoms+num_linkh) = tmp(1:ld,1:qm_natoms+num_linkh)
    
    return

  end subroutine separate_qm_and_linkh

end module at_qmmm_mod
