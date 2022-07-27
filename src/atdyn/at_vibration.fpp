!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_vibration_mod
!> @brief   perform vibrational analysis
!! @authors Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_vibration_mod

  use at_dynvars_str_mod
  use at_enefunc_str_mod
  use at_boundary_str_mod
  use at_output_str_mod
  use at_pairlist_str_mod
  use at_vibration_str_mod
  use at_energy_mod
  use at_input_mod
  use at_output_mod
  use at_qmmm_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_control_mod
  use fileio_minfo_mod
  use select_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_vib_info

    ! For vibrational analysis
    integer          :: runmode             = RunModeHarm
    integer          :: nreplica            = 1
    character(256)   :: vibatm_select_index = ''
    character(256)   :: output_minfo_atm    = ''
    character(256)   :: minfo_folder        = 'minfo.files'
    real(wp)         :: diff_stepsize       = 0.01_wp
    character(256)   :: gridfile            = ''
    character(256)   :: datafile            = ''
!    real(wp)         :: cutoff              = -1.0

  end type s_vib_info

  ! subroutines
  public  :: show_ctrl_vibration
  public  :: read_ctrl_vibration
  public  :: setup_vibration
  private :: setup_vibatoms
  public  :: run_vib
  private :: harmonic_analysis
  private :: calc_hessian
  private :: normal_mode
  private :: match_vibatom
  private :: diagonalize
  private :: print_minfo
  private :: print_minfo_grad
  private :: read_minfo_grad
  private :: calc_gridpoints

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_vibration
  !> @brief        show usage of VIBRATION section
  !! @authors      KY
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !!                          "vib"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_vibration(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('vib')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') 'group1              = atomno:5-8  # Vib group 1'
        write(MsgOut,'(A)') 'group2              = sid:PROA    # Vib group 2'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[VIBRATION]'
        write(MsgOut,'(A)') 'nreplica            = 1           # number of MPI processes' 
        write(MsgOut,'(A)') 'runmode             = HARM        # HARM, QFF, or GRID'
        write(MsgOut,'(A)') 'vibatm_select_index = 1           # atoms subject to vib analysis'
        write(MsgOut,'(A)') 'output_minfo_atm    = 2           # atoms punched to a minfo file'
        write(MsgOut,'(A)') 'minfo_folder        = minfo.files # a folder where minfo files are created'
        write(MsgOut,'(A)') '# diff_stepsize     = 0.01        # displacement for numerical diff.'
        write(MsgOut,'(A)') '# gridfile          = grid.xyz    # the xyz file containing coordinates of grid points'
        write(MsgOut,'(A)') '# datafile          = grid.dat    # the output file for grid data'
        !write(MsgOut,'(A)') 'cutoff              = -1.0      # cutoff distance (in Angs) for Hessian evaluation.'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('vib')

        write(MsgOut,'(A)') '[SELECTION]'
        write(MsgOut,'(A)') 'group1              = atomno:5-8  # Vib group 1'
        write(MsgOut,'(A)') 'group2              = sid:PROA    # Vib group 2'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '[VIBRATION]'
        write(MsgOut,'(A)') 'nreplica            = 1           # number of MPI processes' 
        write(MsgOut,'(A)') 'runmode             = HARM        # HARM, QFF, or GRID'
        write(MsgOut,'(A)') 'vibatm_select_index = 1           # fixed atoms in minimization'
        write(MsgOut,'(A)') 'output_minfo_atm    = 2           # atoms punched to a minfo file'
        write(MsgOut,'(A)') 'minfo_folder        = minfo.files # a folder where minfo files are created'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_vibration
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_vibration
  !> @brief        read VIBRATIOn section in the control file
  !! @authors      KY
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   vib_info : MINIMIZE section in control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_vibration(handle, vib_info)

    ! parameters
    character(*),            parameter     :: Section = 'Vibration'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_vib_info),        intent(inout) :: vib_info

    ! local
    integer :: idx
    logical :: found_error


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    ! Vibrational analysis
    call read_ctrlfile_type   (handle, Section, 'runmode',              &
                               vib_info%runmode, RunModeTypes)
    call read_ctrlfile_integer(handle, Section, 'nreplica',             &
                               vib_info%nreplica)
    call read_ctrlfile_string (handle, Section, 'vibatm_select_index',  &
                               vib_info%vibatm_select_index)
    call read_ctrlfile_string (handle, Section, 'output_minfo_atm',     &
                               vib_info%output_minfo_atm)
    call read_ctrlfile_real   (handle, Section, 'diff_stepsize',        &
                               vib_info%diff_stepsize)
    call read_ctrlfile_string (handle, Section, 'minfo_folder',          &
                               vib_info%minfo_folder)
!    call read_ctrlfile_real   (handle, Section, 'cutoff',               &
!                               vib_info%cutoff)
    call read_ctrlfile_string (handle, Section, 'gridfile',             &
                               vib_info%gridfile)
    call read_ctrlfile_string (handle, Section, 'datafile',             &
                               vib_info%datafile)

    call end_ctrlfile_section(handle)

    if (trim(vib_info%vibatm_select_index) .eq. '') &
      call error_msg('Read_Ctrl_Vibration> No VIB atoms defined')

    select case (vib_info%runmode)
    case (RunModeQFF)

      ! Generate QFF
      if (vib_info%gridfile == '') vib_info%gridfile = 'makeQFF.xyz'

    case (RunModeGRID)
      ! Generate Grid-PES
      if (vib_info%gridfile == '') vib_info%gridfile = 'makeGrid.xyz'
      if (vib_info%datafile == '') then
        idx=index(vib_info%gridfile,'.xyz')
        vib_info%datafile = vib_info%gridfile(1:idx-1)//'.dat'
      end if

    end select

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Vibration> Parameters of VIBRATION'
      write(MsgOut,'(A,A10)')    &
          '  runmode             = ', trim(RunModeTypes(vib_info%runmode))
      write(MsgOut,'(A,I10)')    &
          '  nreplica            = ', vib_info%nreplica
      write(MsgOut,'(A,A)')      &
          '  vibatm_select_index = ', trim(vib_info%vibatm_select_index)
      if (vib_info%output_minfo_atm /= '') write(MsgOut,'(A,A)')      &
          '  output_minfo_atm    = ', trim(vib_info%output_minfo_atm)
      write(MsgOut,'(A,A)')      &
          '  minfo_folder        = ', trim(vib_info%minfo_folder)

      select case (vib_info%runmode)

      case (RunModeHARM)

        ! Harmonic vibrational analysis
        write(MsgOut,'(A,E10.2)')&
          '  diff_stepsize       = ', vib_info%diff_stepsize
        ! write(MsgOut,'(A,F10.5)')       &
        !   '  cutoff              = ', vib_info%cutoff

      case (RunModeQFF)

        write(MsgOut,'(A,A)')    &
          '  gridfile            = ', trim(vib_info%gridfile)

      case (RunModeGRID)

        write(MsgOut,'(A,A)')    &
          '  gridfile            = ', trim(vib_info%gridfile)
        write(MsgOut,'(A,A)')    &
          '  datafile            = ', trim(vib_info%datafile)

      end select

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_vibration

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_vibration
  !> @brief        setup vibration information
  !> @authors      KY
  !! @param[in]    inp_info  : INPUT section in control parameters
  !! @param[in]    vib_info  : VIBRATION section in control parameters
  !! @param[in]    sel_info  : SELECTION section in control parameters
  !! @param[inout] molecule  : molecular information
  !! @param[in]    qmmm      : QM/MM information
  !! @param[out]   vibration : vibration information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_vibration(inp_info, vib_info, sel_info, &
                             molecule, coord, qmmm, vibration)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_vib_info),        intent(in)    :: vib_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(inout) :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_vibration),       intent(inout) :: vibration


    call setup_vibatoms(vib_info, sel_info, molecule, vibration)

    vibration%nreplica         = vib_info%nreplica
    vibration%diff_stepsize    = vib_info%diff_stepsize
    vibration%minfo_folder     = vib_info%minfo_folder

    select case (vib_info%runmode)
    case(RunModeHARM)
      vibration%gengrid          = .false.
      vibration%grid_ene_only    = .false.
      !vibration%cutoff           = vib_info%cutoff

    case(RunModeQFF )
      vibration%gengrid          = .true.
      vibration%grid_ene_only    = .false.

    case(RunModeGRID)
      vibration%gengrid          = .true.
      vibration%grid_ene_only    = .true.

    end select

    vibration%gridfile         = vib_info%gridfile
    vibration%datafile         = vib_info%datafile

    if (inp_info%minfofile /= '') &
      call read_minfo(inp_info%minfofile, vibration%minfo_in, molecule, coord)

    if (qmmm%do_qmmm .and. qmmm%qmmaxtrial > 0 .and. &
      (vib_info%runmode == RunModeHARM .or. vib_info%runmode == RunModeQFF)) then

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Vibration> '
        write(MsgOut,'(A,i4)') &
          '  WARNING: qmmaxtrial in [QMMM] is non-zero: qmmaxtrial = ', &
          qmmm%qmmaxtrial
        write(MsgOut,'(A)') &
          '  WARNING: QM retry causes problems in numerical differentiations, and thus'
        write(MsgOut,'(A)') &
          '  WARNING: not recommended. If you encounter SCF convergence problems,'
        write(MsgOut,'(A)') &
          '  WARNING: try instead with increased SCF maximum iteration in a QM input.'
        write(MsgOut,'(A,i4)') &
          '  WARNING: qmmaxtrial is reset to zero.'
        write(MsgOut,*)
      end if
      qmmm%qmmaxtrial = 0
    end if

    return

  end subroutine setup_vibration

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_vibatoms
  !> @brief        define atoms subject to vibrational analysis
  !! @authors      KY
  !! @param[in]    vib_info    : vibration input information
  !! @param[in]    sel_info    : selector input information
  !! @param[in]    molecule    : molecular information
  !! @param[out]   vibration    : vibration parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_vibatoms(vib_info, sel_info, molecule, vibration)

    ! formal arguments
    type(s_vib_info),        intent(in)    :: vib_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_vibration),       intent(inout) :: vibration

    ! local variables
    integer                :: igroup, ng_vib, ng_minfo, ng_minfo_bk, natom
    integer                :: i, j, nn, offset, temp, iatom
    integer, allocatable   :: glist_vib(:), glist_minfo(:), glist_minfo_bk(:)
    type(s_selatoms), allocatable :: selatoms(:)

    integer, parameter :: max_atm_print = 100


    ! Number of atoms for vibrational analysis
    !
    ng_vib = split_num(trim(vib_info%vibatm_select_index))

    allocate(glist_vib(ng_vib))
    call split(ng_vib, ng_vib, vib_info%vibatm_select_index, glist_vib)

    allocate(selatoms(ng_vib))

    natom = 0
    do i = 1, ng_vib
      igroup = glist_vib(i)
      call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
      natom = natom + size(selatoms(i)%idx)
    end do

    vibration%vib_natoms = natom

    ! List of vibanal atoms
    !
    allocate(vibration%vibatom_id(vibration%vib_natoms))

    offset = 0
    do i = 1, ng_vib
      igroup = glist_vib(i)
      natom = size(selatoms(i)%idx)
      vibration%vibatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
      offset = offset + natom
    end do

    deallocate(selatoms)

    ! sort vibanal atom indices in ascending order
    !
    do i = vibration%vib_natoms, 2, -1
      do j = 1, i - 1
        if (vibration%vibatom_id(j) > vibration%vibatom_id(j+1)) then
          temp = vibration%vibatom_id(j)
          vibration%vibatom_id(j)   = vibration%vibatom_id(j+1)
          vibration%vibatom_id(j+1) = temp
        end if
      end do
    end do

    ! Number of atoms for minfo
    !
    ng_minfo = split_num(trim(vib_info%output_minfo_atm))

    if (ng_minfo > 0) then
      allocate(glist_minfo(ng_minfo))
      call split(ng_minfo, ng_minfo, vib_info%output_minfo_atm, glist_minfo)

      allocate(selatoms(ng_minfo))

      natom = 0
      do i = 1, ng_minfo
        igroup = glist_minfo(i)
        call select_atom(molecule, sel_info%groups(igroup), selatoms(i))
        natom = natom + size(selatoms(i)%idx)
      end do

      vibration%minfo_natoms = natom

      ! List of minfo subatoms
      !
      allocate(vibration%minfoatom_id(vibration%minfo_natoms))

      offset = 0
      do i = 1, ng_minfo
        igroup = glist_minfo(i)
        natom = size(selatoms(i)%idx)
        vibration%minfoatom_id(offset+1:offset+natom) = selatoms(i)%idx(1:natom)
        offset = offset + natom
      end do

      deallocate(selatoms)

      ! sort atom indices in ascending order
      !
      do i = vibration%minfo_natoms, 2, -1
        do j = 1, i - 1
          if (vibration%minfoatom_id(j) > vibration%minfoatom_id(j+1)) then
            temp = vibration%minfoatom_id(j)
            vibration%minfoatom_id(j)   = vibration%minfoatom_id(j+1)
            vibration%minfoatom_id(j+1) = temp
          end if
        end do
      end do

    else
      vibration%minfo_natoms = 0

    end if

    ! Punch out vibatom info.
    !
    if (main_rank) then
      write(MsgOut,'(a)') &
        "Setup_Vibration_Atoms> Atoms subject to vibrational analysis"

      do i = 1, vibration%vib_natoms
        iatom   = vibration%vibatom_id(i)
        write(MsgOut,'(2i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
          i,                             &
          iatom,                         &
          molecule%segment_name(iatom),  &
          molecule%residue_no(iatom),    &
          molecule%residue_name(iatom),  &
          molecule%atom_name(iatom),     &
          molecule%atom_cls_name(iatom)
      end do

      write(MsgOut,'(a,i0)') "  number of VIB atoms = ", vibration%vib_natoms
      write(MsgOut, '(a)') ' '

      if (vibration%minfo_natoms > 0) then
        write(MsgOut,'(a)') &
          "Setup_Vibration_Atoms> Atoms punched to minfo file in addition"

        if (vibration%minfo_natoms < max_atm_print) then
          do i = 1, vibration%minfo_natoms
            iatom   = vibration%minfoatom_id(i)
            write(MsgOut,'(2i6,x,a4,x,i6,x,a4,x,a6,x,a6)') &
              i,                             &
              iatom,                         &
              molecule%segment_name(iatom),  &
              molecule%residue_no(iatom),    &
              molecule%residue_name(iatom),  &
              molecule%atom_name(iatom),     &
              molecule%atom_cls_name(iatom)
          end do
        end if

        write(MsgOut,'(a,i0)') "  number of atoms = ", vibration%minfo_natoms
        write(MsgOut, '(a)') ' '
      end if
    end if

    deallocate(glist_vib)

    return

  end subroutine setup_vibatoms

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_vib
  !> @brief        perform vibrational analysis
  !! @authors      YA, KY
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] vibration   : vibration information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_vib(molecule, enefunc, dynvars, vibration, &
                     output, pairlist, boundary)

    ! formal arguments
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_vibration),        intent(inout) :: vibration
    type(s_output),           intent(inout) :: output
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary


    ! Open output files
    !
    call open_output(output)

    if (.not. vibration%gengrid) then
      ! Harmonic vibrational analysis
      call harmonic_analysis(molecule, enefunc, dynvars, vibration, &
                        output, pairlist, boundary)
    else
      ! Calc ene/grad at grid points
      call calc_gridpoints(molecule, enefunc, dynvars, vibration, &
                        output, pairlist, boundary)
    end if

    ! close output files
    !
    call close_output(output)

    return

  end subroutine run_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    harmonic_analysis
  !> @brief        harmonic vibrational analysis
  !> @authors      KY
  !! @param[inout] molecule  : molecule information
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[inout] dynvars   : dynamic variables information
  !! @param[inout] vibration : vibration information
  !! @param[inout] output    : output information
  !! @param[inout] pairlist  : non-bond pair list information
  !! @param[inout] boundary  : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine harmonic_analysis(molecule, enefunc, dynvars, vibration, &
                          output, pairlist, boundary)

    ! formal arguments
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_vibration), target, intent(inout) :: vibration
    type(s_output),            intent(inout) :: output
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    ! local variables
    integer :: nat, nat3
    real(wp)               :: energy
    real(wp), allocatable  :: grad(:)
    real(wp), allocatable  :: hessian(:,:)
    real(wp), allocatable  :: dipole(:)
    real(wp), allocatable  :: dipole_derivative(:,:)

    integer   :: i, iatom


    ! Print message
    if (main_rank) then
      write(MsgOut,'(''Enter vibrational analysis'')')
      write(MsgOut,*)

      write(MsgOut,'(''  Cartesian coordinates of atoms for vib. analysis'')')
      do i = 1, vibration%vib_natoms
        iatom = vibration%vibatom_id(i)
        write(MsgOut,'(2x,i4,2x,i9,x,a4,x,i6,x,a4,x,a6,4x,3f18.10)') &
          i, iatom, molecule%segment_name(iatom), molecule%residue_no(iatom), &
          molecule%residue_name(iatom), molecule%atom_name(iatom), &
          dynvars%coord(1:3,iatom)
      end do
      write(MsgOut,*)
    end if

    nat  = vibration%vib_natoms
    nat3 = nat*3
    allocate(grad(nat3))
    allocate(hessian(nat3,nat3))
    allocate(dipole(3))
    allocate(dipole_derivative(3,nat3))

    call calc_hessian(molecule, enefunc, dynvars, vibration, &
                      output, pairlist, boundary,            &
                      nat, energy, grad, hessian, dipole, dipole_derivative)

    call normal_mode(molecule, dynvars, vibration, &
                      nat, energy, grad, hessian, dipole, dipole_derivative)

    deallocate(grad, hessian, dipole, dipole_derivative)

    return

  end subroutine harmonic_analysis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_hessian
  !> @brief        Calculate Hessian matrix by numerical differentiations
  !> @authors      KY
  !! @param[inout] molecule  : molecule information
  !! @param[inout] enefunc   : potential energy functions information
  !! @param[inout] dynvars   : dynamic variables information
  !! @param[inout] vibration : vibration information
  !! @param[inout] output    : output information
  !! @param[inout] pairlist  : non-bond pair list information
  !! @param[inout] boundary  : boundary conditions information
  !! @param[in]    nat         : number of atoms
  !! @param[in]    energy      : total energy at the current geometry
  !! @param[in]    grad(nat3)  : gradient
  !! @param[in]    hessian(nat3, nat3) : Hessian matrix
  !! @param[in]    dipole(3)   : dipole moment
  !! @param[in]    dipole_derv(3,nat3) : dipole moment
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_hessian(molecule, enefunc, dynvars, vibration, output,  &
      pairlist, boundary, nat, energy, grad, hessian, dipole, dipole_derv)

    ! formal arguments
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_vibration), target, intent(inout) :: vibration
    type(s_output),            intent(inout) :: output
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    integer  :: nat
    real(wp) :: energy
    real(wp) :: grad(nat*3)
    real(wp) :: hessian(nat*3,nat*3)
    real(wp) :: dipole(3)
    real(wp) :: dipole_derv(3,nat*3)

    ! local variables
    integer :: nreplica, replicaid
    integer :: nat3, nat_calc, ngrid
    integer :: count
    integer :: icyc, istart0, istart, iend
    integer :: icalc, ist1, ist2, isiz1, isiz2
    integer :: iatom, ixyz
    integer :: i, j, ij, n, ierr
    integer :: i1, i2, k1, k2

    real(wp), pointer      :: coord(:,:), mass(:)

    real(wp), allocatable  :: energy_all(:),   recv_energy(:)
    real(wp), allocatable  :: grad_all(:,:),   recv_grad(:,:)
    real(wp), allocatable  :: dipole_all(:,:), recv_dipole(:,:)
    real(wp), allocatable  :: hess_tmp(:,:), dd_tmp(:,:)
    integer,  allocatable  :: vibatom_calc_id(:)
    logical,  allocatable  :: vibatom_calc(:)

    character(1)           :: xyz(3)
    character(MaxFilename) :: fname
    logical                :: ex

    real(wp)  :: rmsg
    real(wp)  :: delta, x0


    ! use pointers
    coord => dynvars%coord
    mass  => molecule%mass

    ! replica 
    nreplica  = vibration%nreplica
    replicaid = my_country_no + 1

    ! setup vibatoms to calc. the gradient
    allocate(vibatom_calc_id(nat), vibatom_calc(nat))
    vibatom_calc_id = vibration%vibatom_id
    vibatom_calc    = .true.
    nat_calc        = nat

    if (vibration%minfo_in%nat > 0) then
      call match_vibatom(vibration%minfo_in, &
                         nat_calc, vibatom_calc_id, hessian, dipole_derv)
      if (nat_calc < nat) then
        vibatom_calc = .false.
        do i = 1, nat_calc
          do j = 1, nat
            if (vibration%vibatom_id(j) == vibatom_calc_id(i)) then
              vibatom_calc(j) = .true.
              exit
            end if
          end do
        end do
        !dbg write(MsgOut,'("nat_calc = ",i5)') nat_calc
        !dbg write(MsgOut,'(i5)') vibatom_calc_id(1:nat_calc)
        !dbg write(MsgOut,'(l5)') vibatom_calc
      end if
    end if

    ! number of grid points
    ngrid           = nat_calc*6 + 1

    nat3 = nat*3
    allocate(energy_all(0:ngrid-1))
    allocate(grad_all(nat3,0:ngrid-1))
    allocate(dipole_all(3,0:ngrid-1))

    ! for print 
    xyz(1) = 'X'
    xyz(2) = 'Y'
    xyz(3) = 'Z'

    ! initialize
    fname        = ""

    ! create a folder to save minfo files
    call system('mkdir -p '//trim(vibration%minfo_folder)//' > /dev/null 2>&1')

    ! read minfo files 
    istart0 = -1
    do icyc = 0, ngrid-1

      call setaddress(icyc, vibration%minfo_folder, iatom, ixyz, fname)
      inquire(file=trim(fname), exist = ex)
      if (ex) then
        call read_minfo_grad(fname, nat, energy_all(icyc), &
                                         grad_all(:,icyc), &
                                         dipole_all(:,icyc))
      else
        istart0 = icyc
        exit
      end if

    end do

    if (istart0 < 0) then
      ! all grid points are calcualted
      goto 1000
    end if

    if (main_rank) then
      write(MsgOut,'(''  Generate Hessian matrix by num. diff. of gradients  '')')
      write(MsgOut,'(''    Loop over atoms'')')
    end if

#ifdef HAVE_MPI_GENESIS
    ! parallelize over replica
    i = ngrid - istart0
    j = mod(i, nreplica)
    if (j == 0) then
      i = i/nreplica
      istart = istart0 + i * (replicaid - 1)
      iend   = istart  + i - 1

    else
      i = (i-j)/nreplica
      j = nreplica - j
      if (replicaid <= j) then
        istart = istart0 + i * (replicaid -1)
        iend   = istart  + i
      else
        istart = istart0 + i*j + (i+1) * (replicaid -j-1)
        iend   = istart  + i
      end if

    end if
#else
    iend = ngrid-1
#endif

    !dbg write(MsgOut,'("replica = ",i4,", istart = ",i4,", iend = ",i4)') &
    !dbg   replicaid, istart, iend

    ! now compute_energy over grid points
    count = 0
    dynvars%step = count
    do icyc = istart, iend

      call setaddress(icyc, vibration%minfo_folder, iatom, ixyz, fname)

      if (icyc /= 0) then
        delta = vibration%diff_stepsize/sqrt(mass(iatom))
        x0 = coord(ixyz,iatom)
        if (mod(icyc, 2) == 1) then
          coord(ixyz,iatom) = x0 + delta
        else
          coord(ixyz,iatom) = x0 - delta
        end if
      end if

      call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                          .false.,                &
                          dynvars%coord,          &
                          dynvars%trans,          &
                          dynvars%coord_pbc,      &
                          dynvars%energy,         &
                          dynvars%temporary,      &
                          dynvars%force,          &
                          dynvars%force_omp,      &
                          dynvars%virial,         &
                          dynvars%virial_extern)

      ! print
      !
      if (replica_main_rank) then
        if (icyc == 0) then
          write(MsgOut,'(6x,"Done for    0   input",7x,"replicaID = ",i5, &
                         6x,"energy = ",f20.6)') &
                replicaid, dynvars%energy%total 
        else
          write(MsgOut,'(6x,"Done for ",i4,x,a6,"+",a1,6x,"replicaID = ",i5)') &
                icyc, molecule%atom_name(iatom), xyz(ixyz), replicaid
        end if
      end if
      call output_vib(output, molecule, enefunc, vibration, boundary, dynvars)

      ! increment
      !
      if (icyc /= 0) coord(ixyz,iatom) = x0
      count = count + 1
      dynvars%step = count

      ! set energy, grad, and dipole, and print minfo
      !
      energy_all(icyc) =  dynvars%energy%total / CONV_UNIT_ENE
      ij = 1
      do i = 1, nat
        do j = 1, 3
          grad_all(ij,icyc) = -dynvars%force(j,vibration%vibatom_id(i)) / CONV_UNIT_FORCE
          ij = ij + 1
        end do
      end do
      dipole_all(:,icyc) = enefunc%qmmm%qm_dipole

      call print_minfo_grad(fname, nat, energy_all(icyc), grad_all(:,icyc), dipole_all(:,icyc))

    end do

#ifdef HAVE_MPI_GENESIS
    icalc = iend - istart + 1
    allocate(recv_energy(0:icalc*nreplica-1))
    allocate(recv_grad(nat3,0:icalc*nreplica-1))
    allocate(recv_dipole(3,0:icalc*nreplica-1))

    call mpi_allgather(energy_all(istart), icalc, mpi_wp_real, &
                           recv_energy(0), icalc, mpi_wp_real, &
                           mpi_comm_airplane, ierr)
    call mpi_allgather(grad_all(1,istart), icalc*nat3, mpi_wp_real, &
                           recv_grad(1,0), icalc*nat3, mpi_wp_real, &
                           mpi_comm_airplane, ierr)
    call mpi_allgather(dipole_all(1,istart), icalc*3, mpi_wp_real, &
                           recv_dipole(1,0), icalc*3, mpi_wp_real, &
                           mpi_comm_airplane, ierr)

    j = mod((ngrid - istart0), nreplica)
    if (j == 0) then
      energy_all(istart0:)   = recv_energy
      grad_all(:,istart0:)   = recv_grad
      dipole_all(:,istart0:) = recv_dipole

    else
      j = nreplica - j

      isiz1 = icalc
      isiz2 = icalc - 1

      ist1 = 0
      ist2 = istart0
      do n = 1, j
        energy_all(ist2:ist2+isiz2-1) = recv_energy(ist1:ist1+isiz1-1)
        grad_all(:,ist2:ist2+isiz2-1) = recv_grad(:,ist1:ist1+isiz1-1)
        dipole_all(:,ist2:ist2+isiz2-1) = recv_dipole(:,ist1:ist1+isiz1-1)
        ist1 = ist1 + isiz1
        ist2 = ist2 + isiz2
      end do 

      isiz2 = icalc
      do n = j + 1, nreplica
        energy_all(ist2:ist2+isiz2-1) = recv_energy(ist1:ist1+isiz1-1)
        grad_all(:,ist2:ist2+isiz2-1) = recv_grad(:,ist1:ist1+isiz1-1)
        dipole_all(:,ist2:ist2+isiz2-1) = recv_dipole(:,ist1:ist1+isiz1-1)
        ist1 = ist1 + isiz1
        ist2 = ist2 + isiz2
      end do 
    end if

    deallocate(recv_energy, recv_grad, recv_dipole)

#endif

    1000 continue

    ! RMS grad at the current geometry
    !
    rmsg = 0.0_wp
    ij   = 0
    do i = 1, nat
      do j = 1,3
         ij = ij + 1
         rmsg = rmsg + grad_all(ij,0)*grad_all(ij,0) * CONV_UNIT_FORCE * CONV_UNIT_FORCE
      end do
    end do
    rmsg = sqrt(rmsg/dble(nat3))

    if (main_rank) then
      write(MsgOut,*)
      write(MsgOut,'(''  RMSD of the gradient at the input geometry = '',d15.6, &
                   & '' [kcal mol-1 Angs-1]'')') rmsg

      if (rmsg > 0.35_wp) then
        write(MsgOut,'(4x,40(''=''))')
        write(MsgOut,'(4x,''!! Warning !!'')')
        write(MsgOut,'(4x,''RMSG is too large for the following vibrational '' ,/, &
                    &  4x,''analysis to be valid. Further miminization until '',/, &
                    &  4x,''RMSG < 0.35 is highly recommended.'')')
        write(MsgOut,'(4x,40(''=''))')
      end if
      write(MsgOut,*)
      write(MsgOut,*)
    end if

    ! Now get all return values
    !
    energy = energy_all(0)
    grad   = grad_all(:,0)
    dipole = dipole_all(:,0)

    allocate(hess_tmp(nat3, nat_calc*3), dd_tmp(3, nat_calc*3))

    do icyc = 1, ngrid-1, 2
      call setaddress(icyc, vibration%minfo_folder, iatom, ixyz, fname)
      delta = vibration%diff_stepsize/sqrt(mass(iatom))

      i = (icyc+1)/2
      hess_tmp(:,i) = (grad_all(:,icyc) - grad_all(:,icyc+1))   &
                      / delta*HALF * CONV_UNIT_LEN
      dd_tmp(:,i)   = (dipole_all(:,icyc) - dipole_all(:,icyc+1)) &
                      / delta*HALF * CONV_UNIT_LEN

    end do

    if (nat == nat_calc) then
      hessian     = hess_tmp
      dipole_derv = dd_tmp

    else

      ! merge the calculataed hessian/dd with the exsiting one
      do i1 = 1, nat
        if (.not. vibatom_calc(i1)) cycle

        do i = 1, nat_calc
          if (vibatom_calc_id(i) == vibration%vibatom_id(i1)) then
            i2 = i
            exit
          end if
        end do
        k1 = (i1-1)*3 + 1
        k2 = (i2-1)*3 + 1

        do j = 1, 3

          hessian(:,k1) = hess_tmp(:,k2)
          dipole_derv(:,k1) = dd_tmp(:,k2)

          k1 = k1 + 1
          k2 = k2 + 1

        end do
      end do

      do i = 1, nat
        if (vibatom_calc(i)) cycle
        k1 = (i-1)*3 + 1

        do j = 1, nat
          if (.not. vibatom_calc(j)) cycle
          k2 = (j-1)*3 + 1

          hessian(k2,k1:k1+2) = hessian(k1:k1+2,k2)
          hessian(k2+1,k1:k1+2) = hessian(k1:k1+2,k2+1)
          hessian(k2+2,k1:k1+2) = hessian(k1:k1+2,k2+2)

        end do
      end do

    end if

    ! Symmetrize the Hessian matrix 
    do k1 = 1, nat3
      do k2 = 1, k1-1
         hessian(k2,k1) = (hessian(k2,k1) + hessian(k1,k2))*HALF
         hessian(k1,k2) = hessian(k2,k1)
      end do
    end do

    deallocate(hess_tmp, dd_tmp)
    deallocate(vibatom_calc_id, vibatom_calc)
    deallocate(energy_all, grad_all, dipole_all)

    !dbg if (main_rank) then
    !dbg   write(MsgOut,'("Hessian")')
    !dbg   do i = 1, nat3
    !dbg     write(MsgOut,*)
    !dbg     write(MsgOut,'(10e15.6)') hessian(:,i)
    !dbg   end do
    !dbg end if

    return

  contains

    subroutine setaddress(icyc, minfo_folder, iatom, ixyz, fname)

      integer :: icyc, iatom, ixyz
      character(MaxFilename) :: minfo_folder, fname

      integer :: m, ivib
      character(MaxFilename) :: basename
      character(10)          :: c_iatom

      if (icyc == 0) then

        ! -------------------------------------------
        ! energy and gradient at the current geometry

        fname = trim(minfo_folder)//'/eq.minfo'

      else 

        ! -------------------------------------------
        ! energy and gradient at +/- delta

        m = mod(icyc, 6)
        if (m == 0) then
          ivib = icyc/6
        else
          ivib = (icyc-m)/6+1
        end if
        iatom = vibatom_calc_id(ivib)
        write(c_iatom,'(i0)') iatom
        basename = trim(minfo_folder)//"/"//trim(c_iatom)//"_"// &
                   trim(molecule%atom_name(iatom))

        if (m == 1 .or. m == 2) then
          ixyz = 1
        else if (m == 3 .or. m == 4) then
          ixyz = 2
        else
          ixyz = 3
        end if

        if (mod(icyc, 2) == 1) then
          ! + delta
          fname = trim(basename)//"+"//xyz(ixyz)//".minfo"
        else
          ! - delta
          fname = trim(basename)//"-"//xyz(ixyz)//".minfo"
        end if

      end if
 
    end subroutine

  end subroutine calc_hessian

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    normal_mode
  !> @brief        Carry out normal mode analysis
  !> @authors      KY
  !! @param[in]    molecule  : molecule information
  !! @param[in]    dynvars   : dynamic variables information
  !! @param[in]    vibration : vibration information
  !! @param[in]    nat         : number of atoms
  !! @param[in]    energy      : total energy at the current geometry
  !! @param[in]    grad(nat3)  : gradient
  !! @param[in]    hessian(nat3, nat3) : Hessian matrix
  !! @param[in]    dipole(3)   : dipole moment
  !! @param[in]    dipole_derv(3,nat3) : dipole moment
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine normal_mode(molecule, dynvars, vibration, &
                         nat, energy, grad, hessian, dipole, dipole_derv)

    ! formal arguments
    type(s_molecule),  intent(inout) :: molecule
    type(s_dynvars),   intent(inout) :: dynvars
    type(s_vibration), intent(inout) :: vibration

    integer  :: nat
    real(wp) :: energy
    real(wp) :: grad(nat*3)
    real(wp) :: hessian(nat*3,nat*3)
    real(wp) :: dipole(3)
    real(wp) :: dipole_derv(3,nat*3)

    ! local variables
    integer  :: nat3
    integer  :: i, j, ij, k, i1, j1, k1, k2
    integer  :: ncolomn, nc, iatom
    real(wp), allocatable :: hess_pack(:), hess_mw_pack(:)
    real(wp), allocatable :: dd_mw(:,:)
    real(wp), allocatable :: sq_mass3(:)
    real(wp), allocatable :: vec(:,:), omega(:)
    real(wp), allocatable :: dd2(:,:), infrared(:)

    character(1)           :: xyz(3)


    nat3 = nat*3

    ! for print 
    xyz(1) = 'X'
    xyz(2) = 'Y'
    xyz(3) = 'Z'

    ! Pack the Hessian matrix 
    allocate(hess_pack(nat3*(nat3+1)/2))
    k = 1
    do k1 = 1, nat3
      do k2 = 1, k1
         hess_pack(k) = hessian(k2,k1)
         k = k + 1
      end do
    end do

    ! Mass-weight the Hessian matrix and dipole derivatives
    allocate(hess_mw_pack(nat3*(nat3+1)/2), dd_mw(3,nat3), sq_mass3(nat3))
    ij = 1
    do i = 1, nat
      do j = 1, 3
         sq_mass3(ij) = sqrt(molecule%mass(vibration%vibatom_id(i))*ELMASS)
         ij = ij + 1
      end do
    end do

    k = 1
    do k1 = 1, nat3
      do k2 = 1, k1
        hess_mw_pack(k) = hess_pack(k) / sq_mass3(k1) / sq_mass3(k2)
        k = k + 1
      end do
    end do

    do k = 1, nat3
      dd_mw(:,k) = dipole_derv(:,k) / sq_mass3(k)
    end do

    allocate(vec(nat3,nat3), omega(nat3))

    call diagonalize(nat3, nat3, hess_mw_pack, vec, omega)

    do k = 1, nat3
      if (omega(k) > 0.0_wp) then
        omega(k) = sqrt(abs(omega(k))) * HARTREE_WAVENUM
      else
        omega(k) =-sqrt(abs(omega(k))) * HARTREE_WAVENUM
      end if
    end do

    ! Transform the dipole_derivatives in terms of normal coordinates
    allocate(dd2(3,nat3), infrared(nat3))
    do i = 1, nat3
      dd2(:,i) = 0.0_wp
      do j = 1, nat3
        dd2(:,i) = dd2(:,i) + dd_mw(:,j)*vec(j,i)
      end do
    end do

    do i = 1, nat3
      infrared(i) = 0.0_wp
      do j = 1, 3
         infrared(i) = infrared(i) + dd2(j,i)*dd2(j,i)
      end do
      infrared(i) = infrared(i) * PI  / 3.0_wp / VLIGHT_IN_AU / VLIGHT_IN_AU &
                  * CONV_UNIT_LEN * 1.0e-13_wp * AVOGADRO
    end do

    ! print normal modes
    ncolomn = 5
    nc = (nat3-mod(nat3,ncolomn))/ncolomn

    if (main_rank) then
      write(MsgOut,'(''  Harmonic frequencies and normal displacement vectors'')')
      k=1
      do i = 1, nc
        write(MsgOut,'(6x,''    mode  '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(i12,$)') j
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''    omega    '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(f12.4,$)') omega(j)
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''   IR int.   '',$)')
        do j = k, k + ncolomn - 1
          write(MsgOut,'(f12.4,$)') infrared(j)
        end do
        write(MsgOut,*)

        k1 = 1
        do i1 = 1, nat
          iatom = vibration%vibatom_id(i1)
          do j1 = 1, 3
            write(MsgOut,'(6x,i4,x,a6,x,a1,$)') i1,molecule%atom_name(iatom),xyz(j1)
            do j = k, k + ncolomn - 1
              write(MsgOut,'(f12.4,$)') vec(k1,j)
            end do
            write(MsgOut,*)
            k1 = k1 + 1
          end do
        end do

        write(MsgOut,*)
        k = k + ncolomn
      end do

      if (mod(nat3,ncolomn) /= 0) then
        write(MsgOut,'(6x,''    mode  '',$)')
        do j = k, nat3
          write(MsgOut,'(i12,$)') j
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''    omega    '',$)')
        do j = k, nat3
          write(MsgOut,'(f12.4,$)') omega(j)
        end do
        write(MsgOut,*)

        write(MsgOut,'(6x,''   IR int.   '',$)')
        do j = k, nat3
          write(MsgOut,'(f12.4,$)') infrared(j)
        end do
        write(MsgOut,*)

        k1 = 1
        do i1 = 1, nat
          iatom = vibration%vibatom_id(i1)
          do j1 = 1, 3
            write(MsgOut,'(6x,i4,x,a6,x,a1,$)') i1,molecule%atom_name(iatom),xyz(j1)
            do j = k, nat3
              write(MsgOut,'(f12.4,$)') vec(k1,j)
            end do
            write(MsgOut,*)
            k1 = k1 + 1
          end do
        end do

        write(MsgOut,*)
      end if

      call print_minfo(vibration, molecule, dynvars, nat, energy, grad, &
                       hess_pack, dipole, dipole_derv, omega, vec)


    end if

    deallocate(hess_pack, hess_mw_pack, sq_mass3)
    deallocate(vec, omega)
    deallocate(dd_mw, dd2, infrared)

    return

  end subroutine normal_mode

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    match_vibatom
  !> @brief        check the match between old and new vibatom
  !! @authors      KY
  !! @param[in]    minfo_in            : data of existing minfo from the input
  !! @param[inout] nat                 : number of atoms
  !! @param[inout] vibatom_id(nat)     : ID of vibatom
  !! @param[inout] hessian(nat3,nat3)  : Hessian matrix
  !! @param[inout] dipole_derv(3,nat3) : Dipole derivatives
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine match_vibatom(minfo_in , nat, vibatom_id, hessian, dipole_derv)

    ! formal arguments
    type(s_minfo), target, intent(in) :: minfo_in
    integer,       intent(inout) :: nat, vibatom_id(nat)
    real(wp),      intent(inout) :: hessian(nat*3, nat*3)
    real(wp),      intent(inout) :: dipole_derv(3, nat*3)

    integer,  pointer  :: vibatom_in_id(:)
    integer :: nat_in, tmp(nat)
    integer :: i, j, k, ncalc
    logical :: found
    real(wp), allocatable :: hess_in(:,:)


    vibatom_in_id => minfo_in%vibatom_id
    nat_in = minfo_in%nat

    allocate(hess_in(nat_in*3, nat_in*3))
    ! unpack Hessian
    k = 0
    do i = 1, nat_in*3
    do j = 1, i
       k = k + 1
       hess_in(j,i) = minfo_in%hessian(k)
       hess_in(i,j) = hess_in(j,i)
    end do
    end do

    tmp = vibatom_id

    ncalc = 0
    do i = 1, nat

      found = .false.
      do j = 1, nat_in
        if (tmp(i) == vibatom_in_id(j)) then
          found = .true.
          hessian((i-1)*3+1:i*3,:) = hess_in((j-1)*3+1:j*3,:)
          hessian(:,(i-1)*3+1:i*3) = hess_in(:,(j-1)*3+1:j*3)
          dipole_derv(:,(i-1)*3+1:i*3) = minfo_in%dipole_derv(:,(j-1)*3+1:j*3)
          exit
        end if 
      end do

      if (.not. found) then
        ncalc = ncalc + 1
        vibatom_id(ncalc) = tmp(i)
      end if

    end do

    nat = ncalc

    return

  end subroutine match_vibatom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    diagonalize
  !> @brief        diagonalize a symmetric matrix using Lapack (dspevx)
  !! @authors      KY
  !! @param[in]    n : dimension of the matrix
  !! @param[in]    m : number of solution
  !! @param[in]    H(n*(n+1)/2) : The matrix to be diagonalized (upper half)
  !! @param[out]   C(n,m)       : The eigenvectors
  !! @param[out]   E(n)         : The eigenvalues
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine diagonalize(n, m, H, C, E)

    integer :: n,m,m1
    real(wp), dimension(n*(n+1)/2) :: H
    real(wp), dimension(n,m) :: C
    real(wp), dimension(n) :: E

    character :: jobz,range,uplo
    real(wp)  :: vl,vu,abstol
    integer   :: il,iu,ldz,info
    integer,  dimension(:), allocatable :: ifail,iwork
    real(wp), dimension(:), allocatable :: work


    !write(6,*) n,m
    !do i = 1, n
    !   write(6,'(11f12.6)') (H(j),j=i*(i-1)/2+1,i*(i+1)/2)
    !dnd do

    allocate(work(10*n),ifail(n),iwork(10*n))

    jobz='V'
    uplo='U'
    vl=0.0_wp
    vu=0.0_wp
    il=0
    iu=0
    if (n==m) then
      range='A'
    else
      range='I'; il=1; iu=m
    end if

    abstol=0.0_wp
    ldz=n

#ifdef LAPACK
    m1=0
    ifail=0; info=0
    call dspevx(jobz,range,uplo,n,H,vl,vu,il,iu,abstol,m1,E, &
                C,ldz,work,iwork,ifail,info)

#else
    call error_msg('diagonalize> ERROR: This subroutine needs LAPACK.')
#endif

    deallocate(work,ifail,iwork)

    if (info==0) return

    if (info<0) then
      write(MsgOut,'(''ERROR IN '',i3,''TH PARAMETER'')') info
    else
      write(MsgOut,'(3x,i3,''EIGENVECTORS FAILED TO CONVERGE'')') info
    end if

    return

  end subroutine diagonalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_minfo
  !> @brief        print the information to minfo format
  !! @authors      KY
  !! @param[in] vibration : vibration information
  !! @param[in] molecule : molecule information
  !! @param[in] dynvars  : dynamic variables information
  !! @param[in] nat      : number of atoms for vibrational analysis
  !! @param[in] energy   : energy
  !! @param[in] grad     : gradient(3*nat)
  !! @param[in] grad     : hessian matrix(3*nat, 3*nat)
  !! @param[in] dipole   : dipole moment (3)
  !! @param[in] dipole_derivative   : dipole moment derivatives (3, nat*3)
  !! @param[in] omega    : harmonic frequency (nat*3)
  !! @param[in] vec      : normal mode displacement vector (nat*3, nat*3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_minfo(vibration, molecule, dynvars, nat, energy, grad, &
                     hess, dipole, dipole_derivative, omega, vec)

    ! formal arguments
    type(s_vibration), target, intent(in)    :: vibration
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_dynvars),   target, intent(in)    :: dynvars

    integer  :: nat
    real(wp) :: energy
    real(wp) :: grad(3*nat), hess(nat*3*(nat*3+1)/2)
    real(wp) :: dipole(3), dipole_derivative(3,nat*3)
    real(wp) :: omega(nat*3), vec(nat*3, nat*3)

    ! local
    type(s_minfo) :: minfo_out


    if (vibration%minfo_natoms > 0) then 
      call init_minfo_atom(vibration%vib_natoms, minfo_out, vibration%minfo_natoms)
      minfo_out%vibatom_id  = vibration%vibatom_id
      minfo_out%subatom_id  = vibration%minfoatom_id
    else
      call init_minfo_atom(vibration%vib_natoms, minfo_out)
      minfo_out%vibatom_id  = vibration%vibatom_id
    end if

    call init_minfo_elec(vibration%vib_natoms, minfo_out)
    minfo_out%energy      = energy
    minfo_out%gradient    = grad
    minfo_out%hessian     = hess
    minfo_out%dipole      = dipole
    minfo_out%dipole_derv = dipole_derivative

    call init_minfo_vib(vibration%vib_natoms*3, minfo_out)
    minfo_out%omega       = omega
    minfo_out%vec         = vec

    call output_minfo(vibration%minfofile, minfo_out, molecule, dynvars%coord)

    return

  end subroutine print_minfo

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_minfo_grad
  !> @brief        print the information to minfo format
  !! @authors      KY
  !! @param[in] filename : filename of the minfo file
  !! @param[in] nat      : number of atoms for vibrational analysis
  !! @param[in] energy   : energy
  !! @param[in] grad     : gradient(3*nat)
  !! @param[in] dipole   : dipole moment (3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine print_minfo_grad(filename, nat, energy, grad, dipole)

    ! formal arguments
    character(MaxFilename) :: filename
    integer                :: nat
    real(wp)               :: energy, grad(nat*3), dipole(3)

    ! local
    type(s_minfo) :: minfo_out


    call init_minfo_elec(nat, minfo_out)
    minfo_out%energy   = energy
    minfo_out%gradient = grad
    minfo_out%dipole   = dipole
    call output_minfo(trim(filename), minfo_out)

    return

  end subroutine print_minfo_grad

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_minfo_grad
  !> @brief        read the information to minfo format
  !! @authors      KY
  !! @param[in]    filename : filename of the minfo file
  !! @param[in]    nat      : number of atoms for vibrational analysis
  !! @param[inout] energy   : energy
  !! @param[inout] grad     : gradient(3*nat)
  !! @param[inout] dipole   : dipole moment (3)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_minfo_grad(filename, nat, energy, grad, dipole)

    ! formal arguments
    character(MaxFilename) :: filename
    integer                :: nat
    real(wp)               :: energy, grad(nat*3), dipole(3)

    ! local
    type(s_minfo)  :: minfo_in
    integer        :: nat_read

    call read_minfo(trim(filename), minfo_in)
    nat_read = size(minfo_in%gradient)/3
    if (nat_read /= nat .and. main_rank) then
      write(MsgOut,'(3x,"Error while reading ",a)') trim(filename)
      write(MsgOut,'(3x,"The number of vibatoms does not match")')
      write(MsgOut,'(3x,"   Nat(input) = ",i8)') nat
      write(MsgOut,'(3x,"   Nat(read ) = ",i8)') nat_read
      write(MsgOut,*)
      call error_msg()
    end if

    energy = minfo_in%energy
    grad   = minfo_in%gradient
    dipole = minfo_in%dipole
 
    return

  end subroutine read_minfo_grad

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_gridpoints
  !> @brief        calculate the energy and gradient at the grid points
  !> @authors      KY
  !! @param[inout] molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] vibration : vibration information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_gridpoints(molecule, enefunc, dynvars, vibration, &
                     output, pairlist, boundary)

    ! formal arguments
    type(s_molecule),  target, intent(inout) :: molecule
    type(s_enefunc),   target, intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_vibration), target, intent(inout) :: vibration
    type(s_output),            intent(inout) :: output
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary

    ! local
    type(s_qmmm), pointer :: qmmm
    integer               :: nat, nat3
    integer               :: i, ia, j, ierr
    real(wp), pointer     :: coord(:,:)
    integer , pointer     :: vibatom_id(:)

    character(MaxFilename) :: folder, minfo_folder, datafile, basename, fname
    character(MaxLine)     :: line
    character(4)   :: fid, fnum, label
    character(7)   :: num
    integer        :: icount, count, replicaid

    integer        :: ifile, idatafile
    logical        :: ex, thisrun_done, restart
    real(wp)       :: energy, grad(3,vibration%vib_natoms), dipole(3)


    nat        =  vibration%vib_natoms
    nat3       =  nat*3
    vibatom_id => vibration%vibatom_id
    coord      => dynvars%coord
    qmmm       => enefunc%qmmm

    ! replica id
    replicaid = my_country_no + 1
    count = 0

    if (qmmm%do_qmmm) then
      icount        = 0
      qmmm%ene_only = vibration%grid_ene_only
    end if

    if (vibration%grid_ene_only) then
      ! set datafile
      write(num,'(i0)') my_country_no
      datafile = trim(vibration%datafile)//'_'//trim(num)

      ! ignore qm error
      qmmm%ignore_qm_error = .true.

      if (main_rank) &
        write(MsgOut,*) 'Compute energy at grid points: &
         &data written to [ '//trim(vibration%datafile)//' ]'

    else
      ! create a folder to save minfo files
      minfo_folder = vibration%minfo_folder
      if (main_rank) then
        call system('mkdir -p '//trim(minfo_folder)//' > /dev/null 2>&1')
      end if

      if (main_rank) &
        write(MsgOut,*) 'Compute energy at grid points: &
         &minfo files created in [ '//trim(minfo_folder)//' ]'

    end if

    call open_file(ifile, trim(vibration%gridfile), IOFileInput)

    restart = .true.
    do while(.true.)
      read(ifile,*,end=20)
      read(ifile,'(a)') basename

      do i = 1, nat
        ia = vibatom_id(i)
        read(ifile,*) label,coord(:,ia)
      end do

      count = count + 1

      if (mod(count,vibration%nreplica) == replicaid-1) then

        ! check for restart
        if (restart) then
          thisrun_done = .false.
          if (vibration%grid_ene_only) then
            fname = trim(datafile)
            inquire(file=trim(datafile),exist=ex)
            if (ex) then
              call open_file(idatafile, trim(fname), IOFileInput)
              do while(.true.)
                read(idatafile,'(a)',end=10) line
                i = index(line,',')
                if (trim(basename) == line(1:i-1)) then
                   thisrun_done = .true.
                   exit
                end if 
              end do
           10 continue
              !close(idatafile)
              call close_file(idatafile)
            end if
    
          else
            fname = trim(minfo_folder)//'/'//trim(basename)
            inquire(file=trim(fname)//'.minfo',exist=thisrun_done)
    
          end if
          if (thisrun_done) then
            cycle
          else
            restart = .false.
          end if

        end if

        call compute_energy(molecule, enefunc, pairlist, boundary,  &
                          .true., &
                          .false.,               &
                          dynvars%coord,         &
                          dynvars%trans,         &
                          dynvars%coord_pbc,     &
                          dynvars%energy,        &
                          dynvars%temporary,     &
                          dynvars%force,         &
                          dynvars%force_omp,     &
                          dynvars%virial,        &
                          dynvars%virial_extern)


        energy =  dynvars%energy%total / CONV_UNIT_ENE
        do i = 1, nat
          do j = 1, 3
            grad(j,i) = -dynvars%force(j,vibatom_id(i)) / CONV_UNIT_FORCE
          end do
        end do
        dipole = qmmm%qm_dipole

        if (replica_main_rank) then
          write(MsgOut,'(6x,"Done for ",a30," :",4x,"replicaID = ",i5)') &
                trim(basename), replicaid

          if (vibration%grid_ene_only) then
            if (qmmm%qmmm_error) then
              energy = 0.0_wp
              dipole = 0.0_wp
            end if

            fname = trim(datafile)
            call open_file(idatafile, trim(fname), IOFileOutputAppend)
            write(idatafile,'(a,'', '',f25.13,'', '',2(es16.9,'',''),es16.9)') &
                  trim(basename), energy, dipole
            call close_file(idatafile)
            
          else
            fname = trim(minfo_folder)//'/'//trim(basename)//'.minfo'
            call print_minfo_grad(fname, nat, energy, grad, dipole)

          end if
        end if

      end if

    end do

 20 continue
    call close_file(ifile)

#ifdef HAVE_MPI_GENESIS
    call mpi_barrier(mpi_comm_world, ierr)
#endif

    ! combine the grid data to one file
    if (vibration%grid_ene_only .and. main_rank) then
      
      call open_file(ifile, trim(vibration%datafile), IOFileOutputAppend)

      do i = 1, vibration%nreplica

        write(num,'(i0)') i-1
        datafile = trim(vibration%datafile)//'_'//trim(num)
        call open_file(idatafile, trim(datafile), IOFileInput)
        do while(.true.)
          read(idatafile,'(a)',end=30) line
          j = index(line,',')
          read(line(j+1:),*) energy
          if (abs(energy) > 1.e-10)  write(ifile,'(a)') trim(line)
        end do
     30 continue

        call close_file(idatafile,'delete')

      end do

      call close_file(ifile)

    end if

    return

  end subroutine calc_gridpoints

end module at_vibration_mod

