!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_experiments_mod
!> @brief   Provides experimental data.
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
  
module cg_experiments_mod

  use molecules_str_mod
  use cg_restraints_str_mod
  use cg_experiments_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use cg_boundary_str_mod
  use fileio_mod
  use fileio_sit_mod
  use fileio_mrc_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use string_mod
  use timers_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
    use omp_lib, only: omp_get_thread_num
#endif

  implicit none
  private

  ! structures
  type, public :: s_exp_info
    logical                         :: emfit           = .false.
    character(MaxFilename)          :: emfit_target    = ''
    real(wp)                        :: emfit_sigma     = 2.5_wp
    real(wp)                        :: emfit_tolerance = 0.001_wp
    integer                         :: emfit_period    = 1
    real(wp)                        :: emfit_zero_threshold = 0.0_wp
  end type s_exp_info

  type(s_experiments), target, save :: experiments
  integer,                     save :: natom_local, emfit_icycle
  real(wp),                    save :: corrcoeff_save
  integer,        allocatable, save :: list_local(:)
  integer,        allocatable, save :: icount(:), domain_index(:)
  integer,        allocatable, save :: ig_min_local(:,:), ig_max_local(:,:)
  integer,        allocatable, save :: num_atoms_local(:)
  integer,        allocatable, save :: ig_lower(:,:), ig_upper(:,:)
  real(wp),       allocatable, save :: dot_exp_drhodx(:), dot_exp_drhody(:), dot_exp_drhodz(:)
  real(wp),       allocatable, save :: erfa(:,:), expa(:,:)
  real(wp),       allocatable, save :: derfa_x(:,:), derfa_y(:,:), derfa_z(:,:)
  real(wp),       allocatable, save :: dexpa_x(:,:), dexpa_y(:,:), dexpa_z(:,:)
  real(wp),       allocatable, save :: dot_sim_drhodx(:), dot_sim_drhody(:), dot_sim_drhodz(:)
  real(wp),       allocatable, save :: emfit_force(:,:)

  ! subroutines
  public  :: show_ctrl_experiments
  public  :: read_ctrl_experiments
  public  :: setup_experiments
  private :: setup_experiments_emfit
  public  :: compute_energy_experimental_restraint
  private :: compute_energy_experimental_restraint_emfit

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_experiments
  !> @brief        Displays EXPERIMENTS section usage.
  !! @authors      TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_experiments(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[EXPERIMENTS]'
        write(MsgOut,'(A)') 'emfit           = NO          # EM fit'
        write(MsgOut,'(A)') '# emfit_target    = sample.sit  # EM data file'
        write(MsgOut,'(A)') '# emfit_sigma     = 2.5         # resolution parameter of the simulated map'
        write(MsgOut,'(A)') '# emfit_tolerance = 0.001       # Tolerance for error'
        write(MsgOut,'(A)') '# emfit_period    = 1           # emfit force update period'
        write(MsgOut,'(A)') ''

      end select

    end if

    return

  end subroutine show_ctrl_experiments
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_experiments
  !> @brief        Reads EXPERIMENTS section from the control file.
  !! @authors      TM
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   exp_info : EXPERIMENTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_experiments(handle, exp_info) 

    ! parameters
    character(*),            parameter     :: Section = 'Experiments'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_exp_info),        intent(inout) :: exp_info 

    ! local variables


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'emfit',           &
                               exp_info%emfit)

    call read_ctrlfile_string (handle, Section, 'emfit_target',    &
                              exp_info%emfit_target)

    call read_ctrlfile_real   (handle, Section, 'emfit_sigma',     &
                              exp_info%emfit_sigma)

    call read_ctrlfile_real   (handle, Section, 'emfit_tolerance', &
                              exp_info%emfit_tolerance)

    call read_ctrlfile_integer(handle, Section, 'emfit_period',    &
                              exp_info%emfit_period)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      if (exp_info%emfit) then

        write(MsgOut,'(a)') 'Read_Ctrl_Experiments: Parameters for experimental data fitting'
        write(MsgOut,'(a20,a)') '  emfit_target    = ', trim(exp_info%emfit_target)
        write(MsgOut,'(a20,F10.4,a20,F10.4)')                      &
              '  emfit_sigma     = ', exp_info%emfit_sigma,        &
              '  emfit_tolerance = ', exp_info%emfit_tolerance
        write(MsgOut,'(a20,I10)')                                  &
              '  emfit_period    = ', exp_info%emfit_period
        write(MsgOut,'(a)') ''

      end if

    end if

    return

  end subroutine read_ctrl_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments
  !> @brief        Sets up experiment information.
  !! @authors      TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameters
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments(exp_info, molecule, restraints, enefunc)

    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i
    logical    :: do_emfit_res

    experiments%do_emfit = .false.

    if (exp_info%emfit) then

      do_emfit_res = .false.
      do i = 1, restraints%nfunctions
        if (enefunc%restraint_kind(i) == RestraintsFuncEM) do_emfit_res = .true.
      end do

      if (.not. do_emfit_res) then
        call error_msg('Setup_Experiments: EM is not defined in [RESTRAINTS].')
      end if

      experiments%do_emfit = .true.
      call setup_experiments_emfit(exp_info, molecule)

    end if

    return

  end subroutine setup_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments_emfit
  !> @brief        Sets up experiment information.
  !! @authors      TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameter
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments_emfit(exp_info, molecule)

    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule

    ! local variables
    type(s_sit)             :: sit
    type(s_mrc)             :: mrc
    integer                 :: i, j, k, n(3), max_nxyz
    real(wp)                :: pi, r, y
    real(wp)                :: x0, y0, z0, dx, dy, dz
    integer                 :: ix, iy, iz, idx
    character(10)           :: file_ext
    logical                 :: format_mrc, format_sit

    ! read mapfile
    !
    format_mrc = .false.
    format_sit = .false.
    file_ext   = get_extension(exp_info%emfit_target)

    if (file_ext == "map" .or. file_ext == "ccp4" .or. file_ext == "mrc") then
      format_mrc = .true.
      call input_mrc(exp_info%emfit_target, mrc)
      n(1) = mrc%nx
      n(2) = mrc%ny
      n(3) = mrc%nz
      x0   = mrc%origin(1)
      y0   = mrc%origin(2)
      z0   = mrc%origin(3)
      dx   = mrc%cella(1)/mrc%nx
      dy   = mrc%cella(2)/mrc%ny
      dz   = mrc%cella(3)/mrc%nz
    else if (file_ext == "sit") then
      format_sit = .true.
      call input_sit(exp_info%emfit_target, sit)
      n(1) = sit%nx
      n(2) = sit%ny
      n(3) = sit%nz
      x0   = sit%x0
      y0   = sit%y0
      z0   = sit%z0
      dx   = sit%dx
      dy   = sit%dx
      dz   = sit%dx
    else
      call error_msg('Setup_Experiments_Emfit: Unrecognized file format '//&
	             'for the EM density map.')
    end if

    ! allocate array
    !
    max_nxyz = maxval(n)
    call alloc_experiments(experiments, ExperimentsEmfit, &
                           n(1), n(2), n(3),              &
                           molecule%num_atoms, max_nxyz)

    experiments%emfit%nx = n(1)
    experiments%emfit%ny = n(2)
    experiments%emfit%nz = n(3)
    experiments%emfit%x0 = x0
    experiments%emfit%y0 = y0
    experiments%emfit%z0 = z0
    experiments%emfit%dx = dx
    experiments%emfit%dy = dy
    experiments%emfit%dz = dz

    experiments%emfit%sigma        = exp_info%emfit_sigma
    experiments%emfit%tolerance    = exp_info%emfit_tolerance
    experiments%emfit%emfit_period = exp_info%emfit_period

    ! set target EM density map
    !
    do iz = 0, n(3) - 1
      do iy = 0, n(2) - 1
        do ix = 0, n(1) - 1
          idx = 1 + ix + iy*n(1) + iz*n(1)*n(2)

          if (format_mrc) then

            if (mrc%map_value(idx) < exp_info%emfit_zero_threshold) then
              experiments%emfit%target_map(ix,iy,iz) = 0.0_wp
            else
              experiments%emfit%target_map(ix,iy,iz) = mrc%map_value(idx)
            end if

          else if (format_sit) then

            if (sit%map_value(idx) < exp_info%emfit_zero_threshold) then
              experiments%emfit%target_map(ix,iy,iz) = 0.0_wp
            else
              experiments%emfit%target_map(ix,iy,iz) = sit%map_value(idx)
            end if

          end if

        end do
      end do
    end do

    ! need to calculate norm for force calculation
    !
    experiments%emfit%norm_exp = 0.0_wp
    do i = 0, n(1) - 1
      do j = 0, n(2) - 1
        do k = 0, n(3) - 1
          experiments%emfit%norm_exp = experiments%emfit%norm_exp &
                                     + experiments%emfit%target_map(i,j,k)**2
        end do
      end do
    end do
    experiments%emfit%norm_exp = sqrt(experiments%emfit%norm_exp)

    ! boundary of voxel i is [bound_x(i) bound_x(i+1)]
    ! voxel from (0,0,0) to (nx-1, ny-1, nz-1)
    !
    do i = 0, n(1)
      experiments%emfit%bound_x(i) = x0 + dx * (dble(i) - 0.5_wp)
    end do
    do i = 0, n(2)
      experiments%emfit%bound_y(i) = y0 + dy * (dble(i) - 0.5_wp)
    end do
    do i = 0, n(3)
      experiments%emfit%bound_z(i) = z0 + dz * (dble(i) - 0.5_wp)
    end do

    ! determine_cutoff
    !
    pi = acos(-1.0_wp)
    r  = 0.0_wp
    y  = 0.0_wp

    do while (1.0_wp - y > exp_info%emfit_tolerance)
      y = erf(sqrt(3.0_wp/2.0_wp)*r) - sqrt(6.0_wp/pi)*r*exp(-3.0_wp*r*r/2.0_wp)
      r = r + 0.01_wp
    end do

    experiments%emfit%n_grid_cut_x = ceiling(r*exp_info%emfit_sigma/dx)
    experiments%emfit%n_grid_cut_y = ceiling(r*exp_info%emfit_sigma/dy)
    experiments%emfit%n_grid_cut_z = ceiling(r*exp_info%emfit_sigma/dz)

    call dealloc_mrc(mrc)
    call dealloc_sit(sit)

    ! for MPI parallelization
    !
    allocate(icount(nproc_city))
    allocate(ig_min_local(nproc_city,3))
    allocate(ig_max_local(nproc_city,3))
    allocate(domain_index(molecule%num_atoms))
    allocate(num_atoms_local(nproc_city))

    ! force update
    !
    allocate(list_local (  molecule%num_atoms))
    allocate(emfit_force(3,molecule%num_atoms))
    emfit_icycle = -1

    ! write summary
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Experiments_Emfit: Sets up variables for EMFIT.'
      write(MsgOut,'(A20,F10.3)') '  radius/sigma    = ', r
      write(MsgOut,'(A20,F10.3)') '  radius          = ', r*experiments%emfit%sigma
      write(MsgOut,'(A20,F10.3)') '  dx              = ', experiments%emfit%dx
      write(MsgOut,'(A20,F10.3)') '  dy              = ', experiments%emfit%dy
      write(MsgOut,'(A20,F10.3)') '  dz              = ', experiments%emfit%dz
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along x = ', &
                              experiments%emfit%n_grid_cut_x
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along y = ', &
                              experiments%emfit%n_grid_cut_y
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along z = ', &
                              experiments%emfit%n_grid_cut_z
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_experiments_emfit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint
  !> @brief        Calculates restraint energy from experimental data.
  !! @authors      TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint(enefunc, coord, inum, &
                                           calc_force, force, virial, eexp, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eexp
    real(wp),                intent(inout) :: cv

    if (experiments%do_emfit) then
      call compute_energy_experimental_restraint_emfit &
             (enefunc, coord, inum, calc_force, force, virial, eexp, cv)
    end if

    return

  end subroutine compute_energy_experimental_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint_emfit
  !> @brief        Calculates restraint energy from experimental data.
  !! @authors      TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint_emfit(inum, calc_force, &
                           domain, enefunc, force, eemfit, emcorr)

    ! formal arguments
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: eemfit
    real(dp),                intent(inout) :: emcorr

    ! local variables
    integer           :: i, j, k, n, m, natom_all, ix
    integer           :: i1, j1, k1, ncell
    integer           :: nlen, ixx, icycle, ncycle
    integer           :: idomain
    integer           :: id
    integer           :: n_grid_cut
    real(wp)          :: yzfactor, zfactor
    real(wp)          :: drho_dx(3), before_allreduce(2), after_allreduce(2)
    real(wp)          :: dcc_dx, dcc_dy, dcc_dz
    real(wp)          :: norm_sim, dot_exp_sim
    real(wp)          :: f, g, d, coeff_rho, coeff_drho
    real(wp)          :: pi, vol, corrcoeff
    real(wp)          :: exp_tmp(3), sim_tmp(3)
    real(wp)          :: weight, inv_dx
    real(wp)          :: fact1, fact2, fact3
    integer           :: i00, i01
    integer           :: natom_domain
    integer           :: irequest6, irequest7
    real(wp)          :: map_ori(3)
#ifdef HAVE_MPI_GENESIS
    integer           :: istatus(mpi_status_size)
#endif
    integer,  pointer :: id_g2l(:)
    integer,  pointer :: natom_cell(:)
    real(wp), pointer :: norm_exp, sigma
    real(wp), pointer :: simulated_map(:,:,:), target_map(:,:,:)
    real(wp), pointer :: bound_x(:), bound_y(:), bound_z(:)

    ! use pointer
    !
    ncell          =  domain%num_cell_local + domain%num_cell_boundary
    natom_cell     => domain%num_atom
    id_g2l         => domain%id_g2l
    natom_all      =  domain%num_atom_all
    sigma          => experiments%emfit%sigma
    map_ori(1)     =  experiments%emfit%x0
    map_ori(2)     =  experiments%emfit%y0
    map_ori(3)     =  experiments%emfit%z0
    bound_x        => experiments%emfit%bound_x
    bound_y        => experiments%emfit%bound_y
    bound_z        => experiments%emfit%bound_z
    norm_exp       => experiments%emfit%norm_exp
    simulated_map  => experiments%emfit%simulated_map
    target_map     => experiments%emfit%target_map
    n_grid_cut     =  experiments%emfit%n_grid_cut_x
    inv_dx         =  1.0_wp/experiments%emfit%dx
    weight         =  enefunc%restraint_const(1,inum)


    ! perform emfit or not
    !
    if (skip_emfit) then
      do i = 1, domain%num_atom_domain
        n = domain%id_l2g(i)
        force(i,1,1) = force(i,1,1) + emfit_force(1,n)
        force(i,2,1) = force(i,2,1) + emfit_force(2,n)
        force(i,3,1) = force(i,3,1) + emfit_force(3,n)
      end do
      corrcoeff = corrcoeff_save
      experiments%emfit%corrcoeff = corrcoeff
      emcorr = corrcoeff
      eemfit = weight * (1.0_wp - corrcoeff)
      return
    end if

    call timer(TimerEmfit, TimerOn)

    ! coefficient for drho/d[xyz] and rho
    !   - sign is missing in drho/dq equstions of BJ 2012
    !
    pi    = acos(-1.0_wp)
    vol   = experiments%emfit%dx * experiments%emfit%dy * experiments%emfit%dz
    coeff_drho = - sigma**2 * pi / 6.0_wp / vol
    coeff_rho  =  (sigma**2 * pi / 6.0_wp)**(3.0_wp/2.0_wp) / vol


    DO idomain = 1, 0, -1

      ! make calculation list in domain00 and domain10
      !
      if (idomain == 0) then

        natom_domain = natom_domain00

        !$omp parallel do private(i, j)
        !
        do i = 1, natom_domain
          ig_domain(1,i) = ig_domain_tmp(1,i)
          ig_domain(2,i) = ig_domain_tmp(2,i)
          ig_domain(3,i) = ig_domain_tmp(3,i)

          emfit_atom_local(i) = .false.
          if (ig_domain(1,i) < ig_min00(1) - n_grid_cut) cycle
          if (ig_domain(1,i) > ig_max00(1) + n_grid_cut) cycle
          if (ig_domain(2,i) < ig_min00(2) - n_grid_cut) cycle
          if (ig_domain(2,i) > ig_max00(2) + n_grid_cut) cycle
          if (ig_domain(3,i) < ig_min00(3) - n_grid_cut) cycle
          if (ig_domain(3,i) > ig_max00(3) + n_grid_cut) cycle
          emfit_atom_local(i) = .true.

          coord_tmp(1:3,i) = coord_domain0(1:3,i)
          do j = 1, 3
            ig_lower(j,i) = ig_domain(j,i) - n_grid_cut
            ig_upper(j,i) = ig_domain(j,i) + n_grid_cut
            if (ig_max00(j) < ig_upper(j,i)) ig_upper(j,i) = ig_max00(j)
            if (ig_min00(j) > ig_lower(j,i)) ig_lower(j,i) = ig_min00(j)
          end do

        end do
        !$omp end parallel do

        simulated_map(ig_min00(1):ig_max00(1), &
                      ig_min00(2):ig_max00(2), &
                      ig_min00(3):ig_max00(3)) = 0.0_wp

      else if (idomain == 1) then

        natom_domain = natom_domain10

        !$omp parallel do private(i, j)
        !
        do i = 1, natom_domain
          ig_domain(1,i) = int((coord_domain1(1,i) - map_ori(1))*inv_dx + 0.5_wp)
          ig_domain(2,i) = int((coord_domain1(2,i) - map_ori(2))*inv_dx + 0.5_wp)
          ig_domain(3,i) = int((coord_domain1(3,i) - map_ori(3))*inv_dx + 0.5_wp)

          emfit_atom_local(i) = .false.
          if (ig_domain(1,i) < ig_min10(1) - n_grid_cut) cycle
          if (ig_domain(1,i) > ig_max10(1) + n_grid_cut) cycle
          if (ig_domain(2,i) < ig_min10(2) - n_grid_cut) cycle
          if (ig_domain(2,i) > ig_max10(2) + n_grid_cut) cycle
          if (ig_domain(3,i) < ig_min10(3) - n_grid_cut) cycle
          if (ig_domain(3,i) > ig_max10(3) + n_grid_cut) cycle
          emfit_atom_local(i) = .true.

          coord_tmp(1:3,i) = coord_domain1(1:3,i)
          do j = 1, 3
            ig_lower(j,i) = ig_domain(j,i) - n_grid_cut
            ig_upper(j,i) = ig_domain(j,i) + n_grid_cut
            if (ig_max10(j) < ig_upper(j,i)) ig_upper(j,i) = ig_max10(j)
            if (ig_min10(j) > ig_lower(j,i)) ig_lower(j,i) = ig_min10(j)
          end do
        end do
        !$omp end parallel do

        simulated_map(ig_min10(1):ig_max10(1), &
                      ig_min10(2):ig_max10(2), &
                      ig_min10(3):ig_max10(3)) = 0.0_wp

      end if


      ! calculate error and exp functions and differences for integral
      ! and save only for the grid within the cutoff
      !   f: coefficient for x for erf calculation
      !   g: coefficient for y for exp calculation
      !
      f = sqrt(3.0_wp / (2.0_wp * sigma**2))
      g = -3.0_wp / (2.0_wp * sigma**2)

      !$omp parallel do &
      !$omp private (i, j, d)
      !
      do i = 1, natom_domain
        if (emfit_atom_local(i)) then

          do j = -n_grid_cut, n_grid_cut+1
            d = bound_x(ig_domain(1,i)+j) - coord_tmp(1,i)
            erfa(j,i) = erf(f*d)
            expa(j,i) = exp(g*d*d)
          end do
          do j = -n_grid_cut, n_grid_cut
            derfa_x(j,i) = erfa(j+1,i) - erfa(j,i)
            dexpa_x(j,i) = expa(j+1,i) - expa(j,i)
          end do

          do j = -n_grid_cut, n_grid_cut+1
            d = bound_y(ig_domain(2,i)+j) - coord_tmp(2,i)
            erfa(j,i) = erf(f*d)
            expa(j,i) = exp(g*d*d)
          end do
          do j = -n_grid_cut, n_grid_cut
            derfa_y(j,i) = erfa(j+1,i) - erfa(j,i)
            dexpa_y(j,i) = expa(j+1,i) - expa(j,i)
          end do

          do j = -n_grid_cut, n_grid_cut+1
            d = bound_z(ig_domain(3,i)+j) - coord_tmp(3,i)
            erfa(j,i) = erf(f*d)
            expa(j,i) = exp(g*d*d)
          end do
          do j = -n_grid_cut, n_grid_cut
            derfa_z(j,i) = erfa(j+1,i) - erfa(j,i)
            dexpa_z(j,i) = expa(j+1,i) - expa(j,i)
          end do

        end if
      end do
      !$omp end parallel do


      ! calculate simulated density map
      !
      !$omp parallel &
      !$omp private(id, m, n, k, j, i, zfactor, yzfactor, i1, j1, k1)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      do n = 1, natom_domain
        if (emfit_atom_local(n)) then
          do k = id+ig_lower(3,n), ig_upper(3,n), nthread
            k1 = k - ig_domain(3,n)
            zfactor = coeff_rho*derfa_z(k1,n)
            do j = ig_lower(2,n), ig_upper(2,n)
              j1 = j - ig_domain(2,n)
              yzfactor = derfa_y(j1,n)*zfactor
              do i = ig_lower(1,n), ig_upper(1,n)
                i1 = i - ig_domain(1,n)
                simulated_map(i,j,k) = simulated_map(i,j,k) + derfa_x(i1,n)*yzfactor
              end do
            end do
          end do
          !$omp barrier
        end if
      end do
      !$omp end parallel


      ! calculate two nominators for each atom first
      !
      dot_exp_drhodx(:,:) = 0.0_wp
      dot_sim_drhodx(:,:) = 0.0_wp

      !$omp parallel do schedule (dynamic)                     &
      !$omp private (n, exp_tmp, sim_tmp, k, j, i,             &
      !$omp          fact1, fact2, fact3, drho_dx, i1, j1, k1)
      !
      do n = 1, natom_domain
        if (emfit_atom_local(n)) then

          exp_tmp(1:3) = 0.0_wp
          sim_tmp(1:3) = 0.0_wp

          do k = ig_lower(3,n), ig_upper(3,n)
            k1 = k - ig_domain(3,n)
            do j = ig_lower(2,n), ig_upper(2,n)
              j1 = j - ig_domain(2,n)

              fact1 = derfa_y(j1,n) * derfa_z(k1,n)
              fact2 = dexpa_y(j1,n) * derfa_z(k1,n)
              fact3 = dexpa_z(k1,n) * derfa_y(j1,n)

              do i = ig_lower(1,n), ig_upper(1,n)
                i1 = i - ig_domain(1,n)
                drho_dx(1) = dexpa_x(i1,n) * fact1
                drho_dx(2) = derfa_x(i1,n) * fact2
                drho_dx(3) = derfa_x(i1,n) * fact3
                exp_tmp(1:3) = exp_tmp(1:3) + target_map   (i,j,k)*drho_dx(1:3)
                sim_tmp(1:3) = sim_tmp(1:3) + simulated_map(i,j,k)*drho_dx(1:3)
              end do

            end do
          end do

          dot_exp_drhodx(1:3,n) = exp_tmp(1:3) * coeff_drho
          dot_sim_drhodx(1:3,n) = sim_tmp(1:3) * coeff_drho
        end if

      end do
      !$omp end parallel do


      if (idomain == 1) then
        dot_exp_drhodx_tmp(1:3,1:natom_domain) = dot_exp_drhodx(1:3,1:natom_domain)
        dot_sim_drhodx_tmp(1:3,1:natom_domain) = dot_sim_drhodx(1:3,1:natom_domain)
      end if

    END DO


    norm_sim    = 0.0_wp
    dot_exp_sim = 0.0_wp

    !$omp parallel do                                   &
    !$omp private(i, j, k)                              &
    !$omp reduction(+:norm_sim) reduction(+:dot_exp_sim)
    !
    do k = ig_min00(3), ig_max00(3)
      do j = ig_min00(2), ig_max00(2)
        do i = ig_min00(1), ig_max00(1)
          norm_sim = norm_sim + simulated_map(i,j,k)**2
          dot_exp_sim = dot_exp_sim + target_map(i,j,k) * simulated_map(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do                                   &
    !$omp private(i, j, k)                              &
    !$omp reduction(+:norm_sim) reduction(+:dot_exp_sim)
    !
    do k = ig_min10(3), ig_max10(3)
      do j = ig_min10(2), ig_max10(2)
        do i = ig_min10(1), ig_max10(1)
          norm_sim = norm_sim + simulated_map(i,j,k)**2
          dot_exp_sim = dot_exp_sim + target_map(i,j,k) * simulated_map(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do


    call timer(TimerEmfit, TimerOff)
    call timer(TimerEnergyEmfit, TimerOff)

#ifdef HAVE_MPI_GENESIS
    before_allreduce(1) = norm_sim
    before_allreduce(2) = dot_exp_sim

    call mpi_allreduce(before_allreduce, after_allreduce, 2,  &
                       mpi_wp_real,  mpi_sum,                 &
                       mpi_comm_city, ierror)

    norm_sim    = after_allreduce(1)
    dot_exp_sim = after_allreduce(2)
#endif

    norm_sim  = sqrt(norm_sim)
    corrcoeff = dot_exp_sim / (norm_exp*norm_sim)
    experiments%emfit%corrcoeff = corrcoeff
    emcorr = corrcoeff
    eemfit = weight * (1.0_wp - corrcoeff)
    corrcoeff_save = corrcoeff

    ! calculate force of atoms
    !
    if (calc_force) then

      call timer(TimerEmfit, TimerOn)
      call timer(TimerEnergyEmfit, TimerOn)

      emfit_force_tmp(:,:) = 0.0_wp
      fact1 = 1.0_wp / norm_exp
      fact2 = 1.0_wp / norm_sim
      fact3 = dot_exp_sim/norm_sim**3

      !$omp parallel do &
      !$omp private (i00, m, i, ix, dcc_dx, dcc_dy, dcc_dz)
      !
      do i00 = 1, natom_domain00
        m  = ematom_l2g(0,i00)
        i  = id_g2l(m) 

        dcc_dx = fact1*(dot_exp_drhodx (1,i00)*fact2 - dot_sim_drhodx (1,i00)*fact3)
        dcc_dy = fact1*(dot_exp_drhodx (2,i00)*fact2 - dot_sim_drhodx (2,i00)*fact3)
        dcc_dz = fact1*(dot_exp_drhodx (3,i00)*fact2 - dot_sim_drhodx (3,i00)*fact3)

        emfit_force_tmp(1,m) = weight*dcc_dx
        emfit_force_tmp(2,m) = weight*dcc_dy
        emfit_force_tmp(3,m) = weight*dcc_dz

        force(i,1,1) = force(i,1,1) + emfit_force_tmp(1,m)
        force(i,2,1) = force(i,2,1) + emfit_force_tmp(2,m)
        force(i,3,1) = force(i,3,1) + emfit_force_tmp(3,m)
      end do
      !$omp end parallel do


      !$omp parallel do &
      !$omp private (j, dcc_dx, dcc_dy, dcc_dz)
      !
      do j = 1, natom_domain10
        dcc_dx = fact1*(dot_exp_drhodx_tmp(1,j)*fact2 - dot_sim_drhodx_tmp(1,j)*fact3)
        dcc_dy = fact1*(dot_exp_drhodx_tmp(2,j)*fact2 - dot_sim_drhodx_tmp(2,j)*fact3)
        dcc_dz = fact1*(dot_exp_drhodx_tmp(3,j)*fact2 - dot_sim_drhodx_tmp(3,j)*fact3)

        force_domain10(1,j) = weight*dcc_dx
        force_domain10(2,j) = weight*dcc_dy
        force_domain10(3,j) = weight*dcc_dz
      end do
      !$omp end parallel do

      call timer(TimerEmfit, TimerOff)
      call timer(TimerEnergyEmfit, TimerOff)

#ifdef HAVE_MPI_GENESIS
      if (my_sorted_rank <= nproc_city/2) then
        call mpi_isend(force_domain10, 3*natom_domain10, mpi_wp_real, &
                       domain_index1,  1, mpi_comm_city, irequest6, ierror)
        call mpi_irecv(force_domain01, 3*natom_domain01, mpi_wp_real, &
                       domain_index1,  1, mpi_comm_city, irequest7, ierror)
      else
        call mpi_irecv(force_domain01, 3*natom_domain01, mpi_wp_real, &
                       domain_index1,  1, mpi_comm_city, irequest6, ierror)
        call mpi_isend(force_domain10, 3*natom_domain10, mpi_wp_real, &
                       domain_index1,  1, mpi_comm_city, irequest7, ierror)
      end if
      call mpi_wait(irequest6, istatus, ierror)
      call mpi_wait(irequest7, istatus, ierror)
#endif

      call timer(TimerEmfit, TimerOn)
      call timer(TimerEnergyEmfit, TimerOn)

      !$omp parallel do &
      !$omp private (i01, m, i, ix)
      !
      do i01 = 1, natom_domain01
        m  = ematom_l2g(1,i01)
        i  = id_g2l(m) 

        force(i,1,1) = force(i,1,1) + force_domain01(1,i01)
        force(i,2,1) = force(i,2,1) + force_domain01(2,i01)
        force(i,3,1) = force(i,3,1) + force_domain01(3,i01)

        emfit_force_tmp(1,m) = emfit_force_tmp(1,m) + force_domain01(1,i01)
        emfit_force_tmp(2,m) = emfit_force_tmp(2,m) + force_domain01(2,i01)
        emfit_force_tmp(3,m) = emfit_force_tmp(3,m) + force_domain01(3,i01)
      end do
      !$omp end parallel do

      call timer(TimerEmfit, TimerOff)
      call timer(TimerEnergyEmfit, TimerOff)


      ! allreduce emfit_force
      !
      if (experiments%emfit%emfit_period >= 2) then
#ifdef HAVE_MPI_GENESIS
        emfit_force(:,:) = 0.0_wp
        ncycle = (natom_all - 1) / mpi_drain + 1
        nlen   = mpi_drain
        ixx    = 1

        do icycle = 1, ncycle
          if (icycle == ncycle) nlen = natom_all - (ncycle-1) * mpi_drain
          call mpi_allreduce(emfit_force_tmp(1,ixx), emfit_force(1,ixx), 3*nlen, &
                          mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
          ixx = ixx + nlen
        end do
#endif
      end if

    end if

    return

  end subroutine compute_energy_experimental_restraint_emfit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quicksort
  !> @brief        Provides a sorting function.
  !! @authors      NT
  !! @param[inout] a     : list of integers
  !! @param[in]    start : start index
  !! @param[in]    end   : end index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  recursive subroutine quicksort(a, start, end)

    ! formal arguments
    integer,                 intent(inout) :: a(*)
    integer,                 intent(in)    :: start
    integer,                 intent(in)    :: end

    ! local variables
    integer                  :: i, j, t, x


    x = a((start + end) / 2)
    i = start
    j = end

    do
      do while (a(i) < x)
        i = i + 1
      end do

      do while (x < a(j))
        j = j - 1
      end do

      if (i >= j) exit

      t = a(i)
      a(i) = a(j)
      a(j) = t
      i = i + 1
      j = j - 1

    end do

    if (start < i - 1) &
      call quicksort(a, start, i - 1)
    if (j + 1 < end) &
      call quicksort(a, j + 1, end)

    return

  end subroutine quicksort

end module cg_experiments_mod
