!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   density_analyze_mod
!> @brief   run density analysis
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module density_analyze_mod

  use density_output_mod
  use density_option_str_mod
  use sa_domain_str_mod
  use sa_domain_mod
  use sa_boundary_str_mod
  use sa_boundary_mod
  use sa_tool_mod
  use sa_ensemble_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use fitting_mod
  use measure_mod
  use select_atoms_mod
  use fitting_str_mod
  use molecules_str_mod
  use select_atoms_str_mod
  use fileio_trj_mod
  use fileio_mod
  use measure_mod
  use atom_libs_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use string_mod
#ifdef DEBUG
  use molecules_mod
  use fileio_pdb_mod
#endif
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: run_density
  private :: setup_density
  private :: check_options
  private :: resetup_boundary

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyze trajectories to caluclate density distribution
  !! @authors      DM, IY
  !! @param[inout] molecule       : molecule information
  !! @param[in]    trj_list       : trajectory file list information
  !! @param[in]    output         : output information
  !! @param[inout] option         : option information
  !! @param[inout] fitting        : fitting information
  !! @param[in]    density_option : density option information
  !! @param[inout] trajectory     : trajectory information
  !! @param[inout] boundary       : boundary information
  !! @param[inout] domain         : domain information
  !! @param[in]    ensemble       : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_density(molecule, trj_list, output, option, fitting, &
                         density_option, trajectory, boundary, domain, ensemble)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option), target,  intent(inout) :: option
    type(s_fitting),         intent(inout) :: fitting
    type(s_density_option),  intent(in)    :: density_option
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble

    ! local variables
    type(s_trj_file)            :: trj_in
    type(s_density)             :: solvent_density
    real(wp)                    :: dist, solute_crd(3), solvent_crd(3)
    real(wp)                    :: box_size_orig(3)
    real(wp)                    :: com_now(3), move(3)
    integer                     :: ifile, istep, totstep, nfile
    integer                     :: alloc_stat, dealloc_stat
    integer                     :: inside, icount
    integer                     :: alloc_step
    integer                     :: icell, ncell
    integer                     :: iatom, jatom, atom_no, ig
    integer                     :: igroup, ngroup
    integer                     :: natom_selec
    integer                     :: gid_solute, gid_solvent
    integer                     :: OW_pos
    character(len=6)            :: atom_type

    integer                     :: min_nx, max_nx
    integer                     :: min_ny, max_ny
    integer                     :: min_nz, max_nz
    integer                     :: ix, iy, iz

#ifdef DEBUG
    type(s_pdb)                 :: check_trj
#endif

    type(s_parray), allocatable :: selec_atom(:)
    real(dp),       allocatable :: sum_value(:)
    real(dp),       allocatable :: density_1d(:)
    real(wp),       allocatable :: atom_charge(:)
    integer,        allocatable :: natom_group(:), offset(:)
    integer,        allocatable :: solvent_in_domain(:), solute_in_bound(:)


    ! check options
    !
    if (main_rank) then
      call check_options(output, boundary, option, density_option)
    end if


    ! re-setup boundary
    !
    call resetup_boundary(density_option, option, molecule, boundary, domain)


    ! keep input values
    !
    box_size_orig(1) = boundary%box_size_x
    box_size_orig(2) = boundary%box_size_y
    box_size_orig(3) = boundary%box_size_z

    nfile       = size(trj_list%md_steps)
    ngroup      = size(option%analysis_atoms)
    gid_solute  = density_option%gid_solute
    gid_solvent = density_option%gid_solvent


    ! allocation
    !
    allocate(selec_atom(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    allocate(natom_group(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    allocate(offset(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc


    ! calculate the number of electrons of each atom
    !
    allocate(atom_charge(molecule%num_atoms), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    do iatom = 1, molecule%num_atoms
      atom_type = adjustl(molecule%atom_cls_name(iatom))

      if ((molecule%residue_name(iatom)(1:3) .eq. 'TIP') .or.  &
          (molecule%residue_name(iatom)(1:3) .eq. 'WAT') .or.  &
          (molecule%residue_name(iatom)(1:3) .eq. 'SPC') .or.  &
          (molecule%residue_name(iatom)(1:3) .eq. 'SOL')) then

        if (atom_type(1:1) .eq. 'O') then
          OW_pos = iatom
        end if

        ! when TIP4P/TIP5P water model is used
        !
        if (atom_type(1:2) .eq. 'LP' .or.  &
            atom_type(1:2) .eq. 'EP' .or.  &
            atom_type(1:2) .eq. 'IW')  then
          atom_charge(iatom) = 0

          ! atom_charge of lone pairs are put on to the water oxygen atom
          !
          atom_charge(OW_pos) = atom_charge(OW_pos) - molecule%charge(iatom)

        else
          atom_charge(iatom) = real(atomic_number_by_name(atom_type), wp) - &
                               molecule%charge(iatom)
        end if

      else
        atom_charge(iatom) = real(atomic_number_by_name(atom_type), wp) - &
                             molecule%charge(iatom)
      end if
    end do

    offset(:) = 0
    natom_selec = 0
    do igroup = 1, ngroup
      selec_atom(igroup)%idx => option%analysis_atoms(igroup)%idx
      natom_group(igroup) = size(selec_atom(igroup)%idx)
      offset(igroup) = natom_selec
      natom_selec = natom_selec + natom_group(igroup)
    end do

    if (main_rank .and. density_option%verbose) then
      write(MsgOut,'(/,(A))') '(information of density output target atoms)'
      write(MsgOut,'(A)')     '    number atom resnam resno  seg atom_num   charge num_elec     mass'
      do ig = 1, natom_group(gid_solvent)
        iatom = selec_atom(gid_solvent)%idx(ig)
        atom_type = adjustl(molecule%atom_cls_name(iatom))

        write(MsgOut,'(i10,1x,a4,1x,a6,1x,i5,1x,a4,1x,i8,2(1x,F8.3))')  &
          molecule%atom_no(iatom),      &
          molecule%atom_name(iatom),    &
          molecule%residue_name(iatom), &
          molecule%residue_no(iatom),   &
          molecule%segment_name(iatom), &
          atomic_number_by_name(atom_type),     &
          molecule%charge(iatom), &
          atom_charge(iatom)
      end do
      write(MsgOut,'(A)') '(end of information)'
      write(MsgOut, '(A)') ''
    end if

    ! centering reference coordinates
    !
    do iatom = 1, molecule%num_atoms
      molecule%atom_coord(1:3,iatom) = molecule%atom_refcoord(1:3,iatom)
    end do

    if (density_option%recenter /= 0) then
      call get_cofm(molecule, selec_atom, density_option%recenter, com_now)
      move(1:3) = -com_now(1:3)
      call shift_molecule(move, molecule)
    end if

    call wrap_molecule(option, trajectory, ensemble, &
                       boundary, molecule)

    do iatom = 1, molecule%num_atoms
      molecule%atom_refcoord(1:3,iatom) = molecule%atom_coord(1:3,iatom)
    end do


    ! begin the analysis
    !
    totstep = 0
    domain%num_atom_selec = natom_selec
    domain%num_group      = ngroup

    ! analysis loop
    !
    do ifile = 1, nfile

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord)
        !
        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(ifile)) == 0) then

          totstep = totstep + 1
          if (main_rank) then
            write(MsgOut,*) '      number of structures = ', totstep
          end if

          do iatom = 1, molecule%num_atoms
            molecule%atom_coord(1:3, iatom) = trajectory%coord(1:3, iatom)
          end do

#ifdef DEBUG
          if (main_rank) then
!            call export_molecules(molecule, pdb = check_trj)
!            call output_pdb('input_snapshot.pdb', check_trj)
          end if
#endif

          ! centering coordinates
          !
          if (density_option%recenter /= 0) then
            call get_cofm(molecule, selec_atom, density_option%recenter,com_now)
            move(1:3) = -com_now(1:3)
            call shift_molecule(move, molecule)
          end if

          ! wrap coordinates
          !
          call wrap_molecule(option, trajectory, ensemble, &
                             boundary, molecule)

          ! fitting
          !
          if (fitting%mass_weight) then
            call run_fitting(fitting,                 &
                             molecule%atom_refcoord, &
                             molecule%atom_coord,     &
                             molecule%atom_coord,     &
                             molecule%mass)

          else
            call run_fitting(fitting,                 &
                             molecule%atom_refcoord, &
                             molecule%atom_coord,     &
                             molecule%atom_coord)
          end if

          ! setup domains
          !
          call setup_atom(option, ensemble, &
                          trajectory, selec_atom, molecule, boundary, domain)

#ifdef DEBUG
          if (main_rank) then
            call export_molecules(molecule, pdb = check_trj)
            call output_pdb('after_wrap.pdb', check_trj)
          end if
#endif

          ! determinie solvent atoms in domain
          !
          do alloc_step = 1, 2
            icount = 0

            do iatom = 1, natom_group(gid_solvent)
              icell = domain%selec_atom2cell(iatom + offset(gid_solvent))

              if (icell > 0) then
                inside = domain%cell_g2l(icell)

                if (inside > 0) then
                  icount = icount + 1

                  if (alloc_step == 2) then
                    atom_no = selec_atom(gid_solvent)%idx(iatom)
                    solvent_in_domain(icount) = atom_no
                  end if
                end if
              end if
            end do

            if (alloc_step == 1) then
              if (allocated(solvent_in_domain)) then
                deallocate(solvent_in_domain, stat = dealloc_stat)
                if (dealloc_stat /= 0) &
                  call error_msg_dealloc

              end if

              allocate(solvent_in_domain(icount), stat = alloc_stat)
              if (alloc_stat /= 0) &
                call error_msg_alloc

              solvent_in_domain(:) = 0
            end if
          end do

#ifdef DEBUG
          if (totstep == 1) then
!            write(MsgOut,*) 'solvent', my_world_rank, size(solvent_in_domain)
          end if
#endif

          ! determinie solute atoms in domain + boundary
          !
          do alloc_step = 1, 2
            icount = 0

            do iatom = 1, natom_group(gid_solute)
              icell = domain%selec_atom2cell(iatom + offset(gid_solute))

              if (icell > 0) then
                inside = domain%cell_g2l(icell) + domain%cell_g2b(icell)

                if (inside > 0) then
                  icount = icount + 1

                  if (alloc_step == 2) then
                    atom_no = selec_atom(gid_solute)%idx(iatom)
                    solute_in_bound(icount) = atom_no
                  end if
                end if
              end if
            end do

            if (alloc_step == 1) then
              if (allocated(solute_in_bound)) then
                deallocate(solute_in_bound, stat = dealloc_stat)
                if (dealloc_stat /= 0) &
                  call error_msg_dealloc

              end if

              allocate(solute_in_bound(icount), stat = alloc_stat)
              if (alloc_stat /= 0) &
                call error_msg_alloc

              solute_in_bound(:) = 0
            end if
          end do

#ifdef DEBUG
!          write(MsgOut,*) 'solute', my_world_rank, size(solute_in_bound)
#endif

          ! prepare voxels for density output
          !
          if (totstep == 1) then
            if (option%determine_box == DetermineBoxManual) then
              boundary%box_size_x = box_size_orig(1)
              boundary%box_size_y = box_size_orig(2)
              boundary%box_size_z = box_size_orig(3)
            end if

            call setup_density(molecule, boundary, molecule%atom_coord, &
                               selec_atom(gid_solute), density_option,  &
                               solvent_density)

            call alloc_density(solvent_density)

            allocate(sum_value(boundary%num_cells_x * &
                               boundary%num_cells_y * &
                               boundary%num_cells_z), &
                               stat = alloc_stat)

            if (alloc_stat /= 0) &
               call error_msg_alloc

            sum_value(:) = 0.0_dp

            allocate(density_1d(size(sum_value)), stat = alloc_stat)
            if (alloc_stat /= 0) &
              call error_msg_alloc

            density_1d(:) = 0.0_dp
          end if

#ifdef DEBUG
          if (main_rank) then
!            call export_molecules(molecule, pdb = check_trj)
!            call output_pdb('after_setup_density.pdb', check_trj)
          end if
#endif

          ! calculate density distribution
          !

          !$omp parallel do default(shared) &
          !$omp private(atom_no, icell, jatom, ig, solvent_crd, solute_crd, &
          !$omp         dist)
          !

          do iatom = 1, size(solvent_in_domain)
            atom_no = solvent_in_domain(iatom)
            solvent_crd(1:3) = molecule%atom_coord(1:3, atom_no)

            icell = domain%selec_atom2cell(domain%id_global2selec(atom_no));
            if (icell == 0) &
              cycle

            do jatom = 1, size(solute_in_bound)
              ig = solute_in_bound(jatom)
              solute_crd(1:3) = molecule%atom_coord(1:3, ig)

              dist = compute_dis(solute_crd, solvent_crd)

              if (dist <= density_option%ana_range) then

                if (density_option%density_type == DensityTypeNumber) then
                  sum_value(icell) = sum_value(icell) + 1.0_dp

                else if (density_option%density_type== DensityTypeElectron) then
                  sum_value(icell) = sum_value(icell) + &
                                                     dble(atom_charge(atom_no))
                end if

                exit
              end if

            end do

          end do

          !
          !$omp end parallel do

        end if
      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    ! reduce
    !
    ncell = size(sum_value)

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(sum_value, density_1d, ncell, &
                    MPI_REAL8, mpi_sum, 0, mpi_comm_country, ierror)
#else

    density_1d(:) = sum_value(:)
#endif

    if (main_rank) then

#ifdef DEBUG
      write(0,'(e15.8)') epsilon(1.0_dp)
#endif

      density_1d(:) = density_1d(:) / dble(totstep) * dble(density_option%magnification)
      solvent_density%value(:,:,:) = 0.0_dp
      solvent_density%value(:,:,:) = reshape(density_1d, shape(solvent_density%value))


      ! extract non-zero density region
      !
      iz_loop_1: do iz = 1, solvent_density%ngrid_z
                   do iy = 1, solvent_density%ngrid_y
                   do ix = 1, solvent_density%ngrid_x
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0d0)) then
                       min_nz = iz
                       exit iz_loop_1
                     end if
                   end do
                   end do
                 end do iz_loop_1

      iz_loop_2: do iz = solvent_density%ngrid_z, 1, -1
                   do iy = solvent_density%ngrid_y, 1, -1
                   do ix = solvent_density%ngrid_x, 1, -1
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0d0)) then
                       max_nz = iz
                       exit iz_loop_2
                     end if
                   end do
                   end do
                 end do iz_loop_2

      iy_loop_1: do iy = 1, solvent_density%ngrid_y
                   do iz = min_nz, max_nz
                   do ix = 1, solvent_density%ngrid_x
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0_dp)) then
                       min_ny = iy
                       exit iy_loop_1
                     end if
                   end do
                   end do
                 end do iy_loop_1

      iy_loop_2: do iy = solvent_density%ngrid_y, 1, -1
                   do iz = max_nz, min_nz, -1
                   do ix = solvent_density%ngrid_x, 1, -1
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0_dp)) then
                       max_ny = iy
                       exit iy_loop_2
                     end if
                   end do
                   end do
                 end do iy_loop_2

      ix_loop_1: do ix = 1, solvent_density%ngrid_x
                   do iz = min_nz, max_nz
                   do iy = min_ny, max_ny
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0_dp)) then
                       min_nx = ix
                       exit ix_loop_1
                     end if
                   end do
                   end do
                 end do ix_loop_1

      ix_loop_2: do ix = solvent_density%ngrid_x, 1, -1
                   do iz = max_nz, min_nz, -1
                   do iy = max_ny, min_ny, -1
                     if (abs(solvent_density%value(ix, iy, iz)) > epsilon(1.0_dp)) then
                       max_nx = ix
                       exit ix_loop_2
                     end if
                   end do
                   end do
                 end do ix_loop_2

      solvent_density%nxmin = min_nx
      solvent_density%nxmax = max_nx
      solvent_density%nymin = min_ny
      solvent_density%nymax = max_ny
      solvent_density%nzmin = min_nz
      solvent_density%nzmax = max_nz

      solvent_density%ngrid_x = max_nx - min_nx + 1
      solvent_density%ngrid_y = max_ny - min_ny + 1
      solvent_density%ngrid_z = max_nz - min_nz + 1

      ! output the result
      !
      do iatom = 1, molecule%num_atoms
        molecule%atom_coord(1:3,iatom) = molecule%atom_refcoord(1:3,iatom)
      end do
      call write_density(output, density_option, solvent_density, molecule)
    end if

    ! deallocation
    !
    call dealloc_density(solvent_density)
    deallocate(sum_value)
    deallocate(atom_charge)
    deallocate(offset)
    deallocate(natom_group)
    deallocate(selec_atom)
    deallocate(solute_in_bound, solvent_in_domain)

    return

  end subroutine run_density

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_density
  !> @brief        setup density information
  !! @authors      DM
  !! @param[in]    molecule      : molecule information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : atom coordinates
  !! @param[in]    selec_atom    : atom selection information
  !! @param[in]    option        : option information
  !! @param[in]    determine_box : method to determine boundary
  !! @param[out]   density       : density information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_density(molecule, boundary, coord, selec_atom, option, &
                           density)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    type(s_parray),          intent(in)    :: selec_atom
    type(s_density_option),  intent(in)    :: option
    type(s_density),         intent(out)   :: density

    ! local variables
    real(wp)                 :: solute_com(3)
    real(wp)                 :: xmin, xmax
    real(wp)                 :: ymin, ymax
    real(wp)                 :: zmin, zmax


    solute_com = compute_com(coord, molecule%mass, selec_atom%idx)

    xmin = solute_com(1) - boundary%box_size_x * half
    xmax = solute_com(1) + boundary%box_size_x * half
    ymin = solute_com(2) - boundary%box_size_y * half
    ymax = solute_com(2) + boundary%box_size_y * half
    zmin = solute_com(3) - boundary%box_size_z * half
    zmax = solute_com(3) + boundary%box_size_z * half

    density%nxmin = nint(xmin / boundary%cell_size_x)
    density%nxmax = nint(xmax / boundary%cell_size_x)-1
    density%nymin = nint(ymin / boundary%cell_size_y)
    density%nymax = nint(ymax / boundary%cell_size_y)-1
    density%nzmin = nint(zmin / boundary%cell_size_z)
    density%nzmax = nint(zmax / boundary%cell_size_z)-1

    density%ngrid_x = boundary%num_cells_x
    density%ngrid_y = boundary%num_cells_y
    density%ngrid_z = boundary%num_cells_z

    density%voxel_size_x = boundary%cell_size_x
    density%voxel_size_y = boundary%cell_size_y
    density%voxel_size_z = boundary%cell_size_z

    return

  end subroutine setup_density 

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_options
  !> @brief        check options for this program
  !! @authors      DM
  !! @param[in]    output         : output information
  !! @param[in]    boundary       : boundary information
  !! @param[in]    option         : option information
  !! @param[in]    density_option : density_option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_options(output, boundary, option, density_option)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_boundary),        intent(in) :: boundary
    type(s_option),          intent(in) :: option
    type(s_density_option),  intent(in) :: density_option

    ! local variable
    logical                  :: file_exist
    logical                  :: find_error


    find_error = .false.

    ! output file check
    !
    inquire(file=output%mapfile, exist=file_exist)
    if (file_exist) then
      write(MsgOut,'(a)') &
        'Check_Options> the output file ' //trim(output%mapfile) // &
        ' already exist.'
      find_error = .true.
    end if

    ! parameter check
    !
    if (boundary%type == BoundaryTypeNOBC) then
      if (option%determine_box == DetermineBoxTrajectory) then
        write(MsgOut,'(a)') &
          'Check_Options> TRAJECTORY option is available only in PBC'
        find_error = .true.
      end if

    else if (boundary%type == BoundaryTypePBC) then
      if (option%determine_box == DetermineBoxMax) then
        write(MsgOut,'(a)') &
          'Check_Options> MAX option is available only in NOBC'
        find_error = .true.
      end if
    end if

    if (density_option%ana_range > option%buffer) then
      write(MsgOut,'(a)') &
        'Check_Options> range for density calculation should be larger than buffer in [SPANA_OPTION]'
      find_error = .true.
    end if

    ! check processor number based on the num of domains
    !
    if (boundary%num_domain(1) /= 0 .and. &
        boundary%num_domain(2) /= 0 .and. &
        boundary%num_domain(3) /= 0) then

      if (product(boundary%num_domain) /= nproc_country) then
        write(MsgOut,'(a)') &
          'Check_Options> # of process is not domain_x * domain_y * domain_z '
        find_error = .true.
      end if
    end if

    if (find_error) then
      write(MsgOut,'(a)')'Check_Options> error detected'
      call error_msg(&
        'Check_Options> analysis terminated because of option error')
    end if

    return

  end subroutine check_options

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    resetup_boundary
  !> @brief        re-setup boundary by voxcel information
  !! @authors      NT
  !! @param[in]    density_option : density_option information
  !! @param[in]    option         : option information
  !! @param[in]    molecule       : molecule information
  !! @param[inout] boundary       : boundary information
  !! @param[inout] domain         : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine resetup_boundary(density_option, option, molecule, &
                              boundary, domain)

    ! formal arguments
    type(s_density_option),  intent(in)    :: density_option
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: num_cells_x
    integer                  :: num_cells_y
    integer                  :: num_cells_z
    integer                  :: remainder


    if (density_option%voxel_size == 0) &
      return

    num_cells_x = int(boundary%box_size_x / density_option%voxel_size)
    if (mod(num_cells_x, boundary%num_domain(1)) /= 0) then
      remainder = mod(num_cells_x, boundary%num_domain(1))
      num_cells_x = num_cells_x + boundary%num_domain(1) - remainder
    end if

    num_cells_y = int(boundary%box_size_y / density_option%voxel_size)
    if (mod(num_cells_y, boundary%num_domain(2)) /= 0) then
      remainder = mod(num_cells_y, boundary%num_domain(2))
      num_cells_y = num_cells_y + boundary%num_domain(2) - remainder
    end if

    num_cells_z = int(boundary%box_size_z / density_option%voxel_size)
    if (mod(num_cells_z, boundary%num_domain(3)) /= 0) then
      remainder = mod(num_cells_z, boundary%num_domain(3))
      num_cells_z = num_cells_z + boundary%num_domain(3) - remainder
    end if

    boundary%num_cells_x = num_cells_x
    boundary%num_cells_y = num_cells_y
    boundary%num_cells_z = num_cells_z
    boundary%cell_size_x = boundary%box_size_x / real(num_cells_x, wp)
    boundary%cell_size_y = boundary%box_size_y / real(num_cells_y, wp)
    boundary%cell_size_z = boundary%box_size_z / real(num_cells_z, wp)

    call dealloc_domain_all(domain)
    call setup_domain(boundary, molecule, domain, option)

    if (main_rank) then
      write(MsgOut,'(a)') 'Re-Setup_Boundary_Cell> '
      write(MsgOut,'(A24,3I10)') &
           '      ncells (x,y,z)  = ', boundary%num_cells_x, &
                                       boundary%num_cells_y, &
                                       boundary%num_cells_z
      write(MsgOut,'(A24,3F10.4)') &
           '   cell size (x,y,z)  = ', boundary%cell_size_x, &
                                       boundary%cell_size_y, &
                                       boundary%cell_size_z
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine resetup_boundary

end module density_analyze_mod
