!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rdf_analyze_mod
!> @brief   run rdf analysis
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rdf_analyze_mod

  use rdf_option_str_mod
  use sa_domain_mod
  use sa_boundary_mod
  use sa_tool_mod
  use sa_domain_str_mod
  use sa_boundary_str_mod
  use sa_ensemble_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use fileio_trj_mod
  use fileio_mod
  use molecules_str_mod
  use measure_mod
  use select_atoms_str_mod
  use sa_timers_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use string_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
  use omp_lib
#endif

  implicit none
  private

  ! subroutines
  public  :: run_rdf
  public  :: check_options
  private :: resetup_boundary

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_rdf
  !> @brief        run analyzing trajectories
  !! @authors      IY
  !! @param[inout] molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[in]    rdf_option : contact_option information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] domain     : domain information
  !! @param[in]    ensemble   : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_rdf(molecule, trj_list, output, option, &
                     rdf_option, trajectory, boundary, domain, ensemble)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option), target,  intent(inout) :: option
    type(s_rdf_option),      intent(in)    :: rdf_option
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble

    ! local variables
    type(s_trj_file)         :: trj_in
    real(wp)                 :: cell_px, cell_py, cell_pz
    real(wp)                 :: mindist, dist
    real(wp)                 :: cofm(3)
    real(wp)                 :: move(3)
    real(wp)                 :: binsize, bindist, bindist_in, bindist_out
    real(wp)                 :: csize_x, csize_y, csize_z, cvolume, acvolume
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: density_1, density_2, binvolume
    real(wp)                 :: avail_volume, ncoordi_atom
    real(wp)                 :: bulk_density_1, bulk_density_2
    real(wp)                 :: norm_density_1, norm_density_2
    real(wp)                 :: p_ix, p_iy, p_iz, p_jx, p_jy, p_jz
    integer                  :: ifile, nfile
    integer                  :: istep, totstep
    integer                  :: iatom, jatom, idx
    integer                  :: txt_out
    integer                  :: inside, alloc_step, icount
    integer                  :: nsolvatm_inbnd
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: icel_global, ncell, icell
    integer                  :: i_x, i_y, i_z
    integer                  :: cell_iposi(3)
    integer                  :: ibin, maxbin, bin_count
    integer                  :: natom_selec
    integer                  :: igroup,  ngroup
    integer                  :: gid_solvent, gid_solute
    integer                  :: offset_solute, offset_solvent
    integer                  :: ig, jg
    integer                  :: ithread
    integer                  :: id, my_id

    type(s_parray),          allocatable   :: selec_atom(:)
    real(wp),                allocatable   :: density_bin(:)
    real(wp),                allocatable   :: ave_density_bin(:)
    integer,                 allocatable   :: atoms_focus(:)
    integer,                 allocatable   :: atoms_solvent(:)
    integer,                 allocatable   :: ncell_bin(:,:)
    integer,                 allocatable   :: ncell_bin_allthread(:)
    integer,                 allocatable   :: ncell_bin_allproc(:)
    integer,                 allocatable   :: ave_ncell_bin(:)
    integer,                 allocatable   :: natom_bin(:,:)
    integer,                 allocatable   :: natom_bin_allthread(:)
    integer,                 allocatable   :: natom_bin_allproc(:)
    integer,                 allocatable   :: ave_natom_bin(:)


    acvolume = 0
    totstep  = 0

    bulk_density_1 = 0.0_wp
    bulk_density_2 = 0.0_wp
    bin_count      = 0

    ! check cntrl parameter
    !
    if (main_rank) then
      call check_options(molecule, boundary, domain, trj_list, &
                         option, rdf_option, totstep)
    end if

    ! re-setup boundary
    !
    call resetup_boundary(rdf_option, option, molecule, boundary, domain)

    nfile        = size(trj_list%md_steps)
    ngroup       = size(option%analysis_atoms)
    gid_solute   = rdf_option%gid_solute
    gid_solvent  = rdf_option%gid_solvent

    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate output arrays
    !
    binsize = rdf_option%binsize
    maxbin  = int(rdf_option%ana_range/rdf_option%binsize)
    allocate(ncell_bin(maxbin,nthread), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(ncell_bin_allthread(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(ave_ncell_bin(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(ncell_bin_allproc(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    allocate(natom_bin(maxbin,nthread), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(natom_bin_allthread(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(ave_natom_bin(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(natom_bin_allproc(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    allocate(density_bin(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc
    allocate(ave_density_bin(maxbin), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ncell_bin(1:maxbin,1:nthread)     = 0
    ncell_bin_allthread(1:maxbin)     = 0
    ave_ncell_bin(1:maxbin)           = 0
    ncell_bin_allproc(1:maxbin)       = 0
    natom_bin(1:maxbin,1:nthread)     = 0
    natom_bin_allthread(1:maxbin)     = 0
    natom_bin_allproc(1:maxbin)       = 0
    ave_natom_bin(1:maxbin)           = 0
    density_bin(1:maxbin)             = 0.0_wp
    ave_density_bin(1:maxbin)         = 0.0_wp

    allocate(selec_atom(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! open output file
    !
    if (main_rank) then
      if (output%txtfile .ne. '') then
        call open_file(txt_out, output%txtfile, IOFileOutputReplace)
      else
        txt_out = MsgOut
      end if
    end if

    natom_selec = 0
    do igroup = 1, ngroup
      selec_atom(igroup)%idx => option%analysis_atoms(igroup)%idx
      natom_selec = natom_selec + size(option%analysis_atoms(igroup)%idx)
    end do

    domain%num_atom_selec   = natom_selec
    domain%num_group        = ngroup

    offset_solute = 0
    if (gid_solute > 1) then
      do igroup = 1, gid_solute-1
        offset_solute = offset_solute + size(selec_atom(igroup)%idx)
      end do
    end if

    offset_solvent = 0
    if (gid_solvent > 1) then
      do igroup = 1, gid_solvent-1
        offset_solvent = offset_solvent + size(selec_atom(igroup)%idx)
      end do
    end if


    ! set timer
    !
    call timer(TimerAnalysis, TimerOn)


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

          ! refine atomcoodinate
          do iatom = 1, molecule%num_atoms
            molecule%atom_coord(1,iatom) = trajectory%coord(1,iatom)
            molecule%atom_coord(2,iatom) = trajectory%coord(2,iatom)
            molecule%atom_coord(3,iatom) = trajectory%coord(3,iatom)
          end do

          if (rdf_option%recenter /= 0) then
            call get_cofm(molecule, selec_atom, rdf_option%recenter, cofm)
            move(1) = -cofm(1)
            move(2) = -cofm(2)
            move(3) = -cofm(3)
            call shift_molecule(move, molecule)
          end if

          ! wrap molecule
          !
          call wrap_molecule(option, trajectory, ensemble, &
                             boundary, molecule)

          ! setup domains
          !
          call setup_atom(option, ensemble, &
                          trajectory, selec_atom, molecule, boundary, domain)

          ! check cntrl parameter not to include image atoms
          ! for NPT trajectory
          if (main_rank) then
            if (boundary%type == BoundaryTypePBC) then
              if (option%determine_box == DetermineBoxTrajectory) then
                if (ensemble%ensemble == EnsembleNPT  .or. &
                    ensemble%ensemble == EnsembleNPAT .or. &
                    ensemble%ensemble == EnsembleNPgT) then
                  if (main_rank) then
                    call check_options(molecule, boundary, domain, trj_list, &
                         option, rdf_option, totstep)
                  end if
                end if
              end if
            end if
          end if

          ! get volume of a cell
          bsize_x = boundary%box_size_x
          bsize_y = boundary%box_size_y
          bsize_z = boundary%box_size_z

          csize_x = bsize_x/boundary%num_cells_x
          csize_y = bsize_y/boundary%num_cells_y
          csize_z = bsize_z/boundary%num_cells_z

          cvolume = csize_x*csize_y*csize_z
          acvolume = acvolume + cvolume


          ! proximal distribution function
          !
          if (rdf_option%rmode == RModeProximal) then

            ! initialize arrays
            !
            ncell_bin(1:maxbin,1:nthread) = 0
            natom_bin(1:maxbin,1:nthread) = 0
            ncell_bin_allthread(1:maxbin) = 0
            natom_bin_allthread(1:maxbin) = 0

            ncell_bin_allproc(1:maxbin)   = 0
            natom_bin_allproc(1:maxbin)   = 0
            density_bin(1:maxbin)         = 0.0_wp

            ! determin the solute atom inside the boundary
            do alloc_step = 1, 2
              icount = 0

              do iatom = 1, size(selec_atom(gid_solute)%idx)
                idx = selec_atom(gid_solute)%idx(iatom)

                icel_global = domain%selec_atom2cell(iatom + offset_solute)

                if(icel_global > 0 )then
                  inside = domain%cell_g2l(icel_global)
                  inside = inside + domain%cell_g2b(icel_global)
                else
                  inside = 0
                end if

                if (inside > 0) then
                  icount = icount +1
                  if (alloc_step == 2) then
                    atoms_focus(icount) = idx
                  end if
                end if

              end do ! iatom

              if (alloc_step ==1) then
                if (allocated(atoms_focus)) then
                  deallocate(atoms_focus, stat=dealloc_stat)
                  if (dealloc_stat /= 0) call error_msg_dealloc
                end if
                allocate(atoms_focus(icount), stat = alloc_stat)
                if (alloc_stat /= 0) call error_msg_alloc
                atoms_focus(1:icount) = 0
              end if

            end do ! alloc_step

            ncell = domain%num_cell_local

            !$omp parallel default(shared)                                    &
            !$omp private(id, icell, my_id, i_x, i_y, i_z,cell_iposi,mindist) &
            !$omp private(ibin, cell_px, cell_py, cell_pz)
            !
#ifdef OMP
            id  = omp_get_thread_num()
#else
            id  = 0
#endif

            my_id = id
            do icell = my_id+1, ncell, nthread
            !do icell = 1, ncell
              i_x = domain%cell_l2gx_orig(icell)
              i_y = domain%cell_l2gy_orig(icell)
              i_z = domain%cell_l2gz_orig(icell)
              cell_iposi(1) =i_x
              cell_iposi(2) =i_y
              cell_iposi(3) =i_z

              !calculate the minimum dist from local cell to solute atom
              call mindist_cell2atom(molecule, boundary, &
                                     cell_iposi, atoms_focus,  mindist)

              ibin = int(mindist/binsize) + 1

              if (ibin <= maxbin) then
                ncell_bin(ibin, id+1) = ncell_bin(ibin,id+1)+1
                natom_bin(ibin, id+1) = natom_bin(ibin,id+1)+ &
                     domain%num_atom_group(gid_solvent,icell)

                cell_px = (real(cell_iposi(1),wp)-0.5_wp)*csize_x
                cell_py = (real(cell_iposi(2),wp)-0.5_wp)*csize_y
                cell_pz = (real(cell_iposi(3),wp)-0.5_wp)*csize_z

              end if

            end do
            !$omp end parallel

            ! combine the natom and ncell for all threads
            !
            do ithread =1, nthread
              do ibin = 1, maxbin
                ncell_bin_allthread(ibin) = ncell_bin_allthread(ibin) &
                                          + ncell_bin(ibin,ithread)

                natom_bin_allthread(ibin) = natom_bin_allthread(ibin) &
                                          + natom_bin(ibin,ithread)
              end do
            end do

            ! comnine the natom and ncell for all proc
            !
#ifdef HAVE_MPI_GENESIS
            call mpi_reduce(ncell_bin_allthread, ncell_bin_allproc, &
                            maxbin, mpi_integer, &
                            mpi_sum, 0, mpi_comm_country, ierror)

            call mpi_reduce(natom_bin_allthread, natom_bin_allproc, &
                            maxbin, mpi_integer, &
                            mpi_sum, 0, mpi_comm_country, ierror)
#else
            do ibin = 1, maxbin
              ncell_bin_allproc(ibin) = ncell_bin_allthread(ibin)
              natom_bin_allproc(ibin) = ncell_bin_allthread(ibin)
            end do
#endif
            ! calculate density profile for one snapshot
            ! and integrate ncell and natom for time average
            ! WARNNG 0 devided 0  ! we need  nsample(ibin)
            do ibin = 1, maxbin
              density_bin(ibin) = real(natom_bin_allproc(ibin),wp)/ &
                   (cvolume*ncell_bin_allproc(ibin))
              ave_density_bin(ibin) = ave_density_bin(ibin) + density_bin(ibin)

              ave_ncell_bin(ibin) = ave_ncell_bin(ibin) +ncell_bin_allproc(ibin)
              ave_natom_bin(ibin) = ave_natom_bin(ibin) +natom_bin_allproc(ibin)
            end do

          ! radial distribution function
          !
          else if (rdf_option%rmode == RModeRadial) then

            ! determin the solute atom inside the local domain
            do alloc_step = 1, 2
              icount = 0

              do iatom = 1, size(selec_atom(gid_solute)%idx)
                idx = selec_atom(gid_solute)%idx(iatom)
                icel_global = domain%selec_atom2cell(iatom + offset_solute)
                inside = domain%cell_g2l(icel_global)

                if (inside > 0) then
                  icount = icount +1
                  if (alloc_step == 2) then
                    atoms_focus(icount) = idx
                  end if
                end if

              end do ! iatom

              if (alloc_step == 1) then
                if(allocated(atoms_focus))then
                  deallocate(atoms_focus, stat=dealloc_stat)
                  if (dealloc_stat /= 0) call error_msg_dealloc
                end if
                allocate(atoms_focus(icount), stat = alloc_stat)
                if (alloc_stat /= 0) call error_msg_alloc
                atoms_focus(1:icount) = 0
              end if

            end do ! alloc_step

            ! determin the solvent atom inside the boundary
            !
            do alloc_step = 1, 2
              icount = 0

              do iatom = 1, size(selec_atom(gid_solvent)%idx)
                idx = selec_atom(gid_solvent)%idx(iatom)
                icel_global = domain%selec_atom2cell(iatom + offset_solvent)
                inside = domain%cell_g2l(icel_global)
                inside = inside + domain%cell_g2b(icel_global)

                if (inside > 0) then
                  icount = icount +1
                  if (alloc_step == 2) then
                    atoms_solvent(icount) = idx
                  end if
                end if

              end do ! iatom

              if (alloc_step ==1) then
                if(allocated(atoms_solvent))then
                  deallocate(atoms_solvent, stat=dealloc_stat)
                  if (dealloc_stat /= 0) call error_msg_dealloc
                end if
                allocate(atoms_solvent(icount), stat = alloc_stat)
                if (alloc_stat /= 0) call error_msg_alloc
                atoms_solvent(1:icount) = 0
              end if

            end do ! alloc_step

            ! calculate the distance between atoms_focus vs solvent
            !
            nsolvatm_inbnd = size(atoms_solvent)

            do iatom = 1, size(atoms_focus)
              !ig = selec_atom(gid_solvent)%idx(iatom)
              ig   = atoms_focus(iatom)
              p_ix = molecule%atom_coord(1,ig)
              p_iy = molecule%atom_coord(2,ig)
              p_iz = molecule%atom_coord(3,ig)

              !$omp parallel default(shared)      &
              !$omp private(id, jatom, my_id, jg) &
              !$omp private(p_jx, p_jy, p_jz, dist, ibin )
              !
#ifdef OMP
              id  = omp_get_thread_num()
#else
              id  = 0
#endif
              my_id = id
              do jatom = my_id+1, nsolvatm_inbnd, nthread
                jg = atoms_solvent(jatom)
                p_jx = molecule%atom_coord(1,jg)
                p_jy = molecule%atom_coord(2,jg)
                p_jz = molecule%atom_coord(3,jg)

                p_jx  =  p_jx  - bsize_x*anint((p_jx-p_ix)/bsize_x)
                p_jy  =  p_jy  - bsize_y*anint((p_jy-p_iy)/bsize_y)
                p_jz  =  p_jz  - bsize_z*anint((p_jz-p_iz)/bsize_z)

                dist = (p_ix - p_jx)*(p_ix - p_jx)+(p_iy - p_jy)*(p_iy - p_jy)+&
                       (p_iz - p_jz)*(p_iz - p_jz)
                dist = sqrt(dist)

                ibin = int(dist/binsize)+1

                if (ibin <= maxbin .and. dist > 0.0_wp) then
                  natom_bin(ibin, id+1) = natom_bin(ibin, id+1)+ 1
                end if

              end do ! jatom
              !$omp end parallel

            end do  ! iatom

          end if ! rdf_option%rmode

        end if ! istep = ana_period

      end do ! istep

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    ! average cell volume
    !
    acvolume = acvolume/totstep

    if (rdf_option%rmode == RModeRadial) then
      do ithread =1, nthread
        do ibin = 1, maxbin
          ncell_bin_allthread(ibin) = ncell_bin_allthread(ibin)&
                                    + ncell_bin(ibin, ithread)

          natom_bin_allthread(ibin) = natom_bin_allthread(ibin)&
                                    + natom_bin(ibin, ithread)
        end do
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(ncell_bin_allthread, ave_ncell_bin, &
                      maxbin, mpi_integer, &
                      mpi_sum, 0, mpi_comm_country, ierror)

      call mpi_reduce(natom_bin_allthread, ave_natom_bin, &
                      maxbin, mpi_integer, &
                      mpi_sum, 0, mpi_comm_country, ierror)
#else
            do ibin = 1, maxbin
              ncell_bin_allproc(ibin) = ncell_bin_allthread(ibin)
              natom_bin_allproc(ibin) = ncell_bin_allthread(ibin)
            end do
#endif

    end if

    if (main_rank) then
      if(output%txtfile .ne. '')then

        ! output for proximal distribution function
        !
        if (rdf_option%rmode == RModeProximal) then
          write(txt_out, '(A)') "# This file was created with rdf_analysis"
          write(txt_out, '(A)') "#"
          write(txt_out, '(A)') "# Proximal Distribution Function"
          write(txt_out, '(A)') "# Column 1: distance r(A) "
          write(txt_out, '(A)') "# Column 2: time averaged available volume <V(r, t)> (A^3)"
          write(txt_out, '(A)') "# Column 3: time averaged coordination number of solvent atom <N(r, t)>"
          write(txt_out, '(A)') "# Column 4: number density profile ρ1(r) calculated by <N(r,t)/V(r,t)> (A^-3)"
          write(txt_out, '(A)') "# Column 5: number density profile ρ2(r) calculated by <N(r,t)>/<V(r,t)> (A^-3)"

          if (rdf_option%bulk_value == 0.0_wp) then
            write(txt_out, '(A)') "# Column 6: normalized ρ1(r) with averaged_bulk_1"
            write(txt_out, '(A)') "# Column 7: normalized ρ2(r) with averaged_bulk_2"
            write(txt_out, '(A,f4.1,A)') &
                 "# averaged_bulk_1: bulk density using the last ", &
                 rdf_option%bulk_region, " % of the profile ρ1(r)"
            write(txt_out, '(A,f4.1,A)') &
                 "# averaged_bulk_2: bulk density using the last ", &
                 rdf_option%bulk_region, " % of the profile ρ2(r)"
          else
            write(txt_out, '(A)') "# Column 6: normalized ρ1(r) with bulk"
            write(txt_out, '(A)') "# Column 7: normalized ρ2(r) with bulk"
            write(txt_out, '(A)') "# bulk: bulk density specified by input"
          end if

          write(txt_out, '(A)') " "

          ! calclate bulk density
          !
          if (rdf_option%bulk_value == 0.0_wp) then
            do ibin = 1, maxbin
              if (real(ibin,wp) > real(maxbin,wp)*(100.0_wp - &
                   rdf_option%bulk_region)/100.0_wp)then
                bin_count = bin_count + 1
                density_1 = ave_density_bin(ibin)/(real(totstep,wp))
                density_2 = real(ave_natom_bin(ibin),wp) &
                          /(real(ave_ncell_bin(ibin),wp)*acvolume)
                bulk_density_1 = bulk_density_1 + density_1
                bulk_density_2 = bulk_density_2 + density_2
              end if
            end do
            bulk_density_1 = bulk_density_1 /real(bin_count,wp)
            bulk_density_2 = bulk_density_2 /real(bin_count,wp)
          else
            bulk_density_1 = rdf_option%bulk_value
            bulk_density_2 = rdf_option%bulk_value
          end if

          if (rdf_option%bulk_value == 0.0_wp) then
            write(txt_out, '(A,ES25.16E3)') "averaged_bulk_1: ", bulk_density_1
            write(txt_out, '(A,ES25.16E3)') "averaged_bulk_2: ", bulk_density_2
          else
            write(txt_out, '(A,ES25.16E3)') "bulk: ", bulk_density_1
          end if

          write(txt_out, '(A)') " "

          ! output values
          !
          do ibin = 1, maxbin
            bindist = (ibin-1)*binsize+ binsize/2

            density_1 = ave_density_bin(ibin)/(real(totstep,wp))
            density_2 = real(ave_natom_bin(ibin),wp) &
                      /(real(ave_ncell_bin(ibin),wp)*acvolume)

            norm_density_1 = density_1/bulk_density_1
            norm_density_2 = density_2/bulk_density_2

            avail_volume = real(ave_ncell_bin(ibin),wp)*acvolume &
                         /(real(totstep,wp))
            ncoordi_atom = real(ave_natom_bin(ibin),wp)/(real(totstep,wp))

            write(txt_out, &
                 '(6(ES25.16E3,A),ES25.16E3)')     &
                 bindist, "  ",avail_volume,"    ", ncoordi_atom,"      ", &
                 density_1, "      ", density_2,"      ", &
                 norm_density_1,"      ",norm_density_2
          end do

        ! output for radial distribution function
        !
        else if (rdf_option%rmode == RModeRadial) then

          write(txt_out, '(A)') "# This file was created with rdf_analysis"
          write(txt_out, '(A)') "#"
          write(txt_out, '(A)') "# Radial Distribution Function"
          write(txt_out, '(A)') "# Column 1: distance r(A)"
          write(txt_out, '(A)') "# Column 2: number density profile ρ(r) (A^-3)"

          if(rdf_option%bulk_value ==0.0_wp)then
            write(txt_out, '(A)') "# Column 3: normalized ρ(r) with averaged_bulk"
            write(txt_out, '(A,f4.1,A)') &
                 "# averaged_bulk: bulk density using the last ", rdf_option%bulk_region,&
                 " % of the profile ρ(r)"
          else
            write(txt_out, '(A)') "# Column 3: normalized ρ(r) with bulk"
            write(txt_out, '(A)') "# bulk: bulk density specified by input"
          end if

          ! calclate bulk density
          !
          if (rdf_option%bulk_value == 0.0_wp) then
            do ibin = 1, maxbin
              if (real(ibin,wp) > real(maxbin,wp)*(100.0_wp - &
                           rdf_option%bulk_region)/100.0_wp) then
                bin_count = bin_count + 1

                bindist    = (ibin-1)*binsize+ binsize/2
                bindist_in = bindist - binsize/2
                bindist_out= bindist + binsize/2
                binvolume  = (4.0_wp/3.0_wp)*PI* &
                            (bindist_out*bindist_out*bindist_out - &
                             bindist_in *bindist_in *bindist_in)
                density_2  = real(ave_natom_bin(ibin),wp)/ &
                            (real(totstep,wp)*binvolume)/ &
                             real(size(selec_atom(gid_solute)%idx),wp)

                bulk_density_2 = bulk_density_2 + density_2
              end if
            end do
            bulk_density_2 = bulk_density_2 /real(bin_count,wp)
            write(txt_out, '(A,f8.6)') "averaged_bulk: ", bulk_density_2
          else
            bulk_density_2 = rdf_option%bulk_value
            write(txt_out, '(A,f8.6)') "bulk: ", bulk_density_2
          end if

          do ibin = 1, maxbin
            bindist    = (ibin-1)*binsize+ binsize/2
            bindist_in = bindist - binsize/2
            bindist_out= bindist + binsize/2
            binvolume  = (4.0_wp/3.0_wp)*PI* &
                        (bindist_out*bindist_out*bindist_out - &
                         bindist_in *bindist_in *bindist_in)
            density_2  = real(ave_natom_bin(ibin),wp)/ &
                        (real(totstep,wp)*binvolume)/ &
                         real(size(selec_atom(gid_solute)%idx),wp)
            norm_density_2 = density_2/bulk_density_2
            write(txt_out, '(ES25.16E3,A,ES25.16E3,A,ES25.16E3)') &
                 bindist, "  ", density_2, "  ", norm_density_2
          end do
        end if

      end if
    end if

    if (main_rank) then
      if (output%txtfile .ne. '') then
        call close_file(txt_out)
      end if
    end if

    if (allocated(atoms_focus)) then
      deallocate(atoms_focus, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if
    if (allocated(atoms_solvent)) then
      deallocate(atoms_solvent, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if
    deallocate(ncell_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(ncell_bin_allthread, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(ave_ncell_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(natom_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(natom_bin_allthread, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(ave_natom_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(selec_atom, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    deallocate(density_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(ave_density_bin, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(ncell_bin_allproc, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc
    deallocate(natom_bin_allproc, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    ! close output file
    !

    call timer(TimerAnalysis, TimerOff)

    return

  end subroutine run_rdf

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_options
  !> @brief        check options in control file
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_options(molecule, boundary, domain, trj_list, option,&
                           rdf_option, totstep)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_option),          intent(in)    :: option
    type(s_rdf_option),      intent(in)    :: rdf_option
    integer,                 intent(in)    :: totstep

    ! local variables
    real(wp)                 :: max_csize
    real(wp)                 :: csize_xyz(3)
    real(wp)                 :: min_bsize
    real(wp)                 :: bsize_xyz(3)
    logical                  :: find_error


    csize_xyz(1)   =  boundary%cell_size_x
    csize_xyz(2)   =  boundary%cell_size_y
    csize_xyz(3)   =  boundary%cell_size_z
    max_csize      =  maxval(csize_xyz)
    bsize_xyz(1)   =  boundary%box_size_x
    bsize_xyz(2)   =  boundary%box_size_y
    bsize_xyz(3)   =  boundary%box_size_z
    min_bsize      =  minval(bsize_xyz)

    find_error  = .false.


    ! check options before the analysis loop
    !
    if (totstep == 0) then
      if (option%buffer < rdf_option%ana_range) then
        write(MsgOut,'(a)') &
          'Check_Options> [SPANA_OPTION]:buffer should be larger than [RDF_OPTION]:range'
        find_error  = .true.
      end if

      if (rdf_option%gid_solute == 0) then
        write(MsgOut,'(a)') &
          'Check_Options> [RDF_OPTION]:solute is not specified'
        find_error  = .true.
      end if

      if (rdf_option%gid_solvent == 0) then
        write(MsgOut,'(a)') &
          'Check_Options> [RDF_OPTION]:solvent is not specified'
        find_error  = .true.
      end if

      if (rdf_option%binsize == 0.0_wp) then
        write(MsgOut,'(a)') &
          'Check_Options> [RDF_OPTION]:binsize is not specified'
        find_error  = .true.
      end if

      if (rdf_option%ana_range == 0.0_wp) then
        write(MsgOut,'(a)') &
          'Check_Options> [RDF_OPTION]:range is not specified'
        find_error  = .true.
      end if

      if (boundary%type == BoundaryTypeNOBC) then
        if (option%determine_box == DetermineBoxTrajectory) then
          write(MsgOut,'(a)') &
            'Check_Options> [SPANA_OPTION]:TRAJECTORY is available only in PBC'
          find_error  = .true.
        end if
      else if (boundary%type == BoundaryTypePBC) then
        if (option%determine_box == DetermineboxMax) then
          write(MsgOut,'(a)') &
            'Check_Options> [SPANA_OPTION]:MAX is available only in NOBC'
          find_error  = .true.
        end if
      end if

      if (option%wrap) then
        if (trj_list%trj_type /= 2) then
          write(MsgOut,'(a)') &
            'Check_Options> cannot wrap atoms when [TRAJECTORY]:trj_type is not COOR+BOX'
          find_error  = .true.
        end if
      end if

      if (boundary%num_domain(1) /= 0 .and. &
          boundary%num_domain(2) /= 0 .and. &
          boundary%num_domain(3) /= 0) then

        if (product(boundary%num_domain) /= nproc_country) then
          write(MsgOut,'(a)') &
            'Check_Options> # of process is not domain_x * domain_y * domain_z '
          find_error  = .true.
        end if
      end if
    end if ! totstep ==0

    ! check buffer length no to include self image
    ! this is better to repeat for every analysis frame in the case of NPT
    !
    if (option%buffer+max_csize >= min_bsize/2.0_wp) then
      write(MsgOut,'(a)') &
        'Check_Options> [SPANA_OPTION]:buffer is too large. There is a risk to include the self image atom'
      find_error  = .true.
    end if

    if (find_error) then
      write(MsgOut,'(a, I8)')'Check_Options> error detected at step ',totstep
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
  !! @param[in]    rdf_option : rdf_option information
  !! @param[in]    option     : option information
  !! @param[in]    molecule   : molecule information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] domain     : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine resetup_boundary(rdf_option, option, molecule, boundary, domain)

    ! formal arguments
    type(s_rdf_option),      intent(in)    :: rdf_option
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain

    ! local variables
    integer                  :: num_cells_x
    integer                  :: num_cells_y
    integer                  :: num_cells_z
    integer                  :: remainder


    if (rdf_option%voxel_size == 0) &
      return

    num_cells_x = int(boundary%box_size_x / rdf_option%voxel_size)
    if (mod(num_cells_x, boundary%num_domain(1)) /= 0) then
      remainder = mod(num_cells_x, boundary%num_domain(1))
      num_cells_x = num_cells_x + boundary%num_domain(1) - remainder
    end if

    num_cells_y = int(boundary%box_size_y / rdf_option%voxel_size)
    if (mod(num_cells_y, boundary%num_domain(2)) /= 0) then
      remainder = mod(num_cells_y, boundary%num_domain(2))
      num_cells_y = num_cells_y + boundary%num_domain(2) - remainder
    end if

    num_cells_z = int(boundary%box_size_z / rdf_option%voxel_size)
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

end module rdf_analyze_mod
