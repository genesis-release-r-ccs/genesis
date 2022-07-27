!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sasa_analyze_mod
!> @brief   run sasa analysis
!! @authors Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sasa_analyze_mod

  use sasa_option_str_mod
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

  ! structure
  type, public :: s_arc
    real(wp)   :: ti
    real(wp)   :: te
  end type s_arc

  ! subroutines
  public  :: run_sasa
  private :: set_radius
  private :: calc_atomic_sasa
  private :: get_zs
  private :: get_aa
  private :: find_overlap
  private :: combine_overlap
  private :: sort_arcs
  private :: check_options

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_sasa
  !> @brief        run analyzing trajectories
  !! @authors      IY
  !! @param[inout] molecule    : molecule information
  !! @param[in]    trj_list    : trajectory file list information
  !! @param[in]    output      : output information
  !! @param[inout] option      : option information
  !! @param[in]    sasa_option : contact_option information
  !! @param[inout] trajectory  : trajectory information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] domain      : domain information
  !! @param[in]    ensemble    : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_sasa(molecule, trj_list, output, option, &
                      sasa_option, trajectory, boundary, domain, ensemble)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option), target,  intent(inout) :: option
    type(s_sasa_option),     intent(in)    :: sasa_option
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble

    ! local variables
    type(s_trj_file)         :: trj_in
    real(wp)                 :: cofm(3)
    real(wp)                 :: move(3)
    real(wp)                 :: csize_x, csize_y, csize_z, cvolume, acvolume
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: total_sasa
    integer                  :: nstru, ifile, istep, totstep, nfile
    integer                  :: iatom
    integer                  :: txt_out
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: icell, ncell
    integer                  :: ig, is
    integer                  :: i
    integer                  :: natom_selec
    integer                  :: igroup, ngroup
    integer                  :: gid_solute
    integer                  :: offset
    integer                  :: id, my_id

    type(s_parray),          allocatable   :: selec_atom(:)
    real(wp),                allocatable   :: selec_atom_radi(:)
    real(wp),                allocatable   :: selec_atom_sasa(:)
    real(wp),                allocatable   :: selec_atom_sasa_allrank(:)
    real(wp),                allocatable   :: ave_selec_atom_sasa(:)


    nstru        = 0
    nfile        = size(trj_list%md_steps)
    ngroup       = size(option%analysis_atoms)
    gid_solute   = sasa_option%gid_solute

    alloc_stat   = 0
    dealloc_stat = 0

    acvolume = 0
    totstep  = 0
    txt_out  = 24

    allocate(selec_atom(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! open output file
    !
    if (main_rank)then
      if (output%txtfile /= '') then
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

    allocate(selec_atom_radi(natom_selec), stat = alloc_stat)
    allocate(selec_atom_sasa(natom_selec), stat = alloc_stat)
    allocate(selec_atom_sasa_allrank(natom_selec), stat = alloc_stat)
    allocate(ave_selec_atom_sasa(natom_selec), stat = alloc_stat)

    ave_selec_atom_sasa(1:natom_selec) = 0.0_wp

    domain%num_atom_selec = natom_selec
    domain%num_group      = ngroup

    call set_radius(molecule, option, sasa_option, &
                    selec_atom, selec_atom_radi)

    ! check cntrl parameter
    !
    if (main_rank) then
      call check_options(molecule, boundary, domain, trj_list, &
                         option, sasa_option, totstep, selec_atom_radi)
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

          ! get volume of a cell
          bsize_x = boundary%box_size_x
          bsize_y = boundary%box_size_y
          bsize_z = boundary%box_size_z

          csize_x = bsize_x/boundary%num_cells_x
          csize_y = bsize_y/boundary%num_cells_y
          csize_z = bsize_z/boundary%num_cells_z

          cvolume = csize_x*csize_y*csize_z
          acvolume = acvolume + cvolume

          ! refine atomcoodinate
          do iatom = 1, molecule%num_atoms
            molecule%atom_coord(1,iatom) = trajectory%coord(1,iatom)
            molecule%atom_coord(2,iatom) = trajectory%coord(2,iatom)
            molecule%atom_coord(3,iatom) = trajectory%coord(3,iatom)
          end do

          if (sasa_option%recenter /= 0) then
            call get_cofm(molecule, selec_atom, sasa_option%recenter, cofm)
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

          ! check cntrl parameter for NPT trajectory
          !
          if (main_rank)then
            if (boundary%type == BoundaryTypePBC) then
              if (option%determine_box == DetermineBoxTrajectory) then
                if (ensemble%ensemble == EnsembleNPT  .or. &
                    ensemble%ensemble == EnsembleNPAT .or. &
                    ensemble%ensemble == EnsembleNPgT ) then
                  if (main_rank) then
                    call check_options(molecule, boundary, domain, trj_list, &
                         option, sasa_option, totstep, selec_atom_radi)
                  end if
                end if
              end if
            end if
          end if

          ncell = domain%num_cell_local

          ! determine the first index of the solute
          offset = 0
          if (gid_solute>1) then
            do igroup = 1, gid_solute-1
              offset = offset + size(selec_atom(igroup)%idx)
            end do
          end if

          selec_atom_sasa(1:natom_selec) = 0.0_wp
          selec_atom_sasa_allrank(1:natom_selec) = 0.0_wp

          !$omp parallel default(shared)  &
          !$omp private(id, icell, my_id)

#ifdef OMP
          id = omp_get_thread_num()
#else
          id = 0
#endif
          my_id = id

          do icell = my_id+1, ncell, nthread

            call calc_atomic_sasa(molecule, boundary, domain, option,  &
                                  sasa_option, selec_atom_radi, icell, &
                                  selec_atom_sasa)

          end do ! icell

          !$omp end parallel

#ifdef HAVE_MPI_GENESIS
          call mpi_reduce(selec_atom_sasa, selec_atom_sasa_allrank, &
                          natom_selec, mpi_real8, &
                          mpi_sum, 0, mpi_comm_country, ierror)
#else
          do i = 1, size(selec_atom_sasa_allrank)          
            selec_atom_sasa_allrank(i) = selec_atom_sasa(i)
          end do
#endif

          if (main_rank) then
            if (output%txtfile .ne. '') then
              if (totstep == 1) then

                write(txt_out,'(A)') &
             "# This file was created with sasa_analysis"
                write(txt_out,'(A)') &
             "#"
                write(txt_out,'(A)') &
             "# Solvent Accessible Surface Area (SASA)(A^2)"

                if (sasa_option%out_style == OutStyleHistory .or. &
                    sasa_option%out_style == OutStyleAtomicHistory) then

                  if(sasa_option%out_style == OutStyleAtomicHistory) then

                    write(txt_out,'(A)') &
             "# Column 1: number of analysis frame"
                    write(txt_out,'(A)') &
             "# Column 2: serial index of atom i in the selected atoms"
                    write(txt_out,'(A)') &
             "# Column 3: global index of atom i in the input file (e.g.,psf)"
                    write(txt_out,'(A)') &
             "# Column 4: SASA of the atom i at analysis frame t: SASA(t,i)"

                  end if

                  write(txt_out,'(A)') &
             "# SASA(t): total SASA for each analysis frame t"

                end if

                write(txt_out,'(A)') &
             ""
              end if

              if (gid_solute == 0) then
                total_sasa = 0
                do i = 1, size(selec_atom_sasa_allrank)
                  ig = domain%id_selec2global(i)

                  if (sasa_option%out_style == OutstyleAtomicHistory) then

                    write(txt_out,'(I8,A,I8,A,I8,A,ES25.16E3)') &
                         totstep, " ", &
                         i,       " ", &
                         ig,      " ", &
                         selec_atom_sasa_allrank(i)

                  else if (sasa_option%out_style == OutstyleHistory) then

                  end if
                  ave_selec_atom_sasa(i) = ave_selec_atom_sasa(i) + &
                                           selec_atom_sasa_allrank(i)
                  total_sasa = total_sasa + selec_atom_sasa_allrank(i)
                end do

                if (sasa_option%out_style == OutStyleHistory .or. &
                    sasa_option%out_style == OutStyleAtomicHistory) then

                  write(txt_out,'(A,I8,A,ES25.16E3)') &
                       "SASA(", totstep, ")= ", total_sasa
                end if

              else
                total_sasa = 0
                do i = 1, size(selec_atom(gid_solute)%idx)
                  is = offset + i
                  ig = domain%id_selec2global(i)

                  if (sasa_option%out_style == OutStyleAtomicHistory) then
                    write(txt_out,'(I8,A,I8,A,I8,A,ES25.16E3)') &
                         totstep, " ", &
                         i,       " ", &
                         ig,      " ", &
                         selec_atom_sasa_allrank(is)

                  else if (sasa_option%out_style == OutstyleHistory) then

                  end if
                  ave_selec_atom_sasa(is) = ave_selec_atom_sasa(is) + &
                                            selec_atom_sasa_allrank(is)
                  total_sasa = total_sasa + selec_atom_sasa_allrank(is)
                end do

                if (sasa_option%out_style == OutStyleHistory .or. &
                    sasa_option%out_style == OutStyleAtomicHistory) then

                  write(txt_out,'(A,I8,A,ES25.16E3)') &
                       "SASA(", totstep, ")= ", total_sasa

                end if

              end if ! gid_solute ==0
            end if
          end if ! main_rank

        end if !if ana_period
      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    !output time-average values
    !
    if (main_rank) then
      if (output%txtfile .ne. '') then
        write(txt_out,'(A)') "# Time averaged SASA"

        if (sasa_option%out_style == OutStyleAtomicHistory .or. &
            sasa_option%out_style == OutStyleAtomic) then

          write(txt_out,'(A)') &
               "# Column 1: serial index of atom i in the selected atoms"
          write(txt_out,'(A)') &
               "# Column 2: global index of atom i in the input file (e.g., psf)"
          write(txt_out,'(A)') &
               "# Column 3: time averaged SASA of atom i: <SASA(i,t)>"
        end if

        write(txt_out,'(A)') &
               "# <SASA>: time averaged total SASA"
        write(txt_out,'(A)') &
               ""

        if (gid_solute == 0) then

          total_sasa = 0

          do i = 1, size(selec_atom_sasa_allrank)
            ig = domain%id_selec2global(i)
            ave_selec_atom_sasa(i) = ave_selec_atom_sasa(i)/real(totstep,wp)

            if (sasa_option%out_style == OutStyleAtomicHistory .or. &
                sasa_option%out_style == OutStyleAtomic) then
              write(txt_out,'(I8,A,I8,A,ES25.16E3)') &
                   i," ",ig," ", ave_selec_atom_sasa(i)
            end if

            total_sasa = total_sasa + ave_selec_atom_sasa(i)

          end do

          write(txt_out,'(A, f12.4)') "<SASA>= ", total_sasa

        else

          total_sasa = 0

          do i = 1, size(selec_atom(gid_solute)%idx)
            is = offset + i
            ig = domain%id_selec2global(i)
            ave_selec_atom_sasa(is) = ave_selec_atom_sasa(is)/real(totstep,wp)

             if (sasa_option%out_style == OutStyleAtomicHistory .or. &
                 sasa_option%out_style == OutStyleAtomic) then
              write(txt_out,'(I8,A,I8,A,ES25.16E3)') &
                   i," ",ig," ", ave_selec_atom_sasa(is)
            end if

            total_sasa = total_sasa + ave_selec_atom_sasa(is)

          end do

          write(txt_out,'(A, f12.4)') "<SASA>= ", total_sasa

        end if ! gid_solute ==0
      end if
    end if ! main_rank

    if (main_rank)then
      if (output%txtfile /= '') then
        call close_file(txt_out)
      end if
    end if

    deallocate(selec_atom, stat = dealloc_stat)
    deallocate(selec_atom_radi, stat = dealloc_stat)
    deallocate(selec_atom_sasa, stat = dealloc_stat)
    deallocate(selec_atom_sasa_allrank, stat = dealloc_stat)
    deallocate(ave_selec_atom_sasa, stat = dealloc_stat)

    if (dealloc_stat /= 0) call error_msg_dealloc

    ! close output file
    !
    call timer(TimerAnalysis, TimerOff)

    return

  end subroutine run_sasa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_radius
  !> @brief        set atomic radius for selected atoms
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_radius(molecule, option, sasa_option, selec_atom, &
                        selec_atom_radi)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_option),          intent(in)    :: option
    type(s_sasa_option),     intent(in)    :: sasa_option
    type(s_parray),          allocatable, intent(in)      :: selec_atom(:)
    real(wp),                allocatable, intent(inout)   :: selec_atom_radi(:)

    ! local variables
    real(wp)                 :: radius
    integer                  :: nline, iline, istep
    integer                  :: alloc_stat, dealloc_stat
    integer                  :: file
    integer                  :: igroup, iatom, serial, ngroup
    integer                  :: id_global
    character(10)            :: atom_cls_name
    logical                  :: match

    real(wp),                allocatable   :: atom_radi_list(:)
    character(10),           allocatable   :: atom_cls_name_list(:)


    alloc_stat   = 0
    dealloc_stat = 0

    ! input radias data
    !
    do istep = 1, 2
      call open_file(file, sasa_option%radi_file, IOFileInput)

      nline = 0
      do
        read(file, *, end=999 ) atom_cls_name, radius
        !write(MsgOut, *) atom_cls_name, radius
        nline = nline + 1

        if (istep == 2) then
          atom_cls_name_list(nline) = atom_cls_name
          atom_radi_list(nline)     = radius
        end if

      end do

      call close_file(file)

999   if (istep == 1) then
        allocate(atom_cls_name_list(nline), stat = alloc_stat)
        allocate(atom_radi_list(nline), stat = alloc_stat)
        call close_file(file)
      end if

    end do ! istep

    ! setup atomic radius for selec_atom
    !
    ngroup = size(option%analysis_atoms)
    serial = 0

    do igroup = 1, ngroup
      do iatom = 1, size(selec_atom(igroup)%idx)
        match     = .false.
        serial    = serial + 1
        id_global = selec_atom(igroup)%idx(iatom)
        atom_cls_name = molecule%atom_cls_name(id_global)

        ! serch the name from the list
        !
        do iline = 1, nline
          if(atom_cls_name == trim(atom_cls_name_list(iline))) then
            selec_atom_radi(serial) = atom_radi_list(iline)
            match = .true.
          end if
        end do

        if (.not. match) then
          selec_atom_radi(serial) = 1.75_wp
          write(Msgout,*) "WARNING: radius for atom class ", atom_cls_name, &
               " is missing in the list. r= 1.75A was assigned"
        end if

      end do ! iatom
    end do ! igroup

    deallocate(atom_cls_name_list, stat = dealloc_stat)
    deallocate(atom_radi_list, stat = dealloc_stat)

    return

  end subroutine set_radius

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_atomic_sasa
  !> @brief        calculate atomic sasa
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_atomic_sasa(molecule, boundary, domain, option, sasa_option, &
                              selec_atom_radi, icell, selec_atom_sasa)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_option),          intent(in)    :: option
    type(s_sasa_option),     intent(in)    :: sasa_option
    real(wp),                intent(in)    :: selec_atom_radi(:)
    integer,                 intent(in)    :: icell
    real(wp), allocatable,   intent(inout) :: selec_atom_sasa(:)

    ! local variables
    real(wp)                 :: wat_radi, aw_radi, z_radi, dz
    real(wp)                 :: azp, czp
    real(wp)                 :: aa, total_aa
    integer                  :: icelln
    integer                  :: natom_celln, natom_all
    integer                  :: natom_center
    integer                  :: i, j, k, ig, is, iatom
    integer                  :: icell_count, iatom_count
    integer                  :: istep, iz
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: i_x, i_y, i_z, n_x, n_y, n_z
    integer                  :: gid_solute, igroup, ngroup
    integer                  :: all_celln(27)
    integer                  :: natom_all_celln(27)

    ! selected atom index for  atoms in all 27 cells
    real(wp),   allocatable  :: zs(:)
    integer,    allocatable  :: is_all(:)


    alloc_stat   = 0
    dealloc_stat = 0

    wat_radi     = sasa_option%probe_radius
    dz           = sasa_option%delta_z

    all_celln(1:27)=0

    gid_solute   = sasa_option%gid_solute
    ngroup       = size(option%analysis_atoms)

    ! determin the neighbour cell around icelg
    ! i_x,y,z = global integer position of icell
    !
    i_x = domain%cell_l2gx_orig(icell)
    i_y = domain%cell_l2gy_orig(icell)
    i_z = domain%cell_l2gz_orig(icell)

    ! n_x,y,z = neighbouring position around icell
    !
    do istep = 1, 2
      icell_count = 0
      iatom_count = 0
      natom_all   = 0
      natom_all_celln(1:27) = 0

      do i = -1, 1
        n_x = i_x + i

        if (boundary%num_domain(1) == 1 .and. n_x == 0) then
          n_x =  boundary%num_cells_x 
        end if

        if (boundary%num_domain(1) == 1 .and. & 
          n_x == boundary%num_cells_x + 1) then
          n_x = 1 
        end if

        do j = -1, 1
          n_y = i_y + j

          if (boundary%num_domain(2) == 1 .and. n_y == 0) then
            n_y = boundary%num_cells_y
          end if

          if (boundary%num_domain(2) == 1 .and. &
            n_y == boundary%num_cells_y + 1) then
            n_y = 1
          end if

          do k = -1, 1
            n_z = i_z + k

            if (boundary%num_domain(3) == 1 .and. n_z == 0) then
              n_z = boundary%num_cells_z
            end if

            if (boundary%num_domain(3) == 1 .and. &
              n_z == boundary%num_cells_z + 1) then
              n_z = 1
            end if

            icell_count = icell_count + 1

            ! should not use wraped position for cell_gxyz2l
            !
            icelln =  domain%cell_gxyz2l(n_x,n_y,n_z)
            all_celln(icell_count) = icelln

            do igroup = 1, ngroup
              natom_celln = 0
              natom_celln = domain%num_atom_group(igroup, icelln)
              natom_all_celln(icell_count) = &
                            natom_all_celln(icell_count) + natom_celln
              natom_all = natom_all + natom_celln
            end do ! igroup

            if (istep == 2) then
              do iatom = 1, natom_all_celln(icell_count)
                iatom_count = iatom_count+1
                is = domain%id_l2s(iatom, icelln)
                is_all(iatom_count) = is

              end do
            end if

          end do
        end do
      end do

      if (istep == 1) then
        deallocate(is_all, stat = dealloc_stat)
        allocate(is_all(natom_all), stat = alloc_stat)
      end if

    end do !istep

    natom_center = natom_all_celln(14)

    ! calculation of atomic sasa for atoms in central
    ! cell starts
    !
    do iatom = 1, natom_center

      ! is: one of the atoms in cental cell
      !
      is = domain%id_l2s(iatom, icell)
      ig = domain%id_selec2global(is)
      aw_radi = selec_atom_radi(is) + wat_radi

      call get_zs(molecule, boundary, domain, sasa_option, &
                  selec_atom_radi, is, zs)

      total_aa = 0
      do iz = 1 , size(zs)

        czp = zs(iz)

        azp    = molecule%atom_coord(3,ig) !z-position of the iatom
        z_radi = aw_radi*aw_radi - (azp-czp)*(azp-czp)
        z_radi = sqrt(z_radi)

        call get_aa(molecule, boundary, domain, sasa_option, &
             is, czp, z_radi, is_all, selec_atom_radi, aa)

        total_aa = total_aa + aa

      end do ! iz

      selec_atom_sasa(is)= aw_radi*dz* total_aa

      if (allocated(zs)) then
        deallocate(zs, stat = dealloc_stat)
      end if

    end do ! iatom

    deallocate(is_all, stat = dealloc_stat)

    return

  end subroutine calc_atomic_sasa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_zs
  !> @brief        get the z-coordinates for zs
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_zs(molecule, boundary, domain, sasa_option, &
                    selec_atom_radi,  is, zs)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_sasa_option),     intent(in)    :: sasa_option
    real(wp),                intent(in)    :: selec_atom_radi(:)
    integer,                 intent(in)    :: is
    real(wp), allocatable,   intent(inout) :: zs(:)

    ! local variables
    real(wp)                 :: z, dz
    real(wp)                 :: wat_radi, atom_radi
    real(wp)                 :: top, bottom
    integer                  :: ig, istep, icount
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! should be get from sasa_option
    !
    dz        = sasa_option%delta_z
    wat_radi  = sasa_option%probe_radius

    atom_radi = selec_atom_radi(is)
    ig        = domain%id_selec2global(is)

    z       = molecule%atom_coord(3,ig)

    do istep = 1, 2

      bottom = z - (atom_radi + wat_radi) + 0.5_wp*dz
      top    = z + (atom_radi + wat_radi)

      icount = 0
      do while (bottom < top)
        icount     = icount +1

        if (istep == 2) then
          zs(icount) = bottom
        end if

        bottom     = bottom + dz
      end do

      if (istep == 1) then
        deallocate(zs, stat = dealloc_stat)
        allocate(zs(icount), stat = alloc_stat)
      end if

    end do ! istep

    return

  end subroutine get_zs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_aa
  !> @brief        calculate total length of accesible arcs
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_aa(molecule, boundary,  domain, sasa_option,  &
                         is, czp, z_radi, is_all, selec_atom_radi, aa)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_sasa_option),     intent(in)    :: sasa_option
    integer,                 intent(in)    :: is
    real(wp),                intent(in)    :: czp
    real(wp),                intent(in)    :: z_radi
    integer,                 intent(in)    :: is_all(:)
    real(wp),                intent(in)    :: selec_atom_radi(:)
    real(wp),                intent(inout) :: aa

    ! local variables
    real(wp)                 :: wat_radi, jaw_radi, iaw_radi
    real(wp)                 :: ix, iy, iz, jx, jy, jz
    real(wp)                 :: jz_radi, iz_radi
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: tti, tte
    integer                  :: j, js, jg, ig
    integer                  :: istep, i
    integer                  :: icount
    integer                  :: alloc_stat, dealloc_stat

    real(wp), allocatable    :: barcs_ti(:)
    real(wp), allocatable    :: barcs_te(:)
    real(wp), allocatable    :: barcs2_ti(:)
    real(wp), allocatable    :: barcs2_te(:)


    alloc_stat   = 0
    dealloc_stat = 0

    bsize_x  = boundary%box_size_x
    bsize_y  = boundary%box_size_y
    bsize_z  = boundary%box_size_z

    iz_radi  = z_radi
    wat_radi = sasa_option%probe_radius
    iaw_radi = wat_radi + selec_atom_radi(is)
    ig       = domain%id_selec2global(is)
    ix       = molecule%atom_coord(1,ig)
    iy       = molecule%atom_coord(2,ig)
    iz       = molecule%atom_coord(3,ig)

    do istep = 1, 2
      icount = 0

      do j = 1, size(is_all)
        js = is_all(j) ! id in selec_atom

        if (is /= js) then ! eliminate itself
          jaw_radi = wat_radi + selec_atom_radi(js)
          jg       = domain%id_selec2global(js)
          jx       = molecule%atom_coord(1,jg)
          jy       = molecule%atom_coord(2,jg)
          jz       = molecule%atom_coord(3,jg)

          ! atomic position should be translate
          !
          jx       = jx - bsize_x*anint((jx-ix)/bsize_x)
          jy       = jy - bsize_y*anint((jy-iy)/bsize_y)
          jz       = jz - bsize_z*anint((jz-iz)/bsize_z)

          jz_radi  = jaw_radi*jaw_radi -(jz-czp)*(jz-czp)

          if (jz_radi > 0.0_wp) then
            icount = icount +1
            jz_radi = sqrt(jz_radi)

            if (istep == 2) then
              call find_overlap(molecule, domain, selec_atom_radi, &
                                ix, iy, iz_radi, jx, jy, jz_radi,  &
                                tti, tte)

              barcs_ti(icount) = tti
              barcs_te(icount) = tte

            end if !step=2
          end if

        end if ! eliminate itself
      end do ! j

      if (istep == 1) then

        if (allocated(barcs_ti)) then
          deallocate(barcs_ti, stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
        allocate(barcs_ti(icount), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

        if (allocated(barcs_te)) then
          deallocate(barcs_te, stat = dealloc_stat)
          if (dealloc_stat /= 0) call error_msg_dealloc
        end if
        allocate(barcs_te(icount), stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc

      end if

    end do ! istep

    call combine_overlap(barcs_ti, barcs_te, barcs2_ti, barcs2_te)

    if (size(barcs2_ti) == 0) then
      aa = 2.0_wp * pi
    else
      aa = 0.0_wp
      aa = aa + barcs2_ti(1)

      do i = 2, size(barcs2_ti)
        aa = aa + (barcs2_ti(i) - barcs2_te(i-1))
      end do

      aa = aa + (2.0_wp*pi - barcs2_te(size(barcs2_te)))
    end if

    if (aa > 2.0_wp*pi) aa = 2.0_wp * pi;

    deallocate(barcs_ti, stat = dealloc_stat)
    deallocate(barcs_te, stat = dealloc_stat)
    deallocate(barcs2_ti, stat = dealloc_stat)
    deallocate(barcs2_te, stat = dealloc_stat)

    return

  end subroutine get_aa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    find_overlap
  !> @brief        find overlap for arcs
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine find_overlap(molecule, domain, selec_atom_radi, ix, iy, iz_radi, &
                          jx, jy, jz_radi, tti, tte)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),          intent(in)    :: domain
    real(wp),                intent(in)    :: selec_atom_radi(:)
    real(wp),                intent(in)    :: ix, iy, iz_radi
    real(wp),                intent(in)    :: jx, jy, jz_radi
    real(wp),                intent(inout) :: tti
    real(wp),                intent(inout) :: tte

    ! local variables
    real(wp)                 :: dij2, dij
    real(wp)                 :: dij_x, dij_y
    real(wp)                 :: mt
    real(wp)                 :: v1, v2, v3, v4, v5
    real(wp)                 :: t1, t2, t3, t4
    real(wp)                 :: mx, my, nx, ny, dm, dn


    dij2 = (ix-jx)*(ix-jx)+(iy-jy)*(iy-jy)
    dij  = sqrt(dij2)

    tti = -999.99_wp
    tte = -999.99_wp

    if ((iz_radi+jz_radi <= dij) .or. (dij + jz_radi <= iz_radi)) then
      tti = -999.99_wp
      tte = -999.99_wp
    else if (dij + iz_radi <= jz_radi) then
      tti = 0.0_wp
      tte = 2.0_wp*pi
    else

      dij_x    = jx - ix
      dij_y    = jy - iy
      mt       = acos(dij_x/dij)
      if (dij_y <0.0_wp) mt = 2.0_wp*pi -mt
      v1       = iz_radi*iz_radi + dij2 -jz_radi*jz_radi
      v2       = 2.0_wp*iz_radi*dij
      v3       = v1/v2
      if (v3 > 1.0_wp) v3 = 1.0_wp
      v5       = acos(v3)
      t1       = mt + v5
      if (t1 > 2.0_wp*pi) t1 = t1 - 2.0_wp*pi
      t2       = mt - v5
      if (t2 < 0.0_wp) t2   = t2 + 2.0_wp*pi;

      if (t1 > t2) then
        v4 = t1
        t1 = t2
        t2 = v4
      end if

      t3  = (t1 +t2)/2.0_wp
      t4  = t3 - pi
      if (t4 < 0.0_wp) t4 = t4 + 2.0_wp*pi
      mx  = iz_radi*cos(t3) + ix
      my  = iz_radi*sin(t3) + iy
      nx  = iz_radi*cos(t4) + ix
      ny  = iz_radi*sin(t4) + iy
      dm  = (mx-jx)*(mx-jx) + (my-jy)*(my-jy);
      dn  = (nx-jx)*(nx-jx) + (ny-jy)*(ny-jy);

      if (dm <= jz_radi*jz_radi) then
        tti = t1
        tte = t2
      else if (dn <= jz_radi*jz_radi) then
        tti = t2
        tte = t1
      else

        if (dm < dn) then
          tti = t1
          tte = t2
        else
          tti = t2
          tte = t1
        end if

      end if

    end if

    return

  end subroutine find_overlap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    combine_overlap
  !> @brief        combine the overlapped arcs
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine combine_overlap(barcs_ti, barcs_te, barcs2_ti, barcs2_te)

    ! formal arguments
    real(wp),                intent(in)    :: barcs_ti(:)
    real(wp),                intent(in)    :: barcs_te(:)
    real(wp), allocatable,   intent(inout) :: barcs2_ti(:)
    real(wp), allocatable,   intent(inout) :: barcs2_te(:)

    ! local variables
    integer                  :: i, j, istep
    integer                  :: icount
    integer                  :: alloc_stat
    integer                  :: dealloc_stat
    integer                  :: max_arcs
    logical                  :: included
    logical                  :: overlapping

    real(wp), allocatable    :: arcs_temp_ti(:)
    real(wp), allocatable    :: arcs_temp_te(:)
    real(wp), allocatable    :: arcs_merged_ti(:)
    real(wp), allocatable    :: arcs_merged_te(:)


    max_arcs     = 100
    alloc_stat   = 0
    dealloc_stat = 0

    !
    do istep = 1, 2
      icount = 0

      do i = 1, size(barcs_ti)
        if (barcs_ti(i) >= 0.0_wp) then
          if (barcs_ti(i) < barcs_te(i)) then
            icount = icount + 1
            if (istep == 2) then
              barcs2_ti(icount) = barcs_ti(i)
              barcs2_te(icount) = barcs_te(i)
            end if
          else
            icount = icount + 1
            if (istep == 2) then
              barcs2_ti(icount) = 0.0_wp
              barcs2_te(icount) = barcs_te(i)
            end if
            icount = icount +1
            if (istep == 2) then
              barcs2_ti(icount) = barcs_ti(i)
              barcs2_te(icount) = 2.0_wp* pi
            end if
          end if ! ti < te
        end if  ! ti .ge. 0.0
      end do ! i

      if (istep == 1) then
        deallocate(barcs2_ti, stat = dealloc_stat)
        allocate(barcs2_ti(icount), stat = alloc_stat)
        deallocate(barcs2_te, stat = dealloc_stat)
        allocate(barcs2_te(icount), stat = alloc_stat)
      end if

    end do ! istep

    !
    do istep = 1, 2
      icount = 0

      do i = 1, size(barcs2_ti)
        included = .false.
        do j = 1, size(barcs2_ti)

          if (i /= j) then

            if ((barcs2_ti(j) <= barcs2_ti(i)) .and. &
                (barcs2_te(i) <= barcs2_te(j))) then

              if((barcs2_ti(j) /= barcs2_ti(i)) .or. &
                 (barcs2_te(i) /= barcs2_te(j))) then
                included = .true.
                exit
              end if

            end if

          end if !i .ne. j

        end do !j

        if (.not. included) then
          icount = icount +1

          if (istep == 2) then
            arcs_temp_ti(icount) = barcs2_ti(i)
            arcs_temp_te(icount) = barcs2_te(i)
          end if

        end if !included = .false.

      end do  !i

      if (istep == 1) then
        deallocate(arcs_temp_ti, stat = dealloc_stat)
        allocate(arcs_temp_ti(icount), stat = alloc_stat)
        deallocate(arcs_temp_te, stat = dealloc_stat)
        allocate(arcs_temp_te(icount), stat = alloc_stat)
      end if

    end do ! istep

    !
    call sort_arcs(arcs_temp_ti, arcs_temp_te)

    deallocate(arcs_merged_ti, stat = dealloc_stat)
    allocate(arcs_merged_ti(max_arcs), stat = alloc_stat)
    deallocate(arcs_merged_te, stat = dealloc_stat)
    allocate(arcs_merged_te(max_arcs), stat = alloc_stat)

    !
    do i = 1, max_arcs
      arcs_merged_ti(i) = -999.99_wp
      arcs_merged_te(i) = -999.99_wp
    end do

    !
    icount = 0 ! number of non-overlapping arcs
    do i = 1, size(arcs_temp_ti)
      overlapping = .false.

      do j = 1, max_arcs
        if (arcs_merged_ti(j)< -999) exit

        if (arcs_temp_ti(i) <= arcs_merged_te(j)) then
          arcs_merged_te(j) = arcs_temp_te(i)
          overlapping = .true.
          exit
        end if

      end do !j

      if (.not. overlapping) then
        icount = icount + 1
        arcs_merged_ti(icount) = arcs_temp_ti(i)
        arcs_merged_te(icount) = arcs_temp_te(i)
      end if

    end do ! i

    !
    deallocate(barcs2_ti, stat = dealloc_stat)
    allocate(barcs2_ti(icount), stat = alloc_stat)
    deallocate(barcs2_te, stat = dealloc_stat)
    allocate(barcs2_te(icount), stat = alloc_stat)

    do i = 1, icount
      barcs2_ti(i) = arcs_merged_ti(i)
      barcs2_te(i) = arcs_merged_te(i)
    end do

    deallocate(arcs_temp_ti,   stat = dealloc_stat)
    deallocate(arcs_merged_ti, stat = dealloc_stat)
    deallocate(arcs_temp_te,   stat = dealloc_stat)
    deallocate(arcs_merged_te, stat = dealloc_stat)

    return

  end subroutine combine_overlap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sort_arcs
  !> @brief        sort arcs
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sort_arcs(arcs_ti, arcs_te)

    ! formal arguments
    real(wp),                intent(inout) :: arcs_ti(:)
    real(wp),                intent(inout) :: arcs_te(:)

    ! local variables
    real(wp)                 :: temp
    integer                  :: i, j
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! perform bubble sort for arcs
    !
    do i = 1, size(arcs_ti)-1
      do j = i, size(arcs_ti)

        if (arcs_ti(i) > arcs_ti(j)) then
          temp = arcs_ti(i)
          arcs_ti(i) = arcs_ti(j)
          arcs_ti(j) = temp

          temp = arcs_te(i)
          arcs_te(i) = arcs_te(j)
          arcs_te(j) = temp
        end if

        if (arcs_ti(i) == arcs_ti(j) .and. &
            arcs_te(i)  > arcs_te(j)) then
          temp = arcs_te(i)
          arcs_te(i) = arcs_te(j)
          arcs_te(j) = temp
        end if

      end do !j
    end do  !i

    return

  end subroutine sort_arcs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_options
  !> @brief        check options in control file
  !! @authors      IY
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_options(molecule, boundary, domain, trj_list, option,&
                           sasa_option, totstep, selec_atom_radi)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_option),          intent(in)    :: option
    type(s_sasa_option),     intent(in)    :: sasa_option
    integer,                 intent(in)    :: totstep
    real(wp), allocatable,   intent(in)    :: selec_atom_radi(:)

    ! local variables
    real(wp)                 :: max_csize
    real(wp)                 :: min_csize
    real(wp)                 :: csize_xyz(3)
    real(wp)                 :: min_bsize
    real(wp)                 :: bsize_xyz(3)
    real(wp)                 :: max_radi
    real(wp)                 :: probe_radi
    logical                  :: find_error


    find_error   = .false.
    probe_radi   =  sasa_option%probe_radius
    csize_xyz(1) =  boundary%cell_size_x
    csize_xyz(2) =  boundary%cell_size_y
    csize_xyz(3) =  boundary%cell_size_z
    max_csize    =  maxval(csize_xyz)
    min_csize    =  minval(csize_xyz)
    bsize_xyz(1) =  boundary%box_size_x
    bsize_xyz(2) =  boundary%box_size_y
    bsize_xyz(3) =  boundary%box_size_z
    min_bsize    =  minval(bsize_xyz)
    max_radi     =  maxval(selec_atom_radi)


    ! check options before the analysis loop
    !
    if (totstep == 0) then
      if (option%buffer < (probe_radi + max_radi)*2.0_wp) then
        write(MsgOut,'(a)') &
          'Check_Options> [SPANA_OPTION]:buffer should be larger than'
        write(MsgOut,'(a)') &
          'Check_Options> (probe_radius + max_atom_radius)*2'
        find_error  = .true.
      end if

      if (boundary%type == BoundaryTypeNOBC) then
        if (option%determine_box == DetermineBoxTrajectory) then
          write(MsgOut,'(a)') &
            'Check_Options> [SPANA_OPTION]:TRAJECTORY is available only in PBC'
          find_error  = .true.
        end if
      else if (boundary%type == BoundaryTypePBC) then
        if (option%determine_box == DetermineBoxMax) then
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

    ! check cell size: cell side should be large enough
    ! to detect overapping circles around the target atom
    !
    if (min_csize < (probe_radi + max_radi)*2.0_wp) then
      write(MsgOut,'(a)') &
        'Check_Options> [BOUNDARY]:num_cells_ is too large '
      write(MsgOut,'(a)') &
        'Check_Options> length of minimum cell side should be larger than '
      write(MsgOut,'(a)') &
        'Check_Options> (probe_radius + max_atom_radius )*2'
      find_error  = .true.
    end if

    ! check buffer length no to include self image
    ! this is better to repeat for every analysis frame in the case of NPT
    !
    if (option%buffer+max_csize >= min_bsize/2.0_wp) then
      write(MsgOut,'(a)') &
        'Check_Options>  [SPANA_OPTION]:buffer is too large. There is a risk to include the self image atom'
      find_error  = .true.
    end if

    if (find_error) then
      write(MsgOut,'(a, I8)')'Check_Options> error detected at step ',totstep
      call error_msg(&
        'Check_Options> analysis terminated because of option error')
    end if

    return

  end subroutine check_options

end module sasa_analyze_mod
