!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   kc_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kc_analyze_mod

  use kc_option_str_mod
  use fitting_mod
  use fileio_trj_mod
  use fitting_str_mod
  use trajectory_str_mod
  use input_str_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use string_mod
  use random_mod
  use select_atoms_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: assign_mass
  private :: get_replicate_name1

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[inout] molecule   : molecule information
  !! @param[in]    input      : input information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, input, trj_list, trajectory, fitting, option, output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_input),           intent(in)    :: input
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_trj_file)         :: trj_in
    type(s_pdb)              :: pdb_out
    real(wp)                 :: rmsd, min_rmsd, convergency, accel, min_convergency
    real(wp)                 :: diff_coord(3)
    integer                  :: i, j, k, idx, nclst, iclst, jclst, iseed
    integer                  :: idx_in, idx_out
    integer                  :: natom, niter, ntraj, nstru
    integer                  :: iatom, iiter, itraj, istep, istru
    integer                  :: rms_out, alloc_stat
    logical                  :: converged
    character(MaxLine)       :: linein

    real(wp),         allocatable :: av0_coord_tmp(:,:), av0_coord_tmp2(:,:)
    real(wp),         allocatable :: ave_coord_tmp(:,:)
    real(wp),         allocatable :: av0_coord(:,:,:)
    real(wp),         allocatable :: ave_coord(:,:,:)
    real(wp),         allocatable :: cnt_coord(:,:,:)
    real(wp),         allocatable :: trj_coord(:,:)
    real(wp),         allocatable :: sqrt_mass(:)
    real(wp),         allocatable :: min_rmsd_clst(:)
    real(wp),         allocatable :: sum_rmsd(:)
    integer,          allocatable :: diff1(:), diff2(:)
    integer,          allocatable :: cluster_index(:)
    integer,          allocatable :: cluster_index_old(:)
    integer,          allocatable :: center_index(:)
    integer,          allocatable :: ndata(:)
    integer,          allocatable :: ndata_old(:)
    logical,          allocatable :: init_cluster(:)
    type(s_trj_file), allocatable :: trj_out(:)

    if (option%check_only) &
      return

    natom   = molecule%num_atoms
    nclst   = option%num_clusters
    iseed   = option%iseed
    ntraj   = size(trj_list%md_steps)

    ! allocate memory
    !
    allocate(sqrt_mass(natom),          &
             av0_coord_tmp(3,natom),    &
             ave_coord_tmp(3,natom),    &
             av0_coord_tmp2(3,natom),   &
             av0_coord (nclst,3,natom), &
             ave_coord (nclst,3,natom), &
             cnt_coord (nclst,3,natom), &
             trj_coord(3,natom),        &
             init_cluster(nclst),       &
             sum_rmsd(nclst),           &
             center_index(nclst),       &
             min_rmsd_clst(nclst),      &
             diff1(nclst),              &
             diff2(nclst),              &
             ndata(nclst), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! setup mass
    !
    if (fitting%mass_weight) then
      call assign_mass(molecule)
    end if

    if (fitting%mass_weight) then
      do iatom = 1, natom
        sqrt_mass(iatom) = sqrt(molecule%mass(iatom))
      end do
    else
      do iatom = 1, natom
        sqrt_mass(iatom) = 1.0_wp
      end do
    end if

    ! count the number of snapshots to be analyzed
    !
    istru = 0
    do itraj = 1, ntraj
      call open_trj(trj_in,                    &
                    trj_list%filenames(itraj), &
                    trj_list%trj_format,       &
                    trj_list%trj_type, IOFileInput)
      do istep = 1, trj_list%md_steps(itraj)
        call read_trj(trj_in, trajectory)
        if (mod(istep, trj_list%ana_periods(itraj)) == 0) then
          istru = istru + 1
        end if
      end do
      call close_trj(trj_in)
    end do
    nstru = istru

    ! setup cluster index
    !
    allocate(cluster_index(nstru),cluster_index_old(nstru))
    cluster_index(1:nstru) = 0

    ! use input cluster index as the initial cluster index
    if (input%indexfile /= '') then
      call open_file(idx_in, input%indexfile, IOFileInput)

      do while (.true.)
        read (idx_in,'(A)',end=10) linein
        read (linein,'(I10)') istru
        if (istru <= nstru) then
          read (linein,'(I10,1x,I10)') istru, cluster_index(istru)
        end if
      end do
10    continue

      call close_file(idx_in)

      do i = 1, nstru
        if (cluster_index(i) <= 0 .or. cluster_index(i) > nclst) then
          iclst = int(random_get_legacy(iseed)*nclst) + 1
          if (iclst <= 0    ) iclst = 1
          if (iclst >  nclst) iclst = nclst
          cluster_index(i) = iclst
        end if
      end do

    else
      do i = 1, nstru
        iclst = int(random_get_legacy(iseed)*nclst) + 1
        if (iclst <= 0    ) iclst = 1
        if (iclst >  nclst) iclst = nclst
        cluster_index(i) = iclst
      end do
    end if
    init_cluster(1:nclst) = .false.

    write(MsgOut,'(A)') 'Analyze> initial cluster index'
    do i = 1, nstru
      write(MsgOut,'(5x,I10,1x,I10)') i, cluster_index(i)
    end do
    write(MsgOut,'(A)') ' '

    ! open output files
    !
    if (output%indexfile /= '') &
      call open_file(idx_out, output%indexfile, IOFileOutputNew)

    if (output%trjfile /= '') then
      allocate(trj_out(nclst))
      do iclst = 1, nclst
        call open_trj(trj_out(iclst),                             &
                      get_replicate_name1(output%trjfile, iclst), &
                      option%trjout_format, &
                      option%trjout_type,   &
                      IOFileOutputNew)
      end do
    end if


    ! analysis loop
    !
    do iiter = 1, option%max_iteration

      write(MsgOut,'(A,i10)') 'Analyze> k-means iteration = ', iiter

      ! store old information
      !
      cluster_index_old(1:nstru) = cluster_index(1:nstru)

      ! iteration to obtain averaged coordinates
      !
      do k = 1, option%num_iterations

        ! reset average coordinates
        !
        do j = 1, nclst
          do iatom = 1, natom
            ave_coord(j,1:3,iatom) = 0.0_wp
          end do
        end do

        istru = 0
        ndata(1:nclst) = 0

        DO itraj = 1, ntraj

          call open_trj(trj_in,                    &
                        trj_list%filenames(itraj), &
                        trj_list%trj_format,       &
                        trj_list%trj_type, IOFileInput)

          do istep = 1, trj_list%md_steps(itraj)

            call read_trj(trj_in, trajectory)

            if (mod(istep, trj_list%ana_periods(itraj)) == 0) then

              istru = istru + 1

              ! Set initial averaged coordinates 
              !
              iclst = cluster_index(istru)

              if (.not. init_cluster(iclst)) then
                do iatom = 1, natom
                  av0_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom)
                  cnt_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom) 
                end do
                center_index(iclst)  = istru
                init_cluster(iclst) = .true.
              end if

              ! fitting (geometrical fitting followed by mass-weighted fitting)
              !
              do iatom = 1, natom
                av0_coord_tmp(1:3,iatom) = av0_coord(iclst,1:3,iatom)
              end do

              call run_fitting(fitting, av0_coord_tmp, trajectory%coord, trajectory%coord)

              do iatom = 1, natom
                av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom) *sqrt_mass(iatom)
                trj_coord(1:3,iatom)      = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
              end do

              call run_fitting(fitting, av0_coord_tmp2, trj_coord, trj_coord)

              ! sum trajectory coordinates
              !
              ndata(iclst) = ndata(iclst) + 1
              do iatom = 1, natom
                ave_coord(iclst,1:3,iatom) = ave_coord(iclst,1:3,iatom) + trj_coord(1:3,iatom)
              end do

            end if

          end do

          call close_trj(trj_in)

        END DO

        ! compute averaged coordinates of each cluster
        !
        do iclst = 1, nclst

          do iatom = 1, natom
            ave_coord(iclst,1:3,iatom) = ave_coord(iclst,1:3,iatom) / real(ndata(iclst), wp)
          end do

          do iatom = 1, natom
            av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom)*sqrt_mass(iatom)
            ave_coord_tmp (1:3,iatom) = ave_coord(iclst,1:3,iatom)
          end do

          call run_fitting(fitting, av0_coord_tmp2, ave_coord_tmp, ave_coord_tmp)

          do iatom = 1, natom
            av0_coord(iclst,1:3,iatom) = ave_coord_tmp(1:3,iatom)
            av0_coord(iclst,1:3,iatom) = av0_coord(iclst,1:3,iatom)/sqrt_mass(iatom)
          end do

        end do

      end do

      ! update cluster index and cluster center coordinates
      !
      istru                  = 0
      sum_rmsd     (1:nclst) = 0.0_wp
      ndata        (1:nclst) = 0
      min_rmsd_clst(1:nclst) = 999999999.9_wp

      DO itraj = 1, ntraj

        call open_trj(trj_in,                    &
                      trj_list%filenames(itraj), &
                      trj_list%trj_format,       &
                      trj_list%trj_type, IOFileInput)

        do istep = 1, trj_list%md_steps(itraj)

          call read_trj(trj_in, trajectory)

          if (mod(istep, trj_list%ana_periods(itraj)) == 0) then

            istru = istru + 1

            ! get my cluster index by comparing RMSD
            !
            min_rmsd = 999999999.99_wp

            do iclst = 1, nclst

              do iatom = 1, natom
                av0_coord_tmp(1:3,iatom) = av0_coord(iclst,1:3,iatom)
              end do
              call run_fitting(fitting, av0_coord_tmp, trajectory%coord, trajectory%coord)

              do iatom = 1, natom
                av0_coord_tmp2(1:3,iatom) = av0_coord(iclst,1:3,iatom) *sqrt_mass(iatom)
                trj_coord(1:3,iatom)      = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
              end do
              call run_fitting(fitting, av0_coord_tmp2, trj_coord, trj_coord)

              ! calculate RMSD for the selected atoms
              !
              rmsd = 0.0_wp
              do iatom = 1, size(option%analysis_atom%idx)
                idx = option%analysis_atom%idx(iatom)
                diff_coord(1:3) = av0_coord_tmp2(1:3, idx) - trj_coord(1:3, idx)
                rmsd = rmsd + dot_product(diff_coord, diff_coord)
              end do

              rmsd = sqrt(rmsd / real(size(option%analysis_atom%idx), wp))

              if (rmsd <= min_rmsd) then
                min_rmsd = rmsd
                cluster_index(istru) = iclst
              end if

            end do

            ! get minimum RMSD and update cluster center information
            !
            if (min_rmsd <= min_rmsd_clst(cluster_index(istru))) then
              min_rmsd_clst(cluster_index(istru)) = min_rmsd
              center_index(cluster_index(istru))  = istru

              iclst = cluster_index(istru)
              do iatom = 1, natom
                cnt_coord(iclst,1:3,iatom) = trajectory%coord(1:3,iatom)
              end do
            end if

            ! accumulate important information
            !
            sum_rmsd(cluster_index(istru)) = sum_rmsd(cluster_index(istru)) + min_rmsd
            ndata   (cluster_index(istru)) = ndata   (cluster_index(istru)) + 1

          end if

        end do

        call close_trj(trj_in)

      END DO

      ! output results
      !
      write(MsgOut,'(A)') '   Cluster index    Cluster center   # of structures    Cluster radius'
      do iclst = 1, nclst
        if (ndata(iclst) /= 0) then
          write(MsgOut,'(I16,I18,I18,F18.5)') iclst, center_index(iclst), ndata(iclst), &
                                              sum_rmsd(iclst) / real(ndata(iclst), wp)
        else
          write(MsgOut,'(I16,I18,I18,F18.5)') iclst, 0, ndata(iclst), 0.0_wp
        end if
      end do
      write(MsgOut,'(A)') ''

      ! check empty cluster
      !
      do iclst = 1, nclst
        if (ndata(iclst) == 0) then
          write(MsgOut,'(A)') ' Warning: Empty cluster was generated'
          write(MsgOut,'(A)') ' Closest structure to the old center was selected as a new center'
          write(MsgOut,'(A)') ''
          do iatom = 1, natom
            av0_coord(iclst,1:3,iatom) = cnt_coord(iclst,1:3,iatom)
          end do
          cluster_index(center_index(iclst)) = iclst
        end if
      end do

      ! check convergence
      !
      if (iiter == 1) then
        diff1(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff1(cluster_index(istru)) = diff1(cluster_index(istru)) + 1
          end if
        end do
        converged = .false.
      else if (iiter == 2) then
        diff2(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff2(cluster_index(istru)) = diff2(cluster_index(istru)) + 1
          end if
        end do
        converged = .false.
      else
        diff1(1:nclst) = diff2(1:nclst)
        diff2(1:nclst) = 0
        do istru = 1, nstru
          if (cluster_index_old(istru) /= cluster_index(istru)) then
            diff2(cluster_index(istru)) = diff2(cluster_index(istru)) + 1
          end if
        end do

        min_convergency = 999999999
        do iclst = 1, nclst
          convergency = 100.0_wp - 100.0_wp * abs(diff2(iclst) - diff1(iclst))/ndata(iclst)
          if (convergency < min_convergency) then
            min_convergency = convergency
          end if
        end do

        converged = .true.
        if (min_convergency < option%stop_threshold) then
          converged = .false.
        end if

        write(MsgOut,'(A,F10.5,A)') '   Convergency = ', min_convergency,' %'
        write(MsgOut,'(A)') ''

      end if

      if (converged) exit

    end do


    ! output index file
    !
    if (output%indexfile /= '') then
      do istru = 1, nstru
        write(idx_out,'(I10,1X,I10)') istru, cluster_index(istru)
      end do
      call close_file(idx_out)
    end if


    ! output PDB file of cluster center and new trajectory files
    !
    if (output%pdbfile /= '' .or. output%trjfile /= '' ) then

      write(MsgOut,'(A)') 'Analyze> output PDB files of the cluster centers'
      write(MsgOut,'(A)') ''

      istru = 0

      DO itraj = 1, ntraj

        call open_trj(trj_in,                    &
                      trj_list%filenames(itraj), &
                      trj_list%trj_format,       &
                      trj_list%trj_type, IOFileInput)

        do istep = 1, trj_list%md_steps(itraj)

          call read_trj(trj_in, trajectory)

          if (mod(istep, trj_list%ana_periods(itraj)) == 0) then

            istru = istru + 1

            do iclst = 1, nclst

              ! output PDB file
              !
              if (output%pdbfile /= '') then
                if (istru == center_index(iclst)) then

                  molecule%atom_coord(1:3,1:natom) = cnt_coord(iclst,1:3,1:natom)

                  call run_fitting(fitting,             &
                                   molecule%atom_coord, &
                                   trajectory%coord,    &
                                   trajectory%coord)

                  molecule%atom_coord(1:3,1:natom) = trajectory%coord(1:3,1:natom)

                  call export_molecules(molecule, option%trjout_atom, pdb_out)
                  write(MsgOut,'(A,I10,A)') '   structure = ',istru, '  >  ' // &
                                             trim(get_replicate_name1(output%pdbfile,iclst))
                  call output_pdb(get_replicate_name1(output%pdbfile,iclst), pdb_out)
                  call dealloc_pdb_all(pdb_out)
                end if
              end if

              ! output TRJ file
              !
              if (output%trjfile /= '') then
                if (iclst == cluster_index(istru)) then
                  molecule%atom_coord(1:3,1:natom) = cnt_coord(iclst,1:3,1:natom)

                  call run_fitting(fitting,             &
                                   molecule%atom_coord, &
                                   trajectory%coord,    &
                                   trajectory%coord)

                  call write_trj(trj_out(iclst), trajectory, option%trjout_atom, molecule)
                end if
              end if

            end do

          end if

        end do

        call close_trj(trj_in)

      END DO


    end if

    write(MsgOut,'(A)') ''


    ! close output file
    !
    call close_file(idx_out)

    if  (output%trjfile /= '') then
      do iclst = 1, nclst
        call close_trj (trj_out(iclst))
      end do
    end if


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [indexfile] ' // trim(output%indexfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Index of the cluster to which the structure belongs'
    write(MsgOut,'(A)') ''


    ! deallocate memory
    !
    deallocate(sqrt_mass, av0_coord_tmp, ave_coord_tmp, av0_coord_tmp2, &
               av0_coord, ave_coord, cnt_coord, trj_coord,              &
               init_cluster, sum_rmsd, center_index, min_rmsd_clst,     &
               diff1, diff2, ndata, cluster_index,cluster_index_old)

    if (output%trjfile /= '') then
      deallocate(trj_out)
    end if


    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_mass
  !> @brief        assign mass
  !! @authors      NT, TM
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_mass(molecule)

    ! parameters
    real(wp),                parameter     :: MassH   =  1.008000_wp
    real(wp),                parameter     :: MassC   = 12.011000_wp
    real(wp),                parameter     :: MassN   = 14.007000_wp
    real(wp),                parameter     :: MassO   = 15.999000_wp
    real(wp),                parameter     :: MassS   = 32.060000_wp
    real(wp),                parameter     :: MassP   = 30.974000_wp
    real(wp),                parameter     :: MassMG  = 24.305000_wp
    real(wp),                parameter     :: MassZN  = 65.370000_wp

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    write(MsgOut,'(a)'),      'WARNING: atom mass is not assigned.'
    write(MsgOut,'(a)'),      '   uses default mass.'
    write(MsgOut,'(a,f9.6)')  '      1) H   : ', MassH
    write(MsgOut,'(a,f9.6)')  '      2) C   : ', MassC
    write(MsgOut,'(a,f9.6)')  '      3) N   : ', MassN
    write(MsgOut,'(a,f9.6)')  '      4) O   : ', MassO
    write(MsgOut,'(a,f9.6)')  '      5) S   : ', MassS
    write(MsgOut,'(a,f9.6)')  '      6) P   : ', MassP
    write(MsgOut,'(a,f9.6)')  '      7) MG  : ', MassMG
    write(MsgOut,'(a,f9.6)')  '      8) ZN  : ', MassZN


    do i = 1, molecule%num_atoms

      if (molecule%atom_name(i)(1:1) == 'H') then
        molecule%mass(i) = MassH
      else if (molecule%atom_name(i)(1:1) == 'C') then
        molecule%mass(i) = MassC
      else if (molecule%atom_name(i)(1:1) == 'N') then
        molecule%mass(i) = MassN
      else if (molecule%atom_name(i)(1:1) == 'O') then
        molecule%mass(i) = MassO
      else if (molecule%atom_name(i)(1:1) == 'S') then
        molecule%mass(i) = MassS
      else if (molecule%atom_name(i)(1:1) == 'P') then
        molecule%mass(i) = MassP
      else if (molecule%atom_name(i)(1:2) == 'MG') then
        molecule%mass(i) = MassMG
      else if (molecule%atom_name(i)(1:2) == 'ZN') then
        molecule%mass(i) = MassZN
      else
        write(MsgOut,'(a,a)') 'Assign_Mass> Unknown atom :', &
             molecule%atom_name(i)
      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_mass

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:)

    return

  end function get_replicate_name1

end module kc_analyze_mod
