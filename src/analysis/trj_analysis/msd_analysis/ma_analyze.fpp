!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_analyze_mod
!> @brief   analyze mean square displacement
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_analyze_mod

  use ma_option_str_mod
  use fileio_trj_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_molecules_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use molecule_manipulate_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: set_analysis_set_molecule_ranges
  private :: get_start_end

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyze mean square displacement
  !! @authors      DS, TM
  !! @param[in]    molecule   : molecule structure
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

#   if defined(_OPENMP) || defined(OMP)
    use omp_lib
#   endif /*_OPENMP || OMP*/

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)                          :: trj_in
    integer                                   :: ifile, istep, num_trjfiles
    integer                                   :: nmol, imol, i, j, iaxis, &
                                                 msd_out, msg_interval, natm, &
                                                 istep_local

    ! number of steps to be analyzed in total
    integer                                   :: md_steps_analysis

    ! sample frames count
    integer, allocatable, dimension(:)        :: nsamples

    ! number of molecules in each analysis set
    integer, allocatable, dimension(:)        :: nmols_aset

    real(wp), dimension(3)                    :: com, com_local
    real(wp), allocatable, dimension(:, :)    :: sd_sums
    real(wp), allocatable, dimension(:, :)    :: coord_prev, msd

    ! array that holds all analysis molecules from all analysis sets
    ! same atoms can apear more than once
    type(s_one_molecule), dimension(:), allocatable :: all_analysis_mols

    ! indexes that show which molecule belongs to which analysis set
    ! in all_analysis_mols
    integer,              dimension(:), allocatable :: imol2iaset


    integer, save                                   :: omp_mol_start, &
                                                       omp_mol_end
    integer, save                                   :: omp_atm_start, &
                                                       omp_atm_end
    integer                                         :: omp_id, omp_size, &
                                                       omp_mol_stride, &
                                                       omp_atm_stride

    ! molecule starts and ends at each analysi set
    integer,  allocatable, dimension(:), save       :: omp_mol_starts
    integer,  allocatable, dimension(:), save       :: omp_mol_ends

    real(wp), allocatable, dimension(:, :, :), save :: sd_sums_local
    real(wp), allocatable, dimension(:, :),    save :: coms, coms_prev, &
                                                       dcoms_sums
    real(wp), allocatable, dimension(:, :, :), save :: dcoms_buf

    ! analysis mol set number
    integer                                         :: iaset

    integer, dimension(2) :: start_end


    !$omp threadprivate(omp_atm_start, omp_atm_end, sd_sums_local)
    !$omp threadprivate(omp_mol_start, omp_mol_end)
    !$omp threadprivate(coms, coms_prev, dcoms_sums, dcoms_buf)
    !$omp threadprivate(omp_mol_starts, omp_mol_ends)

    if (option%check_only) &
      return

    ! count total number of analysis molecules and atoms in them
    !
    allocate(nmols_aset(size(option%analysis_mols)))
    nmols_aset = 0
    nmol = 0
    natm = 0
    do iaset = 1, size(option%analysis_mols)
      nmols_aset(iaset) = size(option%analysis_mols(iaset)%mols)
      nmol = nmol + nmols_aset(iaset)
    end do
    write(MsgOut, '(A, i0)') &
      "Analyze> Total number of analysis molecules in system: ", nmol

    ! assign all analysis molecules to all_analysis_mols
    !
    allocate(all_analysis_mols(nmol), imol2iaset(nmol))
    imol = 1
    do iaset = 1, size(option%analysis_mols)
      all_analysis_mols(imol:imol+nmols_aset(iaset)-1) &
        = option%analysis_mols(iaset)%mols
      imol2iaset(imol:imol+nmols_aset(iaset)-1) = iaset
      imol = imol + nmols_aset(iaset)
    end do
    natm = sum([(all_analysis_mols(i)%num_atoms,i=1,size(all_analysis_mols))])
    write(MsgOut, '(A, i0)') &
      "Analyze> Total number of analysis atoms in system: ", natm

    ! open output file
    !
    if (output%msdfile /= '') &
      call open_file(msd_out, output%msdfile, IOFileOutputNew)

    num_trjfiles = size(trj_list%filenames)

    allocate(sd_sums(size(option%analysis_mols), option%delta), &
             nsamples(option%delta),                            &
             msd(size(option%analysis_mols), option%delta))

    !$omp parallel default(none) &
    !$omp shared(omp_size)
#   if defined(_OPENMP) || defined(OMP)
    !$omp master
    omp_size = omp_get_num_threads()
    !$omp end master
#   else /*!_OPENMP && !OMP*/
    omp_size = 1
#   endif /*_OPENMP || OMP*/
    !$omp end parallel

#   if defined(_OPENMP) || defined(OMP)
    call omp_set_dynamic(.false.)
#   endif /*_OPENMP || OMP*/

    write(MsgOut, '(A, i0)') &
      "Analyze> Number of OpenMP threads: ", omp_size

    !$omp parallel default(none)                                                   &
    !$omp private(omp_id, omp_mol_stride, omp_atm_stride, start_end)               &
    !$omp shared(nmol, omp_size, option, molecule, all_analysis_mols, imol2iaset)

#   if defined(_OPENMP) || defined(OMP)
    omp_id = omp_get_thread_num()
#   else /*!_OPENMP && !OMP*/
    omp_id = 0
#   endif /*_OPENMP || OMP*/


    ! molecules
    start_end     = get_start_end(nmol, omp_size, omp_id)
    omp_mol_start = start_end(1)
    omp_mol_end   = start_end(2)

    ! atoms
    start_end     = get_start_end(molecule%num_atoms, omp_size, omp_id)
    omp_atm_start = start_end(1)
    omp_atm_end   = start_end(2)

    ! molecules for each analysis set
    call set_analysis_set_molecule_ranges(all_analysis_mols, &
      imol2iaset, omp_size, omp_id, omp_mol_starts, omp_mol_ends)

    allocate(sd_sums_local(size(option%analysis_mols), 3, 0:option%delta-1), &
             coms      (omp_mol_start:omp_mol_end, 3),                       &
             coms_prev (omp_mol_start:omp_mol_end, 3),                       &
             dcoms_sums(omp_mol_start:omp_mol_end, 3),                       &
             dcoms_buf (omp_mol_start:omp_mol_end, 3, 0:option%delta-1))
    dcoms_buf     = 0_wp
    coms          = 0_wp
    sd_sums_local = 0_wp
    !$omp end parallel

    md_steps_analysis = sum(ceiling(real(trj_list%md_steps(:), wp) &
      / trj_list%ana_periods(:)))

    ! compute number of samples for each time value
    if (option%oversample) then
      forall (i=1:option%delta)
        nsamples(i) = md_steps_analysis - i
      end forall
    else
      forall (i=1:option%delta)
        nsamples(i) = (md_steps_analysis - 1) / i
      end forall
    end if

    ! open first trajectory to read the first step for smart unwrapping
    ! make sure coord shape matches with other trajectory files
    do ifile = 1, num_trjfiles

      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)
      call read_trj(trj_in, trajectory)
      call close_trj(trj_in)

      if (ifile == 1) then
        ! first unwrap system according to molecules
        call unwrap_molecules(trajectory%pbc_box, option%all_mols, &
          trajectory%coord)

        ! unwrap again for all analysis molecules to align substructures
        call unwrap_molecules(trajectory%pbc_box, all_analysis_mols, &
          trajectory%coord)

        allocate(coord_prev(3, size(trajectory%coord, 2)))
        coord_prev(:, :) = trajectory%coord

      else

        if (any(shape(coord_prev) /= shape(trajectory%coord))) then
          call error_msg("Analyze> Wrong particle number in " &
            // trim(trj_list%filenames(ifile)))
        end if

      end if

    end do


    ! global analysis step number
    istep = 0

    do ifile = 1, num_trjfiles

      write(MsgOut, '(A)') "Analyze> Reading " &
        // trim(trj_list%filenames(ifile))
      msg_interval = max(floor(real(trj_list%md_steps(ifile), wp) &
        / (10 * trj_list%ana_periods(ifile))), 1)

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep_local = 1, trj_list%md_steps(ifile)

        if (mod(istep_local, msg_interval) == 0) then
          write(MsgOut, "(A, i0, A, i0, A, i0, A, i0)") &
            "Analyze> File ", ifile, " of ", size(trj_list%filenames), &
            "; Step ", istep_local, " of ", trj_list%md_steps(ifile)
        end if

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord)
        !
        call read_trj(trj_in, trajectory)

        ! respect analysis periods setting
        if (mod(istep_local-1, trj_list%ana_periods(ifile)) /= 0) cycle

        istep = istep + 1

        com = 0_wp
        !$omp parallel default(none)                               &
        !$omp private(iaxis, com_local)                            &
        !$omp shared(trajectory, coord_prev, molecule, com)        &
        !$omp shared(option, istep, all_analysis_mols, imol2iaset)
        ! unwrap coordinates according to previous coordinates
        call unwrap_according_to_previous(&
          trajectory%pbc_box, &
          coord_prev(:, omp_atm_start:omp_atm_end), &
          trajectory%coord(:, omp_atm_start:omp_atm_end))
        coord_prev(:, omp_atm_start:omp_atm_end) &
            = trajectory%coord(:, omp_atm_start:omp_atm_end)

        !system center of mass
        do iaxis = 1, 3
          com_local(iaxis) = sum(molecule%mass(omp_atm_start:omp_atm_end) &
            * trajectory%coord(iaxis, omp_atm_start:omp_atm_end)) &
            / sum(molecule%mass)
        end do

        !$omp critical
        com = com + com_local
        !$omp end critical
        !$omp barrier

        coms_prev(omp_mol_start:omp_mol_end, :) &
          = coms(omp_mol_start:omp_mol_end, :)
        do imol = omp_mol_start, omp_mol_end
          do iaxis = 1,3
            coms(imol, iaxis) = sum(all_analysis_mols(imol)%mass &
            * trajectory%coord(iaxis, all_analysis_mols(imol)%atom_idx)) &
            / all_analysis_mols(imol)%molecule_mass(1) - com(iaxis)
          end do
        end do

        dcoms_buf(omp_mol_start:omp_mol_end, :,  mod(istep, option%delta)) &
            = coms(omp_mol_start:omp_mol_end, :) &
            - coms_prev(omp_mol_start:omp_mol_end, :)

        dcoms_sums(omp_mol_start:omp_mol_end, :) = 0

        do i = 0, min(istep-2, option%delta-1)
          dcoms_sums(omp_mol_start:omp_mol_end, :) &
            = dcoms_sums(omp_mol_start:omp_mol_end, :) &
              + dcoms_buf(omp_mol_start:omp_mol_end, :, &
                                                modulo(istep-i, option%delta))
          if (.not. option%oversample) then
            if (mod(istep-1, i+1) /= 0) cycle
          end if
          do iaxis = 1,3
            do iaset = 1, size(option%analysis_mols)
              sd_sums_local(iaset,iaxis,i) &
                = sd_sums_local(iaset,iaxis,i) &
                + sum(dcoms_sums(omp_mol_starts(iaset):omp_mol_ends(iaset), &
                                                                      iaxis)**2)
            end do
          end do
        end do
        !$omp end parallel

      end do

      call close_trj(trj_in)

    end do

    sd_sums = 0_wp
    !$omp parallel default(none) &
    !$omp shared(option)         &
    !$omp private(iaset)         &
    !$omp reduction(+:sd_sums)
    do iaset = 1, size(option%analysis_mols)
      sd_sums(iaset, :) = sum(sd_sums_local(iaset, option%axes(iaset)%i, :), 1)
    end do
    !$omp end parallel

    msd = 0

    ! prevent division by zero
    where (nsamples == 0)
      nsamples = 1
    end where
    where (nmols_aset == 0)
      nmols_aset = 1
    end where

    do iaset = 1, size(option%analysis_mols)
      msd(iaset, :) = sd_sums(iaset, :) / nsamples(:) / nmols_aset(iaset)
    end do


    if (output%msdfile /= '') then
      do i = 1, option%delta
        write(msd_out, '(i0)', advance="no") i
        do j = 1, size(option%analysis_mols)
          write(msd_out, '(x, es25.16e3)', advance="no") msd(j, i)
        end do
        write(msd_out, '()')
      end do
    end if

    ! close output file
    !
    if (output%msdfile /= '') call close_file(msd_out)

    ! Output summary
    !
    write(MsgOut,'(A)')      ''
    write(MsgOut,'(A)')      'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)')      ''
    write(MsgOut,'(A)')      '  [msdfile] ' // trim(output%msdfile)
    write(MsgOut,'(A)')      '    Column 1: time'
    if (size(option%analysis_mols)==1) then
      write(MsgOut,'(A)')    '    Column 2: MSD'
    else
      write(MsgOut,'(A,i0,A)') '    Columns 2-', size(option%analysis_mols)+1, ': MSD of each selection'
    end if
    write(MsgOut,'(A)')      '    Time units are set to the trajectory file interval'
    write(MsgOut,'(A)')      '    Distance units are those of the trajectory file'
    write(MsgOut,'(A)')      ''
    write(MsgOut,'(A)')      '    Number of axes used for each MSD:'
    do i  = 1, size(option%axes)
      write(MsgOut,'(A,i0,A,i0)') &
                             '    Column ', i+1, ': ', size(option%axes(i)%i)
    end do
    write(MsgOut, '(A)')     ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_analysis_set_molecule_ranges
  !> @brief        set molecule ranges for each analysis set
  !! @authors      DS
  !! @param[in]    all_analysis_mols : array of all analysis molecules
  !! @param[in]    imol2iaset        : molecule to analysis set mapping
  !! @param[in]    omp_size          : OpenMP size
  !! @param[in]    omp_rank          : OpenMP rank
  !! @param[out]   mol_starts        : molecule starts for each analysis set
  !! @param[out]   mol_ends          : molecule ends for each analysis set
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_analysis_set_molecule_ranges &
       (all_analysis_mols, imol2iaset, omp_size, omp_rank, mol_starts, mol_ends)

    type(s_one_molecule), dimension(:),          intent(in)  :: all_analysis_mols
    integer, dimension(size(all_analysis_mols)), intent(in)  :: imol2iaset
    integer,                                     intent(in)  :: omp_size
    integer,                                     intent(in)  :: omp_rank
    integer,         dimension(:), allocatable,  intent(out) :: mol_starts
    integer,         dimension(:), allocatable,  intent(out) :: mol_ends

    integer                            :: iaset, naset
    integer                            :: nmol, mol_start, mol_end
    integer, dimension(2)              :: mol_start_end
    integer, dimension(:), allocatable :: mol_counts

    naset = maxval(imol2iaset)
    nmol = size(all_analysis_mols)

    allocate(mol_starts(naset), mol_ends(naset), mol_counts(naset))

    mol_start_end = get_start_end(nmol, omp_size, omp_rank)
    mol_start = mol_start_end(1)
    mol_end = mol_start_end(2)

    do iaset = 1, naset
      mol_counts(iaset) = count(imol2iaset(mol_start:mol_end) == iaset)
    end do

    mol_starts(1) = mol_start

    do iaset = 1, naset-1

      mol_ends(iaset) = mol_starts(iaset) + mol_counts(iaset) - 1
      mol_starts(iaset+1) = mol_ends(iaset) + 1

    end do

    mol_ends(naset) = mol_end

    return

  end subroutine set_analysis_set_molecule_ranges

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_start_end
  !> @brief        a simple function to give starting and ending numbers
  !! @authors      DS
  !! @return       an integer array consisting of two elements
  !! @param[in]    nmol              : number of molecules
  !! @param[in]    omp_size          : OpenMP size
  !! @param[in]    omp_rank          : OpenMP rank
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_start_end(nmol, omp_size, omp_rank) result(res)

    integer, intent(in)   :: nmol, omp_size, omp_rank
    integer, dimension(2) :: res
    integer               :: stride, start, end

    stride = nmol / omp_size
    start = omp_rank * stride + 1
    if (omp_rank + 1 == omp_size) then
      end = nmol
    else
      end = start + stride - 1
    end if
    res = [start, end]

    return

  end function get_start_end

end module ma_analyze_mod
