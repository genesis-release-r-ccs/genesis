!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module da_analyze_mod

  use da_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use fileio_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: analyze

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM, DM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: iatom, jatom, mat_out, idx1, idx2
    integer                  :: natoms1, natoms2
    real(wp)                 :: diff(3)
    real(wp)                 :: dist2
    real(wp), allocatable    :: distance(:,:)

    if (option%check_only) &
      return


    ! open output file
    !
    if (output%outfile /= '') &
      call open_file(mat_out, output%outfile, IOFileOutputNew)


    ! allocate distance_matrix
    !
    if (option%calc_mode == intramolecule) then
      natoms1 = size(option%analysis_atom_1%idx)
      allocate(distance(natoms1, natoms1))

    else
      natoms1 = size(option%analysis_atom_1%idx)
      natoms2 = size(option%analysis_atom_2%idx)
      allocate(distance(natoms1, natoms2))

    end if

    distance(:,:) = 0.0_wp


    ! analysis loop
    !
    nstru = 0
    num_trjfiles = size(trj_list%md_steps)

    do ifile = 1, num_trjfiles

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

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru


          ! intra-molecule distance mode

          if (option%calc_mode == intramolecule) then
            do iatom = 1, natoms1-1
              idx1 = option%analysis_atom_1%idx(iatom)

              do jatom = iatom+1, natoms1
                idx2  = option%analysis_atom_1%idx(jatom)

                diff(:) = trajectory%coord(:, idx1) - trajectory%coord(:,idx2)
                dist2   = dot_product(diff, diff)
                distance(iatom,jatom) = distance(iatom,jatom) + sqrt(dist2)

              end do
            end do

          else if (option%calc_mode == intermolecule) then
            do jatom = 1, natoms2
              idx2 = option%analysis_atom_2%idx(jatom)

              do iatom = 1, natoms1
                idx1  = option%analysis_atom_1%idx(iatom)

                diff(:) = trajectory%coord(:, idx1) - trajectory%coord(:,idx2)
                dist2   = dot_product(diff, diff)
                distance(iatom, jatom) = distance(iatom, jatom) + sqrt(dist2)

              end do
            end do
          end if

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! average
    !
    if (nstru /= 0) then
      distance(:,:) = distance(:,:)/real(nstru)
    end if


    ! output results
    !
    if (output%outfile /= '') then

      if (option%matrix_shape == MatrixShapeFull) then

        if (option%calc_mode == intramolecule) then
          do iatom = 1, natoms1 - 1
            do jatom = iatom+1, natoms1
              distance(jatom,iatom) = distance(iatom,jatom)
            end do
          end do

          do iatom = 1, natoms1
            do jatom = 1, natoms1
              write(mat_out,'(i10,i10,f12.5)') iatom, jatom, distance(iatom, jatom)
            end do
            write(mat_out,'(1x)')
          end do

        else if (option%calc_mode == intermolecule) then

          do jatom = 1, natoms2
            do iatom = 1, natoms1
              write(mat_out,'(i10,i10,f12.5)') iatom, jatom, distance(iatom, jatom)
            end do
            write(mat_out,'(1x)')
          end do
        end if

      else if (option%matrix_shape == MatrixShapeHalf) then

        do iatom = 1, natoms1 - 1
          do jatom = iatom+1, natoms1
            write(mat_out,'(i10,i10,f12.5)') iatom, jatom, distance(iatom,jatom)
          end do
        end do

      end if

    end if


    ! close output file
    !
    if (output%outfile /= '') call close_file(mat_out)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
    write(MsgOut,'(A)') '    Column 1: i-th atom'
    write(MsgOut,'(A)') '    Column 2: j-th atom'
    write(MsgOut,'(A)') '    Column 3: averaged distance between i- and j-th atoms (angstrom)'
    write(MsgOut,'(A)') ''


    deallocate(distance)

    return

  end subroutine analyze

end module da_analyze_mod
