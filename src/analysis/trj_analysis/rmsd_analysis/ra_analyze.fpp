!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_analyze_mod

  use ra_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
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
  !! @authors      TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory, fitting)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(in)    :: trj_list
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting

    ! local variables
    type(s_trj_file)         :: trj_in
    integer                  :: nstru, ifile, istep, num_trjfiles
    integer                  :: iatom, rms_out, idx, i
    real(wp)                 :: rmsd, tot_mass

    real(wp), allocatable    :: mass_fitting(:), mass_analysis(:)


    if (option%check_only) &
      return

!    if (fitting%mass_weight) &
!      call error_msg('Analyze> mass weighted is not allowed')

    if (fitting%mass_weight .and. .not. option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted fitting is enable while RMSD is not mass-weighted '
    else if (.not. fitting%mass_weight .and.  option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted RMSD is enable while fitting is not mass-weighted '
    endif

    allocate(mass_fitting(1:size(molecule%mass(:))), &
             mass_analysis(1:size(option%analysis_atom%idx)))

    if (fitting%mass_weight) then
      if (abs(molecule%mass(1)) < 1.0e-05_wp) then
         call error_msg('Analyze> mass is not defined')
      else
         mass_fitting(:) = molecule%mass(:)
      endif
    else
      mass_fitting(:)=1.0_wp
    endif

    if (option%mass_weight) then
      tot_mass = 0.0_wp
      do iatom = 1, size(option%analysis_atom%idx)
        idx = option%analysis_atom%idx(iatom)
        mass_analysis(iatom) = molecule%mass(idx)
        tot_mass = tot_mass + mass_analysis(iatom)
      end do
    else
      mass_analysis(:)=1.0_wp
      tot_mass = real(size(option%analysis_atom%idx),wp)
    endif
    

    ! open output file
    !
    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)


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


          ! fitting
          !
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord, &
                           mass_fitting)


          ! compute RMSD for selected atoms
          !
          rmsd = 0.0_wp
          do iatom = 1, size(option%analysis_atom%idx)

            idx = option%analysis_atom%idx(iatom)
            if (fitting%fitting_method == FittingMethodXYTR_ZROT) then
              rmsd = rmsd + mass_analysis(iatom) * (                                    &
                + (molecule%atom_coord(1,idx) - trajectory%coord(1,idx))**2  &
                + (molecule%atom_coord(2,idx) - trajectory%coord(2,idx))**2)
            else
              rmsd = rmsd + mass_analysis(iatom) * (                                    &
              + (molecule%atom_coord(1,idx) - trajectory%coord(1,idx))**2  &
              + (molecule%atom_coord(2,idx) - trajectory%coord(2,idx))**2  &
              + (molecule%atom_coord(3,idx) - trajectory%coord(3,idx))**2)
            endif

          end do
          rmsd = sqrt(rmsd/tot_mass)


          ! output results
          !
          write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ',rmsd
          write(MsgOut,*) ''

          if (output%rmsfile /= '') &
            write(rms_out, '(i10,1x,f10.5)') nstru, rmsd

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    if (output%rmsfile /= '') call close_file(rms_out)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rmsfile] ' // trim(output%rmsfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Root-mean-square deviation (RMSD)(angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

end module ra_analyze_mod
