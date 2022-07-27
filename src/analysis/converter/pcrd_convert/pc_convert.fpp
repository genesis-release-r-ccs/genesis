!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pc_convert_mod
!> @brief   convert trajectory files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pc_convert_mod

  use pc_option_mod
  use pc_option_str_mod
  use parallel_trj_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use molecules_str_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
 
  implicit none
  private

  ! subroutines
  public :: convert

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert trajectory files
  !! @authors      NT
  !! @param[inout] molecule    molecule information
  !! @param[inout] trj_list    trajectory list information
  !! @param[inout] trajectory  trajectory information
  !! @param[inout] fitting     fitting information
  !! @param[inout] option      option information
  !! @param[inout] output      output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine convert(molecule,   &
                    trj_list,   &
                    trajectory, &
                    fitting,    &
                    option,     &
                    output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: nstru, irun, itrj
    integer                  :: rms_out, trr_out


    ! check-only
    if (option%check_only) &
      return

    ! open output files
    if (output%trjfile /= '') &
      call open_trj (trj_out,              &
                     output%trjfile,       &
                     option%trjout_format, &
                     option%trjout_type,   &
                     IOFileOutputNew)

    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    if (output%trrfile /= '') &
      call open_file(trr_out, output%trrfile, IOFileOutputNew)

    nstru = 0

    do irun = 1, size(trj_list%md_steps)

      call open_parallel_trj(trj_in, trj_list%filenames(irun))

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_parallel_trj(trj_in, trajectory)

        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru          

          ! selection
          !
          call reselect_atom(molecule, &
                             option%trjout_atom_exp, &
                             trajectory%coord, &
                             option%trjout_atom, &
                             option%trjout_atom_trj)

          ! pbc-correct
          !
          call run_pbc_correct(option%pbcc_mode, &
                               molecule, &
                               trajectory)

          ! fitting
          !
          call run_fitting(fitting, &
                           molecule%atom_coord, &
                           trajectory%coord, &
                           trajectory%coord)

          ! write data
          !
          if (output%rmsfile /= '') &
            call out_rmsd (rms_out, nstru, fitting)

          if (output%trrfile /= '') &
            call out_trrot(trr_out, nstru, fitting)

          if (output%trjfile /= '') &
            call write_trj(trj_out, trajectory, option%trjout_atom_trj,molecule)

        end if

      end do

      call close_parallel_trj(trj_in)

    end do

    if (output%trrfile /= '') call close_file(trr_out)
    if (output%rmsfile /= '') call close_file(rms_out)
    if (output%trjfile /= '') call close_trj (trj_out)

    return

  end subroutine convert

end module pc_convert_mod
