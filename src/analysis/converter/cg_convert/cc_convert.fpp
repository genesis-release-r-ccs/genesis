!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cc_convert_mod
!> @brief   convert trajectory files
!! @authors Norio Takase (NT) Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module cc_convert_mod

  use cc_option_mod
  use cc_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_mod
  use molecules_str_mod
  use constants_mod
  use measure_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use string_mod
  use fileio_pdb_mod
 
  implicit none
  private

  ! subroutines
  public  :: convert
  private  :: calculate_com

  real(wp), allocatable, save   :: tmp_coord(:,:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert trajectory files
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !! @param[inout] trajectory_out : trajectory (out) information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine convert(molecule,     &
                    trj_list,     &
                    trajectory,   &
                    option,       &
                    output,       &
                    molecule_out, &
                    trajectory_out)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule_out
    type(s_trajectory),      intent(inout) :: trajectory_out

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: nstru, irun, itrj
    integer                  :: rms_out, trr_out
    integer                  :: i

    type(s_selatoms)              :: seledumm


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

    call alloc_selatoms(seledumm,option%num_groups)
    do i = 1, option%num_groups
     seledumm%idx(i) = i
    end do

    nstru = 0

    do irun = 1, size(trj_list%md_steps)

      call open_trj(trj_in,                   &
                    trj_list%filenames(irun), &
                    trj_list%trj_format,      &
                    trj_list%trj_type,        &
                    IOFileInput)

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_trj(trj_in, trajectory)

        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru          

          ! selection
          !
            call calculate_com(molecule,  option, &
                               trajectory%coord, trajectory_out%coord)
                               

          if (output%trjfile /= '') then
            call write_trj(trj_out, trajectory_out, seledumm, molecule_out)
          end if

        end if

      end do

      call close_trj(trj_in)

    end do

    if (output%trjfile /= '') call close_trj (trj_out)

    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calculate_com
  !> @brief        output the COM of the selected groups
  !! @authors      CK
  !! @param[in]    molecule   : molecule information
  !! @param[in]    coord      : atom coordinates
  !! @param[out]   coord_out  : com coordinates
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calculate_com(molecule, option, coord, coord_out)

  ! formal arguments
  type(s_molecule), intent(in)    :: molecule
  type(s_option),   intent(in)    :: option
  real(wp),         intent(inout) :: coord(:,:)
  real(wp),         intent(out)   :: coord_out(:,:)

  ! local variables
  integer  :: i, nselect
  real(wp) :: com(3)

  do i = 1, option%num_groups
    nselect = option%num_atoms(i)
    com   = compute_com(coord, molecule%mass, option%atomlist(1:nselect,i))
    coord_out(1:3,i) = com(1:3)
  end do

  return

  end subroutine calculate_com

end module cc_convert_mod
