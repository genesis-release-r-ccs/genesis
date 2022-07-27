!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ma_setup_mod
!> @brief   setup variables and structures
!! @authors Donatas Surblys (DS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ma_setup_mod

  use ma_control_mod
  use ma_option_mod
  use ma_option_str_mod
  use trajectory_mod
  use output_mod
  use input_mod
  use trajectory_str_mod
  use output_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use select_molecules_mod

  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, trj_list, trajectory, output, option)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_pdb)                                :: pdb
    type(s_psf)                                :: psf
    type(s_prmtop)                             :: prmtop
    type(s_ambcrd)                             :: ambcrd

    type(s_selmols), dimension(:), allocatable      :: selmols
    type(s_one_molecule), dimension(:), allocatable :: allmols


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
      pdb=pdb,                           &
      psf=psf,                           &
      prmtop=prmtop,                     &
      ambcrd=ambcrd)

    ! define molecules
    !
    call define_molecules(molecule, &
      pdb=pdb,                      &
      psf=psf,                      &
      prmtop=prmtop,                &
      ambcrd=ambcrd)

    call dealloc_pdb_all(pdb)
    call dealloc_psf_all(psf)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, molecule, trj_list, trajectory)


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup molecule selection
    !
    call setup_molselect(ctrl_data%molsel_info, ctrl_data%sel_info, &
      molecule, selmols, allmols=allmols)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, selmols, allmols, molecule, option)


    return

  end subroutine setup

end module ma_setup_mod
