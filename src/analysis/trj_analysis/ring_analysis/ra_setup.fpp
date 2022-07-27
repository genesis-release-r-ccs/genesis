!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_setup_mod
!> @brief   setup variables and structures in TRJ_ANALYSIS
!! @authors Norio Takase (NT), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_setup_mod

  use ra_control_mod
  use ra_option_mod
  use ra_option_str_mod
  use trajectory_mod
  use output_mod
  use input_mod
  use trajectory_str_mod
  use output_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_morph_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_pdb_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, output, option)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd

    integer                  :: i

    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     psf=psf,            &
                     ref=ref,            &
                     prmtop=prmtop,      &
                     ambcrd=ambcrd,      &
                     grotop=grotop,      &
                     grocrd=grocrd)


    ! define molecules
    !
    call define_molecules(molecule, pdb=ref,       &
                                    psf=psf,       &
                                    prmtop=prmtop, &
                                    ambcrd=ambcrd, &
                                    grotop=grotop, &
                                    grocrd=grocrd)

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grocrd_all(grocrd)


    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, molecule,  option)

    return

  end subroutine setup

end module ra_setup_mod
