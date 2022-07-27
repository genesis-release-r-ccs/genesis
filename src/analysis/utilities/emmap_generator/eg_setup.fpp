!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   eg_setup_mod
!> @brief   setup variables and structures
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module eg_setup_mod

  use eg_control_mod
  use eg_option_mod
  use eg_option_str_mod
  use input_mod
  use output_mod
  use output_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_ambcrd_mod
  use fileio_pdb_mod
 
  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, option, output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     pdb=pdb,            &
                     ambcrd=ambcrd,      &
                     grocrd=grocrd)


    ! define molecules
    !
    call define_molecules(molecule, pdb=pdb,       &
                                    ambcrd=ambcrd, &
                                    grocrd=grocrd)

    call dealloc_pdb_all(pdb)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grocrd_all(grocrd)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, option)

    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module eg_setup_mod
