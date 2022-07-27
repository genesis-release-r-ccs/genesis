!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   wa_setup_mod
!> @brief   setup variables and structures in WHAM_ANALYSIS
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module wa_setup_mod

  use wa_control_mod
  use wa_option_mod
  use wa_option_str_mod
  use output_mod
  use input_mod
  use output_str_mod
  use input_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in WHAM_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, option, input, output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     pdb=pdb,            &
                     psf=psf,            &
                     prmtop=prmtop,      &
                     ambcrd=ambcrd,      &
                     grotop=grotop,      &
                     grocrd=grocrd)


    ! define molecules
    !
    if (ctrl_data%inp_info%cvfile == '') then
      call define_molecules(molecule, pdb=pdb,       &
                                      psf=psf,       &
                                      prmtop=prmtop, &
                                      ambcrd=ambcrd, &
                                      grotop=grotop, &
                                      grocrd=grocrd)
    end if

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(pdb)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)


    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)


    return

  end subroutine setup

end module wa_setup_mod
