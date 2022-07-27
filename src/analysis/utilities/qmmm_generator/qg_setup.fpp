!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qg_setup_mod
!> @brief   setup variables and structures in qmmm_generator
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qg_setup_mod

  use qg_control_mod
  use qg_option_mod
  use qg_option_str_mod
  use fitting_mod
  use fitting_str_mod
  use input_mod
  use output_mod
  use output_str_mod
  use trajectory_mod
  use trajectory_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_pdb_mod
 
  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in qmmm_generator
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] ref        : pdb information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data,  &
                   molecule,   &
                   ref,        &
                   trj_list,   &
                   trajectory, &
                   fitting,    &
                   option,     &
                   output)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_pdb),             intent(inout) :: ref 
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_psf)              :: psf ! , psf_out
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd


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

    call duplicate_psf(psf, option%dup_psf)
    if (option%dup_psf%type == PsfTypeXPLOR) then
      option%dup_psf%type = PsfTypeXPLOREXT
    else
      option%dup_psf%type = PsfTypeCHARMMEXT
    end if
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, &
                          molecule, trj_list, trajectory)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup fitting
    !
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      trj_list, molecule, option)

    ! setup output 
    !
    call setup_output(ctrl_data%out_info, output)


    !! export reference molecules (for analysis ...)
    !!
    !if (output%qmmm_pdbfile /= '') then ! qmmm_reffile
!
!      call export_molecules(molecule, option%qmmm_atom, ref_out)
!      call output_pdb(output%qmmm_pdbfile, ref_out) ! qmmm_reffile
!      call dealloc_pdb_all(ref_out)
!
!    end if

!    ! export reference PSFs (currently only psf is allowed) 
!    !
!    if (output%qmmm_psffile /= '') then 
!
!      if(len_trim(ctrl_data%inp_info%psffile) > 0) then
!        !call input_psf(ctrl_data%inp_info%psffile, psf)
!        call export_molecules(molecule, option%qmmm_atom, psf=psf)
!        call output_psf(output%qmmm_psffile, psf)
!        call dealloc_psf_all(psf)
!        !call dealloc_psf_all(psf_out)
!      end if
!
!    endif

    return

  end subroutine setup

end module qg_setup_mod
