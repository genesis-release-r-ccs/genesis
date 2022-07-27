!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sa_setup_mod
!> @brief   setup variables and structures in TRJ_ANALYSIS
!! @authors Norio Takase (NT), Takaharu Mori (TM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module sa_setup_mod

  use sa_control_mod
  use sa_option_mod
  use sa_domain_mod
  use sa_boundary_mod
  use sa_ensemble_mod
  use sa_option_str_mod
  use sa_domain_str_mod
  use sa_boundary_str_mod
  use sa_ensemble_str_mod
  use trajectory_mod
  use fitting_mod
  use output_mod
  use input_mod
  use trajectory_str_mod
  use fitting_str_mod
  use output_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_top_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use mpi_parallel_mod

  implicit none
  private

  ! subroutines
  public  :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      IY
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] boundary   : boundary information
  !! @param[inout] domain     : domain information
  !! @param[inout] ensemble   : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, trj_list, trajectory, fitting, output, &
                   option, boundary, domain, ensemble)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(inout) :: ensemble

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref
    type(s_pdb)              :: pdb
    type(s_top)              :: top
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_ambcrd)           :: ambref
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd
    type(s_grocrd)           :: groref


    ! input files
    !
    call input_files(ctrl_data%inp_info, &
                     psf=psf,            &
                     ref=ref,            &
                     pdb=pdb,            &
                     top=top,            &
                     prmtop=prmtop,      &
                     ambcrd=ambcrd,      &
                     ambref=ambref,      &
                     grotop=grotop,      &
                     groref=groref,      &
                     grocrd=grocrd)

    ! define molecules
    ! defined in /lib/molecule.fpp
    call define_molecules(molecule, ref=ref,       &
                                    psf=psf,       &
                                    pdb=pdb,       &
                                    top=top,       &
                                    prmtop=prmtop, &
                                    ambcrd=ambcrd, &
                                    ambref=ambref, &
                                    grotop=grotop, &
                                    groref=groref, &
                                    grocrd=grocrd)

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(pdb)
    call dealloc_top_all(top)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(groref)
    call dealloc_grocrd_all(grocrd)

    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, molecule, trj_list, trajectory)

    ! setup selection
    ! /lib/select.fpp
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup fitting
    !
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)

    ! setup ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, ensemble)

    ! setup boundary
    ! (IY) some input valuables are fixed to use constant values,
    !
    call setup_boundary(ctrl_data%bound_info,  &
                        .true. ,               &
                        13.5_wp ,              &
                        ensemble%ensemble,     &
                        molecule, boundary,    &
                        trj_list, trajectory,  &
                        option)

    call setup_domain(boundary, molecule, domain, option)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    return

  end subroutine setup

end module sa_setup_mod
