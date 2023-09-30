!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_setup_mod
!> @brief   setup variables and structures in TRJ_ANALYSIS
!! @authors Motoshi Kamiya (MK)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_setup_mod

  use ea_control_mod
  use ea_option_mod
  use ea_option_str_mod
  use trajectory_mod
  use trajectory_str_mod

  use at_input_mod
  use at_enefunc_mod
  use at_enefunc_str_mod
  use at_pairlist_mod
  use at_pairlist_str_mod
  use at_dynvars_mod
  use at_dynvars_str_mod
  use at_boundary_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_setup_mpi_mod
  use at_restraints_mod
  use at_restraints_str_mod
  use at_remd_str_mod
  use at_remd_mod
  use at_energy_mod
  use at_energy_str_mod
  use at_dynamics_str_mod
  use at_dynamics_mod
  use at_ensemble_str_mod
  use at_ensemble_mod
  use at_morph_mod
  use at_qmmm_mod

  use output_mod
  use output_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_table_mod
  use fileio_mode_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_rst_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_crd_mod
  use fileio_eef1_mod
  use fileio_morph_mod
  use fileio_spot_mod
  use messages_mod

  implicit none
  private

  ! subroutines
  public :: setup

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      MK
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] enefunc    : energy function information
  !! @param[inout] pairlist   : nonbonded pairlist information
  !! @param[inout] dynvars    : dynamic variables
  !! @param[inout] boundary   : boundary information
  !! @param[inout] remd       : remd information
  !! @param[inout] ensemble   : ensemble information
  !! @param[inout] dynamics   : dynamics information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, output, molecule, enefunc, pairlist, &
                   dynvars, boundary, remd, ensemble, dynamics,    &
                   trj_list, trajectory, option)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_boundary),        intent(inout) :: boundary
    type(s_remd),            intent(inout) :: remd
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_eef1)             :: eef1
    type(s_table)            :: table
    type(s_morph_in)         :: morph_in
    type(s_spot)             :: spot


    ! disable QM/MM
    !
    enefunc%qmmm%do_qmmm = .false.

    ! setup MPI
    !
    call setup_mpi_md(ctrl_data%ene_info)

    ! input min (not libana/input but is atdyn/at_input)
    !
    call input_min(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref,     &
                   mode, eef1, table, morph_in, spot)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup dynamics
    !
    dynamics%integrator       = ctrl_data%dyn_info%integrator
    dynamics%nsteps           = ctrl_data%dyn_info%nsteps
    dynamics%timestep         = ctrl_data%dyn_info%timestep
    dynamics%crdout_period    = ctrl_data%dyn_info%crdout_period
    dynamics%velout_period    = ctrl_data%dyn_info%velout_period
    dynamics%eneout_period    = ctrl_data%dyn_info%eneout_period
    dynamics%rstout_period    = ctrl_data%dyn_info%rstout_period
    dynamics%stoptr_period    = ctrl_data%dyn_info%stoptr_period
    dynamics%nbupdate_period  = ctrl_data%dyn_info%nbupdate_period
    dynamics%initial_time     = ctrl_data%dyn_info%initial_time
    dynamics%annealing        = ctrl_data%dyn_info%annealing
    dynamics%anneal_period    = ctrl_data%dyn_info%anneal_period
    dynamics%dtemperature     = ctrl_data%dyn_info%dtemperature
    dynamics%istart_step      = 1
    dynamics%iend_step        = dynamics%nsteps
    dynamics%verbose          = ctrl_data%dyn_info%verbose
    dynamics%target_md        = ctrl_data%dyn_info%target_md
    dynamics%steered_md       = ctrl_data%dyn_info%steered_md
    dynamics%initial_value    = ctrl_data%dyn_info%initial_value
    dynamics%final_value      = ctrl_data%dyn_info%final_value
    dynamics%restart          = .false.
    !! dynamics
    dynamics%crdout_period    = dynamics%nsteps
    dynamics%velout_period    = dynamics%nsteps
    dynamics%rstout_period    = dynamics%nsteps

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, molecule, option)

    ! boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist,                &
                        ctrl_data%sel_info,                             &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup solute tempering
    !
    if ( len_trim(option%remfile) /= 0 .or. option%rest_component ) then
      call setup_solute_tempering(ctrl_data%rep_info, &
                                  molecule, restraints, ctrl_data%cons_info)
    end if

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary,  &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, enefunc%qmmm) 

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info,                    &
                        boundary,                              &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)
    call setup_enefunc_morph_in(morph_in, enefunc)

    call dealloc_psf_all(psf)
    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! setup ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, boundary, enefunc, ensemble)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, molecule, trj_list, trajectory)

    ! setup remd
    ! The route cannot be selected.
    if ( len_trim(option%remfile) /= 0 .or. &
         option%rest_component .or. &
         option%mbar ) then
      if ( option%rest_component .or. option%mbar ) then
        rst%rstfile_type = RstfileTypeMd
      end if
      call setup_remd(ctrl_data%rep_info, rst, boundary, dynamics, &
                      molecule, restraints, ensemble, enefunc, remd)
    end if
    call dealloc_rst_all(rst)
    call dealloc_restraints_all(restraints)
    if (remd%rest_mixed_three_atoms) &
        call error_msg('Setup> This combination type of solute/solvent is not allowed in this program. Please execute spdyn with analysis_grest option.')

    return

  end subroutine setup

end module ea_setup_mod
