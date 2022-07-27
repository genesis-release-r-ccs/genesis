!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_setup_atdyn_mod
!> @brief   setup variables and structures in MD simulaton
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_setup_atdyn_mod

  use at_control_mod
  use at_output_mod
  use at_input_mod
  use at_remd_mod
  use at_rpath_mod
  use at_morph_mod
  use at_minimize_mod
  use at_vibration_mod
  use at_dynamics_mod
  use at_dynvars_mod
  use at_ensemble_mod
  use at_enefunc_table_mod
  use at_restraints_mod
  use at_constraints_mod
  use at_boundary_mod
  use at_pairlist_mod
  use at_enefunc_mod
  use at_energy_mod
  use at_output_str_mod
  use at_minimize_str_mod
  use at_vibration_str_mod
  use at_dynamics_str_mod
  use at_dynvars_str_mod
  use at_ensemble_str_mod
  use at_remd_str_mod
  use at_rpath_str_mod
  use at_morph_str_mod
  use at_restraints_str_mod
  use at_constraints_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_qmmm_mod
  use at_experiments_mod
  use select_mod
  use fitting_mod
  use fitting_str_mod
  use at_gamd_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_morph_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_rst_mod
  use fileio_psf_mod
  use fileio_gpr_mod
  use fileio_par_mod
  use fileio_top_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_mode_mod
  use fileio_rstmep_mod
  use fileio_eef1_mod
  use fileio_table_mod
  use fileio_spot_mod
  use structure_check_mod
  use messages_mod
  use mpi_parallel_mod
 
  implicit none
  private

  ! subroutines
  public  :: setup_atdyn_md
  public  :: setup_atdyn_min
  public  :: setup_atdyn_remd
  public  :: setup_atdyn_rpath
  public  :: setup_atdyn_vib
  public  :: setup_atdyn_morph

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_md
  !> @brief        setup variables and structures in MD simulation
  !! @authors      TM, YS
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                            dynvars, dynamics, constraints, ensemble, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary

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


    ! input md
    !
    call input_md(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  mode, eef1, table, morph_in, spot)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! set dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info, ctrl_data%bound_info, &
                        ctrl_data%res_info, ctrl_data%out_info, &
                        rst, molecule, dynamics)

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)

    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars, dynamics, &
                       ctrl_data%ens_info%tpcontrol)

    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary, &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, &
                    enefunc%qmmm)
    call dealloc_psf_all(psf)
    call dealloc_rst_all(rst)

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info, &
                        boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)

    call setup_enefunc_morph_in(morph_in, enefunc)

    call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             molecule, enefunc)

    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_eef1(eef1)

    if (.not. dynamics%target_md) enefunc%target_function = 0
    if (.not. dynamics%steered_md) enefunc%steered_function = 0

    ! setup ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, boundary, enefunc, ensemble)

    ! setup constraints
    !
    call setup_constraints(ctrl_data%cons_info, ctrl_data%sel_info, &
                           dynamics, boundary, molecule, &
                           enefunc, constraints)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ctrl_data%ene_info, molecule, constraints, &
                                enefunc)

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)
    call dealloc_restraints_all(restraints)

    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, molecule, enefunc)

    ! setup output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)') &
        'Setup_Atdyn_Md> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    ! setup ESP/MM-MD
    if (dynamics%esp_mm) &
      call setup_dynamics_espmm(enefunc, molecule, pairlist, dynamics, dynvars, enefunc%qmmm)

    return

  end subroutine setup_atdyn_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_min
  !> @brief        setup variables and structures in minimization
  !! @authors      TM, YS
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[out]   output    : information of output
  !! @param[out]   molecule  : information of molecules
  !! @param[out]   enefunc   : information of energy function
  !! @param[out]   pairlist  : information of nonbonded pairlist
  !! @param[out]   dynvars   : information of dynamic variables
  !! @param[out]   minimize  : information of minimize
  !! @param[out]   boundary  : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_min(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, minimize, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_boundary),        intent(inout) :: boundary

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


    ! input min
    !
    call input_min(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   mode, eef1, table, morph_in, spot)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)
    
    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars)
    

    ! setup structure_check
    !
    call setup_structure_check(molecule, &
                     ctrl_data%ene_info%forcefield_char, dynvars%coord, &
                     ctrl_data%min_info%check_structure, &
                     ctrl_data%min_info%fix_ring_error, &
                     ctrl_data%min_info%fix_chirality_error, &
                     ctrl_data%min_info%exclude_ring_grpid, &
                     ctrl_data%min_info%exclude_chiral_grpid)

    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary, &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, &
                    enefunc%qmmm)
    call dealloc_psf_all(psf)
    call dealloc_rst_all(rst)    

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info, &
                        boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)

    call setup_enefunc_morph_in(morph_in, enefunc)

    call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             molecule, enefunc)

    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_eef1(eef1)

    ! setup minimize
    !
    call setup_minimize(ctrl_data%min_info, ctrl_data%out_info, &
         ctrl_data%sel_info, molecule, enefunc%qmmm, boundary, minimize)

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)
    call dealloc_restraints_all(restraints)

    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)

    dynvars%verbose = minimize%verbose

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)') &
        'Setup_Atdyn_Min> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_atdyn_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_remd
  !> @brief        setup variables and structures in REMD simulation
  !! @authors      TM, YS
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   remd        : information of remd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_remd(ctrl_data, output, molecule, enefunc, pairlist, &
                              dynvars, dynamics, constraints, ensemble, &
                              boundary, remd)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_remd),            intent(inout) :: remd

    ! local variables
    integer                     :: i
    integer                     :: replicaid, parmsetid
    type(s_restraints)          :: restraints
    type(s_top)                 :: top
    type(s_par)                 :: par
    type(s_gpr)                 :: gpr
    type(s_psf)                 :: psf
    type(s_prmtop)              :: prmtop
    type(s_grotop)              :: grotop
    type(s_pdb)                 :: ref
    type(s_pdb)                 :: fit
    type(s_pdb)                 :: pdb
    type(s_crd)                 :: crd
    type(s_ambcrd)              :: ambcrd
    type(s_grocrd)              :: grocrd
    type(s_ambcrd)              :: ambref
    type(s_grocrd)              :: groref
    type(s_rst)                 :: rst
    type(s_mode)                :: mode
    type(s_eef1)                :: eef1
    type(s_table)               :: table
    type(s_morph_in)            :: morph_in
    type(s_spot)                :: spot


    ! input remd
    !
    call input_remd(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                    pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                    mode, eef1, table, morph_in, spot)

    ! define molecule
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info, ctrl_data%bound_info, &
                        ctrl_data%res_info, ctrl_data%out_info, &
                        rst, molecule, dynamics)

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)
    
    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars, dynamics, &
                       ctrl_data%ens_info%tpcontrol)
    
    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup REST (evacuation of water mols)
    !
    call setup_solute_tempering(ctrl_data%rep_info, &
                                molecule, restraints, ctrl_data%cons_info)

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary, &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, &
                    enefunc%qmmm)
    call dealloc_psf_all(psf)

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info, &
                        boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)

    call setup_enefunc_morph_in(morph_in, enefunc)

    call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             molecule, enefunc)

    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_eef1(eef1)


    ! setup ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, boundary, enefunc, ensemble)

    ! setup constraints
    !
    call setup_constraints(ctrl_data%cons_info, ctrl_data%sel_info, &
                           dynamics, boundary, molecule, &
                           enefunc, constraints)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ctrl_data%ene_info, molecule, constraints, &
                                enefunc)

    ! setup remd
    !
    call setup_remd(ctrl_data%rep_info, rst, boundary, dynamics, &
                    molecule, restraints, ensemble, enefunc, remd)
    call dealloc_rst_all(rst) 

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)
    call dealloc_restraints_all(restraints)

    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, molecule, enefunc, remd)

    ! set output
    !
    call setup_output_remd(ctrl_data%out_info, dynamics, remd, output)

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)') &
        'Setup_Atdyn_Remd> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    ! ESP/MM-MD is not available for REMD
    if (dynamics%esp_mm .and. main_rank) &
      call error_msg('Setup_Atdyn_Remd> ESP/MM is not available for REMD.')

    return

  end subroutine setup_atdyn_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_rpath
  !> @brief        setup variables and structures in RPATH simulation
  !! @authors      TM, YS, YK, KY
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   rpath       : information of rpath
  !! @param[out]   minimize    : information of minimize
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_rpath(ctrl_data, output, molecule, enefunc, pairlist, &
                              dynvars, dynamics, constraints, ensemble, &
                              boundary, rpath, minimize)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_rpath),           intent(inout) :: rpath
    type(s_minimize),        intent(inout) :: minimize

    ! local variables
    integer                     :: i
    integer                     :: replicaid, parmsetid
    type(s_restraints)          :: restraints
    type(s_top)                 :: top
    type(s_par)                 :: par
    type(s_gpr)                 :: gpr
    type(s_psf)                 :: psf
    type(s_prmtop)              :: prmtop
    type(s_grotop)              :: grotop
    type(s_pdb)                 :: ref
    type(s_pdb)                 :: fit
    type(s_pdb)                 :: pdb
    type(s_crd)                 :: crd
    type(s_ambcrd)              :: ambcrd
    type(s_grocrd)              :: grocrd
    type(s_ambcrd)              :: ambref
    type(s_grocrd)              :: groref
    type(s_rst)                 :: rst
    type(s_mode)                :: mode
    type(s_rstmep)              :: rstmep
    type(s_eef1)                :: eef1
    type(s_table)               :: table
    type(s_morph_in)            :: morph_in
    type(s_spot)                :: spot


    ! define replica
    !
    call define_nreplica(ctrl_data%rpath_info, rpath)
    my_replica_no = nrep_per_proc*(my_country_no) + 1

    ! input rpath
    !
    call input_rpath(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop,&
                    pdb, crd, ambcrd, grocrd, rst, ref, fit, ambref, &
                    groref, mode, rstmep, eef1, table, morph_in, spot)

    ! define molecule
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup dynamics
    !
    if (ctrl_data%rpath_info%rpathmode == RpathmodeMFEP  .or. &
        ctrl_data%rpath_info%rpathmode == RpathmodeMEPMD .or. &
        ctrl_data%rpath_info%rpathmode == RpathmodeFEP) then
      call setup_dynamics(ctrl_data%dyn_info, ctrl_data%bound_info, &
                          ctrl_data%res_info, ctrl_data%out_info, &
                          rst, molecule, dynamics)
    end if

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)
    
    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars, dynamics, &
                       ctrl_data%ens_info%tpcontrol)


    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary, &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, &
                    enefunc%qmmm)
    call dealloc_psf_all(psf)

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info, &
                        boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)

    call setup_enefunc_morph_in(morph_in, enefunc)

    if (ctrl_data%rpath_info%rpathmode == RpathmodeMFEP) then
      call setup_fitting_atdyn(.true., ctrl_data%fit_info, ctrl_data%sel_info, &
                               molecule, enefunc)

    else
      call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info,&
                               molecule, enefunc)

    end if

    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_eef1(eef1)

    if (ctrl_data%rpath_info%rpathmode == RpathmodeMEP) then
      ! setup minimize
      !
      call setup_minimize(ctrl_data%min_info, ctrl_data%out_info, &
           ctrl_data%sel_info, molecule, enefunc%qmmm, boundary, minimize)

    else if (ctrl_data%rpath_info%rpathmode == RpathmodeMFEP  .or. &
             ctrl_data%rpath_info%rpathmode == RpathmodeMEPMD .or. &
             ctrl_data%rpath_info%rpathmode == RpathmodeFEP) then
      ! setup ensemble
      !
      call setup_ensemble(ctrl_data%ens_info, boundary, enefunc, ensemble)

    end if

    ! setup constraints
    !
    call setup_constraints(ctrl_data%cons_info, ctrl_data%sel_info, &
                           dynamics, boundary, molecule, &
                           enefunc, constraints)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ctrl_data%ene_info, molecule, constraints, &
                                enefunc)

    ! setup rpath
    !
    call setup_rpath(ctrl_data%rpath_info, ctrl_data%sel_info, rst, rstmep, &
                     dynamics, constraints, molecule, restraints, enefunc,  &
                     minimize, dynvars, rpath)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_rst_all(rst)

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    call dealloc_restraints_all(restraints)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    ! set output
    !
    call setup_output_rpath(ctrl_data%out_info, dynamics, rpath, output)

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)') &
        'Setup_Atdyn_Rpath> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif


    return

  end subroutine setup_atdyn_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_vib
  !> @brief        setup variables and structures in vibration
  !! @authors      KY
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[out]   output    : information of output
  !! @param[out]   molecule  : information of molecules
  !! @param[out]   enefunc   : information of energy function
  !! @param[out]   pairlist  : information of nonbonded pairlist
  !! @param[out]   dynvars   : information of dynamic variables
  !! @param[out]   vibration : information of vibration
  !! @param[out]   boundary  : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_vib(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, vibration, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_vibration),       intent(inout) :: vibration
    type(s_boundary),        intent(inout) :: boundary

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


    ! input (use input sequence for minimize)
    !
    call input_min(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   mode, eef1, table, morph_in, spot)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    call dealloc_spot(spot)
    
    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars)

    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! setup qmmm
    !
    call setup_qmmm(ctrl_data%qmmm_info, ctrl_data%sel_info, boundary, &
                    psf, rst, ctrl_data%ene_info%forcefield, molecule, &
                    enefunc%qmmm)
    call dealloc_psf_all(psf)
    call dealloc_rst_all(rst)    

    ! define enefunc
    !
    call define_enefunc(ctrl_data%ene_info, &
                        boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)

    call setup_enefunc_morph_in(morph_in, enefunc)

    call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             molecule, enefunc)

    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_eef1(eef1)
    call dealloc_restraints_all(restraints)

    ! setup vibration
    !
    call setup_vibration(ctrl_data%inp_info, ctrl_data%vib_info, &
                         ctrl_data%sel_info, molecule, dynvars%coord, &
                         enefunc%qmmm, vibration)

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, dynvars%coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! set output
    !
    call setup_output_vib(ctrl_data%out_info, vibration, output)

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)') &
        'Setup_Atdyn_Vib> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_atdyn_vib

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_atdyn_morph
  !> @brief        setup variables and structures in morphing
  !! @authors      CK
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[out]   output    : information of output
  !! @param[out]   molecule  : information of molecules
  !! @param[out]   enefunc   : information of energy function
  !! @param[out]   pairlist  : information of nonbonded pairlist
  !! @param[out]   dynvars   : information of dynamic variables
  !! @param[out]   morph     : information of morph
  !! @param[out]   boundary  : information of boundary condition
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_atdyn_morph(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, morph, boundary)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_morph),           intent(inout) :: morph
    type(s_boundary),        intent(inout) :: boundary

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

    ! input morph
    !
    call input_morph(ctrl_data%inp_info, top, par, gpr, psf, prmtop, grotop, &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, mode, &
                   eef1, table, morph_in)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)
    call dealloc_top_all(top)
    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)
    call dealloc_mode(mode)

    ! setup boundary conditions
    !
    call setup_boundary(ctrl_data%bound_info, ctrl_data%ene_info%table, &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%sel_info, &
                        molecule, rst, spot, boundary)
    
    if (ctrl_data%ene_info%forcefield == ForcefieldRESIDCG) then
      call update_boundary_cg(ctrl_data%ene_info%cg_pairlistdist_ele, &
          ctrl_data%ene_info%cg_pairlistdist_126, &
          ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
          ctrl_data%ene_info%cg_pairlistdist_DNAbp, &
          ctrl_data%ene_info%cg_pairlistdist_exv, &
          boundary)
    else
      call update_boundary(ctrl_data%ene_info%table, &
          ctrl_data%ene_info%pairlistdist, &
          boundary)
    end if

    ! setup dynamic variables
    !
    call setup_dynvars(molecule, rst, dynvars)
    call dealloc_rst_all(rst)    

    ! setup restraints
    !
    call setup_restraints(ctrl_data%res_info, ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup energy
    !
    call setup_energy(restraints, ctrl_data%ene_info, dynvars%energy)

    ! do not qmmm in MORPH
    if (ctrl_data%qmmm_info%do_qmmm) &
      call error_msg('Setup_Atdyn_Moprh> QM/MM is not allowed in Moprh')
    enefunc%qmmm%do_qmmm=.false.

    call dealloc_psf_all(psf)

    ! define enefunc
    !
    enefunc%morph_flag=.true.
    call define_enefunc(ctrl_data%ene_info, boundary, &
                        par, gpr, prmtop, grotop, eef1, table, &
                        molecule, restraints, enefunc)
    call dealloc_par_all(par)
    call dealloc_gpr_all(gpr)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! set parameters for fitting
    !
    call setup_fitting_atdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             molecule, enefunc)

    ! setup morph
    !
    call setup_morph(ctrl_data%morph_info, morph_in, morph, dynvars, &
    molecule, enefunc, restraints)

    call dealloc_morph_in(morph_in)

    call dealloc_restraints_all(restraints)

    if (enefunc%num_morph_bb+enefunc%num_morph_sc <= 0) then
      call error_msg('Setup_Atdyn_Morph> no morphing distance')
    endif

    ! setup pairlist
    !
    call setup_pairlist(enefunc, boundary, molecule%atom_coord, dynvars%trans, &
                        dynvars%coord_pbc, pairlist)

    ! set output
    !
    call setup_output_morph(ctrl_data%out_info, morph, output)

    dynvars%verbose = morph%verbose

    if (enefunc%contact_check) &
      call error_msg('Setup_atdyn_morph> contact_check is not allowed')

    return

  end subroutine setup_atdyn_morph

end module at_setup_atdyn_mod
