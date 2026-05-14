!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_setup_cgdyn
!> @brief   Sets up variables and structures for MD (DD) simulations
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_setup_cgdyn_mod

  use cg_control_mod
  use cg_restart_mod
  use cg_output_mod
  use cg_input_mod
  use cg_minimize_mod
  use cg_dynamics_mod
  use cg_dynvars_mod
  use cg_domain_mod
  use cg_ensemble_mod
  use cg_restraints_mod
  use cg_boundary_mod
  use cg_pairlist_mod
  use cg_enefunc_mod
  use cg_enefunc_fit_mod
  use cg_energy_mod
  use cg_communicate_str_mod
  use cg_communicate_mod
  use cg_output_str_mod
  use cg_minimize_str_mod
  use cg_dynamics_str_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_restraints_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_energy_str_mod
  use cg_domain_str_mod
  use cg_migration_mod
  use molecules_mod
  use molecules_str_mod
  use fitting_mod
  use fitting_str_mod
  use fileio_localres_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_crd_mod
  use fileio_rst_mod
  use fileio_mode_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_cgdyn_md
  public  :: setup_cgdyn_min

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cgdyn_md
  !> @brief        Sets up variables and structures for MD simulation
  !! @authors      JJ
  !! @param[in]    ctrl_data   : control parameters information
  !! @param[out]   output      : output information
  !! @param[out]   molecule    : molecules information
  !! @param[out]   enefunc     : energy function information
  !! @param[out]   pairlist    : nonbonded pairlist information
  !! @param[out]   dynvars     : dynamic variables information
  !! @param[out]   dynamics    : molecular dynamics information
  !! @param[out]   ensemble    : ensemble information
  !! @param[out]   boundary    : boundary condition information
  !! @param[out]   domain      : domain information
  !! @param[out]   comm        : domain communicator
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cgdyn_md(ctrl_data, output, grotop, molecule, enefunc,      &
                            pairlist, dynvars, dynamics, ensemble, boundary,   &
                            domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_grotop),          intent(inout) :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
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
    type(s_localres)         :: localres


    ! read input files
    !
    call input_md(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if


    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,                      &
                        ctrl_data%ene_info,                        &
                        ctrl_data%ens_info%ensemble,               &
                        ctrl_data%ene_info%cg_pairlistdist_ele,    &
                        ctrl_data%ene_info%cg_pairlistdist_126,    &
                        ctrl_data%ene_info%cg_pairlistdist_WF,     &
                        ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
                        ctrl_data%ene_info%cg_pairlistdist_DNAbp,  &
                        ctrl_data%ene_info%cg_pairlistdist_exv,    &
                        ctrl_data%ene_info%cg_cellsize_min,        &
                        molecule, rst, boundary)

    ! set parameters for domain 
    !
    call setup_domain(ctrl_data%ene_info,  &
                      boundary, molecule, enefunc, domain)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup enefunc in each domain
    !
    call dealloc_molecules_bond(molecule)

    call define_enefunc(ctrl_data%ene_info, grotop, &
                        molecule, restraints,       &
                        domain, enefunc, comm)

    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info, molecule, dynamics)
    if (ctrl_data%dyn_info%lbupdate_period == 0) then
      call dealloc_molecules(molecule, MoleculeAtom)
      call dealloc_par_all(par)
      call dealloc_prmtop_all(prmtop)
      call dealloc_grotop_all(grotop)
    end if

    call setup_fitting_cgdyn(ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)
 
    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call update_communicate_size(domain, comm)
    call update_enefunc_contact(domain, enefunc)
    call communicate_contact(domain, comm, enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(boundary, enefunc, domain, pairlist)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, enefunc, ensemble)

    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)


    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if


    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_cgdyn_Md> nonb_limiter : minimum distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    domain%num_deg_freedom = molecule%num_deg_freedom

    return

  end subroutine setup_cgdyn_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_cgdyn_min
  !> @brief        Sets up variables and structures for minimization
  !! @authors      TM, JJ
  !! @param[in]    ctrl_data   : control parameters information
  !! @param[out]   output      : output information
  !! @param[out]   molecule    : molecules information
  !! @param[out]   enefunc     : energy function information
  !! @param[out]   pairlist    : nonbonded pairlist information
  !! @param[out]   dynvars     : dynamic variables information
  !! @param[out]   minimize    : minimize information
  !! @param[out]   boundary    : boundary condition information
  !! @param[out]   domain      : domain information
  !! @param[out]   comm        : domain communicator
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_cgdyn_min(ctrl_data, output, grotop, molecule, enefunc,     &
                             pairlist, dynvars, minimize, boundary, domain,    &
                             comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_grotop),          intent(inout) :: grotop
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
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
    type(s_localres)         :: localres
 
 
    ! read input files
    !
    call input_min(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if


    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,                      &
                        ctrl_data%ene_info,                        &
                        ctrl_data%ens_info%ensemble,               &
                        ctrl_data%ene_info%cg_pairlistdist_ele,    &
                        ctrl_data%ene_info%cg_pairlistdist_126,    &
                        ctrl_data%ene_info%cg_pairlistdist_WF,     &
                        ctrl_data%ene_info%cg_pairlistdist_PWMcos, &
                        ctrl_data%ene_info%cg_pairlistdist_DNAbp,  &
                        ctrl_data%ene_info%cg_pairlistdist_exv,    &
                        ctrl_data%ene_info%cg_cellsize_min,        &
                        molecule, rst, boundary)

    call setup_domain(ctrl_data%ene_info,  &
                      boundary, molecule, enefunc, domain)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup enefunc in each domain
    !
    call dealloc_molecules_bond(molecule)

    call define_enefunc(ctrl_data%ene_info, grotop, &
                        molecule, restraints,       &
                        domain, enefunc, comm)

    call setup_fitting_cgdyn(ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call dealloc_molecules(molecule, MoleculeAtom)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call update_communicate_size(domain, comm)
    call update_enefunc_contact(domain, enefunc)
    call communicate_contact(domain, comm, enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(boundary, enefunc, domain, pairlist)

    ! set parameters for minimize
    !
    call setup_minimize(ctrl_data%min_info, minimize)


    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars)

    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)


    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_cgdyn_Min> nonb_limiter : minimum distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_cgdyn_min

end module cg_setup_cgdyn_mod
