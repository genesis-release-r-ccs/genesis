!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_setup_spdyn
!> @brief   setup variables and structures in MD (DD) simulaton
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_setup_spdyn_mod

  use sp_control_mod
  use sp_restart_mod
  use sp_output_mod
  use sp_input_mod
  use sp_minimize_mod
  use sp_dynamics_mod
  use sp_dynvars_mod
  use sp_domain_mod
  use sp_ensemble_mod
  use sp_enefunc_table_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_boundary_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_energy_mod
  use sp_energy_pme_mod
  use sp_parallel_io_mod
  use sp_update_domain_mod
  use sp_communicate_mod
  use sp_output_str_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_enefunc_fit_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use sp_remd_str_mod
  use sp_rpath_str_mod
  use sp_remd_mod
  use sp_rpath_mod
  use sp_experiments_mod
  use sp_gamd_mod
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
  use structure_check_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_fep_topology_mod
  use sp_alchemy_mod
  use sp_alchemy_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_spdyn_md
  public  :: setup_spdyn_min
  public  :: setup_spdyn_remd
  public  :: setup_spdyn_rpath
  public  :: setup_spdyn_md_pio
  public  :: setup_spdyn_min_pio
  private :: read_parallel_io_rst

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md
  !> @brief        setup variables and structures in MD simulation
  !! @authors      JJ, HO
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_alchemy),optional,intent(inout) :: alchemy

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
    type(s_localres)         :: localres
    logical                  :: use_parallel_io

    ! setup parallel I/O
    !
    use_parallel_io = pio_check_ranked_file(ctrl_data%inp_info%rstfile)
    if (use_parallel_io) then

      call setup_spdyn_md_pio(ctrl_data, output, enefunc, pairlist,     &
                              dynvars, dynamics, constraints, ensemble, &
                              restraints, boundary, domain, comm)
      return

    end if

    ! read input files
    !
    call input_md(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    if (present(alchemy)) then
      ! FEP: define singleA, singleB, dualA, dualB, and preserved regions
      call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
        ctrl_data%alch_info)
    end if

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
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        rst, boundary)
  
    ! set parameters for domain 
    !
    if (present(alchemy)) then
      ! FEP
      call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)
    else
      call setup_domain(ctrl_data%ene_info,  &
                        ctrl_data%cons_info, &
                        boundary, molecule, enefunc, constraints, domain)
    end if

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc, comm)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)

    if (domain%fep_use) then
      ! FEP: setup alchemy
      call setup_alchemy_md(ctrl_data%alch_info, domain, enefunc, alchemy)
    end if

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        if (domain%fep_use) call error_msg('VDWPME is not allowed in FEP')
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        ctrl_data%alch_info,  &
                        molecule, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)


    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    if (domain%fep_use) then
      ! FEP: Some thermostat and barostat are not available in FEP
      if (ensemble%tpcontrol  == TpcontrolBerendsen) then
        call error_msg('Setup_Ensemble> Berendsen is not allowed in FEP')
      else if (ensemble%tpcontrol  == TpcontrolNoseHoover) then
        call error_msg('Setup_Ensemble> NoseHoover is not allowed in FEP')
      else if (ensemble%tpcontrol  == TpcontrolMTK) then
        call error_msg('Setup_Ensemble> MTK is not allowed in FEP')
      end if
    end if

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, domain, &
                           enefunc, constraints)

    if (domain%fep_use) then
      ! FEP: remove degree of freedom of singleB
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)


    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc)


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
        'Setup_Spdyn_Md> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    ! turn off group T/P when vacuum is used
    !
    if (ensemble%group_tp) then
      if (enefunc%vacuum) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
            'Setup_Spdyn_Md> group_tp = No is assigned with vacuum = yes'
        ensemble%group_tp = .false.
      end if
    end if

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Md> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    return

  end subroutine setup_spdyn_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_min
  !> @brief        setup variables and structures in minimization
  !! @authors      TM, JJ, HO 
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   minimize    : information of minimize
  !! @param[out]   constraints : information of constraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_min(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, minimize, constraints, boundary,       &
                             domain, comm, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_alchemy),optional,intent(inout) :: alchemy

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
    type(s_localres)         :: localres
 

    ! setup parallel I/O
    !

    if (pio_check_ranked_file(ctrl_data%inp_info%rstfile)) then

      call setup_spdyn_min_pio(ctrl_data, output, enefunc, pairlist,       &
                               dynvars, minimize, constraints, restraints, &
                               boundary, domain, comm)
      return

    end if


    ! read input files
    !
    call input_min(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    if (present(alchemy)) then
      ! FEP: define singleA, singleB, dualA, dualB, and preserved regions
      call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
        ctrl_data%alch_info)
    end if

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

    ! setup structure_check
    !
    call setup_structure_check(molecule,                                      &
                     ctrl_data%ene_info%forcefield_char, molecule%atom_coord, &
                     ctrl_data%min_info%check_structure,                      &
                     ctrl_data%min_info%fix_ring_error,                       &
                     ctrl_data%min_info%fix_chirality_error,                  &
                     ctrl_data%min_info%exclude_ring_grpid,                   &
                     ctrl_data%min_info%exclude_chiral_grpid)

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        rst, boundary)


    ! set parameters for domain
    !
    if (ctrl_data%cons_info%rigid_bond) then
      ctrl_data%cons_info%rigid_bond = .false.
      if (main_rank) then
        write(Msgout, '(A)') 'Setup_Constraints> WARNING : &
                            & constraints are applied only for water'
        write(Msgout, '(A)') 
      end if
    end if
      
    if (present(alchemy)) then
      ! FEP
      call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)
    else
      call setup_domain(ctrl_data%ene_info,  &
                        ctrl_data%cons_info, &
                        boundary, molecule, enefunc, constraints, domain)
    end if


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc, comm)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call setup_constraints(ctrl_data%cons_info, par, prmtop, &
                           grotop, molecule, domain,         &
                           enefunc, constraints)
    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    if (domain%fep_use) then
      ! FEP: remove degree of freedom of singleB
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_localres(localres, LocalRestraint)
    

    if (domain%fep_use) then
      ! FEP: setup alchemy
      call setup_alchemy_min(ctrl_data%alch_info, domain, enefunc, alchemy)
    end if
    
    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        if (domain%fep_use) call error_msg('VDWPME is not allowed in FEP')
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if

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
        'Setup_Spdyn_Min> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    return

  end subroutine setup_spdyn_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_remd
  !> @brief        setup variables and structures in REMD simulation
  !! @authors      TM, HO
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   remd        : information of remd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_remd(ctrl_data, output, molecule, enefunc, pairlist,  &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, remd, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd
    type(s_alchemy),optional,intent(inout) :: alchemy

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
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
    logical                  :: use_parallel_io


    ! setup parallel I/O
    !
    use_parallel_io = pio_check_ranked_file(ctrl_data%inp_info%rstfile)
    if (use_parallel_io) then
      call setup_spdyn_remd_pio(ctrl_data, output, enefunc, pairlist,     &
                                dynvars, dynamics, constraints, ensemble, &
                                restraints, boundary, domain, comm, remd)
      return
    end if
   
    ! read input files
    !
    call input_remd(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                    pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                    localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    if (present(alchemy)) then
      ! FEP: define singleA, singleB, dualA, dualB, and preserved regions
      call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
        ctrl_data%alch_info)
    end if

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
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        rst, boundary)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    call setup_solute_tempering(ctrl_data%rep_info, &
                                molecule, restraints, ctrl_data%cons_info)

    ! set parameters for domain 
    !
    if (present(alchemy)) then
      ! FEP
      call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)
    else
      call setup_domain(ctrl_data%ene_info,  &
                        ctrl_data%cons_info, &
                        boundary, molecule, enefunc, constraints, domain)
    end if

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc, comm)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    if (domain%fep_use) then
      ! FEP: setup alchemy
      call setup_alchemy_remd(ctrl_data%alch_info, domain, enefunc, alchemy)
    end if

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        if (domain%fep_use) call error_msg('VDWPME is not allowed in FEP')
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if

    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        ctrl_data%alch_info,  &
                        molecule, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    if (domain%fep_use) then
      ! FEP: Some thermostat and barostat are not available in FEP
      if (ensemble%tpcontrol  == TpcontrolBerendsen) then
        call error_msg('Setup_Ensemble> Berendsen is not allowed in FEP')
      else if (ensemble%tpcontrol  == TpcontrolNoseHoover) then
        call error_msg('Setup_Ensemble> NoseHoover is not allowed in FEP')
      else if (ensemble%tpcontrol  == TpcontrolMTK) then
        call error_msg('Setup_Ensemble> MTK is not allowed in FEP')
      end if
    end if

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, domain, &
                           enefunc, constraints)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup remd
    !
    call setup_remd(ctrl_data%rep_info, rst, boundary, dynamics, molecule, &
                    domain, restraints, ensemble, enefunc, remd, alchemy)

    if (domain%fep_use) then
      ! FEP:remove degree of freedom of singleB
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)

    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc, remd)

    ! set output
    !
    call setup_output_remd(ctrl_data%out_info, dynamics, remd, output)

    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    ! turn off group T/P when vacuum is used
    !
    if (ensemble%group_tp) then
      if (enefunc%vacuum) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
            'Setup_Spdyn_Md> group_tp = No is assigned with vacuum = yes'
        ensemble%group_tp = .false.
      end if
    end if

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Remd> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Remd> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    return

  end subroutine setup_spdyn_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_rpath
  !> @brief        setup variables and structures in RPATH simulation
  !! @authors      YK, YM
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   rpath       : information of rpath
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_rpath(ctrl_data, output, molecule, enefunc, pairlist, &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, rpath)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
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


    ! define replica
    !
    call define_nreplica(ctrl_data%rpath_info, rpath)

    ! read input files
    !
    call input_rpath(ctrl_data%inp_info, top, par, psf, prmtop, grotop, &
                    pdb, crd, ambcrd, grocrd, rst, ref, fit, ambref, groref,&
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
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        rst, boundary)

    ! set parameters for domain 
    !
    call setup_domain(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc, comm)

    call setup_fitting_spdyn(.true., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if

    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        ctrl_data%alch_info,  &
                        molecule, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, domain, &
                           enefunc, constraints)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup rpath
    !
    call setup_rpath(ctrl_data%rpath_info, rst, boundary, dynamics, molecule, &
                    domain, restraints, ensemble, enefunc, rpath)
    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)

    ! set output
    !
    call setup_output_rpath(ctrl_data%out_info, dynamics, rpath, output)

    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Rpath> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Rpath> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    return

  end subroutine setup_spdyn_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md_pio
  !> @brief        setup variables and structures in MD simulation
  !! @authors      JJ
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   restraints  : information of restraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md_pio(ctrl_data, output, enefunc, pairlist,     &
                                dynvars, dynamics, constraints, ensemble, &
                                restraints, boundary, domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_restraints),      intent(inout) :: restraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_localres)         :: localres


    ! set parameters for boundary condition
    !
    call setup_boundary_pio(ctrl_data%bound_info, boundary)

    ! read domain restart file (Parallel I/O)
    !
    call read_parallel_io_rst(ctrl_data,   &
                              .false.,     &
                              domain,      &
                              enefunc,     &
                              boundary,    &
                              restraints,  &
                              constraints, &
                              dynvars,     &
                              dynamics)

    ! determine the number of degree of freedom
    !
    domain%num_deg_freedom = 3*domain%num_atom_all

    !
    ! determine restart or not

    dynamics%restart=pio_restart

    ! read input files
    !
    if (ctrl_data%inp_info%local_resfile /= '') &
    call input_localres(ctrl_data%inp_info%local_resfile, localres)

    ! set parameters for domain 
    !
    call setup_domain_pio(ctrl_data%ene_info, ctrl_data%cons_info, &
                          boundary, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! update migration
    !
    call domain_update_pio(domain, enefunc, boundary, constraints, comm)

    ! set parameters for restraints
    ! 
    call setup_restraints_pio(ctrl_data%res_info, &
                              ctrl_data%sel_info, &
                              restraints, enefunc)

    ! setup enefunc in each domain
    !
    call define_enefunc_pio(ctrl_data%ene_info, localres, comm, &
                            constraints, restraints, domain, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics_pio(ctrl_data%dyn_info,   &
                            ctrl_data%bound_info, &
                            ctrl_data%res_info,   &
                            domain, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)


    ! set parameters for constraints
    !
    call setup_constraints_pio(ctrl_data%cons_info, pio_restart, &
                               enefunc, constraints, domain)


    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)

    ! change step number
    !
    dynvars%step = 0

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md_Pio> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Md_Pio> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    return

  end subroutine setup_spdyn_md_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_remd_pio
  !> @brief        setup variables and structures in REMD simulation
  !! @authors      JJ
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   restraints  : information of restraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   remd        : information of REMD
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_remd_pio(ctrl_data, output, enefunc, pairlist,   &
                                dynvars, dynamics, constraints, ensemble, &
                                restraints, boundary, domain, comm, remd)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_restraints),      intent(inout) :: restraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd

    ! local variables
    type(s_localres)         :: localres

    ! set parameters for boundary condition
    !
    call setup_boundary_pio(ctrl_data%bound_info, boundary)

    ! read domain restart file (Parallel I/O)
    !
    call read_parallel_io_rst(ctrl_data,   &
                              .true.,      &
                              domain,      &
                              enefunc,     &
                              boundary,    &
                              restraints,  &
                              constraints, &
                              dynvars,     &
                              dynamics,    &
                              remd)

    ! determine the number of degree of freedom
    !
    domain%num_deg_freedom = 3*domain%num_atom_all

    !
    ! determine restart or not

    dynamics%restart=pio_restart

    ! read input files
    !
    if (ctrl_data%inp_info%local_resfile /= '') &
    call input_localres(ctrl_data%inp_info%local_resfile, localres)

    ! set parameters for domain 
    !
    call setup_domain_pio(ctrl_data%ene_info, ctrl_data%cons_info, &
                          boundary, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! update migration
    !
    call domain_update_pio(domain, enefunc, boundary, constraints, comm)

    ! set parameters for restraints
    ! 
    call setup_restraints_pio(ctrl_data%res_info, &
                              ctrl_data%sel_info, &
                              restraints, enefunc)

    call setup_solute_tempering_pio(ctrl_data%rep_info, restraints, &
                                    ctrl_data%cons_info)

    ! setup enefunc in each domain
    !
    call define_enefunc_pio(ctrl_data%ene_info, localres, comm, &
                            constraints, restraints, domain, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics_pio(ctrl_data%dyn_info,   &
                            ctrl_data%bound_info, &
                            ctrl_data%res_info,   &
                            domain, dynamics)

    ! mass repartitioning
    !
    if (dynamics%hydrogen_mr) &
      call setup_mass_repartitioning(dynamics, constraints, domain, &
                                     enefunc)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)


    ! set parameters for constraints
    !
    call setup_constraints_pio(ctrl_data%cons_info, pio_restart, &
                               enefunc, constraints, domain)

    ! setup remd
    !
    call setup_remd_pio(ctrl_data%rep_info, boundary, dynamics, domain, &
                        restraints, ensemble, enefunc, remd)
    call dealloc_restraints(restraints, RestraintsList)

    ! set output
    !
    call setup_output_remd(ctrl_data%out_info, dynamics, remd, output)

    ! change step number
    !
    dynvars%step = 0

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md_Pio> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    ! pressure degree of freedom for group pressure
    !
    if (ensemble%group_tp) then
      if (constraints%rigid_bond) then
        call compute_group_deg_freedom(domain, constraints)
      else
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Setup_Spdyn_Remd_Pio> group_tp = No is assigned with rigid_bond = No'
        ensemble%group_tp = .false.
      end if
    end if

    return

  end subroutine setup_spdyn_remd_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_min_pio
  !> @brief        setup variables and structures in minimization
  !! @authors      JJ
  !! @param[inout] ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   minimize    : information of minimize
  !! @param[out]   constraints : information of constraints
  !! @param[out]   restraints  : information of restraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_min_pio(ctrl_data, output, enefunc, pairlist,       &
                                 dynvars, minimize, constraints, restraints, &
                                 boundary, domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(inout) :: restraints 
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_localres)         :: localres
 

    ! set parameters for boundary condition
    !
    call setup_boundary_pio(ctrl_data%bound_info, boundary)

 
    ! read domain restart file (Parallel I/O)
    !
    call read_parallel_io_rst(ctrl_data,   &
                              .false.,     &
                              domain,      &
                              enefunc,     &
                              boundary,    &
                              restraints,  &
                              constraints)

    ! read input files
    !
    if (ctrl_data%inp_info%local_resfile /= '') &
    call input_localres(ctrl_data%inp_info%local_resfile, localres)


    ! set parameters for domain 
    !
    call setup_domain_pio(ctrl_data%ene_info, ctrl_data%cons_info, &
                          boundary, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! update migration
    !
    call domain_update_pio(domain, enefunc, boundary, constraints, comm)

    ! set parameters for restraints
    ! 
    call setup_restraints_pio(ctrl_data%res_info, &
                              ctrl_data%sel_info, &
                              restraints, enefunc)

    ! setup enefunc in each domain
    !
    call define_enefunc_pio(ctrl_data%ene_info, localres, comm, &
                            constraints, restraints, domain, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      if (enefunc%vdw == VDWPME) then
        call setup_pme (domain, boundary, enefunc)
        call pme_pre_lj(domain, boundary)
      else
        call setup_pme(domain, boundary, enefunc)
        call pme_pre  (domain, boundary)
      end if
    end if


    ! set parameters for minimize
    !
    call setup_minimize(ctrl_data%min_info, minimize)


    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)

    ! change step number
    !
    dynvars%step = 0

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Min_Pio> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    end if

    return

  end subroutine setup_spdyn_min_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_parallel_io_rst
  !> @brief        read parallel I/O restart file
  !! @authors      JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[in]    use_replica : flag for use replica or not
  !! @param[out]   domain      : information of each domain
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   restraints  : information of restraints
  !! @param[out]   constraints : information of constraints
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   remd        : information of REMD
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_parallel_io_rst(ctrl_data, use_replica, domain, enefunc,    &
                                  boundary, restraints, constraints, dynvars, &
                                  dynamics, remd)

    ! formal arguments
    type(s_ctrl_data),             intent(in)    :: ctrl_data
    logical,                       intent(in)    :: use_replica
    type(s_domain),                intent(inout) :: domain
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_boundary),              intent(inout) :: boundary
    type(s_restraints),            intent(inout) :: restraints
    type(s_constraints),           intent(inout) :: constraints
    type(s_dynvars),   optional,   intent(inout) :: dynvars
    type(s_dynamics),  optional,   intent(inout) :: dynamics
    type(s_remd),      optional,   intent(inout) :: remd    

    ! local variables
    integer                        :: file_num, file_tot_num, ix, iy, iz, ip
    integer                        :: i, iproc(3)
    integer                        :: my_x_rank_pio, my_y_rank_pio
    integer                        :: my_z_rank_pio, my_tot_rank_pio
    integer                        :: nplace, num_pio_tot
    logical                        :: fit_check
    character(MaxMultiFilename)    :: filename
    integer,       parameter :: Places (7) = &
         (/10, 100, 1000, 10000, 100000, 1000000, 10000000/)

    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Parallel_Io_Rst> '
      write(MsgOut,'(A)') ' '
    end if

    if (.not. boundary%multiple_file_read) then
      filename = ctrl_data%inp_info%rstfile
      call include_id_to_filename(filename)
      call pio_read_domain_rst( &
                            pio_get_ranked_filename(  &
                              filename, my_rank_pio), &
                            1, 1,                &
                            boundary,            &
                            domain,              &
                            enefunc,             &
                            constraints,         &
                            dynvars,             &
                            dynamics,            &
                            remd)
    else

      iproc(1) = mod(my_country_rank, boundary%num_domain(1))
      iproc(2) = mod(my_country_rank/boundary%num_domain(1),  &
                     boundary%num_domain(2))
      iproc(3) = my_country_rank &
               / (boundary%num_domain(1)*boundary%num_domain(2))

      file_num = 0
      file_tot_num = boundary%multiple_file(1) &
                    *boundary%multiple_file(2) &
                    *boundary%multiple_file(3)

      num_pio_tot = boundary%num_pio_domain(1) * boundary%num_pio_domain(2) &
                    * boundary%num_pio_domain(3)
      do ip = 5, size(Places)
        if (num_pio_tot < Places(ip)) &
          exit
      end do
      nplace  = ip

      do ix = 1, boundary%multiple_file(1)
        my_x_rank_pio = iproc(1)*boundary%multiple_file(1) + ix - 1
        do iy = 1, boundary%multiple_file(2)
          my_y_rank_pio = iproc(2)*boundary%multiple_file(2) + iy - 1
          do iz = 1, boundary%multiple_file(3)
            my_z_rank_pio = iproc(3)*boundary%multiple_file(3) + iz - 1
            file_num = file_num + 1
            my_tot_rank_pio = my_x_rank_pio                   &
                            + my_y_rank_pio * boundary%num_pio_domain(1) &
                            + my_z_rank_pio * boundary%num_pio_domain(1) &
                                            * boundary%num_pio_domain(2)
            filename = ctrl_data%inp_info%rstfile
            call include_id_to_filename(filename)

            call pio_read_domain_rst( &
                                  pio_get_ranked_filename(      &
                                    filename, my_tot_rank_pio, nplace), &
                                  file_num,            &
                                  file_tot_num,        &
                                  boundary,            &
                                  domain,              &
                                  enefunc,             &
                                  constraints,         &
                                  dynvars,             &
                                  dynamics,            &
                                  remd)
          end do
        end do
      end do

    end if

    if (ctrl_data%inp_info%selfile /= '') then
      filename = ctrl_data%inp_info%selfile
      call pio_read_selection(filename, restraints)
    else
      call pio_read_selection_blank(restraints)
    end if

    if (.not. main_rank) &
      return

    write(MsgOut,'(A,A)') ' Parallel I/O file: ', &
                                             trim(ctrl_data%inp_info%rstfile)
    write(MsgOut,'(A,L)') '          Restart : ', pio_restart
    write(MsgOut,'(A)')   ' '

    ! compatible check for fit

    if (ctrl_data%fit_info%fitting_method /= FittingMethodNO) then
      fit_check=.false.
      do i = 1,enefunc%num_restraintfuncs
        ! restaint posi is not fitted in pio
        if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
            enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
            enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
            fit_check = .true.
            exit
        end if
      end do
      if (fit_check) then
        call error_msg('Setup_Spdyn_Md_Pio> WARNING: Fitting is not allowed')
      else
        if (main_rank) then
          write(MsgOut,*) "Setup_Fitting_Spdyn> NO fitting is applied, skip"
          write(MsgOut,*) 
        end if
        enefunc%fitting_method=FittingMethodNO
      end if
    end if

    return

  end subroutine read_parallel_io_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    include_id_to_filename
  !> @brief        include id to filename
  !! @authors      TM
  !! @param[inout] filename : replicate filename
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine include_id_to_filename(filename)

    ! formal arguments
    character(MaxFilename),  intent(inout) :: filename

    ! local variables
    integer                  :: comp, ndigit, id
    integer                  :: i, j, ci1, ci2, cnumber
    character(MaxFilename)   :: filename_ori
    character(10)            :: frmt, cd, cid


    ! define replica id
    !
    id = my_country_no + 1
    do i = 1, 100
      comp = 10**i
      if (id < comp) then
        ndigit = i
        exit
      end if
    end do

    ! check filename
    !
    filename_ori = filename

    ci1 = scan(filename, '{')
    ci2 = scan(filename, '}')

    if (ci1 == 0 .or. ci2 ==0) &
      return

    if (ci1 > 0 .and. ci2 > ci1) then

      write(cd,'(i10)') ndigit
      frmt = '(i' // trim(adjustl(cd)) // '.' // trim(adjustl(cd)) // ')'
      write(cid,frmt) id

      cnumber = len_trim(filename_ori)
      if (cnumber + ndigit > MaxFilename) &
         call error_msg('Error: too long filename'//filename_ori)

      j = 0
      do i = 1, ci1 - 1
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do
      do i = 1, ndigit
        j = j + 1
        filename(j:j) = cid(i:i)
      end do
      do i = ci2+1, MaxFilename
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do

    end if

    return

  end subroutine include_id_to_filename

end module sp_setup_spdyn_mod
