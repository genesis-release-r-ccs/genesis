#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB1)
else
  INT_LAPACK_LIB=
endif

SUBDIRS	:= ../lib .

if PARALLEL_IO
  PRST_SETUP=prst_setup
else
  PRST_SETUP=
endif

bin_PROGRAMS = spdyn $(PRST_SETUP)

spdyn_SOURCES = \
	sp_domain_str.fpp \
	sp_energy_str.fpp \
	sp_enefunc_str.fpp \
	sp_pairlist_str.fpp \
	sp_boundary_str.fpp \
	sp_constraints_str.fpp \
	sp_experiments_str.fpp \
	sp_restraints_str.fpp \
	sp_ensemble_str.fpp \
	sp_dynvars_str.fpp \
	sp_dynamics_str.fpp \
	sp_minimize_str.fpp \
	sp_alchemy_str.fpp \
	sp_output_str.fpp \
	sp_remd_str.fpp \
	sp_rpath_str.fpp \
	sp_migration.fpp \
	sp_communicate.fpp \
	sp_parallel_io.fpp \
	sp_pairlist_generic.fpp \
	sp_pairlist_fugaku.fpp \
	sp_pairlist_gpu.fpp \
	sp_pairlist.fpp \
	sp_boundary.fpp \
	sp_experiments.fpp \
	sp_energy_efield.fpp \
	sp_energy_bonds.fpp \
	sp_energy_angles.fpp \
	sp_energy_dihedrals.fpp \
	sp_energy_nonbond_generic.fpp \
	sp_energy_nonbond_fugaku.fpp \
	sp_energy_nonbond_intel.fpp \
	sp_energy_nonbond_gpu.fpp \
	sp_energy_nonbond_fep.fpp \
	sp_energy_nonbond_vacuum.fpp \
	sp_energy_table_linear.fpp \
	sp_energy_14_notable.fpp \
	sp_energy_nonbond_pme_table.fpp \
	sp_energy_nonbond_pme_notable.fpp \
	sp_energy_nonbond_ljpme.fpp \
	sp_energy_table_linear_bondcorr.fpp \
	sp_energy_table_cubic.fpp \
	sp_energy_pme_opt_1dalltoall.fpp \
	sp_energy_pme_noopt_1dalltoall.fpp \
	sp_energy_pme_opt_2dalltoall.fpp \
	sp_energy_pme_noopt_2dalltoall.fpp \
	sp_energy_pme_opt_slab.fpp \
	sp_energy_pme_noopt_slab.fpp \
	sp_energy_pme.fpp \
	sp_energy_nonbonds.fpp \
	sp_enefunc_fit.fpp \
	sp_energy_restraints.fpp \
	sp_energy_gamd.fpp \
	sp_energy.fpp \
	sp_grest_energy.fpp \
	sp_fep_energy.fpp \
	sp_enefunc_table.fpp \
	sp_enefunc_restraints.fpp \
	sp_enefunc_localres.fpp \
	sp_enefunc_charmm.fpp \
	sp_enefunc_amber.fpp \
	sp_enefunc_gromacs.fpp \
	sp_enefunc_gamd.fpp \
	sp_enefunc.fpp \
	sp_constraints.fpp \
	sp_restraints.fpp \
	sp_ensemble.fpp \
	sp_domain.fpp \
	sp_dynvars.fpp \
	sp_assign_velocity.fpp \
	sp_update_domain.fpp \
	sp_input.fpp \
	sp_output.fpp \
	sp_md_leapfrog.fpp \
	sp_md_vverlet.fpp \
	sp_md_respa.fpp \
	sp_md_mts.fpp \
	sp_alchemy.fpp \
	sp_remd.fpp \
    sp_fep_topology.fpp \
	sp_rpath.fpp \
	sp_gamd.fpp \
	sp_dynamics.fpp \
	sp_minimize.fpp \
	sp_restart.fpp \
	sp_control.fpp \
	sp_setup_spdyn.fpp \
	sp_setup_mpi.fpp 

if CUDAGPU
spdyn_SOURCES += \
	gpu_sp_energy.cu \
	gpu_common.h \
	spdyn.fpp
else
spdyn_SOURCES += \
	spdyn.fpp
endif

spdyn_ADDLIB = ../lib/lib.a $(INT_LAPACK_LIB)

prst_setup_SOURCES = \
	sp_domain_str.fpp \
	sp_energy_str.fpp \
	sp_enefunc_str.fpp \
	sp_pairlist_str.fpp \
	sp_boundary_str.fpp \
	sp_constraints_str.fpp \
	sp_experiments_str.fpp \
	sp_restraints_str.fpp \
	sp_ensemble_str.fpp \
	sp_dynvars_str.fpp \
	sp_dynamics_str.fpp \
	sp_alchemy_str.fpp \
	sp_minimize_str.fpp \
	sp_output_str.fpp \
	sp_remd_str.fpp \
	sp_rpath_str.fpp \
	sp_migration.fpp \
	sp_communicate.fpp \
	sp_parallel_io.fpp \
	sp_pairlist_generic.fpp \
	sp_pairlist_fugaku.fpp \
	sp_pairlist_gpu.fpp \
	sp_pairlist.fpp \
	sp_boundary.fpp \
	sp_experiments.fpp \
	sp_energy_efield.fpp \
	sp_energy_bonds.fpp \
	sp_energy_angles.fpp \
	sp_energy_dihedrals.fpp \
	sp_energy_nonbond_generic.fpp \
	sp_energy_nonbond_fugaku.fpp \
	sp_energy_nonbond_intel.fpp \
	sp_energy_nonbond_gpu.fpp \
	sp_energy_nonbond_fep.fpp \
	sp_energy_nonbond_vacuum.fpp \
	sp_energy_table_linear.fpp \
	sp_energy_14_notable.fpp \
	sp_energy_nonbond_pme_table.fpp \
	sp_energy_nonbond_pme_notable.fpp \
	sp_energy_nonbond_ljpme.fpp \
	sp_energy_table_linear_bondcorr.fpp \
	sp_energy_table_cubic.fpp \
	sp_energy_pme_opt_1dalltoall.fpp \
	sp_energy_pme_noopt_1dalltoall.fpp \
	sp_energy_pme_opt_2dalltoall.fpp \
	sp_energy_pme_noopt_2dalltoall.fpp \
	sp_energy_pme_opt_slab.fpp \
	sp_energy_pme_noopt_slab.fpp \
	sp_energy_pme.fpp \
	sp_energy_nonbonds.fpp \
	sp_enefunc_fit.fpp \
	sp_energy_restraints.fpp \
	sp_energy_gamd.fpp \
	sp_energy.fpp \
	sp_grest_energy.fpp \
	sp_fep_energy.fpp \
	sp_enefunc_table.fpp \
	sp_enefunc_restraints.fpp \
	sp_enefunc_localres.fpp \
	sp_enefunc_charmm.fpp \
	sp_enefunc_amber.fpp \
	sp_enefunc_gromacs.fpp \
	sp_enefunc_gamd.fpp \
	sp_enefunc.fpp \
	sp_constraints.fpp \
	sp_restraints.fpp \
	sp_ensemble.fpp \
	sp_domain.fpp \
	sp_dynvars.fpp \
	sp_assign_velocity.fpp \
	sp_update_domain.fpp \
	sp_input.fpp \
	sp_output.fpp \
	sp_md_leapfrog.fpp \
	sp_md_vverlet.fpp \
	sp_md_respa.fpp \
	sp_md_mts.fpp \
	sp_alchemy.fpp \
	sp_remd.fpp \
    sp_fep_topology.fpp \
	sp_rpath.fpp \
	sp_gamd.fpp \
	sp_dynamics.fpp \
        pr_fileio.c \
        pr_huge_molecule.fpp \
        pr_pdb2hm.fpp \
        pr_charmm2hm.fpp \
        pr_amber2hm.fpp \
        pr_gromacs2hm.fpp \
        pr_domain_index.fpp \
        pr_control.fpp \
        pr_select_atoms_peer.fpp \
        pr_setup_spdyn_peer.fpp \
        pr_run_prst_setup.fpp

if CUDAGPU
prst_setup_SOURCES += \
	gpu_sp_energy.cu \
	prst_setup.fpp
else
prst_setup_SOURCES += \
	prst_setup.fpp
endif

prst_setup_ADDLIB = ../lib/lib.a $(INT_LAPACK_LIB)

SUFFIXES = .fpp

EXTRA_DIST = \
	Makefile.depends \
	pr_huge_molecule_file.fpp \
	pr_huge_molecule_mem.fpp \
	pr_domain_index_mem.fpp \
	pr_domain_index_file.fpp

DEPENDS	= Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(spdyn_OBJECTS) $(prst_setup_OBJECTS) *.f90 *~ *.mod

sp_pairlist_fugaku.o : sp_pairlist_fugaku.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_KERNELS_PAIRLIST) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
endif

sp_energy_dihedrals.o : sp_energy_dihedrals.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_DIHEDRAL) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
endif

sp_energy_nonbond_fugaku.o : sp_energy_nonbond_fugaku.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_KERNELS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
endif

if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_KERNELS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
endif

sp_energy_pme_noopt_1dalltoall.o : sp_energy_pme_noopt_1dalltoall.fpp 
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_NOOPT) -c $*.f90
endif

sp_energy_pme_noopt_2dalltoall.o : sp_energy_pme_noopt_2dalltoall.fpp 
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_NOOPT) -c $*.f90
endif

sp_energy_pme_noopt_slab.o : sp_energy_pme_noopt_slab.fpp 
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS_NOOPT) -c $*.f90
endif

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../lib $(FCFLAGS) -c $*.f90
endif

if CUDAGPU
.cu.o:
	$(NVCC) $(NVCCFLAG) --ptxas-options=-v -lineinfo $*.cu
endif

spdyn$(EXEEXT) : $(spdyn_ADDLIB) $(spdyn_OBJECTS)
	$(FC) -o spdyn $(spdyn_OBJECTS) \
		$(spdyn_ADDLIB) $(LDFLAGS)

prst_setup$(EXEEXT) : $(spdyn_ADDLIB) $(prst_setup_OBJECTS)
	$(FC) -o prst_setup $(prst_setup_OBJECTS) \
		$(prst_setup_ADDLIB) $(LDFLAGS)

$(spdyn_ADDLIB):
	cd ../lib; ${MAKE}

depend: clean_depend
	python ../../fortdep.py -e ../lib -f ".fpp" *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
