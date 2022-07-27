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

SUBDIRS := ../lib .

bin_PROGRAMS = atdyn

atdyn_SOURCES = \
	at_energy_str.fpp \
	at_enefunc_str.fpp \
	at_pairlist_str.fpp \
	at_boundary_str.fpp \
	at_constraints_str.fpp \
	at_experiments_str.fpp \
	at_restraints_str.fpp \
	at_ensemble_str.fpp \
	at_dynvars_str.fpp \
	at_dynamics_str.fpp \
	at_minimize_str.fpp \
	at_vibration_str.fpp \
	at_output_str.fpp \
	at_remd_str.fpp \
	at_morph_str.fpp \
	at_rpath_str.fpp \
	at_qmmm.fpp \
	at_experiments.fpp \
	at_pairlist.fpp \
	at_boundary.fpp \
	at_energy_morph.fpp \
	at_energy_bonds.fpp \
	at_energy_angles.fpp \
	at_energy_dihedrals.fpp \
	at_energy_cg_nonlocal.fpp \
	at_energy_table_cubic.fpp \
	at_energy_table_linear.fpp \
	at_energy_table_linear_bondcorr.fpp \
	at_energy_pme.fpp \
	at_energy_eef1.fpp \
	at_energy_gbsa.fpp \
	at_energy_nonbonds.fpp \
	at_energy_restraints.fpp \
	at_energy_go.fpp \
	at_energy_soft.fpp \
	at_energy_gamd.fpp \
	at_energy.fpp \
	at_enefunc_restraints.fpp \
	at_enefunc_gamd.fpp \
	at_enefunc_gbsa.fpp \
	at_enefunc_table.fpp \
	at_enefunc_pme.fpp \
	at_enefunc_charmm.fpp \
	at_enefunc_amber.fpp \
	at_enefunc_gromacs.fpp \
	at_enefunc_go.fpp \
	at_enefunc.fpp \
	at_constraints.fpp \
	at_restraints.fpp \
	at_ensemble.fpp \
	at_dynvars.fpp \
	at_input.fpp \
	at_output.fpp \
	at_md_leapfrog.fpp \
	at_md_vverlet.fpp \
	at_md_vverlet_cg.fpp \
	at_minimize.fpp \
	at_vibration.fpp \
	at_remd.fpp \
	at_rpath_mep.fpp \
	at_rpath_fep.fpp \
	at_rpath.fpp \
	at_morph.fpp \
	at_gamd.fpp \
	at_dynamics.fpp \
	at_control.fpp \
	at_setup_atdyn.fpp \
	at_setup_mpi.fpp \
	atdyn.fpp

at_energy_linear_object = \
	at_energy_table_linear.o \
	at_energy_table_linear_bondcorr.o 

if QSIMULATE
  atdyn_ADDLIB = ../lib/lib.a ../lib/qsimulate/libqs.a ../lib/json-fortran/libjsonfortran.a $(INT_LAPACK_LIB)
  INCLUDE_DIR = -I../lib -I../lib/json-fortran
else
  atdyn_ADDLIB = ../lib/lib.a $(INT_LAPACK_LIB)
  INCLUDE_DIR = -I../lib
endif

EXTRA_DIST = \
	Makefile.depends 

SUFFIXES = .fpp
DEPENDS = Makefile.depends

$(at_energy_linear_object) : %.o : %.fpp
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) $(INCLUDE_DIR) $(FCFLAGS_FPP__PGI_NOOPT) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) $(INCLUDE_DIR) $(FCFLAGS_FPP__PGI_NOOPT) -c $*.f90
endif

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) $(INCLUDE_DIR) $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) $(INCLUDE_DIR) $(FCFLAGS) -c $*.f90
endif

atdyn$(EXEEXT) : $(atdyn_ADDLIB) $(atdyn_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(atdyn_OBJECTS) $(atdyn_ADDLIB) $(atdyn_LDADD) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(atdyn_OBJECTS) *.f90 *~ *.mod

depend: clean_depend
	python ../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
