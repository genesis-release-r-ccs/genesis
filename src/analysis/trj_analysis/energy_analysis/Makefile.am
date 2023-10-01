#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : MK
#
#  (c) Copyright 2016-2018 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8
if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana ../../../atdyn .

bin_PROGRAMS = energy_analysis

energy_analysis_SOURCES = \
			ea_option_str.fpp \
			ea_option.fpp \
			ea_control.fpp \
			ea_setup.fpp \
			ea_analyze.fpp \
			ea_main.fpp

EXTRA_DIST =\
			Makefile.depends

if QSIMULATE
energy_analysis_ADDLIB = ../../../lib/lib.a ../../../lib/qsimulate/libqs.a ../../../lib/json-fortran/libjsonfortran.a
else
energy_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
endif

energy_analysis_ADDLIBANA = ../../libana/libana.a
optional_base = at_input.o at_boundary.o at_boundary_str.o \
	at_dynvars.o at_dynvars_str.o at_restraints.o \
	at_restraints_str.o at_remd.o at_remd_str.o \
	at_energy.o at_energy_angles.o at_energy_bonds.o \
	at_energy_dihedrals.o at_energy_go.o at_energy_nonbonds.o \
	at_energy_pme.o at_energy_eef1.o at_energy_gbsa.o at_energy_restraints.o \
	at_energy_morph.o at_energy_soft.o at_energy_cg_nonlocal.o \
	at_energy_str.o at_energy_table_cubic.o \
	at_energy_gamd.o \
	at_energy_table_linear.o at_energy_table_linear_bondcorr.o \
	at_enefunc.o at_enefunc_amber.o at_enefunc_charmm.o \
	at_enefunc_go.o at_enefunc_gromacs.o at_enefunc_pme.o \
	at_enefunc_restraints.o at_enefunc_str.o at_enefunc_gbsa.o at_enefunc_table.o \
	at_enefunc_gamd.o \
	at_experiments_str.o at_experiments.o \
	at_pairlist.o at_pairlist_str.o at_setup_mpi.o \
	at_constraints.o at_constraints_str.o at_output.o \
	at_md_leapfrog.o at_md_vverlet.o at_md_vverlet_cg.o at_qmmm.o \
	at_ensemble.o at_ensemble_str.o \
	at_dynamics.o at_dynamics_str.o \
	at_morph.o at_morph_str.o
optional_objs = $(foreach file, $(optional_base), ../../../atdyn/$(file))

SUFFIXES = .fpp
DEPENDS	= Makefile.depends

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I ../../../atdyn $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana -I ../../../atdyn $(FCFLAGS) -c $*.f90
endif

energy_analysis$(EXEEXT) : $(energy_analysis_OBJECTS) \
		$(energy_analysis_ADDLIBANA) $(energy_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(energy_analysis_OBJECTS) \
		$(optional_objs) \
		$(energy_analysis_ADDLIBANA) $(energy_analysis_ADDLIB) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(energy_analysis_OBJECTS) *.f90 *~ *.mod

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
