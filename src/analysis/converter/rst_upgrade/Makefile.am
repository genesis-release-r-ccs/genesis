#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8
if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif

SUBDIRS := ../../libana .

bin_PROGRAMS = rst_upgrade

rst_upgrade_SOURCES = \
	sp_energy_str.fpp \
	sp_domain_str.fpp \
	sp_enefunc_str.fpp \
	sp_pairlist_str.fpp \
	sp_boundary_str.fpp \
	sp_constraints_str.fpp \
	sp_restraints_str.fpp \
	sp_ensemble_str.fpp \
	sp_dynvars_str.fpp \
	sp_dynamics_str.fpp \
	sp_remd_str.fpp \
	sp_parallel_io.fpp \
	fileio_rst_old.fpp \
	sp_parallel_io_old.fpp \
	ru_control.fpp \
	ru_setup.fpp \
	ru_upgrade.fpp \
	ru_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rst_upgrade_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
rst_upgrade_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends

.NOTPARALLEL:

sp_energy_str.fpp:
	ln -fs ../../../spdyn/sp_energy_str.fpp .

sp_domain_str.fpp:
	ln -fs ../../../spdyn/sp_domain_str.fpp .

sp_enefunc_str.fpp:
	ln -fs ../../../spdyn/sp_enefunc_str.fpp .

sp_pairlist_str.fpp:
	ln -fs ../../../spdyn/sp_pairlist_str.fpp .

sp_boundary_str.fpp:
	ln -fs ../../../spdyn/sp_boundary_str.fpp .

sp_constraints_str.fpp:
	ln -fs ../../../spdyn/sp_constraints_str.fpp .

sp_restraints_str.fpp:
	ln -fs ../../../spdyn/sp_restraints_str.fpp .

sp_ensemble_str.fpp:
	ln -fs ../../../spdyn/sp_ensemble_str.fpp .

sp_dynvars_str.fpp:
	ln -fs ../../../spdyn/sp_dynvars_str.fpp .

sp_dynamics_str.fpp:
	ln -fs ../../../spdyn/sp_dynamics_str.fpp .

sp_remd_str.fpp:
	ln -fs ../../../spdyn/sp_remd_str.fpp .

sp_parallel_io.fpp:
	ln -fs ../../../spdyn/sp_parallel_io.fpp .

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) -I../../../spdyn $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) -I../../../spdyn $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

rst_upgrade$(EXEEXT) : $(rst_upgrade_ADDLIBANA) $(rst_upgrade_ADDLIB) \
						$(rst_upgrade_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(rst_upgrade_OBJECTS) \
		$(rst_upgrade_ADDLIBANA) $(rst_upgrade_ADDLIB) $(LDFLAGS)

clean:
	rm -f $(bin_PROGRAMS) $(rst_upgrade_OBJECTS) *.f90 *~ *.mod \
		sp_energy_str.fpp \
		sp_domain_str.fpp \
		sp_enefunc_str.fpp \
		sp_pairlist_str.fpp \
		sp_boundary_str.fpp \
		sp_constraints_str.fpp \
		sp_restraints_str.fpp \
		sp_ensemble_str.fpp \
		sp_dynvars_str.fpp \
		sp_dynamics_str.fpp \
		sp_remd_str.fpp \
		sp_parallel_io.fpp

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
