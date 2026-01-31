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

SUBDIRS = .
include $(top_srcdir)/src/mk/common-analysis.am

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

rst_upgrade_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
rst_upgrade_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends

.NOTPARALLEL:
spdyndir=$(top_srcdir)/src/spdyn/

sp_energy_str.fpp:
	ln -fs $(spdyndir)/sp_energy_str.fpp .

sp_domain_str.fpp:
	ln -fs $(spdyndir)/sp_domain_str.fpp .

sp_enefunc_str.fpp:
	ln -fs $(spdyndir)/sp_enefunc_str.fpp .

sp_pairlist_str.fpp:
	ln -fs $(spdyndir)/sp_pairlist_str.fpp .

sp_boundary_str.fpp:
	ln -fs $(spdyndir)/sp_boundary_str.fpp .

sp_constraints_str.fpp:
	ln -fs $(spdyndir)/sp_constraints_str.fpp .

sp_restraints_str.fpp:
	ln -fs $(spdyndir)/sp_restraints_str.fpp .

sp_ensemble_str.fpp:
	ln -fs $(spdyndir)/sp_ensemble_str.fpp .

sp_dynvars_str.fpp:
	ln -fs $(spdyndir)/sp_dynvars_str.fpp .

sp_dynamics_str.fpp:
	ln -fs $(spdyndir)/sp_dynamics_str.fpp .

sp_remd_str.fpp:
	ln -fs $(spdyndir)/sp_remd_str.fpp .

sp_parallel_io.fpp:
	ln -fs $(spdyndir)/sp_parallel_io.fpp .

FPP2F90_NVHPC += -I$(spdyndir)
FPP2F90       += -I$(spdyndir)
F90COMPILE    += -I$(spdyndir)

.fpp.o:
if USE_NVHPC
	$(FPP2F90_NVHPC) $*.fpp > $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP2F90) $*.f90
	mv $*.cpp.f90 $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
	$(FPP2F90) $*.fpp $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
endif
endif

rst_upgrade$(EXEEXT) : $(rst_upgrade_ADDLIBANA) $(rst_upgrade_ADDLIB) \
						$(rst_upgrade_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(rst_upgrade_OBJECTS) \
		$(rst_upgrade_ADDLIBANA) $(rst_upgrade_ADDLIB) $(LDFLAGS) \
		$(LIBS)

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
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
