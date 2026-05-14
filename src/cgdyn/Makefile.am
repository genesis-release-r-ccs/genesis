#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS	= .

bin_PROGRAMS = cgdyn

include $(top_srcdir)/src/mk/common-libonly.am

cgdyn_SOURCES = \
	cg_domain_str.fpp \
	cg_energy_str.fpp \
	cg_enefunc_str.fpp \
	cg_pairlist_str.fpp \
	cg_boundary_str.fpp \
	cg_restraints_str.fpp \
	cg_ensemble_str.fpp \
	cg_dynvars_str.fpp \
	cg_dynamics_str.fpp \
	cg_minimize_str.fpp \
	cg_output_str.fpp \
	cg_communicate_str.fpp \
	cg_migration.fpp \
	cg_communicate.fpp \
	cg_pairlist.fpp \
	cg_boundary.fpp \
	cg_energy_bonds.fpp \
	cg_energy_angles.fpp \
	cg_energy_dihedrals.fpp \
	cg_energy_bases.fpp \
	cg_energy_go.fpp \
	cg_energy_nonlocal.fpp \
	cg_energy_martini.fpp \
	cg_enefunc_fit.fpp \
	cg_energy_restraints.fpp \
	cg_energy.fpp \
	cg_enefunc_restraints.fpp \
	cg_enefunc_localres.fpp \
	cg_enefunc_gromacs.fpp \
	cg_enefunc.fpp \
	cg_restraints.fpp \
	cg_ensemble.fpp \
	cg_domain.fpp \
	cg_dynvars.fpp \
	cg_assign_velocity.fpp \
	cg_update_domain.fpp \
	cg_input.fpp \
	cg_output.fpp \
	cg_md_leapfrog.fpp \
	cg_md_vverlet.fpp \
	cg_dynamics.fpp \
	cg_minimize.fpp \
	cg_restart.fpp \
	cg_control.fpp \
	cg_setup_cgdyn.fpp \
	cg_setup_mpi.fpp 

cgdyn_SOURCES += \
	cgdyn.fpp

cgdyn_ADDLIB = $(top_builddir)/src/lib/lib.a

EXTRA_DIST = \
	Makefile.depends 

SUFFIXES = .fpp

DEPENDS	= Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(cgdyn_OBJECTS)  *.f90 *~ *.mod

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

if CUDAGPU
.cu.o:
	$(NVCC) $(NVCCFLAG) --ptxas-options=-v -lineinfo $*.cu
endif

cgdyn$(EXEEXT) : $(cgdyn_ADDLIB) $(cgdyn_OBJECTS)
	$(FC) -o cgdyn $(cgdyn_OBJECTS) \
		$(cgdyn_ADDLIB) $(LDFLAGS)


depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
