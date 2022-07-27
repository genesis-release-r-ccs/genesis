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

noinst_LIBRARIES = libspana.a
libspana_a_SOURCES = \
        sa_setup_mpi.fpp \
        sa_option_str.fpp \
        sa_option.fpp \
        sa_ensemble_str.fpp \
        sa_ensemble.fpp \
        sa_boundary_str.fpp \
        sa_boundary.fpp \
        sa_domain_str.fpp\
        sa_domain.fpp \
        sa_tool.fpp \
        sa_control.fpp \
        sa_timers.fpp \
        sa_setup.fpp

EXTRA_DIST = \
	Makefile.depends 

libspana_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
libspana_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

libspana.a : makelib $(libspana_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(libspana_a_OBJECTS)

clean:
	rm -f $(noinst_LIBRARIES) $(libspana_a_OBJECTS) *.f90 *~ *.mod

makelib:
	cd ../../../lib; make
	cd ../../libana; make

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
