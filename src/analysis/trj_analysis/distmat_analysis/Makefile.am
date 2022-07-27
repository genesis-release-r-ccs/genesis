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

bin_PROGRAMS = distmat_analysis

distmat_analysis_SOURCES = \
	da_option_str.fpp \
	da_option.fpp \
	da_control.fpp \
	da_setup.fpp \
	da_analyze.fpp \
	da_main.fpp

EXTRA_DIST = \
	Makefile.depends 

distmat_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
distmat_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(distmat_analysis_OBJECTS) *.f90 *~ *.mod

.fpp.o:
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) \
		-I../../../lib -I../../libana $(FCFLAGS) -c $*.f90
endif

distmat_analysis$(EXEEXT) : $(distmat_analysis_OBJECTS) \
		$(distmat_analysis_ADDLIBANA) $(distmat_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(distmat_analysis_OBJECTS) \
		$(distmat_analysis_ADDLIBANA) $(distmat_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
