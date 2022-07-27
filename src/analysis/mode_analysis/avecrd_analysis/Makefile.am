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

bin_PROGRAMS = avecrd_analysis

avecrd_analysis_SOURCES = \
	aa_option_str.fpp \
	aa_option.fpp \
	aa_control.fpp \
	aa_setup.fpp \
	aa_analyze.fpp \
	aa_main.fpp

EXTRA_DIST = \
	Makefile.depends 

avecrd_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
avecrd_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(avecrd_analysis_OBJECTS) *.f90 *~ *.mod

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

avecrd_analysis$(EXEEXT) : $(avecrd_analysis_ADDLIBANA) $(avecrd_analysis_ADDLIB) \
							$(avecrd_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(avecrd_analysis_OBJECTS) \
		$(avecrd_analysis_ADDLIBANA) $(avecrd_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
