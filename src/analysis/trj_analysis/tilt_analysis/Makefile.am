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

bin_PROGRAMS = tilt_analysis

tilt_analysis_SOURCES = \
	ta_option_str.fpp \
	ta_option.fpp \
	ta_control.fpp \
	ta_setup.fpp \
	ta_analyze.fpp \
	ta_main.fpp

EXTRA_DIST = \
	Makefile.depends 

tilt_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
tilt_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(tilt_analysis_OBJECTS) *.f90 *~ *.mod

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

tilt_analysis$(EXEEXT) : $(tilt_analysis_OBJECTS) \
		$(tilt_analysis_ADDLIBANA) $(tilt_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(tilt_analysis_OBJECTS) \
		$(tilt_analysis_ADDLIBANA) $(tilt_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
