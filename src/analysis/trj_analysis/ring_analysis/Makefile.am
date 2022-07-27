#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : CK
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

bin_PROGRAMS = ring_analysis

ring_analysis_SOURCES = \
	ra_option_str.fpp \
	ra_option.fpp \
	ra_control.fpp \
	ra_setup.fpp \
	ra_analyze.fpp \
	ra_main.fpp

EXTRA_DIST = \
	Makefile.depends 

ring_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
ring_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(ring_analysis_OBJECTS) *.f90 *~ *.mod

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

ring_analysis$(EXEEXT) : $(ring_analysis_OBJECTS) \
			$(ring_analysis_ADDLIBANA) $(ring_analysis_ADDANA)
	$(FC) -o $(bin_PROGRAMS) $(ring_analysis_OBJECTS) \
		$(ring_analysis_ADDLIBANA) $(ring_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
