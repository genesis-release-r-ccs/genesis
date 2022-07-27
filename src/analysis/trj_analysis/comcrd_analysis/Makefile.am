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

bin_PROGRAMS = comcrd_analysis

comcrd_analysis_SOURCES = \
	ca_option_str.fpp \
	ca_option.fpp \
	ca_control.fpp \
	ca_setup.fpp \
	ca_analyze.fpp \
	ca_main.fpp

EXTRA_DIST = \
	Makefile.depends 

comcrd_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
comcrd_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(comcrd_analysis_OBJECTS) *.f90 *~ *.mod

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

comcrd_analysis$(EXEEXT) : $(comcrd_analysis_OBJECTS) \
				$(comcrd_analysis_ADDLIBANA) $(comcrd_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(comcrd_analysis_OBJECTS) \
		$(comcrd_analysis_ADDLIBANA) $(comcrd_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
