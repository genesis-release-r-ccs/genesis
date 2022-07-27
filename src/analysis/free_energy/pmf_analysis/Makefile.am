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

bin_PROGRAMS = pmf_analysis

pmf_analysis_SOURCES = \
        pm_option_str.fpp \
        pm_option.fpp \
        pm_control.fpp \
        pm_setup.fpp \
        pm_analyze.fpp \
        pm_main.fpp

EXTRA_DIST = \
	Makefile.depends 

pmf_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
pmf_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(pmf_analysis_OBJECTS) *.f90 *~ *.mod

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

pmf_analysis$(EXEEXT) : $(pmf_analysis_ADDLIBANA) $(pmf_analysis_ADDLIB) \
						$(pmf_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(pmf_analysis_OBJECTS) \
		$(pmf_analysis_ADDLIBANA) $(pmf_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
