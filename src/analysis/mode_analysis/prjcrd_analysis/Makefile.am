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

bin_PROGRAMS = prjcrd_analysis

prjcrd_analysis_SOURCES = \
        pa_option_str.fpp \
        pa_option.fpp \
        pa_control.fpp \
        pa_setup.fpp \
        pa_analyze.fpp \
        pa_main.fpp

EXTRA_DIST = \
	Makefile.depends 

prjcrd_analysis_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
prjcrd_analysis_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(prjcrd_analysis_OBJECTS) *.f90 *~ *.mod

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

prjcrd_analysis$(EXEEXT) : $(prjcrd_analysis_OBJECTS) \
			$(prjcrd_analysis_ADDLIBANA) $(prjcrd_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(prjcrd_analysis_OBJECTS) \
		$(prjcrd_analysis_ADDLIBANA) $(prjcrd_analysis_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
