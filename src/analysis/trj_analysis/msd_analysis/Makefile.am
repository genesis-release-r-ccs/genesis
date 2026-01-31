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

bin_PROGRAMS = msd_analysis
include $(top_srcdir)/src/mk/common-analysis.am

msd_analysis_SOURCES = \
	one_molecule_str.fpp \
	molecule_manipulate.fpp \
	molecule_select.fpp \
	ma_option_str.fpp \
	ma_option.fpp \
	ma_control.fpp \
	ma_setup.fpp \
	ma_analyze.fpp \
	ma_main.fpp

EXTRA_DIST = \
	Makefile.depends 

msd_analysis_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
msd_analysis_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(msd_analysis_OBJECTS) *.f90 *~ *.mod

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

msd_analysis$(EXEEXT) : $(msd_analysis_OBJECTS) \
		$(msd_analysis_ADDLIBANA) $(msd_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(msd_analysis_OBJECTS) \
		$(msd_analysis_ADDLIBANA) $(msd_analysis_ADDLIB) $(LDFLAGS) \
		$(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
