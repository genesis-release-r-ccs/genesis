#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : DM
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

if INTERNAL_LAPACK
  INT_LAPACK_LIB=$(LAPACK_LIB3)
else
  INT_LAPACK_LIB=
endif
include $(top_srcdir)/src/mk/common-spana.am

bin_PROGRAMS = hbond_analysis

hbond_analysis_SOURCES = \
	hbond_option_str.fpp \
	hbond_option.fpp \
	hbond_control.fpp \
	hbond_analyze.fpp\
	hbond_main.fpp

EXTRA_DIST = \
	Makefile.depends 

hbond_analysis_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
hbond_analysis_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a
hbond_analysis_ADDLIBSPANA = $(top_builddir)/src/analysis/sp_analysis/libspana/libspana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends

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


hbond_analysis$(EXEEXT) : makelib $(hbond_analysis_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(hbond_analysis_OBJECTS) \
		 $(hbond_analysis_ADDLIBSPANA) $(hbond_analysis_ADDLIBANA) \
		 $(hbond_analysis_ADDLIB) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(bin_PROGRAMS) $(hbond_analysis_OBJECTS) *.f90 *~ *.mod

makelib:
	$(MAKE) -C $(top_srcdir)/src/lib
	$(MAKE) -C $(top_srcdir)/src/analysis/libana
	$(MAKE) -C $(top_srcdir)/src/analysis/sp_analysis/libspana

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
