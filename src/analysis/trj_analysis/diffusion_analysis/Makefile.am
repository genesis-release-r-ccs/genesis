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

bin_PROGRAMS = diffusion_analysis
include $(top_srcdir)/src/mk/common-analysis.am

diffusion_analysis_SOURCES = \
  	da_option_str.fpp \
        da_option.fpp \
        da_control.fpp \
  	da_setup.fpp \
        da_analyze.fpp \
        da_main.fpp

EXTRA_DIST = \
	Makefile.depends 

diffusion_analysis_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
diffusion_analysis_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(diffusion_analysis_OBJECTS) *.f90 *~ *.mod

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

diffusion_analysis$(EXEEXT) : $(diffusion_analysis_OBJECTS) \
		$(diffusion_analysis_ADDLIBANA) $(diffusion_analysis_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(diffusion_analysis_OBJECTS) \
		$(diffusion_analysis_ADDLIBANA) $(diffusion_analysis_ADDLIB) \
		$(LDFLAGS) $(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
