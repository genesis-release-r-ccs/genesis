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
include $(top_srcdir)/src/mk/common-analysis.am

bin_PROGRAMS = rpath_generator

rpath_generator_SOURCES = \
        rg_option_str.fpp \
        rg_option.fpp \
        rg_control.fpp \
        rg_setup.fpp \
        rg_generate.fpp \
        rg_main.fpp

EXTRA_DIST = \
	Makefile.depends 

rpath_generator_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
rpath_generator_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(rpath_generator_OBJECTS) *.f90 *~ *.mod

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

rpath_generator$(EXEEXT) : $(rpath_generator_OBJECTS) \
		$(rpath_generator_ADDLIBANA) $(rpath_generator_ADDLIB)
	$(FC) -o $(bin_PROGRAMS) $(rpath_generator_OBJECTS) \
		$(rpath_generator_ADDLIBANA) $(rpath_generator_ADDLIB) $(LDFLAGS) \
		$(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
