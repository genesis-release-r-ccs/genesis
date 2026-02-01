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

bin_PROGRAMS = pcrd_convert

pcrd_convert_SOURCES = \
	parallel_trj.fpp \
	pc_option_str.fpp \
	pc_option.fpp \
	pc_control.fpp \
	pc_setup.fpp \
	pc_convert.fpp \
	pc_main.fpp

EXTRA_DIST = \
	Makefile.depends 

pcrd_convert_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
pcrd_convert_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(pcrd_convert_OBJECTS) *.f90 *~ *.mod

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

pcrd_convert$(EXEEXT) :	$(pcrd_convert_ADDLIBANA) $(pcrd_convert_ADDLIB) \
						$(pcrd_convert_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(pcrd_convert_OBJECTS) \
		$(pcrd_convert_ADDLIBANA) $(pcrd_convert_ADDLIB) $(LDFLAGS) $(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
