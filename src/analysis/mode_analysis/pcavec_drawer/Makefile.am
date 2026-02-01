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

bin_PROGRAMS = pcavec_drawer
include $(top_srcdir)/src/mk/common-analysis.am

pcavec_drawer_SOURCES = \
        pd_option_str.fpp \
        pd_option.fpp \
        pd_control.fpp \
        pd_setup.fpp \
        pd_draw.fpp \
        pd_main.fpp

EXTRA_DIST = \
	Makefile.depends 

pcavec_drawer_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
pcavec_drawer_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(pcavec_drawer_OBJECTS) *.f90 *~ *.mod

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


pcavec_drawer$(EXEEXT) : $(pcavec_drawer_ADDLIBANA) $(pcavec_drawer_ADDLIB) \
							$(pcavec_drawer_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(pcavec_drawer_OBJECTS) \
		$(pcavec_drawer_ADDLIBANA) $(pcavec_drawer_ADDLIB) $(LDFLAGS) \
		$(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
