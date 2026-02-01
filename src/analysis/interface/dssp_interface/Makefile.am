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

bin_PROGRAMS = dssp_interface

dssp_interface_SOURCES = \
        di_option_str.fpp \
        di_option.fpp \
        di_control.fpp \
        di_setup.fpp \
        di_analyze.fpp \
        di_main.fpp

EXTRA_DIST = \
	Makefile.depends 

dssp_interface_ADDLIB = $(top_builddir)/src/lib/lib.a $(INT_LAPACK_LIB)
dssp_interface_ADDLIBANA = $(top_builddir)/src/analysis/libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(dssp_interface_OBJECTS) *.f90 *~ *.mod

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

dssp_interface$(EXEEXT) : $(dssp_interface_ADDLIBANA) $(dssp_interface_ADDLIB) \
							$(dssp_interface_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(dssp_interface_OBJECTS) \
		$(dssp_interface_ADDLIBANA) $(dssp_interface_ADDLIB) $(LDFLAGS) \
		$(LIBS)

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
