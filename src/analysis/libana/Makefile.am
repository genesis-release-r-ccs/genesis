#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

SUBDIRS = .

noinst_LIBRARIES = libana.a

include $(top_srcdir)/src/mk/common-libonly.am

libana_a_SOURCES = \
	input_str.fpp \
	output_str.fpp \
	trajectory_str.fpp \
	measure.fpp \
	fileio_trj.fpp \
	fileio_atp.fpp \
	fileio_rtp.fpp \
	fileio_mt.fpp \
	input.fpp \
	output.fpp \
	trajectory.fpp \
	pbc_correct.fpp \
	analysis_timers.fpp

EXTRA_DIST = \
	Makefile.depends 

libana_ADDLIB = $(top_builddir)/src/lib/lib.a

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

libana.a : $(libana_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(libana_a_OBJECTS)

clean:
	rm -f $(noinst_LIBRARIES) $(libana_a_OBJECTS) *.f90 *~ *.mod

depend: clean_depend
	python $(top_srcdir)/fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
