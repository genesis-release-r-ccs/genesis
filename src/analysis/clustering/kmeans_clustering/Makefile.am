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

bin_PROGRAMS = kmeans_clustering

kmeans_clustering_SOURCES = \
        kc_option_str.fpp \
        kc_option.fpp \
        kc_control.fpp \
        kc_setup.fpp \
        kc_analyze.fpp \
        kc_main.fpp

EXTRA_DIST = \
	Makefile.depends 

kmeans_clustering_ADDLIB = ../../../lib/lib.a $(INT_LAPACK_LIB)
kmeans_clustering_ADDLIBANA = ../../libana/libana.a

SUFFIXES = .fpp
DEPENDS = Makefile.depends
MAINTAINERCLEANFILES = Makefile.in

CLEANFILES = $(bin_PROGRAMS) $(kmeans_clustering_OBJECTS) *.f90 *~ *.mod

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

kmeans_clustering$(EXEEXT) : $(kmeans_clustering_ADDLIBANA) $(kmeans_clustering_ADDLIB) \
							$(kmeans_clustering_OBJECTS)
	$(FC) -o $(bin_PROGRAMS) $(kmeans_clustering_OBJECTS) \
		$(kmeans_clustering_ADDLIBANA) $(kmeans_clustering_ADDLIB) $(LDFLAGS)

depend: clean_depend
	python ../../../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
