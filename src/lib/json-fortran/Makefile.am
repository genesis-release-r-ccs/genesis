#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : KY
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

noinst_LIBRARIES = libjsonfortran.a

libjsonfortran_a_SOURCES = \
        json_kinds.F90 \
        json_parameters.F90 \
        json_string_utilities.F90 \
        json_value_module.F90 \
        json_file_module.F90 \
        json_module.F90 

EXTRA_DIST = \
	Makefile.depends 

SUFFIXES = .F90
DEPENDS	= Makefile.depends

CLEANFILES = $(noinst_LIBRARIES) $(libjsonfortran_a_OBJECTS) *~ *.mod

.F90.o:
	$(FC)  $(DEFAULT_INCLUDES) $(INCLUDES) -I. $(FCFLAGS) -c -o $@ $*.F90 

libjsonfortran.a : $(libjsonfortran_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(libjsonfortran_a_OBJECTS)

depend: clean_depend
	python ../../../fortdep.py *F90 > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
