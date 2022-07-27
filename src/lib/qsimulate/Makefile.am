noinst_LIBRARIES = libqs.a

libqs_a_SOURCES = interface.cc

EXTRA_DIST = \
        interface.hpp

SUFFIXES = .cc
.cc.o:
	$(CXX) $(CXXFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) -I. -c $*.cc

libqs.a : $(libqs_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(libqs_a_OBJECTS)

clean:
	rm -f $(noinst_LIBRARIES) $(libqs_a_OBJECTS) *~
