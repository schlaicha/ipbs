

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
NORMALTESTS = test_surfacepot 
# list of tests to run
TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

#
## common flags
#

# paranoia
DUNE_EXTRA_CHECKS = -DDUNE_DEVEL_MODE
# output coverage
#COVERAGE = -fprofile-arcs -ftest-coverage
AM_CXXFLAGS = $(COVERAGE)
AM_CPPFLAGS = @AM_CPPFLAGS@ $(DUNE_EXTRA_CHECKS)

# The directory where the example grid files used for testing are located
AM_CPPFLAGS+=-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"$(top_srcdir)/doc/grids/\"

#
## define the programs
#
sphere2d.msh:
	gmsh -2 sphere2d.geo

test_surfacepot_SOURCES = test_surfacepot.cc \
						  ../sysparams.cc ../parser.cc ../boundary.cc
test_surfacepot_CPPFLAGS = $(AM_CPPFLAGS) \
		$(DUNEMPICPPFLAGS) \
		$(UG_CPPFLAGS) \
		$(AMIRAMESH_CPPFLAGS) \
		$(ALBERTA_CPPFLAGS) \
		$(ALUGRID_CPPFLAGS) \
		$(GRIDDIM_CPPFLAGS) \
		-DSURFACE_POT
test_surfacepot_LDADD = \
	 	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
		$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
		$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
		$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
		$(UG_LDFLAGS) $(UG_LIBS) \
		$(DUNEMPILIBS)\
		$(LDADD) 
test_surfacepot_LDFLAGS = $(AM_LDFLAGS) \
	   	$(DUNEMPILDFLAGS) \
		$(UG_LDFLAGS) \
		$(AMIRAMESH_LDFLAGS) \
		$(ALBERTA_LDFLAGS) \
		$(ALUGRID_LDFLAGS) \
		$(DUNE_LDFLAGS)
test_surfacepot_DEPENDENCIES = sphere2d.msh

# distribution tarball
# SOURCES = parser.cc 
# gridcheck not used explicitly, we should still ship it :)
# EXTRA_DIST = $(SOURCES)

CLEANFILES = *.dat *.vtu *.msh
include $(top_srcdir)/am/global-rules
