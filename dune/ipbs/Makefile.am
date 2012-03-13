
SUBDIRS = test
noinst_LTLIBRARIES = libipbs.la

# the ipbs special stuff
libipbs_la_SOURCES = boundary.cc parser.cc sysparams.cc
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libipbs_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libipbs_la_LIBADD = $(DUNE_LIBS)

ipbsincludedir = $(includedir)/dune/ipbs
ipbsinclude_HEADERS = ipbs.hh \
               boundaries.hh \
               ipbsolver.hh \
               p0layout.hh \
               datawriter.hh \
               ipbs_Pk.hh \
	       PBLocalOperator.hh \
				 			 sysparams.hh \
							 parser.hh \
							 boundaries.hh

include $(top_srcdir)/am/global-rules
