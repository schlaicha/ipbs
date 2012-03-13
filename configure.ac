# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/ipbs.cc])
AM_CONFIG_HEADER([config.h])
#GSL_VERSION=1.00 # don#t know if this makes sense here...

AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR(could not find required version of GSLcblas)])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR(could not find required version of GSL)])

# Check for all DUNE dependencies
DUNE_CHECK_ALL

# Check for IPBS specific dependencies
IPBS_CHECKS

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
AC_SUBST([LDADD], '$(top_builddir)/dune/ipbs/libipbs.la $(DUNE_LIBS)')
LIBS="$LIBS $DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/ipbs/Makefile
  dune/ipbs/test/Makefile
  m4/Makefile
  ipbs.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
