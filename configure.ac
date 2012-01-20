# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/iPBS.cc])
AM_CONFIG_HEADER([config.h])
GSL_VERSION=1.00 # don#t know if this makes sense here...

AC_ARG_ENABLE(grid,
              [--enable-grid[=grid type] Specify the default grid manager:
              ug [default]: use UG
              alu: use ALUGRID_SIMPLEX],
              [ if test "x$enableval" = "xug"
                    then AC_DEFINE([UGGRID], [], [default choice replacing NOGRID])
                elif test "x$enableval" = "xalu"
                    then AC_DEFINE([ALUGRID_SIMPLEX], [], [default choice replacing NOGRID])
                else echo "Error!"; echo "Unknown grid manager"; exit -1
                fi])

#AC_ARG_ENABLE(grid,
#              [--enable-grid[=grid type] Specify the default grid manager:
#              ug [default]: use UGGRID
#              alu: use ALUGRID_SIMPLEX],
#              [ if test "x$enableval" = "xug"
#                    then AC_SUBST(MY_GRIDTYPE, [UGGRID])
#                elif test "x$enableval" = "xalu"
#                    then AC_SUBST(MY_GRIDTYPE, [ALUGRID_SIMPLEX])
#                else echo "Error!"; echo "Unknown grid manager"; exit -1
#                fi])


AC_ARG_ENABLE(dim,
              [--enable-dim[=dimension] Specify the default grid dimension:
              2 [default]: use 2 dimensions
              3: use 3 dimensions],
              [ if test "x$enableval" = "x2"
                    then AC_SUBST(MY_GRIDDIM, 2)
                elif test "x$enableval" = "x3"
                    then AC_SUBST(MY_GRIDDIM, 3)
                else echo "No GRRIDDIM specified. Assuming default value GRIDDIM = 2"
                    echo "Use --enable-dim=value to set configure default."
                    AC_SUBST(MY_GRIDDIM, 2)
                fi])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-localfunctions dune-istl dune-pdelab
# this implies checking for [dune-common], [dune-grid], [dune-localfunctions], [dune-istl], [dune-pdelab]
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR(could not find required version of GSLcblas)])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR(could not find required version of GSL)])
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$LIBS $DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/iPBS/Makefile
  m4/Makefile
  iPBS.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
