# checks for the availability of GSL for computing the elliptic integrals

AC_DEFUN([IPBS_CHECK_GSL], [


AC_ARG_WITH(gsl, AC_HELP_STRING([--with-gsl=PATH}],
                                [specify the path to GSL library]),
                                , [with_gsl=/usr/lib])

#DUNE_CHECK_LIB_EXT([$with_gsl], 
#         [gslcblas],
#         [cblas_dgemm],
#         [GSL_CBLAS_LIBS="-L$with_gsl -lgslcblas"],
#         [AC_MSG_ERROR(could not find required version of GSLcblas)])

DUNE_CHECK_LIB_EXT([$with_gsl], 
        [gsl -lgslcblas],
        [gsl_sf_ellint_Kcomp],
        [GSL_LIBS="-L$with_gsl -lgsl -lgslcblas"
         LIBS="$ac_save_LIBS $GSL_LIBS "
         LDFLAGS="$ac_save_LDFLAGS $GSL_LIBS "
        ],[AC_MSG_ERROR(could not find required version of GSL)])

AC_SUBST(GSL_LIBS, $GSL_LIBS)
AC_SUBST(GSL_LDFLAGS, $GSL_LDFLAGS)
AC_SUBST(GSL_CPPFLAGS, $GSL_CPPFLAGS)

AC_MSG_RESULT(ok)

# add to global list
DUNE_PKG_LIBS="$GSL_LIBS $DUNE_PKG_LIBS"
DUNE_PKG_LDFLAGS="$IPBS_LDFLAGS $DUNE_PKG_LDFLAGS"
DUNE_PKG_CPPFLAGS="$IPBS_CPPFLAGS $DUNE_PKG_CPPFLAGS"

DUNE_ADD_SUMMARY_ENTRY([GSL],[$with_gsl])

])
