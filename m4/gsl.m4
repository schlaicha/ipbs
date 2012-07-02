# checks for the availability of GSL for computing the elliptic integrals

AC_DEFUN([_gsl_lib_path], 
    [
	my_include_path=include
	my_lib_path=lib
	my_gsl_found=yes
	if test ! -f "$1/$my_include_path/gsl/$2" ; then
	   #Try to find headers under gsl
	   my_include_path=include
	   if test ! -f "$with_gsl/$my_include_path/$2" ; then
		my_include_path=src
                my_lib_path=""
                if test ! -f "$with_gsl/$my_include_path/$2"; then
                    my_gsl_found=no
                fi
            fi
        fi
    ]
)

AC_DEFUN([_gsl_search_default],
    [
        with_gsl=/usr
        _gsl_lib_path($with_gsl, "gsl_sf_ellint.h")
        
        if test "$my_gsl_found" = "no"; then
            with_gsl=/usr/local
            _gsl_lib_path($with_gsl, "gsl_sf_ellint.h")
        fi
    ]
)

AC_DEFUN([DUNE_PATH_GSL],[
        AC_MSG_CHECKING(for GNU Scientific Library) 
        AC_REQUIRE([AC_PROG_CC])
        AC_REQUIRE([ACX_BLAS])

  #
  # USer hints ...
  #
        my_lib_path=""
        my_include_path=""
        AC_ARG_WITH([gsl],
            [AC_HELP_STRING([--with-gsl],[user defined path to GSL])],
            [dnl
                if test "$withval" != no ; then
                    # get absolute path
                    with_gsl=`eval cd $withval 2>&1 && pwd`
                    if test "$withval" = yes; then
                        # Search in default locations
                        _gsl_search_default
                    else
                        # Search for the headers in the specified location
                        _gsl_lib_path("$with_gsl", "gsl_sf_ellint.h")
                    fi
                fi
                ],
            [dnl
                # Search in default locations
                    _gsl_search_default
                ])
        AC_MSG_RESULT([$with_gsl/$my_include_path])
        AC_ARG_WITH([gsl_lib],
            [AC_HELP_STRING([--with-gsl-lib],[The name of the static GSL to link to. By default the shared library with the name gslcblas.a is tried])],
            [
                if test "$withval" != no ; then
                    with_gsl_lib=$withval
                fi
            ]
        )

  # store old values
        ac_save_LDFLAGS="$LDFLAGS"
        ac_save_CPPFLAGS="$CPPFLAGS"
        ac_save_LIBS="$LIBS"

  ## do nothing if --without-gsl is used
        if test x"$with_gsl" != x"no" ; then

      # defaultpath
            GSL_LIB_PATH="$with_gsl/$my_lib_path"
            GSL_INCLUDE_PATH="$with_gsl/$my_include_path"
	echo $GSL_INCLUDE_PATH

            GSL_LDFLAGS="-L$GSL_LIB_PATH $DUNEMPILDFLAGS"

      # set variables so that tests can use them
            CPPFLAGS="$CPPFLAGS -I$GSL_INCLUDE_PATH $DUNEMPICPPFLAGS"

      # check for central header
            AC_CHECK_HEADER([gsl/gsl_sf_ellint.h],[
                    GSL_CPPFLAGS="$CPPFLAGS"
                    HAVE_GSL="1"],[
                    HAVE_GSL="0"
                    AC_MSG_ERROR([gsl_sf_ellint.h not found in $GSL_INCLUDE_PATH with $CPPFLAGS])]
            )

            GSL_CPPFLAGS="-I$GSL_INCLUDE_PATH $DUNEMPICPPFLAGS"

      # if header is found check for the libs

            if test x$HAVE_GSL = x1 ; then

                # set variables so that tests can use them
                OLDFLAGS="$LDFLAGS"
                LDFLAGS="$LDFLAGS -L$GSL_LIB_PATH $DUNEMPILDFLAGS"
                LIBS="$BLAS_LIBS $LIBS $FLIBS $DUNEMPILIBS"

		AC_CHECK_LIB([gslcblas],[cblas_dgemm],
                    [dnl
                        GSL_LIBS="-lgslcblas"
                        GSL_LDFLAGS="$LDFLAGS"
                        HAVE_GSL="1"
                        AC_MSG_RESULT(yes)
                        ],[dnl
                        HAVE_GSL="0"
                        AC_MSG_ERROR(libgslcblas not found)])

		AC_CHECK_LIB([gsl -lgslcblas],[gsl_sf_ellint_Kcomp],
                    [dnl
                        GSL_LIBS="$GSL_LIBS -lgsl $LIBS"
                        GSL_LDFLAGS="$LDFLAGS"
                        HAVE_GSL="1"
                        AC_MSG_RESULT(yes)
                        ],[dnl
                        HAVE_GSL="0"
                        AC_MSG_ERROR(libgsl not found)])


                if test "$HAVE_GSL" = 0; then
                    #check for the static library
                    if test x$with_gsl_lib = x ; then
                        with_gsl_lib=libsglcblas.a
                    fi
                    AC_MSG_CHECKING([static GSL $with_gsl_lib in "$GSL_LIB_PATH"])

                    if test -f "$GSL_LIB_PATH/$with_gsl_lib"; then
                        LIBS="$GSL_LIB_PATH/$with_gsl_lib $LIBS"
                        LDFLAGS="$OLDFLAGS"
                        AC_CHECK_FUNC(gsl_sf_ellint_Kcomp_e,
                            [
                                GSL_LIBS="$LIBS"
                                GSL_LDFLAGS="$LDFLAGS"
                                HAVE_GSL="1"
                                AC_MSG_RESULT(yes)
                                ],
                            [ 
                                HAVE_GSL="0"
                                AC_MSG_RESULT(failed)
                            ]
                        )
                    else
                        AC_MSG_RESULT(failed)
                        HAVE_GSL="0"
                    fi
                fi
            fi
      # pre-set variable for summary
      #with_gsl="no"

      # did it work?
            AC_MSG_CHECKING([GSL in $with_gsl])
            if test x$HAVE_GSL = x1 ; then
                AC_SUBST(GSL_LDFLAGS, $GSL_LDFLAGS)
                AC_SUBST(GSL_LIBS, $GSL_LIBS)
                AC_SUBST(GSL_CPPFLAGS, $GSL_CPPFLAGS)
                AC_DEFINE(HAVE_GSL, 1, [Define to 1 if GSL is found])
                AC_MSG_RESULT(ok)

    # add to global list
                DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $GSL_LDFLAGS"
                DUNE_PKG_LIBS="$DUNE_PKG_LIBS $GSL_LIBS"
                DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $GSL_CPPFLAGS"

    # re-set variable correctly
                with_gsl="yes"
            else
                with_gsl="no"
                AC_MSG_RESULT(failed)
            fi

  # end of "no --without-gsl"
        else
            with_gsl="no"
        fi

  # tell automake
        AM_CONDITIONAL(GSL, test x$HAVE_GSL = x1)

  # restore variables
        LDFLAGS="$ac_save_LDFLAGS"
        CPPFLAGS="$ac_save_CPPFLAGS"
        LIBS="$ac_save_LIBS"

    DUNE_ADD_SUMMARY_ENTRY([GSL],[$with_gsl])

    ]
)

dnl Local Variables:
dnl mode: shell-script
dnl End:
