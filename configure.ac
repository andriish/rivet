## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[1.0],[],[Rivet])
AC_CONFIG_AUX_DIR([Config])
AM_INIT_AUTOMAKE
AC_LANG(C++)

AC_CONFIG_SRCDIR([Analysis/RivetHandler.h])
AC_CONFIG_HEADERS([Config/config.h])

## Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

## Find CLHEP
AC_SEARCH_CLHEP([AC_TEST_CLHEP], [AC_MSG_ERROR([CLHEP must be installed])])

AM_CPPFLAGS="-I\$(top_builddir)/include"
AM_CXXFLAGS="-ansi -pedantic -Wall"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_EMPTY_SUBST

AC_PATH_PROG(PERL, perl)

AC_CONFIG_FILES([include/Makefile
CLHEPWrap/Makefile
Config/Makefile
Projections/Makefile
Analysis/Makefile
Analysis/Examples/Makefile
lib/Makefile
Makefile],)

AC_OUTPUT
