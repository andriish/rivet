## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[1.0],[],[Rivet])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_LANG(C++)

AC_CONFIG_SRCDIR([src/RivetHandler.cc])
AC_CONFIG_HEADERS([include/Rivet/Config/RivetConfig.h])

## Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

## More comprehensible error messages if possible. See:
## http://www.bdsoft.com/tools/stlfilt.html
AC_PATH_PROG(GFILT, gfilt, $CXX, $PATH:$HOME/bin:$HOME/local/bin)
CXX=$GFILT

## Find CLHEP
AC_SEARCH_CLHEP([AC_TEST_CLHEP], [AC_MSG_ERROR([CLHEP must be installed (for now)])])
AC_CEDAR_LIBRARYANDHEADERS([HepMC], , , [AC_MSG_ERROR([HepMC is required])])

## Logging framework
AC_CEDAR_LIBRARYANDHEADERS([log4cpp], , , [AC_MSG_ERROR([log4cpp is required])])

AM_CPPFLAGS="-I\$(top_builddir)/include -I\$(HEPMCINCPATH) -I\$(LOG4CPPINCPATH) "
AM_CXXFLAGS="-ansi -pedantic -Wall"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_EMPTY_SUBST

AC_PATH_PROG(PERL, perl)

AC_CONFIG_FILES([
include/Makefile
include/Rivet/Makefile
src/Makefile
src/Tools/Makefile
src/Projections/Makefile
src/Analysis/Makefile
src/Test/Makefile
lib/Makefile
Makefile],)

AC_OUTPUT
