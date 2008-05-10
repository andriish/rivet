## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[1.1.0],[rivet@projects.hepforge.org],[Rivet])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/Analysis.cc])
AC_CONFIG_HEADERS([include/Rivet/Config/RivetConfig.hh])

## OS X
AC_CEDAR_OSX

## LCG platform tag
AC_LCG_TAG

## Set default compiler flags
if test "x$CXXFLAGS" == "x"; then CXXFLAGS="-O3"; fi

## Make compiler error messages more readable if gfilt is installed.
dnl AC_CEDAR_CXXFILTER

## Checks for programs.
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

## Work out library suffix for the build
LIB_SUFFIX=\\\"$shrext_cmds\\\"
AC_SUBST([LIB_SUFFIX])

## Find required libraries
AC_CEDAR_LIBRARY([gsl], , , [AC_MSG_ERROR([GNU Scientific Library is required])])
AC_CEDAR_LIBRARY([gslcblas], , , [AC_MSG_ERROR([GSL CBLAS is required])])
AC_CEDAR_LIBRARYANDHEADERS([HepMC], , , [AC_MSG_ERROR([HepMC is required])])
AC_CEDAR_LIBRARYANDHEADERS([fastjet], , , [AC_MSG_ERROR([FastJet is required])])

## Decide/determine whether to build the rivetgun executable
AC_ARG_ENABLE(rivetgun, [AC_HELP_STRING(--enable-rivetgun, [build RivetGun executable @<:@default=yes@:>@])], [], [enable_rivetgun=yes])
AC_CEDAR_LIBRARYANDHEADERS([AGILe], , , [AC_MSG_WARN([AGILe not found])])
AM_CONDITIONAL(ENABLE_RIVETGUN, [test x$enable_rivetgun == xyes && test x$WITH_AGILE != xno])
if test x$enable_rivetgun == xyes && test x$WITH_AGILE == xno; then
  AC_MSG_NOTICE([AGILe not found: not building rivetgun])
fi

## Decide whether to build the Python extension
AC_ARG_ENABLE(pyrivet, [AC_HELP_STRING(--enable-pyrivet, 
  [build pyrivet Python module @<:@default=no@:>@])], [], [enable_pyrivet=no])
AM_CONDITIONAL(ENABLE_PYRIVET, [test x$enable_pyrivet == xyes])
if test x$enable_pyrivet == xyes; then
  AM_PATH_PYTHON(2.3)
  AC_PROG_SWIG(1.3.21)
  SWIG_ENABLE_CXX
  SWIG_PYTHON
fi

## Decide/determine whether to build with FastJet JADE + Durham support
AC_ARG_ENABLE(jade, [AC_HELP_STRING(--enable-jade, 
  [use JADE support from patched FastJet @<:@default=no@:>@])])
AC_MSG_CHECKING([whether to try and use the FastJet JADE + Durham plugin patch])
AC_MSG_RESULT([test x$enable_jade != x])
AM_CONDITIONAL(ENABLE_JADE, [test x$enable_jade != x])

## Optional ROOT analysis framework
AC_ARG_ENABLE([ROOT], [AC_HELP_STRING(--enable-root, 
  [build with ROOT output @<:@default=yes@:>@])], [], [])
if test x$enable_root != xno; then
  AC_SEARCH_ROOT([enable_root=yes; AC_TEST_ROOT], [AC_MSG_WARN([ROOT is not available])])
else
  AC_MSG_NOTICE([Not using ROOT])
fi
AM_CONDITIONAL([WITH_ROOT], [test x$enable_root == xyes])

## Set default build flags
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include"
AM_CPPFLAGS="$AM_CPPFLAGS \$(BOOST_CPPFLAGS)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(HEPMCINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(HEPPIDINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(HEPPDTINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(FASTJETINCPATH)"
AM_CXXFLAGS="$AM_CXXFLAGS -ansi -pedantic -Wall -DLIB_SUFFIX=$LIB_SUFFIX "
if test x$with_root == xyes; then
  AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long "
fi

## Debug flag (default=none)
AC_ARG_ENABLE([debug], [AC_HELP_STRING(--enable-debug, 
  [build with debugging symbols  @<:@default=no@:>@])], [], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_CXXFLAGS="$AM_CXXFLAGS -g"]
fi

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_EMPTY_SUBST
AC_CONFIG_FILES(include/Makefile include/Rivet/Makefile)
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(src/Tools/Makefile)
AC_CONFIG_FILES(src/Projections/Makefile)
AC_CONFIG_FILES(src/Analyses/Makefile)
AC_CONFIG_FILES(src/Test/Makefile)
AC_CONFIG_FILES(pyext/Makefile)
AC_CONFIG_FILES(data/Makefile)
AC_CONFIG_FILES(bin/Makefile bin/rivet-config)
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
