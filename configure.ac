## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[1.0a1],[],[Rivet])

AM_INIT_AUTOMAKE
AC_LANG(C++)

AC_CONFIG_SRCDIR([src/RivetInterface.cc])
AC_CONFIG_HEADERS([include/Rivet/Config/RivetConfig.hh])

## OS X
AC_CHECK_TOOL(SWVERS, sw_vers)
if test x$SWVERS != x; then
  MACOSX_DEPLOYMENT_TARGET=`$SWVERS -productVersion | cut -f 1,2 -d.`
  AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
fi

## Set default compiler flags
if test "x$CXXFLAGS" == "x"; then CXXFLAGS="-O3"; fi

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
AC_CEDAR_LIBRARYANDHEADERS([KtJet], , , [AC_MSG_ERROR([KtJet library is required])])
AC_CEDAR_LIBRARYANDHEADERS([HepMC], , , [AC_MSG_ERROR([HepMC is required])])
AC_CEDAR_LIBRARYANDHEADERS([HepPID], , , [AC_MSG_ERROR([HepPID is required])])
AC_CEDAR_LIBRARYANDHEADERS([HepPDT], , , [AC_MSG_ERROR([HepPDT is required])])
AC_CEDAR_LIBRARYANDHEADERS([fastjet], , , [AC_MSG_WARN([fastjet not found])])

## Conditionally build the test programs
AC_ARG_WITH([testprogs], [AC_HELP_STRING(--with-testprogs, [build internal test programs])], [], [])
AM_CONDITIONAL(WITH_TESTPROGS, [test x$with_testprogs == xyes])

## Find ROOT
AC_ARG_WITH([ROOT], [AC_HELP_STRING(--with-root, [build ROOT persistency])], [], [])
AC_SEARCH_ROOT([with_root=yes
  AC_TEST_ROOT], [AC_MSG_WARN([ROOT functionalities not available])])
## Enable ROOT analysis framework if requested
AM_CONDITIONAL(WITH_ROOT, [test x$with_root == xyes])

## Set default build flags
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include -I\$(CLHEPINCPATH) -I\$(HEPMCINCPATH) -I\$(HEPPIDINCPATH) -I\$(HEPPDTINCPATH) -I\$(KTJETINCPATH)"
AM_CXXFLAGS="$AM_CXXFLAGS -ansi -pedantic -Wall"

## Set fastjet/fastjet include subdirectory in include search path
#if test -d ${FASTJETINCPATH}  
if test "${FASTJETINCPATH}" != ""; then
  AC_MSG_NOTICE([adding subdirectory: ${FASTJETINCPATH}/fastjet into inlcude path])
  AM_CPPFLAGS="$AM_CPPFLAGS -I\$(FASTJETINCPATH)/fastjet"
  AM_CXXFLAGS="$AM_CXXFLAGS -I\$(FASTJETINCPATH)/fastjet"
fi


## Debug flag (default=none)
AC_ARG_ENABLE([debug], [AC_HELP_STRING(--enable-debug, [build with debugging symbols (default=no)])], [], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_CXXFLAGS="$AM_CXXFLAGS -g"]
fi

## KtJet double precision build
## Needed? Does Rivet code actually use this preprocessor directive? Anyway, 
## *if* needed, use AC_ARG_ENABLE instead of asking users to hack this file.
#AM_CPPFLAGS="$AM_CPPFLAGS -DKTDOUBLEPRECISION"

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_EMPTY_SUBST
AC_CONFIG_FILES([
include/Makefile
include/Rivet/Makefile
src/Makefile
src/Tools/Makefile
src/Projections/Makefile
src/Analyses/Makefile
src/Test/Makefile
data/Makefile
bin/Makefile
bin/rivet-config
src/Tools/RivetPaths.cc
Makefile],)

AC_OUTPUT
