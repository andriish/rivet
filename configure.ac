## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Rivet],[1.1.3a0],[rivet@projects.hepforge.org],[Rivet])
AC_CONFIG_SRCDIR([src/Analysis.cc])
AC_CONFIG_HEADERS([include/Rivet/Config/DummyConfig.hh include/Rivet/Config/RivetConfig.hh include/Rivet/Config/BuildOptions.hh])
AM_INIT_AUTOMAKE(dist-bzip2)
AC_CONFIG_MACRO_DIR([m4])
AC_SUBST(LT_OBJDIR)

## Package-specific #defines
AC_DEFINE_UNQUOTED(RIVET_VERSION, "$PACKAGE_VERSION", "Rivet version string")
AC_DEFINE_UNQUOTED(RIVET_NAME, "$PACKAGE_NAME", "Rivet name string")
AC_DEFINE_UNQUOTED(RIVET_STRING, "$PACKAGE_STRING", "Rivet name and version string")
AC_DEFINE_UNQUOTED(RIVET_TARNAME, "$PACKAGE_TARNAME", "Rivet short name string")
AC_DEFINE_UNQUOTED(RIVET_BUGREPORT, "$PACKAGE_BUGREPORT", "Rivet contact email address")

## Set library version info code
## Set for version 1.1.3
VERSIONINFOFLAGS="-version-info 5:0:0"
AC_SUBST(VERSIONINFOFLAGS)

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
RAW_LIB_EXTN=$(echo "$shrext_cmds" | sed -e 's/^\.//')
LIB_SUFFIX=\\\"$shrext_cmds\\\"
AC_SUBST(RAW_LIB_EXTN)
AC_SUBST([LIB_SUFFIX])

## GNU Scientific Library
AC_CEDAR_LIBRARY([gsl], , , [AC_MSG_ERROR([GNU Scientific Library is required])])
AC_CEDAR_LIBRARY([gslcblas], , , [AC_MSG_ERROR([GSL CBLAS is required])])

## Boost C++ header library
AC_CEDAR_HEADERS([Boost], , , [AC_MSG_ERROR([Boost is required])])
oldCPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$BOOSTINCPATH"
BOOSTERRMSG="You need at least version 1.34.0 of Boost for this installation."
BOOSTERRMSG="$BOOSTERRMSG If you want to use Boost from a non-standard location,"
BOOSTERRMSG="$BOOSTERRMSG please specify '--with-boost-incpath=/path/to/boost_1_35/include'"
BOOSTERRMSG="$BOOSTERRMSG as option to the configure script or run 'export"
BOOSTERRMSG="$BOOSTERRMSG BOOSTFLAGS=--with-boost-incpath=/path/to/boost_1_35/include'"
BOOSTERRMSG="$BOOSTERRMSG before you run the rivet-bootstrap script."
AC_CHECK_HEADER([boost/foreach.hpp], [], [AC_MSG_ERROR([Boost foreach not found. $BOOSTERRMSG])])
AC_CHECK_HEADER([boost/smart_ptr.hpp], [], [AC_MSG_ERROR([Boost smart_ptr not found. $BOOSTERRMSG])])
AC_CHECK_HEADER([boost/lexical_cast.hpp], [], [AC_MSG_ERROR([Boost lexical_cast not found. $BOOSTERRMSG])])
AC_CHECK_HEADER([boost/assign.hpp], [], [AC_MSG_ERROR([Boost assign not found. $BOOSTERRMSG])])
AC_CHECK_HEADER([boost/random.hpp], [], [AC_MSG_ERROR([Boost random numbers not found. $BOOSTERRMSG])])
CPPFLAGS=$oldCPPFLAGS


## HepMC library
AC_CEDAR_LIBRARYANDHEADERS([HepMC], , , [AC_MSG_ERROR([HepMC is required])])


## FastJet library
AC_CEDAR_LIBRARYANDHEADERS([fastjet], , , [AC_MSG_ERROR([FastJet is required])])
## Decide/determine whether to build with FastJet JADE + Durham support
AC_ARG_ENABLE(jade, [AC_HELP_STRING(--enable-jade, 
  [use JADE support from patched FastJet @<:@default=no@:>@])], 
  [], [enable_jade=no])
AC_MSG_CHECKING([whether to try and use the FastJet JADE + Durham plugin patch])
AC_MSG_RESULT($enable_jade)
AM_CONDITIONAL(ENABLE_JADE, [test x$enable_jade != xno])
#
AC_PATH_PROG(FJCONFIG, fastjet-config, [], $FASTJETPATH/bin)
if test -f "$FJCONFIG"; then
  FASTJETCONFIGLIBADD="$($FJCONFIG --plugins --shared --libs)";
else
  FASTJETCONFIGLIBADD="-L\$(FASTJETLIBPATH) -l\$(FASTJETLIBNAME)";
fi;
AC_SUBST(FASTJETCONFIGLIBADD)


## Build Doxygen if possible
AC_ARG_ENABLE([doxygen],
  [AC_HELP_STRING(--disable-doxygen, [don't try to make Doxygen documentation])],
  [], [enable_doxygen=yes])
if test x$enable_doxygen == xyes; then
   AC_PATH_PROG(DOXYGEN, doxygen)
fi
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])


## Decide whether to build the Python extension
AZ_PYTHON_DEFAULT
AC_ARG_ENABLE(pyext, [AC_HELP_STRING(--enable-pyext, 
  [build Rivet Python module @<:@default=yes@:>@])], [], [enable_pyext=yes])
AM_CONDITIONAL(ENABLE_PYRIVET, [test x$enable_pyext == xyes])
if test x$enable_pyext == xyes; then
  AZ_PYTHON_PATH
  AZ_PYTHON_VERSION_ENSURE([2.3])
fi


dnl ## Optional ROOT analysis framework
dnl AC_ARG_ENABLE([ROOT], [AC_HELP_STRING(--enable-root, 
dnl   [build with ROOT output @<:@default=yes@:>@])], [], [])
dnl if test x$enable_root != xno; then
dnl   AC_SEARCH_ROOT([enable_root=yes; AC_TEST_ROOT], [AC_MSG_WARN([ROOT is not available])])
dnl else
dnl   AC_MSG_NOTICE([Not using ROOT])
dnl fi
dnl AM_CONDITIONAL([WITH_ROOT], [test x$enable_root == xyes])


## Set default build flags
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(BOOSTINCPATH) \$(BOOST_CPPFLAGS)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(HEPMCINCPATH)"
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(FASTJETINCPATH)"
AM_CXXFLAGS="$AM_CXXFLAGS -DLIB_SUFFIX=$LIB_SUFFIX "
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic "])
AC_CEDAR_CHECKCXXFLAG([-ansi], [AM_CXXFLAGS="$AM_CXXFLAGS -ansi "])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
dnl if test x$with_root == xyes; then
dnl   AC_CEDAR_CHECKCXXFLAG([-Wno-long-long], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long "])
dnl fi


## Debug flag (default=none)
AC_ARG_ENABLE([debug], [AC_HELP_STRING(--enable-debug, 
  [build with debugging symbols  @<:@default=no@:>@])], [], [enable_debug=no])
if test x$enable_debug == xyes; then
  AC_CEDAR_CHECKCXXFLAG([-g], [AM_CXXFLAGS="$AM_CXXFLAGS -g "])
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
AC_CONFIG_FILES(pyext/setup.py)
AC_CONFIG_FILES(data/Makefile)
AC_CONFIG_FILES(data/refdata/Makefile)
AC_CONFIG_FILES(doc/Makefile)
AC_CONFIG_FILES(doc/rivetversion.sty)
AC_CONFIG_FILES(autopackage/default.apspec)
AC_CONFIG_FILES(autopackage/data.apspec)
AC_CONFIG_FILES(bin/Makefile bin/rivet-config)
AC_CONFIG_FILES(plugindemo/Makefile)
AC_CONFIG_FILES(Makefile Doxyfile)

AC_OUTPUT

if test x$enable_pyrivet == xyes; then
   cat <<EOF

************************************************************
RIVET CONFIGURED!

Now build and install (to the $prefix tree) with e.g.
make -j2 && make -j2 install

To use Rivet, we recommend reading HepMC files from a file
or pipe (the latter may be made with mkfifo) using the 
'rivet' executable.

For a more pleasant command line experience, you can include 
the data/rivet-completion file into your .bashrc file,
or your bash_completion.d directory if you have one.
************************************************************
EOF
fi
