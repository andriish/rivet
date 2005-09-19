# Check for CLHEP.
AC_DEFUN([AC_CHECK_CLHEP],
[AC_MSG_CHECKING([if CLHEPPATH is set])
if test -z "$CLHEPPATH"; then
  if test "x$prefix" == "xNONE"; then
    CLHEPPATH=$ac_default_prefix
  else
    CLHEPPATH=$prefix
  fi
  AC_MSG_RESULT([no (using $CLHEPPATH)])
else
  AC_MSG_RESULT([yes ($CLHEPPATH)])
fi
AC_ARG_VAR(CLHEPPATH,[The path to where CLHEP is installed. Default is $prefix.])

AC_MSG_CHECKING([if CLHEPLIB is set])
if test -z "$CLHEPLIB"; then
  CLHEPLIB="-lCLHEP"
  AC_MSG_RESULT([no (using $CLHEPLIB)])
else
  AC_MSG_RESULT([yes ($CLHEPLIB)])
fi
AC_ARG_VAR(CLHEPLIB,[The argument to be used when linking the CLHEP library. Default is "-lCLHEP".])

AC_MSG_CHECKING([if CLHEPINCLUDE is set])
if test -z "$CLHEPINCLUDE"; then
  CLHEPINCLUDE=-I$CLHEPPATH/include
  AC_MSG_RESULT([no (using $CLHEPINCLUDE)])
else
  AC_MSG_RESULT([yes ($CLHEPINCLUDE)])
fi
AC_ARG_VAR(CLHEPINCLUDE,[The argument used when compiling source files which uses CLHEP. Default is "-I$CLHEPPATH/include".])

dnl ###############################
dnl ###############################
dnl ###############################

dnl Now lets see if the libraries work properly
oldLIB="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$LIBS $CLHEPLIB"
LDFLAGS="$LDFLAGS -L$CLHEPPATH/lib"
CPPFLAGS="$CPPFLAGS $CLHEPINCLUDE"

dnl check CLHEP first
AC_MSG_CHECKING([that CLHEP works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <CLHEP/Random/Random.h>
namespace CLHEP {}]], [[using namespace CLHEP; HepRandom r; r.flat();]])],AC_MSG_RESULT(yes),[AC_MSG_RESULT(no) 
AC_MSG_ERROR(CLHEP must be installed to continue.)
AC_SUBST(CLHEPPATH)
AC_SUBST(CLHEPLIB)
AC_SUBST(CLHEPINCLUDE)
])

LIBS="$oldLIB"
LDFLAGS="$oldLDFLAGS"
CPPFLAGS="$oldCPPFLAGS"
])

AC_DEFUN([AC_EMPTY_SUBST],
[EMPTY=""
AC_SUBST(EMPTY)
])

