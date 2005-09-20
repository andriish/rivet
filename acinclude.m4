# Search for CLHEP in standard directories using standard CLHEP names
AC_DEFUN([AC_SEARCH_CLHEP],
[AC_MSG_CHECKING([if CLHEPPATH, CLHEPLIB and CLHEPINCLUDE are set])
notset=""
if test -z "$CLHEPPATH"; then
  notset="true"
  for dirbase in / /usr $ac_default_prefix $prefix; do
    if test -z "$CLHEPLIB"; then
      for filename in $dirbase/lib/libCLHEP-?.?.?.?.so $dirbase/lib/libCLHEP.so; do
        if test -f $filename; then
          CLHEPPATH=$dirbase
          filename=`basename $filename`
          CLHEPLIB=`echo $filename | sed -e 's/^lib/-l/' -e 's/\.so$//'`
        fi
      done
    else
      filename=`echo $CLHEPLIB | sed -e 's/^-l/lib/'`.so
      if test -f $dirbase/lib/$filename; then
	CLHEPPATH=$dirbase
      fi
    fi
  done
else
  if test -z "$CLHEPLIB"; then
    notset="true"
    for filename in $CLHEPPATH/lib/libCLHEP-?.?.?.?.so CLHEPPATH/lib/libCLHEP.so; do
      if test -f $filename; then
        filename=`basename $filename`
        CLHEPLIB=`echo $filename | sed -e 's/^lib/-l/' -e 's/\.so$//'`
      fi
    done
  fi
fi

if test -z"$CLHEPINCLUDE"; then
  notset="true"
  CLHEPINCLUDE=-I$CLHEPPATH/include
fi

if test -z "$notset"; then
  AC_MSG_RESULT([yes ($CLHEPPATH, $CLHEPLIB and $CLHEPINCLUDE)])
else
  AC_MSG_RESULT([no (found $CLHEPPATH, $CLHEPLIB and $CLHEPINCLUDE)])
fi

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

AC_ARG_VAR(CLHEPPATH,[The path to where CLHEP is installed. Default is $prefix.])
AC_ARG_VAR(CLHEPLIB,[The argument to be used when linking the CLHEP library. Default is "-lCLHEP".])

AC_ARG_VAR(CLHEPINCLUDE,[The argument used when compiling source files which uses CLHEP. Default is "-I$CLHEPPATH/include".])

])

AC_DEFUN([AC_EMPTY_SUBST],
[EMPTY=""
AC_SUBST(EMPTY)
])

