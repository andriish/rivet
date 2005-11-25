# Search for CLHEP in standard directories using standard CLHEP names
#AC_SEARCH_CLHEP
#----------------------------------------
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


#AC_CEDAR_CHECK_GENERATOR([Name], [Version])
#----------------------------------------
AC_DEFUN([AC_CEDAR_CHECK_GENERATOR], [
  m4_define([cedar_prettyname], [$1 $2])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1$2], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1$2], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_libname], [lib@&t@cedar_pkgname@&t@.a])dnl
  m4_define([cedar_incname], [cedar_pkgname])dnl

  AC_MSG_RESULT([Checking for cedar_prettyname:]) 
  # Don't know why this isn't working by default:
  test x${prefix} == xNONE && prefix=${ac_default_prefix}
  AC_ARG_WITH([cedar_pkgname],
              AC_HELP_STRING(--with-@&t@cedar_pkgname, path to cedar_prettyname generator),
              [pkgpath=$with_@&t@cedar_pkgname],
              [pkgpath=${prefix}])
  pkglib=yes; pkginc=yes; pkggood=yes

  if test x$pkgpath != xno; then
    AC_CHECK_FILE([$pkgpath/lib/cedar_libname],
                  [AC_MSG_NOTICE([Found cedar_libname])], 
                  [pkglib=no; AC_MSG_NOTICE([trying to find cedar_prettyname library in build directory structure...])])
    if test x$pkglib == xno; then
      AC_CHECK_FILE([$pkgpath/src/cedar_libname],
                    [pkglib=yes; AC_MSG_NOTICE([Found cedar_libname])], 
                    [pkggood=no; AC_MSG_RESULT([cedar_prettyname library "cedar_libname" is not in a standard location])])
    fi
    if test x$pkggood != xno; then
      AC_CHECK_FILE([$pkgpath/include/cedar_incname],
                    [AC_MSG_NOTICE([Found cedar_prettyname header directory])],
                    [pkggood=no; AC_MSG_RESULT([cedar_prettyname header directory is not in a standard location])])
      #              [pkginc=no; AC_MSG_NOTICE([trying to find cedar_prettyname header directory in build directory structure...])])
      #if test x$pkginc == xno; then
      #  AC_CHECK_FILE([$pkgpath/include],
      #                [pkginc=yes; AC_MSG_NOTICE([Found cedar_prettyname header directory])], 
      #                [pkggood=no; AC_MSG_RESULT([cedar_prettyname header directory is not in a standard location])])
      #fi
    fi
  else
    pkggood=no
  fi
  if test x$pkggood != xno; then
    AC_MSG_RESULT([cedar_prettyname paths verified])
  else
    AC_MSG_RESULT([Not building against cedar_prettyname])
  fi
  # Note quoting subtlty on first arg of AM_CONDITIONAL:
  AM_CONDITIONAL(WITH_@&t@cedar_PKGNAME, [test x$pkggood != xno])
  cedar_PKGNAME@&t@PATH=${pkgpath}
  AC_SUBST(cedar_PKGNAME@&t@PATH)
])


#AC_CEDAR_CHECK_GENHEAD([Name], [Version])
#----------------------------------------
AC_DEFUN([AC_CEDAR_CHECK_GENHEAD], [
  m4_define([cedar_prettyname], [$1 $2])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1$2], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1$2], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_libname], [lib@&t@cedar_pkgname@&t@.a])dnl
  m4_define([cedar_incname], [cedar_pkgname])dnl

  AC_MSG_RESULT([Checking for cedar_prettyname:]) 
  # Don't know why this isn't working by default:
  test x${prefix} == xNONE && prefix=${ac_default_prefix}
  AC_ARG_WITH([cedar_pkgname],
              AC_HELP_STRING(--with-@&t@cedar_pkgname, path to cedar_prettyname generator),
              [pkgpath=$with_@&t@cedar_pkgname],
              [pkgpath=${prefix}])
  pkglib=yes; pkginc=yes; pkggood=yes

  if test x$pkgpath != xno; then
    if test x$pkggood != xno; then
      AC_CHECK_FILE([$pkgpath/include/cedar_incname],
                    [AC_MSG_NOTICE([Found cedar_prettyname header directory])],
                    [pkggood=no; AC_MSG_RESULT([cedar_prettyname header directory is not in a standard location])])
      #              [pkginc=no; AC_MSG_NOTICE([trying to find cedar_prettyname header directory in build directory structure...])])
      #if test x$pkginc == xno; then
      #  AC_CHECK_FILE([$pkgpath/include],
      #                [pkginc=yes; AC_MSG_NOTICE([Found cedar_prettyname header directory])], 
      #                [pkggood=no; AC_MSG_RESULT([cedar_prettyname header directory is not in a standard location])])
      #fi
    fi
  else
    pkggood=no
  fi
  if test x$pkggood != xno; then
    AC_MSG_RESULT([cedar_prettyname paths verified])
  else
    AC_MSG_RESULT([Not building against cedar_prettyname])
  fi
  # Note quoting subtlty on first arg of AM_CONDITIONAL:
  AM_CONDITIONAL(WITH_@&t@cedar_PKGNAME, [test x$pkggood != xno])
  cedar_PKGNAME@&t@PATH=${pkgpath}
  AC_SUBST(cedar_PKGNAME@&t@PATH)
])


#AC_CEDAR_CHECK_REQDPKG([Name], [Version])
#--------------------------------------
AC_DEFUN([AC_CEDAR_CHECK_REQDPKG], [
  m4_define([cedar_prettyname], [$1])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_libname], [lib@&t@cedar_pkgname@&t@.a])dnl
  m4_define([cedar_incname], [cedar_pkgname])dnl

  AC_MSG_RESULT([Checking for cedar_prettyname:])
  # Don't know why this isn't working by default:
  test x${prefix} == xNONE && prefix=${ac_default_prefix}
  AC_ARG_WITH([cedar_pkgname],
              AC_HELP_STRING(--with-@&t@cedar_pkgname, path to cedar_prettyname),
              [pkgpath=$with_@&t@cedar_pkgname], 
              [pkgpath=${prefix}])
  pkglib=yes; pkginc=yes; pkggood=yes

  if test x$pkgpath != xno; then
    AC_CHECK_FILE([${pkgpath}/lib/cedar_libname], [],
                  [pkglib=no; AC_MSG_NOTICE([trying to find cedar_prettyname library in build directory structure...])])
    if test x$pkglib == xno; then
      AC_CHECK_FILE([${pkgpath}/src/cedar_libname], [pkglib=yes], 
                    [AC_MSG_ERROR([cedar_prettyname library "cedar_libname" is not in a standard location])])
    fi
    AC_MSG_NOTICE([Found cedar_libname])
    AC_CHECK_FILE([${pkgpath}/include/cedar_incname], [], 
    #              [pkginc=no; AC_MSG_NOTICE([trying to find cedar_prettyname header directory in build directory structure...])])
    #if test x$pkginc == xno; then
    #  AC_CHECK_FILE([${pkgpath}/inc], [pkginc=yes],
                    [AC_MSG_ERROR([cedar_prettyname header directory is not in a standard location])])
    #fi
    AC_MSG_NOTICE([Found cedar_prettyname header directory])
  else
    AC_MSG_ERROR([You've specified "--without-@&t@cedar_pkgname", but cedar_prettyname is required to build ${PACKAGE_NAME}])
  fi
  AC_MSG_RESULT([cedar_prettyname paths verified])
  cedar_PKGNAME@&t@PATH=${pkgpath}
  AC_SUBST(cedar_PKGNAME@&t@PATH)
])


#AC_CEDAR_CHECK_THEPEG()
#--------------------------------------
AC_DEFUN([AC_CEDAR_CHECK_THEPEG], [
  #m4_define([cedar_prettyname], [$1])dnl
  #m4_define([cedar_PKGNAME], [translit([translit([$1], [a-z], [A-Z])], [.])])dnl
  #m4_define([cedar_pkgname], [translit([translit([$1], [A-Z], [a-z])], [.])])dnl
  #m4_define([cedar_libname], [lib@&t@cedar_pkgname@&t@.a])dnl
  #m4_define([cedar_incname], [cedar_pkgname])dnl

  AC_MSG_RESULT([Checking for ThePEG:])
  # Don't know why this isn't working by default:
  test x${prefix} == xNONE && prefix=${ac_default_prefix}
  AC_ARG_WITH(thepeg,
              AC_HELP_STRING(--with-thepeg, path to ThePEG),
              [pkgpath=$with_thepeg], 
              [pkgpath=${prefix}])
  pkglib=yes; pkginc=yes; pkggood=yes

  if test x$pkgpath != xno; then
    AC_CHECK_FILE([${pkgpath}/lib/ThePEG/libThePEG.so], [],
                  [pkglib=no; AC_MSG_ERROR([ThePEG library "libThePEG.so" is not in a standard location])])
    AC_MSG_NOTICE([Found ThePEG libraries])
    AC_CHECK_FILE([${pkgpath}/include/ThePEG], [],
                  [pkginc=no; AC_MSG_ERROR([cedar_prettyname header directory is not in a standard location])])
    AC_MSG_NOTICE([Found ThePEG header directory])
  else
    AC_MSG_ERROR([You've specified "--without-thepeg", but ThePEG is required to build ${PACKAGE_NAME}])
  fi
  AC_MSG_RESULT([ThePEG paths verified])
  THEPEGPATH=${pkgpath}
  AC_SUBST(THEPEGPATH)
])
