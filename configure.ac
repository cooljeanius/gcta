dnl# configure.ac                                  -*- Autoconf -*-
dnl# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.70])
AC_INIT([gcta],[0.93.9],[jian.yang@uq.edu.au])
AC_CONFIG_SRCDIR([cdflib.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_MSG_NOTICE([this configure script is generated from configure.ac by autoconf])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.11.6])
AM_DEP_TRACK
AM_MAINTAINER_MODE
AM_MAKE_INCLUDE
m4_ifdef([AM_WITH_REGEX],[
  AM_WITH_REGEX
])dnl

# Checks for programs.
m4_ifdef([AM_PROG_AR],[
  AM_PROG_AR
  ## ...
],[
  AC_CHECK_TOOL([AR],[ar])dnl
])dnl
m4_ifdef([AC_PROG_AR],[
  if test "x${AR}" = "x" && test "x${ac_cv_prog_AR}" = "x"; then
    test -z "${AR}" || unset AR
    test -z "${ac_cv_prog_AR}" || unset ac_cv_prog_AR
    AC_PROG_AR
  fi
],[
  test -x "`which ar`" && echo "ar is `which ar`"
])dnl
AM_PROG_AS
if test "x${CC}" = "x"; then
  test -z "${CC}" || unset CC
  AC_PROG_CC
else
  test -n "${CC}" && export CC
fi
m4_ifdef([AM_PROG_CC_C_O],[
  AM_PROG_CC_C_O
])dnl
AC_PROG_CPP
AC_SUBST([CPPFLAGS])dnl
AC_PROG_CXX
  ## keep autoscan happy:
if test "x${RANLIB}" = "x"; then
  test -z "${RANLIB}" || unset RANLIB
  AC_PROG_RANLIB
else
  test -n "${RANLIB}" && export RANLIB
  AC_SUBST([RANLIB])dnl
  ## (keep from being empty)
fi
AC_PATH_PROGS([RM_PROG],[rm grm],[rm])
AC_SUBST([RM_PROG])dnl
LT_INIT([win32-dll])
LT_LANG([C++])dnl

# Checks for libraries.
AC_CHECK_LIB([c],[printf])
AC_CHECK_LIB([z],[gzclose])
AC_CHECK_LIB([z],[gzdopen])
AC_CHECK_LIB([z],[gzopen])
AC_CHECK_LIB([z],[gzread])
AC_CHECK_LIB([z],[gzsetparams])
AC_CHECK_LIB([z],[gzwrite])
AC_CHECK_LIB([eigen],[eigenDiagMat],[],[
  AC_CHECK_LIB([eigen],[main])
])
AC_CHECK_LIB([crt0.o],[main])
AC_SEARCH_LIBS([std::locale::locale],[stdc++ c++ msvc++ stdlocale locale intl])
AC_CHECK_LIB([stdc++],[main],[],[
  AC_CHECK_LIB([c++],[main])dnl
])dnl

if test "x${LIBS}" != "x"; then
  test -n "${LIBS}" && export LIBS="`echo ${LIBS} | tr \  \\n | uniq | tr \\n \ `"
fi

# AX_* macros are from the autoconf-archive package
AX_CXX_COMPLEX_MATH_IN_NAMESPACE_STD
AX_CXX_HAVE_NUMERIC_LIMITS
AX_CXX_HAVE_COMPLEX_MATH1
AX_CXX_HAVE_COMPLEX_MATH2
AX_CXX_HAVE_COMPLEX
AX_CXX_HAVE_VECTOR_AT
AX_CXX_HEADER_STDCXX_98
AX_CXX_HEADER_STDCXX_TR1
AX_CHECK_GNU_MAKE

# This is for the Eigen subdir:
CMAKE_FIND_BINARY
AS_IF([test "x${CMAKE_BINARY}" != "x"],[
  test -n "${CMAKE_BINARY}" && test -e "${CMAKE_BINARY}" && test -x "${CMAKE_BINARY}" && echo "${CMAKE_BINARY} seems ok"
  CMAKE_FIND_PACKAGE([EIGEN],[CXX],[GNU],[],[],[
    PKG_CHECK_MODULES([EIGEN],[eigen],[
      if test "x${EIGEN_CFLAGS}" != "x"; then
        test -n "${EIGEN_CFLAGS}" && echo "EIGEN_CFLAGS are ${EIGEN_CFLAGS}"
        if test "x${CFLAGS}" = "x"; then
          test -z "${CFLAGS}" && export CFLAGS="${EIGEN_CFLAGS}"
        elif test "x${CPPFLAGS}" = "x"; then
          test -z "${CPPFLAGS}" && export CPPFLAGS="${EIGEN_CFLAGS}"
        else
          test -n "${CPPFLAGS}" && export CPPFLAGS="${CPPFLAGS} ${EIGEN_CFLAGS}"
        fi
      fi
      if test "x${EIGEN_LIBS}" != "x"; then
        test -n "${EIGEN_LIBS}" && echo "EIGEN_LIBS are ${EIGEN_LIBS}"
        if test "x${LIBS}" = "x"; then
          test -z "${LIBS}" && export LIBS="${EIGEN_LIBS}"
        else
          test -n "${LIBS}" && export LIBS="${LIBS} ${EIGEN_LIBS}"
        fi
      fi
    ],[
      dnl# Avoid eigen2, as that is too old. Just do eigen3 instead:
      PKG_CHECK_MODULES([EIGEN3],[eigen3],[
        if test "x${EIGEN3_CFLAGS}" != "x"; then
          test -n "${EIGEN3_CFLAGS}" && echo "EIGEN3_CFLAGS are ${EIGEN3_CFLAGS}"
          if test "x${CFLAGS}" = "x"; then
            test -z "${CFLAGS}" && export CFLAGS="${EIGEN3_CFLAGS}"
          elif test "x${CPPFLAGS}" = "x"; then
            test -z "${CPPFLAGS}" && export CPPFLAGS="${EIGEN3_CFLAGS}"
          else
            test -n "${CPPFLAGS}" && export CPPFLAGS="${CPPFLAGS} ${EIGEN3_CFLAGS}"
          fi
        fi
        if test "x${EIGEN3_LIBS}" != "x"; then
          test -n "${EIGEN3_LIBS}" && echo "EIGEN3_LIBS are ${EIGEN3_LIBS}"
          if test "x${LIBS}" = "x"; then
            test -z "${LIBS}" && export LIBS="${EIGEN3_LIBS}"
          else
            test -n "${LIBS}" && export LIBS="${LIBS} ${EIGEN3_LIBS}"
          fi
        fi
      ])dnl
    ])dnl
  ])dnl
])dnl

# Checks for header files.
# Autoupdate added the following line to ensure that your configure
# script's behavior did not change.  It is probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_CHECK_HEADERS([Eigen/Dense src/Core/util/DisableStupidWarnings.h \
                  unsupported/Eigen/SparseExtra])dnl

if test "x${ac_cv_header_stdio_h}" = "x"; then
  test -z "${ac_cv_header_stdio_h}" || unset ac_cv_header_stdio_h
  AC_CHECK_HEADERS([stdio.h])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_PROTOTYPES
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([atoi abs pow reverse sqrt seed_str stable_sort str_strm time])dnl

# Output.
AC_MSG_NOTICE([Done configuring; if either of the Makefiles generated by this script (Makefile or Makefile_orig) fail to let you build, try using the hardcoded Makefile instead: Makefile_bak])
AC_SUBST([LDFLAGS])dnl
AC_SUBST([LIB])dnl
AC_SUBST([HDR])dnl

AC_CONFIG_FILES([Makefile Makefile_orig])
AC_OUTPUT
