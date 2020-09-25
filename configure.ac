dnl# configure.ac                                  -*- Autoconf -*-
dnl# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69.185-b502e])
AC_INIT([gcta],[0.93.9],[jian.yang@uq.edu.au])
AC_CONFIG_SRCDIR([cdflib.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.11.6])
AM_DEP_TRACK
AM_MAINTAINER_MODE
AM_MAKE_INCLUDE
m4_ifdef([AM_WITH_REGEX],[
  AM_WITH_REGEX
])dnl

# Checks for programs.
AM_PROG_AS
AC_PROG_CC
m4_ifdef([AM_PROG_CC_C_O],[
  AM_PROG_CC_C_O
])dnl
AC_PROG_CPP
AC_SUBST([CPPFLAGS])dnl
AC_PROG_CXX
AC_PATH_PROG([RM],[rm])
AC_SUBST([RM])dnl
LT_INIT([win32-dll])
LT_LANG([C++])dnl

# Checks for libraries.
AC_CHECK_LIB([c],[printf])
AC_CHECK_LIB([z],[gzread])
# FIXME: Replace `main' with a function in `-leigen':
AC_CHECK_LIB([eigen],[main])
# FIXME: Replace `main' with a function in `-lcrt0.o':
AC_CHECK_LIB([crt0.o],[main])dnl

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
CMAKE_FIND_PACKAGE([EIGEN],[CXX])dnl

# Checks for header files.
# Autoupdate added the following line to ensure that your configure
# script's behavior did not change.  It is probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_CHECK_HEADERS([stdio.h])dnl

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_PROTOTYPES
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([atoi abs pow reverse sqrt seed_str stable_sort str_strm time])

AC_SUBST([LDFLAGS])dnl
AC_SUBST([LIB])dnl
AC_SUBST([HDR])dnl

# Output.
AC_CONFIG_FILES([Makefile Makefile_orig])
AC_OUTPUT
