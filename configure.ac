#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([gcta], [0.93.9], [jian.yang@uq.edu.au])
AC_CONFIG_SRCDIR([cdflib.h])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AM_DEP_TRACK
AM_MAINTAINER_MODE
AM_MAKE_INCLUDE
AM_C_PROTOTYPES
AM_WITH_REGEX

# Checks for programs.
AM_PROG_AS
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_SUBST(CPPFLAGS)
AC_PROG_CXX
AC_PATH_PROG([RM],[rm])
AC_SUBST(RM)
AM_PROG_LIBTOOL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lc':
AC_CHECK_LIB([c], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])
# FIXME: Replace `main' with a function in `-leigen':
AC_CHECK_LIB([eigen], [main])
# FIXME: Replace `main' with a function in `-lcrt0.o':
AC_CHECK_LIB([crt0.o], [main])

# AX_* macros are from the autoconf-archive package
AX_CXX_COMPLEX_MATH_IN_NAMESPACE_STD
AX_CXX_HAVE_NUMERIC_LIMITS
AX_CXX_HAVE_COMPLEX_MATH1
AX_CXX_HAVE_COMPLEX_MATH2
AX_CXX_HAVE_COMPLEX
AX_CXX_HAVE_VECTOR_AT
AX_CXX_HEADER_STDCXX_98
AX_CXX_HEADER_STDCXX_TR1

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([atoi abs pow reverse sqrt seed_str stable_sort str_strm time])

AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([Makefile Makefile.orig])
AC_OUTPUT
