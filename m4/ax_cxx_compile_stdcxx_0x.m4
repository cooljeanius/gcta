# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_0x.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_0X
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++0x
#   standard.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 8

AU_ALIAS([AC_CXX_COMPILE_STDCXX_0X], [AX_CXX_COMPILE_STDCXX_0X])
AC_DEFUN([AX_CXX_COMPILE_STDCXX_0X], [
  AC_CACHE_CHECK(if g++ supports C++0x features without additional flags,
  ax_cv_cxx_compile_cxx0x_native,
  dnl# autoupdate says this is obsolete: Instead of using `AC_LANG',
  dnl# `AC_LANG_SAVE', and `AC_LANG_RESTORE', you should use `AC_LANG_PUSH'
  dnl# and `AC_LANG_POP'.
  AC_LANG_SAVE
  AC_LANG([C++])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);]], [[]])],[ax_cv_cxx_compile_cxx0x_native=yes],[ax_cv_cxx_compile_cxx0x_native=no])
  AC_LANG_POP([])
  ])

  AC_CACHE_CHECK(if g++ supports C++0x features with -std=c++0x,
  ax_cv_cxx_compile_cxx0x_cxx,
  dnl# autoupdate says this is obsolete: Instead of using `AC_LANG',
  dnl# `AC_LANG_SAVE', and `AC_LANG_RESTORE', you should use `AC_LANG_PUSH'
  dnl# and `AC_LANG_POP'
  AC_LANG_SAVE
  AC_LANG([C++])
  ac_save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -std=c++0x"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);]], [[]])],[ax_cv_cxx_compile_cxx0x_cxx=yes],[ax_cv_cxx_compile_cxx0x_cxx=no])
  CXXFLAGS="$ac_save_CXXFLAGS"
  AC_LANG_POP([])
  ])

  AC_CACHE_CHECK(if g++ supports C++0x features with -std=gnu++0x,
  [ax_cv_cxx_compile_cxx0x_gxx],
  [
  dnl# autoupdate says this is obsolete: Instead of using `AC_LANG',
  dnl# `AC_LANG_SAVE', and `AC_LANG_RESTORE', you should use `AC_LANG_PUSH'
  dnl# and `AC_LANG_POP'.
  AC_LANG_SAVE
  AC_LANG([C++])
  ac_save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -std=gnu++0x"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
  template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

    typedef check<check<bool>> right_angle_brackets;

    int a;
    decltype(a) b;

    typedef check<int> check_type;
    check_type c;
    check_type&& cr = static_cast<check_type&&>(c);]], [[]])],[ax_cv_cxx_compile_cxx0x_gxx=yes],[ax_cv_cxx_compile_cxx0x_gxx=no])
  CXXFLAGS="$ac_save_CXXFLAGS"
  AC_LANG_POP([])
  ])

  if test "$ax_cv_cxx_compile_cxx0x_native" = yes ||
     test "$ax_cv_cxx_compile_cxx0x_cxx" = yes ||
     test "$ax_cv_cxx_compile_cxx0x_gxx" = yes; then
    AC_DEFINE(HAVE_STDCXX_0X,,[Define if g++ supports C++0x features. ])
  fi
])
