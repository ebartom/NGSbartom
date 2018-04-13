#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([pear], [0.9.10], [Tomas.Flouri@h-its.org])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/pear-pt.c])
AC_CONFIG_HEADERS([src/config.h])

if test -z "$CFLAGS"; then
    CFLAGS='-O3'
fi

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove pow strtol])


have_pthreads=no
have_bzlib=no
have_zlib=no
AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthreads=yes])

AC_ARG_ENABLE(bzip2, AS_HELP_STRING([--disable-bzip2], [Disable bzip2 support]))
AS_IF([test "x$enable_bzip2" != "xno"], [
  AC_SEARCH_LIBS([BZ2_bzCompress], [bz2], [have_bzlib=yes])
])

AC_ARG_ENABLE(zlib, AS_HELP_STRING([--disable-zlib], [Disable zlib support]))
AS_IF([test "x$enable_zlib" != "xno"], [
  AC_SEARCH_LIBS([zlibVersion], [z], [have_zlib=yes])
])

if test "x${have_pthreads}" = "xyes"; then
  AC_CHECK_HEADERS([pthread.h], [], [have_pthreads=no])
fi

if test "x${have_bzlib}" = "xyes"; then
  AC_CHECK_HEADERS([bzlib.h], [], [have_bzlib=no])
fi

if test "x${have_zlib}" = "xyes"; then
  AC_CHECK_HEADERS([zlib.h], [], [have_zlib=no])
fi

AM_CONDITIONAL(HAVE_BZLIB, test "x${have_bzlib}" = "xyes")
AM_CONDITIONAL(HAVE_ZLIB, test "x${have_zlib}" = "xyes")
AM_CONDITIONAL(HAVE_PTHREADS, test "x${have_pthreads}" = "xyes")
AM_PROG_CC_C_O

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile])

AC_OUTPUT
