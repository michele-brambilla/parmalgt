#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([ParmaLGT], [0.11], [])
AM_INIT_AUTOMAKE([1.11])
AC_CONFIG_SRCDIR([Background.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
AC_LANG_CPLUSPLUS

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lfftw3':
AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lgsl':
# AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
# AC_CHECK_LIB([gslcblas], [main])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h inttypes.h limits.h stdlib.h string.h sys/time.h unistd.h])

AC_ARG_WITH([gtest-srcdir], [AS_HELP_STRING([--with-gtest-srcdir=DIR],
 [Path to Googletest.])], [], [])

AC_SUBST([gtest_srcdir],["$with_gtest_srcdir"])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_UINT32_T
# AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([bzero clock_gettime gettimeofday memset pow sqrt])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_FILES([Makefile testing/Makefile Doxyfile])
# AC_CONFIG_FILES([SU2/makefile
#                  dSFMT/Makefile
#                  makefile
#                  perform/makefile
#                  renorm/makefile
#                  spectr/makefile
#                  tools/makefile])
# AC_CONFIG_SUBDIRS([gtest
#                    testing])
AC_OUTPUT
