#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([PRlgtTests], [0.1], [])
AM_INIT_AUTOMAKE([1.10 no-define foreign])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_LANG_CPLUSPLUS

# Checks for header files.

AC_CHECK_HEADER(array,[AC_DEFINE(HAVE_CXX0X, 1, 
                            [Define to 1 if c++11 is available])])
AC_CHECK_HEADER(tr1/array,[AC_DEFINE(HAVE_TR1, 1, 
                            [Define to 1 if TR1 is available])])


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_LIB(pthread, pthread_create)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
