#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(vcg, 0.1, ponchio@dm.unipi.it)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([../tetra_decimator/main.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lGL':
AC_CHECK_LIB([GL], [main])
# FIXME: Replace `main' with a function in `-lGLU':
AC_CHECK_LIB([GLU], [main])
# FIXME: Replace `main' with a function in `-lSDL':
AC_CHECK_LIB([SDL], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor ftruncate gettimeofday pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		nexus/Makefile
		metro/Makefile
		tri_decimator/Makefile
		shadevis/Makefile])

AC_OUTPUT
