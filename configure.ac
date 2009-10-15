dnl Process this file with autoconf to produce a configure script.

AC_INIT([domemaster], [0.4.1], [reunanen@ftml.net])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([admin])
AC_CANONICAL_SYSTEM 

AM_INIT_AUTOMAKE() 

dnl AC_CONFIG_MACRO_DIR([macros])

AC_PROG_CXX
dnl AC_ISC_POSIX
dnl AC_PROG_CPP
dnl AM_PROG_CC_C_O
AC_PROG_INSTALL
 
dnl ====== FreeImage =================================
AC_CHECK_LIB(freeimage, FreeImage_Load, [LIBS="$LIBS -lfreeimage -lm"],
  [AC_MSG_ERROR([Cannot find libfreeimage])]
)

dnl ====== pthread ===================================
AC_CHECK_LIB(pthread, pthread_create, [LIBS="$LIBS -lpthread"],
  [AC_MSG_ERROR([Cannot find libpthread])]
)

dnl ====== gtk =================================
PKG_CHECK_MODULES(GTK, [gtk+-2.0], [ AC_DEFINE([HAVE_GTK], [1], [Description]) 
 CXXFLAGS="$CXXFLAGS $GTK_CFLAGS"
 LIBS="$LIBS $GTK_LIBS"],
 [AC_MSG_WARN([Cannot use GTK])])

dnl ====== Needed in Solaris/SunOS =============
AC_CHECK_LIB(socket, socket, [LIBS="$LIBS -lsocket"])
AC_CHECK_FUNCS(gethostbyname, , [AC_CHECK_LIB(nsl, gethostbyname), , LIBS="$LIBS -lnsl"])

dnl ====== Output ==============================
AC_OUTPUT(Makefile src/Makefile)

dnl ====== Print configuration =================
echo "

Configuration:

  Detected OS:      ${host_os}
  Software version: ${VERSION}
  Install path:     ${prefix}
  Compiler:         ${CXX}
  Compiler flags:   ${CXXFLAGS}
  Linker flags:     ${LIBS}

"
