# $Id: cm_compilers.m4,v 1.1.1.1 2012-09-11 06:58:21 jwp Exp $

#****f* cm_macros/cm_compilers.m4 *
#
# NAME
#    cm_compilers.m4 - Autoconf / automake macros for various compilers.
#
# SYNOPSIS
#    Place the cm_compilers.m4 file in a subdirectory 'm4' of the top level
#    directory of the project, and add
#
#       ACLOCAL_AMFLAGS  = -I m4
#
#    to the toplevel Makefile.am
#
#    To generate the initial aclocal.m4, run
#
#       aclocal -I m4
#
#    and the usual sequence of automake and autoconf thereafter.
# 
# DESCRIPTION
#    The file cm_compilers.m4 provides a small number of useful autoconf /
#    automake macros for working with Fortran 9x projects.
#
# NOTES
#    cm_compilers.m4 requires (at least) autoconf 2.59 and a patched version
#    of automake 1.8.2.
#
# SEE ALSO
#    cm_prog_cc_vendor
#    cm_prog_cxx_vendor
#    cm_prog_f77_vendor
#    cm_prog_fc_vendor
#
#    cm_fc_srcext
#
#    cm_fortran.m4
#    cm_libraries.m4
#
# AUTHOR
#    C. Marquardt, West Hill, UK        <christian@marquardt.fsnet.co.uk>
#
#****


#****m* cm_compilers.m4/cm_prog_cc_vendor *
#
# NAME
#    CM_PROG_CC_VENDOR - C compiler's vendor.
#
# SYNOPSIS
#    CM_PROG_CC_VENDOR(vendor list)
# 
# DESCRIPTION
#    This macro tries to find the vendor of the C (CC) compiler by comparing
#    it with a (given or defaulted) list of vendors, and sets a shell variable
#    accordingly.
#
# INPUTS
#    vendor list     List of vendors (seperated by blanks). The default list is
#                       'Intel Fujitsu MIPS AIX GCC'
#
# OUTPUT
#    CM_CC_VENDOR   Shell variable set to the matching entry in <vendor list>.
#
# NOTES
#    The macro relies on version information given by the compiler; some compilers
#    do not report a company name, but some other information. This is a list of
#    known identification strings that can be used, but do not agree with common
#    vendor names:
#
#       GCC               GNU compilers
#       MIPS, MIPSPro     Native compilers on SGI platforms
#       AIX, VisualAge    Native compilers on IBM platforms (xlC)
#
# EXAMPLE
#    To see if the C compiler is in the default list of compilers, use
#
#       CM_PROG_CC_VENDOR()
#
#    To check if the compiler is from SGI or IBM, use
#
#       CM_PROG_CC_VENDOR(MIPS AIX)
#
# SEE ALSO
#    cm_prog_cxx_vendor
#    cm_prog_f77_vendor
#    cm_prog_fc_vendor
#
#    cm_fc_srcext
#
#    cm_libraries.m4
#    cm_compilers.m4
#
# BUGS
#    If the compiler does not provide some vendor / version information when
#    called with -V, -v, -version, --version, -qversion or +version, it's vendor cannot
#    be determined with the current version of the macro.
#
#    At least one C compiler (under HP-UX) does not provide any version
#    information if called without a file to compile. Currently, this macro
#    cannot detect this compiler, as no compilation is undertaken for
#    performance reasons.
#
#****

# 1. CM_PROG_CC_VENDOR(VENDOR ...)
# -------------------------------

AC_DEFUN([CM_PROG_CC_VENDOR],
[AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_CACHE_CHECK([for C compiler vendor],
               [CM_CC_VENDOR],
[dnl
if test x"$1" = x ; then
   vendors="Intel Fujitsu MIPS AIX GCC PathScale Portland Cray Sun SX-6 NEC IBM Absoft Open64 GNU"
else
   vendors="$1"
fi

CM_CC_VENDOR="Unknown"

for ac_version in -qversion -V -v -version --version +version
do
   case $host in
   *-ibm-aix*)
     $CC $ac_version 2>&1 > conftest.txt
     cm_cc_v_output=`head -c 2048 conftest.txt`
     rm -f conftest.txt
     ;;
   *)
     cm_cc_v_output=`$CC $ac_version 2>&1`
     ;;
   esac
   for cm_cc_vendor_string in $vendors
   do
      cm_cc_v_string=`echo $cm_cc_v_output | grep $cm_cc_vendor_string`
      if test x"$cm_cc_v_string" != x ; then
         CM_CC_VENDOR=$cm_cc_vendor_string
         break 2
      fi
   done
done
])
])


#****m* cm_compilers.m4/cm_prog_cxx_vendor *
#
# NAME
#    CM_PROG_CXX_VENDOR - C++ compiler's vendor.
#
# SYNOPSIS
#    CM_PROG_CXX_VENDOR(vendor list)
# 
# DESCRIPTION
#    This macro tries to find the vendor of the C++ (CXX) compiler by comparing
#    it with a (given or defaulted) list of vendors, and sets a shell variable
#    accordingly.
#
# INPUTS
#    vendor list     List of vendors (seperated by blanks). The default list is
#                       'Intel Fujitsu MIPS AIX GCC'
#
# OUTPUT
#    CM_CXX_VENDOR   Shell variable set to the matching entry in <vendor list>.
#
# NOTES
#    The macro relies on version information given by the compiler; some compilers
#    do not report a company name, but some other information. This is a list of
#    known identification strings that can be used, but do not agree with common
#    vendor names:
#
#       GCC               GNU compilers
#       MIPS, MIPSPro     Native compilers on SGI platforms
#       AIX, VisualAge    Native compilers on IBM platforms (xlC)
#
# EXAMPLE
#    To see if the C++ compiler is in the default list of compilers, use
#
#       CM_PROG_CXX_VENDOR()
#
#    To check if the compiler is from SGI or IBM, use
#
#       CM_PROG_CXX_VENDOR(MIPS AIX)
#
# SEE ALSO
#    cm_prog_cc_vendor
#    cm_prog_f77_vendor
#    cm_prog_fc_vendor
#
#    cm_fc_srcext
#
#    cm_libraries.m4
#    cm_compilers.m4
#
# BUGS
#    If the compiler does not provide some vendor / version information when
#    called with -V, -v, -version, --version or +version, it's vendor cannot
#    be determined with the current version of the macro.
#
#    I have not tested this macro properly (I'm not a C++ programmer, really).
#
#****

# 2. CM_PROG_CXX_VENDOR(VENDOR ...)
# --------------------------------

AC_DEFUN([CM_PROG_CXX_VENDOR],
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_CACHE_CHECK([for C++ compiler vendor],
               [CM_CXX_VENDOR],
[dnl
if test x"$1" = x ; then
   vendors="Intel Fujitsu MIPS AIX GCC PathScale Portland Cray Sun SX-6 NEC IBM Absoft Open64 GNU"
else
   vendors="$1"
fi

CM_CXX_VENDOR="Unknown"

for ac_version in -qversion -V -v -version --version +version
do
   case $host in
   *-ibm-aix*)
     $CXX $ac_version 2>&1 > conftest.txt
     cm_cxx_v_output=`head -c 2048 conftest.txt`
     rm -f conftest.txt
     ;;
   *)
     cm_cxx_v_output=`$CXX $ac_version 2>&1`
     ;;
   esac
   for cm_cxx_vendor_string in $vendors
   do
      cm_cxx_v_string=`echo $cm_cxx_v_output | grep $cm_cxx_vendor_string`
      if test x"$cm_cxx_v_string" != x ; then
         CM_CXX_VENDOR=$cm_cxx_vendor_string
         break 2
      fi
   done
done
])
])


#****m* cm_compilers.m4/cm_prog_f77_vendor *
#
# NAME
#    CM_PROG_F77_VENDOR - Fortran compiler's vendor.
#
# SYNOPSIS
#    CM_PROG_F77_VENDOR(vendor list)
# 
# DESCRIPTION
#    This macro tries to find the vendor of the Fortran 77 (F77) compiler by
#    comparing it with a (given or defaulted) list of vendors, and sets a shell
#    variable accordingly.
#
# INPUTS
#    vendor list     List of vendors (seperated by blanks). The default list is
#                       'Intel NAG Fujitsu Portland HP MIPS AIX GCC'
#
# OUTPUT
#    CM_F77_VENDOR   Shell variable set to the matching entry in <vendor list>.
#
# NOTES
#    The macro relies on version information given by the compiler; some compilers
#    do not report a company name, but some other information. This is a list of
#    known identification strings that can be used, but do not agree with common
#    vendor names:
#
#       GCC               GNU compilers
#       MIPS, MIPSPro     Native compilers on SGI platforms
#       AIX, XL           Native compilers on IBM platforms (xlf77, xlf)
#
# EXAMPLE
#    To see if the Fortran 77 compiler is in the default list of compilers, use
#
#       CM_PROG_F77_VENDOR()
#
#    To check if the compiler is from SGI or IBM, use
#
#       CM_PROG_F77_VENDOR(MIPS AIX)
#
# SEE ALSO
#    cm_prog_fc_vendor
#    cm_prog_cc_vendor
#    cm_prog_cxx_vendor
#
#    cm_fc_srcext
#
#    cm_libraries.m4
#    cm_compilers.m4
#
# BUGS
#    If the compiler does not provide some vendor / version information when
#    called with -V, -v, -version, --version or +version, it's vendor cannot
#    be determined with the current version of the macro.
#
#****

# 3. CM_PROG_F77_VENDOR(VENDOR ...)
# --------------------------------

AC_DEFUN([CM_PROG_F77_VENDOR],
[AC_REQUIRE([AC_PROG_F77])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_CACHE_CHECK([for Fortran 77 compiler vendor],
               [CM_F77_VENDOR],
[dnl
if test x"$1" = x ; then
   vendors="Intel NAG Fujitsu MIPS AIX GCC PathScale Portland Cray Sun SX-6 NEC IBM Absoft Open64 HP GNU"
else
   vendors="$1"
fi

CM_F77_VENDOR="Unknown"

for ac_version in -qversion -V -v -version --version +version
do
   case $host in
   *-ibm-aix*)
     $F77 $ac_version 2>&1 > conftest.txt
     cm_f77_v_output=`head -c 2048 conftest.txt`
     rm -f conftest.txt
     ;;
   *)
     cm_f77_v_output=`$F77 $ac_version 2>&1`
     ;;
   esac
   for cm_f77_vendor_string in $vendors
   do
      cm_f77_v_string=`echo $cm_f77_v_output | grep $cm_f77_vendor_string`
      if test x"$cm_f77_v_string" != x ; then
         CM_F77_VENDOR=$cm_f77_vendor_string
         break 2
      fi
   done
done
])
])


#****m* cm_compilers.m4/cm_prog_fc_vendor *
#
# NAME
#    CM_PROG_FC_VENDOR - Fortran compiler's vendor.
#
# SYNOPSIS
#    CM_PROG_FC_VENDOR(vendor list)
# 
# DESCRIPTION
#    This macro tries to find the vendor of the Fortran (FC) compiler by comparing
#    it with a (given or defaulted) list of vendors, and sets a shell variable
#    accordingly.
#
# INPUTS
#    vendor list     List of vendors (seperated by blanks). The default list is
#                       'Intel NAG Fujitsu Portland HP MIPS AIX GCC'
#
# OUTPUT
#    CM_FC_VENDOR    Shell variable set to the matching entry in <vendor list>.
#
# NOTES
#    The macro relies on version information given by the compiler; some compilers
#    do not report a company name, but some other information. This is a list of
#    known identification strings that can be used, but do not agree with common
#    vendor names:
#
#       GCC               GNU compilers
#       MIPS, MIPSPro     Native compilers on SGI platforms
#       AIX, XL           Native compilers on IBM platforms (xlf95, xlf90, xlf)
#
# EXAMPLE
#    To see if the Fortran compiler is in the default list of compilers, use
#
#       CM_PROG_FC_VENDOR()
#
#    To check if the compiler is from SGI or IBM, use
#
#       CM_PROG_FC_VENDOR(MIPS AIX)
#
# SEE ALSO
#    cm_prog_f77_vendor
#    cm_prog_cc_vendor
#    cm_prog_cxx_vendor
#
#    cm_fc_srcext
#
#    cm_libraries.m4
#    cm_compilers.m4
#
# BUGS
#    If the compiler does not provide some vendor / version information when
#    called with -V, -v, -version, --version or +version, it's vendor cannot
#    be determined with the current version of the macro.
#
#****

# 4. CM_PROG_FC_VENDOR(VENDOR ...)
# --------------------------------

AC_DEFUN([CM_PROG_FC_VENDOR],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_CACHE_CHECK([for Fortran compiler vendor],
               [CM_FC_VENDOR],
[dnl
if test x"$1" = x ; then
   vendors="Intel NAG Fujitsu MIPS AIX GCC PathScale Portland Cray Sun SX-6 NEC IBM Absoft Open64 HP GNU"
else
   vendors="$1"
fi

CM_FC_VENDOR="Unknown"

for ac_version in -qversion -V -v -version --version +version
do
   case $host in
   *-ibm-aix*)
     $FC $ac_version 2>&1 > conftest.txt
     cm_fc_v_output=`head -c 2048 conftest.txt`
     rm -f conftest.txt
     ;;
   *)
     cm_fc_v_output=`$FC $ac_version 2>&1`
     ;;
   esac
   for cm_fc_vendor_string in $vendors
   do
      cm_fc_v_string=`echo $cm_fc_v_output | grep $cm_fc_vendor_string`
      if test x"$cm_fc_v_string" != x ; then
         CM_FC_VENDOR=$cm_fc_vendor_string
         break 2
      fi
   done
done
])
])


#****m* cm_compilers.m4/cm_fc_srcext *
#
# NAME
#    CM_FC_SRCEXT - Alternative version of AC_FC_SRCEXT.
#
# SYNOPSIS
#    CM_FC_SRCEXT(ext, [code], [action-if-success], [action-if-failure])
# 
# DESCRIPTION
#    This macro provides an extension of the AC_FC_SRCEXT macro. It adds the
#    optional CODE which allows to check if the compiler is able to handle a
#    certain language construct. Examples are preprocessor directives or
#    dialect specific features.
#
# INPUTS
#    ext
#    code
#    action-if-success
#    action-if-failure
#
# OUTPUT
#    FCFLAGS_<ext>   Compiler options required to compile .<ext> files.
#
# NOTES
#    For ordinary extensions like f90 the modified FCFLAGS are currently
#    needed for IBM's xlf* and Intel's ifc.  Unfortunately,  xlf* will only
#    take flags to recognize one extension at a time, so if the user wants
#    to compile multiple extensions (.f90 and .f95, say), she will need to
#    use the FCFLAGS_f90 and FCFLAGS_f95 individually rather than just adding
#    them all to FCFLAGS. Also, for Intel's ifc compiler (which does not
#    accept .f95 by default in some versions), the $FCFLAGS_<ext> variable
#    must go immediately before the source file on the command line, unlike
#    other $FCFLAGS.
#
#    Note that the usual syntax -DDEFINITION=<something> does not work with
#    IBM's xlf* compilers; therefore, the Makefile variable DEFS needs to
#    be changed if it is to be applied with preprocessed Fortran source code.
#    See the first example below.
#
#    Under Irix, the compiler does not by default expand preprocessor macro
#    definitions. If this is required in the application, the source code
#    for the test should be written in a way that it relies on macro
#    expansion; see the second example.
#
# EXAMPLE
#    In order to find out if the compiler requires an additional flag for .f90
#    or .f95 files, try
#
#       CM_FC_SRCEXT(f90)
#       CM_FC_SRCEXT(f95)
#
#    To check if the compiler supports a Fortran 95 construct which is not
#    part of Fortran 90, try
#
#       CM_FC_SRCEXT(f95, [real, dimension(:), pointer :: ptr = null()])
#
#    IBM's xlf* family of compilers require a different set of command line
#    options for standard and preprocessed options. To get the correct option
#    to compile Fortran code which requires preprocessing
#
#       CM_FC_SRCEXT(F90, [#define choke me])
#
#    As noted above, the xlf family of compilers requires a change in the
#    definition of preprocessor variables. With GNU make, this can be achived
#    by using a target specific variable definition for those objects that
#    depend on preprocessed Fortran files. In the corresponding Makefile.am,
#    this requires something like
#
#       if AIX_XL
#          sub.o : DEFS = $(DEFS:-D%=-WF,-D%)
#       endif
#
#    where AIX_XL has been defined as an automake conditional in configure.ac,
#    e.g. via
#
#       AM_CONDITIONAL(AIX_XL, test x$FCFLAGS_F90 = x-qsuffix=cpp=F90)
#
#    Note that this solution probably requires the use of GNU make, as the
#    native make may note support target specific variables and/or pattern
#    substitution rules.
#
#    If the expansion of macros is required, the CM_FC_SRCEXT test should
#    be written as something like
#
#       CM_FC_SRCEXT(F90, [
#          implicit none
#          real :: x
#       #  define choke_me 1.0
#          x = choke_me
#       ])
#
#    to enforce checking that the expansion is actually done at configuration
#    time. This is required as some compiler (notably the native Irix one) do
#    not expand preprocessor macros by default.
#
# SEE ALSO
#    cm_prog_f77_vendor
#    cm_prog_fc_vendor
#    cm_prog_cc_vendor
#    cm_prog_cxx_vendor
#
#    cm_libraries.m4
#    cm_compilers.m4
#
#****

# 5. CM_FC_SRCEXT(EXT, [CODE], [ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
# ----------------------------------------------------------------------

AC_DEFUN([CM_FC_SRCEXT],
[AC_LANG_ASSERT(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag to compile .$1 files],
                cm_cv_fc_srcext_$1,
[ac_ext=$1
cm_fc_srcext_FCFLAGS_SRCEXT_save=$FCFLAGS_SRCEXT
FCFLAGS_SRCEXT=""
cm_cv_fc_srcext_$1=unknown
for cm_flag in none -Tf -qsuffix=f=$1 -qsuffix=cpp=$1 -macro_expand +cpp=yes; do
  test "x$cm_flag" != xnone && FCFLAGS_SRCEXT="$cm_flag"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[$2])], [cm_cv_fc_srcext_$1=$cm_flag; break])
done
rm -f conftest.$ac_objext conftest.$1
FCFLAGS_SRCEXT=$cm_fc_srcext_FCFLAGS_SRCEXT_save
])
if test "x$cm_cv_fc_srcext_$1" = xunknown; then
  m4_default([$4],[AC_MSG_ERROR([Fortran could not compile .$1 files])])
else
  FC_SRCEXT=$1
  if test "x$cm_cv_fc_srcext_$1" = xnone; then
    FCFLAGS_SRCEXT=""
    FCFLAGS_[]$1[]=""
  else
    FCFLAGS_SRCEXT=$cm_cv_fc_srcext_$1
    FCFLAGS_[]$1[]=$cm_cv_fc_srcext_$1
  fi
  AC_SUBST(FCFLAGS_[]$1)
  $3
fi
])# CM_FC_SRCEXT
