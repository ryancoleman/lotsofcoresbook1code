# Copyright (c) 2014, Per Berg and Jacob Weismann Poulsen, DMI
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

AC_PREREQ(2.59)
AC_INIT([hbm], [hbm-micro-tflow-1.6], jwp@dmi.dk)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src/main.f90])
CM_PROG_FC_VENDOR

BYTESWAPIO=""
EXTRA_CPPFLAGS=""

AC_ARG_ENABLE(openmp,
  AC_HELP_STRING([--enable-openmp], [try to build with openMP support]),
  [case $enableval in
    yes | true) want_openmp=yes ;;
    *)          want_openmp=no ;;
  esac],
  [want_openmp=no]
)

AC_ARG_ENABLE(mpi,
  AC_HELP_STRING([--enable-mpi], [try to build with MPI support]),
  [case $enableval in
    yes | true) want_mpi=yes ;;
    *)          want_mpi=no ;;
  esac],
  [want_mpi=no]
)

AC_ARG_ENABLE(contiguous,
  AC_HELP_STRING([--enable-contiguous], [try to build with F2008 contiguous attribute support]),
  [case $enableval in
    yes | true) want_contiguous=yes ;;
    *)          want_contiguous=no ;;
  esac],
  [want_contiguous=no]
)
if test x"$want_contiguous" == x"yes"; then
  cp src/optional/contiguous/*.inc src/
else
  cp src/optional/no-contiguous/*.inc src/
fi

if test x"$want_mpi" == x"yes"; then
  AX_PROG_FC_MPI([],[EXTRA_CPPFLAGS="-DMPI"],[AC_MSG_ERROR([MPI support is missing on the local system and this is required when using option --enable-mpi])])
fi

AC_ARG_ENABLE(openacc,
  AC_HELP_STRING([--enable-openacc], [try to build with openacc support]),
  [case $enableval in
    yes | true) want_openacc=yes ;;
    *)          want_openacc=no ;;
  esac],
  [want_openacc=no]
)

# test exit
AC_LANG_PUSH([Fortran])
AC_MSG_CHECKING(for external exit() function)
AC_LINK_IFELSE([AC_LANG_SOURCE([
      program foo
      call exit(0)
      end program foo])], [ac_have_exit=yes], [ac_have_exit=no])
AC_MSG_RESULT($ac_have_exit)
AC_LANG_POP([Fortran])
if test x"$ac_have_exit" == x"yes"; then
  cp src/optional/exits_exit_external.f90 src/exits_exit.f90
else
  cp src/optional/exits_noexit.f90 src/exits_exit.f90
fi

# test flush
cp src/optional/io_subs_noflush.f90 src/io_subs_flush.f90
if test x"$ac_have_flush2003" == x"no"; then
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING(for external flush() function)
  AC_LINK_IFELSE([AC_LANG_SOURCE([
      program foo
      call flush(100)
      end program foo])], [ac_have_flush=yes], [ac_have_flush=no])
  AC_MSG_RESULT($ac_have_flush)
  AC_LANG_POP([Fortran])
fi
AC_LANG_PUSH([Fortran])
AC_MSG_CHECKING(for fortran2003 flush() function)
AC_LINK_IFELSE([AC_LANG_SOURCE([
      program foo
      integer :: lun
      flush(lun)
      end program foo])], [ac_have_flush2003=yes], [ac_have_flush2003=no])
AC_MSG_RESULT($ac_have_flush2003)
AC_LANG_POP([Fortran])
if test x"$ac_have_flush" == x"yes"; then
	cp src/optional/io_subs_flush_external.f90 src/io_subs_flush.f90
fi
if test x"$ac_have_flush2003" == x"yes"; then
  cp src/optional/io_subs_flush_2003.f90 src/io_subs_flush.f90
fi

if test x"$want_openmp" == x"yes"; then
  if test x"$want_openacc" == x"yes"; then
	AC_MSG_ERROR([We support openMP and openACC but you cannot combine the two])
  fi
fi

EXTRA_FFLAGS="$EXTRA_CPPFLAGS"

case x"$CM_FC_VENDOR" in
    xCray)
        FCCPP=""
		AUTODBL="-sreal64 -DREAL_8"
        BYTESWAPIO="-hbyteswapio"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -O omp"
		else
			EXTRA_FFLAGS="$EXTRA_FFLAGS -O noomp"
		fi
		if test x"$want_openacc" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -O acc"
		fi
        FCFLAGS="$EXTRA_FFLAGS $FCFLAGS -em"
        ;;
	xAIX|xIBM)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP=""
		EXTRA_FFLAGS=""
		AUTODBL="-qautodbl=dbl4 -qdpc=e"
        BYTESWAPIO=""
		if test x"$want_openmp" == x"yes"; then
            EXTRA_FFLAGS="$EXTRA_FFLAGS -qsmp=omp"
        fi
		if test x"$want_mpi" == x"yes"; then
            EXTRA_FFLAGS="$EXTRA_FFLAGS -WF,-DMPI"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''-q64'' ''$EXTRA_FFLAGS'' $FCFLAGS)
#		AC_CHECK_LIB(mass, sin, , AC_MSG_WARN(libmass missing))
		;;
	xNAG)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-fpp"
		AUTODBL="-r8 -DREAL_8"
        BYTESWAPIO="-convert=BIG_ENDIAN"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -D_OPENMP -openmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''-kind=byte'' ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xGCC|xGNU)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-cpp"
		AUTODBL="-fdefault-real-8 -fdefault-double-8 -DREAL_8"
		if test x"$want_endian" == x"yes"; then
            BYTESWAPIO="-fconvert=big-endian"
		fi
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -fopenmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xIntel)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-cpp"
		AUTODBL="-r8"
        BYTESWAPIO="-convert big_endian"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -openmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xPortland)
        FCCPP="-Mpreprocess"
		AUTODBL="-r8"
        BYTESWAPIO="-byteswapio"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -mp"
		fi
		if test x"$want_openacc" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -ta=nvidia"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xPathScale)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-ftpp"
		AUTODBL="-r8"
        BYTESWAPIO="-byteswapio"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -openmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xOpen64)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-ftpp"
		AUTODBL="-r8"
        BYTESWAPIO="-byteswapio"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -openmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS,  ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
	xSun)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-fpp"
		AUTODBL="-xtypemap=real:64,double:64"
        BYTESWAPIO="-xfilebyteorder=big8:%all"
		if test x"$want_openmp" == x"yes"; then
			EXTRA_FFLAGS="$EXTRA_FFLAGS -xopenmp"
		fi
		AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
		;;
    xFujitsu)
		if test x"$want_openacc" == x"yes"; then
		    AC_MSG_WARN(This compiler does not support openACC)
		fi
        FCCPP="-Cpp"
		AUTODBL="--dbl"
        BYTESWAPIO="-byteswapio"
        if test x"$want_openmp" == x"yes"; then
            EXTRA_FFLAGS="$EXTRA_FFLAGS --openmp"
        fi
        AC_UNIQIFY_LIST(FCFLAGS, ''$EXTRA_FFLAGS'' $FCFLAGS)
        ;;
	*)
		AC_MSG_ERROR(Compiler is unknown to the hbm configuration - please consult
		the author for information on required compiler flags)
		;;
esac

AC_SUBST(BYTESWAPIO)
AC_SUBST(FCCPP)
AC_SUBST(FCFLAGS)
AC_SUBST(ac_have_exit)
AC_SUBST(ac_have_flush2003)
AC_SUBST(ac_have_flush)
AC_SUBST(ac_want_netcdf_op)
AC_SUBST(want_openacc)
AC_SUBST(want_dynvertical)
AC_SUBST(want_dmiergom9)
AC_SUBST(want_dmiergom10)
AC_SUBST(want_dmiergom12)
AC_SUBST(want_bench)
AC_SUBST(want_toy)
AC_CONFIG_FILES([Makefile 
				Makefile.include
				Makefile.manual
                src/.objs.spec
                src/.dependencies
				 ])
AC_OUTPUT
