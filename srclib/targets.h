# "$Id$"
PVM_PATH = ${PVM_ROOT}/lib/${PVM_ARCH}

ifdef PVM_ROOT
   pvm_linker_line=-L${PVM_PATH} -lfpvm3 -lgpvm3 -lpvm3
   pvm_message=Building program with PVM support
else
   pvm_linker_line=
   pvm_message=Building program without PVM support
   # SnoopMLSL2.o: SnoopMLSL2.f90 $(OBJS) $(libdir)/libmls.a
endif

ifdef BLAS
ifdef LIB_BLAS
    blas_linker_line := $(shell ${REECHO} -dir ${BLAS} -prefixn=-l -lib ${LIB_BLAS})
    blas_linker_line := -L${BLAS} ${blas_linker_line}
else
   blas_linker_line=-L${BLAS} -lblas
endif
   blas_message=Building program with externally-supplied BLAS
else
   blas_linker_line=
   blas_message=Building program with internally-supplied BLAS
endif

ifdef LAPACK
ifdef LIB_LAPACK
    lapack_linker_line := $(shell ${REECHO} -dir ${LAPACK} -prefixn=-l -lib ${LIB_LAPACK})
    lapack_linker_line := -L${LAPACK} ${lapack_linker_line}
else
   lapack_linker_line=-L${LAPACK} -llapack
endif
   lapack_message=Building program with externally-supplied LAPACK
else
   lapack_linker_line=
   lapack_message=Building program with internally-supplied LAPACK
endif

ifdef PGSTK
   tk_message=Building mlslib with toolkit
   tk_linker_line=-L${PGSTK} -lPGSTK
else
   tk_message=Building mlslib without toolkit; using SDPToolkitSubstitute
   tk_linker_line=
endif

fftw_linker_line = $(shell ${UTILDIR}/which_fftw.sh ${FFTW_ROOT} ${FFTW_PREC})

mlspack_dir=$(CONFDIR)/blas/$(MLSCONFG)

mlspack_linker_line=-L${mlspack_dir} -lmlspack
mlspack_message=Building program with mlspack

utctotai_linker_line=-L${INSTALLDIR} -lutctotai
utctotai_message=Building program with toolkitless utc to tai conversion

# $Log$
# Revision 1.6  2010/01/29 18:14:11  pwagner
# Added -dl to Lahey link for xml compatibility
#
# Revision 1.5  2008/07/11 23:54:42  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 1.4  2007/04/16 23:04:20  pwagner
# Changes in line with Intel v9.1
#
# Revision 1.3  2005/08/17 17:45:05  pwagner
# Moved more definitions to srclib/targets.h
#
# Revision 1.2  2005/03/07 17:38:17  pwagner
# Fixed bug in linking non-static LF95
#
# Revision 1.1  2005/03/04 19:00:19  pwagner
# First commit
#
