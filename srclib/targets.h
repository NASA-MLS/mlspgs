# "$Id$"
STATIC_LINK := $(shell echo $(LDOPTS) | grep -e "-static")

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
#   blas_linker_line=-L${BLAS} -l${LIB_BLAS}
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
#   lapack_linker_line=-L${LAPACK} -l${LIB_LAPACK}
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

ifdef HDF5
ifeq ($(MLSPLAT),Sun)
   hdf5_linker_line=-L${HDF5} -lhdf5_fortran -lhdf5 -R${HDF5}
else
#   hdf5_linker_line=-L${HDF5} -lhdf5_fortran -lhdf5 -Xlinker -rpath=${HDF5} \
#    -lssl -lz
   hdf5_linker_line=-L${HDF5} -lhdf5_fortran -lhdf5 
   l_specials := $(l_specials) -lz
endif
ifdef HDFEOS5
   hdfeos5_linker_line=-L${HDFEOS5} -lhe5_hdfeos
ifeq ($(MLSF95),LF95)
   gcc_version=$(shell ${UTILDIR}/which_fftw.sh -gcc)
   gcc_libloc=$(shell ${REECHO} -d /usr/lib/gcc*/i386-redhat-linux/${gcc_version})
   static_object=$(shell ${UTILDIR}/which_fftw.sh -static)
ifdef STATIC_LINK
   l_specials=-lz -L${gcc_libloc} -lgcc \
     ${static_object}
else
   l_specials=-lz -L${gcc_libloc} -lgcc
endif
endif
ifeq ($(MLSF95),IFC)
   gcc_version=$(shell ${UTILDIR}/which_fftw.sh -gcc)
   gcc_libloc=$(shell ${REECHO} -d /usr/lib/gcc*/i386-redhat-linux/${gcc_version})
   l_specials=-lz -L${gcc_libloc} -lgcc
endif
else
   hdfeos5_linker_line=
endif
ifdef HDFVERSIONS
   he5lib_linker_line=-L$(CONFDIR)/he5lib/$(MLSCONFG) -lhe5
   he5lib_message=mlsl2 will support files with either hdf version
   mlsl2_prereqs := $(mlsl2_prereqs) $(he5libdir)/libhe5.a
endif
   hdf5_message=Building program with hdf5
else
   hdf5_linker_line=
   hdfeos5_linker_line=
   hdf5_message=Building program without hdf5
endif

fftw_linker_line = $(shell ${UTILDIR}/which_fftw.sh ${FFTW_ROOT} ${FFTW_PREC})
# $Log$
# Revision 1.1  2005/03/04 19:00:19  pwagner
# First commit
#
