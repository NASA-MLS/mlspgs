# $Id$
# $RCSfile$

# Notice that SDPToolkitSubstitute is not in SRCS or OBJS.
# The intent is that you will put SDPToolkitSubstitute.o into
# your build explicitly, so that when libmls.a is searched, the
# entry points in SDPToolkit will already be satisfied, so
# SDPToolkit will not be included from libmls.a

SRCS =	Allocate_Deallocate.f90 GriddedData.f90 HDF-EOS.f90 Hdf.f90 \
        L1BData.f90 MLSCF.f90 MLSCommon.f90 MLSMessageModule.f90 \
        MLSNumerics.f90 MLSPCF.f90 MLSSignalNomenclature.f90 MLSStrings.f90 \
        MatrixModule_0.f90 SDPToolkit.f90 Units.f90 \
        VerticalCoordinate.f90 dates_module.f90 dump_0.f90 \
        gridded_data_module.f90 heights_module.f90 machine.f90 swapi.f90 \
        swapi_character_array.f90 swapi_character_scalar.f90 \
        swapi_double.f90 swapi_integer.f90 swapi_real.f90


OBJS =	Allocate_Deallocate.o GriddedData.o HDF-EOS.o Hdf.o \
        L1BData.o MLSCF.o MLSCommon.o MLSMessageModule.o \
        MLSNumerics.o MLSPCF.o MLSSignalNomenclature.o MLSStrings.o \
        MatrixModule_0.o SDPToolkit.o Units.o \
        VerticalCoordinate.o dates_module.o dump_0.o \
        gridded_data_module.o heights_module.o machine.o swapi.o \
        swapi_character_array.o swapi_character_scalar.o \
        swapi_double.o swapi_integer.o swapi_real.o

all: libmls.a SDPToolkitSubstitute.a

libmls.a: $(OBJS)
	ar ru $@ $?

SDPToolkitSubstitute.a: SDPToolkitSubstitute.o
	ar ru $@ $?

clean:
	rm -f *.o *.mod *.a

.SUFFIXES: .o .f90

#.f90.o:
%.o: %.f90
	$(F90) $(F90FLAGS) -c $I$S $I$C $I$T $<

Allocate_Deallocate.o: MLSMessageModule.o
GriddedData.o: MLSCommon.o MLSMessageModule.o MLSStrings.o
Hdf.o: MLSCommon.o
L1BData.o: MLSCommon.o MLSMessageModule.o
	$(F90) $(F90FLAGS) -c $I$S $I$C $I$T $(DUSTY) $S/L1BData.f90
MatrixModule_0.o: Allocate_Deallocate.o dump_0.o \
          MLSCommon.o MLSMessageModule.o $C/output_m.o
MLSCF.o: Hdf.o MLSCommon.o MLSMessageModule.o MLSStrings.o SDPToolkit.o
MLSMessageModule.o: SDPToolkit.o
MLSNumerics.o: MLSCommon.o MLSMessageModule.o MLSStrings.o
MLSSignalNomenclature.o: MLSCommon.o MLSMessageModule.o MLSStrings.o
MLSStrings.o: MLSMessageModule.o
SDPToolkit.o: MLSCommon.o $T/PGS_CBP.f $T/PGS_CBP_6.f $T/PGS_CSC_4.f \
        $T/PGS_EPH_5.f $T/PGS_IO.f $T/PGS_IO_1.f $T/PGS_MEM_7.f \
	$T/PGS_MET.f $T/PGS_MET_13.f $T/PGS_PC.f $T/PGS_PC_9.f \
	$T/PGS_SMF.f $T/PGS_TD.f $T/PGS_TD_3.f
Units.o: MLSCommon.o MLSMessageModule.o MLSStrings.o
VerticalCoordinate.o: MLSCommon.o MLSStrings.o Units.o
dates_module.o: MLSMessageModule.o MLSStrings.o
dump_0.o: 
gridded_data_module.o: MLSCommon.o MLSMessageModule.o MLSStrings.o \
	dates_module.o
heights_module.o: MLSCommon.o MLSMessageModule.o
swapi.o: swapi_character_array.o swapi_character_scalar.o swapi_double.o \
	swapi_integer.o swapi_real.o

# ========================================================================
# Test harnesses

test_multiply.o: $S/test_multiply.f90 MatrixModule_0.o
	$(F90) -c $(FOPTS) $I$C $S/test_multiply.f90

test_multiply: test_multiply.o
	$(F90) -o test_multiply test_multiply.o -L. -lmls -L$C -lCF \
                 -L$(HOME)/mlspgs/SDPToolkit -lPGSTK \

test_sparsify.o: $S/test_sparsify.f90 MatrixModule_0.o
	$(F90) -c $(FOPTS) $I$C $S/test_sparsify.f90

test_sparsify: test_sparsify.o
	$(F90) -o test_sparsify test_sparsify.o -L. -lmls -L$C -lCF \
                 -L$(HOME)/mlspgs/SDPToolkit -lPGSTK \

# ========================================================================

# $Log$
# Revision 2.6  2000/10/05 01:35:34  vsnyder
# Revised to account for moving test_multiply from mlspgs/l2.  Also,
# mlslib.a was deleted, because libraries are required to begin with
# lib, so we have libmls.a instead.
#
# Revision 2.5  2000/10/04 20:24:21  vsnyder
# Account for moving Allocate_Deallocate, MatrixModule_0 and test_sparsify
# from mlspgs/l2
#
# Revision 2.4  2000/10/04 18:03:06  vsnyder
# Added a rule to make libmls.a from mlslib.a
#
# Revision 2.3  2000/10/04 01:25:11  vsnyder
# Added an explicit action for L1BData, to specify the DUSTY option,
# because sfrdata is called with inconsistent argument lists.
#
# Revision 2.2  2000/09/29 18:05:32  vsnyder
# Add CVS comments.  Replace inference rule ".f90.o" with "%.o: %.f90"
# because the former sometimes didn't "make" out-of-date files, but the
# latter seems always to do so.
#
