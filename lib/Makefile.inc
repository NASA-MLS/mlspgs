# $Id$
# $RCSfile$

SRCS =	GriddedData.f90 HDF-EOS.f90 Hdf.f90 L1BData.f90 MLSCF.f90 \
	MLSCommon.f90 MLSMessageModule.f90 MLSNumerics.f90 MLSPCF.f90 \
	MLSSignalNomenclature.f90 MLSStrings.f90 SDPToolkit.f90 \
	SDPToolkitSubstitute.f90 Units.f90 VerticalCoordinate.f90 \
	dates_module.f90 dump_0.f90 gridded_data_module.f90 \
	heights_module.f90 swapi.f90 swapi_character_array.f90 \
	swapi_character_scalar.f90 swapi_double.f90 swapi_integer.f90 \
	swapi_real.f90

OBJS =	GriddedData.o HDF-EOS.o Hdf.o L1BData.o MLSCF.o MLSCommon.o \
	MLSMessageModule.o MLSNumerics.o MLSPCF.o MLSSignalNomenclature.o \
	MLSStrings.o SDPToolkit.o SDPToolkitSubstitute.o Units.o \
	VerticalCoordinate.o dates_module.o dump_0.o gridded_data_module.o \
	heights_module.o swapi.o swapi_character_array.o \
	swapi_character_scalar.o swapi_double.o swapi_integer.o swapi_real.o

LIBS =	

all: mlslib.a SDPToolkitSubstitute.a

mlslib.a: $(OBJS)
	ar ru $@ $?

SDPToolkitSubstitute.a: SDPToolkitSubstitute.o
	ar ru $@ $?

clean:
	rm -f *.o *.mod *.a

.SUFFIXES: .o .f90

#.f90.o:
%.o: %.f90
	$(F90) $(F90FLAGS) -c $I$S $I$C $I$T $<

GriddedData.o: MLSCommon.o MLSMessageModule.o MLSStrings.o
Hdf.o: MLSCommon.o
L1BData.o: MLSCommon.o MLSMessageModule.o
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

# $Log$
