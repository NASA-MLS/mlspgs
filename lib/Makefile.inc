# $Id$

# This is a makefile for the utils section of the mlspgs suite

# Define some suffix rules

.SUFFIXES: .o .dbg.o .f90
.f90.o:
	${F90} ${F90FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -c \
               $(INCLUDE) $<

.f90.dbg.o:
	${F90} ${F90DEBUGFLAGS} ${MODULEPATHPREFIX}${MODULEPATH} \
		-c -o $*.dbg.o $(INCLUDE) $<

MLSPGS_LIB_PATH=.

#OBJS=	MLSStrings.o MLSSignalNomenclature.o VerticalCoordinate.o \
#	QuantityTemplates.o MLSNumerics.o MLSMessageModule.o \
#	SDPToolkit.o L2GPData.o L2AUXData.o VectorsModule.o Units.o \
#	VGrid.o HGrid.o L1BData.o ConstructQuantityTemplates.o\
#	ConstructVectorTemplates.o gridded_data_module.o \
#	dates_module.o MLSPCF.o MLSCF.o GriddedData.o HDF-EOS.o

OBJS=	dump_0.o GriddedData.o HDF-EOS.o L1BData.o L2AUXData.o L2GPData.o  \
        MLSCommon.o MLSMessageModule.o MLSNumerics.o MLSPCF.o \
        MLSSignalNomenclature.o MLSStrings.o SDPToolkit.o \
        swapi.o swapi_character_array.o swapi_character_scalar.o \
        swapi_double.o swapi_integer.o swapi_real.o Units.o \
        VerticalCoordinate.o dates_module.o gridded_data_module.o

OBJS_DEBUG= 	MLSStrings.dbg.o MLSSignalNomenclature.dbg.o \
		VerticalCoordinate.dbg.o QuantityTemplates.dbg.o \
		MLSNumerics.dbg.o MLSMessageModule.dbg.o \
		SDPToolkit.dbg.o L2GPData.dbg.o L2AUXData.dbg.o \
		VectorsModule.dbg.o Units.dbg.o VGrid.dbg.o \
		HGrid.dbg.o L1BData.o ConstructQuantityTemplates.dbg.o \
		ConstructVectorTemplates.dbg.o \
		gridded_data_module.dbg.o dates_module.dbg.o \
		MLSPCF.dbg.o MLSCF.dbg.o GriddedData.dbg.o

Libpath=.
MODULEPATH=.

all:	mlslib.a SDPToolkitSubstitute.a
debug:	mlslib.a_debug SDPToolkitSubstitute.a_debug
both: all debug

mlslib.a: ${OBJS}
	ar ru $@ $?
mlslib.a_debug: ${OBJS_DEBUG}
	ar ru $@ $?

SDPToolkitSubstitute.a: SDPToolkitSubstitute.o
	ar ru $@ $?
SDPToolkitSubstitute.a_debug: SDPToolkitSubstitute.dbg.o
	ar ru $@ $?

clean:
	-rm -f *.o *.a *.mod *.M *.dbg.o *.a_debug
#
# Here we list manually dependencies present through USE statements.
# Later on we may try to adapt Hugh's script to deal with these
#

dump_0.o: $S/dump_0.f90 $C/output_m.o
	${F90} ${F90FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -c \
               $(INCLUDE) ${MODULEPATHPREFIX}$C $S/dump_0.f90
ConstructQuantityTemplates.o: MLSCommon.o MLSMessageModule.o \
	QuantityTemplates.o VGrid.o HGrid.o VerticalCoordinate.o \
	Units.o MLSSignalNomenclature.o MLSCF.o
ConstructVectorTemplates.o: MLSCommon.o MLSMessageModule.o VectorsModule.o \
	QuantityTemplates.o MLSCF.o
gridded_data_module.o: MLSCommon.o MLSStrings.o MLSMessageModule.o \
	dates_module.o
HGrid.o: MLSCommon.o MLSMessageModule.o MLSStrings.o L1BData.o MLSNumerics.o MLSCF.o
L1BData.o: Hdf.o SDPToolkit.o MLSCommon.o
	${F90} ${F90FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -c \
               $(INCLUDE) $(DUSTY) ../L1BData.f90
MLSCF.o: SDPToolkit.o Hdf.o MLSMessageModule.o
MLSMessageModule.o: SDPToolkit.o
MLSNumerics.o: MLSMessageModule.o 
MLSSignalNomenclature.o: MLSStrings.o MLSMessageModule.o MLSCommon.o
MLSStrings.o: MLSMessageModule.o
QuantityTemplates.o: MLSSignalNomenclature.o MLSStrings.o MLSMessageModule.o \
	VGrid.o MLSCF.o
swapi.o: swapi_character_array.o swapi_character_scalar.o swapi_double.o \
        swapi_integer.o swapi_real.o
VerticalCoordinate.o: MLSCommon.o MLSStrings.o Units.o
VGrid.o: MLSCommon.o VerticalCoordinate.o MLSStrings.o MLSCF.o

#
# Now the debug dependencies
#

ConstructQuantityTemplates.dbg.o: MLSCommon.dbg.o MLSMessageModule.dbg.o \
	QuantityTemplates.dbg.o VGrid.dbg.o HGrid.dbg.o VerticalCoordinate.dbg.o \
	Units.dbg.o MLSSignalNomenclature.dbg.o
ConstructVectorTemplates.dbg.o: MLSCommon.dbg.o MLSMessageModule.dbg.o VectorsModule.dbg.o \
	QuantityTemplates.dbg.o
gridded_data_module.dbg.o: MLSCommon.dbg.o MLSStrings.dbg.o MLSMessageModule.dbg.o \
	dates_module.dbg.o
HGrid.dbg.o: MLSCommon.dbg.o MLSMessageModule.dbg.o MLSStrings.dbg.o L1BData.dbg.o MLSNumerics.dbg.o
L1BData.dbg.o: Hdf.dbg.o SDPToolkit.dbg.o MLSCommon.dbg.o
MLSCF.dbg.o: SDPToolkit.dbg.o Hdf.dbg.o MLSMessageModule.dbg.o
MLSMessageModule.dbg.o: SDPToolkit.dbg.o
MLSNumerics.dbg.o: MLSMessageModule.dbg.o
MLSSignalNomenclature.dbg.o: MLSStrings.dbg.o MLSMessageModule.dbg.o MLSCommon.dbg.o
MLSStrings.dbg.o: MLSMessageModule.dbg.o
QuantityTemplates.dbg.o: MLSSignalNomenclature.dbg.o MLSStrings.dbg.o MLSMessageModule.dbg.o \
	VGrid.dbg.o 
VerticalCoordinate.dbg.o: MLSCommon.dbg.o MLSStrings.dbg.o Units.dbg.o
VGrid.dbg.o: MLSCommon.dbg.o VerticalCoordinate.dbg.o MLSStrings.dbg.o \
	 MLSCF.dbg.o

# $Log$
# Revision 2.1  2000/09/13 20:38:23  vsnyder
# Support to compile dump_0.
#
# Revision 2.0  2000/09/05 17:41:06  dcuddy
# Change revision to 2.0
#
# Revision 1.1  2000/09/02 01:51:18  vsnyder
# Initial entry
#
# Revision 1.26  2000/02/23 20:33:48  perun
# Changed .o_debug to .dbg.o for SunOS
#
# Revision 1.25  2000/02/05 00:30:14  livesey
# Removed temporary types stuff
#
# Revision 1.24  2000/01/21 21:00:52  livesey
# No longer necessary
#
# Revision 1.23  2000/01/11 22:51:34  livesey
# Dealt with ramifications of change from read_parse_l2cf to MLSCF
#
# Revision 1.22  2000/01/07 23:53:35  livesey
# Nearly integrated, just a few tweaks.
#
# Revision 1.21  1999/12/20 14:50:23  pumphrey
# Added gridded_data_module to Makefile
#
# Revision 1.20  1999/12/18 01:07:13  livesey
# Regular commit
#
# Revision 1.19  1999/12/17 00:59:32  livesey
# Nightly checkin
#
# Revision 1.18  1999/12/16 23:45:15  livesey
# Regular commit
#
# Revision 1.18  1999/12/16 23:45:15  livesey
# Regular commit
#
# Revision 1.17  1999/12/16 18:23:31  livesey
# Routine checkin
#
# Revision 1.16  1999/12/16 01:27:45  livesey
# Routine checkin
#
# Revision 1.15  1999/12/14 22:56:02  livesey
# Regular commit
#
# Revision 1.14  1999/12/04 00:27:43  livesey
# Regular end of day checkin
#
# Revision 1.13  1999/12/03 00:14:18  livesey
# Nightly checkin
#
# Revision 1.12  1999/12/02 01:56:11  livesey
# Changed filenames to mixed case, redid makefiles.
#
# Revision 1.11  1999/12/01 23:14:06  livesey
# Renamed files to upper/lower case
#
# Revision 1.10  1999/12/01 23:01:39  livesey
# Before renaming things to upper/lower case
#
# Revision 1.9  1999/12/01 05:03:55  livesey
# Nightly checkin
#
# Revision 1.8  1999/11/30 04:04:17  livesey
# Nightly checkin
#
# Revision 1.7  1999/11/24 23:07:20  livesey
# Regular commit
#
# Revision 1.6  1999/11/10 01:52:00  livesey
# Transfer to home
#
# Revision 1.5  1999/11/04 21:59:40  livesey
# Move contents of util to lib and renamed to mlslib.a
#
# Revision 1.4  1999/11/04 03:06:49  livesey
# Got the make stuff working I think
#
