# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make.inc
MLSPGS_LIB_PATH=.

OBJS=	MLSStrings.o MLSSignalNomenclature.o VerticalCoordinate.o \
	QuantityTemplates.o MLSNumerics.o MLSMessageModule.o \
	SDPToolkit.o L2GPData.o L2AUXData.o VectorsModule.o Units.o \
	VGrid.o HGrid.o L1BData.o ConstructQuantityTemplates.o\
	Temporary_types.o ConstructVectorTemplates.o gridded_data_module.o \
	dates_module.o read_parse_l2cf.o MLSPCF.o

OBJS_DEBUG= 	MLSStrings.o_debug MLSSignalNomenclature.o_debug \
		VerticalCoordinate.o_debug QuantityTemplates.o_debug \
		MLSNumerics.o_debug MLSMessageModule.o_debug \
		SDPToolkit.o_debug L2GPData.o_debug L2AUXData.o_debug \
		VectorsModule.o_debug Units.o_debug VGrid.o_debug \
		HGrid.o_debug L1BData.o ConstructQuantityTemplates.o_debug \
		Temporary_types.o_debug ConstructVectorTemplates.o_debug \
		gridded_data_module.o_debug dates_module.o_debug \
		read_parse_l2cf.o MLSPCF.o_debug

LIBPATH=.
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
SDPToolkitSubstitute.a_debug: SDPToolkitSubstitute.o_debug
	ar ru $@ $?

clean:
	-rm -f *.o *.a *.mod *.M *.o_debug *.a_debug
#
# Here we list manually dependencies present through USE statements.
# Later on we may try to adapt Hugh's script to deal with these
#

read_parse_l2cf.o: SDPToolkit.o Hdf.o MLSMessageModule.o
gridded_data_module.o: MLSCommon.o MLSStrings.o MLSMessageModule.o \
	dates_module.o
MLSMessageModule.o: SDPToolkit.o
MLSStrings.o: MLSMessageModule.o
MLSSignalNomenclature.o: MLSStrings.o MLSMessageModule.o MLSCommon.o
QuantityTemplates.o: MLSSignalNomenclature.o MLSStrings.o MLSMessageModule.o \
	VGrid.o Temporary_types.o read_parse_l2cf.o
MLSNumerics.o: MLSMessageModule.o 
L1BData.o: Hdf.o SDPToolkit.o MLSCommon.o
HGrid.o: MLSCommon.o MLSMessageModule.o MLSStrings.o L1BData.o MLSNumerics.o read_parse_l2cf.o
VGrid.o: MLSCommon.o VerticalCoordinate.o MLSStrings.o read_parse_l2cf.o
ConstructQuantityTemplates.o: MLSCommon.o MLSMessageModule.o \
	QuantityTemplates.o VGrid.o HGrid.o VerticalCoordinate.o \
	Units.o MLSSignalNomenclature.o read_parse_l2cf.o
VerticalCoordinate.o: MLSCommon.o MLSStrings.o Units.o
ConstructVectorTemplates.o: MLSCommon.o MLSMessageModule.o VectorsModule.o \
	QuantityTemplates.o read_parse_l2cf.o

#
# Now the debug dependencies
#

gridded_data_module.o_debug: MLSCommon.o_debug MLSStrings.o_debug MLSMessageModule.o_debug \
	dates_module.o_debug
MLSMessageModule.o_debug: SDPToolkit.o_debug
MLSStrings.o_debug: MLSMessageModule.o_debug
MLSSignalNomenclature.o_debug: MLSStrings.o_debug MLSMessageModule.o_debug MLSCommon.o_debug
QuantityTemplates.o_debug: MLSSignalNomenclature.o_debug MLSStrings.o_debug MLSMessageModule.o_debug \
	VGrid.o_debug Temporary_types.o_debug
MLSNumerics.o_debug: MLSMessageModule.o_debug
L1BData.o_debug: Hdf.o_debug SDPToolkit.o_debug MLSCommon.o_debug
HGrid.o_debug: MLSCommon.o_debug MLSMessageModule.o_debug MLSStrings.o_debug L1BData.o_debug MLSNumerics.o_debug
VGrid.o_debug: MLSCommon.o_debug VerticalCoordinate.o_debug MLSStrings.o_debug Temporary_types.o_debug
ConstructQuantityTemplates.o_debug: MLSCommon.o_debug MLSMessageModule.o_debug \
	QuantityTemplates.o_debug VGrid.o_debug HGrid.o_debug VerticalCoordinate.o_debug \
	Units.o_debug MLSSignalNomenclature.o_debug
VerticalCoordinate.o_debug: MLSCommon.o_debug MLSStrings.o_debug Units.o_debug
ConstructVectorTemplates.o_debug: MLSCommon.o_debug MLSMessageModule.o_debug VectorsModule.o_debug \
	QuantityTemplates.o_debug Temporary_types.o_debug

# $Log$
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
