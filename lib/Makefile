# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make.inc

OBJS=	MLSStrings.o MLSSignalNomenclature.o VerticalCoordinate.o \
	VectorQuantities.o MLSNumerics.o MLSMessageModule.o

OBJS_DEBUG= 	MLSStrings.o_debug MLSSignalNomenclature.o_debug \
		VerticalCoordinate.o_debug VectorQuantities.o_debug \
		MLSNumerics.o_debug MLSMessageModule.o_debug

TKS_OBJS= SDPToolkitSubstitute.o
TKS_OBJS_DEBUG= SDPToolkitSubstitute.o_debug

MODULEPATH=.

both: all debug
all: mlslib.a SDPToolkitSubstitute.a
debug: mlslib.a_debug SDPToolkitSubstitute.a_debug

mlslib.a: ${OBJS}
	ar ru $@ $?

mlslib.a_debug: ${OBJS_DEBUG}
	ar ru $@ $?

SDPToolkitSubstitute.a: ${TKS_OBJS}
	ar ru SDPToolkitSubstitute.a $?
SDPToolkitSubstitute.a_debug: ${TKS_OBJS}
	ar ru SDPToolkitSubstitute.a $?

clean:
	-rm -f *.o *.a *.mod *.M *.o_debug *.a_debug

MLSStrings.o: MLSMessageModule.o
MLSSignalNomenclature.o: MLSStrings.o MLSMessageModule.o
VectorQuantities.o: MLSSignalNomenclature.o MLSStrings.o MLSMessageModule.o
MLSNumerics.o: MLSMessageModule.o

MLSStrings.o_debug: MLSMessageModule.o_debug
MLSSignalNomenclature.o_debug: MLSStrings.o_debug MLSMessageModule.o_debug
VectorQuantities.o_debug: MLSSignalNomenclature.o_debug MLSStrings.o_debug \
		MLSMessageModule.o_debug
MLSNumerics.o_debug: MLSMessageModule.o_debug


#
# $Log$
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
