# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make.inc

OBJS=	mlsstrings.o mlssignalnomenclature.o verticalcoordinates.o \
	vectorquantities.o mlsnumerics.o mlsmessagemodule.o
OBJS_DEBUG= mlsstrings.o_debug mlssignalnomenclature.o_debug \
	verticalcoordinates.o_debug vectorquantities.o_debug \
	mlsnumerics.o_debug mlsmessagemodule.o_debug
MODULEPATH=.

both: all debug
all: mlslib.a SDPToolkitSubstitute.a
debug: mlslib.a_debug SDPToolkitSubstitute.a_debug

mlslib.a: ${OBJS}
	ar ru mlslib.a $?

mlslib.a_debug: ${OBJS_DEBUG}
	ar ru mlslib.a_debug $?

SDPToolkitSubstitute.a: ${TKS_OBJS}
	ar ru SDPToolkitSubstitute.a $?
SDPToolkitSubstitute.a_debug: ${TKS_OBJS}
	ar ru SDPToolkitSubstitute.a $?

clean:
	-rm -f *.o *.a *.mod *.M *.o_debug *.a_debug

mlsstrings.o: mlsmessagemodule.o
mlssignalnomenclature.o: mlsstrings.o mlsmessagemodule.o
vectorquantities.o: mlssignalnomenclature.o mlsstrings.o mlsmessagemodule.o
mlsnumerics.o: mlsmessagemodule.o

mlsstrings.o_debug: mlsmessagemodule.o_debug
mlssignalnomenclature.o_debug: mlsstrings.o_debug mlsmessagemodule.o
vectorquantities.o_debug: mlssignalnomenclature.o_debug mlsstrings.o_debug \
	mlsmessagemodule.o_debug
mlsnumerics.o_debug: mlsmessagemodule.o_debug

#
# $Log$
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
