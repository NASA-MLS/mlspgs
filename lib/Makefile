# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make.inc

OBJS=	mlsstrings.o mlssignalnomenclature.o verticalcoordinates.o \
	vectorquantities.o
OBJS_DEBUG= mlsstrings.o_debug mlssignalnomenclature.o_debug \
	verticalcoordiantes.o_debug vectorquantities.o_debug
MODULEPATH=.

all: mlslib.a
debug: mlslib.a_debug

mlslib.a: ${OBJS}
	ar ru mlslib.a $?

mlslib.a_debug: ${OBJS_DEBUG}
	ar ru mlslib.a_debug $?

clean:
	-rm -f *.o *.a *.mod *.M *.o_debug *.a_debug

mlssignalnomenclature.o: mlsstrings.o
vectorquantities.o: mlssignalnomenclature.o mlsstrings.o

mlssignalnomenclature.o_debug: mlsstrings.o_debug
vectorquantities.o_debug: mlssignalnomenclature.o_debug mlsstrings.o_debug

#
# $Log$
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
