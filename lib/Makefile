# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make_${OSTYPE}.inc

OBJS=mlsstrings.o signalsfile.o

all: util.a

util.a: ${OBJS}
	ar ru util.a $?

clean:
	-rm -f *.o *.a *.mod *.M
# Revision 1.4  1999/11/04 03:06:49  livesey
# Got the make stuff working I think
#
