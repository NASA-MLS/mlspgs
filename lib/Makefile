# $Id$
#
# This is a Makefile for my test programs

OBJ=mlsstrings.o signalsfile.o

all: util.a

util.a: ${OBJ}
	ar ru util.a $?

clean:
	-rm -f *.o *.a *.mod *.M
# Revision 1.4  1999/11/04 03:06:49  livesey
# Got the make stuff working I think
#
