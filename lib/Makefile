# $Id$
#
# This is a makefile for the utils section of the mlspgs suite
#
include make.inc

OBJS= mlsstrings.o signalsfile.o
OBJS_DEBUG= mlsstrings.o_debug signalsfile.o_debug
MODULEPATH=.

all: util.a
debug: util.a_debug

util.a: ${OBJS}
	ar ru util.a $?

util.a_debug: ${OBJS_DEBUG}
	ar ru util.a_debug $?

clean:
	-rm -f *.o *.a *.mod *.M *.o_debug *.a_debug
#
# $Log$
# Revision 1.4  1999/11/04 03:06:49  livesey
# Got the make stuff working I think
#
