# $Id$
#
# This is a Makefile for my test programs
#

include ../../lib/make.inc

MODULEPATH=../../lib/
LIB=../../lib/mlslib.a
LIB_DEBUG=../../lib/mlslib.a_debug

OBJS=test.o
OBJS_DEBUG=test.o_debug

all:	test
debug:	test_debug
both:	all debug

lib:
	cd ../../lib ; ${MAKE} mlslib.a
lib_debug:
	cd ../../lib ; ${MAKE} mlslib.a_debug


test:	lib ${OBJS}
	${F90} -o test ${OBJS} ${LIB}

test_debug:	lib_debug ${OBJS_DEBUG}
	${F90} -g -o test_debug ${OBJS_DEBUG} ${LIB_DEBUG}

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean



