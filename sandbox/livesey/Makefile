# $Id$
#
# This is a Makefile for my test programs
#

include ../../util/make.inc

MODULEPATH=../../util/
LIB=../../util/util.a
LIB_DEBUG=../../util/util.a_debug

OBJS=test.o
OBJS_DEBUG=test.o_debug

all:	test

lib:
	cd ../../util ; ${MAKE} util.a
lib_debug:
	cd ../../util ; ${MAKE} util.a_debug


test:	lib ${OBJS}
	${F90} -o test ${OBJS} ${LIB}

test_debug:	lib_debug ${OBJS_DEBUG}
	${F90} -o test_debug ${OBJS_DEBUG} ${LIB_DEBUG}

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../util ; ${MAKE} clean



