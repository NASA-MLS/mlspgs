# $Id$
#
# This is a Makefile for my test programs

F90FLAGS= -I../../util
LNKFLAGS= -I../../util
LIBS=../../util/util.a

all:	test

util:
	cd ../../util ; ${MAKE} all

test:	util test.o
	${F90} ${LNKFLAGS} -o test test.o ${LIBS}

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../util ; ${MAKE} clean
