# $Id$
#
# This is a Makefile for my test programs
#

include ../../lib/make.inc

MODULEPATH=../../lib/
LIB=../../lib/mlslib.a
LIB_DEBUG=../../lib/mlslib.a_debug

TKSLIB=../../lib/SDPToolkitSubstitute.a

all:	testnomenclature testhunt testinterpolate
debug:	testnomenclature_debug testhunt_debug testinterpolate_debug
both:	all debug

../../lib/SDPToolkitSubstitute.o:
	cd ../../lib ; ${MAKE} SDPToolkitSubstitute.o
../../lib/SDPToolkitSubstitute.o_debug:
	cd ../../lib ; ${MAKE} SDPToolkitSubstitute.o_debug

lib:	../../lib/mlslib.a
	cd ../../lib ; ${MAKE} mlslib.a
lib_debug: ../../lib/mlslib.a_debug
	cd ../../lib ; ${MAKE} mlslib.a_debug


testnomenclature: lib testnomenclature.o 
	${F90} -o $@ testnomenclature.o ${LIB} ${TKSLIB}

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean
