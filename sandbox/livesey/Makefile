# $Id$
#
# This is a Makefile for my test programs
#
MODULEPATH=../../lib/
MLSPGS_LIB_PATH=../../lib
TOOLKIT_LIB=../../lib/SDPToolkitSubstitute.a

include ${MLSPGS_LIB_PATH}/make.inc

all:	testnomenclature
debug:	testnomenclature_debug
both:	all debug

testnomenclature: testnomenclature.o mlslib #SDPToolkitSubstitue
	${F90} -o $@ testnomenclature.o \
	${MLSPGS_LIB_PATH}/mlslib.a \
	${TOOLKIT_LIB}

testnomenclature.o: mlslib

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean

include ${MLSPGS_LIB_PATH}/makelibs.inc




