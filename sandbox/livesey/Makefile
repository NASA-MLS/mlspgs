# $Id$
#
# This is a Makefile for my test programs
#
MODULEPATH=../../lib/
MLSPGS_LIB_PATH=../../lib
TOOLKIT_LIB=../../lib/SDPToolkitSubstitute.a

include ${MLSPGS_LIB_PATH}/make.inc

all:	testnomenclature testunits
debug:	testnomenclature_debug testunits_debug
both:	all debug

# testnomenclature: testnomenclature.o mlslib #SDPToolkitSubstitue
# 	${F90} -o $@ testnomenclature.o \
# 	${MLSPGS_LIB_PATH}/mlslib.a \
# 	${TOOLKIT_LIB}

# testnomenclature.o: mlslib

testunits: testunits.o mlslib
	${F90} -o $@ testunits.o \
	${MLSPGS_LIB_PATH}/mlslib.a \
	${TOOLKIT_LIB}

testunits.o: mlslib

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean

include ${MLSPGS_LIB_PATH}/makelibs.inc




