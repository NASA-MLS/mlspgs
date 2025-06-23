# $Id$
#
# This is a Makefile for my test programs
#
MODULEPATH=../../lib/
MLSPGS_LIB_PATH=../../lib
TOOLKIT_LIB=../../lib/SDPToolkitSubstitute.a

include ${MLSPGS_LIB_PATH}/make.inc

all:	testnomenclature testunits testvert
debug:	testnomenclature_debug testunits_debug testvert_debug
both:	all debug

testnomenclature: testnomenclature.o mlslib #SDPToolkitSubstitue
	${F90} -o $@ testnomenclature.o \
	${MLSPGS_LIB_PATH}/mlslib.a \
	${TOOLKIT_LIB}

testnomenclature.o: mlslib

# testunits: testunits.o mlslib
# 	${F90} -o $@ testunits.o \
# 	${MLSPGS_LIB_PATH}/mlslib.a \
# 	${TOOLKIT_LIB}

# testunits.o: mlslib

# testvert: testvert.o mlslib
# 	${F90} -o $@ testvert.o \
# 	${MLSPGS_LIB_PATH}/mlslib.a \
# 	${TOOLKIT_LIB}

# testvert.o: mlslib

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean

include ${MLSPGS_LIB_PATH}/makelibs.inc




