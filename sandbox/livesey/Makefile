# $Id$
#
# This is a Makefile for my test programs
#

include ../../lib/make.inc

MODULEPATH=../../lib/
LIB=../../lib/mlslib.a
LIB_DEBUG=../../lib/mlslib.a_debug

all:	testnomenclature testhunt testinterpolate
debug:	testnomenclature_debug testhunt_debug testinterpolate_debug
both:	all debug

lib:	../../lib/mlslib.a
	cd ../../lib ; ${MAKE} mlslib.a
lib_debug: ../../lib/mlslib.a_debug
	cd ../../lib ; ${MAKE} mlslib.a_debug


# testnomenclature: lib testnomenclature.o
# 	${F90} -o $@ testnomenclature.o ${LIB} 
# testnomenclature_debug: lib_debug testnomenclature.o_debug
# 	${F90} -o $@ testnomenclature.o_debug ${LIB_DEBUG}

# testhunt: lib testhunt.o
# 	${F90} -o $@ testhunt.o ${LIB} 
# testhunt_debug: lib_debug testhunt.o_debug
# 	${F90} -o $@ testhunt.o_debug ${LIB_DEBUG}

testinterpolate: lib testinterpolate.o
	${F90} -o $@ testinterpolate.o ${LIB} 
testinterpolate_debug: lib_debug testinterpolate.o_debug
	${F90} -o $@ testinterpolate.o_debug ${LIB_DEBUG}

clean:
	-rm *.o *.mod *.M

allclean: clean
	cd ../../lib ; ${MAKE} clean
