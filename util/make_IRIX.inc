# $Id$
#
# This is a makefile include for IRIX
#
F90=f90
F90FLAGS=
F90DEBUGFLAGS= -g -C
#
.SUFFIXES: .o .o_debug .f90
.f90.o:
	${F90} ${F90FLAGS} ${MODPATHFLAG} -c $<

.f90.o_debug:
	${F90} ${F90DEBUGFLAGS} ${MODPATHFLAG} -c -o $*.o_debug $<
