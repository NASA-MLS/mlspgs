# $Id$
#
# This is the makefile for the EOS MLS Level 1 software
#
LDFLAGS= -lm -lnsl
MODULEPATH=../lib/
MLSPGS_LIB_PATH=../lib
TOOLKIT_LIB=${PGSLIB}/libPGSTK.a $
HDF_LIB=/usr/local/HDF/lib/libmfhdf.a /usr/local/HDF/lib/libdf.a \
	/usr/local/HDF/lib/libz.a

# Note we define OBJS to be all the files EXCEPT MLSL1.o as this allows
# the dependencies to be defined more simply.  Where we need to, we add
# MLSL1.o in.

OBJS= MLSL1Common.o MLSL1Rad.o OutputL1B.o open_init.o sort_qualify.o calib.o \
	close_files.o Scan.o Orbit.o TkL1B.o
OBJS_DEBUG= MLSL1Common.o_debug MLSL1Rad.o_debug OutputL1B.o_debug \
	open_init.o_debug sort_qualify.o_debug calib.o_debug \
	close_files.o_debug Scan.o_debug Orbit.o TkL1B.o_debug

include ${MLSPGS_LIB_PATH}/make.inc

all: MLSL1
debug: MLSL1

MLSL1: $(OBJS) MLSL1.o mlslib
	${F90} ${F90FLAGS} ${LDFLAGS} PGS_TD_ADEOSIItoTAI.obj -o $@ ${OBJS} \
	MLSL1.o ${MLSPGS_LIB_PATH}/mlslib.a ${TOOLKIT_LIB} ${HDF_LIB}
MLSL1_debug: $(OBJS_DEBUG) MLSL1.o_debug mlslib_debug
	${F90} ${F90DEBUGFLAGS} ${LDFLAGS} -o $@ ${OBJS} MLSL1.o_debug  \
	${MLSPGS_LIB_PATH}/mlslib.a_debug ${TOOLKIT_LIB} ${HDF_LIB}

MLSL1.o: $(OBJS)

clean:
	-rm -f *.o *.mod *.M

# Now give the additional dependencies through the use statement

include ${MLSPGS_LIB_PATH}/makelibs.inc

#
# $Log$
# Revision 1.3  2000/02/15 19:05:41  nakamura
# Removed Mc module.
#
# Revision 1.2  2000/02/14 17:27:56  perun
# Added modules for L1BOA records and L1B data
#
# Revision 1.1  2000/02/08 21:45:16  perun
# Initial version
#
#
