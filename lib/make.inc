# $Id$
#
# This is an include file for makefiles that is intended to be as
# platform independent as possible.  Where things are platform dependent,
# an if statement is used.
#
F90=f95
#
# First define some suffix rules
#
.SUFFIXES: .o .dbg.o .f90
.f90.o:
	${F90} ${F90FLAGS}${INCLUDEPATHPREFIX}  ${MODULEPATHPREFIX}${MODULEPATH} -c -C -g  -mismatch $<

.f90.dbg.o:
	${F90} ${F90DEBUGFLAGS} ${MODULEPATHPREFIX}${MODULEPATH} ${INCLUDEPATHPREFIX} \
		-c -g $*.dbg.o $<
#
# Now include the nasty platform dependent stuff
#
ostype = $(OSTYPE)
#
# Check for Solaris OS which is the same as SunOS
#
ifeq ($(ostype),solaris)
ostype = SunOS
endif

# First the linux dependent stuff, this assumes absoft Fortran for the moment
#
ifeq ($(ostype),linux)
MODULEPATHPREFIX=-p
F90DEBUGFLAGS=-g
INCLUDEPATHPREFIX=-I/software/SRC/TOOLKIT-5.2.6v1.00.NAG/include
endif
#
ifeq ($(ostype),IRIX)
MODULEPATHPREFIX=-I
F90DEBUGFLAGS=-g -C
INCLUDEPATHPREFIX=-I/software/SRC/TOOLKIT-5.2.6v1.00.NAG/include
endif
#
ifeq ($(ostype),IRIX64)
MODULEPATHPREFIX=-I
F90DEBUGFLAGS=-g -C
INCLUDEPATHPREFIX=-I/software/SRC/TOOLKIT-5.2.6v1.00.NAG/include
endif
#
ifeq ($(ostype),SunOS)
MODULEPATHPREFIX=-I
INCLUDEPATHPREFIX=-I/software/SRC/TOOLKIT-5.2.6v1.00.NAG/include
F90FLAGS=-g 
F90DEBUGFLAGS=-g
endif

#
# $Log$
# Revision 1.13  2000/06/30 22:55:24  lungu
# updated for NAG.
#
# Revision 1.12  2000/06/20 00:01:20  lungu
# Added include files path /usr/local/TOOLKIT/include.
#
# Revision 1.11  2000/05/18 00:12:06  lungu
# added INCLUDEPATHPREFIX=-I/usr/local/TOOLKIT/include
#
# Revision 1.10  2000/05/17 19:38:52  livesey
# Added IRIX64 option for abacus.
#
# Revision 1.9  2000/02/23 20:34:16  perun
# Changed .o_debug to .dbg.o for SunOS
#
# Revision 1.8  2000/02/22 19:12:45  perun
# Added F90FLAGS for SunOS
#
# Revision 1.7  2000/01/18 15:19:04  perun
# Make solaris same as SunOS
#
# Revision 1.6  1999/12/02 01:56:11  livesey
# Changed filenames to mixed case, redid makefiles.
#
# Revision 1.5  1999/11/04 23:31:27  livesey
# Made the SunOS stuff correct!
#
# Revision 1.4  1999/11/04 23:29:56  livesey
# Added the stuff for SunOS
#
# Revision 1.3  1999/11/04 03:06:49  livesey
# Got the make stuff working I think
#
# Revision 1.2  1999/11/04 02:51:15  livesey
# Fixed bug
#
# Revision 1.1  1999/11/04 02:46:22  livesey
# Added make.inc
#
#
