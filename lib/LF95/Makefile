# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# Places to look for sources without being explicitly told to do so:
VPATH = ..

# End of the path to the directory containing this Makefile:
ME = LF95

# Fortran compiler:
F90 = lf95

# Compiler options:
FOPTS = -g --chk --info

# DUSTY is used to allow picky compilers to handle references to "sfrdata",
# which has an inconsistent interface in that it can have any type
# object for its fifth argument.

DUSTY =

# Loader options:
LDOPTS = -g --chk

# Option to specify path to .mod and include files:
I = -I

# Path to source files
S = ..

# Beginning of path to toolkit libraries:
TOOLKIT = /nas/software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Platform (used within directory paths in ../Makefile.inc):
PLATFORM = linux

# Path to libPGSTK.a
PGSTK = $(HOME)/mlspgs/SDPToolkit/$(ME)
#PGSTK = $(TOOLKIT)/lib/$(PLATFORM)

include $S/Makefile.inc

machine.o: machine.f90
	$(F90) -c $(FOPTS) machine.f90

#
# $Log$
# Revision 2.5  2000/10/11 18:23:13  vsnyder
# Insert copyright notice
#
# Revision 2.4  2000/10/11 18:04:57  vsnyder
# Put rule for machine.o back in -- it was erroneously deleted
#
# Revision 2.3  2000/10/10 22:44:25  vsnyder
# Correct actions for test_sparsify and test_multiply
#
# Revision 2.2  2000/10/10 00:12:47  vsnyder
# Corrected an incorrect path for CF
#
# Revision 2.1  2000/09/13 20:41:07  vsnyder
# Support to compile dump_0.o
#
# Revision 2.0  2000/09/05 17:41:31  dcuddy
# Change revision to 2.0
#
# Revision 1.1  2000/09/02 01:50:26  vsnyder
# Initial installation
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
