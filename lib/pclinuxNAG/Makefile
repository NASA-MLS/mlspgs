# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# Places to look for sources without being explicitly told to do so:
VPATH = ..

# End of the path to the directory containing this Makefile:
ME = pclinuxNAG

# Fortran compiler:
F90 = f95

# Compiler options:
FOPTS = -g

# DUSTY is used to allow picky compilers to handle references to "sfrdata",
# which has an inconsistent interface in that it can have any type
# object for its fifth argument.

DUSTY = -dusty

# Loader options:
LDOPTS = -g

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
# Revision 2.5  2000/10/12 22:05:28  vsnyder
# Remove CF variable, introduce ME variable.
#
# Revision 2.4  2000/10/11 18:23:15  vsnyder
# Insert copyright notice
#
# Revision 2.3  2000/10/11 18:05:09  vsnyder
# Put rule for machine.o back in -- it was erroneously deleted
#
# Revision 2.2  2000/10/10 22:44:26  vsnyder
# Correct actions for test_sparsify and test_multiply
#
# Revision 2.1  2000/10/10 00:13:04  vsnyder
# Corrected an incorrect path for CF
#
# Revision 2.0  2000/09/13 23:07:59  ahanzel
# Removed old logs and changed file revision to 2.0.
#
# Revision 1.2  2000/09/13 20:41:10  vsnyder
# Support to compile dump_0.o
#
