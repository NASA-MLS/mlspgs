# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# End of path to this Makefile:
ME = LF95

# Fortran compiler:
F90 = lf95

# Compiler options:
FOPTS =- g --chk --info

# Loader options:
LDOPTS = -g --chk

# Option to specify path to .mod and include files:
I = -I

# Path to source files for test harness:
S = ..

# Path to init_tables_module.f90:
A = $S/../../l2

# Path to MLS library sources:
MLS = $S/..

# Path to MLS library .o and .mod files:
M = $(MLS)/$(ME)

# Beginning of path to toolkit libraries:
TOOLKIT = /nas/software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Platform (used within directory paths in ../Makefile.inc):
PLATFORM = linux

# Path to libPGSTK.a
PGSTK = $(HOME)/mlspgs/SDPToolkit/$(ME)
#PGSTK = $(TOOLKIT)/lib/$(PLATFORM)

include $S/Makefile.inc

# $Log$
# Revision 1.1  2000/10/11 20:56:40  vsnyder
# Initial entry
#
