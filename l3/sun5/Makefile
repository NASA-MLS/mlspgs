# Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# End of path to this Makefile:
ME = sun5

# Fortran compiler:
F90 = f95

# Compiler options:
FOPTS = -g

# DUSTY is used to allow picky compilers to handle references to "sfrdata",
# which has an inconsistent interface in that it can have any type
# object for its fifth argument.

DUSTY=-dusty

# Loader options:
LDOPTS = -g

# Option to specify path to .mod and include files:
I = -I

# Path to source files
S=..

# Path to init_tables_module.f90:
A = ..

# Path to MLS library sources:
MLS = $S/../lib

# Path to MLS library .o and .mod files:
M = $(MLS)/$(ME)

# Path to MLS source library:
MLSS = $S/../srclib

# Beginning of path to toolkit libraries:
TOOLKIT = /software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Beginning of path to HDFEOS libraries:
HDFEOS = $(TOOLKIT)

# Platform (used within directory paths in ../Makefile.inc):
PLATFORM = sun5

# Path to libPGSTK.a
PGSTK = $(TOOLKIT)/lib/$(PLATFORM)

include $S/Makefile.inc

#
# $Log$
# Revision 1.1  2000/10/24 19:49:44  nakamura
# Customized for sun5/NAG make of L3.
#
#
