# $Id$

# Places to look for sources without being explicitly told to do so:
VPATH=..

# Fortran compiler:
F90=f95

# Compiler options:
FOPTS=-g

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

# Path to cf_parser
CF=$S/cf_parser/NAG

# Beginning of path to toolkit libraries:
TOOLKIT = /nas/software/SRC/TOOLKIT-5.2.6v1.00.NAG

include ../Makefile.inc

#
# $Log$
# Revision 2.1  2000/10/10 00:13:04  vsnyder
# Corrected an incorrect path for CF
#
# Revision 2.0  2000/09/13 23:07:59  ahanzel
# Removed old logs and changed file revision to 2.0.
#
# Revision 1.2  2000/09/13 20:41:10  vsnyder
# Support to compile dump_0.o
#
