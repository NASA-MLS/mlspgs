# $Id$
#
# This is an include file for makefiles that is intended to be as
# platform independent as possible.  Where things are platform dependent,
# an if statement is used.
#

# DUSTY is used to allow picky compilers to handle references to "sfrdata",
# which has an inconsistent interface in that it can have any type
# object for its fifth argument.

DUSTY=-dusty
F90=f95
F90DEBUGFLAGS=-g
INCLUDE=-I/usr/local/TOOLKIT/include
MODULEPATHPREFIX=-I
VPATH=..

# Platform path
P=NAG

# Path to source files
S=..

# Path to cf_parser
C=$S/cf_parser/$P

include ../Makefile.inc

#
# $Log$
# Revision 2.0  2000/09/13 23:07:59  ahanzel
# Removed old logs and changed file revision to 2.0.
#
# Revision 1.2  2000/09/13 20:41:10  vsnyder
# Support to compile dump_0.o
#
