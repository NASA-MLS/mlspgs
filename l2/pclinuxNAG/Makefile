# $Id$

# End of the path to this Makefile:
ME = pclinuxNAG

# Fortran compiler:
F90 = f95

# Compiler options:
FOPTS = -g -C=all
#FOPTS = -O

# Path to source files:
S = ..

# Path to MLS library source files:
MLS = $S/../lib

# Path to .mod files:
M = $(MLS)/$(ME)

# Option to specify path to .mod files:
I = -I

# Beginning of path to libraries:
TOOLKIT = /software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Beginning of path to HDFEOS libraries:
HDFEOS = /van1/opt/TOOLKIT
#HDFEOS = $(TOOLKIT)

# Platform (used in a directory path in ../Makefile.inc):
PLATFORM = linux

# Path to libPGSTK
PGSTK = $(HOME)/mlspgs/SDPToolkit/$(ME)
#PGSTK = $(TOOLKIT)/lib/$(PLATFORM)

include $S/Makefile.inc

# $Log$
# Revision 2.2  2000/10/12 00:38:10  vsnyder
# Account for lib/cf_parser moving to lib
#
# Revision 2.1  2000/10/10 20:33:48  vsnyder
# Repairs to get rid of WVS pathnames
#
