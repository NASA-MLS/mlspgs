# $Id$

# End of the path to this Makefile:
ME = LF95

# Fortran compiler:
F90 = lf95

# Compiler options:
FOPTS = -g --chk --info
#FOPTS = -O

# Path to source files:
S = ..

# Path to MLS library source files:
MLS = $S/../lib

# Path to .mod files:
M = $(MLS)/$(ME)

# Option to specify path to .mod files:
I = -I

# Beginning of path to toolkit libraries:
TOOLKIT = /software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Beginning of path to HDFEOS libraries:
HDFEOS = /van1/opt/TOOLKIT
#HDFEOS = $(TOOLKIT)

# Platform (used within directory paths in ../Makefile.inc):
PLATFORM = linux

# Path to libPGSTK.a
PGSTK = $(HOME)/mlspgs/SDPToolkit/$(ME)
#PGSTK = $(TOOLKIT)/lib/$(PLATFORM)

include $S/Makefile.inc

# $Log$
# Revision 2.1  2000/10/10 20:33:32  vsnyder
# Repairs to get rid of WVS pathnames
#
