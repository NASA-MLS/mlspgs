# $Id$

# Fortran compiler:
FC = lf95

# Compiler options:
FOPTS = -g --chk

# Loader options:
LDOPTS = -g --chk

# Path to source files:
S = ..

# Option to specify path to .mod files:
I = -I

# Path to .mod files for MLS library:
M = $(HOME)/mlspgs/lib/LF95

# SDPToolkit or SDPToolkitSubstitute?
#SDP = 
SDP = Substitute

machine.o: machine.f90
	$(FC) -c $(FOPTS) machine.f90

include ../Makefile.inc

# $Log$
# Revision 2.0  2000/09/05 17:42:01  dcuddy
# Change revision to 2.0
#
# Revision 1.2  2000/07/06 18:04:05  vsnyder
# Insert CVS "variables"
#
