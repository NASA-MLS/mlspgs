# $Id$

# Fortran compiler:
FC = f95

# Compiler options:
FOPTS = -g

# Loader options:
LDOPTS = -g

# Path to source files:
S = ..

# Option to specify path to .mod files:
I = -I

# Path to .mod files for MLS library:
M = $(HOME)/mlspgs/lib/NAG

# SDPToolkit or SDPToolkitSubstitute?
#SDP = 
SDP = Substitute

machine.o: machine.f90
	$(FC) -c $(FOPTS) machine.f90

include ../Makefile.inc

# $Log$
# Revision 2.2  2000/10/09 23:25:07  vsnyder
# Repair SDP variable
#
# Revision 2.1  2000/10/09 23:17:41  vsnyder
# Search more paths for modules, because second-order USEs need to be
# accomodated by specifying paths for them.
#
# Revision 2.0  2000/09/05 17:42:08  dcuddy
# Change revision to 2.0
#
# Revision 1.2  2000/07/06 18:04:08  vsnyder
# Insert CVS "variables"
#
