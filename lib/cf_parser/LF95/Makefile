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

machine.o: machine.f90
	$(FC) -c $(FOPTS) machine.f90

include ../Makefile.inc

# $Log$
