# Fortran compiler:
FC = lf95

# Compiler options:
FOPTS = -g --chk

# Path to source files:
S = ..

# Path to .mod files:
M = $(HOME)/mlspgs/lib/LF95

# Option to specify path to .mod files:
I = -I

# Beginning of path to libraries:
L = /nas/software/SRC/TOOLKIT-5.2.6v1.00
#L = /nas/software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Beginning of path to HDFEOS  libraries:
E = /van1/opt/TOOLKIT

# Path to Configuration parser
C = $(HOME)/mlspgs/lib/cf_parser/LF95

# Platform (used in a directory path in ../Makefile.inc):
P = linux

machine.o: machine.f90
	$(FC) -c $(FOPTS) machine.f90

include ../Makefile.inc
