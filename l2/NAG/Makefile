# Fortran compiler:
FC = f95

# Compiler options:
FOPTS = -g -C=all

# Path to source files:
S = ..

# Path to .mod files:
M = $(HOME)/mlspgs/lib/NAG

# Option to specify path to .mod files:
I = -I

# Beginning of path to libraries:
#L = /nas/software/SRC/TOOLKIT-5.2.6v1.00
L = /nas/software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Beginning of path to HDFEOS  libraries:
E = /van1/opt/TOOLKIT

# Path to Configuration parser
C = $(HOME)/mlspgs/lib/cf_parser/NAG

# Platform (used in a directory path in ../Makefile.inc):
P = linux

machine.o: machine.f90
	$(FC) -c $(FOPTS) machine.f90

include ../Makefile.inc
