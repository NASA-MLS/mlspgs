# $Id$

# This is the main makefile for Lahey/Fujitsu Fortran 95.

F90=lf95
F90FLAGS = -O --f95
F90DEBUGFLAGS=-g
#MODULEPATH="-I$(HOME)/mlspgs/lib/LF95 -I."
MODULEPATH=-I.
#INCLUDE=-I/usr/local/TOOLKIT/include
INCLUDE=
VPATH=..

include ../Makefile.inc

machine.o: machine.f90
	$(F90) -c machine.f90

# $Log$
# Revision 1.3  2000/05/05 17:19:39  vsnyder
# Replaced because CVS lost v1.1
#
# Revision 1.1  2000/05/04 18:12:13  vsnyder
# Initial conversion to Fortran 90
#
