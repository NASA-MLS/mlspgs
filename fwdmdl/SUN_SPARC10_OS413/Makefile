# $Id$

# This is the main makefile for the SUN Compiler

F90=f90
F90DEBUGFLAGS=-g
#MODULEPATH="-I$(HOME)/mlspgs/lib/SUN_SPARC10_OS413 -I."
MODULEPATH=-M.
#INCLUDE=-I/usr/local/TOOLKIT/include
INCLUDE=
VPATH=..

include ../Makefile.inc

machine.o: machine.f90
	$(F90) -c machine.f90

# $Log$