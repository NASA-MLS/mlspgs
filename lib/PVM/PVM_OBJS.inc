# PVM_OBJS.inc
# A file of Makefile statements that will be included
# in the Makefile if the PVM library is to be optionally linked

# "$Id$"

PVM_OBJS=PVMWrappers.o

vpath %.c ../ $(srclib)

PVMWrappers.o: $(S)/PVMWrappers.c
	cc -c -I${PGSINC} -DNAGf90Fortran $(S)/PVMWrappers.c

# $Log$
