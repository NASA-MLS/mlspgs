#!/bin/sh
    #makemakedep.sh
#Creates file Makefile.dep of dependencies to be included by Makefile
#to compile a Fortran 9x program or library
#uses makedepf90 to trace dependencies based on USEs and INCLUDEs
#
#
# "$Id$"
if [ -f Makefile.dep ]
then
	echo "Renaming older Makefile.dep Make.dep.n"
	name=Make.dep.`ls -1 Make*.dep* | wc -l`
	cname=`echo $name | sed 's/ //'g`
	mv Makefile.dep $cname
fi
echo "#Makefile.dep -- a file to be included by a Makefile" > Makefile.dep
echo "#to compile a Fortran 9x program or library" >> Makefile.dep
if [ -f Makefile.mac ]
then
	echo "Including macro definitions in Makefile.mac"
	echo "include ../Makefile.mac"  >> Makefile.dep
fi
echo "OBJS = \\"  >> Makefile.dep
(ls -C *.f90 | sed 's/.f90/.o  /g; s/$/\\/') >> Makefile.dep
echo " "  >> Makefile.dep
makedepf90 *.f90 | sed 's/^makedepf90/#makedepf90/' >> Makefile.dep
echo " "  >> Makefile.dep
echo "#End of Makefile.dep" >> Makefile.dep
echo " "  >> Makefile.dep

# $Log$

