#!/bin/sh
    #makemakedep.sh
#Creates file Makefile.dep of dependencies to be included by Makefile
#to compile a Fortran 9x program or library
#this version uses either of the two depmakers
#to trace dependencies based on USEs and INCLUDEs:
#(1)  makedepf90 (a compiled program)
#(2)  f90makedep.pl (a perl script)
#
#which of (1) or (2) to use is set by the following line
DEPMAKER=2
#        ^  -- set this to 1 for makedepf90, 2 for f90makedep.pl
#
#unless you have obtained and compiled makedepf90, you will almost
#certainly need to set DEPMAKER to 2 in the above line
#If you prefer to use makedepf90 (perhaps because you lack perl)
#it is available from
#        http://www.helsinki.fi/~eedelman/makedepf90.html
#as a compiled binary for Linux and as sources for other platforms
#After compiling makedepf90, you will need to place it in your PATH,
#perhaps by moving a copy to $(HOME)/bin if it exists,
#or else replace the line
# 	makedepf90 | sed ...
# below with
# 	/dir1/dir2/../dirn/makedepf90 | sed ...
#where /dir1.. is directory where you created makedpf90 (use 'pwd')
#
#if you use f90makedep.pl, you may have a problem if the path
#your copy of perl is different from the one in its 1st line
#compare 'which perl' with 'sed -n "1 p" f90makedep.pl
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
#
if [ $DEPMAKER = "1" ]
then
	#
	# use makedepf90 to calculate dependencies
	echo " using makedepf90 to calculate dependencies "
	echo "# using makedepf90 to calculate dependencies "  >> Makefile.dep
	#
	#makedepf90 *.f90 | sed 's/^makedepf90/#makedepf90/' >> Makefile.dep
	#
	#Pipe through sed to overcome an apparent bug in makedepf90
	makedepf90 *.f90 | sed 's/\.f90\.o/\.o/g' >> Makefile.dep
else
	#
	# use f90makedep.pl to calculate dependencies
	echo " using f90makedep.pl to calculate dependencies "
	echo "# using f90makedep.pl to calculate dependencies "  >> Makefile.dep
	#
	# Prefix f90makedep.pl with the path to util
        # which is assumed to be the same as the path to this script
	the_DEPMAKER="`echo $0 | sed 's/makemakedep.sh/f90makedep.pl/'`"
	echo " Your perl is `which perl` "
	echo " f90makedep.pl is looking for it at `sed -n '1 p' $the_DEPMAKER`"
#	f90makedep.pl >> Makefile.dep
	$the_DEPMAKER >> Makefile.dep
fi
echo " "  >> Makefile.dep
echo "#End of Makefile.dep" >> Makefile.dep
echo " "  >> Makefile.dep

# $Log$
# Revision 1.2  2000/10/18 20:31:58  pwagner
# Fix to overcome apparent bug in makedepf90
#
# Revision 1.1  2000/10/16 18:31:49  pwagner
# first commit
#

