#!/bin/sh
    #makemakedep.sh
#Creates file Makefile.dep of dependencies to be included by Makefile
#to compile a Fortran 9x program or library among files =~ the pattern %.f90
#
#All files and Makefile.dep are in the current working directory
#except in case of (2) (see below) any args are passed through to
#f90makedep.pl where they will be interpreted as additional directories
#in which to search for additional prerequisite modules that may be among
#the lists of dependencies (please see f90makedep.pl for details)
#except for any args preceded by the "-d" flag (all of which should
#come first, before any of the directory args)
#Arguments preceded by the "-d" flag, each of which must appear
#separately and with its own flag, are Fortran 9x files in the
#current working directory to be excluded from the list of dependencies
#(often because they are special-purpose, "test" programs, or some such)
#
#Usage:
#makemakedep.sh [-f77] [-d file1.f90 -d file2.f90] [arg1 arg2 ..]
#
#Note:
#The option, -f77, if present must be the first option on the command line
#Result:
#creates Makefile.dep in current working directory
#Makefile.dep may then be inserted in a generic Makefile which
#thereby is able to keep current with a "moving target" in case
#the source files are changed in a way that alters their dependencies
#E.g. examine the following mock session
#  mlspgs/l2%ls
#source1.f90  source2.f90  source3.f90
#  mlspgs/l2%makemakedep.sh
#  mlspgs/l2%ls
#Makefile.dep  source1.f90  source2.f90  source3.f90
#  mlspgs/l2%cat Makefile.dep
##Makefile.dep -- a file to be included by a Makefile
##to compile a Fortran 9x program or library
#OBJS = \
#source1.o               source2.o          source3.o  \
#
#source1.o: source1.f90 source2.o source3.o
#source2.o: source2.f90 source3.o
#source3.o: source3.f90
# 
##End of Makefile.dep
# 
#Options
# -f77   include files with .f extensions as well as .f90
#         (otherwise only .f90 files considered)
# -d file_name  exclude file_name
# arg1   search directory named by arg1 as well as cwd for files
#
# Function to prompt for user response
#

UserPrompt()
{
    if [ "$BRAND" = "linux" ] ; then
	/bin/echo -n "$* " > /dev/tty
    else
	echo "$* \\c" > /dev/tty
    fi
    read user_response
}

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
#where /dir1.. is directory where you created makedepf90 (use 'pwd')
#
#if you use f90makedep.pl, you may have a problem if the path where you keep
#your copy of perl is different from the one in its 1st line
#compare 'which perl' with 'sed -n "1 p" f90makedep.pl
#The script will attempt ot anticipate this problem and ask you if it
#should fix f90makedep.pl
#
#makemakedep.sh may act courteously in the following sense:
#if makemakedep.sh finds there is already a file named Makefile.dep
#it will attempt to rename the older file rather than deleting it
#
ACT_COURTEOUS=1
#             ^  -- set this to 1 to rename older Makefile.dep, 0 deletes it
#
TRY_CLEANUP=1
#           ^  -- set this to 1 to try cleaning up from a prior faulty run
#
PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience
#
EDIT_GB_PERL_PATH=1
#                 ^  -- set this to 1 to change path to perl in ghostbuster.sh
#
# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#----------------------- Implementation -----------------------
#
#           How to rename or hide excluded files so they !~= %.f90
dsuffix=".xui"
#         ^^^----- this is the suffix stuck onto any excluded files
#                   or else the name of a temp directory hiding them
# Do we have write permission in the current working directory
if [ -w "`pwd`" ]
then
	f90suffix=".f90"
	if [ -d "$dsuffix" ] ; then
		if [ $TRY_CLEANUP = "1" ] ; then
			mv "$dsuffix"/* .
			rmdir "$dsuffix"
		else
			echo "Sorry--$dsuffix already exists"
			echo "Aborting new Makefile.dep"
                exit
      fi
    fi
    mkdir "$dsuffix"
else
	f90suffix="$dsuffix"
fi

#   Rename older Makefile.dep if one already exists
if [ -f Makefile.dep ]
then
	if [ $ACT_COURTEOUS = "1" ]
	then
		if [ $PRINT_TOO_MUCH = "1" ]
		then
			echo "Renaming older Makefile.dep Make.dep.n"
		fi
		name=Make.dep.`ls -1 Make*.dep* | wc -l`
		cname=`echo $name | sed 's/ //'g`
		mv Makefile.dep $cname
        else
                chmod u+w Makefile.dep
                rm -f Makefile.dep
	fi
fi
#

if [ "$1" = "-f77" ] ; then
   include_f77="yes"
   shift
else
   include_f77="no"
fi

#                  Rename or hide excluded files so they !~= %.f90
wrong_list=""
while [ "$1" = "-d" ] ; do
	if [ -f "$2" ]
        then
		if [ "$f90suffix" != "$dsuffix" ]
        	then
#        		mv "$2" "`echo $2 | sed 's/'$f90suffix'/'$dsuffix'/'`"
			mv "$2" "$dsuffix"
		else
			echo "$2 wrongly added to dependency lists"
                        wrong_list="$wrong_list $2"
	        fi
        fi
        shift
	shift

done
#
#                Create Makefile.dep; write 1st line
echo "#Makefile.dep -- a file to be included by a Makefile" > Makefile.dep
echo "#to compile a Fortran 9x program or library" >> Makefile.dep
#
#Warn of files wrongly added to dependency lists
if [ "$wrong_list" != "" ]
then
	echo "#Delete the following files from the dependency lists:" >> Makefile.dep
        echo "#$wrong_list" >> Makefile.dep
fi
#The following may be useful in rare circumstances to include special macros
#but violates assertion that Makefile.dep includes only dependency info
#Therefore uniformity and simplicity suggest against relying on it
#e.g., the mlspgs software does not require such a file "Makefile.mac"
#
if [ -f Makefile.mac ]
then
	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo "Including special macro definitions from Makefile.mac"
	fi
	echo "include ../Makefile.mac"  >> Makefile.dep
fi
#
echo "OBJS = \\"  >> Makefile.dep

# Note:
# The following won't work if there are no files matching the specified
# extensions in the directory;
# instead the echo will just produce a bogus "*.f90" string
if [ "$include_f77" = "yes" ] ; then
   (echo *.f90 | sed 's/\.f90/.o  /g; s/$/\\/') >> Makefile.dep
   (echo *.f | sed 's/\.f/.o  /g; s/$/\\/') >> Makefile.dep
else
   (ls -C *.f90 | sed 's/\.f90/.o  /g; s/$/\\/') >> Makefile.dep
fi
echo " "  >> Makefile.dep
#
if [ $DEPMAKER = "1" ]
then
	#
	# use makedepf90 to calculate dependencies
	# this assumes makedepf90 both exists && is in PATH
	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo " using makedepf90 to calculate dependencies "
	fi
	echo "# using makedepf90 to calculate dependencies "  >> Makefile.dep
	#
	#makedepf90 *.f90 | sed 's/^makedepf90/#makedepf90/' >> Makefile.dep
	#
	#Pipe through sed to overcome an apparent bug in makedepf90
	makedepf90 *.f90 | sed 's/\.f90\.o/\.o/g' >> Makefile.dep
else
	#
	# use f90makedep.pl to calculate dependencies
	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo " using f90makedep.pl to calculate dependencies "
	fi
	echo "# using f90makedep.pl to calculate dependencies "  >> Makefile.dep
	#
	# Prefix f90makedep.pl with the path to util
        # this assumes f90makedep.pl is in same dir as this script
	the_DEPMAKER="`echo $0 | sed 's/makemakedep.sh/f90makedep.pl/'`"
	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo " Your perl is `which perl` "
		echo " f90makedep.pl is looking for it at `sed -n '1 p' $the_DEPMAKER`"
	fi

# Check whether script is looking for perl in right place
# && give user a chance to redirect it if it is not
	script_perl=`sed -n '1 p' $the_DEPMAKER`
	your_perl='#!'`which perl` 
   if [ "$script_perl" != "$your_perl" ]
        then
		#Warn user that perl script may need to be changed

		echo " *** Warning: f90makedep.pl may need to be changed"
		echo " ***          to point to where your perl actually is"
		echo " ***          "
		echo " ***  (unless you know of a compelling reason to do otherwise"
		echo " ***  you probably want to answer 'yes' to the following)"
		echo " ***          "
 		UserPrompt "Change f90makedep.pl to look in [$your_perl](yes) or no?"
   		if [ "$user_response" != "" ] ; then
			case  "$user_response" in
	    		n* | N* )
				echo "Continuing to use $script_perl"
              		  	change_perl=
	   		 ;;
	    		y* | Y* )
				echo "Changing to use $your_perl"
                		change_perl="$your_perl"
	    		;;
	   		 * )
				echo "Changing to use $user_response"
                		change_perl="$user_response"
	    		;;
			esac
			else
				echo "Changing to use $your_perl"
                		change_perl="$your_perl"
			fi
   		if [ "$change_perl" != "" ] ; then
            sed -n "1 s%$script_perl%$change_perl%p;2,$ p" $the_DEPMAKER > temp.pl
				chmod u+w "$the_DEPMAKER"
         	mv temp.pl "$the_DEPMAKER"
				chmod a+x "$the_DEPMAKER"
				echo "*** You have fixed f90makedep.pl to look for $your_perl"
   		   if [ "$EDIT_GB_PERL_PATH" = "1" ] ; then
	            #
	            # Prefix f90GhostFiles.pl with the path to util
               # this assumes f90GhostFiles.pl is in same dir as this script
	            the_GHOSTFINDER="`echo $0 | sed 's/makemakedep.sh/f90GhostFiles.pl/'`"
            	script_perl=`sed -n '1 p' $the_GHOSTFINDER`
               if [ "$script_perl" != "$your_perl" ] ; then
                  sed -n "1 s%$script_perl%$change_perl%p;2,$ p" $the_GHOSTFINDER > temp.pl
		      		chmod u+w "$the_GHOSTFINDER"
         	      mv temp.pl "$the_GHOSTFINDER"
		      		chmod a+x "$the_GHOSTFINDER"
				      echo "*** You have also fixed f90GhostFiles.pl to look for $your_perl"
               fi
            fi
         fi
       fi
#	f90makedep.pl >> Makefile.dep
	$the_DEPMAKER "$@" >> Makefile.dep
fi
echo " "  >> Makefile.dep
echo "#End of Makefile.dep" >> Makefile.dep
echo " "  >> Makefile.dep

#         clean up
# renamed or hidden files
if [ -w "$dsuffix" ]
then
	moved_files_list="$dsuffix"/*
	moved_files=`echo $moved_files_list`
#	the above will expand the wild card * if dsuffix is non-empty
	if [ "$moved_files" != "$dsuffix/*" ]
	then
        	mv "$dsuffix"/* .
	fi
        rmdir "$dsuffix"
fi
exit
# $Log$
# Revision 1.9  2001/05/01 17:04:32  pwagner
# Also edits perl path of f90GhostFiles.pl
#
# Revision 1.8  2001/03/23 00:41:23  pwagner
# Prints less
#
# Revision 1.7  2001/02/26 00:25:40  pwagner
# More informative message before fixing perl path; cleanup after abort
#
# Revision 1.6  2000/11/21 00:47:34  pwagner
# Warns user if script's perl not which's; allows change
#
# Revision 1.5  2000/11/14 21:50:14  pwagner
# Can exclude specified files from dependencies
#
# Revision 1.4  2000/11/02 23:22:38  pwagner
# Dependencies may cross directories
#
# Revision 1.3  2000/10/27 23:22:37  pwagner
# works with f90makedep.pl
#
# Revision 1.2  2000/10/18 20:31:58  pwagner
# Fix to overcome apparent bug in makedepf90
#
# Revision 1.1  2000/10/16 18:31:49  pwagner
# first commit
#

