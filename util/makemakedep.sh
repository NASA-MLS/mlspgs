#!/bin/sh
#makemakedep.sh

# --------------- makemakedep.sh help
#Creates file Makefile.dep of dependencies to be included by Makefile
#to compile a Fortran 9x program or library among files.
#First a list of OBJS is assembled: this is the list of
#objects that need to be built. By default this list is based on
#all files in the current directory matching the pattern %.f90,
#and/or, depending on options set on the command-line,
#files matching patterns %.f and %.c. In other words, source files.
#Second, dependencies are calculated for each of the OBJS,
#searching for prerequisites in, by default, the current directory,
#and, depending on command line args, extra search directories
#
#A separate possible use: assembling lists of dependencies
#to "compile" source files with suffixes matching one pattern 
#to create "objects" files with suffixes matching another
#Under this use, the orthodox patterns "%.f90", etc. will be ignored
#Usage:
#makemakedep.sh [opt1] [opt2] ..  [dir1 dir2 ..]
#
#    O p t i o n s
# -[n]f90       [don't] include files with .f90 extensions
#                (default is to include them)
# -[n]f77       [don't] include files with .f extensions, too
#                (default is to include them)
# -[n]c         [don't] include files with .c extensions, too
#                (default is to include them)
# -d file_name  exclude file_name
# -s pattern    match source file suffixes against "%.pattern"
# -o pattern    match object file suffixes against "%.pattern"
# -h[elp]       print brief help message; exit
# -nodep        skip tracing dependencies; stop after list of objects
# dir1          search directory named by dir1 as well as cwd for files
#
#
#Note:
#The option(s) marked with "-", if present,
#must precede the extra search directories on the command line
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
# --------------- End makemakedep.sh help
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#----------------------- UserPrompt -----------------------

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

#------------------------------- extant_files ------------
#
# Function to return only those files among the args
# that actually exist
# Useful when passed something like *.f which may 
# (1) expand to list of files, returned as extant_files_result, or
# (2) stay *.f, in which case a blank is returned as extant_files_result 
#     (unless you have perversely named a file '*.f')
# usage: extant_files arg1 [arg2] ..

extant_files()
{
   extant_files_result=
   # Trivial case ($# = 0)
   if [ "$1" != "" ]
   then
      for file
      do
         if [ -f "$file" ]
         then
               extant_files_result="$extant_files_result $file"
         fi
      done
   fi
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
#this version uses either of the two depmakers
#to trace dependencies based on USEs and INCLUDEs (or none at all):
#(1)  makedepf90 (a compiled program)
#(2)  f90makedep.pl (a perl script)
#(3)  none (skips tracing dependencies)
#
#which of (1) or (2) to use is set by the following line
DEPMAKER=2
#        ^  -- set this to 1 for makedepf90, 2 for f90makedep.pl
# (option (3) is selected by 
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
ACT_COURTEOUS=0
#             ^  -- set this to 1 to rename older Makefile.dep, 0 deletes it
#
TRY_CLEANUP=1
#           ^  -- set this to 1 to try cleaning up from a prior faulty run
#
PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience of users
#
EDIT_GB_PERL_PATH=1
#                 ^  -- set this to 1 to change path to perl in ghostbuster.sh
#

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
         extant_files "$dsuffix"/*
         if [ "$extant_files_result" != "" ]
         then
			   mv "$dsuffix"/* .
         fi
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
# The initial settings are
# include *.f90   yes
# include *.f     yes
# include *.c     yes
include_f90="yes"
include_f77="yes"
include_c="yes"
orthodox="yes"
s_pattern="f90"
o_pattern="o"
me="$0"
my_name=makemakedep.sh

wrong_list=""
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -f77 )
       include_f77="yes"
       shift
       ;;
    -nf77 )
       include_f77="no"
       shift
       ;;
    -f90 )
       include_f90="yes"
       shift
       ;;
    -nf90 )
       include_f90="no"
       shift
       ;;
    -c )
       include_c="yes"
       shift
       ;;
    -nc )
       include_c="no"
       shift
       ;;
#                  Rename or hide excluded files so they !~= %.f90
    -d )
	    if [ -f "$2" ]
       then
	       if [ "$f90suffix" != "$dsuffix" ]
          then
#            mv "$2" "`echo $2 | sed 's/'$f90suffix'/'$dsuffix'/'`"
	          mv "$2" "$dsuffix"
	       else
	          echo "$2 wrongly added to dependency lists"
                            wrong_list="$wrong_list $2"
	       fi
       fi
       shift
	    shift
       ;;
    -s )
	    s_pattern="$2"
       shift
	    shift
       orthodox="no"
       include_f90="no"
       include_f77="no"
       include_c="no"
       ;;
    -nodep )
	    DEPMAKER=0
       shift
       ;;
    -o )
	    o_pattern="$2"
       shift
	    shift
       orthodox="no"
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;

    * )
       more_opts="no"
       ;;
    esac
done

if [ $PRINT_TOO_MUCH = "1" ]                            
then                                                    
   echo " Summary of options to $me "  
   echo " include_f90: $include_f90 "  
   echo " include_f77: $include_f77 "  
   echo " include_c: $include_c "  
   echo " orthodox: $orthodox "  
   echo " wrong list: $wrong_list "  
   echo " s_pattern: $s_pattern "  
   echo " o_pattern: $o_pattern "  
   echo " ACT_COURTEOUS: $ACT_COURTEOUS "  
   echo " TRY_CLEANUP: $TRY_CLEANUP "  
   echo " DEPMAKER: $DEPMAKER "  
   echo " EDIT_GB_PERL_PATH: $EDIT_GB_PERL_PATH "  
fi                                                      

#                Create Makefile.dep; write 1st line
echo "#Makefile.dep -- a file to be included by a Makefile" > Makefile.dep
if [ "$orthodox" = "yes" ]
then
   echo "#to compile a Fortran 9x program or library" >> Makefile.dep
else
   echo "#to hold automatically calculated dependencies" >> Makefile.dep
fi
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

if [ "$include_f90" = "yes" ] ; then
   extant_files *.f90
   (echo $extant_files_result | sed 's/\.f90/.o  /g; s/$/\\/') >> Makefile.dep
fi

if [ "$include_f77" = "yes" ] ; then
   extant_files *.f
   (echo $extant_files_result | sed 's/\.f/.o  /g; s/$/\\/') >> Makefile.dep
fi

if [ "$include_c" = "yes" ] ; then
   extant_files *.c
   (echo $extant_files_result | sed 's/\.c/.o  /g; s/$/\\/') >> Makefile.dep
fi

if [ "$orthodox" != "yes" ] ; then
   extant_files *.$s_pattern
   (echo $extant_files_result | sed 's/\.'$s_pattern'/.'$o_pattern'  /g; s/$/\\/') >> Makefile.dep
fi

echo " "  >> Makefile.dep
#
if [ $DEPMAKER = "0" ]
then
  echo "#(Not tracing dependencies)" >> Makefile.dep
  echo "#End of simplified Makefile.dep" >> Makefile.dep
elif [ $DEPMAKER = "1" ]
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
   if [ "$orthodox" = "yes" ]
   then
	   if [ $PRINT_TOO_MUCH = "1" ]
	   then
         $the_DEPMAKER "$@"
      fi
      $the_DEPMAKER "$@" >> Makefile.dep
   else
      echo $the_DEPMAKER -s "$s_pattern" -o "$o_pattern"
	   if [ $PRINT_TOO_MUCH = "1" ]
	   then
         $the_DEPMAKER -s "$s_pattern" -o "$o_pattern"
      fi
      $the_DEPMAKER -s "$s_pattern" -o "$o_pattern" >> Makefile.dep
   fi
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
# Revision 1.16  2002/04/09 19:38:00  pwagner
# Changes allow it to compute m4 dependencies to make l2cf
#
# Revision 1.15  2001/11/06 00:19:26  pwagner
# Turned off ACT_COURTEOUS
#
# Revision 1.14  2001/08/14 16:05:08  pwagner
# Added options; defaults to include_ = yes
#
# Revision 1.13  2001/08/10 23:48:27  pwagner
# Cosmetic changes only
#
# Revision 1.12  2001/08/10 17:42:33  pwagner
# Added -h(elp) option; general housekeeping
#
# Revision 1.11  2001/08/09 23:35:11  pwagner
# Omits bogus *.o lines from Makefile.dep
#
# Revision 1.10  2001/05/30 17:37:51  pwagner
# Added include_f77
#
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

