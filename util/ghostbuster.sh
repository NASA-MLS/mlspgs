#!/bin/sh
    #ghostbuster.sh
# --------------- ghostbuster.sh help
#rm ghosts from the binary directory (the 1st arg)
#where ghosts are defined as .o and .mod-suffixed
#files for which no corresponding source file currently
#exists in any of the source directories (arg# 2 ..)
#The most likely reason being that the original source file
#has been moved or deleted, yet the persistent .o or .mod
#file still "haunts" the binary-directory.
#
#Usage:
#ghostbuster.sh [-lib] arg1 arg2 [arg3 ..]
#
#where argn is a path, either absolute or relative to current working directory
#Result:
#rm the list of such files, unless it is empty in which case does nothing
#
#
#    O p t i o n s
# -d file_name  ignore file_name (from current directory)
#                in effect, acting as though it is already dead
#                (this would be more useful if we could ignore
#                 files from each of the source directories)
# -h[elp]       print brief help message; exit
# -lib          if this flag is present, 
#                executes "rm -f *.a" if any ghosts found
# -perl program  
#               use program instead of perl to run f90GhostFiles.pl
#Note:
#The option(s) marked with "-", if present,
#must precede the extra search directories on the command line
# 
# --------------- End ghostbuster.sh help
# Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# Function to prompt for user response
#

UserPrompt()
{
    if [ "$AUTO_REPLY" = 1 ] ; then
      user_response='y'
    elif [ "$BRAND" = "linux" ] ; then
	/bin/echo -n "$* " > /dev/tty
    read user_response
    else
	echo "$* \\c" > /dev/tty
    read user_response
    fi
}

#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in pace of "."

get_unique_name()
{

   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
      fi
   # Is $HOST defined?
      if [ "$HOST" != "" ]
      then
         our_host_name="$HOST"
      elif [ "$HOSTNAME" != "" ]
      then
         our_host_name="$HOSTNAME"
      else
         our_host_name="host"
      fi
    #  echo $our_host_name
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
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
# Bug: -d file_name only ignores files in current directory
#         instead of every directory searched for sources
# Fix: Form list of bad files from reading command-line args
# then loop over extra search directories, mv-ing files as necessary
# then bust ghosts
# then loop again over extra search directories, cleaning up this time

#

DEBUG=0
#     ^  -- set this to 1 to echo rather than execute rm command

PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience
#
#
#----------------------- Implementation -----------------------

TRY_CLEANUP=1
#           ^  -- set this to 1 to try cleaning up from a prior faulty run
#
#           How to rename or hide excluded files so they !~= %.f90
#dsuffix=".xug"
#          ^^^----- this is the suffix stuck onto any excluded files
#                   or else the name of a temp directory hiding them
dsuffix=".`get_unique_name g`"
if [ $PRINT_TOO_MUCH = "1" ]
then
  echo "Hiding excluded files behind temp directory $dsuffix"
fi
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

me="$0"
my_name=ghostbuster.sh
PERL=perl

wrong_list=""
rm_any_libs="no"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -d )
	    if [ -f "$2" ]
       then
	       if [ "$f90suffix" != "$dsuffix" ]
          then
	          mv "$2" "$dsuffix"
	       else
	          echo "$2 wrongly added to dependency lists"
                            wrong_list="$wrong_list $2"
	       fi
       fi
       shift
	    shift
       ;;
    -lib )
       rm_any_libs="yes"
       shift
	    ;;
    -perl )
       PERL="$2"
       shift
	    shift
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

	if [ $# -lt "2" ]
	then
		echo " Usage: ghostbuster.sh arg1 arg2 [arg3 ..] "
		echo " rm ghosts from the binary directory (arg1)"
		echo " where ghosts are defined as .o and .mod-suffixed"
		echo " files for which no corresponding source file currently"
		echo " exists in any of the source directories (arg2 ..)"
      exit
	fi

	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo " using f90GhostFiles.pl to find ghosts "
	fi

	#
	# Prefix f90GhostFiles.pl with the path to util
        # this assumes f90GhostFiles.pl is in same dir as this script
	the_GHOSTFINDER="`echo $0 | sed 's/ghostbuster.sh/f90GhostFiles.pl/'`"
	if [ $PRINT_TOO_MUCH = "1" ]
	then
		echo " Your perl is $PERL "
	fi

	if [ "$DEBUG" = "1" ]
	then
      echo "About to call $the_GHOSTFINDER $@"
      echo "from directory `pwd`"
   fi
	the_ghosts=`$PERL $the_GHOSTFINDER "$@"`

   from_here=`pwd`
	if [ "$DEBUG" = "1" ]
	then
   	if [ "$the_ghosts" != "" ]
	   then
         	echo cd "$1"
         	echo "rm -f $the_ghosts"
	        if [ "$rm_any_libs" = "yes" ]
	        then
             echo	"rm -f *.a"
   	    fi
          echo cd "$from_here"
      else
         	echo "Sorry, the ghosts are away--may I take a message?"
	   fi
      exit
	fi

	if [ "$the_ghosts" != "" ]
	then
        	cd "$1"
        	rm -f $the_ghosts
	   if [ "$rm_any_libs" = "yes" ]
	   then
        	rm -f *.a
   	fi
      cd "$from_here"
	fi
#         clean up
# renamed or hidden files
if [ -w "$dsuffix" ]
then
         extant_files "$dsuffix"/*
   if [ "$extant_files_result" != "" ]
	then
        	mv "$dsuffix"/* .
	fi
   rmdir "$dsuffix"
fi
exit
# $Log$
# Revision 1.10  2004/10/27 22:34:08  pwagner
# Set AUTO_REPLY=1 to automate pointing to wherever perl relocates
#
# Revision 1.9  2003/06/21 00:34:17  pwagner
# Added fail-safe measures and MAY_EDIT_PERL flag
#
# Revision 1.8  2002/07/25 19:54:37  pwagner
# Comments on hiding excluded files only if PRINT_TOO_MUCH
#
# Revision 1.7  2002/07/22 22:08:44  pwagner
# Uses get_unique_name for temp file names
#
# Revision 1.6  2002/01/29 23:43:13  pwagner
# Cleans up excluded files from .xug dir before exiting
#
# Revision 1.5  2001/11/21 00:26:35  pwagner
# Fixed bug in ignored files
#
# Revision 1.4  2001/11/14 19:42:08  pwagner
# Added the stuff left out: extant_files code and -lib implementation
#
# Revision 1.3  2001/11/13 23:21:44  pwagner
# option -d bad_file_name added
#
# Revision 1.2  2001/05/31 20:16:26  pwagner
# Added -lib option to rm any haunted libs
#
# Revision 1.1  2001/05/01 17:03:31  pwagner
# First commit
#
#

