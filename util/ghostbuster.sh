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
# -d file_name  exclude file_name
# -h[elp]       print brief help message; exit
# -lib   if this flag is present, executes "rm -f *.a" if any ghosts found
#Note:
#The option(s) marked with "-", if present,
#must precede the extra search directories on the command line
# 
# --------------- End ghostbuster.sh help
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

#

DEBUG=0
#     ^  -- set this to 1 to echo rather than execute rm command

PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience
#
# Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#----------------------- Implementation -----------------------
#
#if [ "$1" = "-lib" ] ; then
#   rm_any_libs="yes"
#   shift
#else
#   rm_any_libs="no"
#fi

TRY_CLEANUP=1
#           ^  -- set this to 1 to try cleaning up from a prior faulty run
#
#           How to rename or hide excluded files so they !~= %.f90
dsuffix=".xug"
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

me="$0"
my_name=ghostbuster.sh

wrong_list=""
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
		echo " Your perl is `which perl` "
		echo " f90GhostFiles.pl is looking for it at `sed -n '1 p' $the_GHOSTFINDER`"
	fi

# Check whether script is looking for perl in right place
# && give user a chance to redirect it if it is not
	script_perl=`sed -n '1 p' $the_GHOSTFINDER`
	your_perl='#!'`which perl` 
   if [ "$script_perl" != "$your_perl" ]
        then
		#Warn user that perl script may need to be changed

		echo " *** Warning: f90GhostFiles.pl may need to be changed"
		echo " ***          to point to where your perl actually is"
		echo " ***          "
		echo " ***  (unless you know of a compelling reason to do otherwise"
		echo " ***  you probably want to answer 'yes' to the following)"
		echo " ***          "
 		UserPrompt "Change f90GhostFiles.pl to look in [$your_perl](yes) or no?"
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
            sed -n "1 s%$script_perl%$change_perl%p;2,$ p" $the_GHOSTFINDER > temp.pl
				chmod u+w "$the_GHOSTFINDER"
         	mv temp.pl "$the_GHOSTFINDER"
				chmod a+x "$the_GHOSTFINDER"
				echo "*** You have fixed f90GhostFiles.pl to look for $your_perl"
         fi
       fi
	the_ghosts=`$the_GHOSTFINDER "$@"`

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
	fi
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
exit
# $Log$
# Revision 1.2  2001/05/31 20:16:26  pwagner
# Added -lib option to rm any haunted libs
#
# Revision 1.1  2001/05/01 17:03:31  pwagner
# First commit
#
#

