#!/bin/sh
    #ghostbuster.sh
#rm ghosts from the binary directory (the 1st arg)
#where ghosts are defined as .o and .mod-suffixed
#files for which no corresponding source file currently
#exists in any of the source directories (arg# 2 ..)
#The most likely reason being that the original source file
#has been moved or deleted, yet the persistent .o or .mod
#file still "haunts" the binary-directory.
#
#Usage:
#ghostbuster.sh arg1 arg2 [arg3 ..]
#where argn is a path, either absolute or relative to current working directory
#Result:
#rm the list of such files, unless it is empty in which case does nothing
# 
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

#

DEBUG=0
#     ^  -- set this to 1 to echo rather than execute rm command

PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience
#
# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#----------------------- Implementation -----------------------
#
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
      else
         	echo "Sorry, the ghosts are away--may I take a message?"
	   fi
      exit
	fi

	if [ "$the_ghosts" != "" ]
	then
        	cd "$1"
        	rm -f $the_ghosts
	fi
exit
# $Log$
#

