#!/bin/sh
#identl2cf.sh

# --------------- identl2cf.sh help
# Print to stdout idents for the l2cf fragments included in
# an l2cf file.  Takes the l2cf template as an arg. 

# Method:
# (1) Assemble list of l2cf fragments to be included
# (2) Extract '$Id' strings from fragments found in standard parts directory
# (3) Print strings to stdout

# Usage:
# identl2cf.sh [opt1] [opt2] ..  template
#
#     O p t i o n s
#    -dryrun              Merely echo the commands that would be executed
#    -I dir               Use dir as the standard parts directory
#                           default is ${HOME}/mlspgs/l2/l2cf/lib
#    -h[elp]              Print help; quit
#
# Bugs and limitations
# (1) It is assumed that all the l2cf fragments are included directly
#      from the l2cf fragment; it will miss any that are included indirectly
#      although we make a brave effort (see below)
# (2) All the l2cf fragments must be found in the standard parts directory
# (3) All the parts are included using syntax like
#     !include(name.l2cf)..
#     with possible exception of m4defs.l2cf
# (4) l2cf-formatted calibration files 
#     (like MLS-Aura_L2Cal-Tsys_v2-0-0_0000d000.l2cf)
#     are ignored (not explicitly, but because 
#     (a) they don't obey the syntax above, nor 
#     (b) they are not found in the standard parts directory
# --------------- End identl2cf.sh help
# Copyright 2009, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

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
#
files='m4defs.l2cf'
debug=0
#     ^  -- set this to 1 if worried
me="$0"
my_name=identl2cf.sh
I=identl2cf
dryrun="no"
stdparts="${HOME}/mlspgs/l2/l2cf/lib"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -I )
	    stdparts="$2"
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

template="$1"
if [ "$debug" = 1 ]
then
  echo "template $template"
  echo "stdparts $stdparts"
fi
if [ ! -f "$template" ]
then
  echo "Sorry--the l2cf template $template not found"
  exit 1
elif [ ! -d "$stdparts" ]
then
  echo "Sorry--the standard parts directory $stdparts not found"
  exit 1
fi

# Extract names of files to be included
incfiles=`grep '\!include' $template | sed -n 's/\(.*\)\!include(\(.*\.l2cf\))\(?*\)/\2/p' | sed 's/}.*//'`

# Append any known to be missed by this syntax rule
# ---------------------------------------------------------------
# We can find some of these missed that are themselves !include-d
# by their siblings in the standard parts directory by doing this:
# grep '\!include' ${HOME}/mlspgs/l2/l2cf/lib/*.l2cf | \
#   sed -n 's/\(.*\)\!include(\(.*\.l2cf\))\(?*\)/\2/p' | \
#   sed 's/}.*//' | sort | uniq
files="filllesserstate.l2cf writeapriori.l2cf constructstandardgrids.l2cf copystdprods.l2cf $files"
# ---------------------------------------------------------------
files="$incfiles $files"

if [ "$debug" = 1 ]
then
  echo "files $files"
fi

cd $stdparts

for i in $files
do
  # echo $i
  extant_files $i
  # echo $extant_files_result
  if [ "$extant_files_result" != "" ]
  then
    # echo grep -i '\$Id' $i
    grep -i '\$Id' $i
  fi
done
exit
# $Log$
# Revision 1.2  2009/11/06 18:17:16  pwagner
# Fixed bug in statement defining incfiles
#
# Revision 1.1  2009/06/16 22:35:53  pwagner
# First commit
#
