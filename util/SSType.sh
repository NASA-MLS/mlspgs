#!/bin/sh
# SSType.sh
# usage: SSType.sh [option_1] [option_2] ..

# Copyright 2010, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.
# "$Id$"

# --------------- SSType.sh help
# returns the type if Intel SSE (streaming microarchitecture instruction set)
# or Advanced Vector Extensions (AVX, also known as Sandy Bridge New Extensions)
# supported by the host where the command is executed
# E.g., on riverrun would return
# SSSE3
#Usage:
#SSType.sh [options]
#
#    O p t i o n s
# -a            show all SSE enabled on your host; e.g., "SSE2 SSE3 SSSE3"; 
#                 by default shows only the highest
# -d dir        use dir to search for tool hasSSE2, etc.
# -p "p1 p2 .." check on programs p1, p2, .. instead of default hasSSE2, ..
# -v            turn on verbose mode
# -h[elp]       print brief help message; exit
#Notes:
#(1) Proper operation assumes that the following programs
#    are in one of: your path, where MLSTOOLS is defined, or set using -d option
#    hasSSE2  hasSSE3  hasSSSE3 .. hasAVX2 (but see -p option)
#(2) The pattern followed by the program names is the program to check on
#    SSxxx is named "hasSSxxx" (but see below)
#(3) We are assuming that these programs were built for the same bitness
#    as the host where the commands are executed; e.g. 64-bit OS
#(4) To add a new SSE or AV, build /users/pwagner/docs/fortran/helloworld.f90
#    enabling the appropriate option to be tested, 
#    and rename the a.out "hasSSxxx" or "hasAVxxx"
# --------------- End SSType.sh help

NORMAL_STATUS=0
# The checking programs has*
# were compiled using -x * with ifort 17.0.0 20160721
# Note: the following are arranged in order of ascending powers
# i.e., chips that support hasAVX2 must also support hasAVX-I, hasAVX, ..
# program    compiler option
# -------        -------
# hasSSE2       -xSSe2
# hasSSE3       -xSSe3
# hasSSSE3      -xSSSe3
# hasSSE4.1     -xSSe4.1
# hasSSE4.2     -xSSe4.2
# hasAVX        -xAVX
# hasAVX-I      -xCORE-AVX-I
# hasAVX2       -xCORE-AVX2
# Use the following line to add extra programs
PROGS="hasSSE2  hasSSE3  hasSSE4.1  hasSSE4.2  hasSSSE3  hasAVX  hasAVX-I hasAVX2"
me="$0"
my_name=SSType.sh
all="no"
verbose="no"
dir="$MLSTOOLS"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;
    -a )
       all="yes"
       shift
	;;
    -v )
       verbose="yes"
       shift
	;;
    -d )
       dir="$2"
       shift
       shift
	;;
    -p )
       PROGS="$2"
       shift
       shift
	;;
    * )
       more_opts="no"
       ;;
    esac
done

if [ "$verbose" = "yes" ]
then
  echo "Will search for progs in mls tools directory $dir"
  echo "programs to check $PROGS"
fi
#Is MLSTOOLS defined and is it a permissible directory?
if [ -d "$dir" ]
then
  cd $dir
fi
response=""

for PROG in $PROGS
do
  if [ ! -e "$PROG" ]
  then
    echo "$PROG is not executable; did you define MLSTOOLS? use -d option?"
  else
    ./$PROG > /dev/null 2>&1
    return_status=`expr $?`
    SSstring=`echo $PROG | sed 's/has//'`
    if [ "$return_status" = "$NORMAL_STATUS" ]
    then
      if [ "$all" = "yes" ]
      then
        response="$response $SSstring"
      else
        response="$SSstring"
      fi
    fi
    if [ "$verbose" = "yes" ]
    then
      echo "$PROG returns $return_status"
    fi
  fi
done
echo $response
# $Log$
# Revision 1.4  2015/01/22 01:11:58  pwagner
# Expanded default PROGS; added new commandline opts
#
# Revision 1.3  2010/03/11 21:56:47  pwagner
# May specify tools directory with -d option
#
# Revision 1.2  2010/03/06 01:15:50  pwagner
# Was not hiding stderr properly
#
# Revision 1.1  2010/03/06 00:59:57  pwagner
# First commit
#
