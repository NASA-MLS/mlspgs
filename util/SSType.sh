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
# supported by the host where the command is executed
# E.g., on riverrun would return
# SSSE3
#Usage:
#SSType.sh [options]
#
#    O p t i o n s
# -a            show all SSE; e.g., "SSE2 SSE3 SSSE3"; 
#                 by default shows only highest
# -h[elp]       print brief help message; exit
#Notes:
#(1) Proper operation assumes that either the following programs
#    are in your path or that MLSTOOLS is defined and that can be found there
#    hasSSE2  hasSSE3  hasSSSE3
#(2) We are assuming that these programs were built for the same bitness
#    as the host where the commands are executed; e.g. 64-bit OS
# --------------- End SSType.sh help


NORMAL_STATUS=0
# The checking programs has*
# were compiled using -x * with ifort v11.1.56
# program    compiler option
# -------        -------
# hasSSE2       -x -xSSe2
# hasSSE3       -x -xSSe3
# hasSSE4.1     -x -xSSe4.1
# hasSSE4.2     -x -xSSe4.2
# hasSSSE3      -x -xSSSe3
# Use the following line to add extra options to MLSPROG
PROGS="hasSSE2  hasSSE3  hasSSSE3"
me="$0"
my_name=SSType.sh
all="no"
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
    * )
       more_opts="no"
       ;;
    esac
done

#Is MLSTOOLS defined and is it a permissible directory?
if [ -d "$MLSTOOLS" ]
then
  cd $MLSTOOLS
fi
response=""

for PROG in $PROGS
do
  $PROG > /dev/null
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
done
echo $response
# $Log$
