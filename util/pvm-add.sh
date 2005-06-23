#!/bin/bash

# --------------- pvm-add help
# Script to add nodes to pvm and optionally set them to useable by
# an already-running l2q
# Usage:
# pvm-add.sh [opt] node1 [node2] .. [noden]
#
#    O p t i o n s
# -h[elp]       print brief help message; exit
# -l2q          tell l2q to resurrect nodes, too
#
# Note:
# (1) pvm (and l2q) must be in your PATH and runnable
# (2) passwordless ssh must be enabled
# (3) nodes must be visible and reachable by ssh
# 
# Result:
# pvm will be started on nodes
# Optionally, already-running l2q will be able to grant nodes to master tasks
# --------------- End pvm-add help
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

# ************
# Main Program
# ************
# 
me="$0"
my_name=pvm-add
I=pvm-add
resurrect="no"
list=""
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -l2q )
	    shift
       resurrect="yes"
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

list="$@"
for node
do
pvm <<EOF
add $node
quit
EOF
done

if [ "$resurrect" = "no" ]
then
  exit 0
fi

comma_list=`echo $list | sed 's/ /,/g'`
l2q -c check $comma_list

exit 0
# $Log$
# Revision 1.1  2005/03/26 00:13:25  pwagner
# First commit
#
