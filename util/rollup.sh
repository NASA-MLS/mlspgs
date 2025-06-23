#!/bin/sh
# --------------- rollup.sh help
# rollup.sh
# l2gpcat all the files in dir1 for each day, 
# saving them in dir2 under the same day
#
# rollup.sh dir1 dir2

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -nrt2         Choose files for the set of products in nrt2
#                   (default is to just cat the O3 and Temperature files)
#    -b command    Use command instead of batch_l2gpcat.sh

# Bugs and limitations
# (1) batch_l2gpcat.sh is assumed to exist, to be an executable, and to be in
#     $HOME/mlspgs/util
# (2) mkpath is assumed to exist, to be an executable, and to be in
#     your path

# usage: see (1) above

# --------------- End rollup.sh help
# Copyright 2012, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
#
debug=1
#     ^  -- set this to 1 if worried
me="$0"
my_name=rollup.sh
I=rollup
l2gpcat=$HOME/mlspgs/util/batch_l2gpcat.sh
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
dryrun="no"
nrtversion="nrt"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -nrt2 )
	    nrtversion="nrt2"
	    shift
       ;;
    -b )
	    l2gpcat="$2"
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

dir1="$1"
dir2="$2"

if [ "$debug" = 1 ]
then
  echo "batch_l2gpcat.sh $l2gpcat"
  echo "dir1 $dir1"
  echo "dir2 $dir2"
fi
if [ "$dir1" = "$dir2" ]
then
  echo "Sorry--dir1 and dir2 must be different"
elif [ ! -d "$dir1" ]
then
  echo "Sorry--dir1 must be a directory"
elif [ ! -x "$l2gpcat" ]
then
  echo "Sorry--batch_l2gpcat.sh must exist and be executable"
  exit
fi

days=`/bin/ls $dir1`
for day in $days
do
  if [ ! -d "$dir1/$day" ]
  then
    echo "Sorry--$dir1/$day is not a subdirectory containing l2gp files"
  elif [ "$dryrun" = "yes" ]
  then
    echo $l2gpcat -dryrun $nrtversion - $dir1/$day $dir2/$day
    $l2gpcat -dryrun $nrtversion - $dir1/$day $dir2/$day
  else
    if [ ! -d $dir2/$day ]
    then
      mkpath $dir2/$day
    fi
    $l2gpcat $nrtversion - $dir1/$day $dir2/$day
  fi
done

# $Log$
