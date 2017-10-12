#!/bin/sh
# postl2script.sh
# script to be run after (a) and (b) master tasks
# by mlsnrt-dual-l2.sh
# its name is passed by the POSTL2SCRIPT environment variables
#
# Copyright 2017, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# Now if the tool h5repack in LEVEL1_BINARY_DIR
# and if the current working directory houses the product files,
# then as a final step repack them

# usage: postl2script.sh c_dir a_dir b_dir
#   arg         meaning
#   ---         -------
#  c_dir        where to put combined std prods
#  a_dir        where to find a task's std prods
#  b_dir        where to find b task's std prods

# Notes:
# It is assumed that l2gpcat and l2auxcat are executable and in your path
# Possible alternatives are
# (1) Define MLSTOOLS and put them there
# (2) Define env variables L2GPCAT and l2AUXCAT with their locations
# directory as this script
# 
#------------------------------- executable_or_exit ------------
#
# Check that a named file is executable;
# if not, then print an error message and quit

executable_or_exit()
{
  a=`which $2`
  if [ ! -x "$a"  ]
  then
    echo "Sorry, $1 not executable or not in your PATH"
    exit 1
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

I=postl2script
# Is there an env file we are supposed to source?
ENVFILE=job.env
if [ -f "$ENVFILE" ]
then
 . ./"$ENVFILE"
fi

# Now check on possible candidates for l2gpcat and l2auxcat
l2gpcat=`which l2gpcat 2>/dev/null`
if [ ! -x "$l2gpcat" ]
then
  l2gpcat=$MLSTOOLS/l2gpcat
fi
l2auxcat=`which l2auxcat 2>/dev/null`
if [ ! -x "$l2auxcat" ]
then
  l2auxcat=$MLSTOOLS/l2auxcat
fi
if [ "$L2GPCAT" != "" ]
then
  l2gpcat=$L2GPCAT
fi
if [ "$L2AUXCAT" != "" ]
then
  l2auxcat=$L2AUXCAT
fi

executable_or_exit l2gpcat "$l2gpcat"
executable_or_exit l2auxcat "$l2auxcat"

list='CO H2O HNO3 N2O O3 SO2 Temperature DGG'
echo "Launching postl2script with args $@"
c_dir=$1
a_dir=$2
b_dir=$3
pwd
cwd=`pwd`
for i in $list
do
  echo "cwd/a_dir $cwd/$a_dir"
  #ls $cwd/$a_dir
  #ls $cwd/$a_dir | grep "L2GP-${i}_" | grep -e 'he5$'
  if [ "$#" -lt "4" ]
  then
    echo No PCF supplied, so we will look in $a_dir
    name1=`ls $a_dir | grep "L2GP-${i}_" | grep -e 'he5$'`
    if [ ! -f "$a_dir/$name1" ]
    then
      echo "Sorry: unable to find $i as $name1"
      exit 1
    fi
  else
    echo Hope to find name of this file in PCF $4
    grep "L2GP-${i}_" $4 | awk -F'|' '{print $2}'
    name1=`grep "L2GP-${i}_" $4 | awk -F'|' '{print $2}'`
  fi
  echo "name1 $name1"
  a_file=$a_dir/$name1
  b_file=$b_dir/$name1
  c_file=$c_dir/$name1
  echo "a,b,c_file $a_file $b_file $c_file"
  # Now there are 4 possibilities as to whether a_file and b_file exist
  if [ -f "$a_file" -a -f "$b_file" ]
  then
    # (1) Both exist :: prefer b over a
    echo "$i: both exist"
    echo $l2gpcat -v -nodup -o $c_file $b_file $a_file
    $l2gpcat -nodup -v -o $c_file $b_file $a_file
    cp $b_file.xml $c_dir
  elif [ -f "$a_file" ]
  then
    # (2) Only a_file exists
    echo "$i: a exists"
    cp $a_file* $c_dir
  elif [ -f "$b_file" ]
  then
    # (3) Only b_file exists
    echo "$i: b exists"
    cp $b_file* $c_dir
  else
    # (4) Neither exists
    echo "Sorry unable to find $i in $a_dir $b_dir"
  fi
done
# Now the DGM files
name1=`ls $a_dir | grep "L2AUX-DGM_" | grep -e 'h5$'`
echo "name1 $name1"
a_file=$a_dir/$name1
b_file=$b_dir/$name1
c_file=$c_dir/$name1
echo "a,b,c_file $a_file $b_file $c_file"
a_file=$a_dir/$name1
echo $l2auxcat -nodup -v -o $c_file -g $b_file $b_file $a_file
$l2auxcat -nodup -v -o $c_file -g $b_file $b_file $a_file
cp $b_file.xml $c_dir
# $Log$
# Revision 1.2  2017/07/13 17:42:31  pwagner
# Fixed error when species in b but not in a
#
# Revision 1.1  2017/05/17 22:26:44  pwagner
# First commit
#
