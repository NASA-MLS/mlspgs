#!/bin/sh
# mlsxxyyzz.sh
# run a program specified as the variable MLSPROG
# assuming that it's in directory MLSBIN
# Then return an exit status
# of 1 if program's exit status is different from
# the variable specified as TARGET_STATUS
# else 0
#
# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# usage: mlsxxyyzz.sh

# Used by mlspgs/MakeFC to install script versions of each mls program
# MLSPROG. The script version returns a status of 0 only if 
# the MLSPROG goes past some "finish" line in the code
# causing it to exit with status=TARGET_STATUS

TARGET_STATUS=2

MLSPROG=mlsxxyyzz
# This directory may be a relative path or an absolute one
MLSBIN=mlsbbiinn
# If relative, it must be relative to the following absolute path
MLSHOME=mlshhoommee

# Does MLSBIN start with "/" or not?
# (in other words is it absolute or relative?)

is_absolute=`echo "$MLSBIN" | grep '^\/'`

if [ "$is_absolute" = "" ]
then
   echo $MLSHOME/$MLSBIN/$MLSPROG "$@"
   $MLSHOME/$MLSBIN/$MLSPROG "$@"
else
   echo $MLSBIN/$MLSPROG "$@"
   $MLSBIN/$MLSPROG "$@"
fi

return_status=`expr $?`

if [ $return_status != $TARGET_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
