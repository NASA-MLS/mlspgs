#!/bin/sh
# mlsxxyyzz.sh
# usage: mlsxxyyzz.sh [option_1] [option_2] ..

# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# run separate programs specified as the variables
# MLSPROG_1, MLSPROG_2, MLSPROG_3
# assuming that they're in directory MLSBIN
# Then return an exit status of ( for each program):
# 1 if program's exit status is different from
# the variable specified as NORMAL_STATUS; otherwise
# 0
# Then assemble the sequence of status for the three programs as follows:
# MLSPROG_1  MLSPROG_2  MLSPROG_3
#  (0 or 1)   (0 or 1)   (0 or 1)
# in a binary string of length three; e.g. '001' means
# MLSPROG_1 and _2 had normal status, but _3 did not
# Finally return a final exit status based on this binary string:
Successful_codes="000 001"
# Possible and sensible values for the above would be "000" or "000 001"
# with the following effects
# "000"     ==> returns error-free status only if all three programs do
# "000 001" ==> returns error-free status even if 3rd program fails
# We're assuming MLSPROG_2 only runs if MLSPROG_1 succeeds
# SImilarly, MLSPROG_3 only runs if MLSPROG_2 succeeds
# Used by mlspgs/MakeFC to install script versions of each mls program
# MLSPROG. The script version returns a status of 0 only if 
# the MLSPROG goes past some "finish" line in the code
# causing it to exit with status=NORMAL_STATUS
# MakeFC sed's this file to replace xxyyzz, hhoommee, etc. as appropriate

NORMAL_STATUS=2
# Use the following line to add extra options to MLSPROG
EXTRA_OPTIONS=mlseexxttrraa

MLSPROG_1=mlsxxyyzz_1
MLSPROG_2=mlsxxyyzz_2
MLSPROG_3=mlsxxyyzz_3
# This directory may be a relative path or an absolute one
MLSBIN=mlsbbiinn
# If relative, it must be relative to the following absolute path
MLSHOME=mlshhoommee

# Does MLSBIN start with "/" or not?
# (in other words is it absolute or relative?)

is_absolute=`echo "$MLSBIN" | grep '^\/'`

if [ "$is_absolute" = "" ]
then
   echo $MLSHOME/$MLSBIN/$MLSPROG_1 $EXTRA_OPTIONS "$@"
   $MLSHOME/$MLSBIN/$MLSPROG_1 $EXTRA_OPTIONS "$@"
else
   echo $MLSBIN/$MLSPROG_1 $EXTRA_OPTIONS "$@"
   $MLSBIN/$MLSPROG_1 $EXTRA_OPTIONS "$@"
fi

return_status_1=`expr $?`

if [ $return_status_1 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_1=0
fi

if [ "$is_absolute" = "" ]
then
   echo $MLSHOME/$MLSBIN/$MLSPROG_2 $EXTRA_OPTIONS "$@"
   $MLSHOME/$MLSBIN/$MLSPROG_2 $EXTRA_OPTIONS "$@"
else
   echo $MLSBIN/$MLSPROG_2 $EXTRA_OPTIONS "$@"
   $MLSBIN/$MLSPROG_2 $EXTRA_OPTIONS "$@"
fi

return_status_2=`expr $?`

if [ $return_status_2 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_2=0
fi

if [ "$is_absolute" = "" ]
then
   echo $MLSHOME/$MLSBIN/$MLSPROG_3 $EXTRA_OPTIONS "$@"
   $MLSHOME/$MLSBIN/$MLSPROG_3 $EXTRA_OPTIONS "$@"
else
   echo $MLSBIN/$MLSPROG_3 $EXTRA_OPTIONS "$@"
   $MLSBIN/$MLSPROG_3 $EXTRA_OPTIONS "$@"
fi

return_status_3=`expr $?`

if [ $return_status_3 = $NORMAL_STATUS ]
then
   return_status_3=0
else
   return_status_3=1
fi

binary_string="${return_status_1}${return_status_2}${return_status_3}"

return_status=1
for code in $Successful_codes
do
  if [ $binary_string = $code ]
  then
    return_status=$NORMAL_STATUS
  fi
done

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
