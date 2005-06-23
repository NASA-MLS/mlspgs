#!/bin/sh
# mlsxxyyzz.sh
# run a program specified as the variable MLSPROG
# assuming that it's in directory MLSBIN
# Then return an exit status of:
# 1 if program's exit status is different from
# the variable specified as NORMAL_STATUS; otherwise
# 0
#
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# usage: mlsxxyyzz.sh [option_1] [option_2] ..

# Used by mlspgs/MakeFC to install script versions of each mls program
# MLSPROG. The script version returns a status of 0 only if 
# the MLSPROG goes past some "finish" line in the code
# causing it to exit with status=NORMAL_STATUS
# MakeFC sed's this file to replace xxyyzz, hhoommee, etc. as appropriate

NORMAL_STATUS=2
# Use the following line to add extra options to MLSPROG
EXTRA_OPTIONS=mlseexxttrraa

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
   echo $MLSHOME/$MLSBIN/$MLSPROG $EXTRA_OPTIONS "$@"
   $MLSHOME/$MLSBIN/$MLSPROG $EXTRA_OPTIONS "$@"
else
   echo $MLSBIN/$MLSPROG $EXTRA_OPTIONS "$@"
   $MLSBIN/$MLSPROG $EXTRA_OPTIONS "$@"
fi

return_status=`expr $?`

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.2  2002/02/21 21:57:40  pwagner
# EXTRA_OPTIONS now settable
#
# Revision 1.1  2001/08/07 20:57:47  pwagner
# First commit
#
