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

# Now if the tool h5repack in the same directory as MLSPROG
# and if the current working directory houses the product files
# then as a final step repack them

# For this final step to operate we must be able to recognize the product files
# A regular expression rreeggeexx will be sed-ed by the install script
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
   echo $extant_files_result
}

run_prog()
{
  echo $BIN_DIR/$MLSPROG $EXTRA_OPTIONS "$@"
  if [ -f $BIN_DIR/license.txt ]
  then
    cat $BIN_DIR/license.txt
  fi
  if [ "$STDERRFILE" != "" ]
  then
    $BIN_DIR/$MLSPROG $EXTRA_OPTIONS "$@" 2> "$STDERRFILE"
  else
    $BIN_DIR/$MLSPROG $EXTRA_OPTIONS "$@"
  fi
  return_status=`expr $?`
  H5REPACK=$BIN_DIR/h5repack
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

echo "We were called as $0 $@"
NORMAL_STATUS=2

# If the first arguments is "-Ef" then source the file of environment settings"
if [ "$1" = "-Ef" ]
then
  shift
  . $1
  shift
fi
# Use the following line to add extra options to MLSPROG
EXTRA_OPTIONS="$OTHEROPTS mlseexxttrraa"

MLSPROG=mlsxxyyzz
# This directory may be a relative path or an absolute one
MLSBIN=mlsbbiinn
# If relative, it must be relative to the following absolute path
MLSHOME=mlshhoommee

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

# The following environmental variable may already have been set
if [ "$PGSMEM_USESHM" = "" ]
then
  PGSMEM_USESHM=NO
fi
ulimit -s unlimited
ulimit -a
export FLIB_DVT_BUFFER=0
# Does MLSBIN start with "/" or not?
# (in other words is it absolute or relative?)

is_absolute=`echo "$MLSBIN" | grep '^\/'`

if [ "$PGE_BINARY_DIR" != "" ]
then
  BIN_DIR=$PGE_BINARY_DIR
elif [ "$is_absolute" = "" ]
then
  BIN_DIR=$MLSHOME/$MLSBIN
else
  BIN_DIR=$MLSBIN
fi

run_prog

# Last chance to find h5repack
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi

# repack level 2 product files to speed things up
if [ -x "$H5REPACK" ]
then
  files=`extant_files rreeggeexx`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files rreeggeexx`
    fi
  fi
  for file in $files
  do
    if [ -w "$file" ]
    then
      packed="$file".p
      if [ "$GZIPLEVEL" != "" ] 
      then
        filter="-f GZIP=$GZIPLEVEL"
      else
        filter=""
      fi
      echo "Packing $file into $packed"
      echo $H5REPACK -i "$file" -o "$packed" $filter
      $H5REPACK -i "$file" -o "$packed" $filter
      # Here we could insert some check involving h5diff if we were dubious
      mv "$packed" "$file"
    fi
  done
fi

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.10  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.9  2009/04/18 00:51:24  pwagner
# May use environment variable OTHEROPTS in command line options
#
# Revision 1.8  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.7  2009/01/16 01:51:39  pwagner
# Takes -Ef settings.env as optional args
#
# Revision 1.6  2008/09/09 16:55:48  pwagner
# May h5repack even products in outputs subdirectory
#
# Revision 1.5  2006/04/12 22:39:24  pwagner
# Needed to define H5REPACK before its use
#
# Revision 1.4  2006/04/03 23:09:09  pwagner
# Added repacking
#
# Revision 1.3  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.2  2002/02/21 21:57:40  pwagner
# EXTRA_OPTIONS now settable
#
# Revision 1.1  2001/08/07 20:57:47  pwagner
# First commit
#
