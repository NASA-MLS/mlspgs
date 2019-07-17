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
  if [ "$PGE_BINARY" = "" ]
  then
    PGE_BINARY=$BIN_DIR/$MLSPROG
  fi
  echo $PGE_BINARY $EXTRA_OPTIONS "$@"
  
  if [ ! -x "$PGE_BINARY"  ]
  then
    echo "************************************************************************"
    echo "$PGE_BINARY doesn't exist!"
    echo "Check for typos in its path and program names"
    echo "************************************************************************"
    exit 1
  fi
  if [ -f $BIN_DIR/license.txt ]
  then
    cat $BIN_DIR/license.txt
  fi
  echo $PGE_BINARY $EXTRA_OPTIONS "$@"
  if [ "$CAPTURE_MT" = "yes" -a "$STDERRFILE" = "" ]
  then
    STDERRFILE=$MLSPROG.stderr
  fi
  if [ "$CAPTURE_MT" = "yes" ]
  then
    /usr/bin/time -f 'M: %M t: %e' $PGE_BINARY $EXTRA_OPTIONS "$@" 2> "$STDERRFILE"
  elif [ "$STDERRFILE" != "" ]
  then
    $PGE_BINARY $EXTRA_OPTIONS "$@" 2> "$STDERRFILE"
  else
    $PGE_BINARY $EXTRA_OPTIONS "$@"
  fi
  return_status=`expr $?`
  H5REPACK=$BIN_DIR/h5repack
  MISALIGNMENT=$BIN_DIR/misalignment
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

# We print the file license.txt to stdout
if [ -f $PGE_BINARY_DIR/license.txt ]
then
  cat $PGE_BINARY_DIR/license.txt
fi

# If the first arguments is "-Ef" then source the file of environment settings"
if [ "$1" = "-Ef" ]
then
  shift
  . $1
  shift
elif [ -r "$JOBDIR/job.env"  ]
then
  . $JOBDIR/job.env
elif [ -r "$PGE_ROOT/science_env.sh"  ]
then
  . ${PGE_ROOT}/science_env.sh
elif [ -r "$PGSBIN/pgs-env.ksh" ]
then
  # Oops, this stomps on any PCF we might have selected
  # so save it to be restored
  PCF=$PGS_PC_INFO_FILE
  . $PGSBIN/pgs-env.ksh
  export PGS_PC_INFO_FILE=$PCF
fi
# Use the following line to add extra options to MLSPROG
EXTRA_OPTIONS="$OTHEROPTS mlseexxttrraa"

MLSPROG=mlsxxyyzz

# So far only mlsl2 can accept cmdline options
# If 
# (1) we are invoked as mlsl2; and if 
# (2) HOSTNAME is defined 
# it will be stored in metadata as ProductionLocation
if [ "$MLSPROG" = "mlsl2"  -a "$HOSTNAME" != "" ]
then
  EXTRA_OPTIONS="--loc $HOSTNAME $EXTRA_OPTIONS"
fi
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

run_prog $@

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

# Check products for misaligned geolocations
# If they are found to be misaligned, set return_status to 99
if [ -x "$MISALIGNMENT" ]
then
  files=`extant_files *L2GP-DGG_*.he5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2GP-DGG_*.he5`
    fi
  fi
  a=""
  if [ "$files" != "" ]
  then
    echo $MISALIGNMENT -silent *L2GP-DGG_*.he5
    $MISALIGNMENT -silent *L2GP-DGG_*.he5
    a=`$MISALIGNMENT -silent *L2GP-DGG_*.he5`
  fi
  if [ "$a" != "" ]
  then
    echo $a
    return_status=99
  fi
fi

echo "return_status $return_status"

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.16  2019/01/31 19:17:11  pwagner
# Pass HOSTNAME on commandline to mlsl2
#
# Revision 1.15  2017/02/09 23:33:17  pwagner
# Avoid stomping on any already-selected PCF
#
# Revision 1.14  2016/05/17 17:07:14  pwagner
# 'dot' job.env if found
#
# Revision 1.13  2016/05/12 17:00:14  pwagner
# Obey CAPTURE_MT by capturing time, mmory footpint to stderr
#
# Revision 1.12  2015/11/02 23:20:47  pwagner
# Now checks for mialigned geolocations
#
# Revision 1.11  2015/10/07 22:58:42  pwagner
# Automtically stores stderr to STDERRFILE; housekeeping
#
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
