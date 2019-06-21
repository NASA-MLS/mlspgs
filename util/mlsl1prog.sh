#!/bin/sh
# mlsxxyyzz.sh
# usage: mlsxxyyzz.sh [option_1] [option_2] ..

# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# run separate programs specified as the variables
# MLSPROG_0, MLSPROG_1, MLSPROG_2, MLSPROG_3
# assuming that they're in directory MLSBIN
# Then return an exit status of ( for each program):
# 1 if program's exit status is different from
# the variable specified as NORMAL_STATUS; otherwise
# 0
# Now if MLSPROG_0 fails, let's bail out immediately--
# no sense trying with any of the others
# Then assemble the sequence of status for three of the programs as follows:
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
# Similarly, MLSPROG_3 only runs if MLSPROG_2 succeeds
# Used by mlspgs/MakeFC to install script versions of each mls program
# MLSPROG. The script version returns a status of 0 only if 
# the MLSPROG goes past some "finish" line in the code
# causing it to exit with status=NORMAL_STATUS
# Typing 
#    make install
# sed's this file to replace xxyyzz, hhoommee, etc. as appropriate

# Now if the tool h5repack in the same directory as the three level 1 programs
# and if the current working directory houses the l1b files created by
# the level 1 programs, then as a final step repack the l1b files

NORMAL_STATUS=2
# Use the following line to add extra options to MLSPROG
EXTRA_OPTIONS=mlseexxttrraa

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MLSPROG_0=mlsxxyyzz_0
MLSPROG_1=mlsxxyyzz_1
MLSPROG_2=mlsxxyyzz_2
MLSPROG_3=mlsxxyyzz_3
# This directory may be a relative path or an absolute one
MLSBIN=mlsbbiinn
# If relative, it must be relative to the following absolute path
MLSHOME=mlshhoommee

# Did we create an env file to be sourced by this job?
JOBENV=job.env
if [ -f "$JOBENV" ]
then
  . ./$JOBENV
fi

if [ "$PGE_BINARY_DIR" = "" ]
then
  is_absolute=`echo "$MLSBIN" | grep '^\/'`
  if [ "$is_absolute" = "" ]
  then
    PGE_BINARY_DIR=$MLSBIN
  else
    PGE_BINARY_DIR=$MLSHOME/$MLSBIN
  fi
fi

ulimit -s unlimited
ulimit -a

# We print the file license.txt to stdout
echo $PGE_BINARY_DIR/$MLSPROG_0 $EXTRA_OPTIONS "$@"
if [ -f $PGE_BINARY_DIR/license.txt ]
then
  cat $PGE_BINARY_DIR/license.txt
fi
$PGE_BINARY_DIR/$MLSPROG_0 $EXTRA_OPTIONS "$@"

return_status_0=`expr $?`

if [ $return_status_0 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_0=0
fi

echo $PGE_BINARY_DIR/$MLSPROG_1 $EXTRA_OPTIONS "$@"
$PGE_BINARY_DIR/$MLSPROG_1 $EXTRA_OPTIONS "$@"
return_status_1=`expr $?`
H5REPACK=$PGE_BINARY_DIR/h5repack

if [ $return_status_1 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_1=0
fi

echo $PGE_BINARY_DIR/$MLSPROG_2 $EXTRA_OPTIONS "$@"
$PGE_BINARY_DIR/$MLSPROG_2 $EXTRA_OPTIONS "$@"

return_status_2=`expr $?`

if [ $return_status_2 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_2=0
fi

echo $PGE_BINARY_DIR/$MLSPROG_3 $EXTRA_OPTIONS "$@"
$PGE_BINARY_DIR/$MLSPROG_3 $EXTRA_OPTIONS "$@"

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

# Last chance to find h5repack
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi

# repack level 1 files to speed things up
if [ -x "$H5REPACK" ]
then
  files=`echo *L1*.h5`
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
      # Here we could insert some check involving l1bdiff if we were dubious
      mv "$packed" "$file"
    fi
  done
fi

# Exit with status according to whether we succeeded or failed
if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.13  2016/09/07 00:33:30  pwagner
# Should work better when running goldbrick with prebuilt binaries
#
# Revision 1.12  2015/01/16 17:16:17  pwagner
# stack size limit lifted
#
# Revision 1.11  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.10  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.9  2006/04/03 23:11:07  pwagner
# Repack only files with write access
#
# Revision 1.8  2006/04/03 22:14:41  pwagner
# Fixed serious syntax errors, oyther bugs
#
# Revision 1.7  2006/03/27 19:21:15  pwagner
# Added new pge mlsl0sn to store switch settings
#
# Revision 1.6  2006/03/23 19:21:18  pwagner
# Add gzip compression when repacking hdf5 files
#
# Revision 1.5  2005/12/22 19:07:09  pwagner
# Restored idea of h5repacking l1b files
#
# Revision 1.4  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.3  2005/04/19 21:02:38  pwagner
# Previous committal an oversight--undoing it
#
# Revision 1.2  2005/04/18 16:27:18  pwagner
# Mistakenly deallocated timings before possibly needing to use it--fixed
#
# Revision 1.1  2003/02/03 23:57:06  pwagner
# First commit
#
