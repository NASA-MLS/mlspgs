#!/bin/sh
# mlsnrt.sh
# wrapper script that combines the tasks of running two kinds of pges:
# (1) mls level 1
# (2) master task for mlsl2
#
# Assumes that:

# (0) It has been called from PGS_PC_Shell.sh as the (pge) in
#  PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v
# (1) $LEVEL1_BINARY_DIR contains mlsl0sn, mlsl1log, and mlsl1g
# (2) $LEVEL2_BINARY_DIR contains mlsl2
# (3) JOBDIR is defined as an environment variable
#     It should be the path where the job is run
# (4) PGE_ROOT is defined as an environment variable
#     It should be the path where the science_env.sh script is kept
# (5) OTHEROPTS is defined as an environment variable
#     It would contain other meaningful runtimeoptions, 
#       e.g. OTHEROPTS="--skipRetrieval"
#
# Copyright 2008, by the California Institute of Technology. ALL
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

# usage: see (1) above
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

#------------------------------- hide_files ------------
#
# Hide files when something goes awry
# usage: hide_files arg1 [arg2] ..

hide_files()
{
   hide_files_result=
   # Trivial case ($# = 0)
   if [ "$1" != "" ]
   then
      for file
      do
         if [ -f "$file" ]
         then
               hide_files_result="$hide_files_result $file"
         fi
      done
   fi
   echo $hide_files_result
   if [ ! -d hidden ]
   then
     mkdir hidden
   fi
   mv $hide_files_result hidden
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)


REECHO="`echo $0 | sed 's/mlsnrt/reecho/'`"
# In addition to whatever options and switches may be set by the environment
# variable OTHEROPTS, the following are set:
# -g       trace path of execution through code sections
# --wall   show timing in wall clock times
# chu      show chunk divisions
# opt1     show command line options
# log      copy any log file messages to stdout
# pro      announce input files at opening, output files at creation
# time     summarize time consumed by each code  section, phase, etc.
otheropts="$OTHEROPTS -g --wall -S'chu,opt1,log,pro,time'"

# Check that assumptions are valid
if [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo 'PGS_PC_INFO_FILE undefined'
  echo 'usage:'
  echo 'PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v'
  exit 1
elif [ "$JOBDIR" = "" ]
then
  echo 'JOBDIR undefined'
  echo 'It should be the path where the job is run'
  exit 1
elif [ "$PGE_ROOT" = "" ]
then
  echo 'PGE_ROOT undefined'
  echo 'It should be the path where the science_env.sh script is kept'
  exit 1
fi

# In case you used the ep=(LEVEL2_BINARY_DIR)
# in the host file
if [ "$LEVEL2_BINARY_DIR" = "" ]
then
  LEVEL2_BINARY_DIR=$HOME/pvm3/bin/LINUX
fi
LEVEL2_BINARY=$LEVEL2_BINARY_DIR/mlsl2

if [ ! -x "$LEVEL2_BINARY"  ]
then
  echo "$LEVEL2_BINARY doesn't exist!"
  exit 1
fi

masterlog="${JOBDIR}/exec_log/process.stdout"
if [ "$MASTERLOG" != "" ]
then
  masterlog="$MASTERLOG"
fi

# The following environmental variable may already have been set
if [ "$PGSMEM_USESHM" = "" ]
then
  PGSMEM_USESHM=NO
fi
if [ "$LOCOUNT" = "" ]
then
  LOCOUNT=0
fi
if [ "$HICOUNT" = "" ]
then
  HICOUNT=4000
fi

export FLIB_DVT_BUFFER=0

NORMAL_STATUS=2

# We print the file license.txt to stdout
if [ -f $LEVEL2_BINARY_DIR/license.txt ]
then
  cat $LEVEL2_BINARY_DIR/license.txt
fi

env
ulimit -s unlimited
ulimit -a


GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MLSPROG_0=mlsl0sn
MLSPROG_1=mlsl1log
MLSPROG_2=mlsl1g
MLSPROG_3=mlsl1t
echo $LEVEL1_BINARY_DIR/$MLSPROG_0 $EXTRA_OPTIONS "$@"
$LEVEL1_BINARY_DIR/$MLSPROG_0 $EXTRA_OPTIONS "$@"

return_status_0=`expr $?`

if [ $return_status_0 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_0=0
fi

echo $LEVEL1_BINARY_DIR/$MLSPROG_1 $EXTRA_OPTIONS "$@"
$LEVEL1_BINARY_DIR/$MLSPROG_1 $EXTRA_OPTIONS "$@"
return_status_1=`expr $?`

if [ $return_status_1 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_1=0
fi

H5REPACK=$LEVEL1_BINARY_DIR/h5repack
NETCDFAUGMENT=$LEVEL1_BINARY_DIR/aug_hdfeos5
L2GPDUMP=$LEVEL1_BINARY_DIR/l2gpdump
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$MLSTOOLS/H5REPACK
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_hdfeos5
fi
# Last chance to find h5repack
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi

# We are going to insist that both H5REPACK and NETCDFAUGMENT
# be defined before going any further. To override this, set
# the environment variable OKTONOTAUGMENT to "yes"
if [ "$OKTONOTAUGMENT" = "" ]
then
  if [ ! -x "$NETCDFAUGMENT" ]
  then
    echo "NETCDFAUGMENT not defined"
    exit 1
  elif [ ! -x "$H5REPACK" ]
  then
    echo "H5REPACK not defined"
    exit 1
  fi
fi

if [ ! -x "$L2GPDUMP" ]
then
  L2GPDUMP=$MLSTOOLS/l2gpdump
fi

echo $LEVEL1_BINARY_DIR/$MLSPROG_2 $EXTRA_OPTIONS "$@"
$LEVEL1_BINARY_DIR/$MLSPROG_2 $EXTRA_OPTIONS "$@"

return_status_2=`expr $?`

if [ $return_status_2 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_2=0
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

# Now we launch the master task itself to set everything in motion
echo $LEVEL2_BINARY --sharedPCF --tk $otheropts
$LEVEL2_BINARY --sharedPCF --tk $otheropts

# Save return status
return_status=`expr $?`

# repack level 2 product files to speed things up
if [ -x "$H5REPACK" ]
then
  files=`echo *L2FWM*.h5 *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5 *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`$REECHO *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
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

# augment level 2 product files to make them netcdf-compatible
# (must not reverse order of repack, augment)
if [ -x "$NETCDFAUGMENT" ]
then
  files=`$REECHO *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`$REECHO *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
    fi
  fi
  for file in $files
  do
    if [ -w "$file" ]
    then
      echo $NETCDFAUGMENT $file
      $NETCDFAUGMENT $file
    fi
  done
fi

echo "Check that the number of profiles is within range"
#pwd
#ls
if [ -d outputs ]
then
#  ls outputs/*
  cd outputs
else
  echo "Separate outputs directory not found"
fi
#pwd
#ls
$L2GPDUMP -status *L2GP-Temper*.he5
files=`extant_files *L2GP-Temp*.he5`
echo "$files"
if [ -f "$files" ]
then
  echo "Checking $L2GPDUMP -status $files"
  count=`$L2GPDUMP -status "$files" \
    | grep 'valid data co' | awk '{print $4}'`
  if [ "$count" -lt "$LOCOUNT" -o "$count" -gt "$HICOUNT" ]
  then
    echo "Too few or too many profiles; number was $count"
    echo hide_files *.he5 *.met *.xml *.h5
    hide_files *.he5 *.met *.xml *.h5
    exit 1
  fi
fi

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.9  2016/04/06 21:10:24  pwagner
# Check on return status immediately on return from mlsl1log
#
# Revision 1.8  2013/11/23 00:59:41  pwagner
# Hide product files if number of profiles too many or too few
#
# Revision 1.7  2013/09/06 16:39:15  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.6  2013/07/02 22:49:09  pwagner
# LOCOUNT and HICOUNT set acceptable range for count of profiles to catch bogus geodetic angles
#
# Revision 1.5  2012/08/10 20:09:26  pwagner
# Some changes to accommodate goldbrick
#
# Revision 1.4  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.3  2010/01/28 01:12:10  pwagner
# May augment level 2 product files to be netcdf-compatible
#
# Revision 1.2  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.1  2008/01/25 18:24:55  pwagner
# First commit
#
