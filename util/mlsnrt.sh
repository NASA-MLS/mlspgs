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

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)


REECHO="`echo $0 | sed 's/mlsnrt/reecho/'`"
# In addition to whatever options and switches may be set by the environment
# variable OTHEROPTS, the following are set:
# -g       trace path of execution through code sections
# --wall   show timing in wall clock times
# --cat    have master task catenate any split dgg/dgm files
# chu      show chunk divisions
# opt1     show command line options
# log      copy any log file messages to stdout
# pro      announce input files at opening, output files at creation
# time     summarize time consumed by each code  section, phase, etc.
otheropts="$OTHEROPTS -g --wall --cat -S'chu,opt1,log,pro,time'"

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
H5REPACK=$LEVEL1_BINARY_DIR/h5repack
NETCDFAUGMENT=$LEVEL1_BINARY_DIR/aug_hdfeos5
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$MLSTOOLS/H5REPACK
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_hdfeos5
fi

if [ $return_status_1 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_1=0
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

# Now we launch the master task itself to set everything in motion
echo $LEVEL2_BINARY --sharedPCF --tk $otheropts
$LEVEL2_BINARY --sharedPCF --tk $otheropts

# Save return status
return_status=`expr $?`

# repack level 2 product files to speed things up
if [ -x "$H5REPACK" ]
then
  files=`echo *L2FWM*.h5 *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5 *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
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

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
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
