#!/bin/sh
# mlsnrtp.sh
# wrapper script that combines the tasks of running two kinds of pges:
# (1) mls level 1
# (2) master task for mlsl2
# It does so in a pvm-mediated parallel environment, communicating
# with the l2q queue manager
#
# Assumes that:

# (0) It has been called as
#  (pge) 0111 (PCF_file) 25 -v
# (1) $LEVEL1_BINARY_DIR contains mlsl0sn, mlsl1log, and mlsl1g
# (2) $LEVEL2_BINARY_DIR contains mlsl2
# (3) JOBDIR is defined as an environment variable
#     It should be the path where the job is run
# (4) PGE_ROOT is defined as an environment variable
#     It should be the path where the pgs-env.ksh script is kept
# (5) OTHEROPTS is defined as an environment variable
#     It would contain other meaningful runtimeoptions, 
#       e.g. OTHEROPTS="--skipRetrieval"
# (6) l2q and pvm are both runnung
# (7) $MLSTOOLS is defined and ontains the following binary executables
#     and scripts
#     Spartacus
#     ronin.sh
#     set_read_env.sh
# (2) $PGE_SCRIPT_DIR contains
#       slavetmplt.sh and jobstat-sips.sh
#     (or else define an enviromental variable (PVM_EP) and put them there)
# (3) PVM_HOSTS_INFO is defined as an environment variable
#     It should be the path and name of the host file,
#     a text file containing the hosts available
#     for running the slave tasks, one host per line
# (4) The toolkit environment is correctly set and that PGS_PC_Shell.sh
#     is in your PATH and can be used to run each binary executable correctly
#
# Copyright 2012, by the California Institute of Technology. ALL
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
#
GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

# In addition to whatever options and switches may be set by the environment
# variable OTHEROPTS, the following are set:
# -g       trace path of execution through code sections
# --wall   show timing in wall clock times
# chu      show chunk divisions
# opt1     show command line options
# log      copy any log file messages to stdout
# pro      announce input files at opening, output files at creation
# time     summarize time consumed by each code  section, phase, etc.
#EXTRA_OPTIONS="$@"
echo "Launching mlsnrtp with args $@"
PCF=$2
otheropts="$OTHEROPTS --sharedPCF -g --wall --submit l2q --delay 20000 -S'l2q,glob,mas,chu,opt1,log,pro,time'"

# Define the tools we will need
SPARTACUS=$MLSTOOLS/Spartacus
RONIN=$MLSTOOLS/ronin.sh
SETREADENV=$MLSTOOLS/set_read_env.sh
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
if [ ! -x "$L2GPDUMP" ]
then
  L2GPDUMP=$MLSTOOLS/l2gpdump
fi

# Last chance to find h5repack
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi
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
  echo 'It should be the path where the pgs-env.ksh script is kept'
  exit 1
elif [ "$MLSTOOLS" = "" ]
then
  echo 'MLSTOOLS undefined'
  echo 'It should be the path where the other stuff is kept'
  exit 1
elif [ ! -x "$SPARTACUS" ]
then
  echo 'Spartacus not found'
  echo "It should be in $MLSTOOLS"
  exit 1
elif [ ! -x "$RONIN" ]
then
  echo 'ronin.sh not found'
  echo "It should be in $MLSTOOLS"
  exit 1
elif [ ! -x "$SETREADENV" ]
then
  echo 'set_read_env.sh not found'
  echo "It should be in $MLSTOOLS"
  exit 1
fi

# In case you used the ep=(LEVEL2_BINARY_DIR)
# in the host file
if [ "$LEVEL2_BINARY_DIR" = "" ]
then
  LEVEL2_BINARY_DIR=$HOME/pvm3/bin/LINUX
fi
LEVEL2_BINARY=$LEVEL2_BINARY_DIR/mlsl2p.sh

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

# Unlike most mls wrapper scripts, ours must use 0 for normal status
# because Spartacus will never return "2"
NORMAL_STATUS=0

# We print the file license.txt to stdout
if [ -f $LEVEL2_BINARY_DIR/license.txt ]
then
  cat $LEVEL2_BINARY_DIR/license.txt
fi

#env
ulimit -s unlimited
ulimit -a


GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MLSPROG_0=mlsl0sn
#MLSPROG_0=showme.sh
MLSPROG_1=mlsl1log
MLSPROG_2=mlsl1g
MLSPROG_3=mlsl1t

# Here's how we'll use the l2q queue manager:

# We'll start by creating an environment script
# so that each job can inherit our settings

JOBENV=job.env
echo "#!/bin/sh" > $JOBENV
# echo "export PGS_PC_INFO_FILE=$PCF" >> $JOBENV
if [ -f "$MLSTOOLS/tkreset.sh" ]
then
  echo ". $MLSTOOLS/tkreset.sh" >> $JOBENV
fi
echo ". $PGE_ROOT/pgs-env.ksh" >> $JOBENV
echo "export JOBDIR=$JOBDIR" >> $JOBENV
echo "export MLSTOOLS=$MLSTOOLS" >> $JOBENV
echo "export PGSMEM_USESHM=$PGSMEM_USESHM" >> $JOBENV
echo "export FLIB_DVT_BUFFER=$FLIB_DVT_BUFFER" >> $JOBENV
echo "export PGE_BINARY_DIR=$PGE_BINARY_DIR" >> $JOBENV
echo "export PGE_ROOT=$PGE_ROOT" >> $JOBENV
echo "export PVM_HOSTS_INFO=$PVM_HOSTS_INFO" >> $JOBENV
#echo "export OTHEROPTS=$otheropts" >> $JOBENV

# For each of the level 1 jobs, we'll have Spartacus request a host
# from l2q
# For level 2, we'll imitate the mlsl2p.sh script

# (1) Create 3 level 1 jobs and run them
JOB1SCRIPT="job1script.sh"
echo "#!/bin/sh" > $JOB1SCRIPT
echo "$RONIN `pwd` . $JOBDIR/$JOBENV; env; PGS_PC_Shell.sh $LEVEL1_BINARY_DIR/$MLSPROG_0 $EXTRA_OPTIONS $@" >> $JOB1SCRIPT
chmod a+x $JOB1SCRIPT
echo $SPARTACUS $JOB1SCRIPT
$SPARTACUS $JOB1SCRIPT

return_status_0=`expr $?`

if [ $return_status_0 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_0=0
fi

JOB2SCRIPT="job2script.sh"
echo "#!/bin/sh" > $JOB2SCRIPT
echo "$RONIN `pwd` . $JOBDIR/$JOBENV; PGS_PC_Shell.sh $LEVEL1_BINARY_DIR/$MLSPROG_1 $EXTRA_OPTIONS $@" >> $JOB2SCRIPT
chmod a+x $JOB2SCRIPT
echo $SPARTACUS $JOB2SCRIPT
$SPARTACUS $JOB2SCRIPT

return_status_1=`expr $?`

if [ $return_status_1 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_1=0
fi

JOB3SCRIPT="job3script.sh"
echo "#!/bin/sh" > $JOB3SCRIPT
echo "$RONIN `pwd` . $JOBDIR/$JOBENV; PGS_PC_Shell.sh $LEVEL1_BINARY_DIR/$MLSPROG_2 $EXTRA_OPTIONS $@" >> $JOB3SCRIPT
chmod a+x $JOB3SCRIPT
echo $SPARTACUS $JOB3SCRIPT
$SPARTACUS $JOB3SCRIPT

return_status_2=`expr $?`

if [ $return_status_2 != $NORMAL_STATUS ]
then
   exit 1
else
   return_status_2=0
fi

# (2) repack level 1 files to speed things up
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

# (3) Create the level 2 job and run it
# We try to reuse the already existing mlsl2p.sh script, presumably in 
#    $LEVEL2_BINARY_DIR
#. $PGE_ROOT/pgs-env.ksh
#. $JOBDIR/$JOBENV
export OTHEROPTS="$otheropts"
MASTERSCRIPT="masterscript.sh"
echo "#!/bin/sh" > $MASTERSCRIPT
echo ". $JOBDIR/$JOBENV" >> $MASTERSCRIPT
#echo "export PGS_PC_INFO_FILE=$PCF" >> $MASTERSCRIPT
echo "env" >> $MASTERSCRIPT
echo "$LEVEL2_BINARY" >> $MASTERSCRIPT
export JOBDIR="`pwd`"
/bin/rm -fr job1logs
mv pvmlog job1logs
# Now we launch the master task itself to set everything in motion
chmod a+x $MASTERSCRIPT
MASTERSCRIPT="$LEVEL2_BINARY"
echo PGS_PC_Shell.sh $MASTERSCRIPT $@
which PGS_PC_Shell.sh
PGS_PC_Shell.sh $MASTERSCRIPT $@

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

# $Log$
# Revision 1.6  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.5  2013/07/03 17:53:39  pwagner
# LOCOUNT and HICOUNT set acceptable range for count of profiles to catch bogus geodetic angles
#
# Revision 1.4  2012/08/10 20:10:42  pwagner
# Some changes to accommodate goldbrick
#
# Revision 1.3  2012/07/10 22:31:48  pwagner
# Changes needed to work at the sips
#
# Revision 1.2  2012/06/07 21:16:00  pwagner
# Many changes; appears to work now
#
# Revision 1.1  2012/03/23 00:55:20  pwagner
# Firsst commit
#
