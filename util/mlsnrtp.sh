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
BREAKER_PY=checknrtgranuleforbreaks.py
BREAKER_SAV=checknrtgranuleforbreaks.sav
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

if [ "$JOBDIR" = "" ]
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
# But don't clobber if it exists already
JOBENV=job.env
if [ ! -f "$JOBENV" ]
then
  echo "Must build $JOBENV"
  echo "#!/bin/sh" > $JOBENV
  echo "export JOBDIR=$JOBDIR" >> $JOBENV
  echo "export MLSTOOLS=$MLSTOOLS" >> $JOBENV
  echo "export PGE_BINARY_DIR=$PGE_BINARY_DIR" >> $JOBENV
  echo "export PGE_ROOT=$PGE_ROOT" >> $JOBENV
  echo "export PVM_HOSTS_INFO=$PVM_HOSTS_INFO" >> $JOBENV
else
  echo "Will reuse $JOBENV"
fi
if [ -f "$MLSTOOLS/tkreset.sh" ]
then
  echo ". $MLSTOOLS/tkreset.sh" >> $JOBENV
fi
echo ". $PGE_ROOT/pgs-env.ksh" >> $JOBENV
echo "export PGSMEM_USESHM=$PGSMEM_USESHM" >> $JOBENV
echo "export FLIB_DVT_BUFFER=$FLIB_DVT_BUFFER" >> $JOBENV
echo "export PGS_PC_INFO_FILE=$2" >> $JOBENV
echo "export UsingPCF=yes" >> $JOBENV

# For the level 1 jobs, we'll have Spartacus request a host
# from l2q
# For level 2, we'll imitate the mlsl2p.sh script

. $JOBDIR/$JOBENV

# Do we have an outputs subdirectory of JOBDIR for std prods?
if [ -d "$JOBDIR/outputs" ]
then
  STDPRODDIR="$JOBDIR/outputs"
else
  STDPRODDIR="$JOBDIR"
fi

if [ "$MUSTHAVEBREAKER" = "yes" ]
then
  cp $MLSTOOLS/$BREAKER_PY $MLSTOOLS/$BREAKER_SAV $STDPRODDIR
  echo cp $MLSTOOLS/$BREAKER_PY $MLSTOOLS/$BREAKER_SAV $STDPRODDIR
  if [ ! -x "$STDPRODDIR/$BREAKER_PY" ]
  then
    echo "$BREAKER_PY not found"
    echo "It should have been in $MLSTOOLS"
    exit 1
  fi
fi

# (1) Create a single level 1 job and run it
JOB1SCRIPT="job1script.sh"
JOB1STDERR="job1script.stderr"
echo "#!/bin/sh" > $JOB1SCRIPT
if [ "$CAPTURE_MT" = "yes" ]
then
  echo "$RONIN `pwd` . $JOBDIR/$JOBENV; env; /usr/bin/time -f 'M: %M t: %e' PGS_PC_Shell.sh $LEVEL1_BINARY_DIR/mlsl1.sh $EXTRA_OPTIONS $@ 2> $JOB1STDERR" >> $JOB1SCRIPT
else
  echo "$RONIN `pwd` . $JOBDIR/$JOBENV; env; PGS_PC_Shell.sh $LEVEL1_BINARY_DIR/mlsl1.sh $EXTRA_OPTIONS $@" >> $JOB1SCRIPT
fi
chmod a+x $JOB1SCRIPT
echo $SPARTACUS $JOB1SCRIPT
$SPARTACUS $JOB1SCRIPT

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
cd $STDPRODDIR
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
  # Also check for chunk breaks, i.e. anomalously large
  # derivatives [ds / dt] where 
  # ds   distance between successive profiles
  # dt   time between successive profiles
  if [ "$MUSTHAVEBREAKER" = "yes" ]
  then
    echo "Check for anomalously large gaps between profiles"
    file=*L2GP-Temper*.he5
    # Unfortunately, BREAKER_PY expects the file to have a string like
    # _yyyydDoythhmm.he5 in its name. If it dosn't the script treats
    # the name as if it were a directory and gets sore when it
    # turns out not to be a directory.
    # So we'll check if the file name instead looks like _yyyydDoy.he5
    # and if it does, we'll copy it a file with a more complaisant name
    file2=`echo $file | sed '/_20[0-9][0-9]d[0-9][0-9][0-9]\.he5/ s/\.he5/t0000.he5/'`
    if [ "$file2" != "$file" ]
    then
      cp $file $file2
    fi
    echo ./$BREAKER_PY --verbose $file2
    echo ./$BREAKER_PY --verbose $file2 > $BREAKER_PY.out
    ./$BREAKER_PY --verbose $file2 >> $BREAKER_PY.out
    return_status=`expr $?`
    cat $BREAKER_PY.out
    if [ "$return_status" != 0 ]
    then
      echo "Break detected in $file, status $return_status"
      hide_files *.he5 *.met *.xml *.h5
      exit 1
    fi
  fi
fi

# $Log$
# Revision 1.16  2018/03/01 00:22:38  pwagner
# Now sets UsingPCF for use with new mlsl2p.sh
#
# Revision 1.15  2017/05/19 20:48:10  pwagner
# repaired another error
#
# Revision 1.14  2017/05/19 20:31:11  pwagner
# Should now work both at scf and at sips
#
# Revision 1.13  2017/02/16 23:52:57  pwagner
# Added gap detector
#
# Revision 1.12  2016/05/19 19:48:46  pwagner
# Hope we fixed the plast fix
#
# Revision 1.11  2016/05/19 17:32:49  pwagner
# Restore PCF file name for master level 2
#
# Revision 1.10  2016/05/13 21:58:42  pwagner
# Obey CAPTURE_MT by capturing time, mmory footpint to stderr
#
# Revision 1.9  2016/05/13 00:36:53  pwagner
# Dont clobber an existing JOBENV
#
# Revision 1.8  2016/02/12 20:14:38  pwagner
# Use mlsl1.sh wrapper script to stop with error status when appropriate
#
# Revision 1.7  2013/11/23 00:59:48  pwagner
# Hide product files if number of profiles too many or too few
#
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
