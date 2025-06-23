#!/bin/sh
# mlsl2pntk.sh
# runs the master task for mlsl2
#
# Assumes that:

# (1) It has been called to work directly w/o the toolkit panoply
#       with a single command-line arg: the l2cf file
# (2) $PGE_BINARY_DIR contains mlsl2 and $PGE_SCRIPT_DIR contains
#       slavetmpltntk.sh and jobstat-sips.sh
#     (or else define an enviromental variable (PVM_EP) and put them there)
# (3) PVM_HOSTS_INFO is defined as an environment variable
#     It should be the path and name of the host file,
#     a text file containing the hosts available
#     for running the slave tasks, one host per line
# (4) JOBDIR is defined as an environment variable
#     It should be the path where the job is run
# (5) PGE_ROOT is defined as an environment variable
#     It should be the path where the science_env.sh script is kept
# (6) OTHEROPTS is defined as an environment variable
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

# Now if the tool h5repack in PGE_BINARY_DIR
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
CHECKPATHS="yes"
#           ^^^---- "yes" if extra preflight non-retrieval run to check paths

CHECKIDENTS="yes"
#            ^^^---- "yes" if each slave to check its ident against master

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MLSL2PLOG="$JOBDIR/mlsl2p.log"
#             ^^^---- "yes" if progress of chunks thru phases recorded

if [ -f "$MLSL2PLOG" ]
then
  mv $MLSL2PLOG $MLSL2PLOG.1
fi

echo "$MLSL2PLOG"
echo "$MLSL2PLOG" > "$MLSL2PLOG"
env | sort >> "$MLSL2PLOG"

SAVEJOBSTATS="no"
#             ^^^---- "yes" if progress of chunks thru phases recorded

# In addition to whatever options and switches may be set by the environment
# variable OTHEROPTS, the following are set:
# -g       trace path of execution through code sections
# --wall   show timing in wall clock times
# slv      insert slave outputs among master's stdout
# mas      show master's pvm activities as they occur
# chu      show chunk divisions
# opt1     show command line options
# log      copy any log file messages to stdout
# pro      announce input files at opening, output files at creation
# time     summarize time consumed by each code  section, phase, etc.
otheropts="$OTHEROPTS -g --wall -S'mas,chu,opt1,log,pro,time'"
l2cf="$1"
echo "We were called as $0 $@"
# Check that assumptions are valid
if [ ! -f "$l2cf" ]
then
  echo 'Usage: mlsl2pntk.sh l2cf_file'
  exit 1
elif [ "$PVM_HOSTS_INFO" = "" ]
then
  echo 'PVM_HOSTS_INFO undefined' >> "$MLSL2PLOG"
  echo 'PVM_HOSTS_INFO undefined'
  echo 'It should be the path and name of the host file'
  echo 'a text file containing the hosts available'
  echo 'for running the slave tasks, one host per line'
  exit 1
elif [ "$JOBDIR" = "" ]
then
  echo 'JOBDIR undefined' >> "$MLSL2PLOG"
  echo 'JOBDIR undefined'
  echo 'It should be the path where the job is run'
  exit 1
elif [ "$PGE_ROOT" = "" ]
then
  echo 'PGE_ROOT undefined' >> "$MLSL2PLOG"
  echo 'PGE_ROOT undefined'
  echo 'It should be the path where the science_env.sh script is kept'
  exit 1
fi

# In case you used the ep=(PGE_BINARY_DIR)
# in the host file
if [ "$PGE_BINARY_DIR" = "" ]
then
  PGE_BINARY_DIR=$HOME/pvm3/bin/LINUX
fi

if [ "$PGE_BINARY" = "" ]
then
  PGE_BINARY=$PGE_BINARY_DIR/mlsl2
fi

if [ ! -x "$PGE_BINARY"  ]
then
  echo "$PGE_BINARY doesn't exist!" >> "$MLSL2PLOG"
  echo "$PGE_BINARY doesn't exist!"
  exit 1
elif [ ! -r "$PGE_SCRIPT_DIR/slavetmpltntk.sh"  ]
then
  echo "slavetmpltntk.sh not in $PGE_SCRIPT_DIR" >> "$MLSL2PLOG"
  echo "slavetmpltntk.sh not in $PGE_SCRIPT_DIR"
  exit 1
fi

H5REPACK=$PGE_BINARY_DIR/h5repack
MISALIGNMENT=$PGE_BINARY_DIR/misalignment
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$MLSTOOLS/h5repack
fi
if [ ! -x "$MISALIGNMENT" ]
then
  MISALIGNMENT=$MLSTOOLS/misalignment
fi
masterlog="${JOBDIR}/exec_log/process.stdout"
if [ "$MASTERLOG" != "" ]
then
  masterlog="$MASTERLOG"
fi

# We print the file license.txt to stdout
if [ -f $PGE_BINARY_DIR/license.txt ]
then
  cat $PGE_BINARY_DIR/license.txt >> "$MLSL2PLOG"
  cat $PGE_BINARY_DIR/license.txt
fi

# The logs will be written as separate files into ${JOBDIR}/pvmlog
# before being catenated at end of run
# Make certain this directory exists and that it is empty
if [ ! -d ${JOBDIR}/pvmlog ]
then
  mkdir ${JOBDIR}/pvmlog
else
  rm -f ${JOBDIR}/pvmlog/*.log
fi

# The following environmental variable may already have been set
if [ "$PGSMEM_USESHM" = "" ]
then
  PGSMEM_USESHM=NO
fi

export FLIB_DVT_BUFFER=0

# Use sed to convert slavetmpltntk.sh into an executable script
# The resulting script sets some toolkit-savvy environment variables
# and then launches the regular mlsl2 binary when summoned to do so.
SLV_SUF=slave
slave_script=$JOBDIR/mlsl2.$SLV_SUF
rm -f $slave_script
sed "s=ppggssmmeemmuusseesshhmm=${PGSMEM_USESHM}=;
s=ssllaavveessccrriipptt=${slave_script}=;
s=ootthheerrooppttss=\"${otheropts}\"=;
s=ppggeebbiinnaarryy=${PGE_BINARY}=;
s=ppggssbbiinn=${PGSBIN}=;
s=jjoobbddiirr=${JOBDIR}=;
s=ppggeerroott=${PGE_ROOT}=;" $PGE_SCRIPT_DIR/slavetmpltntk.sh > $slave_script
chmod a+x $slave_script

NORMAL_STATUS=2

env
env | sort >> "$MLSL2PLOG"
ulimit -s unlimited
ulimit -a
#ulimit -a >> ${JOBDIR}/mlsl2pntk.env
# First a pre-flight run to check paths
# If a problem, disclosed by exit status, exit before starting big run
if [ "$CHECKPATHS" = "yes" ]
then
  $PGE_BINARY --checkPaths --ntk $otheropts $l2cf
  return_status=`expr $?`
  if [ $return_status != $NORMAL_STATUS ]
  then
     echo "Preflight checkPaths run ended badly"
     echo "Possibly an error in pathnames; please check your l2cf"
     exit 1
  fi
fi

# Next, create a fresh ident based on master task's binary executable
if [ "$CHECKIDENTS" = "yes" ]
then
  IDENTFILE="$JOBDIR/master.ident"
  rm -f "$IDENTFILE"
  ident $PGE_BINARY > "$IDENTFILE"
else
  IDENTFILE="none"
fi

# Now we launch the master task itself to set everything in motion
$PGE_BINARY --pge $slave_script --ntk --master $PVM_HOSTS_INFO \
  --idents "$IDENTFILE" $otheropts $l2cf 2> $JOBDIR/master.l2cfname

# Save return status
return_status=`expr $?`

echo SAVEJOBSTATS $SAVEJOBSTATS
echo masterlog $masterlog
echo PGE_SCRIPT_DIR/jobstat-sips.sh $PGE_SCRIPT_DIR/jobstat-sips.sh
# Save record of progress thru phases
if [ "$SAVEJOBSTATS" = "yes" -a -x "$PGE_SCRIPT_DIR/jobstat-sips.sh" ]
then
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  JOBSTATSFILE="$JOBDIR/phases.stats"
  l2cf=`grep -i 'Level 2 configuration file' $masterlog | head -1 | awk '{print $11}'`
  $PGE_SCRIPT_DIR/jobstat-sips.sh -S $PGE_BINARY_DIR/mlsqlog-scan-sips.py \
    -t $PGE_SCRIPT_DIR/split_path.sh ${JOBDIR}/pvmlog "$l2cf" "$masterlog" \
    > "$JOBSTATSFILE"
else
  JOBSTATSFILE="none"
fi

# catenate each slave's log to a log file
LOGFILE="${JOBDIR}/pvmlog/mlsl2.log"
if [ ! -w "$LOGFILE"  ]
then
  echo "catenating slave logs in $LOGFILE"
  echo "catenating slave logs in $LOGFILE" >> "$MLSL2PLOG"
  echo "catenating slave logs in $LOGFILE" > "$LOGFILE".1
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  cat ${JOBDIR}/pvmlog/*.log >> "$LOGFILE".1
  rm -f ${JOBDIR}/pvmlog/*.log
  mv "$LOGFILE".1 "$LOGFILE"
fi

if [ -f "$LOGFILE" -a -f "$JOBSTATSFILE" ]
then
  cat "$LOGFILE" "$JOBSTATSFILE" > "$LOGFILE".1
  mv "$LOGFILE".1 "$LOGFILE"
fi
# Last chance to find h5repack
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi

# repack level 2 product files to speed things up
if [ -x "$H5REPACK" ]
then
  files=`extant_files *L2FWM*.h5 *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5 *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2FWM*.h5 *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5 *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
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
      echo "Packing $file into $packed" >> "$MLSL2PLOG"
      echo $H5REPACK -i "$file" -o "$packed" $filter
      echo $H5REPACK -i "$file" -o "$packed" $filter >> "$MLSL2PLOG"
      $H5REPACK -i "$file" -o "$packed" $filter
      # Here we could insert some check involving h5diff if we were dubious
      mv "$packed" "$file"
    fi
  done
fi

# Check products for misaligned geolocations
# If they are found to be misaligned, set return_status to 99
file=`extant_files *L2GP-DGG_*.he5`
if [ -x "$MISALIGNMENT" -a "$file" != "" ]
then
  a=`$MISALIGNMENT -silent $file`
  if [ "$a" != "" ]
  then
    echo $a
    echo "Misalignment detected"
    echo hide_files *.he5 *.met *.xml *.h5
    hide_files *.he5 *.met *.xml *.h5
    return_status=99
  fi
fi

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.13  2018/02/28 21:04:13  pwagner
# Made more similar to toolkit-dependent mlsl2p.sh
#
# Revision 1.12  2016/01/05 00:45:08  pwagner
# Correctly captures l2cf name to create job stats
#
# Revision 1.11  2015/12/09 01:32:41  pwagner
# Tries harder to find h5repack and misalignment; fix small bug
#
# Revision 1.10  2015/11/02 23:21:43  pwagner
# Now checks for mialigned geolocations
#
# Revision 1.9  2015/10/07 22:56:29  pwagner
# Automtically writes out l2cf name to master.l2cfname
#
# Revision 1.8  2013/10/29 16:58:37  pwagner
# May set environment variable PGE_BINARY
#
# Revision 1.7  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.6  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.5  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.4  2008/07/11 20:20:18  pwagner
# Fixed bugs regarding products in outputs subdirectory
#
# Revision 1.3  2008/07/10 17:40:06  pwagner
# Handle case when product files are in outputs subdirectory
#
# Revision 1.2  2008/06/17 00:04:20  pwagner
# Fixed, small non-consequential bugs
#
# Revision 1.1  2008/04/04 20:51:22  pwagner
# first commit
#
