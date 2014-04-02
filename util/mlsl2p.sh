#!/bin/sh
# mlsl2p.sh
# runs the master task for mlsl2
#
# Assumes that:

# (1) It has been called from PGS_PC_Shell.sh as the (pge) in
#  PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v
# (2) $PGE_BINARY_DIR contains mlsl2 and $PGE_SCRIPT_DIR contains
#       slavetmplt.sh and jobstat-sips.sh
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
# Copyright 2005, by the California Institute of Technology. ALL
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

SAVEJOBSTATS="yes"
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
verbose=`echo "$otheropts" | grep -i verbose`

# Check that assumptions are valid
if [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo 'PGS_PC_INFO_FILE undefined'
  echo 'usage:'
  echo 'PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v'
  exit 1
elif [ "$PVM_HOSTS_INFO" = "" ]
then
  echo 'PVM_HOSTS_INFO undefined'
  echo 'It should be the path and name of the host file'
  echo 'a text file containing the hosts available'
  echo 'for running the slave tasks, one host per line'
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
  echo "$PGE_BINARY doesn't exist!"
  exit 1
elif [ ! -r "$PGE_SCRIPT_DIR/slavetmplt.sh"  ]
then
  echo "slavetmplt.sh not in $PGE_SCRIPT_DIR"
  exit 1
fi

H5REPACK=$PGE_BINARY_DIR/h5repack
NETCDFAUGMENT=$PGE_BINARY_DIR/aug_hdfeos5
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$MLSTOOLS/H5REPACK
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_hdfeos5
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_eos5
fi
masterlog="${JOBDIR}/exec_log/process.stdout"
if [ "$MASTERLOG" != "" ]
then
  masterlog="$MASTERLOG"
fi

# We print the file license.txt to stdout
if [ -f $PGE_BINARY_DIR/license.txt ]
then
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

# Use sed to convert slavetmplt.sh into an executable script
# The resulting script sets some toolkit-savvy environment variables
# and then launches the regular mlsl2 binary when summoned to do so.
SLV_SUF=slave
slave_script=$JOBDIR/mlsl2.$SLV_SUF
rm -f $slave_script
sed "s=ppccff=${PGS_PC_INFO_FILE}=;
s=ppggssmmeemmuusseesshhmm=${PGSMEM_USESHM}=;
s=ssllaavveessccrriipptt=${slave_script}=;
s=ootthheerrooppttss=\"${otheropts}\"=;
s=ppggeebbiinnaarryy=${PGE_BINARY}=;
s=ppggssbbiinn=${PGSBIN}=;
s=jjoobbddiirr=${JOBDIR}=;
s=ppggeerroott=${PGE_ROOT}=;" $PGE_SCRIPT_DIR/slavetmplt.sh > $slave_script
chmod a+x $slave_script

NORMAL_STATUS=2

env
#env | sort > ${JOBDIR}/mlsl2p.env
ulimit -s unlimited
ulimit -a
#ulimit -a >> ${JOBDIR}/mlsl2p.env
# First a pre-flight run to check paths
# If a problem, disclosed by exit status, exit before starting big run
if [ "$CHECKPATHS" = "yes" ]
then
  $PGE_BINARY --checkPaths --tk $otheropts
  return_status=`expr $?`
  if [ $return_status != $NORMAL_STATUS ]
  then
     echo "Preflight checkPaths run ended badly"
     echo "Possibly an error in pathnames; please check your PCF"
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

# HOSTNAME will be stored in metadata as ProductionLocation
# to trace provenance
if [ "$HOSTNAME" = "" ]
then
  echo "Need HOSTNAME to trace provenance as ProductionLocation"
  exit 1
fi

# Now we launch the master task itself to set everything in motion
$PGE_BINARY --pge $slave_script --tk --master $PVM_HOSTS_INFO \
  --idents "$IDENTFILE" --loc "$HOSTNAME" $otheropts

# Save return status
return_status=`expr $?`

echo SAVEJOBSTATS $SAVEJOBSTATS
echo masterlog $masterlog
echo PGE_SCRIPT_DIR/jobstat-sips.sh $PGE_SCRIPT_DIR/jobstat-sips.sh
# Save record of progress thru phases
#l2cf=`grep -i l2cf $masterlog | head -1 | awk '{print $9}'`
l2cf=`grep -i 'Level 2 configuration file' $masterlog | head -1 | awk '{print $13}'`
if [ "$SAVEJOBSTATS" = "yes" -a -x "$PGE_SCRIPT_DIR/jobstat-sips.sh" -a -f "$l2cf" ]
then
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  JOBSTATSFILE="$JOBDIR/phases.stats"
  JOBSOPTS="-S $PGE_BINARY_DIR/mlsqlog-scan-sips.py -t $PGE_SCRIPT_DIR/split_path.sh ${JOBDIR}/pvmlog $l2cf $masterlog"
  if [ "$verbose" != "" ]
  then
    JOBSOPTS="-v $JOBSOPTS"
  fi
  echo "$PGE_SCRIPT_DIR/jobstat-sips.sh $JOBSOPTS"
  $PGE_SCRIPT_DIR/jobstat-sips.sh $JOBSOPTS > "$JOBSTATSFILE"
else
  JOBSTATSFILE="none"
fi

# catenate each slave's log to a log file
LOGFILE="${JOBDIR}/pvmlog/mlsl2.log"
if [ ! -w "$LOGFILE"  ]
then
  echo "catenating slave logs in $LOGFILE"
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
  files=`extant_files *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
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
# Revision 1.27  2013/12/05 00:25:24  pwagner
# Fix latest problems with creating JOBSTATSFILE
#
# Revision 1.26  2013/10/29 16:58:20  pwagner
# May set environment variable PGE_BINARY
#
# Revision 1.25  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.24  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.23  2009/12/10 18:54:21  pwagner
# Must repack before augmenting, not after
#
# Revision 1.22  2009/10/27 21:09:12  pwagner
# Added NETCDFAUGMENTing hdfeos std products
#
# Revision 1.21  2009/05/01 22:10:38  honghanh
# Change the initialization of l2cf back to the old way because environment
# variables defined in .env file is not avaiable to mlsl2p.sh.
#
# Revision 1.20  2009/05/01 20:40:22  honghanh
# Use L2CFVERSION to supply l2cf variable
#
# Revision 1.19  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.18  2008/07/11 20:20:18  pwagner
# Fixed bugs regarding products in outputs subdirectory
#
# Revision 1.17  2008/07/10 17:39:50  pwagner
# Handle case when product files are in outputs subdirectory
#
# Revision 1.16  2007/05/10 23:40:34  pwagner
# Used ulimit to increase tiny stacksize so Intel-built mlsl2 can finish
#
# Revision 1.15  2007/03/23 00:32:00  pwagner
# By default no longer dumps slave stdouts into masters
#
# Revision 1.14  2007/02/09 21:31:19  pwagner
# Sets environmental variable as work-around to Lahey 6.2 bug
#
# Revision 1.13  2006/10/19 18:31:10  pwagner
# Now saves chunk stats at end of output
#
# Revision 1.12  2006/04/03 23:10:01  pwagner
# Fixed serious and various bugs
#
# Revision 1.11  2006/03/23 19:22:42  pwagner
# repack with gzip compression all hdf5/hdfeos5 product files
#
# Revision 1.10  2005/12/22 19:07:58  pwagner
# Imitated l1b repacking to defragment forward model files
#
# Revision 1.9  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.8  2004/03/05 19:08:04  pwagner
# Made --cat the default option
#
# Revision 1.7  2004/01/07 17:26:09  pwagner
# Merged in sips-friendly changes
#
# Revision 1.6  2003/12/11 23:07:57  pwagner
# May check each slaves ident against master to verify pge versions are the same
#
# Revision 1.5  2003/11/15 00:45:08  pwagner
# Uses --checkPaths to perform extra preflight checks
#
# Revision 1.4  2003/10/22 23:00:17  pwagner
# Catenates each slaves log file to mlsl2.log at end
#
# Revision 1.3  2003/10/15 17:01:18  pwagner
# Added slv option
#
# Revision 1.2  2003/09/11 20:15:57  pwagner
# Allow OTHEROPTS environmental variable pass-through to mlsl2
#
# Revision 1.1  2003/08/01 16:46:31  pwagner
# First commit
#
