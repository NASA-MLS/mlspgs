#!/bin/sh
# mlsl2p.sh
# runs the master task for mlsl2
#
# Assumes that:

# (1) It has been called from PGS_PC_Shell.sh as the (pge) in
#  PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v
# (2) $PGE_BINARY_DIR contains mlsl2 and $PGE_SCRIPT_DIR slavetmplt.sh
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
# Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# usage: see (1) above

CHECKPATHS="yes"
#           ^^^---- "yes" if extra preflight non-retrieval run to check paths

CHECKIDENTS="yes"
#            ^^^---- "yes" if each slave to check its ident against master

otheropts="$OTHEROPTS -g --wall -S'slv,mas,chu,opt1,log,pro,time'"

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
PGE_BINARY=$PGE_BINARY_DIR/mlsl2

if [ ! -x "$PGE_BINARY"  ]
then
  echo "$PGE_BINARY doesn't exist!"
  exit 1
elif [ ! -r "$PGE_SCRIPT_DIR/slavetmplt.sh"  ]
then
  echo "slavetmplt.sh not in $PGE_SCRIPT_DIR"
  exit 1
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
s=jjoobbddiirr=${JOBDIR}=;
s=ppggeerroott=${PGE_ROOT}=;" $PGE_SCRIPT_DIR/slavetmplt.sh > $slave_script
chmod a+x $slave_script

NORMAL_STATUS=2

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

# Now we launch the master task itself to set everything in motion
$PGE_BINARY --pge $slave_script --tk --master $PVM_HOSTS_INFO \
  --idents "$IDENTFILE" $otheropts

# Save return status
return_status=`expr $?`

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

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
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
