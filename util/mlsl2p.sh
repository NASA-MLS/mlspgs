#!/bin/sh
# mlsl2p.sh
# runs the master task for mlsl2
#
# Assumes that:

# (1) It has been called from PGS_PC_Shell.sh as the (pge) in
#  PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v
# (2) (HOME)/pvm3/bin/LINUX/ contains both mlsl2 and slavetmplt.sh
#     (or else define an enviromental variable (PVM_EP) and put them there)
# (3) PVM_HOSTS_INFO is defined as an environment variable
#     It should be the path and name of the host file,
#     a text file containing the hosts available
#     for running the slave tasks, one host per line
#
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# usage: see (1) above

CHECKPATHS="yes"
#           ^^^---- "yes" if extra prefilght non-retrieval run to check paths

otheropts="$OTHEROPTS -g -S'slv,mas,chu,opt1,log,pro,time'"

# Check that assumptions are valid
if [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo 'PGS_PC_INFO_FILE undefined'
  echo 'usage:'
  echo 'PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v'
elif [ "$PVM_HOSTS_INFO" = "" ]
then
  echo 'PVM_HOSTS_INFO undefined'
  echo 'It should be the path and name of the host file'
  echo 'a text file containing the hosts available'
  echo 'for running the slave tasks, one host per line'
fi

# In case you used the ep=(PVM_EP)
# in the host file
if [ "$PVM_EP" = "" ]
then
  PVM_EP=$HOME/pvm3/bin/LINUX
fi

if [ ! -x "$PVM_EP/mlsl2"  ]
then
  echo "mlsl2 not in $PVM_EP"
elif [ ! -r "$PVM_EP/slavetmplt.sh"  ]
then
  echo "slavetmplt.sh not in $PVM_EP"
fi

# Use sed to convert slavetmplt.sh into an executable script
# The resulting script sets some toolkt-savvy environment variables
# and then launches the regular mlsl2 binary when summoned to do so.
SLV_SUF=slave
rm -f $PVM_EP/mlsl2.$SLV_SUF
sed "s=ssllaavvee=$SLV_SUF=; s=ppggssbbiinn=$PGSBIN=; \
s=ppccff=$PGS_PC_INFO_FILE=; s=ppvvmmbbiinn=$PVM_EP=" \
"$PVM_EP/slavetmplt.sh" > $PVM_EP/mlsl2.$SLV_SUF
chmod a+x $PVM_EP/mlsl2.$SLV_SUF

NORMAL_STATUS=2

# First a pre-flight run to check paths
# If a problem, disclosed by exit status, exit before starting big run
if [ "$CHECKPATHS" = "yes" ]
then
  $PVM_EP/mlsl2 --checkPaths --tk $otheropts
  return_status=`expr $?`
  if [ $return_status != $NORMAL_STATUS ]
  then
     echo "Preflight checkPaths run ended badly"
     echo "Possibly an error in pathnames; please check your PCF"
     exit 1
  fi
fi

# Now we launch the master task itself to set everything in motion
#$PVM_EP/mlsl2 --pge mlsl2.$SLV_SUF --tk --master $PVM_HOSTS_INFO -g -S'mas,chu,pro,log,opt1,pcf,time'
$PVM_EP/mlsl2 --pge mlsl2.$SLV_SUF --tk --master $PVM_HOSTS_INFO $otheropts

# Save return status
return_status=`expr $?`

# catenate each slave's log to a log file
LOGFILE="$PVM_EP/mlsl2.log"
if [ ! -w "$LOGFILE"  ]
then
  echo "catenating slave logs in $LOGFILE"
  echo "catenating slave logs in $LOGFILE" > "$LOGFILE".1
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  cat $PVM_EP/*.log >> "$LOGFILE".1
  rm -f $PVM_EP/*.log
  mv "$LOGFILE".1 "$LOGFILE"
fi
if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
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
