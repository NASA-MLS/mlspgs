#!/bin/sh
# slavetmplt.sh
# runs a slave task when mlsl2 is in parallel mode
#
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# usage: not called directly, but as mlsl2.slave after being seded

# Used by mlsl2 master task mlsl2p.sh to launch slave tasks in toolkit
# environment when running mlsl2 in parallel mode
# mlsl2p.sh sed's this file to replace ssllaavvee, ppggssbbiinn, etc.
# as appropriate

# The resulting script will be called by mlsl2 master task
# It will attempt to set some toolkt-savvy environment variables
# and then launch the regular mlsl2 binary

#---------------------------- do_the_call
#
# Put after a trip through the main program simply to
# separate the various options which may have been conglomerated
# before being passed as args to this script

do_the_call()
{
# Do we need MLSHOME?
# MLSHOME=mlshhoommee
# MLSUTIL=$MLSHOME/util
# PVM_BIN=$HOME/pvm3/bin/LINUX
PVM_BIN=ppvvmmbbiinn
PGSMEM_USESHM=NO
PGS_PC_INFO_FILE=ppccff

# The following choice puts the outputs of all the slaves into a single file
LOGFILE="$PVM_BIN/mlsl2.log"

# The next choice, in contrast, puts each slave's output into its own unique file
#temp_file_name=`unique_name.sh mlsl2`
#LOGFILE="$PVM_BIN/$temp_file_name"

if [ ! -w "$LOGFILE" ]
then
  echo "#mlsl2.log" > "$LOGFILE"
fi
#echo "Called as $@" >> "$LOGFILE"
#echo "command is $0" >> "$LOGFILE"
#echo "1st arg is $1" >> "$LOGFILE"
#echo "2nd arg is $2" >> "$LOGFILE"
#echo "3rd arg is $3" >> "$LOGFILE"

SLVPROG=mlsl2.ssllaavvee
# It's possible that $1 is the command name--in which case we
# need to do a shift to get the actual args
if [ "$1" = "$SLVPROG" ]
then
  shift
  echo "Need to shift because 1st arg is command name" >> "$LOGFILE"
fi

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in
    --chunk )
       shift
       shift
       ;;
    --slave )
       masterTid="$2"
       shift
       shift
       ;;
    --* )
       shift
       ;;
    -* )
       shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

echo "PGS_PC_INFO_FILE: $PGS_PC_INFO_FILE" >> "$LOGFILE"
echo "masterTid: $masterTid" >> "$LOGFILE"
echo "executable: $PVM_BIN/mlsl2" >> "$LOGFILE"

export PGSMEM_USESHM
export PGS_PC_INFO_FILE

if [ "$masterTid" = "" ]
then
  echo "masterTid undefined" >> "$LOGFILE"
  exit 1
elif [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo "PCF undefined" >> "$LOGFILE"
  exit 1
elif [ ! -x "$PVM_BIN/mlsl2" ]
then
  echo "$PVM_BIN/mlsl2 not found/executable" >> "$LOGFILE"
  exit 1
fi

$PVM_BIN/mlsl2 --tk -m --slave $masterTid -g -S'slv,opt1,log,pro,time' >> "$LOGFILE"

}
      
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

PGSBIN=ppggssbbiinn
. $PGSBIN/pgs-env.ksh

all_my_opts=$@
do_the_call $all_my_opts
exit 0

# $Log$
# Revision 1.1  2003/08/01 16:46:31  pwagner
# First commit
#
