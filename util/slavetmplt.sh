#!/bin/sh
# slavetmplt.sh
# run a program specified as the variable MLSPROG
# assuming that it's in directory MLSBIN
# Then return an exit status of:
# 1 if program's exit status is different from
# the variable specified as NORMAL_STATUS; otherwise
# 0
#
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# usage: slavetmplt.sh [option_1] [option_2] ..

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

if [ ! -w "$PVM_BIN/mlsl2.log" ]
then
  echo "#mlsl2.log" > "$PVM_BIN/mlsl2.log"
fi
# echo "Called as $@" >> "$PVM_BIN/mlsl2.log"
# echo "command is $0" >> "$PVM_BIN/mlsl2.log"
# echo "1st arg is $1" >> "$PVM_BIN/mlsl2.log"
# echo "2nd arg is $2" >> "$PVM_BIN/mlsl2.log"
# echo "3rd arg is $3" >> "$PVM_BIN/mlsl2.log"

SLVPROG=mlsl2.ssllaavvee
# It's possible that $1 is the command name--in which case we
# need to do a shift to get the actual args
if [ "$1" = "$SLVPROG" ]
then
  shift
  # echo "Need to shift because 1st arg is command name" >> "$PVM_BIN/mlsl2.log"
fi

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in
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

echo "PGS_PC_INFO_FILE: $PGS_PC_INFO_FILE" >> "$PVM_BIN/mlsl2.log"
echo "masterTid: $masterTid" >> "$PVM_BIN/mlsl2.log"
echo "executable: $PVM_BIN/mlsl2" >> "$PVM_BIN/mlsl2.log"

export PGSMEM_USESHM
export PGS_PC_INFO_FILE

if [ "$masterTid" = "" ]
then
  echo "masterTid undefined" >> "$PVM_BIN/mlsl2.log"
  exit 1
elif [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo "PCF undefined" >> "$PVM_BIN/mlsl2.log"
  exit 1
elif [ ! -x "$PVM_BIN/mlsl2" ]
then
  echo "$PVM_BIN/mlsl2 not found/executable" >> "$PVM_BIN/mlsl2.log"
  exit 1
fi

$PVM_BIN/mlsl2 --tk --slave $masterTid -g -S'slv,opt1,log,pro,time'  >> "$PVM_BIN/mlsl2.log"

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
