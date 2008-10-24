#!/bin/sh
# ronin.sh
# runs a task (e.g. mlsl2) in masterless mode
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

# usage:
# ronin.sh dir pge opts
# (or else ronin.sh dir pge -dryrun opts
#  to merely echo the commands that would be executed)
# where we will 
# (1) cd to dir
# (2) . the file ronin.sh to establish environment settings (if avialable)
# (3) launch "pge opts"

# Used by Spartacus to launch mlsl2 in serial mode
# It will attempt to set some toolkt-savvy environment variables
# and then launch the regular mlsl2 binary

#---------------------------- get_unique_name
get_unique_name()
{

      me="$0"
      my_name=unique_name.sh
   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        mayday=`echo "$1" | grep '\-h'`
        if [ "$mayday" != "" ]
        then
          sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
              | sed -n 's/^.//p' | sed '1 d; $ d'
          exit
        fi
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
      fi
   # Was second arg "-reverse?" If so,
   # that's a reverse part oder in name flag, not a redef of pt
   if [ "$pt" = "-reverse" ] 
   then
     reverse="yes"
     pt="."
   else
     reverse="no"
   fi
   # Is $HOST defined?
   if [ "$HOST" != "" ]
   then
      our_host_name="$HOST"
   elif [ "$HOSTNAME" != "" ]
   then
      our_host_name="$HOSTNAME"
   else
      our_host_name="host"
   fi
   # if in form host.moon.planet.star.. extract host
   our_host_name=`echo $our_host_name | sed 's/\./,/g'`
   our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
   if [ "$reverse" = "yes" ]
   then
     echo $our_host_name.$$.$temp
   else
     echo $temp${pt}$our_host_name${pt}$$
   fi
}
      
#---------------------------- log_failed_run
log_failed_run()
{
  FAILLOG=$HOME/failedjobs.log
  if [ ! -f "$FAILLOG" ]
  then
    echo "log of failed jobs" > "$FAILLOG"
  fi
  echo "`date` $rcmd" >> "$FAILLOG"
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

ENV_SCRIPT="./ronin.env"
DRYRUN="no"
rcmd="$@"

if [ ! -d "$1" ]
then
  temp_file_name=`get_unique_name log -reverse`
  JOBDIR=`pwd`
  LOGFILE="${JOBDIR}/$temp_file_name"
  echo "Sorry--can not cd to $1" > $LOGFILE
  log_failed_run
fi
cd $1
shift
PGE=$1
shift

if [ "$1" = "-dryrun" ]
then
  DRYRUN="yes"
  shift
fi
if [ -r "${HOME}/mlspgs/util/tkreset.sh"  ]
then
. ${HOME}/mlspgs/util/tkreset.sh
fi
if [ -r "$PGE_ROOT/science_env.sh"  ]
then
. ${PGE_ROOT}/science_env.sh
elif [ -r "$PGSBIN/pgs-env.ksh" ]
then
. $PGSBIN/pgs-env.ksh
fi
export
PGSMEM_USESHM=NO
JOBDIR=`pwd`
export JOBDIR PGS_PC_INFO_FILE PGSMEM_USESHM
export FLIB_DVT_BUFFER=0
ulimit -s unlimited

if [ -f "$ENV_SCRIPT" ]
then
  . "$ENV_SCRIPT"
fi

# The next choice, in contrast, puts each output into its own unique file
temp_file_name=`get_unique_name log -reverse`

if [ ! -d "${JOBDIR}/pvmlog" ]
then
  mkdir "${JOBDIR}/pvmlog"
fi
LOGFILE="${JOBDIR}/pvmlog/$temp_file_name"

if [ ! -w "$LOGFILE" ]
then
  echo "#$LOGFILE mlsl2.log" > "$LOGFILE"
fi
echo "Masterless task ronin.sh started with arguments: $PGE $@" >> $LOGFILE

if [ "$DRYRUN" = "yes" ]
then
  echo $PGE $@ >> $LOGFILE
else
  eval $PGE $@ >> $LOGFILE
  # Save return status
  return_status=`expr $?`
  if [ "$return_status" = 1 ]
  then
    log_failed_run
  fi
fi
exit 0

# $Log$
# Revision 1.1  2008/04/09 17:13:54  pwagner
# First commit
#
