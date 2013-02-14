#!/bin/sh
# slavetmplt.sh
# runs a slave task when mlsl2 is in parallel mode
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

# usage: not called directly, but as mlsl2.slave after being seded

# Used by mlsl2 master task mlsl2p.sh to launch slave tasks in toolkit
# environment when running mlsl2 in parallel mode
# mlsl2p.sh sed's this file to replace ssllaavvee, ppggssbbiinn, etc.
# as appropriate

# The resulting script will be called by mlsl2 master task
# It will attempt to set some toolkt-savvy environment variables
# and then launch the regular mlsl2 binary

#---------------------------- add_option
# To prevent adding the same option if it's already part
# usage:
# instead of: otheropts="--skipRetrieval $otheropts"
# try: otheropts=`add_option "$otheropts" --skipRetrieval`
# By default the new options get added on at the end, each separated by a space
# To add them on at the front, supply a third arg
# e.g.: otheropts=`add_option "$otheropts" --skipRetrieval to_front`
add_option()
{
   # How many args?
      if [ $# -gt 1 ]
      then
      # Does the option to be added duplicate one already present?
        mayday=`echo "$1" | grep -e "$2"`
        if [  "$1" = "" ]
        then
          echo $2
        elif [  "$mayday" != "" ]
        then
          echo $1
        elif [ $# -gt 2 ]
        then
          echo "$2 $1"
        else
          echo "$1 $2"
        fi
      fi
}

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
      
#---------------------------- do_the_call
#
# Put after a trip through the main program simply to
# separate the various options which may have been conglomerated
# before being passed as args to this script.
# It's possible that "--slave {task ID}" can come through as one
# argument containing a space, rather than as 2 arguments.
# The same thing will happen with --idents {masterIdent}.

do_the_call()
{

# Not sure why this isn't done automatically
if [ -r "$HOME/.bashrc" ]
then
. $HOME/.bashrc
fi

JOBDIR=jjoobbddiirr
PGSMEM_USESHM=ppggssmmeemmuusseesshhmm
#LOGFILE=${HOME}/slave.log
PGE_ROOT=ppggeerroott
#echo $PGE_ROOT > $LOGFILE
#env >> $LOGFILE
PGSBIN=ppggssbbiinn
#echo $PGSBIN >> $LOGFILE
if [ -r "$JOBDIR/job.env"  ]
then
  . $JOBDIR/job.env
elif [ -r "$PGE_ROOT/science_env.sh"  ]
then
  . ${PGE_ROOT}/science_env.sh
elif [ -r "$PGSBIN/pgs-env.ksh" ]
then
  . $PGSBIN/pgs-env.ksh
fi
PGS_PC_INFO_FILE=ppccff
SLVPROG=ssllaavveessccrriipptt
OTHEROPTS=ootthheerrooppttss
PGE_BINARY=ppggeebbiinnaarryy

export PGS_PC_INFO_FILE PGSMEM_USESHM
export FLIB_DVT_BUFFER=0

# Puts each slave's output into its own unique file
temp_file_name=`get_unique_name log -reverse`
#OLDLOGFILE=$LOGFILE
LOGFILE="${JOBDIR}/pvmlog/$temp_file_name"
UNBUFFERED="${LOGFILE}.u"
#ENVSETTINGS="${LOGFILE}.env"
#echo $LOGFILE >> $OLDLOGFILE

#env  2>&1 | tee "$ENVSETTINGS"

if [ ! -w "$LOGFILE" ]
then
  echo "#$LOGFILE mlsl2.log" 2>&1 | tee "$LOGFILE"
fi
echo "Slave task slavetmplt.sh started with arguments: $@" >> $LOGFILE

# It's possible that $1 is the command name--in which case we
# need to do a shift to get the actual args
if [ "$1" = "$SLVPROG" ]
then
  shift
  echo "Need to shift because 1st arg is command name" 2>&1 | tee -a "$LOGFILE"
fi

masterIdent="none"
otheropts="-g"
switches="-S'slv,opt1,log,pro,time,glob1'"
# otheropts="-g -S'slv,opt1,log,pro,time,glob1'"
# otheropts="$OTHEROPTS"
echo "otheropts starting as $otheropts" 2>&1 | tee -a "$LOGFILE"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    echo "Considering argument: $1" >> $LOGFILE
    case "$1" in
    --chunk )
       echo "Skipping chunk-setting arguments: $1 $2" >> $LOGFILE
       shift
       shift
       ;;
    --crash* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding argument to crash on any error: $1" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --delay )
       echo "Skipping argument to change master delay: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --idents )
       masterIdent="$2"
       shift
       shift
       ;;
    --skipDir* )
       otheropts=`add_option "$otheropts" --skipRetrieval`
       echo "Adding argument to skip direct writes: $1" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --skipR* )
       otheropts=`add_option "$otheropts" --skipRetrieval`
       echo "Adding argument to skip retrieval: $1" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --skipS* )
       otheropts=`add_option "$otheropts" --skipSections`
       otheropts=`add_option "$otheropts" $2`
       echo "Adding argument to skip certain sections: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --slave* )
	   masterTid="$2"
	   shift
	   shift
       ;;
    --state* )
       otheropts=`add_option "$otheropts" $1`
       otheropts=`add_option "$otheropts" $2`
       echo "Adding argument to fill skipped retrievals: $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --stdout )
       otheropts=`add_option "$otheropts" $1`
       # Note that we can't have all the slaves and masters directing
       # unbuffered stdout to $2
       otheropts=`add_option "$otheropts" $UNBUFFERED`
       echo "Adding argument to buffer stdout to $UNBUFFERED" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --share* )
       otheropts=`add_option "$otheropts" --sharedPCF`
       shift
       ;;
    --submit )
       echo "Skipping argument to submit to a special l2 queue manager: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --tk )
       otheropts=`add_option "$otheropts" --tk`
       shift
       ;;
    --cat )
       otheropts=`add_option "$otheropts" --cat`
       echo "Adding argument to catenate split dgg/dgm files: $1" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --wall )
       otheropts=`add_option "$otheropts" --wall`
       echo "Adding argument to use wall clock: $1" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --patch* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding argument to patch existing directWrite files" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       ;;
    --maxFailuresPerCh* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding argument to set max falures per chunk" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       otheropts=`add_option "$otheropts" $1`
       shift
       ;;
    --maxFailuresPerM* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding argument to set max falures per machine" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       otheropts=`add_option "$otheropts" $1`
       shift
       ;;
    --lac* )
       otheropts=`add_option "$otheropts" $1`
       otheropts=`add_option "$otheropts" $2`
       echo "Adding arguments to reduce output: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding unrecognized option: $1" >> $LOGFILE
       shift
       ;;
#   --* )
#      echo "Skipping unknown argument: $1" >> $LOGFILE
#      shift
#      ;;
#    -S* )
#       echo "Skipping switch-setting argument: $1" >> $LOGFILE
#       shift
#       ;;
    -S* )
       switches="$switches $1"
       echo "Appending switch-setting argument: $1" >> $LOGFILE
       shift
       ;;
    -f* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding argument to trace forward model: $1" >> $LOGFILE
       shift
       ;;
    -g* )
       echo "Skipping superfluous argument to trace execution: $1" >> $LOGFILE
       shift
       ;;
    -* )
       echo "Skipping unknown argument: $1" >> $LOGFILE
       shift
       ;;
    * )
       echo "Finished considering arguments" >> $LOGFILE
       more_opts="no"
       ;;
    esac
done

otheropts="$otheropts $switches"

echo "PGS_PC_INFO_FILE: $PGS_PC_INFO_FILE" 2>&1 | tee -a "$LOGFILE"
echo "masterTid: $masterTid" 2>&1 | tee -a "$LOGFILE"
echo "masterIdent file: $masterIdent" 2>&1 | tee -a "$LOGFILE"
echo "executable: $PGE_BINARY" 2>&1 | tee -a "$LOGFILE"
echo "otheropts $otheropts" 2>&1 | tee -a "$LOGFILE"

if [ "$masterTid" = "" ]
then
  echo "masterTid undefined" 2>&1 | tee -a "$LOGFILE"
  exit 1
elif [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo "PCF undefined" 2>&1 | tee -a "$LOGFILE"
  exit 1
elif [ ! -x "$PGE_BINARY" ]
then
  echo "$PGE_BINARY not found/executable" 2>&1 | tee -a "$LOGFILE"
  exit 1
fi

# Compare master task's ident with this slave's
# If different, quit with error message
if [ -f "$masterIdent" ]
then
  the_diff=`ident "$PGE_BINARY" | diff - "$masterIdent"`
  if [ ! "$the_diff" = "" ]
  then
     echo "master task ident differs from slave"
     echo "Possibly an error in paths; please check PVM_BIN, PVM_EP, ~/bin"
     exit 1
  fi
fi

#env | sort > "$ENVSETTINGS"
ulimit -s unlimited
#ulimit -a >> "$ENVSETTINGS"

echo $PGE_BINARY --tk -m --slave $masterTid $otheropts 2>&1 | tee -a "$LOGFILE"
$PGE_BINARY --tk -m --slave $masterTid $otheropts 2>&1 | tee -a "$LOGFILE"

}
      
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

# This just separates any args that have suffered conglomeration
# e.g., "arg1 arg2" being passed as a single space-containing arg
all_my_opts=$@
do_the_call $all_my_opts
exit 0

# $Log$
# Revision 1.23  2012/08/10 20:08:41  pwagner
# Some changes to accommodate goldbrick
#
# Revision 1.22  2012/07/02 23:07:15  pwagner
# Fixed long-standing bug
#
# Revision 1.21  2010/04/23 23:18:42  pwagner
# Removed argumentless 'export' which was killing pvmd
#
# Revision 1.20  2009/05/26 20:04:21  pwagner
# Can pass 3 more options from master
#
# Revision 1.19  2008/07/31 23:57:10  pwagner
# Pass --skipDirectWrite option to slave tasks
#
# Revision 1.18  2008/04/22 18:00:41  pwagner
# Removed things causing more harm than good
#
# Revision 1.17  2008/01/08 18:44:12  pwagner
# Pass --sharedPCF option to slaves
#
# Revision 1.16  2007/08/31 00:06:43  pwagner
# Passes -f, --crash options to slaves
#
# Revision 1.15  2007/08/17 00:42:02  pwagner
# Passes switches from master task to slave
#
# Revision 1.14  2007/05/10 23:40:12  pwagner
# Used ulimit to increase tiny stacksize so Intel-built mlsl2 can finish
#
# Revision 1.13  2007/02/09 21:31:07  pwagner
# Sets environmental variable as work-around to Lahey 6.2 bug; accepts --delay --stdout opts
#
# Revision 1.12  2006/10/05 23:41:32  pwagner
# Needed for latest options, to work at scf (needs sips testing)
#
# Revision 1.11  2006/04/21 23:58:39  pwagner
# Ugly LASTDITCHPGSBIN set to overcome unknown problem with some scf hosts; remove later
#
# Revision 1.10  2005/06/23 22:20:46  pwagner
# Reworded Copyright statement
#
# Revision 1.9  2004/03/05 19:09:33  pwagner
# Passes --cat option to mlsl2
#
# Revision 1.8  2004/01/07 17:26:09  pwagner
# Merged in sips-friendly changes
#
# Revision 1.7  2003/12/11 23:07:57  pwagner
# May check each slaves ident against master to verify pge versions are the same
#
# Revision 1.6  2003/10/22 23:01:51  pwagner
# Changed each slaves temp log file name to bzzz.host.log
#
# Revision 1.5  2003/10/20 22:22:00  pwagner
# Use tee to let slaves output to stdout--which pvm reechoes to master
#
# Revision 1.4  2003/10/15 20:56:27  pwagner
# Each slave sends its stdout to a unique file
#
# Revision 1.3  2003/09/11 20:17:13  pwagner
# Passes --skipRetr[] option to mlsl2
#
# Revision 1.2  2003/09/05 23:45:51  pwagner
# Made name of log file a variable: LOGFILE; tweaked initial comments
#
# Revision 1.1  2003/08/01 16:46:31  pwagner
# First commit
#

