#!/bin/sh
# slavetmpltntk.sh
# runs a slave task when mlsl2 is in parallel mode w/o toolkit
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

# usage: not called directly, but as mlsl2.slave after being seded

# Used by mlsl2 master task mlsl2p.sh to launch slave tasks w/o toolkit
# environment when running mlsl2 in parallel mode
# mlsl2p.sh seds this file to replace ssllaavvee, ppggssbbiinn, etc.
# as appropriate

# The resulting script will be called by mlsl2 master task
# It will attempt to set some toolkt-savvy environment variables
# and then launch the regular mlsl2 binary

# Note: Why don't you look into combining this with its
# toolkit-aware counterpart slavetmplt.sh so you don't have
# to maintain the two separate files?

#---------------------------- add_option
# Accumulate a list of commandline options one-by-one
# Useful to prevent adding the same option if it's already present
# usage:
# instead of: otheropts="--skipRetrieval $otheropts"
# e.g.: otheropts=`add_option "$otheropts" --skipRetrieval`

# By default the new options get added on at the end, each separated by a space
# To add them on at the front, use prepend_option

# Be especially careful if you plan to supply an auxiliary arg,
# lest it duplicate an existing one
# e.g.: otheropts=`add_option "--maxFailuresPerChunk 1" "$otheropts" `
# e.g.: otheropts=`add_option "$otheropts"  "--maxFailuresPerChunk 1"`
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
        else
          echo "$1 $2"
        fi
      fi
}

# Exit the entire script
# If a separate stderr file exists, append it to our stdout
Exit_With_Status()
{
   echo "Exit with status $1" >> "$LOGFILE"
   echo "LOGFILE $LOGFILE" >> "$LOGFILE"
   echo "STDERRFILE $STDERRFILE" >> "$LOGFILE"
   # How many args?
   if [ $# -lt 1 ]
   then
     exit 0
   fi
   if [ -f "$STDERRFILE" ]
   then
     echo cat "$STDERRFILE" >> "$LOGFILE"
     cat "$STDERRFILE" >> "$LOGFILE"
   fi
   exit $1
}

#---------------------------- prepend_option
# Accumulate a list of commandline options one-by-one
# Useful to prevent adding the same option if it's already present
# usage:
# instead of: otheropts="--skipRetrieval $otheropts"
# e.g.: otheropts=`prepend_option "$otheropts" --skipRetrieval`

# By default the new options get added on at the beginning, each separated by a space
# To add them on at the end, use add_option

# Be especially careful if you plan to supply an auxiliary arg,
# lest it duplicate an existing one
# e.g.: otheropts=`prepend_option "--maxFailuresPerChunk 1" "$otheropts" `
# e.g.: otheropts=`prepend_option "$otheropts"  "--maxFailuresPerChunk 1"`
prepend_option()
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
        else
          echo "$2 $1"
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

#LOGFILE=${HOME}/slave.log
PGE_ROOT=ppggeerroott
#echo $PGE_ROOT > $LOGFILE
#env >> $LOGFILE
PGSBIN=ppggssbbiinn
#echo $PGSBIN >> $LOGFILE
# **************************************************
# These 2 sleep times can be reset here by hand 
# or by one of the .env files below
# pgestartdelay  how long to wait after starting pge to grep its tid
# pgekilldelay   how long to wait per cycle before killing it
#                  Note: this must be < WaitBeforeKillingSlaves in L2Parallel
pgestartdelay=20
pgekilldelay=120
# **************************************************
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

PGSMEM_USESHM=ppggssmmeemmuusseesshhmm
SLVPROG=ssllaavveessccrriipptt
OTHEROPTS=ootthheerrooppttss
PGE_BINARY=ppggeebbiinnaarryy
JOBDIR=jjoobbddiirr
export PGSMEM_USESHM
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
echo "Slave task slavetmpltntk.sh started with arguments: $@" >> $LOGFILE

# It's possible that $1 is the command name--in which case we
# need to do a shift to get the actual args
if [ "$1" = "$SLVPROG" ]
then
  shift
  echo "Need to shift because 1st arg is command name" 2>&1 | tee -a "$LOGFILE"
fi

masterIdent="none"
runinbackground="no"
otheropts="-g"
switches="-S'slv,opt,log,pro,time,glob'"
OPTSFILE="${JOBDIR}/slave.opts"
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
       otheropts=`add_option "$otheropts" --skipDirect`
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
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Adding argument to fill skipped retrievals: $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --stdout )
       otheropts=`add_option "$otheropts" "$1 $UNBUFFERED"`
       # Note that we can't have all the slaves and masters directing
       # unbuffered stdout to $2
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
    --backg* )
       otheropts=`add_option "$otheropts" $1`
       echo "Skipping argument to run slave pge in background" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       runinbackground="yes"
       shift
       ;;
    --tk )
       otheropts=`add_option "$otheropts" --tk`
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
    --maxChun* )
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Adding arguments to set max chunk size: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --maxFailuresPerCh* )
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Adding arguments to set max failures per chunk: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --maxFailuresPerM* )
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Adding arguments to set max failures per machine: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --lac* )
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Adding arguments to reduce output: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --setf* )
       # Read opts from a file
       if [ ! -f "$OPTSFILE" ]
       then
         sed 's/chunk=/#chunk=/; s/submit=/#submit=/' "$2" > "$OPTSFILE"
       fi
       otheropts=`add_option "$otheropts" "$1 $OPTSFILE"`
       echo "Adding arguments to read opts from file: $1 $OPTSFILE" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       # Are we setting to run in background?
       a=`grep '^backg=true' $OPTSFILE`
       if [ "$a" != "" ]
       then
         runinbackground="yes"
       fi
       shift
       shift
       # Are we setting to capture stderr?
       a=`grep '^stderr=true' $OPTSFILE`
       if [ "$a" != "" ]
       then
         CAPTURE_STDERR="yes"
       fi
       shift
       shift
       ;;
    --set* )
       # set one cmdline opt
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Setting one cmdline opt: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --versid* )
       # set current version id
       otheropts=`add_option "$otheropts" "$1 $2"`
       echo "Setting current version id: $1 $2" >> $LOGFILE
       echo "$otheropts" >> $LOGFILE
       shift
       shift
       ;;
    --* )
       otheropts=`add_option "$otheropts" $1`
       echo "Adding unrecognized option: $1" >> $LOGFILE
       shift
       ;;
    -[DRSV]* )
       switches="$switches $1"
       echo "Appending one of the switch-setting arguments: $1" >> $LOGFILE
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
l2cf=$1

otheropts="$otheropts $switches"

if [ "$STDERRFILE" = "" ]
then
  if [ "$CAPTURE_MT" = "yes" -o "$CAPTURE_STDERR" = "yes"  ]
  then
    STDERRFILE="$LOGFILE.stderr"
  fi
fi

echo "(running without toolkit panoply)" 2>&1 | tee -a "$LOGFILE"
echo "masterTid: $masterTid" 2>&1 | tee -a "$LOGFILE"
echo "masterIdent file: $masterIdent" 2>&1 | tee -a "$LOGFILE"
echo "executable: $PGE_BINARY" 2>&1 | tee -a "$LOGFILE"
echo "otheropts $otheropts" 2>&1 | tee -a "$LOGFILE"
echo "l2cf $l2cf" 2>&1 | tee -a "$LOGFILE"

if [ "$masterTid" = "" ]
then
  echo "masterTid undefined" 2>&1 | tee -a "$LOGFILE"
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

echo $PGE_BINARY --ntk -m --slave $masterTid $otheropts $l2cf 2>&1 >> "$LOGFILE"
if [ "$runinbackground" != "yes" ]
then
  echo "Must run $PGE_BINARY in foreground" >> "$LOGFILE"
  # Run pge in foreground
  if [ "$CAPTURE_MT" = "yes" ]
  then
    /echo usr/bin/time -f 'M: %M t: %e' $PGE_BINARY --ntk -m --slave $masterTid $otheropts $l2cf  "$LOGFILE"  "$STDERRFILE" >> "$LOGFILE"
    /usr/bin/time -f 'M: %M t: %e' $PGE_BINARY --ntk -m --slave $masterTid $otheropts $l2cf  1>> "$LOGFILE" 2> "$STDERRFILE"
  elif [ "$STDERRFILE" != "" ]
  then
    $PGE_BINARY --ntk -m --slave $masterTid $otheropts $l2cf 1>> "$LOGFILE" 2> "$STDERRFILE"
  else
    $PGE_BINARY --ntk -m --slave $masterTid $otheropts $l2cf 2>&1 >> "$LOGFILE"
  fi
  Exit_With_Status 0
fi

# Run pge in background
echo "Must run $PGE_BINARY in background" >> "$LOGFILE"
NOTEFILE=`echo "$LOGFILE" | sed 's/.log$/.note/'`
notdone="true"

if [ "$CAPTURE_MT" = "yes" ]
then
  /usr/bin/time -f 'M: %M t: %e' $PGE_BINARY --ntk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts $l2cf  \
    1>> "$LOGFILE" 2> "$STDERRFILE" &
elif [ "$STDERRFILE" != "" ]
then
  $PGE_BINARY --ntk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts $l2cf  \
    1>> "$LOGFILE" 2> "$STDERRFILE" &
else
  $PGE_BINARY --ntk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts $l2cf  \
    >> "$LOGFILE" &
fi
pgepid=$!
echo "pge pid: $pgepid" "$LOGFILE"
sleep $pgestartdelay

# Did the launch fail immediately?
if [ "$pgepid" = "" ]
then
  echo "Failed to launch $PGE_BINARY in background" "$LOGFILE"
  exit 1
fi

# Write this pid to a uniquely-named file
echo "$pgepid" > "$NOTEFILE"
echo "$pgepid echoed to $NOTEFILE" >> "$LOGFILE"
# Now cycle endlessly until we learn the pge we launched is finished
# We do this so that we can
# (1) append the stderr file onto the end of the logfile
# (2) kill -9 in case it was left only "mostly dead"
while [ "$notdone" = "true" ]
do
  sleep $pgekilldelay
  # Either of two signs that we finished:
  # the pge wrote Finished to the note file
  a=`grep Finished "$NOTEFILE"`
  if [ "$a" != "" ]
  then
    notdone=finished
  fi
  # its pid disappears
  a=`ps p $pgepid | grep "$pgepid"`
  if [ "$a" = "" ]
  then
    notdone=disappeared
  fi
done
echo "$notdone; Returned from $PGE_BINARY with status $?" >> "$LOGFILE"
if [ -f "$STDERRFILE" ]
then
  echo cat "$STDERRFILE" >> "$LOGFILE"
  cat "$STDERRFILE" >> "$LOGFILE"
fi
echo "cat $NOTEFILE" >> "$LOGFILE"
cat $NOTEFILE >> "$LOGFILE"
# For good measure, we alo kill the pge's own pid (n case it was left hanging)
echo "killing $pgepid" >> "$LOGFILE"
kill -9 "$pgepid"
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
Exit_With_Status 0

# $Log$
# Revision 1.20  2018/02/09 17:41:28  pwagner
# Correct spelling error
#
# Revision 1.19  2018/01/26 17:43:33  pwagner
# Tried to fix hanging when running in foreground
#
# Revision 1.18  2017/12/22 00:55:42  pwagner
# Correct some bugs in saving stderrfile
#
# Revision 1.17  2016/10/20 23:26:31  pwagner
# Append chunk stderr to chunk log
#
# Revision 1.16  2016/09/23 00:13:07  pwagner
# Reduce default switches for opt and glob
#
# Revision 1.15  2016/05/12 17:02:05  pwagner
# Obey CAPTURE_MT by capturing time, mmory footpint to stderr
#
# Revision 1.14  2015/09/25 00:12:52  pwagner
# Added --maxChunkSize option
#
# Revision 1.13  2014/11/06 01:52:18  pwagner
# Now add_option with auxiliary args pairwise
#
# Revision 1.12  2014/09/29 22:34:29  pwagner
# Should run in background correctly now
#
# Revision 1.11  2014/06/25 23:08:27  pwagner
# Handle --set, --setf, --versid
#
# Revision 1.10  2013/11/14 23:58:15  pwagner
# Treats options -D, -V, -R, -S equally
#
# Revision 1.9  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.8  2013/02/14 19:06:53  pwagner
# Defaults to adding instead of skipping unknown args
#
# Revision 1.7  2012/09/12 19:48:37  pwagner
# Skips '--submit l2q' args
#
# Revision 1.6  2012/02/15 18:13:47  pwagner
# echo name of l2cf to log file
#
# Revision 1.5  2010/04/23 23:19:42  pwagner
# Removed argumentless 'export' which was killing pvmd
#
# Revision 1.4  2009/05/26 20:04:21  pwagner
# Can pass 3 more options from master
#
# Revision 1.3  2008/07/31 23:57:10  pwagner
# Pass --skipDirectWrite option to slave tasks
#
# Revision 1.2  2008/04/22 18:01:02  pwagner
# Removed things causing more harm than good
#
# Revision 1.1  2008/04/04 20:51:22  pwagner
# first commit
#
