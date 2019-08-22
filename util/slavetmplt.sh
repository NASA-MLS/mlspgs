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
# mlsl2p.sh seds this file to replace ssllaavvee, ppggssbbiinn, etc.
# as appropriate

# The resulting script will be called by mlsl2 master task
# It will attempt to set some toolkt-savvy environment variables
# and then launch the regular mlsl2 binary

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

#-------------- mycp ----------------
# read each line of stdin
# and then echo it to stdout
# Why?
# Lines like
#   b=$a
# will be echoed as
#   b=100
# if the env variable a has been set to 100
mycp()
{
  set -o noglob
  while read line; do
    # Does line begin with a '#'?
    acomment=`echo $line | grep '^#'`
    if [ "$acomment" != "" ]
    then
      # Dont do anything special--just a comment
      echo $line
    else
      # May need to evaluate twice if line contains a $variable
      eval echo $line
    fi
  done
  set +o noglob
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
      
#---------------------------- launch_mlsl2_ntk
#
# launch the pge without the toolkit panoply
launch_mlsl2_ntk()
{
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
}

#---------------------------- launch_mlsl2_tk
#
# launch the pge assuming you have the toolkit panoply
launch_mlsl2_tk()
{
echo $PGE_BINARY --tk -m --slave $masterTid $otheropts 2>&1 >> "$LOGFILE"
if [ "$PGS_PC_INFO_FILE" = "" ]
then
  echo "PCF undefined" 2>&1 | tee -a "$LOGFILE"
  exit 1
fi
if [ "$runinbackground" != "yes" ]
then
  echo "Must run $PGE_BINARY in foreground" >> "$LOGFILE"
  # Run pge in foreground
  if [ "$CAPTURE_MT" = "yes" ]
  then
    /usr/bin/time -f 'M: %M t: %e' \
      $PGE_BINARY --tk -m --slave $masterTid $otheropts \
       1>> "$LOGFILE" 2> "$STDERRFILE"
  elif [ "$STDERRFILE" != "" ]
  then
    $PGE_BINARY --tk -m --slave $masterTid $otheropts 1>> "$LOGFILE" 2> "$STDERRFILE"
  else
    $PGE_BINARY --tk -m --slave $masterTid $otheropts 2>&1 >> "$LOGFILE"
  fi
  echo "Returned from $PGE_BINARY with status $?" "$LOGFILE"
  Exit_With_Status 0
fi

# Run pge in background
echo "Must run $PGE_BINARY in background" >> "$LOGFILE"
NOTEFILE=`echo "$LOGFILE" | sed 's/.log$/.note/'`
notdone="true"

if [ "$CAPTURE_MT" = "yes" ]
then
  echo   /usr/bin/time -f 'M: %M t: %e' \
    $PGE_BINARY --tk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts  \
    "$LOGFILE" "$STDERRFILE" >> "$LOGFILE"
  /usr/bin/time -f 'M: %M t: %e' \
    $PGE_BINARY --tk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts  \
    1>> "$LOGFILE" 2> "$STDERRFILE" &
elif [ "$STDERRFILE" != "" ]
then
  echo "Writing standard error to $STDERRFILE" >> "$LOGFILE"
  $PGE_BINARY --tk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts \
    1>> "$LOGFILE" 2> "$STDERRFILE" &
else
  $PGE_BINARY --tk -m --slave $masterTid --pidf "$NOTEFILE" $otheropts  \
    >> "$LOGFILE" &
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
  # SETREAD=`which set_read_env.sh`
  . $JOBDIR/job.env
  # . $SETREAD < $JOBDIR/job.env
elif [ -r "$PGE_ROOT/science_env.sh"  ]
then
  . ${PGE_ROOT}/science_env.sh
elif [ -r "$PGSBIN/pgs-env.ksh" ]
then
  # Oops, this stomps on any PCF we might have selected
  # so save it to be restored
  PCF=$PGS_PC_INFO_FILE
  . $PGSBIN/pgs-env.ksh
  export PGS_PC_INFO_FILE=$PCF
fi
PGS_PC_INFO_FILE=ppccff
SLVPROG=ssllaavveessccrriipptt
OTHEROPTS=ootthheerrooppttss
PGE_BINARY=ppggeebbiinnaarryy
UsingPCF=UUssiinnggPPCCFF
export PGS_PC_INFO_FILE PGSMEM_USESHM PGSHOME
export FLIB_DVT_BUFFER=0

# Puts each slave's output into its own unique file
temp_file_name=`get_unique_name log -reverse`
#OLDLOGFILE=$LOGFILE
LOGFILE="${JOBDIR}/pvmlog/$temp_file_name"
UNBUFFERED="${LOGFILE}.u"
pid="$$"
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
runinbackground="no"
CAPTURE_STDERR="yes"
otheropts="-g --uid $pid"
switches="--stdout out -S'slv,opt,log,pro,time,glob'"
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
       # Are we setting to capture stderr?
       a=`grep '^stderr=true' $OPTSFILE`
       if [ "$a" != "" ]
       then
         CAPTURE_STDERR="yes"
       fi
       # Are we pre-processing the opts using environment variables?
       a=`grep '^USEOPTSENV' $OPTSFILE`
       echo "a is $a" >> "$LOGFILE"
       ls $JOBDIR/job.env  >> "$LOGFILE"
       if [ "$a" != "" ]
       then
         mv $OPTSFILE $OPTSFILE.1
         mycp < $OPTSFILE.1 > $OPTSFILE
         echo "Planning to use these opts" >> "$LOGFILE"
         mycp < $OPTSFILE.1 >> "$LOGFILE"
         rm $OPTSFILE.1
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

otheropts="$otheropts $switches"

if [ "$STDERRFILE" = "" ]
then
  if [ "$CAPTURE_MT" = "yes" -o "$CAPTURE_STDERR" = "yes"  ]
  then
    STDERRFILE="$LOGFILE.stderr"
  fi
fi

echo "PGS_PC_INFO_FILE: $PGS_PC_INFO_FILE" 2>&1 | tee -a "$LOGFILE"
echo "masterTid: $masterTid" 2>&1 | tee -a "$LOGFILE"
echo "masterIdent file: $masterIdent" 2>&1 | tee -a "$LOGFILE"
echo "executable: $PGE_BINARY" 2>&1 | tee -a "$LOGFILE"
echo "otheropts $otheropts" 2>&1 | tee -a "$LOGFILE"
echo "STDERRFILE $STDERRFILE" 2>&1 | tee -a "$LOGFILE"

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
if [ "$UsingPCF" != "" ]
then
  launch_mlsl2_tk
else
  l2cf=$1
  launch_mlsl2_ntk
fi
pgepid=$!
echo "pge pid: $pgepid"  >> "$LOGFILE"
sleep $pgestartdelay
echo "Find the pge's pid; call it pgepid2" >> "$LOGFILE"
ps aux | grep "uid"  >> "$LOGFILE"
ps aux | grep "uid $pid " | egrep -v '(grep|time)' >> "$LOGFILE"
pgepid2=`ps aux | grep "uid $pid " | egrep -v '(grep|/time)' | awk '{print $2}'`
if [ "$pgepid" != "$pgepid2" -a "$CAPTURE_MT" != "yes" ]
then
  echo "Warning-- $pgepid and $pgepid2 differ" >> "$LOGFILE"
  ps -l -p "$pgepid" >> "$LOGFILE"
  ps -l -p "$pgepid2" >> "$LOGFILE"
  kill -9 $pgepid
  kill -9 $pgepid2
  Exit_With_Status 1
else
  echo "pgepid and pgepid2 are both $pgepid2" >> "$LOGFILE"
fi
pgepid=$pgepid2
# Did the launch fail immediately?
if [ "$pgepid" = "" ]
then
  echo "Failed to launch $PGE_BINARY in background" >> "$LOGFILE"
  Exit_With_Status 1
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
# Revision 1.47  2019/04/18 16:19:09  pwagner
# May evaluate variables in opts file if USEOPTSENV is set
#
# Revision 1.46  2018/02/28 21:05:30  pwagner
# This version handles both --tk and --ntk cases
#
# Revision 1.45  2018/02/21 21:20:23  pwagner
# Remove unised exits
#
# Revision 1.44  2018/02/09 17:41:10  pwagner
# Correct spelling error
#
# Revision 1.43  2018/01/26 01:29:58  pwagner
# Tried to fix hanging when running in foreground
#
# Revision 1.42  2017/12/22 00:55:42  pwagner
# Correct some bugs in saving stderrfile
#
# Revision 1.41  2017/08/02 22:33:04  pwagner
# Fixed an error in defining pgepid
#
# Revision 1.40  2016/11/16 19:26:50  pwagner
# Avoids stomping on an already-selected PCF
#
# Revision 1.39  2016/10/20 23:26:31  pwagner
# Append chunk stderr to chunk log
#
# Revision 1.38  2016/09/23 00:13:07  pwagner
# Reduce default switches for opt and glob
#
# Revision 1.37  2016/06/04 00:25:12  pwagner
# Tried to prevent another cause of slave deaths
#
# Revision 1.36  2016/05/12 17:01:28  pwagner
# Obey CAPTURE_MT by capturing time, mmory footpint to stderr
#
# Revision 1.35  2015/09/25 00:12:52  pwagner
# Added --maxChunkSize option
#
# Revision 1.34  2014/11/06 01:52:18  pwagner
# Now add_option with auxiliary args pairwise
#
# Revision 1.33  2014/09/29 22:32:00  pwagner
# pgepid obtained from $! mechanism
#
# Revision 1.32  2014/09/12 22:24:48  pwagner
# Enhanced diagnostic logging
#
# Revision 1.31  2014/09/02 17:54:48  pwagner
# Corrected bugs related to running in background
#
# Revision 1.30  2014/08/14 00:49:58  pwagner
# Should correctly run PGE_BINARY as backg task
#
# Revision 1.29  2014/06/25 23:08:27  pwagner
# Handle --set, --setf, --versid
#
# Revision 1.28  2013/11/14 23:57:50  pwagner
# Treats options -D, -V, -R, -S equally
#
# Revision 1.27  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.26  2013/08/30 17:22:46  pwagner
# Trying to stop slaves from writing to toolkit log
#
# Revision 1.25  2013/05/09 18:04:08  pwagner
# Log return status
#
# Revision 1.24  2013/02/14 19:06:36  pwagner
# Defaults to adding instead of skipping unknown args
#
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

