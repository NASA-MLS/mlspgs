#!/bin/bash
# --------------- sipsl2 help
# Show level 2 jobs running on lightspeed, scramjet
# Usage:
# sipsl2.sh [options]
#    O p t i o n s
# -h[elp]     print brief help message; exit
# -l          show lightspeed jobs only
# -s          show scramjet jobs only
# -c          convert dates from yyyydoy to yyyy Month day
# -d          show directory names
# -x          show chunks, nodes that died
#   The following options carry an automatic restriction to running jobs
# -r          show how many chunks running
# -full       show full stats: running, completed, crashed, on deck, total
# -t          show date, time the run started
#   The following options print one line per machine
# -fail       show machines that jobs could not be spawned on
#
# Note:
# (1) you must be logged into a host that can see /workops/jobs/science
# (2) jobs should still be putting stdout into JOBDIR/exec_log/process.stdout
# (3) dateconverter must be in your path to use the -c option
#     (or else /home/pwagner/bin/dateconverter must be visible)
# 
# Result:
# a table of jobs, machine, date, status, etc.
# --------------- End sipsl2 help
# Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

# "$Id$"

note_failures()
  {
  # How many args?
  if [ $# -lt 2 ]
  then
    echo "note_failures needs exactly two arguments: line_numbers and file_name"
    exit
  elif [ ! -r "$2" ]
  then
    echo "$2 not found"
    exit
  fi
  echo -e "machine \t failure with"
  for line_number in $1
  do
    #echo "line number $line_number of $2"
    #sed -n ''$line_number' p' "$2"
    machine_info=`sed -n ''$line_number' p' "$2" | awk '{print $9, $11}'`
    #echo "machine_info: $machine_info"
    machine=`echo $machine_info | awk '{print $1}'`
    info=`echo $machine_info | awk '{print $2}'`
    case "$info" in
      -6)
        note="$machine \t pvm demon"
        ;;
      -7)
        note="$machine \t file system"
        ;;
      *)
        note="$machine \t ($info: unknown)"
        ;;
    esac
    echo -e $note
  done
}

# ************
# Main Program
# ************
# 
debug="no"
me="$0"
my_name=sipsl2
I=sipsl2
convert="no"
lightspeed="yes"
scramjet="yes"
dirnames="no"
died="no"
restrict="no"
fail="no"
running="no"
full="no"
time="no"
more_opts="yes"
DATECONVERTER=`which dateconverter 2>/dev/null`
if [ ! -r "$DATECONVERTER" ]
then
  DATECONVERTER="/home/pwagner/bin/dateconverter"
fi 
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -c )
	    shift
       convert="yes"
       ;;
    -d )
	    shift
       dirnames="yes"
       ;;
    -x )
	    shift
       died="yes"
       ;;
    -r )
	    shift
       running="yes"
       restrict="yes"
       ;;
    -fail )
	    shift
       fail="yes"
       restrict="yes"
       ;;
    -full )
	    shift
       full="yes"
       restrict="yes"
       ;;
    -t )
	    shift
       time="yes"
       restrict="yes"
       ;;
    -l )
	    shift
       scramjet="no"
       ;;
    -s )
	    shift
       lightspeed="no"
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;

    * )
       more_opts="no"
       ;;
    esac
done

dirsplus=`echo /workops/jobs/science/*/mlsl2.slave`
dirs=`echo $dirsplus | sed "s:mlsl2.slave::g"`
if [ "$debug" = "yes" ]
then
  echo "dirs: $dirs"
fi
if [ "$restrict" = "yes" ]
then
  list='data date \t machine'
else
  list='data date \t machine \t status'
fi
if [ "$dirnames" = "yes" ]
then
  list="$list \t directory"
fi
if [ "$died" = "yes" ]
then
  list="$list \t dead chunks \t killer nodes"
fi
if [ "$running" = "yes" ]
then
  list="$list \t chunks running"
fi
if [ "$full" = "yes" ]
then
  echo "Stats order: running complete crashed waiting total"
  list="$list \t full stats "
fi
if [ "$time" = "yes" ]
then
  list="$list \t run start (date time)"
fi
if [ "$fail" = "no" ]
then
  echo -e $list
fi
for dir in $dirs
do
  # These are pretty quick (at least compared with grepping the whole
  # file as later tests require)
  testl=`head $dir/exec_log/process.stdout | grep -i lightspeed`
  tests=`head $dir/exec_log/process.stdout | grep -i scramjet`
  statbad=`tail $dir/exec_log/process.stdout | grep -i "ended badly"`
  statgood=`tail $dir/exec_log/process.stdout | grep -i "catenating slave"`
  l1boa=`echo $dir/*L1BOA_*`
  date=`echo $l1boa | sed "s/_/\n/g;s/.h5//" | tail -1`
  list=""
  if [ "$lightspeed" = "yes" -a "$testl" != "" ]
  then
    list="lightspeed"
    machine="lightspeed"
  fi
  if [ "$scramjet" = "yes" -a "$tests" != "" ]
  then
    list="scramjet"
    machine="scramjet"
  fi
  if [ "$debug" = "yes" ]
  then
    echo "testl: $testl"
    echo "tests: $tests"
    echo "statbad: $statbad"
    echo "statgood: $statgood"
    echo "l1boa: $l1boa"
    echo "date: $date"
    echo "list: $list"
  fi
  if [ "$list" != "" ]
  then
    skipifrestricting="$restrict"
    if [ "$statbad" != "" ]
    then
      newlist="$list \t ended badly"
    elif [ "$statgood" != "" ]
    then
      newlist="$list \t completed normally"
    else
      newlist="$list \t  still running"
      skipifrestricting="no"
    fi
    if [ "$restrict" = "no" ]
    then
      list="$newlist"
    fi
    if [ "$convert" = "yes" ]
    then
      list="`$DATECONVERTER $date` \t  $list"
    else
      list="$date \t $list"
    fi
    if [ "$dirnames" = "yes" ]
    then
      list="$list \t $dir"
    fi
    # These other test are slower, so skip them if we won't print the results
    if [ "$died" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep -i 'died,' $dir/exec_log/process.stdout | awk '{print $7}'`
	   nodes=`grep -i 'died,' $dir/exec_log/process.stdout | awk '{print $9}'`
      list="$list \t $chunks \t $nodes"
    fi
    if [ "$running" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep -i completed $dir/exec_log/process.stdout | tail -1 | \
        awk '{print $8}'`
      if [ "$chunks" = "" ]
      then
        chunks="(unknown)"
      fi
      list="$list \t $chunks"
    fi
    if [ "$full" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep -i completed $dir/exec_log/process.stdout | tail -1 | \
        awk '{print $8, $3, $10, $12, $5}'`
      if [ "$chunks" = "" ]
      then
        chunks="(unknown)"
      fi
      list="$list \t $chunks"
    fi
    if [ "$time" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep -i 'starting mlsl2' $dir/exec_log/process.stdout | head -1 | \
        awk '{print $3, $4}'`
      list="$list \t $chunks"
    fi
    if [ "$fail" = "yes" -a "$skipifrestricting" = "no" ]
    then
      echo "Nodes on $machine on which a job could not be spawned"
      a=`grep -n 'Unable to start slave task' $dir/exec_log/process.stdout | \
       awk '{print $1}' | sed 's/://'`
      #echo "a: $a"
      if [ "$a" = "" ]
      then
        echo "All attempts at spawning jobs succeeded"
        exit
      fi
      note_failures "$a" "$dir/exec_log/process.stdout"
    fi
    if [ "$skipifrestricting" = "no" -a $fail = "no" ]
    then
      echo -e $list
    fi
  fi
done
exit 0
# $Log$
# Revision 1.2  2005/04/06 22:13:58  pwagner
# Fixed bug; optionally show unspawnable nodes
#
# Revision 1.1  2005/04/01 00:13:15  pwagner
# First commit
#
