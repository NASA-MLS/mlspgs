#!/bin/bash
# --------------- sipsl2 help
# Show level 2 jobs running on lightspeed, scramjet
# Usage:
# sipsl2.sh [options]
#    O p t i o n s
# -h[elp]     print brief help message; exit
# -m cluster  show jobs only for cluster named cluster
# -Cf file    choose clusters named in file
# -cf file    get expected compute times from file
# -vn version show jobs only for version (e.g., V01-51)
# -bug        show chunks lost to level 2 bug "list out of order in hunt"
# -debug      print lots of extra debugging inof
# -c          convert dates from yyyydoy to yyyy Month day
# -d          show directory names
# -x          show chunks, nodes that died
# -D d1,d2,   ignore jobs in directories d1, d2 (declared "legally dead")
# -Df file,   jobs in directories named in file are declared "legally dead"
#   The following options carry an automatic restriction to running jobs
# -r          show how many chunks running
# -full       show full stats: running, completed, crashed, on deck, total
# -finish     show date, time the run should complete (approximate)
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
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

#---------------------------- cat_args
#
# Function to catenate args with optional third arg separating them,
# else a space
# Treat special cases where 1st or second args are blanks sensibly

cat_args()
{

   # Do we have enough args?
      separator=" "
      if [ $# -lt 2 ]
      then
         echo "Usage: cat_args arg1 arg2 [arg3]"
         exit 1
      elif [ $# -eq 3 ]
      then
         separator="$3"
      fi
      if [ "$1" = "" ]
      then
        echo $2
      elif [ "$2" = "" ]
      then
        echo $1
      else
        echo "$1${separator}$2"
      fi
}
      
#---------------------------- how_many_array_elements
#
# Function to return the number of elements in a shell array
# where the array is the lone arg
# eg given
# how_many_array_elements 'a b c d'
# writes '4' to standard output
# (Note: we have thus corrected for perl's c-centric array origin)

how_many_array_elements()
{

   # Do we have enough args?
      if [ $# -lt 1 ]
      then
         echo "Usage: how_many_array_elements 'a b c ..'"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); print eval($#parts + 1)' "$1"
}
      
#---------------------------- note_failures
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
    machine_info=`sed -n ''$line_number' p' "$2" | awk '{print $9, $11}'`
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

#---------------------------- print_nth_array_element
#
# Function to print the nth element in a space-delimited list
# where n is the first arg, the list is the second arg, 
# print_nth_array_element 3 'a b c d'
# writes 'c' to standard output

print_nth_array_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: print_nth_array_element n 'a b c ..'"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); print $parts[$ARGV[1]-1]' "$2" "$1"
}

#---------------------------- print_hash_element
#
# Function to print value corresponding to key given two space-delimited arrays
# where the first array are the keys and the second the values
# print_hash_element c 'a b c d' '1 2 3 4'
# writes '3' to standard output
print_hash_element()
{
   # Do we have enough args?
      if [ $# -lt 3 ]
      then
         echo "Usage: print_hash_element c 'a b c ..' '1 2 3 4 ..'"
         exit 1
      fi
      arg_num=`what_array_element $1 "$2"`
      print_nth_array_element $arg_num "$3"
}

#---------------------------- read_file_into_array
#
# read each line of stdin
# catenating them into an array which we will return
# Ignore any entries beginning with '#' character
# In fact, only first entry in each line is kept
# Possible improvements:
#   Other comment signifiers
#   Choose field number other than 1

read_file_into_array()
{
  array_result=''
  while read line; do
    element=`echo $line | awk '$1 !~ /^#/ {print $1}'`
    if [ "$element" != "" ]
    then
      array_result="$array_result $element"
    fi
  done
  echo $array_result
}
      
#---------------------------- what_array_element
#
# Function to show what element in a space-delimited list has been supplied
# in the first arg, the list is the second arg, 
# what_array_element c 'a b c d'
# writes '3' to standard output
what_array_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: what_array_element c 'a b c ..'"
         exit 1
      fi
      n=0
      count_up=1
      for i in $2
      do
        # echo $i $1
        if [ "$i" = "$1" ]
        then
          # echo "really, $i $1"
          n=$count_up
        fi
        count_up=`expr $count_up + 1`
      done
      echo $n
}

# ************
# Main Program
# ************
# 
bug="no"
# This is a shell array of all cluster names
# In our nomenclature, array elements are separated by spaces
# List elements are separated by commas
clusternames="lightspeed scramjet speedracer"
computetimes=""
convert="no"
corpses=""
debug="no"
died="no"
dirnames="no"
dontprintemptydates="yes"
#                      ^------ skip printing from where sips removed l1b files
I=sipsl2
fail="no"
finish="no"
full="no"
maxobituaries="24"
#               ^------ stop printing dead chunks, nodes past this
me="$0"
my_name=sipsl2
restrict="no"
running="no"
time="no"
version=""
DATECONVERTER=`which dateconverter 2>/dev/null`
REECHO=`which reecho.sh 2>/dev/null`
if [ ! -r "$DATECONVERTER" ]
then
  DATECONVERTER="/home/pwagner/bin/dateconverter"
fi 
if [ ! -r "$REECHO" ]
then
  REECHO="/home/pwagner/bin/reecho.sh"
fi 
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -bug )
	    shift
       bug="yes"
       ;;
    -debug )
	    shift
       debug="yes"
       ;;
    -c )
	    shift
       convert="yes"
       ;;
    -Cf )
	    shift
       clusternames=`cat $1 | uniq | read_file_into_array`
       shift
       ;;
    -cf )
	    shift
       computetimes=`cat $1 | uniq | read_file_into_array`
       shift
       ;;
    -d )
	    shift
       dirnames="yes"
       ;;
    -D )
	    shift
       corpses="$1"
       shift
       ;;
    -Df )
	    shift
       corpses=`cat $1 | uniq | read_file_into_array`
       shift
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
    -finish )
	    shift
       finish="yes"
       restrict="yes"
       ;;
    -full )
	    shift
       full="yes"
       restrict="yes"
       ;;
    -m )
	    shift
       clusternames="$1"
	    shift
       ;;
    -vn )
	    shift
       version="$1"
	    shift
       ;;
    -t )
	    shift
       time="yes"
       restrict="yes"
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
nclusters=`how_many_array_elements "$clusternames"`
dirsplus=`echo /workops/jobs/science/*/mlsl2.slave`
dirs=`echo $dirsplus | sed "s:mlsl2.slave::g"`
if [ "$debug" = "yes" ]
then
  echo "dirs: $dirs"
fi
if [ "$nclusters" = "1" ]
then
  echo "status for cluster $clusternames only"
  list="date of data"
else
  list='date of data \t machine'
fi
if [ "$restrict" != "yes" ]
then
  list=`cat_args "$list" "---- status ----" "\t"`
fi
if [ "$dirnames" = "yes" ]
then
  list=`cat_args "$list" "directory" "\t"`
fi
if [ "$bug" = "yes" ]
then
  list=`cat_args "$list" "lost to bug" "\t"`
fi
if [ "$died" = "yes" ]
then
  list=`cat_args "$list" "dead chunks \t killer nodes" "\t"`
fi
if [ "$running" = "yes" ]
then
  list=`cat_args "$list" "chunks running" "\t"`
fi
if [ "$full" = "yes" ]
then
  echo "Stats order: running complete crashed waiting total"
  list=`cat_args "$list" "full stats " "\t"`
fi
if [ "$time" = "yes" ]
then
  list=`cat_args "$list" "run start (date time)" "\t"`
fi
if [ "$finish" = "yes" ]
then
  list=`cat_args "$list" "run finish (date time)" "\t"`
fi
if [ "$fail" = "no" ]
then
  echo -e $list
fi

#
# for machine in lightspeed scramjet speedracer unknown
# do
# computetime=`print_hash_element $machine "$clusternames" "$computetimes"`
# echo "$machine computetime $computetime"
# done
# exit

afterfirst="no"
for dir in $dirs
do
  # These are pretty quick (at least compared with grepping the whole
  # file as later tests require)
  machine=""
  for name in $clusternames
  do
    # testl=`grep -i executing $dir/exec_log/process.stdout | grep -i $name`
    testl=`head -40 $dir/exec_log/process.stdout | grep -n Executing  | grep -i $name`
    if [ "$testl" != "" ]
    then
      machine="$name"
    fi
  done
  theversion=`grep -i pgeversion $dir/job.PCF | awk -F"|" '{print $3}'`
  case $theversion in
    V01-5*)
      computetime=24
      ;;
    V02-20)
      computetime=30
      ;;
    V02-21)
      computetime=33
      ;;
    *)
      computetime="depends"
      ;;
  esac
  statbad=`tail $dir/exec_log/process.stdout | grep -i "ended badly"`
  # Some other signs that the job ended badly are that it was killed or
  # the network went down
  if [ "$statbad" = "" ]
  then
    statbad=`tail $dir/exec_log/process.stderr | \
      grep -i "/home/sips/ops/bin/dispatch_l2" | grep -i killed`
  fi
  if [ "$statbad" = "" ]
  then
    statbad=`tail $dir/exec_log/process.stderr | grep -i "connection lost"`
  fi
  if [ "$statbad" = "" ]
  then
    statbad=`tail $dir/exec_log/process.stderr | grep -i "terminated"`
  fi
  if [ "$statbad" = "" -a "$corpses" != "" ]
  then
    statbad=`echo "$corpses" | grep -i "$dir"`
    if [ "$statbad" != "" ]
    then
      statbad="legally dead"
    fi
  fi
  statnochunks=`tail $dir/exec_log/process.stdout | grep -i "No chunks were processed"`
  statpvmtrouble=`tail $dir/exec_log/process.stdout | grep -i "probably pvm trouble"`
  statgood=`tail $dir/exec_log/process.stdout | grep -i "catenating chunk"`
  l1boa=`$REECHO $dir/*L1BOA_*`
  date=`echo $l1boa | sed "s/_/\n/g;s/.h5//" | tail -1`
  bugs=`grep -ic 'list out of order' $dir/pvmlog/mlsl2.log 2>/dev/null`
  if [ "$bugs" = "" ]
  then
    bugs=0
  fi
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
  if [ "$machine" != "" -a "$nclusters" != "1" ]
  then
    list="$machine"
  fi
  if [ "$debug" = "yes" ]
  then
    echo "machine: $machine"
    echo "version sought: $version"
    echo "version found: $theversion"
    echo "statbad: $statbad"
    echo "statnochunks: $statnochunks"
    echo "statpvmtrouble: $statpvmtrouble"
    echo "statgood: $statgood"
    echo "l1boa: $l1boa"
    echo "date: $date"
    echo "bugs: $bugs"
    echo "list: $list"
    echo "computetime: $computetime"
  fi
  if [ "$dontprintemptydates" = "yes" -a "$date" = "" ]
  then
    machine=""
  fi
  if [ "$version" != "" -a "$version" != "$theversion" ]
  then
    machine=""
  fi
  if [ "$machine" != "" ]
  then
    skipifrestricting="$restrict"
    if [ "$statbad" != "" ]
    then
      statkilled=`echo "$statbad" | grep -i killed`
      statterminated=`echo "$statbad" | grep -i terminated`
      statlegallydead=`echo "$statbad" | grep -i legally`
      if [ "$statkilled" != "" ]
      then
        newlist=`cat_args "$list" "killed" "\t"`
      elif [ "$statterminated" != "" ]
      then
        newlist=`cat_args "$list" "terminated" "\t"`
      elif [ "$statlegallydead" != "" ]
      then
        newlist=`cat_args "$list" "legally dead" "\t"`
      else
        newlist=`cat_args "$list" "ended badly" "\t"`
      fi
    elif [ "$statnochunks" != "" ]
    then
      newlist=`cat_args "$list" "failed (no chunks)" "\t"`
    elif [ "$statpvmtrouble" != "" ]
    then
      newlist=`cat_args "$list" "failed (pvm)" "\t"`
    elif [ "$statgood" != "" ]
    then
      newlist=`cat_args "$list" "completed normally" "\t"`
    else
      newlist=`cat_args "$list" "still running" "\t"`
      skipifrestricting="no"
    fi
    if [ "$restrict" = "no" ]
    then
      list="$newlist"
    fi
    if [ "$date" = '' ]
    then
      list="(empty date) \t  $list"
    elif [ "$convert" = "yes" ]
    then
      list="`$DATECONVERTER $date` \t  $list"
    else
      list="$date \t $list"
    fi
    if [ "$dirnames" = "yes" ]
    then
      list="$list \t $dir"
    fi
    if [ "$bug" = "yes" ]
    then
      list="$list \t $bugs"
    fi
    if [ "$debug" = "yes" ]
    then
      echo "skipifrestricting: $skipifrestricting"
      echo "list: $list"
      echo "newlist: $newlist"
    fi
    # These other test are slower, so skip them if we won't print the results
    if [ "$died" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep 'died,' $dir/exec_log/process.stdout | awk '{print $7}'`
	   nodes=`grep 'died,' $dir/exec_log/process.stdout | awk '{print $9}'`
      ndied=`echo $chunks | wc | awk '{print $2}'`
      if [ "$ndied" -gt "$maxobituaries" ]
      then
        list="$list \t ($ndied chunks) \t ($ndied nodes)"
      else
        list="$list \t $chunks \t $nodes"
      fi
    fi
    if [ "$running" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep completed $dir/exec_log/process.stdout | tail -1 | \
        awk '{print $8}'`
      if [ "$chunks" = "" ]
      then
        chunks="(unknown)"
      fi
      list="$list \t $chunks"
    fi
    if [ "$full" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep completed $dir/exec_log/process.stdout | tail -1 | \
        awk '{print $8, $3, $10, $12, $5}'`
      if [ "$chunks" = "" ]
      then
        chunks="(unknown)"
      fi
      list="$list \t $chunks"
    fi
    if [ "$time" = "yes" -a "$skipifrestricting" = "no" ]
    then
	   chunks=`grep 'starting mlsl2' $dir/exec_log/process.stdout | head -1 | \
        awk '{print $3, $4}'`
      list="$list \t $chunks"
    fi
    if [ "$finish" = "yes" -a "$skipifrestricting" = "no" ]
    then
      numChunks=`grep -i 'nochunks:' $dir/exec_log/process.stdout | tail -1 | awk '{print $4}'`
      # numLaunched=`grep -i launched $dir/exec_log/process.stdout | wc -l`
      numLaunched=` grep Launched $dir/exec_log/process.stdout | tail -1 | sed -n 's/.*chunk//p' | awk '{print $1}'`
      if [ "$numChunks" = "" ]
      then
        # Not all chunks launched yet, so can't tell when last one will finish
        chunks="(unknown)"
      elif [ "$numLaunched" = "" ]
      then
        # Not all chunks launched yet, so can't tell when last one will finish
        # chunks="0 / $numChunks"
        chunks=" "
      elif [ "$numLaunched" -lt "$numChunks" ]
      then
        # Not all chunks launched yet, so can't tell when last one will finish
        # chunks="$numLaunched / $numChunks"
        chunks=" "
      else
        # Try to find when when last chunk started
        lastline=`grep -h 'starting mlsl2' $dir/pvmlog/* | sort  | tail -1`
        lastdate=`echo $lastline | awk '{print $1}'`
        lasttime=`echo $lastline | awk '{print $2}'`
        # And add an average compute time to it
        if [ "$computetime" = "depends" ]
        then
          computetime=`print_hash_element $machine "$clusternames" "$computetimes"`
        fi
        chunks=`$DATECONVERTER +H $computetime -t $lasttime -o 'M dd' $lastdate`
      fi
      list="$list \t $chunks"
    fi
    if [ "$fail" = "yes" -a "$skipifrestricting" = "no" ]
    then
      if [ "$afterfirst" = "yes" ]
      then
        echo ""
      fi
      afterfirst="yes"
      echo "Nodes on $machine on which a job could not be spawned"
      a=`grep -n 'Unable to start slave task' $dir/exec_log/process.stdout | \
       awk '{print $1}' | sed 's/://'`
      if [ "$a" = "" ]
      then
        echo "All attempts at spawning jobs succeeded"
      else
        note_failures "$a" "$dir/exec_log/process.stdout"
      fi
    fi
    if [ "$skipifrestricting" = "no" -a $fail = "no" ]
    then
      echo -e $list
    fi
  fi
done
exit 0
# $Log$
# Revision 1.17  2007/01/22 17:48:04  pwagner
# Repaired bug affecting runs just starting out
#
# Revision 1.16  2007/01/18 23:30:12  pwagner
# Added -finish option to display predicted finish time
#
# Revision 1.15  2006/08/21 21:56:09  pwagner
# Works better with v2.1 versions
#
# Revision 1.14  2006/07/13 18:12:18  pwagner
# Accepts that certain jobs (corpses) may be legally declared dead
#
# Revision 1.13  2006/06/28 00:03:40  pwagner
# Notes when jobs are terminated rather than completed normally
#
# Revision 1.12  2006/05/19 19:56:25  pwagner
# Catches sign that a job was killed
#
# Revision 1.11  2006/04/17 23:22:31  pwagner
# Tests for another condition indicating ended badly
#
# Revision 1.10  2006/03/23 19:18:01  pwagner
# Handles multiple pge versions explicitly
#
# Revision 1.9  2006/03/18 00:06:45  pwagner
# Avoids printing empty dates, limits printing dead chunks
#
# Revision 1.8  2005/07/15 22:42:43  pwagner
# Will handle empty date fields more gracefully
#
# Revision 1.7  2005/06/23 22:20:46  pwagner
# Reworded Copyright statement
#
# Revision 1.6  2005/05/20 23:10:27  pwagner
# Show how many chunks succumbed to 'list out of order' bug
#
# Revision 1.5  2005/04/28 18:42:32  pwagner
# Fixed bug preventing more than one machine  from reporting pvm failures
#
# Revision 1.4  2005/04/21 20:38:02  pwagner
# Replace -l, -s options with -m; added speedracer name failed status
#
# Revision 1.3  2005/04/12 20:33:26  pwagner
# Tightened criteria for dead chunks, killer nodes
#
# Revision 1.2  2005/04/06 22:13:58  pwagner
# Fixed bug; optionally show unspawnable nodes
#
# Revision 1.1  2005/04/01 00:13:15  pwagner
# First commit
#
