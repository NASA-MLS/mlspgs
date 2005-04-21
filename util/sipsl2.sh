#!/bin/bash
# --------------- sipsl2 help
# Show level 2 jobs running on lightspeed, scramjet
# Usage:
# sipsl2.sh [options]
#    O p t i o n s
# -h[elp]     print brief help message; exit
# -m cluster  show jobs only for cluster named cluster
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
      
#---------------------------- change_nth_array_element
#
# Function to change the nth element in a comma-delimited list
# where n is the first arg, the list is the second arg, 
# and the new value is the third
# eg given
# change_nth_list_element 3 'a b c d' e
# writes 'a b e d' to standard output

change_nth_array_element()
{

   # Do we have enough args?
      if [ $# -lt 3 ]
      then
         echo "Usage: change_nth_array_element n 'a b c ..' newvalue"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); $parts[$ARGV[1]-1]=$ARGV[2]; print join(" ",@parts)' "$2" "$1" "$3"
}
      
#---------------------------- change_nth_list_element
#
# Function to change the nth element in a comma-delimited list
# where n is the first arg, the list is the second arg, 
# and the new value is the third
# eg given
# change_nth_list_element 3 'a,b,c,d' e
# writes 'a,b,e,d' to standard output

change_nth_list_element()
{

   # Do we have enough args?
      if [ $# -lt 3 ]
      then
         echo "Usage: change_nth_list_element n 'a,b,c,..' newvalue"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[0]"); $parts[$ARGV[1]-1]=$ARGV[2]; print join(",",@parts)' "$2" "$1" "$3"
}
      
#---------------------------- convert_array_to_list
#
# Function to convert a space-separated array into a comma-delimited list
# eg given
# convert_array_to_list 'a b c d'
# writes 'a,b,e,d' to standard output
convert_array_to_list()
{

   # Do we have enough args?
      if [ $# -lt 1 ]
      then
         echo "Usage: convert_array_to_list 'a b c d ..'"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); print join(",",@parts)' "$1"
}
      
#---------------------------- convert_list_to_array
#
# Function to convert a comma-delimited list into a space-separated array
# eg given
# convert_list_to_array 'a,b,e,d'
# writes 'a b c d' to standard output
convert_list_to_array()
{

   # Do we have enough args?
      if [ $# -lt 1 ]
      then
         echo "Usage: convert_list_to_array 'a,b,c,..'"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[0]"); print join(" ",@parts)' "$1"
}

#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in pace of "."

get_unique_name()
{

   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
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
      echo $temp${pt}$our_host_name${pt}$$
}
      
#---------------------------- how_many_array_elements
#
# Function to return the number of elements in a shell array
# where the array is the lone arg
# eg given
# how_many_list_elements 'a b c d'
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
      
#---------------------------- how_many_list_elements
#
# Function to return the number of elements in a comma-delimited list
# where the list is the lone arg
# eg given
# how_many_list_elements 'a,b,c,d'
# writes '4' to standard output
# (Note: we have thus corrected for perl's c-centric array origin)

how_many_list_elements()
{

   # Do we have enough args?
      if [ $# -lt 1 ]
      then
         echo "Usage: how_many_list_elements 'a,b,c,..'"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[0]"); print eval($#parts + 1)' "$1"
}
      
#---------------------------- lhs_eq_rhs
#
# Function to assign second arg to first: lhs=rhs
# where (lhs rhs) = ($1 $2)
# if given optional third arg "n"
# assigns $$..$lhs = $$..$rhs

lhs_eq_rhs()
{

   # Do we have write permission in ./?
      if [ ! -w "./" ]
      then
         echo "Sorry--need write permission in ./ to operate"
         exit 1
      fi
      
      lhs_unique_name=`get_unique_name lhs`
      rm -f $lhs_unique_name
      if [ $# -lt 3 -a "$1" != "" ]
      then
         echo "$1"="'"$2"'" > $lhs_unique_name
         . $lhs_unique_name
      else
         mega_buck $3 $1
         lhs=$mega_buck_result
         if [ "$lhs" != "" ] ; then
            mega_buck $3 $2
            rhs=$mega_buck_result
            echo "$lhs"="'"$rhs"'" > $lhs_unique_name
            . $lhs_unique_name
         fi
      fi
     rm -f $lhs_unique_name
}
      
#------------------------------- mega_buck ------------
#
# Function to force multiple-evaluation of its second arg
# i.e., $$..$color (where the number of $ signs is n)
# usage: mega_buck n color

# (uses PID to generate unique name) 

mega_buck()
{

   # Trivial case (n is 0)
   if [ "$1" -lt 1 ]
      then
         mega_buck_result="$2"
   elif [ "$2" = "" ]
      then
         mega_buck_result="$2"
   # Do we have write permission in ./?
   elif [ ! -w "./" ]
      then
         echo "Sorry--need write permission in ./ to operate"
         exit 1
   else
      
      unique_name=`get_unique_name`
      rm -f $unique_name

      number=0
      mega_buck_result=$2
      while [ "$number" -lt "$1" ]
      do
         echo "echo \$arg" | sed 's/arg/'$mega_buck_result'/' > $unique_name
         mega_buck_result=`. $unique_name`
         rm -f $unique_name
         number=`expr $number + 1`
      done
   fi
}

#---------------------------- nth_arg
#
# Function to return the nth (n > 0) arg of the std args after n
# eg given
# nth_arg 3 a b 'c or d' d e ..
# writes 'c or d' to standard output (w/o ' marks)

nth_arg()
{

   # Do we have enough args? No error if not, just blank output
     # echo "Entering nth_arg with args $@"
      the_arg=
      if [ $# -gt 1 ]
      then
         count_up=1
         n=$1
         shift
         while [ $count_up -lt $n ]
         do
           shift
           count_up=`expr $count_up + 1`
         done
         the_arg="$1"
      fi
     # echo "$the_arg"
}
      
#---------------------------- nth_list_element
#
# Function to return the nth (n > 0) element of a comma-delimited list
# where n is arg 1 and the list arg 2
# eg given
# nth_list_element 3 'a,b,c,d'
# writes 'c' to standard output (w/o ' marks)
# (Note: we have thus corrected for perl's c-centric array origin)

nth_list_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: nth_list_element n 'a,b,c,..'"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[1]"); print $parts[$ARGV[0]-1]' $1 "$2"
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
# This is a shell array of all cluster names
# In our nomenclature, array elements are separated by spaces
# List elements are separated by commas
clusternames="lightspeed scramjet speedracer"
shows="yes yes yes"
debug="no"
me="$0"
my_name=sipsl2
I=sipsl2
convert="no"
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
    -m )
	    shift
       clusternames="$1"
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
  list="data date"
else
  list='data date \t machine'
fi
if [ "$restrict" != "yes" ]
then
  list=`cat_args "$list" "status" "\t"`
fi
if [ "$dirnames" = "yes" ]
then
  list=`cat_args "$list" "directory" "\t"`
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
if [ "$fail" = "no" ]
then
  echo -e $list
fi
afterfirst="no"
for dir in $dirs
do
  # These are pretty quick (at least compared with grepping the whole
  # file as later tests require)
  #testl=`head $dir/exec_log/process.stdout | grep -i lightspeed`
  #tests=`head $dir/exec_log/process.stdout | grep -i scramjet`
  machine=""
  for name in $clusternames
  do
    testl=`head $dir/exec_log/process.stdout | grep -i $name`
    if [ "$testl" != "" ]
    then
      machine="$name"
    fi
  done
  statbad=`tail $dir/exec_log/process.stdout | grep -i "ended badly"`
  statnochunks=`tail $dir/exec_log/process.stdout | grep -i "No chunks were processed"`
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
  if [ "$machine" != "" -a "$nclusters" != "1" ]
  then
    list="$machine"
  fi
  if [ "$debug" = "yes" ]
  then
    echo "machine: $machine"
    echo "statbad: $statbad"
    echo "statnochunks: $statnochunks"
    echo "statgood: $statgood"
    echo "l1boa: $l1boa"
    echo "date: $date"
    echo "list: $list"
  fi
  if [ "$machine" != "" ]
  then
    skipifrestricting="$restrict"
    if [ "$statbad" != "" ]
    then
      newlist=`cat_args "$list" "ended badly" "\t"`
    elif [ "$statnochunks" != "" ]
    then
      newlist=`cat_args "$list" "failed" "\t"`
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
# Revision 1.3  2005/04/12 20:33:26  pwagner
# Tightened criteria for dead chunks, killer nodes
#
# Revision 1.2  2005/04/06 22:13:58  pwagner
# Fixed bug; optionally show unspawnable nodes
#
# Revision 1.1  2005/04/01 00:13:15  pwagner
# First commit
#
