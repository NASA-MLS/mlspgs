#!/bin/bash
# --------------- nodestatus help
# Check on the number of processors on a node at the clusters at the sips
# determining whether the number is as expected (usu. 2)
# If not, pvm delete the node
# Optionally, mail this fact to me
# Usage:
# nodestatus.sh [options]  [nodes] ..
#    O p t i o n s
# -dryrun     don't run the actual commands
# -nodelete   don't delete the deficient nodes
# -v          verbose--show commands, results even when normal
# -S          silent--run (more) silently
# -H HOSTER   use HOSTER instead of default pvmhost.sh to delete host(s)
# -F address  show mail as having come from address
# -M MAILER   use MAILER instead of default mailtome.sh to mail report(s)
# -R list     change RECIPIENTS to list
# -mail       mail report(s) to RECIPIENTS
# -n number   change expected number of processors to number
# -f file     read node names from file; e.g., PVM_ALL_SLAVES
# -h[elp]     print brief help message; exit
#
# Note:
# (1) executable HOSTER and MAILER scripts must exist and be in your path
# (2) passwordless ssh must be set up
# (3) the number of processors must be findable from /proc/cpuinfo
# (4) -v and -S options uneasily coexist--
#     best to treat them as mutually exclusive
# (5) if the file $HOME/ncpus.txt exists it will be read to change the
#      expected number of processors (may be overridden by -n option)
# 
# Result:
# Won't try to spawn jobs to deficient nodes
# --------------- End nodestatus help
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

#
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
    #  echo $our_host_name
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
}
      
#
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
      

# ************
# Main Program
# ************
# 
# return values from HOSTER regarding pvmhost's presence in pvmtable
PRESENT=0
ABSENT=1
BUNDLE="yes"
SORTUNIQ="yes"
# For debugging, only
# godzilla="compute-0-181"
# godzilla="c0-181"

me="$0"
my_name=nodestatus
I=nodestatus

# Initial settings--may be overridden by command-line options
FROMME="Paul Wagner <paul.a.wagner@jpl.nasa.gov>"
RECIPIENTS="pwagner sysadmin@sdsio.jpl.nasa.gov"
# RECIPIENTS="pwagner"
# MAILER="/home/pwagner/bin/mailtome.sh"
# HOSTER="/home/pwagner/bin/pvmhost.sh"
MAILER="mailtome.sh"
HOSTER="pvmhost.sh"
delete="yes"
dryrun="no"
expected="2"
# Is there a file to read to change the expected number of processors?
file=$HOME/ncpus.txt
if [ -f "$file" ]
then
  expected=`cat $file | read_file_into_array`
fi
file=""
mail="no"
silent="no"
verbose="no"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -nodelete )
	    shift
       delete="no"
       ;;
    -dryrun )
	    shift
       dryrun="yes"
       ;;
    -v )
	    shift
       verbose="yes"
       ;;
    -S )
	    shift
       silent="yes"
       ;;
    -f )
	    shift
       file="$1"
       shift
       ;;
    -F )
	    shift
       FROMME="$1"
       shift
       ;;
    -R )
	    shift
       RECIPIENTS="$1"
       shift
       ;;
    -M )
	    shift
       MAILER="$1"
       shift
       ;;
    -H )
	    shift
       HOSTER="$1"
       shift
       ;;
    -mail )
	    shift
       mail="yes"
       ;;
    -n )
	    shift
       expected="$1"
       shift
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

TODAY=`date +"%Y-%m-%d"`
MAILBUNDLE=`get_unique_name $TODAY`
username=`whoami`
pvmrunning=`ps -auxxxx | grep -i pvmd3 | grep $username | grep -v grep`

# Where can we write the MAILBUNDLE? Current directory? HOME?
if [ -w "`pwd`" ]
then
  BUNDLEDIR=`pwd`
elif [ -w "$HOME" ]
then
  BUNDLEDIR="$HOME"
else
  # Oops--we'll need to mail reports one-by-one afterall
  BUNDLE="no"
fi
MAILBUNDLE="$BUNDLEDIR/$MAILBUNDLE"

if [ "$file" = "" ]
then
  hosts="$@"
elif [ "$SORTUNIQ" = "yes" ]
then
  hosts=`cat $file | sort | uniq | read_file_into_array`
else
  hosts=`read_file_into_array < $file`
fi

if [ "$verbose" = "yes" ]
then
  echo $TODAY
  echo "userame $username"
  echo "pvmrunning $pvmrunning"
  if [ "$mail" = "yes" -a "$BUNDLE" = "yes" ]
  then
    echo "MAILBUNDLE file $MAILBUNDLE"
  fi
  if [ "$file" != "" ]
  then
    echo "pvmhosts file $file"
  fi
  echo "pvmhosts $hosts"
fi

if [ "$mail" = "yes" -a "$BUNDLE" = "yes" ]
then
  echo "$HOSTNAME nodestatus $TODAY" > "$MAILBUNDLE"
fi

for pvmhost in $hosts
do
  #Check whether the number of processors matches the expected number
  actual=`ssh -a -x $pvmhost 'grep -i processor /proc/cpuinfo | wc -l' | \
    awk '{print $1}'`
  # This is just for debugging purposes
  # So don't forget to disable it before turning the creature loose
  # on the world, or you'll end up with a scenario like ..
  # "Godzilla is marching through Tokyo, leaving a path of
  # destruction, downed powerlines, and general mayhem, and it's
  # all your fault, Mr. Wagner!"
  # if [ "$pvmhost" = "$godzilla" ]
  # then
  #   actual=999
  # fi

  if [ "$actual" != "$expected" ]
  then
    if [ "$dryrun" = "yes" ]
    then
       #   d r y r u n   s e c t i o n
      echo "Expected $expected, got $actual on $pvmhost"
      if [ "$mail" = "yes" ]
      then
        echo $MAILER -s "$HOSTNAME nodestatus $TODAY" \
        -f "$FROMME" -R "$RECIPIENTS" "$OUTPUT"
      fi
      if [ "$pvmrunning" != "" ]
      then
        hoststatus=`$HOSTER -s $pvmhost`
        if [ "$hoststatus" = "$PRESENT" ]
        then
          echo $HOSTER -delete $pvmhost
        else
          echo $pvmhost already deleted or never added
        fi
      else
        echo pvm not currently running
      fi
    else
       #   l i v e f i r e   s e c t i o n
      if [ "$silent" != "yes" ]
      then
        echo "Expected $expected, got $actual on $pvmhost"
      fi
      if [ "$pvmrunning" != "" ]
      then
        # 1st--check if the deficient node needs to be deleted from pvmtable
        # so we won't try to delete the same host more than once
        # nor send multiple emails about it
        hoststatus=`$HOSTER -s $pvmhost`
        if [ "$hoststatus" = "$PRESENT" -a "$delete" = "yes" ]
        then
          $HOSTER -delete $pvmhost
        fi
      else
        hoststatus="$ABSENT"
      fi
      if [ "$mail" = "yes" -a "$hoststatus" = "$PRESENT" ]
      then
        if [ "$BUNDLE" = "yes" ]
        then
          # Accumulate the deficient node reports in MAILBUNDLE
          echo "Expected $expected, got $actual on $pvmhost" >> "$MAILBUNDLE" 
        else
          # Mail the reports one-by-one
          echo "Expected $expected, got $actual on $pvmhost" | \
          $MAILER -s "$HOSTNAME nodestatus $TODAY" \
          -f "$FROMME" -R "$RECIPIENTS" "$OUTPUT"
        fi
      fi
    fi
  elif [ "$verbose" = "yes" ]
  then
    echo "Expected $expected, got $actual on $pvmhost"
  fi
done

# Do we have a MAILBUNDLE to mail off?
if [ "$mail" = "yes" -a "$BUNDLE" = "yes" ]
then
  # Does it have any lines besides the 1st one? (Which repeats the subject)
  actual=`wc -l "$MAILBUNDLE" | awk '{print $1}'`
  if [ "$actual" -gt 1 ]
  then
    cat "$MAILBUNDLE" | \
    $MAILER -s "$HOSTNAME nodestatus $TODAY" \
    -f "$FROMME" -R "$RECIPIENTS" "$OUTPUT"
  fi
  /bin/rm -f "$MAILBUNDLE"
fi

exit 0
# $Log$
# Revision 1.1  2005/10/01 00:36:28  pwagner
# First commit
#
