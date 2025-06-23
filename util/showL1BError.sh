#!/bin/sh
# --------------- showL1BErrors.sh help
# showL1BErrors.sh
# Shows the L1B errors in the LOGSTATUS file for each day in a given year 

#
# usage: showL1BErrors.sh [options] [years]

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -silent       Prevent most output during execution
#    -Ef env_file  Get options by sourcing env_file

# You might run it like this
# /users/pwagner/mlspgs/util/showL1BError.sh -Ef v5.02.env 2005 2006 2009 2012 2013 2018 2019 2020 2021 2022 2023 2024

# where the file v5.02.env would contain
##------env file for showL1BErrors.sh----
#dryrun=yes
#debug=1
#silent=no
#version=v05.02
#L1BDir=/data/emls/l1b/$version
# (without the leading #'s)

# --------------- End showL1BErrors.sh help
# Copyright 2024, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

#---------------------------- logit
# Either echo input string(s), or else append them to the master log file, 
# or both
# Which to do depends on $MUSTLOG

logit()
{
   if [ "$dryrun" = "yes" ]
   then
     echo "$@"
   else
     case "$MUSTLOG" in
     no )
       echo "$@"
       ;;
     yes )
       echo "$@" >> "$MASTERLOG"
       ;;
     both )
       echo "$@" | tee -a "$MASTERLOG"
       ;;
     esac
   fi
}

#---------------------------- magnify
# Returns the file name with highest cycle number
#

magnify()
{
   echo *$1*.he5 | awk '{print $NF}'
}
      
#---------------------------- commando
# Either execute a command with all its args, or merely echo the command that
# would be executed
# Which to do depends on $dryrun

commando()
{
   if [ "$dryrun" = "yes" ]
   then
     logit $@
   else
     $@
   fi
}
      
#---------------------------- showIt
# Add the L1BDay to the list of dates when the L1B filee contains
# the error code.
# Unless you sent silent to yes, print the messages

showIt()
{
   if [ "$silent" != "yes" ]
   then
     grep "$ERRORCODE" "$L1B"
   fi
   b=`grep "$ERRORCODE" "$L1B"`
   if [ "$b" != "" ]
   then
     list="$list $L1BDay"
   fi
}
      
# ------------------ one_day ----------------------
# ----------------------------------------------------------------
# Process each day in one day's worth of files.
one_day()
{
  # Showing an entire day of files
  files=`/bin/ls $L1BDay`
  for file in $files
  do
    L1B=$L1BDay/$file
    a=`echo $file | grep 'LOGSTATUS'`
    if [ "$a" != "" ]
    then
      if [ "$debug" = "1" ]
      then
        echo $file
      fi
      showIt
      # repack_files
      # augment_files
    fi
  done
}

# ------------------ one_year ----------------------
# ----------------------------------------------------------------
# Process each day in one year's worth of files.
one_year()
{
  # Repairing an entire year of files
  days=`/bin/ls $L1BYear`
  for day in $days
  do
    if [ ! -d "$L1BYear/$day" ]
    then
      echo "Sorry--$L1BYear/$day is not a subdirectory containing L1B files"
    else
      # Get files with highest cycle numbers (except for L1B)
      thisDir=`pwd`
      L1BDay=$L1BYear/$day
      one_day
    fi
  done
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script    *
#****************************************************************
#
debug=1
#     ^  -- set this to 1 if worried

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MUSTLOG="no"
MASTERLOG=/dev/null

me="$0"
my_name=showL1BErrors.sh
I=showL1BErrors
insertoptions="-v"
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
dryrun="no"
more_opts="yes"
L1BDir=""
L1B=""
ERRORCODE=PGSEPH_E_BAD_EPHEM_FILE_HDR
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -silent )
	    insertoptions="-silent"
	    shift
       ;;
    -Dd )
	    L1BDir="$2"
	    shift
	    shift
       ;;
    -Df )
	    L1B="$2"
	    shift
	    shift
       ;;
    -Ef )
	    . "$2"
	    shift
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

years="$@"
if [ "$debug" = 1 ]
then
  echo "dryrun $dryrun"
  echo "L1BDir $L1BDir"
  echo "L1B file $L1B"
  echo "years $years"
fi

#echo "Past (1), anyway"
if [ "$L1BDir" != "" -a ! -d "$L1BDir" ]
then
  echo "Sorry--L1BDir must be a directory"
  exit
fi

list=""
#echo "Past (2), anyway"
if [ "$L1BDir" = "" ]
then
  # Repairing a single file
  echo "Repairing $L1B"
  # Can't repair a non-existent L1B file
  if [ ! -f "$L1B" ]
  then
    echo "No standard product L1B file defined; did you use -Df option?"
    exit
  else
    showIt
  fi
elif [ "$years" = "" ]
then
  # L1BDir contain the year datum;
  # you're doing just one year
  L1BYear=$L1BDir
  one_year
else
  # L1BDir are just the version;
  # you're doing $years
  for year in $years
  do
    logit "year     $year"
    L1BYear=$L1BDir/$year
    one_year
  done
fi
for dates in $list
do
  echo $dates
done
# $Log$
