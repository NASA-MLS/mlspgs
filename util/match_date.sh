#!/bin/sh
# match_date.sh
#
# --------------- match_date.sh help
# Use 1:
# choose the date that best matches the target
# given a sorted set of choices
# E.g., if the target is 2018d015
# and the choices are 
# 2016d364 2018d016 2018364
# return 2016d364
#
# Usage:
# match_date.sh day choice_1 [choice_2] ..
#
# 
# Result:
# choice_n such that
# choice_n < day < choice_(n+1)
#
# Use 2:
# Check if the date matches any of the dates read
# from a file named on the commandline. The file is a text formatted
# like the following
#
# # file_name       # ignore comments
# 2018d002          # a bare date
# 2018d006          # another bare date
# 2018d217-2018d224 # a range of dates
#
# Usage:
# match_date.sh [options] -Df file_name day
# 
# Result:
# Returns day if match, " " if not
#
## Note: We rely on one of the following (in descending priority)
# (1) Define the env variable DATECONVERTER
# (2) dateconverter being in your PATH
# (3) MLSTOOLS defined as an env variable, and
#     $MLSTOOLS/dateconverter exists and you have execute permisiion for it
# (4) options, if any, must precede the other args
#
# Bugs and limitations:
# (1) We fail silently if DATECONVERTER can't be found or executed
# (2) All dates must be in a format DATECONVERTER groks, e.g. yyyyDddd
# (3) All dates must be within the range 1993d001 .. 2999d365
# --------------- End match_date.sh help
# Copyright 2018, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

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

# -------------------- Choose_match -------------
Choose_match()
{
  # Method: dateconverter will compare two dates: date1 and date2
  # if date1 is after date2 it returns '+'
  # otherwise '-' or '='
  # So we'll iterate over each choice_n until
  # we find the first that produces '+'
  # The one before that one is what we want
  date=$1
  shift
  # iterate over choice_n
  choose=$1
  for choice
  do
    #echo $DATECONVERTER -c $date $choice
    #$DATECONVERTER -c $date $choice
    if [ "`$DATECONVERTER -c $date $choice`" = '+' ]
    then
      # went too far, man
      too_far=1
    else
      choose=$choice
    fi
  done
  echo $choose
}

# -------------------- Read_match -------------
Read_match()
{
  # Method: dateconverter will compare two dates: date1 and date2
  # if date1 is after date2 it returns '+'
  # otherwise '-' or '='
  # So we'll iterate over each choice_n looking for one
  # that produces '='
  # If we find one, we echo it; otherwise simply exit
  date=$1
  choices=`read_file_into_array < $2`
  # iterate over choice_n
  choose=$1
  if [ "$debug" = "yes" ]
  then
    echo "Read_match $1 $2"
    echo "choose $choose"
    echo "choices $choices"
  fi
  for choice in $choices
  do
    # Does choice represent a single date or a range?
    # We do this by checking for a string like
    #  2018d217-2018d224 # a range of dates
    #         ^^^
    a=`echo $choice | grep '[0-9]-[12]'`
    if [ "$a" = "" ]
    then
      # single date
      if [ "`$DATECONVERTER -c $date $choice`" = '=' ]
      then
        # matched!
        echo $choice
        exit 0
      fi
    else
      # a range; find start and end
      startdate=`echo $choice | awk -F"-" '{print $1}'`
      enddate=`echo $choice | awk -F"-" '{print $2}'`
      if [ "$debug" = "yes" ]
      then
        echo "startdate $startdate"
        echo "enddate $enddate"
        $DATECONVERTER -c $date $startdate
        $DATECONVERTER -c $date $enddate
      fi
      # Is date prior to startdate? If so, we don't care
      if [ "`$DATECONVERTER -c $date $startdate`" = '+' ]
      then
        echo " "
      # Is date after enddate? If so we don't care
      elif [ "`$DATECONVERTER -c $date $enddate`" = '-' ]
      then
        echo " "
      # Aha, we know it must be matching somewhere in this range
      else
        echo $date
        exit 0
      fi
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
cmdline=`echo $0 $@`
me="$0"
my_name="match_date.sh"
debug="no"
dryrun="no"
verbose="no"
datesfile="default.txt"

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    # echo "option: $1"
    case "$1" in

    -Df )
       shift
       datesfile="$1"
       shift
       ;;
    -dryrun )
       dryrun="yes"
       shift
       ;;
    -d )
       debug="yes"
       shift
       ;;
    -v )
       verbose="yes"
       shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       rm -f $settings_file
       exit
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

if [ "$DATECONVERTER" = "" ]
then
  DATECONVERTER=`which dateconverter`
fi
if [ "$DATECONVERTER" = "" ]
then
  DATECONVERTER=$MLSTOOLS/dateconverter
fi

# Which use case are we?
# Do we have a datesfile?
if [ ! -f "$datesfile" ]
then
  Choose_match $@
else
  Read_match $1 $datesfile
fi
exit 0
# $Log$
# Revision 1.1  2018/10/11 00:08:50  pwagner
# First commit
#
