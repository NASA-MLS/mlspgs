#!/bin/sh
# match_date.sh
#
# --------------- match_date.sh help
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
# Note: We rely on one of the following (in descending priority)
# (1) Define the env variable DATECONVERTER
# (2) dateconverter being in your PATH
# (3) MLSTOOLS dfined as an env variable, and
#     $MLSTOOLS/dateconverter exists and you have execute permisiion for it
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

  # Method: dateconverter will compare two dates: date1 and date2
  # if date1 is after date2 it returns '+'
  # otherwise '-' or '='
  # So we'll iterate over each choice_n until
  # we find the first that produces '+'
  # The one before that one is what we want
  if [ "$DATECONVERTER" = "" ]
  then
    DATECONVERTER=`which dateconverter`
  fi
  if [ "$DATECONVERTER" = "" ]
  then
    DATECONVERTER=$MLSTOOLS/dateconverter
  fi
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
  exit 0
# $Log$
