#!/bin/sh
# tkreset.sh
# Usage:
# (1)
# . tkreset.sh path/to/new/tk/pgs-env.csh
# This will replace previous toolkit-related paths
# with those specified by the file path/to/new/tk/pgs-env.csh
# or else
# (2)
# . tkreset.sh
# This will remove previous toolkit-related paths
# It also snips duplicate paths
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
echo "We were called as $0"
standsalone=`echo $0 | grep -i tkreset`
echo "tkreset.sh called with $@"
if [ "$PGSHOME" = "" ]
  then
  # No previous toolkit paths
  echo "No previous PGSHOME defined"
  if [ -f "$1" -a "$standsalone" != "" ]
    then
     . $1
  fi
  echo "new PATH is $PATH"
  echo "new PGSHOME is $PGSHOME"
else
  oldPATH=`echo $PATH | sed 's^:^ ^g'`
  oldhome="$PGSHOME"
  echo "oldPath is $oldPATH"
  echo "oldhome is $oldhome"
  usage="2"
  if [ -f "$1" -a "$standsalone" != "" ]
    then
    usage="1"
     . $1
  fi
  echo "usage is $usage"
  PATH=""
  checkdpths=""
  for arg in $oldPATH
  do
    echo "arg is $arg"
    repl=""
    repeat=""
    match=`echo $arg | /bin/grep "$oldhome"`
    if [ "$checkdpths" != "" ]
    then
      for oldarg in $checkdpths
      do
        if [ "$arg" = "$oldarg" ]
        then
          repeat="yes"
        fi
      done
      checkdpths="$checkdpths $arg"
    else
      checkdpths=$arg
    fi
    if [ "$repeat" = "yes" ]
    then
      echo "Snipping duplicate path $arg"
    elif [ "$usage" = "2" ]
    then
      # Usage (2)
      if [ "$match" = "" ]
      then
        if [ "$PATH" = "" ]
        then
          PATH=$arg
        else
          PATH="${PATH}:$arg"
        fi
      fi
    else
      # Usage (1)
      if [ "$match" = "" ]
      then
        if [ "$PATH" = "" ]
        then
          PATH=$arg
        else
          PATH="${PATH}:$arg"
        fi
      else
        repl=`echo $arg | /bin/sed "s^$oldhome^$PGSHOME^"`
        if [ "$PATH" = "" ]
        then
          PATH=$repl
        else
          PATH="${PATH}:$repl"
        fi
      fi
    fi
  done
  if [ "$usage" = "2" ]
  then
    PGSHOME=
  else
    echo "new PGSHOME is $PGSHOME"
  fi
  echo "new path is $PATH"
  export PATH
  export PGSHOME
  # housekeeping (because we may "dot" this script)
  unset oldPATH
  unset oldhome
  unset usage
  unset checkdpths
  unset arg
  unset repl
  unset repeat
  unset match
  unset oldarg
fi
# $Log$
# Revision 1.3  2005/06/23 22:20:46  pwagner
# Reworded Copyright statement
#
# Revision 1.2  2003/05/08 20:09:07  pwagner
# Removed exit commands where possible
#
# Revision 1.1  2003/05/06 20:40:52  pwagner
# First commit
#
