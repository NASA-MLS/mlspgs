#!/bin/csh
# tkreset.csh
# Usage:
# (1)
# source tkreset.csh path/to/new/tk/pgs-env.csh
# This will replace previous toolkit-related paths
# with those specified by the file path/to/new/tk/pgs-env.csh
# or else
# (2)
# source tkreset.csh
# This will remove previous toolkit-related paths
# It also snips duplicate paths
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
if ("$?PGSHOME" == 0) then
  # No previous toolkit paths
  echo "No previous PGSHOME defined"
  if (-e "$1") then
     source $1
  endif
  echo "new count for path is $#path"
  echo "new path is $path"
  echo "new PGSHOME is $PGSHOME"
else
  set oldpath = ('' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '')
  set count = 0
  if ($#path > 21) then
    echo "path too long to reset toolkit portion; sorry"
    exit 1
  endif
  foreach arg ($path[*])
    @ count++
    set oldpath[$count] = $path[$count]
  end
  # Now, if usage (1), will need to learn old, new PGSHOME values
  set oldhome = "$PGSHOME"
  set usage = 2
  if (-e "$1") then
       source $1
       set usage = 1
  endif
  # 
  set oldcount = 0
  set newcount = 0
  setenv PATH ""
  # echo "usage is $usage"
  # echo "count is $count"
  # echo "Looping over $oldpath[1-$count]"
  foreach arg ($oldpath[1-$count])
    # echo "arg is $arg"
    set repl = ""
    set repeat = ""
    set match = `echo $arg | /bin/grep "$oldhome"`
    # echo "match is $match"
    # Check for duplicates among 1st oldcount-1
    if ($oldcount > 0) then
      foreach oldarg ($oldpath[1-$oldcount])
        if ($arg == $oldarg) then
          set repeat = "yes"
        endif
      end
    endif
    #
    @ oldcount++
    if ($repeat == "yes") then
      echo "Snipping duplicate path $arg"
    else if ($usage == 2) then
      # Usage (2)
      if ("$match" == "") then
        @ newcount++
        if ("$newcount" == 1) then
          setenv PATH "$arg"
        else
          setenv PATH "${PATH}:$arg"
        endif
      endif
    else
      # Usage (1)
      @ newcount++
      if ("$match" == "") then
        # echo "newcount is $newcount"
        if ("$newcount" == 1) then
          setenv PATH "$arg"
        else
          setenv PATH "${PATH}:$arg"
        endif
      else
        set repl = `echo $arg | /bin/sed "s^$oldhome^$PGSHOME^"`
        if ("$newcount" == 1) then
          setenv PATH "$repl"
        else
          setenv PATH "${PATH}:$repl"
        endif
      endif
    endif
  end
  # echo "new count for path is $newcount"
  echo "new path is $path"
  # echo "with count $#path"
  if ($usage == 2) then
    unset PGSHOME
  else
    echo "new PGSHOME is $PGSHOME"
  endif
  # housekeeping (because we will source this script)
  unset arg
  unset count
  unset newcount
  unset oldcount
  unset match
  unset usage
  unset repl 
  unset repeat 
  unset oldhome
  unset oldpath
endif
# $Log$
# Revision 1.1  2003/05/06 20:40:52  pwagner
# First commit
#
