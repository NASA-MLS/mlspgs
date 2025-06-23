#!/bin/sh
# mkpath.sh

# Usage: mkpath.sh path
# mkdir each step or subpath along a path

# given a path like p1/p2/p3/../pn
# mkdir p1 unless true =  [ -d p1 ]
# mkdir p2 unless true =  [ -d p2 ]
# mkdir p3 unless true =  [ -d p3 ]
#   . . .
# mkdir pn unless true =  [ -d pn ]
# i.e., don't leave before p1/p2/p3/../pn is a valid path

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

if [ $# -lt 1 ]
then
   echo "Usage:  mkpath.sh path"
   exit 1
fi

# Does path start with "/" or not?
# (in other words is it absolute or relative?)

is_absolute=`echo "$1" | grep '^\/'`

the_pn=`echo $1 | sed 's;/; ;g'`

if [ "$is_absolute" = "" ]
then
   the_path=
else
   the_path="/"
fi

for p in $the_pn
do 
   the_path="${the_path}$p/"
   if [ ! -d "$the_path" ]
   then
      mkdir "$the_path"
#      echo mkdir "$the_path"
   fi
done

exit

# $Log$
# Revision 1.1  2003/06/03 20:36:58  pwagner
# First commit under this name
#
