#!/bin/sh
# mkdir_cascade.sh

# Usage: mkdir_cascade.sh path
# mkdir each step or subpath along a path

# given a path like p1/p2/p3/../pn
# mkdir p1 unless true =  [ -d p1 ]
# mkdir p2 unless true =  [ -d p2 ]
# mkdir p3 unless true =  [ -d p3 ]
#   . . .
# mkdir pn unless true =  [ -d pn ]
# i.e., don't leave before p1/p2/p3/../pn is a valid path

# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

if [ $# -lt 1 ]
then
   echo "Usage:  mkdir_cascade.sh path"
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
