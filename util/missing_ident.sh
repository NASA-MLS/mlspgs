#!/bin/sh
# missing_ident.sh

# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
# --------------- missing_ident.sh help
#Print list of files lacking $id and $RCS lines
#List compiled from args on command-line
#or else automatically computed to match pattern
#  *{$suffix} in supplied directory names
#
#Usage:
#missing_ident.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#
#    O p t i o n s
# -d path       path of directory if automatic
# -suf suffix   file name suffix if automatic; e.g. ".f90"
# -h[elp]       print brief help message; exit
#Notes:
#(1)If option -suf present but -d absent, will search current working directory
#(2)The options -suf and -d may be repeated; e.g.
#   missing_ident -d dir1 -d dir2 .. -suf sufx1 -suf sufx2 ..
# --------------- End missing_ident.sh help
#------------------------------- extant_files ------------
#
# Function to return only those files among the args
# that actually exist
# Useful when passed something like *.f which may 
# (1) expand to list of files, returned as extant_files_result, or
# (2) stay *.f, in which case a blank is returned as extant_files_result 
#     (unless you have perversely named a file '*.f')
# usage: extant_files arg1 [arg2] ..

extant_files()
{
   extant_files_result=
   # Trivial case ($# = 0)
   if [ "$1" != "" ]
   then
      for file
      do
         if [ -f "$file" ]
         then
               extant_files_result="$extant_files_result $file"
         fi
      done
   fi
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
id_list=
RCS_list=
echo "Who or what is $@"
suffixes=""
directories=""
me="$0"
my_name=missing_ident.sh
#
# Get arguments from command line
#

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    case "$1" in
	-d )
	    directories="$directories $2"
	    shift
       shift
	;;
	-suf )
	    suffixes="$suffixes $2"
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
extant_files "$@"
arglist="$extant_files_result"
if [ "$suffixes" != "" ]
then
  # Assemble automatic list of files
  if [ "$directories" = "" ]
  then
    directories="."
  fi
  args_dir=`pwd`
  for dir in $directories
  do
    cd $dir
    for suffix in $suffixes
    do
      extant_files *"$suffix"
      if [ "$extant_files_result" != "" ]
      then
        for file in $extant_files_result
        do
          arglist="$arglist $dir/$file"
        done
      fi
    done
    cd $args_dir
  done
fi
for file in $arglist
do
   # echo file is $file
   test=`grep -i '\$id' $file`
      if [ "$test" = "" ]
      then
         id_list="$id_list $file"
      fi
   test=`grep -i '\$RCSFI' $file`
      if [ "$test" = "" ]
      then
         RCS_list="$RCS_list $file"
      fi
   
done
echo "These files lack a \$id line"
echo "$id_list"
echo "These files lack a \$RCS line"
echo "$RCS_list"

exit 0
# $Log$
