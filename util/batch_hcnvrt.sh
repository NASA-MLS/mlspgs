#!/bin/sh
# batch_hcnvrt.sh
# --------------- batch_hcnvrt.sh help
#batch convert l2gp files from version 2 hdfeos to version 5.
#By default uses heconvert (for which see below).
#
#Usage:
#batch_hcnvrt.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [file1] .. [filen]
#
#Result:
#heconverts list of files made up of 
#(a) all files matching the calculated pattern in the source directory,
#(b) any extra arguments file1 .. filen that follow the last option 
#The resulting files will go in the destination directory
#  C a l c u l a t e d   p a t t e r n s
#Given character strings for the prefix and suffix,
#the pattern to be matched is [prefix]*[suffix]
#the resulting file names will be the same as source file names
#unless differing pre(suf)fixes are supplied via -npre(suf) options
#
#    O p t i o n s
# -c command    name of command if different from heconvert
# -s path       path of source directory
# -d path       path of destination directory (created if not already exists)
# -pre prefix   character string prefix in calculated pattern
# -suf suffix   character string suffix in calculated pattern
# -npre prefix  prefix for the resulting file if different from source
# -nsuf suffix  suffix for the resulting file if different from source
# -o output     send stdin and stderr to file output (defaults to your screen)
# -v            print verbose messages while heconverting
# -h[elp]       print brief help message; exit
#Notes:
#(1)The options -s and -d are mandatory (except if lone -h option)
#(2)All options except for -v and -h require a following arg
#(3)If -pre(suf) omitted, only extra arguments will be converted
#Thus, even if no files are found in the source directory matching the
#calculated pattern, or if no options -pre(suf) are submitted, a list
#of files to be converted can be entered directly
# 
#  h e c o n v e r t
#heconvert is an executable binary that must be built and installed
#Try the command 'which heconvert' to see if heconvert exists
#and is within your execution PATH
#If not, the command 'make heconvert' will build and install it
# --------------- End batch_hcnvrt.sh help
#
# Function to mkdir each step or subpath along a path
#

mkdir_cascade()
{

   # Does path start with "/" or not?
   # (in other words is it absolute or relative?)

   is_absolute=`echo "$*" | grep '^\/'`

   the_pn=`echo $* | sed 's;/; ;g'`

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

}

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

DEBUG=0
#     ^  -- set this to 1 to echo rather than execute all commands
if [ "$DEBUG" = "1" ]
then
  run_mode="echo"
else
  run_mode=""
fi

PRINT_TOO_MUCH=1
#              ^  -- set this to 1 if patient enough to print more
if [ "$PRINT_TOO_MUCH" = "1" ]
then
  echo "$0 called with args $@"
fi

#
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# variable        meaning            default
# command     name of command       "heconvert"
# prefix   in calculated pattern    ""
# suffix   in calculated pattern    ""
# spath      source directory       "wrong"
# dpath    destination directory    "wrong"
# me                                "$0"
# my_name                           batch_hcnvrt.sh
# result_pattern   same or diff     "same"
# vopt            -v or blank       ""
# output    output file or blank    ""
# arglist  list of source files     ""
#----------------------- Implementation -----------------------

command="heconvert"
prefix=""
suffix=""
result_pattern="same"
spath="wrong"
dpath="wrong"
me="$0"
my_name=batch_hcnvrt.sh
vopt=""
output=""
arglist=""
#
# Get arguments from command line
#

while [ "$1" != "" ] ; do

    case "$1" in
	-c )
	    command=$2
	    shift
	;;
	-s )
	    spath=$2
	    shift
	;;
	-d )
	    dpath=$2
	    shift
	;;
	-pre )
	    prefix=$2
	    shift
	;;
	-suf )
	    suffix=$2
	    shift
	;;
	-npre )
	    nprefix=$2
       result_pattern="diff"
	    shift
	;;
	-nsuf )
	    nsuffix=$2
       result_pattern="diff"
	    shift
	;;
	-o )
	    output=$2
	    shift
	;;
	-v )
	    vopt="-v"
	;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;
	* )
      arglist="$arglist $1"
	;;
    esac
    shift

done

# Error checks

# Mandatory options missing?
if [ "$spath" = "wrong" -o "$dpath" = "wrong" ]
then
   echo "You must set both options -s and -d"
   exit 1
fi

# Does the destination directory exist?
if [ ! -d "$dpath" ]
then
   mkdir_cascade $dpath
fi
# Do we have write permission to it?
if [ ! -w "$dpath" ]
then
   echo "Cannot write to destination path $dpath"
   exit 1
fi

args_dir=`pwd`
cd "$spath"
if [ "$prefix" = "" -a "$suffix" = "" ]
then
   echo "No calculated pattern matching; converting only supplied list"
else
   extant_files "$prefix"*"$suffix"
   arglist="$arglist $extant_files_result"
fi

if [ "$arglist" = "" ]
then
   echo "No files in $spath match the calculated pattern or supplied list"
   exit 1
fi

if [ "$PRINT_TOO_MUCH" = "1" ]
then
  echo "list of files: $arglist"
fi

cd $args_dir
if [ "$output" != "" ]      
then     
  rm -f $output
  rm -f "$output".tmp
  echo "Output from batch_hcnvrt.sh" > "$output"
fi              

for file in $arglist
do
   if [ -r $spath/$file ]
   then

#     Target name
      if [ "$result_pattern" = "same" ]
      then
         target=$file
      elif [ "$prefix" = "" ] 
      then
          target="`echo $file | sed 's/'$suffix'/'$nsuffix'/'`"
      elif [ "$suffix" = "" ] 
      then
          target="`echo $file | sed 's/'$prefix'/'$nprefix'/'`"
      else
          target="`echo $file | sed 's/'$prefix'/'$nprefix'/'`"
          target="`echo $target | sed 's/'$suffix'/'$nsuffix'/'`"
      fi

      if [ "$output" = "" ]
      then
        $run_mode $command "$vopt" -i $spath/$file -o $dpath/$target 
      else
        $run_mode $command "$vopt" -i $spath/$file -o $dpath/$target > "$output".tmp 2>&1
        cat "$output".tmp >> "$output"
        rm -f "$output".tmp
      fi
   fi
done

exit 0

# $Log$
