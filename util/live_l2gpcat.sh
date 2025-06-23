#!/bin/sh
# --------------- live_l2gpcat.sh help
# live_l2gpcat.sh
# a wrapper for l2gpcat
# assuring that the files to l2gpcatted are "pure"
# i.e., not in the middle of being written to
#
# live_l2gpcat.sh [opt1] [opt2] ..  -o output-file file1 [file2] ..
# Special notes
# (1) l2gpcat must be in your PATH
# (2) h5debug must be in your PATH
# (3) split_path.sh must have the same path prefix as live_l2gpcat.sh
# (4) A sub-directory, live_l2gpcat, will be created in your home directory
#     to hold copies of the files being copied; it will be removed afterwards
# i.e.,  
# BrO CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI Temperature

#     O p t i o n s
#    -dryrun              Merely echo the commands that would be executed
#    -l2gpcat command     Use command instead of l2gpcat
#    -h5debug command     Use command instead of h5debug

# Bugs and limitations
# (1) All the input files are assumed to exist and, potentially, be valid hdf5
# (2) We can create ${HOME}/live_l2gpcat

# usage: see (1) above

# --------------- End live_l2gpcat.sh help
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

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
   echo $extant_files_result
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
#
debug=0
#     ^  -- set this to 1 if worried
cleanup=1
#       ^  -- set this to 1 to clean up afterwards
echo $@
me="$0"
my_name=live_l2gpcat.sh
I=live_l2gpcat
#L2GPCAT=$HOME/bin/l2gpcat
#H5DEBUG=$HOME/bin/h5debug
L2GPCAT=l2gpcat
H5DEBUG=h5debug
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
l2gpcat_opts=""
list=""
live=""
dryrun="no"
outfile=""
live_dir="$HOME/live_l2gpcat"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -o )
	    outfile="$2"
	    shift
	    shift
       ;;
    -l2gpcat )
	    L2GPCAT="$2"
	    shift
	    shift
       ;;
    -h5debug )
	    H5DEBUG="$2"
	    shift
	    shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;
    -* )
       l2gpcat_opts="$l2gpcat_opts $1"
	    shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

echo $@
list=`extant_files $@`

if [ "$debug" = 1 ]
then
  echo "l2gpcat $L2GPCAT"
  echo "h5debug $H5DEBUG"
  echo "l2gpcat_opts $l2gpcat_opts"
  echo "list $list"
fi

if [ "$list" = "" ]
then
  sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
      | sed -n 's/^.//p' | sed '1 d; $ d'
  exit
fi

if [ -d "$live_dir" ]
then
  echo "Error--you must be able to create $live_dir"
  exit 1
fi

mkdir "$live_dir"

for i in $list
do
  name=`$the_splitter $i`
  returnstatus=1
  while [ "$returnstatus" != 0 ]
  do
    if [ "$dryrun" = "yes" ]
    then
      echo cp "$i" "$live_dir"
      echo $H5DEBUG "$live_dir/$name"
    else
      cp "$i" "$live_dir"
      $H5DEBUG "$live_dir/$name" > /dev/null
    fi
    returnstatus="$?"
    if [ "$returnstatus" != 0 ]
    then
      echo "Retry with $i"
      sleep 10
    fi
  done
  live="$live $live_dir/$name"
done
if [ "$dryrun" = "yes" ]
then
  echo $L2GPCAT $l2gpcat_opts -o $outfile $live
else
  $L2GPCAT $l2gpcat_opts -o $outfile $live
  returnstatus="$?"
fi
if [ "$cleanup" = 1 ]
then
  /bin/rm $live_dir/*
  rmdir $live_dir
fi
exit 0

# $Log$
