#!/bin/sh
#batch_l2gpdiff.sh

# --------------- batch_l2gpdiff.sh help
# Performs a l2gpdiff on files from two separate runs for various species
# saving the results to one file per species
# Usage:
# batch_l2gpdiff.sh [opt1] [opt2] ..  str1 [str2] .. - dir1 dir2
# where each string is a valid species name; e.g., H2O
# and dir1, dir2 are the directories holding the l2gp files from the two runs
#
#     O p t i o n s
#    -append              Append the results to an existing file
#                          (Otherwise replace any existing file)
#    -unique              Save results to a uniquely-named new file
#    -l2gpdiff command    Use command instead of l2gpdiff
#  (and all the normal l2gpdiff options; e.g. -rms -ignore are wise choices)
#
# Note:
# (1) Each option is preceded by a dash
# (2) The options must come before the strings
# (3) The bare "-", must be present 
#      (to separate the strings from the directories on the command line)
# (4) The files being l2gpdiffed must be in different directories 
# (5) Each comparison will be saved in the current working directory as str.diff
#
# Bugs and limitations
# (1) l2gpdiff is assumed to exist, to be an executable, and to be in
#     $HOME/mlspgs/bin/LF95.Linux
# (2) The l2gp file names are assumed to match the pattern
#      MLS-Aura_L2GP-xxxx_*.he5
# (3) It would be convenient to have some helpful short string that
#      is translated into all the standard species
# (4) Why don't you just make "-rms -ignore" the default options, either
#      here or in l2gpdiff itself?
# --------------- End batch_l2gpdiff.sh help
# Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in pace of "."

get_unique_name()
{

   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
      fi
   # Is $HOST defined?
      if [ "$HOST" != "" ]
      then
         our_host_name="$HOST"
      elif [ "$HOSTNAME" != "" ]
      then
         our_host_name="$HOSTNAME"
      else
         our_host_name="host"
      fi
    #  echo $our_host_name
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
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
#
debug=1
#     ^  -- set this to 1 if worried
PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience of users
me="$0"
my_name=batch_l2gpdiff.sh
I=batch_l2gpdiff
L2GPDIFF=~/mlspgs/bin/LF95.Linux/l2gpdiff
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
l2gpdiff_opts=""
list=""
append="no"
unique="no"
more_opts="yes"
more_strs="yes"
while [ "$more_strs" = "yes" ] ; do

    case "$1" in

    - )
	    shift
       more_opts="no"
       more_strs="no"
       ;;
    -append )
	    append="yes"
	    shift
       ;;
    -unique )
	    unique="yes"
	    shift
       ;;
    -l2gpdiff )
	    L2GPDIFF="$2"
	    shift
	    shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;
    -* )
       if [ "$l2gpdiff_opts" = "" ]
       then
	       l2gpdiff_opts="$1"
       else
	       l2gpdiff_opts="$l2gpdiff_opts $1"
       fi
       if [ "$1" = "-d" -o "$1" = "-f" ]
       then
	       l2gpdiff_opts="$l2gpdiff_opts $2"
          shift
       fi
	    shift
       ;;
    * )
       more_opts="no"
       if [ "$list" != "" ]
       then
         list="$list $1"
       else
         list="$1"
       fi
	    shift
       ;;
    esac
done

if [ "$list" = "" ]
then
  sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
      | sed -n 's/^.//p' | sed '1 d; $ d'
  exit
fi
dir1="$1"
dir2="$2"
if [ "$debug" = 1 ]
then
  echo "L2GPDIFF $L2GPDIFF"
  echo "l2gpdiff_opts $l2gpdiff_opts"
  echo "list $list"
  echo "unique $unique"
  echo "append $append"
  echo "dir1 $dir1"
  echo "dir2 $dir2"
fi
if [ "$dir1" = "$dir2" ]
then
  echo "Sorry--dir1 and dir2 must be different"
elif [ ! -d "$dir1" ]
then
  echo "Sorry--dir1 must be a directory"
elif [ ! -d "$dir2" ]
then
  echo "Sorry--dir2 must also be a directory"
elif [ ! -x "$L2GPDIFF" ]
then
  echo "Sorry--L2GPDIFF must exist and be executable"
fi
for i in $list
do
  extant_files $dir1/MLS-Aura_L2GP-${i}_*.he5
  file1="$extant_files_result"
  extant_files $dir2/MLS-Aura_L2GP-${i}_*.he5
  file2="$extant_files_result"
  if [ "$debug" = 1 ]
  then
    echo "species $i"
    echo "file1 $file1"
    echo "file2 $file2"
  fi
  if [ "$unique" = "yes" ]
  then
     outfile=`get_unique_name $i`
  else
    outfile="$i.diff"
  fi
  if [ ! -f "$outfile" ]
  then
    echo "$file1 $file2 $outfile" > "$outfile"
  elif [ "$append" = "yes" ]
  then
     echo "$file1 $file2 $outfile" >> "$outfile"
  else
    rm -f "$outfile"
    echo "$file1 $file2 $outfile" > "$outfile"
  fi
  if [ "$file1" = "" ]
  then
    echo "No $i file found in $dir1"
  elif [ "$file2" = "" ]
  then
    echo "No $i file found in $dir2"
  elif [ "$l2gpdiff_opts" = "" ]
  then
    $L2GPDIFF $file1 $file2 >> "$outfile"
  else
    $L2GPDIFF $l2gpdiff_opts $file1 $file2 >> "$outfile"
  fi
done
exit
# $Log$
# Revision 1.1  2004/12/29 22:26:03  pwagner
# First commit
#
