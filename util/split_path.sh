#!/bin/sh
# split_path.sh
# --------------- split_path.sh help
# split arg (e.g., path/file_name) into path file_name
# Usage: split_path [option] arg
#
#    O p t i o n s
# -h            print help; exit
# -p            print path only
# -f            print file_name only
# -s char       separate path and file_name with char instead of space
#
# Notes:
# option must appear before arg
# options -s, -p and -f are mutually exclusive
# w/o option, prints path then file_name separated by a space
# so that in Bourne shell script
# arg_list=`split_path.sh path/file_name`
# for arg in $arg_list
#     do
#       echo $arg
#     done
# will print out
#     path
#     file_name
# * * *   s p e c i a l   c a s e s   * * *
# (in the following any opt matches {-f -p -s})
# if given arg ends with "/", then
#   `split_path.sh [any opt but -f] arg` returns arg, wile opt=-f returns null
# if given arg matches "/file_name", then
#   `split_path.sh [-f] arg` returns file_name
#   `split_path.sh [-p] arg` returns /
#   `split_path.sh arg` returns / file_name
# if given arg lacks "/" symbol (thus is not splittable), then
#   `split_path.sh [any opt] arg` returns arg no matter what opt
# --------------- End split_path.sh help

# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#---------------------------- nth_arg
#
# Function to return the nth (n > 0) arg of the std args after n
# eg given
# nth_arg 3 a b 'c or d' d e ..
# writes 'c or d' to standard output (w/o ' marks)

nth_arg()
{

   # Do we have enough args? No error if not, just blank output
     # echo "Entering nth_arg with args $@"
      the_arg=
      if [ $# -gt 1 ]
      then
         count_up=1
         n=$1
         shift
         while [ $count_up -lt $n ]
         do
           shift
           count_up=`expr $count_up + 1`
         done
         the_arg="$1"
      fi
     echo "$the_arg"
}
      
#----------------------- reverse -----------------------

# Function to echo reverse(arg)
#

reverse()
{
    perl -e '$parts=reverse("$ARGV[0]"); print $parts' "$1"
}


#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
me="$0"
my_name=split_path.sh
which_part="both"
sep=" "
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;
    -p )
        which_part="path"
	;;
    -f )
        which_part="file"
	;;
    -s )
        sep="$2"
        shift
	;;
    * )
       more_opts="no"
       arg="$1"
       ;;
    esac
    shift
done

# Does path contain "/" or not?

is_splittable=`echo "$arg" | grep '\/'`
if [ "$is_splittable" = "" ]
then
   echo $arg
   exit 0
fi

# Does path end with "/" or not?
# (in other words is it a pure path?)

is_pure_path=`echo "$arg" | grep '\/$'`
if [ "$is_pure_path" != "" ]
then
  if [ "$which_part" != "file" ]
  then
     echo $arg
  fi
  exit 0
fi

# Does path start with "/" or not?
# (in other words is it absolute or relative?)

is_absolute=`echo "$arg" | grep '^\/'`

#echo "arg is $arg"
gra=`reverse $arg`
#echo "reverse arg is $gra"
the_strap=`echo $gra | sed 's;/; ;'`
#echo "reverse the_parts is $the_strap"
the_parts=`reverse "$the_strap"`
#echo "the_parts is $the_parts"
#for arg in $the_parts
#do
#  echo $arg
#done

path=`nth_arg 1 $the_parts`
file_name=`nth_arg 2 $the_parts`

if [ "$file_name" = "" -a "$is_absolute" != "" ]
then
  # format of arg was /file_name
  # so path = "/" and file_name = file_name
  file_name="$path"
  path="/"
fi

if [ "$which_part" = "path" ]
then
  echo $path
elif [ "$which_part" = "file" ]
then
  echo $file_name
else
  echo $path${sep}$file_name
fi
exit 0

# $Log$
