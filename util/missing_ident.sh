#!/bin/sh
# missing_ident.sh

# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
# --------------- missing_ident.sh help
#Use(1) (default)
#Print list of files lacking $id and $RCS lines
#List compiled from args on command-line
#or else automatically computed to match pattern
#  *{$suffix} in supplied directory names
#
#Use(2)
#Print list of files whose names do not appear among idents actually
#incorporated into executable
#Usage:
#missing_ident.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#
#    O p t i o n s
# -d path       path of directory if automatic
# -suf suffix   file name suffix if automatic; e.g. ".f90"
# -h[elp]       print brief help message; exit
# -x executable the executable, switching to use(2)
#Notes:
#(1)If option -suf present but -d absent, will search current working directory
#(2)The options -suf and -d may be repeated; e.g.
#   missing_ident -d dir1 -d dir2 .. -suf sufx1 -suf sufx2 ..
# --------------- End missing_ident.sh help
#Another example of my regrettable tendency to bundle a number of
#different functionalities within a single multi-use file
#Isn't it cleaner to devote a single file to a single use?

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
# echo "All the args: $@"
suffixes=""
directories=""
me="$0"
my_name=missing_ident.sh
ex_name=
# $the_splitter is split_path with me's path (if any) prepended
the_splitter="`echo $0 | sed 's/missing_ident/split_path/'`"
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
	-x )
	    ex_name="$2"
	    shift
       shift
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

temp1=`get_unique_name id`
temp2=`get_unique_name rcs`
if [ "$ex_name" != "" ]
then
  rm -f "$temp1"
  rm -f "$temp2"
  ident "$ex_name" | grep -i '\$id' > "$temp1"
  ident "$ex_name" | grep -i '\$rcs' > "$temp2"
  #ls -tl "$temp1" "$temp2"
fi

for file in $arglist
do
   # echo file is $file
   file_name=`$the_splitter -f $file`
   if [ "$ex_name" != "" ]
   then
     test=`grep "$file_name" "$temp1"`
   else
     test=`grep -i '\$id' $file`
   fi
   #echo "id test"
   #echo "file_name $file_name"
   #echo "test $test"
   if [ "$test" = "" ]           
   then                          
      id_list="$id_list $file"   
   fi                            
   if [ "$ex_name" != "" ]
   then
     test=`grep "$file_name" "$temp2"`
   else
     test=`grep -i '\$RCSFI' $file`
   fi
   #echo "rcs test"
   #echo "file_name $file_name"
   #echo "test $test"
   if [ "$test" = "" ]             
   then                            
      RCS_list="$RCS_list $file"   
   fi                              
   
done
if [ "$ex_name" != "" ]   
then
  echo "These files do not appear as \$id lines in $ex_name"
  echo "$id_list"
  echo "These files do not appear as \$RCS lines in $ex_name"
  echo "$RCS_list"
  rm -f "$temp1" "$temp2"
else                      
  echo "These files lack a \$id line"
  echo "$id_list"
  echo "These files lack a \$RCS line"
  echo "$RCS_list"
fi

exit 0
# $Log$
# Revision 1.1  2002/10/08 16:57:33  pwagner
# First commit
#
