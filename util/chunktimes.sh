#!/bin/sh
#chunktimes.sh

# --------------- chunktimes.sh help
#Creates a table of chunk timings similar to the following
#chunk phase_1 phase_2 ..
#99 788.45 1288.13
#276 656.02 1326.15
#10 637.48 1267.51
#100 604.73 1393.60
#    .   .   .
#Usage:
#chunktimes.sh [opt1] [opt2] ..  mlsl2.log
#where mlsl2.log is the catenation of all the chunks' logs from a run
#
#    O p t i o n s
# -head list    creates table headed by comma-separated list of phases
#                 e.g., -head "core,corumba,corcorvado,papaya"
# -r            reverse sense of sort
# -sort k       where k is a number: sort by column k
#                 k=1: by chunk; k=2: by 1st phase; k=n: by (n-1)st phase
# -h[elp]       print brief help message; exit
#
#Note:
#The option(s) marked with "-", if present,
#must precede the file mlsl2.log on the command line
# --------------- End chunktimes.sh help
# Bugs and limitations
#(1) See chunktimes.pl
#(2) If sorting, must have write access to cwd
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
PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if willing to try patience of users
me="$0"
my_name=chunktimes.sh
I=chunktimes
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
# $the_perl_script is me but 'sh' -> 'pl'
the_perl_script="`echo $0 | sed 's/'$my_name'/chunktimes.pl/'`"
list=""
reverse="no"
sort="no"
s_column="0"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -head )
	    list="$2"
       shift
	    shift
       ;;
    -r )
	    shift
       reverse="yes"
       ;;
    -sort )
	    s_column="$2"
       shift
	    shift
       sort="yes"
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;

    -* )
	    k_column=`expr 0 - $2`
       shift
	    shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

if [ "$list" != "" ]
then
  extra_args="-head $list $1"
else
  extra_args="$1"
fi
temp_file1=`get_unique_name ct1`
temp_file2=`get_unique_name ct2`
if [ $PRINT_TOO_MUCH = "1" ]                            
then                                                    
   echo " Summary of options to $me "  
   echo " sort?: $sort "  
   echo " sort column: $s_column "  
   echo " column number: $k_column "  
   echo " list of phases: $list "  
   echo " temp_file1: $temp_file1 "  
   echo " temp_file2: $temp_file2 "  
   echo " extra args: $extra_args "  
   echo " perl script: $the_perl_script "  
fi                                                      
rm -f $temp_file1 $temp_file2
if [ "$sort" = "no" ]
then
  # echo $the_perl_script $extra_args
  $the_perl_script $extra_args > $temp_file1
  # $the_perl_script $extra_args | head
elif [ "$reverse" = "yes" ]
then
  $the_perl_script -headonly $extra_args > $temp_file1
  $the_perl_script -nohead $extra_args > $temp_file2
  sort -r +$s_column $temp_file2 >> $temp_file1
else
  $the_perl_script -headonly $extra_args > $temp_file1
  $the_perl_script -nohead $extra_args > $temp_file2
  sort +$s_column $temp_file2 >> $temp_file1
fi
cat $temp_file1
# head $temp_file1
rm -f $temp_file1 $temp_file2
exit
# $Log$
