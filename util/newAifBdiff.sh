#!/bin/sh
#newAifBdiff.sh

# --------------- newAifBdiff.sh help
# Executes supplied command 
# which presumably creates a file newA & possibly newB;
# if new B same as old B, keeps old A;
# else deletes old A and B, replacing with new A and B
# (works even if A and B name the same files
#   meaning simply keep old A, and its date, unless new A different)
# leaves record file that names A and tells whether it is new or old
# if record file already exists, appends new line to it with above info
#Usage:
#newAifBdiff.sh [options] old_A old_B command [arg1] [arg2] ..
#
#    O p t i o n s
# -a            no B; usage newAifBdiff.sh -a old_A command [args] ..
#                 meaning keep old A (and its date) unless new A different
# -k            if keeping old A, keep old B, too
#                (as no diff betw. old and new B, means keep date of old B)
# -nr           leave no record file
# -x n          ignore actual error status from command: exit with "n"
# -h[elp]       print brief help message; exit
# old_A         the file to be conditionally replaced
# old_B         the file to be compared with new_B
# command       the command to be executed (producing new_A and _B)
# arg1 ..       optional arguments passed to command
#Result:
#file(s) named the same as old_A and old_B
#record file named newAifBdiff.out
#and anything else done by command as a side-effect
#Useful as an adjunct to Makefile for conditionally
#compiling according to whether .mod file of prerequisite has changed;
#this reduces the cascade of massive recompilations that results
#when the prerequisite changes
#
#Notes:
#(1) The option(s) marked with "-", if present,
#    must precede the command and args on the command line
#(2) Because the option, "-a", is legal to gnu diff, but not to Sun's diff
#    we have ceased relying on it; instead we use the octal dump routine od
#    if your platform/os doesn't support it please let me know
#(3) Unless you specify otherwise, using the -x n option,
#    the script will attempt to pass the exit status of command back to
#    whoever called it
#
#Record file
#file name: newAifBdiff.out
#format: ascii text, two columns, 1st is file names, 2nd new/old
#  file_name_of_1st_A   new
#  file_name_of_2nd_A   old
#      .  .  .
#  file_name_of_this_A  new
# --------------- End newAifBdiff.sh help
#Bugs and limitations:
#(1) Why have two ways to do essentially the same thing: newAifAdiff.sh?
#(2) Why such a lousy name?
#(3) How about optionally using a user-supplied program instead of diff?
#      that way we could get rid of relying on octal dump routine
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#------------------------------- what_diff_opt ------------
#
# Function to determine whether a diff_opt of "-a" is
# warranted by the type of file passed as argument
# i.e., if ascii, then diff_opt is ""
# but if data, then "-a"
# usage: what_diff_opt arg1
# returns result as diff_opt
# Superseded: see diff_fun

what_diff_opt()
{
   diff_opt=
   # Trivial case ($# = 0)
   if [ "$1" != ""  -a -f "$1" ]
   then
      file_type=`file "$1" | sed 's/^.*://'`
#      echo $file_type
      reduced_file_type=`echo $file_type | grep -i text`
      if [ "$reduced_file_type" = "" ]
      then
         diff_opt="-a"
      fi
   fi
}

#------------------------------- diff_fun ------------
#
# Function to determine whether two files are different
# w/o relying on Gnu's diff with its handy "-a" option
# (Inferior standard versions of diff, like Sun's, lack this option)
# usage: extant_files arg1 [arg2] ..
# returns (number) result as the_diff: 0 if none, else non-zero

diff_fun()
{
   the_diff="0"
	if [ $# -gt "1" ]
	then
      rm -f temp1 temp2
      od -t c "$1" > temp1
      od -t c "$2" > temp2
      the_diff=`diff temp1 temp2 | wc -l`
      rm -f temp1 temp2
   fi
}

#
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
# The initial settings are
# A, B same files              no  
# old A exists                 yes 
# old B exists                 yes 
# keeping old unchanged B      no
# keeping records of A         yes
# diff_opt                     -a
ASameAsB="no"
oldAExists="yes"
oldBExists="yes"
keepUnchangedB="no"
keepRecords="yes"
me="$0"
my_name=newAifBdiff.sh
record_file=newAifBdiff.out
the_status="undefined"
DEEBUG="no"
NORMAL_STATUS=0
MUST_EXIT_STATUS=""
return_status=0

wrong_list=""
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;
    -a )
       ASameAsB="yes"
       shift
	;;
    -k )
       keepUnchangedB="yes"
       shift
	;;
    -nr )
       keepRecords="no"
       shift
	;;
    -x )
       shift
       MUST_EXIT_STATUS=`expr $1`
       shift
	;;
    * )
       more_opts="no"
       ;;
    esac
done

old_A="$1"
shift
if [ "$ASameAsB" = "no" ]; then
  old_B="$1"
  shift
else
  old_B="$old_A"
fi
the_command="$1"
shift
if [ "$old_A" = "$old_B" ] ; then
  ASameAsB="yes"
fi
if [ ! -f "$old_A" ] ; then
  oldAExists="no"
fi
if [ ! -f "$old_B" ] ; then
  oldBExists="no"
fi
if [ "$DEEBUG" = "yes" ]; then
  echo "old_A: $old_A"
  echo "old_B: $old_B"
  echo "MUST_EXIT_STATUS: $MUST_EXIT_STATUS"
  echo "the_command: $the_command"
  echo "the_args: $@"
  echo "ASameAsB: $ASameAsB"
  echo "oldAExists: $oldAExists"
  echo "oldBExists: $oldBExists"
fi

if [ "$keepRecords" = "yes" -a ! -f "$record_file" ]; then    
  echo "#filename   status" > "$record_file"
fi
if [ "$oldBExists" = "no" ]; then
# There's no old B => automatically new and old B differ, so need new A
    message="There is no old B => automatically need new A"
  "$the_command" "$@"
  return_status=`expr $?`
  the_status="new"
elif [ "$ASameAsB" = "yes" ]; then
# A and B the same, so check if new A diff from old
  mv "$old_A" "$old_A.1"
  "$the_command" "$@"
  return_status=`expr $?`
#  the_diff=`diff -a $old_A $old_A.1 | wc -l`
#  what_diff_opt $old_A
#  echo diff $diff_opt $old_A $old_A.1 | wc -l
#  the_diff=`diff $diff_opt $old_A $old_A.1 | wc -l`
  diff_fun $old_A $old_A.1
  if [ "$the_diff" -gt 0 ] ; then
    rm -f "$old_A.1"
      message="($the_diff) A, B same; old A, new A differ =>need new A"
      the_status="new"
  else
    mv "$old_A.1" "$old_A"
      message="($the_diff) A, B same; old A, new A same =>keep old A"
      the_status="old"
  fi
elif [ "$oldAExists" = "no" ]; then
# There's no old A so need new A
  "$the_command" "$@"
  return_status=`expr $?`
    message="There is no old A => automatically need new A"
    the_status="new"
else
# Most general case
  mv "$old_A" "$old_A.1"
  mv "$old_B" "$old_B.1"
  "$the_command" "$@"
  return_status=`expr $?`
#  the_diff=`diff -a $old_B $old_B.1 | wc -l`
#  what_diff_opt $old_B
#  echo diff $diff_opt $old_B $old_B.1 | wc -l
#  the_diff=`diff $diff_opt $old_B $old_B.1 | wc -l`
  diff_fun $old_B $old_B.1
  if [ "$the_diff" -gt 0 ] ; then
    rm -f "$old_A.1" "$old_B.1"
      message="($the_diff) old B, new B differ =>need new A"
     the_status="new"
  else
    mv "$old_A.1" "$old_A"
      message="($the_diff) old B, new B same =>keep old A"
      the_status="old"
    if [ "$keepUnchangedB" = "yes" ] ; then
      mv "$old_B.1" "$old_B"
    else
      rm -f "$old_B.1"
    fi
  fi
fi
if [ "$DEEBUG" = "yes" ]; then    
  echo "$message"                 
fi                                
if [ "$keepRecords" = "yes" ]; then    
  echo "$old_A   $the_status" >> "$record_file"
fi

# What status must we exit with?
#echo "MUST_EXIT_STATUS: $MUST_EXIT_STATUS"
#echo "return_status: $return_status"
#echo "NORMAL_STATUS: $NORMAL_STATUS"
if [ "$MUST_EXIT_STATUS" != "" ]; then    
#   echo "exiting with status $MUST_EXIT_STATUS (special)"
   exit "$MUST_EXIT_STATUS"
elif [ "$return_status" != "$NORMAL_STATUS" ]; then
#   echo "exiting with status 1"
   exit 1
else
#   echo "exiting with status 0"
   exit 0
fi
# $Log$
# Revision 1.4  2002/06/25 18:04:17  pwagner
# Relies on octal dump routine instead of -a option to diff
#
# Revision 1.3  2002/06/21 20:32:06  pwagner
# Conditionally passes option -a to diff
#
# Revision 1.2  2002/06/21 00:09:32  pwagner
# Added -nr option
#
# Revision 1.1  2002/05/22 00:36:28  pwagner
# First commit
#
