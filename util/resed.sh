#!/bin/sh
#resed.sh

# --------------- resed.sh help
# Runs the command line arguments through sed, replacing
# (depending on option(s) selected)
# them with the changes sed brings about
# Useful for automating scripts that modify text files
# (1) comment/uncomment lines
# (2) add/remove line(s)
# (3) change hard-coded paths (e.g. in perl scripts)
#
# Usage:
# resed.sh [opt] ..  file1 [file2 ..]
#
#    O p t i o n s
# -o options    pass options to sed (default is none)
# -c command    command to pass (surrounded by 's) to sed
# -f file       file of commands to pass through to sed
# -h[elp]       print brief help message; exit
# -n            rename new file "old_name"+suffix, keeping old as old_name
# -o            rename new file "old_name", keeping old as old_name+suffix
# -suffix=xxx   suffix to use when renaming either new or old file
#               (w/o any separator; e.g., to add ".bak" use -suffix=.bak)
# filen         file name (with path)
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede any filen on the command line
# (2) The options -c and -f are mutually exclusive
#     i.e., you must pass commands to sed either by file or by command-line
# (3) The options -n and -o are mutually exclusive
#     if you choose -n but omit suffix, the new files will be "old_name.new"
#     if you choose -o but omit suffix, the old files will be "old_name.old"
#     if you choose neither the new file will simply replace the old one
#     (leaving you with no backup)
# 
# Result:
# Files in the list will be modified according to
# the conditions set up by the options
# --------------- End resed.sh help
# Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
#
#   Notes
# Bugs: 
# (1) If the old_name+suffix exceeds the longest possible file name
#     the script does not fail gracefully
# (2) If filen+suffix is the same as filem for some pair (n,m)
#     one will replace the other w/o any warning of this possibility
# Unimplemented improvements: 
# (1) Why not allow the user to input original and replacement strings
#     via -os "string1" -rs "string2"
# (2) Also allow inserting a block of text stored in a file after  
#     line nnn via -l nnn -b block_file
# (3) Instead of messing with suffixes, let the script edit files
#     from one dir, saving the modified versions in another
#     via -d1 d_orig -d2 d_mod
#     (you will have to disable reecho part if you do this)
me="$0"
my_name=resed.sh
I=resed
NORMAL_STATUS=0
return_status=0
# $unique_name is unique_name with me's path prepended
unique_name="`echo $0 | sed 's/'$I'/unique_name/'`"
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
DEEBUG=on
NO_REPLACE=off
if [ $DEEBUG = "on" ]
then
   echo "Called me as $0"
   echo "with args $@"
fi

command_file=""
the_command=""
the_opt=""
more_opts="yes"
rename_which="neither"
the_suffix=""
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -o )
       the_opt="$2"
       shift
       shift
       ;;
    -c )
       the_command="$2"
       shift
       shift
       ;;
    -f )
       command_file="$2"
       shift
       shift
       ;;
    -n )
       rename_which="new"
       shift
       ;;
    -o )
       rename_which="old"
       shift
       ;;
    -suffix=* )
       the_suffix=`echo $1 | sed 's/-suffix=//'`
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

if [ $DEEBUG = "on" ]
then
   echo "the_opt $the_opt"
   echo "command_file $command_file"
   echo "the_command $the_command"
   echo "the_suffix $the_suffix"
   echo "rename_which? $rename_which"
   echo "remaining args $@"
fi
# Check for forbidden arguments
if [ "$the_command" != "" -a "$command_file" != "" ]
then
  echo 'You cannot have both a command file and command-line'
elif [ "$command_file" != "" ]
then
  the_opt="$the_opt -f $command_file"
elif [ "$the_command" = "" ]
then
#
  echo 'You must have either a command file or a command-line'
fi
if [ "$rename_which" = "new" -a "$the_suffix" = "" ]
then
  the_suffix=.new
elif [ "$rename_which" = "old" -a "$the_suffix" = "" ]
then
  the_suffix=.old
fi

the_list=`$reecho "$@"`
if [ "$the_list" = "" ]
then
  echo 'No valid files were found among the command line arguments'
  echo 'Make sure you entered their path/names correctly'
fi

for file in $the_list                                      
do
  file_name=`$the_splitter -f $file`
  if [ "$file" = "$file_name" ]
  then
    # No path supplied
    path=""
  else
    path=`$the_splitter -p $file`
  fi
  temp_name=`$unique_name resed`
  if [ "$path" != "" ]
  then
    temp_name="$path"/"$temp_name"
  fi
  if [ "$the_command" != "" ]
  then
    # echo sed $the_opt "'"$the_command"'" $file
    sed $the_opt "$the_command" $file > "$temp_name"
  else
    # echo sed $the_opt $file
    sed $the_opt $file > "$temp_name"
  fi
  # If sed failed, give up right now
  return_status=`expr $?`
  if [ "$return_status" != "$NORMAL_STATUS" ]; then       
     echo "Sorry--sed returned an error"   
    if [ "$the_command" != "" ]
    then
     echo "Possibly a syntax error in your command: $the_command"   
    else
     echo "Possibly a syntax error in your command-file: $command_file"   
    fi
    rm "$temp_name"
    exit 1                                               
  elif [ "$NO_REPLACE" = "on" ]
  then
    exit 0
  fi                                                      
  # Now, which file do we rename? (The other retains the original name)
  case $rename_which in
    new)
      mv "$temp_name" "$file${the_suffix}"
      ;;
    old)
      mv "$file" "$file${the_suffix}"
      mv "$temp_name" "$file"
      ;;
    *)
      mv "$temp_name" "$file"
      ;;
  esac
done                                                       
exit
# $Log$
# Revision 1.1  2002/10/29 00:56:57  pwagner
# First commit
#
