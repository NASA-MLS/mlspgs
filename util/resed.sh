#!/bin/sh
#resed.sh
# --------------- resed.sh help
# Runs the command line arguments through sed, possibly
# (depending on option(s) selected)
# replacing them with the changes sed brings about
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
# -example      print brief example; exit
# -rn           rename new file "old_name"+suffix, keeping old as old_name
# -ro           rename new file "old_name", keeping old as old_name+suffix
# -rd           replace the file only if the editing commands change it
# -print        instead of creating new file(s), print results of sed to stdout
# -suffix=xxx   suffix to use when renaming either new or old file
#               (w/o any separator; e.g., to add ".bak" use -suffix=.bak)
# -sed mysed    run mysed instead of sed (e.g. /opt/enhanced/bin/mightysed)
# filen         file name (with path)
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede any filen on the command line
# (2) The options -c and -f are mutually exclusive, but one is mandatory
#     i.e., you must pass commands to sed either by file or by command-line
# (3) The options -rd, -rn and -ro are mutually exclusive
#     if you choose -rn but omit suffix, the new files will be "old_name.new"
#     if you choose -ro but omit suffix, the old files will be "old_name.old"
#     if you choose -rd and the new file is different, or
#     if you choose none the new file will simply replace the old one
#     (leaving you with no backup--so test with -print first)
# 
# Result:
# Files in the list may be modified or new files created
# according to the the options you supply
# --------------- End resed.sh help
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
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
#     though you will have to check whether the delimiters "/", ":", etc.
#     are present among the strings
# (2) Also allow inserting a block of text stored in a file after  
#     line nnn via -l nnn -b block_file
# (3) Instead of messing with suffixes, let the script edit files
#     from one dir, saving the modified versions in another
#     via -d1 d_orig -d2 d_mod
#     (you will have to disable reecho part if you do this)
#****************************************************************
# --------------- resed.sh example
# Example:
# I have a set of files, in each of which I wish to change one line
# in which a variable "init_m_dir" is currently being set to a value "l1":
# ../mlspgs:76% grep 'init_m_dir = l1' tests/fwdmdl/platforms/[AD-Z]* tests/cloudfwdm/platforms/[AD-Z]*
# tests/fwdmdl/platforms/LF95.Linux:init_m_dir = l1
# tests/fwdmdl/platforms/NAG.Linux:init_m_dir = l1
# tests/fwdmdl/platforms/NAG.SGI:init_m_dir = l1
# tests/fwdmdl/platforms/NAG.Sun:init_m_dir = l1
# tests/fwdmdl/platforms/Sun.Sun:init_m_dir = l1
# tests/fwdmdl/platforms/Unknown.None:init_m_dir = l1
# tests/cloudfwdm/platforms/LF95.Linux:init_m_dir = l1
# tests/cloudfwdm/platforms/NAG.Linux:init_m_dir = l1
# tests/cloudfwdm/platforms/NAG.SGI:init_m_dir = l1
# tests/cloudfwdm/platforms/NAG.Sun:init_m_dir = l1
# tests/cloudfwdm/platforms/Sun.Sun:init_m_dir = l1
# tests/cloudfwdm/platforms/Unknown.None:init_m_dir = l1
# 
# So to set that variable instead to the value "l2" I would first test
# (to safely assure that my changes accomplish what I want):
# ../mlspgs:77% util/resed.sh -print -c 's/init_m_dir = l1/init_m_dir = l2/' tests/fwdmdl/platforms/[AD-Z]* tests/cloudfwdm/platforms/[AD-Z]* | grep 'init_m_dir = l2'
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
# init_m_dir = l2
#
# and then reenter the last command w/o the -print option
# --------------- End resed.sh example
#****************************************************************
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
DEEBUG=off
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
# Possible values for rename_which: {new, old, neither, diff}
rename_which="neither"
the_suffix=""
print_to_stdout="no"
SED="sed"
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
    -sed )
       SED="$2"
       shift
       shift
       ;;
    -rn )
       rename_which="new"
       shift
       ;;
    -ro )
       rename_which="old"
       shift
       ;;
    -rd )
       rename_which="diff"
       shift
       ;;
    -print )
       print_to_stdout="yes"
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
    -example )
       sed -n '/'$my_name' example/,/End '$my_name' example/ p' $me \
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
   echo "print_to_stdout? $print_to_stdout"
   echo "sed command to use $SED"
   echo "remaining args $@"
fi
# Check for forbidden arguments
if [ "$the_command" != "" -a "$command_file" != "" ]
then
  echo 'You cannot have both a command file and command-line'
  exit
elif [ "$command_file" != "" ]
then
  the_opt="$the_opt -f $command_file"
elif [ "$the_command" = "" ]
then
#
  echo 'You must have either a command file or a command-line'
  exit
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
  exit
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
    # sed $the_opt "$the_command" $file > "$temp_name"
    $SED $the_opt "$the_command" $file > "$temp_name"
  else
    # echo sed $the_opt $file
    # sed $the_opt $file > "$temp_name"
    $SED $the_opt $file > "$temp_name"
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
  elif [ "$print_to_stdout" = "yes" ]
  then
    echo "== $file =="  
    cat "$temp_name"    
    rm  "$temp_name"    
    # exit 0
  elif [ "$NO_REPLACE" = "on" ]
  then
    echo "Created $temp_name from $file"
  else
    # Now, which file do we rename? (The other retains the original name)
    case $rename_which in
      new)
        mv "$temp_name" "$file${the_suffix}"
        ;;
      old)
        mv "$file" "$file${the_suffix}"
        mv "$temp_name" "$file"
        ;;
      diff)
        if [ "`diff $file $temp_name`" != "" ]
        then
          mv "$temp_name" "$file"
        fi
        ;;
      *)
        mv "$temp_name" "$file"
        ;;
    esac
  fi                                                      
done                                                       
exit
# $Log$
# Revision 1.2  2002/10/30 01:11:21  pwagner
# Actually works now
#
# Revision 1.1  2002/10/29 00:56:57  pwagner
# First commit
#
