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
# (3) create new files, based on contents of existing ones
#
# Usage:
# resed.sh [opt] ..  [file1] [file2 ..]
#
#    O p t i o n s
# -o options    pass options to sed (default is none)
# -c command    command to pass (surrounded by ') to sed
#                as an alternative to "-c command", you may use the next pair
# -os oldstring string to be replaced
# -rs oldstring string to replace it with
# -Of os_file   file with list of old strings to be replaced
# -Rf rs_file   file with list of replacement strings
# -g            "globally" replace old string with new
# -f file       file of commands to pass through to sed
# -d1 dir1      operate on every file in dir1
# -d2 dir2      storing results in dir2 (if dir2 doesn't exist, it creates it)
#                (if omitted, replace each file in dir1 with its result)
# -[n]grep text restrict sed to only those files [not] containing text
# -h[elp]       print brief help message; exit
# -example      print brief example; exit
# -name=xxx     name to use when renaming either new or old file
# -rn           rename new file, letting old file retain old_name
# -ro           rename old file, letting new file inherit old_name
# -rd           replace old file with new only if the editing commands change it
# -dryrun       instead of creating new file(s), print results to stdout
# -suffix=xxx   suffix to use when renaming either new or old file
#               (w/o any separator; e.g., to add ".bak" use -suffix=.bak)
# -sed mysed    run mysed instead of sed (e.g. /opt/enhanced/bin/mightysed)
# filen         file name (with path)
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede any file on the command line
# (2) If you supply either -suffix=xxx or -name=xxx then
#     you must also choose one of -rn, -ro, or -rd
# (3) The options -c and -f are mutually exclusive, but one is mandatory
#     i.e., you must pass commands to sed either by file or by command-line
# (4) The options -rd, -rn and -ro are mutually exclusive
#     if you omit suffix=xxx and name=xxx, then with 
#     (option)    old file becomes       new file becomes
#      -rn            old_name           "old_name.new"
#      -ro          "old_name.old"          old_name
#      -rd            (deleted)             old_name
#      (none)         (deleted)             old_name
#     To repeat, if you choose -rd and the new file is different, or
#     if you choose none the new file will simply replace the old one
#     (leaving you with no backup--so test with -dryrun first)
# (5) The options -name=xxx and -sufffix=xxx are mutually exclusive
# (6) With the option -d1 dir1 don't name any files on the command line
#     -- they will be ignored
# (7) Unless you supply one of options -d2, -rn, or -ro resed will overwrite
#     each filen it modifies. This could be dangerous. You have been warned.
# (8) Using -Of and -Rf lets you you replace not just a single string
#     with another in each filen, but as many as you have in os_file and rs_file
# 
# Result:
# Files in the list may be modified or new files created
# according to the the options you supply
# --------------- End resed.sh help
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

#---------------------------- read_file_into_array
#
# read each line of stdin
# catenating them into an array which we will return
# Ignore any entries beginning with '#' character
# In fact, only first entry in each line is kept
# Possible improvements:
#   Other comment signifiers
#   Choose field number other than 1

read_file_into_array()
{
  array_result=''
  while read line; do
    element=`echo $line | awk '$1 !~ /^#/ {print $1}'`
    if [ "$element" != "" ]
    then
      array_result="$array_result $element"
    fi
  done
  echo $array_result
}

#---------------------------- print_nth_array_element
#
# Function to print the nth element in a space-delimited list
# where n is the first arg, the list is the second arg, 
# print_nth_array_element 3 'a b c d'
# writes 'c' to standard output

print_nth_array_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: print_nth_array_element n 'a b c ..'"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); print $parts[$ARGV[1]-1]' "$2" "$1"
}

#
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script    *
#****************************************************************
#
#   Notes
# Bugs: 
# (1) If the old_name+suffix exceeds the longest possible file name
#     the script does not fail gracefully
# (2) If filen+suffix is the same as filem for some pair (n,m)
#     one will replace the other w/o any warning of this possibility
# (3) The renaming mechanism suffix=xxx and name=xxx has not been
#     implemented for the option -rd yet; maybe it makes no sense to do so
# (4) When using -os "string1" -rs "string2"
#     you must not have both delimiters "/", ":" present among the strings
# Unimplemented improvements: 
# (1) Also allow inserting a block of text to be stored in a file after  
#     line nnn via -l nnn -b block_file
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
# ../mlspgs:77% util/resed.sh -dryrun -c 's/init_m_dir = l1/init_m_dir = l2/' tests/fwdmdl/platforms/[AD-Z]* tests/cloudfwdm/platforms/[AD-Z]* | grep 'init_m_dir = l2'
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
# and then reenter the last command w/o the -dryrun option
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
# $mkpath is mkpath with me's path prepended
mkpath="`echo $0 | sed 's/'$I'/mkpath/'`"
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
DEEBUG=off
NO_REPLACE=off
if [ $DEEBUG = "on" ]
then
   echo "Called me as $0"
   echo "with args $@"
fi

the_text=""
reecho_opt=""
command_file=""
the_command=""
the_opt=""
all_the_opts="$@"
more_opts="yes"
# Possible values for rename_which: {new, old, neither, diff}
rename_which="neither"
the_suffix=""
dir1=""
dir2=""
new_name=""
old_string=""
new_string=""
os_file=""
rs_file=""
dryrun="no"
globally="no"
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
    -os )
       old_string="$2"
       shift
       shift
       ;;
    -rs )
       new_string="$2"
       shift
       shift
       ;;
    -Of )
       os_file="$2"
       shift
       shift
       ;;
    -Rf )
       rs_file="$2"
       shift
       shift
       ;;
    -d1 )
       dir1="$2"
       shift
       shift
       ;;
    -d2 )
       dir2="$2"
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
    -dryrun )
       dryrun="yes"
       shift
       ;;
    -g )
       globally="yes"
       shift
       ;;
    -grep )
       shift
       the_text="$1"
       reecho_opt="-grep"
       shift
       ;;
    -ngrep )
       shift
       the_text="$1"
       reecho_opt="-ngrep"
       shift
       ;;
    -name=* )
       new_name=`echo $1 | sed 's/-name=//'`
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
   echo "old_string $old_string"
   echo "new_string $new_string"
   echo "os_file $os_file"
   echo "rs_file $rs_file"
   echo "the_suffix $the_suffix"
   echo "new_name $new_name"
   echo "rename_which? $rename_which"
   echo "dryrun? $dryrun"
   echo "sed command to use $SED"
   echo "remaining args $@"
fi
# Check for forbidden arguments or inconsistencies
if [ "$the_command" != "" -a "$command_file" != "" ]
then
  echo 'You cannot have both a command file and command-line'
  exit
elif [ "$command_file" != "" ]
then
  the_opt="$the_opt -f $command_file"
elif [ "$old_string" != "" -a "$os_file" != "" ]
then
  echo 'You cannot have both an old string and a file of them'
  exit
elif [ "$new_string" != "" -a "$rs_file" != "" ]
then
  echo 'You cannot have both a replacement string and a file of them'
  exit
elif [ "$rs_file" != "" ]
then
  # We'll check later
  echo "$rs_file" > /dev/null
elif [ "$old_string" != "" ]
then
  # Which delimiter shall we use? / or : ?
  old_slashed=`echo "$old_string" | grep '/'`
  new_slashed=`echo "$new_string" | grep '/'`
  if [ "$old_slashed" = "" -a "$new_slashed" = "" ]
  then
    the_command='s/'$old_string'/'$new_string'/'
    # the_opt='s/'$old_string'/'$new_string'/'
  else
    the_command='s:'$old_string':'$new_string':'
    # the_opt='s:'$old_string':'$new_string':'
  fi
  if [ "$globally" = "yes" ]
  then
    the_command="$the_command"g
  fi
elif [ "$the_command" = "" ]
then
#
  echo 'You must have either a command file or a command-line'
  exit
fi

if [ "$new_name" != "" -a "$the_suffix" != "" ]
then
  echo 'You cannot specify both a new_name and a suffix'
  exit
elif [ "$rename_which" = "neither" ]
then
  if [ "$new_name" != "" -o "$the_suffix" != "" ]
  then
    echo 'You must specify whether to rename old or new'
    exit
  fi
fi

# Were we given -Of and -Rf?
# ------------------------------------------------------------
# Method: Create a new commandline for each old string and its replacement
# Then call ourself ($0) once for each new commandline
if [ "$os_file" != "" -a "$rs_file" != "" ]
then
  # Create arrays of old and replacement strings
  os_array=`read_file_into_array < "$os_file"`
  rs_array=`read_file_into_array < "$rs_file"`
  # Snip -Of and -Rf from the commandline to pass to $0
  temp_opts=$all_the_opts
  all_the_opts=`echo "$temp_opts" | sed "s/-Of $os_file//"`
  temp_opts=$all_the_opts
  all_the_opts=`echo "$temp_opts" | sed "s/-Rf $rs_file//"`
  # Iteratively call $0 with each old string and its replacement
  n=1
  for old_string in $os_array
  do
    new_string=`print_nth_array_element $n "$rs_array"`
    if [ $DEEBUG = "on" ]
    then
      echo "old_string $old_string"
      echo "new_string $new_string"
      echo "$0 -os $old_string -rs $new_string $all_the_opts"
    fi
    $0 -os $old_string -rs $new_string $all_the_opts
    n=`expr $n + 1`
  done
  # That's it
  exit 0
elif [ "$os_file" != "" ]
then
  echo "you must specify both -Of and -Rf files or neither"
  exit 1
elif [ "$rs_file" != "" ]
then
  echo "you must specify both -Of and -Rf files or neither"
  exit 1
fi
# ------------------------------------------------------------

if [ "$rename_which" = "new" -a "$the_suffix" = "" ]
then
  the_suffix=.new
elif [ "$rename_which" = "old" -a "$the_suffix" = "" ]
then
  the_suffix=.old
fi

if [ "$dir1" = "" ]
then
  the_list=`$reecho $reecho_opt $the_text "$@"`
else
  the_list=`$reecho -dir "$dir1" $reecho_opt $the_text '*'`
fi

if [ $DEEBUG = "on" ]
then
   echo "the_command $the_command"
   echo "the_suffix $the_suffix"
   echo "the_list $the_list"
fi

if [ "$the_list" = "" ]
then
  echo 'No valid files were found among the command line arguments'
  echo 'Make sure you entered their path/names correctly  (resed.sh)'
  exit
fi

for anyfile in $the_list                                      
do
  file="$anyfile"
  file_name=`$the_splitter -f $file`
  if [ "$dir1" != "" ]
  then
    # The path is the arg
    path="$dir1"
    file="$dir1/$anyfile"
  elif [ "$file" = "$file_name" ]
  then
    # No path supplied
    path=""
  else
    path=`$the_splitter -p $file`
  fi
  temp_name=`$unique_name resed`
  if [ "$path" != "" ]
  then
    if [ -w "$path" ]
    then
      temp_name="$path"/"$temp_name"
    fi
  fi
  if [ "$the_command" != "" ]
  then
    # echo sed $the_opt "'"$the_command"'" $file
    # sed $the_opt "$the_command" $file > "$temp_name"
    $SED $the_opt "$the_command" $file > "$temp_name"
  else
    # echo sed $the_opt $file
    # sed $the_opt $file > "$temp_name"
    # the_opt='s/'$old_string'/'$new_string'/'
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
  elif [ "$dryrun" = "yes" ]
  then
    echo "== $file =="  
    cat "$temp_name"    
    rm  "$temp_name"    
    # exit 0
  elif [ "$NO_REPLACE" = "on" ]
  then
    echo "Created $temp_name from $file"
  elif [ "$dir1" != "" -a "$dir2" != "" ]
  then
    if [ ! -d "$dir2" ]
    then
      $mkpath "$dir2"
    fi
    mv "$temp_name" "$dir2/$anyfile"
  elif [ "$dir1" != "" ]
  then
    mv "$temp_name" "$file"
  else
    # 1st--what do we call the renamed file
    if [ "$new_name" != "" ]
    then
      new_file_name=`$the_splitter -f $new_name`
      if [ "$new_name" != "$new_file_name" ]
      then
        new_path=`$the_splitter -p $new_name`
      elif [ -w "$path" ]
      then
        # No path supplied for new name--assume same as old path
        new_path="$path"
      else
        new_path=""
      fi
      call_it="$new_name"
      if [ "$new_path" != "" ]
      then
        call_it="$new_path"/"$new_name"
      fi
    else
      call_it="$file${the_suffix}"
    fi
    # Now, which file do we rename? (The other retains the original name)
    case $rename_which in
      new)
        mv "$temp_name" "$call_it"
        ;;
      old)
        mv "$file" "call_it"
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
# Revision 1.11  2009/12/10 18:50:16  pwagner
# Added debugging
#
# Revision 1.10  2009/07/28 21:57:23  pwagner
# print origin of error mesg
#
# Revision 1.9  2008/10/13 23:30:03  pwagner
# Added the -g option to globally do the subsitution old to new
#
# Revision 1.8  2006/06/07 17:06:49  pwagner
# Did this fix a longstanding bug, or create one?
#
# Revision 1.7  2005/06/23 22:20:46  pwagner
# Reworded Copyright statement
#
# Revision 1.6  2005/04/27 23:32:48  pwagner
# Many new options
#
# Revision 1.5  2003/05/22 23:35:41  pwagner
# Tried to fix problem with writing to read-only dirs
#
# Revision 1.4  2003/04/01 23:14:25  pwagner
# Added -name=xxx option as alternative to suffixing
#
# Revision 1.3  2003/02/28 19:14:09  pwagner
# Many new options, helpful example added
#
# Revision 1.2  2002/10/30 01:11:21  pwagner
# Actually works now
#
# Revision 1.1  2002/10/29 00:56:57  pwagner
# First commit
#
