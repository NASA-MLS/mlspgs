#!/bin/sh
# Copyright 2008, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# $Id$

# --------------- expandl2cf.sh help
# This shell script expands an l2cf template supplied as an arg

# Used when expanding l2cf files
#
# Usage:
# expandl2cf.sh [opt1] ..  [optn] template
#
#    O p t i o n s
# -I path             pass -I path to m4
# -Dmacro             pass -Dmacro to m4; may be repeated
# -Df file            pass macro definitions found in file, too
# -E envstg           set environment according to envstg
#                        e.g., "-E TEST=GC-01" means "export TEST=GC-01
# -Ef file            set environment according to definitions found in file
# -dot file          "dot" file before executing m4;
#                        it may contain env settings and much more besides
# -dryrun             merely echo the command that would be executed
# -m4 cmd             use cmd instead of m4
# -o file             store expanded l2cf in file instead of stdout
# -v                  verbose; prints handy summary at end of l2cf
# -w                  wrap lines in l2cf
# -example            print brief example of how to use; exit
# -h[elp]             print brief help message; exit
#
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede the template on the command line
# (2) The default M4PATH of ~/mlspgs/l2/l2cf/lib will be used
#     unless M4PATH is already defined or the -I option supplied
# (3) The -w option assumes that wrapLines exists, is executable,
#     and is in your PATH
# (4) -Dmacro definitions, plus any in the -Df file, will be combined
#     and passed to m4
# (5) If the line
#      M4=/some/path/to/m4
#     appears in the macros file, it will have the same effect as if
#       -m4 /some/path/to/m4
#     was among the command-line options
# Result:
# An expanded l2cf is written using appropriate macros
# --------------- End expandl2cf.sh help
# --------------- expandl2cf.sh example
# Example:
# Assume macros.txt contains the following (without the '#')
## start of macros.txt
# M4=$HOME/bin/m4
# machine=me 
# day=2008d037
# flagCopyStandardProducts 
# flagWriteAPrioriToDGG 
# flagWriteL2CFToDGM 
# flagUncompressRadiances
# V2ID=$TEST-$SUBTEST
## end of macros.txt
# and that env.txt contains the following (without the '#')
## start of env.txt
# TEST=GC-01
# SUBTEST=061
## end of env.txt
#
# Then if you execute this script as
#  expandl2cf.sh -Ef env.txt -Df macros.txt -o out.l2cf l2cf.m4
# the following commands will be executed
#  export TEST=GC-01
#  export SUBTEST=061
#  $HOME/bin/m4 -Dmachine=me -Dday=2008d037 \
#  -DflagCopyStandardProducts -DflagWriteAPrioriToDGG \
#  -DflagWriteL2CFToDGM -DflagUncompressRadiances \
#  -DV2ID=$TEST-$SUBTEST \
#  l2cf.m4 > out.l2cf
# --------------- End expandl2cf.sh example

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
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
}
      
#---------------------------- read_file_into_array
#
# read each line of stdin
# catenating them into an array which we will return
# Ignore any entries beginning with '#' character
# In fact, only first entry in each line is kept
# Possible improvements:
#   Other comment signifiers
#   Choose field number other than 1
# Notes and bugs:
# We treat lines with ' and " characters specially:
# on such lines we replace every space with a '&' character
# You will probably want to undo that replacement when you
# actually use the line later

read_file_into_array()
{
  array_result=''
  while read line; do
    element=""
    acomment=`echo $line | grep '^#'`
    a=`echo $line | grep [\'\"]`
    if [ "$acomment" != "" ]
    then
      # Dont do anything--just a comment
      acomment="$acomment"
    elif [ "$a" != "" ]
    then
      element=`echo "$line" | sed 's/ /\&/g'`
    else
      element=`echo $line | awk '$1 !~ /^#/ {print $1}'`
    fi
    if [ "$element" != "" ]
    then
      array_result="$array_result $element"
    fi
  done
  echo $array_result
}
      
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

cmdline="$0 $@"
debug="no"
dotfile=""
dryrun="no"
envfile=""
I=expandl2cf
l2cf="STDOUT"
M4=m4
macros=""
macrofile=""
me="$0"
my_name=expandl2cf.sh
myPATH=""
stempl2cf=$HOME/`get_unique_name l2cf2`
templ2cf=$HOME/`get_unique_name l2cf1`
verbose="no"
wrap="no"

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    # echo "option: $1"
    case "$1" in

    -Df )
       shift
       macrofile="$1"
       shift
       ;;
    -D* )
       macros="$macros $1"
       shift
       ;;
    -d )
       debug="yes"
       verbose="yes"
       shift
       ;;
    -dryrun )
       dryrun="yes"
       shift
       ;;
    -Ef )
       shift
       envfile="$1"
       shift
       ;;
    -E )
       shift
       eval export `echo $1`
       shift
       ;;
    -E* )
       a=`echo $b | sed 's/^-E//'`
       eval export `echo $a`
       shift
       ;;
    -dot )
       shift
       dotfile="$1"
       shift
       ;;
    -example )
       sed -n '/'$my_name' example/,/End '$my_name' example/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       rm -f $settings_file
       exit
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       rm -f $settings_file
       exit
       ;;
    -I )
       shift
       myPATH="$1"
       shift
       ;;
    -m4 )
       shift
       M4="$1"
       shift
       ;;
    -v )
       verbose="yes"
       shift
       ;;
    -w )
       wrap="yes"
       shift
       ;;
    -o )
       shift
       l2cf="$1"
       shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

if [ "$debug" = "yes" ]
then
  echo "dotfile: $dotfile"
  echo "envfile: $envfile"
  echo "macros: $macros"
  echo "macrofile: $macrofile"
  echo "dryrun: $dryrun"
  echo "mypath: $mypath"
  echo "M4: $M4"
  which m4
  echo "wrap: $wrap"
  echo "template: $1"
  echo "l2cf: $l2cf"
  echo "templ2cf: $templ2cf"
fi

if [ "$1" = "" ]
then
  echo "Sorry--no template found among args"
  exit 1
fi

if [ -f "$dotfile" ]
then
  . "$dotfile"
fi

if [ -f "$envfile" ]
then
  envlines=`cat $envfile | uniq | read_file_into_array`
  for linenosp in $envlines
  do
    line=`echo $linenosp | sed 's/\&/ /g'`
    eval export `echo $line`
  done
fi

if [ -f "$macrofile" ]
then
  l2cflines=`cat $macrofile | uniq | read_file_into_array`
  for linenosp in $l2cflines
  do
    line=`echo $linenosp | sed 's/\&/ /g'`
    a=`echo $line | grep 'M4='`
    # echo "line: $line"
    # echo "a: $a"
    if [ "$a" != "" ]
    then
      eval $line
      # echo $M4
      # exit
    else
      macros="-D${line} $macros"
    fi
  done
fi

if [ "$mypath" = "" -a "$M4PATH" = "" ]
then
  if [ -d "$HOME/mlspgs/l2/l2cf/lib" ]
  then
    mypath="$HOME/mlspgs/l2/l2cf/lib"
  fi
fi

if [ "$mypath" != "" ]
then
  ALLOPTS="-I $mypath $macros"
else
  ALLOPTS="$macros"
fi

if [ "$dryrun" = "yes" ]
then
  echo "$M4 $ALLOPTS $1 > $templ2cf"
  exit 0
else
  eval $M4 $ALLOPTS $1 > $templ2cf
fi

if [ "$wrap" = "yes" ]
then
  mv $templ2cf $stempl2cf
  wrapLines -v -blank 1 < $stempl2cf > $templ2cf
  rm $stempl2cf
fi

if [ "$l2cf" = "STDOUT" ]
then
  cat $templ2cf
  rm $templ2cf
  exit 0
else
  mv $templ2cf $l2cf
fi

if [ "$verbose" != "yes" ]
then
  exit 0
fi

echo ";;; Expanded `date` by expandl2cf.sh" >> $l2cf
echo ";;; $cmdline" >> $l2cf

if [ -f "$dotfile" ]
then
  dotlines=`cat $dotfile | uniq | read_file_into_array`
  echo " " >> $l2cf
  echo ";;; ---------- contents of $dotfile -----------" >> $l2cf
  for linenosp in $dotlines
  do
    line=`echo $linenosp | sed 's/\&/ /g'`
    echo ";;; ${line}" >> $l2cf
  done
fi

if [ -f "$envfile" ]
then
  echo " " >> $l2cf
  echo ";;; ---------- contents of $envfile -----------" >> $l2cf
  for linenosp in $envlines
  do
    line=`echo $linenosp | sed 's/\&/ /g'`
    echo ";;; ${line}" >> $l2cf
  done
fi

if [ -f "$macrofile" ]
then
  echo " " >> $l2cf
  echo ";;; ---------- contents of $macrofile -----------" >> $l2cf
  for linenosp in $l2cflines
  do
    line=`echo $linenosp | sed 's/\&/ /g'`
    echo ";;; ${line}" >> $l2cf
  done
fi

exit 0
# $Log$
# Revision 1.1  2008/07/16 21:03:04  pwagner
# First commit
#
