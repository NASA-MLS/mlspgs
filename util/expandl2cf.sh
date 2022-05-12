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
# -Cf file            use combined definitions found in file for environment 
#                       and macros
# -I path             pass -I path to m4
# -Dmacro             pass -Dmacro to m4; may be repeated
# -Df file            pass macro definitions found in file, too
# -E envstg           set environment according to envstg
#                        e.g., "-E TEST=GC-01" means "export TEST=GC-01
# -Ef file            set environment according to definitions found in file
# -nmacros            don't "dot" macros file
# -dot file           "dot" file before executing m4;
#                        it may contain env settings and much more besides
# -dryrun             merely echo the command that would be executed
# -i                  ignore any TEMPLATE defs in env or macros files; use
#                        final arg on commandline instead
# -ident              print idents of l2cf fragments at end of file
# -m4 cmd             use cmd instead of m4
# -o file             store expanded l2cf in file instead of stdout
# -v                  verbose; prints handy summary at end of l2cf
# -d                  debug; even more verbose
# -w                  wrap lines in l2cf
# -example            print brief example of how to use; exit
# -example--dotfile   print brief example of how to use dotfile; exit
# -h[elp]             print brief help message; exit
#
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede the template on the command line
# (2) The default M4PATH of ~/mlspgs/l2/l2cf/lib will be used
#     unless M4PATH is already defined or the -I option supplied
# (3) -Dmacro definitions, plus any in the -Df file, will be combined
#     and passed to m4
# (4) If the line
#      M4=/some/path/to/m4
#     appears in the macros file, it will have the same effect as if
#       -m4 /some/path/to/m4
#     was among the command-line options
# (5) If the line
#      TEMPLATE=/some/path/to/template
#     appears in either the macros file or the env settings file,
#     and "-i" is not among the command line options,
#     it will have the same effect as if
#       /some/path/to/template.m4
#     was the template in command-line arguments
# (6) Use a dot file in place of an env file
#     (or in addition)
#     if you need to use shell control structures or other features
#     For an example of this see -example--dotfile
# Result:
#     An expanded l2cf is written using appropriate macros
# Other settings in the .macros or .env file and their effect:
#   WIDTH          maximum line length if wrapping lines
#   IDENTMAKER     use this instead of $HOME/mlspgs/util/identl2cf.sh
#   SETREAD        use this instead of `which set_read_env.sh`
#   WRAPLINES      use this instead of `which wrapLines`
# Other uses:
# (PCF)
#     expandl2cf.sh can also be used to expand a PCF template
# (no overridepaths.sh)
#     append to the macros file all the paths normally found in the
#     overridepaths.sh file
#
# Bugs and limitations:
# (1) set_read_env.sh must be both executable and in your path
#     (unless you set SETREAD in your .env file)
# (2) The -w option assumes that wrapLines exists, is executable,
#     and is in your PATH
# (3) Nobody uses the "dotfile"--should we remove it as an option?
# (4) Nobody uses overridepaths.sh--should we omit mentioning it?
# (5) How about substituting a variable for the ';' comment character?
# --------------- End expandl2cf.sh help
# --------------- expandl2cf.sh example
# Example:
# Assume macros.txt contains the following
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

# --------------- expandl2cf.sh dotfile example
# Example:
# Assume SO.macros contains the following (without the '#')
## start of SO.macros
# file of m4 macros
# each line contains a macro definition that will be passed to m4
# e.g., a line consisting of (without the #)
#category=noGoodBand17
# would add the command line option "-Dcategory=noGoodBand17"
# to the invocation of m4
# M4=$HOME/bin/m4
# EXPANDL2CF=${HOME}/mlspgs/util/expandl2cf.sh
# machine=scramjet
# day=$day
# flagComputeAvks
# l2pcVersion=$l2pcVersion
# V2ID=$TEST-$SUBTEST
# L1BSIDSVERSION=s7-f52
# l2pcversion=$l2pcVersion
# pfaversion=FS-03
# outpathl2mtx=$JOBDIR/outputs
# inpathleapsec=/science/pge/v0223/toolkit5.2.14/database/common/TD
## end of SO.macros
#
# and that env.txt contains the following
## start of SO.sh
#!/bin/sh
# dot file to establish environmental settings
# executed prior to the invocation of m4
# These will be expanded in your overridepaths.sh file
# as well as in your file of m4 macros (if any)
# These are not m4 macros, however, and so they
# will not be automatically expanded in your l2cf template
#
# Use dot file instead of env file if you wish to use
# shell control structures or commands
#
# Don't forget to export at the end
#
# Change these with each test
# day=2008d037
# CHUNK=96
# TEST=SO-01
# SUBTEST=a0096
# L2CFVERSION=v3-09-so3001
# TEMPLATE=v3-09-michaelavgkrnl.m4
# l2pcVersion=v3-00-FS-05
# case "$l2pcVersion" in
#   "v3-00-HO-01")
#     INPATHL2PC=/data1/pwagner/l2pc_30H1
#     ;;
#   "v3-00-FS-05")
#     INPATHL2PC=/data1/pwagner/l2pc_305
#     ;;
# esac
# PGE_BINARY_DIR=/home/pwagner/shrnkwrp-v3/pge/bin/IFC.Linux
# PGE_SCRIPT_DIR=${HOME}/shrnkwrp-v3/pge/util
# PGE_ROOT=/data1/pwagner_home/toolkits/toolkit5.2.14/bin/linux
# PGSHOME=/home/pwagner/shrnkwrp-v3/toolkit5.2.14
# JOBDIR="/data1/pwagner/l2tests/avgkrnls/$TEST/$SUBTEST"
# # Retrieval version
# OTHEROPTS="--crash --ntk -g --chunk $CHUNK --skipDirec -S'l2q,glob,mas,chu,opt1,log,pro1,time,apr'"
# #
# # Don't forget to export at the end
# export day CHUNK TEST SUBTEST L2CFVERSION TEMPLATE l2pcVersion INPATHL2PC
# export PGE_BINARY_DIR PGE_SCRIPT_DIR PGE_ROOT PGSHOME JOBDIR OTHEROPTS
## end of SO.sh
#
# Then you would execute this script as
#  expandl2cf.sh -dot SO.sh -Df SO.macros -o out.l2cf
# --------------- End expandl2cf.sh dotfile example

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
      
#---------------------------- write_file_to_stdout
#
# read each line of stdin
# writing it to stdout,
# optionally prefixing each line with $1, e.g. comment characters

write_file_to_stdout()
{
  while read line; do
     echo "$1 $line"
  done
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
combfile=""
debug="no"
dotfile=""
dryrun="no"
envfile=""
expandmacros="yes"
I=expandl2cf
ident="no"
#IDENTMAKER=identl2cf.sh
IDENTMAKER=""
ignore="no"
l2cf="STDOUT"
M4=m4
macros=""
macrofile=""
me="$0"
my_name=expandl2cf.sh
mypath=""
stempl2cf=$HOME/`get_unique_name l2cf2`
templ2cf=$HOME/`get_unique_name l2cf1`
TEMPLATE=""
verbose="no"
wrap="no"

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    # echo "option: $1"
    case "$1" in

    -Cf )
       shift
       combfile="$1"
       shift
       ;;
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
    -i )
       ignore="yes"
       shift
       ;;
    -ident* )
       ident="yes"
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
    -example--dotfile )
       sed -n '/'$my_name' dotfile example/,/End '$my_name' dotfile example/ p' $me \
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
       mypath="$1"
       shift
       ;;
    -m4 )
       shift
       M4="$1"
       shift
       ;;
    -nmacros )
       expandmacros="no"
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

if [ -f "./$dotfile" ]
then
  . "./$dotfile"
elif [ -f "$dotfile" ]
then
  . "$dotfile"
fi

# Did we use a single file combining env and macros?
# If so, disgorge the separate env and macros files
if [ -f "$combfile" ]
then
  envfile="$combfile".env
  macrofile="$combfile".macros
  sed -n '/- e n v/,/- m 4   m/ p' $combfile | sed '$ d' > $envfile
  sed -n '/- m 4   m/,$ p' $combfile | sed '1 d' > $macrofile
fi

# Are there any special settings we should be aware of in the env file?
# Namely, do we set alternatives for SETREAD, etc.?
# If so, we may need to set them now, expecially SETREAD
# because it's used to set the others.
if [ -f "$envfile" ]
then
  SETREAD=`grep SETREAD $envfile | sed 's/.*=//'`
  WRAPLINES=`grep WRAPLINES $envfile | sed 's/.*=//'`
fi

if [ "$SETREAD" = "" ]
then
  SETREAD=`which set_read_env.sh`
fi

if [ "$WRAPLINES" = "" ]
then
  WRAPLINES=`which wrapLines`
fi

if [ -f "$envfile" ]
then
    . $SETREAD < $envfile
else
    echo ";;; Warning--no envfile"
fi

if [ -f "$macrofile" ]
then
  if [ "$expandmacros" = "yes" ]
  then
    . $SETREAD < $macrofile
  fi
  l2cflines=`cat $macrofile | uniq | read_file_into_array`
  for linenosp in $l2cflines
  do
    line=`eval echo $linenosp | sed 's/\&/ /g'`
    a=`echo $line | grep 'M4='`
    assignment=`echo $line | grep '='`
    atemplate=`echo $line | grep 'TEMPLATE='`
    # echo "line: $line"
    # eval echo "line: $line"
    # echo "a: $a"
    if [ "$a" != "" ]
    then
      eval $line
      # echo $M4
      # exit
    elif [ "$atemplate" != "" ]
    then
      eval $line
      # echo $TEMPLATE
      # exit
    else
      macros="$macros -D${line}"
    fi
  done
else
    echo ";;; Warning--no macrofile"
fi

if [ "$IDENTMAKER" = "" ]
then
  IDENTMAKER=$HOME/mlspgs/util/identl2cf.sh
fi

if [ "$debug" = "yes" ]
then
  echo "dotfile: $dotfile"
  echo "dryrun: $dryrun"
  echo "envfile: $envfile"
  echo "macros: $macros"
  echo "macrofile: $macrofile"
  echo "dryrun: $dryrun"
  echo "expandmacros: $expandmacros"
  echo "mypath: $mypath"
  echo "M4: $M4"
  which m4
  echo "wrap: $wrap"
  echo "template: $1"
  echo "l2cf: $l2cf"
  echo "templ2cf: $templ2cf"
  echo "PATH $PATH"
  echo "M4PATH $M4PATH"
fi

if [ "$1" = "" -a "$TEMPLATE" = "" ]
then
  echo "Sorry--no template found among args"
  exit 1
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

# Did we supply TEMPLATE?
if [ "$TEMPLATE" = "" -o "$ignore" = "yes" ]
then
  TEMPLATE="$1"
fi

if [ "$dryrun" = "yes" ]
then
  echo "$M4 $ALLOPTS $TEMPLATE > $templ2cf"
  exit 0
elif [ "$debug" = "yes" ]
then
  echo "$M4 $ALLOPTS $TEMPLATE  $templ2cf"
  eval $M4 $ALLOPTS $TEMPLATE > $templ2cf
else
  eval $M4 $ALLOPTS $TEMPLATE > $templ2cf
fi

#exit 0

if [ "$wrap" = "yes" ]
then
  mv $templ2cf $stempl2cf
  # WRAPOPTS="-v -blank 1"
  WRAPOPTS="-blank 1"
  if [ "$WIDTH" != "" ]
  then
    WRAPOPTS="$WRAPOPTS -width $WIDTH"
  fi
  if [ "$verbose" = "yes" ]
  then
    WRAPOPTS="-v $WRAPOPTS"
    echo " ; --- expandl2cf.sh settings ---"      > $templ2cf
    echo " ; dotfile: $dotfile"                  >> $templ2cf
    echo " ; dryrun: $dryrun"                    >> $templ2cf
    echo " ; envfile: $envfile"                  >> $templ2cf
    # echo " ; macros: $macros"                    >> $templ2cf
    echo " ; macrofile: $macrofile"              >> $templ2cf
    echo " ; dryrun: $dryrun"                    >> $templ2cf
    echo " ; expandmacros: $expandmacros"        >> $templ2cf
    echo " ; mypath: $mypath"                    >> $templ2cf
    echo " ; M4: $M4"                            >> $templ2cf
    echo " ; --- End expandl2cf.sh settings ---" >> $templ2cf
    $WRAPLINES $WRAPOPTS -i $stempl2cf            >> $templ2cf
  else
    $WRAPLINES $WRAPOPTS -i $stempl2cf            > $templ2cf
  fi
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

if [ "$ident" = "yes" ]
then
echo ";;; ----- idents -----" >> $l2cf
echo ";;; Note: any of the following files may differ from" >> $l2cf
echo ";;; the cvs repository contents if they have been locally modified" >> $l2cf
  if [ "$mypath" != "" ]
  then
    echo ";;; (from $mypath)" >> $l2cf
    $IDENTMAKER -I $mypath $TEMPLATE >> $l2cf
  else
    echo ";;; (from ${HOME}/mlspgs/l2/l2cf/lib)" >> $l2cf
    $IDENTMAKER $TEMPLATE >> $l2cf
  fi
fi

if [ "$verbose" != "yes" ]
then
  exit 0
fi

echo ";;; Expanded `date` by expandl2cf.sh" >> $l2cf
echo ";;; $cmdline" >> $l2cf

if [ -f "$dotfile" ]
then
  echo " " >> $l2cf
  echo ";;; ---------- contents of $dotfile -----------" >> $l2cf
  cat $dotfile | write_file_to_stdout ";;;" >> $l2cf
fi

if [ -f "$envfile" ]
then
  echo " " >> $l2cf
  echo ";;; ---------- contents of $envfile -----------" >> $l2cf
  cat $envfile | write_file_to_stdout ";;;" >> $l2cf
fi

if [ -f "$macrofile" ]
then
  echo " " >> $l2cf
  echo ";;; ---------- contents of $macrofile -----------" >> $l2cf
  cat $macrofile | write_file_to_stdout ";;;" >> $l2cf
fi

exit 0
# $Log$
# Revision 1.14  2017/12/07 00:36:24  pwagner
# Add a new comdline option -Cf combining -Df and -Ef in a single file
#
# Revision 1.13  2017/03/28 20:34:50  pwagner
# Handle case where . is not in PATH
#
# Revision 1.12  2016/12/16 22:00:01  pwagner
# verbose setting affects us and wrapLines, too
#
# Revision 1.11  2016/06/01 16:34:20  pwagner
# Expand l2cf even if debug; fallback for IDENTMAKER now in util
#
# Revision 1.10  2013/12/04 21:36:58  pwagner
# env variable WIDTH can set width if wrapping lines
#
# Revision 1.9  2013/05/30 20:38:52  pwagner
# Fixed error where set_read_env.sh not local
#
# Revision 1.8  2013/01/09 18:50:07  pwagner
# Cant recall why, but made this change
#
# Revision 1.7  2009/06/16 22:37:08  pwagner
# Added -ident option to store version id for l2cf fragments
#
# Revision 1.6  2009/04/13 20:43:17  pwagner
# Fixed a bug preventing macros file from using its own macros properly
#
# Revision 1.5  2009/03/26 20:24:08  honghanh
# Change myPATH to mypath so the -I option can work
#
# Revision 1.4  2009/03/13 17:26:30  pwagner
# Many improvements; tested with no overridepaths.sh file
#
# Revision 1.3  2008/08/18 17:35:00  pwagner
# Added -i, -example--dotfile options
#
# Revision 1.2  2008/07/16 23:37:01  pwagner
# Improved example; consistent with email announcement
#
# Revision 1.1  2008/07/16 21:03:04  pwagner
# First commit
#
