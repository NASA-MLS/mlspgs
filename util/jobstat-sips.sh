#!/bin/sh
#	jobstat-sips.sh				Author P.A. Wagner (JPL)
# --------------- jobstat help
# Prepares a progress table of the chunks for a sips-like run of mlsl2
# during or after run
# Usage:
# jobstat [options] log_dir l2cf masterlog
# where log_dir is the directory holding all the log files from the run
# l2cf is the l2cf controlling the run
# and masterlog is the log file from the master task
#
#     O p t i o n s
#  -S scan-sript     Use scan-sript instead of mlsqlog-scan-sips.py
#  -t split-sript    Use split-sript instead of split_path.sh
#  -dryrun           Don't execute commands, just echo them
#  -v                verbose; print extra
#  -h[elp]           Show this help message
#  -vers number      The mlsl2 was version number (default is v2.2)

# Notes:
# (1) Requires the python script mlsqlog-scan-sips.py be in your PATH
#   (unless you override it by a command-line option)
# (2) if the log files have already been catenated into a single file
#    we will split them (temporarily)
#    (a) Assumes the files all started with the same or a similar line; e.g.
#     "#/users/pwagner/l2tests/v2.1/2006d121/pvmlog/riverrun.11864.log mlsl2.log"
#    (b) Requires the script split_path.sh be in your PATH
#      (unless you override it by a command-line option)
# (3) for mlsl2 versions prior to 2.2, a different format of the masterlog
#     is assumed to extract the name of the l2cf; version 2.2 or later share
#     a common format as far as extracting that name is concerned
#
# --------------- End jobstat help
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

#---------------------------- uncat_logs
#
# Function reverses the catenation performed on log files
uncat_logs()
{
#tempdir=$HOME/uncattemp
tempdir=$HOME/`get_unique_name uncat`
#the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
# $mkpath is mkpath with me's path prepended
# mkpath="`echo $0 | sed 's/'$I'/mkpath/'`"
# $reecho is reecho with me's path prepended
# reecho="`echo $0 | sed 's/'$I'/reecho/'`"
DEEBUG=off

if [ "$tempdir" = "$HOME/" -o "$tempdir" = "/" ]
then
  echo "Problem forming tempdir"
  exit 1
fi

mkdir $tempdir
#echo "args: $@"
#echo sed -n '2,$ p' $1 | grep 'mlsl2.log' | awk '{print $1}'
p=`sed -n '2,$ p' $1 | grep 'mlsl2.log' | awk '{print $1}'`
lastq=`sed -n '2,$ p' $1 | grep 'mlsl2.log' | tail -1 | awk '{print $1}'`
prevq=""
for q in $p
do
  # echo $q
  qname=`$the_splitter -f $q`
  # echo $qname
  if [ "$prevq" != "" ]
  then
    sed -n '/'$prevq'/,/'$qname'/ p' $1 | sed '$ d' > $tempdir/$prevq
  fi
  prevq=$qname
done
sed -n '/'$prevq'/,$ p' $1 > $tempdir/$prevq
echo $tempdir

}      

#---------------------------- execute
#
# Function either executes or else echoes a command
# depending on value of dryrun
execute()
{
  if [ "$dryrun" = "yes" ]
  then
    echo $@
  elif [ "$verbose" = "yes" ]
  then
    echo $@
    $@
  else
    $@
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
me="$0"
my_name=jobstat
I=uncat
TIDYUPAFTER=1
#           ^----- Set this to 1 to rm any temp files we create
dryrun="no"
verbose="no"
SCANNER="mlsqlog-scan-sips.py"
the_splitter=split_path.sh
version=2.2

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -S )
       shift
       SCANNER="$1"
       shift
       ;;
    -t )
       shift
       the_splitter="$1"
       shift
       ;;
    -v )
       verbose="yes"
       shift
       ;;
    -vers )
       shift
       version="$1"
       shift
       ;;
    -dryrun )
       dryrun="yes"
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

# Do we need to split the log files?
nlogs=`echo $1/*.log | wc -w`
if [ "$verbose" = "yes" ]
then
  echo "number of log files $nlogs"
fi

if [ "$nlogs" = 1 -a "$1/*.log" = "$1/mlsl2.log" ]
then
  echo "Splitting already-catenated logfiles $1/*.log"
  ls $1
  time usedir=`uncat_logs $1/*.log`
  execute $SCANNER $usedir xxx $2 $3 $version
  if [ "$TIDYUPAFTER" = "1" ]
  then
    echo /bin/rm -fr $usedir
    /bin/rm -fr $usedir
  fi
else
  execute $SCANNER $1 xxx $2 $3 $version
fi
exit 0
# $Log$
# Revision 1.3  2013/12/05 00:24:03  pwagner
# Added verbose option; explain why mlsl2 version matters
#
# Revision 1.2  2012/07/03 15:20:41  pwagner
# Fixed bug confusing single chunk with result of cat
#
# Revision 1.1  2006/10/19 18:32:03  pwagner
# First commit
#
