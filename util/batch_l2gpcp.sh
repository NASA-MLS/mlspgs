#!/bin/sh
#batch_l2gpcp.sh

# --------------- batch_l2gpcp.sh help
# Performs a l2gpcat on a DGG file splitting it
# into standard product files for various species
# saving the results to one file per species
# Usage:
# batch_l2gpcp.sh [opt1] [opt2] ..  str1 [str2] .. - file
# where each string is a valid species name; e.g., H2O
# and file is the DGG file
#       * Special value expanded into all standard species *
# If str1 is StdProd it is expanded into the current standard products
# i.e.,  
# BrO CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI SO2 Temperature
# If str1 is APriori it is expanded into the current standard products but
# the swaths copied from the DGG file are the a priori ones, not the std prods
#
#     O p t i o n s
#    -dryrun              Merely echo the commands that would be executed
#    -profiles m n        Keep only profiles in range m n
#    -l2gpcat  command    Use command instead of l2gpcat
#  (and all the normal l2gpcat options; e.g. -v are wise choices)
#
# Note:
# (1) Each option is preceded by a dash
# (2) The options must come before the strings
# (3) The bare "-", must be present 
#      (to separate the strings from the file name on the command line)
#
# Bugs and limitations
# (1) l2gpcat is assumed to exist, to be an executable, and to be either in
#     your PATH
#     $MLSTOOLS/l2gpcat
# (2) The l2gp file names are assumed to match the pattern
#      MLS-Aura_L2GP-xxxx_*.he5
# (3) If multiple matches, we try always to pick out the last
#      in alphabetical order (which may or may not have the highest cycle num:
#      which is why we have the -c1 and -c2)
# --------------- End batch_l2gpcp.sh help
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
stdprods='BrO CH3Cl CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI SO2 Temperature'
debug=1
#     ^  -- set this to 1 if worried
keep=1
#    ^  -- set this to 1 to keep temp files (else delete)
me="$0"
my_name=batch_l2gpcp.sh
I=batch_l2gpcp
L2GPCAT=`which l2gpcat 2>/dev/null`
if [ ! -x "$L2GPCAT" ]
then
  L2GPCAT=$MLSTOOLS/l2gpcat
fi
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
l2gpcat_opts=""
list=""
profile1=""
profile2=""
apriori="no"
dryrun="no"
unique="no"
more_opts="yes"
more_strs="yes"
while [ "$more_strs" = "yes" ] ; do

    case "$1" in

    - )
	    shift
       more_opts="no"
       more_strs="no"
       ;;
    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -prof* )
	    profile1="$2"
	    shift
	    profile2="$2"
	    shift
	    shift
       ;;
    -l2gpcat )
	    L2GPCAT="$2"
	    shift
	    shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;
    -* )
       if [ "$l2gpcat_opts" = "" ]
       then
	       l2gpcat_opts="$1"
       else
	       l2gpcat_opts="$l2gpcat_opts $1"
       fi
       opt_takes_args=`echo "-d,-f,-chunks,-fields,-s1,-s2" | grep -e "$1"`
       # if [ "$1" = "-d" -o "$1" = "-f" ]
       echo opt_takes_args: $opt_takes_args
       if [ "$opt_takes_args" != "" ]
       then
	       l2gpcat_opts="$l2gpcat_opts $2"
          shift
       fi
	    shift
       ;;
    * )
       more_opts="no"
       aptest=`echo $1 | grep -i apr`
       sptest=`echo $1 | grep -i std`
       if [ "$sptest" != "" ]
       then
         list="$stdprods"
       elif [ "$aptest" != "" ]
       then
         list="$stdprods"
         apriori="yes"
       elif [ "$list" != "" ]
       then
         list="$list $1"
       else
         list="$1"
       fi
	    shift
       ;;
    esac
done

if [ "$list" = "" ]
then
  sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
      | sed -n 's/^.//p' | sed '1 d; $ d'
  exit
fi
dgg="$1"
if [ "$debug" = 1 ]
then
  echo "L2GPCAT $L2GPCAT"
  echo "l2gpcat_opts $l2gpcat_opts"
  echo "list $list"
  echo "unique $unique"
  echo "profile1 $profile1"
  echo "profile2 $profile2"
  echo "dgg $dgg"
fi
if [ ! -x "$L2GPCAT" ]
then
  echo "Sorry--l2gpcat must exist and be executable"
fi
for i in $list
do
  file1=`echo $dgg | sed -n "s/DGG/$i/ p"`
  if [ "$apriori" = "yes" ]
  then
    tail="APriori"
  else
    tail="StdProd"
  fi
  oldswath1=${i}-${tail}
  oldswath2="${i}-${tail} column"
  newswath1=${i}
  newswath2="${i} column"
  case "$i" in
  IWC )
    oldswaths="$oldswath1,IWP-${tail}"
    newswaths="$newswath1,IWP"
    ;;
  O3 )
    oldswaths="O3-${tail} column-MLS,O3-${tail} column-GEOS5,O3-${tail}"
    newswaths="O3 column-MLS,O3 column-GEOS5,O3"
    ;;
  Temperature )
    oldswaths="$oldswath1,WMOTPPressure-${tail}-MLS,WMOTPPressure-${tail}-GEOS5"
    newswaths="$newswath1,WMOTPPressure-MLS,WMOTPPressure-GEOS5"
    ;;
  RHI )
    oldswaths="$oldswath1"
    newswaths="$newswath1"
    ;;
  H2O )
    oldswaths="$oldswath1"
    newswaths="$newswath1"
    if [ "$apriori" = "yes" ]
    then
      oldswaths=H2O_LR-APriori
    fi
    ;;
  * )
    oldswaths="$oldswath1"
    newswaths="$newswath1"
    ;;
  esac
  
  if [ "$dryrun" = "yes" ]
  then
    echo "$L2GPCAT $l2gpcat_opts -s $oldswaths -r $newswaths -o $file1 $dgg"
  elif [ "$l2gpcat_opts" = "" ]
  then
    $L2GPCAT $l2gpcat_opts -s "$oldswaths" -r "$newswaths" -o $file1 $dgg
  else
    $L2GPCAT -s "$oldswaths" -r "$newswaths" -o $file1 $dgg
  fi
done
exit
# $Log$
# Revision 1.7  2013/05/30 20:37:03  pwagner
# Corrected erroneous comments about default location of l2gpcat
#
# Revision 1.6  2009/12/10 18:52:58  pwagner
# Attempts to find command first in PATH, then inder MLSTOOLS
#
# Revision 1.5  2009/08/18 21:05:09  pwagner
# Added new CH3Cl stdprod
#
# Revision 1.4  2008/10/13 23:30:56  pwagner
# Fixed typo
#
# Revision 1.3  2007/03/06 21:31:06  pwagner
# Added magic APriori string
#
# Revision 1.2  2007/02/06 23:17:31  pwagner
# Swath names fitted to v2.21
#
# Revision 1.1  2006/09/29 00:32:41  pwagner
# First commit
#
