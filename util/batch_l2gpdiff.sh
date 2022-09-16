#!/bin/sh
#batch_l2gpdiff.sh

# --------------- batch_l2gpdiff.sh help
# Performs a l2gpdiff on files from two separate runs for various species
# saving the results to one file per species
# Usage:
# batch_l2gpdiff.sh [opt1] [opt2] ..  str1 [str2] .. - dir1 dir2
# where each string is a valid species name; e.g., H2O
# and dir1, dir2 are the directories holding the l2gp files from the two runs
#       * Special value expanded into all standard species *
# If str1 is StdProd it is expanded into the current standard products
# i.e.,  
# BrO CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI Temperature
# If the -Ef file mechanism is used, and the variable PRODUCTS is defined
# its value will override the other species
#
#     O p t i o n s
#    -Ef file             set environment according to definitions found in file
#    -append              Append the results to an existing file
#                          (Otherwise replace any existing file)
#    -dryrun              Merely echo the commands that would be executed
#    -unique              Save results to a uniquely-named new file
#    -c1 cycle1           Insist on cycle number c"cycle1" in dir1
#                         (if cycle1 is "none" then file name may lack cycle number)
#    -c2 cycle2           Insist on cycle number c"cycle2" in dir2
#                         (if cycle2 is "none" then file name may lack cycle number)
#    -profiles m n        Keep only profiles in range m n
#    -l2gpdiff command    Use command instead of l2gpdiff
#  (and all the normal l2gpdiff options; e.g. -au -rms -ignore are wise choices)
#
# Note:
# (1) Each option is preceded by a dash
# (2) The options must come before the strings
# (3) The bare "-", must be present 
#      (to separate the strings from the directories on the command line)
# (4) The files being l2gpdiffed must be in different directories 
# (5) Each comparison will be saved in the current working directory as str.diff
#
# Bugs and limitations
# (1) l2gpdiff is assumed to exist, to be an executable, and (in desc order)
#     (a) defined by the env variables L2GPDIFF
#     (b) in your PATH
#     (c) /software/toolkit/mlstools/l2gpdiff
# (2) The l2gp file names are assumed to match the pattern
#      MLS-Aura_L2GP-xxxx_*.he5
# (3) If multiple matches, we choose the last
#      in alphabetical order (which may or may not have the highest cycle num:
#      that's why we have the -c1 and -c2)
# --------------- End batch_l2gpdiff.sh help
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
debug=0
#     ^  -- set this to 1 if worried
keep=0
#    ^  -- set this to 1 to keep temp files (else delete)
me="$0"
my_name=batch_l2gpdiff.sh
I=batch_l2gpdiff
if [ "$L2GPDIFF" = "" ]
then
  L2GPDIFF=`which l2gpdiff 2> /dev/null`
fi
if [ ! -x "$L2GPDIFF" ]
then
  L2GPDIFF=$MLSTOOLS/l2gpdiff
fi
# $reecho is reecho with me's path prepended
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
envfile="default-env.txt"
l2gpdiff_opts=""
list=""
PRODUCTS=""
append="no"
cycle1=""
cycle2=""
profile1=""
profile2=""
dryrun="no"
matchTimes="no"
unique="no"
more_opts="yes"
more_strs="yes"
while [ "$more_strs" = "yes" ] ; do

    if [ "$debug" = "1" ]
    then
      echo "arg $1"
    fi
    case "$1" in

    - )
	    shift
       more_opts="no"
       more_strs="no"
       ;;
    -append )
	    append="yes"
	    shift
       ;;
    -c1 )
	    cycle1="$2"
	    shift
	    shift
       ;;
    -c2 )
	    cycle2="$2"
	    shift
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
            if [ -f "$envfile" ]
            then
              . $envfile
            fi
       ;;
    -prof* )
	    profile1="$2"
	    shift
	    profile2="$2"
	    shift
	    shift
       ;;
    -unique )
	    unique="yes"
	    shift
       ;;
    -l2gpdiff )
	    L2GPDIFF="$2"
	    shift
	    shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;
    -* )
       if [ "$l2gpdiff_opts" = "" ]
       then
	       l2gpdiff_opts="$1"
       else
	       l2gpdiff_opts="$l2gpdiff_opts $1"
       fi
       opt_takes_args=`echo "-d,-f,-chunks,-fields,-pressures,-s1,-s2" | grep -e "$1"`
       # if [ "$1" = "-d" -o "$1" = "-f" ]
       if [ "$debug" = 1 ]
       then
         echo opt_takes_args: $opt_takes_args
       fi
       if [ "$opt_takes_args" != "" -a "$1" != "-s" ]
       then
	       l2gpdiff_opts="$l2gpdiff_opts $2"
          shift
       fi
	    shift
       ;;
    * )
       more_opts="no"
       sptest=`echo $1 | grep -i st`
       if [ "$sptest" != "" ]
       then
         list="$stdprods"
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

if [ "$PRODUCTS" != "" ]
then
  list="$PRODUCTS"
elif [ "$list" = "" ]
then
  sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
      | sed -n 's/^.//p' | sed '1 d; $ d'
  exit
fi
dir1="$1"
dir2="$2"
if [ "$profile1" != "" -a "$profile2" != "" ]
then
  frstdir=~/`get_unique_name frst`
  scnddir=~/`get_unique_name scnd`
  mkdir $frstdir
  mkdir $scnddir
fi
if [ "$debug" = 1 ]
then
  echo "L2GPDIFF $L2GPDIFF"
  echo "l2gpdiff_opts $l2gpdiff_opts"
  echo "list $list"
  echo "unique $unique"
  echo "append $append"
  echo "cycle1 $cycle1"
  echo "cycle2 $cycle2"
  echo "profile1 $profile1"
  echo "profile2 $profile2"
  echo "dir1 $dir1"
  echo "dir2 $dir2"
  echo "frstdir $frstdir"
  echo "scnddir $scnddir"
fi

if [ "$dir1" = "$dir2" ]
then
  echo "Sorry--dir1 and dir2 must be different"
elif [ ! -d "$dir1" ]
then
  echo "Sorry--dir1 must be a directory"
elif [ ! -d "$dir2" ]
then
  echo "Sorry--dir2 must also be a directory"
elif [ ! -x "$L2GPDIFF" ]
then
  echo "Sorry--L2GPDIFF must exist and be executable"
fi
for i in $list
do
  # We'll use this bit of trickery to skip any products that weren't created
  # e.g., "OH" on dates when the THz radiometer is off
  # of products not created by the nrt process
  file1="skip"
  file2="skip"
  extant_files $dir1/MLS-Aura_L2GP-${i}_*c${cycle1}*.he5
  if [ "$cycle1" = "none" ]
  then
    extant_files $dir1/MLS-Aura_L2GP-${i}_*.he5
  fi
  
  echo dir1/MLS-Aura_L2GP-{i}_ $dir1/MLS-Aura_L2GP-${i}_*.he5
  echo extant_files_result $extant_files_result

  nfiles=`echo "$extant_files_result" | wc | awk '{print $2}'`
  if [ "$nfiles" -gt 1 ]
  then
    file1=`echo $extant_files_result | sed 's/ /\n/g' | tail -1`
  elif [ "$nfiles" -lt 1 ]
  then
    echo "Sorry--no matching files for species $i in $dir1"
  else
    file1="$extant_files_result"
  fi

  extant_files $dir2/MLS-Aura_L2GP-${i}_*c${cycle2}*.he5
  if [ "$cycle2" = "none" ]
  then
    extant_files $dir2/MLS-Aura_L2GP-${i}_*.he5
  fi
  
  echo dir2/MLS-Aura_L2GP-{i}_ $dir2/MLS-Aura_L2GP-${i}_*.he5
  echo extant_files_result $extant_files_result

  nfiles=`echo "$extant_files_result" | wc | awk '{print $2}'`
  if [ "$nfiles" -gt 1 ]
  then
    file2=`echo $extant_files_result | sed 's/ /\n/g' | tail -1`
  elif [ "$nfiles" -lt 1 ]
  then
    echo "Sorry--no matching files for species $i in $dir2"
  else
    file2="$extant_files_result"
  fi

  if [ "$debug" = 1 ]
  then
    echo "species $i"
    echo "file1 $file1"
    echo "file2 $file2"
  fi
  a=`echo "$file1 $file2" | grep skip`
  # Run the diffs tool only if both files are present
  # i.e., neither is named "skip"
  if [ "$a" = "" ]
  then
    if [ "$profile1" != "" -a "$profile2" != "" ]
    then
      fname1=`perl -e '$reverse=reverse("$ARGV[0]"); @parts=split("/",$reverse); $reverse=$parts[0]; $reverse=reverse($reverse); print $reverse' $file1`
      fname2=`perl -e '$reverse=reverse("$ARGV[0]"); @parts=split("/",$reverse); $reverse=$parts[0]; $reverse=reverse($reverse); print $reverse' $file2`
      newfile1=$frstdir/$fname1
      newfile2=$scnddir/$fname2
      if [ "$dryrun" = "yes" ]
      then
        echo "l2gpcat -v -o $newfile1 -profiles $profile1 $profile2 $file1"
        echo "l2gpcat -v -o $newfile2 -profiles $profile1 $profile2 $file2"
      else
        l2gpcat -v -o $newfile1 -profiles $profile1 $profile2 $file1
        l2gpcat -v -o $newfile2 -profiles $profile1 $profile2 $file2
      fi
      file1=$newfile1
      file2=$newfile2
    fi
    if [ "$unique" = "yes" ]
    then
       outfile=`get_unique_name $i`
    else
      outfile="$i.diff"
    fi
    if [ "$dryrun" = "yes" ]
    then
      echo "l2gpdiff $file1 $file2 to $outfile"
    elif [ ! -f "$outfile" ]
    then
      echo "$file1 $file2 $outfile" > "$outfile"
    elif [ "$append" = "yes" ]
    then
       echo "$file1 $file2 $outfile" >> "$outfile"
    else
      rm -f "$outfile"
      echo "$file1 $file2 $outfile" > "$outfile"
    fi
    if [ "$file1" = "" ]
    then
      echo "No $i file found in $dir1"
    elif [ "$file2" = "" ]
    then
      echo "No $i file found in $dir2"
    elif [ "$dryrun" = "yes" ]
    then
      echo "$L2GPDIFF $l2gpdiff_opts $file1 $file2"
    elif [ "$l2gpdiff_opts" = "" ]
    then
      $L2GPDIFF $file1 $file2 >> "$outfile"
    else
      $L2GPDIFF $l2gpdiff_opts $file1 $file2 >> "$outfile"
    fi
  fi
done
if [ "$profile1" != "" -a "$profile2" != "" -a "$keep" != 1 ]
then
  /bin/rm -fr $frstdir
  /bin/rm -fr $scnddir
fi
exit
# $Log$
# Revision 1.19  2018/03/12 21:45:25  pwagner
# Respect the env variable L2GPDIFF as the preferred l2gpdiff
#
# Revision 1.18  2012/07/02 23:06:34  pwagner
# The environment variable PRODUCTS now determines what is diffed
#
# Revision 1.17  2012/04/20 23:16:19  pwagner
# Disabled more debugging messages
#
# Revision 1.16  2012/04/20 20:18:30  pwagner
# Turn off extra debugging
#
# Revision 1.15  2012/02/13 23:43:59  pwagner
# Moved where envfile gets sourced to better location
#
# Revision 1.14  2011/06/16 23:19:09  pwagner
# Fixed bug when finding L2GPDIFF
#
# Revision 1.13  2011/05/26 20:46:31  pwagner
# May use file of environment settings
#
# Revision 1.12  2009/12/10 18:52:58  pwagner
# Attempts to find command first in PATH, then inder MLSTOOLS
#
# Revision 1.11  2009/10/27 21:14:00  pwagner
# Corrected bug in options; l2gpdiff now assumed to be in PATH
#
# Revision 1.10  2009/08/18 21:04:52  pwagner
# Added new CH3Cl stdprod
#
# Revision 1.9  2007/10/09 18:24:05  pwagner
# Passes -pressures 'p1,p2,..' to l2gpdiff
#
# Revision 1.8  2007/06/01 16:48:58  pwagner
# Tiny changes, defaults to not keeping temp files
#
# Revision 1.7  2006/11/22 18:31:13  pwagner
# Special cycle number 'none' for file names lacking 'cnn' string
#
# Revision 1.6  2006/05/24 22:22:47  pwagner
# Added -profiles option
#
# Revision 1.5  2006/03/23 19:15:37  pwagner
# Better handling of multiple cycle numbers
#
# Revision 1.4  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.3  2005/04/01 00:06:47  pwagner
# StdProd species name expands to all stndard products
#
# Revision 1.2  2005/01/03 18:57:55  pwagner
# Fixed suntax errors in if--elif--fi sequence
#
# Revision 1.1  2004/12/29 22:26:03  pwagner
# First commit
#
