#!/bin/sh
# --------------- batch_l2gpcat.sh help
# batch_l2gpcat.sh
# l2gpcat all the files in dir1, saving them in dir2
#
# batch_l2gpcat.sh [opt1] [opt2] ..  str1 [str2] .. - dir1 dir2
# where each string is a valid species name; e.g., H2O
# and dir1, dir2 are the directories holding the l2gp files from the two runs
#       * Special value expanded into all standard species *
# Special notes on species names:
# (1) If str1 is StdProd it is expanded into the current standard products
#   BrO CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI Temperature
# (2) If str1 is nrt2 it is expanded into the current nrt standard products
#   CO H2O HNO3 N2O O3 SO2 Temperature
# (3) If str1 is nrt it is expanded into the former nrt standard products
#   O3 Temperature
# (4) If str1 is L2FWM-.. it is expanded into the appropriate l2fwm files
# (5) If str1 is L2AUX-.. it is expanded into the appropriate l2aux files

#     O p t i o n s
#    -dryrun              Merely echo the commands that would be executed
#    -c cycle             Insist on cycle number c"cycle"
#    -l2gpcat command     Use command instead of l2gpcat

# Bugs and limitations
# (1) l2gpcat is assumed to exist, to be an executable, and to be in
#     $HOME/pvm3/bin/LINUX
# (2) The l2gp file names are assumed to match the pattern
#      MLS-Aura_L2GP-xxxx_*.he5 ; l2fwm, l2aux as indicated above
# (3) If multiple matches, we try always to pick out the last
#      in alphabetical order (which may or may not have the highest cycle num:
#      which is why we have the -c1 and -c2)

# usage: see (1) above

# --------------- End batch_l2gpcat.sh help
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

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
GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)


nrt2prods='CO H2O HNO3 N2O O3 SO2 Temperature'
nrtdprods='O3 Temperature'
stdprods='BrO CH3CN ClO CO GPH H2O HCl HCN HNO3 HO2 HOCl IWC N2O O3 OH RHI Temperature'
debug=1
#     ^  -- set this to 1 if worried
me="$0"
my_name=batch_l2gpcat.sh
I=batch_l2gpcat
l2gpcat=`which l2gpcat 2>/dev/null`
if [ ! -x "$l2gpcat" ]
then
  l2gpcat=$MLSTOOLS/l2gpcat
fi
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
l2gpcat_opts="-cat"
list=""
cycle=""
dryrun="no"
more_opts="yes"
more_strs="yes"
while [ "$more_strs" = "yes" ] ; do

    case "$1" in

    - )
	    shift
       more_opts="no"
       more_strs="no"
       ;;
    -c )
	    cycle="$2"
	    shift
	    shift
       ;;
    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -l2gpcat )
	    l2gpcat="$2"
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
       if [ "$1" = "-d" -o "$1" = "-f" ]
       then
	       l2gpcat_opts="$l2gpcat_opts $2"
          shift
       fi
	    shift
       ;;
    * )
       more_opts="no"
       nrt2test=`echo $1 | grep -i nrt2`
       nrttest=`echo $1 | grep -i nrt`
       sptest=`echo $1 | grep -i st`
       if [ "$nrt2test" != "" ]
       then
         list="nrt2prods"
       elif [ "$nrttest" != "" ]
       then
         list="$nrtdprods"
       elif [ "$sptest" != "" ]
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

if [ "$list" = "" ]
then
  sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
      | sed -n 's/^.//p' | sed '1 d; $ d'
  exit
fi
dir1="$1"
dir2="$2"

if [ "$debug" = 1 ]
then
  echo "l2gpcat $l2gpcat"
  echo "l2gpcat_opts $l2gpcat_opts"
  echo "list $list"
  echo "cycle $cycle"
  echo "dir1 $dir1"
  echo "dir2 $dir2"
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
elif [ ! -x "$l2gpcat" ]
then
  echo "Sorry--l2gpcat must exist and be executable"
  exit
fi

for i in $list
do
  l2aux=`echo $i | grep -i l2aux`
  l2fwm=`echo $i | grep -i l2fwm`
  # echo "l2aux $l2aux"
  # echo "l2fwm $l2fwm"
  if [ "$l2fwm" != "" ]
  then
    # echo $dir1/MLS-Aura_${i}_*c${cycle}*.h5
    extant_files $dir1/MLS-Aura_${i}_*c${cycle}*.h5
  elif [ "$l2aux" != "" ]
  then
    # echo $dir1/MLS-Aura_${i}_*c${cycle}*.h5
    extant_files $dir1/MLS-Aura_${i}_*c${cycle}*.h5
  else
    extant_files $dir1/MLS-Aura_L2GP-${i}_*c${cycle}*.he5
  fi
  nfiles=`echo "$extant_files_result" | wc | awk '{print $2}'`
  # echo "files $extant_files_result"
  # echo "nfiles $nfiles"
  if [ "$nfiles" -gt 1 ]
  then
    file1=`echo $extant_files_result | sed 's/ /\n/g' | tail -1`
  elif [ "$nfiles" -lt 1 ]
  then
    echo "Sorry--no matching files for species $i in $dir1"
  else
    file1="$extant_files_result"
  fi
  # file2=`echo "$file1" | sed 's:'$dir1':'$dir2':'`
  name1=`$split_path -f $file1`
  name2_head=`echo $name1 | cut -d_ -f 1,2,3`
  name2_tail=`echo $name1 | cut -d_ -f 4 | cut -c 1-8,14-17`
  name2=${name2_head}_${name2_tail}
  file2=$dir2/$name2

  if [ "$debug" = 1 ]
  then
    echo "species $i"
    echo "file1 $file1"
    echo "file2 $file2"
  fi
  if [ "$dryrun" = "yes" ]
  then
    echo $l2gpcat $l2gpcat_opts -o $file2 $extant_files_result
  else
    echo "cating $extant_files_result into $file2"
    rm -f $file2
    $l2gpcat $l2gpcat_opts -o $file2 $extant_files_result
  fi
done

# $Log$
# Revision 1.5  2012/04/05 20:11:07  pwagner
# nrt2 option cats the extra nrt standard prods from 2nd generation
#
# Revision 1.4  2009/12/10 18:52:58  pwagner
# Attempts to find command first in PATH, then inder MLSTOOLS
#
# Revision 1.3  2009/04/18 00:50:18  pwagner
# Reverted last buggy commit
#
# Revision 1.2  2009/04/13 20:43:17  pwagner
# Fixed a bug preventing macros file from using its own macros properly
#
# Revision 1.1  2008/02/28 01:37:24  pwagner
# 1st commit
#
