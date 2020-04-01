#!/bin/sh
# Copyright 2015, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# $Id$

# A chameleon script:
# The same file exists as two different links; the action it
# performs depends on the name by which it is summoned
# In particular it can check for either of the following problems
# name                      problem it checks for
# zeros.sh               0s in the Pricision field
# misalignment.sh        misaligned geolocations

# --------------- zeros.sh help
# zeros.sh
# checks for 0s in the Precision field of ClO
#
# usage:
# zeros.sh [-Ef envsettings-file] [options] [files]
# zeros.sh [-Ef envsettings-file] [options] [years]
#
# where envsettings-file is a file of environmental settings
#
#    O p t i o n s
# -Ef file            set environment according to definitions found in file
# -dryrun             merely echo the command that would be executed
# -v                  verbose
# -h[elp]             print brief help message; exit
#
#
# will run l2gpdump to discover the dates when there 0 values in the Precision field
#
# You can use the envsettings file to change any of
# species             defaults to ClO
# version             defaults to v04
# root                defaults to /data/emls/l2gp
# --------------- End zeros.sh help

# --------------- misalignment.sh help
# misalignment.sh
# checks if the horizontal geolocations are misalgned
#
# usage:
# misalignment.sh [-Ef envsettings-file] [options] [files]
# misalignment.sh [-Ef envsettings-file] [options] [years]
#
# where envsettings-file is a file of environmental settings
#
#    O p t i o n s
# -Ef file            set environment according to definitions found in file
# -dryrun             merely echo the command that would be executed
# -v                  verbose
# -h[elp]             print brief help message; exit
#
#
# will run misalignment to discover the dates when it occurs
#
# You can use the envsettings file to change any of
# swath1              defaults to 'refGPH-InitPtan'
# swath2              defaults to 'O3-StdProd'
# version             defaults to v04
# root                defaults to /data/emls/l2gp
# --------------- End misalignment.sh help

#------------------------------- print_if_zeros ------------
print_if_zeros()
{
  if [ "$dryrun" = "yes" ]
  then
    echo "~/mlspgs/bin/IFC.Linux/l2gpdump -s $species -prec 0 -rel statuseven $1 | awk '{print $6}'"
  else
    a=`~/mlspgs/bin/IFC.Linux/l2gpdump -s $species -prec 0 -rel statuseven $1 | awk '{print $6}'`
    if [ "$a" != "0" ]
    then
      echo $1
    fi
  fi
}

#------------------------------- print_if_misaligned ------------
print_if_misaligned()
{
  for k
  do
    yyyyddd=`echo $k | awk -F_ '{print $4}' | sed 's/\.he5//'`
    yyyymmdd=`dateconverter $yyyyddd`
    if [ "$dryrun" = "yes" ]
    then
      echo "/software/toolkit/mlstools/misalignment -silent $k "
    else
      a=`/software/toolkit/mlstools/misalignment -silent $k`
      if [ "$a" != "" ]
      then
        echo "$k  ($yyyymmdd)"
      fi
    fi
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
cmdline=`echo $0 $@`
#my_name=zeros.sh
me=$0
root="/data/emls/l2gp"
# Were we called with a relative or absolute path prefix?
# i.e., does dollar-0 begin with a slash?
slashtest=`echo $0 | grep '^/.*'`
if [ "$slashtest" != "" ]
then
  I=`echo $0 | sed 's:/.*/:/:g' | sed 's:/::' | sed 's/\.sh//'`
else
  I=`echo $0 | sed 's:.*/::g' | sed 's/\.sh//'`
fi
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
my_name=` $the_splitter -f $me`
species=ClO
dryrun="no"
reuse="no"
verbose="no"
envfile="default-env.txt"
years=""
version=v04.2
swath1='refGPH-InitPtan'
swath2='O3-StdProd'
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    # echo "option: $1"
    case "$1" in

    -dryrun )
       dryrun="yes"
       shift
       ;;
    -Ef )
       shift
       envfile="$1"
       shift
       ;;
    -h | -help )
       echo "my_name $my_name"
       echo "me $me"
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       rm -f $settings_file
       exit
       ;;
    -v )
       verbose="yes"
       shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done
if [ -f "$envfile" ]
then
  . "$envfile"
fi

if [ "$verbose" = "yes" ]
then
  echo $cmdline
fi

if [ "$my_name" = "zeros.sh" ]
then
  command=print_if_zeros
else
  command=print_if_misaligned
  species=DGG
fi
# Do we loop over years? Or do we loop over file names?
# If the first arg is a file name, then presumably they all are
# If the 1st arg is not a filename, we interpret it as a year number, e.g. "2004"
if [ ! -f "$1" ]
then
  # Year number
  for v in $root/$version*
  do
    # echo $v
    for i
    do
      # echo $v/$i
      if [ -d $v/$i ]
      then
        cd $v/$i
        # echo $v/$i
        # exit
        for j in *
        do
          $command $j/*$species*.he5
        done
      fi
    done
  done
else
  # File name
  for i
  do
    $command $i
  done
fi
# $Log$
# Revision 1.1  2015/11/03 17:33:42  pwagner
# First commit
#
