#!/bin/sh
# --------------- updatemf.sh help
#                 updatemf.sh    updates the Makefile in the MLSCONFG
#
# Called by make as a subtask of the update target
#
#        O p t i o n s
#  -h[elp]              Help (a work in progress)
#  -cf confg_file       configuration file (with path if needed)
#                        (defaults to ../.configure)
#  -pc plat_dir         Directory where platforms directory to be found
#                        (defaults to ./)
#  -machine             Check that machine files up-to-date
#  -reverse             Reverse hiding the file
#  -d                   Turn on debugging
#  -v                   Verbose
#
# --------------- End updatemf.sh help
#
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
#----------------------- EchoVerb -----------------------

# Function to echo its argument only if $verbose is true
#

EchoVerb()
{
   if [ "$verbose" = "true" -o "$debugging" = "true" ]
   then
     echo  "$1"
   fi
}

#----------------------- Exclude -----------------------

# Function to temporarily hide a file to exclude it from calculations
# or to reverse the hiding after the calculations
#

Exclude()
{
   if [ "$reverse" = "true" -a -f "$1-x" ]
   then
     #echo mv "$1-x" "$1"
     mv "$1-x" "$1"
   elif [  "$reverse" = "false" -a -f "$1" ]
   then
     #echo mv "$1" "$1-x"
     mv "$1" "$1-x"
   fi
}

#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in place of "."

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

#---------------------------- my_dot
#
# Function first replaces the rhs in assignment statements like
# a=b c d
# with
# a="b c d"
# so that they can be safely "dot"ted
# (otherwise the shell will complain: "c": command not found)
# (if a line is already of the form a="b c d", it will be left unchanged)
# use in place of "."

my_dot()
{
   temp1=`get_unique_name updatemf`
   sed 's/\(^[^ ][^ ]*\)=\(.*\)/\1="\2"/; s/\"\"\(.*\)\"\"/\"\1\"/' \
      "$1" > $temp1
   . ./"$temp1"
   rm "$temp1"
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

# Set the following to "on" to print degugging info
debugging="off"
me="$0"
me_simple=updatemf.sh
confg_file=../.configure
REECHO="`echo $0 | sed 's/updatemf\.sh/reecho.sh/'`"
MYPATH="`echo $0 | sed 's/updatemf\.sh//'`"
# The following preset some holdovers from when I was part of mlsconfigure
del_exist_dir="false"
mcfg_dir="./"
plat_dir="./"
conf_dir="$MYPATH/.."
m_update="false"
reverse="false"
verbose="false"

while [ "$1" != "" ] ; do

    case "$1" in

	-h | -help )
#	    Help
       sed -n '/'$me_simple' help/,/End '$me_simple' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	;;
	-cf )
	    confg_file=$2
	    shift
	;;
	-d )
	    debugging="true"
	;;
	-pc )
	    plat_dir=$2
	    shift
	;;
	-machine )
	    m_update="true"
	;;
	-reverse )
	    reverse="true"
	;;
	-v )
	    verbose="true"
	;;
    esac
    shift
done
if [ "$debugging" = "on" ] ; then
	echo "configuration file: $confg_file"
	echo "Reconfiguring Makefile in `pwd`"
fi
if [ ! -f "$confg_file" ]
then
	echo "configuration file $confg_file not found"
   exit 1
fi

my_dot "$confg_file"

machines_root=$conf_dir/lib/machines

# Have we been tasked with reversing a previous call where we hid
# one or more files temporarily via Exclude commands?
if [ "$reverse" = "true" ]
then
  if [ "$REECHOMACHOPTS" != "" -a "$MLSCONFG" = "$MLSF95.$MLSPLAT" ]
  then
    the_files=`$REECHO -dir $machines_root/"$MLSCONFG" -escape \*-x`
    if [ "$the_files" = "" ]
    then
      exit 0
    fi
    EchoVerb "About to restore $the_files"
    EchoVerb "in $machines_root/$MLSCONFG"
    for file in $the_files
    do
      argmx=`echo $file | sed 's/-x$//'`
      var=`Exclude $machines_root/$MLSCONFG/$argmx`
      # echo $var
    done
  fi
  exit 0
fi
EchoVerb "Re-create object-file directory for this $MLSCONFG"
EchoVerb "del_exist_dir $del_exist_dir"
#Re-create object-file directory for this MLSCONFG
#unless $del_exist_dir reset to false
#(Because mlsconfigure may be called simply to update $MLSCONFG/Makefile)
made_new_dir="false"

if [ "$MLSCONFG" = "." -o "$MLSCONFG" = "" ] ; then
	echo "MLSCONFG not yet defined--still $MLSCONFG"
	echo "To repair this, enter 'make configure'"
   exit 1
elif [ "$del_exist_dir" != "true" ] ; then
  if [ ! -d $mcfg_dir/$MLSCONFG ] ; then
      mkdir $mcfg_dir/$MLSCONFG
		made_new_dir="true"
   #else
      #rm -f $mcfg_dir/$MLSCONFG/Makefile
  fi
elif [ -d $mcfg_dir/$MLSCONFG ] ; then
	chmod u+w $mcfg_dir/$MLSCONFG/*
        if [ -f $mcfg_dir/$MLSCONFG/Makefile ]
        then
          EchoVerb "Preserving Makefile"
          ls -l $mcfg_dir/$MLSCONFG/Makefile
          mv $mcfg_dir/$MLSCONFG/Makefile $mcfg_dir/$MLSCONFG/.Makefile
	  rm -f $mcfg_dir/$MLSCONFG/*
          mv $mcfg_dir/$MLSCONFG/.Makefile $mcfg_dir/$MLSCONFG/Makefile
        else
          EchoVerb "Makefile not found; so creating it"
	  rm -f $mcfg_dir/$MLSCONFG/*
        fi
else
   mkdir $mcfg_dir/$MLSCONFG   
	made_new_dir="true"        
fi

if [ ! -f $conf_dir/srclib/$MLSF95.$MLSPLAT ] ; then
	echo "File $MLSF95.$MLSPLAT not found in $conf_dir/srclib"
	exit
fi

if [ ! -f $plat_dir/platforms/targets ] ; then
 echo "File targets not found in $plat_dir/platforms"
 exit
fi

# cat $confg_file $plat_dir/platforms/common \
#   $plat_dir/platforms/targets > $mcfg_dir/$MLSCONFG/Makefile

if [ -f $mcfg_dir/$MLSCONFG/Makefile ]
then
  mv $mcfg_dir/$MLSCONFG/Makefile $mcfg_dir/$MLSCONFG/Makefile.1
  cat $confg_file $plat_dir/platforms/targets > $mcfg_dir/$MLSCONFG/Makefile
  n=`diff $mcfg_dir/$MLSCONFG/Makefile $mcfg_dir/$MLSCONFG/Makefile.1 | wc -w`
  EchoVerb "n $n"
  if [ "$n" -lt 1 ]
  then
    mv $mcfg_dir/$MLSCONFG/Makefile.1 $mcfg_dir/$MLSCONFG/Makefile
  fi
else
  cat $confg_file $plat_dir/platforms/targets > $mcfg_dir/$MLSCONFG/Makefile
fi

if [ "$made_new_dir" = "true" ] ; then
	echo "Configured mls to use $mcfg_dir/$MLSCONFG "
	echo "(This is where your object files, .mod files, libraries,"
	echo "and executable programs will be created)"
	echo "To build the software, enter the 'make' command again"
	echo ""
	echo "If you later wish to rebuild for a different compiler"
	echo "or platform, enter 'make configure'"
   if [ "$MLSCONFG" != "$MLSF95.$MLSPLAT" ] ; then
     # Custom configuration name -- create a machines sub-directory for it
     # We will assume that the machines sub-directories live under the path
     # $conf_dir/lib/machines/$MLSCONFG
     if [ ! -d "$machines_root" ] ; then
        echo "Unable to locate root machines directory $machines_root"
        exit
     fi
     machines_base=$machines_root/$MLSF95.$MLSPLAT
     if [ ! -d "$machines_base" ] ; then
        echo "Unable to locate base machines directory $machines_base"
        exit
     fi
     if [ ! -d "$machines_root/$MLSCONFG" ] ; then
        mkdir "$machines_root/$MLSCONFG"
        cp $machines_base/*.f9* $machines_root/$MLSCONFG
        echo "New directory $machines_root/$MLSCONFG based on $machines_base"
     fi
   fi
fi

machines_root=$conf_dir/lib/machines            
machines_base=$machines_root/$MLSF95.$MLSPLAT   
if [ "$MLSCONFG" != "$MLSF95.$MLSPLAT" -a "$m_update" = "true" ] ; then
# In case of a custom configuration name, make sure machine files are up-to-date
  the_files=`$REECHO -dirn $machines_root/"$MLSF95.$MLSPLAT"`
  EchoVerb "About to make sure $the_files"
  EchoVerb "are up-to-date in $machines_root/$MLSCONFG"
  for file in $the_files
  do
    make -f $MYPATH/uptodate.make \
      SOURCE=$machines_base/$file TARGET=$machines_root/$MLSCONFG/$file \
        $machines_root/$MLSCONFG/$file
  done
fi

# In case you're using the REECHOMACHOPTS mechanism
# to exclude certain files from the machines subdirectory, basically
# rename them or remove them
# echo "m_update $m_update REECHOMACHOPTS $REECHOMACHOPTS"
if [ "$m_update" = "true" -a "$REECHOMACHOPTS" != "" ] ; then
  the_files=`$REECHO -dir $machines_root/$MLSCONFG fuzzy $REECHOMACHOPTS`
  if [ "$the_files" = "" ]
  then
    exit 0
  fi
  EchoVerb "About to rename $the_files"
  EchoVerb "in $machines_root/$MLSCONFG"
  for file in $the_files
  do
    var=`Exclude $machines_root/$MLSCONFG/$file`
    # echo $var
  done
fi


exit 0

# $Log$
# Revision 1.10  2014/01/29 21:03:03  pwagner
# Added -d and -v commandline opts
#
# Revision 1.9  2014/01/28 19:39:00  pwagner
# Create new Makefile only if different
#
# Revision 1.8  2013/06/05 18:42:40  pwagner
# Added ability to reverse hiding of Excluded files
#
# Revision 1.7  2009/12/01 21:32:47  pwagner
# May use the REECHOMACHOPTS mechanism to exclude certain files from machines
#
# Revision 1.6  2007/05/29 17:40:29  pwagner
# compiler-specific settings now in srclib, not platforms
#
# Revision 1.5  2005/06/23 22:20:46  pwagner
# Reworded Copyright statement
#
# Revision 1.4  2004/03/26 23:31:09  pwagner
# A quickie--but fails if we cvs remove from parent
#
# Revision 1.3  2002/09/24 18:10:28  pwagner
# Fixed -h option
#
# Revision 1.2  2002/09/04 18:11:36  pwagner
# my_dot repairs problem of multiple LIB_BLAS libraries
#
# Revision 1.1  2002/07/29 22:59:53  pwagner
# First commit
#
