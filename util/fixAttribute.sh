#!/bin/sh
# --------------- fixAttribute.sh help
# fixAttribute.sh
# Checks the files in dir1 for granule, 
# (A) echo the InstrumentName global attribute and change it to "MLS UARS"
#     if the file name contains "MUSTARD-UMLS"
# (B) Correct the metatdata field named LocalGranuleID 
# (C) Fix PGEVersion attribute amd metadata
# (D) Apply h5repack and utilize compression
# (E) Apply aug_hdfeos5 if the file is an hdfeos
#
# fixAttribute.sh dir1

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -e command    Use command instead of h5edit
#    -f command    Use command instead of fixAttribute
#    -m path       Use path instead of MLSTOOLS
#    -v            operate verbosely

# Bugs and limitations
# (1) MLSTOOLS is assumed to an environment variable, and to hold the files
#     fixAttribute h5repack aug_hdfeos5
# (2) h5edit is assumed to exist, to be an executable, and to be equal to
#     the H5EDIT environment variable; or you may use -e option
# (3) You have write permission to dir1

# ~/mlspgs/util/fixAttribute.sh -v \
#   -e /software/toolkit/hdftools/h5edit \
#   -f /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test \
#   ~/mlspgs/tmp
# --------------- End fixAttribute.sh help
# Copyright 2018, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

#---------------------------- bogify
# Create a file name with a bogus cycle number in place of its current one
#

bogify()
{
   file=$1
   nn=`echo $file | sed 's/.*c//; s/_.*//'`
   bogusnn=00
   echo $file | sed "s/c$nn/c$bogusnn/"
}
      
#---------------------------- magnify
# Returns the file name with highest cycle number
#

magnify()
{
   echo *$1*.he5 | awk '{print $NF}'
}
      
#---------------------------- commando
# Either execute a command with all its args, or merely echo the command that
# would be executed
# Which to do depends on $dryrun

commando()
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
#	The entry point where control is given to the script    *
#****************************************************************
#
debug=1
#     ^  -- set this to 1 if worried

# Some annoying warnings come from H5EDIT unless 
# we set the following environment variable
HDF5_DISABLE_VERSION_CHECK=2
export HDF5_DISABLE_VERSION_CHECK

me="$0"
my_name=fixAttribute.sh
I=fixAttribute
dryrun="no"
verbose="no"

more_opts="yes"
undo="no"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -e )
	    H5EDIT="$2"
	    shift
	    shift
       ;;
    -f )
	    FIXER="$2"
	    shift
	    shift
       ;;
    -m )
	    MLSTOOLS="$2"
	    shift
	    shift
       ;;
    -v )
	    verbose="yes"
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

if [ "$#" -lt 1 ]
then
  echo "Sorry, you must supply a starting directory"
  echo "E.g., /data/mustard/umls/l2gp/vMUS01.50/1993"
  exit 1
fi

if [ "$MLSTOOLS" = "" ]
then
  MLSTOOLS=/software/toolkit/mlstools
fi

if [ ! -x "$FIXER" ]
then
  FIXER=$MLSTOOLS/fixAttribute
fi

H5REPACK=$MLSTOOLS/h5repack
NETCDFAUGMENT=$MLSTOOLS/aug_hdfeos5
GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)
if [ "$GZIPLEVEL" != "" ] 
then
  filter="-f GZIP=$GZIPLEVEL"
else
  filter=""
fi

DUMPER=`which l2auxdump`

if [ ! -x "$DUMPER" ]
then
  DUMPER=$MLSTOOLS/l2auxdump
fi

if [ ! -x "$FIXER" ]
then
  echo "Sorry--unable to find fixAttribute"
  exit 1
elif [ ! -x "$H5EDIT" ]
then
  echo "Sorry--unable to find h5edit"
  exit 1
elif [ ! -x "$H5REPACK" ]
then
  echo "Sorry--unable to find h5repack"
  exit 1
elif [ ! -x "$NETCDFAUGMENT" ]
then
  echo "Sorry--unable to find aug_hdfeos5"
  exit 1
fi

dir1="$1"

if [ "$debug" = 1 ]
then
  echo "dryrun $dryrun"
  echo "verbose $verbose"
  echo "H5REPACK $H5REPACK"
  echo "H5EDIT $H5EDIT"
  echo "NETCDFAUGMENT $NETCDFAUGMENT"
  echo "DUMPER $DUMPER"
  echo "FIXER $FIXER"
  echo "dir1 $dir1"
  # ls
fi

days=`/bin/ls $dir1`
#days=001
for day in $days
do
  if [ ! -d "$dir1/$day" ]
  then
    echo "Sorry--$dir1/$day is not a subdirectory containing emls or umls files"
  else
    cd $dir1/$day
    # pwd
    # ls
    TFiles=`/bin/ls MUSTARD-*MLS_L2*.h*5`
    # cd $JOBDIR
    for file in $TFiles
    do
      if [ "$verbose" = "yes" ]
      then
        echo $file
        if [ -x "$DUMPER" ]
        then
          # Is file a DGM or a DGG?
          a=`$DUMPER -r /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES -a InstrumentName $file 2> /dev/null`
          if [ "`echo $a | grep 'Unable to open'`" != "" ]
          then
            # DGM
            ftype='DGM'
            $DUMPER -a InstrumentName $file
          else
            # DGG
            ftype='DGG'
            echo $a
          fi
        fi
      fi
      # Is instrument umls or emls?
      # echo "$file"
      # echo "$file" | grep 'MUSTARD-UMLS'
      a=`echo "$file" | grep 'MUSTARD-UMLS'`
      if [ "$a" = "" ]
      then
        instrument='emls'
      else
        instrument='umls'
      fi
      echo "instrument $instrument"
      # exit 0
      
      #
      # Our agenda:
      # (0) Fix the InstrumentName and PGEVersion attributes
      # (1) Extract the 2 metadata datasets, core and xml, into text files
      # (2) Edit the text files as necessary
      # (3) Overwrite the 2 metadadata sets in the hdf file
      #
      #              (0)
      # Fix the InstrumentName attribute, but only if necessary
      if [ "$instrument" = 'umls' ]
      then
        $FIXER -A -a 'MLS UARS' $file
      fi
      
      vold="MUS01-50"
      vnew="MUS01-60"
      # Fix the PGEVersion attribute
      $FIXER -A -n PGEVersion -a V$vnew $file

      #              (1)
      $FIXER -X $file > $file.xml
      $FIXER -C $file > $file.met
      
      #              (2)
      # local version id contains the vold string
      # so does the PGEVersion
      # They must be changed to vnew
      # Do this first for the odl-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $file.met temp.met
      sed 's/'"$vold"'/'"$vnew"'/g' $file.met > temp.met
      mv temp.met $file.met
      
      # Do this next for the xml-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $file.xml temp.xml
      sed 's/'"$vold"'/'"$vnew"'/g' $file.xml > temp.xml
      mv temp.xml $file.xml
      
    
      #              (3)
      $FIXER -c $file.met $file
      $FIXER -x $file.xml $file
      # We are unable to use commando on the following command because
      # it loses the double quotes surrounding the inline command
      # A recourse would be to set the --command-file option to H5EDIT
      # and read in the single command
      if [ "$miscnotes" != "" -a "$dryrun" != "yes" ]
      then
        $H5EDIT -c \
        "CREATE /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/AddedNotes {DATATYPE  H5T_STRING {STRSIZE 21} DATASPACE  SCALAR DATA {'Repaired faulty swath, metadata, attribute'}};"  \
        $file
      fi
      # repack level 2 product files to speed things up
      if [ -x "$H5REPACK" ]
      then
        pfile="$file".p
        $H5REPACK -i "$file" -o "$pfile" $filter
        mv "$pfile" "$file"
      fi
      # augment level 2 product files to make them netcdf-compatible
      # (must not reverse order of repack, augment)
      if [ -x "$NETCDFAUGMENT" -a "$ftype" = 'DGG' ]
      then
        commando $NETCDFAUGMENT $file
      fi

    done
  fi
done

exit 0
# $Log$
