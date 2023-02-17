#!/bin/sh
# --------------- repairDMP.sh help
# repairDMP.sh
# Repairs the DMP files for each day in a given year 
# (1) rename the existing l2GP-DMP file
# (2) insertl2gpvalues creating the new DMP with the original name
#     using dataset values from the original DMP file
# for the swaths
#   Altitude
#   DTDZ
#   DynTropAltitude
#   DynTropDTDZ
#   DynTropFlag
#   DynTropPressure
#   DynTropStaticStability
#   DynTropTemperature
#   DynTropTheta
#   GeopotentialHeight
#   GlobalAverageOfPVGrad
#   HorizontalPVGradient
#   HorizontalTemperatureGradient
#   MeridionalCompOfPVGrad
#   MeridionalWind
#   MontgomeryStreamFunction
#   PVEquivalentLatitude
#   PotentialVorticity
#   RelativeVorticity
#   ScaledPV
#   StaticStability
#   Temperature
#   Theta
#   VortexEdgeCenter
#   VortexEdgeInner
#   VortexEdgeOuter
#   WMOTropAltitude
#   WMOTropPV
#   WMOTropPressure
#   WMOTropStaticStability
#   WMOTropTemperature
#   WMOTropTheta
#   ZonalWind

#
# usage: repairDMP.sh [options] [years]

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -silent       Prevent most output during execution
#    -b command    Use command instead of insertl2gpvalues
#    -e command    Use command (e.g., h5edit) to add an AddedNotes attribute
#    -m "whatever" Modify value of AddedNotes attribute to "whatever"
#    -Ad Attr_dir  Root directory for legitimate L2GP files, e.g. L2GP-CO
#    -Dd DMP_dir   Root directory for DMP files
#    -Af Attr_file Legitimate L2GP file, e.g. L2GP-CO
#    -Df DMP_file  Existing DMP file
#    -Ef env_file  Get options by sourcing env_file
#    -undo         Undo effects of a previous repair, restoring the old files

# Effect:
# The DMP file will be replaced with a genuine hdfeos formatted DMP file

# Bugs and limitations
# (1) insertl2gpvalues is assumed to exist, to be an executable, and to be in
#     your path; otherwise use the -b option
# (2) the h5edit command is assumed to  be an executable if you use the -e option; 
#     e.g., /software/toolkit/hdftools/h5edit
# (3) the backup versions of the original files still exist and were
#     not moved or renamed if you use the -undo option
# (4) If -Ad and -Dd are both left undefined, you must
#     supply -Af and -Df instead. Perhaps most useful when testing.
# (5) If the years args are supplied: 
#     Attr_dir and DMP_dir include only paths up to
#     the version string; e.g.,
#       /data/emls/dmp/l2edmp/v05.01
#     and you will repair all the files for each day in all the years supplied.
# (6) If, however, the years args are absent:
#     Attr_dir and DMP_dir must include the year, too; e.g.
#       /data/emls/dmp/l2edmp/v05.01/2023
#     and you will repair only the files for each day in that year.

# --------------- End repairDMP.sh help
# Copyright 2023, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

#---------------------------- logit
# Either echo input string(s), or else append them to the master log file, 
# or both
# Which to do depends on $MUSTLOG

logit()
{
   if [ "$dryrun" = "yes" ]
   then
     echo "$@"
   else
     case "$MUSTLOG" in
     no )
       echo "$@"
       ;;
     yes )
       echo "$@" >> "$MASTERLOG"
       ;;
     both )
       echo "$@" | tee -a "$MASTERLOG"
       ;;
     esac
   fi
}

#---------------------------- bogify
# Create a file name with a bogus file extension in place of its current one
# i.e., file.he5 -> file.bogus

bogify()
{
   file=$1
#    nn=`echo $file | sed 's/.*c//; s/_.*//'`
#    bogusnn=00
#    echo $file | sed "s/c$nn/c$bogusnn/"
   echo $file | sed 's/\.he5/.bogus/'
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
     logit $@
   else
     $@
   fi
}
      
#---------------------------- undoit
# Undo the effects of a previous repair by restoring the old DMP.

undoit()
{
      commando mv $bogusDMP $DMP
}
      
#---------------------------- applyit
# Apply the repair by using va;ues from the old DMP and attributes
# copied from the legitimate L2GP-CO.

applyit()
{
    # Have we repaired this already?
#   echo "DMP $DMP"
#   echo "bogusDMP $bogusDMP"
#   echo "AttrFile $AttrFile"
    if [ -f "$bogusDMP" ]
    then
      echo "We already fixed this DMP $DMP"
    else
      commando mv $DMP $bogusDMP
      # Loop over all the swaths
      for swath in $swaths
      do
        # Now we do the overwriting: copy the DMP-190 swath values
        commando $insertl2gpvalues $insertoptions -Af $AttrFile -Vf $bogusDMP \
          -d $swath -dmp $DMP
      done
      # Have we been asked to add an extra Global Attribute
      # stating that we have monkeyed with the DMP file?
      # it loses the double quotes surrounding the inline command
      # A recourse would be to set the --command-file option to h5edit
      # and read in the single command
      if [ "$h5edit" != "" -a "$dryrun" != "yes" ]
      then
        $h5edit -c \
        "CREATE /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/AddedNotes {DATATYPE  H5T_STRING {STRSIZE 21} DATASPACE  SCALAR DATA {'Repaired faulty DMP file'}};"  \
        $DMP
      fi
    fi
}
      
# ------------------ repack ----------------------
# ----------------------------------------------------------------
# repack level 2 product files to speed things up
# Last chance to find h5repack
repack_files()
{
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$HDFTOOLS/h5repack
fi

if [ -x "$H5REPACK" ]
then
  for file in $DMP
  do
    if [ -w "$file" -o "$dryrun" = "yes" ]
    then
      packed="$file".p
      if [ "$GZIPLEVEL" != "" ] 
      then
        filter="-f GZIP=$GZIPLEVEL"
      else
        filter=""
      fi
      #logit "Packing $file into $packed"
      #logit $H5REPACK -i "$file" -o "$packed" $filter
      commando $H5REPACK -i "$file" -o "$packed" $filter
      # Here we could insert some check involving h5diff if we were dubious
      commando mv "$packed" "$file"
    else
      logit "No write permission to repack $file"
    fi
  done
fi
# ----------------------------------------------------------------
}

# ------------------ augment ----------------------
# ----------------------------------------------------------------
# augment level 2 product files to make them netcdf-compatible
# (must not reverse order of repack, augment)
# After we become more comfortable with l2gp2nc4 actually
# converting the hdfeos formats to NetCDF4 we may delete
# this section.
augment_files()
{
if [ -x "$NETCDFAUGMENT" ]
then
  for file in $DMP
  do
    if [ -w "$file" -o "$dryrun" = "yes" ]
    then
      # logit $NETCDFAUGMENT $file
      commando $NETCDFAUGMENT $file
    else
      logit "No write permission to augment $file"
    fi
  done
fi
}

# ------------------ one_day ----------------------
# ----------------------------------------------------------------
# Process each day in one day's worth of files.
one_day()
{
  # Repairing an entire day of files
  files=`/bin/ls $DMPDay`
  for file in $files
  do
    echo $file
    DMP=$DMPDay/$file
    bogusDMP=`bogify $DMP`
    # bogusDMP=$DMPDay/$bogusDMP
    AttrFile=$AttrDay/$AttrFile
    # Can't repair a non-existent DMP file
    if [ ! -f "$DMP" ]
    then
      echo "No standard product DMP file in $DMPYear/$day"
    # Have we been asked to undo an earlier repair?
    elif [ "$undo" = "yes" ]
    then
      undoit
    else
      applyit
      repack_files
      augment_files
    fi
  done
}

# ------------------ one_year ----------------------
# ----------------------------------------------------------------
# Process each day in one year's worth of files.
one_year()
{
  # Repairing an entire year of files
  days=`/bin/ls $DMPYear`
  for day in $days
  do
    if [ ! -d "$DMPYear/$day" ]
    then
      echo "Sorry--$DMPYear/$day is not a subdirectory containing DMP files"
    else
      # Get files with highest cycle numbers (except for DMP)
      thisDir=`pwd`
      DMPDay=$DMPYear/$day
      AttrDay=$AttrYear/$day
      cd $AttrDay
      AttrFile=`magnify L2GP-CO`
      # cd $DMPDay
      # DMP=`magnify DMP`
      # Now prefix with appropriate paths
      cd $thisDir
      # pwd
      # echo "DMPDay $DMPDay"
      one_day
    fi
  done
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

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

# Some annoying warnings come from h5edit unless 
# we set the following environment variable
HDF5_DISABLE_VERSION_CHECK=2
export HDF5_DISABLE_VERSION_CHECK

# But where is MLSTOOLS defined? Do you still rely on env variables?
# We should accept as a command line variable
# a file of settings using the
#   -Ef envfile
# mechanism we know and practice.
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=$MLSTOOLS/h5repack
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_hdfeos5
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=$MLSTOOLS/aug_eos5
fi

MUSTLOG="no"
MASTERLOG=/dev/null

# Last ditch definitions--good only for the scf filesystem
if [ ! -x "$H5REPACK" ]
then
  H5REPACK=/software/toolkit/hdftools/h5repack
fi
if [ ! -x "$NETCDFAUGMENT" ]
then
  NETCDFAUGMENT=/software/toolkit/mlstools/aug_eos5
fi

swaths="Altitude DTDZ DynTropAltitude DynTropDTDZ DynTropFlag DynTropPressure DynTropStaticStability DynTropTemperature DynTropTheta GeopotentialHeight GlobalAverageOfPVGrad HorizontalPVGradient HorizontalTemperatureGradient MeridionalCompOfPVGrad MeridionalWind MontgomeryStreamFunction PVEquivalentLatitude PotentialVorticity RelativeVorticity ScaledPV StaticStability Temperature Theta VortexEdgeCenter VortexEdgeInner VortexEdgeOuter WMOTropAltitude WMOTropPV WMOTropPressure WMOTropStaticStability WMOTropTemperature WMOTropTheta ZonalWind"

me="$0"
my_name=repairDMP.sh
I=repairDMP
h5edit=""
insertl2gpvalues=insertl2gpvalues
insertoptions="-v"
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
dryrun="no"
miscnotes=""
more_opts="yes"
undo="no"
AttrDir=""
DMPDir=""
AttrFile=""
DMP=""
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    dryrun="yes"
	    shift
       ;;
    -silent )
	    insertoptions="-silent"
	    shift
       ;;
    -undo )
	    undo="yes"
	    shift
       ;;
    -b )
	    insertl2gpvalues="$2"
	    shift
	    shift
       ;;
    -e )
	    h5edit="$2"
	    shift
	    shift
       ;;
    -m )
	    miscnotes="$2"
	    shift
	    shift
       ;;
    -Ad )
	    AttrDir="$2"
	    shift
	    shift
       ;;
    -Dd )
	    DMPDir="$2"
	    shift
	    shift
       ;;
    -Af )
	    AttrFile="$2"
	    shift
	    shift
       ;;
    -Df )
	    DMP="$2"
	    shift
	    shift
       ;;
    -Ef )
	    . "$2"
	    shift
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

years="$@"
if [ "$debug" = 1 ]
then
  echo "dryrun $dryrun"
  echo "undo $undo"
  echo "h5edit $h5edit"
  echo "insertl2gpvalues $insertl2gpvalues"
  echo "AttrDir $AttrDir"
  echo "DMPDir $DMPDir"
  echo "AttrFile $AttrFile"
  echo "DMP file $DMP"
  echo "miscnotes $miscnotes"
  echo "H5REPACK $H5REPACK"
  echo "NETCDFAUGMENT $NETCDFAUGMENT"
  echo "years $years"
fi

#echo "Past (1), anyway"
if [ "$DMPDir" != "" -a ! -d "$DMPDir" ]
then
  echo "Sorry--DMPDir must be a directory"
  exit
elif [ "$AttrDir" != "" -a ! -d "$AttrDir" ]
then
  echo "Sorry--AttrDir must be a directory"
  exit
elif [ ! -x "$insertl2gpvalues" ]
then
  echo "Sorry--insertl2gpvalues must exist and be executable"
  exit
fi

#echo "Past (2), anyway"
if [ "$AttrDir" = "" -a "$DMPDir" = "" ]
then
  # Repairing a single file
  bogusDMP=`bogify $DMP`
  echo "Repairing $DMP"
  echo "AttrFile $AttrFile"
  # Can't repair a non-existent DMP file
  if [ ! -f "$DMP" ]
  then
    echo "No standard product DMP file defined; did you use -Df option?"
    exit
  elif [ "$undo" = "yes" ]
  then
    undoit
  else
    applyit
    repack_files
    augment_files
  fi
elif [ "$years" = "" ]
then
  # DMPDir and AttrDir contain the year datum;
  # you're doing just one year
  DMPYear=$DMPDir
  AttrYear=$AttrDir
  one_year
else
  # DMPDir and AttrDir are jusst the version;
  # you're doing $years
  for year in $years
  do
    logit "year     $year"
    DMPYear=$DMPDir/$year
    AttrYear=$AttrDir/$year
    one_year
  done
fi
# $Log$
# Revision 1.1  2023/02/02 23:08:23  pwagner
# First commit
#
