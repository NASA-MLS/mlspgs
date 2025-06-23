#!/bin/sh
# --------------- repairDMP.sh help
# repairDMP.sh
# Repairs the DMP files for each day in a given year 
# (1) rename the existing l2GP-DMP file
# (2) insertl2gpvalues creating the new DMP with the original name
#     using dataset values from the original DMP file
# (3) repair the DOI value
# (4) write the .met and .xml files as hdf 
#
# Note that (3) does what the standalone script repairDOI.sh once did
# Note that (2) is done
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
#    -redo         Redo effects of a previous, presumably faulty, repair,
#                    assuming the .bogus files were left in place
#    -undo         Undo effects of a previous repair, restoring the old files
#    -refreesh     Refresh the contents of DMPDir from SourceDMPDir

# Effect:
# The DMP file will be replaced with a genuine hdfeos formatted DMP file

# Bugs and limitations
# (1) insertl2gpvalues is assumed to exist, to be an executable, and to be in
#     your path; otherwise use the -b option
# (2) the h5edit command is assumed to  be an executable if you use the -e option; 
#     e.g., /software/toolkit/hdftools/h5edit
# (3) the backup versions of the original files still exist and were
#     not moved or renamed if you use the -undo or -redo options
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
# Apply the repair by using values from the old DMP and attributes
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
        echo $insertl2gpvalues $insertoptions -Af $AttrFile -Vf $bogusDMP \
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
      # Now for my next trick, repair the DOI value(s)
      # Use the default DOIs which depend solely on the file name:
      # Does the file name contain the string MERRA?
      a=`echo "$DMP" | grep MERRA`
      if [ "$a" = "" ]
      then
        # non-MERRA
        commando $resetl2gpDOIs -v -a identifier_product_doi \
          -V "10.5067/AURA/MLS/DATA2525" "$DMP"
        LongName="MLS L2 EDMP"
        ShortName=ML2EDMP
        ShortDOI=DATA2525
      else
        # MERRA
        commando $resetl2gpDOIs -v -a identifier_product_doi \
          -V "10.5067/AURA/MLS/DATA2524" "$DMP"
        ShortName=ML2EDMP
        ShortDOI=DATA2524
        LongName="MLS L2 EDMP"
      fi
      # Last but not least, the metadata must be edited
      # Because we stopped archiving .met files we must extract it instead
      # cp $AttrFile.met $DMPmet
      $L2AUXDUMP -d "/HDFEOS INFORMATION/coremetadata.0" $AttrFile > $DMPmet
      cp $AttrFile.xml $DMPxml

      #              (1)
      vold=DATA2506
      vnew=$ShortDOI
      # DOI contains the vold string
      # They must be changed to vnew
      # Do this first for the odl-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet $DMPmet.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet > $DMPmet.temp
      mv $DMPmet.temp $DMPmet
      
      # Do this next for the xml-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml $DMPxml.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml > $DMPxml.temp
      mv $DMPxml.temp $DMPxml

      #              (2)
      vold=$AttrGranID
      vnew=$DMPGranID
      # local gran id contains the vold string
      # They must be changed to vnew
      # Do this first for the odl-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet $DMPmet.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet > $DMPmet.temp
      mv $DMPmet.temp $DMPmet
      
      # Do this next for the xml-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml $DMPxml.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml > $DMPxml.temp
      mv $DMPxml.temp $DMPxml
      

      #              (3)
      vold=ML2CO
      vnew=$ShortName
      # short name contains the vold string
      # They must be changed to vnew
      # Do this first for the odl-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet $DMPmet.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPmet > $DMPmet.temp
      mv $DMPmet.temp $DMPmet
      
      # Do this next for the xml-formatted metadata
      echo sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml $DMPxml.temp
      sed 's/'"$vold"'/'"$vnew"'/g' $DMPxml > $DMPxml.temp
      mv $DMPxml.temp $DMPxml

      #              (4)
      vnew=$LongName
      # These lines are not yet in the metadata files so they must be added
      # OBJECT                 = LONGNAME
      #   NUM_VAL              = 1
      #   VALUE                = "$LONGNAME"
      # END_OBJECT             = LONGNAME
      sed -n '1,/END_OBJECT             = SHORTNAME/ p' $DMPmet > $DMPmet.tmp.1  
      sed -n '/END_GROUP              = COLLECTIONDESCRIPTIONCLASS/,$ p' $DMPmet > $DMPmet.tmp.2  
      cat $DMPmet.tmp.1 > $DMPmet.tmp
      echo "" >> $DMPmet.tmp
      echo "    OBJECT                 = LONGNAME  " >> $DMPmet.tmp
      echo "      NUM_VAL              = 1          " >> $DMPmet.tmp
      echo "      VALUE                = "$vnew >> $DMPmet.tmp
      echo "    END_OBJECT             = LONGNAME  " >> $DMPmet.tmp
      echo "" >> $DMPmet.tmp
      cat $DMPmet.tmp.2 >> $DMPmet.tmp
      mv $DMPmet.temp $DMPmet

      # Do this next for the xml-formatted metadata
      # These lines are not yet in the metadata files so they must be added
      # <LongName>$LongName</LongName>
      sed -n '1,/<ShortName>/ p' $DMPxml > $DMPxml.tmp.1  
      sed -n '/<VersionID>/,$ p' $DMPxml > $DMPxml.tmp.2  
      cat $DMPxml.tmp.1 > $DMPxml.tmp
      echo "      <ShortName>$LongName</ShortName>" >> $DMPxml.tmp
      cat $DMPxml.tmp.2 >> $DMPxml.tmp
      mv $DMPxml.temp $DMPxml

      # Now insert the two kinds of metadata
      echo "Now insert the two kinds of metadata"
      echo $FIXER -c $DMPmet $DMP
      echo $FIXER -x $DMPxml $DMP
      $FIXER -c $DMPmet $DMP
      $FIXER -x $DMPxml $DMP
      
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
  AttrGranID=$AttrFile
  AttrFile=$AttrDay/$AttrFile
  if [ "$refresh" = "yes" ]
  then
    /bin/rm -f $DMPDay/*
    echo cp $SourceDay/* $DMPDay
    cp $SourceDay/* $DMPDay
  else
    # Repairing an entire day of files
    files=`/bin/ls $DMPDay`
    for file in $files
    do
      echo $file
      DMPGranID=$file
      DMP=$DMPDay/$file
      DMPmet=$DMP.met
      DMPxml=$DMP.xml
      bogusDMP=`bogify $DMP`
      # bogusDMP=$DMPDay/$bogusDMP
      a=`echo "$file" | grep '\.met'`
      b=`echo "$file" | grep '\.xml'`
      echo "a means .met $a"
      echo "b means .xml $b"
      if [ "$a" != "" ]
      then
        echo "metadata file will not be converted"
      elif [ "$b" != "" ]
      then
        echo "metadata file will not be converted"
      # Can't repair a non-existent DMP file
      elif [ ! -f "$DMP" ]
      then
        echo "No standard product DMP file in $DMPYear/$day"
      # Have we been asked to undo an earlier repair?
      elif [ "$redo" = "yes" ]
      then
        # redo means first undo, then do it over again
        # with whatever corrected tools are needed
        echo "First, undo the faulty conversion"
        undoit
        echo "Second, do the new conversion"
        applyit
        echo "Third, repack the files"
        repack_files
        echo "Fourth, augment the files"
        augment_files
      elif [ "$undo" = "yes" ]
      then
        undoit
      else
        applyit
        repack_files
        augment_files
      fi
    done
  fi
}

# ------------------ one_year ----------------------
# ----------------------------------------------------------------
# Process each day in one year's worth of files.
one_year()
{
  # Repairing an entire year of files
  days=`/bin/ls $DMPYear`
  #  days=038
  for day in $days
  do
    if [ ! -d "$DMPYear/$day" ]
    then
      echo "Sorry--$DMPYear/$day is not a subdirectory containing DMP files"
    else
      # Get files with highest cycle numbers (except for DMP)
      thisDir=`pwd`
      DMPDay=$DMPYear/$day
      SourceDay=$SourceYear/$day
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

if [ ! -x "$L2AUXDUMP" ]
then
  L2AUXDUMP=$MLSTOOLS/l2auxdump
fi

if [ ! -x "$FIXER" ]
then
  FIXER=$MLSTOOLS/fixAttribute
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
#swaths="Altitude DTDZ ZonalWind"
#swaths="Altitude DynTropDTDZ"

me="$0"
my_name=repairDMP.sh
I=repairDMP
h5edit=""
insertl2gpvalues=insertl2gpvalues
resetl2gpDOIs=resetl2gpDOIs
insertoptions="-v"
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
dryrun="no"
miscnotes=""
more_opts="yes"
redo="no"
refresh="no"
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
    -redo )
	    redo="yes"
	    shift
       ;;
    -refresh )
	    refresh="yes"
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
  echo "redo $redo"
  echo "refresh $refresh"
  echo "undo $undo"
  echo "h5edit $h5edit"
  echo "insertl2gpvalues $insertl2gpvalues"
  echo "FIXER $FIXER"
  echo "AttrDir $AttrDir"
  echo "DMPDir $DMPDir"
  echo "AttrFile $AttrFile"
  echo "DMP file $DMP"
  echo "miscnotes $miscnotes"
  echo "H5REPACK $H5REPACK"
  echo "L2AUXDUMP $L2AUXDUMP"
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
  elif [ "$redo" = "yes" ]
  then
    # redo means first undo, then do it over again
    # with whatever corrected tools are needed
    echo "First, undo the faulty conversion"
    undoit
    echo "Second, do the new conversion"
    applyit
    echo "Third, repack the files"
    repack_files
    echo "Fourth, augment the files"
    augment_files
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
  # DMPDir and AttrDir already contain the year datum;
  # you're doing just one year
  DMPYear=$DMPDir
  AttrYear=$AttrDir
  SourceYear=$SourceDMPDir
  one_year
else
  # DMPDir and AttrDir are jusst the version;
  # you're doing $years
  for year in $years
  do
    logit "year     $year"
    DMPYear=$DMPDir/$year
    AttrYear=$AttrDir/$year
    SourceYear=$SourceDMPDir/$year
    one_year
  done
fi
# $Log$
# Revision 1.5  2024/08/01 20:48:40  pwagner
# Copies metadata, too; adds function of repairing doi
#
# Revision 1.4  2024/07/19 16:43:22  pwagner
# Fixed bug with multiple files
#
# Revision 1.3  2023/05/08 21:57:02  pwagner
# Added -redo option
#
# Revision 1.2  2023/02/17 21:12:14  pwagner
# Now can repair every file in a days directory
#
# Revision 1.1  2023/02/02 23:08:23  pwagner
# First commit
#
