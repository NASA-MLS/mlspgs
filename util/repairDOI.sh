#!/bin/sh
# --------------- repairDOI.sh help
# repairDOI.sh
# (Superseded by repairDMP.sh which does all this does and more now)
# Repairs the DOI attribute in DMP files for each day in a given year 
# * resetl2gpDOIs to the correct value

#
# usage: repairDOI.sh [options] [years]

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -silent       Prevent most output during execution
#    -b command    Use command instead of resetl2gpDOIs
#    -Dd DMP_dir   Root directory for DMP files
#    -Df DMP_file  Existing DMP file
#    -Ef env_file  Get options by sourcing env_file

# Effect:
# The DMP file's DOI will be replaced with its correct value
# i.e.,
#    "10.5067/AURA/MLS/DATA2523" for non-MERRA
#    "10.5067/AURA/MLS/DATA2524" for MERRA

# Bugs and limitations
# (1) resetl2gpDOIs is assumed to exist, to be an executable, and to be in
#     your path; otherwise use the -b option
# (2) If -Dd is left undefined, you must
#     supply -Df instead. Perhaps most useful when testing.
# (3) If the years args are supplied: 
#     DMP_dir should include only paths up to
#     the version string; e.g.,
#       /data/emls/dmp/l2edmp/v05.01
#     and you will repair all the files for each day in all the years supplied.
# (4) If, however, the years args are absent:
#     DMP_dir must include the year, too; e.g.
#       /data/emls/dmp/l2edmp/v05.01/2023
#     and you will repair only the files for each day in that year.

# --------------- End repairDOI.sh help
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
      
#---------------------------- applyit
# Either set the file attribute to a value when you supplied
# both its name and its value explicitly (probably through a .env file)
# or else do the default which is to correct any faulty DOI attribute values.
# Which to do depends on whether the default blank values for
# attrName and attrValue have been overridden.

applyit()
{
   if [ "$attrName" != "" -a "$attrValue" != "" ]
   then
     commando $resetl2gpDOIs -v -a "$attrName" -V "$attrValue" "$DMP"
   else
     # Use the default DOIs which depend solely on the file name:
     # Does the file name contain the string MERRA?
     a=`echo "$DMP" | grep MERRA`
     if [ "$a" = "" ]
     then
       # non-MERRA
       commando $resetl2gpDOIs -v -a identifier_product_doi \
         -V "10.5067/AURA/MLS/DATA2525" "$DMP"
     else
       # MERRA
       commando $resetl2gpDOIs -v -a identifier_product_doi \
         -V "10.5067/AURA/MLS/DATA2524" "$DMP"
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
    # bogusDMP=`bogify $DMP`
    # bogusDMP=$DMPDay/$bogusDMP
    # Can't repair a bugus DMP file
    a=`echo $file | grep '\.bogus'`
    if [ "$a" != "" ]
    then
      echo "Skipping bogus file $file in $DMPYear/$day"
    # Can't repair a non-existent DMP file
    elif [ ! -f "$DMP" ]
    then
      echo "No standard product DMP file in $DMPYear/$day"
    # Have we been asked to undo an earlier repair?
    # elif [ "$undo" = "yes" ]
    # then
    #  undoit
    else
      applyit
      # repack_files
      # augment_files
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
RUNANYWAY="no"
#           ^  -- set this to yes if you're really, really sure
if [ "$RUNANYWAY" != "yes" ]
then
  echo "repairDOI.sh has been superseded by repairDMP.sh"
  echo "so run that instead"
  echo "For more info contact Paul Wagner"
  exit 1
fi
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

if [ ! -x "$FIXER" ]
then
  FIXER=$MLSTOOLS/fixAttribute
fi

me="$0"
my_name=repairDOI.sh
I=repairDOI
h5edit=""
resetl2gpDOIs=resetl2gpDOIs
insertoptions="-v"
split_path="`echo $0 | sed 's/'$I'/split_path/'`"
dryrun="no"
miscnotes=""
more_opts="yes"
DMPDir=""
DMP=""
attrName=""
attrValue=""
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
    -b )
	    resetl2gpDOIs="$2"
	    shift
	    shift
       ;;
    -Dd )
	    DMPDir="$2"
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
  echo "resetl2gpDOIs $resetl2gpDOIs"
  echo "DMPDir $DMPDir"
  echo "DMP file $DMP"
#  echo "H5REPACK $H5REPACK"
#  echo "NETCDFAUGMENT $NETCDFAUGMENT"
  echo "years $years"
fi

#echo "Past (1), anyway"
if [ "$DMPDir" != "" -a ! -d "$DMPDir" ]
then
  echo "Sorry--DMPDir must be a directory"
  exit
elif [ ! -x "$resetl2gpDOIs" ]
then
  echo "Sorry--resetl2gpDOIs must exist and be executable"
  exit
fi

#echo "Past (2), anyway"
if [ "$DMPDir" = "" ]
then
  # Repairing a single file
  echo "Repairing $DMP"
  # Can't repair a non-existent DMP file
  if [ ! -f "$DMP" ]
  then
    echo "No standard product DMP file defined; did you use -Df option?"
    exit
  else
    applyit
    # repack_files
    # augment_files
  fi
elif [ "$years" = "" ]
then
  # DMPDir contain the year datum;
  # you're doing just one year
  DMPYear=$DMPDir
  one_year
else
  # DMPDir are just the version;
  # you're doing $years
  for year in $years
  do
    logit "year     $year"
    DMPYear=$DMPDir/$year
    one_year
  done
fi
# $Log$
# Revision 1.3  2024/07/19 16:45:48  pwagner
# Avoid collision between non-Merra and CloudTopPressure dois
#
# Revision 1.2  2023/03/30 21:35:09  pwagner
# Fixed many bugs
#
# Revision 1.1  2023/03/16 16:24:29  pwagner
# First commit
#
