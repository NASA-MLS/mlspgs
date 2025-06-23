#!/bin/sh
# Generate the key=value file from the attributes of the hdf data file.
#    Some parameters will be hardcoded because they are not the hdf data file.
#    This version is specified for UARS level 1B
#

# Copyright 2014, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

#------------------------------- zero_pad ------------
#
# Function to pad a number with a zero if its
# length is less than 2
# usage: zero_pad n

zero_pad()
{
   # Count chars in n (-n  means skip newline)
   if [ `echo -n "$1" | wc -m` -lt 2 ]
   then
     echo "0$1"
   else
     echo $1
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

# ------------------ ccsds ascii date formats ------------------
# There aree 2 possible date formats we may impose on the metadata
# (a) yyyy-mm-dd
# (b) yyyy-Doy
#
# Let the next line assign either 'a' or 'b' as the format you prefer
dateformat='b'
# ----------------------------------------------------------------

# Assume the 1st argument is the hdf filename
filein=$1

# Generate all the known attributes from L1b
echo "LocalGranuleID=$filein"
l2auxdump -a PGEVersion  $filein | perl -pe 's/Dump of (.*)\n/$1=/'
if [ "$dateformat" = "b" ]
then
  # The file-level attributes were already formatted like 'b'
  l2auxdump -a StartUTC    $filein | perl -ne 'print if $. == 2' | perl -pe 's/(.*)T(.*)\n/RangeBeginningDate=$1\nRangeBeginningTime=$2.00000\nEquatorCrossingLongitude=0.0\nEquatorCrossingDate=$1\nEquatorCrossingTime=$2\n/' 
  l2auxdump -a EndUTC      $filein | perl -ne 'print if $. == 2' | perl -pe 's/(.*)T(.*)\n/RangeEndingDate=$1\nRangeEndingTime=$2.00000\n/' 
else
  GranuleYear=`l2auxdump -a GranuleYear $filein | awk '{print $3}'`
  GranuleMonth=`l2auxdump -a GranuleMonth $filein | awk '{print $3}'`
  GranuleDay=`l2auxdump -a GranuleDay $filein | awk '{print $3}'`
  # In case they are single digits, we pad them with zeros
  GranuleMonth=`zero_pad $GranuleMonth`
  GranuleDay=`zero_pad $GranuleDay`
  StartUTC=`l2auxdump -a StartUTC    $filein`
  EndUTC=`l2auxdump -a EndUTC    $filein`
  StartTime=`echo $StartUTC | awk '{print $4}' | awk -FT '{print $2}'`
  EndTime=`echo $EndUTC | awk '{print $4}' | awk -FT '{print $2}'`
  echo "RangeBeginningDate=$GranuleYear-$GranuleMonth-$GranuleDay"
  echo "RangeBeginningTime=$StartTime"
  echo "EquatorCrossingLongitude=0.0"
  echo "EquatorCrossingDate=0.0"
  echo "EquatorCrossingTime=$StartTime"
  echo "RangeEndingDate=$GranuleYear-$GranuleMonth-$GranuleDay"
  echo "RangeEndingTime=$EndTime"
fi
I=gen_uars_l1b_odl_txt
split_path="`echo $0 | sed 's/'$I'/split_path/'`"

# Guess from the filename 
# Must split path/filename
filename=`$split_path -f $filein`
filepath=`$split_path -p $filein`

LocalVersionID=`echo $filename | perl -nle 'print $1 if /-([a-z0-9]+)_/' `
#echo $filename | perl -nle 'print $1 if /-([a-z0-9]+)_/' 

# L1BOA or L1BRAD?
L1BOA=`echo $filename | grep L1BOA`
if [ "$L1BOA" != "" ]
then
  Identifier_Product_DOI=10.5067/AURA/MLS/DATA1001
  AAII=MUSUL1BOA_V1.html
else
  Identifier_Product_DOI=10.5067/AURA/MLS/DATA1002
  AAII=MUSUL1BRAD_V1.html
fi
PGEVersion=`echo $filename | awk -F_ '{print $3}' | awk -F"-" '{print $1"-"$2}'`
ProductionDateTime=`date +"%C%y-%m%d %T"`

# Give default values to the "stuttered" attributes
TTIITTLLEE=$filename
MMIIIIDD=$filename
AABBSSTTRRAACCTT=$filename
PPUURRPPOOSSEE=$filename
FFIILLEEIIDD=$filename
DDSS=$filename
MMUUSSTTAARRDD=$filename
DDQQDD=$filename
EEXXTTEENNTT=$filename

# May reset the stuttered values (or any others for that matter)
if [ -f "$filepath/attributes.txt" ]
then
  . "$filepath/attributes.txt"
fi

if [ -f "$filepath/abstract.txt" ]
then
  AABBSSTTRRAACCTT=`cat $filepath/abstract.txt`
  AABBSSTTRRAACCTT=`echo $AABBSSTTRRAACCTT`
fi

cat << EOF
StartOrbitNumber=-9999
StopOrbitNumber=-9999
LocalVersionID=$LocalVersionID
PGEVersion=$PGEVersion
ProductionDateTime=$ProductionDateTime
EOF

# Hardcode for L1b for missing additional information
cat  << EOF
DayNightFlag=Both
QAPercentMissingData=0
QAPercentOutofBoundsData=0
QAPercentInterpolatedData=0
ParameterName=Orbit/attitude and tangent point
ReprocessingPlanned=further update anticipated using enhanced PGE
InputPointer=Found at /PCF
VerticalSpatialDomainType=Atmosphere Layer
VerticalSpatialDomainValue=Brightness Temperature
ZoneIdentifier=Other Grid System
WestBoundingCoordinate=-180.0
NorthBoundingCoordinate=90.0
EastBoundingCoordinate=180.0
SouthBoundingCoordinate=-90.0
LocalityValue=Limb
DataProducer=MLS_SIPS
Identifier_Product_DOI=$Identifier_Product_DOI
title=$TTIITTLLEE
MI_Identifier=$MMIIIIDD
edition=v01-00
abstract=$AABBSSTTRRAACCTT
purpose=$PPUURRPPOOSSEE
MI_AcquisitionInformation=http://disc.sci.gsfc.nasa.gov/datacollection/$AAII
MI_CoverageDescription=+-180 degrees longitude, +-80 degrees latitude, day and night
MD_Keywords=atmospheric chemistry, climate
MI_Operation=Raytheon Corp.
CI_ResponsibleParty=GES-DAAC; another DAAC?
MD_GridSpatialRepresentation=grid
fileIdentifier=$FFIILLEEIIDD
DS_Series=$DDSS
LE_Lineage=Algorithm Theoretical Basis Document : $MMUUSSTTAARRDD
DQ_Element=MUSTARD Data Quality Document $DDQQDD
MD_BrowseGraphic=none
EX_Extent=$EEXXTTEENNTT
EOF



# $Log$
# Revision 1.4  2015/04/16 22:33:32  pwagner
# Reduce number of warnings logged (but why are there any?)
#
# Revision 1.3  2015/01/06 18:08:53  pwagner
# Revised default attribute values; added abstract.txt and attributes.txt mechanisms
#
# Revision 1.2  2014/11/25 01:06:28  pwagner
# Correctly set missing or erroneous metadata values
#
# Revision 1.1  2014/06/12 22:54:49  quyen
# tool to extract attributes from UARS L1B file
#
