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

# Assume the 1st argument is the hdf filename
filein=$1

# Generate the known attributes from L1b
echo "LocalGranuleID=$filein"
l2auxdump -a PGEVersion  $filein | perl -pe 's/Dump of (.*)\n/$1=/'
l2auxdump -a StartUTC    $filein | perl -ne 'print if $. == 2' | perl -pe 's/(.*)T(.*)\n/RangeBeginningDate=$1\nRangeBeginningTime=$2.00000\nEquatorCrossingLongitude=0.0\nEquatorCrossingDate=$1\nEquatorCrossingTime=$2\n/' 
l2auxdump -a EndUTC      $filein | perl -ne 'print if $. == 2' | perl -pe 's/(.*)T(.*)\n/RangeEndingDate=$1\nRangeEndingTime=$2.00000\n/' 

# Guess from the filename 
#LocalVersionID=`echo $filein | perl -nle 'print $1 if /-([a-z0-9]+)_/' `
# Hardcode because LocalVersionID can not guess from the filename
LocalVersionID=c02
ProductionDateTime=`date +"%C%y-%m%d %T"`
cat << EOF
StartOrbitNumber=-9999
StopOrbitNumber=-9999
LocalVersionID=$LocalVersionID
ProductionDateTime=$ProductionDateTime
EOF

# Hardcode for L1b for missing additional information
cat  << EOF
DayNightFlag=Both
QAPercentMissingData=0
QAPercentOutofBoundsData=0
QAPercentInterpolatedData=0
AutomaticQualityFlag=Passed
AutomaticQualityFlagExplanation=Validated
OperationalQualityFlag=Not Investigated
OperationalQualityFlagExplanation=Not Investigated
ScienceQualityFlag=Inferred Passed
ScienceQualityFlagExplanation=Validated, see http://disc.gsfc.nasa.gov/Aura/MLS/ for quality document
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
title=TBD
MI_Identifier=TBD
edition=TBD
abstract=TBD
purpose=TBD
MI_AcquisitionInformation=TBD
MI_CoverageDescription=TBD
MD_Keywords=TBD
MI_Operation=TBD
CI_ResponsibleParty=TBD
MD_GridSpatialRepresentation=TBD
fileIdentifier=TBD
DS_Series=TBD
LE_Lineage=TBD
DQ_Element=TBD
MD_BrowseGraphic=TBD
EX_Extent=TBD
EOF



# $Log$
