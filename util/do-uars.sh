#!/bin/sh
#Run the executable using the PCF file and whole toolkit panoply
# usage: do-uars.sh conv_uars_PGE genmet_PGE gen_l1b_odl_script PCF "
#          conv_uars_PGE          converting UARS to L1B program"
#          genmet_PGE             generate .met and .xml program"
#          gen_l1b_odl_script     script extract/fill the ODL text file"
#          PCF                    PCF file"
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

# Update PCF giving LID and filename
function replace_LID {

   filein=`basename $3`
   dir_filein=`dirname $3`

   # Search LID in PCF file
   looking_LID=`grep "^$2" $1`
   if [ -z "$looking_LID" ]; then
      #echo "Can not find the LID $2 in $1 ; going to add LID at the end"
      eval "perl -pi -e 's#\?   SUPPORT INPUT FILES#$2|$filein|$dir_filein||||1\n?   SUPPORT INPUT FILES#' $1"
   else
      eval "perl -pi -e 's#^(.*)\$#$2|$filein|$dir_filein||||1# if /^$2\|/ ' $1 "
   fi

}

# Read filename in PCF given LID
function readfname_LID {
   # Search LID in PCF file
   looking_LID=`grep "^$2" $1`
   if [ -z "$looking_LID" ] ; then
      #echo "Can not find the LID $2 in $1; return null"
      local fname=""
   else
      local fname=`echo $looking_LID | awk -F"|"  '{print $3 "/" $2}' `
   fi
   echo $fname
}

# Convert the binary-formatted uars data file to hdf 
# suitable for the mlspgs software
call_convert(){
   # Setup the toolkit environment
   . $PGSHOME/bin/linux/pgs-dev-env.ksh 

   # export variables for executing the conv_uars tool
   export FLIB_DVT_BUFFER=0
   export PGSMEM_USESHM=NO
   export PGE_BINARY_DIR
   export PGE_SCRIPT_DIR
   export PGE_ROOT
   $conv_uars_PGE $UARS_File -o $TEMP_DIR   > $CONV_UARS_LOG
}

# Generate the xml and xml files from executing the genmet routine
# Usage: call_genmet genmet PCF hd_file odl_txt mcf_file
call_genmet(){
   # Setup the toolkit environment
   . $PGSHOME/bin/linux/pgs-dev-env.ksh 

   # export variables for executing genmet tool
   export FLIB_DVT_BUFFER=0
   export PGSMEM_USESHM=NO
   export PGE_BINARY_DIR
   export PGE_SCRIPT_DIR
   export PGE_ROOT

   # Replace  LID 40000 with input hdf filename
   replace_LID $2 40000 $3 
   # Replace the LID 40001 with the odl text file
   replace_LID $2 40001 $4
   # Replace the LID 40002 with the MCF file
   replace_LID $2 40002 $5

   PGS_PC_Shell.sh $1 0111 $2 25 -v
}


# Main script specially for UARS files
echo "do-uars.sh $1 $2 $3 $4"
if [ "$#" -lt 4 ]
then
  echo "usage: do-uars.sh conv_uars_PGE genmet_PGE gen_l1b_odl_script PCF "
  echo "   conv_uars_PGE          converting UARS to L1B program"
  echo "   genmet_PGE             generate .met and .xml program"
  echo "   gen_l1b_odl_script     script extract/fill the ODL text file"
  echo "   PCF                    PCF file"
  exit 0
fi

# Assume in PCF,
#   LID 40010 UARS input file (as input)
#   LID 40011 L1BOA output file (as output directory)
#   LID 40012 L1BRAD output file (as output directory)
#   LID 40013 temporary working directory 
#
#   LID 40021 L1BOA  MCF file (fixed)
#   LID 40022 L1BRAD MCF file (fixed)
#
# Temporary LIDs
#   LID 40000 used in the genmet_PGE program
#   LID 40001 used in the genmet_PGE program
#   LID 40002 used in the genmet_PGE program
#

# Make the names to be easier to read
conv_uars_PGE=$1
genmet_PGE=$2
gen_uars_l1b_odl_script=$3
PCF=$4

# Extract the input UARS file name
UARS_File=$(readfname_LID $PCF 40010)

# Extract the output directory name
OUTPUT_L1BOA_File=$(readfname_LID $PCF 40011)
OUTPUT_L1BRAD_File=$(readfname_LID $PCF 40012)

# Extract the temporary working directory for conv_uars program
TEMP_LOC=$(readfname_LID $PCF 40013)

# Creating the temporary working directory (Append unique ID to make sure doesn't collide)
TEMP_DIR=${TEMP_LOC}_$$
if [ -d $TEMP_DIR ]; then
   # Clean out any h5 files
   /bin/rm -f ${TEMP_DIR}/*.h5
else
   # Create the temporary working directory
   /bin/mkdir -p $TEMP_DIR
fi
if [ ! -d $TEMP_DIR ]; then
  echo "ERROR: Can not access the working directory $TEMP_DIR"
  exit 1
fi


# Execute convert UARS program to generate L1BOA and L1BRAD
CONV_UARS_LOG=./conv_uars_log_$$.txt
/bin/rm -f $CONV_UARS_LOG

# The bare pge command silently failed for strange reasons,
# yielding product files with erroneous data
#$conv_uars_PGE $UARS_File -o $TEMP_DIR   > $CONV_UARS_LOG
#/users/pwagner/l1tests/conv_uars/do-testnopcf.sh $conv_uars_PGE $UARS_File -o $TEMP_DIR   > $CONV_UARS_LOG
call_convert
conv_stat=$?
if [ $conv_stat != 0 ] ; then
   echo "ERROR: the conv_uars program failed with error code is $conv_stat"
   echo "       check the output log file $CONV_UARS_LOG "
   exit 1
fi

# Update the PCF file with the correct L1OA filename and L1BRAD filename based on the conv_uars program
TEMP_L1BRAD_File=` grep '^rad file:' ./$CONV_UARS_LOG | awk '{print $3}'`
TEMP_L1BOA_File=` grep '^OA file:' ./$CONV_UARS_LOG | awk '{print $3}'`

echo "The conv_uars program generate two files:"
echo $TEMP_L1BOA_File
echo $TEMP_L1BRAD_File

# Move the output of the conv_uars program to the output specified in the PCF
/bin/mv -f $TEMP_L1BOA_File $OUTPUT_L1BOA_File
/bin/mv -f $TEMP_L1BRAD_File $OUTPUT_L1BRAD_File

# Execute genmet tool to generate .met and .xml files for L1BOA file
L1BOA=$(readfname_LID $PCF 40011)
ODL_L1BOA_File=./odl_input_L1BOA_$$.txt
MCF_L1BOA_File=$(readfname_LID $PCF 40021)
$gen_uars_l1b_odl_script $L1BOA > $ODL_L1BOA_File
call_genmet $genmet_PGE $PCF $L1BOA $ODL_L1BOA_File $MCF_L1BOA_File


# Execute genmet tool to generate .met and .xml files for L1BRAD file
L1BRAD=$(readfname_LID $PCF 40012)
ODL_L1BRAD_File=./odl_input_L1BRAD_$$.txt
MCF_L1BRAD_File=$(readfname_LID $PCF 40022)
$gen_uars_l1b_odl_script $L1BRAD > $ODL_L1BRAD_File
call_genmet $genmet_PGE $PCF $L1BRAD $ODL_L1BRAD_File $MCF_L1BRAD_File


# Cleanup the temporary file and directory
/bin/rm -f $ODL_L1BOA_File $ODL_L1BRAD_File $CONV_UARS_LOG
/bin/rm -rf $TEMP_DIR

# $Log$
# Revision 1.3  2015/05/27 17:41:38  pwagner
# Fixed obscure bug by setting environment variables
#
# Revision 1.2  2014/11/25 01:07:17  pwagner
#  Fixed errors in obtaining names of temporary l1b files
#
# Revision 1.1  2014/06/24 20:16:41  quyen
# tool to execute conv_uars program and generate the corresponding .met/.xml files
#
