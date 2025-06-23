#!/bin/sh
# --------------- fixAttribute_batch.sh help
# fixAttribute_batch.sh
# Checks the files in dir1 for granule, 
# (A) echo the InstrumentName global attribute and change it to "MLS UARS"
#     if the file name contains "MUSTARD-UMLS"
# (B) Correct the metatdata field named LocalGranuleID 
# (C) Fix PGEVersion attribute amd metadata
# (D) Apply h5repack and utilize compression
# (E) Apply aug_hdfeos5 if the file is an hdfeos
#
# fixAttribute_batch.sh dir1

#     O p t i o n s
#    -dryrun       Merely echo the commands that would be executed
#    -e command    Use command instead of h5edit
#    -f command    Use command instead of fixAttribute_batch
#    -m path       Use path instead of MLSTOOLS
#    -v            operate verbosely

# Bugs and limitations
# (1) MLSTOOLS is assumed to an environment variable, and to hold the files
#     fixAttribute_batch h5repack aug_hdfeos5
# (2) h5edit is assumed to exist, to be an executable, and to be equal to
#     the H5EDIT environment variable; or you may use -e option
# (3) You have write permission to dir1
#
# cd /nas/testing/workspace/pwagner/mustard/umls/l2aux/vMUS01.60
# ~/mlspgs/util/fixAttribute_batch.sh -v \
#   -e /software/toolkit/hdftools/h5edit \
#   -f /users/pwagner/mlspgs/tests/lib/NAG.Linux-6.2-cool/test \
#   
# --------------- End fixAttribute_batch.sh help
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
my_name=fixAttribute_batch.sh
I=fixAttribute_batch
dryrun="no"
verbose="no"

old_dir=`pwd`
years=`/bin/ls `
#days=001
for year in $years
do
  ~/mlspgs/util/fixAttribute.sh $@ $old_dir/$year
done
exit
