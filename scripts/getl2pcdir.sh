#! /bin/sh
# This shell script prints the name of the l2pc directory, which is either
# /scratch/livesey or if /scratch is not available, is /bigdata/livesey/auraintercomparison
if [ -d /scratch/livesey ]; then
    echo /scratch/livesey
else
    echo /bigdata/livesey/v1.0.1
fi
