#!/bin/sh
# postl2script.sh
# script to be run after (a) and (b) master tasks
# by mlsnrt-dual-l2.sh
# its name is passed by the POSTL2SCRIPT environment variables
#
# Copyright 2017, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# Now if the tool h5repack in LEVEL1_BINARY_DIR
# and if the current working directory houses the product files,
# then as a final step repack them

# usage: postl2script.sh c_dir a_dir b_dir
#   arg         meaning
#   ---         -------
#  c_dir        where to put combined std prods
#  a_dir        where to find a task's std prods
#  b_dir        where to find b task's std prods

# Notes:
# It is assumed that l2gpcat and l2auxcat are executable and in your path
# Possible alternatives are
# (1) Define MLSTOOLS and put them there
# (2) Define env variables L2GPCAT and l2AUXCAT with their locations
# directory as this script
# 

#---------------------------- commando
# May
# * echo a command with all its args;
# * execute a command with all its args, or 
# * merely echo the command that would be executed
# Which to do depends on $dryrun and $verbose

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
      
#------------------------------- executable_or_exit ------------
#
# Check that a named file is executable;
# if not, then print an error message and quit

executable_or_exit()
{
  a=`which $2`
  if [ ! -x "$a"  ]
  then
    echo "Sorry, $1 not executable or not in your PATH"
    exit 1
  fi
}

#------------------------------- h2o_nn_retrieval ------------
#
# Perform a separate retrieval for H2O
# using a trained n-n algorithm
# Then overwrite the swaths in the std. prod. "file" and in the "DGG"
# Args
# L1BRADG L1BRADD L1BOA Weights_File Prediction_File 
#
# In effect we are running
#   python3 $H2ONNSCRIPT $L1BRADG $L1BRADD $L1BOA $Weights_File $Prediction_File
# Afterwards we 
# * use insertl2gpvalues to overwrite the relevant levels in the h2o swath
# * use l2gpcat to copy the changed swath to the DGG from the std. prod. file
# *** Scratch that last point ***
# We could have used a tool like l2gpcp instead of l2gpcat. But, since we
# have insertl2gpvalues handy, why not use insertl2gpvalues on the DGG file, too?
# ***                         ***

h2o_nn_retrieval()
{
  pwd
  ls
  swath="H2O"
  echo $PYTHON3 $H2ONNSCRIPT $@
  $PYTHON3 $H2ONNSCRIPT $@
  file=`echo *L2GP-$swath*.he5`
  DGG=`echo *L2GP-DGG_*.he5`
  echo $insertl2gpvalues $file
  $insertl2gpvalues -s $swath -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    $file
  # $l2gpcat -nodup -s $swath -r $swath-StdProd -o $DGG $file 
  # l2gpcat doesn't act like we thought it did, so we'll just 
  # reuse insertl2gpvalues to handle the DGG file
  $insertl2gpvalues -s $swath-StdProd -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    $DGG
}

#------------------------------- swath_nn_retrieval ------------
#
# Perform a separate retrieval for a general swath named "$swath"
# using a trained n-n algorithm
# Then overwrite the swaths in the std. prod. "file" and in the "DGG"
# Args
# swath SWATHNNSCRIPT L1BRADG L1BRADD L1BOA Weights_File Prediction_File 
#
# In effect we are running
#   python3 $SWATHNNSCRIPT $L1BRADG $L1BRADD $L1BOA $Weights_File $Prediction_File
# Afterwards we 
# * use insertl2gpvalues to overwrite the relevant levels in the swath
# ***                         ***

swath_nn_retrieval()
{
  pwd
  ls
  swath="$1"
  SWATHNNSCRIPT="$2"
  shift
  shift
  echo $PYTHON3 $SWATHNNSCRIPT $@
  $PYTHON3 $SWATHNNSCRIPT $@
  file=`echo *L2GP-$swath*.he5`
  DGG=`echo *L2GP-DGG_*.he5`
  echo $insertl2gpvalues $file
  $insertl2gpvalues -s $swath -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    $file
  $insertl2gpvalues -s $swath-StdProd -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    $DGG
}

#------------------------------- merge_files ----------------------------
# Mrge the DGG and DGM files produced in the a and b subdirectories
# by the two independent nrt mlsl2 runs
# ***                         ***

merge_files()
{
for i in $list
do
  echo "cwd/a_dir $cwd/$a_dir"
  #ls $cwd/$a_dir
  #ls $cwd/$a_dir | grep "L2GP-${i}_" | grep -e 'he5$'
  if [ "$#" -lt "4" ]
  then
    echo No PCF supplied, so we will look in $a_dir
    name1=`ls $a_dir | grep "L2GP-${i}_" | grep -e 'he5$'`
    if [ ! -f "$a_dir/$name1" ]
    then
      echo "Sorry: unable to find $i as $name1"
      exit 1
    fi
  else
    echo Hope to find name of this file in PCF $4
    grep "L2GP-${i}_" $4 | awk -F'|' '{print $2}'
    name1=`grep "L2GP-${i}_" $4 | awk -F'|' '{print $2}'`
  fi
  echo "name1 $name1"
  a_file=$a_dir/$name1
  b_file=$b_dir/$name1
  c_file=$c_dir/$name1
  echo "a,b,c_file $a_file $b_file $c_file"
  # Now there are 4 possibilities as to whether a_file and b_file exist
  if [ -f "$a_file" -a -f "$b_file" ]
  then
    # (1) Both exist :: prefer b over a
    echo "$i: both exist"
    echo $l2gpcat -v -nodup -o $c_file $b_file $a_file
    $l2gpcat -nodup -v -o $c_file $b_file $a_file
    cp $b_file.xml $c_dir
  elif [ -f "$a_file" ]
  then
    # (2) Only a_file exists
    echo "$i: a exists"
    cp $a_file* $c_dir
  elif [ -f "$b_file" ]
  then
    # (3) Only b_file exists
    echo "$i: b exists"
    cp $b_file* $c_dir
  else
    # (4) Neither exists
    echo "Sorry unable to find $i in $a_dir $b_dir"
  fi
done
# Now the DGM files
name1=`ls $a_dir | grep "L2AUX-DGM_" | grep -e 'h5$'`
echo "name1 $name1"
a_file=$a_dir/$name1
b_file=$b_dir/$name1
c_file=$c_dir/$name1
echo "a,b,c_file $a_file $b_file $c_file"
a_file=$a_dir/$name1
echo $l2auxcat -nodup -v -o $c_file -g $b_file $b_file $a_file
$l2auxcat -nodup -v -o $c_file -g $b_file $b_file $a_file
cp $b_file.xml $c_dir
}
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
#
verbose="yes"
I=postl2script
# Is there an env file we are supposed to source?
ENVFILE=job.env
if [ -f "$ENVFILE" ]
then
 . ./"$ENVFILE"
fi

# Now check on possible candidates for l2gpcat, l2auxcat, etc.
l2gpcat=`which l2gpcat 2>/dev/null`
if [ ! -x "$l2gpcat" ]
then
  l2gpcat=$MLSTOOLS/l2gpcat
fi
l2auxcat=`which l2auxcat 2>/dev/null`
if [ ! -x "$l2auxcat" ]
then
  l2auxcat=$MLSTOOLS/l2auxcat
fi
if [ "$L2GPCAT" != "" ]
then
  l2gpcat=$L2GPCAT
fi
if [ "$L2AUXCAT" != "" ]
then
  l2auxcat=$L2AUXCAT
fi

insertl2gpvalues=`which insertl2gpvalues 2>/dev/null`
if [ ! -x "$insertl2gpvalues" ]
then
  insertl2gpvalues=$MLSTOOLS/insertl2gpvalues
fi
if [ "$INSERTL2GPVALUES" != "" ]
then
  insertl2gpvalues=$INSERTL2GPVALUES
fi

executable_or_exit l2gpcat "$l2gpcat"
executable_or_exit l2auxcat "$l2auxcat"

list='CO H2O HNO3 N2O O3 SO2 Temperature DGG'
echo "Launching postl2script with args $@"
c_dir=$1
a_dir=$2
b_dir=$3
pwd
cwd=`pwd`

# Merge the product files from a_dir and b_dir
# but only if they exist
if [ -d "$a_dir" -a -d "$b_dir" ]
then
  merge_files $@
fi

if [ "$H2ONNSCRIPT" = '' ]
then
  H2ONNSCRIPT=$MLSTOOLS/h2o_prediction.py
fi

if [ "$O3NNSCRIPT" = '' ]
then
  O3NNSCRIPT=$MLSTOOLS/o3_prediction.py
fi

if [ "$CONNSCRIPT" = '' ]
then
  CONNSCRIPT=$MLSTOOLS/co_prediction.py
fi

if [ "$TEMPNNSCRIPT" = '' ]
then
  TEMPNNSCRIPT=$MLSTOOLS/temp_prediction.py
fi

echo "PYTHON3 $PYTHON3"
echo "H2ONNSCRIPT $H2ONNSCRIPT"

# -------------------------- Uncomment the following ------------
if [ ! -f "$PYTHON3" ]
then
  echo "PYTHON3 not defined so skipping n-n retrievals"
elif [ ! -f "$insertl2gpvalues" ]
then
  echo "insertl2gpvalues not found so skipping n-n retrievals"
else
  cd $c_dir
  radg=*RADG*.h5
  radd=*RADD*.h5
  boa=*L1BOA*.h5
  if [ -f "$H2ONNSCRIPT" ]
  then
    pred=h2o_prediction.h5
    commando swath_nn_retrieval H2O $H2ONNSCRIPT $radg $radd $boa $H2OANNWEIGHTS $pred
  fi
  if [ -f "$O3NNSCRIPT" ]
  then
    pred=o3_prediction.h5
    commando swath_nn_retrieval O3 $O3NNSCRIPT $radg $radd $boa $O3ANNWEIGHTS $pred
  fi
  if [ -f "$CONNSCRIPT" ]
  then
    pred=co_prediction.h5
    commando swath_nn_retrieval CO $CONNSCRIPT $radg $radd $boa $COANNWEIGHTS $pred
  fi
  if [ -f "$TEMPNNSCRIPT" ]
  then
    pred=temp_prediction.h5
    commando swath_nn_retrieval Temperature $TEMPNNSCRIPT $radg $radd $boa $TEMPANNWEIGHTS $pred
  fi
  cd $cwd
fi
# -------------------------- -------------------------- ------------
# $Log$
# Revision 1.6  2022/05/12 22:36:46  pwagner
# Corrected buggy use of l2gpcat
#
# Revision 1.5  2022/05/11 23:41:44  pwagner
# Fixed bug in args to insertl2gpvalues
#
# Revision 1.4  2022/05/05 21:48:58  pwagner
# Added internal function h2o_nn_retrieval
#
# Revision 1.3  2017/10/12 00:00:49  pwagner
# Source job.env if it exists
#
# Revision 1.2  2017/07/13 17:42:31  pwagner
# Fixed error when species in b but not in a
#
# Revision 1.1  2017/05/17 22:26:44  pwagner
# First commit
#
