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
#     the same directory as this script
# (3) The most recent nrt versions rely entirely on neural networks and
#     don't require the simultaneous level 2 tasks runing in a_dir, b_dir;
#     eventually we should strip out the now-obsolete option, leaving
#     this script leaner and cleaner
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
      
#---------------------------- print_nth_array_element
#
# Function to print the nth element in a space-delimited list
# where n is the first arg, the list is the second arg, 
# print_nth_list_element 3 'a b c d'
# writes 'c' to standard output

print_nth_array_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: print_nth_array_element n 'a b c ..'"
         exit 1
      fi
      
      perl -e '@parts=split(" ","$ARGV[0]"); print $parts[$ARGV[1]-1]' "$2" "$1"
}

#---------------------------- what_array_element
#
# Function to show what element in a space-delimited list has been supplied
# in the first arg, the list is the second arg, 
# what_array_element c 'a b c d'
# writes '3' to standard output
what_array_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: what_array_element c 'a b c ..'"
         exit 1
      fi
      n=0
      count_up=1
      for i in $2
      do
        # echo $i $1
        if [ "$i" = "$1" ]
        then
          # echo "really, $i $1"
          n=$count_up
        fi
        count_up=`expr $count_up + 1`
      done
      echo $n
}

#---------------------------- print_hash_element
#
# Function to print value corresponding to key given two space-delimited arrays
# where the first array are the keys and the second the values
# print_hash_element c 'a b c d' '1 2 3 4'
# writes '3' to standard output
print_hash_element()
{
   # Do we have enough args?
      if [ $# -lt 3 ]
      then
         echo "Usage: print_hash_element c 'a b c ..' '1 2 3 4 ..'"
         exit 1
      fi
      arg_num=`what_array_element $1 "$2"`
      print_nth_array_element $arg_num "$3"
}

#------------------------------- mega_buck ------------
#
# Function to force multiple-evaluation of its second arg
# i.e., $$..$color (where the number of $ signs is n)
# usage: mega_buck n color

mega_buck()
{
   # Trivial case (n is 0)
   if [ "$1" -lt 1 ]
   then
     mega_buck_result="$2"
   elif [ "$2" = "" ]
   then
     mega_buck_result="$2"
   else
     mega_buck_result=`pr_mega_buck $1 $2`
   fi
}
      
pr_double_buck()
{
      eval echo $`echo $1`
}
      
pr_mega_buck()
{
      number=0
      mega_buck_temp=$2
      while [ "$number" -lt "$1" ]
      do
        mega_buck_temp=`pr_double_buck $mega_buck_temp`
        number=`expr $number + 1`
      done
      echo $mega_buck_temp
}

      
#---------------------------- lhs_eq_rhs
#
# Function to assign second arg to first: lhs=rhs
# where (lhs rhs) = ($1 $2)
# if given optional third arg "n"
# assigns $$..$lhs = $$..$rhs

lhs_eq_rhs()
{

      if [ $# -lt 3 ]
      then
         eval "$1"="'"$2"'"
      else
         mega_buck $3 $1
         lhs=$mega_buck_result
         mega_buck $3 $2
         rhs=$mega_buck_result
         eval "$lhs"="'"$rhs"'"
      fi
}

#---------------------------- set_if_not_def
#
# Function to assign second arg to first only if first not already defined
set_if_not_def()
{
     # echo $@
     mega_buck 1 $1
     # echo $mega_buck_result
     if [ "$mega_buck_result" = "" -o "$mega_buck_result" = '$' ]
     then
       lhs_eq_rhs $@
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
#  pwd
#  ls
  swath="$1"
  SWATHNNSCRIPT="$2"
  shift
  shift
  echo $PYTHON3 $SWATHNNSCRIPT $@
  $PYTHON3 $SWATHNNSCRIPT $@
  file=`echo *L2GP-$swath*.he5`
  DGG=`echo *L2GP-DGG_*.he5`
  # echo $insertl2gpvalues $file
  commando $insertl2gpvalues -s $swath -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    -conv 1 \
    -qual 0 \
    -stat 68 \
    $file
  commando $insertl2gpvalues -s $swath-StdProd -d ANN_Prediction -p ANN_Precision \
    -Lf $3 \
    -Vf $5 \
    -conv 1 \
    -qual 0 \
    -stat 68 \
    $DGG
  # If createattributes is defined and is an executable, use it
  # to add job.env as a file of attribute name=value pairs
  echo "createattributes $createattributes"
  echo "ATTRFILE $ATTRFILE"
  if [ -x "$createattributes" -a -f "$ATTRFILE" ]
  then
    commando $createattributes -Vf "$ATTRFILE" -l2gp $file
    commando $createattributes -Vf "$ATTRFILE" -l2gp $DGG
  fi
}

#------------------------------- merge_files ----------------------------
#------------------------------- (Obsolete) ----------------------------
# Merge the DGG and DGM files produced in the a and b subdirectories
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
#pwd
#ls
# Are there any env files we are supposed to source?
ENVFILE=job.env
if [ -f "$ENVFILE" ]
then
  . ./"$ENVFILE"
fi

# Did the ENVFILE define an ATTRFILE?
if [ -f "$ATTRFILE" ]
then
  echo "dotting $ATTRFILE"
  . "$ATTRFILE"
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

createattributes=`which createattributes 2>/dev/null`
if [ ! -x "$createattributes" ]
then
  createattributes=$MLSTOOLS/createattributes
fi
if [ "$CREATEATTRIBUTES" != "" ]
then
  createattributes=$CREATEATTRIBUTES
fi

executable_or_exit l2gpcat "$l2gpcat"
executable_or_exit l2auxcat "$l2auxcat"

list='CO H2O HNO3 N2O O3 SO2 Temperature DGG'
echo "Launching postl2script with args $@"
c_dir=$1
a_dir=$2
b_dir=$3
#pwd
cwd=`pwd`

# Merge the product files from a_dir and b_dir
# but only if they exist
if [ -d "$a_dir" -a -d "$b_dir" ]
then
  merge_files $@
fi

set_if_not_def H2ONNSCRIPT $MLSTOOLS/h2o_prediction.py
set_if_not_def O3NNSCRIPT $MLSTOOLS/o3_prediction.py
set_if_not_def CONNSCRIPT $MLSTOOLS/co_prediction.py
set_if_not_def SO2NNSCRIPT $MLSTOOLS/so2_prediction.py
set_if_not_def TEMPNNSCRIPT $MLSTOOLS/temp_prediction.py
set_if_not_def N2ONNSCRIPT $MLSTOOLS/n2o_prediction.py
set_if_not_def HNO3NNSCRIPT $MLSTOOLS/hno3_prediction.py
set_if_not_def CLOUDNNSCRIPT $MLSTOOLS/cloud_prediction.py

echo "PYTHON3 $PYTHON3"
echo "H2ONNSCRIPT $H2ONNSCRIPT"

if [ ! -f "$PYTHON3" ]
then
  echo "PYTHON3 not defined so skipping n-n retrievals"
elif [ ! -f "$insertl2gpvalues" ]
then
  echo "insertl2gpvalues not found so skipping n-n retrievals"
else
  cd $c_dir
  radg=`ls *RADG*.h5`
  radd=`ls *RADD*.h5`
  boa=`ls *L1BOA*.h5`
  # It's possible we'll have a PCF defined but the above files don't exist
  if [ ! -f "$boa" ]
  then
    line=`sed -n '/^21050/p' ${PGS_PC_INFO_FILE}`
    a=`echo $line | awk -F'|' '{print $2}'`
    b=`echo $line | awk -F'|' '{print $3}'`
    radd=$b/$a
    line=`sed -n '/^21051/p' ${PGS_PC_INFO_FILE}`
    a=`echo $line | awk -F'|' '{print $2}'`
    b=`echo $line | awk -F'|' '{print $3}'`
    radg=$b/$a
    line=`sed -n '/^21110/p' ${PGS_PC_INFO_FILE}`
    a=`echo $line | awk -F'|' '{print $2}'`
    b=`echo $line | awk -F'|' '{print $3}'`
    boa=$b/$a
  fi
  echo "boa $boa"
  echo "radg $radg"
  echo "radd $radd"
  # The following arrays hold corresponding names for
  # molecules  swaths  pythonscripts  weightsfiles
  mols="h2o o3 co so2 temp n2o hno3"
  swaths="H2O O3 CO SO2 Temperature N2O HNO3"
  scripts="$H2ONNSCRIPT $O3NNSCRIPT $CONNSCRIPT $SO2NNSCRIPT $TEMPNNSCRIPT $N2ONNSCRIPT $HNO3NNSCRIPT"
  weights="$H2OANNWEIGHTS $O3ANNWEIGHTS $COANNWEIGHTS $SO2ANNWEIGHTS $TEMPANNWEIGHTS $N2OANNWEIGHTS $HNO3ANNWEIGHTS"
  # We may be called to predict only cloud tops, instead of the standard molecules
  if [ ! -f "$H2OANNWEIGHTS" ]
  then
    mols="cloud"
    swaths="CloudTopPressure"
    scripts="$CLOUDNNSCRIPT"
    weights="$CLOUDANNWEIGHTS"
  fi
  echo "weights: $weights"
  for mol in $mols
  do
    # mol_num is the index number in mols corresponding to mol
    # i.e., 1 2 .. n
    mol_num=`what_array_element $mol "$mols"`
    swath=`print_nth_array_element $mol_num "$swaths"`
    script=`print_nth_array_element $mol_num "$scripts"`
    weight=`print_nth_array_element $mol_num "$weights"`
    pred=${mol}_prediction.h5
    commando swath_nn_retrieval $swath $script $radg $radd $boa $weight $pred
  done
  cd $cwd
fi
# -------------------------- -------------------------- ------------
# $Log$
# Revision 1.12  2024/02/02 21:44:56  pwagner
# Cope with being used for CloudTopPressure by v6 level 2
#
# Revision 1.11  2023/09/01 18:38:43  pwagner
# Added feature to create attributes in product files
#
# Revision 1.10  2022/10/13 22:19:48  pwagner
# Added hno3 and n2o swaths; using arrays like mols instead of multiple ifs
#
# Revision 1.9  2022/07/13 22:53:49  pwagner
# Use new insertl2gpvalues cmdline opts
#
# Revision 1.8  2022/06/02 21:39:35  pwagner
# Added SO2 to mols predicted by python script
#
# Revision 1.7  2022/05/27 20:59:32  pwagner
# Can now use python scripts to predict CO, O3, and Temperature; can skip merger of a, b
#
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
