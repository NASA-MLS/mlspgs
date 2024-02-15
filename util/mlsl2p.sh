#!/bin/sh
# mlsl2p.sh
# runs the master task for mlsl2
#
# Assumes that:

# (1) Either
#  ----
# (a) It has been called from PGS_PC_Shell.sh as the (pge) in
#  PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v; or
# (b) It has been called to work directly w/o the toolkit panoply
#       with a single command-line arg: the l2cf file
# (a) and (b) are distinguished by the env variable UsingPCF
# If it's defined, and "yes" (or in fact just non-blank)
# then we'll assume we have case(a). Otherwise, we have case (b)
#  ----
# (2) $PGE_BINARY_DIR contains mlsl2 and $PGE_SCRIPT_DIR contains
#       slavetmplt.sh and jobstat-sips.sh
#     (or else define an enviromental variable (PVM_EP) and put them there)
# (3) PVM_HOSTS_INFO is defined as an environment variable
#     It should be the path and name of the host file,
#     a text file containing the hosts available
#     for running the slave tasks, one host per line
# (4) JOBDIR is defined as an environment variable
#     It should be the path where the job is run
# (5) PGE_ROOT is defined as an environment variable
#     It should be the path where the science_env.sh script is kept
# (6) OTHEROPTS is defined as an environment variable
#     It would contain other meaningful runtime options, 
#       e.g. OTHEROPTS="--skipRetrieval"
# (7) POSTL2SCRIPT (just like mlsnrt) is defined
# (8) ATTRFILE (just like mlsnrt) is defined
#
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# At the end of the master task, if the tools
#  h5repack 
#  aug_hdfeos5
# are in PGE_BINARY_DIR
# and if either
#   the current working directory ($cwd) or
#   $cwd/outputs
# houses the product files,
# then as a final step repack and augment them with netCDF attributes.

# As an alternative, if NETCDFCONVERT is defined and executable,
# use that instead to create an explicitly NetCDF-formatted file.
# Additionally, if NCCOPY is defined and executable,
# use that tool to compress the NetCDF file.

# usage: see (1) above

#------------------------------- extant_files ------------
#
# Function to return only those files among the args
# that actually exist
# Useful when passed something like *.f which may 
# (1) expand to list of files, returned as extant_files_result, or
# (2) stay *.f, in which case a blank is returned as extant_files_result 
#     (unless you have perversely named a file '*.f')
# usage: extant_files arg1 [arg2] ..

extant_files()
{
   extant_files_result=
   # Trivial case ($# = 0)
   if [ "$1" != "" ]
   then
      for file
      do
         if [ -f "$file" ]
         then
               extant_files_result="$extant_files_result $file"
         fi
      done
   fi
   echo $extant_files_result
}

#------------------------------- hide_files ------------
#
# Hide files when something goes awry
# usage: hide_files arg1 [arg2] ..

hide_files()
{
   hide_files_result=
   # Trivial case ($# = 0)
   if [ "$1" != "" ]
   then
      for file
      do
         if [ -f "$file" ]
         then
               hide_files_result="$hide_files_result $file"
         fi
      done
   fi
   echo $hide_files_result
   if [ ! -d hidden ]
   then
     mkdir hidden
   fi
   mv $hide_files_result hidden
}

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
       echo "$@" >> "$MLSL2PLOG"
       ;;
     both )
       echo "$@" | tee -a "$MLSL2PLOG"
       ;;
     esac
   fi
}
      
#-------------- mycp ----------------
# read each line of stdin
# and then echo it to stdout
# Why?
# Lines like
#   b=$a
# will be echoed as
#   b=100
# if the env variable a has been set to 100
mycp()
{
  set -o noglob
  while read line; do
    # Does line begin with a '#'?
    acomment=`echo $line | grep '^#'`
    if [ "$acomment" != "" ]
    then
      # Dont do anything special--just a comment
      echo $line
    else
      # May need to evaluate twice if line contains a $variable
      eval echo $line
    fi
  done
  set +o noglob
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

#---------------------------- rename_all_nc
#
# Function to rename all the .nc-marked l2gp files to .he5
rename_all_nc()
{
  files=`extant_files *L2GP-[A-CE-Z]*.nc *L2GP-DGG_*.nc`
  for file in $files
  do
    hname=`echo $file | sed 's/\.nc/\.he5/'`
    mv $file $hname
  done
}

#---------------------------- create_nc_metafiles
#
# Function to create metadata files for
# all the .nc-marked l2gp files
create_nc_metafiles()
{
  files=`extant_files *L2GP-[A-CE-Z]*.nc *L2GP-DGG_*.nc`
  for file in $files
  do
    hname=`echo $file | sed 's/\.nc/\.he5/'`
    # 1st: the .met files
    if [ -f "$hname.met" ]
    then
      cp $hname.met $file.met.tmp
      sed 's/\.he5/\.nc/' $file.met.tmp > $file.met
      rm $file.met.tmp
    fi
    # 2nd: the .xml files
    if [ -f "$hname.xml" ]
    then
      cp $hname.xml $file.xml.tmp
      sed 's/\.he5/\.nc/' $file.xml.tmp > $file.xml
      rm $file.xml.tmp
    fi
  done
}

# ------------------ check for misalignment ----------------------
# ----------------------------------------------------------------
# Check products for misaligned geolocations
# If they are found to be misaligned, set return_status to 99
check_for_misalignment()
{
if [ -x "$MISALIGNMENT" -a "$file" != "" ]
then
  logit "Checking for misalignment"
  a=`$MISALIGNMENT -silent $file`
  if [ "$a" != "" ]
  then
    logit $a
    logit "Misalignment detected"
    logit hide_files *.he5 *.met *.xml *.h5
    hide_files *.he5 *.met *.xml *.h5
    return_status=99
  fi
fi

}

# ------------------ convert the product file format ----------------------
# ----------------------------------------------------------------
# Convert hdfeos product files to netcdf
# A trick: the hdfeos file names already end in '.nc' so
# we'll need to mv them so they end in .he5
# When we're done we'll have 2 versions of each std. prod. file:
#   product.he5      hdfeos
#   product.nc       NetCDF4
# 
# Possible improvements:
# Option to end file names with ".nc4" instead of ".nc"
# Option to hide the hdfeos versions so they aren't copied
# to the scf or ingested into the database
convert_file_formats()
{
logit "NETCDFCONVERT: $NETCDFCONVERT"
logit "NCCOPY: $NCCOPY"
if [ -x "$NETCDFCONVERT" ]
then
  logit "Converting file format to NetCDF4"
  files=`extant_files *L2GP-*.he5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2GP-*.he5`
    fi
  fi
  for file in $files
  do
    ncname=`echo $file | sed 's/\.he5/\.nc/'`
    logit "Converting $file to $ncname"
    $NETCDFCONVERT $file
    # The tool insists on naming its output file "something.nc4"
    # So we must rename it yet again to "something.nc"
    # What a shameful mess we've made
    nc4name=`echo $file | sed 's/\.he5/\.nc4/'`
    logit "mving $nc4name to $ncname"
    mv $nc4name $ncname
  done
fi
}

# ------------------ nccopy_files ----------------------
# ----------------------------------------------------------------
# h5repack doesn't work properly with NetCDF files
# Instead we must use nccopy
nccopy_files()
{
if [ -x "$NCCOPY" ]
then
  logit "Doing nccopy on files"
  files=`extant_files *.nc`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *.nc`
    fi
  fi
  for file in $files
  do
    if [ -w "$file" ]
    then
      packed="$file".p
      if [ "$GZIPLEVEL" != "" ] 
      then
        filter="-d $GZIPLEVEL"
      else
        filter=""
      fi
      logit "nccopying $file into $packed"
      logit $NCCOPY  $filter "$file" "$packed"
      $NCCOPY  $filter "$file" "$packed"
      # Here we could insert some check involving h5diff if we were dubious
      mv "$packed" "$file"
    fi
  done
fi
# ----------------------------------------------------------------
}
# ----------------------------------------------------------------

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
  logit "Doing h5repack on files"
  files=`extant_files *L2FWM*.h5 *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5 *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2FWM*.h5 *L2GP-*.he5 *L2GP-*.nc *L2AUX-[A-C]*.h5 *L2AUX-DGM_*.h5`
    fi
  fi
  for file in $files
  do
    if [ -w "$file" ]
    then
      packed="$file".p
      if [ "$GZIPLEVEL" != "" ] 
      then
        filter="-f GZIP=$GZIPLEVEL"
      else
        filter=""
      fi
      logit "Packing $file into $packed"
      logit $H5REPACK -i "$file" -o "$packed" $filter
      $H5REPACK -i "$file" -o "$packed" $filter
      # Here we could insert some check involving h5diff if we were dubious
      mv "$packed" "$file"
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
  logit "Doing augment on files"
  files=`extant_files *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
  if [ "$files" = "" ]
  then
    if [ -d "outputs" ]
    then
      cd "outputs"
      files=`extant_files *L2GP-[A-CE-Z]*.he5 *L2GP-DGG_*.he5`
    fi
  fi
  logit "files: $files"
  for file in $files
  do
    if [ -w "$file" ]
    then
      logit $NETCDFAUGMENT $file
      $NETCDFAUGMENT $file
    fi
  done
fi
# ----------------------------------------------------------------
}

# ----------------------reformat_and_repack-----------------
# ----------------------------------------------------------------
reformat_and_repack()
{

  # Will we use python scripts to predict values for any other
  # products, e.g. CloudTopHeight?
  if [ -f "$POSTL2SCRIPT" -a -f "$ATTRFILE" ]
  then
      $POSTL2SCRIPT "$OUTPUTS_C" OUTPUTS_A OUTPUTS_B PCF_A PCF_B
  fi
  convert_file_formats
  
  # Do we have metadata files for NetCDF-formatted product files yet?
  ncmetadata=`extant_files *L2GP-DGG_*.nc.met`
  if [ "$ncmetadata" = "" ]
  then
    create_nc_metafiles
  fi


  nccopy_files

  repack_files

  augment_files

}

# ----------------------wait_for_slaves_to_finish-----------------
# ----------------------------------------------------------------
wait_for_slaves_to_finish()
{
logit SAVEJOBSTATS $SAVEJOBSTATS
logit masterlog $masterlog
logit PGE_SCRIPT_DIR/jobstat-sips.sh $PGE_SCRIPT_DIR/jobstat-sips.sh
# Save record of progress thru phases
if [ -f "$masterlog" ]
then
  l2cf=`grep -i 'Level 2 configuration file' $masterlog | head -1 | awk '{print $13}'`
fi
if [ "$SAVEJOBSTATS" = "yes" -a -x "$PGE_SCRIPT_DIR/jobstat-sips.sh" -a -f "$l2cf" ]
then
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  JOBSTATSFILE="$JOBDIR/phases.stats"
  JOBSOPTS="-S $PGE_BINARY_DIR/mlsqlog-scan-sips.py -t $PGE_SCRIPT_DIR/split_path.sh ${JOBDIR}/pvmlog $l2cf $masterlog"
  if [ "$verbose" != "" ]
  then
    JOBSOPTS="-v $JOBSOPTS"
  fi
  logit "$PGE_SCRIPT_DIR/jobstat-sips.sh $JOBSOPTS"
  $PGE_SCRIPT_DIR/jobstat-sips.sh $JOBSOPTS > "$JOBSTATSFILE"
else
  JOBSTATSFILE="none"
fi

# catenate each slave's log to a log file
LOGFILE="${JOBDIR}/pvmlog/mlsl2.log"
if [ ! -w "$LOGFILE"  ]
then
  logit "catenating slave logs in $LOGFILE"
  echo "catenating slave logs in $LOGFILE" > "$LOGFILE".1
  # This sleep is to give slave tasks extra time to complete stdout
  sleep 20
  cat ${JOBDIR}/pvmlog/*.log >> "$LOGFILE".1
  rm -f ${JOBDIR}/pvmlog/*.log
  mv "$LOGFILE".1 "$LOGFILE"
fi

if [ -f "$LOGFILE" -a -f "$JOBSTATSFILE" ]
then
  cat "$LOGFILE" "$JOBSTATSFILE" > "$LOGFILE".1
  mv "$LOGFILE".1 "$LOGFILE"
fi
}

#------------------------------- run_mlsl2_ntk ------------
#
# Function to do one level 2 run after expanding an l2cf template
# so: no toolkit is required
# usage:
# run_mlsl2_ntk [options] l2cf
# args:
# options        recognized cmdline args for mlsl2
# l2cf           level 2 config file
#

run_mlsl2_ntk()
{
  logit "Running mlsl2 master w/o toolkit"
  logit "args $@"
  # First a pre-flight run to check paths                               
  # If a problem, disclosed by exit status, exit before starting big run
  if [ "$CHECKPATHS" = "yes" ]                                          
  then                                                                  
    $PGE_BINARY --checkPaths --ntk $otheropts $l2cf
    return_status=`expr $?`                                             
    if [ $return_status != $NORMAL_STATUS ]                             
    then                                                                
                                                                        
       logit "Preflight checkPaths run ended badly"                      
       logit "Possibly an error in pathnames; please check your l2cf"    
       exit 1                                                           
    fi                                                                  
  fi                                                                    
                                                                        
  # Next, create a fresh ident based on master task's binary executable 
  if [ "$CHECKIDENTS" = "yes" ]                                         
  then                                                                  
    IDENTFILE="$JOBDIR/master.ident"                                    
    rm -f "$IDENTFILE"                                                  
    ident $PGE_BINARY > "$IDENTFILE"                                    
  else                                                                  
    IDENTFILE="none"                                                    
  fi         

  # HOSTNAME will be stored in metadata as ProductionLocation
  # to trace provenance
  if [ "$HOSTNAME" = "" ]
  then
    HOSTNAME=local
  fi

  # Now we launch the master task itself to set everything in motion  
  $PGE_BINARY --pge $slave_script --ntk --master $PVM_HOSTS_INFO \
    --idents "$IDENTFILE" --loc "$HOSTNAME" \
    $otheropts $l2cf 2> $JOBDIR/master.l2cfname
  return_status=`expr $?`
  wait_for_slaves_to_finish
}

#------------------------------- run_mlsl2_tk ------------
#
# Function to do one level 2 run after expanding a pcf template
# so this run requires the toolkit panoply
# usage:
# run_mlsl2_tk [options]
# args:
# options        recognized cmdline args for mlsl2
#

run_mlsl2_tk()
{
  # First a pre-flight run to check paths
  # If a problem, disclosed by exit status, exit before starting big run
  if [ "$CHECKPATHS" = "yes" ]
  then
    $PGE_BINARY --checkPaths --tk $otheropts 2> $JOBDIR/master.l2cfname
    return_status=`expr $?`
    if [ $return_status != $NORMAL_STATUS ]
    then
       cat $JOBDIR/master.l2cfname
       logit "Preflight checkPaths run ended badly"
       logit "Possibly an error in pathnames; please check your PCF"
       cat $JOBDIR/master.l2cfname >> "$MLSL2PLOG"
       exit 1
    fi
  fi

  # Next, create a fresh ident based on master task's binary executable
  if [ "$CHECKIDENTS" = "yes" ]
  then
    IDENTFILE="$JOBDIR/master.ident"
    rm -f "$IDENTFILE"
    ident $PGE_BINARY > "$IDENTFILE"
  else
    IDENTFILE="none"
  fi

  # HOSTNAME will be stored in metadata as ProductionLocation
  # to trace provenance
  if [ "$HOSTNAME" = "" ]
  then
    logit "Need HOSTNAME to trace provenance as ProductionLocation"
    exit 1
  fi

    if [ "$CAPTURE_MT" = "yes" -a "$STDERRFILE" = "" ]
    then
      STDERRFILE=mlsl2p.stderr
    fi

  # Now we launch the master task itself to set everything in motion
  if [ "$CAPTURE_MT" = "yes" ]
  then
    /usr/bin/time -f 'M: %M t: %e' $PGE_BINARY --pge $slave_script --tk --master $PVM_HOSTS_INFO \
      --idents "$IDENTFILE" --loc "$HOSTNAME" $otheropts 2> "$STDERRFILE"
  elif [ "$STDERRFILE" != "" ]
  then
    $PGE_BINARY --pge $slave_script --tk --master $PVM_HOSTS_INFO \
      --idents "$IDENTFILE" --loc "$HOSTNAME" $otheropts 2> "$STDERRFILE"
  else
    $PGE_BINARY --pge $slave_script --tk --master $PVM_HOSTS_INFO \
      --idents "$IDENTFILE" --loc "$HOSTNAME" $otheropts
  fi
  # Save return status
  return_status=`expr $?`
  wait_for_slaves_to_finish
}

# ------------------------------- run_serially
run_serially()
{
  if [ "$PGE_BINARY" = "" ]
  then
    PGE_BINARY=$BIN_DIR/$MLSPROG
  fi
  echo $PGE_BINARY $EXTRA_OPTIONS "$@"
  
  if [ ! -x "$PGE_BINARY"  ]
  then
    echo "************************************************************************"
    echo "$PGE_BINARY doesn't exist!"
    echo "Check for typos in its path and program names"
    echo "************************************************************************"
    exit 1
  fi
  if [ -f $BIN_DIR/license.txt ]
  then
    cat $BIN_DIR/license.txt
  fi
  echo $PGE_BINARY $otheropts "$@"
  if [ "$CAPTURE_MT" = "yes" -a "$STDERRFILE" = "" ]
  then
    STDERRFILE=$MLSPROG.stderr
  fi
  if [ "$UsingPCF" != "" ]
  then
    otheropts="--tk $otheropts"
  fi
  if [ "$CAPTURE_MT" = "yes" ]
  then
    /usr/bin/time -f 'M: %M t: %e' $PGE_BINARY $otheropts "$@" 2> "$STDERRFILE"
  elif [ "$STDERRFILE" != "" ]
  then
    $PGE_BINARY $otheropts "$@" 2> "$STDERRFILE"
  else
    $PGE_BINARY $otheropts "$@"
  fi
  return_status=`expr $?`
  H5REPACK=$BIN_DIR/h5repack
  MISALIGNMENT=$BIN_DIR/misalignment
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
CHECKPATHS="yes"
#           ^^^---- "yes" if extra preflight non-retrieval run to check paths

CHECKIDENTS="yes"
#            ^^^---- "yes" if each slave to check its ident against master

GZIPLEVEL="1"
#          ^^^---- compression level ("" means none)

MLSL2PLOG="$JOBDIR/mlsl2p.log"
#             ^^^---- "yes" if progress of chunks thru phases recorded

MUSTLOG="both"
if [ -f "$MLSL2PLOG" ]
then
  mv $MLSL2PLOG $MLSL2PLOG.1
fi
echo "Entering mlsl2p.sh with args $@"
pwd
echo "MLSL2PLOG " > "$MLSL2PLOG"

logit "$MLSL2PLOG"
env | sort >> "$MLSL2PLOG"

SAVEJOBSTATS="no"
#             ^^^---- "yes" if recording progress of chunks through each phase

# In addition to whatever options and switches may be set by the environment
# variable OTHEROPTS, the following are set:
# -g       trace path of execution through code sections
# --wall   show timing in wall clock times
# slv      insert slave outputs among master's stdout
# mas      show master's pvm activities as they occur
# chu      show chunk divisions
# opt1     show command line options
# log      copy any log file messages to stdout
# pro      announce input files at opening, output files at creation
# time     summarize time consumed by each code  section, phase, etc.
otheropts="$OTHEROPTS -g --wall -S'mas,chu,opt1,log,pro,time'"
verbose=`echo "$otheropts" | grep -i verbose`

# Check that assumptions are valid
if [ "$PGS_PC_INFO_FILE" = "" -a "$UsingPCF" != "" ]
then
  logit "************************************************************************"
  logit 'PGS_PC_INFO_FILE undefined'
  logit 'usage:'
  logit 'PGS_PC_Shell.sh (pge) 0111 (PCF_file) 25 -v'
  logit "************************************************************************"
  exit 1
elif [ "$PVM_HOSTS_INFO" = "" ]
then
  logit "************************************************************************"
  logit 'PVM_HOSTS_INFO undefined' 
  logit 'It should be the path and name of the host file'
  logit 'a text file containing the hosts available'
  logit 'for running the slave tasks, one host per line'
  logit "************************************************************************"
  exit 1
elif [ "$JOBDIR" = "" ]
then
  logit "************************************************************************"
  logit 'JOBDIR undefined' 
  logit 'It should be the path where the job is run'
  logit "************************************************************************"
  exit 1
elif [ "$PGE_ROOT" = "" ]
then
  logit "************************************************************************"
  logit 'PGE_ROOT undefined' 
  logit 'It should be the path where the science_env.sh script is kept'
  logit "************************************************************************"
  exit 1
fi

# In case you used the ep=(PGE_BINARY_DIR)
# in the host file
set_if_not_def PGE_BINARY_DIR  $HOME/pvm3/bin/LINUX
set_if_not_def PGE_BINARY      $PGE_BINARY_DIR/mlsl2

if [ ! -x "$PGE_BINARY"  ]
then
  logit "************************************************************************"
  logit "$PGE_BINARY doesn't exist!"
  logit "Check for typos in its path and program names"
  logit "************************************************************************"
  exit 1
elif [ ! -r "$PGE_SCRIPT_DIR/slavetmplt.sh"  ]
then
  logit "************************************************************************"
  logit "slavetmplt.sh not in $PGE_SCRIPT_DIR"
  logit "Check for typos in its path and program names"
  logit "************************************************************************"
  exit 1
fi

set_if_not_def H5REPACK       $PGE_BINARY_DIR/h5repack
set_if_not_def NETCDFAUGMENT  $PGE_BINARY_DIR/aug_hdfeos5
set_if_not_def NETCDFCONVERT  $PGE_BINARY_DIR/l2gp2nc4
set_if_not_def MISALIGNMENT   $PGE_BINARY_DIR/misalignment
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

# We are going to insist that both H5REPACK and NETCDFAUGMENT
# be defined before going any further. To override this, set
# the environment variable OKTONOTAUGMENT to "yes"
if [ "$OKTONOTAUGMENT" = "" ]
then
  if [ ! -x "$NETCDFAUGMENT" ]
  then
    echo "NETCDFAUGMENT not defined"
    exit 1
  elif [ ! -x "$H5REPACK" ]
  then
    echo "H5REPACK not defined"
    exit 1
  fi
fi

if [ ! -x "$NETCDFCONVERT" ]
then
  NETCDFCONVERT=$MLSTOOLS/l2gp2nc4
fi
if [ ! -x "$MISALIGNMENT" ]
then
  MISALIGNMENT=$MLSTOOLS/misalignment
fi

masterlog="${JOBDIR}/exec_log/process.stdout"
if [ "$MASTERLOG" != "" ]
then
  masterlog="$MASTERLOG"
fi

# We print the file license.txt to stdout
if [ -f $PGE_BINARY_DIR/license.txt ]
then
  cat $PGE_BINARY_DIR/license.txt
  cat $PGE_BINARY_DIR/license.txt >> "$MLSL2PLOG"
fi

# Oops, this stomps on any PCF we might have selected
# so save it to be restored
logit "Original PCF: $PGS_PC_INFO_FILE"
PCF=$PGS_PC_INFO_FILE
if [ -r "$JOBDIR/job.env"  ]
then
  . $JOBDIR/job.env
  PCF=$PGS_PC_INFO_FILE
elif [ -r "$PGE_ROOT/science_env.sh"  ]
then
  . ${PGE_ROOT}/science_env.sh
elif [ -r "$PGSBIN/pgs-env.ksh" ]
then
  . $PGSBIN/pgs-env.ksh
fi
PGS_PC_INFO_FILE=$PCF
export PGS_PC_INFO_FILE
logit "new PCF: $PGS_PC_INFO_FILE" 

# The logs will be written as separate files into ${JOBDIR}/pvmlog
# before being catenated at end of run
# Make certain this directory exists and that it is empty
if [ ! -d ${JOBDIR}/pvmlog ]
then
  mkdir ${JOBDIR}/pvmlog
else
  rm -f ${JOBDIR}/pvmlog/*.log
fi

# The following environmental variable may already have been set
if [ "$PGSMEM_USESHM" = "" ]
then
  PGSMEM_USESHM=NO
fi

export FLIB_DVT_BUFFER=0

if [ "$MASTEROPTSFILE" != "" ]
then
  mycp < $MASTEROPTSFILE > ${JOBDIR}/master.opts
  otheropts="-g --wall -S'mas,chu,opt1,log,pro,time' --setf ${JOBDIR}/master.opts"
fi
env
env | sort >> "$MLSL2PLOG"
ulimit -s unlimited
ulimit -a
NORMAL_STATUS=2

# ------------------ launch the level 2 run ----------------------
# ----------------------------------------------------------------
if [ "$SERIALRUN" = "yes" ]
then
  run_serially
  reformat_and_repack
elif [ "$UsingPCF" != "" ]
then
  # Use sed to convert slavetmplt.sh into an executable script
  # The resulting script sets some toolkit-savvy environment variables
  # and then launches the regular mlsl2 binary when summoned to do so.
  SLV_SUF=slave
  slave_script=$JOBDIR/mlsl2.$SLV_SUF
  rm -f $slave_script
  sed "s=ppccff=${PGS_PC_INFO_FILE}=;
  s=UUssiinnggPPCCFF=yes=;
  s=ppggssmmeemmuusseesshhmm=${PGSMEM_USESHM}=;
  s=ssllaavveessccrriipptt=${slave_script}=;
  s=ootthheerrooppttss=\"${otheropts}\"=;
  s=ppggeebbiinnaarryy=${PGE_BINARY}=;
  s=ppggssbbiinn=${PGSBIN}=;
  s=jjoobbddiirr=${JOBDIR}=;
  s=ppggeerroott=${PGE_ROOT}=;" $PGE_SCRIPT_DIR/slavetmplt.sh > $slave_script
  chmod a+x $slave_script

  run_mlsl2_tk
  l2cf=$L2CF
  
  # Do we have an outputs subdirectory of JOBDIR for std prods?
  if [ -d "$JOBDIR/outputs" ]
  then
    STDPRODDIR="$JOBDIR/outputs"
    OUTPUTS_C=outputs
  else
    STDPRODDIR="$JOBDIR"
    OUTPUTS_C=./
  fi
  
  cd $STDPRODDIR
  OUTPUTS_C=`pwd`
  # Did we name the dgg file *L2GP-DGG*.nc in the PCF?
  # if so, rename all the l2gp files so they end with .he5
  dggfilenc=`extant_files *L2GP-DGG_*.nc`
  if [ "$dggfilenc" != "" ]
  then
    rename_all_nc
  fi

  check_for_misalignment
  reformat_and_repack
  
else
  # Use sed to convert slavetmplt.sh into an executable script
  # The resulting script sets some toolkitless environment variables
  # and then launches the regular mlsl2 binary when summoned to do so.
  l2cf=$1
  # Check that we were called properly
  if [ ! -f "$l2cf" ]
  then
    logit "************************************************************************"
    logit 'Usage (w/o toolkit): mlsl2p.sh l2cf_file'
    logit 'Check that you really intended to run w/o toolkit'
    logit 'To run with a PCF you must set the environment variable UsingPCF=yes'
    logit "************************************************************************"
    exit 1
  fi
  SLV_SUF=slave
  slave_script=$JOBDIR/mlsl2.$SLV_SUF
  rm -f $slave_script
  sed "s=ppggssmmeemmuusseesshhmm=${PGSMEM_USESHM}=;
  s=UUssiinnggPPCCFF==;
  s=ssllaavveessccrriipptt=${slave_script}=;
  s=ootthheerrooppttss=\"${otheropts}\"=;
  s=ppggeebbiinnaarryy=${PGE_BINARY}=;
  s=ppggssbbiinn=${PGSBIN}=;
  s=jjoobbddiirr=${JOBDIR}=;
  s=ppggeerroott=${PGE_ROOT}=;" $PGE_SCRIPT_DIR/slavetmplt.sh > $slave_script
  chmod a+x $slave_script

  run_mlsl2_ntk
  check_for_misalignment
fi

# ----------------------------------------------------------------

# End toolkitless run's log file with a grep-able sign-off message
if [ "$UsingPCF" = "" ]
then
  logit " "
  logit "mlsl2p.sh: return from PGE: $return_status"
fi

if [ $return_status != $NORMAL_STATUS ]
then
   exit 1
else
   exit 0
fi

# $Log$
# Revision 1.46  2023/05/11 22:42:13  pwagner
# Now insists that H5REPACK and NETCDFAUGMENT be defined before running
#
# Revision 1.45  2022/09/28 21:41:52  pwagner
# Add comments regarding NCCOPY variable
#
# Revision 1.44  2022/05/27 21:00:41  pwagner
# Printed a more helpful error message when run w/o toolkit
#
# Revision 1.43  2022/02/03 18:45:10  pwagner
# Move augment, convert, etc. into internal functions
#
# Revision 1.42  2022/01/20 22:01:30  pwagner
# Cant SAVEJOBSTATS until we debug mlsqlog-scan-sips.py
#
# Revision 1.41  2020/04/22 16:45:40  pwagner
# May use NETCDFCONVERT instead of augmenting hdfeos files
#
# Revision 1.40  2019/07/17 20:23:04  pwagner
# Try to makee error mesgs more helpful and conspicuous
#
# Revision 1.39  2019/06/27 21:07:05  pwagner
# Print sign-off at end of toolkitless run; use logit mechanism
#
# Revision 1.38  2019/04/18 16:22:08  pwagner
# May evaluate variables in opts file if USEOPTSENV is set; drop slavetmpltntk.sh because slavetmplt.sh now does double-duty
#
# Revision 1.37  2018/03/01 00:21:24  pwagner
# Now handles both --tk and --ntk cases
#
# Revision 1.36  2017/12/07 22:49:20  pwagner
# Try harder not to stomp on choice of PCF
#
# Revision 1.35  2016/11/16 19:26:50  pwagner
# Avoids stomping on an already-selected PCF
#
# Revision 1.34  2016/10/20 23:24:55  pwagner
# cat master.l2cfname to stdout and to master stdout
#
# Revision 1.33  2016/05/17 17:07:52  pwagner
# 'dot' job.env if found; Obey CAPTURE_MT
#
# Revision 1.32  2015/12/09 01:30:59  pwagner
# Fixed some typos
#
# Revision 1.31  2015/11/02 23:21:18  pwagner
# Now checks for mialigned geolocations
#
# Revision 1.30  2015/10/07 22:56:07  pwagner
# Automtically writes out l2cf name to master.l2cfname
#
# Revision 1.29  2014/08/14 00:51:17  pwagner
# Creates mlsl2p log file
#
# Revision 1.28  2014/04/02 23:07:10  pwagner
# Pass HOSTNAME on commandline to use in metadata
#
# Revision 1.27  2013/12/05 00:25:24  pwagner
# Fix latest problems with creating JOBSTATSFILE
#
# Revision 1.26  2013/10/29 16:58:20  pwagner
# May set environment variable PGE_BINARY
#
# Revision 1.25  2013/09/04 17:44:45  pwagner
# Replaced '--cat' cmdline option; 'Catenate' now an Output section command
#
# Revision 1.24  2012/02/15 18:12:06  pwagner
# Offer last chance to find h5repack in HDFTOOLS directory
#
# Revision 1.23  2009/12/10 18:54:21  pwagner
# Must repack before augmenting, not after
#
# Revision 1.22  2009/10/27 21:09:12  pwagner
# Added NETCDFAUGMENTing hdfeos std products
#
# Revision 1.21  2009/05/01 22:10:38  honghanh
# Change the initialization of l2cf back to the old way because environment
# variables defined in .env file is not avaiable to mlsl2p.sh.
#
# Revision 1.20  2009/05/01 20:40:22  honghanh
# Use L2CFVERSION to supply l2cf variable
#
# Revision 1.19  2009/02/13 17:37:05  pwagner
# Running mlspgs automatically prints license text
#
# Revision 1.18  2008/07/11 20:20:18  pwagner
# Fixed bugs regarding products in outputs subdirectory
#
# Revision 1.17  2008/07/10 17:39:50  pwagner
# Handle case when product files are in outputs subdirectory
#
# Revision 1.16  2007/05/10 23:40:34  pwagner
# Used ulimit to increase tiny stacksize so Intel-built mlsl2 can finish
#
# Revision 1.15  2007/03/23 00:32:00  pwagner
# By default no longer dumps slave stdouts into masters
#
# Revision 1.14  2007/02/09 21:31:19  pwagner
# Sets environmental variable as work-around to Lahey 6.2 bug
#
# Revision 1.13  2006/10/19 18:31:10  pwagner
# Now saves chunk stats at end of output
#
# Revision 1.12  2006/04/03 23:10:01  pwagner
# Fixed serious and various bugs
#
# Revision 1.11  2006/03/23 19:22:42  pwagner
# repack with gzip compression all hdf5/hdfeos5 product files
#
# Revision 1.10  2005/12/22 19:07:58  pwagner
# Imitated l1b repacking to defragment forward model files
#
# Revision 1.9  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.8  2004/03/05 19:08:04  pwagner
# Made --cat the default option
#
# Revision 1.7  2004/01/07 17:26:09  pwagner
# Merged in sips-friendly changes
#
# Revision 1.6  2003/12/11 23:07:57  pwagner
# May check each slaves ident against master to verify pge versions are the same
#
# Revision 1.5  2003/11/15 00:45:08  pwagner
# Uses --checkPaths to perform extra preflight checks
#
# Revision 1.4  2003/10/22 23:00:17  pwagner
# Catenates each slaves log file to mlsl2.log at end
#
# Revision 1.3  2003/10/15 17:01:18  pwagner
# Added slv option
#
# Revision 1.2  2003/09/11 20:15:57  pwagner
# Allow OTHEROPTS environmental variable pass-through to mlsl2
#
# Revision 1.1  2003/08/01 16:46:31  pwagner
# First commit
#
