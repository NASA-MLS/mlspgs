#!/bin/sh
# mark_as_uptodate.sh

#
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# --------------- mark_as_uptodate.sh help
#Main use:
#mark targets as uptodate in current working directory
#and optionally in any prerequisite directories.
#
#Usage:
#mark_as_uptodate.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [pd1] .. [pdn]
#
#Result:
#make will show targets uptodate with prerequisites
#so subsequent makes should do nothing
#
#Special touch use (invoked by -t option):
#if target exists--just touch target
#else try "rm -f 1st_arg; $(MYMAKE) 1st_arg"
#Example of special touch usage:
#mark_as_uptodate.sh -t -T target_name alternate_target_name
#
#Special compile-touch use (invoked by -c option):
#infer object_name--then touch it
#Example of special compile usage:
#mark_as_uptodate.sh -c file_name
#
#Result:
#target_name will be touched or rebuilt if necessary, and thus brought uptodate
#also exit status will be passed to calling routine
#
#    O p t i o n s
# -h[elp]            brief summary of usage and options; then exit
# -T target_name     main target_name of make
#                      (may be repeated for multiple main targets)
# -M MYMAKE          use MYMAKE for make command instead of make
# -cg MLSCONFG       use arg for MLSCONFG instead of any current setting
# -C MLSCFILE        use arg for MLSCFILE instead of .configure
# -fyes              begin with MARK_ALL_AS_UPTODATE set to "yes"
# -fno               begin with MARK_ALL_AS_UPTODATE set to "no"
# -mod ext           use .ext instead of .mod for module name extension
# -syes              skip going through prerequisite directories
# -sno               don't skip going through prerequisite directories
# -t                 special touch use: touch target_name 
# -c                 special compile use: touch inferred object_name
# -o target_name     special link use: touch target
# -v                 verbose--note marking up to date
#Notes:
#(1)The options -x must appear before the prerequisite directories pdn
#(2)The options -[T M cg C mod] require an arg
#(3)The -fxxx are mutually exclusive; if neither, or -f"", 
#   MARK_ALL_AS_UPTODATE="" (the default)
#(4)The -sxxx are mutually exclusive; if neither, or -s"", 
#   SKIP_PDS="no" (the default)
#(5)The prerequisite directories pdn must be ordered so that they obey:
#   pd1 << pd2 << .. << pdn << (current directory)
#   where the symbol << means "is a prerequisite of"
#(6)Marking is conditional: at least one of the following must be met
#    in order for marking to actually take place:
#   (a) MARK_ALL_AS_UPTODATE=yes
#   (b) a file named newAifBdiff.out exists in the target's directory
#(7)If -T target_name is supplied, target_name will be touched; otherwise not
# --------------- End mark_as_uptodate.sh help
# Bugs and limitations:
# Assumes the script split_path.sh exists in the same directory
# 
# Purpose: A detailed explanation
# (main usage)
# make calls me
# First, in order to forestall triggering
# cascades of recompilation when modules are linked in an intricate
# tree of "USE"-based interdependencies, object files may be rebuilt
# without rebuilding .mod files. Afterward, to bring everything forward
# in time, some targets may need to marked as up to date. That's my job.

# (special touch use)
# build commands for .o targets with .mod prerequisites; e.g.,
# tree_checker.o: tree_checker.mod tree_checker.f90
#	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
#   -T tree_checker.o tree_checker.mod 
# Solves problem of teaching make that .o depends on .mod
# w/o requiring a duplicate recompilation to build it

# (special compile touch use)
# invoked when I call make myself to bring entire directory
# full of targets up to date. I do it with the following command
#   make FC=$me
# so all the build commands that look like
#   $(FC) -c $(DUSTY) $(INC_PATHS) file_name.f90
# get passed back to me as if I were a compiler;
# what I do wearing my compiler-hat is to touch the already-built target
# (or exit with an error if the target is non-existant)

# (special link use)
# invoked when I call make (like special compile touch use)
# when make reaches the link statement that looks like
#   $(FC) $(LDOPTS) -o target $(ALLOBJS) \
# simply touches target
# (or exit with an error if the target is non-existant)

#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in pace of "."

get_unique_name()
{

   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
      fi
   # Is $HOST defined?
      if [ "$HOST" != "" ]
      then
         our_host_name="$HOST"
      elif [ "$HOSTNAME" != "" ]
      then
         our_host_name="$HOSTNAME"
      else
         our_host_name="host"
      fi
    #  echo $our_host_name
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
}

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
}

#------------------------------- mark_files ------------
#
# Function to mark files
# usage: mark_files arg1 [arg2] ..

# The main doing procedure in maras_uptodate
# After (optionally) touching modules and objects
#   ( a preparatory bit of throat-clearing fallen out of favor )
# call make
# Once we trusted that make with the -t option would do the job, but
# now we rely on our own strength and resolve
# Note the settings in the calls to make below:
# (1) MARK_ALL_AS_UPTODATE=no
#     so make won't call me back to go through this again (no infinite loop)
# (2) FC="$me"
#     so compiler and linker commands get passed back to me to handle
# (3) LDOPTS=
#     so the linker command is simplified to $(FC) -o target [ignored stuff]
# (4) FAFTER= and LAFTER=
#     so Van's thing with FAFTER=2>&1|hl|grep -v "Id is set but" won't bomb

# Anomalies and bugs:
# Why the ugly DOUBLE_MAKE?
# Sometimes a single invocation left a prerequisite subdirectory
# not uptodate (Why?--this needs study)

mark_files()
{
    TOUCH_OBJECTS="no"
    TOUCH_MODULES="no"
    DOUBLE_MAKE="yes"
    PRINT_STDOUT="no"
    temp1=`get_unique_name mark`
    if [ "$PRINT_STDOUT" != "yes" ]            
    then  
      echo "mark_files record $temp1"
      echo "mark_files record" > $temp1
    fi
#   extant_files *.o *.mod
#   echo *.${ext}
   extant_files *.${ext}
#   echo "extant_files: $extant_files_result"
   if [ "$extant_files_result" != "" -a "$TOUCH_MODULES" = "yes" ]            
   then  
     touch $extant_files_result
     if [ "$verbose" != "" ]            
     then  
        echo "touched: $extant_files_result"                      
     fi                                 
   fi
   extant_files *.o
   if [ "$extant_files_result" != "" -a "$TOUCH_OBJECTS" = "yes" ]            
   then  
     touch $extant_files_result
     if [ "$verbose" != "" ]            
     then  
        echo "touched: $extant_files_result"                      
        echo "Now will $MYMAKE in `pwd`"                      
     fi                                 
   fi
   if [ "$verbose" != "" ]            
   then  
      echo "Now will $MYMAKE in `pwd`"                        
   fi                                 
   if [ "$DOUBLE_MAKE" = "yes" ]            
   then  
      make_times="1 2"                        
   else                                 
      make_times="1"                        
   fi                                 
   for the_time in $make_times
   do
     if [ "$MLSCONFG" != "" -a "$MLSCFILE" != "" ]
     then
#       $MYMAKE -t MLSCONFG="$MLSCONFG" MLSCFILE="$MLSCFILE" \
#           MARK_ALL_AS_UPTODATE=no
        if [ "$PRINT_STDOUT" = "yes" ]
        then
          $MYMAKE MLSCONFG="$MLSCONFG" MLSCFILE="$MLSCFILE" \
             MARK_ALL_AS_UPTODATE=no FC="$me" LDOPTS=
        else
          $MYMAKE MLSCONFG="$MLSCONFG" MLSCFILE="$MLSCFILE" \
             MARK_ALL_AS_UPTODATE=no FC="$me" LDOPTS= \
             FAFTER= LAFTER= >> $temp1
        fi
     elif [ "$MLSCONFG" != "" ]
     then
#        $MYMAKE -t MLSCONFG="$MLSCONFG" MARK_ALL_AS_UPTODATE=no
        if [ "$PRINT_STDOUT" = "yes" ]
        then
          $MYMAKE MLSCONFG="$MLSCONFG" MARK_ALL_AS_UPTODATE=no \
           FC="$me" LDOPTS=
        else
          $MYMAKE MLSCONFG="$MLSCONFG" MARK_ALL_AS_UPTODATE=no \
           FC="$me" LDOPTS= FAFTER= LAFTER= >> $temp1
        fi
     elif [ "$MLSCFILE" != "" ]
     then
#       $MYMAKE -t MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no
        if [ "$PRINT_STDOUT" = "yes" ]
        then
          $MYMAKE MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no \
           FC="$me" LDOPTS=
        else
          $MYMAKE MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no \
           FC="$me" LDOPTS= FAFTER= LAFTER= >> $temp1
        fi
     else
#       $MYMAKE -t MARK_ALL_AS_UPTODATE=no
#        echo $MYMAKE MARK_ALL_AS_UPTODATE=no FC="$me" LDOPTS=
        if [ "$PRINT_STDOUT" = "yes" ]
        then
          $MYMAKE MARK_ALL_AS_UPTODATE=no FC="$me" LDOPTS=
        else
          $MYMAKE MARK_ALL_AS_UPTODATE=no FC="$me" LDOPTS= \
           FAFTER= LAFTER= >> $temp1
        fi
     fi
     return_status=`expr $?`
     if [ "$return_status" != "$NORMAL_STATUS" ]; then
        echo "Error in mark_files; see $temp1 for details"
        exit 1
     elif [ "$PRINT_STDOUT" != "yes" ]           
     then  
       rm -f "$temp1"
     fi
   done
}

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

DEBUG=0
#     ^  -- set this to 1 to echo rather than execute all commands
if [ "$DEBUG" = "1" ]
then
  run_mode="echo"
else
  run_mode=""
fi

PRINT_TOO_MUCH=0
#              ^  -- set this to 1 if patient enough to print more
if [ "$PRINT_TOO_MUCH" = "1" ]
then
  echo "$0 called with args $@"
fi

# variable        meaning            default
# target_name    of make            ""
# me                                "$0"
# my_name                           mark_as_uptodate.sh
# MYMAKE                            make
# MLSCONFG                          (as found in Makefile)
# MLSCFILE                          (as found in Makefile)
# MARK_ALL_AS_UPTODATE              ""
# special_use                       no (special touch use)
# compile_use                       no (special compile use)
# verbose                           no
# SKIP_PDS                          no
# ext                               mod
#----------------------- Implementation -----------------------

target_name=""
last_target_name=""
record_file="newAifBdiff.out"
me="$0"
my_name=mark_as_uptodate.sh
# $the_splitter is split_path with me's path prepended
the_splitter="`echo $0 | sed 's/mark_as_uptodate/split_path/'`"
MYMAKE=make
MLSCONFG=""
MLSCFILE=""
MARK_ALL_AS_UPTODATE=""
special_use=no
compile_use=no
link_use=no
verbose=no
ext=mod
SKIP_PDS=no
NORMAL_STATUS=0
return_status=0
#
# Get arguments from command line
#

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    case "$1" in
	-T )
	    last_target_name="$2"
	    target_name="$2 $target_name"
	    shift
	    shift
	;;
	-o )
	    last_target_name="$2"
	    target_name=$2
	    special_use=yes
	    link_use=yes
	    shift
	    shift
	;;
	-M )
	    MYMAKE=$2
	    shift
	    shift
	;;
	-C )
	    MLSCFILE=$2
	    shift
	    shift
	;;
	-cg )
	    MLSCONFG=$2
	    shift
	    shift
	;;
	-c )
	    compile_use=yes
	    shift
	;;
	-fyes )
	    MARK_ALL_AS_UPTODATE=yes
	    shift
	;;
	-fno )
	    MARK_ALL_AS_UPTODATE=no
	    shift
	;;
	-mod )
	    ext=$2
	    shift
	    shift
	;;
	-syes )
	    SKIP_PDS=yes
	    shift
	;;
	-s* )
	    SKIP_PDS=no
	    shift
	;;
	-f* )
	    MARK_ALL_AS_UPTODATE=
	    shift
	;;
	-t )
	    special_use=yes
	    shift
	;;
	-v )
	    verbose=yes
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

if [ "$PRINT_TOO_MUCH" = "1" ]
then
  verbose=yes
  echo "MYMAKE $MYMAKE"
  echo "MLSCONFG $MLSCONFG"
  echo "MLSCFILE $MLSCFILE"
  echo "MARK_ALL_AS_UPTODATE? $MARK_ALL_AS_UPTODATE"
  echo "mod file name ext $ext"
  echo "target_name(s) $target_name"
  echo "last target_name $last_target_name"
  echo "special touch use? $special_use"
  echo "special compile use? $compile_use"
  echo "skip prerequisite directories? $SKIP_PDS"
  echo "remaining args (prereq dirs): $@"
fi
if [ "$verbose" != "" ]
then
  if [ "$special_use" = "yes" ]
  then
    echo "Marking $target_name as up to date in `pwd`"
  elif [ "$compile_use" = "yes" ]
  then
    echo "Marking inferred object and module names as up to date in `pwd`"
  else
    echo "Marking OBJS and MODULES as up to date in `pwd`"
  fi
fi

if [ "$compile_use" = "yes" ]
then
	source_name=""
   while [ "$1" != "" ] ; do
   #  if [ "$verbose" != "" ]
   #  then
   #    echo "Processing option $1"
   #  fi
    case "$1" in
	*.f90 )
	    source_name=$1
       break
	;;
#	* )
#       do_I_care="no"
#	;;
    esac
     shift
   done
   if [ "$source_name" = "" ]
   then
     echo "source_name not found in mark_as_uptodate"
     exit 1
   fi
   # Now extract just file_name
#   file_name=`$the_splitter -f $source_name`
   file_name=`perl -e '$reverse=reverse("$ARGV[0]"); @parts=split("/",$reverse); $reverse=$parts[0]; $reverse=reverse($reverse); print $reverse' $source_name`
   object_name="`echo $file_name | sed 's/\.f90/\.o/'`"
   if [ ! -f "$object_name" ]
   then
     echo "object $object_name not found in mark_as_uptodate"
     exit 1
   fi
   touch "$object_name"
   if [ "$verbose" != "" ]                     
   then                                        
     echo "Touching object name $object_name"  
   fi                                          
   exit 0
elif [ "$special_use" = "yes" ]
then
   if [ -f "$last_target_name" ]
   then
      if [ "$verbose" != "" -a "$special_use" = "yes" ]
      then
         echo "touching $target_name"
      fi
      touch $target_name
   elif [ "$link_use" = "yes" ]
   then
      echo "link target $target_name not found in mark_as_uptodate.sh"
      exit 1
   else
      if [ "$verbose" != "" -a "$special_use" = "yes" ]
      then
         echo "rebuilding $target_name via $1"
      fi
      if [ -f "$1" ]
      then
         rm -f "$1"
      fi
      if [ "$MLSCONFG" != "" -a "$MLSCFILE" != "" ]
      then
         $MYMAKE $1 MLSCONFG="$MLSCONFG" MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no
      elif [ "$MLSCONFG" != "" ]
      then
         $MYMAKE $1 MLSCONFG="$MLSCONFG" MARK_ALL_AS_UPTODATE=no
      elif [ "$MLSCFILE" != "" ]
      then
         $MYMAKE $1 MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no
      else
         $MYMAKE $1 MARK_ALL_AS_UPTODATE=no
      fi
   fi
   return_status=`expr $?`
   if [ "$return_status" != "$NORMAL_STATUS" ]; then
   #   echo "exiting with status 1"
      exit 1
   else
   #   echo "exiting with status 0"
      exit 0
   fi
fi
#
# Process and mark pdn (if any) as arguments from command line
# (unless $SKIP_PDS is yes)

if [ "$SKIP_PDS" != "yes" ]
then
   while [ "$1" != "" ] ; do
     if [ "$verbose" != "" ]
     then
       echo "Marking prerequisite directory $1 as up to date"
     fi
     if [ -f "$1/$record_file" ]; then                                                
       MARK_ALL_AS_UPTODATE=yes                                                       
       rm -f "$1/$record_file"                                                        
     fi                                                                               
#    Note the logic in the following argument assignments:
#    *  MARK_ALL_AS_UPTODATE=$MARK_ALL_AS_UPTODATE  *
#    Once we know we must mark 
#    (either by command-line arg or from a prerequisite dir)
#    we continue in all dependent directories
#
#    *  SKIP_PREREQ_DIRS=yes  *
#    in calls to this shell command resulting from prereq dirs,
#    don't go down list of {\em their} prereq dirs
#
     $MYMAKE -C $1 mark_all_as_uptodate \
        MARK_ALL_AS_UPTODATE=$MARK_ALL_AS_UPTODATE \
        SKIP_PREREQ_DIRS=yes   
     shift
   done
fi

if [ -f "$record_file" ]
then
    MARK_ALL_AS_UPTODATE=yes                                                          
fi
if [ "$MARK_ALL_AS_UPTODATE" != "yes" ]
then
   if [ "$verbose" != "" ]            
   then                               
     echo "Targets already up to date--no marking necessary"
   fi                                 
   exit 0
fi
mark_files
rm -f "$record_file"                                                           
if [ -f "$last_target_name" ]
then
   if [ "$verbose" != "" ]
   then
     echo "Marking $target_name as up to date"
   fi
   touch $target_name
fi

if [ "$verbose" != "" ]              
then                                 
  echo "Targets marked up to date"   
fi                                   
exit 0

# $Log$
# Revision 1.6  2002/08/07 17:03:19  pwagner
# Fixes error when FAFTER or LAFTER redirect output
#
# Revision 1.5  2002/07/26 23:49:59  pwagner
# Faster at marking files uptodate (but still slow)
#
# Revision 1.4  2002/07/25 20:58:08  pwagner
# Improved marking up to date
#
# Revision 1.3  2002/07/23 23:18:01  pwagner
# Added -mod ext; double makes in mark_files
#
# Revision 1.2  2002/07/01 20:52:04  pwagner
# Added special usage via -t option
#
# Revision 1.1  2002/06/21 00:10:13  pwagner
# First commit
#
