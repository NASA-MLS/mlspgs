#!/bin/sh
# mark_as_uptodate.sh
# --------------- mark_as_uptodate.sh help
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
#    O p t i o n s
# -h[elp]            brief summary of usage and options; then exit
# -T target_name     main target_name of make
# -M MYMAKE          use MYMAKE for make command instead of make
# -c MLSCONFG        use arg for MLSCONFG instead of any current setting
# -C MLSCFILE        use arg for MLSCFILE instead of .configure
# -fyes              begin with MARK_ALL_AS_UPTODATE set to "yes"
# -fno               begin with MARK_ALL_AS_UPTODATE set to "no"
# -syes              skip going through prerequisite directories
# -sno               don't skip going through prerequisite directories
# -v                 verbose--note marking up to date
#Notes:
#(1)The options -x must appear before the prerequisite directories pdn
#(2)All options except for -fxxx, -sxxx, -n and -h require a following arg
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

mark_files()
{
   extant_files *.o *.mod
   touch $extant_files_result
   if [ "$verbose" != "" ]            
   then  
      echo "touched: $extant_files_result"                        
      echo "Now will $MYMAKE in `pwd`"                        
   fi                                 
   if [ "$MLSCONFG" != "" -a "$MLSCFILE" != "" ]
   then
      $MYMAKE -t MLSCONFG="$MLSCONFG" MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no
   elif [ "$MLSCONFG" != "" ]
   then
      $MYMAKE -t MLSCONFG="$MLSCONFG" MARK_ALL_AS_UPTODATE=no
   elif [ "$MLSCFILE" != "" ]
   then
      $MYMAKE -t MLSCFILE="$MLSCFILE" MARK_ALL_AS_UPTODATE=no
   else
      $MYMAKE -t MARK_ALL_AS_UPTODATE=no
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

#
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# variable        meaning            default
# target_name    of make            ""
# me                                "$0"
# my_name                           mark_as_uptodate.sh
# MYMAKE                            make
# MLSCONFG                          (as found in Makefile)
# MLSCFILE                          (as found in Makefile)
# MARK_ALL_AS_UPTODATE              ""
# verbose                           no
# SKIP_PDS                          no
#----------------------- Implementation -----------------------

target_name=""
record_file="newAifBdiff.out"
me="$0"
my_name=mark_as_uptodate.sh
MYMAKE=make
MLSCONFG=""
MLSCFILE=""
MARK_ALL_AS_UPTODATE=""
verbose=no
SKIP_PDS=no
#
# Get arguments from command line
#

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    case "$1" in
	-T )
	    target_name=$2
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
	-c )
	    MLSCONFG=$2
	    shift
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
  echo "target_name $target_name"
  echo "skip prerequisite directories? $SKIP_PDS"
  echo "remaining args (prereq dirs): $@"
fi
if [ "$verbose" != "" ]
then
  echo "Marking OBJS and MODULES as up to date in `pwd`"
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
if [ -f "$target_name" ]
then
   if [ "$verbose" != "" ]
   then
     echo "Marking $target_name as up to date"
   fi
   touch "$target_name"
fi

if [ "$verbose" != "" ]              
then                                 
  echo "Targets marked up to date"   
fi                                   
exit 0

# $Log$
