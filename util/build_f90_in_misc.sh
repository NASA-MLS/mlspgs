#!/bin/sh
# build_f90_in_misc.sh

# --------------- build_f90_in_misc.sh help
# Usage: build_f90_in_misc.sh [options] file [files]
# given fortran file(s) named as arg(s) on command line
# builds an executable by
# (1) hiding any .f90, .f, .h, .f9h or .c curently in tests/misc
# (2) copying the arg(s) to tests/misc
# (3) cd tests/misc; make depends; make
# (4) mv \$MLSCONFG/test ../../bin/\$MLSCONFG/test

# options
# -h[elp]            brief summary of usage and options; then exit
# -p prog_name       name resulting program prog_name instead of test
# -d prog_path       install resulting program in prog_path
# -m test_dir_name   use tests/test_dir_name instead of tests/misc
# -M MYMAKE          use MYMAKE for make command instead of make
# -t test_dir_path   use test_dir_path/misc instead of tests/misc
# -c MLSCONFG        use arg for MLSCONFG instead of any current setting
# -C MLSCFILE        use arg for MLSCFILE instead of .configure
# -CC compiler       use arg for c compiler instead of cc
# -FC compiler       use arg for Fortran compiler instead of MLSCONFG
# -i INC_PATHS       let arg override usual value for INC_PATHS
#                     used for search paths when compiling
#                     e.g., -i "../ /software/SRC..."
#                     these values can be relative or absolute
#                     relative means w.r.t the current working directory
#                     e.g., '../other/directory'
#                     absolute means the path begins with "/"
#                     e.g., '/usr/lib/gl'
# -I PATH            add "-I PATH" to EXTRA_PATHS
#                     you may have multiple occurrences of '-I PATH'
#                     that assemble a longer EXTRA_PATHS
#                     e.g., '-I path1 -I path2 ..'
#                     however '-i INC_PATHS' would override these
#                     therefore don't use the -i  and -I options simultaneously
# --------------- End build_f90_in_misc.sh help

# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# Known bugs and limitations
# (1) Doesn't pick up MLSCONFG from .configure if run directly from command-line
# (2) Doesn't allow you to override other .configure or make variables
#----------------------- Implementation -----------------------
#

#****************************************************************
#                                                               *
#                   * * * Functions  * * *                      *
#                                                               *
#****************************************************************

#
# Function to mkdir each step or subpath along a path
#

mkdir_cascade()
{

   # Does path start with "/" or not?
   # (in other words is it absolute or relative?)

   is_absolute=`echo "$*" | grep '^\/'`

   the_pn=`echo $* | sed 's;/; ;g'`

   if [ "$is_absolute" = "" ]
   then
      the_path=
   else
      the_path="/"
   fi

   for p in $the_pn
   do
      the_path="${the_path}$p/"
      if [ ! -d "$the_path" ]
      then
         mkdir "$the_path"
   #      echo mkdir "$the_path"
      fi
   done

}

#
# Function to repair each relative subpath among the args
#

repair_subpaths()
{

   repaired_subpaths=
   # 1st delete all -I prefixes
   unprefixed=`$REECHO -d $@`

   # Trivial case ($# = 0)
   if [ "$unprefixed" != "" ]
   then
      for file in $unprefixed
      do

         # Does path start with "/" or not?
         # (in other words is it absolute or relative?)

         is_absolute=`echo "$file" | grep '^\/'`

         the_pn=`echo $file | sed 's;/; ;g'`

         if [ "$is_absolute" = "" ]
         then
         # Relative--Need to repair path
            repaired_subpaths="`pwd`/$file $repaired_subpaths"
        else
         # Absolute path--no need to repair
            repaired_subpaths="$file $repaired_subpaths"
         fi

      done
   # Now that the sub_paths are repaired, re-prefix them
   repaired_subpaths=`$REECHO -d -prefix=-I $repaired_subpaths`
   fi

}

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#****************************************************************

if [ $# -lt 1 ]
then
   echo "Usage: build_f90_in_misc.sh [options] file [files]"
   exit 1
elif [ "$MLSCONFG" = "" ]
then
#   echo "Sorry--you must first configure the mls software"
#   echo "(or at least define MLSCONFG)"
#   exit 1
   NEED_MLSCONFG="yes"
else
   NEED_MLSCONFG="no"
fi

# Initialize settings to defaults
DEEBUG=on
BUILD=on
prog_name=test
test_dir_name=misc
hidden_dir_name=hideme
# The following paths work well starting from, say, mlspgs/l2
test_dir_path=../tests
prog_path=../bin/$MLSCONFG
override_INC_PATHS="no"
INC_PATHS=""
EXTRA_PATHS=""
MYMAKE=make
MYMAKEOPTS=

me="$0"
my_name=build_f90_in_misc.sh
args_dir=`pwd`
REECHO="`echo $0 | sed 's/build_f90_in_misc.sh/reecho.sh/'`"

arglist=

#
# Get arguments from command line
#

while [ "$1" != "" ] ; do

    case "$1" in
	-p )
	    prog_name=$2
	    shift
	;;
	-d )
	    prog_path=$2
	    shift
	;;
	-m )
	    test_dir_name=$2
	    shift
	;;
	-M )
	    MYMAKE=$2
	    shift
	;;
	-t )
	    test_dir_path=$2
	    shift
	;;
	-c )
	    MLSCONFG=$2
       NEED_MLSCONFG="no"
	    shift
	;;
	-C )
	    MYMAKEOPTS="MLSCFILE=$2 $MYMAKEOPTS"
	    shift
	;;
	-CC )
	    MYMAKEOPTS="CC=$2 $MYMAKEOPTS"
	    shift
	;;
	-FC )
	    MYMAKEOPTS="FC=$2 $MYMAKEOPTS"
	    shift
	;;
	-I )
	    EXTRA_PATHS="-I $2 $EXTRA_PATHS"
	    shift
	;;
	-i )
	    INC_PATHS=$2
       override_INC_PATHS="yes"
	    shift
	;;
	-h | -help )
	   sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
   		| sed -n 's/^..//p' | sed '1 d; $ d'
        exit
	;;
	* )
      arglist="$arglist $1"
	;;
    esac
    shift

done

if [ $DEEBUG = "on" ]
then
   echo "Done with parsing control line"
   echo "prog_name: $prog_name"
   echo "prog_path: $prog_path"
   echo "test_dir_name: $test_dir_name"
   echo "test_dir_path: $test_dir_path"
   echo "make command: $MYMAKE"
   echo "NEED_MLSCONFG: $NEED_MLSCONFG"
   echo "MLSCONFG: $MLSCONFG"
   echo "override_INC_PATHS: $override_INC_PATHS"
   echo "INC_PATHS: $INC_PATHS"
   echo "EXTRA_PATHS: $EXTRA_PATHS"
   echo "extra make options: $MYMAKEOPTS"
   echo "arglist: $arglist"
fi

# Check on args--need any more? Are they self-consistent?
if [ "$NEED_MLSCONFG" = "yes" ]
then
   echo "Sorry--MLSCONFG must be defined"
   echo "(Try running with -h option for help)"
   exit 1
elif [ "override_INC_PATHS" = "yes" -a "$EXTRA_PATHS" != "" ]
then
   echo "Warning--attempt to specify both -i and -I arguments"
   exit 1
fi

# Check whether misc directory exists and that we have write permission
if [ ! -w "$test_dir_path/$test_dir_name" ]
then
   echo "Sorry--you must have write permission to $test_dir_path/$test_dir_name"
   exit 1
elif [ -d "$prog_path" -a ! -w "$prog_path" ]
then
   echo "Sorry--you must have write permission to $test_dir_path/$test_dir_name"
   exit 1
elif [ "$arglist" = "" ]
then
   echo "Sorry--you must have at least one source file named as an arg"
   exit 1
fi

# Check whether a hidden directory exists yet; if not make one
if [ ! -d "$test_dir_path/$test_dir_name/$hidden_dir_name" ]
then
   if [ "$DEEBUG" = "on" ]
   then
      echo mkdir $test_dir_path/$test_dir_name/$hidden_dir_name
   fi
   if [ "$BUILD" = "on" ]
   then
      mkdir $test_dir_path/$test_dir_name/$hidden_dir_name
   fi
fi

cd $test_dir_path/$test_dir_name

# (1) Hide any pre-existing source files
if [ $DEEBUG = "on" ]
then
   echo "Hiding any pre-existing source files"
fi
for source_suffix in c f f90 h f9h
do
   pre_files=*.$source_suffix
   for file in $pre_files
   do
      if [ -f $file ]
      then
         if [ "$BUILD" = "on" ]
         then
            mv $file $hidden_dir_name
         fi
         if [ "$DEEBUG" = "on" ]
         then
            echo mv $file $hidden_dir_name
         fi
      fi
   done
done

cd $args_dir

# (2) Copy files to misc directory
if [ $DEEBUG = "on" ]
then
   echo "Copying $arglist to misc"
fi
for file in $arglist
do
   if [ -r $file ]
   then
      if [ "$BUILD" = "on" ]
      then
         cp $file $test_dir_path/$test_dir_name
      fi
      if [ "$DEEBUG" = "on" ]
      then
         echo cp $file $test_dir_path/$test_dir_name
      fi
   fi
done

# (3) Build the executable
if [ $DEEBUG = "on" ]
then
   echo "Building the executable"
fi

if [ "$EXTRA_PATHS" != "" ]
then
   repair_subpaths $EXTRA_PATHS
   EXTRA_PATHS=$repaired_subpaths
   if [ "$BUILD" = "on" ]
   then
      cd $test_dir_path/$test_dir_name
#      make depends ghostbuster
#      make EXTRA_PATHS="$EXTRA_PATHS"
      $MYMAKE depends ghostbuster $MYMAKEOPTS
      $MYMAKE EXTRA_PATHS="$EXTRA_PATHS" $MYMAKEOPTS
   fi
   if [ "$DEEBUG" = "on" ]
   then
      echo $MYMAKE EXTRA_PATHS="$EXTRA_PATHS" $MYMAKEOPTS
   fi
elif [ "$override_INC_PATHS" = "no" ]
then
   if [ "$BUILD" = "on" ]
   then
      cd $test_dir_path/$test_dir_name
#      make depends ghostbuster
#      make
      $MYMAKE depends ghostbuster $MYMAKEOPTS
      $MYMAKE $MYMAKEOPTS
   fi
else
   repair_subpaths $INC_PATHS
   INC_PATHS=$repaired_subpaths
   if [ "$BUILD" = "on" ]
   then
      cd $test_dir_path/$test_dir_name
#      make depends ghostbuster
#      make INC_PATHS="$INC_PATHS"
      $MYMAKE depends ghostbuster $MYMAKEOPTS
      $MYMAKE INC_PATHS="$INC_PATHS" $MYMAKEOPTS
   fi
   if [ "$DEEBUG" = "on" ]
   then
      echo $MYMAKE INC_PATHS="$INC_PATHS" $MYMAKEOPTS
   fi
fi

# (4) Install
if [ $DEEBUG = "on" ]
then
   echo "Installing $prog_name in $prog_path"
fi
cd $args_dir

# Check whether the final directory exists yet; if not make it
if [ ! -d "$prog_path" -a "$BUILD" = "on" ]
then
   mkdir_cascade $prog_path
fi

if [ -x "$test_dir_path/$test_dir_name/$MLSCONFG/test" -a "$BUILD" = "on" ]
then
   mv $test_dir_path/$test_dir_name/$MLSCONFG/test $prog_path/$prog_name
fi

exit 0

# $Log$
# Revision 1.8  2001/08/28 18:26:32  pwagner
# Added -M option; no more bare make
#
# Revision 1.7  2001/08/17 23:15:26  pwagner
# EXTRA_PATHS now implemented; relative paths w.r.t. cwd
#
# Revision 1.6  2001/08/14 23:45:41  pwagner
# CHanged file requirement before hiding from -w to just -r
#
# Revision 1.5  2001/08/13 23:57:53  pwagner
# Forgot ot turn DEEBUG off
#
# Revision 1.4  2001/08/13 23:25:25  pwagner
# Added rules for compiling .c files, too
#
# Revision 1.3  2001/08/10 23:51:17  pwagner
# Cosmetic changes only
#
# Revision 1.2  2001/08/10 17:43:19  pwagner
# Fixed -h(elp) option; general housekeeping
#
# Revision 1.1  2001/07/26 22:49:56  pwagner
# First commit
#
