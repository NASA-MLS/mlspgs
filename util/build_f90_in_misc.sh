#!/bin/sh
# build_f90_in_misc.sh

# --------------- build_f90_in_misc.sh help
# Usage: build_f90_in_misc.sh [options] file [files]
# given fortran file(s) named as arg(s) on command line
# builds an executable by
# (1) hiding any .f90, .f, or .c curently in tests/misc
# (2) copying the arg(s) to tests/misc
# (3) cd tests/misc; make depends; make
# (4) mv \$MLSCONFG/test ../../bin/\$MLSCONFG/test

# options
# -h                 helpful summary of usage and options
# -p prog_name       name resulting program prog_name instead of test
# -d prog_path       install resulting program in prog_path
# -m test_dir_name   use tests/test_dir_name instead of tests/misc
# -t test_dir_path   use test_dir_path/misc instead of tests/misc
# --------------- End build_f90_in_misc.sh help

# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#
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
   echo "Sorry--you must first configure the mls software"
   echo "(or at least define MLSCONFG)"
   exit 1
fi

# Initialize settings to defaults
DEEBUG=off
prog_name=test
test_dir_name=misc
hidden_dir_name=hideme
# The following paths work well starting from, say, mlspgs/l2
test_dir_path=../tests
prog_path=../bin/$MLSCONFG

me="$0"
my_name=build_f90_in_misc.sh
args_dir=`pwd`

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

	-t )
	    test_dir_path=$2
	    shift
	;;

	-h | -help )
	   sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
   		| sed -n 's/^..//p' | sed '$ d'
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
   echo "arglist: $arglist"
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
   mkdir $test_dir_path/$test_dir_name/$hidden_dir_name
fi

cd $test_dir_path/$test_dir_name

# (1) Hide any pre-existing source files
if [ $DEEBUG = "on" ]
then
   echo "Hiding any pre-existing source files"
fi
for source_suffix in c f f90
do
   pre_files=*.$source_suffix
   for file in $pre_files
   do
      if [ -w $file ]
      then
         mv $file $hidden_dir_name
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
   if [ -w $file ]
   then
      cp $file $test_dir_path/$test_dir_name
   fi
done

# (3) Build the executable
if [ $DEEBUG = "on" ]
then
   echo "Building the executable"
fi
cd $test_dir_path/$test_dir_name
make depends ghostbuster
make

# (4) Install
if [ $DEEBUG = "on" ]
then
   echo "Installing $prog_name in $prog_path"
fi
cd $args_dir

# Check whether the final directory exists yet; if not make it
if [ ! -d "$prog_path" ]
then
   mkdir_cascade $prog_path
fi

mv $test_dir_path/$test_dir_name/$MLSCONFG/test $prog_path/$prog_name

exit 0

# $Log$
# Revision 1.1  2001/07/26 22:49:56  pwagner
# First commit
#
