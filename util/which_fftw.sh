#!/bin/sh
#which_fftw.sh

# --------------- which_fftw.sh help
# Use (1):
# Creates fftw_link_line to build mlspgs level 3d program (default)
# Use (2):
# Echoes version number of gcc compiler; e.g. '2.96' or '3.2' (-gcc option)
# Use (3):
# object needed to build static executable; e.g. './temp/ctype-info.o' (-static option)
# 
# Usage:
# which_fftw.sh [options] FFTW_ROOT [FFTW_PREC]
#
#    O p t i o n s   a n d   a r g s
# -h[elp]       print brief help message; exit
# -gcc          Use (2) above
# -static       Use (3) above
# FFTW_ROOT     a directory path where the fftw libraries might be found
#                These two libraries must have names matching the patterns
#                1st library: lib[ sd]rfftw.a
#                2nd library: lib[ sd]fftw.a
# FFTW_PREC     if present, either 's' or 'd' specifying which precision
#                library to choose
#                (if absent, the order of preference will be {none, d, s}
#                in other words if librfftw.a exists, use it
#                elseif libdrfftw.a exists, use it
#                elseif libsrfftw.a exists, use it
#                else return nothing

# Result:
# Writes fftw_link_line to stdout; similar to the following
#   -L${FFTW_ROOT} -lrfftw -lfftw
# Also creates a file containing fftw_link_message in the cwd;
#  similar to the following
#    Building mlsl3d with default-precision fftw

# Note:
# (1) The option(s) marked with "-", if present,
#     must precede the args on the command line
# (2) It is assumed that the shell script reecho.sh exists
#      and is in the same directory as which_fftw.sh
#
# --------------- End which_fftw.sh help
# Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
me="$0"
my_name=which_fftw.sh
# What to echo as need object for building static executable
# (if Lahey compiler < v6.2 and Red Hat > v8 )
my_stobject="/users/pwagner/lib/lf95/ctype-info.o"
DEEBUG=off
USE="1"

# The following is echoed if we can't any fftw libraries in FFTW_ROOT
nothing='-L${FFTW_ROOT} -ldrfftw -ldfftw'

if [ $DEEBUG = "on" ]
then
   echo "Called me as $0"
   echo "with args $@"
fi

REECHO="`echo $0 | sed 's/which_fftw.sh/reecho.sh/'`"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
       ;;
    -gcc )
       USE="2"
       shift
       ;;
    -static )
       USE="3"
       shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done

if [ $USE != "1" ]
then
   stobject=''
#   - - - U s e  ( 2 ) o r ( 3 ) - - -
   test_295=`gcc -v 2>&1 | grep -i '2\.95'`
   test_296=`gcc -v 2>&1 | grep -i '2\.96'`
   # test_31=`gcc -v 2>&1 | grep -i '3\.1'`
   # test_32=`gcc -v 2>&1 | grep -i '3\.2'`
   # test_322=`gcc -v 2>&1 | grep -i '3\.2\.2'`
   # test_323=`gcc -v 2>&1 | grep -i '3\.2\.3'`
   gcc_version=
   if [ "$test_295" != "" ]
   then
     gcc_version="2.95"
   elif [ "$test_296" != "" ]
   then
     gcc_version="2.96"
  #elif [ "$test_31" != "" ]
  #then
  #  gcc_version="3.1"
  #elif [ "$test_322" != "" ]
  #then
  #  gcc_version="3.2.2"
  #  stobject="$my_stobject"
  #elif [ "$test_323" != "" ]
  #then
  #  gcc_version="3.2.3"
  #  stobject="$my_stobject"
  #elif [ "$test_32" != "" ]
  #then
  #  gcc_version="3.2"
   else
     gcc_version=`gcc -dumpversion`
     #exit 1
     case "$gcc_version" in
        3.1)
        stobject=
        ;;
        3.2)
        stobject=
        ;;
        3.2.2)
        stobject="$my_stobject"
        ;;
        3.2.3)
        stobject="$my_stobject"
        ;;
        *)
        stobject="$my_stobject"
        ;;
      esac
   fi
   if [ "$USE" = "2" ]
   then
     echo $gcc_version
   else
     echo $stobject
   fi
   exit 0
fi

#   - - - U s e  ( 1 ) - - -
FFTW_ROOT=$1
shift

FFTW_PREC=$1

if [ $DEEBUG = "on" ]
then
   echo "FFTW_ROOT: $FFTW_ROOT"
   echo "FFTW_PREC: $FFTW_PREC"
fi

rm -f fftw_link_message
# No specified precision
test=`$REECHO $FFTW_ROOT/librfftw.a $FFTW_ROOT/libfftw.a | wc -w`
if [ $DEEBUG = "on" ]
then
   echo "First test: $test"
fi
if [ `expr "$test"` = 2 -a "$FFTW_PREC" = "" ]
then
  echo "Building mlsl3[dm] with default-precision fftw" > fftw_link_message
  echo '-L${FFTW_ROOT} -lrfftw -lfftw'
  exit
fi

# double precision
test=`$REECHO $FFTW_ROOT/libdrfftw.a $FFTW_ROOT/libdfftw.a | wc -w`
if [ $DEEBUG = "on" ]
then
   echo "dp test: $test"
fi
if [ `expr "$test"` = 2 -a "$FFTW_PREC" != "s" ]
then
  echo "Building mlsl3[dm] with double-precision fftw" > fftw_link_message
  echo '-L${FFTW_ROOT} -ldrfftw -ldfftw'
  exit
fi

# single precision (and last chance)
test=`$REECHO $FFTW_ROOT/libsrfftw.a $FFTW_ROOT/libsfftw.a | wc -w`
if [ $DEEBUG = "on" ]
then
   echo "sp test: $test"
fi
if [ `expr "$test"` = 2 -a "$FFTW_PREC" != "d" ]
then
  echo "Building mlsl3[dm] with single-precision fftw" > fftw_link_message
  echo '-L${FFTW_ROOT} -lsrfftw -lsfftw'
  exit
fi

# uh-oh, no qualified libraries found, so echo $nothing (unless debugging)
if [ $DEEBUG = "on" ]
then
   echo "uh-oh, no qualified libraries found, so echo nothing"
   echo "No fftw libraries were found in $FFTW_ROOT"         
   echo "(assuming your precision $FFTW_PREC)"               
   echo "You probably have to reset FFTW_ROOT in .configure" 
   echo "Do that by 'make configure_pvm'"                    
else
   echo $nothing
fi
echo "No fftw libraries were found in $FFTW_ROOT" > fftw_link_message           
echo "(assuming your precision $FFTW_PREC)" >> fftw_link_message                 
echo "You probably have to reset FFTW_ROOT in .configure" >> fftw_link_message   
echo "Do that by 'make configure_pvm'" >> fftw_link_message                      
exit
# $Log$
# Revision 1.6  2004/03/17 00:19:35  pwagner
# Use(3) to fix building static Lahey executables when Red Hat>8
#
# Revision 1.5  2004/02/06 18:31:34  pwagner
# Added gcc version 3.2.3
#
# Revision 1.4  2003/05/12 22:08:46  jonathan
# added v3.2.2
#
# Revision 1.3  2002/10/03 23:06:10  pwagner
# Added Use(2) to return gcc version
#
# Revision 1.2  2001/10/25 17:53:43  pwagner
# Slightly more robust given wrong FFTW_ROOT
#
# Revision 1.1  2001/10/09 20:51:22  pwagner
# First commit
#
