#!/bin/sh
#which_fftw.sh

# --------------- which_fftw.sh help
# Creates fftw_link_line to build mlspgs level 3d program
# 
# Usage:
# which_fftw.sh [options] FFTW_ROOT [FFTW_PREC]
#
#    O p t i o n s   a n d   a r g s
# -h[elp]       print brief help message; exit
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
# Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
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
DEEBUG=off

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
    * )
       more_opts="no"
       ;;
    esac
done

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
# Revision 1.1  2001/10/09 20:51:22  pwagner
# First commit
#
