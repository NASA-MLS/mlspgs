#!/bin/sh
#reecho.sh

# --------------- reecho.sh help
# (Re)echoes the command line arguments, minus any that
# (depending on option(s) selected
# are not actual files (with special properties); or
# are not directories; or
# contain the glob character '*'
# Useful to weed out formulas like *.c that
# have no matching files and so stay just '*.f'
# 'no stars => reecho (night)' is the mnemonic
#
# Usage:
# reecho.sh [opt] ..  [arg1 arg2 ..]
#
#    O p t i o n s
# -[n]f         [don't] reecho args that are ordinary files and exist
#                (default is -f)
# -[n]d         [don't] reecho args that are directories and exist
# -[n]r         [don't] reecho args that exist and you have read permission
# -[n]w         [don't] reecho args that exist and you have write permission
# -[n]x         [don't] reecho args that exist and you have execute permission
# -[n]glob      [don't] reecho any arg containing the glob character '*'
# -h[elp]       print brief help message; exit
# arg1          an arg that may or may not be reechoed
#
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede the args on the command line
# (2) The options are mutually exclusive (for the present)
# 
# Result:
# A filtered list of args, e.g. arg3 arg6 .. , that satisfy
# the conditions set up by the options
# Example:
# Say you have a directory containingg the following files
#    a.b c.d/ e.d/ f.7*   
# (where /~a directory, *~ an executable; none of which are parts of the names)
# then
#    `reecho *.b *.c`       returns "a.b"
#    `reecho *.d`           returns ""
#    `reecho -nf *.d`       returns "c.d e.d"
#    `reecho -x *`          returns "f.7"
#    `reecho -glob *.7`     returns "" 
#    `reecho -glob *.8`     returns "*.8" 
# 
# --------------- End reecho.sh help
# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

#

#------------------------------- extant_files ------------
#
# Function to return only those files among the args
# that satisfy the_opt
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
         if [ $the_opt = "-glob" ]
         then
            check=`echo $file | grep -i '\*'`
            if [ $the_sense = "yes" -a $check != "" ]
            then
                  extant_files_result="$extant_files_result $file"
            elif [ $the_sense = "no" -a $check = "" ]
            then
                  extant_files_result="$extant_files_result $file"
            fi
         elif [ $the_sense = "yes" ]
         then
            if [ $the_opt "$file" ]
            then
                  extant_files_result="$extant_files_result $file"
            fi
         else
            if [ ! $the_opt "$file" ]
            then
                  extant_files_result="$extant_files_result $file"
            fi
         fi
      done
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
me="$0"
my_name=reecho.sh
DEEBUG=off
if [ $DEEBUG = "on" ]
then
   echo "Called me as $0"
   echo "with args $@"
fi

the_opt="-f"
the_sense="yes"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -f )
       the_opt="-f"
       the_sense="yes"
       shift
       ;;
    -nf )
       the_opt="-f"
       the_sense="no"
       shift
       ;;
    -d )
       the_opt="-d"
       the_sense="yes"
       shift
       ;;
    -nd )
       the_opt="-d"
       the_sense="no"
       shift
       ;;
    -r )
       the_opt="-r"
       the_sense="yes"
       shift
       ;;
    -nr )
       the_opt="-r"
       the_sense="no"
       shift
       ;;
    -w )
       the_opt="-w"
       the_sense="yes"
       shift
       ;;
    -nw )
       the_opt="-w"
       the_sense="no"
       shift
       ;;
    -x )
       the_opt="-x"
       the_sense="yes"
       shift
       ;;
    -nx )
       the_opt="-x"
       the_sense="no"
       shift
       ;;
    -glob )
       the_opt="-glob"
       the_sense="yes"
       shift
       ;;
    -nglob )
       the_opt="-glob"
       the_sense="no"
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

if [ $DEEBUG = "on" ]
then
   echo "the_opt $the_opt"
   echo "the_sense $the_sense"
   echo "remaining args $@"
fi

extant_files "$@"
echo $extant_files_result
exit
# $Log$
