#!/bin/sh
#reecho.sh

# --------------- reecho.sh help
# (Re)echoes the command line arguments, minus any that
# (depending on option(s) selected)
# are not actual files (with special properties); or
# are not directories; or
# contain the glob character '*'
# Useful 
# (1) to weed out formulas like *.c that
# have no matching files and so stay just '*.c'
# rather than expanding into a list of matching files
# (2) to construct link or include lines inside Makefiles
#
# Optionally pre- or suffixes each with a specified pre- or suffix; e.g.
# Given 'reecho -prefix=-I *-dir'
# might echo '-I a-dir -I b-dir -I c-dir ..'
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
# -dir "dir"    cd to "dir" before filtering
# -dirn "dir"   ignore args; cd to "dir" then filter all files/directories there
# -lib          filter each arg "x" based on the file "libx.a" exists or ...
# -h[elp]       print brief help message; exit
# -prefix=xxx   reecho xxx in front of each arg (separated by a space)
#               e.g. '-I a-dir -I b-dir -I c-dir ..'
# -prefixn=xxx  reecho xxx in front of each arg (without a separating space)
#               e.g. '-Ia-dir -Ib-dir -Ic-dir ..'
# -suffix=xxx   reecho xxx after each arg (separated by a space)
#               e.g. 'a plop b plop c plop ...'
# -suffixn=xxx  reecho xxx after each arg (without a separating space)
#               e.g. 'Xshabam Yshabam Zshabam ..'
# -excl "bad"   exclude any arg named "bad" (before any pre- or suffixes)
#               may be repeated; e.g. -excl bad1 -excl bad2 excludes both
# arg1          an arg that may or may not be reechoed
#
# Note:
# (1) The option(s) marked with "-", if present,
#     must precede the args on the command line
# (2) The options [n](f d r w x glob) are mutually exclusive (for the present)
#      i.e. don't try 'reecho.sh -d -w *' to get all directories you have 
#      write permission to (though that is a logical improvement to make)
# 
# Result:
# A filtered list of args, e.g. arg3 arg6 .. , that satisfy
# the conditions set up by the options
# Example:
# Say you have a directory containing the following files
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
# More practical examples:
#  67%reecho.sh -lib -dir /usr/lib -prefixn=-l dfftw drfftw
#     -ldfftw -ldrfftw
#  69%util/reecho.sh -lib -dir /usr/lib -prefixn=lib -suffixn=.a dfftw drfftw
#     libdfftw.a libdrfftw.a
# --------------- End reecho.sh help
# Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
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
      for arg
      do
         arg_not_bad="yes"
         for file in $bad_args
         do
            if [ "$arg" = "$file" ]
            then
               arg_not_bad="no"
            fi
         done

         if [ "$as_lib" = "yes" ]
         then
            file="lib${arg}.a"
         else
            file=$arg
         fi
         
         if [ "$arg_not_bad" != "yes" ]
         then
#           arg a bad one--automatically excluded--no operation needed
            file=$file
         elif [ $the_opt = "-glob" ]
         then
            check=`echo $file | grep -i '\*'`
            if [ $the_sense = "yes" -a "$check" != "" ]
            then
                  extant_files_result="$extant_files_result $arg"
            elif [ $the_sense = "no" -a "$check" = "" ]
            then
                  extant_files_result="$extant_files_result $arg"
            fi
         elif [ $the_sense = "yes" ]
         then
            if [ $the_opt "$file" ]
            then
                  extant_files_result="$extant_files_result $arg"
            fi
         else
            if [ ! $the_opt "$file" ]
            then
                  extant_files_result="$extant_files_result $arg"
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
#
#   Notes
#  A logical improvement would be to allow multiple options
# among the set {[n]f [n]d [n]w [n]r [n]x}
# which should be doable by forming a list of the_opts and
# looping over them with calls to extant_files
# steadily narrowing down the files that survive being refiltered
me="$0"
my_name=reecho.sh
DEEBUG=off
if [ $DEEBUG = "on" ]
then
   echo "Called me as $0"
   echo "with args $@"
fi

new_dir=""
new_list="no"
as_lib="no"
the_opt="-f"
the_sense="yes"
more_opts="yes"
the_prefix=""
the_suffix=""
separate_prefix="yes"
separate_suffix="yes"
bad_args=""
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
    -prefix=* )
       the_prefix=`echo $1 | sed 's/-prefix=//'`
       shift
       ;;
    -prefixn=* )
       the_prefix=`echo $1 | sed 's/-prefixn=//'`
       separate_prefix="no"
       shift
       ;;
    -suffix=* )
       the_suffix=`echo $1 | sed 's/-suffix=//'`
       shift
       ;;
    -suffixn=* )
       the_suffix=`echo $1 | sed 's/-suffixn=//'`
       separate_suffix="no"
       shift
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
       ;;
    -dir )
       shift
       new_dir="$1"
       shift
       ;;
    -dirn )
       shift
       new_dir="$1"
       new_list="yes"
       shift
       ;;
    -lib )
       as_lib="yes"
       shift
       ;;
    -excl )
       shift
       bad_args="$1 $bad_args"
       shift
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
   echo "the_prefix $the_prefix"
   echo "separate_prefix? $separate_prefix"
   echo "the_suffix $the_suffix"
   echo "separate_suffix? $separate_suffix"
   echo "as lib? $as_lib"
   echo "new_dir $new_dir"
   echo "remaining args $@"
fi

if [ "$new_dir" = "" ]
then
   extant_files "$@"
elif [ "$new_list" = "yes" ]
then
   old_dir=`pwd`
   cd "$new_dir"
   my_files=`echo *`
   set `echo $my_files`
   extant_files "$@"
   cd "$old_dir"
else
   old_dir=`pwd`
   cd "$new_dir"
#  Since we changed directories, args which failed to glob
#  in old_dir may glob successfully in new_dir
   extant_files "$@"
#   Ahh, but this didn't work; not sure why. Worry about it later
#   new_args=`echo $@`
#   if [ $DEEBUG = "on" ]
#   then
#      pwd
#      echo "new_args $new_args"
#   fi
#   extant_files "$new_args"
   cd "$old_dir"
fi

if [ "$separate_suffix" = "yes" ]   
then                                
   the_suffix=" $the_suffix"        
fi                                  
if [ "$separate_prefix" = "yes" ]   
then                                
   the_prefix="$the_prefix "        
fi                                  

if [ $DEEBUG = "on" ]          
then                           
   echo "Before pre-suffixing: $extant_files_result"   
fi                             

if [ "$the_prefix" = "" -a "$the_suffix" = "" ]
then
   # No prefix--just echo filtered args
   echo $extant_files_result
elif [ "$the_suffix" = "" ]
then
   # prefix each filtered arg
   prefixed_result=
   for file in $extant_files_result
   do
      prefixed_result="$prefixed_result ${the_prefix}$file"
   done
   echo $prefixed_result
elif [ "$the_prefix" = "" ]
then
   # suffix each filtered arg
   prefixed_result=
   for file in $extant_files_result
   do
      prefixed_result="$prefixed_result $file${the_suffix}"
   done
   echo $prefixed_result
else
   # prefix and suffix each filtered arg
   prefixed_result=
   for file in $extant_files_result
   do
      prefixed_result="$prefixed_result ${the_prefix}$file${the_suffix}"
   done
   echo $prefixed_result
fi
exit
# $Log$
# Revision 1.3  2001/11/07 00:20:40  pwagner
# New -lib -dir and -suffix options
#
# Revision 1.2  2001/08/16 18:14:35  pwagner
# Added -prefix=xxx stuff
#
# Revision 1.1  2001/08/13 16:42:10  pwagner
# First commit under this name; previous name was foolish
#
