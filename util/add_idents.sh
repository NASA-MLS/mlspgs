#!/bin/sh
# add_idents.sh

# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
# --------------- add_idents.sh help
#add idents to list of files having $id and $RCS lines
#List compiled from args on command-line
#or else automatically computed to match pattern
#  *{$suffix} in supplied directory names
#The fix involves adding the following lines to each file:
#    "$RCSfile$"
#  >  private :: not_used_here
#   . . .
#  >  logical function not_used_here()
#  >    not_used_here = (id(1:1) == ModuleName(1:1))
#  >  end function not_used_here
#end module L1BData
#(where the added lines have been marked with the ">")
#
#Usage:
#add_idents.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#
#    O p t i o n s
# -d path       path of directory if automatic
# -suf suffix   file name suffix if automatic; e.g. ".f90"
# -h[elp]       print brief help message; exit
#Notes:
#(1)If option -suf present but -d absent, will search current working directory
#(2)The options -suf and -d may be repeated; e.g.
#   add_idents -d dir1 -d dir2 .. -suf sufx1 -suf sufx2 ..
# --------------- End add_idents.sh help
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

#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************
id_list=
RCS_list=
echo "Who or what is $@"
suffixes=""
directories=""
me="$0"
my_name=add_idents.sh
#
# Get arguments from command line
#

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do
    case "$1" in
	-d )
	    directories="$directories $2"
	    shift
       shift
	;;
	-suf )
	    suffixes="$suffixes $2"
	    shift
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
extant_files "$@"
arglist="$extant_files_result"
if [ "$suffixes" != "" ]
then
  # Assemble automatic list of files
  if [ "$directories" = "" ]
  then
    directories="."
  fi
  args_dir=`pwd`
  for dir in $directories
  do
    cd $dir
    for suffix in $suffixes
    do
      extant_files *"$suffix"
      if [ "$extant_files_result" != "" ]
      then
        for file in $extant_files_result
        do
          arglist="$arglist $dir/$file"
        done
      fi
    done
    cd $args_dir
  done
fi
rm -f temp
for file in $arglist
do
  echo file is $file                   
  test=`grep -i 'not_used_here' $file` 
  if [ "$test" = "" ]
  then
# > >     test1=`grep -i 'END MODULE' $file`
# > >     if [ "$test1" != "" ]
# > >     then
# > >       sed 's/END MODULE/end module/' $file > temp
# > >       mv temp $file
# > >     fi
    test1=`grep -i '^ *contains' $file`
    if [ "$test1" = "" ]
    then
      sed '/\$RCSfile/ a\
  private \:\: not_used_here ' $file | \
      sed '/[eE][nN][dD] [mM][oO][dD][uU][lL][eE]/ i\
contains \
  logical function not_used_here()\
    not_used_here = (id(1:1) == ModuleName(1:1))\
  end function not_used_here\
'     > temp
    else
      sed '/\$RCSfile/ a\
  private \:\: not_used_here ' $file | \
      sed '/[eE][nN][dD] [mM][oO][dD][uU][lL][eE]/ i\
  logical function not_used_here()\
    not_used_here = (id(1:1) == ModuleName(1:1))\
  end function not_used_here\
'     > temp
    fi
    mv temp $file
  else
    echo "$file already has ident added in this way"
  fi
done
exit 0
# $Log$
