#!/bin/sh
# edit_l2cfs.sh

# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
# --------------- edit_l2cfs.sh help
#Use(1) (the default)
#edit l2cf files to comment-out any retrieve sections
# i.e., any sections bounded by
# begin Retrieve       ;;; begin Retrieve
#    .  .  .        => ;;;    .  .  .
# end Retrieve         ;;; end Retrieve
# 
#Use(2) (the reverse)
#uncomment any retrieve sections
#Usage:
#edit_l2cfs.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#Result: each file replaced with its edited copy (unless -print option)
#
#    O p t i o n s
# -1            use(1)
# -2            use(2)
# -s mySection  instead of Retrieve, comment out section mySection
# -print        instead of replacing each file, print result to stdout
# -diff         instead of replacing each file, print diff between
#                  result and the original to stdout
# -h[elp]       print brief help message; exit
#Notes:
#(1)Options -1, -2
#(2)The default use is Use(1)
# --------------- End edit_l2cfs.sh help
#Another example of my regrettable tendency to bundle a number of
#different functionalities within a single multi-use file
#Isn't it cleaner to devote a single file to a single use?

#------------------------------- print_the_lines ------------
#
# Function to print the lines from a file bookended by 
# e.g. where object="INPUTPOINTER"
#'OBJECT     = INPUTPOINTER'
# . . .
#'END_OBJECT             = INPUTPOINTER'
#
# Optionally end each line (except the last) with a "\" 
# so it can serve as replacement text in a sed script
# usage: print_the_lines object file [opt]

print_the_lines()
{
object="$1"
file="$2"
if [ "$usage" = "1" ]
then
  sed '/begin *'$object'/,/end *'$object'/ s/^./;;; &/' $file
  #sed -n \
  #'/begin *'$object'/,/end *'$object'/ p' $file 
else
  sed '/;;; begin *'$object'/,/;;; end *'$object'/ s/^;;; //' $file
  #sed -n \
  #'/begin *'$object'/,/end *'$object'/ p' $file \
  #| sed 's/.*$/&\\/' | sed '$ s/\\$//'

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
# echo "All the args: $@"
suffixes=""
directories=""
me="$0"
my_name=edit_l2cfs.sh
I=edit_l2cfs
unique_name="`echo $0 | sed 's/'$I'/unique_name/'`"
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
resed="`echo $0 | sed 's/'$I'/resed/'`"

print_to_stdout="no"
more_opts="yes"
usage="1"
mySection="Retrieve"

#
# Get arguments from command line
#

while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -1 )
       usage="1"
       shift
       ;;
    -2 )
       usage="2"
       shift
       ;;
    -diff )
       print_to_stdout="diff"
       shift
       ;;
    -print )
       print_to_stdout="yes"
       shift
       ;;
    -s )
       mySection="$2"
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
diff_list=""
nodiff_list=""
rm -f temp1.sed
for the_file
do
  if [ "$print_to_stdout" = "yes" ]
  then
    print_the_lines "$mySection" "$the_file"
  elif [ "$print_to_stdout" = "diff" ]
  then
    print_the_lines "$mySection" "$the_file" > temp1.sed
    how_many=`diff "$the_file" temp1.sed | wc -l`
    if [ "$how_many" -gt 0 ]
    then
      diff_list="$diff_list $the_file"
    else
      nodiff_list="$nodiff_list $the_file"
    fi
    rm -f temp1.sed
  else
    print_the_lines "$mySection" "$the_file" > temp1.sed
    mv temp1.sed "$the_file"
  fi
done
if [ "$print_to_stdout" = "diff" ]
then
  if [ "$diff_list" = "" ]
  then
    echo "None of the files would be changed"
  elif [ "$nodiff_list" = "" ]
  then
    echo "All of the files pass would be changed"
  else
    echo "The following files would be changed"
    echo "$diff_list"
    echo "The following files would be left the same"
    echo "$nodiff_list"
  fi
fi
exit 0
# $Log$
