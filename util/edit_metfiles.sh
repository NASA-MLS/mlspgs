#!/bin/sh
# edit_metfiles.sh

# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
# --------------- edit_metfiles.sh help
#Use(1) (the default)
#edit .met files so that their INPUTPOINTER values
#match what we have decided level 1 hdf5 files should have
#The fix involves changing lines to each file to the following:
#  >  OBJECT                 = INPUTPOINTER
#  >    NUM_VAL              = 1
#  >    VALUE                = "Found at /PCF"
#  >  END_OBJECT             = INPUTPOINTER
#(where the changed lines have been marked with the ">")
#
#Use(2)
#As in use (1) except to
#match what we have decided level 2 hdfeos5 files should have
#The fix involves changing lines to each file to the following:
#  >  OBJECT                 = INPUTPOINTER
#  >    NUM_VAL              = 1
#  >    VALUE                = "Found at /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/PCF"
#  >  END_OBJECT             = INPUTPOINTER
#(where the changed lines have been marked with the ">")
#
#Use(3)
#More general than (1) and (2); given myOBJECT and myVALUE
#The fix involves changing one line in each file as the following illustrates:
#   OBJECT                 = myOBJECT
#    .  .  .
#  >    VALUE                = "myVALUE"
#   END_OBJECT             = myOBJECT
#(where only the changed line has been marked with the ">")
#Usage:
#edit_metfiles.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#
#    O p t i o n s
# -1            use(1)
# -2            use(2)
# -obj myOBJECT instead of INPUTPOINTER, change value of myOBJECT
# -val myVALUE  instead of use(1) or (2), change value of OBJECT to myVALUE
# -print        instead of creating new file(s), print results of sed to stdout
# -test         instead of creating new file(s), confirm that files already 
#                  match output of running script; i.e., diff would be null
# -h[elp]       print brief help message; exit
#Notes:
#(1)Options -1, -2, -val are mutually exclusive
#(2)If -obj is supplied, then -val should be, too
#(3)The default use is Use(1)
# --------------- End edit_metfiles.sh help
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
if [ "$3" = "" ]
then
  sed -n \
  '/OBJECT * = *'$object'/,/END_OBJECT * = *'$object'/ p' $file 
else
  sed -n \
  '/OBJECT * = *'$object'/,/END_OBJECT * = *'$object'/ p' $file \
  | sed 's/.*$/&\\/' | sed '$ s/\\$//'

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
my_name=edit_metfiles.sh
I=edit_metfiles
unique_name="`echo $0 | sed 's/'$I'/unique_name/'`"
the_splitter="`echo $0 | sed 's/'$I'/split_path/'`"
reecho="`echo $0 | sed 's/'$I'/reecho/'`"
resed="`echo $0 | sed 's/'$I'/resed/'`"

print_to_stdout="no"
more_opts="yes"
usage="1"
myOBJECT=""
myVALUE=""

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
    -obj )
       myOBJECT="$2"
       shift
       shift
       ;;
    -val )
       myVALUE="$2"
       shift
       shift
       ;;
    -print )
       print_to_stdout="yes"
       shift
       ;;
    -test )
       print_to_stdout="diff"
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
# Check for constraints on options
if [ "$myOBJECT" != "" ]
then
  if [ "$myVALUE" = "" ]
  then
    echo "You must supply -val myVALUE to go with -obj myOBJECT"
    exit 1
  fi
fi
if [ "$myVALUE" != "" ]
then
  usage="3"
fi
case "$usage" in
   1)
     myOBJECT="INPUTPOINTER"
     myVALUE="Found at /PCF"
     ;;
   2)
     myOBJECT="INPUTPOINTER"
     myVALUE="Found at /HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/PCF"
     ;;
   *)
     ;;
esac
# Prepare a file containing sed instructions to carry out
# the editing task
example_file="$1"
rm -f temp.sed
echo '/OBJECT * = *'$myOBJECT'/,/END_OBJECT * = *'$myOBJECT'/ c \' > temp.sed
print_the_lines "$myOBJECT" "$example_file" "yes" >> temp.sed
rm -f temp1.sed
sed '/VALUE/ s^".*"^"'"$myVALUE"'"^' temp.sed > temp1.sed
mv temp1.sed temp.sed
diff_list=""
nodiff_list=""
# Now use the file of sed instructions plus our resed script
# to act upon the args (which we hope are the filenames)
#if [ ! -x "$resed" ]
#then
#  chmod u+x "$resed"
#fi
#Hmmmm-skip using resed for now
#until problem of creating temp files in read-only directory is solved
#if [ "$print_to_stdout" = "yes" ]
#then
#  $resed -print -f temp.sed $@
#else
#  $resed -f temp.sed $@
#fi
for the_file
do
  if [ "$print_to_stdout" = "yes" ]
  then
    sed -f temp.sed "$the_file"
  elif [ "$print_to_stdout" = "diff" ]
  then
    sed -f temp.sed "$the_file" > temp1.sed
    how_many=`diff "$the_file" temp1.sed | wc -l`
    if [ "$how_many" -gt 0 ]
    then
      diff_list="$diff_list $the_file"
    else
      nodiff_list="$nodiff_list $the_file"
    fi
    rm -f temp1.sed
  else
    sed -f temp.sed "$the_file" > temp1.sed
    mv temp1.sed "$the_file"
  fi
done
rm temp.sed
if [ "$print_to_stdout" = "diff" ]
then
  if [ "$diff_list" = "" ]
  then
    echo "All the files pass the test"
  elif [ "$nodiff_list" = "" ]
  then
    echo "Sorry--none of the files pass the test"
  else
    echo "The following files failed the test"
    echo "$diff_list"
    echo "The following files passed the test"
    echo "$nodiff_list"
  fi
fi
exit 0
# $Log$
# Revision 1.3  2003/05/28 20:25:46  pwagner
# Fixed embarrassing error where usage(2) no different from (1)
#
# Revision 1.2  2003/05/28 16:20:34  pwagner
# Small fix--not yet tested with resed script
#
# Revision 1.1  2003/05/22 20:52:19  pwagner
# First commit
#
