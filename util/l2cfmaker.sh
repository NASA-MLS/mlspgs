#!/bin/sh
# "$Id$"
# --------------- l2cfmaker.sh help
# l2cfmaker.sh    preprocessor acting on a section template of the l2cf
# to be processed later by the m4 macro processor
# Writes its output to stdout, error msgs to stderr
#
# Typical usage:
# l2cfmaker.sh [options] template_name [lists_of_args] > section_name.m4
# (this will be made automatic by the suffix rules in our Makefile)
# where template_name is a file laid out like the following:
# ;; (Beginning of template file)
# ;; --------------- args  (each line should be a comma-delimited list)
# a1,a2,a3,...
# b1,b2,b3,...
# ;; --------------- End args
# ;; --------------- foreword
# ;; this is a normal comment, probably describing the section
# ;; The preamble appears just once, at beginning of section
#    One or more (optional) m4 macro definition(s) or invocation(s)
#    One or more (optional) l2cf line(s)
# ;; --------------- End foreword
# ;; --------------- Repeat
# ;; The Repeat gets repeated, once for each arg in list
#    One or more (optional) m4 macro invocation(s)
#    One or more (optional) l2cf line(s)
# ;; --------------- End Repeat
# ;; --------------- afterword
# ;; The afterword appears just once, at end of section
#    One or more (optional) m4 macro definition(s) or invocation(s)
#    One or more (optional) l2cf line(s)
# ;; --------------- End afterword
# ;; (End of template file)
#
# Command-line options control:
# (1) Whether to use the lists of args embedded in the template file,
#     other lists found similarly embedded in a different file, or to use
#     lists supplied as arguments on the command line
# a list_of_args should look like 'a1,a2,a3,a4,..,an'
# Similarly, if extra_args is present, it might be 'b1,b2,b3,..,bn'
# and, if needed, followed by 'c1,c2,..,cn', .., 'z1,z2,..,zn'
# Then in the jth repetition of the Repeat portion the 
#      occurrence of     gets replaced by
#           #1                   a[j]
#           #2                   b[j]
#           ..                  ...
#           #k                   k[j]
#
# In contrast, any occurrence of #1 in the foreword or afterword
# portions gets replaced by the entire list of args 'a1,a2,a3,a4,..,an'
# (and similarly for #2 with the 'b' list, and so on)
# Would an example help? Maybe not, but here's one anyway
# 1% cat myTemplate
# ;; --------------- args  (each line should be a comma-delimited list)
# bird,cloud,snowflake
# plane,ufo,bridge
# superman,Zarathustra,J.L. Seagull
# ;; --------------- End args
# ;; --------------- foreword
# !define(nietzsche,{lookupinthesky: itsa_$1, itsa_$2, no_its_$3})
#     Enter: phone_booth
# ;; --------------- End foreword
# ;; --------------- Repeat
#   !nietzsche(#1,#2,#3)
# ;; --------------- End Repeat
# ;; --------------- afterword
#     Apotheosize: #1
# ;; --------------- End afterword
# 2% l2cfmaker.sh myTemplate 
#     Enter: phone_booth
#   lookupinthesky: itsa_bird, itsa_plane, no_its_superman
#   lookupinthesky: itsa_cloud, itsa_ufo, no_its_Zarathustra
#   lookupinthesky: itsa_snowflake, itsa_bridge, no_its_J.L. Seagull
#     Apotheosize: superman,Zarathustra,J.L. Seagull
#
#      Options
# -h           Print this help message and quit
# -A           Use the lists of args supplied as command-line arguments
# -a file_name Use the lists of args embedded within file_name
# -t file_name Use file_name for a template file
#
# Notes:
# If no options are supplied, then lists of args and template file are
# both assumed to be the first command-line argument; i.e.
# (command file_name) = (command -t file_name -a file_name)
# The hyphenated options must precede all other command line arguments
# If lists of args are supplied as command-line arguments, they must
# be the final command-line arguments
# You can't use both -A and -a simultaneously
#
# --------------- End l2cfmaker.sh help
# Notes:
# Requires perl
# An empty entry within a ist_of_args is denoted by '...,,..'
# e.g., in 'a1,a2,,a4' the 3rd element is ''
# The first list_of_args must not have any empty spaces
# i.e., no occurrences of ',,'

# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

#---------------------------- extract_loa
#
# Function to return the args embedded in the file named as $1
# as '"'-surrounded space-delimited things
# E.g. extract_loa myFile
# where myFile is a succession of <cr>-delimited lines
# bookended by lines with "--- args" and "--- End args"
# returns "line_1" "line_2" .. "line_n"

extract_loa()
{

   # Do we have enough args? No error if not, just blank output
      if [ $# -gt 0 ]
      then
         sed -n '/-- args/,/-- End args/ p' "$1" \
          | sed '1 d; $ d' | \
          perl -e "while (<STDIN>) {chop; print '\''; print; print '\' '}"
#          perl -e 'while (<STDIN>) {chop; print "\""; print; print "\" "}'
      fi
}
      
#---------------------------- nth_arg
#
# Function to return the nth (n > 0) arg of the std args after n
# eg given
# nth_arg 3 a b 'c or d' d e ..
# writes 'c or d' to standard output (w/o ' marks)

nth_arg()
{

   # Do we have enough args? No error if not, just blank output
     # echo "Entering nth_arg with args $@"
      the_arg=
      if [ $# -gt 1 ]
      then
         count_up=1
         n=$1
         shift
         while [ $count_up -lt $n ]
         do
           shift
           count_up=`expr $count_up + 1`
         done
         the_arg="$1"
      fi
     # echo "$the_arg"
}
      
#---------------------------- nth_list_element
#
# Function to return the nth (n > 0) element of a comma-delimited list
# where n is arg 1 and the list arg 2
# eg given
# nth_list_element 3 'a,b,c,d'
# writes 'c' to standard output (w/o ' marks)
# (Note: we have thus corrected for perl's c-centric array origin)

nth_list_element()
{

   # Do we have enough args?
      if [ $# -lt 2 ]
      then
         echo "Usage: nth_list_element n 'a,b,c,..'"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[1]"); print $parts[$ARGV[0]-1]' $1 "$2"
}
      
#---------------------------- how_many_list_elements
#
# Function to return the number of elements in a comma-delimited list
# where the list is the lone arg
# eg given
# how_many_list_elements 'a,b,c,d'
# writes '4' to standard output
# (Note: we have thus corrected for perl's c-centric array origin)

how_many_list_elements()
{

   # Do we have enough args?
      if [ $# -lt 1 ]
      then
         echo "Usage: how_many_list_elements 'a,b,c,..'"
         exit 1
      fi
      
      perl -e '@parts=split(",","$ARGV[0]"); print eval($#parts + 1)' "$1"
}
      
#------------------------------- Main Program ------------

#****************************************************************
#                                                               *
#                  * * * Main Program  * * *                    *
#                                                               *
#                                                               *
#	The entry point where control is given to the script         *
#****************************************************************

DEEBUG="false"
SPLITLONGLINES="true"
LOCALEPROBLEMS="true"
if [ "$LOCALEPROBLEMS" = "true" ]
then
   LC_ALL="en_US"
   export LC_ALL
   LOCALE="en_US"
   export LOCALE
fi
#
# Get arguments from command line
#

me="$0"
me_simple=l2cfmaker.sh
MYPATH="`echo $0 | sed 's/'$me_simple'//'`"
PERLUTIL="`echo $0 | sed 's/'$me_simple'/f90makedep.pl/'`"
REECHO="`echo $0 | sed 's/'$me_simple'/reecho.sh/'`"
templateFile=undefined
loaFile=undefined
loa_on_command_line="no"
loa_in_template_file="yes"

# Check that there are enough args
if [ $# -lt 1 ] ; then
   echo "Sorry--not enough args" 1>&2
   sed -n '/'$me_simple' help/,/End '$me_simple' help/ p' $me \
  | sed -n 's/^..//p' | sed '1 d; $ d'  1>&2
   exit                                                            
fi

more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

   case "$1" in                                                     

   -A )                                                     
      loa_on_command_line="yes"
      loa_in_template_file="no"
      shift
   ;;                                                               
   -a )                                                     
      loa_in_template_file="no"
      shift
      loaFile="$1"
      shift
   ;;                                                               
   -h | -help )                                                     
      sed -n '/'$me_simple' help/,/End '$me_simple' help/ p' $me \
     | sed -n 's/^..//p' | sed '1 d; $ d'                           
      exit                                                          
   ;;                                                               
   -t )                                                     
      shift
      templateFile="$1"
      shift
   ;;                                                               
   * )                                                              
      more_opts="no"                                             
   ;;                                                               

   esac                                                             
done
# Check whether we have template file yet
if [ "$templateFile" = "undefined" ] ; then
   templateFile="$1"
	shift
fi
# Check whether we have lists of args file yet
if [ "$loaFile" = "undefined" -a "$loa_in_template_file" = "yes" ] ; then
   loaFile="$templateFile"
fi
# Check that template file exists and that you may read it
if [ ! -r "$templateFile" ] ; then
   echo "Sorry--you cannot read template file $templateFile" 1>&2
	exit
# Do we have write permission in ./?                         
elif [ ! -w "./" ]                                           
   then                                                      
      echo "Sorry--need write permission in ./ to operate"    1>&2
      exit 1                                                 
fi

number_of_args=$#
how_many=`how_many_list_elements $1`
# if we are using lists of args embedded within a file,
if [ "$loa_on_command_line" = "no" ] ; then
# Check that lists of args file exists and that you may read it
   if [ ! -r "$loaFile" ] ; then
      echo "Sorry--you cannot read file $loaFile for lists of args" 1>&2
      exit
   fi
   if [ "$DEEBUG" = "true" ]
   then
      echo "About to extract list of args from $loaFile"
      sed -n '/-- args/,/-- End args/ p' "$loaFile" \
             | sed '1 d; $ d'
      sed -n '/-- args/,/-- End args/ p' "$loaFile" \
             | sed '1 d; $ d' | \
             perl -e 'while (<STDIN>) {chop; print "\""; print; print "\" "}'
      echo `extract_loa $loaFile`
      echo set `extract_loa $loaFile`
   fi
   # Unfortunately, the following "set" command fails
   # set `extract_loa $loaFile`
   # so instead we are forced to a more roundabout way
   arg_list=`extract_loa $loaFile`
   how_many=`how_many_list_elements $arg_list`
   # 1st check that we have any args:
   if [ "$arg_list" = "" ]
   then
     # No args--just print file
     cat "$templateFile"
     exit
   else
     # Need to know how many args we have
     number_of_args=`sed -n '/-- args/,/-- End args/ p' "$loaFile" \
         | sed '1 d; $ d'  | wc -l`
   fi
fi

# write foreword portion
rm -f temp.fwd 
sed -n '/-- foreword/,/-- End foreword/ p' "$templateFile" \
  | sed '1 d; $ d' > temp.fwd                     

arg_number=1
if [ "$DEEBUG" = "true" ]
then
   echo "Entering foreword while loop"
   echo "Number of args is $number_of_args"
fi
while [ $arg_number -le $number_of_args ]
do
   nth_arg $arg_number "$@"
   arg_list="$the_arg"
   if [ "$loa_on_command_line" = "no" ] ; then                   
     next_arg=`expr $arg_number + 1`                             
     arg_list=`sed -n '/-- args/,/-- End args/ p' "$loaFile" \
       | sed -n "$next_arg p"`                                   
   fi                                                            
   mv temp.fwd tempp.fwd
   sed 's/#'$arg_number'/'"$arg_list"'/g' tempp.fwd > temp.fwd
   arg_number=`expr $arg_number + 1`
done
if [ "$DEEBUG" = "true" ]   
then                        
   echo "Done with foreword portion"
   cat temp.fwd
fi

# write Repeat portions
rep=1
if [ "$DEEBUG" = "true" ]
then
   echo "Entering outer Repeat while loop"
fi
while [ $rep -le $how_many ]
do
   sed -n '/-- Repeat/,/-- End Repeat/ p' "$templateFile" \
     | sed '1 d; $ d' >> temp.fwd
   arg_number=1
   if [ "$DEEBUG" = "true" ]
   then
      echo "Entering inner Repeat while loop with $rep repetition"
   fi
   while [ $arg_number -le $number_of_args ]
   do
      nth_arg $arg_number "$@"
      arg_list="$the_arg"
      if [ "$loa_on_command_line" = "no" ] ; then
        next_arg=`expr $arg_number + 1`
        arg_list=`sed -n '/-- args/,/-- End args/ p' "$loaFile" \
          | sed -n "$next_arg p"`
      fi
      element=`nth_list_element $rep "$arg_list"`
      if [ "$DEEBUG" = "true" ]
      then
         echo "element $element of $arg_list"
      fi
      mv temp.fwd tempp.fwd
      sed 's/#'$arg_number'/'"$element"'/g' tempp.fwd > temp.fwd
      arg_number=`expr $arg_number + 1`
   done
   rep=`expr $rep + 1`
done

if [ "$DEEBUG" = "true" ]   
then                        
   echo "Done with Repeat portion"
   cat temp.fwd
fi

# write afterword portion
sed -n '/-- afterword/,/-- End afterword/ p' "$templateFile" \
  | sed '1 d; $ d' >> temp.fwd                     

arg_number=1
if [ "$DEEBUG" = "true" ]
then
   echo "Entering afterword while loop"
fi
while [ $arg_number -le $number_of_args ]
do
   nth_arg $arg_number "$@"   
   arg_list="$the_arg"        
   if [ "$loa_on_command_line" = "no" ] ; then                   
     next_arg=`expr $arg_number + 1`                             
     arg_list=`sed -n '/-- args/,/-- End args/ p' "$loaFile" \
       | sed -n "$next_arg p"`                                   
   fi                                                            
   mv temp.fwd tempp.fwd
   sed 's/#'$arg_number'/'"$arg_list"'/g' tempp.fwd > temp.fwd
   arg_number=`expr $arg_number + 1`
done
if [ "$DEEBUG" = "true" ]   
then                        
   echo "Done with afterword portion"
fi
if [ "$SPLITLONGLINES" = "true" ]
then
   $PERLUTIL -f temp.fwd -d "," -t '$' -c ';'
else
   cat temp.fwd
fi
rm -f temp.fwd tempp.fwd

# Done at last
exit
# $Log$
