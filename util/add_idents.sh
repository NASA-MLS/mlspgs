#!/bin/sh
# add_idents.sh

# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"
# --------------- add_idents.sh help
#Use(1) (default)
#add idents to list of files having $id and $RCS lines
#List compiled from args on command-line
#or else automatically computed to match pattern
#  *{$suffix} in supplied directory names
#The fix involves adding the following two blocs to each file:
#       (bloc 1)
# !---------------------------- RCS Ident Info -------------------------------
#   character (len=*), private, parameter :: ModuleName= "$RCSfile$"
#   private :: not_used_here
# !---------------------------------------------------------------------------
#   . . .
#       (bloc 2)
# !--------------------------- end bloc --------------------------------------
#   logical function not_used_here()
#   character (len=*), parameter :: IdParm = &
#        "$Id$"
#   character (len=len(idParm)) :: Id = idParm
#     not_used_here = (id(1:1) == ModuleName(1:1))
#     print *, Id ! .mod files sometimes change if PRINT is added
#   end function not_used_here
# !---------------------------------------------------------------------------
#end module L1BData
#(where the added lines have been marked with the ">")
# Note that we assume the older style of rcs info is already present
# In the present incarnation, 
# we also replace the older wording of the copyright statement
#
#Use(2)
#Print toc or api from list of files
#List compiled from args on command-line
#or else automatically computed to match pattern
#  *{$suffix} in supplied directory names
#Usage:
#add_idents.sh [opt1 [arg1]] [opt2 [arg2]] .. [optm [argm]] [filelist]
#
#    O p t i o n s
# -h[elp]       print brief help message; exit
# -d path       path of directory if automatic
# -suf suffix   file name suffix if automatic; e.g. ".f90"
#       with use (1) only
# -dryrun       print diff between result and original only, don't replace
# -dtxt path    where to find RCSIdent.txt, RCSModule.txt, and COPYRIGHT.txt
# -f n          forcibly insert bloc n 
#                 where  n       bloc               where
#                        1    COPYRIGHT.txt       beginning of file
#                        2    RCSModule.txt       before 1st "contains"
#                        3    endbloc.txt        before "end module"
# -rep script   use script instead of replacetext.sh
#       with use (2) only
# -api          printing api
# -toc          printing toc
# -id           printing $id line, too
# -rcs          printing $rcs line, too
# -n            print only if there is a (api,toc) section
#Notes:
#(1)If option -suf present but -d absent, will search current working directory
#(2)The options -suf and -d may be repeated; e.g.
#   add_idents -d dir1 -d dir2 .. -suf sufx1 -suf sufx2 ..
#(3)The default use is Use(1); -api and -toc are mutually exclusive
#(4)use(1) requires 
#   (a) replacetext.sh to be in the path (unless overridden by -rep)
#   (b) endbloc.txt, RCSModule.txt, and COPYRIGHT.txt
#       to be in the current working directory (unless overridden by -dtxt)
#(5)Not tested with c, f77, shell scripts, Makefiles, or include files
# --------------- End add_idents.sh help
#Another example of my regrettable tendency to bundle a number of
#different functionalities within a single multi-use file
#Isn't it cleaner to devote a single file to a single use?

#---------------------------- get_unique_name
#
# Function returns a unique name based on arg, PID and HOSTNAME
# e.g.,
#           temp_file_name=`get_unique_name foo`
#           echo $temp_file_name
# might print foo.colossus.21455
# if no arg, defaults to "temp" (very original name)
# if two args present, assumes second is punctuation to
# use in pace of "."

get_unique_name()
{

   # How many args?
      if [ $# -gt 1 ]
      then
        pt="$2"
        temp="$1"
      elif [ $# -gt 0 ]
      then
        pt="."
        temp="$1"
      else
        pt="."
        temp="temp"
      fi
   # Is $HOST defined?
      if [ "$HOST" != "" ]
      then
         our_host_name="$HOST"
      elif [ "$HOSTNAME" != "" ]
      then
         our_host_name="$HOSTNAME"
      else
         our_host_name="host"
      fi
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
}
      
#------------------------------- mv_if_diff ------------
#
# replace 2nd arg with 1st if diff
# unless dryrun is "yes", in which case just print diff
# usage: mv_if_diff file1 file2

mv_if_diff()
{
  test=`diff "$1" "$2"`
  if [ "$1" = "$2" ]
  then
    echo "$1 and $2 are the same files"
  elif [ "$test" = "" ]
  then
    echo "$1 and $2 are the same"
    rm "$1"
  elif [ "$dryrun" = "yes" ]
  then
    echo "$1 and $2 diff:"
    diff "$1" "$2"
    rm "$1"
  else
    mv "$1" "$2"
  fi
}

#------------------------------- add_the_lines ------------
#
# Function to add the lines needed to assure that after compilation
# a file's $id and $RCS idents will be found in the executable
# Relies upon a trick to circumvent Intel's zealous optimizer
# RCS id no longer split into two blocs; who cares about NAG's propensity
# to launch cascade of recompilation when date field in ID changes
# at each cvs commit?
# Upgraded to:
# (1) Replace the copyright statement with its new wording
# (2) Rejoin $id and $RCS identifiers,
# (3) Replace the not_used_here private function code bloc with new end bloc
# usage: add_the_lines file

add_the_lines()
{
file="$1"
tempfile=`get_unique_name addlines`
# has the file already had the new copyright statement added?
test=`grep -i 'Office of Technology Transfer' $file`                     
if [ "$test" = "" ]                                      
then                                                     
  # First: replace the old Copyright
  test1=`echo $force_blocs | grep 1`
  $REPLACER -i "Copyright" "Government Sponsorship" \
    "$file" $dtxt/COPYRIGHT.txt > "$tempfile"
  mv_if_diff "$tempfile" "$file"
  testC=`grep -i copyright "$file"`
  if [ "$test1" != "" -a "$testC" = "" ]
  then
    cat $dtxt/COPYRIGHT.txt "$file" > "$tempfile"
    mv_if_diff "$tempfile" "$file"
  fi
fi
# Second: replace the old RCS bloc
test2=`echo $force_blocs | grep 2`
$REPLACER -i "-- RCS Module Info --" "!--------------------------------------------" \
  $file $dtxt/RCSModule.txt > "$tempfile"
mv_if_diff "$tempfile" "$file"
testC=`grep -i "rcs module" "$file"`
if [ "$test2" != "" -a "$testC" = "" ]
then
  $REPLACER -i -b "contains" "contains" \
  $file $dtxt/RCSModule.txt > "$tempfile"
  mv_if_diff "$tempfile" "$file"
fi

# Third: replace the old not_used_here bloc
test3=`echo $force_blocs | grep 3`
test=`grep -ie '-- end bloc --' $file`                     
if [ "$test" != "" ]                                      
then                                                     
  # echo Must replace the old end bloc
$REPLACER -i "-- end bloc --" "!--------------------------------------------" \
    $file $dtxt/endbloc.txt > "$tempfile"
  mv_if_diff "$tempfile" "$file"
else
  $REPLACER -i "logical function not_used_here" "end function not_used_here" \
    $file $dtxt/endbloc.txt > "$tempfile"
  mv_if_diff "$tempfile" "$file"
  testC=`grep -i "end function not_used_here" "$file"`
  if [ "$test3" != "" -a "$testC" = "" ]
  then
    $REPLACER -i -b "end module" "end module" \
    $file $dtxt/endbloc.txt > "$tempfile"
    mv_if_diff "$tempfile" "$file"
  fi
fi
}

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

#------------------------------- print_the_lines ------------
#
# Function to print the lines from a file bookended by e.g. where section="toc"
#'=== (start of toc) ==='
# and
#'=== (end of toc) ==='
# 
# usage: print_the_lines section file

print_the_lines()
{
section="$1"
file="$2"
#echo "Trying to print lines in $file, section $section"
#echo '/=== (start of '$section') ===/,/=== (end of '$section') ===/ p' $file
sed -n \
  '/=== (start of '$section') ===/,/=== (end of '$section') ===/ p' $file \
           | sed -n 's/^.//p' | sed '1 d; $ d'
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
my_name=add_idents.sh
# The following takes possible values "add", "api", "toc"
# the default is "add"
my_use="add"
dryrun="no"
print_id="no"
print_rcs="no"
skip_if_none="no"
dtxt=`pwd`
REPLACER="replacetext.sh"
force_blocs=""
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
	-dryrun )
	    dryrun="yes"
	    shift
	;;
	-dtxt )
	    dtxt="$2"
	    shift
       shift
	;;
	-rep )
	    REPLACER="$2"
	    shift
       shift
	;;
	-f )
	    force_blocs="$2 $force_blocs"
	    shift
       shift
	;;
	-api )
	    my_use="api"
	    shift
	;;
	-toc )
	    my_use="toc"
	    shift
	;;
	-n )
	    skip_if_none="yes"
	    shift
	;;
	-id )
	    print_id="yes"
	    shift
	;;
	-rcs )
	    print_rcs="yes"
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
for the_file in $arglist
do
  if [ "$my_use" = "add" ]
  then
    echo Adding idents to $the_file > temp
  else
    echo '******** ' "$my_use for $the_file" ' ********' > temp
  fi
  if [ "$print_id" = "yes" ]
  then
    grep -i '\$id' $the_file >> temp
  fi
  if [ "$print_rcs" = "yes" ]
  then
    grep -i '\$rcs' $the_file >> temp
  fi
  case "$my_use" in
    add)
      cat temp
      add_the_lines $the_file
      ;;
    api | toc)
      print_the_lines "$my_use" $the_file > temp1
      how_many=`cat temp1 | wc -l`
      if [ "$how_many" -gt 0 -o "$skip_if_none" = "no" ]
      then
        echo ""
        cat temp temp1
      fi
      rm -f temp1
      ;;
    *)
      echo "Programming error in add_idents.sh"
      ;;
  esac
  rm -f temp
done
exit 0
# $Log$
# Revision 1.5  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.4  2005/06/22 22:43:21  pwagner
# Changes to add new Copyright statement, rcs blocs
#
# Revision 1.3  2002/10/11 23:01:05  pwagner
# Added -n option with use(2)
#
# Revision 1.2  2002/10/11 19:37:42  pwagner
# Added use(2): print toc/api sections from source files
#
# Revision 1.1  2002/10/08 16:57:33  pwagner
# First commit
#
