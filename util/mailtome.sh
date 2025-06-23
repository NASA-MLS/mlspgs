#!/bin/sh
# --------------- mailtome help
# mail files to list of recipients
# Usage:
# mailtome.sh [options] file-list
#    O p t i o n s
# -h[elp]     print brief help message; exit
# -dryrun     don't actually mail the files--just echo commands
# -F          set mail Return-path to (hard wired) pwagner@mail.jpl.nasa.gov
# -f address  show mail as having come from address
# -r dir      keep record of mailing in directory dir
# -s "subject line"  use "subject line"
# -R "recipients list"  use "recipients list"
# -Rf file    change RECIPIENTS to list of addresses found in file
#
# Note:
# (1) you must be allowed to mail from this machine
# (2) files must exist
# (3) if file-list not supplied, it will mail stdin
# 
# Result:
# files mailed; optionally a log of the transfer created
# --------------- End mailtome help
# Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

# "$Id$"

#---------------------------- read_file_into_array
#
# read each line of stdin
# catenating them into an array which we will return
# Ignore any entries beginning with '#' character
# In fact, only first entry in each line is kept
# Possible improvements:
#   Other comment signifiers
#   Choose field number other than 1

read_file_into_array()
{
  array_result=''
  while read line; do
    element=`echo $line | awk '$1 !~ /^#/ {print $1}'`
    if [ "$element" != "" ]
    then
      array_result="$array_result $element"
    fi
  done
  echo $array_result
}
      
#
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
    #  echo $our_host_name
   # if in form host.moon.planet.star.. extract host
      our_host_name=`echo $our_host_name | sed 's/\./,/g'`
      our_host_name=`perl -e '@parts=split(",","$ARGV[0]"); print $parts[0]' $our_host_name`
      echo $temp${pt}$our_host_name${pt}$$
}
      
# ************
# Main Program
# ************
# 
# First define where we're going to do all our work
TODAY=`date +"%Y-%m-%d"`
SOMETHING=`get_unique_name $TODAY`
STDINSUB=`get_unique_name $TODAY`
TEMPLATE=`get_unique_name sndmltmp`

# Where can we write the MAILBUNDLE? Current directory? HOME?
if [ -w "`pwd`" ]
then
  record_dir=`pwd`
elif [ -w "$HOME" ]
then
  record_dir="$HOME"
fi
debug="no"
forgereturn="no"
forge=""
me="$0"
my_name=mailtome
I=mailtome
dryrun="no"
record="no"
subject="$SOMETHING"
RECIPIENTS="pwagner"
more_opts="yes"
/bin/rm -f ~/$SOMETHING
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    shift
       dryrun="yes"
       ;;
    -F )
	    shift
       forgereturn="yes"
       ;;
    -f )
	    shift
       forge="$1"
       shift
       ;;
    -r )
	    shift
       record="yes"
       record_dir="$1"
       shift
       ;;
    -R )
	    shift
       RECIPIENTS="$1"
       shift
       ;;
    -Rf )
	    shift
       RECIPIENTS=`cat $1 | uniq | read_file_into_array`
       shift
       ;;
    -s )
	    shift
       subject="$1"
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

files="$@"
SOMETHING="$record_dir/$SOMETHING"
STDINSUB="$record_dir/$STDINSUB"
TEMPLATE="$record_dir/$TEMPLATE"

if [ "$files" = "" ]
then
  cat > "$STDINSUB"
  files="$STDINSUB"
fi

#Do we need to forge headers?
if [ "$forge" != "" ]
then
  cat > "$TEMPLATE" << EOF
To: $RECIPIENTS
From: $forge
Subject: $subject
EOF
fi

SENDMAILOPTS="-t"
if [ "$forgereturn" = "yes" ]
then
  SENDMAILOPTS="-t -fpwagner@mail"
fi

if [ $dryrun = "no" ]
then
  #echo "About to email $files to $RECIPIENTS"
  if [ "$forge" != "" ]
  then
    cat "$TEMPLATE" $files | /usr/sbin/sendmail $SENDMAILOPTS
    # cat "$TEMPLATE" $files | /usr/sbin/sendmail -t
  else
    cat $files | mail -s"$subject" $RECIPIENTS
  fi
elif [ $dryrun = "yes" ]
then
  if [ "$forge" != "" ]
  then
    cat "$TEMPLATE"
    echo 'cat '"$TEMPLATE"' '$files' | /usr/sbin/sendmail ' $SENDMAILOPTS
    # echo 'cat '"$TEMPLATE"' '$files' | /usr/sbin/sendmail -t'
  else
    echo 'cat '"$files"' | mail -s"'$subject'"' $RECIPIENTS
  fi
fi
if [ $record = "no" ]
then
  rm -f "$TEMPLATE"
  if [ "$files" = "$STDINSUB" ]
  then
    rm "$files"
  fi
  exit
fi
# ------------------------------------------summary of the email
if [ "$forge" != "" ]
then
  cat "$TEMPLATE" > "$SOMETHING"
else
  echo "mail -s$subject $RECIPIENTS" > "$SOMETHING"
fi
if [ "$files" != "$STDINSUB" ]
then
  echo "files: $files" >> "$SOMETHING"
fi
echo '<<contents>>' >> "$SOMETHING"
cat $files >> "$SOMETHING"
if [ "$files" = "$STDINSUB" ]
then
  rm "$files"
fi
rm -f "$TEMPLATE"
exit
# $Log$
# Revision 1.2  2005/10/29 00:22:17  pwagner
# May set Return-path to pwagner@mail.jpl.nasa.gov
#
# Revision 1.1  2005/10/01 00:36:28  pwagner
# First commit
#
