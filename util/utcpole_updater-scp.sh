#!/bin/sh
# utcpole_updater-scp.sh
# Plagiarised Raytheon's update_utcpole.sh
# distributed with the Toolkit
# but I've removed all the email and user-input stuff
# Should be useable by a cron job
# or you can use it directly
#
# REQUIRES:
# PGSHOME must have been defined (by e.g., pgs-env.csh)
# write access to $PGSHOME/database/common/CSC
 
# Copyright 2023, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

# Check that PGSHOME is defined
if [ "$PGSHOME" = "" ]
then
  echo "PGSHOME is not yet defined"
  echo "You probably need to run one of the scripts"
  echo "pgs-env.csh or pgs-env.ksh"
  exit 1
fi

# Check that you have write access permission to the file
if [ ! -w $PGSHOME/database/common/CSC/utcpole.dat ]
then
  echo "You lack write access permission to the file utcpole.dat"
  echo "Likely causes are"
  echo "(1) It is a directory belonging to someone else"
  echo "(2) The name or location of the directory has been changed"
  echo "(3) A prior attempt to update the file badly misfired"
  exit 1
fi

# change directory to the utcpole.dat directory:

cd $PGSHOME/database/common/CSC
# if original file is missing abort with a message
if [ ! -s utcpole.dat ] ; then


   exit 20
fi

# get file of recent and predicted data from Naval Observatory 
#echo "Are you a DAAC-location [yes/no]"
#read user_response
user_response="No"
case  "$user_response" in
    y* | Y* )
echo "Please enter your DAAC-specific ftp machine:"
read ftp_machine
echo "Please enter your user name:"
read ftp_user
ftp -n <<GG>holder
#open p0fwi09
#user cmts2
open $ftp_machine
user $ftp_user
quote site maia.usno.navy.mil 
user anonymous EOSuser
prompt
cd ser7
dir
get finals.data
bye
GG
    ;;
    * )
#scp pwagner@fox.jpl.nasa.gov:/science/pge/v1300-nrt/toolkit5.2.18/database/common/CSC/utcpole.dat ./
scp pwagner@kestrel.jpl.nasa.gov:/science/pge/v0502-l2/toolkit5.2.18/database/common/CSC/utcpole.dat ./
    ;;
    esac

#abort and return status if ftp failed

state=$?

if [ $state != 0 ];  then 

exit $state
fi

exit 0
# $Log$
