#!/bin/sh
# utcpole_updater.sh
# Plagiarised Raytheon's update_utcpole.sh
# distributed with the Toolkit
# but I've removed all the email and user-input stuff
# Should be useable by a cron job
# or you can use it directly
#
# REQUIRES:
# PGSHOME must have been defined (by e.g., pgs-env.csh)
# write access to $PGSHOME/database/common/CSC
 
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
echo "Please wait ...."
ftp -n <<RRR>holder
open maia.usno.navy.mil
user anonymous EOSuser
prompt
cd ser7
dir
get finals.data
bye
RRR
    ;;
    esac

#abort and return status if ftp failed

state=$?

if [ $state != 0 ];  then 

exit $state
fi

#abort if file is missing or short (less than about a year's data)

if [ ! -s finals.data ] || [ `wc -c < finals.data` -lt  50000 ]
then
echo "Failure to get data file from maia.usno.navy.mil"
exit 22
 
fi
 
# parse the results of ftp for the date, to be used in the new header

grep "finals.data" holder | awk '{print $5 " " $6 " " $7 " " $8}'  > func_input

# overwrite all the data in the old file whose dates are duplicated in new data
# extend file length as needed to accommodate new data
# growth is a line a day. The "func_input" temp file has the header data

PGS_CSC_UT1_update < func_input

state=$?

echo "Status of PGS_CSC_UT1_update call was ($state)"

# test return value.  

    case $state in
      0)
 
#  return value OK but make one further check to be sure file 
#  reaches 30 days or more into the future

DAYTOTL=`date '+%m %Y' | awk '{ sss = 30*($1-1) + 365.2425*($2-1972) + 30;  print(sss) }'`
DAYTOTL=`echo ${DAYTOTL} | awk -F"." '{print $1}'`

   new_length=` cat utcpole.dat.NEW | wc -l `  ;
   grew=`expr $new_length - $DAYTOTL `  ;
   if [ $grew -ge 0 ] ; then

#      all is A-OK

       /bin/mv -f utcpole.dat.NEW utcpole.dat

   mv_state=$?
   echo "Status of MOVE command was ($mv_state)"

#BEGIN nested if - move file in if allowed - branch as error if mv fails

       if [ $mv_state -eq 0 ] ; then

       /bin/rm func_input holder
       /bin/rm finals.data
       exit 0

      else
#error branch within nested if - used if "mv" fails

state=21

#      clean up
       /bin/rm func_input holder
       /bin/rm finals.data
       /bin/rm utcpole.dat.NEW

        fi
#EXITED nested if

   else


#      clean up
       /bin/rm func_input holder
       /bin/rm finals.data
       /bin/rm utcpole.dat.NEW
       exit 23
   fi

;;

      2)

;;

      3)

;;


      4)


;;

      5)

;;

      6)

;;

      7)

;;

      8)

;;

      9)

;;

      10)

;;

      11)

;;

      12)

;;

      13)

;;

       14)

;;

       15)

;;

      *)

;;

    esac

#   clean up files 
   if [ -f func_input ] ; then /bin/rm func_input ; fi
   if [ -f holder ] ; then /bin/rm holder ; fi
   if [ -f finals.data ] ; then /bin/rm finals.data ; fi
   if [ -f utcpole.dat.NEW ] ; then /bin/rm utcpole.dat.NEW ; fi

   exit $state
# $Log$
