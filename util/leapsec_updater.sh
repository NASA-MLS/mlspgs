#!/bin/sh
# leapsec_updater.sh
# Plagiarised Raytheon's update_leapsec.sh
# distributed with the Toolkit
# but I've removed all the email and user-input stuff
# Should be useable by a cron job
# or you can use it directly
#
# REQUIRES:
# PGSHOME must have been defined (by e.g., pgs-env.csh)
# write access to $PGSHOME/database/common/TD
 
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
if [ ! -w $PGSHOME/database/common/TD/leapsec.dat ]
then
  echo "You lack write access permission to the file leapsec.dat"
  echo "Likely causes are"
  echo "(1) It is a directory belonging to someone else"
  echo "(2) The name or location of the directory has been changed"
  echo "(3) A prior attempt to update the file badly misfired"
  exit 1
fi

# change directory to the leapsec.dat directory:

cd $PGSHOME/database/common/TD

#
#ftp to USNO  (maia.usno.navy.mil 192.5.41.22 ) and get the current tai-utc.dat in ser7

#echo "Are you a DAAC-location [yes/no]"
#read user_response
user_response="No"
case  "$user_response" in
    y* | Y* )
echo "Please enter your DAAC-specific ftp machine"
read ftp_machine

echo "Please enter your user name"
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
get tai-utc.dat
bye
GG
    ;;
    * )
echo "Please wait ...."
ftp -n <<GG>holder
open maia.usno.navy.mil
user anonymous EOSuser
prompt
cd ser7
dir
get tai-utc.dat
bye
GG
    ;;
    esac


#parse the results of ftp for the date, to be used in the new header

grep "tai-utc.dat" holder | awk '{print  " " $6 " " $7 " " $8}'  > func_input

#clean up a bit
/bin/rm holder

#run update C program PGS_TD_NewLeap:
#PGS_TD_NewLeap.c assumed compiled and in $PGSBIN (done by the INSTALL script)
# if the executable PGS_TD_NewLeap has been lost or corrupted, "cd" to the
# directory $PGSSRC/TD and at the keyboard enter "make utilities"

   PGS_TD_NewLeap < func_input

   state=$?

echo "Status of PGS_TD_NewLeap call was ($state)"

#test return value

  case $state in

       0)

#All is OK  - replace file

   /bin/mv -f leapsec.dat.new $PGSHOME/database/common/TD/leapsec.dat

   mv_state=$?

#report results of the mv

   echo "Status of MOVE command was ($mv_state)"
 
 if [ -f holder ]; then /bin/rm holder ; fi
 if [ -f tai-utc.dat ]; then /bin/rm tai-utc.dat ; fi
 if [ -f func_input ]; then /bin/rm func_input ; fi

       if [ $mv_state -eq 0 ] ; then

       exit 0
 
      else

   /bin/rm -f leapsec.dat.new
 
#error branch within if - used if "mv" fails
 
 
       exit 22

fi
#EXITED if test on bin/mv
 
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

       9)


;;

       12)


;;
       13)


;;
       14)
 
 
;;

       15)


;;
       18)


;;

      *)
 
;;

esac

#clean up files

 if [ -f tai-utc.dat ]; then /bin/rm tai-utc.dat ; fi
 if [ -f func_input ]; then /bin/rm func_input ; fi
 if [ -f leapsec.dat.new ]; then /bin/rm leapsec.dat.new ; fi
 if [ -f holder ]; then /bin/rm holder ; fi

exit $state
# $Log$
