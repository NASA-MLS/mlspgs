#!/bin/bash
# --------------- clusterstatus help
# Create status files for level 2 jobs running on lightspeed, scramjet
# Optionally, mail or scp it to me
# Usage:
# clusterstatus.sh [options]
#    O p t i o n s
# -dryrun     don't run the sipsl2.sh script
# -mail       mail the file to RECIPIENTS
# -scp        scp the file to me
# -h[elp]     print brief help message; exit
#
# Note:
# (1) see notes for /home/pwagner/bin/sipsl2.sh
# (2) scp option does not seem to work
# (3) so use mail instead
# 
# Result:
# a file showing jobs, machine, date, status, etc.
# --------------- End clusterstatus help
# Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

# "$Id$"

# ************
# Main Program
# ************
# 
RECIPIENTS="pwagner David.T.Cuddy@jpl.nasa.gov ahanzel@mls.jpl.nasa.gov pzimdars@sdsio-mail.jpl.nasa.gov sneely@sdsio.jpl.nasa.gov eparaiso@sdsio-mail.jpl.nasa.gov bsaha@sdsio-mail.jpl.nasa.gov dromo@sdsio.jpl.nasa.gov cvuu@mls.jpl.nasa.gov"
#RECIPIENTS="cvuu@mls.jpl.nasa.gov paul.a.wagner@jpl.nasa.gov"
#RECIPIENTS="pwagner"
MAILER="/home/pwagner/bin/mailtome.sh"
debug="no"
me="$0"
my_name=clusterstatus
I=clusterstatus
dryrun=""
mail="no"
scp="no"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    shift
       dryrun="yes$dryrun"
       ;;
    -mail )
	    shift
       mail="yes"
       ;;
    -scp )
	    shift
       scp="yes"
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

TODAY=`date +"%Y-%m-%d"`
OUTPUT=/home/pwagner/clusterstatus/"$TODAY".txt
mustrun="yes"
if [ -f "$OUTPUT" ]
then
  if [ "$dryrun" != "" ]
  then
    mustrun="no"
  else
    mv "$OUTPUT" "$OUTPUT.1"
  fi
fi
if [ "$mustrun" = "yes" ]
then
  echo "Cluster status for $TODAY" > "$OUTPUT"
  /home/pwagner/bin/sipsl2.sh -c -x >> "$OUTPUT"
  echo "" >> "$OUTPUT"
  echo "Status of running jobs:" >> "$OUTPUT"
  /home/pwagner/bin/sipsl2.sh -c -full -t >> "$OUTPUT"
  echo "" >> "$OUTPUT"
  echo "Nodes showing pvm failures:" >> "$OUTPUT"
  /home/pwagner/bin/sipsl2.sh -fail >> "$OUTPUT"
else
  echo "Cluster status for $TODAY $OUTPUT"
  echo "/home/pwagner/bin/sipsl2.sh -c -x"
  echo ""
  echo "Status of running jobs: >> $OUTPUT"
  echo "/home/pwagner/bin/sipsl2.sh -c -full -t"
  echo ""
  echo "Nodes showing pvm failures: >> $OUTPUT"
  echo "/home/pwagner/bin/sipsl2.sh -fail"
fi
if [ "$mail" = "yes" -a "$dryrun" != "yesyes" ]
then
  $MAILER -r /home/pwagner/maillogs -s "clusterstatus $TODAY" \
   -R "$RECIPIENTS" "$OUTPUT"
elif [ "$mail" = "yes" ]
then
  $MAILER -r /home/pwagner/maillogs -s "clusterstatus $TODAY" \
   -R "$RECIPIENTS" -dryrun "$OUTPUT"
elif [ "$scp" = "no" ]
then
  exit 0
elif [ "$dryrun" = "no" ]
then
  echo 'Warning: scp not yet properly tested'
  # /home/pwagner/private/scptome.sh -r /home/pwagner/scplogs "$OUTPUT" 
else
  /home/pwagner/private/scptome.sh -dryrun "$OUTPUT" 
fi
exit 0
# $Log$
