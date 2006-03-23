#!/bin/bash
# --------------- clusterstatus help
# Create status files for level 2 jobs running on mls clusters at the sips
# Optionally, mail or scp it to me
# Usage:
# clusterstatus.sh [options]
#    O p t i o n s
# -debug        print lots of extra stuff
# -dryrun       don't run the sipsl2.sh script
# -vn versions  show separate listings for versions
#                (e.g., "V01-51,V01-52")
# -R list       change RECIPIENTS to list
# -mail         mail the file to RECIPIENTS
# -scp          scp the file to me
# -[n]sort      [don't] sort the initial table according to machine
#                (sort by default)
# -temp         don't save the resulting clusterstatus file
# -h[elp]       print brief help message; exit
#
# Note:
# (1) see notes for /home/pwagner/bin/sipsl2.sh
# (2) scp option does not seem to work
# (3) so use mail instead
# 
# Result:
# a file showing jobs, machine, date, status, etc.
# --------------- End clusterstatus help
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

# ************
# Main Program
# ************
# 
#RECIPIENTS="pwagner David.T.Cuddy@jpl.nasa.gov ahanzel@mls.jpl.nasa.gov sysadmin@sdsio.jpl.nasa.gov eparaiso@sdsio-mail.jpl.nasa.gov bsaha@sdsio-mail.jpl.nasa.gov dromo@sdsio.jpl.nasa.gov cvuu@mls.jpl.nasa.gov Brian.W.Knosp@jpl.nasa.gov"
RECIPIENTS="pwagner David.T.Cuddy@jpl.nasa.gov ahanzel@mls.jpl.nasa.gov sysadmin@sdsio.jpl.nasa.gov mliukis@sdsio-mail.jpl.nasa.gov bsaha@sdsio-mail.jpl.nasa.gov dromo@sdsio.jpl.nasa.gov cvuu@mls.jpl.nasa.gov Brian.W.Knosp@jpl.nasa.gov"
MAILER="/home/pwagner/bin/mailtome.sh"
# clusternames="lightspeed scramjet speedracer"
clusternames="lightspeed speedracer"
debug="no"
me="$0"
my_name=clusterstatus
I=clusterstatus
dryrun=""
mail="no"
scp="no"
sort="yes"
temp="no"
versions="(none)"
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    shift
       dryrun="yes$dryrun"
       ;;
    -debug )
	    shift
       debug="yes"
       ;;
    -R )
	    shift
       RECIPIENTS="$1"
       shift
       ;;
    -vn )
	    shift
       versions="$1"
       shift
       ;;
    -mail )
	    shift
       mail="yes"
       ;;
    -scp )
	    shift
       scp="yes"
       ;;
    -nsort )
	    shift
       sort="no"
       ;;
    -sort )
	    shift
       sort="yes"
       ;;
    -temp )
	    shift
       temp="yes"
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
if [ "$temp" = "yes" ]
then
  OUTPUT=/home/pwagner/clusterstatus/"$TODAY".tmp
else
  OUTPUT=/home/pwagner/clusterstatus/"$TODAY".txt
fi

versions=`echo $versions | sed 's/,/ /g'`
if [ "$dryrun" != "" ]
then
  mustrun="no"
else
  mustrun="yes"
fi
if [ -f "$OUTPUT" ]
then
  if [ "$dryrun" != "" ]
  then
    mustrun="no"
  else
    mv "$OUTPUT" "$OUTPUT.1"
  fi
fi
options="-c -bug -x"

if [ "$debug" = "yes" ]
then
  echo "dryrun $dryrun"
  echo "mail $mail"
  echo "sort $sort"
  echo "TODAY $TODAY"
  echo "OUTPUT $OUTPUT"
  echo "mustrun $mustrun"
  echo "options $options"
  echo "versions $versions"
  echo "clusternames $clusternames"
  echo "MAILER $MAILER"
  echo "RECIPIENTS $RECIPIENTS"
fi

if [ "$mustrun" = "yes" ]
then
  echo "Cluster status for $TODAY" > "$OUTPUT"
fi

for version in $versions
do
  sipsoptions="$options"

  if [ "$mustrun" = "yes" ]
  then
    if [ "$version" != "(none)" ]
    then
      echo "(version $version)" >> "$OUTPUT"
      sipsoptions="-vn $version $options"
    fi
    if [ "$sort" = "yes" ]
    then
      for name in $clusternames
      do
        machoptions="-m $name $sipsoptions"
        # echo "/home/pwagner/bin/sipsl2.sh $machoptions"
        /home/pwagner/bin/sipsl2.sh $machoptions >> "$OUTPUT"
        # /home/pwagner/bin/sipsl2.sh -c -bug -x -m $name >> "$OUTPUT"
        echo "" >> "$OUTPUT"
      done
    else
      # echo "/home/pwagner/bin/sipsl2.sh $sipsoptions"
      /home/pwagner/bin/sipsl2.sh $sipsoptions >> "$OUTPUT"
      # /home/pwagner/bin/sipsl2.sh -c -bug -x >> "$OUTPUT"
      echo "" >> "$OUTPUT"
    fi
  else
    echo "Cluster status for $TODAY $OUTPUT"
    if [ "$version" != "(none)" ]
    then
      echo "(version $version)"
      sipsoptions="-vn $version $sipsoptions"
    fi
    if [ "$sort" = "yes" ]
    then
      for name in $clusternames
      do
        machoptions="-m $name $sipsoptions"
        echo "/home/pwagner/bin/sipsl2.sh $machoptions"
        # echo "/home/pwagner/bin/sipsl2.sh -c -bug -x $name"
        echo ""
      done
    else
      echo "/home/pwagner/bin/sipsl2.sh $sipsoptions"
      # echo "/home/pwagner/bin/sipsl2.sh -c -bug -x"
      echo ""
    fi
  fi

done

if [ "$mustrun" = "yes" ]
then
  echo "Status of running jobs:" >> "$OUTPUT"
  # echo "/home/pwagner/bin/sipsl2.sh -c -full -t >> $OUTPUT"
  /home/pwagner/bin/sipsl2.sh -c -full -t >> "$OUTPUT"
  echo "" >> "$OUTPUT"
  echo "Nodes showing pvm failures:" >> "$OUTPUT"
  # echo "/home/pwagner/bin/sipsl2.sh -fail >> $OUTPUT"
  /home/pwagner/bin/sipsl2.sh -fail >> "$OUTPUT"
else
  echo "Status of running jobs: >> $OUTPUT"
  echo "/home/pwagner/bin/sipsl2.sh -c -full -t"
  echo ""
  echo "Nodes showing pvm failures: >> $OUTPUT"
  echo "/home/pwagner/bin/sipsl2.sh -fail"
fi

if [ "$mail" = "yes" -a "$dryrun" != "yesyes" ]
then
  $MAILER -r /home/pwagner/maillogs -s "clusterstatus $TODAY" \
   -f "Paul Wagner <paul.a.wagner@jpl.nasa.gov>" -F -R "$RECIPIENTS" "$OUTPUT"
elif [ "$mail" = "yes" ]
then
  $MAILER -r /home/pwagner/maillogs -s "clusterstatus $TODAY" \
   -f "Paul Wagner <paul.a.wagner@jpl.nasa.gov>" -F -R "$RECIPIENTS" -dryrun "$OUTPUT"
elif [ "$scp" = "no" ]
then
  exit 0
elif [ "$dryrun" = "no" ]
then
  echo 'Warning: scp not yet properly tested'
else
  /home/pwagner/private/scptome.sh -dryrun "$OUTPUT" 
fi
exit 0
# $Log$
# Revision 1.6  2005/10/29 00:19:48  pwagner
# Changes to have mesgs show Return-path we want
#
# Revision 1.5  2005/09/09 16:43:20  pwagner
# Changed args to MAILER to forge return address to jpl address
#
# Revision 1.4  2005/06/23 22:20:45  pwagner
# Reworded Copyright statement
#
# Revision 1.3  2005/05/20 23:10:49  pwagner
# Show how many chunks succumbed to 'list out of order' bug
#
# Revision 1.2  2005/04/21 20:39:48  pwagner
# Improved readability, added speedracer name
#
# Revision 1.1  2005/04/06 22:16:22  pwagner
# First commit
#
