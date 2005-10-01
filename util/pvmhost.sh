#!/bin/sh
# --------------- pvmhost help
#pvmhost.sh
#add/delete the pvmhost[s] to/from the currently available ones
#
#Usage:
#pvmhost.sh option1 [option2] .. pvmhost1 [pvmhost2] ..
#
#    O p t i o n s
# -dryrun     don't run the actual pvm commands
# -add        add the pvmhosts to list
# -delete     delete the pvmhosts from list
# -mstat      print the whole status line of the pvmhost; e,g,
#              c0-181  No such host
# -s          show the status of the pvmhost [0 is ok, 1 is not]
#              (whose value can be used in tests conveniently)
# -h[elp]     print brief help message; exit
#
# Note:
# 
# Result:
# Won't try to spawn jobs to deficient nodes
# --------------- End pvmhost help
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
me="$0"
my_name=pvmhost
I=pvmhost
OK=0
notOK=1
dryrun="no"
command=""
more_opts="yes"
while [ "$more_opts" = "yes" ] ; do

    case "$1" in

    -dryrun )
	    shift
       dryrun="yes"
       ;;
    -h | -help )
       sed -n '/'$my_name' help/,/End '$my_name' help/ p' $me \
           | sed -n 's/^.//p' | sed '1 d; $ d'
       exit
	     ;;
    -* )
       command=`echo "$1" | sed 's/-//'`
       shift
       ;;
    * )
       more_opts="no"
       ;;
    esac
done
# if [ $# -lt 1 ]
# then
#   echo "Usage: pvmhost.sh pvmhost1 [pvmhost2] .."
#   exit 1
# fi
for pvmhost
do
  if [ "$dryrun" = "yes" ]
  then
    echo "$command from pvm table $pvmhost"
  else
    case "$command" in
      delete | add | mstat )
      echo "$command from pvm table $pvmhost"
      pvm << EOF
      $command $pvmhost
      quit
EOF
    ;;
      s )
      result=`pvm << EOFs
      mstat $pvmhost
      quit
EOFs
`
      s=`echo $result | grep ' ok '`
      if [ "$s" = "" ]
      then
        echo $notOK
      else
        echo $OK
      fi
    ;;
      * )
      echo "$command from pvm table $pvmhost (not recognized)"
      pvm << EOFt
      $command $pvmhost
      quit
EOFt
    ;;
    esac
  fi
done
exit 0
# $Log$
