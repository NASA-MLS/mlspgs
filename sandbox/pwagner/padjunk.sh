#!/bin/sh
# padjunk.sh
# pads a file by putting the contents of a file named "junk" ahead
# of the original file's contents
#
# Usage
# padjunk.sh -nnn file_name
# nnn the number of times to pad file_name by
# file_name  the file to be padded
#
# restrictions: both junk and file_name must exist and junk
# must be in the current working directory

case "$1" in
	-[0-9]*)
        	count=`expr 0 - $1`
                shift
         	;;
          *)
          	count=1
                ;;
esac

while [ "$count" -gt 0 ]
do
	cat junk "$1" > "temp.$1"
        mv "temp.$1" "$1"
	count=`expr $count - 1`
done
exit
