#!/bin/sh
#     Comment-out.sh
# Comments out the lines beginning with FOPTS =
# in the list of files among the args to the routine
#
if [ $# = 0 ] ; then
	echo "Usage: Comment-out.sh list_of_file_names"
	exit
fi

rm -f temp

debugging="off"

replacing="on"

for arg
do
	echo "$arg"
	if [ -f "$arg" ] ; then
		if [ "$debugging" = "on" ] ; then
			echo "sed 's/^FOPTS =/#FOPTS =/' $arg"
			sed 's/^FOPTS =/#FOPTS =/' "$arg"
		fi

		if [ "$replacing" = "on" ] ; then
			sed 's/^FOPTS =/#FOPTS =/' "$arg" > temp
			chmod u+w "$arg"
			mv temp "$arg"
		fi
	fi
done
