#!/bin/sh
#     Suffix_to_line.sh
# Suffixes more text at ends of lines containing (FC) -c
# in the list of files among the args to the routine
#
if [ $# = 0 ] ; then
	echo "Usage: Suffix_to_line-out.sh list_of_file_names"
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
			echo "sed 's/(FC) -c.*$/& \$(FAFTER)/' $arg"
			sed -n 's/(FC) -c.*$/& $(FAFTER)/p' "$arg"
#			echo "sed 's/\$(l_specials)/\$(l_specials) \$(LAFTER)/' $arg"
#			sed -n 's/$(l_specials)/$(l_specials) $(LAFTER)/p' "$arg"
		fi

		if [ "$replacing" = "on" ] ; then
			sed 's/(FC) -c.*$/& $(FAFTER)/' "$arg" > temp
#			sed 's/$(l_specials)/$(l_specials) $(LAFTER)/' "$arg" > temp
			chmod u+w "$arg"
			mv temp "$arg"
		fi
	fi
done
