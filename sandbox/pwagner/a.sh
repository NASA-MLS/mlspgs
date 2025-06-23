#!/bin/sh
# a silly test script that either 
# (depending on whether it is called with an argument or not)
# (1) sets a to the value given as arg[1]; or
# (2) asks the user for a value, and sets a to that
#
echo "Currently, a is $a"
if [ $# = 0 ] ; then
	/bin/echo -n "Please enter a value for a"
	read user_response
	a="$user_response"
else
	a="$1"
fi
echo "a is now $a"

conf_dir="."
rm -f $conf_dir/.configure
echo "#	.configure for mls" > $conf_dir/.configure
echo 'a="'$a'"' >> $conf_dir/.configure
echo "export a" >> $conf_dir/.configure

