#! /usr/bin/python3
# $Id$
# remove repeated molecules in the lbl molecule list used in the PFA forward model
# handle (R1A,R1B) and (R5H,R5V) duplicate definitions
# Converted to python3 (2021-12-22)

import sys
from pyutils import uniq 

# read the args
s = sys.argv[1:]

#print 'xxx'+str(s)+'xxx'

# force upper case
u = str(s).upper()

# replace 'NONE' with ''
if u.find('NONE') >= 0 : u = u.replace('NONE','')

# if both _R1A and _R1B are present then return only _R1A
if u.find('_R1A') >= 0 & u.find('_R1B') >= 0 : u = u.replace('_R1B','_R1A')

# if both _R5H and R5V are present then return only _R5H
if u.find('_R5H') >= 0 & u.find('_R5V') >= 0 : u = u.replace('_R5V','_R5H')

# print uniquified string for m4 to grab
print(uniq(u), end=" ")
# $Log$
