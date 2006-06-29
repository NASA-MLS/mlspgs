#! /usr/bin/python
# $Id$

import sys
from pyutils import diff

# read the args
l = sys.argv[1]
u = sys.argv[2]

# remove unwanted chars
l = l.replace("'",'')
l = l.replace(' ','')
l = l.replace(',',' ')  
l = l.replace('[','')  
l = l.replace(']','')

u = u.replace("'",'')
u = u.replace(' ','')
u = u.replace(',',' ')  
u = u.replace('[','')  
u = u.replace(']','')

# convert to lists
l = l.split(' ')
u = u.split(' ')

# print string for m4 to grab
print diff(l,u),

