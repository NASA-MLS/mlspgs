#! /usr/bin/python
# Copyright 2010, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"
# Usage
# (1) pydiff.py a b
#     remove string b from string a
#     e.g., if a is "p,q,r"
#              b is "q"
#     return        "p,r"
# (2) pydiff.py a b s
#     treat s as the separator in (1) instead of ","
# (3) pydiff.py -s b a
#     treat s as the separator in (1) instead of ","
#     (Note the "-" in front of s and reversal of a and b)
import sys
from pyutils import diff

# Find which usage (1) - (3)
if len(sys.argv) > 3:
  if sys.argv[1][0] == "-":
# usage (3)
    s = sys.argv[1][1]
    l = sys.argv[3]
    u = sys.argv[2]
  else:
# usage (2)
    s = sys.argv[3][0]
    l = sys.argv[1]
    u = sys.argv[2]
else:
# usage (1)
  s = ","
  # read the args
  l = sys.argv[1]
  u = sys.argv[2]

# remove unwanted chars
l = l.replace("'",'')
l = l.replace(' ','')
l = l.replace(s,' ')  
l = l.replace('[','')  
l = l.replace(']','')

u = u.replace("'",'')
u = u.replace(' ','')
u = u.replace(s,' ')  
u = u.replace('[','')  
u = u.replace(']','')

# convert to lists
# (instead of converting s to ' ' above, why not just split on s?)
l = l.split(' ')
u = u.split(' ')

# print string for m4 to grab
print diff(l,u),

# $Log$
