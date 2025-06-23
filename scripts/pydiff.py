#! /usr/bin/python3
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
#from __future__ import print_function
import sys
from pyutils import diff, mrClean

# Find which usage (1) - (3)
if len(sys.argv) > 3:
  if sys.argv[1][0] == "-":
# usage (3)
    s = sys.argv[1][1]
    l = mrClean(sys.argv[3],s)
    u = mrClean(sys.argv[2],s)
  else:
# usage (2)
    s = sys.argv[3][0]
    l = mrClean(sys.argv[1],s)
    u = mrClean(sys.argv[2],s)
else:
# usage (1)
  s = ","
  # read the args
  l = mrClean(sys.argv[1],s)
  u = mrClean(sys.argv[2],s)

# convert to lists
# (instead of converting s to ' ' above, why not just split on s?)
l = l.split(' ')
u = u.split(' ')

# print string for m4 to grab
# print("l", file=sys.stderr, end="=")
# print(l, file=sys.stderr)
# print("u", file=sys.stderr, end="=")
# print(u, file=sys.stderr)
print(diff(l,u), end=" ")

# $Log$
# Revision 1.3  2010/05/21 23:54:31  pwagner
# Use mrClean
#
# Revision 1.2  2010/05/15 00:30:45  pwagner
# Added usages (2) and (3)
#
