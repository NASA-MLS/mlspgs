#! /usr/bin/python3
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
# Return the input list with duplicates removed
# Usage
# (1) pyuniq.py a
#     remove duplicates from string a
#     e.g., if a is "p,q,r,p"
#     return        "p,q,r"
# (2) pyuniq.py -s a
#     treat s as the separator in (1) instead of ","
# (3) pyuniq.py -c -s a
#     treat s as the separator in (1) instead of ","
#     and treat list a as case insensitive
import sys
from pyutils import uniq

casesense = True
# Find which usage (1) - (3)
if len(sys.argv) > 3:
# usage (3)
    s = sys.argv[2][1]
    a = sys.argv[3]
    a = a.replace(s,',')
    casesense = False
elif len(sys.argv) > 2:
# usage (2)
    s = sys.argv[1][1]
    a = sys.argv[2]
    a = a.replace(s,',')
else:
# usage (1)
  s = ","
  a = sys.argv[1]

if casesense == False:
  a = a.upper()

# print string for m4 to grab
print(uniq(a), end=" ")

# $Log$
# Revision 1.1  2010/05/15 00:31:11  pwagner
# First commit
#
