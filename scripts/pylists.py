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
# Perform one of several possible list-oriented operations
# Usage
# (1) pylists.py -1 -s n a
#     return item number n from list a
#     e.g., if a is "p,q,r"
#              n is "2"
#     return        "q"
# (2) pylists.py -2 -s c a
#     return the item number which c occupies in list a
#     e.g., if a is "p,q,r"
#              c is "q"
#     return        "2"
# (3) pylists.py -3 -s key keys values
#     return values[key] where {keys,values} is a dict (associative array)
#     e.g., if keys is "p,q,r" and values is "jack,queen,king"
#              key is "q"
#     return          "queen"
# Notes:
# (a) we treat s as the separator between list items
# (b) (1) returns nothing if fewer than n items in list
# (c) (2) and (3) return nothing if 3rd arg not in list
# (d) we return nothing if too few args or if first arg is not one of "-1,-2,-3"

import sys
from pyutils import diff, mrClean

# Find which usage (1) - (3)
if len(sys.argv) > 4:
  s = sys.argv[2][1]
  # print "s: ", s
  # print "Usage: ", sys.argv[1][1]
  if sys.argv[1] == "-1":
# usage (1)
    n = int(mrClean(sys.argv[3],s))
    a = sys.argv[4]
    a = mrClean(a,s)
    a = a.split(' ')
    # print "n: ", n
    # print "a: ", a
# print string for m4 to grab
    if  len(a) > n - 1:
      print(a[n-1], end=" ")
  elif sys.argv[1] == "-2":
# usage (2)
    c = sys.argv[3]
    a = sys.argv[4]
    a = mrClean(a,s)
    a = a.split(' ')
    # print "c: ", c
    # print "a: ", a
    if c in a:
      print(a.index(c) + 1, end=" ")
  elif sys.argv[1] == "-3":
# usage (3)
    key = sys.argv[3]
    keys = mrClean(sys.argv[4],s)
    values = mrClean(sys.argv[5], s)
    keys = keys.split(' ')
    values = values.split(' ')
    # print "key: ", key
    # print "keys: ", keys
    # print "values: ", values
    if key in keys:
      n = keys.index(key)
      if len(values) > n:
        print(values[n], end=" ")

# $Log$
# Revision 1.1  2010/05/22 00:01:39  pwagner
# First commit
#
# Revision 1.2  2010/05/15 00:30:45  pwagner
# Added usages (2) and (3)
#
