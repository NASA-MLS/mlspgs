#! /usr/bin/python
import sys
import string
from normalDate import NormalDate

a = sys.argv[1]

year = string.atoi (  a[0:4] )
day = string.atoi ( a[5:8] )
date = NormalDate()
date.setNormalDate( ( year, 1, 1 ) )
date = date + day - 1
date = date + string.atoi( sys.argv[2] )

print "%04d-%03d" % ( date.year(), date.dayOfYear() )
    
