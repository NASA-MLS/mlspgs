#! /usr/bin/python
import sys
import string
from normalDate import NormalDate

a = sys.argv[1]

year = string.atoi (  a[0:4] )
month = string.atoi ( a[5:7] )
day = string.atoi ( a[8:10] )
date = NormalDate()
date.setNormalDate( ( year, month, day ) )
date = date + string.atoi( sys.argv[2] )

print "%04d-%02d-%02d" % ( date.year(), date.month(), date.day() )
    
