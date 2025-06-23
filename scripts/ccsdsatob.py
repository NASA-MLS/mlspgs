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

print "%04d-%03d" % ( date.year(), date.dayOfYear() )

    
