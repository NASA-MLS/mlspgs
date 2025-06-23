#! /usr/bin/python3

# calculate magnetic field vector from the
# components fieldStrength, fieldAzimuth, fieldElevation
# Not debugged under python3 yet (2021-12-16)
import sys
from math import sin, cos, pi

# read the args
fieldStrength = float(sys.argv[1])
fieldAzimuth = float(sys.argv[2])
fieldElevation = float(sys.argv[3])

dtor = pi/180.

B  = [fieldStrength*sin(fieldElevation*dtor)*cos(fieldAzimuth*dtor), fieldStrength*sin(fieldElevation*dtor)*sin(fieldAzimuth*dtor), fieldStrength*cos(fieldElevation*dtor)]

print( ('%15.6E, %15.6E, %15.6E' % ( B[0], B[1], B[2] ) ), end=" ")
# $Log$
