#!/bin/env python

################################################################
## readl2gp.py
################################################################
## 
## Author: Brian Mills
## Date: 2010 08 19
## Based on readl2gp.pro by Nathaniel Livesey 2009 09 16
##
## Depends on:
## * Python 2.4+ http://www.python.org
## * PyTables http://www.pytables.org
## * Numpy http://www.numpy.org
##
## Copyright 2010, by the California Institute of
## Technology. ALL RIGHTS RESERVED. United States Government
## Sponsorship acknowledged. Any commercial use must be
## negotiated with the Office of Technology Transfer at the
## California Institute of Technology. 
## 
## This software may be subject to U.S. export control
## laws. By accepting this software, the user agrees to
## comply with all applicable U.S. export laws and
## regulations. User has the responsibility to obtain export
## licenses, or other export authority as may be required
## before exporting such information to foreign countries or
## providing access to foreign persons. 
## 
################################################################


# REQUIREMENT: PyTables for HDF5 reading capabilities
# http://www.pytables.org
import tables
# Note that numpy isn't directly required for this file, but
# results from PyTables are in numpy array formats, so it's required
# to do anything with the results (and probably required for PyTables...)

def readl2gp(filename,
              swathName=None,
              variableName=None,
              precisionName=None,
              noError=False):
    """
Read a MLS L2GP swath from the given file and return a dictionary containing
a structured version of that information.

Unless the swath is specified, only the first one is read and returned--but this
is usually the only swath in the file.  Either way, the name of the swath read
is part of the return.

TODO: Deeper error checking is forthcoming.

PyTables and Numpy are both requirements to run this function.

ARGUMENTS:
* filename
    The full or relative path to the L2GP HDF-5 file to be read
* swathName

    The name of the HDF-EOS swath.  If blank, the first swath in the
    file will be read by default. For files with multiple swaths (most
    files for v4, some for v3) you'll have to pass the swath name to
    get anything but the first swath

* variableName
    The name of the regular HDF variable entry.  "L2gpValue" by default.
* precisionName
    The name of the regular HDF precision entry.  "L2gpPrecision" by default.
* noError
    If set, exceptions will not be raised on error and instead we will
    silently return a None.

* returns: A numpy ND variable or None, if a failure occurred

* Note: One of the variable attributes contains a variable of type
  HE5_REFERENCE. PyTables doesn't like this, so it complains. But this
  only happens on the first call and it is not fatal. The data is
  still returned and valid. 

    """

    result = {}

    # Set up defaults
    if variableName is None:
        variableName = 'L2gpValue'
        precisionName = 'L2gpPrecision'
    elif precisionName is None:
        precisionName = variableName + 'Precision'

    #DESIGN_DECISION: We're assuming the file is HDF-5 only.
    fHandle = tables.openFile(filename, mode='r')

    # Query swaths
    swathGroup = fHandle.root.HDFEOS.SWATHS
    nSwaths = swathGroup._v_nchildren
    if nSwaths == 0:
        # No swaths in this file!!
        fHandle.close()
        if noError: return None
        else: raise IOError("No swaths in file")

    # If we specified a swath, go right there
    # Otherwise, take the first one
    swathNode = None
    if swathName is not None:
        result['swathName'] = swathName
        if swathName in swathGroup:
            swathNode = swathGroup._f_getChild(swathName)
        else:
            # Desired swath not here
            fHandle.close()
            if noError: return None
            else: raise IOError("Swath %s not in file!"%swathName)
    else:
        swathNode = swathGroup._f_listNodes()[0]
        result['swathName'] = swathNode._v_name

    # Get the geolocation fields
    geoLocNode = swathNode._f_getChild('Geolocation Fields')
    result['time'] = geoLocNode._f_getChild('Time').read()
    result['latitude'] = geoLocNode._f_getChild('Latitude').read()
    result['longitude'] = geoLocNode._f_getChild('Longitude').read()

    # Sort out time dimensions
    result['nTimes'] = len(result['time'])
    
    # These might not be in all files; tread carefully
    for (hdfName, dictName) in (
        ('ChunkNumber','chunkNumber'),
        ('LineOfSightAngle', 'lineOfSightAngle'),
        ('LocalSolarTime', 'localSolarTime'),
        ('OrbitGeodeticAngle', 'orbitGeodeticAngle'),
        ('SolarZenithAngle', 'solarZenithAngle')):
        
        if hdfName in geoLocNode:
            result[dictName] = geoLocNode._f_getChild(hdfName).read()

    # Sort out dimensions of fields that don't already exist
    for (hdfName, dictName, dimName) in (
        ('Pressure', 'pressure', 'nLevels'),
        ('Frequency', 'frequency', 'nFreqs'),
        ('Theta', 'theta', 'nTheta')):
        
        if hdfName in geoLocNode:
            result[dictName] = geoLocNode._f_getChild(hdfName).read()
            result[dimName] = len(result[dictName])
    
    # Now, get the data fields
    dataNode = swathNode._f_getChild('Data Fields')
    if variableName in dataNode:
        result['l2gpValue'] = dataNode._f_getChild(variableName).read()

        # Also create an attributes dict
        attrs = {}
        attr_node = dataNode._f_getChild(variableName)._v_attrs
        attrs['MissingValue'] = attr_node.MissingValue[0]
        attrs['Title'] = attr_node.Title
        attrs['UniqueFieldDefinition'] = attr_node.UniqueFieldDefinition
        attrs['Units'] = attr_node.Units
        attrs['_FillValue'] = attr_node._FillValue[0]
        result['attrs'] = attrs


    if precisionName in dataNode:
        result['l2gpPrecision'] = dataNode._f_getChild(precisionName).read()

    # If precision wasn't in the geoloc data, see if it's in data
    if 'pressure' not in result and 'Pressure' in dataNode:
        result['pressure'] = dataNode._f_getChild('Pressure').read()
        result['nLevels'] = len(result['pressure'])

    # Get status and quality if it's there
    for (hdfName, dictName) in (
        ('Status', 'status'),
        ('Quality', 'quality'),
        ('Convergence', 'convergence')):

        if hdfName in dataNode:
            result[dictName] = dataNode._f_getChild(hdfName).read()

    # Close the file, return the final result dictionary!
    fHandle.close()
    return result




#
# $Id$
#
# $Name$
#
# Modifications:
# $Log$
#                

