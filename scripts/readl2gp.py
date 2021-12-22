#!/bin/env python3

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

# -------------------------
# Not tested under python3 yet
# Who uses this script?
# If you use it, and it fails (!)
# please contact Paul Wagner
# or, even better, repair it
# and then commit the repaired
# version to the cvs repository
# -------------------------

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

PyTables and Numpy are both requirements to run this function.

ARGUMENTS:
* filename
    The full or relative path to the L2GP HDF-5 file to be read (required!)

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

* returns: Failure: None. If `filename' doesn't exist, an exception
                    will be thrown and a traceback printed
 
           Success: A dictionary various tags, the most important of
           which are...

          'nTimes': the number of profiles
          'nLevels': the number of surfaces

          'l2gpValue': the actual data, a (nTimes by nLevels) Numpy ND
                       array

          'l2gpPrecision': The precision for the data, a (nTimes by
                            nLevels) Numpy ND Array

          'pressure' : the pressure surfaces (the nLevels dimension)
                        as Numpy array

          'time' : time as TAI93 (nTimes ), as Numpy array
          'longitude' : longitude (nTimes)
          'latitude': latitude (nTimes)

          'status' : a 'nTimes' byte array that gives status
          information. In general. See the quality document for how to
          interpret this field.
          
          'attrs': the attributes. Contains the keys 'Units', 'Title',
                   '_FillValue' and 'MissingValue'. The last two are
                   the value of data to avoid.



* Note: One of the attributes for l2gpValue contains a variable of
        type HE5_REFERENCE. PyTables doesn't like this, so it
        complains. But this only happens on the first call and it is
        not fatal. The data is still returned and is valid.

* Example

    % ipython # start interactive python interpreter
    # Issue commands to it.

    In [1]: from readl2gp import *
    In [2]: data=readl2gp('/path/to/your/data/file.he5')

    In [3]: value     = data['l2gpValue']
    In [4]: precision = data['l2gpPrecision']
    In [5]: time      = data['time']
    In [6]: latitude  = data['latitude']
    In [7]: longitude = data['longitude']
    In [8]: pressure  = data['pressure']
    In [9]: status    = data['status']

    In [10]: baddata=data['attrs']['MissingData']

    # use status and baddata to quality check your data and 
    # proceed from there.

    To read this documentation, do

    In [11]: help(readl2gp)


    Calling it from a python script is simply a matter of issuing the
    import statement, followed by the read (statements 1 and 2 above)
    in your script, and pulling out whatever data you wish from the
    returned dictionary.
  
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
# Revision 1.3  2015/03/17 20:20:40  whdaffer
# Small mods to documentation. Add a readme file
#
# Revision 1.2  2015/02/13 18:17:15  whdaffer
# Added comments, removed commented out code
#
#                

