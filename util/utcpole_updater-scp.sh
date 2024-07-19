#!/bin/sh
# utcpole_updater-scp.sh
# simplified: 
# (1) merely saves utcpole.dat in the current working directory
#     thereby eliminating all the PGSHOME stuff
# (2) Since we'll never use ftp again, eliminate that stuff, too
# (3) However, you must now manually cp the utcpole.dat to whereever it belongs
#
# Copyright 2023, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

#scp pwagner@fox.jpl.nasa.gov:/science/pge/v1300-nrt/toolkit5.2.18/database/common/CSC/utcpole.dat ./
scp pwagner@jackal.jpl.nasa.gov:/science/pge/v1800-nrt/toolkit5.2.18/database/common/CSC/utcpole.dat ./

state=$?

if [ $state != 0 ];  then 

exit $state
fi

exit 0
# $Log$
# Revision 1.2  2024/01/04 22:21:23  pwagner
# Simplified--no longer uses PGSHOME
#
# Revision 1.1  2023/05/25 21:53:24  pwagner
# ftp no longer supported; scp from sips machine instead
#
