# $Id$
#
# This make include file defines how to make the MLSPGS libraries
#
mlslib:
	cd ${MLSPGS_LIB_PATH} ; ${MAKE} mlslib.a
mlslib_debug:
	cd ${MLSPGS_LIB_PATH} ; ${MAKE} mlslib.a_debug
#
# $Log$
# Revision 1.1  1999/12/04 00:27:07  livesey
# First version
#
#
