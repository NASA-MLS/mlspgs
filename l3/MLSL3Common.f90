
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE MLSL3Common
!===============================================================================

   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Remarks:  This module defines quantities common to multiple modules of L3.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD1 = 'Latitude'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD2 = 'Longitude'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD3 = 'Time'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD4 = 'LocalSolarTime'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD5 = 'SolarZenithAngle'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD6 = 'LineOfSightAngle'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD7 = 'OrbitGeodeticAngle'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD8 = 'ChunkNumber'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD9 = 'Pressure'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD10 = 'Frequency'

   CHARACTER (LEN=*), PARAMETER :: DIM_NAME1 = 'nTimes'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME2 = 'nLevels'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME3 = 'nFreqs'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME12 = 'nLevels,nTimes'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME123 = 'nFreqs,nLevels,nTimes'

   CHARACTER (LEN=*), PARAMETER :: DIML_NAME = 'nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMT_NAME = 'TDim'
   CHARACTER (LEN=*), PARAMETER :: DIMX_NAME = 'XDim'
   CHARACTER (LEN=*), PARAMETER :: DIMY_NAME = 'YDim'
   CHARACTER (LEN=*), PARAMETER :: DIMZ_NAME = 'ZDim'
   CHARACTER (LEN=*), PARAMETER :: DIMLL_NAME = 'nLevels,nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMXYZ_NAME = 'ZDim,YDim,XDim'

   CHARACTER (LEN=*), PARAMETER :: DAT_ERR = 'Failed to define data field '
   CHARACTER (LEN=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
   CHARACTER (LEN=*), PARAMETER :: GEO_ERR = 'Failed to define geolocation &
                                             &field '
   CHARACTER (LEN=*), PARAMETER :: GD_ERR = 'Failed to detach from grid '
   CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = 'Error writing metadata &
                                                &attribute '
   CHARACTER (LEN=*), PARAMETER :: NOOUT_ERR = ' data expected but not found &
                                                &for output.'
   CHARACTER (LEN=*), PARAMETER :: TAI2A_ERR = 'Error converting time from &
                                               &TAI to UTC.'
   CHARACTER (LEN=*), PARAMETER :: WR_ERR = 'Failed to write field '
   CHARACTER (LEN=*), PARAMETER :: SZ_ERR = 'Failed to get size of dimension '

   INTEGER, PARAMETER :: CCSDS_LEN = 27
   INTEGER, PARAMETER :: CCSDSB_LEN = 25
   INTEGER, PARAMETER :: GCTP_GEO = 0
   INTEGER, PARAMETER :: HDFE_NOMERGE = 0
   INTEGER, PARAMETER :: INVENTORYMETADATA = 2
   INTEGER, PARAMETER :: GridNameLen = 64
   INTEGER, PARAMETER :: maxGridPoints = 500
   INTEGER, PARAMETER :: maxWindow = 30

!=====================
END MODULE MLSL3Common
!=====================

!# $Log$
!# Revision 1.5  2001/03/27 19:32:12  nakamura
!# Added some grid parameters.
!#
!# Revision 1.4  2001/02/21 20:52:21  nakamura
!# Added SZ_ERR parameter.
!#
!# Revision 1.3  2001/02/09 19:18:33  nakamura
!# Added some dimension parameters.
!#
!# Revision 1.2  2001/01/16 17:46:18  nakamura
!# Added parameter for metadata writing error.
!#
!# Revision 1.1  2000/12/29 20:52:48  nakamura
!# Module for parameters common across L3 routines.
!#
