
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

   CHARACTER (LEN=*), PARAMETER :: DIML_NAME = 'nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMLL_NAME = 'nLevels,nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMT_NAME = 'TDim'

   CHARACTER (LEN=*), PARAMETER :: DAT_ERR = 'Failed to define data field '
   CHARACTER (LEN=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
   CHARACTER (LEN=*), PARAMETER :: GEO_ERR = 'Failed to define geolocation &
                                             &field '
   CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = 'Error writing metadata &
                                                &attribute '
   CHARACTER (LEN=*), PARAMETER :: TAI2A_ERR = 'Error converting time from &
                                               &TAI to UTC.'
   CHARACTER (LEN=*), PARAMETER :: WR_ERR = 'Failed to write field '
   CHARACTER (LEN=*), PARAMETER :: SZ_ERR = 'Failed to get size of dimension '

   INTEGER, PARAMETER :: CCSDS_LEN = 27
   INTEGER, PARAMETER :: CCSDSB_LEN = 25
   INTEGER, PARAMETER :: INVENTORYMETADATA = 2
   INTEGER, PARAMETER :: GridNameLen = 64
   INTEGER, PARAMETER :: maxWindow = 30

!=====================
END MODULE MLSL3Common
!=====================

!# $Log$
!# Revision 1.3  2001/02/09 19:18:33  nakamura
!# Added some dimension parameters.
!#
!# Revision 1.2  2001/01/16 17:46:18  nakamura
!# Added parameter for metadata writing error.
!#
!# Revision 1.1  2000/12/29 20:52:48  nakamura
!# Module for parameters common across L3 routines.
!#
!#
