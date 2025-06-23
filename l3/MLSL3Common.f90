
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
MODULE MLSL3Common
!==============================================================================

   USE MLSCommon, ONLY: FileNameLen
   IMPLICIT NONE
   PUBLIC


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

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
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD11 = 'Date'

   CHARACTER (LEN=*), PARAMETER :: DG_FIELD = 'GRss'
   CHARACTER (LEN=*), PARAMETER :: MD_FIELD = 'MaxDiff'
   CHARACTER (LEN=*), PARAMETER :: MDT_FIELD = 'MaxDiffTime'
   CHARACTER (LEN=*), PARAMETER :: DG_FIELD1 = 'LatRss'
   CHARACTER (LEN=*), PARAMETER :: DG_FIELD2 = 'PerMisPoints'

   CHARACTER (LEN=*), PARAMETER :: DIM_NAME1 = 'nTimes'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME2 = 'nLevels'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME3 = 'nFreqs'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME12 = 'nLevels,nTimes'
   CHARACTER (LEN=*), PARAMETER :: DIM_NAME123 = 'nFreqs,nLevels,nTimes'

   CHARACTER (LEN=*), PARAMETER :: DIML_NAME = 'nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMN_NAME = 'N'
   CHARACTER (LEN=*), PARAMETER :: DIMR_NAME = 'RDim'
   CHARACTER (LEN=*), PARAMETER :: DIMT_NAME = 'TDim'
   CHARACTER (LEN=*), PARAMETER :: DIMX_NAME = 'XDim'
   CHARACTER (LEN=*), PARAMETER :: DIMY_NAME = 'YDim'
   CHARACTER (LEN=*), PARAMETER :: DIMZ_NAME = 'ZDim'
   CHARACTER (LEN=*), PARAMETER :: DIMLL_NAME = 'nLevels,nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMRL_NAME = 'RDim,nLats'
   CHARACTER (LEN=*), PARAMETER :: DIMXYZ_NAME = 'ZDim,YDim,XDim'
   CHARACTER (LEN=*), PARAMETER :: DIMZYX_NAME = 'XDim,YDim,ZDim'
   CHARACTER (LEN=*), PARAMETER :: PROJ_NAME = 'Projection'
   CHARACTER (LEN=*), PARAMETER :: NAMEPROJ = 'Simple Cylindrical'
                                                                            
   CHARACTER (LEN=*), PARAMETER :: FILEATTR_ERR = 'Failed to write attribute'
   CHARACTER (LEN=*), PARAMETER :: DAT_ERR = 'Failed to define data field '
   CHARACTER (LEN=*), PARAMETER :: GRID_ORIGIN = 'GridOrigin'
   CHARACTER (LEN=*), PARAMETER :: GRID_NAME = 'Center'
   CHARACTER (LEN=*), PARAMETER :: GRID_SPACING = 'GridSpacing'
   CHARACTER (LEN=*), PARAMETER :: GSPACING_VALUE = '(4,2)'
   CHARACTER (LEN=*), PARAMETER :: GRID_SPACING_UNIT = 'GridSpacingUnit'
   CHARACTER (LEN=*), PARAMETER :: GSPACING_UNIT = 'Degree'
   CHARACTER (LEN=*), PARAMETER :: GRID_SPAN = 'GridSpan'
   CHARACTER (LEN=*), PARAMETER :: GSPAN_VALUE = '(0,360,-82,+82)'
   CHARACTER (LEN=*), PARAMETER :: GRID_SPAN_UNIT = 'GridSpanUnit'
   CHARACTER (LEN=*), PARAMETER :: GSPAN_UNIT = 'Degree'
   CHARACTER (LEN=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
   CHARACTER (LEN=*), PARAMETER :: GEO_ERR = & 
        & 'Failed to define geolocation field '
   CHARACTER (LEN=*), PARAMETER :: GD_ERR = 'Failed to detach from grid '
   CHARACTER (LEN=*), PARAMETER :: SW_ERR = 'Failed to detach from swath '
   CHARACTER (LEN=*), PARAMETER :: METAWR_ERR = & 
        & 'Error writing metadata attribute '
   CHARACTER (LEN=*), PARAMETER :: NOOUT_ERR = & 
        & ' data expected but not found for output.'
   CHARACTER (LEN=*), PARAMETER :: TAI2A_ERR = & 
        & 'Error converting time from TAI to UTC.'
   CHARACTER (LEN=*), PARAMETER :: WR_ERR = 'Failed to write field '
   CHARACTER (LEN=*), PARAMETER :: SZ_ERR = 'Failed to get size of dimension '

   INTEGER, PARAMETER :: CCSDS_LEN = 27
   INTEGER, PARAMETER :: CCSDSB_LEN = 25
   INTEGER, PARAMETER :: DATE_LEN = 8
   INTEGER, PARAMETER :: GCTP_GEO = 0
   INTEGER, PARAMETER :: HDFE_NOMERGE = 0
   INTEGER, PARAMETER :: HE5_HDFE_NOMERGE = 0
   INTEGER, PARAMETER :: INVENTORYMETADATA = 2
   INTEGER, PARAMETER :: GridNameLen = 64
   INTEGER, PARAMETER :: maxDgProds = 42
   INTEGER, PARAMETER :: maxGridPoints = 500
   INTEGER, PARAMETER :: maxWindow = 31
   INTEGER, PARAMETER :: maxMisDays = 1000 ! 300
   INTEGER, PARAMETER :: MIN_MAX = 2

! This data type is used to store the names/dates of output files 
! actually created.

   TYPE OutputFiles_T

     CHARACTER (LEN=FileNameLen) :: name(maxwindow)
	! array of names of the created files

     CHARACTER (LEN=8) :: date(maxWindow)       ! CCSDS B format dates of files

     INTEGER :: nFiles		! number of distinct output files created

   END TYPE OutputFiles_T

contains
!=====================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE MLSL3Common
!=====================

!# $Log$
!# Revision 1.17  2006/02/28 17:56:56  cvuu
!# V2.00 commit
!#
!# Revision 1.16  2005/06/23 19:07:39  pwagner
!# Reworded Copyright statement, moved rcs id
!#
!# Revision 1.15  2004/01/07 21:43:18  cvuu
!# version 1.4 commit
!#
!# Revision 1.14  2003/03/22 02:44:16  jdone
!# use only, indentation added
!#
!# Revision 1.13  2001/12/13 20:49:21  nakamura
!# Added field name for MaxDiffTime.
!#
!# Revision 1.12  2001/12/10 18:33:42  nakamura
!# Moved DM dg fields here.
!#
!# Revision 1.11  2001/11/26 19:25:51  nakamura
!# Added some dg fields.
!#
!# Revision 1.10  2001/10/04 18:23:31  nakamura
!# Removed lev as dim for local solar fields.
!#
!# Revision 1.9  2001/09/27 17:47:08  nakamura
!# Added stuff for local solar ancillary fields.
!#
!# Revision 1.8  2001/07/18 15:56:29  nakamura
!# Added a dg parameter; expanded maxWindow for calendar month.
!#
!# Revision 1.7  2001/05/04 18:37:51  nakamura
!# Added date-related params, OutputFiles_T.
!#
!# Revision 1.6  2001/04/24 19:41:13  nakamura
!# Added parameters that L2 wishes to keep PRIVATE.
!#
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
