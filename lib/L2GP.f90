
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L2GP
!===============================================================================

   USE L1BData
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

! Contents:

! Definition -- L2GP_T
! Subroutines -- L2GP_createFile
!                L2GP_writeGeo
!                DeallocateL2GP

! Remarks:  This is a prototype module containing parameters, a derived type
!           definition, and subroutines used in producing an L2GP output file.

! Parameters

   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD1 = 'latitude'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD2 = 'longitude'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD3 = 'time'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD4 = 'ccsdsTime'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD5 = 'solarTime'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD6 = 'solarZenith'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD7 = 'losAngle'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD8 = 'geodAngle'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD9 = 'chunkNumber'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD10 = 'pressure'
   CHARACTER (LEN=*), PARAMETER :: GEO_FIELD11 = 'frequency'

   INTEGER, PARAMETER :: CCSDS_LEN = 27		! len of CCSDS time string
   INTEGER, PARAMETER :: HDFE_AUTOMERGE = 1	! merge fields with share dim
   INTEGER, PARAMETER :: HDFE_NOMERGE = 0	! don't merge

! This data type is used to store all or part of an L2GP quantity prior to
! output.

   TYPE L2GP_T

      CHARACTER (LEN=NAME_LEN) :: name  ! Name for quantity to be output

      INTEGER :: noProfs	! Total number of profiles
      INTEGER :: noSurfs	! Total number of surfaces
      INTEGER :: noFreqs	! Number of frequencies in breakdown

      INTEGER :: noProfsStored		! # of profiles stored in this variable
      INTEGER :: firstProfStored	! First profile in this variable

! Now we store the geolocation fields, first the vertical one:

      REAL, DIMENSION(:), POINTER :: pressure
					! vertical coordinates (noSurfs)

! Now the horizontal geolocation information,
! dimensioned (firstProfStored:firstProfStored+noProfsStored-1)

      DOUBLE PRECISION, DIMENSION(:), POINTER :: time
      INTEGER, DIMENSION(:), POINTER :: chunkNumber
      REAL, DIMENSION(:), POINTER :: latitude, longitude, solarTime, &
                                     solarZenith, losAngle, geodAngle

! dimensioned (CCSDS_LEN,firstProfStored:firstProfStored+noProfsStored-1)

      CHARACTER, DIMENSION(:,:), POINTER :: ccsdsTime

! Now we store the "frequency" geolocation field

      DOUBLE PRECISION, DIMENSION(:), POINTER :: freq
							! dimensioned (noFreqs)
! Finally we store the data fields

      REAL, DIMENSION(:,:,:), POINTER :: l2gpValue
      REAL, DIMENSION(:,:,:), POINTER :: l2gpPrecision
		! dimensioned (noFreqs, noSurfs,
		!	firstProfStored:firstProfStored+noProfsStored-1)

      CHARACTER (LEN=1), DIMENSION(:), POINTER :: l2gpStatus
					! ( status is a reserved word in F90)
      REAL, DIMENSION(:), POINTER :: quality
		! both of the above are dimensioned:
		! 	(firstProfStored:firstProfStored+noProfsStored-1)
   END TYPE L2GP_T

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE L2GP_createFile (L2FileHandle, name, nProf, nSurf, nFreq, flag)
!----------------------------------------------------------------------------

! Brief description of program
! This program contains prototype code for setting up the structural
! definitions in an empty L2GP file.

! Arguments

      CHARACTER (LEN=NAME_LEN), INTENT(IN) :: name

      INTEGER, INTENT(IN) :: L2FileHandle, nProf, nSurf, nFreq

      INTEGER, INTENT(OUT) :: flag

! Parameters

      CHARACTER (LEN=*), PARAMETER :: DATA_FIELD1 = 'l2gpValue'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELD2 = 'l2gpPrecision'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELD3 = 'l2gpStatus'
      CHARACTER (LEN=*), PARAMETER :: DATA_FIELD4 = 'quality'

      CHARACTER (LEN=*), PARAMETER :: DIM_NAME1 = 'noProfs'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME2 = 'noSurfs'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME3 = 'noFreqs'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME4 = 'CCSDSLen'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME12 = 'noProfs,noSurfs'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME123 = 'noProfs,noSurfs,noFreqs'
      CHARACTER (LEN=*), PARAMETER :: DIM_NAME41 = 'CCSDSLen,noProfs'

      CHARACTER (LEN=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
      CHARACTER (LEN=*), PARAMETER :: GEO_ERR = 'Failed to define geolocation &
                                                &field '
      CHARACTER (LEN=*), PARAMETER :: DAT_ERR = 'Failed to define data field '  

! Functions

      INTEGER :: Pgs_smf_generateStatusReport
      INTEGER :: swcreate, swdefdfld, swdefdim, swdefgfld, swdetach

! Variables

      CHARACTER (LEN=480) :: msr
   
      INTEGER :: i, returnStatus, swid, status

      flag = 0

! Create the swath within the file

      swid = swcreate(L2FileHandle, name)
      IF (swid == -1) THEN
         msr = 'Failed to create swath ' // name
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

! Define dimensions

      status = swdefdim(swid, DIM_NAME1, nProf)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME1
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      IF (nSurf > 0) THEN
         status = swdefdim(swid, DIM_NAME2, nSurf)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME2
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF
      ENDIF

      IF (nFreq > 0) THEN
         status = swdefdim(swid, DIM_NAME3, nFreq)
         IF (status == -1) THEN
            msr = DIM_ERR // DIM_NAME3
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF
      ENDIF

      status = swdefdim(swid, DIM_NAME4, CCSDS_LEN)
      IF (status == -1) THEN
         msr = DIM_ERR // DIM_NAME4
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

! Define horizontal geolocation fields using above dimensions

      status = swdefgfld(swid, GEO_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD1
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD2
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD3, DIM_NAME1, DFNT_FLOAT64, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD3
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD4, DIM_NAME41, DFNT_CHAR8, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD4
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD5
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD6, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD6
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD7, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD7
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD8, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD8
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefgfld(swid, GEO_FIELD9, DIM_NAME1, DFNT_INT32, HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = GEO_ERR // GEO_FIELD9
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      IF (nSurf > 0) THEN
         status = swdefgfld(swid, GEO_FIELD10, DIM_NAME2, DFNT_FLOAT32, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD10
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF
      ENDIF

      IF (nFreq > 0) THEN
         status = swdefgfld(swid, GEO_FIELD11, DIM_NAME3, DFNT_FLOAT64, &
                            HDFE_NOMERGE)
         IF (status == -1) THEN
            msr = GEO_ERR // GEO_FIELD11
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF
      ENDIF

! Define data fields using above dimensions

      IF (nFreq > 0) THEN

         status = swdefdfld(swid, DATA_FIELD1, DIM_NAME123, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         status = swdefdfld(swid, DATA_FIELD2, DIM_NAME123, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD1 // ' or ' // DATA_FIELD2 // ' for 3D &
                 &quantity.'
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF

      ELSE IF (nSurf > 0) THEN

         status = swdefdfld(swid, DATA_FIELD1, DIM_NAME12, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         status = swdefdfld(swid, DATA_FIELD2, DIM_NAME12, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD1 // ' or ' // DATA_FIELD2 // ' for 2D &
                 &quantity.'
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF

      ELSE

         status = swdefdfld(swid, DATA_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         status = swdefdfld(swid, DATA_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
                            HDFE_NOMERGE)

         IF (status == -1) THEN
            msr = DAT_ERR // DATA_FIELD1 // ' or ' // DATA_FIELD2 // ' for 1D &
                 &quantity.'
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF

      ENDIF

      status = swdefdfld(swid, DATA_FIELD3, DIM_NAME1, DFNT_CHAR8, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELD3
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swdefdfld(swid, DATA_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
                         HDFE_NOMERGE)
      IF (status == -1) THEN
         msr = DAT_ERR // DATA_FIELD4
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

! Detach from the swath interface.  This stores the swath info within the file
! and must be done before writing or reading data to or from the swath.

      status = swdetach(swid)
      IF (status == -1) THEN
         returnStatus = Pgs_smf_generateStatusReport('Failed to detach from &
                                         &swath interface after definition.')
         flag = -1

      ENDIF

!--------------------------------
   END SUBROUTINE L2GP_createFile
!--------------------------------

!------------------------------------------------
   SUBROUTINE L2GP_writeGeo (swid, l2gpGeo, flag)
!------------------------------------------------

! Brief description of subroutine
! This program contains prototype code for writing the geolocation fields to an
! L2GP output file.

! Arguments

      TYPE( L2GP_T ), INTENT(IN) :: l2gpGeo

      INTEGER, INTENT(IN) :: swid

      INTEGER, INTENT(OUT) :: flag

! Parameters

      CHARACTER (LEN=*), PARAMETER :: WR_ERR = 'Failed to write geolocation &
                                               &field '

! Functions

      INTEGER :: Pgs_smf_generateStatusReport, swwrfld

! Variables

      CHARACTER (LEN=480) :: msr
   
      INTEGER :: i, returnStatus, status
      INTEGER :: start(2), stride(2), edge(2)

      flag = 0

! Write data to the fields

      stride = 1
      start(1) = 0
      start(2) = l2gpGeo%firstProfStored
      edge(1) = CCSDS_LEN
      edge(2) = l2gpGeo%noProfsStored

      status = swwrfld(swid, GEO_FIELD1, start(2), stride(2), edge(2), &
                       l2gpGeo%latitude)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD1
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD2, start(2), stride(2), edge(2), &
                       l2gpGeo%longitude)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD2
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD3, start(2), stride(2), edge(2), &
                       l2gpGeo%time)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD3
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld( swid, GEO_FIELD4, start, stride, edge, &
                        l2gpGeo%ccsdsTime)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD4
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD5, start(2), stride(2), edge(2), &
                       l2gpGeo%solarTime)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD5
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD6, start(2), stride(2), edge(2), &
                       l2gpGeo%solarZenith)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD6
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD7, start(2), stride(2), edge(2), &
                       l2gpGeo%losAngle)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD7
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF
 
      status = swwrfld(swid, GEO_FIELD8, start(2), stride(2), edge(2), &
                       l2gpGeo%geodAngle)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD8
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      status = swwrfld(swid, GEO_FIELD9, start(2), stride(2), edge(2), &
                       l2gpGeo%chunkNumber)
      IF (status == -1) THEN
         msr = WR_ERR // GEO_FIELD9
         returnStatus = Pgs_smf_generateStatusReport(msr)
         flag = -1
         RETURN
      ENDIF

      IF (l2gpGeo%noSurfs > 0) THEN

         start(1) = 0
         edge(1) = l2gpGeo%noSurfs

         status = swwrfld(swid, GEO_FIELD10, start(1), stride(1), edge(1), &
                          l2gpGeo%pressure)
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD10
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF

      ENDIF

      IF (l2gpGeo%noFreqs > 0) THEN

         start(1) = 0
         edge(1) = l2gpGeo%noFreqs

         status = swwrfld(swid, GEO_FIELD11, start(1), stride(1), edge(1), &
                          l2gpGeo%freq)
         IF (status == -1) THEN
            msr = WR_ERR // GEO_FIELD11
            returnStatus = Pgs_smf_generateStatusReport(msr)
            flag = -1
            RETURN
         ENDIF

      ENDIF

!------------------------------
   END SUBROUTINE L2GP_writeGeo
!------------------------------

!----------------------------------------
   SUBROUTINE DeallocateL2GP (l2gp, flag)
!----------------------------------------

! Brief description of subroutine
! This subroutine deallocates the internal field pointers of the L2GP_T
! derived type, after the calling program has finished with the data.

! Arguments

      TYPE( L2GP_T ), INTENT(IN) :: l2gp

      INTEGER, INTENT(OUT) :: flag

! Parameters

! Functions

      INTEGER :: Pgs_smf_generateStatusReport

! Variables

      INTEGER :: dealloc_err, returnStatus

      flag = 0

! Horizontal geolocation fields

      IF ( ASSOCIATED(l2gp%time) ) DEALLOCATE (l2gp%time, STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%chunkNumber) ) DEALLOCATE (l2gp%chunkNumber, &
                                                      STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%latitude) ) DEALLOCATE (l2gp%latitude, &
                                                   STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%longitude) ) DEALLOCATE (l2gp%longitude, &
                                                    STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%solarTime) ) DEALLOCATE (l2gp%solarTime, &
                                                    STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%solarZenith) ) DEALLOCATE (l2gp%solarZenith, &
                                                      STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%losAngle) ) DEALLOCATE (l2gp%losAngle, &
                                                   STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%geodAngle) ) DEALLOCATE (l2gp%geodAngle, &
                                                    STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%ccsdsTime) ) DEALLOCATE (l2gp%ccsdsTime, &
                                                    STAT=dealloc_err)

! Vertical geolocation field

      IF ( ASSOCIATED(l2gp%pressure) ) DEALLOCATE (l2gp%pressure, &
                                                   STAT=dealloc_err)
! Frequency "geolocation field"

      IF ( ASSOCIATED(l2gp%freq) ) DEALLOCATE (l2gp%freq, STAT=dealloc_err)

! Data fields

      IF ( ASSOCIATED(l2gp%l2gpValue) ) DEALLOCATE (l2gp%l2gpValue, &
                                                    STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%l2gpPrecision) ) DEALLOCATE (l2gp%l2gpPrecision, &
                                                        STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%l2gpStatus) ) DEALLOCATE (l2gp%l2gpStatus, &
                                                     STAT=dealloc_err)
      IF ( ASSOCIATED(l2gp%quality) ) DEALLOCATE (l2gp%quality, &
                                                  STAT=dealloc_err)

! Error checking

      IF ( dealloc_err /= 0 ) THEN
         returnStatus = Pgs_smf_generateStatusReport('Failed deallocation of &
                                                     &L2GP pointers.')
         flag = -1
      ENDIF

!-------------------------------
   END SUBROUTINE DeallocateL2GP
!-------------------------------

!==============
END MODULE L2GP
!==============

!# $Log$
!# Revision 1.3  1999/11/18 16:02:36  nakamura
!# Changed name of first subroutine to L2GP_createFile; added subroutines L2GP_writeGeo and DeallocateL2GP.
!#
!# Revision 1.2  1999/11/15 19:09:54  nakamura
!# Changed dimensions to start from 1, instead of 0.
!#
!# Revision 1.1  1999/11/05 19:52:41  nakamura
!# Initial revision
!#
