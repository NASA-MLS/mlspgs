
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Temporary_types
!===============================================================================

   USE L1BData
   USE L2GP
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

! Contents:

! ProfileGrid_T
! L2cf
! L2cfEntry
! L2cfCell

! Remarks:  This module contains derived types which are used by L2 but are
!           currently not attached to any module.

! Parameters used by derived types

   INTEGER, PARAMETER :: L2cfEntryLen = 64
   INTEGER, PARAMETER :: MaxKeyLen = 64
   INTEGER, PARAMETER :: MaxCharValueLen = 64
   INTEGER, PARAMETER :: UnitsLen = 64

! This Fortran 90 data type defines a horizontal profile grid upon which the
! retrievals are to be performed.

   TYPE ProfileGrid_T
      CHARACTER (LEN=NAME_LEN) :: name	! Name for this profile grid
      INTEGER :: chunkNo		! Chunk that grid refers to
      INTEGER :: l2gpStart		! Profile # of chunk start in L2GP data
      INTEGER :: noProfs		! # of profiles in the chunk
      INTEGER :: firstNonOverlap	! Index of first non-overlapped profile
      INTEGER :: lastNonOverlap		! Index of last non-overlapped profile

! Geolocation information, dimensioned (noProfs)

      DOUBLE PRECISION, DIMENSION(:), POINTER :: latitude, longitude, time, &
                         solarTime, solarZenith, losAngle, geodAngle

      INTEGER, DIMENSION(:), POINTER :: profMIF

! dimensioned (CCSDS_LEN,noProfs)

      CHARACTER, DIMENSION(:,:), POINTER :: ccsdsTime

   END TYPE ProfileGrid_T

! This datatype defines the l2cf and its sections

   TYPE L2cf

! For each section:

      INTEGER :: NoGlobalSettingsEntries
		! Number of entries within section

      TYPE (L2cfEntry), DIMENSION(:), POINTER :: GlobalSettings
		! The actual entries, dimensioned (NoGlobalSettingsEntries)

      INTEGER :: NoReadAprioriEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: ReadApriori

      INTEGER :: NoMergeAprioriEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: MergeApriori

      INTEGER :: NoChunkDivideEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: ChunkDivide

      INTEGER :: NoProfileLayoutEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: ProfileLayout

      INTEGER :: NoConstructEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: Construct

      INTEGER :: NoFillEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: Fill

      INTEGER :: NoJoinEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: Join

      INTEGER :: NoOutputEntries
      TYPE (L2cfEntry), DIMENSION(:), POINTER :: Output

   END TYPE l2cf

! This datatype defines an L2cf entry, i.e. a name followed by a list of 
! "keyword=value" cells

   TYPE L2cfEntry

! Entry name, e.g. aprioriTemp

      CHARACTER (len=L2cfEntryLen) :: L2cfEntryName

! Number of "keyword=value" cells

      INTEGER :: L2cfEntryNoKeys

! Cell holder, dimensioned (L2cfEntryNoKeys)

      TYPE (L2cfCell), DIMENSION(:), POINTER :: Cells

   END TYPE L2cfEntry

! This datatype defines a L2cf "cell", i.e. a contruction of type:
! "keyword=value" or "keyword=LowerBound..UpperBound[units]"

   TYPE L2cfCell

! First the "keyword" holder:

      CHARACTER (len=MaxKeyLen) :: Keyword

! Value holder for "real" type

      DOUBLE PRECISION :: RealValue

! Value holder for integer type

      INTEGER :: IntValue

! Value holder for "character" type

      CHARACTER (len= MaxCharValueLen) :: CharValue

! Value holder for "units" field

      CHARACTER (len=UnitsLen) :: Units

! Value holder for upper and lower bounds fields

      DOUBLE PRECISION :: RangeLowerBound
      DOUBLE PRECISION :: RangeUpperBound

   END TYPE L2cfCell

!=========================
END MODULE Temporary_types
!=========================

!# $Log$
