! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! --------------------------------------------------
module l2cf         ! global declarations
! --------------------------------------------------

type l2cfKey
  character (len=MaxKeyLen) :: Keyword
  character (len=MaxTypeLen):: Type
  integer :: Keylen
end type l2cfKey

type (l2cfKey), dimension(MaxNoL2cfKeys) :: L2cfTable




TYPE L2cfCell

   ! This datatype defines a L2cf "cell", i.e. a contruction of type:
   ! "keyword=value" or "keyword=LowerBound..UpperBound[units]"

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








TYPE L2cfEntry

   ! This datatype defines am L2cf entry, i.e. a name followed by
   ! a list of "keyword=value" cells
 
   ! Entry name, e.g. aprioriTemp

   CHARACTER (len=L2cfEntryLen) L2cfEntryName

   ! Number of "keyword=value" cells

   INTEGER L2cfEntryNoKeys

   !Cells holder 

   TYPE (L2cfCell), DIMENSION(MaxNoKeysPerEntry) :: Cells

END TYPE L2cfEntry


TYPE L2cfSection


   ! For each section:
   ! Name, Number of entries within section
   CHARACTER (len=L2cfEntryLen) L2cfSectionName
   INTEGER NoSectionEntries

   ! The actual entries:

   TYPE (L2cfEntry), DIMENSION(MaxNoEntriesPerSection) :: Entries ! (NoSectionEntries)


END TYPE L2cfSection




TYPE L2cf

! This datatype defines the l2cf as an array of L2cfSections

  INTEGER NoSections

  TYPE (L2cfSection), POINTER, DIMENSION (:) :: Sections ! (NoSections) 

END TYPE l2cf

end module l2cf

