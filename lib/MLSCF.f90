! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSCF                    ! MLSCF stuff excluding reading and parsing
!=============================================================================

  ! This is a grossly stripped-down instance of the version 0.1 MLSCF
  ! module.  All it does is provide some parameters, some types and an
  ! allocator. It is only used in the test harness.  It is not used
  ! in the revised MLSCF parser.

  implicit NONE
  public

  !------------------- RCS Ident Info -----------------------
  character(len=130), private :: Id = &                                                    
    & "$Id$"
  character (len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !----------------------------------------------------------

  integer, parameter :: MaxKeyLen = 25
  integer, parameter :: MaxTypeLen = 6
  integer, parameter :: MaxNOMlscfKeys = 20
  integer, parameter :: MaxCharValueLen = 132
  integer, parameter :: ShortCharValueLen = 30 ! for upper end of string range
  integer, parameter :: UnitsLen = 10
  integer, parameter :: MlscfEntryLen = 25
  integer, parameter :: MaxNoKeysPerEntry = 60
  integer, parameter :: MaxNoDefsPerSection = 100
  integer, parameter :: MaxNoEntriesPerSection = 100

  type MlscfKey_T
    character (len=MaxKeyLen) :: Keyword
    character (len=MaxTypeLen):: TYPE
    integer :: Keylen
  end type MlscfKey_T


  type MlscfCell_T

    ! This datatype defines a Mlscf "cell", i.e. a contruction of type:
    ! "keyword=value" or "keyword=LowerBound..UpperBound[units]"

    ! First the "keyword" holder:

    character (len=MaxKeyLen) :: Keyword

    ! Type: Label of a specification, String, Real Value, Range, etc.
    ! See the DECLARATION_TABLE module.
    integer :: TYPE

    ! Value holder for "real" type

    double precision :: RealValue

    ! Value holder for "character" type

    character (len= MaxCharValueLen) :: CharValue

    ! Value holder for "units" field.  See the UNITS module.

    integer :: Units

    ! Value holder for numerical upper bound.  Lower bound is in RealValue.

    double precision :: RangeUpperBound

    ! Value holder for character upper bound.  Lower bound is CharValue

    character ( len= ShortCharValueLen) :: CharRangeUpperBound

    ! How many more CELLS are part of the same datum.  The following
    ! ones have a blank keyword.

    integer :: More

  end type MlscfCell_T


  type MlscfEntry_T

    ! This datatype defines am Mlscf entry, i.e. a name followed by
    ! a list of "keyword=value" cells

    ! Entry label, i.e. the name before the colon before the entry,
    ! if any.

    character (len=MlscfEntryLen) :: MlscfLabelName

    ! Entry name, e.g. aprioriTemp

    character (len=MlscfEntryLen) :: MlscfEntryName

    ! Number of "keyword=value" cells

    integer :: MlscfEntryNoKeys

    !Cells holder 

    type (MlscfCell_T), dimension(MaxNoKeysPerEntry) :: Cells

  end type MlscfEntry_T


  type MlscfSection_T

    ! For each section:
    ! Name, Number of "X=..." definitions, the definitions,
    ! Number Y,... entries, the entries.
    character (len=MlscfEntryLen) MlscfSectionName
    integer NoSectionDefs, NoSectionEntries

    ! The definitions:

    type (MlscfCell_T), dimension(MaxNoDefsPerSection) :: Cells

    ! The entries:

    type (MlscfEntry_T), dimension(MaxNoEntriesPerSection) :: Entries ! (NoSectionEntries)

  end type MlscfSection_T


  type Mlscf_T

    ! This datatype defines the mlscf as an array of MlscfSections

    integer :: NoSections

    type (MlscfSection_T), POINTER, DIMENSION (:) :: Sections ! (NoSections) 

  end type Mlscf_T



! =====     Public Procedures     ======================================

contains

! -----------------------------------------------  ALLOCATE_MLSCF  -----
  subroutine ALLOCATE_MLSCF ( MLSCF_DATA, HOW_MANY_SECTIONS )
  ! Allocate the "mlscf_data % sections" field with "how_many_sections"
  ! elements

    type(mlscf_t), intent(out) :: MLSCF_DATA
    integer, intent(in) :: HOW_MANY_SECTIONS
    integer :: STATUS
    allocate ( mlscf_data % sections ( how_many_sections ), stat = status )
    mlscf_data%NoSections = how_many_sections
    if ( status == 0 ) return
!   call MLSMessage ( MLSMSG_Error, ModuleName, &
!     & 'Unable to allocate "sections" field of "mlscf_data"' )
    stop
  end subroutine ALLOCATE_MLSCF
!=======================
end module MLSCF
!=======================

! $Log$
! Revision 2.3  2001/02/22 01:54:41  vsnyder
! Periodic commit
!
! Revision 2.2  2000/10/19 22:05:02  vsnyder
! Try to get the new one to stick
!
! Revision 2.0  2000/09/05 17:41:49  dcuddy
! Change revision to 2.0
!
! Revision 1.2  2000/09/01 21:42:15  vsnyder
! Add "more" field to MLSCF_CELL type
!
