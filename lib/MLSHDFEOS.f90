! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSHDFEOS

  ! This module contains MLS specific routines to do common HDFEOS
  ! tasks not specifically found in he5_swapi or he5_gdapi.

  use MLSStrings, only: GetStringElement, NumStringElements
  use HE5_SWAPI, only: HE5_SWSETFILL, HE5_SWWRATTR, HE5_SWWRLATTR
  use HE5_SWAPI_CHARACTER_ARRAY, only: HE5_EHWRGLATT_CHARACTER_ARRAY
  use HE5_SWAPI_CHARACTER_SCALAR, only: HE5_EHWRGLATT_CHARACTER_SCALAR
  use HE5_SWAPI_DOUBLE, only: HE5_EHWRGLATT_DOUBLE
  use HE5_SWAPI_INTEGER, only: HE5_EHWRGLATT_INTEGER
  use HE5_SWAPI_REAL, only: HE5_EHWRGLATT_REAL

  implicit NONE
  private

  public :: HE5_EHWRGLATT, MLS_SWSETFILL
  logical, parameter :: HE5_SWSETFILL_BROKEN = .true.
  character(len=*), parameter :: SETFILLTITLE = '_FillValue'

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! HE5_EHWRGLATT     Sets the global attributes at file level
! MLS_SWSETFILL     Sets the fill value for one or more fields
! === (end of toc) ===

! === (start of api) ===
! int HE5_EHWRGLATT (int fileID, char* attrName, int datatype, int count, value)
!     value can be one of:
!    {char* value, char* value(:), int value(:), r4 value(:), r8 value(:)}
! int MLS_SWSETFILL (int swathID, char* names, int datatype, value) 
!     value can be one of:
!    {char* value, int value, r4 value, r8 value}
! === (end of api) ===
  interface HE5_EHWRGLATT   ! From its name, might better be in he5_ehapi.f90
    module procedure HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHWRGLATT_DOUBLE, &
    HE5_EHWRGLATT_INTEGER, HE5_EHWRGLATT_REAL, HE5_EHWRGLATT_CHARACTER_ARRAY
  end interface

  interface MLS_SWSETFILL
    module procedure MLS_SWSETFILL_DOUBLE, &
      & MLS_SWSETFILL_INTEGER, MLS_SWSETFILL_REAL
  end interface

contains ! ======================= Public Procedures =========================

  integer function MLS_SWSETFILL_DOUBLE ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    double precision, intent(in) :: FILLVALUE

    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    mls_swsetfill_double = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields < 1 ) return
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_double = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_double = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_double == -1 ) return
    enddo

  end function MLS_SWSETFILL_DOUBLE

  integer function MLS_SWSETFILL_INTEGER ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    integer, intent(in) :: FILLVALUE

    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    mls_swsetfill_integer = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields < 1 ) return
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_integer = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_integer = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_integer == -1 ) return
    enddo

  end function MLS_SWSETFILL_INTEGER

  integer function MLS_SWSETFILL_REAL ( SWATHID, FIELDNAMES, DATATYPE, &
    & FILLVALUE )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAMES     ! Field names
    integer, intent(in) :: DATATYPE
    real, intent(in) :: FILLVALUE
    integer, external :: HE5_SWsetfill
    integer :: Field
    integer :: numFields
    character(len=len(FIELDNAMES)) :: FIELDNAME
    print *, 'SWATHID: ', SWATHID
    print *, 'FIELDNAMES: ', FIELDNAMES
    mls_swsetfill_real = 0
    numFields = NumStringElements(fieldnames, .false.)
    if ( numFields == -1 ) return
    print *, 'numFields: ', numFields
    do Field=1, numFields
      call GetStringElement(fieldnames, fieldname, Field, .true.)
      print *, 'trim(fieldname): ', trim(fieldname)
      if ( HE5_SWSETFILL_BROKEN ) then
        mls_swsetfill_real = he5_swwrlattr(swathid, trim(fieldname), &
          & trim(SETFILLTITLE), &
          & DATATYPE, 1, (/ FILLVALUE /) )
      else
        mls_swsetfill_real = HE5_SWsetfill(swathid, trim(fieldname), &
        & datatype, fillvalue)
      endif
      if ( mls_swsetfill_real == -1 ) return
    enddo

  end function MLS_SWSETFILL_REAL


! ======================= Private Procedures =========================  

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSHDFEOS

! $Log$
! Revision 2.2  2003/04/15 21:58:54  pwagner
! Now sets _FillValue attribute because swsetfill seems broken
!
! Revision 2.1  2003/04/11 23:28:44  pwagner
! FIrst commit
!
