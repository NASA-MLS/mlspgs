! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HE5_SWAPI_CHARACTER_SCALAR

  public :: HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHRDGLATT_CHARACTER_SCALAR, &
    & HE5_SWRDFLD_CHARACTER_SCALAR, HE5_SWWRFLD_CHARACTER_SCALAR, &
    & HE5_SWWRATTR_CHARACTER_SCALAR, HE5_SWWRLATTR_CHARACTER_SCALAR, &
    & HE5_SWRDATTR_CHARACTER_SCALAR, HE5_SWRDLATTR_CHARACTER_SCALAR
  private

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  integer function HE5_EHWRGLATT_CHARACTER_SCALAR ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_character_scalar = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_CHARACTER_SCALAR

  integer function HE5_EHRDGLATT_CHARACTER_SCALAR ( FILEID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER   ! Buffer for read

    integer, external :: HE5_EHRDGLATT

    he5_ehrdglatt_CHARACTER_SCALAR = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_CHARACTER_SCALAR

  integer function HE5_SWRDFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_character_scalar = HE5_SWrdfld(swathid, fieldname, starts, &
         strides,edges, buffer )
  end function HE5_SWRDFLD_CHARACTER_SCALAR

  integer function HE5_SWWRFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_character_scalar = HE5_SWwrfld(swathid, fieldname, starts, &
         strides,edges, buffer )
  end function HE5_SWWRFLD_CHARACTER_SCALAR

  integer function HE5_SWWRATTR_CHARACTER_SCALAR ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_character_scalar = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_CHARACTER_SCALAR

  integer function HE5_SWWRLATTR_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_character_scalar = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_CHARACTER_SCALAR

  integer function HE5_SWRDATTR_CHARACTER_SCALAR ( SWATHID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER       ! Buffer for read

    integer, external :: HE5_SWRDATTR
    buffer = ''
    HE5_SWRDATTR_CHARACTER_SCALAR = he5_swrdattr(swathid, &
         & attrname, buffer )
  end function HE5_SWRDATTR_CHARACTER_SCALAR

  integer function HE5_SWRDLATTR_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER       ! Buffer for read

    integer, external :: HE5_SWRDLATTR
    buffer = ''
    HE5_SWRDLATTR_CHARACTER_SCALAR = he5_swrdlattr(swathid, fieldname, &
         & attrname, buffer )
  end function HE5_SWRDLATTR_CHARACTER_SCALAR

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI_CHARACTER_SCALAR

! $Log$
! Revision 2.8  2004/03/24 23:53:02  pwagner
! Switched from HE5T_NATIVE_SCHAR to MLS_CHARTYPE
!
! Revision 2.7  2004/02/26 21:57:25  pwagner
! Takes care that attribute return values nulled out before being read
!
! Revision 2.6  2004/02/13 00:16:02  pwagner
! New stuff for reading swath attributes
!
! Revision 2.5  2003/10/28 00:27:46  pwagner
! Corrected some comments
!
! Revision 2.4  2003/04/11 23:32:23  pwagner
! Moved he5_swsetfill he5_ehwrglatt interfaces
!
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
