! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HE5_SWAPI_CHARACTER_ARRAY

  public :: HE5_EHWRGLATT_CHARACTER_ARRAY, HE5_EHRDGLATT_CHARACTER_ARRAY, &
    & HE5_SWRDFLD_CHARACTER_ARRAY, HE5_SWWRFLD_CHARACTER_ARRAY, &
    & HE5_SWWRATTR_CHARACTER_ARRAY, HE5_SWWRLATTR_CHARACTER_ARRAY, &
    & HE5_SWRDATTR_CHARACTER_ARRAY, HE5_SWRDLATTR_CHARACTER_ARRAY
  private

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  integer function HE5_EHWRGLATT_CHARACTER_ARRAY ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many to write
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_character_ARRAY = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_CHARACTER_ARRAY

  integer function HE5_EHRDGLATT_CHARACTER_ARRAY ( FILEID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_EHRDGLATT

    he5_ehrdglatt_CHARACTER_ARRAY = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_CHARACTER_ARRAY

  integer function HE5_SWRDFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER(:) ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_character_array = HE5_SWrdfld(swathid, fieldname, starts, &
         strides, edges, buffer )
  end function HE5_SWRDFLD_CHARACTER_ARRAY

  integer function HE5_SWWRFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_character_array=HE5_SWwrfld(swathid, fieldname, starts, &
         strides, edges, buffer )
  end function HE5_SWWRFLD_CHARACTER_ARRAY

  integer function HE5_SWWRATTR_CHARACTER_ARRAY ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many to write
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_character_ARRAY = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_CHARACTER_ARRAY

  integer function HE5_SWWRLATTR_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer, intent(in) :: COUNT   ! How many to write
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_character_ARRAY = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_CHARACTER_ARRAY

  integer function HE5_SWRDATTR_CHARACTER_ARRAY ( SWATHID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDATTR
    buffer = ''
    HE5_SWRDATTR_CHARACTER_ARRAY = he5_swrdattr(swathid, &
         & attrname, buffer )
  end function HE5_SWRDATTR_CHARACTER_ARRAY

  integer function HE5_SWRDLATTR_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    character(len=*), intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDLATTR
    buffer = ''
    HE5_SWRDLATTR_CHARACTER_ARRAY = he5_swrdlattr(swathid, fieldname, &
         & attrname, buffer )
  end function HE5_SWRDLATTR_CHARACTER_ARRAY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI_CHARACTER_ARRAY

