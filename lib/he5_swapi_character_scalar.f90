! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HE5_SWAPI_CHARACTER_SCALAR

  public

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
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_character_scalar = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_CHARACTER_SCALAR

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
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_character_scalar = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_CHARACTER_SCALAR

  integer function HE5_SWWRLATTR_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_character_scalar = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_CHARACTER_SCALAR

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI_CHARACTER_SCALAR

! $Log$
! Revision 2.4  2003/04/11 23:32:23  pwagner
! Moved he5_swsetfill he5_ehwrglatt interfaces
!
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
