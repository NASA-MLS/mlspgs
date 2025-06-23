! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HE5_SWAPI_CHARACTER_ARRAY

  public :: HE5_EHWRGLATT_CHARACTER_ARRAY, HE5_EHRDGLATT_CHARACTER_ARRAY, &
    & HE5_SWRDFLD_CHARACTER_ARRAY, HE5_SWWRFLD_CHARACTER_ARRAY, &
    & HE5_SWWRATTR_CHARACTER_ARRAY, HE5_SWWRLATTR_CHARACTER_ARRAY, &
    & HE5_SWRDATTR_CHARACTER_ARRAY, HE5_SWRDLATTR_CHARACTER_ARRAY
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  integer function HE5_EHWRGLATT_CHARACTER_ARRAY ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many to write
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
    BUFFER = ' ' ! Why is this necessary?

    he5_ehrdglatt_CHARACTER_ARRAY = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_CHARACTER_ARRAY

  integer function HE5_SWRDFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER(:) ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_character_array = HE5_SWrdfld(swathid, fieldname, starts, &
         strides, edges, buffer )
  end function HE5_SWRDFLD_CHARACTER_ARRAY

  integer function HE5_SWWRFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer(kind=size_t) :: COUNTS(1)
    integer, external :: HE5_SWWRFLD

    counts = max(edges(1), int(size(buffer), size_t) )

    HE5_SWwrfld_character_array=HE5_SWwrfld(swathid, fieldname, starts, &
         strides, counts, buffer )
  end function HE5_SWWRFLD_CHARACTER_ARRAY

  integer function HE5_SWWRATTR_CHARACTER_ARRAY ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many to write
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_character_ARRAY = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_CHARACTER_ARRAY

  integer function HE5_SWWRLATTR_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many to write
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HE5_SWAPI_CHARACTER_ARRAY

