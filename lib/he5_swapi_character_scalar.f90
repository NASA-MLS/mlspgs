! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HE5_SWAPI_CHARACTER_SCALAR

  public :: HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHRDGLATT_CHARACTER_SCALAR, &
    & HE5_SWRDFLD_CHARACTER_SCALAR, HE5_SWWRFLD_CHARACTER_SCALAR, &
    & HE5_SWWRATTR_CHARACTER_SCALAR, HE5_SWWRLATTR_CHARACTER_SCALAR, &
    & HE5_SWRDATTR_CHARACTER_SCALAR, HE5_SWRDLATTR_CHARACTER_SCALAR
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  integer function HE5_EHWRGLATT_CHARACTER_SCALAR ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many
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
    BUFFER = ' ' ! Why is this necessary?

    he5_ehrdglatt_CHARACTER_SCALAR = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_CHARACTER_SCALAR

  integer function HE5_SWRDFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_character_scalar = HE5_SWrdfld(swathid, fieldname, starts, &
         strides,edges, buffer )
  end function HE5_SWRDFLD_CHARACTER_SCALAR

  integer function HE5_SWWRFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_character_scalar = HE5_SWwrfld(swathid, fieldname, starts, &
         strides,edges, buffer )
  end function HE5_SWWRFLD_CHARACTER_SCALAR

  integer function HE5_SWWRATTR_CHARACTER_SCALAR ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_character_scalar = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_CHARACTER_SCALAR

  integer function HE5_SWWRLATTR_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_charType
    integer(kind=size_t), intent(in) :: COUNT   ! How many
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HE5_SWAPI_CHARACTER_SCALAR

! $Log$
! Revision 2.13  2022/09/30 20:48:08  pwagner
! Fixed strange bug in HE5_EHRDGLATT_CHARACTER_SCALAR
!
! Revision 2.12  2009/10/05 23:37:05  pwagner
! Moved use hdf5 statements from module scope to speedup Lahey; this is the last time we do taht
!
! Revision 2.11  2009/09/29 23:34:38  pwagner
! Changes needed by 64-bit build
!
! Revision 2.10  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
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
