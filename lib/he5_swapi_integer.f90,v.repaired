! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HE5_SWAPI_INTEGER

  public :: HE5_SWRDFLD_INTEGER, HE5_SWRDFLD_INTEGER_2D, HE5_SWRDFLD_INTEGER_3D, &
    & HE5_SWWRFLD_INTEGER, HE5_SWWRFLD_INTEGER_2D, HE5_SWWRFLD_INTEGER_3D, &
    & HE5_EHWRGLATT_INTEGER, HE5_EHRDGLATT_INTEGER, &
    & HE5_SWWRATTR_INTEGER, HE5_SWWRLATTR_INTEGER, &
    & HE5_SWRDATTR_INTEGER, HE5_SWRDLATTR_INTEGER
  public :: HE5_SWSETFILL_INTEGER
  public :: HE5_EHRDGLATT
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface HE5_EHRDGLATT
    module procedure HE5_EHRDGLATT_INTEGER_ARR, HE5_EHRDGLATT_INTEGER_SCA
  end interface
contains

  integer function HE5_EHWRGLATT_INTEGER ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    integer, intent(in) :: BUFFER(:)   ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_INTEGER = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_INTEGER

  integer function HE5_EHRDGLATT_INTEGER_ARR ( FILEID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(out) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_EHRDGLATT

    he5_ehrdglatt_INTEGER_ARR = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_INTEGER_ARR

  integer function HE5_EHRDGLATT_INTEGER_SCA ( FILEID, &
    & ATTRNAME, VALUE )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(out) :: VALUE

    integer, external :: HE5_EHRDGLATT
    integer, dimension(1) :: buffer

    he5_ehrdglatt_INTEGER_SCA = he5_ehrdglatt(fileID, &
         & attrname, buffer )
    value = buffer(1)
  end function HE5_EHRDGLATT_INTEGER_SCA

  integer function HE5_SWSETFILL_INTEGER (SWATHID, FIELDNAME, NUMBERTYPE, &
    & FILLVALUE)
    integer,intent(in)::SWATHID
    character(len=*),intent(IN)::FIELDNAME
    integer, intent(in) :: NUMBERTYPE
    integer, intent(in) :: FILLVALUE
    integer, external :: HE5_SWSETFILL
    HE5_SWSETFILL_INTEGER = &
      & HE5_SWSETFILL (SWATHID, FIELDNAME, NUMBERTYPE, FILLVALUE)
  end function HE5_SWSETFILL_INTEGER

  integer function HE5_SWRDFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER

  integer function HE5_SWRDFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer_2d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER_2D

  integer function HE5_SWRDFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:,:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer_3d =HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER_3D

  integer function HE5_SWWRFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:)    ! Buffer for write

    integer(kind=size_t) :: COUNTS(1)
    integer, external :: HE5_SWWRFLD

    counts = max(edges(1), int(size(buffer), size_t) )

    HE5_SWwrfld_integer = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & counts, buffer )
  end function HE5_SWWRFLD_INTEGER

  integer function HE5_SWWRFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:,:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_integer_2d=HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_INTEGER_2D

  integer function HE5_SWWRFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:,:,:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_integer_3d =HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_INTEGER_3D

  integer function HE5_SWWRATTR_INTEGER ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    integer, intent(in) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_INTEGER = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_INTEGER

  integer function HE5_SWWRLATTR_INTEGER ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
  use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    integer, intent(in) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_INTEGER = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_INTEGER

  integer function HE5_SWRDATTR_INTEGER ( SWATHID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDATTR

    HE5_SWRDATTR_INTEGER = he5_swrdattr(swathid, &
         & attrname, buffer )
  end function HE5_SWRDATTR_INTEGER

  integer function HE5_SWRDLATTR_INTEGER ( SWATHID, FIELDNAME, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDLATTR

    HE5_SWRDLATTR_INTEGER = he5_swrdlattr(swathid, fieldname, &
         & attrname, buffer )
  end function HE5_SWRDLATTR_INTEGER

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HE5_SWAPI_INTEGER

! $Log$
! Revision 2.11  2010/03/25 18:41:46  pwagner
! args to max function now the same integer type
!
! Revision 2.10  2010/01/14 23:18:37  pwagner
! Separate scalar and array generic forms for HE5_EHRDGLATT_INTEGER
!
! Revision 2.9  2009/10/05 23:37:05  pwagner
! Moved use hdf5 statements from module scope to speedup Lahey; this is the last time we do taht
!
! Revision 2.8  2009/09/29 23:34:38  pwagner
! Changes needed by 64-bit build
!
! Revision 2.7  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2004/02/13 00:16:02  pwagner
! New stuff for reading swath attributes
!
! Revision 2.4  2003/04/11 23:32:23  pwagner
! Moved he5_swsetfill he5_ehwrglatt interfaces
!
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
