! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HE5_SWAPI_REAL

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  integer function HE5_SWRDFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:)      ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_real  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL

  integer function HE5_SWRDFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:)    ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_real_2d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL_2D

  integer function HE5_SWRDFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_real_3d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL_3D

  integer function HE5_SWWRFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL

  integer function HE5_SWWRFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:)     ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real_2d  = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL_2D

  integer function HE5_SWWRFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:,:)   ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real_3d  = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL_3D

  integer function HE5_EHWRGLATT_REAL ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_REAL = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_REAL

  integer function HE5_SWWRATTR_REAL ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_REAL = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_REAL

  integer function HE5_SWWRLATTR_REAL ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE    ! Start array
    integer, intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_REAL = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_REAL

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI_REAL

! $Log$
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
