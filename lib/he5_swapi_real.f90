! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HE5_SWAPI_REAL

  public :: HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, HE5_SWRDFLD_REAL_3D, &
    & HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, HE5_SWWRFLD_REAL_3D, &
    & HE5_EHWRGLATT_REAL, HE5_EHRDGLATT_REAL, &
    & HE5_SWWRATTR_REAL, HE5_SWWRLATTR_REAL, &
    & HE5_SWRDATTR_REAL, HE5_SWRDLATTR_REAL
  public :: HE5_SWSETFILL_REAL
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  integer function HE5_EHWRGLATT_REAL ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_EHWRGLATT

    he5_ehwrglatt_REAL = he5_ehwrglatt(fileID, &
         & attrname, datatype, count, buffer )
  end function HE5_EHWRGLATT_REAL

  integer function HE5_EHRDGLATT_REAL ( FILEID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    real, intent(out) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_EHRDGLATT

    he5_ehrdglatt_REAL = he5_ehrdglatt(fileID, &
         & attrname, buffer )
  end function HE5_EHRDGLATT_REAL

  integer function HE5_SWSETFILL_REAL (SWATHID, FIELDNAME, NUMBERTYPE, &
    & FILLVALUE)
    integer,intent(in)::SWATHID
    character(len=*),intent(IN)::FIELDNAME
    integer, intent(in) :: NUMBERTYPE
    real, intent(in) :: FILLVALUE
    integer, external :: HE5_SWSETFILL
    HE5_SWSETFILL_REAL = &
      & HE5_SWSETFILL (SWATHID, FIELDNAME, NUMBERTYPE, FILLVALUE)
  end function HE5_SWSETFILL_REAL

  integer function HE5_SWRDFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:)      ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_real  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL

  integer function HE5_SWRDFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:)    ! Buffer for read

    integer, external :: HE5_SWRDFLD
    HE5_SWrdfld_real_2d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL_2D

  integer function HE5_SWRDFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_real_3d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_REAL_3D

  integer function HE5_SWWRFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    ! integer(kind=size_t) :: COUNTS(1)
    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL

  integer function HE5_SWWRFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:)     ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real_2d  = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL_2D

  integer function HE5_SWWRFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer(kind=size_t), intent(in) :: STARTS(*)    ! Start array
    integer(kind=size_t), intent(in) :: STRIDES(*)   ! Stride array
    integer(kind=size_t), intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:,:)   ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_real_3d  = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_REAL_3D

  integer function HE5_SWWRATTR_REAL ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_SWWRATTR

    he5_swwrattr_REAL = he5_swwrattr(swathid, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRATTR_REAL

  integer function HE5_SWWRLATTR_REAL ( SWATHID, FIELDNAME, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    use hdf5, only: size_t
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    integer, intent(in) :: DATATYPE
    integer(kind=size_t), intent(in) :: COUNT   ! Stride array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: HE5_SWWRLATTR

    he5_swwrlattr_REAL = he5_swwrlattr(swathid, fieldname, &
         & attrname, datatype, count, buffer )
  end function HE5_SWWRLATTR_REAL

  integer function HE5_SWRDATTR_REAL ( SWATHID, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    real, intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDATTR

    HE5_SWRDATTR_REAL = he5_swrdattr(swathid, &
         & attrname, buffer )
  end function HE5_SWRDATTR_REAL

  integer function HE5_SWRDLATTR_REAL ( SWATHID, FIELDNAME, &
    & ATTRNAME, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Field name
    real, intent(out) :: BUFFER(:)       ! Buffer for read

    integer, external :: HE5_SWRDLATTR

    HE5_SWRDLATTR_REAL = he5_swrdlattr(swathid, fieldname, &
         & attrname, buffer )
  end function HE5_SWRDLATTR_REAL

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HE5_SWAPI_REAL

! $Log$
! Revision 2.10  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
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
