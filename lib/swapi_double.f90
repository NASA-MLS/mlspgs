module SWAPI_DOUBLE

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function SWRDFLD_DOUBLE ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(out) :: BUFFER(*)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_double  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_DOUBLE

  integer function SWRDFLD_DOUBLE_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(out) :: BUFFER(1,*)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_double_2d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_DOUBLE_2D

  integer function SWRDFLD_DOUBLE_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(out) :: BUFFER(1,1,*)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_double_3d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_DOUBLE_3D

  integer function SWWRFLD_DOUBLE ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(in) :: BUFFER(*)   ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_double = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_DOUBLE

  integer function SWWRFLD_DOUBLE_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(in) :: BUFFER(1,*)  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_double_2d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_DOUBLE_2D

  integer function SWWRFLD_DOUBLE_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    double precision, intent(in) :: BUFFER(1,1,*)  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_double_3d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_DOUBLE_3D

end module SWAPI_DOUBLE

! $Log$
! Revision 1.1  2000/09/02 01:56:36  vsnyder
! Initial code
!
