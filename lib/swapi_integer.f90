module SWAPI_INTEGER

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function SWRDFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(*)   ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_integer  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_INTEGER

  integer function SWRDFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(1,*)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_integer_2d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_INTEGER_2D

  integer function SWRDFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(1,1,*)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_integer_3d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_INTEGER_3D

  integer function SWWRFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(*)    ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_integer = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_INTEGER

  integer function SWWRFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(1,*)  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_integer_2d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_INTEGER_2D

  integer function SWWRFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(1,1,*)  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_integer_3d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_INTEGER_3D

end module SWAPI_INTEGER

! $Log$
