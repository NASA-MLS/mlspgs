module HE5_SWAPI_INTEGER

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function HE5_SWRDFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:)   ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER

  integer function HE5_SWRDFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer_2d  = HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER_2D

  integer function HE5_SWRDFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(out) :: BUFFER(:,:,:)  ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_integer_3d =HE5_SWrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWRDFLD_INTEGER_3D

  integer function HE5_SWWRFLD_INTEGER ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:)    ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_integer = HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_INTEGER

  integer function HE5_SWWRFLD_INTEGER_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:,:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_integer_2d=HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_INTEGER_2D

  integer function HE5_SWWRFLD_INTEGER_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    integer, intent(in) :: BUFFER(:,:,:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_integer_3d =HE5_SWwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function HE5_SWWRFLD_INTEGER_3D

end module HE5_SWAPI_INTEGER

