module SWAPI_REAL

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function SWRDFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:)      ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_real  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWRDFLD_REAL

  integer function SWRDFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:)    ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_real_2d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWRDFLD_REAL_2D

  integer function SWRDFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(out) :: BUFFER(:,:,:)  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_real_3d  = swrdfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWRDFLD_REAL_3D

  integer function SWWRFLD_REAL ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:)       ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_real = swwrfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWWRFLD_REAL

  integer function SWWRFLD_REAL_2D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:)     ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_real_2d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWWRFLD_REAL_2D

  integer function SWWRFLD_REAL_3D ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    real, intent(in) :: BUFFER(:,:,:)   ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_real_3d  = swwrfld(swathid, fieldname, starts, strides, &
      & edges, reshape(buffer, (/size(buffer)/) ) )
  end function SWWRFLD_REAL_3D

end module SWAPI_REAL

! $Log$
