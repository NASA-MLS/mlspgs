module HE5_SWAPI_CHARACTER_SCALAR

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

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

end module HE5_SWAPI_CHARACTER_SCALAR

