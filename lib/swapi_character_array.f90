module SWAPI_CHARACTER_ARRAY

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function SWRDFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character, intent(out) :: BUFFER(*) ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_character_array = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_CHARACTER_ARRAY

  integer function SWWRFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character, intent(in) :: BUFFER(*)  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_character_array = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_CHARACTER_ARRAY

end module SWAPI_CHARACTER_ARRAY

! $Log$
