module SWAPI_CHARACTER_SCALAR

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  integer function SWRDFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER  ! Buffer for read

    integer, external :: SWRDFLD

    swrdfld_character_scalar = swrdfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWRDFLD_CHARACTER_SCALAR

  integer function SWWRFLD_CHARACTER_SCALAR ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    integer, external :: SWWRFLD

    swwrfld_character_scalar = swwrfld(swathid, fieldname, starts, strides, &
      & edges, buffer )
  end function SWWRFLD_CHARACTER_SCALAR

end module SWAPI_CHARACTER_SCALAR

! $Log$
! Revision 1.1  2000/09/02 01:56:36  vsnyder
! Initial code
!
