module HE5_SWAPI_CHARACTER_ARRAY

  public

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  integer function HE5_SWRDFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(out) :: BUFFER(:) ! Buffer for read

    integer, external :: HE5_SWRDFLD

    HE5_SWrdfld_character_array = HE5_SWrdfld(swathid, fieldname, starts, &
         strides, edges, buffer )
  end function HE5_SWRDFLD_CHARACTER_ARRAY

  integer function HE5_SWWRFLD_CHARACTER_ARRAY ( SWATHID, FIELDNAME, &
    & STARTS, STRIDES, EDGES, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, intent(in) :: STARTS(*)    ! Start array
    integer, intent(in) :: STRIDES(*)   ! Stride array
    integer, intent(in) :: EDGES(*)     ! Edge array
    character(len=*), intent(in) :: BUFFER(:)  ! Buffer for write

    integer, external :: HE5_SWWRFLD

    HE5_SWwrfld_character_array=HE5_SWwrfld(swathid, fieldname, starts, &
         strides, edges, buffer )
  end function HE5_SWWRFLD_CHARACTER_ARRAY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI_CHARACTER_ARRAY

