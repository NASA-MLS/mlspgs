module SWAPI

  use SWAPI_CHARACTER_SCALAR, only: SWRDFLD_CHARACTER_SCALAR, &
    &                               SWWRFLD_CHARACTER_SCALAR
  use SWAPI_CHARACTER_ARRAY, only: SWRDFLD_CHARACTER_ARRAY, &
                                   SWWRFLD_CHARACTER_ARRAY
  use SWAPI_DOUBLE, only: SWRDFLD_DOUBLE, SWRDFLD_DOUBLE_2D, SWRDFLD_DOUBLE_3D, &
    &                     SWWRFLD_DOUBLE, SWWRFLD_DOUBLE_2D, SWWRFLD_DOUBLE_3D
  use SWAPI_INTEGER, only: SWRDFLD_INTEGER, SWRDFLD_INTEGER_2D, SWRDFLD_INTEGER_3D, &
    &                      SWWRFLD_INTEGER, SWWRFLD_INTEGER_2D, SWWRFLD_INTEGER_3D
  use SWAPI_REAL, only: SWRDFLD_REAL, SWRDFLD_REAL_2D, SWRDFLD_REAL_3D, &
    &                   SWWRFLD_REAL, SWWRFLD_REAL_2D, SWWRFLD_REAL_3D

  private
  public :: SWRDFLD, SWWRFLD

  interface SWRDFLD
    module procedure SWRDFLD_CHARACTER_SCALAR, SWRDFLD_CHARACTER_ARRAY
    module procedure SWRDFLD_DOUBLE, SWRDFLD_DOUBLE_2D, SWRDFLD_DOUBLE_3D
    module procedure SWRDFLD_INTEGER, SWRDFLD_INTEGER_2D, SWRDFLD_INTEGER_3D
    module procedure SWRDFLD_REAL, SWRDFLD_REAL_2D, SWRDFLD_REAL_3D
  end interface

  interface SWWRFLD
    module procedure SWWRFLD_CHARACTER_SCALAR, SWWRFLD_CHARACTER_ARRAY
    module procedure SWWRFLD_DOUBLE, SWWRFLD_DOUBLE_2D, SWWRFLD_DOUBLE_3D
    module procedure SWWRFLD_INTEGER, SWWRFLD_INTEGER_2D, SWWRFLD_INTEGER_3D
    module procedure SWWRFLD_REAL, SWWRFLD_REAL_2D, SWWRFLD_REAL_3D
  end interface

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

end module SWAPI

! $Log$
