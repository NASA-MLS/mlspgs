module HE5_SWAPI

  use HE5_SWAPI_CHARACTER_SCALAR, only: HE5_SWRDFLD_CHARACTER_SCALAR, &
    &                               HE5_SWWRFLD_CHARACTER_SCALAR
  use HE5_SWAPI_CHARACTER_ARRAY, only: HE5_SWRDFLD_CHARACTER_ARRAY, &
                                   HE5_SWWRFLD_CHARACTER_ARRAY
  use HE5_SWAPI_DOUBLE, only: HE5_SWRDFLD_DOUBLE, HE5_SWRDFLD_DOUBLE_2D, &
       HE5_SWRDFLD_DOUBLE_3D, HE5_SWWRFLD_DOUBLE, HE5_SWWRFLD_DOUBLE_2D, &
       HE5_SWWRFLD_DOUBLE_3D
  use HE5_SWAPI_INTEGER, only: HE5_SWRDFLD_INTEGER, HE5_SWRDFLD_INTEGER_2D, &
       HE5_SWRDFLD_INTEGER_3D, HE5_SWWRFLD_INTEGER, HE5_SWWRFLD_INTEGER_2D, &
       HE5_SWWRFLD_INTEGER_3D
  use HE5_SWAPI_REAL, only: HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, &
       HE5_SWRDFLD_REAL_3D, HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, &
       HE5_SWWRFLD_REAL_3D

  private
  public :: HE5_SWRDFLD, HE5_SWWRFLD

  interface HE5_SWRDFLD
    module procedure HE5_SWRDFLD_CHARACTER_SCALAR, HE5_SWRDFLD_CHARACTER_ARRAY
    module procedure HE5_SWRDFLD_DOUBLE, HE5_SWRDFLD_DOUBLE_2D, &
         HE5_SWRDFLD_DOUBLE_3D
    module procedure HE5_SWRDFLD_INTEGER, HE5_SWRDFLD_INTEGER_2D,&
         HE5_SWRDFLD_INTEGER_3D
    module procedure HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, HE5_SWRDFLD_REAL_3D
  end interface

  interface HE5_SWWRFLD
    module procedure HE5_SWWRFLD_CHARACTER_SCALAR, HE5_SWWRFLD_CHARACTER_ARRAY
    module procedure HE5_SWWRFLD_DOUBLE, HE5_SWWRFLD_DOUBLE_2D,&
         HE5_SWWRFLD_DOUBLE_3D
    module procedure HE5_SWWRFLD_INTEGER, HE5_SWWRFLD_INTEGER_2D, &
         HE5_SWWRFLD_INTEGER_3D
    module procedure HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, HE5_SWWRFLD_REAL_3D
  end interface

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

end module HE5_SWAPI

