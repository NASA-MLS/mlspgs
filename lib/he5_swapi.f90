! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HE5_SWAPI

  use HE5_SWAPI_CHARACTER_SCALAR, only: HE5_SWRDFLD_CHARACTER_SCALAR, &
    &               HE5_SWWRFLD_CHARACTER_SCALAR, &
    & HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_SWWRATTR_CHARACTER_SCALAR, &
    & HE5_SWWRLATTR_CHARACTER_SCALAR
  use HE5_SWAPI_CHARACTER_ARRAY, only: HE5_SWRDFLD_CHARACTER_ARRAY, &
                                   HE5_SWWRFLD_CHARACTER_ARRAY, &
    & HE5_EHWRGLATT_CHARACTER_ARRAY, HE5_SWWRATTR_CHARACTER_ARRAY, &
    & HE5_SWWRLATTR_CHARACTER_ARRAY
  use HE5_SWAPI_DOUBLE, only: HE5_SWRDFLD_DOUBLE, HE5_SWRDFLD_DOUBLE_2D, &
       HE5_SWRDFLD_DOUBLE_3D, HE5_SWWRFLD_DOUBLE, HE5_SWWRFLD_DOUBLE_2D, &
       HE5_SWWRFLD_DOUBLE_3D, &
       HE5_EHWRGLATT_DOUBLE, HE5_SWWRATTR_DOUBLE, HE5_SWWRLATTR_DOUBLE
  use HE5_SWAPI_INTEGER, only: HE5_SWRDFLD_INTEGER, HE5_SWRDFLD_INTEGER_2D, &
       HE5_SWRDFLD_INTEGER_3D, HE5_SWWRFLD_INTEGER, HE5_SWWRFLD_INTEGER_2D, &
       HE5_SWWRFLD_INTEGER_3D, &
       HE5_EHWRGLATT_INTEGER, HE5_SWWRATTR_INTEGER, HE5_SWWRLATTR_INTEGER
  use HE5_SWAPI_REAL, only: HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, &
       HE5_SWRDFLD_REAL_3D, HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, &
       HE5_SWWRFLD_REAL_3D, &
       HE5_EHWRGLATT_REAL, HE5_SWWRATTR_REAL, HE5_SWWRLATTR_REAL

  private
  public :: HE5_SWRDFLD, HE5_SWWRFLD, HE5_EHWRGLATT, HE5_SWWRATTR, HE5_SWWRLATTR

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

  interface HE5_EHWRGLATT   ! From its name, might better be in he5_ehapi.f90
    module procedure HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHWRGLATT_DOUBLE, &
    HE5_EHWRGLATT_INTEGER, HE5_EHWRGLATT_REAL, HE5_EHWRGLATT_CHARACTER_ARRAY
  end interface

  interface HE5_SWWRATTR
    module procedure HE5_SWWRATTR_CHARACTER_SCALAR, HE5_SWWRATTR_DOUBLE, &
    HE5_SWWRATTR_INTEGER, HE5_SWWRATTR_REAL, HE5_SWWRATTR_CHARACTER_ARRAY
  end interface

  interface HE5_SWWRLATTR
    module procedure HE5_SWWRLATTR_CHARACTER_SCALAR, HE5_SWWRLATTR_DOUBLE, &
    HE5_SWWRLATTR_INTEGER, HE5_SWWRLATTR_REAL, HE5_SWWRLATTR_CHARACTER_ARRAY
  end interface

  !------------------------------- RCS Ident Info ----------------------------
  character(len=256), private :: Id = &
    & "$Id$"
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HE5_SWAPI

! $Log$
! Revision 2.4  2003/02/10 22:05:52  pwagner
! New attributes writing for char arrays
!
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
