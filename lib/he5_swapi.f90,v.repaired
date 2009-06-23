! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HE5_SWAPI

  use HE5_SWAPI_CHARACTER_SCALAR, only: HE5_SWRDFLD_CHARACTER_SCALAR, &
    &               HE5_SWWRFLD_CHARACTER_SCALAR, &
    & HE5_SWWRATTR_CHARACTER_SCALAR, HE5_SWWRLATTR_CHARACTER_SCALAR, &
    & HE5_SWRDATTR_CHARACTER_SCALAR, HE5_SWRDLATTR_CHARACTER_SCALAR
  use HE5_SWAPI_CHARACTER_ARRAY, only: HE5_SWRDFLD_CHARACTER_ARRAY, &
                                   HE5_SWWRFLD_CHARACTER_ARRAY, &
    & HE5_SWWRATTR_CHARACTER_ARRAY, HE5_SWWRLATTR_CHARACTER_ARRAY, &
    & HE5_SWRDATTR_CHARACTER_ARRAY, HE5_SWRDLATTR_CHARACTER_ARRAY
  use HE5_SWAPI_DOUBLE, only: HE5_SWRDFLD_DOUBLE, HE5_SWRDFLD_DOUBLE_2D, &
    & HE5_SWRDFLD_DOUBLE_3D, HE5_SWWRFLD_DOUBLE, HE5_SWWRFLD_DOUBLE_2D, &
    & HE5_SWWRFLD_DOUBLE_3D, &
    & HE5_SWWRATTR_DOUBLE, HE5_SWWRLATTR_DOUBLE, HE5_SWSETFILL_DOUBLE, &
    & HE5_SWRDATTR_DOUBLE, HE5_SWRDLATTR_DOUBLE
  use HE5_SWAPI_INTEGER, only: HE5_SWRDFLD_INTEGER, HE5_SWRDFLD_INTEGER_2D, &
    & HE5_SWRDFLD_INTEGER_3D, HE5_SWWRFLD_INTEGER, HE5_SWWRFLD_INTEGER_2D, &
    & HE5_SWWRFLD_INTEGER_3D, &
    & HE5_SWWRATTR_INTEGER, HE5_SWWRLATTR_INTEGER, HE5_SWSETFILL_INTEGER, &
    & HE5_SWRDATTR_INTEGER, HE5_SWRDLATTR_INTEGER
  use HE5_SWAPI_REAL, only: HE5_SWRDFLD_REAL, HE5_SWRDFLD_REAL_2D, &
    & HE5_SWRDFLD_REAL_3D, HE5_SWWRFLD_REAL, HE5_SWWRFLD_REAL_2D, &
    & HE5_SWWRFLD_REAL_3D, &
    & HE5_SWWRATTR_REAL, HE5_SWWRLATTR_REAL, HE5_SWSETFILL_REAL, &
    & HE5_SWRDATTR_REAL, HE5_SWRDLATTR_REAL

  private
  public :: HE5_SWRDFLD, HE5_SWWRFLD, HE5_EHWRGLATT, HE5_SWWRATTR, &
    & HE5_SWWRLATTR, HE5_SWSETFILL, HE5_SWRDATTR, HE5_SWRDLATTR

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

  interface HE5_SWWRATTR
    module procedure HE5_SWWRATTR_CHARACTER_SCALAR, HE5_SWWRATTR_DOUBLE, &
    HE5_SWWRATTR_INTEGER, HE5_SWWRATTR_REAL, HE5_SWWRATTR_CHARACTER_ARRAY
  end interface

  interface HE5_SWWRLATTR
    module procedure HE5_SWWRLATTR_CHARACTER_SCALAR, HE5_SWWRLATTR_DOUBLE, &
    HE5_SWWRLATTR_INTEGER, HE5_SWWRLATTR_REAL, HE5_SWWRLATTR_CHARACTER_ARRAY
  end interface

  interface HE5_SWRDATTR
    module procedure HE5_SWRDATTR_CHARACTER_SCALAR, HE5_SWRDATTR_DOUBLE, &
    HE5_SWRDATTR_INTEGER, HE5_SWRDATTR_REAL, HE5_SWRDATTR_CHARACTER_ARRAY
  end interface

  interface HE5_SWRDLATTR
    module procedure HE5_SWRDLATTR_CHARACTER_SCALAR, HE5_SWRDLATTR_DOUBLE, &
    HE5_SWRDLATTR_INTEGER, HE5_SWRDLATTR_REAL, HE5_SWRDLATTR_CHARACTER_ARRAY
  end interface

  interface HE5_SWSETFILL
    module procedure HE5_SWSETFILL_DOUBLE, &
    HE5_SWSETFILL_INTEGER, HE5_SWSETFILL_REAL
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains 
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HE5_SWAPI

! $Log$
! Revision 2.8  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.6  2004/02/13 00:16:02  pwagner
! New stuff for reading swath attributes
!
! Revision 2.5  2003/04/11 23:32:23  pwagner
! Moved he5_swsetfill he5_ehwrglatt interfaces
!
! Revision 2.4  2003/02/10 22:05:52  pwagner
! New attributes writing for char arrays
!
! Revision 2.3  2003/02/08 00:32:54  pwagner
! New attribute-related interfaces
!
