! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

end module SWAPI

! $Log$
! Revision 2.3  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/08 00:09:14  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/09/02 01:56:36  vsnyder
! Initial code
!
