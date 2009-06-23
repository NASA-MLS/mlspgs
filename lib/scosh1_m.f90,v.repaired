! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SCOSH1_M

!--S replaces "?": ?COSH1, ?COSH1_M

  implicit NONE
  private
  public :: COSH1, SCOSH1

  interface COSH1; module procedure SCOSH1; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  elemental function SCOSH1 ( X ) result ( COSH1 )
    integer, parameter :: RK = kind(0.0e0)
    include 'cosh1.f9h'
  end function SCOSH1

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module SCOSH1_M

! $Log$
! Revision 2.4  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/11 20:36:26  vsnyder
! Make it elemental
!
! Revision 2.1  2002/10/11 19:44:00  vsnyder
! Initial commit
!
