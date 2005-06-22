! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DCOSH1_M

!--D replaces "?": ?COSH1, ?COSH1_M

  implicit NONE
  private
  public :: COSH1, DCOSH1

  interface COSH1; module procedure DCOSH1; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  elemental function DCOSH1 ( X ) result ( COSH1 )
    integer, parameter :: RK = kind(0.0d0)
    include 'cosh1.f9h'
  end function DCOSH1

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module DCOSH1_M

! $Log$
! Revision 2.3  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/11 20:36:26  vsnyder
! Make it elemental
!
! Revision 2.1  2002/10/11 19:44:00  vsnyder
! Initial commit
!
