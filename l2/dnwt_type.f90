! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DNWT_TYPE

! This is a separate module only because the parameter RK is needed in
! the interface body for FANDJ

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, public, parameter :: RK = kind(0.0d0)

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

end module DNWT_TYPE

! $Log$
! Revision 2.7  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2002/10/07 23:43:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.3  2001/02/06 23:25:29  vsnyder
! Initial commit
!
! Revision 2.1  2001/02/06 23:23:53  vsnyder
! Initial commit
