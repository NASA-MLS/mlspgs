! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DNWT_TYPE

! This is a separate module only because the parameter RK is needed in
! the interface body for FANDJ

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256),private :: Id = &
       "$Id$"
  CHARACTER (LEN=*), private, PARAMETER :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, public, parameter :: RK = kind(0.0d0)

contains 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DNWT_TYPE

! $Log$
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
