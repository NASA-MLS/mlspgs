module DNWT_TYPE

! This is a separate module only because the parameter RK is needed in
! the interface body for FANDJ

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, public, parameter :: RK = kind(0.0d0)

end module DNWT_TYPE

! $Log$
! Revision 2.3  2001/02/06 23:25:29  vsnyder
! Initial commit
!
! Revision 2.1  2001/02/06 23:23:53  vsnyder
! Initial commit
