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

!$ Log: $
