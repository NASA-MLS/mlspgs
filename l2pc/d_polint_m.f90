module D_POLINT_M
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public D_POLINT, POLINT
  interface POLINT; module procedure D_POLINT; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine D_POLINT (XA,YA,N,X,Y)
    include 'polint.f9h'
  end subroutine D_POLINT
end module D_POLINT_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90

