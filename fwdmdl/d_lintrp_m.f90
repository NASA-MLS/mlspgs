module D_LINTRP_M
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public D_LINTRP, LINTRP
  interface LINTRP; module procedure D_LINTRP; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = kind(0.0d0)

contains

  subroutine D_LINTRP ( XIN, XOUT, YIN, YOUT, NIN, NOUT )
    include 'lintrp.f9h'
  end subroutine D_LINTRP
end module D_LINTRP_M

! $Log$
