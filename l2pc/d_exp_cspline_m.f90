module D_EXP_CSPLINE_M
  use D_CSPLINE_M, only: CSPLINE
  implicit NONE
  private
  public D_EXP_CSPLINE, EXP_CSPLINE
  interface EXP_CSPLINE; module procedure D_EXP_CSPLINE; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
   "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
   "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine d_exp_cspline (XIN,XOUT,YIN,YOUT,NIN,NOUT)
    include 'exp_cspline.f9h'
  end subroutine d_exp_cspline
end module D_EXP_CSPLINE_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
