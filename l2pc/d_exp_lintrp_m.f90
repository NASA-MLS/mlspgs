module D_EXP_LINTRP_M
  implicit NONE
  private
  public D_EXP_LINTRP, EXP_LINTRP
  interface EXP_LINTRP; module procedure D_EXP_LINTRP; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
   "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
   "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine d_exp_lintrp (XIN,XOUT,YIN,YOUT,NIN,NOUT)
    include 'exp_lintrp.f9h'
  end subroutine d_exp_lintrp
end module D_EXP_LINTRP_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
!
