module D_QSORT2_M
  implicit NONE
  private
  public D_QSORT2, QSORT2
  interface QSORT2; module procedure D_QSORT2; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine d_qsort2 (N,RA,RB)
    include 'qsort2.f9h'
  end subroutine d_qsort2
end module D_QSORT2_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
!
