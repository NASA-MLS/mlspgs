module D_QSORT_M
  implicit NONE
  private
  public D_QSORT, QSORT
  interface QSORT; module procedure D_QSORT; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  subroutine d_qsort (N,RA)
    include 'qsort.f9h'
  end subroutine d_qsort
end module D_QSORT_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
!
