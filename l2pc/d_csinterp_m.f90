module D_CSINTERP_M
  use MLSCommon, only: I4
  implicit NONE
  private
  public :: D_CSINTERP, CSINTERP
  interface CSINTERP; module procedure D_CSINTERP; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = kind(0.0d0)
contains
  function D_CSINTERP ( x_in, x_out, y_in, y_out,        &
 &                      len_in, len_out, max_len_in, max_len_out, &
 &                      a, b, c, d) result (RESULT)
    include 'csinterp.f9h'
  end function D_CSINTERP
end module D_CSINTERP_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
