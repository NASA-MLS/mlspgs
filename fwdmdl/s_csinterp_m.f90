module S_CSINTERP_M
  use MLSCommon, only: I4
  implicit NONE
  private
  public :: S_CSINTERP, CSINTERP
  interface CSINTERP; module procedure S_CSINTERP; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter :: RK = kind(0.0e0)
contains

  function S_CSINTERP ( x_in, x_out, y_in, y_out,        &
 &                      len_in, len_out, max_len_in, max_len_out, &
 &                      a, b, c, d) result (RESULT)
    include 'csinterp.f9h'
  end function S_CSINTERP
end module S_CSINTERP_M

! $Log$
