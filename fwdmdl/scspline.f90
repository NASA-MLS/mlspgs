! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SCSPLINE_DER_M

  implicit NONE
  private
  public :: SCSPLINE_DER, CSPLINE_DER

  interface CSPLINE_DER; module procedure SCSPLINE_DER; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_useS_here 
!---------------------------------------------------------------------------

contains
  subroutine SCSPLINE_DER ( XIN, XOUT, YIN, YOUT, DYOUT, NIN, NOUT, YMIN, YMAX)
    use S_HUNT_M, only: HUNT
    use S_PCSPL_M, only: PCSPL
    integer, parameter :: RK = kind(0.0e0)
    include 'cspline_der.f9h'
  end subroutine SCSPLINE_DER
  logical function not_useS_here()
    not_useS_here = (id(1:1) == ModuleName(1:1))
  end function not_useS_here

end module SCSPLINE_DER_M

! $Log$
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
