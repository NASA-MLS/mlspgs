! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DREXP_M

! ----------------------------------------------------------------------
!            EVALUATE THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------

  implicit NONE
  private
  public :: REXP, DREXP

  interface REXP; module procedure DREXP; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

  function DREXP ( X ) result ( REXP )
    integer, parameter :: RK = kind(0.0d0)
    include 'rexp.f9h'
  end function DREXP

end module DREXP_M
