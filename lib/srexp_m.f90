! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SREXP_M

! ----------------------------------------------------------------------
!            EVALUATE THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------

  implicit NONE
  private
  public :: REXP, SREXP

  interface REXP; module procedure SREXP; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

  function SREXP ( X ) result ( REXP )
    integer, parameter :: RK = kind(0.0e0)
    include 'rexp.f9h'
  end function SREXP

end module SREXP_M

! $Log$
! Revision 2.2  2002/10/10 20:36:50  vsnyder
! Darn, I forgot the CVS log comment
!
