! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SCOSH1_M

!--S replaces "?": ?COSH1, ?COSH1_M

  implicit NONE
  private
  public :: COSH1, SCOSH1

  interface COSH1; module procedure SCOSH1; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

  elemental function SCOSH1 ( X ) result ( COSH1 )
    integer, parameter :: RK = kind(0.0e0)
    include 'cosh1.f9h'
  end function SCOSH1
end module SCOSH1_M

! $Log$
! Revision 2.2  2002/10/11 20:36:26  vsnyder
! Make it elemental
!
! Revision 2.1  2002/10/11 19:44:00  vsnyder
! Initial commit
!
