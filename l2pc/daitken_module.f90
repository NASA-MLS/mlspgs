module DAITKEN_MODULE
!--D replaces "?": ?AITKEN_MODULE ?AITKEN
  implicit NONE
  private
  public :: AITKEN, DAITKEN, RK
  integer, parameter :: RK=kind(0.0d0)
  interface AITKEN; module procedure DAitken; end interface
  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------
contains
! ================================================     DAitken     =====
! pure real(rk) Function DAitken ( XN, XNP1, XNP2 ) result (AITKEN)
  real(rk) Function DAitken ( XN, XNP1, XNP2 ) result (AITKEN)
    include 'aitken.f9h'
  End Function DAitken
end module DAITKEN_MODULE
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
