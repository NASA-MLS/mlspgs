module SAITKEN_MODULE
!--S replaces "?": ?AITKEN_MODULE ?AITKEN
  implicit NONE
  private
  public :: AITKEN, SAITKEN, RK
  integer, parameter :: RK=kind(0.0e0)

  interface AITKEN; module procedure SAitken; end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

contains
! ================================================     DAitken     =====
  pure real(rk) Function SAitken ( XN, XNP1, XNP2 ) result (AITKEN)
    include 'aitken.f9h'
  End Function SAitken
end module SAITKEN_MODULE

! $Log$
