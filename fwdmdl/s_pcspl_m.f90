module S_PCSPL_M
  implicit NONE
  private
  public :: S_PCSPL, PCSPL
  interface PCSPL; module procedure S_PCSPL; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, parameter :: RK = kind(0.0e0)

contains

  subroutine S_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    include 'pcspl.f9h'
  end subroutine S_PCSPL
end module S_PCSPL_M

! $Log$
