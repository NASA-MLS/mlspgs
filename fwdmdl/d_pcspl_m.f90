module D_PCSPL_M
  implicit NONE
  private
  public :: D_PCSPL, PCSPL
  interface PCSPL; module procedure D_PCSPL; end interface

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  integer, parameter :: RK = kind(0.0d0)

contains

  subroutine D_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    include 'pcspl.f9h'
  end subroutine D_PCSPL
end module D_PCSPL_M

! $Log$
