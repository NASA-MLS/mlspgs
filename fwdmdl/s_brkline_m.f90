module S_BRKLINE_M
  implicit NONE
  private
  public :: BRKLINE   ! generic
  public :: S_BRKLINE ! specific

  interface BRKLINE; module procedure S_BRKLINE; end interface

  integer, private, parameter :: RK = kind(0.0e0)

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine S_BRKLINE ( LINE, V, NV, IO )
    include 'brkline.f9h'
  end subroutine S_BRKLINE
end module S_BRKLINE_M

! $Log$
