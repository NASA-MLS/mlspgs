module D_BRKLINE_M
  implicit NONE
  private
  public :: BRKLINE   ! generic
  public :: D_BRKLINE ! specific
  interface BRKLINE; module procedure D_BRKLINE; end interface
  integer, private, parameter :: RK = kind(0.0d0)
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  subroutine D_BRKLINE ( LINE, V, NV, IO )
    include 'brkline.f9h'
  end subroutine D_BRKLINE
end module D_BRKLINE_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
