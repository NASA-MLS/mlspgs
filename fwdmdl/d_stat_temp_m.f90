module D_STAT_TEMP_M
  use MLSCommon, only: R8
  use PHYSICS, only: H_OVER_K
  implicit NONE
  private
  public :: D_STAT_TEMP, STAT_TEMP

  interface STAT_TEMP; module procedure D_STAT_TEMP; end interface

!{ This module calculates the energy emitted at the input frequency(MHz)
! from a black body source with temperature TEMP kelvins. The Returned
! energy is the statistical temperature equivalent in kelvins:
!
! \begin{equation*}
!  B = \frac{h \nu}{k \left( \exp(\frac{h \nu}{k T}) -1 \right )}
!  \end{equation*}
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = r8
contains
! pure function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
  function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
    include 'stat_temp.f9h'
  End function D_STAT_TEMP
end module D_STAT_TEMP_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
