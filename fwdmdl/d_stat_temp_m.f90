! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_STAT_TEMP_M

  implicit NONE
  private
  public :: D_STAT_TEMP, STAT_TEMP

  interface STAT_TEMP; module procedure D_STAT_TEMP; end interface

!{ This module calculates the energy emitted at the input frequency(MHz)
! from a black body source with temperature TEMP kelvins. The Returned
! energy is the statistical temperature equivalent in kelvins:

!{ \begin{equation*}
!  B = \frac{h \nu}{k \left( \exp(\frac{h \nu}{k T}) -1 \right )}
!  \end{equation*}

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! pure function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
  elemental function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
    use MLSCommon, ONLY: R8
    use PHYSICS, only: H_OVER_K
    integer, parameter :: RK = r8
    include 'stat_temp.f9h'
  end function D_STAT_TEMP
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module D_STAT_TEMP_M

! $Log$
! Revision 2.1  2002/09/27 00:10:57  vsnyder
! Move USEs from module scope to procedure scope, cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.5.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
