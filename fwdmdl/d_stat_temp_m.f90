! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! pure function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
  elemental function D_STAT_TEMP ( TEMP, FREQ ) result ( STAT_TEMP )
    use MLSKinds, ONLY: R8
    use PHYSICS, only: H_OVER_K
    integer, parameter :: RK = r8
    include 'stat_temp.f9h'
  end function D_STAT_TEMP

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module D_STAT_TEMP_M

! $Log$
! Revision 2.4  2008/07/17 21:03:58  vsnyder
! Get real kind from MLSKinds instead of MLSCommon
!
! Revision 2.3  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
