module ELLIPSE
  use MLSCommon, only: R8
  implicit NONE
  public
  real(r8) :: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, HT2, RR, &
 &            PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC
  logical :: EARTHX
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
end module ELLIPSE
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
