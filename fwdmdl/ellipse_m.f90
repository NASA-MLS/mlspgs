module ELLIPSE_M
  use MLSCommon, only: R8
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
! This structure contains the "ELLIPSE parameters", containing various
! variables needed for computing various elliptical geometry entities.
  Type ELLIPSE
    Logical  :: EARTHX
    Real(r8) :: A2
    Real(r8) :: C2
    Real(r8) :: C2OA2
    Real(r8) :: CPT
    Real(r8) :: SPT
    Real(r8) :: CPS
    Real(r8) :: SPS
    Real(r8) :: CPTS
    Real(r8) :: SPTS
    Real(r8) :: HT
    Real(r8) :: HT2
    Real(r8) :: RR
    Real(r8) :: PHI_TAN
    Real(r8) :: NPHI_TAN
    Real(r8) :: PHI_S
    Real(r8) :: NPHI_S
    Real(r8) :: PS
    Real(r8) :: ROC
    Real(r8) :: XOC
    Real(r8) :: YOC
  End Type ELLIPSE
end module ELLIPSE_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
