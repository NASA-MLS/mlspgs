module EARTH_INTERSECTION_M
  use ELLIPSE, only: A2, CPS, CPTS, C2, C2OA2, EARTHX, HT, NPHI_S, &
                     PHI_S, PHI_TAN, ROC, SPS, SPTS
  use MLSCommon, only: R8
  implicit NONE
  private
  public :: Earth_Intersection

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
! This routine solves for the Phi_s and Rs for Earth intersection case.
!  ** Note: This routine is using The Equivalent Circle concept
!
  Subroutine Earth_Intersection ( Rs )
    real(r8), intent(out) :: RS
!
    real(r8) :: V

    EarthX = .true.            ! Set Earth Intersecting Ray flag
!
    Phi_s = Phi_tan - Acos( (ht+RoC)/RoC )
    sps = Sin(Phi_s)
    cps = Cos(Phi_s)
    NPhi_s = a2 / Sqrt(a2*cps*cps + c2*sps*sps)
    v = 2.0_r8 * Phi_s - Phi_tan
    spts = Sin(v)
    cpts = Cos(v)
!
    v = c2oa2 * sps
    Rs = NPhi_s * Sqrt(cps*cps + v*v)
!
    Return
  End Subroutine Earth_Intersection
end module EARTH_INTERSECTION_M

! $Log$
