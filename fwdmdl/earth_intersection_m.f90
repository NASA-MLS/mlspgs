module EARTH_INTERSECTION_M
  use MLSCommon, only: R8
  use ELLIPSE_M, only: ELLIPSE
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
  Subroutine Earth_Intersection (elvar,Rs)

    Type(ELLIPSE), intent(IN OUT) :: elvar
    real(r8), intent(out) :: RS
!
    real(r8) :: r,V
!
    elvar%EarthX = .true.            ! Set Earth Intersecting Ray flag
    v = (elvar%ht + elvar%RoC)/elvar%RoC
    elvar%Phi_s = elvar%Phi_tan - Acos(v)
    elvar%sps = Sin(elvar%Phi_s)
    elvar%cps = Cos(elvar%Phi_s)
    r = elvar%a2*elvar%cps*elvar%cps + elvar%c2*elvar%sps*elvar%sps
    elvar%NPhi_s = elvar%a2 / Sqrt(r)
    v = 2.0_r8 * elvar%Phi_s - elvar%Phi_tan
    elvar%spts = Sin(v)
    elvar%cpts = Cos(v)
!
    v = elvar%c2oa2 * elvar%sps
    Rs = elvar%NPhi_s * Sqrt(elvar%cps*elvar%cps + v*v)
!
    Return
  End Subroutine Earth_Intersection
end module EARTH_INTERSECTION_M
! $Log$
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
