module GEOC_GEOD_CONV_M
  use MLSCommon, only: R8
  use Geometry, only: GeodToGeocLat
  use STRINGS, only: STRLWR
  use UNITS, only: DEG2RAD
  use L2PC_PFA_STRUCTURES, only: EARTH_MAJOR, EARTH_MINOR
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, &
                     NPHI_TAN, ROC, XOC, YOC

  implicit NONE
  private
  public :: GEOC_GEOD_CONV

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

  ! ---------------------------------------------  GEOC_GEOD_CONV  -----
  subroutine GEOC_GEOD_CONV ( BETA_INC, PHI_TAN, GEOD_LAT, GEOC_LAT, Rp)

    real(r8), intent(in) :: BETA_INC    ! In Degrees
    real(r8), intent(in) :: PHI_TAN     ! In Radians
    real(r8), intent(in) :: GEOD_LAT    ! In Radians

    real(r8), intent(out) :: GEOC_LAT   ! In Radians

    real(r8), intent(out) :: RP         ! In Kilometers

    real(r8) :: q, r, incl, cw, sw, b2

    q = real(earth_major,r8)
    a2 = q * q

    r = real(earth_minor,r8)
    b2 = r * r

    incl = (beta_inc - 90.0_r8) * deg2rad
    q = tan(incl)
    r = q * q

! c is the Minor axis for 2D ellipse, c2 = c * c

    c2 = (1.0_r8+r)*a2*b2/(a2+b2*r)

    c2oa2 = c2 / a2
!
! Get the Geocentric Lat. (geoc_lat) from the Geodetic Lat.

    geoc_lat = GeodToGeocLat(geod_lat)
!
    spt = sin(Phi_tan)
    cpt = cos(Phi_tan)
    sw = spt * spt
    cw = 1.0_r8 - sw

    nphi_tan = a2 / sqrt(c2-(c2-a2)*cw)

!  Compute Radius of Curvature Circle (RoC) and its center coordinates:

    r = c2 * cpt / a2
    q = spt * spt + r * r
    RoC = nphi_tan * sqrt(q)
    XoC = (nphi_tan - RoC) * cpt
    YoC = (c2oa2 * nphi_tan - RoC) * spt

!  Compute Earth Radius (Elliptical)

    q = ((a2*a2)*cw+(b2*b2)*sw)/(a2*cw+b2*sw)
    Rp = sqrt(q)

    return

  end subroutine GEOC_GEOD_CONV

end module GEOC_GEOD_CONV_M

! $Log$
! Revision 1.3  2001/03/20 23:22:40  zvi
! Change to new geoc_geod routine..
!
! Revision 1.2  2001/03/20 11:03:15  zvi
! Fixing code for "real" data run, increase dim. etc.
!
! Revision 1.1  2001/01/31 22:40:12  zvi
! Add new version
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
