module GEOC_GEOD_CONV_M
  use MLSCommon, only: R8
  use Geometry  !  , only: GeodToGeocLat
  use STRINGS, only: STRLWR
  use L2PC_FILE_PARAMETERS, only: DEG2RAD
  use L2PC_PFA_STRUCTURES, only: EARTH_MAJOR, EARTH_MINOR
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, &
                     NPHI_TAN, ROC, XOC, YOC
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!---------------------------------------------------------------------

SUBROUTINE geoc_geod_conv(beta_inc,phi_tan,geod_lat,geoc_lat,Rp)

Real(r8), INTENT(IN) :: beta_inc                ! In Degrees
Real(r8), INTENT(IN) :: phi_tan                 ! In Radians
Real(r8), INTENT(IN) :: geod_lat                ! In Radians

Real(r8), INTENT(OUT) :: geoc_lat               ! In Radians

Real(r8), INTENT(OUT) :: rp

Real(r8) :: q,r,incl,cw,sw,b2

  q = Real(earth_major,r8)
  a2 = q * q

  r = Real(earth_minor,r8)
  b2 = r * r

  incl = (beta_inc - 90.0D0) * deg2rad
  q = TAN(incl)
  r = q * q

! c is the Minor axis for 2D ellipse, c2 = c * c

  c2 = (1.0D0+r)*a2*b2/(a2+b2*r)

  c2oa2 = c2 / a2
!
! Get the Geocentric Lat. (geoc_lat) from the Geodetic Lat.

  geoc_lat = GeodToGeocLat(geod_lat)
!
  spt = SIN(Phi_tan)
  cpt = COS(Phi_tan)
  sw = spt * spt
  cw = 1.0d0 - sw

  nphi_tan = a2 / SQRT(c2-(c2-a2)*cw)

!  Compute Radius of Curvature Circle (RoC) and its center coordinates:

  r = c2 * cpt / a2
  q = spt * spt + r * r
  RoC = nphi_tan * SQRT(q)
  XoC = (nphi_tan - RoC) * cpt
  YoC = (c2oa2 * nphi_tan - RoC) * spt

!  Compure Earth Radius (Elliptical)

  q = ((a2*a2)*cw+(b2*b2)*sw)/(a2*cw+b2*sw)
  Rp = SQRT(q)

  Return

 END SUBROUTINE geoc_geod_conv

End module GEOC_GEOD_CONV_M
! $Log$
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
