module GEOC_GEOD_CONV_M
  use MLSCommon, only: R8
  use Geometry, only: GeodToGeocLat
  use UNITS, only: DEG2RAD, RAD2DEG
  use L2PC_PFA_STRUCTURES, only: EARTH_MAJOR, EARTH_MINOR
  use ELLIPSE_M, only: ELLIPSE
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

  ! ------------------------------------- --------  GEOC_GEOD_CONV  -----
  subroutine GEOC_GEOD_CONV (elvar,BETA_INC,PHI_TAN,GEOD_LAT,GEOC_LAT,Rp)

    real(r8), intent(in) :: BETA_INC    ! In Degrees
    real(r8), intent(in) :: PHI_TAN     ! In Radians
    real(r8), intent(in) :: GEOD_LAT    ! In Radians
    
    type(ELLIPSE), intent(in out) :: elvar

    real(r8), intent(out) :: GEOC_LAT   ! In Radians

    real(r8), intent(out) :: RP         ! In Kilometers

    real(r8) :: CW, INCL, SW, r
!   real(r8), parameter :: B2 = real(earth_minor,r8) ** 2 ! Doesn't work
    real(r8), parameter :: B2 = (earth_minor + 0.0_r8) ** 2

! Fill in some variables in the Ellipse module

    elvar%a2 = real(earth_major,r8) ** 2

    incl = (beta_inc - 90.0_r8) * deg2rad
!   r = tan(incl) ** 2

! c is the Minor axis for 2D ellipse, c2 = c * c
! c2oa2 is c**2 / a**2

    elvar%c2oa2 = b2/(elvar%a2 * cos(incl)**2 + b2 * sin(incl)**2 )
!   c2oa2 = (1.0_r8+r)*b2/(a2+b2*r)

    elvar%c2 = elvar%a2 * elvar%c2oa2

    elvar%phi_tan = Phi_tan
    elvar%spt = sin(Phi_tan)
    elvar%cpt = cos(Phi_tan)
    sw = elvar%spt * elvar%spt
    cw = elvar%cpt * elvar%cpt

    elvar%nphi_tan = elvar%a2 / sqrt(elvar%c2*sw + elvar%a2*cw)

!  Compute Radius of Curvature Circle (RoC) and its center coordinates:

    elvar%RoC = elvar%nphi_tan * sqrt(sw + elvar%c2oa2**2 * cw)
    elvar%XoC = (elvar%nphi_tan - elvar%RoC) * elvar%cpt
    elvar%YoC = (elvar%c2oa2 * elvar%nphi_tan - elvar%RoC) * elvar%spt

! Get the Geocentric Lat. (geoc_lat) from the Geodetic Lat.

    geoc_lat = GeodToGeocLat ( geod_lat * rad2deg )

!  Compute Earth Radius (Elliptical)

    r = elvar%a2*cw + b2*sw
    Rp = sqrt( ((elvar%a2*elvar%a2)*cw + (b2*b2)*sw) / r)

    return

  end subroutine GEOC_GEOD_CONV

end module GEOC_GEOD_CONV_M

! $Log$
! Revision 1.6  2001/03/28 19:55:26  vsnyder
! Revised some computations to make them more stable.
! Corrected units on call to GeodToGeocLat.
!
! Revision 1.4  2001/03/27 19:48:00  vsnyder
! Get Deg2Rad from Units, add some comments, make everything but
! geoc_geod_conv private, make constants _r8 instead of d0, cosmetic changes
!
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
