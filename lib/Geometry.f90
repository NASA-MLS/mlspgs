! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Geometry

  ! This module contains some geometry routines common to the forward model and
  ! the scan model.

  use MLSCommon, only: R8, RP

  implicit none
  private
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------


  public :: EarthRadA, EarthRadB, GM, J2, J4, W, PI, LN10, GeodToGeocLat

  ! Now parameters etc

  real(r8), parameter :: EarthRadA=6378137.0_r8    ! Major axis radius in m 
  real(r8), parameter :: EarthRadB=6356752.3141_r8 ! Minor axis radius in m 

  ! These are the 1980 reference geoid values.

  real(rp), parameter :: GM = 3.98600436e14_rp ! m^3/sec^2
  real(rp), parameter :: J2 = 0.0010826256_rp
  real(rp), parameter :: J4 = -.0000023709122_rp

  ! earth rotational velocity

  real(rp), parameter :: W = 7.292115e-05_rp ! rad/sec

  ! Just in case these constants aren't defined everywhere else.
  
  real(r8), parameter :: PI=3.14159265358979323846264338328_r8   
  real(r8), parameter :: LN10=2.30258509299404568401799145468_r8 

contains ! ------------------------------- Subroutines and functions ----

  ! This function converts a geodetic latitude (IN DEGREES!) into a geocentric
  ! one (IN RADIANS!)
  
  real(r8) elemental function GeodToGeocLat(geodLat)

    ! Arguments
    real (r8), intent(IN) :: geodLat

    ! Executable code, use special method for high latitudes.
    if (abs(geodLat).gt.89.0) then
       if (geodLat.gt.0.0) then
          geodtoGeocLat=PI/2.0+((earthRadA/earthRadB)**2)*&
               (geodLat-90.0)*PI/180.0
       else
          geodToGeocLat=-PI/2.0+((earthRadA/earthRadB)**2)*&
               (geodLat+90.0)*PI/180.0
       endif
    else
       geodToGeocLat=atan(((earthRadB/earthRadA)**2)*tan(geodLat*PI/180.0))
    endif
  end function GeodToGeocLat
  
end module Geometry

! $Log$
! Revision 2.5  2002/09/26 16:27:14  livesey
! Changes from Van, new constants etc.
!
! Revision 2.4  2001/03/28 19:50:19  vsnyder
! Change constants from d0 to _r8
!
! Revision 2.3  2001/03/27 18:09:11  vsnyder
! Revised CVS stuff to use CHARACTER(len=*) parameter
!
! Revision 2.2  2001/03/26 13:23:34  livesey
! Added CVS stuff
!
