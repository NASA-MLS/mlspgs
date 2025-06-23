! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

MODULE oa_file_contents

  use MLSKinds, only: R8

  IMPLICIT NONE

  TYPE emls_oa_t
     sequence

     ! for satellite:

     real :: sc_orbincl(32)
     real :: sc_geoclat(32)         ! Geocentric latitude, degrees
     real :: sc_geodlat(32)         ! Geodetic latitude, degrees
     real :: sc_lon(32)             ! Longitude, degrees
     real(r8) :: sc_geocalt(32)     ! Geocentric altitude, meters
     real(r8) :: sc_geodalt(32)     ! Geodetic height, meters
     real(r8) :: sc_MIF_TAI(32)
     real(r8) :: sc_ECI(3,32)       ! ECI position, meters
     real(r8) :: sc_ECR(3,32)       ! ECR position, meters
     real(r8) :: sc_VelECI(3,32)    ! ECI velocity, meters/second
     real(r8) :: sc_VelECR(3,32)    ! ECR velocity, meters/second
     real(r8) :: sc_ypr(3,32)
     real(r8) :: sc_ypr_rate(3,32)

     ! for instrument:

     real(r8) :: geoc_alt(32)       ! Geocentric altitude, meters
     real(r8) :: geod_alt(32)       ! Geodetic height, meters
     real :: azimAngle(32)
     real :: geoc_lat(32)           ! Geocentric latitude, degrees
     real :: geod_lat(32)           ! Geodetic latitude, degrees
     real :: geod_ang(32)           ! Orbit angle, degrees
     real :: lon(32)                ! longitude, degrees
     real :: scanAngle(32)
     real :: solartime(32)
     real :: solarzenith(32)
     integer :: bo_stat(32)
     real(r8) :: MAFStartTimeTAI
     real(r8) :: ECI(3,32)          ! ECI position, meters
     real(r8) :: ECR(3,32)          ! ECR position, meters
     real(r8) :: ECRtoFOV(9,32)
     real(r8) :: LosAngle(32)
     real(r8) :: LosVel(32)         ! Line-of-sight velocity, meters/second
     real(r8) :: OrbY(32)

  END TYPE

END MODULE oa_file_contents

! $Log$
! Revision 1.2  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
