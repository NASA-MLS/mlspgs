MODULE oa_file_contents

IMPLICIT NONE

TYPE emls_oa_t
   sequence
   
   ! for satellite:

   real :: sc_orbincl(32)
   real :: sc_geoclat(32)
   real :: sc_geodlat(32)
   real :: sc_lon(32)
   real*8 :: sc_geocalt(32)
   real*8 :: sc_geodalt(32)
   real*8 :: sc_MIF_TAI(32)
   real*8 :: sc_ECI(3,32)
   real*8 :: sc_ECR(3,32)
   real*8 :: sc_VelECI(3,32)
   real*8 :: sc_VelECR(3,32)
   real*8 :: sc_ypr(3,32)
   real*8 :: sc_ypr_rate(3,32)

   ! for instrument:

   real*8 :: geoc_alt(32)
   real*8 :: geod_alt(32)
   real :: azimAngle(32)
   real :: geoc_lat(32)
   real :: geod_lat(32)
   real :: geod_ang(32)
   real :: lon(32)
   real :: scanAngle(32)
   real :: solartime(32)
   real :: solarzenith(32)
   integer :: bo_stat(32)
   real*8 :: MAFStartTimeTAI
   real*8 :: ECI(3,32)
   real*8 :: ECR(3,32)
   real*8 :: ECRtoFOV(9,32)
   real*8 :: LosAngle(32)
   real*8 :: LosVel(32)
   real*8 :: OrbY(32)

END TYPE

TYPE (emls_oa_t) :: emls_oa

END MODULE oa_file_contents
