program test_geom

use MLSCommon
use heights_module
real(kind=r8)::geom,gph,geom2
!real(kind=r8)::latc,latd
integer::i

do i=0,100
  geom=i*1000 ! all heights are in metres
  gph=geom_to_gph(geom,89.9_r8)
  geom2=gph_to_geom(gph,89.9_r8)

  print "(5g16.8)",geom,gph,geom2,geom2-geom,geom-gph
end do

!do i=-90,90
!  latd=i
!  latc=lat_geod_to_geoc(latd)
!  print*,latd,latc,latd-latc
!end do


end program test_geom
