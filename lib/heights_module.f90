module heights_module

  ! Copyright (c) © 1999 the University of Edinburgh. All rights reserved 
  ! U.K. Government funding under NERC contract F14/6/33 acknowledged. 
  ! U.S. Government Sponsorship under NASA Contract NAS7-1407 also 
  ! acknowledged.

  ! This module provides conversions from geodetic (geometric height above
  ! earth surface)  to geopotential height (above geoid). 
  ! Formulae taken from  "Theoretical basis for EOS-MLS geopotential 
  ! and scan calculations" by N. J. Livesey

  !  use kinds_module
  use MLSCommon  !Only the kinds info is used
  !
  use MLSMessageModule

  implicit none
  private
  public:: geom_to_gph,lat_geod_to_geoc,gph_to_geom
  ! Geodetic parameters used by more than one routine. 
  ! Some of these should be public if they don't interfere with anyone elses
  ! code
  real(kind=r8),private,parameter::   g0 = 9.80665_r8,&
       dpi                       =  3.141592653589793115997963_r8,&
       earth_major_axis         = 6378137.0_r8, & 
       earth_minor_axis         = 6356752.3141_r8, &
       earth_axis_ratio_squared = (earth_major_axis*earth_major_axis / &
       (earth_minor_axis*earth_minor_axis))

  character(len=*),parameter::ModuleName=&
       "$RCSfile$"

contains

  function geom_to_gph(geom_height,lat_geoc,lat_geod) result (gph)
    ! Converts geometric height above the geoid to  geopotential height
    ! This applies Equation 40 in "Theoretical basis for EOS-MLS geopotential 
    ! and scan calculations" by N. J. Livesey
    ! If I have got this right, the difference between the two heights is 
    ! piffling at ground level, going to about 1km at 80km. 

    !----Argument----!
    real(kind=r8),intent(in)::geom_height
    real(kind=r8),intent(in),optional::lat_geoc,lat_geod
    !----Function result----!
    real(kind=r8)::gph
    !---Local variabules---!
    real(kind=r8),parameter::omega=7.292115e-5_r8,      gm=3.986005e14_r8, &
         j2=              0.0010826256,       j4=-0.0000023709122_r8
    real(kind=r8)::p2,p4,sinlat,h_inf,ahrat,coslat,h0,hc,lat
    !---------Executable statements-------!

    if(present(lat_geoc) .and. present(lat_geod)) then
       call MLSMessage(MLSMSG_Warning,ModuleName,&
            "in function geom_to_gph,"// &
            "Both geometric and geodetic lats supplied: using geometric")
       lat=lat_geoc
    else if(present(lat_geoc) )then
       lat=lat_geoc
    else if(present(lat_geod)) then
       lat=lat_geod_to_geoc(lat_geod)
    else
       call MLSMessage(MLSMSG_Warning,ModuleName,&
            "in function geom_to_gph, Called with no latitude - assuming 45deg N")
       lat=45.0_r8
    endif

    sinlat=sin(lat_geoc*dpi/180.0_r8)
    coslat=cos(lat_geoc*dpi/180.0_r8)

    ! Calculate radius of earth at given latitude

    h0=sqrt( (earth_major_axis*coslat)**2 + (earth_minor_axis*sinlat)**2)
    hc=h0+geom_height
    p2=(3*sinlat*sinlat -1)/2
    p4=(35*(sinlat**4) - 30*sinlat*sinlat +3)/8

    ahrat=earth_major_axis/hc
    ahrat=ahrat*ahrat

    h_inf=(gm/(g0*hc))*(1-j2*p2*ahrat - j4*p4*ahrat*ahrat) + &
         (omega*omega*hc*hc*coslat*coslat)/(2*g0)

    gph=6387182.265_r8-h_inf

  end function geom_to_gph

  function lat_geod_to_geoc(lat_geod) result(lat_geoc)
    ! Converts geodetic latitude to geometric (geocentric) latitude
    ! If I have got this right, the difference between the two latitudes 
    ! is 0 at equator and pole with a maximum of 0.2 degrees at midlatitudes
    real(kind=r8),intent(in)::lat_geod
    real(kind=r8)::lat_geoc
    lat_geoc=(180.0_r8/dpi)*atan(tan(lat_geod*dpi/180.0_r8)&
         /earth_axis_ratio_squared)
  end function lat_geod_to_geoc

  function gph_to_geom(gph,lat_geoc,lat_geod) result (geom_height)
    ! This converts Geopotential height (GPH) back to geometric height
    ! We use a yucky secant-method iteration mainly out of idleness -- 
    ! we can't be bothered to invert the function analytically. Of course, 
    ! the present hack has the advantage that gph_to_geom and geom_to_gph 
    ! are always consistent with each other.
    ! Tests suggest that the secant method converges after either 3 or 4
    ! iterations over the 0-100km range.

    !----Argument----!
    real(kind=r8),intent(in)::gph
    real(kind=r8),intent(in),optional::lat_geoc,lat_geod
    !----Function result----!
    real(kind=r8)::geom_height
    !---Local variabules---!
    real(kind=r8),parameter:: accuracy=0.001 
    real(kind=r8)::geom1,geom2,newgeom,dgph1,dgph2,lat
    !---------Executable statements-------!

    if(present(lat_geoc) .and. present(lat_geod)) then
       call MLSMessage(MLSMSG_Warning,ModuleName,&
            "in function gph_to_geom,"// &
            "Both geometric and geodetic lats supplied: using geometric")
       lat=lat_geoc
    else if(present(lat_geoc) )then
       lat=lat_geoc
    else if(present(lat_geod)) then
       lat=lat_geod_to_geoc(lat_geod)
    else
       call MLSMessage(MLSMSG_Warning,ModuleName,&
            "in function gph_to_geom, heights_module"// &
            "Called with no latitude - assuming 45deg N")
       lat=45.0_r8
    endif

    geom1=gph
    geom2=gph+1.0
    dgph1=geom_to_gph(geom1,lat)-gph
    secantloop:do 
      dgph2=geom_to_gph(geom2,lat)-gph
      newgeom= geom1 - dgph1*(geom2-geom1)/(dgph2-dgph1)
      !     print*,geom1,geom2,newgeom
      if (abs(newgeom-geom2) < accuracy) then
         exit secantloop
      endif
      geom1=geom2
      dgph1=dgph2
      geom2=newgeom

    enddo secantloop
    geom_height=newgeom

  end function gph_to_geom

end module heights_module
