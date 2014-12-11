! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module uars_to_emls_oa_m

  implicit NONE
  private

  public :: UARS_to_EMLS_OA

  logical, parameter, private :: Deebug = .false.

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public procedures ===================================

  subroutine uars_to_emls_oa ( n_days, limb_oa, emls_oa, Angle )

    use Calc_GeodAngle_m, only: Calc_GeodAngle
    use Constants, only: Deg2Rad, Rad2Deg, Pi
    use Cross_m, only: Cross
    use Dates_Module, only: SecondsBetween2UTCs
    use Geometry, only: GeocToGeodLat, Omega => W
    use MLSKinds, only: R8
    use OA_File_Contents, only: emls_oa_t
    use Rad_File_Contents, only: limb_oa_t
    use Rotation_m, only: Rotate_3d
    use SDPToolkit, only: PGS_S_SUCCESS

    ! Args
    integer, intent(in) :: n_days
    type(limb_oa_t), intent(in) :: limb_oa
    type(emls_oa_t), intent(inout) :: emls_oa
    real, intent(in) :: Angle ! Between first MIFs of consecutive MAFs, radians

    integer, parameter :: noMIFs = 32
    integer :: n_years
    integer :: i, m, mif_range(2), ms, hrs, mins, yrdoy, stat, mif1_ms, year
    real(r8) :: forgedsecs
    real(r8) :: secs, cosalf, sinalf, coseps, sineps, cosphi, sinphi, costheta, &
         sintheta
    real(r8) :: fcol1(3), fcol2(3), fcol3(3), fmaf(3,3), flosf(3,3), fecr(3), &
         ft(3), ipf2eci(3,3,noMIFs), ecrtofov(3,3,noMIFs), phi, theta
    integer, parameter :: moon_in_fov = int(z'04000400')  ! bits to set for EMLS
    real, parameter :: orbincl = 57.0  ! orbit incline
    real, parameter :: mif_inc = 65.536 / noMIFs   ! increment per MIF in seconds
    real(r8) :: eciV(6,noMIFs), ecrV(6,noMIFs), los_vec(3,noMIFs), &
         offsets(noMIFs), orb_norm(3), posecr(3), rotated(3), sectai93, &
         tngtvel(3)
    real(r8) :: Delta, xi, cosxi, sinxi
    character(len=25) :: asciiutc, forgedasciiutc, err, msg
    character (len=*), parameter :: earthellipstag = 'WGS84', nl = char(10)
    integer, external :: pgs_td_utctotai, pgs_csc_geotoecr, pgs_smf_getmsg, &
         pgs_csc_ecrtoeci, pgs_csc_ecitoecr

  ! Initialize offsets per mif in seconds:

    do m = 1, noMIFs
      offsets(m) = (m-1) * mif_inc
    end do

    emls_oa%sc_geocalt = limb_oa%sat_gcrad * 1000.0    ! array copy
   !emls_oa%sc_orbincl = orbincl                       ! array copy

    emls_oa%geod_alt = limb_oa%tngt_geod_alt * 1000.0  ! meters
    emls_oa%geod_lat = limb_oa%tngt_geod_lat
    emls_oa%lon = limb_oa%tngt_long

  ! bright object status:

    emls_oa%bo_stat = 0  !clear to start
    if (limb_oa%ptg_fov_bo_diag_mmif_num /= 0) then  ! moon in FOV
      mif_range(1) = limb_oa%ptg_fov_bo_diag_mmif_fst - 1
      mif_range(2) = limb_oa%ptg_fov_bo_diag_mmif_lst - 1
      emls_oa%bo_stat(mif_range(1):mif_range(2)) = moon_in_fov
    end if

  ! Use same value for all MIFs:

    emls_oa%solartime = limb_oa%ref_solar_time
    emls_oa%solarzenith = limb_oa%ref_solar_zen
    do m = 1, noMIFs   ! same for every MIF
      emls_oa%sc_ypr(:,m) = limb_oa%ypr
      emls_oa%sc_ypr_rate(:,m) = limb_oa%ypr_rate
    end do

  ! Calculate TAI time at start of frame:

    mif1_ms = (limb_oa%ref_mmif - 1) * 2048   ! millisecs of MIF 1
    yrdoy = limb_oa%ref_time(1)  ! - n_days      ! year plus day of year
    ms = limb_oa%ref_time(2)     ! millisecs of day for ref_mmif
    ms = ms - mif1_ms            ! millisecs at start of MAF
    if (ms < 0) then             ! in previous day!
      ms = ms + 86400000
      yrdoy = yrdoy - 1
    end if
    hrs = mod (ms/3600000, 24)
    mins = mod (ms/60000, 60)
    secs = amod (float(ms), 60000.0) /1000.0
    year = yrdoy / 1000
    if (year > 99) then      ! take care of years 2000 or greater
       year = mod(year, 100) + 2000
    else                     ! years in the 1900's
       year = year + 1900
    end if

    write (asciiutc, fmt= &
     '(i4, "-", i3.3, "T", i2.2, ":", i2.2, ":", f9.6, "Z", TL10, i2.2)') &
     year, mod(yrdoy, 1000), hrs, mins, secs, int(secs)  ! force leading 0's
    if ( deebug ) print *, 'asciiutc: ', trim(asciiutc)
    stat = pgs_td_utctotai (asciiutc, sectai93)

    emls_oa%MAFStartTimeTAI= sectai93
    ! Adjust for possible backdating
    if ( n_days > 999 ) then
      ! Account for our convention that a "year" has 1000 days
      n_years = (n_days/1000)
      write (forgedasciiutc, fmt= &
       '(i4, "-", i3.3, "T", i2.2, ":", i2.2, ":", f9.6, "Z", TL10, i2.2)') &
       year-n_years, mod(yrdoy, 1000), hrs, mins, secs, int(secs)  ! force leading 0's
      forgedsecs = secondsbetween2utcs ( forgedasciiutc, asciiutc )
      emls_oa%MAFStartTimeTAI = sectai93 - forgedsecs
    elseif ( n_days > 0 ) then
      forgedsecs = n_days * 24 * 3600
      emls_oa%MAFStartTimeTAI = sectai93 - forgedsecs
    else
      forgedsecs = 0.d0
    end if

  ! spacecraft MIF TAI:

    do m = 1, noMIFs
       emls_oa%sc_MIF_TAI(m) = sectai93 + (m-1) * mif_inc - forgedsecs
    end do

  ! spacecraft ECR positions, at first the same for all MIFs:

    stat = pgs_csc_geotoecr ( limb_oa%sat_long*Deg2Rad, &
                              limb_oa%sat_geod_lat*Deg2Rad, &
                              limb_oa%sat_geod_alt*1000.0d0, &
                              earthellipstag, posecr )
    do m = 1, noMIFs
      emls_oa%sc_ECR(:,m) = posecr ! Same ECR for every MIF
    end do
    emls_oa%sc_geoclat(:) = ATAN2 ( posecr(3), &
                                  & SQRT (posecr(1)**2 + posecr(2)**2)) * rad2deg

  ! Convert spacecraft ECR position to ECI:

    ecrV(1:3,:) = emls_oa%sc_ECR
    ecrV(4:6,:) = 0.0

    stat = pgs_csc_ecrtoeci ( noMIFs, asciiUTC, offsets, ecrV, eciV )
    emls_oa%sc_ECI = eciV(1:3,:)

  ! Get SC velocity (in meters/sec) in ECI:

    do m = 1, noMIFs
       eciV(4:6,m) = limb_oa%sat_vel * 1000.0 ! only need one set of values
    end do

  ! Convert back to ECR (to get SC velocity in ECR):

    stat = pgs_csc_ecitoecr ( noMIFs, asciiUTC, offsets, eciV, ecrV )

  ! Rotate ECR positions and velocities after the first MIF about the normal
  ! to the orbit plane by Angle/noMIFs * ( MIF# - 1 ).

    orb_norm = cross ( ecrV(1:3,1), ecrV(4:6,1), norm=.true. )
    delta = angle / noMIFs
    do m = 2, noMIFs
      call rotate_3d ( ecrV(1:3,m), delta * ( m-1 ), orb_norm, rotated )
      ecrV(1:3,m) = rotated
      call rotate_3d ( ecrV(4:6,m), delta * ( m-1 ), orb_norm, rotated )
      ecrV(4:6,m) = rotated
    end do
    emls_oa%sc_ECR = ecrV(1:3,:)
    emls_oa%sc_VelECR = ecrV(4:6,:)

  ! Now convert rotated positions and velocities back to ECI.
    stat = pgs_csc_ecrtoeci ( noMIFs, asciiUTC, offsets, ecrV, eciV )
    emls_oa%sc_ECI = eciV(1:3,:)
    emls_oa%sc_VelECI = eciV(4:6,:)

  ! Spacecraft geolocation
    emls_oa%sc_lon = atan2 ( ecrV(2,:), ecrV(1,:) ) * rad2deg
    emls_oa%sc_geoclat = atan2 ( ecrV(3,:), sqrt(ecrV(1,:)**2 + ecrV(2,:)**2) ) * &
                       & rad2deg

  ! Remap sc/Lon
    where (emls_oa%sc_lon > 180.0)
       emls_oa%sc_lon = emls_oa%sc_lon - 360.0
    end where
  ! Get SC geodetic latitude from its geocentric latitude  
    emls_oa%sc_geodlat = geocToGeodLat ( emls_oa%sc_geoclat )
  ! SC Geodetic altitude is the same for every MIF; hope that's OK
    emls_oa%sc_geodalt = limb_oa%sat_geod_alt * 1000.0  ! meters

  ! GHz ECR and geoc_alt and geoc_lat:

    do m = 1, noMIFs
      stat = pgs_csc_geotoecr (limb_oa%tngt_long(m)*Deg2Rad, &
           limb_oa%tngt_geod_lat(m)*Deg2Rad, limb_oa%tngt_geod_alt(m)*1000.0d0, &
           earthellipstag, posecr) !emls_oa%ECR(:,m))
      emls_oa%ECR(:,m) = posecr
      emls_oa%geoc_alt(m) = SQRT (emls_oa%ECR(1,m)**2 + emls_oa%ECR(2,m)**2 + &
           emls_oa%ECR(3,m)**2)
      emls_oa%geoc_lat(m) = Rad2Deg * ATAN2 (posecr(3) , &
           SQRT (posecr(1)**2 + posecr(2)**2))
    end do

  ! Calculate ECRtoFOV:

    do m = 1, noMIFs
      emls_oa%azimAngle(m) = limb_oa%ptg_fov_azim_offset(m) + 90.0
      cosalf = COS (emls_oa%azimAngle(m)*Deg2Rad)
      sinalf = SIN (emls_oa%azimAngle(m)*Deg2Rad)

      emls_oa%scanAngle(m) = limb_oa%ptg_fov_elev_offset(m) + 23.3
      coseps = COS (emls_oa%scanAngle(m)*Deg2Rad)
      sineps = SIN (emls_oa%scanAngle(m)*Deg2Rad)

      ! Changed to account for rotation of axes during ascan
      xi = 113.6 + emls_oa%scanAngle(m) - 23.3 ! p.30 of UARS MLS cal.report, B1 is pointing reference
      cosxi = COS (xi*Deg2Rad)
      sinxi = SIN (xi*Deg2Rad)

      fcol1 = (/ cosalf*sineps*sinxi+sinalf*cosxi, sinalf*sineps*sinxi-cosalf*cosxi, -coseps*sinxi /)
      fcol2 = (/ cosalf*sineps*cosxi-sinalf*sinxi, sinalf*sineps*cosxi+cosalf*sinxi, -coseps*cosxi /)
      fcol3 = (/ cosalf*coseps, sinalf*coseps, sineps /)

      fmaf = RESHAPE((/ fcol1, fcol2, fcol3 /), (/ 3, 3 /))
      ipf2eci(:,:,m) = TRANSPOSE (MATMUL (limb_oa%trans_inst2eci, fmaf))
      flosf =  TRANSPOSE (RESHAPE ((/ (/ 0.0d0, 0.0d0, -1.0d0 /), &
           (/ -sinalf, cosalf, 0.0d0 /), &
           (/ cosalf, sinalf, 0.0d0 /) /), (/ 3, 3 /)))
      fmaf = TRANSPOSE (fmaf)
      do i = 1, 3
         eciV(1:3,m) = ipf2eci(i,:,m)
         eciV(4:6,m) = 0.0d0
         stat = pgs_csc_ecitoecr (1, asciiUTC, offsets(m), eciV(:,m), ecrV(:,m))
         ecrtofov(:,i,m) = ecrV(1:3,m)
      end do
      ecrtofov(:,:,m) = TRANSPOSE (ecrtofov(:,:,m))    ! Put in "right" order
      emls_oa%ECRtoFOV(:,m) = RESHAPE (ecrtofov(:,:,m), (/ 9 /))  ! make 2-d

    ! Calculate GeodAngle (using REC's algorithm):

      call Calc_GeodAngle ( m, fmaf, flosf, emls_oa )

    end do

  ! Convert GHz ECR to ECI:

    ecrV(1:3,:) = emls_oa%ECR
    ecrV(4:6,:) = 0.0

    stat = pgs_csc_ecrtoeci (noMIFs, asciiUTC, offsets, ecrV, eciV)
    emls_oa%ECI = eciV(1:3,:)

  ! Calculate LosAngle

    do m = 1, noMIFs
      phi = emls_oa%lon(m) * Deg2Rad
      cosphi = cos (phi)
      sinphi = sin (phi)
      theta = (90.0 - emls_oa%geod_lat(m)) * Deg2Rad
      costheta = cos (theta)
      sintheta = sin (theta)

      fcol1 = (/ -sinphi, cosphi, 0.0d0 /)
      fcol2 = (/ -costheta*cosphi, -costheta*sinphi, sintheta /)
      fcol3 = (/ sintheta*cosphi, sintheta*sinphi, costheta /)
      fmaf = RESHAPE((/ fcol1, fcol2, fcol3 /), (/ 3, 3 /))
      fecr = [emls_oa%ecrtofov(3,m), emls_oa%ecrtofov(6,m), emls_oa%ecrtofov(9,m)]
      ft = matmul (fecr, fmaf)
      emls_oa%LosAngle(m) = mod ( 90.0d0-atan2 (ft(2), ft(1))*Rad2Deg+360.0, &
                                & 360.0d0 )

    end do

  ! Calculate LosVel

    do m = 1, noMIFs
      tngtVel = omega * (/ -emls_oa%ECI(2,m), emls_oa%ECI(1,m), 0.0d0 /)
      los_vec(:,m) = emls_oa%ECI(:,m) - emls_oa%sc_ECI(:,m)
      los_vec(:,m) = los_vec(:,m) / SQRT (SUM (los_vec(:,m)**2))
      emls_oa%LosVel(m) = DOT_PRODUCT (tngtVel, los_vec(:,m)) - &
           DOT_PRODUCT (limb_oa%sat_vel*1000.0d0, los_vec(:,m))
    end do

  ! Calculate OrbY

    do m = 1, noMIFs
      orb_norm = cross ( emls_oa%sc_ECI(:,m), emls_oa%sc_VelECI(:,m), norm=.true. )
      emls_oa%OrbY(m) = DOT_PRODUCT (-emls_oa%ECI(:,m), orb_norm)
    end do

    ! Remap GHz/Lon
    where (emls_oa%lon > 180.0) emls_oa%lon = emls_oa%lon - 360.0

  !!$write (30, *) limb_oa%ptg_fov_azim_offset
  !!$write (31, *) limb_oa%ptg_fov_elev_offset
  !!$write (32, *) limb_oa%sat_geod_alt
  !!$write (33, *) limb_oa%sat_geod_lat
  !!$write (34, *) limb_oa%sat_long
  !!$write (35, *) limb_oa%tngt_geod_alt
  !!$write (36, *) limb_oa%tngt_geod_lat
  !!$write (37, *) limb_oa%tngt_long
  !!$write (40, *) emls_oa%ECR
  !!$write (81, *) emls_oa%sc_ECI
  !!$write (82, *) emls_oa%sc_VelECI
  !!$write (83, *) emls_oa%ECI

  !write (3x, *) limb_oa%

  end subroutine uars_to_emls_oa

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module uars_to_emls_oa_m

! $Log$
