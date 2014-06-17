subroutine uars_to_emls_oa

  USE rad_file_contents, ONLY: limb_oa
  USE oa_file_contents, ONLY: emls_oa
  USE SDPToolkit, ONLY: PGS_S_SUCCESS
  USE Constants, ONLY: Deg2Rad, Rad2Deg, Pi
  USE Geometry, ONLY: Omega => W

  implicit none

  integer, parameter :: nomifs = 32
  integer :: i, m, mif_range(2), ms, hrs, mins, yrdoy, stat, mif1_ms, year
  real*8 :: secs, cosalf, sinalf, coseps, sineps, cosphi, sinphi, costheta, &
       sintheta
  real*8 :: fcol1(3), fcol2(3), fcol3(3), fmaf(3,3), flosf(3,3), fecr(3), &
       ft(3), ipf2eci(3,3,nomifs), ecrtofov(3,3,nomifs), phi, theta
  integer, parameter :: moon_in_fov = x'04000400'  ! bits to set for EMLS
  real, parameter :: orbincl = 57.0  ! orbit incline
  real, parameter :: mif_inc = 65.536 / 32   ! increment per MIF in seconds
  REAL*8 :: sectai93, eciV(6,nomifs), ecrV(6,nomifs), offsets(nomifs), &
       tngtvel(3), los_vec(3,nomifs), posecr(3), rv(3), orb_norm(3)
  character*25 asciiutc, err, msg
  CHARACTER (LEN=*), PARAMETER :: earthellipstag = 'WGS84', nl = char(10)
  integer, external :: pgs_td_utctotai, pgs_csc_geotoecr, pgs_smf_getmsg, &
       pgs_csc_ecrtoeci, pgs_csc_ecitoecr

! Initialize offsets per mif in seconds:

  do m = 1, nomifs
     offsets(m) = (m-1) * mif_inc
  enddo

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
  endif

! Use same value for all MIFs:

  emls_oa%solartime = limb_oa%ref_solar_time
  emls_oa%solarzenith = limb_oa%ref_solar_zen
  emls_oa%sc_lon = limb_oa%sat_long
  where (emls_oa%sc_lon > 180.0)
     emls_oa%sc_lon = emls_oa%sc_lon - 360.0
  end where
  emls_oa%sc_geodlat = limb_oa%sat_geod_lat
  emls_oa%sc_geodalt = limb_oa%sat_geod_alt * 1000.0  ! meters
  do m = 1, nomifs   ! same for every MIF
     emls_oa%sc_ypr(:,m) = limb_oa%ypr
     emls_oa%sc_ypr_rate(:,m) = limb_oa%ypr_rate
  enddo

! Calculate TAI time at start of frame:

  mif1_ms = (limb_oa%ref_mmif - 1) * 2048   ! millisecs of MIF 1
  yrdoy = limb_oa%ref_time(1)  ! year plus day of year
  ms = limb_oa%ref_time(2)     ! millisecs of day for ref_mmif
  ms = ms - mif1_ms            ! millisecs at start of MAF
  if (ms < 0) then             ! in previous day!
     ms = ms + 86400000
     yrdoy = yrdoy - 1
  endif
  hrs = mod (ms/3600000, 24)
  mins = mod (ms/60000, 60)
  secs = amod (float(ms), 60000.0) /1000.0
  year = yrdoy / 1000
  if (year > 99) then      ! take care of years 2000 or greater
     year = mod(year, 100) + 2000
  else                     ! years in the 1900's
     year = year + 1900
  endif

  write (asciiutc, fmt= &
   '(i4, "-", i3.3, "T", i2.2, ":", i2.2, ":", f9.6, "Z", TL10, i2.2)') &
   year, mod(yrdoy, 1000), hrs, mins, secs, int(secs)  ! force leading 0's
  stat = pgs_td_utctotai (asciiutc, sectai93)

  emls_oa%MAFStartTimeTAI= sectai93

! spacecraft MIF TAI:

  do m = 1, nomifs
     emls_oa%sc_MIF_TAI(m) = sectai93 + (m-1) * mif_inc
  enddo

! spacecraft ECR:

  do m = 1, nomifs
     stat = pgs_csc_geotoecr (limb_oa%sat_long*Deg2Rad, &
          limb_oa%sat_geod_lat*Deg2Rad, limb_oa%sat_geod_alt*1000.0d0, &
          earthellipstag, posecr) !emls_oa%sc_ECR(:,m))
     emls_oa%sc_ECR(:,m) = posecr
     emls_oa%sc_geoclat(m) = Rad2Deg * ATAN (posecr(3) / &
          SQRT (posecr(1)**2 + posecr(2)**2))
  enddo

! Convert spacecraft ECR to ECI:

  ecrV(1:3,:) = emls_oa%sc_ECR
  ecrV(4:6,:) = 0.0

  stat = pgs_csc_ecrtoeci (nomifs, asciiUTC, offsets, ecrV, eciV)
  emls_oa%sc_ECI = eciV(1:3,:)

! Convert back to ECR and set velocity (in meters/sec):

  do m = 1, nomifs
     eciV(4:6,m) = limb_oa%sat_vel * 1000.0 ! only need one set of values
  enddo

  stat = pgs_csc_ecitoecr (nomifs, asciiUTC, offsets, eciV, ecrV)
  emls_oa%sc_VelECI = eciV(4:6,:)
  emls_oa%sc_VelECR = ecrV(4:6,:)

! GHz ECR and geoc_alt and geoc_lat:

  do m = 1, nomifs
     stat = pgs_csc_geotoecr (limb_oa%tngt_long(m)*Deg2Rad, &
          limb_oa%tngt_geod_lat(m)*Deg2Rad, limb_oa%tngt_geod_alt(m)*1000.0d0, &
          earthellipstag, posecr) !emls_oa%ECR(:,m))
     emls_oa%ECR(:,m) = posecr
     emls_oa%geoc_alt(m) = SQRT (emls_oa%ECR(1,m)**2 + emls_oa%ECR(2,m)**2 + &
          emls_oa%ECR(3,m)**2)
     emls_oa%geoc_lat(m) = Rad2Deg * ATAN (posecr(3) / &
          SQRT (posecr(1)**2 + posecr(2)**2))
  enddo

! Calculate ECRtoFOV:

  do m = 1, nomifs
     emls_oa%azimAngle(m) = limb_oa%ptg_fov_azim_offset(m) + 90.0
     cosalf = COS (emls_oa%azimAngle(m)*Deg2Rad)
     sinalf = SIN (emls_oa%azimAngle(m)*Deg2Rad)

     emls_oa%scanAngle(m) = limb_oa%ptg_fov_elev_offset(m) + 23.3
     coseps = COS (emls_oa%scanAngle(m)*Deg2Rad)
     sineps = SIN (emls_oa%scanAngle(m)*Deg2Rad)

     fcol1 = (/ cosalf*sineps, sinalf*sineps, -coseps /)
     fcol2 = (/ -sinalf, cosalf, 0.0d0 /)
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
     enddo
     ecrtofov(:,:,m) = TRANSPOSE (ecrtofov(:,:,m))    ! Put in "right" order
     emls_oa%ECRtoFOV(:,m) = RESHAPE (ecrtofov(:,:,m), (/ 9 /))  ! make 2-d

! Calculate GeodAngle (using REC's algorithm):

     call Calc_GeodAngle (m, fmaf, flosf)

  enddo

! Convert GHz ECR to ECI:

  ecrV(1:3,:) = emls_oa%ECR
  ecrV(4:6,:) = 0.0

  stat = pgs_csc_ecrtoeci (nomifs, asciiUTC, offsets, ecrV, eciV)
  emls_oa%ECI = eciV(1:3,:)

! Calculate LosAngle

  DO m = 1, nomifs
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
     emls_oa%LosAngle(m) = mod (90.0d00-atan2 (ft(2), ft(1))*Rad2Deg+360.0 , 360.0)

  ENDDO

! Calculate LosVel

  DO m = 1, nomifs
     tngtVel = omega * (/ -emls_oa%ECI(2,m), emls_oa%ECI(1,m), 0.0d0 /)
     los_vec(:,m) = emls_oa%ECI(:,m) - emls_oa%sc_ECI(:,m)
     los_vec(:,m) = los_vec(:,m) / SQRT (SUM (los_vec(:,m)**2))
     emls_oa%LosVel(m) = DOT_PRODUCT (tngtVel, los_vec(:,m)) - &
          DOT_PRODUCT (limb_oa%sat_vel*1000.0d0, los_vec(:,m))
  ENDDO

! Calculate OrbY

  DO m = 1, nomifs
     call cross_product (emls_oa%sc_ECI(:,m), emls_oa%sc_VelECI(:,m), rv)
     orb_norm = rv / sqrt (rv(1)*rv(1) + rv(2)*rv(2) + rv(3)*rv(3))
     emls_oa%OrbY(m) = DOT_PRODUCT (-emls_oa%ECI(:,m), orb_norm)
  ENDDO

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
