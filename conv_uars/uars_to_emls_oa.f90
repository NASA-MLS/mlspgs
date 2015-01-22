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

  public :: Frame_TAI, UARS_to_EMLS_OA

  character (len=*), parameter, public :: EarthEllipseTag = 'WGS84'
  integer, parameter, public :: noMIFs = 32

  logical, parameter, private :: Deebug = .false.

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public procedures ===================================

  subroutine uars_to_emls_oa ( n_days, limb_oa, emls_oa, dV, dVN, dVV, Sat_Vel, &
    & In_posECR, Vel_ECI, dVI, dVIN, In_AsciiUTC )

    use Calc_GeodAngle_m, only: Calc_GeodAngle
    use Constants, only: Deg2Rad, Rad2Deg, Pi
    use Cross_m, only: Cross
    use Geometry, only: GeocToGeodLat, Omega => W, XYZ_to_Geod
    use MLSKinds, only: R8
    use OA_File_Contents, only: emls_oa_t
    use PGS_Interfaces, only: pgs_csc_ecitoecr, pgs_csc_ecrtoeci, pgs_csc_geotoecr
    use Rad_File_Contents, only: limb_oa_t
    use Rotation_m, only: Rotate_3d
    use SDPToolkit, only: PGS_S_SUCCESS

    ! Args
    integer, intent(in) :: n_days
    type(limb_oa_t), intent(in) :: limb_oa
    type(emls_oa_t), intent(inout) :: emls_oa
    real(r8), intent(in) :: dV     ! delta V angle per MIF, Radians, ECR
    real(r8), intent(in) :: dVN    ! relative delta V length per MIF (tiny) km/s
    real(r8), intent(in) :: dVV(3) ! delta V per MAF, in ECR, km/s
    real(r8), intent(in) :: Sat_Vel(3,2)   ! Satellite velocities, in ECR, m/s
    real(r8), intent(in) :: In_posECR(3,2) ! Satellite positions, ECR, meters
    real(r8), intent(in) :: Vel_ECI(3,2)   ! Satellite velocities, ECI, km/s
    real(r8), intent(in) :: dVI    ! delta V angle per MIF, ECI
    real(r8), intent(in) :: dVIN   ! relative delta V length per MIF, ECI
    character(*), intent(in) :: In_AsciiUTC(2)

    integer :: i, m, mif_range(2), stat
    real(r8) :: forgedsecs
    real(r8) :: cosalf, sinalf, coseps, sineps, cosphi, sinphi, costheta, &
         sintheta
    real(r8) :: fmaf(3,3), flosf(3,3), fecr(3), &
         ft(3), ipf2eci(3,3), ecrtofov(3,3,noMIFs), phi, theta
    integer, parameter :: moon_in_fov = int(z'04000400')  ! bits to set for EMLS
    real, parameter :: orbincl = 57.0  ! orbit incline
    real, parameter :: mif_inc = 65.536 / noMIFs   ! increment per MIF in seconds
    real(r8), parameter :: Offsets(noMIFs) = [ ( (m-1) * mif_inc, m = 1, noMIFs ) ]
    real(r8) :: eciV(6,noMIFs), ecrV(6,noMIFs), Geod(3), los_vec(3), &
         orb_norm(3), posecr(3), sectai93, tngtvel(3)
    real(r8) :: Delta, xi, cosxi, sinxi
    real(r8) :: dVVR(3) ! delta V, in ECR, per MIF
    real(r8) :: PosECI(3,2)
    character(len=25) :: asciiutc

    emls_oa%sc_geocalt = limb_oa%sat_gcrad * 1000.0   ! MIF resolved, meters
   !emls_oa%sc_orbincl = orbincl                      ! MIF resolved

    emls_oa%geod_alt = limb_oa%tngt_geod_alt * 1000.0 ! MIF resolved, meters
    emls_oa%geod_lat = limb_oa%tngt_geod_lat          ! MIF resolved
    emls_oa%lon = limb_oa%tngt_long                   ! MIF resolved

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

    call frame_TAI ( limb_oa, secTAI93, asciiUTC, n_days, forgedSecs )

    emls_oa%MAFStartTimeTAI= sectai93 - forgedSecs

  ! spacecraft MIF TAI:

    do m = 1, noMIFs
       emls_oa%sc_MIF_TAI(m) = sectai93 + (m-1) * mif_inc - forgedsecs
    end do

  ! spacecraft ECR position at the first MIF:

    stat = pgs_csc_geotoecr ( limb_oa%sat_long*Deg2Rad, &
                              limb_oa%sat_geod_lat*Deg2Rad, &
                              limb_oa%sat_geod_alt*1000.0d0, &
                              earthEllipseTag, posecr )

  ! Convert spacecraft ECR position in meters (and zero velocity) at the first
  ! MIF to ECI:

    ecrV(1:3,1) = posECR
    ecrV(4:6,1) = 0
write ( 42, '(2x,2a)' ) '-------- In_PosECR -------- --------- PosECR ----------', &
' --------- Vel ECI ---------'
write ( 42, '(2x,5(3f9.2:","))' ) in_posECR(:,1)/1000, posECR/1000, limb_oa%sat_vel*1000
    stat = pgs_csc_ecrtoeci ( 1, asciiUTC, offsets, ecrV, eciV )

  ! Get SC velocity (in meters/sec) at the first MIF in ECI:

    eciV(4:6,1) = limb_oa%sat_vel * 1000.0

  ! Convert position and velocity back to ECR (to get SC velocity in ECR):

    stat = pgs_csc_ecitoecr ( 1, asciiUTC, offsets, eciV, ecrV )
write ( 42,'(2x,2a)' ) ' ------------------------ ECIV ------------------------', &
& ' ------------------------ ECRV ------------------------'
write ( 42, '(2x,4(3f9.2:","))' ) eciV(1:3,1)/1000,eciV(4:6,1), ecrV(1:3,1)/1000,ecrV(4:6,1)
write ( 42, '(2x,3a)' ) '------- Vel ECR here ------ --------- Vel ECR 1 -------', &
& ' --------- Vel ECR 2 ------- --------- dVV ECR ---------', &
& ' -------- Delta Vel -------- --------- Vel ECI ---------'
write ( 42, '(2x,6(3f9.2:","))' ) ecrV(4:6,1), sat_vel, dVV, sat_vel(:,2)-sat_vel(:,1), limb_oa%sat_vel*1000
    dVVR = dVV / noMIFs

  ! Normal to the orbit plane in ECR, at the first MIF, is the cross
  ! product of the spacecraft position and velocity vectors
    orb_norm = cross ( ecrV(1:3,1), ecrV(4:6,1) )

  ! Normal to the orbit plane in ECI, at the first MIF, is the cross
  ! product of the spacecraft position and velocity vectors
    orb_norm = cross ( eciV(1:3,1), eciV(4:6,1) )

  !{ Angle between consecutive MIFs = $|V|\, T\, \sin(\theta) / |H|$, where
  !  $H$ is the spacecraft position vector, $V$ is the spacecraft velocity
  !  vector, $\sin\theta = |H \times V|\, /\, |H|\, |V|$, and $T$ is the MIF
  !  duration.
    delta = mif_inc * norm2 ( orb_norm ) / emls_oa%sc_geocalt(1)**2

!   ! Rotate ECR positions and velocities after the first MIF about the normal
!   ! to the orbit plane by Delta * ( MIF# - 1 ).  Adjust the velocity vector
!   ! length at each MIF according to the total relative change during the MAF.
! write ( 42, '(2x,2a,3x,a)' ) '------- Vel ECR here ------ --------- dVV ECR ---------', &
! & ' --- |Vel|', asciiUTC
! write ( 42, '(i2, 5(3f9.2:","))' ) 1, ecrV(4:6,1), dVVR*noMIFs, norm2(ecrV(4:6,1))
!     do m = 2, noMIFs
!       call rotate_3d ( ecrV(1:3,1), delta * ( m-1 ), orb_norm, ecrV(1:3,m) )
!       call rotate_3d ( ecrV(4:6,1), dV * ( m-1 ), orb_norm, ecrV(4:6,m) )
!       ecrV(4:6,m) = ecrV(4:6,m) * ( ( m-1) * dVN + 1.0 )
! write ( 42, '(i2, 5(3f9.2:","))' ) m, ecrV(4:6,m), dVVR*(m-1), norm2(ecrV(4:6,m))
!     end do
!     ! Set the length of the SC ECR position vector to the SC Geocentric Altitude.
!     ! The ECR length gotten using PGS_CSC_GeoToECR's computation of SC position
!     ! from geodetic coordinates is different from the SC Geocentric Altitude by
!     ! +/- 300 meters, with a period of about 120 MAFs.  We choose to believe the
!     ! SC Geocentric Altitude.
!     do m = 1, noMIFs
!       ecrV(1:3,m) = ( ecrV(1:3,m) / norm2(ecrV(1:3,m)) ) * emls_oa%sc_geocalt(m)
!     end do

  ! Rotate ECI positions and velocities after the first MIF about the normal
  ! to the orbit plane by Delta * ( MIF# - 1 ).  Adjust the velocity vector
  ! length at each MIF according to the total relative change during the MAF.
write ( 42, '(2x,2a,3x,a)' ) '------- Vel ECI here ------ --------- dVV ECI ---------', &
& ' --- |Vel|', asciiUTC
write ( 42, '(i2, 5(3f9.2:","))' ) 1, eciV(4:6,1), dVV*noMIFs, norm2(eciV(4:6,1))
    do m = 2, noMIFs
      call rotate_3d ( eciV(1:3,1), delta * ( m-1 ), orb_norm, eciV(1:3,m) )
      call rotate_3d ( eciV(4:6,1), dVI * ( m-1 ), orb_norm, eciV(4:6,m) )
      eciV(4:6,m) = eciV(4:6,m) * ( ( m-1) * dVIN + 1.0 )
write ( 42, '(i2, 5(3f9.2:","))' ) m, eciV(4:6,m), dVV*(m-1), norm2(eciV(4:6,m))
    end do
    ! Set the length of the SC ECR position vector to the SC Geocentric Altitude.
    ! The ECR length gotten using PGS_CSC_GeoToECR's computation of SC position
    ! from geodetic coordinates is different from the SC Geocentric Altitude by
    ! +/- 300 meters, with a period of about 120 MAFs.  We choose to believe the
    ! SC Geocentric Altitude.
    do m = 1, noMIFs
      eciV(1:3,m) = ( eciV(1:3,m) / norm2(eciV(1:3,m)) ) * emls_oa%sc_geocalt(m)
    end do

    emls_oa%sc_ECI = eciV(1:3,:)
    emls_oa%sc_VelECI = eciV(4:6,:)
  ! Now convert rotated positions and velocities to ECR.
    stat = pgs_csc_ecitoecr ( noMIFs, asciiUTC, offsets, eciV, ecrV )

    emls_oa%sc_ECR = ecrV(1:3,:)
    emls_oa%sc_VelECR = ecrV(4:6,:)

!   ! Now convert rotated positions and velocities back to ECI.
!     stat = pgs_csc_ecrtoeci ( noMIFs, asciiUTC, offsets, ecrV, eciV )
!     emls_oa%sc_ECI = eciV(1:3,:)
!     emls_oa%sc_VelECI = eciV(4:6,:)

  ! Spacecraft geocentric geolocation
    emls_oa%sc_lon = atan2 ( ecrV(2,:), ecrV(1,:) ) * rad2deg
    emls_oa%sc_geoclat = atan2 ( ecrV(3,:), norm2(ecrV(1:2,:),1) ) * rad2deg

  ! Get SC geodetic latitude and altitude from ECR.  We already have longitude.

    do m = 1, noMIFs
      geod = xyz_to_geod ( ecrV(1:3,m) )
      emls_oa%sc_geodlat(m) = geod(1)*rad2deg
      emls_oa%sc_geodalt(m) = geod(3)
    end do

  ! Get GHz ECR, geoc_alt, geoc_lat:

    do m = 1, noMIFs
      stat = pgs_csc_geotoecr (limb_oa%tngt_long(m)*Deg2Rad, &
           limb_oa%tngt_geod_lat(m)*Deg2Rad, limb_oa%tngt_geod_alt(m)*1000.0d0, &
           earthEllipseTag, posecr) !emls_oa%ECR(:,m))
      emls_oa%ECR(:,m) = posecr
      emls_oa%geoc_alt(m) = SQRT (emls_oa%ECR(1,m)**2 + emls_oa%ECR(2,m)**2 + &
           emls_oa%ECR(3,m)**2)
      emls_oa%geoc_lat(m) = Rad2Deg * ATAN2 (posecr(3) , &
           SQRT (posecr(1)**2 + posecr(2)**2))
    end do

  !{ Calculate ECRtoFOV.  Let $\alpha$ be the azimuth angle, $\epsilon$ be
  !  the scan angle, and $\xi$ be the scan angle as rotated during the scan.
  !  Then
  !  \begin{equation*}
  !  F_{\text{MAF}} = \left[ \begin{array}{ccccc}
  !    \cos\alpha \sin\epsilon \sin\xi + \sin\alpha \cos\xi &&
  !    \cos\alpha \sin\epsilon \cos\xi - \sin\alpha \sin\xi &&
  !    \cos\alpha \cos\epsilon \\
  !    \sin\alpha \sin\epsilon \sin\xi - \cos\alpha \cos\xi &&
  !    \sin\alpha \sin\epsilon \cos\xi - \cos\alpha \sin\xi &&
  !    \sin\alpha \cos\epsilon \\
  !    -\cos\epsilon \sin\xi && -\cos\epsilon \cos\xi && \sin\epsilon \\
  !    \end{array} \right]
  !  \end{equation*}

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

      fmaf(:,1) = (/ cosalf*sineps*sinxi+sinalf*cosxi, sinalf*sineps*sinxi-cosalf*cosxi, -coseps*sinxi /)
      fmaf(:,2) = (/ cosalf*sineps*cosxi-sinalf*sinxi, sinalf*sineps*cosxi+cosalf*sinxi, -coseps*cosxi /)
      fmaf(:,3) = (/ cosalf*coseps,                    sinalf*coseps,                     sineps /)

      ipf2eci = TRANSPOSE (MATMUL (limb_oa%trans_inst2eci, fmaf))
      fmaf = TRANSPOSE (fmaf)
      do i = 1, 3
         eciV(1:3,m) = ipf2eci(i,:)
         eciV(4:6,m) = 0.0d0
         stat = pgs_csc_ecitoecr (1, asciiUTC, offsets(m), eciV(:,m), ecrV(:,m))
         ecrtofov(:,i,m) = ecrV(1:3,m)
      end do
      ecrtofov(:,:,m) = TRANSPOSE (ecrtofov(:,:,m))    ! Put in "right" order
      emls_oa%ECRtoFOV(:,m) = RESHAPE (ecrtofov(:,:,m), (/ 9 /))  ! make 2-d

  !{ Calculate GeodAngle (using REC's algorithm):
  !   Let $L =$ \text{limb_oa\%trans_inst2eci}.  Then
  !   ECI = $L\, F_{\text{MAF}}$.
  !
  !  \begin{equation*}
  !  F_{\text{LOS}} = \left[ \begin{array}{ccc}
  !  0           & 0          & -1 \\
  !  -\sin\alpha & \cos\alpha & 0  \\
  !  \cos\alpha  & \sin\alpha & 0  \\
  !    \end{array} \right]
  !  \end{equation*}

      flosf =  TRANSPOSE (RESHAPE ((/ &
           (/ 0.0d0,   0.0d0, -1.0d0 /), &
           (/ -sinalf, cosalf, 0.0d0 /), &
           (/  cosalf, sinalf, 0.0d0 /) /), (/ 3, 3 /)))
      call Calc_GeodAngle ( m, fmaf, flosf, emls_oa )

    end do

  ! Convert GHz ECR to ECI:

    ecrV(1:3,:) = emls_oa%ECR
    ecrV(4:6,:) = 0.0

    stat = pgs_csc_ecrtoeci (noMIFs, asciiUTC, offsets, ecrV, eciV)
    emls_oa%ECI = eciV(1:3,:)

  !{ Calculate LosAngle.  Let $\phi$ = longitude and $\theta$ = geodetic
  ! colatitude (90 degrees $-$ geodetic latitude).  Then
  !  \begin{equation*}
  !  F_{\text{MAF}} = \left[ \begin{array}{ccccc}
  !    -\sin\phi && -\cos\theta \cos\phi && \sin\theta \cos\phi \\
  !    \cos\phi  && -\cos\theta \sin\phi && \sin\theta \sin\phi \\
  !    0         && \sin\theta           && \cos\theta \\
  !    \end{array} \right]
  !  \end{equation*}
  !
  !  $F_{\text{ECR}_i} = \text{ECR to FOV}_{3,i}$, $i = 1 \dots 3$.
  !
  !  $F_t = F_{\text{ECR}}^T \, F_{\text{MAF}}$.
  !
  !  LosAngle = $\mod ( 90 - \tan^{-1} \frac{F_{t_2}}{F_{t_1}} + 360,\, 360)$.

    do m = 1, noMIFs
      phi = emls_oa%lon(m) * Deg2Rad
      cosphi = cos (phi)
      sinphi = sin (phi)
      theta = (90.0 - emls_oa%geod_lat(m)) * Deg2Rad
      costheta = cos (theta)
      sintheta = sin (theta)

      fmaf(:,1) = (/ -sinphi,           cosphi,          0.0_r8 /)
      fmaf(:,2) = (/ -costheta*cosphi, -costheta*sinphi, sintheta /)
      fmaf(:,3) = (/  sintheta*cosphi,  sintheta*sinphi, costheta /)
      fecr = [emls_oa%ecrtofov(3,m), emls_oa%ecrtofov(6,m), emls_oa%ecrtofov(9,m)]
      ft = matmul (fecr, fmaf)
      emls_oa%LosAngle(m) = mod ( 270.0_r8-atan2(ft(2), ft(1))*Rad2Deg, &
                                & 360.0_r8 )

    end do

  !{ Calculate LosVel.  Let $T$ be the tangent position in ECI, and $\omega$
  !  the Earth's angular rotation velocity. Then the tangent velocity, in ECI,
  !  $V_t = \omega [ -T_y, T_x, 0 ]^T$.
  !  Let $S$ be the spacecraft position in ECI.  Then the line-of-sight vector
  !  in ECI,
  !  $L = T - S$.  Let $V_s$ be the spacecraft velocity in ECI.  Then the
  !  line-of-sight velocity,
  !  $V = ( ( V_t  - V_s ) \cdot L ) / |L|$.

    do m = 1, noMIFs
      tngtVel = omega * (/ -emls_oa%ECI(2,m), emls_oa%ECI(1,m), 0.0_r8 /)
      los_vec = emls_oa%ECI(:,m) - emls_oa%sc_ECI(:,m)
      ! Unlike limb_oa%sat_vel, emls_oa%sc_VelECI is now MIF resolved
      emls_oa%LosVel(m) = DOT_PRODUCT ( tngtVel - emls_oa%sc_VelECI(:,m), &
                                      & los_vec ) / norm2 ( los_vec )
    end do

  !{ Calculate OrbY.  Let $\hat{n} = S \times V_s\, /\, | S \times V_s |$ be
  !  the unit normal to the orbit plane in ECI.  Then
  !  $Y_{\text{ORB}} = -T \cdot \hat{n}$.  Let $\gamma$ be the angle between
  !  the normal to the orbit plane and the tangent-point position vector, and
  !  $H_t = |T|$ be the tangent altitude in geocentric coordinates.  Then
  !  $Y_{\text{ORB}} = H_t \cos\gamma$.

    do m = 1, noMIFs
    ! Unit normal to the orbit plane, in ECI
      orb_norm = cross ( emls_oa%sc_ECI(:,m), emls_oa%sc_VelECI(:,m), norm=.true. )
      emls_oa%OrbY(m) = - DOT_PRODUCT ( emls_oa%ECI(:,m), orb_norm )
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

  subroutine Frame_TAI ( limb_oa, SecTAI93, AsciiUTC, n_days, ForgedSecs )

  ! Calculate TAI time at start of frame:

    use Dates_Module, only: SecondsBetween2UTCs
    use MLSKinds, only: R8
    use Rad_File_Contents, only: limb_oa_t

    type(limb_oa_t), intent(in) :: limb_oa
    real(r8), intent(out) :: SecTAI93
    character(len=25), intent(out), optional :: AsciiUTC
    integer, intent(in), optional :: n_days
    real(r8), intent(out), optional :: ForgedSecs

    integer, external :: pgs_td_utctotai

    character(len=25) :: MyAsciiUTC, ForgedAsciiUTC

    integer :: Hrs
    integer :: MIF1_ms
    integer :: Mins
    integer :: Ms
    integer :: N_Years
    real(r8) :: Secs
    integer :: Stat
    integer :: Year
    integer :: YrDoy

    mif1_ms = (limb_oa%ref_mmif - 1) * 2048   ! millisecs of MIF 1
    yrdoy = limb_oa%ref_time(1)  ! - n_days   ! 1000*year plus day of year
    ms = limb_oa%ref_time(2)     ! millisecs of day for ref_mmif
    ms = ms - mif1_ms            ! millisecs at start of MAF
    if (ms < 0) then             ! in previous day!
      ms = ms + 86400000
      yrdoy = yrdoy - 1
    end if
    hrs = mod (ms/3600000, 24)
    mins = mod (ms/60000, 60)
    secs = mod (real(ms), 60000.0) /1000.0
    year = yrdoy / 1000
    if (year > 99) then      ! take care of years 2000 or greater
       year = mod(year, 100) + 2000
    else                     ! years in the 1900's
       year = year + 1900
    end if

    write (myAsciiUTC, fmt= &
     '(i4, "-", i3.3, "T", i2.2, ":", i2.2, ":", f9.6, "Z", TL10, i2.2)') &
     year, mod(yrdoy, 1000), hrs, mins, secs, int(secs)  ! force leading 0's
    if ( deebug ) print *, 'asciiutc: ', trim(myAsciiUTC)
    stat = pgs_td_utctotai (myAsciiUTC, sectai93)

    if ( present(n_days) ) then
      ! Adjust for possible backdating
      if ( n_days > 999 ) then
        ! Account for our convention that a "year" has 1000 days
        n_years = (n_days/1000)
        write (forgedasciiutc, fmt= &
         '(i4, "-", i3.3, "T", i2.2, ":", i2.2, ":", f9.6, "Z", TL10, i2.2)') &
         year-n_years, mod(yrdoy, 1000), hrs, mins, secs, int(secs)  ! force leading 0's
        forgedsecs = secondsBetween2UTCs ( forgedasciiutc, myAsciiUTC )
      else if ( n_days > 0 ) then
        forgedsecs = n_days * 24 * 3600
      else
        forgedsecs = 0.d0
      end if
    end if

    if ( present(asciiUTC) ) asciiUTC = myAsciiUTC

  end subroutine Frame_TAI

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
! Revision 1.5  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
