! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

MODULE TkL1B

  USE Geometry, ONLY: Omega => W
  USE MLSCommon, ONLY: R8, DEFAULTUNDEFINEDVALUE
  USE MLSL1Common
  USE MLSMessageModule, ONLY: MLSMESSAGE, MLSMSG_Error, MLSMSG_Warning, &
       MLSMSG_Allocate
  USE OUTPUT_M, ONLY: BLANKS, OUTPUT
  USE OutputL1B_DataTypes, ONLY: L1BOAsc_T, L1BOATP_T, L1BOAINDEX_T, LENCOORD, &
       LENG, LENT
  USE OutputL1B, ONLY: OUTPUTL1B_THZ, OUTPUTL1B_SC, OUTPUTL1B_INDEX, &
       OUTPUTL1B_GHZ
  USE Scan, ONLY : Scan_guess, Scan_start
  USE SDPToolkit
  USE MLSL1Utils, ONLY: Finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L1BOA_MAF, Flag_Bright_Objects, LOG_ARR1_PTR_T, GHz_GeodAlt, &
       GHz_GeodLat, GHz_GeodAngle

  TYPE LOG_ARR1_PTR_T
     LOGICAL, DIMENSION(:), POINTER :: ptr
  END TYPE LOG_ARR1_PTR_T

  REAL :: GHz_GeodAlt(LENG), GHz_GeodLat(LENG), GHz_GeodAngle(LENG)

  LOGICAL, PARAMETER :: ORBINCLINE_IS_CONSTANT = .FALSE.
  REAL, PARAMETER ::    UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE ! -999.99
  REAL, PARAMETER ::    HUGE_F = HUGE (1.0)

  !------------------- RCS Ident Info -----------------------
  CHARACTER(len=*), PARAMETER :: IdParm = &
    & "$Id$"
  CHARACTER(len=LEN(idParm)) :: Id = idParm
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  ! This module contains subroutines for producing the L1BOA records on
  ! a MAF by MAF basis.

CONTAINS

!=============================================================================
  SUBROUTINE Init_L1BOAsc (sc)
!=============================================================================

    TYPE (L1BOAsc_T) :: sc

    sc%scECI = HUGE_F
    sc%scECR = HUGE_F
    sc%scGeocAlt = HUGE_F
    sc%scGeodAlt = HUGE_F
    sc%scGeocLat = HUGE_F
    sc%scGeodLat = HUGE_F
    sc%scLon = HUGE_F
    sc%scGeodAngle = HUGE_F
    sc%scOrbIncl = HUGE_F
    sc%scVelECI = HUGE_F
    sc%scVelECR = HUGE_F
    sc%ypr = HUGE_F
    sc%yprRate = HUGE_F

  END  SUBROUTINE Init_L1BOAsc

!=============================================================================
  SUBROUTINE Init_L1BOAtp (tp, MIFbad)
!=============================================================================

    TYPE (L1BOAtp_T) :: tp
    LOGICAL, OPTIONAL, INTENT(IN) :: MIFbad(:)

    LOGICAL :: badindx (SIZE(tp%tpLon))  ! Indexes for "bad" MIFs

    IF (PRESENT (MIFbad)) THEN
       badindx = MIFbad
    ELSE
       badindx = .TRUE.
    ENDIF

    WHERE (badindx)
       tp%tpGeodAlt = UNDEFINED_VALUE
       tp%tpGeocAlt = UNDEFINED_VALUE
       tp%tpOrbY = UNDEFINED_VALUE
       tp%tpGeocLat = UNDEFINED_VALUE
       tp%tpGeocAltRate = UNDEFINED_VALUE
       tp%tpGeodLat = UNDEFINED_VALUE
       tp%tpGeodAltRate = UNDEFINED_VALUE
       tp%tpLon = UNDEFINED_VALUE
       tp%tpGeodAngle = UNDEFINED_VALUE
       tp%tpSolarTime = UNDEFINED_VALUE
       tp%tpSolarZenith = UNDEFINED_VALUE
       tp%tpLosAngle = UNDEFINED_VALUE
       tp%tpLosVel = UNDEFINED_VALUE
    ENDWHERE

  END SUBROUTINE Init_L1BOAtp

  !------------------------------------------ TkL1B_sc ---------
  SUBROUTINE TkL1B_sc (numValues, offsets, asciiUTC, mafTAI, sc, ecrtosc, &
       oastat)
    ! This subroutine contains prototype code for creating the desired s/c
    ! record from the EPHEMATTIT output.

    ! Arguments
    TYPE (L1BOAsc_T) :: sc
    CHARACTER (LEN=27), INTENT(IN) :: asciiUTC
    INTEGER, INTENT(IN) :: numValues
    REAL(r8), INTENT(IN) :: offsets(numValues), mafTAI
    REAL, INTENT(OUT) :: ecrtosc(3,3,numValues)
    INTEGER, INTENT(OUT) :: oastat

    ! Functions
    INTEGER :: Pgs_csc_eciToECR, Pgs_csc_ecrToGEO, Pgs_eph_ephemAttit

    ! Variables
    CHARACTER (LEN=32) :: mnemonic
    CHARACTER (LEN=480) :: msg, msr
    INTEGER :: i, j, returnStatus
    INTEGER :: qualityFlags(2, numValues)
    REAL(r8) :: ecrVec(3), eulerangles(3,numValues)
    REAL(r8) :: radC(numValues), radD(numValues), radL(numValues)
    REAL(r8), TARGET :: attitQuat(4,numValues)
    REAL(r8) :: eciV(6,numValues), ecrV(6,numValues)
    REAL(r8) :: sctoeci(6,3*numValues), sctoecr(6,3*numValues) 
    REAL(r8) , POINTER :: w(:), x(:), y(:), z(:)
    REAL(r8), PARAMETER :: SCtoGHz(3,3) = RESHAPE ((/ &
         -0.0000086, -0.4260405, 0.9047041, &
          1.0000000,  0.0000000, 0.0000095, &
         -0.0000041,  0.9047041, 0.4260405 /), (/ 3, 3 /))

    ! Executable code

    ! Read oa data

    returnStatus = Pgs_eph_ephemAttit (spacecraftId, numValues, asciiUTC,  &
         offsets, pgs_true, pgs_true, qualityFlags, sc%scECI, sc%scVelECI, &
         eulerangles, sc%yprRate, attitQuat)

! save YPR:

    sc%ypr(1,:) = eulerangles(1,:)
    sc%ypr(2,:) = eulerangles(3,:)
    sc%ypr(3,:) = eulerangles(2,:)

! Trap any bad attitQuat data

    IF (returnStatus == PGS_S_SUCCESS) THEN
       outer: DO j = 1, numValues
          DO i = 1, 4
             IF (.NOT. Finite (REAL(attitQuat(i,j)))) THEN
                returnStatus = PGSEPH_W_BAD_EPHEM_VALUE
                EXIT outer
             ENDIF
          ENDDO
       ENDDO outer
    ENDIF

    oastat = returnStatus

    IF (returnStatus /= PGS_S_SUCCESS) THEN
       CALL Pgs_smf_getMsg (returnStatus, mnemonic, msg)
       msr = 'Routine ephemAttit, ' // mnemonic // ':  ' // msg
       CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
       ! Initial out values:

       CALL Init_L1BOAsc (sc)
       ecrtosc = HUGE_F

       RETURN
    ENDIF

    ! attitQuat is spacecraft to eci rotation quaternion
    w => attitquat(1,:)
    x => attitquat(2,:)
    y => attitquat(3,:)
    z => attitquat(4,:)

    sctoeci = 0.0

    ! Convert quaternion to matrix (optimized for readability)

    sctoeci(1, 1::3) = w**2 + x**2 - y**2 -z**2
    sctoeci(2, 2::3) = w**2 - x**2 + y**2 -z**2
    sctoeci(3, 3::3) = w**2 - x**2 - y**2 +z**2
    sctoeci(1, 2::3) = 2*x*y - 2*w*z
    sctoeci(2, 1::3) = 2*x*y + 2*w*z
    sctoeci(1, 3::3) = 2*x*z + 2*w*y
    sctoeci(3, 1::3) = 2*x*z - 2*w*y
    sctoeci(2, 3::3) = 2*y*z - 2*w*x
    sctoeci(3, 2::3) = 2*y*z + 2*w*x

    ! rotate the columns of the sctoeci matrix to form ecrtosc matrix columns

    returnStatus = PGS_CSC_ECItoECR (3*numValues, asciiUTC, &
         PACK(SPREAD(offsets,1,3), .TRUE.), sctoeci, sctoecr)

    ! we really prefer the other index to be exposed for matrix multiplication
    ! now we have [ECR x (SC*numValues)] which requires a loop for rotation
    ! SC-->GHz reference (R3,B8)
    
    DO i = 1, numValues    
       ecrtosc(:,:,i) = TRANSPOSE (sctoecr(1:3,3*i-2:3*i))
    ENDDO 

    ! Convert scECI to scECR
    eciV(1:3,:) = sc%scECI
    eciV(4:6,:) = sc%scVelECI    ! was 0.0
    returnStatus = Pgs_csc_eciToECR (numValues, asciiUTC, offsets, eciV, ecrV)
    sc%scECR = ecrV(1:3,:)
    sc%scVelECR = ecrV(4:6,:)

    ! Calculate geocentric/geodetic altitude, latitude & longitude from scECR
    DO i = 1, numValues
       sc%scGeocAlt(i) = SQRT (ecrV(1,i)**2 + ecrV(2,i)**2 + ecrV(3,i)**2)
       radC(i) = ATAN (ecrV(3,i) / SQRT (ecrV(1,i)**2 + ecrV(2,i)**2))
       sc%scGeocLat(i) = Rad2Deg * radC(i)
       ecrVec = ecrV(1:3,i)
       returnStatus = Pgs_csc_ecrToGEO (ecrVec, earthModel, radL(i), &
            radD(i), sc%scGeodAlt(i))
       sc%scGeodLat(i) = Rad2Deg * radD(i)
       sc%scLon(i) = Rad2Deg * radL(i)
       IF (ORBINCLINE_IS_CONSTANT) THEN
          sc%scOrbIncl(i) = orbInclineCrossProd (sc%scECI(:,i), &
               sc%scVelECI(:,i))
       ELSE
          sc%scOrbIncl(i) = orbInclineCalculated (sc%scECR(:,i), &
               sc%scVelECR(:,i), sc%scGeocLat(i), sc%scLon(i))
       ENDIF
       sc%MIF_TAI(i) = mafTAI + offsets(i)
    ENDDO

  END SUBROUTINE TkL1B_sc

  !------------------------------------------------------TkL1B_tp ----
  SUBROUTINE TkL1B_tp (asciiTAI, asciiUTC, lenG, numValues, offsets, posECR, &
    posECI, velECI, scAngle, encoderAngle, tp, ecrtosc, GroundToFlightMounts, &
    ScToGroundMounts, gtindx)
    ! This subroutine fills the tangent point record.

    USE FOV, ONLY: CalcMountsToFOV
    USE SDPToolkit, ONLY: PGS_S_SUCCESS, PGSCSC_W_LOOK_AWAY, PGSCSC_W_HIT_EARTH 

    ! Arguments
    TYPE (L1BOAtp_T) :: tp
    CHARACTER (LEN=27), INTENT(IN) :: asciiUTC
    INTEGER, INTENT(IN) :: lenG, numValues, gtindx
    REAL, INTENT(IN) :: scAngle(numValues), encoderAngle(numValues)
    REAL(r8), INTENT(IN) :: asciiTAI
    REAL(r8), INTENT(IN) :: offsets(numValues)
    REAL(r8), INTENT(IN) :: posECR(3,lenG), posECI(3,lenG), velECI(3,lenG)
    REAL, INTENT(IN) :: ecrtosc(3,3,numValues)

    REAL(r8), INTENT(IN) :: GroundToFlightMounts(3,3)
    REAL(r8), INTENT(IN) :: ScToGroundMounts(3,3)

    ! Functions
    INTEGER :: Pgs_cbp_sat_cb_vector, Pgs_cbp_solarTimeCoords
    INTEGER :: Pgs_csc_scToOrb, Pgs_csc_scToECI, Pgs_csc_eciToECR
    INTEGER :: Pgs_csc_grazingRay, Pgs_csc_ecrToECI, Pgs_csc_eciToOrb
    INTEGER :: Pgs_td_taiToUTC

    ! Variables
    CHARACTER (LEN=27) :: time
    INTEGER :: flag, flagQ, i, returnStatus
    LOGICAL :: MIFbad(lenG)
    REAL(r8) :: declination, deltaAlt, deltaLat, deltaLon
    REAL(r8) :: greenwich, localApparent, rightAscension, tai
    REAL(r8) :: dot(lenG), latD(lenG), localMean(lenG), lon(lenG), los(lenG)
    REAL(r8) :: ecr_sign(lenG), slantRange(lenG)
    REAL(r8) :: angleRad(numValues)
    REAL(r8) :: eci(3,lenG), ecr(3,lenG), hECR(3,lenG), nts(3,lenG)
    REAL(r8) :: posSurf(3,lenG), sc_frame_vector(3,lenG), sc_sun(3,lenG)
    REAL(r8) :: sc_tp(3,lenG), tp_sun(3,lenG), tpOrb(3,lenG), unitAlt(3,lenG)
    REAL(r8) :: unitLat(3,lenG), unitLon(3,lenG), vECR(3,lenG)
    REAL(r8) :: fov_sc(3,numValues), fov_orb(3,numValues)
    REAL(r8) :: eciV(6,lenG), ecrV(6,lenG)
    REAL(r8) :: tngtVel(3), los_vec(3)
    REAL(r8) :: MountsToFOV(3,3), ECRtoFOV(3,3)
    REAL(r8) :: GroundMountsToFOV(3,3), ScToFOV(3,3)
    CHARACTER (LEN=32) :: mnemonic
    CHARACTER (LEN=480) :: msg, msr

    ! Executable code

    deltaLat = 0.1
    deltaLon = 0.1
    deltaAlt = 1000.0
    tp%encoderAngle = encoderAngle

    ! Determine ECR to FOV

    DO i = 1, numValues

       CALL CalcMountsToFOV (scAngle(i), gtindx, MountsToFOV)

       GroundMountsToFOV = MATMUL (MountsToFOV, GroundToFlightMounts)
       ScToFOV = MATMUL (GroundMountsToFOV, ScToGroundMounts)
       ECRtoFOV = MATMUL (ScToFOV, ecrtosc(:,:,i))

       IF (i <= lenG) tp%tpECRtoFOV(:,i) = RESHAPE (ECRtoFOV, (/ 9 /))
       
       fov_sc(1,i) = ScToFOV(3,1)
       fov_sc(2,i) = ScToFOV(3,2)
       fov_sc(3,i) = ScToFOV(3,3)

    ! Put sc angle

       tp%scAngle(i) = 90.0 - &
            ACOS (MAX (MIN(fov_sc(3,i), 1.0d0), -1.0d0)) * Rad2Deg
    ENDDO

    ! Convert s/c vector to Orb vector/angle/degrees

    returnStatus = Pgs_csc_scToOrb (spacecraftId, numValues, asciiUTC, &
      offsets, fov_sc, fov_orb)
    IF (returnStatus /= PGS_S_SUCCESS) THEN
       CALL Pgs_smf_getMsg (returnStatus, mnemonic, msg)
       msr = 'Routine scToOrb, ' // mnemonic // ':  ' // msg
       CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
 
       ! Initial out values:

       CALL Init_L1BOAtp (tp)
       RETURN
    ENDIF

    tp%scanAngle = 90.0 - ACOS (MAX (MIN (fov_orb(3,:), 1.0d0), -1.0d0)) * &
         Rad2Deg
    tp%scanRate(1) = 0.0
    tp%scanRate(2:) = ABS(tp%scanAngle(2:) - tp%scanAngle(1:)) / offsets(2)
    tp%azimAngle = ATAN2 (fov_orb(2,:), fov_orb(1,:)) * Rad2Deg

    ! Convert s/c vector to ECR
    returnStatus = Pgs_csc_scToECI (spacecraftId, lenG, asciiUTC, &
      offsets(1:lenG), fov_sc(:,1:lenG), eci)
    eciV(1:3,:) = eci
    eciV(4:6,:) = 0.0
    returnStatus = Pgs_csc_eciToECR (lenG, asciiUTC, offsets(1:lenG), eciV, &
         ecrV)
    ecr = ecrV(1:3,:)

    ! For each scanning MIF

    MIFbad = .FALSE.
    DO i = 1, lenG

      ! Calculate tangent point (geodetic & ECR)
      returnStatus = Pgs_csc_grazingRay (earthModel, posECR(:,i), ecr(:,i), &
           latD(i), lon(i), tp%tpGeodAlt(i), &
           slantRange(i), tp%tpECR(:,i), posSurf(:,i))

      IF (returnStatus /= PGS_S_SUCCESS .AND. &
           returnStatus /= PGSCSC_W_HIT_EARTH) THEN  ! success or "hit" earth
         latD(i) = 0.0 ! HUGE_F
         lon(i) = 0.0 ! HUGE_F
         tp%tpGeodAlt(i) = 0.0 ! HUGE_F
         slantRange(i) = 0.0 ! HUGE_F
         tp%tpECR(:,i) = 0.0 ! HUGE_F
         posSurf(:,i) = 0.0 ! HUGE_F
         MIFbad(i) = .TRUE.
      ENDIF

      ! Create ECR unit vector quantities -- lat=1, lon=2, alt=3
      flagQ = 2
      CALL Tp_unit (flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLon, &
           unitLon(:,i), flag)
      flagQ = 1
      CALL Tp_unit (flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLat, &
           unitLat(:,i), flag)
      flagQ = 3
      CALL Tp_unit (flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaAlt, &
           unitAlt(:,i), flag)

      ! Get local mean solar time from Toolkit
      tai = asciiTAI + (i-1)*offsets(2)
      returnStatus = Pgs_td_taiToUTC (tai, time)
      returnStatus = Pgs_cbp_solarTimeCoords (time, lon(i), greenwich, &
        localMean(i), localApparent, rightAscension, declination)
    ENDDO

    tp%tpGeodLat = Rad2Deg * latD
    tp%tpLon = Rad2Deg * lon
    tp%tpSolarTime = localMean / 3600.0

    ! Calculate solarZenith

    returnStatus = Pgs_cbp_sat_cb_vector (spacecraftId, lenG, asciiUTC, &
         offsets(1:lenG), PGSd_SUN, sc_frame_vector)
    returnStatus = Pgs_csc_scToECI (spacecraftId, lenG, asciiUTC, &
         offsets(1:lenG), sc_frame_vector, eci)
    eciV(1:3,:) = eci
    eciV(4:6,:) = 0.0
    returnStatus = Pgs_csc_eciToECR (lenG, asciiUTC, offsets(1:lenG), eciV, &
         ecrV)

    sc_sun = ecrV(1:3,:)
    DO i = 1, lenG
       sc_tp(:,i) = ecr(:,i) * slantRange(i)
    ENDDO
    tp_sun = sc_sun - sc_tp
    DO i = 1, lenG
       nts(:,i) = tp_sun(:,i) / SQRT (tp_sun(1,i)**2 + tp_sun(2,i)**2 + &
            tp_sun(3,i)**2)
    ENDDO
    dot = nts(1,:)*unitAlt(1,:) + nts(2,:)*unitAlt(2,:) + &
         nts(3,:)*unitAlt(3,:)
    tp%tpSolarZenith = ACOS(dot) * Rad2Deg

    ! Calculate losAngle
    DO i = 1, lenG
       vECR(:,i) = ecr(:,i) - (ecr(1,i)*unitAlt(1,i) + &
            ecr(2,i)*unitAlt(2,i) + ecr(3,i)*unitAlt(3,i))
       hECR(:,i) = vECR(:,i) / SQRT(vECR(1,i)**2 + vECR(2,i)**2 + &
            vECR(3,i)**2)
    ENDDO
    los = ACOS(hECR(1,:)*unitLat(1,:) + hECR(2,:)*unitLat(2,:) + &
         hECR(3,:)*unitLat(3,:))
    ecr_sign = ecr(1,:)*unitLon(1,:) + ecr(2,:)*unitLon(2,:) + &
      ecr(3,:)*unitLon(3,:)
    DO i = 1, lenG
       IF (ecr_sign(i) < 0 ) los(i) = 2*PI - los(i)
    ENDDO
    tp%tpLosAngle = los * Rad2Deg

    ! Convert tpECR to tpECI
    ecrV(1:3,:) = tp%tpECR
    ecrV(4:6,:) = 0.0
    returnStatus = Pgs_csc_ecrToECI (lenG, asciiUTC, offsets(1:lenG), ecrV, &
      eciV)
    tp%tpECI = eciV(1:3,:)

    ! Calculate losVel
    DO i = 1, lenG
      tngtVel = omega * (/ -tp%tpECI(2,i), tp%tpECI(1,i), 0.0_r8 /)
      los_vec = tp%tpECI(:,i) - posECI(:,i)
      los_vec = los_vec / SQRT (SUM (los_vec**2))
      tp%tpLosVel(i) = DOT_PRODUCT (tngtVel, los_vec) - &
           DOT_PRODUCT (velECI(:,i), los_vec)
    ENDDO

    ! Convert tpECI to tpOrb
    returnStatus = Pgs_csc_eciToOrb (spacecraftId, lenG, asciiUTC, &
         offsets(1:lenG), tp%tpECI, tpOrb)
    tp%tpOrbY = tpOrb(2,:)

    ! Calculate tp geocentric coordinates
    DO i = 1, lenG
      tp%tpGeocAlt(i) = SQRT( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2 + &
           tp%tpECR(3,i)**2 )
      tp%tpGeocLat(i) = Rad2Deg * ATAN (tp%tpECR(3,i) &
           / SQRT( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2))
    ENDDO

    ! Calculate dummy values

    DO i = 2, lenG
      tp%tpGeocAltRate(i-1) = (tp%tpGeocAlt(i) - tp%tpGeocAlt(i-1)) &
           / offsets(2)
      tp%tpGeodAltRate(i-1) = (tp%tpGeodAlt(i) - tp%tpGeodAlt(i-1)) &
           / offsets(2)
    ENDDO
    tp%tpGeocAltRate(lenG) = tp%tpGeocAltRate(lenG-1)
    tp%tpGeodAltRate(lenG) = tp%tpGeodAltRate(lenG-1)

    IF (ANY (MIFbad)) THEN
       CALL Init_L1BOAtp (tp, MIFbad)
       print *, 'some bad MIFs...'
    ENDIF

  END SUBROUTINE TkL1B_tp

  !------------------------------------------------- Tp_unit ---------
  SUBROUTINE Tp_unit (flagQ, lon, lat, alt, delta, unitQ, flag)
    ! This subroutine creates unit ECR vector quantities used by TkL1B_tp to
    ! calculate solarZenith and losAngle.

    ! Arguments
    INTEGER, INTENT(IN) :: flagQ
    REAL(r8), INTENT(IN) :: alt, delta, lat, lon
    INTEGER, INTENT(OUT) :: flag
    REAL(r8), INTENT(OUT) :: unitQ(3)

    ! Functions
    INTEGER :: Pgs_csc_geoToECR

    ! Variables
    CHARACTER (LEN=32) :: mnemonic
    CHARACTER (LEN=480) :: msg, msr
    INTEGER :: returnStatus
    REAL(r8) :: del, geoMinus, geoPlus
    REAL(r8) :: ecrMinus(3), ecrPlus(3), vec(3), vecSqrtSum

    ! Exectuable code

    flag = 0

    ! Create unit vector
    IF (flagQ == 2) THEN
      ! longitude
      del = delta * Deg2Rad
      geoPlus  = (lon + del/2)
      geoMinus = (lon - del/2)
      IF (geoPlus  >  PI) geoPlus  = geoPlus  - 2*PI
      IF (geoMinus < -PI) geoMinus = geoMinus + 2*PI
      returnStatus = Pgs_csc_geoToECR (geoPlus, lat, alt, earthModel, &
           ecrPlus)
      returnStatus = Pgs_csc_geoToECR (geoMinus, lat, alt, earthModel, &
           ecrMinus)
    ELSE IF (flagQ == 1) THEN
      ! latitude
      del = delta * Deg2Rad
      geoPlus  = (lat + del/2)
      geoMinus = (lat - del/2)
      returnStatus = Pgs_csc_geoToECR (lon, geoPlus, alt, earthModel, &
           ecrPlus)
      returnStatus = Pgs_csc_geoToECR (lon, geoMinus, alt, earthModel, &
           ecrMinus)
    ELSE
      geoPlus  = (alt + delta/2)
      geoMinus = (alt - delta/2)
      returnStatus = Pgs_csc_geoToECR (lon, lat, geoPlus, earthModel, &
           ecrPlus)
      returnStatus = Pgs_csc_geoToECR (lon, lat, geoMinus, earthModel, &
           ecrMinus)
    ENDIF

    IF (returnStatus /= PGS_S_SUCCESS) THEN
       CALL Pgs_smf_getMsg (returnStatus, mnemonic, msg)
       msr = mnemonic // ':  ' // msg
       CALL MLSMessage (MLSMSG_Warning, ModuleName, msr)
       flag = -1
    ENDIF

    vec = ecrPlus - ecrMinus
    vecSqrtSum = SQRT (vec(1)**2 + vec(2)**2 + vec(3)**2)
    IF (vecSqrtSum > 0.0) THEN
       unitQ = vec / vecSqrtSum
    ELSE
       unitQ = 0.0
    ENDIF

  END SUBROUTINE Tp_unit

  !-------------------------------------------------- L1BOA_MAF ----------------
  SUBROUTINE L1BOA_MAF (altG, altT, ascTAI, counterMAF, dscTAI, L1FileHandle, &
       MAFinfo, noMAF, MIFsPerMAF, numOrb, scAngleG, scAngleT, encAngleG, &
       encAngleT)

    ! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
    ! HDF output file.

    USE MLSL1Config, ONLY: L1Config
    USE FOV, ONLY: GroundToFlightMountsGHz, GroundToFlightMountsTHz, &
         ScToGroundMountsGHz,ScToGroundMountsTHz

    ! Arguments
    TYPE (MAFinfo_T) :: MAFinfo
    INTEGER, INTENT(IN) :: L1FileHandle, counterMAF, noMAF, numOrb, MIFsPerMAF
    REAL, INTENT(INOUT) :: scAngleG(:), scAngleT(:)
    REAL, INTENT(IN) :: encAngleG(:), encAngleT(:)  ! encoder Angles
    REAL(r8), INTENT(IN) :: altG, altT
    REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:)

    ! Functions
    INTEGER :: Pgs_td_taiToUTC

    ! Variables
    TYPE( L1BOAindex_T ) :: index
    TYPE( L1BOAsc_T ) :: sc
    TYPE( L1BOAtp_T ) :: tp
    CHARACTER (LEN=27) :: mafTime
    CHARACTER (LEN=480) :: msr
    INTEGER :: error, i, nV, returnStatus, oastat, gtindx
    REAL(r8) :: mafTAI
    REAL(r8) :: initRay(3), q(3,MAFinfo%MIFsPerMAF)
    REAL(r8) :: offsets(MAFinfo%MIFsPerMAF)
    REAL :: angle_del
    REAL :: ecrtosc(3,3,MAFinfo%MIFsPerMAF)

    ! Calculate time (secs) offsets array

    mafTAI = MAFinfo%startTAI
    nV = MAFinfo%MIFsPerMAF
    DO i = 1, nV 
      offsets(i) = (i-1)*MAFinfo%MIF_dur
    ENDDO

    ! Write "index" information

    returnStatus = Pgs_td_taiToUTC (mafTAI, mafTime)
    index%MAFStartTimeUTC = mafTime
    index%MAFStartTimeTAI = mafTAI
    index%noMIFs = MIFsPerMAF
    index%counterMAF = counterMAF

    CALL OutputL1B_index (noMAF, L1FileHandle, index)

    ! Allocate the MIF variables in the output structures

    IF (.NOT. ASSOCIATED(sc%MIF_TAI)) THEN
       ALLOCATE (sc%MIF_TAI(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  MIF_TAI quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scGeocAlt)) THEN
       ALLOCATE (sc%scGeocAlt(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c  GeocAlt quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scGeocLat)) THEN
       ALLOCATE (sc%scGeocLat(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c  GeocLat quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scGeodAlt)) THEN
       ALLOCATE (sc%scGeodAlt(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c  GeodAlt quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scGeodLat)) THEN
       ALLOCATE (sc%scGeodLat(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c GeodLat  quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scLon)) THEN
       ALLOCATE (sc%scLon(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c Lon quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scGeodAngle)) THEN
       ALLOCATE (sc%scGeodAngle(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c GeodAngle quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scOrbIncl)) THEN
       ALLOCATE (sc%scOrbIncl(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c OrbIncl quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%encoderAngle)) THEN
       ALLOCATE (tp%encoderAngle(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  tp encoder angle quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%scAngle)) THEN
       ALLOCATE (tp%scAngle(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  tp angle quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%scanAngle)) THEN
       ALLOCATE (tp%scanAngle(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  tp scan angle quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%azimAngle)) THEN
       ALLOCATE (tp%azimAngle(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  tp azim angle quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%scanRate)) THEN
       ALLOCATE (tp%scanRate(nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  tp scan rate quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scECI)) THEN
       ALLOCATE (sc%scECI(lenCoord,nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c ECI quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scECR)) THEN
       ALLOCATE (sc%scECR(lenCoord,nV), STAT=error) 
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c ECR quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scVelECI)) THEN
       ALLOCATE (sc%scVelECI(lenCoord,nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c VelECI quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%scVelECR)) THEN
       ALLOCATE (sc%scVelECR(lenCoord,nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c VelECR quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%ypr)) THEN
       ALLOCATE (sc%ypr(lenCoord,nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c ypr quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(sc%yprRate)) THEN
       ALLOCATE (sc%yprRate(lenCoord,nV), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  s/c ypr rate quantities.'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    ! Get oa data

    CALL TkL1B_sc (nV, offsets, mafTime, mafTAI, sc, ecrtosc, oastat)

    IF (oastat == PGS_S_SUCCESS) THEN

    ! Get s/c master coordinate

       CALL Mc_aux (mafTime, offsets, sc%scECR, q)

       CALL TkL1B_mc (ascTAI, dscTAI, sc%scECR, nV, numOrb, &
            q, mafTAI, offsets, sc%scGeodAngle, sc%scOrbIncl)

    ENDIF

    ! Write s/c information

    CALL OutputL1B_sc (noMAF, L1FileHandle, sc)

    ! Allocate the output structure

    IF (.NOT. ASSOCIATED(tp%tpECI)) THEN
       ALLOCATE (tp%tpECI(lenCoord,lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: ECI'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpECR)) THEN
       ALLOCATE (tp%tpECR(lenCoord,lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: ECR'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpECRtoFOV)) THEN
       ALLOCATE (tp%tpECRtoFOV(lenCoord*lenCoord,lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: ECRtoFOV'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpOrbY)) THEN
       ALLOCATE (tp%tpOrbY(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: OrbY'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeocAlt)) THEN
       ALLOCATE (tp%tpGeocAlt(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeocAlt'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeocLat)) THEN
       ALLOCATE (tp%tpGeocLat(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeocLat'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeocAltRate)) THEN
       ALLOCATE (tp%tpGeocAltRate(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeocAltRate'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeodAlt)) THEN
       ALLOCATE (tp%tpGeodAlt(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeodAlt'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeodLat)) THEN
       ALLOCATE (tp%tpGeodLat(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeodLat'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeodAltRate)) THEN
       ALLOCATE (tp%tpGeodAltRate(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeodAltRate'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpLon)) THEN
       ALLOCATE (tp%tpLon(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: Lon'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpGeodAngle)) THEN
       ALLOCATE (tp%tpGeodAngle(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: GeodAngle'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpSolarTime)) THEN
       ALLOCATE (tp%tpSolarTime(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: SolarTime'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpSolarZenith)) THEN
       ALLOCATE (tp%tpSolarZenith(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: SolarZenith'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpLosAngle)) THEN
       ALLOCATE (tp%tpLosAngle(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: LosAngle'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (.NOT. ASSOCIATED(tp%tpLosVel)) THEN
       ALLOCATE (tp%tpLosVel(lenG), STAT=error)
       IF (error /= 0) THEN
          msr = MLSMSG_Allocate // '  GHz tp quantities: LosVel'
          CALL MLSMessage (MLSMSG_Error, ModuleName, msr)
       ENDIF
    ENDIF

    IF (oastat == PGS_S_SUCCESS) THEN

    ! Calculate initial guess for look vector in ECR

    CALL Scan_guess (mafTime, initRay)

    ! Find angle, tan pt for start of GHZ scan

    CALL Scan_start (altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

    ! Calculate GHZ tan pt record

    gtindx = 1
    CALL TkL1B_tp (mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenG), &
         sc%scECI(:,1:lenG), sc%scVelECI(:,1:lenG), scAngleG, encAngleG, tp, &
         ecrtosc, GroundToFlightMountsGHz, ScToGroundMountsGHz, gtindx)
    IF (L1Config%Globals%SimOA) THEN   ! correct nominal scan angles for sim
       angle_del = tp%tpGeodAlt(1) / 5200.0 * 0.1
       scAngleG = scAngleG + angle_del
       CALL TkL1B_tp (mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenG), &
            sc%scECI(:,1:lenG), sc%scVelECI(:,1:lenG), scAngleG, encAngleG, &
            tp, ecrtosc, GroundToFlightMountsGHz, ScToGroundMountsGHz, gtindx)
       angle_del = tp%tpGeodAlt(1) / 5200.0 * 0.1
       scAngleG = scAngleG + angle_del
       CALL TkL1B_tp (mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenG), &
            sc%scECI(:,1:lenG), sc%scVelECI(:,1:lenG), scAngleG, encAngleG, &
            tp, ecrtosc, GroundToFlightMountsGHz, ScToGroundMountsGHz, gtindx)
    ENDIF

    ! Compute GHz master coordinate

    CALL TkL1B_mc (ascTAI, dscTAI, tp%tpECR, lenG, numOrb, &
         q, mafTAI, offsets(1:lenG), tp%tpGeodAngle, sc%scOrbIncl)

    ! Save info for baseline corrections

    GHz_GeodAlt = tp%tpGeodAlt
    GHz_GeodLat = tp%tpGeodLat
    GHz_GeodAngle = tp%tpGeodAngle

    ELSE

       CALL Init_L1BOAtp (tp)

    ENDIF

    ! Write GHz information

    CALL OutputL1B_GHz (noMAF, L1FileHandle, tp)

    IF (oastat == PGS_S_SUCCESS) THEN

    ! Find angle, tan pt for start of THz scan

    CALL Scan_start (altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

    ! Calculate THz tan pt record

    gtindx = 2
    CALL TkL1B_tp (mafTAI, mafTime, lenT, nV, offsets, sc%scECR(:,1:lenT), &
         sc%scECI(:,1:lenG), sc%scVelECI(:,1:lenG), scAngleT, encAngleT, tp, &
         ecrtosc, GroundToFlightMountsTHz, ScToGroundMountsTHz, gtindx)
    IF (L1Config%Globals%SimOA) THEN   ! correct nominal scan angles for sim
       angle_del = tp%tpGeodAlt(1) / 5200.0 * 0.1
       scAngleT = scAngleT + angle_del
       CALL TkL1B_tp (mafTAI, mafTime, lenT, nV, offsets, sc%scECR(:,1:lenT), &
            sc%scECI(:,1:lenG), sc%scVelECI(:,1:lenG), scAngleT, encAngleT, &
            tp, ecrtosc, GroundToFlightMountsTHz, ScToGroundMountsTHz, gtindx)
       angle_del = tp%tpGeodAlt(1) / 5200.0 * 0.1
       scAngleT = scAngleT + angle_del
       CALL TkL1B_tp (mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenT), &
            sc%scECI(:,1:lenT), sc%scVelECI(:,1:lenT), scAngleT, encAngleT, &
            tp, ecrtosc, GroundToFlightMountsTHz, ScToGroundMountsTHz, gtindx)
    ENDIF

    ! Compute THz master coordinate

    CALL TkL1B_mc (ascTAI, dscTAI, tp%tpECR, lenT, numOrb, &
         q, mafTAI, offsets(1:lenT), tp%tpGeodAngle, sc%scOrbIncl)

    ! Write THZ information

    ELSE

       CALL Init_L1BOAtp (tp)

    ENDIF

    CALL OutputL1B_THz (noMAF, L1FileHandle, tp)

  END SUBROUTINE L1boa_MAF

  !------------------------------------------- Mc_Aux --------
  SUBROUTINE Mc_aux (asciiUTC, offsets, scECR, q)
    ! This subroutine computes q, an auxilliary vector used in the calculation
    ! of the master coordinate.  Q is a vector that points from the center of
    ! the Earth to the ascending node of the orbit in ECI coordinates.  Thus any
    ! point dotted with q can give you the master coordinate.
    ! Arguments
    CHARACTER (LEN=27), INTENT(IN) :: asciiUTC
    REAL(r8), INTENT(in) :: OFFSETS(:)
    REAL(r8), INTENT(IN) :: scECR(:,:)
    REAL(r8), INTENT(OUT) :: q(:,:)

    ! Parameters

    ! Functions
    INTEGER :: Pgs_csc_orbToECI, PGS_CSC_ECItoECR

    ! Variables
    INTEGER :: returnStatus
    INTEGER :: nV
    REAL(r8), DIMENSION(SIZE(offsets)) :: l
    REAL(r8), DIMENSION(SIZE(offsets)) :: DIST1
    REAL(r8), DIMENSION(SIZE(offsets)) :: DIST2
    REAL(r8), DIMENSION(3,SIZE(offsets)) :: AECR
    REAL(r8), DIMENSION(3,SIZE(offsets)) :: AUX1
    REAL(r8), DIMENSION(3,SIZE(offsets)) :: AUX2
    REAL(r8), DIMENSION(6,SIZE(offsets)) :: AUX1ECI
    REAL(r8), DIMENSION(6,SIZE(offsets)) :: AUX2ECI
    REAL(r8), DIMENSION(6,SIZE(offsets)) :: AUX1ECR
    REAL(r8), DIMENSION(6,SIZE(offsets)) :: AUX2ECR
    REAL(r8), DIMENSION(SIZE(offsets)) :: QSIZE
    LOGICAL, DIMENSION(SIZE(offsets)) :: SMALL1
    LOGICAL, DIMENSION(SIZE(offsets)) :: SMALL2

    ! Executable code
    nV = SIZE(offsets)

    ! Construct two auxilliary vectors -- the first pointing directly ahead 
    ! along the s/c orbit, and the second pointing 45 degrees downward.
    aux1(1,:) = 1.0
    aux1(2,:) = 0.0
    aux1(3,:) = 0.0
    aux2(1,:) = 1.0
    aux2(2,:) = 0.0
    aux2(3,:) = 1.0

    ! Transform the auxilliary vectors from orbital to ECI coordinates
    returnStatus = Pgs_csc_orbToECI (spacecraftId, nV, asciiUTC, &
      offsets, aux1, aux1ECI(1:3,:) )
    returnStatus = Pgs_csc_orbToECI (spacecraftId, nV, asciiUTC, &
      offsets, aux2, aux2ECI(1:3,:) )

    ! Transform again to ECI coordinates
    aux1ECI(4:6,:) = 0.0_r8
    aux2ECI(4:6,:) = 0.0_r8
    returnStatus = Pgs_csc_ECItoECR (nV, asciiUTC, offsets, aux1ECI, aux1ECR)
    returnStatus = Pgs_csc_ECItoECR (nV, asciiUTC, offsets, aux2ECI, aux2ECR)

    ! Define l & a, where l is the distance to the equator along the direction
    ! of one of the auxilliary vectors from the s/c, and a is the auxilliary
    ! vector for which l was defined.

    small1 = ABS (aux1ECR(3,:)) < SQRT (TINY(0.0_r8))
    small2 = ABS (aux2ECR(3,:)) < SQRT (TINY(0.0_r8))
    IF (ANY (small1 .AND. small2)) THEN
      CALL MLSMessage (MLSMSG_Error, ModuleName, &
           'Problem computing auxilliary vector for master coordinate')
    END IF

    WHERE ( small1 )
      ! If aux1 has a zero z-component, set l & a for aux2, if its z-component
      !  /= 0
      l = -scECR(3,:)/aux2ECR(3,:)
      aECR(1,:) = aux2ECR(1,:)
      aECR(2,:) = aux2ECR(2,:)
      aECR(3,:) = aux2ECR(3,:)
    ELSEWHERE ( small2 )
      ! If aux1(3) is not zero, but aux2(3) is, then set l & a for aux1
      l = -scECR(3,:)/aux1ECR(3,:)
      aECR(1,:) = aux1ECR(1,:)
      aECR(2,:) = aux1ECR(2,:)
      aECR(3,:) = aux1ECR(3,:)
    ELSEWHERE
      ! If both aux1(3) and aux2(3) are non-zero, choose the one which gives the
      ! minimum absolute value of l.
      dist1 = -scECR(3,:)/aux1ECR(3,:)
      dist2 = -scECR(3,:)/aux2ECR(3,:)
      WHERE ( ABS(dist1) < ABS(dist2) )
        l = dist1
        aECR(1,:) = aux1ECR(1,:)
        aECR(2,:) = aux1ECR(2,:)
        aECR(3,:) = aux1ECR(3,:)
      ELSEWHERE
        l = dist2
        aECR(1,:) = aux2ECR(1,:)
        aECR(2,:) = aux2ECR(2,:)
        aECR(3,:) = aux2ECR(3,:)
      END WHERE
    END WHERE

    ! Define the vector q = scECR + la, such that its z-component = 0
    q(1,:) = scECR(1,:) + l(:)*aECR(1,:)
    q(2,:) = scECR(2,:) + l(:)*aECR(2,:)
    q(3,:) = 0.0

    ! Modify q, depending on which hemisphere the s/c is in, and the sign of l
    WHERE (((scECR(3,:) < 0.0) .AND. (l < 0.0)) .OR. &
         ((scECR(3,:) >= 0.0) .AND. (l >= 0.0)))
      q(1,:) = -q(1,:)
      q(2,:) = -q(2,:)
    END WHERE 

    ! Normalize q
    qSize = SQRT (q(1,:)**2 + q(2,:)**2) 

    q(1,:) = q(1,:) / qSize(:)
    q(2,:) = q(2,:) / qSize(:)

  END SUBROUTINE Mc_aux

  !---------------------------------------orbInclineCalculated -----------------
  FUNCTION orbInclineCalculated (scECR, scVelECR , lambda, mu)
    ! This function computes the orbital inclination angle beta' in degrees
    ! where 90 would mean a perfectly polar orbit in ECR coordinates
    ! Method: let [r] be the vector of the s/c position (in ECR coords)
    ! i.e., [r] = (x, y, z)
    ! also [v] its instantaneous velocity, 
    ! and scalar values include lambda its geocentric latitude
    ! and mu its longitude
    ! then                         (vy cos[mu] - vx sin[mu]
    ! then sin_beta' = cos[lambda] ------------------------
    !                                         |v|
    ! Arguments
    REAL(r8), INTENT(IN) :: scECR(3)             ! s/c pos.
    REAL(r8), INTENT(IN) :: scVelECR(3)          ! s/c vel.
    REAL, INTENT(IN) ::     lambda               ! s/c geocentric latitude
    REAL, INTENT(IN) ::     mu                   ! s/c longitude
    REAL(r8) ::             orbInclineCalculated

    ! Variables
    LOGICAL, PARAMETER :: DEBUG = .FALSE.
    INTEGER, SAVE :: HOWMANYSOFAR=0
    REAL(r8) :: vUnrotated(3)          ! s/c vel.
    REAL(r8) :: v
    REAL :: muRad, lamRad

    ! Executable code
    HOWMANYSOFAR = HOWMANYSOFAR + 1
    vUnrotated(1) = scVelECR(1) - Omega*scECR(2)
    vUnrotated(2) = scVelECR(2) + Omega*scECR(1)
    vUnrotated(3) = scVelECR(3)
    v = SQRT( vUnrotated(1)**2 + vUnrotated(2)**2 + vUnrotated(3)**2 )
    lamRad = (Pi/180)*lambda
    muRad = (Pi/180)*mu
    IF (v == 0.d0) THEN
       orbInclineCalculated = UNDEFINED_VALUE
       RETURN
    ENDIF
    orbInclineCalculated = (180/Pi) * ASIN( &
        COS(lamRad) * (vUnrotated(2)*COS(muRad) - vUnrotated(1)*SIN(muRad)) / v)
!    sc_velv = [sc_vel(i,j)       - 7.27221d-08*datascecr(1,i,j), $
!               sc_vel(i+mmifs,j) + 7.27221d-08*datascecr(0,i,j), $
!               sc_vel(i+2*mmifs,j)]
    ! Now added contraints: 90 < beta < 180
    orbInclineCalculated = ABS(orbInclineCalculated) + 90
    IF (orbInclineCalculated < 90.d0) THEN
      orbInclineCalculated = 180. - orbInclineCalculated
    ELSEIF (orbInclineCalculated > 180.d0) THEN
      orbInclineCalculated = 360. - orbInclineCalculated
    ENDIF
    
    IF (DEBUG) THEN
      CALL output ('vx, vy, vz ', advance='no')
      CALL blanks (3, advance='no')
      CALL output (scVelECR, advance='yes')
      CALL output ('lambda ', advance='no')
      CALL blanks (3, advance='no')
      CALL output (lambda, advance='no')
      CALL blanks (3, advance='no')
      CALL output ('mu ', advance='no')
      CALL blanks (3, advance='no')
      CALL output (mu, advance='yes')
      CALL output ('orbital inclination ', advance='no')
      CALL blanks (3, advance='no')
      CALL output ((180/Pi) * ASIN( &
           COS(lamRad) * (vUnrotated(2)*COS(muRad) - vUnrotated(1)*SIN(muRad)) &
           / v), advance='no')
      CALL blanks (3, advance='no')
      CALL output (orbInclineCalculated, advance='yes')
      IF (HOWMANYSOFAR > 40 ) STOP
    ENDIF

  END FUNCTION orbInclineCalculated

  !----------------------------------------orbInclineCrossProd -----------------
  FUNCTION orbInclineCrossProd  (scECI, scVelECI)
    ! This function computes the orbital inclination angle beta in degrees
    ! where 90 would mean a perfectly polar orbit
    ! Method: let [r] be the vector of the s/c position (in ECI coords)
    ! i.e., [r] = (x, y, z)
    ! also [v] its instantaneous velocity, [omega] its orbital frequency,
    ! we'll calculate its orbital moment [m]: using 'x' as the cross-product
    ! Then [v] = [omega] x [r]
    ! [p] = [r] x [v] = [omega] r^2 - [r] [omega] . [r]
    ! and if [r] . [omega] = 0
    ! [omega] = [p] / r^2 = omega (sin_beta cos_alfa, sin_beta sin_alfa, 
    ! cos_beta)
    ! Arguments
    REAL(r8), INTENT(IN) :: scECI(3)             ! s/c pos.
    REAL(r8), INTENT(IN) :: scVelECI(3)          ! s/c vel.
    REAL(r8) :: orbInclineCrossProd

    ! Variables
    LOGICAL, PARAMETER :: DEBUG = .FALSE.
    INTEGER, SAVE :: HOWMANYSOFAR=0
    REAL(r8) :: orbMoment(3), OMagnitude

    ! Executable code
    HOWMANYSOFAR = HOWMANYSOFAR + 1
    orbMoment(1) = scECI(2)*scVelECI(3) - scECI(3)*scVelECI(2)
    orbMoment(2) = scECI(3)*scVelECI(1) - scECI(1)*scVelECI(3)
    orbMoment(3) = scECI(1)*scVelECI(2) - scECI(2)*scVelECI(1)
    OMagnitude = SQRT (orbMoment(1)**2 + orbMoment(2)**2 + orbMoment(3)**2)
    IF (OMagnitude == 0.d0) THEN
       orbInclineCrossProd = UNDEFINED_VALUE
       RETURN
    ENDIF
    orbInclineCrossProd = (180/Pi) * ACOS (orbMoment(3) / OMagnitude)

    ! Now added contraints: 90 < beta < 180
    orbInclineCrossProd = ABS(orbInclineCrossProd)
    IF (orbInclineCrossProd < 90.d0) THEN
       orbInclineCrossProd = 180. - orbInclineCrossProd
    ELSEIF (orbInclineCrossProd > 180.d0) THEN
       orbInclineCrossProd = 360. - orbInclineCrossProd
    ENDIF
    
    IF (DEBUG) THEN
       CALL output ('rx, ry, rz ', advance='no')
       CALL blanks (3, advance='no')
       CALL output (scECI, advance='yes')
       CALL output ('vx, vy, vz ', advance='no')
       CALL blanks (3, advance='no')
       CALL output (scVelECI, advance='yes')
       CALL output ('px, py, pz ', advance='no')
       CALL blanks (3, advance='no')
       CALL output (orbMoment, advance='yes')
       CALL output ('orbital inclination ', advance='no')
       CALL blanks (3, advance='no')
       CALL output ((180/Pi) * ACOS( orbMoment(3) / OMagnitude ), advance='no')
       CALL blanks (3, advance='no')
       CALL output (orbInclineCrossProd, advance='yes')
       IF (HOWMANYSOFAR > 40 ) STOP
    ENDIF

  END FUNCTION orbInclineCrossProd

  !---------------------------------------------------TkL1B_mc -----------------
  SUBROUTINE TkL1B_mc (ascTAI, dscTAI, dotVec, nV, numOrb, &
    q, timeTAI, time_offset, geodAngle, scOrbIncl)
    ! This subroutine computes phi, the master coordinate for the spacecraft and
    ! tangent point records.
    ! Arguments
    INTEGER, INTENT(IN) :: nV, numOrb
    REAL(r8), INTENT(IN) :: timeTAI
    REAL(r8), INTENT(IN) :: q(3,nV)
    REAL(r8), INTENT(IN) :: time_offset(nV)
    REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:)
    REAL(r8), INTENT(IN) :: dotVec(3,nV)
    REAL, INTENT(OUT) ::    geodAngle(nV)
    REAL, INTENT(in) ::     scOrbIncl(nV)
!    real(r8), intent(IN) :: scECR(3,nV)             ! s/c pos.
!    real(r8), intent(IN) :: scVelECR(3,nV)          ! s/c vel.

    ! Functions
    INTEGER :: Pgs_csc_getEarthFigure

    ! Variables
    INTEGER :: i, j, returnStatus, scOrb

    REAL(r8) :: a, asciiTAI, b, cSq, equatRad_a, orbRad, polarRad_c
    REAL(r8) :: cosPhi(nV), gamma(nV), phi(nV), sinPhi(nV)
    REAL(r8) :: s(3,nV)
    REAL(r8) :: orbInclineNow, DvecSqrtSum

    ! Executable code

    ! Read a & b from earthfigure.dat
    returnStatus = Pgs_csc_getEarthFigure (earthModel, equatRad_a, polarRad_c)
    a = equatRad_a/1000
    b = polarRad_c/1000

    ! Set s = normalized dotVec
    DO i = 1, nV

       orbInclineNow = scOrbIncl(i)
       IF (orbInclineNow == UNDEFINED_VALUE) THEN
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               'Error in calculating orbital inclination angle')
       ENDIF

       ! Get a unit ECR vector to the point.
       DvecSqrtSum = SQRT(dotVec(1,i)**2 + dotVec(2,i)**2 + &
            dotVec(3,i)**2)
       IF (DvecSqrtSum == 0.0) THEN
          phi(i) = UNDEFINED_VALUE / Rad2Deg
          CYCLE
       ENDIF
       s(:,i) = dotVec(:,i) / DvecSqrtSum

       ! Calculate the geocentric angle as a number of radians between 0 and PI
       gamma(i) = ACOS( q(1,i)*s(1,i) + q(2,i)*s(2,i) )

       ! Place angle between PI and 2*PI, if in Southern Hemisphere
       IF (dotVec(3,i) < 0.0 ) gamma(i) = 2*PI - gamma(i)

       ! Going to convert this to geodetic, calculate som parameters.
       orbRad= Deg2Rad * (orbInclineNow - 90)
       cSq = (1 + (TAN(orbRad)**2)) * (a**2)*(b**2) / &
            (a**2 + (b**2)*(TAN(orbRad)**2))

       ! If |gamma| <= 45 deg of the equator, calculate phi using the SIN !
       ! equation
       IF ((gamma(i) <= PI/4 ) .OR. ( (gamma(i) >= 3*PI/4) .AND. &
            (gamma(i) <= 5*PI/4) ) .OR. (gamma(i) >= 7*PI/4)) THEN
          sinPhi(i) = SQRT( (a**4)*(SIN(gamma(i))**2)/( (cSq**2)* &
               (COS(gamma(i))**2) + (a**4)*(SIN(gamma(i))**2) ) )
          phi(i) = ASIN(sinPhi(i))
       ELSE
          ! If gamma is within 45 deg of a pole, calculate phi using the COS 
          ! equation
          cosPhi(i) = SQRT( (cSq**2)*(COS(gamma(i))**2)/( (cSq**2)* &
               COS(gamma(i))**2 + (a**4)*(SIN(gamma(i))**2) ) )
          phi(i) = ACOS(cosPhi(i))
       ENDIF

       ! Place phi in same quadrant as gamma
       IF ((gamma(i) > PI/2) .AND. (gamma(i) <= PI)) phi(i) = PI - phi(i)
       IF ((gamma(i) > PI) .AND. (gamma(i) <= 3*PI/2)) phi(i) = phi(i) + PI
       IF (gamma(i) > 3*PI/2) phi(i) = 2*PI - phi(i)

       ! Make phi cumulative over orbits
       ! First what is the time?
       asciiTAI = timeTAI + time_offset(i)
       scOrb = 0

       ! Now always make our calculation based on some threshold way back in
       ! time, to avoid any ambiguity
       IF ((phi(i) < pi/2) .OR. (phi(i) >= 3*pi/2)) THEN
          ! If in the ascending part of the orbit, base correction
          ! on time of last descending node
          ! MAKE THIS USE HUNT LATER!
          DO j = 1, numOrb
             IF (asciiTAI > dscTAI(j)) scOrb = j
          ENDDO

          ! If we're in the NH we've begun a new orbit, so add one
          IF (phi(i) < pi) scOrb = scOrb + 1
       ELSE
          ! If in descending part of the orbit, base it on the time
          ! of the last ascending node
          ! MAKE THIS USE HUNT LATER!
          DO j = 1, numOrb
             IF (asciiTAI > ascTAI(j)) scOrb = j
          ENDDO

       END IF
       phi(i) = (scOrb-2)*2*PI + phi(i)

    ENDDO

    ! Convert to degrees for output
    geodAngle = Rad2Deg * phi

  END SUBROUTINE TkL1B_mc

!=============================================================================
  SUBROUTINE Flag_Bright_Objects (TAI, ScAngle, MoonLimbTol, LimbView, &
       SpaceView)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    REAL(r8) :: TAI(:)
    REAL :: ScAngle(0:), MoonLimbTol
    TYPE (LOG_ARR1_PTR_T) :: LimbView(:)
    TYPE (LOG_ARR1_PTR_T), OPTIONAL :: SpaceView(:)

    CHARACTER (LEN=27) :: asciiUTC
    INTEGER :: i, MIF, returnStatus
    LOGICAL :: HasSpaceView
    REAL :: limb_angle, space_angle
    REAL(r8) :: offset (SIZE(TAI))
    REAL(r8) :: sc_frame_vector(3,0:(lenG-1))  ! start at MIF 0
    REAL(r8) :: sc_unit_vector(3)

    INTEGER, PARAMETER :: BO_defs(2) = (/ PGSd_Moon, PGSd_Venus /)
    REAL, PARAMETER :: VenusInSpace = 1.0
    REAL, PARAMETER :: VenusInLimb = 1.0
    REAL :: space_tol(2) = (/ 0.0, VenusInSpace /) ! tolerance for space port
    REAL :: limb_tol(2) = (/ 0.0, VenusInLimb /)   ! tolerance for limb port

    ! Functions

    INTEGER, EXTERNAL :: PGS_TD_taiToUTC, PGS_CBP_Sat_CB_Vector

    HasSpaceView = PRESENT (SpaceView)

! Get Moon tolerances from CF inputs:

    space_tol(1) = L1Config%Calib%MoonToSpaceAngle
!    limb_tol(1) = L1Config%Calib%MoonToLimbAngle
    limb_tol(1) = MoonLimbTol

    returnStatus = PGS_TD_taiToUTC (TAI(1), asciiUTC)
    offset = TAI - TAI(1)   ! offset (secs) from start TAI

    DO i = 1, 2
       IF (HasSpaceView) SpaceView(i)%ptr = .FALSE.
       LimbView(i)%ptr = .FALSE.
       returnStatus = PGS_CBP_Sat_CB_Vector (spacecraftId, lenG, asciiUTC, &
            offset(1:lenG), BO_defs(i), sc_frame_vector)
       IF (returnStatus /= 0) CYCLE

       DO MIF = 0, (lenG - 1)
          sc_unit_vector = sc_frame_vector(:,MIF) / &
               SQRT (sc_frame_vector(1,MIF)**2 + sc_frame_vector(2,MIF)**2 + &
               sc_frame_vector(3,MIF)**2)

! Space port angle check

          IF (HasSpaceView) THEN
             space_angle = ACOS (sc_unit_vector(2)) * Rad2Deg  ! Y vector
             SpaceView(i)%ptr(MIF) = (space_angle < space_tol(i))
          ENDIF

! Limb port angle check

          limb_angle = ACOS ( DOT_PRODUCT( &
               (/ COS (scAngle(MIF) * Deg2Rad), 0.0, &
               SIN (scAngle(MIF) * Deg2Rad) /), &
               sc_unit_vector) ) * Rad2Deg
          LimbView(i)%ptr(MIF) = (limb_angle < limb_tol(i))

       ENDDO

    ENDDO

  END SUBROUTINE Flag_Bright_Objects

END MODULE TkL1B

! $Log$
! Revision 2.21  2005/01/28 17:04:05  perun
! Pass in MoonToLimb tolerance instead of getting it from L1Config
!
! Revision 2.20  2005/01/25 18:01:06  perun
! Calculate tangent point scan rates (deg/sec)
!
! Revision 2.19  2004/11/10 15:32:05  perun
! Add encoder values; latest FOVs; correct YPR values order
!
! Revision 2.18  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.17  2004/08/03 20:41:14  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.16  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.15  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.14  2003/09/15 21:50:04  pwagner
! Removed illegal midline continuation--Lahey disapproves
! 
! Revision 2.13  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.12  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.11  2002/11/07 21:56:20  jdone
! Added Level 1 output datatypes.
!
! Revision 2.10  2002/09/26 20:52:26  vsnyder
! Get Omega from Geometry instead of Units
!
! Revision 2.9  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.8  2001/12/14 01:43:46  livesey
! Working version with ECR based master coordinate
!
! Revision 2.7  2001/12/12 19:05:44  livesey
! Fixed bug with master coordinate for very low latitudes.
! However, this version it turns out is basing master coordinate
! on the ECI orbital node, not ECR which I think is what we really
! want. So expect another version soon!
!
! Revision 2.6  2001/12/12 03:07:32  livesey
! Interim version needs debugging more
!
! Revision 2.5  2001/12/11 00:55:51  livesey
! Slightly kludgy fix for problem where occasionally gets
! phi 180 degrees out.
!
! Revision 2.4  2001/12/07 00:51:44  pwagner
! Finally calculates scOrbIncl correctly
!
! Revision 2.3  2001/12/06 01:03:46  pwagner
! Now writes orbit incline angle in ECR
!
! Revision 2.2  2001/10/12 22:11:05  livesey
! Tidied things up a bit, added scVelECR, but not filled yet
!
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
! Revision 2.0  2000/09/05 18:55:15  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/15 18:48:37  nakamura
! Absorbed module Mc; moved _init subroutine to Orbit; account for parametrization of lenG & lenT.
!
