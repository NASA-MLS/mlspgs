! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TkL1B

!!!!!! THIS LINE DOESN'T SEEM TO WORK!!!!!  use MLSNumerics, only: Hunt
  use Dump_0, only: DUMP
  use MLSCommon, only: R8
  use MLSL1Common
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_Error
  use OUTPUT_M, only: BLANKS, OUTPUT
  use OutputL1B, only: L1BOAsc_T, L1BOATP_T, L1BOAINDEX_T, LENCOORD, &
    OUTPUTL1B_THZ, OUTPUTL1B_SC, OUTPUTL1B_INDEX, OUTPUTL1B_GHZ, LENG, LENT
  use Scan
  use SDPToolkit
  use units, only: Omega
  implicit none
  private

  public :: TKL1B_SC, TKL1B_TP, L1BOA_MAF, MC_AUX, TKL1B_MC
  logical, public, parameter :: ORBINCLINE_IS_CONSTANT = .FALSE.
  real, parameter ::    UNDEFINED_VALUE = -999.99

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  character (LEN=*), parameter :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

  ! This module contains subroutines for producing the L1BOA records on
  ! a MAF by MAF basis.

contains

  !------------------------------------------ TkL1B_sc ---------
  subroutine TkL1B_sc(numValues, offsets, asciiUTC, sc)
    ! This subroutine contains prototype code for creating the desired s/c record
    ! from the EPHEMATTIT output.

    ! Arguments
    type( L1BOAsc_T ) :: sc
    character (LEN=27), intent(IN) :: asciiUTC
    integer, intent(IN) :: numValues
    real(r8), intent(IN) :: offsets(numValues)

    ! Functions
    integer :: Pgs_csc_eciToECR, Pgs_csc_ecrToGEO, Pgs_eph_ephemAttit

    ! Variables
    character (LEN=32) :: mnemonic
    character (LEN=480) :: msg, msr
    integer :: i, returnStatus
    integer :: qualityFlags(2, numValues)
    real(r8) :: ecrVec(3)
    real(r8) :: radC(numValues), radD(numValues), radL(numValues)
    real(r8) :: attitQuat(4,numValues)
    real(r8) :: eciV(6,numValues), ecrV(6,numValues)

    ! Executable code

    ! Read oa data
    returnStatus = Pgs_eph_ephemAttit (spacecraftId, numValues, asciiUTC,  &
      offsets, pgs_true, pgs_true, qualityFlags, &
      sc%scECI, sc%scVelECI, sc%ypr, sc%yprRate, attitQuat)
    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      msr = 'Routine ephemAttit, ' // mnemonic // ':  ' // msg
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    ! Convert scECI to scECR
    eciV(1:3,:) = sc%scECI
    eciV(4:6,:) = sc%scVelECI    ! was 0.0
    returnStatus = Pgs_csc_eciToECR (numValues, asciiUTC, offsets, eciV, ecrV)
    sc%scECR = ecrV(1:3,:)
    sc%scVelECR = ecrV(4:6,:)

    ! Calculate geocentric/geodetic altitude, latitude & longitude from scECR
    do i = 1, numValues
      sc%scGeocAlt(i) = sqrt( ecrV(1,i)**2 + ecrV(2,i)**2 + ecrV(3,i)**2 )
      radC(i) = atan( ecrV(3,i) / sqrt(ecrV(1,i)**2 + ecrV(2,i)**2) )
      sc%scGeocLat(i) = Rad2Deg * radC(i)
      ecrVec = ecrV(1:3, i)
      returnStatus = Pgs_csc_ecrToGEO (ecrVec, earthModel, radL(i), &
        radD(i), sc%scGeodAlt(i))
      sc%scGeodLat(i) = Rad2Deg * radD(i)
      sc%scLon(i) = Rad2Deg * radL(i)
      if ( ORBINCLINE_IS_CONSTANT ) then
        sc%scOrbIncl(i) = orbInclineCrossProd(sc%scECI(:,i), sc%scVelECI(:,i))
      else
        sc%scOrbIncl(i) = orbInclineCalculated(sc%scECR(:,i), sc%scVelECR(:,i), &
        & sc%scGeocLat(i), sc%scLon(i) )
      endif
    enddo


  end subroutine TkL1B_sc

  !------------------------------------------------------TkL1B_tp ----
  subroutine TkL1B_tp(asciiTAI, asciiUTC, lenG, numValues, offsets, posECR, &
    scRate, startAngle, tp)
    ! This subroutine fills the tangent point record.

    ! Arguments
    type( L1BOAtp_T ) :: tp
    character (LEN=27), intent(IN) :: asciiUTC
    integer, intent(IN) :: lenG, numValues
    real, intent(IN) :: scRate(lenG)
    real(r8), intent(IN) :: asciiTAI, startAngle
    real(r8), intent(IN) :: offsets(numValues)
    real(r8), intent(IN) :: posECR(3,lenG)

    ! Functions
    integer :: Pgs_cbp_sat_cb_vector, Pgs_cbp_solarTimeCoords
    integer :: Pgs_csc_scToOrb, Pgs_csc_scToECI, Pgs_csc_eciToECR
    integer :: Pgs_csc_grazingRay, Pgs_csc_ecrToECI, Pgs_csc_eciToOrb
    integer :: Pgs_td_taiToUTC

    ! Variables
    character (LEN=27) :: time
    integer :: flag, flagQ, i, returnStatus
    real(r8) :: declination, delAngle, delTime, deltaAlt, deltaLat, deltaLon
    real(r8) :: greenwich, localApparent, rightAscension, tai
    real(r8) :: dot(lenG), latD(lenG), localMean(lenG), lon(lenG), los(lenG)
    real(r8) :: sign(lenG), slantRange(lenG)
    real(r8) :: angleRad(numValues)
    real(r8) :: eci(3,lenG), ecr(3,lenG), hECR(3,lenG), nts(3,lenG)
    real(r8) :: posSurf(3,lenG), sc_frame_vector(3,lenG), sc_sun(3,lenG)
    real(r8) :: sc_tp(3,lenG), tp_sun(3,lenG), tpOrb(3,lenG), unitAlt(3,lenG)
    real(r8) :: unitLat(3,lenG), unitLon(3,lenG), vECR(3,lenG)
    real(r8) :: angleSc(3,numValues), angleOrb(3,numValues)
    real(r8) :: eciV(6,lenG), ecrV(6,lenG)
    ! Executable code

    deltaLat = 0.1
    deltaLon = 0.1
    deltaAlt = 1000.0

    ! Calculate MIF scan angles (in DEGREES)
    tp%scAngle(1) = startAngle
    do i = 2, lenG
      tp%scAngle(i) = tp%scAngle(i-1) - scRate(i)*offsets(2)
    enddo

    ! Calculate retrace angle, rate
    delAngle = tp%scAngle(lenG) - tp%scAngle(1)
    delTime = (numValues - lenG)*offsets(2)
    tp%scanRate(1:lenG) = scRate
    tp%scanRate( (lenG+1):numValues ) = delAngle/delTime
    do i = lenG+1, numvalues
      tp%scAngle(i) = tp%scAngle(i-1) - tp%scanRate(i)*offsets(2)
    enddo

    ! Put angle in s/c coordinates
    angleRad = Deg2Rad * tp%scAngle
    angleSc(1,:) = cos(angleRad)
    angleSc(2,:) = 0.0
    angleSc(3,:) = sin(angleRad)

    ! Convert s/c vector to Orb vector/angle/degrees
    returnStatus = Pgs_csc_scToOrb(spacecraftId, numValues, asciiUTC, &
      offsets, angleSc, angleOrb)

    tp%scanAngle = Rad2Deg * acos( angleOrb(1,:) )

    ! Convert s/c vector to ECR
    returnStatus = Pgs_csc_scToECI(spacecraftId, lenG, asciiUTC, &
      offsets(1:lenG), angleSc(:,1:lenG), eci)
    eciV(1:3,:) = eci
    eciV(4:6,:) = 0.0
    returnStatus = Pgs_csc_eciToECR(lenG, asciiUTC, offsets(1:lenG), eciV, &
      ecrV)
    ecr = ecrV(1:3,:)

    ! For each scanning MIF,
    do i = 1, lenG

      ! Calculate tangent point (geodetic & ECR)
      returnStatus = Pgs_csc_grazingRay(earthModel, posECR(:,i), ecr(:,i), &
        latD(i), lon(i), tp%tpGeodAlt(i), &
        slantRange(i), tp%tpECR(:,i), posSurf(:,i) )

      ! Create ECR unit vector quantities -- lat=1, lon=2, alt=3
      flagQ = 2
      call Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLon, &
        unitLon(:,i), flag)
      flagQ = 1
      call Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLat, &
        unitLat(:,i), flag)
      flagQ = 3
      call Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaAlt, &
        unitAlt(:,i), flag)

      ! Get local mean solar time from Toolkit
      tai = asciiTAI + (i-1)*offsets(2)
      returnStatus = Pgs_td_taiToUTC(tai, time)
      returnStatus = Pgs_cbp_solarTimeCoords(time, lon(i), greenwich, &
        localMean(i), localApparent, rightAscension, declination)
    enddo

    tp%tpGeodLat = Rad2Deg * latD
    tp%tpLon = Rad2Deg * lon
    tp%tpSolarTime = localMean/3600.0

    ! Calculate solarZenith
    returnStatus = Pgs_cbp_sat_cb_vector(spacecraftId, lenG, asciiUTC, &
      offsets(1:lenG), PGSd_SUN, sc_frame_vector)
    returnStatus = Pgs_csc_scToECI(spacecraftId, lenG, asciiUTC, &
      offsets(1:lenG), sc_frame_vector, eci)
    eciV(1:3,:) = eci
    eciV(4:6,:) = 0.0
    returnStatus = Pgs_csc_eciToECR(lenG, asciiUTC, offsets(1:lenG), eciV, &
      ecrV)
    sc_sun = ecrV(1:3,:)
    do i = 1, lenG
      sc_tp(:,i) = ecr(:,i) * slantRange(i)
    enddo
    tp_sun = sc_sun - sc_tp
    do i = 1, lenG
      nts(:,i) = tp_sun(:,i) / sqrt( tp_sun(1,i)**2 + tp_sun(2,i)**2 + &
        &tp_sun(3,i)**2 )
    enddo
    dot = nts(1,:)*unitAlt(1,:) + nts(2,:)*unitAlt(2,:) + &
      nts(3,:)*unitAlt(3,:)
    tp%tpSolarZenith = acos(dot) * Rad2Deg

    ! Calculate losAngle
    do i = 1, lenG
      vECR(:,i) = ecr(:,i) - ( ecr(1,i)*unitAlt(1,i) + &
        &ecr(2,i)*unitAlt(2,i) + ecr(3,i)*unitAlt(3,i) )
      hECR(:,i) = vECR(:,i) / sqrt(vECR(1,i)**2 + vECR(2,i)**2 + &
        &vECR(3,i)**2)
    enddo
    los = acos( hECR(1,:)*unitLat(1,:) + hECR(2,:)*unitLat(2,:) + &
      &hECR(3,:)*unitLat(3,:) )
    sign = ecr(1,:)*unitLon(1,:) + ecr(2,:)*unitLon(2,:) + &
      ecr(3,:)*unitLon(3,:)
    do i = 1, lenG
      if (sign(i) < 0 ) los(i) = 2*PI - los(i)
    enddo
    tp%tpLosAngle = los * Rad2Deg

    ! Convert tpECR to tpECI
    ecrV(1:3,:) = tp%tpECR
    ecrV(4:6,:) = 0.0
    returnStatus = Pgs_csc_ecrToECI (lenG, asciiUTC, offsets(1:lenG), ecrV, &
      eciV)
    tp%tpECI = eciV(1:3,:)

    ! Convert tpECI to tpOrb
    returnStatus = Pgs_csc_eciToOrb(spacecraftId, lenG, asciiUTC, &
      offsets(1:lenG), tp%tpECI, tpOrb)
    tp%tpOrbY = tpOrb(2,:)

    ! Calculate tp geocentric coordinates
    do i = 1, lenG
      tp%tpGeocAlt(i) = sqrt( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2 + &
        &tp%tpECR(3,i)**2 )
      tp%tpGeocLat(i) = Rad2Deg * atan( tp%tpECR(3,i) &
        &/ sqrt( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2 ) )
    enddo

    ! Calculate dummy values
    tp%encoderAngle = tp%scAngle
    do i = 2, lenG
      tp%tpGeocAltRate(i-1) = ( tp%tpGeocAlt(i) - tp%tpGeocAlt(i-1) )&
        &/offsets(2)
      tp%tpGeodAltRate(i-1) = ( tp%tpGeodAlt(i) - tp%tpGeodAlt(i-1) )&
        &/offsets(2)
    enddo
    tp%tpGeocAltRate(lenG) = tp%tpGeocAltRate(lenG-1)
    tp%tpGeodAltRate(lenG) = tp%tpGeodAltRate(lenG-1)

  end subroutine TkL1B_tp

  !------------------------------------------------- Tp_unit ---------
  subroutine Tp_unit (flagQ, lon, lat, alt, delta, unitQ, flag)
    ! This subroutine creates unit ECR vector quantities used by TkL1B_tp to
    ! calculate solarZenith and losAngle.
    ! Arguments
    integer, intent(IN) :: flagQ
    real(r8), intent(IN) :: alt, delta, lat, lon
    integer, intent(OUT) :: flag
    real(r8), intent(OUT) :: unitQ(3)

    ! Functions
    integer :: Pgs_csc_geoToECR

    ! Variables
    character (LEN=32) :: mnemonic
    character (LEN=480) :: msg, msr
    integer :: returnStatus
    real(r8) :: del, geoMinus, geoPlus
    real(r8) :: ecrMinus(3), ecrPlus(3), vec(3)

    ! Exectuable code

    flag = 0

    ! Create unit vector
    if (flagQ == 2) then
      ! longitude
      del = delta * Deg2Rad
      geoPlus  = (lon + del/2)
      geoMinus = (lon - del/2)
      if (geoPlus  >  PI) geoPlus  = geoPlus  - 2*PI
      if (geoMinus < -PI) geoMinus = geoMinus + 2*PI
      returnStatus = Pgs_csc_geoToECR(geoPlus, lat, alt, earthModel, &
        ecrPlus)
      returnStatus = Pgs_csc_geoToECR(geoMinus, lat, alt, earthModel, &
        ecrMinus)
    else if (flagQ == 1) then
      ! latitude
      del = delta * Deg2Rad
      geoPlus  = (lat + del/2)
      geoMinus = (lat - del/2)
      returnStatus = Pgs_csc_geoToECR(lon, geoPlus, alt, earthModel, &
        ecrPlus)
      returnStatus = Pgs_csc_geoToECR(lon, geoMinus, alt, earthModel, &
        ecrMinus)
    else
      geoPlus  = (alt + delta/2)
      geoMinus = (alt - delta/2)
      returnStatus = Pgs_csc_geoToECR(lon, lat, geoPlus, earthModel, &
        ecrPlus)
      returnStatus = Pgs_csc_geoToECR(lon, lat, geoMinus, earthModel, &
        ecrMinus)
    endif

    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      msr = mnemonic // ':  ' // msg
      call MLSMessage(MLSMSG_Warning, ModuleName, msr)
      flag = -1
    endif

    vec = ecrPlus - ecrMinus
    unitQ = vec / sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

  end subroutine Tp_unit

  !--------------------------------------------------- L1BOA_MAF ----------------
  subroutine L1boa_MAF(altG, altT, ascTAI, counterMAF, dscTAI, L1FileHandle, &
    MAFinfo, noMAF, numOrb, orbIncline, orbitNumber, &
    scRate, scRateT)
    ! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
    ! HDF output file.

    ! Arguments
    type (MAFinfo_T) :: MAFinfo
    integer, intent(IN) :: L1FileHandle, counterMAF, noMAF, numOrb
    integer, intent(IN) :: orbitNumber(:)
    real, intent(IN) :: scRate(:)
    real, intent(IN) :: scRateT(:)
    real(r8), intent(IN) :: altG, altT, orbIncline
    real(r8), intent(IN) :: ascTAI(:), dscTAI(:)

    ! Functions
    integer :: Pgs_td_taiToUTC

    ! Variables
    type( L1BOAindex_T ) :: index
    type( L1BOAsc_T ) :: sc
    type( L1BOAtp_T ) :: tp
    character (LEN=27) :: mafTime
    character (LEN=480) :: msr
    integer :: error, i, nV, returnStatus
    real(r8) :: mafTAI
    real(r8) :: initRay(3), q(3,MAFinfo%MIFsPerMAF)
    real(r8) :: offsets(MAFinfo%MIFsPerMAF)

    ! Calculate offsets array
    mafTAI = MAFinfo%startTAI
    nV = MAFinfo%MIFsPerMAF
    do i = 1, nV 
      offsets(i) = (i-1)*MAFinfo%integTime
    enddo

    ! Write "index" information
    returnStatus = Pgs_td_taiToUTC(mafTAI, mafTime)
    index%MAFStartTimeUTC = mafTime
    index%MAFStartTimeTAI = mafTAI
    index%noMIFs = nV
    index%counterMAF = counterMAF
    call OutputL1B_index(noMAF, L1FileHandle, index)

    ! Allocate the MIF variables in the output structures
    allocate(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
      sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
      sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), sc%scOrbIncl(nV), &
      sc%scVelECI(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
      sc%scVelECR(lenCoord,nV), &
      tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
      tp%scanRate(nV), STAT=error)
    if ( error /= 0 ) then
      msr = MLSMSG_Allocate // '  s/c quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    ! Get oa data
    call TkL1B_sc(nV, offsets, mafTime, sc)

    ! Get s/c master coordinate
    call Mc_aux(mafTime, offsets, sc%scECR, q )

    call TkL1B_mc(ascTAI, dscTAI, sc%scECR, nV, numOrb, &
      & orbIncline, orbitNumber, q, mafTAI, offsets, sc%scGeodAngle, &
      & sc%scOrbIncl)

    ! Write s/c information
    call OutputL1B_sc(noMAF, L1FileHandle, sc)

    ! Calculate initial guess for look vector in ECR
    call Scan_guess(mafTime, initRay)

    ! Allocate the output structure

    allocate(tp%tpECI(lenCoord,lenG), tp%tpECR(lenCoord,lenG), &
      tp%tpOrbY(lenG), tp%tpGeocAlt(lenG), tp%tpGeocLat(lenG), &
      tp%tpGeocAltRate(lenG), tp%tpGeodAlt(lenG), &
      tp%tpGeodLat(lenG), tp%tpGeodAltRate(lenG), tp%tpLon(lenG), &
      tp%tpGeodAngle(lenG), tp%tpSolarTime(lenG), &
      tp%tpSolarZenith(lenG), tp%tpLosAngle(lenG), STAT=error)
    if ( error /= 0 ) then
      msr = MLSMSG_Allocate // '  GHz tp quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    ! Find angle, tan pt for start of GHZ scan
    call Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

    ! Calculate GHZ tan pt record
    call TkL1B_tp(mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenG), &
      scRate, tp%scAngle(1), tp)

    ! Compute GHz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECR, lenG, numOrb, &
      & orbIncline, orbitNumber, q, mafTAI, offsets(1:lenG), tp%tpGeodAngle, &
      & sc%scOrbIncl)

    ! Write GHz information
    call OutputL1B_GHz(noMAF, L1FileHandle, tp)

    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=error)
    if ( error /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
      &deallocation of GHz quantities.')

    ! Find angle, tan pt for start of THz scan
    allocate(tp%tpECI(lenCoord,lenT), tp%tpECR(lenCoord,lenT), &
      tp%tpOrbY(lenT), tp%tpGeocAlt(lenT), tp%tpGeocLat(lenT), &
      tp%tpGeocAltRate(lenT), tp%tpGeodAlt(lenT), &
      tp%tpGeodLat(lenT), tp%tpGeodAltRate(lenT), tp%tpLon(lenT), &
      tp%tpGeodAngle(lenT), tp%tpSolarTime(lenT), &
      tp%tpSolarZenith(lenT), tp%tpLosAngle(lenT), STAT=error)
    if ( error /= 0 ) then
      msr = MLSMSG_Allocate // '  THz tp quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    ! Find angle, tan pt for start of THz scan
    call Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

    ! Calculate THz tan pt record
    call TkL1B_tp(mafTAI, mafTime, lenT, nV, offsets, sc%scECR(:,1:lenT), &
      scRateT, tp%scAngle(1), tp)

    ! Compute THz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECR, lenT, numOrb, &
      & orbIncline, orbitNumber, q, mafTAI, offsets(1:lenT), tp%tpGeodAngle, &
      & sc%scOrbIncl)

    ! Write THZ information
    call OutputL1B_THz(noMAF, L1FileHandle, tp)

    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=error)
    if ( error /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed deallocation of THz quantities.')

    ! Deallocate the MIF quantities
    deallocate(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
      sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle, sc%scVelECI, &
      sc%scOrbIncl, sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, &
      tp%scanAngle, tp%scanRate, STAT=error)
    if ( error /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed deallocation of MIF quantities.')

  end subroutine L1boa_MAF

  !------------------------------------------- Mc_Aux --------
  subroutine Mc_aux (asciiUTC, offsets, scECR, q)
    ! This subroutine computes q, an auxilliary vector used in the calculation of
    ! the master coordinate.  Q is a vector that points from the center of the Earth
    ! to the ascending node of the orbit in ECI coordinates.  Thus any point
    ! dotted with q can give you the master coordinate.
    ! Arguments
    character (LEN=27), intent(IN) :: asciiUTC
    real(r8), intent(in) :: OFFSETS(:)
    real(r8), intent(IN) :: scECR(:,:)
    real(r8), intent(OUT) :: q(:,:)

    ! Parameters

    ! Functions
    integer :: Pgs_csc_orbToECI, PGS_CSC_ECItoECR

    ! Variables
    integer :: returnStatus
    integer :: nV
    real(r8), dimension(size(offsets)) :: l
    real(r8), dimension(size(offsets)) :: DIST1
    real(r8), dimension(size(offsets)) :: DIST2
    real(r8), dimension(3,size(offsets)) :: AECR
    real(r8), dimension(3,size(offsets)) :: AUX1
    real(r8), dimension(3,size(offsets)) :: AUX2
    real(r8), dimension(6,size(offsets)) :: AUX1ECI
    real(r8), dimension(6,size(offsets)) :: AUX2ECI
    real(r8), dimension(6,size(offsets)) :: AUX1ECR
    real(r8), dimension(6,size(offsets)) :: AUX2ECR
    real(r8), dimension(size(offsets)) :: QSIZE
    logical, dimension(size(offsets)) :: SMALL1
    logical, dimension(size(offsets)) :: SMALL2

    ! Executable code
    nV = size(offsets)

    ! Construct two auxilliary vectors -- the first pointing directly ahead 
    ! along the s/c orbit, and the second pointing 45 degrees downward.
    aux1(1,:) = 1.0
    aux1(2,:) = 0.0
    aux1(3,:) = 0.0
    aux2(1,:) = 1.0
    aux2(2,:) = 0.0
    aux2(3,:) = 1.0

    ! Transform the auxilliary vectors from orbital to ECI coordinates
    returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
      offsets, aux1, aux1ECI(1:3,:) )
    returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
      offsets, aux2, aux2ECI(1:3,:) )

    ! Transform again to ECI coordinates
    aux1ECI(4:6,:) = 0.0_r8
    aux2ECI(4:6,:) = 0.0_r8
    returnStatus = Pgs_csc_ECItoECR(nV, asciiUTC, offsets, &
      & aux1ECI, aux1ECR )
    returnStatus = Pgs_csc_ECItoECR(nV, asciiUTC, offsets, &
      & aux2ECI, aux2ECR )

    ! Define l & a, where l is the distance to the equator along the direction of
    ! one of the auxilliary vectors from the s/c, and a is the auxilliary vector
    ! for which l was defined.

    small1 = abs(aux1ECR(3,:)) < sqrt( tiny(0.0_r8) )
    small2 = abs(aux2ECR(3,:)) < sqrt( tiny(0.0_r8) )
    if ( any ( small1 .and. small2 ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Problem computing auxilliary vector for master coordinate' )
    end if

    where ( small1 )
      ! If aux1 has a zero z-component, set l & a for aux2, if its z-component /= 0
      l = -scECR(3,:)/aux2ECR(3,:)
      aECR(1,:) = aux2ECR(1,:)
      aECR(2,:) = aux2ECR(2,:)
      aECR(3,:) = aux2ECR(3,:)
    elsewhere ( small2 )
      ! If aux1(3) is not zero, but aux2(3) is, then set l & a for aux1
      l = -scECR(3,:)/aux1ECR(3,:)
      aECR(1,:) = aux1ECR(1,:)
      aECR(2,:) = aux1ECR(2,:)
      aECR(3,:) = aux1ECR(3,:)
    elsewhere
      ! If both aux1(3) and aux2(3) are non-zero, choose the one which gives the
      ! minimum absolute value of l.
      dist1 = -scECR(3,:)/aux1ECR(3,:)
      dist2 = -scECR(3,:)/aux2ECR(3,:)
      where ( abs(dist1) < abs(dist2) )
        l = dist1
        aECR(1,:) = aux1ECR(1,:)
        aECR(2,:) = aux1ECR(2,:)
        aECR(3,:) = aux1ECR(3,:)
      elsewhere
        l = dist2
        aECR(1,:) = aux2ECR(1,:)
        aECR(2,:) = aux2ECR(2,:)
        aECR(3,:) = aux2ECR(3,:)
      end where
    end where

    ! Define the vector q = scECR + la, such that its z-component = 0
    q(1,:) = scECR(1,:) + l(:)*aECR(1,:)
    q(2,:) = scECR(2,:) + l(:)*aECR(2,:)
    q(3,:) = 0.0

    ! Modify q, depending on which hemisphere the s/c is in, and the sign of l
    where ( ( (scECR(3,:) < 0.0) .and. (l < 0.0) ) .or. &
      &     ( (scECR(3,:) >= 0.0) .and. (l >= 0.0) ) )
      q(1,:) = -q(1,:)
      q(2,:) = -q(2,:)
    end where 

    ! Normalize q
    qSize = sqrt( q(1,:)**2 + q(2,:)**2 ) 

    q(1,:) = q(1,:) / qSize(:)
    q(2,:) = q(2,:) / qSize(:)

  end subroutine Mc_aux

  !----------------------------------------------------orbInclineCalculated -----------------
  function orbInclineCalculated( scECR, scVelECR , lambda, mu)
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
    real(r8), intent(IN) :: scECR(3)             ! s/c pos.
    real(r8), intent(IN) :: scVelECR(3)          ! s/c vel.
    real, intent(IN) ::     lambda               ! s/c geocentric latitude
    real, intent(IN) ::     mu                   ! s/c longitude
    real(r8) ::             orbInclineCalculated

    ! Variables
    logical, parameter :: DEBUG = .FALSE.
    integer, save :: HOWMANYSOFAR=0
    real(r8) :: vUnrotated(3)          ! s/c vel.
    real(r8) :: v
    real :: muRad, lamRad

    ! Executable code
    HOWMANYSOFAR = HOWMANYSOFAR + 1
    vUnrotated(1) = scVelECR(1) - Omega*scECR(2)
    vUnrotated(2) = scVelECR(2) + Omega*scECR(1)
    vUnrotated(3) = scVelECR(3)
    v = sqrt( vUnrotated(1)**2 + vUnrotated(2)**2 + vUnrotated(3)**2 )
    lamRad = (Pi/180)*lambda
    muRad = (Pi/180)*mu
    if ( v == 0.d0 ) then
      orbInclineCalculated = UNDEFINED_VALUE
      return
    endif
    orbInclineCalculated = (180/Pi) * asin( &
    & cos(lamRad) * (vUnrotated(2)*cos(muRad) - vUnrotated(1)*sin(muRad)) &
    & / &
    & v &
    & )
!    sc_velv = [sc_vel(i,j)       - 7.27221d-08*datascecr(1,i,j), $
!               sc_vel(i+mmifs,j) + 7.27221d-08*datascecr(0,i,j), $
!               sc_vel(i+2*mmifs,j)]
    ! Now added contraints: 90 < beta < 180
    orbInclineCalculated = abs(orbInclineCalculated) + 90
    if ( orbInclineCalculated < 90.d0 ) then
      orbInclineCalculated = 180. - orbInclineCalculated
    elseif (orbInclineCalculated > 180.d0 ) then
      orbInclineCalculated = 360. - orbInclineCalculated
    endif
    
    if ( DEBUG ) then
      call output('vx, vy, vz ', advance='no')
      call blanks(3, advance='no')
      call output(scVelECR, advance='yes')
      call output('lambda ', advance='no')
      call blanks(3, advance='no')
      call output(lambda, advance='no')
      call blanks(3, advance='no')
      call output('mu ', advance='no')
      call blanks(3, advance='no')
      call output(mu, advance='yes')
      call output('orbital inclination ', advance='no')
      call blanks(3, advance='no')
      call output((180/Pi) * asin( &
       & cos(lamRad) * (vUnrotated(2)*cos(muRad) - vUnrotated(1)*sin(muRad)) &
       & / &
       & v &
       & ), advance='no')
      call blanks(3, advance='no')
      call output(orbInclineCalculated, advance='yes')
      if ( HOWMANYSOFAR > 40 ) stop
    endif

  end function orbInclineCalculated

  !----------------------------------------------------orbInclineCrossProd -----------------
  function orbInclineCrossProd( scECI, scVelECI )
    ! This function computes the orbital inclination angle beta in degrees
    ! where 90 would mean a perfectly polar orbit
    ! Method: let [r] be the vector of the s/c position (in ECI coords)
    ! i.e., [r] = (x, y, z)
    ! also [v] its instantaneous velocity, [omega] its orbital frequency,
    ! we'll calculate its orbital moment [m]: using 'x' as the cross-product
    ! Then [v] = [omega] x [r]
    ! [p] = [r] x [v] = [omega] r^2 - [r] [omega] . [r]
    ! and if [r] . [omega] = 0
    ! [omega] = [p] / r^2 = omega (sin_beta cos_alfa, sin_beta sin_alfa, cos_beta)
    ! Arguments
    real(r8), intent(IN) :: scECI(3)             ! s/c pos.
    real(r8), intent(IN) :: scVelECI(3)          ! s/c vel.
    real(r8) :: orbInclineCrossProd

    ! Variables
    logical, parameter :: DEBUG = .FALSE.
    integer, save :: HOWMANYSOFAR=0
    real(r8) :: orbMoment(3), OMagnitude

    ! Executable code
    HOWMANYSOFAR = HOWMANYSOFAR + 1
    orbMoment(1) = scECI(2)*scVelECI(3) - scECI(3)*scVelECI(2)
    orbMoment(2) = scECI(3)*scVelECI(1) - scECI(1)*scVelECI(3)
    orbMoment(3) = scECI(1)*scVelECI(2) - scECI(2)*scVelECI(1)
    OMagnitude = sqrt( orbMoment(1)**2 + orbMoment(2)**2 + orbMoment(3)**2 )
    if ( OMagnitude == 0.d0 ) then
      orbInclineCrossProd = UNDEFINED_VALUE
      return
    endif
    orbInclineCrossProd = (180/Pi) * acos( orbMoment(3) / OMagnitude )

    ! Now added contraints: 90 < beta < 180
    orbInclineCrossProd = abs(orbInclineCrossProd)
    if ( orbInclineCrossProd < 90.d0 ) then
      orbInclineCrossProd = 180. - orbInclineCrossProd
    elseif (orbInclineCrossProd > 180.d0 ) then
      orbInclineCrossProd = 360. - orbInclineCrossProd
    endif
    
    if ( DEBUG ) then
      call output('rx, ry, rz ', advance='no')
      call blanks(3, advance='no')
      call output(scECI, advance='yes')
      call output('vx, vy, vz ', advance='no')
      call blanks(3, advance='no')
      call output(scVelECI, advance='yes')
      call output('px, py, pz ', advance='no')
      call blanks(3, advance='no')
      call output(orbMoment, advance='yes')
      call output('orbital inclination ', advance='no')
      call blanks(3, advance='no')
      call output((180/Pi) * acos( orbMoment(3) / OMagnitude ), advance='no')
      call blanks(3, advance='no')
      call output(orbInclineCrossProd, advance='yes')
      if ( HOWMANYSOFAR > 40 ) stop
    endif
  end function orbInclineCrossProd

  !----------------------------------------------------TkL1B_mc -----------------
  subroutine TkL1B_mc(ascTAI, dscTAI, dotVec, nV, numOrb, &
    orbIncline, orbitNumber, q, timeTAI, time_offset, geodAngle, &
    & scOrbIncl)
    ! This subroutine computes phi, the master coordinate for the spacecraft and
    ! tangent point records.
    ! Arguments
    integer, intent(IN) :: nV, numOrb
    integer, intent(IN) :: orbitNumber(:)
    real(r8), intent(IN) :: orbIncline, timeTAI
    real(r8), intent(IN) :: q(3,nV)
    real(r8), intent(IN) :: time_offset(nV)
    real(r8), intent(IN) :: ascTAI(:), dscTAI(:)
    real(r8), intent(IN) :: dotVec(3,nV)
    real, intent(OUT) ::    geodAngle(nV)
    real, intent(in) ::     scOrbIncl(nV)
!    real(r8), intent(IN) :: scECR(3,nV)             ! s/c pos.
!    real(r8), intent(IN) :: scVelECR(3,nV)          ! s/c vel.

    ! Functions
    integer :: Pgs_csc_getEarthFigure

    ! Variables
    integer :: i, j, returnStatus, scOrb

    real(r8) :: a, asciiTAI, b, cSq, equatRad_a, orbRad, phiMin, polarRad_c
    real(r8) :: cosPhi(nV), gamma(nV), phi(nV), sinPhi(nV)
    real(r8) :: s(3,nV)
    real(r8) :: orbInclineNow

    ! Executable code

    ! Read a & b from earthfigure.dat
    returnStatus = Pgs_csc_getEarthFigure(earthModel, equatRad_a, polarRad_c)
    a = equatRad_a/1000
    b = polarRad_c/1000

    ! Set s = normalized dotVec
    do i = 1, nV

      orbInclineNow = scOrbIncl(i)
      if ( orbInclineNow == UNDEFINED_VALUE ) then
        call MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Error in calculating orbital inclination angle')
      endif

      ! Get a unit ECR vector to the point.
      s(:,i) = dotVec(:,i) / sqrt(dotVec(1,i)**2 + dotVec(2,i)**2 + &
        & dotVec(3,i)**2)
      
      ! Calculate the geocentric angle as a number of radians between 0 and PI
      gamma(i) = acos( q(1,i)*s(1,i) + q(2,i)*s(2,i) )
      
      ! Place angle between PI and 2*PI, if in Southern Hemisphere
      if ( dotVec(3,i) < 0.0 ) gamma(i) = 2*PI - gamma(i)

      ! Going to convert this to geodetic, calculate som parameters.
      orbRad= Deg2Rad * (orbInclineNow - 90)
      cSq = (1 + (tan(orbRad)**2)) * (a**2)*(b**2) / &
        & ( a**2 + (b**2)*(tan(orbRad)**2) )

      ! If |gamma| <= 45 deg of the equator, calculate phi using the SIN equation
      if ( (gamma(i) <= PI/4 ) .or. ( (gamma(i) >= 3*PI/4) .and. &
        & (gamma(i) <= 5*PI/4) ) .or. (gamma(i) >= 7*PI/4) ) then
        sinPhi(i) = sqrt( (a**4)*(sin(gamma(i))**2)/( (cSq**2)*&
          &(cos(gamma(i))**2) + (a**4)*(sin(gamma(i))**2) ) )
        phi(i) = asin(sinPhi(i))
      else
        ! If gamma is within 45 deg of a pole, calculate phi using the COS equation
        cosPhi(i) = sqrt( (cSq**2)*(cos(gamma(i))**2)/( (cSq**2)*&
          &cos(gamma(i))**2 + (a**4)*(sin(gamma(i))**2) ) )
        phi(i) = acos(cosPhi(i))
      endif

      ! Place phi in same quadrant as gamma
      if ( (gamma(i) > PI/2) .and. (gamma(i) <= PI) ) phi(i) = PI - phi(i)
      if ( (gamma(i) > PI) .and. (gamma(i) <= 3*PI/2) ) phi(i) = phi(i) + PI
      if ( gamma(i) > 3*PI/2 ) phi(i) = 2*PI - phi(i)

      ! Make phi cumulative over orbits
      ! First what is the time?
      asciiTAI = timeTAI + time_offset(i)
      scOrb = 0

      ! Now always make our calculation based on some threshold way back in
      ! time, to avoid any ambiguity
      if ( ( phi(i) < pi/2 ) .or. ( phi(i) >= 3*pi/2 ) ) then
        ! If in the ascending part of the orbit, base correction
        ! on time of last descending node
        ! MAKE THIS USE HUNT LATER!
        do j = 1, numOrb
          if ( asciiTAI > dscTAI(j) ) scOrb = j
        enddo
        ! If we're in the NH we've begun a new orbit, so add one
        if ( phi(i) < pi ) scOrb = scOrb + 1
      else
        ! If in descending part of the orbit, base it on the time
        ! of the last ascending node
        ! MAKE THIS USE HUNT LATER!
        do j = 1, numOrb
          if ( asciiTAI > ascTAI(j) ) scOrb = j
        enddo
      end if
      phi(i) = (scOrb-2)*2*PI + phi(i)

    enddo

    ! Convert to degrees for output
    geodAngle = Rad2Deg * phi

  end subroutine TkL1B_mc

end module TkL1B

! $Log$
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

