
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE TkL1B
!===============================================================================

   USE MLSCommon
   USE MLSL1Common
   USE MLSMessageModule
   USE OutputL1B, ONLY: L1BOAsc_T, L1BOATP_T, L1BOAINDEX_T, lenCoord, &
        OutputL1B_THz, OutputL1B_SC, OutputL1B_Index, OutputL1B_GHz, lenG, lenT
   USE Scan
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- TkL1B_sc
!                TkL1B_tp
!                Tp_unit
!                L1boa_MAF
!                Mc_aux
!                TkL1B_mc

! Remarks:  This module contains subroutines for producing the L1BOA records on
!           a MAF by MAF basis.

CONTAINS

!-------------------------------------------------------
   SUBROUTINE TkL1B_sc(numValues, offsets, asciiUTC, sc)
!-------------------------------------------------------

! Brief description of subroutine
! This subroutine contains prototype code for creating the desired s/c record
! from the EPHEMATTIT output.

! Arguments

      TYPE( L1BOAsc_T ) :: sc

      CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

      INTEGER, INTENT(IN) :: numValues

      REAL(r8), INTENT(IN) :: offsets(numValues)

! Functions

      INTEGER :: Pgs_csc_eciToECR, Pgs_csc_ecrToGEO, Pgs_eph_ephemAttit

! Variables

      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: i, returnStatus
      INTEGER :: qualityFlags(2, numValues)

      REAL(r8) :: ecrVec(3)
      REAL(r8) :: radC(numValues), radD(numValues), radL(numValues)
      REAL(r8) :: attitQuat(4,numValues)
      REAL(r8) :: eciV(6,numValues), ecrV(6,numValues)

! Read oa data

      returnStatus = Pgs_eph_ephemAttit (spacecraftId, numValues, asciiUTC,  &
                             offsets, pgs_true, pgs_true, qualityFlags, &
                             sc%scECI, sc%scVel, sc%ypr, sc%yprRate, attitQuat)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = 'Routine ephemAttit, ' // mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Convert scECI to scECR

      eciV(1:3,:) = sc%scECI
      eciV(4:6,:) = 0.0

      returnStatus = Pgs_csc_eciToECR (numValues, asciiUTC, offsets, eciV, ecrV)

      sc%scECR = ecrV(1:3,:)
   
! Calculate geocentric/geodetic altitude, latitude & longitude from scECR

      DO i = 1, numValues

         sc%scGeocAlt(i) = SQRT( ecrV(1,i)**2 + ecrV(2,i)**2 + ecrV(3,i)**2 )

         radC(i) = ATAN( ecrV(3,i) / SQRT(ecrV(1,i)**2 + ecrV(2,i)**2) )
         sc%scGeocLat(i) = Rad2Deg * radC(i)

         ecrVec = ecrV(1:3, i)

         returnStatus = Pgs_csc_ecrToGEO (ecrVec, earthModel, radL(i), &
                                          radD(i), sc%scGeodAlt(i))

         sc%scGeodLat(i) = Rad2Deg * radD(i)
         sc%scLon(i) = Rad2Deg * radL(i)

      ENDDO

!-------------------------
   END SUBROUTINE TkL1B_sc
!-------------------------

!-----------------------------------------------------------------------------
   SUBROUTINE TkL1B_tp(asciiTAI, asciiUTC, lenG, numValues, offsets, posECR, &
                       scRate, startAngle, tp)
!-----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine fills the tangent point record.

! Arguments

      TYPE( L1BOAtp_T ) :: tp

      CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

      INTEGER, INTENT(IN) :: lenG, numValues

      REAL, INTENT(IN) :: scRate(lenG)

      REAL(r8), INTENT(IN) :: asciiTAI, startAngle
      REAL(r8), INTENT(IN) :: offsets(numValues)
      REAL(r8), INTENT(IN) :: posECR(3,lenG)

! Parameters

! Functions

      INTEGER :: Pgs_cbp_sat_cb_vector, Pgs_cbp_solarTimeCoords
      INTEGER :: Pgs_csc_scToOrb, Pgs_csc_scToECI, Pgs_csc_eciToECR
      INTEGER :: Pgs_csc_grazingRay, Pgs_csc_ecrToECI, Pgs_csc_eciToOrb
      INTEGER :: Pgs_td_taiToUTC

! Variables

      CHARACTER (LEN=27) :: time

      INTEGER :: error, flag, flagQ, i, returnStatus

      REAL(r8) :: declination, delAngle, delTime, deltaAlt, deltaLat, deltaLon
      REAL(r8) :: greenwich, localApparent, rightAscension, tai
      REAL(r8) :: dot(lenG), latD(lenG), localMean(lenG), lon(lenG), los(lenG)
      REAL(r8) :: sign(lenG), slantRange(lenG)
      REAL(r8) :: angleRad(numValues)
      REAL(r8) :: eci(3,lenG), ecr(3,lenG), hECR(3,lenG), nts(3,lenG)
      REAL(r8) :: posSurf(3,lenG), sc_frame_vector(3,lenG), sc_sun(3,lenG)
      REAL(r8) :: sc_tp(3,lenG), tp_sun(3,lenG), tpOrb(3,lenG), unitAlt(3,lenG)
      REAL(r8) :: unitLat(3,lenG), unitLon(3,lenG), vECR(3,lenG)
      REAL(r8) :: angleSc(3,numValues), angleOrb(3,numValues)
      REAL(r8) :: eciV(6,lenG), ecrV(6,lenG)

      deltaLat = 0.1
      deltaLon = 0.1
      deltaAlt = 1000.0

! Calculate MIF scan angles (in DEGREES)

      tp%scAngle(1) = startAngle

      DO i = 2, lenG
         tp%scAngle(i) = tp%scAngle(i-1) - scRate(i)*offsets(2)
      ENDDO

! Calculate retrace angle, rate

      delAngle = tp%scAngle(lenG) - tp%scAngle(1)

      delTime = (numValues - lenG)*offsets(2)

      tp%scanRate(1:lenG) = scRate
      tp%scanRate( (lenG+1):numValues ) = delAngle/delTime

      DO i = lenG+1, numvalues
         tp%scAngle(i) = tp%scAngle(i-1) - tp%scanRate(i)*offsets(2)
      ENDDO
      
! Put angle in s/c coordinates

      angleRad = Deg2Rad * tp%scAngle

      angleSc(1,:) = COS(angleRad)
      angleSc(2,:) = 0.0
      angleSc(3,:) = SIN(angleRad)

! Convert s/c vector to Orb vector/angle/degrees

      returnStatus = Pgs_csc_scToOrb(spacecraftId, numValues, asciiUTC, &
                                     offsets, angleSc, angleOrb)
   
      tp%scanAngle = Rad2Deg * ACOS( angleOrb(1,:) )

! Convert s/c vector to ECR

      returnStatus = Pgs_csc_scToECI(spacecraftId, lenG, asciiUTC, &
                                     offsets(1:lenG), angleSc(:,1:lenG), eci)
  
      eciV(1:3,:) = eci
      eciV(4:6,:) = 0.0

      returnStatus = Pgs_csc_eciToECR(lenG, asciiUTC, offsets(1:lenG), eciV, &
                                      ecrV)

      ecr = ecrV(1:3,:)

! For each scanning MIF,

      DO i = 1, lenG

! Calculate tangent point (geodetic & ECR)

         returnStatus = Pgs_csc_grazingRay(earthModel, posECR(:,i), ecr(:,i), &
                                latD(i), lon(i), tp%tpGeodAlt(i), &
                                slantRange(i), tp%tpECR(:,i), posSurf(:,i) )

! Create ECR unit vector quantities -- lat=1, lon=2, alt=3

         flagQ = 2
         CALL Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLon, &
                      unitLon(:,i), flag)

         flagQ = 1
         CALL Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaLat, &
                      unitLat(:,i), flag)

         flagQ = 3
         CALL Tp_unit(flagQ, lon(i), latD(i), tp%tpGeodAlt(i), deltaAlt, &
                      unitAlt(:,i), flag)

! Get local mean solar time from Toolkit

         tai = asciiTAI + (i-1)*offsets(2)

         returnStatus = Pgs_td_taiToUTC(tai, time)

         returnStatus = Pgs_cbp_solarTimeCoords(time, lon(i), greenwich, &
                      localMean(i), localApparent, rightAscension, declination)

      ENDDO

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
      DO i = 1, lenG
         sc_tp(:,i) = ecr(:,i) * slantRange(i)
      ENDDO
      tp_sun = sc_sun - sc_tp

      DO i = 1, lenG
         nts(:,i) = tp_sun(:,i) / SQRT( tp_sun(1,i)**2 + tp_sun(2,i)**2 + &
                                       &tp_sun(3,i)**2 )
      ENDDO

      dot = nts(1,:)*unitAlt(1,:) + nts(2,:)*unitAlt(2,:) + &
            nts(3,:)*unitAlt(3,:)

      tp%tpSolarZenith = ACOS(dot) * Rad2Deg

! Calculate losAngle

      DO i = 1, lenG
         vECR(:,i) = ecr(:,i) - ( ecr(1,i)*unitAlt(1,i) + &
                               &ecr(2,i)*unitAlt(2,i) + ecr(3,i)*unitAlt(3,i) )
         hECR(:,i) = vECR(:,i) / SQRT(vECR(1,i)**2 + vECR(2,i)**2 + &
                                     &vECR(3,i)**2)
      ENDDO

      los = ACOS( hECR(1,:)*unitLat(1,:) + hECR(2,:)*unitLat(2,:) + &
                 &hECR(3,:)*unitLat(3,:) )

      sign = ecr(1,:)*unitLon(1,:) + ecr(2,:)*unitLon(2,:) + &
             ecr(3,:)*unitLon(3,:)

      DO i = 1, lenG
         IF (sign(i) < 0 ) los(i) = 2*PI - los(i)
      ENDDO

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

      DO i = 1, lenG

         tp%tpGeocAlt(i) = SQRT( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2 + &
                                &tp%tpECR(3,i)**2 )

         tp%tpGeocLat(i) = Rad2Deg * ATAN( tp%tpECR(3,i) &
                          &/ SQRT( tp%tpECR(1,i)**2 + tp%tpECR(2,i)**2 ) )

      ENDDO

! Calculate dummy values

      tp%encoderAngle = tp%scAngle

      DO i = 2, lenG
         tp%tpGeocAltRate(i-1) = ( tp%tpGeocAlt(i) - tp%tpGeocAlt(i-1) )&
                                &/offsets(2)
         tp%tpGeodAltRate(i-1) = ( tp%tpGeodAlt(i) - tp%tpGeodAlt(i-1) )&
                                &/offsets(2)
      ENDDO
 
      tp%tpGeocAltRate(lenG) = tp%tpGeocAltRate(lenG-1)
      tp%tpGeodAltRate(lenG) = tp%tpGeodAltRate(lenG-1)

!-------------------------
   END SUBROUTINE TkL1B_tp
!-------------------------

!---------------------------------------------------------------
   SUBROUTINE Tp_unit (flagQ, lon, lat, alt, delta, unitQ, flag)
!---------------------------------------------------------------

      IMPLICIT NONE

! Brief description of subroutine
! This subroutine creates unit ECR vector quantities used by TkL1B_tp to
! calculate solarZenith and losAngle.

! Arguments

      INTEGER, INTENT(IN) :: flagQ

      REAL(r8), INTENT(IN) :: alt, delta, lat, lon
   
      INTEGER, INTENT(OUT) :: flag

      REAL(r8), INTENT(OUT) :: unitQ(3)

! Parameters

! Functions

      INTEGER :: Pgs_csc_geoToECR

! Variables

      CHARACTER (LEN=32) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: returnStatus

      REAL(r8) :: del, geoMinus, geoPlus
      REAL(r8) :: ecrMinus(3), ecrPlus(3), vec(3)

      flag = 0

! Create unit vector

      IF (flagQ == 2) THEN

! longitude

         del = delta * Deg2Rad

         geoPlus  = (lon + del/2)
         geoMinus = (lon - del/2)

         IF (geoPlus  >  PI) geoPlus  = geoPlus  - 2*PI
         IF (geoMinus < -PI) geoMinus = geoMinus + 2*PI
 
         returnStatus = Pgs_csc_geoToECR(geoPlus, lat, alt, earthModel, &
                                         ecrPlus)

         returnStatus = Pgs_csc_geoToECR(geoMinus, lat, alt, earthModel, &
                                         ecrMinus)

      ELSE IF (flagQ == 1) THEN

! latitude

         del = delta * Deg2Rad

         geoPlus  = (lat + del/2)
         geoMinus = (lat - del/2)

         returnStatus = Pgs_csc_geoToECR(lon, geoPlus, alt, earthModel, &
                                         ecrPlus)
  
         returnStatus = Pgs_csc_geoToECR(lon, geoMinus, alt, earthModel, &
                                         ecrMinus)

      ELSE

         geoPlus  = (alt + delta/2)
         geoMinus = (alt - delta/2)

         returnStatus = Pgs_csc_geoToECR(lon, lat, geoPlus, earthModel, &
                                         ecrPlus)
  
         returnStatus = Pgs_csc_geoToECR(lon, lat, geoMinus, earthModel, &
                                         ecrMinus)

      ENDIF

      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         flag = -1
      ENDIF

      vec = ecrPlus - ecrMinus

      unitQ = vec / SQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )

!------------------------
   END SUBROUTINE Tp_unit
!------------------------

!------------------------------------------------------------------------------
   SUBROUTINE L1boa_MAF(altG, altT, ascTAI, counterMAF, dscTAI, L1FileHandle, &
                        MAFinfo, noMAF, numOrb, orbIncline, orbitNumber, &
                        scRate, scRateT)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
! HDF output file.

! Arguments

      TYPE (MAFinfo_T) :: MAFinfo

      INTEGER, INTENT(IN) :: L1FileHandle, counterMAF, noMAF, numOrb
      INTEGER, INTENT(IN) :: orbitNumber(:)

      REAL, INTENT(IN) :: scRate(:)
      REAL, INTENT(IN) :: scRateT(:)

      REAL(r8), INTENT(IN) :: altG, altT, orbIncline
      REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:)

! Parameters

! Functions

      INTEGER :: Pgs_td_taiToUTC

! Variables

      TYPE( L1BOAindex_T ) :: index
      TYPE( L1BOAsc_T ) :: sc
      TYPE( L1BOAtp_T ) :: tp

      CHARACTER (LEN=27) :: mafTime
      CHARACTER (LEN=480) :: msr

      INTEGER :: error, i, nV, returnStatus

      REAL(r8) :: mafTAI
      REAL(r8) :: initRay(3), q(3)
      REAL(r8) :: offsets(MAFinfo%MIFsPerMAF)

! Calculate offsets array

      mafTAI = MAFinfo%startTAI
      nV = MAFinfo%MIFsPerMAF

      DO i = 1, nV 
         offsets(i) = (i-1)*MAFinfo%integTime
      ENDDO

! Write "index" information

      returnStatus = Pgs_td_taiToUTC(mafTAI, mafTime)
      
      index%MAFStartTimeUTC = mafTime
      index%MAFStartTimeTAI = mafTAI
      index%noMIFs = nV
      index%counterMAF = counterMAF

      CALL OutputL1B_index(noMAF, L1FileHandle, index)

! Allocate the MIF variables in the output structures

      ALLOCATE(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
       sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
       sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), &
       sc%scVel(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
       tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
       tp%scanRate(nV), STAT=error)
      IF ( error /= 0 ) THEN
         msr = MLSMSG_Allocate // '  s/c quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Get oa data

      CALL TkL1B_sc(nV, offsets, mafTime, sc)

! Get s/c master coordinate

      CALL Mc_aux(mafTime, sc%scECI(:,1), sc%scGeocLat(1), q)

      CALL TkL1B_mc(ascTAI, dscTAI, sc%scECI, sc%scGeocLat, nV, numOrb, &
                  orbIncline, orbitNumber, q, mafTAI, offsets, sc%scGeodAngle)

! Write s/c information

      CALL OutputL1B_sc(noMAF, L1FileHandle, sc)

! Calculate initial guess for look vector in ECR

      CALL Scan_guess(mafTime, initRay)

! Allocate the output structure

      ALLOCATE(tp%tpECI(lenCoord,lenG), tp%tpECR(lenCoord,lenG), &
               tp%tpOrbY(lenG), tp%tpGeocAlt(lenG), tp%tpGeocLat(lenG), &
               tp%tpGeocAltRate(lenG), tp%tpGeodAlt(lenG), &
               tp%tpGeodLat(lenG), tp%tpGeodAltRate(lenG), tp%tpLon(lenG), &
               tp%tpGeodAngle(lenG), tp%tpSolarTime(lenG), &
               tp%tpSolarZenith(lenG), tp%tpLosAngle(lenG), STAT=error)
      IF ( error /= 0 ) THEN
         msr = MLSMSG_Allocate // '  GHz tp quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Find angle, tan pt for start of GHZ scan

      CALL Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

! Calculate GHZ tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, lenG, nV, offsets, sc%scECR(:,1:lenG), &
                    scRate, tp%scAngle(1), tp)

! Compute GHz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, lenG, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:lenG), tp%tpGeodAngle)

! Write GHz information

      CALL OutputL1B_GHz(noMAF, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=error)
      IF ( error /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                             &deallocation of GHz quantities.')

! Find angle, tan pt for start of THz scan

      ALLOCATE(tp%tpECI(lenCoord,lenT), tp%tpECR(lenCoord,lenT), &
               tp%tpOrbY(lenT), tp%tpGeocAlt(lenT), tp%tpGeocLat(lenT), &
               tp%tpGeocAltRate(lenT), tp%tpGeodAlt(lenT), &
               tp%tpGeodLat(lenT), tp%tpGeodAltRate(lenT), tp%tpLon(lenT), &
               tp%tpGeodAngle(lenT), tp%tpSolarTime(lenT), &
               tp%tpSolarZenith(lenT), tp%tpLosAngle(lenT), STAT=error)
      IF ( error /= 0 ) THEN
         msr = MLSMSG_Allocate // '  THz tp quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

! Calculate THz tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, lenT, nV, offsets, sc%scECR(:,1:lenT), &
                    scRateT, tp%scAngle(1), tp)

! Compute THz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, lenT, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:lenT), tp%tpGeodAngle)

! Write THZ information

      CALL OutputL1B_THz(noMAF, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=error)
      IF ( error /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                         &deallocation of THz quantities.')

! Deallocate the MIF quantities

      DEALLOCATE(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
       sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle,  sc%scVel, &
       sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, tp%scanAngle, &
       tp%scanRate, STAT=error)
      IF ( error /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed &
                                             &deallocation of MIF quantities.')

!--------------------------
   END SUBROUTINE L1boa_MAF
!--------------------------

!---------------------------------------------------
   SUBROUTINE Mc_aux (asciiUTC, scECI, scGeocLat, q)
!---------------------------------------------------

! Brief description of subroutine
! This subroutine computes q, an auxilliary vector used in the calculation of
! the master coordinate.

! Arguments

      CHARACTER (LEN=27), INTENT(IN) :: asciiUTC

      REAL, INTENT(IN) :: scGeocLat

      REAL(r8), INTENT(IN) :: scECI(3)

      REAL(r8), INTENT(OUT) :: q(3)

! Parameters

      INTEGER, PARAMETER :: nV = 1

      REAL(r8), PARAMETER :: time_offset = 0.0

! Functions

      INTEGER :: Pgs_csc_orbToECI

! Variables

      INTEGER :: i, returnStatus

      REAL(r8) :: l
      REAL(r8) :: dist(2)
      REAL(r8) :: aECI(3), aux1(3), aux1ECI(3), aux2(3), aux2ECI(3)

! Construct two auxilliary vectors -- the first pointing directly ahead 
! along the s/c orbit, and the second pointing 45 degrees downward.

      aux1(1) = 1.0
      aux1(2) = 0.0
      aux1(3) = 0.0

      aux2(1) = 1.0
      aux2(2) = 0.0
      aux2(3) = 1.0

! Transform the auxilliary vectors from orbital to ECI coordinates

      returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
                                      time_offset, aux1, aux1ECI)

      returnStatus = Pgs_csc_orbToECI(spacecraftId, nV, asciiUTC, &
                                      time_offset, aux2, aux2ECI)

! Define l & a, where l is the distance to the equator along the direction of
! one of the auxilliary vectors from the s/c, and a is the auxilliary vector
! for which l was defined.

      IF ( aux1ECI(3) == 0 ) THEN

! If aux1 has a zero z-component, set l & a for aux2, if its z-component /= 0

         IF ( aux2ECI(3) /= 0 ) THEN
            l = -scECI(3)/aux2ECI(3)
            aECI = aux2ECI
         ENDIF

      ELSE

! If aux1(3) is not zero, but aux2(3) is, then set l & a for aux1

         IF ( aux2ECI(3) == 0 ) THEN

            l = -scECI(3)/aux1ECI(3)
            aECI = aux1ECI

         ELSE

! If both aux1(3) and aux2(3) are non-zero, choose the one which gives the
! minimum absolute value of l.

            dist(1) = -scECI(3)/aux1ECI(3)
            dist(2) = -scECI(3)/aux2ECI(3)

            IF ( ABS(dist(1)) < ABS(dist(2)) ) THEN
               l = dist(1)
               aECI = aux1ECI
            ELSE
               l = dist(2)
               aECI = aux2ECI
            ENDIF

         ENDIF

      ENDIF
   
! Define the vector q = scECI + la, such that its z-component = 0

      q(1:2) = scECI(1:2) + l*aECI(1:2)
      q(3) = 0.0

! Modify q, depending on which hemisphere the s/c is in, and the sign of l

      IF ( scGeocLat < 0.0 ) THEN
         IF ( l < 0 ) q = -q
      ELSE
         IF ( l > 0 ) q = -q
      ENDIF

! Normalize q

      q = q / SQRT( q(1)**2 + q(2)**2 )
 
!-----------------------
   END SUBROUTINE Mc_aux
!-----------------------

!------------------------------------------------------------------------------
   SUBROUTINE TkL1B_mc(ascTAI, dscTAI, dotVec, geocLat, nV, numOrb, &
                   orbIncline, orbitNumber, q, timeTAI, time_offset, geodAngle)
!------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine computes phi, the master coordinate for the spacecraft and
! tangent point records.

! Arguments

      INTEGER, INTENT(IN) :: nV, numOrb
      INTEGER, INTENT(IN) :: orbitNumber(:)

      REAL, INTENT(IN) :: geocLat(nV)

      REAL(r8), INTENT(IN) :: orbIncline, timeTAI
      REAL(r8), INTENT(IN) :: q(3)
      REAL(r8), INTENT(IN) :: time_offset(nV)
      REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:)
      REAL(r8), INTENT(IN) :: dotVec(3,nV)

      REAL, INTENT(OUT) :: geodAngle(nV)

! Parameters

! Functions

      INTEGER :: Pgs_csc_getEarthFigure

! Variables

      INTEGER :: i, j, returnStatus, scOrb

      REAL(r8) :: a, asciiTAI, b, cSq, equatRad_a, orbRad, phiMin, polarRad_c
      REAL(r8) :: cosPhi(nV), gamma(nV), phi(nV), sinPhi(nV)
      REAL(r8) :: s(3,nV)

! Read a & b from earthfigure.dat

      returnStatus = Pgs_csc_getEarthFigure(earthModel, equatRad_a, polarRad_c)

      a = equatRad_a/1000
      b = polarRad_c/1000

! Calculate C-squared

      orbRad= Deg2Rad * (orbIncline - 90)

      cSq = (1 + (TAN(orbRad)**2)) * (a**2)*(b**2)/(a**2 + &
                                                  &(b**2)*(TAN(orbRad)**2))

! Set s = normalized dotVec

      DO i = 1, nV

         s(:,i) = dotVec(:,i) / SQRT(dotVec(1,i)**2 + dotVec(2,i)**2 + &
                                                     &dotVec(3,i)**2)

! Calculate the geocentric angle as a number of radians between 0 and PI

         gamma(i) = ACOS( q(1)*s(1,i) + q(2)*s(2,i) )

! Place angle between PI and 2*PI, if lat indicates Southern Hemisphere

         IF ( geocLat(i) < 0.0 ) gamma(i) = 2*PI - gamma(i)

! If |gamma| <= 45 deg of the equator, calculate phi using the SIN equation

         IF ( (gamma(i) <= PI/4 ) .OR. ( (gamma(i) >= 3*PI/4) .AND. &
              &(gamma(i) <= 5*PI/4) ) .OR. (gamma(i) >= 7*PI/4) ) THEN

            sinPhi(i) = SQRT( (a**4)*(SIN(gamma(i))**2)/( (cSq**2)*&
                             &(COS(gamma(i))**2) + (a**4)*(SIN(gamma(i))**2) ) )

            phi(i) = ASIN(sinPhi(i))

         ELSE

! If gamma is within 45 deg of a pole, calculate phi using the COS equation

            cosPhi(i) = SQRT( (cSq**2)*(COS(gamma(i))**2)/( (cSq**2)*&
                              &COS(gamma(i))**2 + (a**4)*(SIN(gamma(i))**2) ) )

            phi(i) = ACOS(cosPhi(i))

         ENDIF

! Place phi in same quadrant as gamma

         IF ( (gamma(i) > PI/2) .AND. (gamma(i) < PI) ) phi(i) = PI - phi(i)
         IF ( (gamma(i) > PI) .AND. (gamma(i) < 3*PI/2) ) phi(i) = phi(i) + PI
         IF ( gamma(i) > 3*PI/2 ) phi(i) = 2*PI - phi(i)

! Make phi cumulative over orbits

         asciiTAI = timeTAI + time_offset(i)

         DO j = 1, numOrb
            IF ( asciiTAI > ascTAI(j) ) scOrb = orbitNumber(j)
         ENDDO

         phi(i) = (scOrb-1)*2*PI + phi(i)

         IF ( asciiTAI > dscTAI(scOrb+1) ) THEN
            phiMin = scOrb*2*PI - PI
         ELSE
            phiMin = (scOrb-1)*2*PI - PI
         ENDIF

         IF ( phi(i) < phiMin ) phi(i) = phi(i) + 2*PI

      ENDDO

! Convert to degrees for output

      geodAngle = Rad2Deg * phi

!----------------------
END SUBROUTINE TkL1B_mc
!----------------------

!===============
END MODULE TkL1B
!===============

! $Log$
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
! Revision 2.0  2000/09/05 18:55:15  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/15 18:48:37  nakamura
! Absorbed module Mc; moved _init subroutine to Orbit; account for parametrization of lenG & lenT.
!

