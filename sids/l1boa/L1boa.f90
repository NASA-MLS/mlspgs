
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE L1boa
!===============================================================================

   USE MLSCommon
   USE MLSMessageModule
   USE OutputL1B
   USE Scan
   USE Sd
   USE TkL1B
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- L1boa_nofill
!                L1boa_fill

! Remarks:  This module contains the subroutines needed to write the SIDS L1BOA
!           data to SD arrays without re-setting the fill mode/values.

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE L1boa_nofill(altG, altT, ascTAI, dscTAI, L1FileHandle, mifG, &
                           mifT, index, noMAF, numOrb, offsets, orbIncline, &
                           orbitNumber, scRate, scRateT)
!----------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
! HDF output file.

! Arguments

      TYPE( L1BOAindex_T), INTENT(IN) :: index

      INTEGER, INTENT(IN) :: L1FileHandle, mifG, mifT, noMAF, numOrb
      INTEGER, INTENT(IN) :: orbitNumber(:)

      REAL, INTENT(IN) :: scRate(mifG)
      REAL, INTENT(IN) :: scRateT(mifT)

      REAL(r8), INTENT(IN) :: altG, altT, orbIncline
      REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:), offsets(:)

! Parameters
 
! Functions

! Variables
 
      TYPE( L1BOAsc_T ) :: sc
      TYPE( L1BOAtp_T ) :: tp

      CHARACTER (LEN=27) :: mafTime
      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, dealloc_err, nV

      REAL(r8) :: mafTAI
      REAL(r8) :: initRay(3), q(3)

! Write "index" information

      CALL OutputL1B_index(noMAF, L1FileHandle, index)

      mafTime = index%MAFStartTimeUTC
      mafTAI = index%MAFStartTimeTAI
      nV = index%noMIFs

! Get oa data

      ALLOCATE(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
       sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
       sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), &
       sc%scVel(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
       tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
       tp%scanRate(nV), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  s/c quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL TkL1B_sc(nV, offsets(1:nV), mafTime, sc)

! Get s/c master coordinate

      CALL Mc_aux(mafTime, sc%scECI(:,1), sc%scGeocLat(1), q)

      CALL TkL1B_mc(ascTAI, dscTAI, sc%scECI, sc%scGeocLat, nV, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:nV), sc%scGeodAngle)

! Write s/c information

      CALL OutputL1B_sc(noMAF, L1FileHandle, sc)

! Calculate initial guess for look vector in ECR

      CALL Scan_guess(mafTime, initRay)

! Find angle, tan pt for start of GHZ scan

      ALLOCATE(tp%tpECI(lenCoord,mifG), tp%tpECR(lenCoord,mifG), &
               tp%tpOrbY(mifG), tp%tpGeocAlt(mifG), tp%tpGeocLat(mifG), &
               tp%tpGeocAltRate(mifG), tp%tpGeodAlt(mifG), &
               tp%tpGeodLat(mifG), tp%tpGeodAltRate(mifG), tp%tpLon(mifG), &
               tp%tpGeodAngle(mifG), tp%tpSolarTime(mifG), &
               tp%tpSolarZenith(mifG), tp%tpLosAngle(mifG), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
          msr = MLSMSG_Allocate // '  GHz quantites.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

! Calculate GHZ tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, mifG, nV, offsets(1:nV), &
                    sc%scECR(:,1:mifG), scRate, tp%scAngle(1), tp)

! Compute GHz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifG, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:mifG), tp%tpGeodAngle)

! Write GHz information

      CALL OutputL1B_GHz(noMAF, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of GHz quantities.')

! Find angle, tan pt for start of THz scan

      ALLOCATE(tp%tpECI(lenCoord,mifT), tp%tpECR(lenCoord,mifT), &
               tp%tpOrbY(mifT), tp%tpGeocAlt(mifT), tp%tpGeocLat(mifT), &
               tp%tpGeocAltRate(mifT), tp%tpGeodAlt(mifT), &
               tp%tpGeodLat(mifT), tp%tpGeodAltRate(mifT), tp%tpLon(mifT), &
               tp%tpGeodAngle(mifT), tp%tpSolarTime(mifT), &
               tp%tpSolarZenith(mifT), tp%tpLosAngle(mifT), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  THz quantites.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

! Calculate THz tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, mifT, nV, offsets(1:nV), &
                    sc%scECR(:,1:mifT), scRateT, tp%scAngle(1), tp)

! Compute THz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifT, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:mifT), tp%tpGeodAngle)

! Write THZ information
                     
      CALL OutputL1B_THz(noMAF, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of THz quantites.')

! Deallocate the MIF quantities

      DEALLOCATE(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
       sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle,  sc%scVel, &
       sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, tp%scanAngle, &
       tp%scanRate, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of MIF quantities.')

!-----------------------------
   END SUBROUTINE L1boa_nofill
!-----------------------------

!--------------------------------------------------------------------------
   SUBROUTINE L1boa_fill(altG, altT, ascTAI, dscTAI, L1FileHandle, mifG, &
                         mifT, index, noMAF, numOrb, offsets, orbIncline, &
                         orbitNumber, scRate, scRateT)
!--------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
! HDF output file.

! Arguments

      TYPE( L1BOAindex_T), INTENT(IN) :: index

      INTEGER, INTENT(IN) :: L1FileHandle, mifG, mifT, noMAF, numOrb
      INTEGER, INTENT(IN) :: orbitNumber(:)

      REAL, INTENT(IN) :: scRate(mifG)
      REAL, INTENT(IN) :: scRateT(mifT)

      REAL(r8), INTENT(IN) :: altG, altT, orbIncline
      REAL(r8), INTENT(IN) :: ascTAI(:), dscTAI(:), offsets(:)

! Parameters
 
! Functions

! Variables

      TYPE( L1BOAsc_T ) :: sc
      TYPE( L1BOAtp_T ) :: tp

      CHARACTER(LEN=27) :: mafTime
      CHARACTER (LEN=480) :: msr

      INTEGER :: alloc_err, dealloc_err, nV

      REAL(r8) :: mafTAI
      REAL(r8) :: initRay(3), q(3)

! Write "index" information

      CALL Sd_index(noMAF, L1FileHandle, index)
      
      mafTime = index%MAFStartTimeUTC
      mafTAI = index%MAFStartTimeTAI
      nV = index%noMIFs

! Get oa data

      ALLOCATE(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
       sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
       sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), &
       sc%scVel(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
       tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
       tp%scanRate(nV), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  s/c quantities.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL TkL1B_sc(nV, offsets(1:nV), mafTime, sc)

! Get s/c master coordinate

      CALL Mc_aux(mafTime, sc%scECI(:,1), sc%scGeocLat(1), q)

      CALL TkL1B_mc(ascTAI, dscTAI, sc%scECI, sc%scGeocLat, nV, numOrb, &
            orbIncline, orbitNumber, q, mafTAI, offsets(1:nV), sc%scGeodAngle)

! Write s/c information

      CALL Sd_sc(noMAF, nV, L1FileHandle, sc)

! Calculate initial guess for look vector in ECR

      CALL Scan_guess(mafTime, initRay)

! Find angle, tan pt for start of GHZ scan

      ALLOCATE(tp%tpECI(lenCoord,mifG), tp%tpECR(lenCoord,mifG), &
               tp%tpOrbY(mifG), tp%tpGeocAlt(mifG), tp%tpGeocLat(mifG), &
               tp%tpGeocAltRate(mifG), tp%tpGeodAlt(mifG), &
               tp%tpGeodLat(mifG), tp%tpGeodAltRate(mifG), tp%tpLon(mifG), &
               tp%tpGeodAngle(mifG), tp%tpSolarTime(mifG), &
               tp%tpSolarZenith(mifG), tp%tpLosAngle(mifG), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
          msr = MLSMSG_Allocate // '  GHz quantites.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

! Calculate GHZ tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, mifG, nV, offsets(1:nV), &
                    sc%scECR(:,1:mifG), scRate, tp%scAngle(1), tp)

! Compute GHz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifG, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:mifG), tp%tpGeodAngle)

! Write GHz information

      CALL Sd_GHz(mifG, noMAF, nV, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of GHz quantities.')

! Find angle, tan pt for start of THz scan

      ALLOCATE(tp%tpECI(lenCoord,mifT), tp%tpECR(lenCoord,mifT), &
               tp%tpOrbY(mifT), tp%tpGeocAlt(mifT), tp%tpGeocLat(mifT), &
               tp%tpGeocAltRate(mifT), tp%tpGeodAlt(mifT), &
               tp%tpGeodLat(mifT), tp%tpGeodAltRate(mifT), tp%tpLon(mifT), &
               tp%tpGeodAngle(mifT), tp%tpSolarTime(mifT), &
               tp%tpSolarZenith(mifT), tp%tpLosAngle(mifT), STAT=alloc_err)
      IF ( alloc_err /= 0 ) THEN
         msr = MLSMSG_Allocate // '  THz quantites.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

      CALL Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

! Calculate THz tan pt record

      CALL TkL1B_tp(mafTAI, mafTime, mifT, nV, offsets(1:nV), &
                    sc%scECR(:,1:mifT), scRateT, tp%scAngle(1), tp)

! Compute THz master coordinate

      CALL TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifT, numOrb, &
           orbIncline, orbitNumber, q, mafTAI, offsets(1:mifT), tp%tpGeodAngle)

! Write THZ information
                     
      CALL Sd_THz(mifT, noMAF, nV, L1FileHandle, tp)

      DEALLOCATE(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
                 tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
                 tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
                 tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of THz quantites.')

! Deallocate the MIF quantities

      DEALLOCATE(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
       sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle,  sc%scVel, &
       sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, tp%scanAngle, &
       tp%scanRate, STAT=dealloc_err)
      IF ( dealloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
                                      'Failed deallocation of MIF quantities.')

!---------------------------
   END SUBROUTINE L1boa_fill
!---------------------------

!===============
END MODULE L1boa
!===============

!# $Log$
!#
