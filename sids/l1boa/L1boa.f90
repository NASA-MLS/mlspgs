! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L1boa

  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR
  use OutputL1B, only: OUTPUTL1B_INDEX, L1BOAINDEX_T, L1BOASC_T, L1BOATP_T, LENCOORD, &
    & OUTPUTL1B_GHZ, OUTPUTL1B_THZ, OUTPUTL1B_SC
  use Scan, only: SCAN_GUESS, SCAN_START
  use Sd, only: SD_SC, SD_GHZ, SD_THZ, SD_INDEX
  use TkL1B, only: TKL1B_MC, TKL1B_TP, TKL1B_SC, MC_AUX
  implicit none
  private

  public :: L1BOA_NOFILL, L1BOA_FILL

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &                                                    
    "$Id$"
  character (LEN=*), parameter :: ModuleName="$RCSfile$"
  !----------------------------------------------------------

contains

  !---------------------------------------------- L1BOA_NOFILL -----
  subroutine L1boa_nofill(altG, altT, ascTAI, dscTAI, L1FileHandle, mifG, &
    mifT, index, noMAF, numOrb, offsets, orbIncline, &
    orbitNumber, scRate, scRateT)
    ! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
    ! HDF output file.

    ! Arguments
    type( L1BOAindex_T), intent(IN) :: index
    integer, intent(IN) :: L1FileHandle, mifG, mifT, noMAF, numOrb
    integer, intent(IN) :: orbitNumber(:)
    real, intent(IN) :: scRate(mifG)
    real, intent(IN) :: scRateT(mifT)
    real(r8), intent(IN) :: altG, altT, orbIncline
    real(r8), intent(IN) :: ascTAI(:), dscTAI(:), offsets(:)

    ! Variables
    type( L1BOAsc_T ) :: sc
    type( L1BOAtp_T ) :: tp
    character (LEN=27) :: mafTime
    character (LEN=480) :: msr
    integer :: alloc_err, dealloc_err, nV
    real(r8) :: mafTAI
    real(r8) :: initRay(3), q(3)

    ! Executable code

    ! Write "index" information
    call OutputL1B_index(noMAF, L1FileHandle, index)

    mafTime = index%MAFStartTimeUTC
    mafTAI = index%MAFStartTimeTAI
    nV = index%noMIFs

    ! Get oa data
    allocate(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
      sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
      sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), &
      sc%scVel(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
      tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
      tp%scanRate(nV), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  s/c quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call TkL1B_sc(nV, offsets(1:nV), mafTime, sc)

    ! Get s/c master coordinate
    call Mc_aux(mafTime, sc%scECI(:,1), sc%scGeocLat(1), q)

    call TkL1B_mc(ascTAI, dscTAI, sc%scECI, sc%scGeocLat, nV, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:nV), sc%scGeodAngle)

    ! Write s/c information
    call OutputL1B_sc(noMAF, L1FileHandle, sc)

    ! Calculate initial guess for look vector in ECR
    call Scan_guess(mafTime, initRay)

    ! Find angle, tan pt for start of GHZ scan
    allocate(tp%tpECI(lenCoord,mifG), tp%tpECR(lenCoord,mifG), &
      tp%tpOrbY(mifG), tp%tpGeocAlt(mifG), tp%tpGeocLat(mifG), &
      tp%tpGeocAltRate(mifG), tp%tpGeodAlt(mifG), &
      tp%tpGeodLat(mifG), tp%tpGeodAltRate(mifG), tp%tpLon(mifG), &
      tp%tpGeodAngle(mifG), tp%tpSolarTime(mifG), &
      tp%tpSolarZenith(mifG), tp%tpLosAngle(mifG), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  GHz quantites.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

    ! Calculate GHZ tan pt record
    call TkL1B_tp(mafTAI, mafTime, mifG, nV, offsets(1:nV), &
      sc%scECR(:,1:mifG), scRate, tp%scAngle(1), tp)

    ! Compute GHz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifG, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:mifG), tp%tpGeodAngle)

    ! Write GHz information
    call OutputL1B_GHz(noMAF, L1FileHandle, tp)

    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of GHz quantities.')

    ! Find angle, tan pt for start of THz scan
    allocate(tp%tpECI(lenCoord,mifT), tp%tpECR(lenCoord,mifT), &
      tp%tpOrbY(mifT), tp%tpGeocAlt(mifT), tp%tpGeocLat(mifT), &
      tp%tpGeocAltRate(mifT), tp%tpGeodAlt(mifT), &
      tp%tpGeodLat(mifT), tp%tpGeodAltRate(mifT), tp%tpLon(mifT), &
      tp%tpGeodAngle(mifT), tp%tpSolarTime(mifT), &
      tp%tpSolarZenith(mifT), tp%tpLosAngle(mifT), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  THz quantites.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

    ! Calculate THz tan pt record
    call TkL1B_tp(mafTAI, mafTime, mifT, nV, offsets(1:nV), &
      sc%scECR(:,1:mifT), scRateT, tp%scAngle(1), tp)

    ! Compute THz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifT, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:mifT), tp%tpGeodAngle)

    ! Write THZ information
    call OutputL1B_THz(noMAF, L1FileHandle, tp)

    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of THz quantites.')

    ! Deallocate the MIF quantities
    deallocate(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
      sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle,  sc%scVel, &
      sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, tp%scanAngle, &
      tp%scanRate, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of MIF quantities.')

  end subroutine L1boa_nofill

  !--------------------------------------------- L1BOA_FILL -----
  subroutine L1boa_fill(altG, altT, ascTAI, dscTAI, L1FileHandle, mifG, &
    mifT, index, noMAF, numOrb, offsets, orbIncline, &
    orbitNumber, scRate, scRateT)
    ! This subroutine creates the SIDS L1BOA MAF records, and writes them to an
    ! HDF output file.

    ! Arguments
    type( L1BOAindex_T), intent(IN) :: index
    integer, intent(IN) :: L1FileHandle, mifG, mifT, noMAF, numOrb
    integer, intent(IN) :: orbitNumber(:)
    real, intent(IN) :: scRate(mifG)
    real, intent(IN) :: scRateT(mifT)
    real(r8), intent(IN) :: altG, altT, orbIncline
    real(r8), intent(IN) :: ascTAI(:), dscTAI(:), offsets(:)

    ! Variables
    type( L1BOAsc_T ) :: sc
    type( L1BOAtp_T ) :: tp
    character(LEN=27) :: mafTime
    character (LEN=480) :: msr
    integer :: alloc_err, dealloc_err, nV
    real(r8) :: mafTAI
    real(r8) :: initRay(3), q(3)

    ! Executable code

    call Sd_index(noMAF, L1FileHandle, index)
    mafTime = index%MAFStartTimeUTC
    mafTAI = index%MAFStartTimeTAI
    nV = index%noMIFs

    ! Get oa data
    allocate(sc%scECI(lenCoord,nV), sc%scECR(lenCoord,nV), &
      sc%scGeocAlt(nV), sc%scGeocLat(nV), sc%scGeodAlt(nV), &
      sc%scGeodLat(nV), sc%scLon(nV), sc%scGeodAngle(nV), &
      sc%scVel(lenCoord,nV), sc%ypr(lenCoord,nV), sc%yprRate(lenCoord,nV), &
      tp%encoderAngle(nV), tp%scAngle(nV), tp%scanAngle(nV), &
      tp%scanRate(nV), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  s/c quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call TkL1B_sc(nV, offsets(1:nV), mafTime, sc)

    ! Get s/c master coordinate
    call Mc_aux(mafTime, sc%scECI(:,1), sc%scGeocLat(1), q)
    call TkL1B_mc(ascTAI, dscTAI, sc%scECI, sc%scGeocLat, nV, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:nV), sc%scGeodAngle)

    ! Write s/c information
    call Sd_sc(noMAF, nV, L1FileHandle, sc)

    ! Calculate initial guess for look vector in ECR
    call Scan_guess(mafTime, initRay)

    ! Find angle, tan pt for start of GHZ scan
    allocate(tp%tpECI(lenCoord,mifG), tp%tpECR(lenCoord,mifG), &
      tp%tpOrbY(mifG), tp%tpGeocAlt(mifG), tp%tpGeocLat(mifG), &
      tp%tpGeocAltRate(mifG), tp%tpGeodAlt(mifG), &
      tp%tpGeodLat(mifG), tp%tpGeodAltRate(mifG), tp%tpLon(mifG), &
      tp%tpGeodAngle(mifG), tp%tpSolarTime(mifG), &
      tp%tpSolarZenith(mifG), tp%tpLosAngle(mifG), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  GHz quantites.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call Scan_start( altG, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1) )

    ! Calculate GHZ tan pt record
    call TkL1B_tp(mafTAI, mafTime, mifG, nV, offsets(1:nV), &
      sc%scECR(:,1:mifG), scRate, tp%scAngle(1), tp)

    ! Compute GHz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifG, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:mifG), tp%tpGeodAngle)

    ! Write GHz information

    call Sd_GHz(mifG, noMAF, nV, L1FileHandle, tp)
    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of GHz quantities.')

    ! Find angle, tan pt for start of THz scan
    allocate(tp%tpECI(lenCoord,mifT), tp%tpECR(lenCoord,mifT), &
      tp%tpOrbY(mifT), tp%tpGeocAlt(mifT), tp%tpGeocLat(mifT), &
      tp%tpGeocAltRate(mifT), tp%tpGeodAlt(mifT), &
      tp%tpGeodLat(mifT), tp%tpGeodAltRate(mifT), tp%tpLon(mifT), &
      tp%tpGeodAngle(mifT), tp%tpSolarTime(mifT), &
      tp%tpSolarZenith(mifT), tp%tpLosAngle(mifT), STAT=alloc_err)
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  THz quantites.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call Scan_start(altT, sc%scECR(:,1), mafTime, initRay, tp%scAngle(1))

    ! Calculate THz tan pt record
    call TkL1B_tp(mafTAI, mafTime, mifT, nV, offsets(1:nV), &
      sc%scECR(:,1:mifT), scRateT, tp%scAngle(1), tp)

    ! Compute THz master coordinate
    call TkL1B_mc(ascTAI, dscTAI, tp%tpECI, tp%tpGeocLat, mifT, numOrb, &
      orbIncline, orbitNumber, q, mafTAI, offsets(1:mifT), tp%tpGeodAngle)

    ! Write THZ information
    call Sd_THz(mifT, noMAF, nV, L1FileHandle, tp)

    deallocate(tp%tpECI, tp%tpECR, tp%tpOrbY, tp%tpGeocAlt, tp%tpGeocLat, &
      tp%tpGeocAltRate, tp%tpGeodAlt, tp%tpGeodLat, &
      tp%tpGeodAltRate, tp%tpLon, tp%tpGeodAngle, tp%tpSolarTime, &
      tp%tpSolarZenith, tp%tpLosAngle, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of THz quantites.')

    ! Deallocate the MIF quantities

    deallocate(sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
      sc%scGeodAlt, sc%scGeodLat, sc%scLon, sc%scGeodAngle,  sc%scVel, &
      sc%ypr, sc%yprRate, tp%encoderAngle, tp%scAngle, tp%scanAngle, &
      tp%scanRate, STAT=dealloc_err)
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      'Failed deallocation of MIF quantities.')

  end subroutine L1boa_fill

end module L1boa

! $Log$
! Revision 1.3  2001/10/11 23:27:23  livesey
! Tried to change the chmod stuff
!
! Revision 1.2  2001/10/11 23:25:29  livesey
! Tidied up a bit
!
! Revision 1.1  2000/11/30 16:26:56  nakamura
! Module for writing SIDS L1BOA data without re-setting fill mode/values.
!
!
