
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Time

  use MLSCommon, only: R8, TAI93_RANGE_T
  use MLSMessageModule, only: MLSMSG_ALLOCATE, MLSMSG_DeAllocate, &
    & MLSMSG_ERROR, MLSMSG_WARNING, MLSMESSAGE
  use SDPToolkit, only: PI, EARTHMODEL
  use OutputL1B, only: L1BOASC_T, LENCOORD
  use Scan, only: SCAN_GUESS
  use TkL1B, only: TKL1B_SC
  implicit none
  private

  public :: TIME_LAT, TIME_SEARCH, TIME_MAF, TIME_PRE, TIME_DAY, TIME_POST

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  character (LEN=*), parameter :: ModuleName= "$RCSfile$"
  !----------------------------------------------------------

  ! This module contains subroutines needed to calculate time
  ! information for the MLS scan programs.

  ! Parameters
  integer, parameter :: N_RECORDS = 148

contains

  !--------------------------------------------------------- Time_Lat ----
  subroutine Time_lat(asciiUTC, homeAlt, tpLat, converge)

    ! This subroutine takes an initial-guess view vector and iteratively runs the
    ! toolkit GrazingRay routine to find the actual view vector for a given
    ! tangent height.  The returned value is the latitude for the view vector.

    ! Arguments
    character(LEN=27), intent(IN) :: asciiUTC
    real(r8), intent(IN) :: homeAlt
    integer, intent(OUT) :: converge
    real(r8), intent(OUT) :: tpLat

    ! Parameters
    integer, parameter :: numValues = 1

    ! Functions
    integer :: Pgs_csc_geoToECR, Pgs_csc_grazingRay

    ! Variables
    type( L1BOAsc_T ) :: sc
    character (LEN=480) :: msr
    integer :: error, i, returnStatus
    real(r8) :: longitude, missAltitude, slantRange
    real(r8) :: ecr(3), initGuess(3), posNear(3), posSurf(3), ray(3)
    real(r8) :: time_offset(numValues)

    ! Executable code

    time_offset(numValues) = 0.0

    ! Get initial guess look-vector, initial s/c position, initial tp info
    call Scan_guess(asciiUTC, initGuess)

    allocate(sc%scECI(lenCoord,numValues), sc%scECR(lenCoord,numValues), &
      sc%scGeocAlt(numValues), sc%scGeocLat(numValues), &
      sc%scGeodAlt(numValues), sc%scGeodLat(numValues), sc%scLon(numValues), &
      sc%scGeodAngle(numValues), sc%scVelECI(lenCoord,numValues), &
      sc%scVelECR(lenCoord,numValues), sc%scOrbIncl(numValues), &
      sc%ypr(lenCoord,numValues), sc%yprRate(lenCoord,numValues), STAT=error)
    if ( error /= 0 ) then
      msr = MLSMSG_Allocate // '  s/c quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call TkL1B_sc(numValues, time_offset, asciiUTC, sc)

    returnStatus = Pgs_csc_grazingRay(earthModel, sc%scECR, initGuess, &
      tpLat, longitude, missAltitude, slantRange, posNear, posSurf)

    ! Iterations on view vector to see given tangent height
    converge = -1
    do i = 1, 8
      if ( abs(missAltitude - homeAlt) < 1.0 ) then
        converge = 1
        exit
      endif
      returnStatus = Pgs_csc_geoToECR(longitude, tpLat, homeAlt, &
        earthModel, ecr)
      ray = ecr - sc%scECR(:,1)
      returnStatus = Pgs_csc_grazingRay(earthModel, sc%scECR, ray, tpLat, &
        longitude, missAltitude, slantRange, posNear, posSurf)
    enddo

    ! Not sure why these weren't deallocated
    ! but better late than never (pw)
    deallocate( sc%scECI, sc%scECR, sc%scGeocAlt, sc%scGeocLat, &
      sc%scGeodAlt, sc%scGeodLat, sc%scLon, &
      sc%scGeodAngle, sc%scOrbIncl, sc%scVelECI, &
      sc%scVelECR, &
      sc%ypr, sc%yprRate, &
      & STAT=error )
    if ( error /= 0 ) then
      msr = MLSMSG_DeAllocate // ' Failed to deallocate s/c quantities.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

  end subroutine Time_lat

  !--------------------------------------------- Time_Search -------------
  subroutine Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
    finalTime, finalTAI, finalLat)
    ! This program calculates the MIF starting times for the orbits.

    ! Arguments
    character (LEN=27), intent(IN) :: estTime
    integer, intent(IN) :: scansPerOrb
    real(r8), intent(IN) :: homeAlt, homeLat, spm
    character (LEN=27), intent(OUT) :: finalTime
    real(r8), intent(OUT) :: finalLat, finalTAI

    ! Parameters
    character (LEN=*), parameter :: GR_ERR = 'Iterations timed out before &
      &height converged for '

    ! Functions
    integer :: Pgs_td_utcToTAI, Pgs_td_taiToUTC

    ! Variables
    character (LEN=27) :: timeF
    character (LEN=480) :: msg

    integer :: back, converge, i, returnStatus

    real(r8) :: estTAI93, latE, latF, tai

    ! Exectuable code

    finalTime = estTime

    ! Get tp info for estTime
    call Time_lat(estTime, homeAlt, latE, converge)
    if ( converge /= 1 ) then
      msg = GR_ERR // estTime
      call MLSMessage(MLSMSG_Warning, ModuleName, msg)
    endif

    returnStatus = Pgs_td_utcToTAI (estTime, estTAI93)

    ! Set search direction
    if (latE < homeLat) then
      back = -1
    else
      back = 1
    endif

    ! Check determined direction in MIF steps
    do i = 1, scansPerOrb*N_RECORDS
      tai = estTAI93 - back * spm
      returnStatus = Pgs_td_taiToUTC(tai, timeF)
      call Time_lat(timeF, homeAlt, latF, converge)
      if ( converge /= 1 ) then
        msg = GR_ERR // timeF
        call MLSMessage(MLSMSG_Warning, ModuleName, msg)
        exit
      endif
      if ( abs(latF - homeLat) >= abs(latE - homeLat) ) then
        exit
      else
        latE = latF
        estTAI93 = tai
        finalTime = timeF
      endif
    enddo

    ! Add an MIF, if closest lat found is below homeLat
    if (latE < homeLat) then
      tai = estTAI93 + spm
      returnStatus = Pgs_td_taiToUTC(tai, timeF)
      call Time_lat(timeF, homeAlt, latF, converge)
      if ( converge /= 1 ) then
        msg = GR_ERR // timeF
        call MLSMessage(MLSMSG_Warning, ModuleName, msg)
      endif
      latE = latF
      estTAI93 = tai
      finalTime = timeF
    endif
    finalLat = latE
    finalTAI = estTAI93

  end subroutine Time_search

  !-------------------------------------------------- Time_MAF ------------------
  subroutine Time_MAF(startTime, times, homeAlt, homeLat, mifRate, &
    scansPerOrb, spm, scanTimes, scanTAI, numMIFs, orbPerDay)
    ! This program calculates the starting times for the orbits and their scans,
    ! and the number of MIFs per MAF.

    ! Arguments
    type( TAI93_Range_T ), intent(IN) :: times
    character (LEN=27), intent(IN) :: startTime
    integer, intent(IN) :: mifRate, scansPerOrb
    real(r8), intent(IN) :: homeAlt, homeLat, spm
    character (LEN=27), intent(OUT) :: scanTimes(:,:)
    integer, intent(OUT) :: orbPerDay
    integer, intent(OUT) :: numMIFs(:,:)
    real(r8), intent(OUT) :: scanTAI(:,:)

    ! Parameters
    character (LEN=*), parameter ::GR_ERR = 'Iterations timed out before &
      &height converged for '
    real(r8), parameter :: minPerOrb = 98.9

    ! Functions
    integer :: Pgs_td_taiToUTC

    ! Variables
    character (LEN=27) :: estTime, timeF
    character (LEN=27), allocatable :: orbTime(:)
    character (LEN=480) :: msg, msr
    integer :: alloc_err, ascend, converge, dealloc_err, delMIF, delOrb, i, j
    integer :: k, opd, returnStatus
    real(r8) :: delAngle, delMAF, delTime, estLat, latF, orbRate, startLat
    real(r8) :: tai
    real(r8), allocatable :: orbTAI(:)

    ! Exectuable code

    orbRate = minPerOrb * 60 / (2*PI)
    delOrb = aint(minPerOrb * 60 * mifRate)
    opd = size(scanTimes,2)
    allocate( orbTime(opd), orbTAI(opd), STAT=alloc_err )
    if ( alloc_err /= 0 ) then
      msr = MLSMSG_Allocate // '  orbit variables.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    ! Get tp info for startTime
    call Time_lat(startTime, homeAlt, startLat, converge)
    if ( converge /= 1 ) then
      msg = GR_ERR // startTime
      call MLSMessage(MLSMSG_Warning, ModuleName, msg)
    endif

    ! Find whether ascending or descending at startTime
    ascend = -1
    tai = times%startTime + spm
    returnStatus = Pgs_td_taiToUTC(tai, timeF)
    call Time_lat(timeF, homeAlt, latF, converge)
    if ( converge /= 1 ) then
      msg = GR_ERR // timeF
      call MLSMessage(MLSMSG_Warning, ModuleName, msg)
    endif
    if ( (latF - startLat) > 0 ) ascend = 1

    ! Estimate MIF time of nearest ascending homeLat crossing
    if (ascend == 1) then
      if ( homeLat >= startLat ) then
        delAngle = homeLat - startLat
      else
        delAngle = 2*PI + homeLat - startLat
      endif
    else
      delAngle = PI + homeLat + startLat
    endif

    delTime = delAngle * orbRate
    delMIF = aint(delTime * mifRate)
    tai = times%startTime + delMIF * spm
    returnStatus = Pgs_td_taiToUTC(tai, estTime)

    ! Find MIF time of nearest ascending homeLat crossing
    orbTAI = 0.0

    call Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
      orbTime(1), orbTAI(1), estLat)

    ! Estimate time of remaining orbit crossings
    do i = 2, opd
      tai = orbTAI(i-1) + delOrb * spm
      if ( tai > times%endTime) then
        orbTAI(i) = 0.0
        orbPerDay = i-1
        exit
      endif
      returnStatus = Pgs_td_taiToUTC( tai, estTime)
      call Time_search(estTime, homeAlt, homeLat, scansPerOrb, spm, &
        orbTime(i), orbTAI(i), estLat)
      if ( orbTAI(i) > times%endTime) then
        orbTAI(i) = 0.0
        orbPerDay = i-1
        exit
      endif
    enddo

    ! Divide orbits into 240 MAFs
    scanTAI = 0.0
    do i = 1, orbPerDay
      if ( orbTAI(i+1) /= 0.0 ) then
        delMAF = ( orbTAI(i+1) - orbTAI(i) ) / scansPerOrb
        do j = 1, scansPerOrb
          delMIF = anint( (j-1) * delMAF * mifRate )
          scanTAI(j,i) = orbTAI(i) + delMIF * spm
          returnStatus = Pgs_td_taiToUTC( scanTAI(j,i), scanTimes(j,i) )
        enddo
      else
        scanTAI(1,i) = orbTAI(i)
        scanTimes(1,i) = orbTime(i)
        exit
      endif
    enddo

    ! Calculate number of MIFs
    do i = 1, orbPerDay-1
      do j = 1, scansPerOrb
        if ( j /= scansPerOrb) then
          do k = 1, N_RECORDS
            if ( scanTAI(j,i) + k*spm > scanTAI(j+1,i) ) exit
            numMIFs(j,i) = k
          enddo
        else 
          do k = 1, N_RECORDS
            if ( ( scanTAI(j,i) + k*spm ) > scanTAI(1,i+1) ) exit
            numMIFs(j,i) = k
          enddo
        endif
      enddo
    enddo

    deallocate( orbTime, orbTAI, STAT=dealloc_err )
    if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Warning, ModuleName, &
      'Failed deallocation of orbit variables.')

  end subroutine Time_MAF

  !-------------------------------------------- Time_Pre -----------------
  subroutine Time_pre(numValues, orbTime, startTAI, scansPerOrb, spm, &
    preTimes, preTAI, numMIFs, preMAF)

    ! This subroutine fills in the time gap for a partial orbit at the beginning of
    ! a file.  It duplicates the scan pattern for the first MAF of the first full
    ! orbit.

    ! Arguments
    character(LEN=27), intent(IN) :: orbTime
    integer, intent(IN) :: scansPerOrb
    integer, intent(IN) :: numValues(scansPerOrb)
    real(r8), intent(IN) :: spm, startTAI
    character(LEN=27), intent(OUT) :: preTimes(scansPerOrb)
    integer, intent(OUT) :: preMAF
    integer, intent(OUT) :: numMIFs(scansPerOrb)
    real(r8), intent(OUT) :: preTAI(scansPerOrb)

    ! Functions
    integer :: Pgs_td_utcToTAI, Pgs_td_taiToUTC

    ! Variables
    character(LEN=27) :: backTimes(scansPerOrb)
    integer :: i, returnStatus
    real(r8) :: mafTAI, orbTAI
    real(r8) :: backTAI(scansPerOrb)

    ! Executable code

    preMAF = 0

    ! Calculate MAF times prior to first full orbit
    returnStatus = Pgs_td_utcToTAI (orbTime, orbTAI)
    mafTAI = orbTAI
    do i = 1, scansPerOrb
      mafTAI = mafTAI - numValues(scansPerOrb - (i-1))*spm
      if ( mafTAI < startTAI ) exit
      returnStatus = Pgs_td_taiToUTC( mafTAI, backTimes(i) )
      backTAI(i) = mafTAI
      preMAF = i
    enddo

    ! Reverse numbering of MAF times, numValues
    do i = 1, preMAF
      preTimes(i) = backTimes( preMAF - (i-1) )
      preTAI(i)   = backTAI (preMAF - (i-1) )
      numMIFs(i) = numValues( scansPerOrb - (preMAF - i) )
    enddo

  end subroutine Time_pre

  !--------------------------------------------Time_Pos ------------------
  subroutine Time_post(numValues, orbTime, endTAI, scansPerOrb, spm, &
    postTimes, postTAI, postMAF)
    ! This subroutine fills in the time gap for a partial orbit at the end of a
    ! file.  It duplicates the scan pattern for the last MAF of the last full
    ! orbit.

    ! Arguments
    character(LEN=27), intent(IN) :: orbTime
    integer, intent(IN) :: scansPerOrb
    integer, intent(IN) :: numValues(scansPerOrb)
    real(r8), intent(IN) :: endTAI, spm
    character(LEN=27), intent(OUT) :: postTimes(scansPerOrb)
    integer, intent(OUT) :: postMAF
    real(r8), intent(OUT) :: postTAI(scansPerOrb)

    ! Functions
    integer :: Pgs_td_utcToTAI, Pgs_td_taiToUTC

    ! Variables
    integer :: i, returnStatus
    real(r8) :: orbTAI, mafTAI

    ! Executable code

    postMAF = 0

    ! Calculate MAF times after the last full orbit
    returnStatus = Pgs_td_utcToTAI (orbTime, orbTAI)
    mafTAI = orbTAI
    do i = 1, scansPerOrb
      mafTAI = mafTAI + numValues(i)*spm
      if ( mafTAI > endTAI ) exit
      returnStatus = Pgs_td_taiToUTC( mafTAI, postTimes(i+1) )
      postTAI(i+1) = mafTAI
      postMAF = i + 1
    enddo

    postTimes(1) = orbTime
    postTAI(1) = orbTAI

  end subroutine Time_post

  !------------------------------------ Time_Day -----------------------------
  subroutine Time_day(lenMAF, numMIFs, numValues, orbPerDay, postMAF, &
    postTAI, postTimes, preMAF, preTAI, preTimes, &
    scansPerOrb, scanTAI, scanTimes, mafTAI, mafTime, nV)
    ! This program reorganizes the various time arrays into a single array for the
    ! entire day/file.

    ! Arguments
    integer, intent(IN) :: lenMAF, postMAF, preMAF, orbPerDay, scansPerOrb
    character (LEN=27), intent(IN) :: postTimes(postMAF)
    character (LEN=27), intent(IN) :: preTimes(preMAF)
    character (LEN=27), intent(IN) :: scanTimes(scansPerOrb,orbPerDay)
    integer, intent(IN) :: numMIFs(preMAF)
    integer, intent(IN) :: numValues(scansPerOrb,orbPerDay)
    real(r8), intent(IN) :: postTAI(postMAF)
    real(r8), intent(IN) :: preTAI(preMAF)
    real(r8), intent(IN) :: scanTAI(scansPerOrb,orbPerDay)
    character (LEN=27), intent(OUT) :: mafTime(lenMAF)
    integer, intent(OUT) :: nV(lenMAF)
    real(r8), intent(OUT) :: mafTAI(lenMAF)

    ! Variables
    integer :: i, j, k

    ! Executable code

    ! MAFs prior to first full orbit
    do i = 1, preMAF
      if ( i > lenMAF ) exit
      nV(i) = numMIFs(i)
      mafTime(i) = preTimes(i)
      mafTAI(i) = preTAI(i)
    enddo

    ! MAFs for full orbits
    do i = 1, orbPerDay-1
      do j = 1, scansPerOrb

        k = j + (i-1)*scansPerOrb + preMAF
        if ( k > lenMAF ) exit

        nV(k) = numValues(j,i)
        mafTime(k) = scanTimes(j,i)
        mafTAI(k) = scanTAI(j,i)

      enddo
    enddo

    ! MAFs after last full orbit
    do i = 1, postMAF
      j = i + (orbPerDay-1)*scansPerOrb + preMAF
      if ( j > lenMAF ) exit

      nV(j) = numValues(i,orbPerDay-1)
      mafTime(j) = postTimes(i)
      mafTAI(j) = postTAI(i)
    enddo

  end subroutine Time_day

end module Time

! $Log$
! Revision 1.3  2001/12/04 00:25:25  pwagner
! sc%scvelECI is new name of scvel
!
! Revision 1.2  2001/10/11 23:27:23  livesey
! Tried to change the chmod stuff
!
! Revision 1.1  2000/11/30 16:31:13  nakamura
! Module for calculating time information for the MLS scan programs.
