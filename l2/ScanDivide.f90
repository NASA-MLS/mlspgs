
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================
module ScanDivide
!=======================================================================

  use EXPR_M, only: EXPR   
  use INIT_TABLES_MODULE, only:   FIRST_PARM, L_BOTH, L_EITHER, L_GHZ, &
    L_NONE, L_THZ, LAST_PARM, P_CRITICAL_BANDS, P_CRITICAL_SCANNING_MODULES, &
    P_HOME_GEOD_ANGLE, P_HOME_MODULE, P_IDEAL_LENGTH, P_MAX_GAP, P_OVERLAP, &
    P_SCAN_LOWER_LIMIT, P_SCAN_UPPER_LIMIT, PARM_INDICES, &
    PHYQ_INVALID, PHYQ_LENGTH, PHYQ_MAFS, PHYQ_TIME
  use L1BData, only: deallocateL1BDATA, L1BDATA_T, NAME_LEN, READL1BDATA
  use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use MLSNumerics, only: HUNT, R8
  use MLSStrings, only: MLSMSG_L1BRead
  use SDPToolkit, only: MAX_ORBITS
  use STRING_TABLE, only: GET_STRING, STRING_LENGTH
  use TREE, only: DECORATION, NSONS, SUBTREE
  implicit none
  private

  public :: DestroyChunkDatabase, ScanAndDivide

  private :: ID, ModuleName

!------------------- RCS Ident Info ------------------------------------
  character(len=130) :: Id = &
  "$Id$"
  character(len=*), parameter :: ModuleName= "$RCSfile$"
!------------------------------------------------------------------------

! Contents:

! Subroutines -- DestroyChunkDatabase             Public
!                ScanDivide_firstGuess            Private
!                ScanDivide_wall                  Private
!                ScanDivide_range                 Private
!                ScanDivide_checkReorganizeRange  Private
!                ScanDivide_mlscf                 Private
!                ScanAndDivide                    Public

! Remarks:  This prototype module contains subroutines for the Scan/Divide task
!           of the L2 software.

contains ! =====     Public Procedures     =============================

!------------------------------------------  DestroyChunkDatabase  -----
  subroutine DestroyChunkDatabase ( chunks )
    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS
    integer :: STATUS ! From deallocate

    deallocate ( chunks, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_DeAllocate // "Chunks" )
  end subroutine DestroyChunkDatabase

!-------------------------------------------------  ScanAndDivide  -----
  subroutine ScanAndDivide ( root, processingRange, l1bInfo, chunks )
!-----------------------------------------------------------------------

! Brief description of subroutine 
! This subroutine is called by the main L2 program to divide the input L1B
! dataset into chunks.

! Arguments

    integer, intent(in) :: ROOT    ! Root of the L2CF tree for ChunkDivide

    type( L1BInfo_T ), intent(in) :: L1BINFO

    type( TAI93_Range_T ), intent(in) :: PROCESSINGRANGE

    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS

! Parameters

    character (len=*), parameter :: MLSMSG_DeallocateL1b = &
                                      'Deallocation failed:  l1b pointers.'
    character (len=*), parameter :: MLSMSG_SUB = &
                                      'Error return from subroutine '

! Functions

! Variables

    type( L1BData_T ) :: DATA

    integer :: MODHOME
    integer :: MODCRITICAL
    character (len=480) :: MSR
    character (len=name_len) :: QUANTITY
    integer :: UNITSGAP
    integer :: BANDS

    real(r8) :: DIFF, HOME, LLB, LUB, MAX, MAXGAP, MIN, ORBLEN
    real(r8) :: ULB, UUB
    real(r8), allocatable :: PHI(:), TIME(:)

    integer :: ERROR, FIRSTMAF, FLAG, FOUNDEND, FOUNDSTART, I, J
    integer :: LASTMAF, LOOP, NOMAFS, NUMBAD, NUMCHUNKS, NUMFILE, NUMOA
    integer :: NUMSCAN, OVERLAP
    integer, allocatable :: BADMAF(:), COUNTERMAF(:), FIRSTOA(:)
    integer, allocatable :: FIRSTSCAN(:), LASTOA(:), LASTSCAN(:)

! Check that input times are reasonable

    if (processingRange%startTime > processingRange%endTime) &
       call MLSMessage(MLSMSG_Error, ModuleName, 'Input start time &
                                                 &greater than end time.')

! Read the L1 data for MAF times

    quantity = 'MAFStartTimeTAI'
    call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
    if (flag /= 0) then
       msr = MLSMSG_L1BRead // quantity
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

! Allocate arrays which will hold re-used L1BOA quantities for the entire file

    allocate ( time(noMAFs), counterMAF(noMAFs), phi(noMAFs), STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Allocate // ' time, phi, counterMAF.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    time = data%dpField(1,1,:)
    counterMAF = data%counterMAF

    call DeallocateL1BData(data, flag)
    if (flag /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
                                   MLSMSG_DeallocateL1b)

! Check that input times fall within file

    if ( (processingRange%startTime >= time(noMAFs)) .OR. &
         (processingRange%endTime <= time(1)) ) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
                      'Input times out of range of L1 input file.')
   
! Get MLSCF values

    call ScanDivide_mlscf ( root, bands, home, llb, lub, &
          maxGap, modCritical, modHome, orbLen, overlap, ulb, uub, unitsGap )

! Read phi for correct module

    if (modHome == l_thz) then
      quantity = 'THz.tpGeodAngle'
    else
      quantity = 'GHz.tpGeodAngle'
    end if

    call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
    if (flag /= 0) then
      msr = MLSMSG_L1BRead // quantity
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    phi = data%dpField(1,1,:)

    call DeallocateL1BData(data, flag)
    if (flag /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
                                   MLSMSG_DeallocateL1b)

! Find firstMAF for data set

    if ( processingRange%startTime < time(1) ) then
      firstMAF = 0
    else
      foundStart = 0
      do i = 2, noMAFs
        if (foundStart == 1) exit
        if ( time(i) >= processingRange%startTime ) then
          firstMAF = i-1
          foundStart = 1
        end if
      end do
    end if

! Find lastMAF for data set

    numFile = noMAFs-1

    if ( processingRange%endTime >= time(noMAFs) ) then
      lastMAF = numFile
    else
      foundEnd = 0
      do i = 1, numFile
        if (foundEnd == 1) exit
        if ( time(noMAFs-i) < processingRange%endTime ) then
          lastMAF = noMAFs-i-2
          foundEnd = 1
        end if
      end do
    end if

! Calculate a first guess for chunk boundaries

    call ScanDivide_firstGuess(l1bInfo%L1BOAId, firstMAF, home, lastMAF, &
                               numFile, orbLen, overlap, phi, chunks, flag)
    if (flag /= 0) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_SUB // 'ScanDivide_firstGuess')

! If necessary, check for completely missing telemetry

    if (maxGap == -1.0) return

    allocate ( badMAF(noMAFs), STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Allocate // ' holder for bad MAF numbers.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    badMAF = 0
    numBad = 0
    numChunks = SIZE(chunks)
    do i = chunks(1)%firstMAFIndex+2, chunks(numChunks)%lastMAFIndex

      select case ( unitsGap )
      case ( phyq_mafs )
        diff = counterMAF(i) - counterMAF(i-1)
      case ( phyq_time )
        diff = time(i) - time(i-1)
      case default
        diff = phi(i) - phi(i-1)
      end select

      if (diff > maxGap) then
        badMAF(i) = 1
        numBad = numBad + 1
      end if

    end do

    if (numBad > 0) then

      call ScanDivide_wall ( chunks, badMAF, flag )
      if (flag /= 0) then
        msr = MLSMSG_SUB // 'ScanDivide_wall'
        call MLSMessage(MLSMSG_Warning, ModuleName, msr)
      end if

    end if
   
! Scan scGeodAngle for value flagging bad toolkit returns

    quantity = 'scGeodAngle'
    call ReadL1BData ( l1bInfo%L1BOAId, quantity, data, noMAFs, flag )
    if (flag /= 0) then
      msr = MLSMSG_L1BRead // quantity
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

! Check new set of chunks for unusable oa data

    allocate ( firstOa(noMAFs), lastOa(noMAFs), STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Allocate // ' first and last bad Oa MAF holders.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    firstOa = -1
    lastOa = -1
    numOa = 0
    numBad = 1
    numChunks = SIZE(chunks)
    do i = chunks(1)%firstMAFIndex, chunks(numChunks)%lastMAFIndex

      if ( data%dpField(1,1,i+1) <= -999.9 ) then

        if ( firstOa(numBad) == -1 ) firstOa(numBad) = i
        lastOa(numBad) = i
        numOa = numBad

      else if ( i > 0 ) then
        if ( data%dpField(1,1,i) <= -999.9 ) then
          numBad = numBad + 1
        end if
      end if

    end do

    call DeallocateL1BData(data, flag)
    if (flag /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
                                   MLSMSG_DeallocateL1b)

! If gaps are found, check whether they (along with any previous adjacent
! ones) exceed MaxGap; re-organize chunks, accordingly

    if ( numOa > 0 ) then
      call ScanDivide_checkReorganizeRange(counterMAF, firstOa(1:numOa), &
                       lastOa(1:numOa), maxGap, phi, time, unitsGap, chunks)
    end if

! If necessary, check for valid scans

    if ( modCritical /= l_none ) then

      badMAF = 0
      numBad = 0
      numChunks = SIZE(chunks)
   
! Unless THz alone specified as critical module, read GHz first

      if ( modCritical /= l_thz ) then

        quantity = 'GHz.tpGeodAlt'
        call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
        if (flag /= 0) then
          msr = MLSMSG_L1BRead // quantity
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if

! Check GHz for bad scans

        do i = chunks(1)%firstMAFIndex+1, chunks(numChunks)%lastMAFIndex+1

          min = MINVAL(data%dpField(1,:,i))
          max = MAXVAL(data%dpField(1,:,i))

          if ( (min < llb) .OR. (min > lub) .OR. &
               (max < ulb) .OR. (max > uub) ) then

            badMAF(i) = 1
            numBad = numBad + 1

          end if

        end do

        call DeallocateL1BData(data, flag)
        if (flag /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
                                       MLSMSG_DeallocateL1b)

      end if

! If necessary, read THz data

      if ( (modCritical == l_thz) .OR. (modCritical == l_both) .OR. &
           ((modCritical == l_either) .AND. (numBad > 0)) ) then

        quantity = 'THz.tpGeodAlt'
        call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
        if (flag /= 0) then
          msr = MLSMSG_L1BRead // quantity
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if

! If EITHER is critical, check only bad GHz scans; if THz okay,
! unmark MAF as bad.

        if ( modCritical == l_either ) then

          do i = chunks(1)%firstMAFIndex+1,chunks(numChunks)%lastMAFIndex+1

            if (badMAF(i) == 1) then

              min = MINVAL(data%dpField(1,:,i))
              max = MAXVAL(data%dpField(1,:,i))

              if ( (min >= llb) .AND. (min <= lub) .AND. &
                   (max >= ulb) .AND. (max <= uub) ) then

                badMAF(i) = 0
                numBad = numBad - 1

              end if

            end if

          end do

! If THz is critical, check entire data set for bad THz scans

        else

          do i = chunks(1)%firstMAFIndex+1,chunks(numChunks)%lastMAFIndex+1

            min = MINVAL(data%dpField(1,:,i))
            max = MAXVAL(data%dpField(1,:,i))
            if ( (min < llb) .OR. (min > lub) .OR. &
                 (max < ulb) .OR. (max > uub) ) then

              badMAF(i) = 1
              numBad = numBad + 1

            end if

          end do

        end if

        call DeallocateL1BData ( data, flag )
        if (flag /= 0) call MLSMessage(MLSMSG_Warning, ModuleName, &
                                       MLSMSG_DeallocateL1b)

      end if

! If bad scans found, add "bad" markings for any MAFs previously discarded

      if (numBad > 0) then

        do i = 1, numOa
          do j = firstOa(i), lastOa(i)
            badMAF(j+1) = 1
          end do
        end do

! Set first & last badMAF for all gaps found thus far

        allocate ( firstScan(numBad+numOa), lastScan(numBad+numOa), &
                 STAT=error )
        if (error /= 0) then
          msr = MLSMSG_Allocate // ' first and last bad scan MAF holders.'
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if

        firstScan = -1
        lastScan = -1
        loop = 0
        numScan = 1

        do i = chunks(1)%firstMAFIndex+1, chunks(numChunks)%lastMAFIndex+1

          if ( badMAF(i) == 1) then

            if (firstScan(numScan) == -1) firstScan(numScan) = i-1
            lastScan(numScan) = i-1
            loop = numScan

          else if ( (i > 1) .AND. (lastScan(numScan) == i-2) ) then

            numScan = numScan + 1

          end if

        end do

! Check whether these (including MAF "walls" previously found) exceed MaxGap;
! re-organize chunks, if necessary

        call ScanDivide_checkReorganizeRange(counterMAF, &
                              firstScan(1:loop), lastScan(1:loop), &
                              maxGap, phi, time, unitsGap, chunks)

        deallocate ( firstScan, lastScan, STAT=error )
        if (error /= 0) then
          msr = MLSMSG_Deallocate // ' first & lastScan local variables.'
          call MLSMessage(MLSMSG_Warning, ModuleName, msr)
        end if

      end if

    end if

! If necessary, check for critical bands in the L1bRad file
    numChunks = size(chunks)
    if ( bands /= 0 ) then
    end if

! Deallocations

    deallocate ( badMAF, counterMAF, firstOa, lastOa, phi, time, STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Deallocate // ' local variables.'
      call MLSMessage(MLSMSG_Warning, ModuleName, msr)
    end if

!------------------------------
  end subroutine ScanAndDivide
!------------------------------
! =====     Private Procedures     =====================================
!-----------------------------------------  ScanDivide_firstGuess  -----
  subroutine ScanDivide_firstGuess ( L1boaFileHandle, firstMAF, home, &
                                     lastMAF, numFile, orbLen, overlap, &
                                     phi, chunks, flag)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine creates a set of data structures defining a "first guess" for
! MAF chunk boundaries.

! Arguments

    integer, intent(in) :: FIRSTMAF, LASTMAF, NUMFILE, OVERLAP
    integer, intent(in) :: L1BOAFILEHANDLE

    real(r8), intent(in) :: HOME, ORBLEN
    real(r8), intent(in) :: PHI(:)

    integer, intent(out) :: FLAG

    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS

! Parameters

! Functions

! Variables

    character (LEN=480) :: MSR

    integer :: ABOVE, ALLOC_ERR, BDRYMAX, BELOW, DEALLOC_ERR, I
    integer :: NUMBDRY, NUMCHUNKS
    integer, allocatable :: BDRYMAF(:), INDICES(:)

    real(r8) :: MAX, MIN, NEWHOME, PHILEN
    real(r8), allocatable :: ABOVEHOME(:), BELOWHOME(:), BDRYPHI(:)

! Find max & min phi for data set

    max = MAXVAL( phi(firstMAF+1:lastMAF+1) )
    min = MINVAL( phi(firstMAF+1:lastMAF+1) )

! Convert from "fraction of orbit" to "degrees of phi"

    phiLen = orbLen

! Find phi boundaries above & below HomeGeodAngle

    bdryMax = ANINT(max_orbits * 360.0 / orbLen)

    allocate (  aboveHome(bdryMax), belowHome(bdryMax), STAT=alloc_err )
    if ( alloc_err /= 0 ) then
       msr = MLSMSG_Allocate // ' arrays for phi bdries above & below Home.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

! If min above home, reset HomeGeodAngle to closest ideal bdry below min

    if (min >= home) then

       below = 0

       do i = 1, bdryMax
          newHome = home + i*phiLen
          if (newHome > min) exit
       end do

       newHome = home + (i-1)*phiLen

! Calculate boundaries above home phi

       do i = 1, bdryMax
          aboveHome(i) = newHome + i*phiLen
          if (aboveHome(i) >= max) exit
       end do

       above = i-1

    else
      
! If min below home, check whether max also is

        if (max <= home) then

! If max is also below, reset HomeGeodAngle to closest ideal bdry above max

          above = 0

          do i = 1, bdryMax
             newHome = home - i*phiLen
             if (newHome <= max) exit
          end do

          newHome = home - (i-1)*phiLen

       else

! If max is above, calculate boundaries >= HomeGeodAngle

          do i = 1, bdryMax
             aboveHome(i) = home + (i-1)*phiLen
             if (aboveHome(i) >= max) exit
          end do

          above = i-1

          newHome = home

       end if

! Calculate boundaries below home phi

       do i = 1, bdryMax
          belowHome(i) = newHome - i*phiLen
          if (belowHome(i) <= min) exit
       end do

       below = i-1

    end if

! Put phi boundaries into a single array in ascending order, from min to max

    numBdry = above + below + 2

    allocate (  bdryPhi(numBdry), bdryMAF(numBdry), indices(numBdry), &
              STAT=alloc_err )
    if ( alloc_err/= 0 ) then
       msr = MLSMSG_Allocate // ' bdry arrays.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    bdryPhi(1) = min

    do i = 1, below
      bdryPhi(i+1) = belowHome(below-(i-1))
    end do

    do i = 1, above
      bdryPhi(below+1+i) = aboveHome(i)
    end do

    bdryPhi(above+below+2) = max

    deallocate ( aboveHome, belowHome, STAT=dealloc_err)
    if (dealloc_err /= 0) then
       msr = MLSMSG_Deallocate // ' arrays for phi bdries above & below Home.'
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
       flag = -1
    end if

! Find MAFs containing the phi boundaries

    call Hunt(phi, bdryPhi, indices, allowTopValue=.TRUE.)

! Convert MAF bdry indices to begin at 0 (like file), rather than 1 (like phi)

    bdryMAF = indices - 1

    deallocate ( bdryPhi, indices, STAT=dealloc_err)
    if (dealloc_err /= 0) then
       msr = MLSMSG_Deallocate // ' phi, index bdries.'
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
       flag = -1
    end if

! Fill data structure for each chunk

    numChunks = numBdry-1
    allocate (  chunks(numChunks), STAT=alloc_err )
    if ( alloc_err/= 0 ) call MLSMessage  (MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // ' number of chunks.' )

! Calculate firstMAFIndex, subtracting overlap from MAF boundaries

    chunks%firstMAFIndex = bdryMAF(1:numChunks) - overlap
    if (chunks(1)%firstMAFIndex < 0 ) chunks(1)%firstMAFIndex = 0

! Calculate lastMAFIndex, adding overlap to MAF boundaries

    chunks(1:numChunks-1)%lastMAFIndex = bdryMAF(2:numBdry-1) - 1 + overlap
    chunks(numChunks)%lastMAFIndex = bdryMAF(numBdry) + overlap

    if (chunks(numChunks)%lastMAFIndex > numFile ) &
        chunks(numChunks)%lastMAFIndex = numFile

! Set lower overlap to MLSCF input; consider first MAF separately

    chunks(2:numChunks)%noMAFsLowerOverlap = overlap

    chunks(1)%noMAFsLowerOverlap = bdryMAF(1) - chunks(1)%firstMAFIndex

! Set upper overlap to MLSCF input; consider last MAF separately

    chunks(1:numChunks-1)%noMAFsUpperOverlap = overlap

    chunks(numChunks)%noMAFsUpperOverlap = chunks(numChunks)%lastMAFIndex - &
                                          &bdryMAF(numBdry)

! Calculate accumulatedMAFs

    chunks(1)%accumulatedMAFs = 0

    do i = 2, numChunks
       chunks(i)%accumulatedMAFs = chunks(i-1)%accumulatedMAFs + &
        &chunks(i-1)%lastMAFIndex - chunks(i-1)%firstMAFIndex + 1 - &
        &chunks(i-1)%noMAFsLowerOverlap - chunks(i-1)%noMAFsUpperOverlap
    end do

    deallocate ( bdryMAF, STAT=dealloc_err)
    if (dealloc_err /= 0) then
       msr = MLSMSG_Deallocate // ' MAF bdries.'
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
       flag = -1
    end if

!--------------------------------------
  end subroutine ScanDivide_firstGuess
!--------------------------------------

!-----------------------------------------------  ScanDivide_wall  -----
  subroutine ScanDivide_wall(chunks, bad, flag)
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reorganizes chunks around "walls."

! Arguments

    integer, intent(in) :: BAD(:)

    integer, intent(out) :: FLAG

    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS

! Parameters

! Functions

! Variables

    character (len=480) :: MSR

    integer :: ERROR, I, J, MOVE, NUMBAD, NUMCHUNKS, SCALAR, WALLIND
    integer, allocatable :: BADMAF(:)

    type( MLSChunk_T ), dimension(:), pointer  :: WALLEDCHUNKS

! If input is a single MAF #, change its format to that of an array,
! indexed such that badMAF(input) = 1

    numBad = SIZE(bad)

    if (numBad == 1) then

       scalar = bad(1)

       allocate (  badMAF(scalar), STAT=error )
       if (error /= 0) then
          msr = MLSMSG_Allocate // ' badMAF array.'
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       end if

       badMAF = 0
       badMAF(scalar) = 1

    else

       allocate (  badMAF(numBad), STAT=error )
       if (error /= 0) then
          msr = MLSMSG_Allocate // ' badMAF array.'
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       end if

       badMAF = bad

    end if

! Initialize walledChunks array

    numChunks = SIZE(chunks)

    allocate ( walledChunks(numChunks+numBad), STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Allocate // ' walled chunk arrays.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    walledChunks(1) = chunks(1)

! If wall in first chunk splits off a portion containing only lower overlap,
! jettison overlap portion

    if (chunks(1)%noMAFsLowerOverlap > 0) then

       do i = 1, chunks(1)%noMAFsLowerOverlap

          if ( badMAF(chunks(1)%firstMAFIndex+&
                      chunks(1)%noMAFsLowerOverlap-i+2) == 1) then
             walledChunks(1)%firstMAFIndex = chunks(1)%firstMAFIndex+&
                                           chunks(1)%noMAFsLowerOverlap-(i-1)
             walledChunks(1)%noMAFsLowerOverlap = i-1
             exit
          end if

       end do

    end if

! If wall in non-overlapped portion, split chunk into 2, with no overlaps at
! new bdry.

    move = 0
    wallInd = 1

    do i = chunks(1)%firstMAFIndex+chunks(1)%noMAFsLowerOverlap+1, &
           chunks(2)%firstMAFIndex-1

       if (badMAF(i+1) == 1) then
          walledChunks(wallInd)%lastMAFIndex = i-1
          walledChunks(wallInd)%noMAFsUpperOverlap = 0
          walledChunks(wallInd+1)%firstMAFIndex = i
          walledChunks(wallInd+1)%noMAFsLowerOverlap = 0
          walledChunks(wallInd+1)%accumulatedMAFs = i
          walledChunks(wallInd+1)%lastMAFIndex = chunks(1)%lastMAFIndex
          walledChunks(wallInd+1)%noMAFsUpperOverlap = &
                        chunks(1)%noMAFsUpperOverlap
          wallInd = wallInd + 1
       end if

    end do

! If in upper overlap, re-set bdries without overlaps

    do i = chunks(2)%firstMAFIndex, chunks(1)%lastMAFIndex+1

       if (badMAF(i+1) == 1) then
          walledChunks(wallInd)%lastMAFIndex = i-1
          walledChunks(wallInd)%noMAFsUpperOverlap = 0
          walledChunks(wallInd+1)%firstMAFIndex = i
          walledChunks(wallInd+1)%noMAFsLowerOverlap = 0
          walledChunks(wallInd+1)%accumulatedMAFs = i
          walledChunks(wallInd+1)%lastMAFIndex = chunks(2)%lastMAFIndex
          walledChunks(wallInd+1)%noMAFsUpperOverlap = &
                        chunks(2)%noMAFsUpperOverlap
          wallInd = wallInd + 1
          move = 1
       end if

    end do

! Re-organize rest of chunks, except for last

    do i = 2, numChunks-1

       if (move == 1) then
          move = 0
       else
          wallInd = wallInd + 1
          walledChunks(wallInd) = chunks(i)
       end if

       do j = chunks(i-1)%lastMAFIndex+2, chunks(i+1)%firstMAFIndex-1

          if (badMAF(j+1) == 1) then
             walledChunks(wallInd)%lastMAFIndex = j-1
             walledChunks(wallInd)%noMAFsUpperOverlap = 0
             walledChunks(wallInd+1)%accumulatedMAFs = j
             walledChunks(wallInd+1)%firstMAFIndex = j
             walledChunks(wallInd+1)%noMAFsLowerOverlap = 0
             walledChunks(wallInd+1)%lastMAFIndex = chunks(i)%lastMAFIndex
             walledChunks(wallInd+1)%noMAFsUpperOverlap = &
                           chunks(i)%noMAFsUpperOverlap
             wallInd = wallInd + 1
          end if

       end do

       do j = chunks(i+1)%firstMAFIndex, chunks(i)%lastMAFIndex+1

          if (badMAF(j+1) == 1) then
             walledChunks(wallInd)%lastMAFIndex = j-1
             walledChunks(wallInd)%noMAFsUpperOverlap = 0
             walledChunks(wallInd+1)%accumulatedMAFs = j
             walledChunks(wallInd+1)%firstMAFIndex = j
             walledChunks(wallInd+1)%noMAFsLowerOverlap = 0
             walledChunks(wallInd+1)%lastMAFIndex = chunks(i+1)%lastMAFIndex
             walledChunks(wallInd+1)%noMAFsUpperOverlap = &
                         chunks(i+1)%noMAFsUpperOverlap
             wallInd = wallInd + 1
             move = 1
          end if

       end do

    end do

! Re-organize non-overlapped portion of last chunk

    if (move == 0) then
       wallInd = wallInd + 1
       walledChunks(wallInd) = chunks(numChunks)
    end if

    do i = chunks(numChunks-1)%lastMAFIndex+2, &
        chunks(numChunks)%lastMAFIndex-chunks(numChunks)%noMAFsUpperOverlap

       if (badMAF(i+1) == 1) then
          walledChunks(wallInd)%lastMAFIndex = i-1
          walledChunks(wallInd)%noMAFsUpperOverlap = 0
          walledChunks(wallInd+1)%accumulatedMAFs = i
          walledChunks(wallInd+1)%firstMAFIndex = i
          walledChunks(wallInd+1)%noMAFsLowerOverlap = 0
          walledChunks(wallInd+1)%lastMAFIndex = &
                chunks(numChunks)%lastMAFIndex
          walledChunks(wallInd+1)%noMAFsUpperOverlap = &
                chunks(numChunks)%noMAFsUpperOverlap
          wallInd = wallInd +1
       end if

    end do

! Re-organize upper overlap portion (if it exists); truncate any resulting
! chunk that is purely overlap
 
    do i = chunks(numChunks)%lastMAFIndex-&
        chunks(numChunks)%noMAFsUpperOverlap+1,chunks(numChunks)%lastMAFIndex

       if (badMAF(i+1) == 1) then

         walledChunks(wallInd)%lastMAFIndex = i-1
         walledChunks(wallInd)%noMAFsUpperOverlap = &
          walledChunks(wallInd)%lastMAFIndex-chunks(numChunks)%lastMAFIndex+&
          chunks(numChunks)%noMAFsUpperOverlap

         if (walledChunks(wallInd)%noMAFsUpperOverlap < 0) then
            walledChunks(wallInd)%noMAFsUpperOverlap = 0
         end if

         exit

       end if

    end do

! Deallocations

    deallocate ( badMAF, STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Deallocate // ' badMAF array.'
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
       flag = -1
    end if

    deallocate ( chunks, STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Deallocate // ' chunks array.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

! Reset input chunks array

    allocate ( chunks(wallInd), STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Allocate // ' new output chunks array.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    chunks = walledChunks(1:wallInd)

    deallocate ( walledChunks, STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Deallocate // ' walledChunks array.'
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
       flag = -1
    end if

!--------------------------------
  end subroutine ScanDivide_wall
!--------------------------------

!----------------------------------------------  ScanDivide_range  -----
  subroutine ScanDivide_range ( chunks, firstBadMAFIndex, lastBadMAFIndex, &
                                flag )
!------------------------------------------------------------------------------

! Brief description of subroutine
! This subroutine reorganizes chunks around periods of bad data.

! Arguments

    integer, intent(in) :: firstBadMAFIndex, lastBadMAFIndex

    integer, intent(out) :: flag

    type( MLSChunk_T ), dimension(:), pointer  :: chunks

! Parameters

! Functions

! Variables

    character (len=480) :: MSR

    integer :: CHUNKIND, ERROR, I, J, NEWIND, NUMCHUNKS

    type( MLSChunk_T ), dimension(:), pointer  :: NEWCHUNKS

! Initialize expanded newChunks array

    numChunks = SIZE(chunks)

    allocate ( newChunks(numChunks+1), STAT=error )
    if (error /= 0) then
       msr = MLSMSG_Allocate // ' newChunks array.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

! Find the chunk where the gap begins

    newInd = 1
    chunkInd = 1
    do i = 1, numChunks

      if ( (firstBadMAFIndex <= chunks(chunkInd)%lastMAFIndex) .AND. &
           (firstBadMAFIndex >= chunks(chunkInd)%firstMAFIndex) ) then

! If it begins above the lower overlap,

        if (firstBadMAFIndex > chunks(chunkInd)%firstMAFIndex+&
                               chunks(chunkInd)%NoMAFsLowerOverlap ) then

! It breaks off a chunk below the gap; set info for this new chunk

          newChunks(newInd)%accumulatedMAFs = &
                            chunks(chunkInd)%accumulatedMAFs
          newChunks(newInd)%firstMAFIndex = chunks(chunkInd)%firstMAFIndex
          newChunks(newInd)%noMAFsLowerOverlap = &
                            chunks(chunkInd)%noMAFsLowerOverlap
          newChunks(newInd)%lastMAFIndex = firstBadMAFIndex - 1

          if ( firstBadMAFIndex > chunks(numChunks)%lastMAFIndex-&
                              chunks(numChunks)%noMAFsUpperOverlap+1) then
            newChunks(newInd)%noMAFsUpperOverlap = &
                              newChunks(newInd)%lastMAFIndex-&
                              chunks(numChunks)%lastMAFIndex+&
                              chunks(numChunks)%noMAFsUpperOverlap
          else
            newChunks(newInd)%noMAFsUpperOverlap = 0
          end if

          newInd = newInd + 1

       end if

! If the gap ends in the upper overlap of the last chunk, exit

       if ( lastBadMAFIndex >= chunks(numChunks)%lastMAFIndex-&
                               chunks(numChunks)%noMAFsUpperOverlap ) exit

! Otherwise, find the chunk where the gap ends

       do j = 1, numChunks

         chunkInd = chunkInd + (j-1)

         if ( lastBadMAFIndex <= chunks(chunkInd)%lastMAFIndex ) then

! Set start info for a new chunk beginning after the gap

           if (newInd == 1) then
              newChunks(newInd)%accumulatedMAFs = 0
           else
              newChunks(newInd)%accumulatedMAFs = &
                                newChunks(newInd-1)%lastMAFIndex + 1
           end if

           newChunks(newInd)%firstMAFIndex = lastBadMAFIndex + 1

           newChunks(newInd)%noMAFsLowerOverlap = &
                             chunks(chunkInd)%firstMAFIndex + &
                             chunks(chunkInd)%noMAFsLowerOverlap - &
                             newChunks(newInd)%firstMAFIndex
           if ( newChunks(newInd)%noMAFsLowerOverlap < 0 ) &
                newChunks(newInd)%noMAFsLowerOverlap = 0

! If gap ends below the upper overlap, set end info to the current old chunk

           if (lastBadMAFIndex < chunks(chunkInd)%lastMAFIndex-&
                              chunks(chunkInd)%noMAFsUpperOverlap ) then

              newChunks(newInd)%lastMAFIndex = &
                                chunks(chunkInd)%lastMAFIndex
              newChunks(newInd)%noMAFsUpperOverlap = &
                                chunks(chunkInd)%NoMAFsUpperOverlap
              chunkInd = chunkInd + 1

           else

! If not, set end info to next old chunk

             newChunks(newInd)%lastMAFIndex = &
                               chunks(chunkInd+1)%lastMAFIndex
             newChunks(newInd)%noMAFsUpperOverlap = &
                               chunks(chunkInd+1)%NoMAFsUpperOverlap
             chunkInd = chunkInd + 2

           end if

           newInd = newInd + 1
           exit

         end if

       end do

     else

! This chunk has no gaps in it; retain old chunk info

       newChunks(newInd) = chunks(chunkInd)
       newInd = newInd + 1
       chunkInd = chunkInd + 1

     end if

     if (chunkInd > numChunks) exit

   end do

! Reset input chunks array

   deallocate ( chunks, STAT=error )
   if (error /= 0) then
     msr = MLSMSG_Deallocate // ' chunks array.'
     call MLSMessage(MLSMSG_Error, ModuleName, msr)
   end if

   allocate ( chunks(newInd-1), STAT=error )
   if (error /= 0) then
     msr = MLSMSG_Allocate // ' new output chunks array.'
     call MLSMessage(MLSMSG_Error, ModuleName, msr)
   end if

   chunks = newChunks(1:newInd-1)

   deallocate ( newChunks, STAT=error )
   if (error /= 0) then
     msr = MLSMSG_Deallocate // ' newChunks array.'
     call MLSMessage(MLSMSG_Warning, ModuleName, msr)
     flag = -1
   end if

!---------------------------------
  end subroutine ScanDivide_range
!---------------------------------

!-------------------------------  ScanDivide_checkReorganizeRange  -----
  subroutine ScanDivide_checkReorganizeRange ( counterMAF, firstBad, &
                                               lastBad, maxGap, phi, &
                                               time, units, chunks )
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine checks data gaps against a maximum and calls the "range
! re-organization routine, if MaxGap is exceeded.

! Arguments

    integer, intent(in) :: units

    integer, intent(in) :: counterMAF(:), firstBad(:), lastBad(:)

    real(r8), intent(in) :: maxGap
    real(r8), intent(in) :: phi(:), time(:)

    type( MLSChunk_T ), dimension(:), pointer  :: chunks
   
! Parameters

! Functions

! Variables

    real(r8) :: DIFF

    integer :: END, FLAG, I, NUMCHUNKS, NUMGAPS, START

! Add completely missing telemetry to bad MAFs, and check against MaxGap

    numGaps = SIZE(firstBad)
    numChunks = SIZE(chunks)

    do i = 1, numGaps

      if ( firstBad(i) == chunks(1)%firstMAFIndex) then
        start = chunks(1)%firstMAFIndex + 1
      else
        start = firstBad(i)
      end if

      if ( lastBad(i) == chunks(numChunks)%lastMAFIndex ) then
        end = chunks(numChunks)%lastMAFIndex + 1
      else
        end = lastBad(i) + 2
      end if

      select case ( units )
      case ( phyq_mafs )
        diff = counterMAF(end) - counterMAF(start)
      case ( phyq_time )
        diff = time(end) - time(start)
      case default
        diff = phi(end) - phi(start)
      end select

! If diff exceeds MaxGap, re-organize chunks

      if (diff > maxGap) &
        call ScanDivide_range(chunks, firstBad(i), lastBad(i), flag)

    end do

!------------------------------------------------
  end subroutine ScanDivide_checkReorganizeRange
!------------------------------------------------

!----------------------------------------------  ScanDivide_mlscf  -----
  subroutine ScanDivide_mlscf ( root, bands, home, llb, lub, maxGap, &
                   modCritical, modHome, orbLen, overlap, ulb, uub, unitsGap )
!-----------------------------------------------------------------------

! Brief description of subroutine
! This subroutine identifies, separates, and checks values from the section
! of the MLSCF (ChunkDivide) passed to Scan/Divide.

! Arguments

    integer, intent(in) :: ROOT    ! Root of the ChunkDivide section of the
                                   ! MLSCF abstract syntax tree
    integer, intent(out) :: modHome     ! Home module
    integer, intent(out) :: modCritical ! Critical scanning module
    integer, intent(out) :: bands       ! String index
    integer, intent(out) :: unitsGap    ! PHYQ_... from Init_Tables_Module

    integer, intent(out) :: overlap

    real(r8), intent(out) :: home, llb, lub, maxGap, orbLen, ulb, uub

! Parameters

    character (len=*), parameter :: MLSMSG_CF = &
                                    ' not specified in the configuration file.'
! Functions

! Variables

    logical :: GOT(first_parm:last_parm) = .false.
    integer :: I         ! Loop inductor
    integer :: KEY       ! A P_... parameter from Init_Tables_Module
    character (len=480) :: MSR
    integer :: SON       ! A son of the ChunkDivide section node
    integer :: UNITS(2)  ! Units of expression
    double precision :: VALUE(2)   ! Value of expression

! Initialize variables to 'unfound' values

    modHome = l_ghz
    llb = 999.9
    lub = 999.9
    ulb = -999.9
    uub = -999.9
    bands = 0               ! String index, = 0 means no decl specified
    maxGap = -1.0
    unitsGap = phyq_invalid

! Loop through MLSCF section, identifying keywords and setting variables

    do i = 2, nsons(root)-1 ! Skip the section identifiers
      son = subtree(i,root)
      key = decoration(subtree(1,son)) ! P_... index from Init_Tables_Module
      got(key) = .true.
      call expr ( subtree(2,son), units, value )

      select case ( key )
      case ( p_ideal_length )
        orbLen = value(1)
      case ( p_overlap )
        overlap = value(1)
      case ( p_home_geod_angle )
        home = value(1)
      case ( p_home_module )
        modHome = value(1)
      case ( p_scan_lower_limit )
        if ( units(1) /= phyq_Length) then
          call MLSMessage( MLSMSG_Error, ModuleName, &
            &              'ScanLowerLimit not specified as a length.')
        end if
        llb = value(1)
        lub = value(2)
      case ( p_scan_upper_limit )
        if ( units(1) /= phyq_Length) then
          call MLSMessage( MLSMSG_Error, ModuleName, &
            &              'ScanUpperLimit not specified as a length.')
        end if
        ulb = value(1)
        uub = value(2)
      case ( p_critical_scanning_modules )
        modCritical = value(1)
      case ( p_critical_bands )
        bands = value(1)
      case ( p_max_gap )
        maxGap = value(1)
        unitsGap = units(1)
      end select

    end do

! Check for missing values that would cause error exits

    do i = first_parm, last_parm
      select case ( i )
      case ( p_ideal_length, p_overlap, p_home_module, &
             p_critical_scanning_modules )
        if ( .not. got(i) ) then
          call get_string ( parm_indices(i), msr )
          msr(string_length(parm_indices(i))+1:) = MLSMSG_CF
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if
      end select
    end do

! Check for values that depend on the presence of other quantities

    if ( .not. got(p_scan_lower_limit) .AND. (modCritical /= l_none) ) then
      msr = 'ScanLowerLimit' // MLSMSG_CF
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    if ( .not. got(p_scan_upper_limit) .AND. (modCritical /= l_none) ) then
      msr = 'ScanUpperLimit' // MLSMSG_CF
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

!--------------------------------
  end subroutine ScanDivide_mlscf
!--------------------------------

!====================
end module ScanDivide
!====================

!# $Log$
!# Revision 2.2  2000/09/11 20:00:20  ahanzel
!# Removed old log entries in file.
!#
!# Revision 2.1  2000/09/08 22:55:56  vsnyder
!# Revised to use the tree output by the parser
!#
