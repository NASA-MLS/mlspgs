! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================
module ScanDivide
  !=======================================================================

  use Dumper, only : DUMP
  use EXPR_M, only: EXPR   
  use INIT_TABLES_MODULE, only:  FIRST_PARM, L_BOTH, L_EITHER, L_GHZ, &
    L_NONE, L_THZ, LAST_PARM, L_TRUE
  use INIT_TABLES_MODULE, only: P_CRITICAL_BANDS, P_CRITICAL_SCANNING_MODULES, &
    P_HOME_GEOD_ANGLE, P_HOME_MODULE, P_IDEAL_LENGTH, P_IGNOREL1B, P_MAX_GAP, &
    P_NOCHUNKS, P_OVERLAP, P_SCAN_LOWER_LIMIT, P_SCAN_UPPER_LIMIT, S_TIME
  use INIT_TABLES_MODULE, only: PHYQ_ANGLE, PHYQ_INVALID, PHYQ_LENGTH, &
    PHYQ_MAFS, PHYQ_TIME
  use Intrinsic, only: PARM_INDICES
  use L1BData, only: deallocateL1BDATA, L1BDATA_T, NAME_LEN, READL1BDATA
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning, MLSMSG_L1BRead
  use MLSNumerics, only: HUNT, R8
  !  use MLSStrings, only: MLSMSG_L1BRead
  use MoreTree, only: Get_Spec_ID
  use Output_M, only: Output
  use SDPToolkit, only: MAX_ORBITS
  use STRING_TABLE, only: Display_String
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, Node_ID, NSONS, Source_Ref, SUBTREE
  use Tree_Types, only: N_Equal, N_Named

  implicit none

  private

  public :: DestroyChunkDatabase, ScanAndDivide

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

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

  ! Some private types

  type ScanDivideConfig_T
    logical   :: ignoreL1B              ! Ignore l1 file completely
    integer   :: noChunks               ! Force this many chunks (0=ignore)
    real(r8)  :: idealLength            ! Ideal length of chunk
    integer   :: idealLength_family     ! Time, MAFs, angle
    real(r8)  :: homeGeodAngle          ! Try to have a chunk start here
    integer   :: homeModule             ! Which key for chunk division
    real(r8)  :: overlap                ! Length for overlaps
    integer   :: overlap_family         ! Time, MAFs, angle
    logical   :: scanLLSet              ! True if scan lower limit should be used
    logical   :: scanULSet              ! True if scan upper limit should be used
    real(r8), dimension(2) :: scanLowerLimit ! Range for bottom of scan
    real(r8), dimension(2) :: scanUpperLimit ! Range for top of scan
    integer   :: criticalModules        ! Which modules must be scanning
    real(r8)  :: maxGap                 ! Length of time/MAFs/orbits allowed for gap
    integer   :: maxGap_family          ! Time, MAFs, angle
    ! Think later about bands etc.
  end type ScanDivideConfig_T

  logical :: Timing
  real :: T1

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

    ! Local Variables
    type( L1BData_T ) :: data

    character (len=480) :: MSR
    character (len=name_len) :: QUANTITY

    real(r8) :: DIFF                    ! Difference
    real(r8) :: MAXV                    ! Result of MINVAL
    real(r8) :: MINV                    ! Result of MAXVAL
    real(r8), pointer, dimension(:) :: PHI ! Array of geodAngles
    real(r8), pointer, dimension(:) :: TIME ! Array of times

    integer :: ERROR                    ! Error flag
    integer :: FIRSTMAF                 ! First MAF to consider
    integer :: FLAG                     ! For L1B reads etc.
    integer :: I                        ! Loop counter
    integer :: J                        ! Loop counter
    integer :: LASTMAF                  ! Last MAF to consider
    integer :: LOOP                     ! Major frame index
    integer :: NOMAFS                   ! Total no major frames
    integer :: NUMBAD                   ! Number bad
    integer :: NUMCHUNKS                ! Number of chunks so far
    integer :: NUMFILE                  ! Last MAF in file
    integer :: NUMOA                    ! Frame counter
    integer :: NUMSCAN                  ! Frame index of some kind

    logical :: FOUNDEND                 ! Flag
    logical :: FOUNDSTART               ! Flag

    integer, pointer, dimension(:) :: BADMAF
    integer, pointer, dimension(:) :: COUNTERMAF
    integer, pointer, dimension(:) :: FIRSTOA
    integer, pointer, dimension(:) :: FIRSTSCAN
    integer, pointer, dimension(:) :: LASTOA
    integer, pointer, dimension(:) :: LASTSCAN

    type (ScanDivideConfig_T) :: config

    ! Executable code

    if ( toggle(gen) ) call trace_begin ("ScanDivide", root )

    ! Get MLSCF values for configuration of this bit
    call ScanDivide_mlscf ( root, config )

    ! Check that input times are reasonable
    if (processingRange%startTime > processingRange%endTime) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Input start time greater than end time.')

    ! Now if the ignoreL1B flag is set we need to just make things up!
    if (config%ignoreL1B) then
      allocate ( chunks(config%noChunks), STAT=error )
      if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'chunks')

      if ( config%overlap >= config%idealLength ) call MLSMessage(MLSMSG_Error, &
        & ModuleName, "Chunks are all overlap!")

      chunks%noMAFsLowerOverlap = config%overlap
      chunks%noMAFsUpperOverlap = config%overlap
      do i = 1, size(chunks)
        chunks(i)%firstMAFIndex = (i-1)*config%idealLength+1
        chunks(i)%lastMAFIndex = i*config%idealLength
        chunks(i)%accumulatedMAFs = (i-1) * &
          & (config%idealLength - 2*config%overlap ) + 1
      enddo
      if ( toggle(gen) ) call trace_end ( "ScanDivide" )
      return
    endif

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
    if ( (processingRange%startTime >= time(noMAFs)) .or. &
      (processingRange%endTime <= time(1)) ) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
      'Input times out of range of L1 input file.')


    ! Read phi for correct module
    if (config%homeModule == l_thz) then
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
      foundStart = .false.
      do i = 2, noMAFs
        if (foundStart) exit
        if ( time(i) >= processingRange%startTime ) then
          firstMAF = i-1
          foundStart = .true.
        end if
      end do
    end if

    ! Find lastMAF for data set
    numFile = noMAFs-1

    if ( processingRange%endTime >= time(noMAFs) ) then
      lastMAF = numFile
    else
      foundEnd = .false.
      do i = 1, numFile
        if (foundEnd) exit
        if ( time(noMAFs-i) < processingRange%endTime ) then
          lastMAF = noMAFs-i-2
          foundEnd = .true.
        end if
      end do
    end if

    ! Calculate a first guess for chunk boundaries
    call ScanDivide_firstGuess(l1bInfo%L1BOAId, firstMAF, &
      & config%homeGeodAngle, lastMAF, &
      numFile, config%IdealLength, nint(config%overlap), phi, chunks, flag)
    if (flag /= 0) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_SUB // 'ScanDivide_firstGuess')
    ! If necessary, check for completely missing telemetry
    if (config%maxGap == -1.0) return

    allocate ( badMAF(noMAFs), STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Allocate // ' holder for bad MAF numbers.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    badMAF = 0
    numBad = 0
    numChunks = size(chunks)
    do i = chunks(1)%firstMAFIndex+2, chunks(numChunks)%lastMAFIndex
      select case ( config%maxGap_family )
      case ( phyq_mafs )
        diff = counterMAF(i) - counterMAF(i-1)
      case ( phyq_time )
        diff = time(i) - time(i-1)
      case default
        diff = phi(i) - phi(i-1)
      end select
      if (diff > config%maxGap) then
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
    numChunks = size(chunks)
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
        lastOa(1:numOa), config%maxGap, phi, time, config%maxGap_family, chunks)
    end if

    ! If necessary, check for valid scans
    if ( config%criticalModules /= l_none ) then
      badMAF = 0
      numBad = 0
      numChunks = size(chunks)

      ! Unless THz alone specified as critical module, read GHz first
      if ( config%criticalModules /= l_thz ) then
        quantity = 'GHz.tpGeodAlt'
        call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
        if (flag /= 0) then
          msr = MLSMSG_L1BRead // quantity
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if
        ! Check GHz for bad scans
        do i = chunks(1)%firstMAFIndex+1, chunks(numChunks)%lastMAFIndex+1
          minv = minval(data%dpField(1,:,i))
          maxv = maxval(data%dpField(1,:,i))
          if ( (minv < config%scanLowerLimit(1)) .or. &
               (minv > config%scanLowerLimit(2)) .or. &
               (maxv < config%scanUpperLimit(1)) .or. &
               (maxv > config%scanUpperLimit(2)) ) then
            badMAF(i) = 1
            numBad = numBad + 1
          end if
        end do

        call DeallocateL1BData(data, flag)
        if (flag /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
          MLSMSG_DeallocateL1b)
      end if

      ! If necessary, read THz data

      if ( (config%criticalModules == l_thz) .or. (config%criticalModules == l_both) .or. &
        ((config%criticalModules == l_either) .and. (numBad > 0)) ) then
        quantity = 'THz.tpGeodAlt'
        call ReadL1BData (l1bInfo%L1BOAId, quantity, data, noMAFs, flag)
        if (flag /= 0) then
          msr = MLSMSG_L1BRead // quantity
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
        end if

        ! If EITHER is critical, check only bad GHz scans; if THz okay,
        ! unmark MAF as bad.
        if ( config%criticalModules == l_either ) then
          do i = chunks(1)%firstMAFIndex+1,chunks(numChunks)%lastMAFIndex+1
            if (badMAF(i) == 1) then
              minv = minval(data%dpField(1,:,i))
              maxv = maxval(data%dpField(1,:,i))
              if ( (minv >= config%scanLowerLimit(1)) .and. &
                &  (minv <= config%scanLowerLimit(2)) .and. &
                &  (maxv >= config%scanUpperLimit(1)) .and. &
                &  (maxv <= config%scanUpperLimit(2)) ) then
                badMAF(i) = 0
                numBad = numBad - 1
              end if
            end if
          end do
          ! If THz is critical, check entire data set for bad THz scans
        else
          do i = chunks(1)%firstMAFIndex+1,chunks(numChunks)%lastMAFIndex+1
            minv = minval(data%dpField(1,:,i))
            maxv = maxval(data%dpField(1,:,i))
            if ( (minv < config%scanLowerLimit(1)) .and. &
              &  (minv > config%scanLowerLimit(2)) .and. &
              &  (maxv < config%scanUpperLimit(1)) .and. &
              &  (maxv > config%scanUpperLimit(2)) ) then
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
          else if ( (i > 1) .and. (lastScan(numScan) == i-2) ) then
            numScan = numScan + 1
          end if
        end do

        ! Check whether these (including MAF "walls" previously found) exceed MaxGap;
        ! re-organize chunks, if necessary
        call ScanDivide_checkReorganizeRange(counterMAF, &
          firstScan(1:loop), lastScan(1:loop), &
          config%maxGap, phi, time, config%maxGap_family, chunks)
        deallocate ( firstScan, lastScan, STAT=error )
        if (error /= 0) then
          msr = MLSMSG_Deallocate // ' first & lastScan local variables.'
          call MLSMessage(MLSMSG_Warning, ModuleName, msr)
        end if
      end if
    end if
    ! If necessary, check for critical bands in the L1bRad file
    numChunks = size(chunks)
    !if ( bands /= 0 ) then
    !end if

    ! Deallocations
    deallocate ( badMAF, counterMAF, firstOa, lastOa, phi, time, STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Deallocate // ' local variables.'
      call MLSMessage(MLSMSG_Warning, ModuleName, msr)
    end if
    if ( timing ) call sayTime
    if ( toggle(gen) ) call trace_end ( "ScanDivide" )

  end subroutine ScanAndDivide

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

    flag = 0        ! Assume no errors occur -- but make sure FLAG is defined!

    ! Find max & min phi for data set

    max = maxval( phi(firstMAF+1:lastMAF+1) )
    min = minval( phi(firstMAF+1:lastMAF+1) )

    ! Convert from "fraction of orbit" to "degrees of phi"

    phiLen = orbLen

    ! Find phi boundaries above & below HomeGeodAngle

    bdryMax = anint(max_orbits * 360.0 / orbLen)

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

    call Hunt(phi, bdryPhi, indices, allowTopValue=.true.)

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
    ! Think there is a problem here, get 1 when would expect 0, NJL.

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

    numBad = size(bad)

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

    numChunks = size(chunks)

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

    numChunks = size(chunks)

    allocate ( newChunks(numChunks+1), STAT=error )
    if (error /= 0) then
      msr = MLSMSG_Allocate // ' newChunks array.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    end if

    ! Find the chunk where the gap begins

    newInd = 1
    chunkInd = 1
    do i = 1, numChunks

      if ( (firstBadMAFIndex <= chunks(chunkInd)%lastMAFIndex) .and. &
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

    integer :: FINISH, FLAG, I, NUMCHUNKS, NUMGAPS, START

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
        finish = chunks(numChunks)%lastMAFIndex + 1
      else
        finish = lastBad(i) + 2
      end if

      select case ( units )
      case ( phyq_mafs )
        diff = counterMAF(finish) - counterMAF(start)
      case ( phyq_time )
        diff = time(finish) - time(start)
      case default
        diff = phi(finish) - phi(start)
      end select

      ! If diff exceeds MaxGap, re-organize chunks

      if (diff > maxGap) &
        call ScanDivide_range(chunks, firstBad(i), lastBad(i), flag)

    end do

  end subroutine ScanDivide_checkReorganizeRange
!------------------------------------------------

!----------------------------------------------  ScanDivide_mlscf  -----
subroutine ScanDivide_mlscf ( root, config )
  !-----------------------------------------------------------------------

  ! Brief description of subroutine
  ! This subroutine identifies, separates, and checks values from the section
  ! of the MLSCF (ChunkDivide) passed to Scan/Divide.

  ! Arguments

  integer, intent(in) :: ROOT    ! Root of the ChunkDivide section of the
  ! MLSCF abstract syntax tree
  type (ScanDivideConfig_T), intent(out) :: CONFIG ! Result of operation

  ! Parameters

  ! For announce_error:
  integer, parameter :: BadUnitsForIgnore = 1
  integer, parameter :: NotLength = BadUnitsForIgnore + 1
  integer, parameter :: NotSpecified = notLength + 1
  integer, parameter :: NotAngle = NotSpecified + 1

  ! Functions

  ! Variables
  integer :: Error     ! Error level
  logical :: GOT(first_parm:last_parm) = .false.
  integer :: I         ! Loop inductor
  integer :: KEY       ! A P_... parameter from Init_Tables_Module
  integer :: SON       ! A son of the ChunkDivide section node
  integer :: UNITS(2)  ! Units of expression
  real(r8) :: VALUE(2)   ! Value of expression

  ! Initialize variables to 'unfound' values

  config%ignoreL1B = .false.
  config%noChunks = 0
  config%idealLength_family = PHYQ_Invalid
  config%homeGeodAngle = 0.0
  config%homeModule = L_None
  config%overlap = 0.0
  config%overlap_family = PHYQ_Invalid
  config%scanLLSet = .false.
  config%scanULSet = .false.
  config%criticalModules = l_none
  config%maxGap_family = PHYQ_Invalid

  error = 0

  ! Loop through MLSCF section, identifying keywords and setting variables

  do i = 2, nsons(root)-1 ! Skip the section identifiers
    son = subtree(i,root)
    if ( node_id(son) == n_equal ) then
      key = decoration(subtree(1,son)) ! P_... index from Init_Tables_Module
      got(key) = .true.
      call expr ( subtree(2,son), units, value )

      select case ( key )
      case ( p_ignoreL1B )
        config%ignoreL1B = value(1) == l_true
      case ( p_noChunks )
        config%noChunks = nint(value(1))
      case ( p_ideal_length )
        config%idealLength = value(1)
        config%idealLength_family = units(1)
      case ( p_home_geod_angle )
        config%homeGeodAngle = value(1)
        if ( units(1) /= PHYQ_Angle ) call announce_error ( son, notAngle, key )
      case ( p_home_module )
        config%homeModule = value(1)
      case ( p_overlap )
        config%overlap = value(1)
        config%overlap_family = units(1)
      case ( p_scan_lower_limit )
        if ( units(1) /= phyq_Length) &
          & call announce_error ( son, notLength, key )
        config%scanLLSet = .true.
        config%scanLowerLimit = value
      case ( p_scan_upper_limit )
        if ( units(1) /= phyq_Length) &
          & call announce_error ( son, notLength, key )
        config%scanULSet = .true.
        config%scanUpperLimit = value
      case ( p_critical_scanning_modules )
        config%criticalModules = value(1)
      case ( p_max_gap )
        config%maxGap = value(1)
        config%maxGap_family = units(1)
      case ( p_critical_bands )
        !bands = value(1)
      end select
    else
      if ( node_id(son) == n_named ) son = subtree(2,son) ! ignore label
      select case ( get_spec_id(son) )
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      case default
      end select
    end if

  end do

  ! Check for missing values that would cause error exits

  do i = first_parm, last_parm
    select case ( i )
    case ( p_ideal_length, p_overlap )
      if ( .not. got(i) ) call announce_error ( root, notSpecified, i )
    case ( p_home_module, p_critical_scanning_modules )    
      if ( .not. got(i) .and. .not. config%ignoreL1B) &
        & call announce_error ( root, notSpecified, i )
    case ( p_noChunks) 
      if ( .not. got(i) .and. config%ignoreL1B) &
        & call announce_error ( root, notSpecified, i )
    end select
  end do

  if (config%ignoreL1B) then
    if (config%idealLength_family /= phyq_MAFs ) &
      & call announce_error ( root, badUnitsForIgnore, p_ideal_length )
    if (config%overlap_family /= phyq_MAFs ) &
      & call announce_error ( root, badUnitsForIgnore, p_overlap )
  endif    

  ! Check for values that depend on the presence of other quantities

  if ( config%criticalModules /= l_none ) then
    if ( .not. got(p_scan_lower_limit) ) &
      & call announce_error ( root, notSpecified, p_scan_lower_limit )
    if ( .not. got(p_scan_upper_limit) ) &
      & call announce_error ( root, notSpecified, p_scan_upper_limit )
  end if

  if ( error > 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
    & 'Errors in specification prevent processing.' )

contains

  subroutine Announce_Error ( where, Code, Param )
    integer, intent(in) :: where, Code, Param

    error = max(error,1)
    call print_source ( source_ref(where) )
    call output ( ' ScanDivide complained: ' )
    select case ( code )
    case ( BadUnitsForIgnore )
      call output ( ' In the IgnoreL1B case ' )
      call display_string ( parm_indices(param) )
      call output (' must have units of MAFs' )
    case ( notLength )
      call output ( ' The value of ' )
      call display_string ( parm_indices(param) )
      call output ( ' does not have units of Length.', advance='yes' )
    case ( notAngle )
      call output ( ' The value of ' )
      call display_string ( parm_indices(param) )
      call output ( ' does not have units of Angle.', advance='yes' )
    case ( notSpecified )
      call output ( ' The parameter ' )
      call display_string ( parm_indices(param) )
      call output ( ' is required but not specified.', advance='yes' )
    end select
  end subroutine Announce_Error
!--------------------------------
end subroutine ScanDivide_mlscf
!--------------------------------

  ! ....................................................  SayTime  .....
  subroutine SayTime
    real :: T2
    call cpu_time ( t2 )
    call output ( 'Timing for ScanDivide = ' )
    call output ( dble(t2-t1), advance='yes' )
    timing = .false.
  end subroutine SayTime

!====================
end module ScanDivide
!====================

!# $Log$
!# Revision 2.12  2001/04/26 02:44:17  vsnyder
!# Moved *_indices declarations from init_tables_module to intrinsic
!#
!# Revision 2.11  2001/04/24 23:12:12  vsnyder
!# Add timing
!#
!# Revision 2.10  2001/04/24 19:08:58  livesey
!# Interim version, needs changes later on.
!#
!# Revision 2.6  2001/03/08 23:36:19  vsnyder
!# Make sure FLAG is defined in ScanDivide_firstGuess
!#
!# Revision 2.5  2001/03/02 19:32:06  pwagner
!# Gets MLSMSG_L1BRead from MLSMessageModule, not MLSStrings
!#
!# Revision 2.4  2001/02/13 00:09:02  vsnyder
!# Simplify and improve MLSCF error messages
!#
!# Revision 2.3  2001/02/12 20:29:34  livesey
!# Flagged a possible error region
!#
!# Revision 2.2  2000/09/11 20:00:20  ahanzel
!# Removed old log entries in file.
!#
!# Revision 2.1  2000/09/08 22:55:56  vsnyder
!# Revised to use the tree output by the parser
!#
