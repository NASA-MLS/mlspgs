! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ChunkDivide_m

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use EXPR_M, only: EXPR
  use Dumper, only: DUMP
  use Dump_0, only: DUMP
  use MLSCommon, only: R8, RP, L1BINFO_T, MLSCHUNK_T, TAI93_Range_T, FINDFIRST
  use MLSNumerics, only: Hunt
  use Intrinsic, only: L_NONE, FIELD_INDICES, LIT_INDICES
  use Tree, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUBTREE
  use Lexer_core, only: PRINT_SOURCE
  use Tree_types, only: N_EQUAL, N_NAMED
  use Init_Tables_Module, only: F_OVERLAP, F_MAXLENGTH, F_NOCHUNKS, &
    & F_METHOD, F_HOMEMODULE, F_CRITICALMODULES, F_HOMEGEODANGLE, F_SCANLOWERLIMIT, &
    & F_SCANUPPERLIMIT, F_NOSLAVES, FIELD_FIRST, FIELD_LAST, L_EVEN, &
    & L_FIXED, F_MAXGAP, L_ORBITAL, S_CHUNKDIVIDE, L_BOTH, L_EITHER
  use L1BData, only: DEALLOCATEL1BDATA, L1BDATA_T, NAME_LEN, READL1BDATA
  use Units, only: PHYQ_INVALID, PHYQ_LENGTH, PHYQ_MAFS, PHYQ_TIME, PHYQ_ANGLE, &
    & PHYQ_LENGTH, PHYQ_DIMENSIONLESS
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, MLSMSG_WARNING
  use MLSSignals_m, only: MODULES, ISMODULESPACECRAFT
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
  use Output_M, only: BLANKS, OUTPUT
  use String_table, only: GET_STRING, DISPLAY_STRING
  use Time_M, only: Time_Now
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END

  ! This module replaces the old ScanDivide, and is a new approach to dividing
  ! the data into chunks.

  implicit none
  private

  public :: DestroyChunkDatabase, ChunkDivide

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This type is filled by the l2cf and describes the configuration of the
  ! chunk division process.
  type ChunkDivideConfig_T
    integer :: method = l_none          ! See below.
    real(rp) :: maxLength = 0           ! Maximum length of chunks
    integer :: maxLengthFamily = PHYQ_Invalid ! PHYQ_Angle etc.
    integer :: noChunks = 0             ! Number of chunks for [fixed]
    real(rp) :: overlap = 0.0           ! Desired length of overlaps
    integer :: overlapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
    integer :: noSlaves = 0             ! Number of slave nodes [even]
    integer :: homeModule = l_none      ! Which module to consider [orbital]
    real(rp) :: homeGeodAngle = 0.0     ! Aim for one chunk to start here [orbital]
    logical   :: scanLLSet = .false.    ! True if scan lower limit should be used
    logical   :: scanULSet = .false.    ! True if scan upper limit should be used
    real(rp), dimension(2) :: scanLowerLimit ! Range for bottom of scan
    real(rp), dimension(2) :: scanUpperLimit ! Range for top of scan
    integer   :: criticalModules = l_none ! Which modules must be scanning
    real(rp)  :: maxGap = 0.0           ! Length of time/MAFs/orbits allowed for gap
    integer   :: maxGapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
  end type ChunkDivideConfig_T
  ! The chunk divide methods are:
  !
  ! Fixed - Ignore the L1B file, just give a fixed set of chunks as described
  !
  ! Orbital - Chunks are some ideal fraction of an orbit.  The algorithm
  !           desires to keep the chunk boundaries at the same positions each
  !           orbit where possible.
  !
  ! Even - Hope to have chunks all about the same length, as quoted. 

  ! This type describes obstructions in the Level 1 data which will affect the
  ! selection of chunk divisions.
  type Obstruction_T
    logical :: range                    ! If set is a range not a wall
    integer, dimension(2) :: MAFS       ! Affected MAF or MAF range
  end type Obstruction_T

  logical :: Timing
  real :: T1

  interface dump
    module procedure DUMP_OBSTRUCTIONS
  end interface

contains ! =================================== Public Procedures==============

  !------------------------------------------  DestroyChunkDatabase  -----
  subroutine DestroyChunkDatabase ( chunks )
    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS
    integer :: STATUS ! From deallocate

    deallocate ( chunks, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_DeAllocate // "Chunks" )
  end subroutine DestroyChunkDatabase

  ! ----------------------------------------  Chunk Divide --------------
  subroutine ChunkDivide ( root, processingRange, l1bInfo, chunks )
    integer, intent(in) :: ROOT    ! Root of the L2CF tree for ChunkDivide
    type( L1BInfo_T ), intent(in) :: L1BINFO
    type( TAI93_Range_T ), intent(in) :: PROCESSINGRANGE
    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS

    ! Local variables
    type (ChunkDivideConfig_T) :: CONFIG ! Configuration
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    integer :: STATUS                   ! From deallocate
    integer, dimension(2) :: MAFRANGE   ! Processing range in MAFs

    ! Executable code
    
    if ( toggle(gen) ) call trace_begin ("ScanDivide", root )

    timing = section_times
    if ( timing ) call time_now ( t1 )

    ! First decode the l2cf instructions
    call ChunkDivideL2CF ( root, config )

    ! For methods other than fixed, we want to survey the L1 data and note the
    ! location of obstructions
    nullify ( obstructions )
    if ( config%method /= l_fixed ) &
      & call SurveyL1BData ( processingRange, l1bInfo, config, mafRange,&
      & obstructions )

    ! Now go place the chunks.
    select case ( config%method )
    case ( l_fixed )
      call ChunkDivide_Fixed ( config, chunks )
    case ( l_orbital )
      call ChunkDivide_Orbital ( config, mafRange, l1bInfo, &
        & obstructions, chunks )
    case ( l_even )
      call ChunkDivide_Even ( config, mafRange, l1bInfo, &
        & obstructions, chunks )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unexpected problem with ChunkDivide' )
    end select

    ! Tidy up
    if ( associated(obstructions) ) then
      deallocate ( obstructions, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'obstructions' )
    end if

  end subroutine ChunkDivide

  ! ============================== Private Procedures ====================

  !---------------------------------------- Add obstruction to database --
  subroutine AddObstructionToDatabase ( database, item )

    ! Dummy arguments
    type (Obstruction_T), dimension(:), pointer :: DATABASE
    type (Obstruction_T), intent(in) :: ITEM

    ! Local variables
    type (Obstruction_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"

  end subroutine AddObstructionToDatabase

  !---------------------------------------- Add chunk to database --
  subroutine AddChunkToDatabase ( database, item )

    ! Dummy arguments
    type (MLSChunk_T), dimension(:), pointer :: DATABASE
    type (MLSChunk_T), intent(in) :: ITEM

    ! Local variables
    type (MLSChunk_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"

  end subroutine AddChunkToDatabase

  !----------------------------------------- ChunkDivide_Even --------
  subroutine ChunkDivide_Even ( config, mafRange, l1bInfo, &
    & obstructions, chunks )
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    integer, dimension(2), intent(in) :: MAFRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local variables
    
    ! Exectuable code
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "I'm sorry, we haven't written the code for even chunk divides yet" )
  end subroutine ChunkDivide_Even

  !----------------------------------------- ChunkDivide_Even --------
  subroutine ChunkDivide_Fixed ( config, chunks )
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local variables
    integer :: I                        ! Loop inductor
    integer :: STATUS                   ! From allocate
    integer :: MAXLENGTH                ! nint(config%maxLength)
    integer :: OVERLAP                  ! nint(config%overlap)
    integer :: NONONOVERLAP             ! maxLength-2*overlap
    
    ! Exectuable code
    allocate ( chunks(config%noChunks), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'chunks' )

    maxLength = nint ( config%maxLength )
    overlap = nint ( config%overlap )
    noNonOverlap = maxLength - 2 * overlap
    do i = 1, config%noChunks
      chunks(i)%firstMAFIndex = max ( (i-1)*maxLength - overlap, 0 )
      chunks(i)%lastMAFIndex = i*maxLength + overlap - 1
      chunks(i)%noMAFsUpperOverlap = overlap
      if ( i > 1 ) then
        chunks(i)%noMAFsLowerOverlap = overlap
      else
        chunks(i)%noMAFsLowerOverlap = 0
      end if
      chunks(i)%accumulatedMAFs = chunks(i)%firstMAFIndex
    end do

  end subroutine ChunkDivide_Fixed

  !----------------------------------------- ChunkDivide_Orbital --------
  subroutine ChunkDivide_Orbital ( config, mafRange, l1bInfo, &
    & obstructions, chunks )
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    integer, dimension(2), intent(in) :: MAFRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local parameters
    real(r8), parameter :: HOMEACCURACY = 3.0 ! Try to hit homeGeodAngle within this

    ! Local variables
    type (L1BData_T) :: TAITIME         ! From L1BOA
    type (L1BData_T) :: TPGEODANGLE     ! From L1BOA

    character(len=10) :: MODNAMESTR     ! Home module name as string

    integer :: CHUNK                    ! Loop counter
    integer :: FLAG                     ! From ReadL1B
    integer :: HOMEMAF                  ! first MAF after homeGeodAngle
    integer :: HOME                     ! Index of home MAF in array
    integer :: M1, M2                   ! MafRange + 1
    integer :: NOCHUNKSBELOWHOME        ! Used for placing chunks
    integer :: NOMAFSATORABOVEHOME      ! Fairly self descriptive
    integer :: NOMAFS                   ! Number of MAFs to consider
    integer :: NOMAFSREAD               ! From ReadL1B
    integer :: ORBIT                    ! Used to locate homeMAF
    integer :: STATUS                   ! From allocate etc.
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: OVERLAPS                 ! Overlaps as integer (MAFs)
    integer :: MAXLENGTH                ! Max length as integer (MAFs)

    integer, dimension(:), pointer :: NEWFIRSTMAFS ! For thinking about overlaps
    integer, dimension(:), pointer :: NEWLASTMAFS ! For thinking about overlaps

    real(r8) :: ANGLEINCREMENT          ! Increment in hunt for homeMAF
    real(r8) :: MAXANGLE                ! Of range in data
    real(r8) :: MAXTIME                 ! Time range in data
    real(r8) :: MAXV                    ! Either minTime or minAngle
    real(r8) :: MINANGLE                ! Of range in data
    real(r8) :: MINTIME                 ! Time range in data
    real(r8) :: MINV                    ! Either minTime or minAngle
    real(r8) :: TESTANGLE               ! Angle to check for
    real(r8) :: HOMEV                   ! Value of angle/time at home

    real(r8), dimension(:), pointer :: BOUNDARIES ! Used in placing chunks
    real(r8), dimension(:), pointer :: FIELD ! Used in placing chunks

    ! Exectuable code
    
    ! Read in the data we're going to need
    call get_string ( lit_indices(config%homeModule), modNameStr, strip=.true. )
    call ReadL1BData ( l1bInfo%l1bOAId, trim(modNameStr)//'.tpGeodAngle', &
      & tpGeodAngle, noMAFsRead, flag )
    call ReadL1BData ( l1bInfo%l1bOAId, 'MAFStartTimeTAI', &
      & taiTime, noMAFsRead, flag )
    noMAFs = mafRange(2) - mafRange(1) + 1
    m1 = mafRange(1) + 1
    m2 = mafRange(2) + 1

    minAngle = minval ( tpGeodAngle%dpField(1,1,m1:m2) )
    maxAngle = maxval ( tpGeodAngle%dpField(1,1,m1:m2) )
    minTime = minval ( taiTime%dpField(1,1,m1:m2) )
    maxTime = maxval ( taiTime%dpField(1,1,m1:m2) )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! First try to locate the last MAF before the homeGeodAngle
    orbit = int ( tpGeodAngle%dpField(1,1,m1)/360.0 )
    if ( tpGeodAngle%dpField(1,1,mafRange(1)) < 0.0 ) orbit = orbit - 1
    testAngle = config%homeGeodAngle + orbit*360.0
    if ( config%maxLengthFamily == PHYQ_Angle ) then
      angleIncrement = config%maxLength
    else
      angleIncrement = 360.0
    end if

    homeHuntLoop: do
      if ( testAngle < minAngle ) then
        testAngle = testAngle + angleIncrement
        cycle
      endif
      if ( testAngle > maxAngle ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to establish a home major frame, using the first' )
        home = 1
        exit homeHuntLoop
      end if
      ! Find MAF which starts before this test angle
      call Hunt ( tpGeodAngle%dpField(1,1,:), testAngle, home, nearest=.true.,&
        & allowTopValue = .true. )
      ! Now if this is close enough, accept it
      if ( abs ( tpGeodAngle%dpField(1,1,home) - &
        & testAngle ) < HomeAccuracy ) exit homeHuntLoop
      ! Otherwise, keep looking
      testAngle = testAngle + angleIncrement
    end do homeHuntLoop
    homeMAF = home + m1 - 1

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! OK, now we have a home MAF, get a first cut for the chunks
    ! We work out the chunk ends for each chunk according to how the
    ! maxLength field is specified.
    if ( config%maxLengthFamily == PHYQ_MAFs ) then
      maxLength = nint ( config%maxLength )
      noChunks = home / maxLength
      if ( mod ( home, maxLength ) /= 0 ) noChunks = noChunks + 1
      noChunksBelowHome = noChunks
      noMAFsAtOrAboveHome = noMAFs - home + 1
      noChunks = noChunks + noMAFsAtOrAboveHome / maxLength
      if ( mod ( noMAFsAtOrAboveHome, maxLength ) /= 0 ) &
        & noChunks = noChunks + 1
      
      ! Allocate the chunks
      allocate ( chunks(noChunks), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'chunks (maxLength/MAFs)' )
      
      ! Work out their positions
      do chunk = 1, noChunks
        chunks(chunk)%lastMAFIndex = homeMAF + &
          & ( chunk - noChunksBelowHome ) * maxLength
      end do
    else
      ! For angle and time, they are similar enough we'll just do some stuff
      ! with pointers to allow us to use common code to sort them out
      select case ( config%maxLengthFamily )
      case ( PHYQ_Angle )
        field => tpGeodAngle%dpField(1,1,:)
        minV = minAngle
        maxV = maxAngle
      case ( PHYQ_Time )
        field => taiTime%dpField(1,1,:)
        minV = minTime
        maxV = maxTime
      case ( PHYQ_MAFs)
      end select
      homeV = field(home)

      noChunks = int ( ( homeV - minV ) / config%maxLength )
      if ( homeV > minV ) noChunks = noChunks + 1
      noChunksBelowHome = noChunks
      noChunks = noChunks + int ( ( maxV - homeV ) / config%maxLength )
      if ( homeV + config%maxLength * ( noChunks - noChunksBelowHome ) < maxV ) &
        & noChunks = noChunks + 1

      ! Allocate the chunks
      allocate ( chunks(noChunks), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'chunks (maxLength/time/angle)' )

      ! Work out their positions
      ! Boundaries are the angles/times at the end of the chunks
      nullify ( boundaries )
      call Allocate_test ( boundaries, noChunks, 'boundaries', ModuleName )
      do chunk = 1, noChunks
        boundaries(chunk) = homeV + ( chunk - noChunksBelowHome ) * config%maxLength
      end do
      boundaries = min ( boundaries, maxV )
      boundaries = max ( boundaries, minV )
      call Hunt ( field, boundaries, chunks%lastMAFIndex, &
        & allowTopValue=.true., nearest=.true. )
      call Deallocate_test ( boundaries, 'boundaries', ModuleName )
    end if

    ! Now deduce the chunk starts from the ends of their predecessors
    chunks(2:noChunks)%firstMAFIndex = chunks(1:noChunks-1)%lastMAFIndex + 1
    chunks(1)%firstMAFIndex = 1
    
    ! Now offset these to the index in the file not the array
    chunks%firstMAFIndex = chunks%firstMAFIndex + mafRange(1) - 1
    chunks%lastMAFIndex = chunks%lastMAFIndex + mafRange(1) - 1 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Think about overlaps
    nullify ( newFirstMAFs, newLastMAFs )
    call Allocate_test ( newFirstMAFs, noChunks, 'newFirstMAFs', ModuleName )
    call Allocate_test ( newLastMAFs, noChunks, 'newLastMAFs', ModuleName )
    if ( config%overlapFamily == PHYQ_MAFs ) then
      newFirstMAFs = max(chunks%firstMAFIndex - nint(config%overlap), mafRange(1) )
      newLastMAFs = min(chunks%lastMAFIndex + nint(config%overlap), mafRange(2) )
    else
      ! For angle and time, they are similar enough we'll just do some stuff
      ! with pointers to allow us to use common code to sort them out
      ! Note that before we search only over the range in mafRange, here we
      ! search over the whole file as we can spill over the processing range
      ! for overlaps.
      select case ( config%maxLengthFamily )
      case ( PHYQ_Angle )
        field => tpGeodAngle%dpField(1,1,:)
        minV = minAngle
      case ( PHYQ_Time )
        field => taiTime%dpField(1,1,:)
        minV = minTime
      case ( PHYQ_MAFs)
      end select
      call Hunt ( field, field(chunks%firstMAFIndex-1) + &
        & config%overlap, newFirstMAFs, allowTopValue=.true. )
      call Hunt ( field, field(chunks%firstMAFIndex-1) - &
        & config%overlap, newLastMAFs, allowTopValue=.true. )
      ! Correct this to be real MAF indicies (starting from zero)
      newFirstMAFs = newFirstMAFs - 1
      newLastMAFs = newLastMAFs - 1
    end if
    chunks%noMAFsLowerOverlap = chunks%firstMAFIndex - newFirstMAFs
    chunks%noMAFsUpperOverlap = newLastMAFs - chunks%lastMAFIndex
    chunks%firstMAFIndex = newFirstMAFs
    chunks%lastMAFIndex = newLastMAFs
    call Deallocate_test ( newFirstMAFs, 'newFirstMAFs', ModuleName )
    call Deallocate_test ( newLastMAFs, 'newLastMAFs', ModuleName )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Now think about the obstructions

    call DealWithObstructions ( chunks, obstructions )

    ! Tidy up
    call DeallocateL1BData ( tpGeodAngle )
    call DeallocateL1BData ( taiTime )

  end subroutine ChunkDivide_Orbital

  !---------------------------------------------- ChunkDivideL2CF -----
  subroutine ChunkDivideL2CF ( sectionRoot, config )
    ! This subroutine identifies, separates, and checks values from the section
    ! of the MLSCF (ChunkDivide) passed to Scan/Divide.
    integer, intent(in) :: SECTIONROOT    ! Root of the ChunkDivide section of the
    ! MLSCF abstract syntax tree
    type (ChunkDivideConfig_T), intent(out) :: CONFIG ! Result of operation

    ! For announce_error:
    integer, parameter :: BadUnits = 1
    integer, parameter :: NotSpecified = BadUnits + 1
    integer, parameter :: Unnecessary = NotSpecified + 1

    integer, target, dimension(3) :: NeededForFixed = &
      & (/ f_noChunks, f_maxLength, f_overlap /)
    integer, target, dimension(7) :: NotWantedForFixed = &
      & (/ f_noSlaves, f_homeModule, f_homeGeodAngle, f_scanLowerLimit, &
      &    f_scanUpperLimit, f_criticalModules, f_maxGap /)

    integer, target, dimension(6) :: NeededForOrbital = &
      & (/ f_maxLength, f_overlap, f_homeModule, f_homeGeodAngle, &
      &    f_criticalModules, f_maxGap /)
    integer, target, dimension(2) :: NotWantedForOrbital = &
      & (/ f_noChunks, f_noSlaves /)

    integer, target, dimension(6) :: NeededForEven = &
      & (/ f_maxLength, f_overlap, f_maxLength, f_noSlaves, f_maxGap, &
      &    f_criticalModules /)
    integer, target, dimension(3) :: NotWantedForEven = &
      & (/ f_noChunks, f_homeModule, f_homeGeodAngle /)

    ! Local variables
    integer :: ERROR                    ! Error level
    integer :: FIELDINDEX               ! Tree type
    integer :: GSON                     ! A son of son the ChunkDivide section node
    integer :: I                        ! Loop inductor
    integer :: ROOT                     ! Root of ChunkDivide command
    integer :: SON                      ! A son of the ChunkDivide section node
    integer :: UNITS(2)                 ! Units of expression
    integer, dimension(:), pointer :: NEEDED ! Which fields are needed
    integer, dimension(:), pointer :: NOTWANTED ! Which fields are not wanted
    logical :: GOT(field_first:field_last) = .false.
    real(rp) :: VALUE(2)                ! Value of expression

    ! Executable code
    error = 0

    ! Eventually the ChunkDivide command will be free floating, in the meantime
    ! find it within the section
    ! WE CAN GET RID OF THIS BIT WHEN THE COMMAND FLOATS FREE LATER
    do i = 2, nsons(sectionRoot)-1      ! Skip the begin/end section
      root = subtree(i,sectionRoot)
      if ( node_id(root) /= n_named ) cycle
      if ( get_spec_id(root) == s_chunkDivide ) exit
    end do

    got = .false.

    ! Loop through the command identifying parameters.
    do i = 2, nsons(root) ! Skip the command
      son = subtree(i,root)
      fieldIndex = get_field_id(son)
      got(fieldIndex) = .true.
      if (nsons(son) > 1 ) then
        gson = subtree(2,son)
        call expr ( gson, units, value )
      end if
      ! Get value for this field if appropriate
      select case ( fieldIndex )
      case ( f_method )
        config%method = decoration ( gson )
      case ( f_noChunks )
        config%noChunks = nint ( value(1) )
      case ( f_maxLength )
        config%maxLength = value(1)
        config%maxLengthFamily = units(1)
      case ( f_overlap )
        config%overlap = value(1)
        config%overlapFamily = units(1)
      case ( f_noSlaves )
        config%noSlaves = value(1)
        if ( units(1) /= PHYQ_DimensionLess ) &
          & call AnnounceError ( root, BadUnits, fieldIndex )
      case ( f_homeModule )
        config%homeModule = decoration ( gson )
      case ( f_homeGeodAngle )
        config%homeGeodAngle = value(1)
        if ( units(1) /= PHYQ_Angle ) &
          & call AnnounceError ( root, BadUnits, fieldIndex )
      case ( f_scanLowerLimit )
        if ( any ( units /= PHYQ_Dimensionless .and. units /= PHYQ_Length ) &
          & .or. .not. any ( units == PHYQ_Length ) ) &
          & call AnnounceError ( root, BadUnits, fieldIndex )
        config%scanLowerLimit = value
        config%scanLLSet = .true.
      case ( f_scanUpperLimit )
        if ( any ( units /= PHYQ_Dimensionless .and. units /= PHYQ_Length ) &
          & .or. .not. any ( units == PHYQ_Length ) ) &
          & call AnnounceError ( root, BadUnits, fieldIndex )
        config%scanUpperLimit = value
        config%scanULSet = .true.
      case ( f_criticalModules )
        config%criticalModules = decoration ( gson )
      case ( f_maxGap )
        config%maxGap = value(1)
        config%maxGapFamily = units(1)
      end select
    end do

    ! Now check the sanity of what we've been given, this varies a bit
    ! depending on the method
    select case ( config%method )
    case ( l_fixed )
      needed => NeededForFixed
      notWanted => NotWantedForFixed
    case ( l_orbital )
      needed => NeededForOrbital
      notWanted => NotWantedForOrbital
    case ( l_even )
      needed => NeededForEven
      notWanted => NotWantedForEven
    end select

    ! Check we've got all the arguments we need
    do i = 1, size(needed)
      if ( .not. got(needed(i) ) ) &
        & call AnnounceError ( root, notSpecified, needed(i) )
    end do
    ! Check we don't have unnecessary ones
    do i = 1, size(notWanted)
      if ( got(notWanted(i) ) ) &
        & call AnnounceError ( root, unnecessary, notWanted(i) )
    end do

    ! Make other checks of parameters
    if ( config%criticalModules /= l_none ) then
      if ( .not. got(f_scanLowerLimit) ) &
        & call AnnounceError ( root, notSpecified, f_scanLowerLimit )
      if ( .not. got (f_scanUpperLimit) ) &
        & call AnnounceError ( root, notSpecified, f_scanUpperLimit )
    else
      if ( got(f_scanLowerLimit) ) &
        & call AnnounceError ( root, unnecessary, f_scanLowerLimit )
      if ( got (f_scanUpperLimit) ) &
        & call AnnounceError ( root, unnecessary, f_scanUpperLimit )
    end if

    ! Now check the units for various cases
    if ( got(f_maxgap) .and. all ( config%maxGapFamily /= &
      & (/ PHYQ_MAFs, PHYQ_Angle, PHYQ_Time /))) &
      & call AnnounceError ( root, badUnits, f_maxGap )
    if ( config%method == l_orbital ) then
      if (all(config%maxLengthFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
        & call AnnounceError ( root, badUnits, f_maxLength )
      if (all(config%overlapFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
        & call AnnounceError ( root, badUnits, f_overlap )
    else
      if ( config%maxLengthFamily /= PHYQ_MAFs ) &
        & call AnnounceError ( root, badUnits, f_maxLength )
      if ( config%overlapFamily /= PHYQ_MAFs ) &
        & call AnnounceError ( root, badUnits, f_overlap )
    end if

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Problem with ChunkDivide command' )

    ! That's it, we're all valid now.

  contains ! - - - - - - - - - - - - - - - - - - - -

    subroutine AnnounceError ( where, Code, field )
      integer, intent(in) :: where, Code, field

      error = max(error,1)
      call print_source ( source_ref(where) )
      call output ( ' ChunkDivide complained: ' )
      select case ( code )
      case ( BadUnits )
        call output ( ' The field ' )
        call display_string ( field_indices(field) )
        call output (' has inappropriate units', advance='yes' )
      case ( notSpecified )
        call output ( ' The parameter ' )
        call display_string ( field_indices(field) )
        call output ( ' is required but not specified.', advance='yes' )
      case ( unnecessary )
        call output ( ' The parameter ' )
        call display_string ( field_indices(field) )
        call output ( ' is specified but not appropriate.', advance='yes' )
      end select
    end subroutine AnnounceError

  end subroutine ChunkDivideL2CF

  ! ------------------------------------------ ConvertFlagsToObstructions --
  subroutine ConvertFlagsToObstructions ( valid, mafRange, obstructions )
    ! This routine takes an array of logicals indicating good/bad data
    ! and converts it into obstruction information.
    logical, dimension(:), intent(in) :: VALID
    integer, dimension(:), intent(in) :: MAFRANGE
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

    ! Local variables
    integer :: MAF                      ! Loop counter
    logical :: LASTONEVALID             ! Flag
    type (Obstruction_T) :: NEWOBSTRUCTION ! In progrss

    ! Executable code
    lastOneValid = .true.
    do maf = 1, size(valid)
      if ( valid(maf) .neqv. lastOneValid ) then
        ! A transition either from good to bad or bad to good
        if ( .not. valid(maf) ) then
          ! From good to bad
          newObstruction%range = .true.
          newObstruction%mafs(1) = maf + mafRange(1) - 1
        else
          newObstruction%mafs(2) = maf + mafRange(1) - 2
          call AddObstructionToDatabase ( obstructions, newObstruction )
        end if
      end if
      lastOneValid = valid(maf)
    end do
    
    ! Make sure any range at the end gets added
    if ( .not. lastOneValid ) then
      newObstruction%mafs(2) = mafRange(2)
      call AddObstructionToDatabase ( obstructions, newObstruction )
    end if
  end subroutine ConvertFlagsToObstructions

  ! ----------------------------------- DealWithObstructions ------------
  subroutine DealWithObstructions ( chunks, obstructions )
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    ! This routine modifies the chunks according to the information
    ! given in the obstructions.

    ! Local variables
    type (MLSChunk_T) :: NEWCHUNK       ! A chunk to create
    integer :: CHUNK                    ! Index of chunk under some consideration
    integer :: DEADCHUNK                ! Index of chunk to kill
    integer :: FIRSTMAF                 ! Index of first MAF in range
    integer :: LASTMAF                  ! Index of last MAF in range
    integer :: MAF                      ! Index of MAF for wall
    integer :: OBSTRUCTION              ! Loop counter

    ! Executable code

    do obstruction = 1, size ( obstructions )
      ! Find chunk where this obstruction is/starts
      if ( obstructions(obstruction)%range ) then
        ! A range obstruction

        ! first identify the chunks whose non overlapped portion falls completely
        ! within the range and delete them
        firstMAF = obstructions(obstruction)%mafs(1)
        lastMAF = obstructions(obstruction)%mafs(2)
        insideRange: do
          deadChunk = FindFirst ( &
            & chunks%firstMAFIndex + chunks%noMAFsLowerOverlap >= firstMAF .and. &
            & chunks%lastMAFIndex - chunks%noMAFsUpperOverlap <= lastMAF )
          if ( deadChunk == 0 ) exit insideRange
          call DeleteChunk ( chunks, deadChunk )
        end do insideRange

        ! Now think about chunks who's non overlapped part completely
        ! encompas the range, and split them.
        chunk = FindFirst ( &
          & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap <= firstMAF ) .and. &
          & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap >= lastMAF ) )
        if ( chunk /= 0 ) then
          ! Create two new chunks with the wall between, first the lower portion
          newChunk = chunks ( chunk )
          newChunk%lastMAFIndex = firstMAF - 1
          newChunk%noMAFsUpperOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now the upper portion
          newChunk = chunks ( chunk )
          newChunk%firstMAFIndex = lastMAF + 1
          newChunk%noMAFsLowerOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Delete the old chunk
          call DeleteChunk ( chunks, chunk )
        endif

        ! So the cases we have left do not involve deleting or creating new
        ! chunks, just modifying existing ones
        ! Look for chunks where the range starts in the chunk
        chunk = FindFirst ( chunks%lastMAFIndex >= firstMAF )
        if ( chunk /= 0 ) then
          chunks(chunk)%noMAFsUpperOverlap = max ( 0, &
            & chunks(chunk)%noMAFsUpperOverlap - &
            &    ( chunks(chunk)%lastMAFIndex - (firstMAF-1) ) )
          chunks(chunk)%lastMAFIndex = firstMAF - 1
        end if
        ! Look for chunks where the range ends in the chunk
        chunk = FindFirst ( chunks%firstMAFIndex <= lastMAF )
        if ( chunk /= 0 ) then
          chunks(chunk)%noMAFsLowerOverlap = max ( 0, &
            & chunks(chunk)%noMAFsLowerOverlap - &
            &    ( (lastMAF+1) - chunks(chunk)%firstMafIndex ) )
          chunks(chunk)%firstMAFIndex = lastMAF + 1
        end if

      else ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! A wall obstruction
        maf = obstructions(obstruction)%mafs(1)

        chunk = FindFirst ( &
          & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap <= maf ) .and. &
          & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap >= maf ) )

        if ( chunk /= 0 ) then
          ! Create two new chunks with the wall between, first the lower portion
          newChunk = chunks ( chunk )
          newChunk%lastMAFIndex = maf - 1
          newChunk%noMAFsUpperOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now the upper portion
          newChunk = chunks ( chunk )
          newChunk%firstMAFIndex = maf
          newChunk%noMAFsLowerOverlap = 0
          if ( newChunk%lastMAFIndex /= newChunk%firstMAFIndex ) &
            & call AddChunkToDatabase ( chunks, newChunk )
          ! Now delete the original chunk
          call DeleteChunk ( chunks, chunk )
        endif

        ! For chunks where the wall is in the overlap, just make the 
        ! overlap shorter 
        ! First the lower overlap
        chunk = FindFirst ( &
          & ( chunks%firstMAFIndex <= maf ) .and. &
          & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap - 1 >= maf ) )
        if ( chunk /= 0 ) then
          chunks(chunk)%noMAFsLowerOverlap = &
          & max ( chunks(chunk)%noMAFsLowerOverlap - &
          & ( maf - chunks(chunk)%firstMAFIndex ), 0 )
          chunks(chunk)%firstMAFIndex = maf
        end if
        ! Now the upper overlap
        chunk = FindFirst ( &
          & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap + 1 <= maf ) .and. &
          & ( chunks%lastMAFIndex >= maf ) )
        if ( chunk /= 0 ) then
          chunks(chunk)%noMAFsUpperOverlap = &
          & max ( chunks(chunk)%noMAFsUpperOverlap - &
          & ( chunks(chunk)%lastMAFIndex - (maf-1) ), 0 )
          chunks(chunk)%lastMAFIndex = maf - 1
        end if
      end if                            ! Wall obstructions
    end do                              ! Loop over obstructions

    ! Sort the chunks back into order
    call SortChunks ( chunks )

  end subroutine DealWithObstructions

  ! ------------------------------------------- DeleteObstruction -------
  subroutine DeleteObstruction ( obstructions, index )
    ! Dummy arguments
    type (Obstruction_T), pointer, dimension(:) :: OBSTRUCTIONS
    integer, intent(in) :: INDEX

    ! Local variables
    type (Obstruction_T), pointer, dimension(:) :: TEMP
    integer :: STATUS                   ! From allocate

    ! Executable code
    allocate ( temp ( size(obstructions) - 1 ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'temp' )

    if ( index > 1 ) temp(1:index-1) = obstructions(1:index-1)
    if ( index < size(obstructions) .and. size(obstructions) > 1 ) &
      & temp(index:) = obstructions(index+1:)

    deallocate ( obstructions, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'obstructions' )

    obstructions => temp

  end subroutine DeleteObstruction

  ! ------------------------------------------- DeleteChunk -------
  subroutine DeleteChunk ( chunks, index )
    ! Dummy arguments
    type (MLSChunk_T), pointer, dimension(:) :: CHUNKS
    integer, intent(in) :: INDEX

    ! Local variables
    type (MLSChunk_T), pointer, dimension(:) :: TEMP
    integer :: STATUS                   ! From allocate

    ! Executable code
    allocate ( temp ( size(chunks) - 1 ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'temp' )

    if ( index > 1 ) temp(1:index-1) = chunks(1:index-1)
    if ( index < size(chunks) .and. size(chunks) > 1 ) &
      & temp(index:) = chunks(index+1:)

    deallocate ( chunks, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'chunks' )

    chunks => temp

  end subroutine DeleteChunk

  ! --------------------------------------- Dump Obstructions ------------
  subroutine Dump_Obstructions ( obstructions ) 
    type (Obstruction_T), dimension(:), pointer :: obstructions
    ! Local variables
    integer :: i                        ! Loop counter
    ! Executable code
    if ( associated ( obstructions ) ) then
      if ( size(obstructions) == 0 ) then
        call output ( 'Obstructions is a zero size array.', advance='yes' )
      else
        call output ( 'Dumping ' )
        call output ( size(obstructions) )
        call output ( ' obstructions:', advance='yes' )
        do i = 1, size(obstructions)
          call output ( i )
          if ( obstructions(i)%range ) then
            call output ( ' : Range [ ' )
            call output ( obstructions(i)%mafs(1) )
            call output ( ' : ' )
            call output ( obstructions(i)%mafs(2) )
            call output ( ' ]', advance='yes' )
          else
            call output ( ' : Wall [ ' )
            call output ( obstructions(i)%mafs(1) )
            call output ( ' ]', advance='yes' )
          end if
        end do
      end if
    else
      call output ( 'Obstructions is not associated.', advance='yes')
    end if
  end subroutine Dump_Obstructions

  ! --------------------------------------- Prune Obstructions ----------
  subroutine PruneObstructions ( obstructions )
    ! This routine merges overlapping range obstructions and deletes
    ! wall obstructions inside ranges.  The job is made easier
    ! by sorting the obstructions into order
    type(Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

    ! Local variables
    integer :: I,J                      ! Loop counters
    type (Obstruction_T) :: newObs      ! New Obstruction
    logical :: FOUNDONE                 ! Found at least one
    integer :: STATUS                   ! Flag from allocate

    ! Executable code
    ! If no obstructions make sure allocate to size zero, not just unassociated pointer
    if ( .not. associated(obstructions) ) then
      allocate ( obstructions(0), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'obstructions(0)' )
      return
    end if

    ! Otherwise, do the tidying up
    outerLoop: do
      foundOne = .false.
      i = 0
      call SortObstructions(obstructions)
      middleLoop: do
        i = i + 1
        if ( i >= size(obstructions) ) exit middleLoop
        j = i
        innerLoop: do
          j = j + 1
          if ( j > size(obstructions) ) exit innerLoop
          if ( all ( obstructions((/i,j/))%range ) ) then
            ! --------------------------- ( Range, range )
            if ( obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) ) then
              ! Combine overlapping range obstructions
              newObs%range = .true.
              newObs%mafs(1) = obstructions(i)%mafs(1)
              newObs%mafs(2) = &
                & max ( obstructions(i)%mafs(2), obstructions(j)%mafs(2) )
              call DeleteObstruction ( obstructions, i )
              call DeleteObstruction ( obstructions, j )
              call AddObstructionToDatabase ( obstructions, newObs )
              call SortObstructions ( obstructions )
              foundOne = .true.
              exit middleLoop
            end if
          else if ( obstructions(i)%range .and. .not. obstructions(j)%range ) then
            ! --------------------------- ( Range, wall )
            if ( obstructions(j)%mafs(1) >= obstructions(i)%mafs(1) .and. &
              &  obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) ) then
              ! Delete wall obstruction inside range
              call DeleteObstruction ( obstructions, j )
              foundOne = .true.
              exit middleLoop
            end if
          else
            ! --------------------------- ( Wall, range ) or ( Wall, wall )
            ! Becuase the obstructions are in order, we know in the wall, range
            ! case that the wall must be at the start of the range, not inside it.
            if ( obstructions(i)%mafs(1) == obstructions(j)%mafs(1) ) then
              ! Delete wall obstruction at start of a range or at another wall
              call DeleteObstruction ( obstructions, i )
              foundOne = .true.
              exit MiddleLoop
            end if
          end if
          ! I'm pretty sure this covers all the possibilities.  It might seem
          ! not at first glance, but I think the fact that I always re-sort the
          ! obstructions into order means that the above code does catch everything.

        end do innerLoop
      end do middleLoop 
      if ( .not. foundOne ) exit outerLoop
    end do outerLoop

  end subroutine PruneObstructions
  
  ! ----------------------------------------- SortObstructions ----------
  subroutine SortObstructions ( obstructions )
    ! Sort the obstructions into order of increasing
    ! mafs(1) (start/wall MAF index)
    type (Obstruction_T), dimension(:), intent(inout) :: OBSTRUCTIONS

    ! Local variables
    type (Obstruction_T) :: TEMP
    integer :: I                        ! Loop counters
    integer, dimension(1) :: TOSWAP     ! Index

    ! Executable code
    do i = 1, size(obstructions) - 1
      toSwap = minloc ( obstructions(i:)%mafs(1) ) + (/ i-1 /)
      if ( toSwap(1) /= i ) then
        temp = obstructions(i)
        obstructions(i) = obstructions(toSwap(1))
        obstructions(toSwap(1)) = temp
      end if
    end do
  end subroutine SortObstructions

  ! ----------------------------------------- SortChunks ----------
  subroutine SortChunks ( chunks )
    ! Sort the chunks into order of increasing firstMAFIndex
    type (MLSChunk_T), dimension(:), intent(inout) :: CHUNKS

    ! Local variables
    type (MLSChunk_T) :: TEMP
    integer :: I                        ! Loop counters
    integer, dimension(1) :: TOSWAP     ! Index

    ! Executable code
    do i = 1, size(chunks) - 1
      toSwap = minloc ( chunks(i:)%firstMAFIndex ) + (/ i-1 /)
      if ( toSwap(1) /= i ) then
        temp = chunks(i)
        chunks(i) = chunks(toSwap(1))
        chunks(toSwap(1)) = temp
      end if
    end do
  end subroutine SortChunks

  ! ------------------------------------------- SayTime -----------------
  subroutine SayTime
    real :: T2
    call time_now ( t2 )
    if ( total_times ) then
      call output ( "Total time = " )
      call output ( dble(t2), advance = 'no' )
      call blanks ( 4, advance = 'no' )
    endif
    call output ( 'Timing for ScanDivide = ' )
    call output ( dble(t2-t1), advance='yes' )
    timing = .false.
  end subroutine SayTime

  ! ------------------------------------------ SurveyL1BData -----------
  subroutine SurveyL1BData ( processingRange, l1bInfo, config, mafRange,&
    & obstructions )
    ! This goes through the L1B data files and trys to spot possible
    ! obstructions.
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    integer, dimension(2), intent(out) :: MAFRANGE   ! Processing range in MAFs
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

    ! Local variables
    type (L1BData_T) :: TAITIME         ! Read from L1BOA file
    type (L1BData_T) :: TPGEODALT       ! Read from L1BOA file
    type (Obstruction_T) :: NEWOBSTRUCTION ! A single obstruction

    character(len=10) :: MODNAMESTR     ! Module name

    integer :: FLAG                     ! From L1B
    integer :: MAF                      ! Loop inductor
    integer :: MOD                      ! Loop inductor
    integer :: NOMAFS                   ! Number of MAFs in processing range
    integer :: NOMAFSREAD               ! From L1B

    logical :: LASTONEVALID             ! To run through valid
    logical :: THISONEVALID             ! To go into valid
    logical, dimension(:), pointer :: VALID ! Flag for each MAF

    real(r8) :: SCANMAX                 ! Range of scan each maf
    real(r8) :: SCANMIN                 ! Range of scan each mif

    ! Executable code
    
    ! Read time from the L1BOA file
    call ReadL1BData ( l1bInfo%l1boaId, 'MAFStartTimeTAI', taiTime, noMAFsRead, flag )

    ! Deduce the first and last MAFs to consider
    call Hunt ( taiTime%dpField(1,1,:), &
      & (/ processingRange%startTime, processingRange%endTime /), &
      & mafRange, allowTopValue=.true., allowBelowValue=.true. )

    ! Check the validity of the MAF range returned
    if ( mafRange(2) == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L1B data starts after requested processing range' )
    if ( mafRange(1) == taiTime%noMAFs ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'L1B data ends before requested processing range' )
    mafRange = mafRange - 1             ! Index from zero
    noMAFS = mafRange(2) - mafRange(1) + 1

    ! At this point we'd look through the L1B data for data gaps.
    ! I want to defer this for now until I've reached agreement with VSP, RFJ
    ! and RRL on how L1 would report these. NJL

    ! Now look through the L1B data first look for scan problems
    if ( config%criticalModules /= l_none ) then
      nullify ( valid )
      call Allocate_test ( valid, noMAFs, 'valid', ModuleName )
      if ( config%criticalModules == l_both ) then
        valid = .true.
      else
        valid = .false.
      endif
      do mod = 1, size(modules)
        call get_string ( modules(mod)%name, modNameStr, strip=.true. )
        if ( .not. modules(mod)%spacecraft .and. &
          & any ( config%criticalModules == &
          &     (/ modules(mod)%name, l_either, l_both /) ) ) then
          ! Read the tangent point altitude
          call ReadL1BData ( l1bInfo%l1boaID, trim(modNameStr)//'.tpGeodAlt', &
            & tpGeodAlt, noMAFsRead, flag )
          ! Consider the scan range in each MAF in turn
          do maf = 1, noMAFs
            scanMax = maxval ( tpGeodAlt%dpField(1,:,maf) )
            scanMin = minval ( tpGeodAlt%dpField(1,:,maf) )
            thisOneValid = ( scanMin >= config%scanLowerLimit(1) .and. &
              &              scanMin <= config%scanLowerLimit(2) ) .and. &
              &            ( scanMax >= config%scanUpperLimit(1) .and. &
              &              scanMax <= config%scanUpperLimit(2) )
            if ( config%criticalModules == l_both ) then
              valid(maf) = valid(maf) .and. thisOneValid
            else
              valid(maf) = valid(maf) .or. thisOneValid
            end if
          end do                        ! Maf loop
          call DeallocateL1BData ( taiTime )
        end if                          ! Consider this module
      end do                            ! Module Loop
      ! Convert this information into obstructions and tidy up.
      call ConvertFlagsToObstructions ( valid, mafRange, obstructions )
      call Deallocate_test ( valid, 'valid', ModuleName )
    end if                              ! Consider scan issues

    ! Here we'll eventually look at radiances and switch changes.  For the
    ! moment I'm deferring this one too. NJL.

    ! Sort the obstructions into order and prune them of repeats, overlaps etc.
    call PruneObstructions ( obstructions ) 

    ! Tidy up
    call DeallocateL1BData ( taiTime )

  end subroutine SurveyL1BData

end module ChunkDivide_m

! $Log$
! Revision 2.15  2001/11/19 23:53:50  livesey
! Works better, fixed bug with cases with all missing data, also handles
! ends of processing range properly.
!
! Revision 2.14  2001/11/19 23:33:45  livesey
! Interim version, some bugs to track down
!
! Revision 2.13  2001/11/16 20:38:42  livesey
! Nullified some pointers I missed
!
! Revision 2.12  2001/11/15 17:43:46  livesey
! Tidied up, more comments etc.
!
! Revision 2.11  2001/11/14 22:33:40  livesey
! This version seems pretty good.
!
! Revision 2.10  2001/11/14 01:49:12  livesey
! Improvements, getting closer
!
! Revision 2.9  2001/11/12 21:15:34  livesey
! More bug fixes
!
! Revision 2.8  2001/11/12 20:31:23  livesey
! Bug fix check in
!
! Revision 2.7  2001/11/10 01:09:41  livesey
! Bug fixes
!
! Revision 2.6  2001/11/10 01:01:03  livesey
! Just tidied up
!
! Revision 2.5  2001/11/10 00:56:24  livesey
! OK, this is getting close to being ready for testing.
!
! Revision 2.4  2001/11/10 00:03:44  livesey
! More work, still got issues!
!
! Revision 2.3  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.2  2001/11/09 06:34:39  livesey
! Minor bug fix, and added Log stuff
!
