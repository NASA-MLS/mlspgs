! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ChunkDivide_m

  ! This module replaces the old ScanDivide, and is a new approach to dividing
  ! the data into chunks.

  use Intrinsic, only: L_NONE, PHYQ_INVALID
  use MLSCommon, only: RP

  implicit none
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

! ChunkDivide              Divide MAFs in processing range among chunks  
! DestroyChunkDatabase     Deallocate memory taken by chunk database     
! ReduceChunkDatabase      Reduce chunk database to [first, last] chunks
! === (end of toc) ===                                                   
! === (start of api) ===
! ChunkDivide (int root, TAI93_Range_T processingRange,
!    L1BInfo_T l1bInfo, *mlSChunk_T Chunks(:) )      
! DestroyChunkDatabase (*mlSChunk_T Chunks(:) )      
! === (end of api) ===
  public :: DestroyChunkDatabase, ChunkDivide, ReduceChunkDatabase

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! This type is filled by the l2cf and describes the configuration of the
  ! chunk division process.
  type ChunkDivideConfig_T
    integer :: method = l_none          ! See below.
    real(rp) :: maxLength = 0           ! Maximum length of chunks
    integer :: maxLengthFamily = PHYQ_Invalid ! PHYQ_Angle etc.
    integer :: noChunks = 0             ! Number of chunks for [fixed]
    real(rp) :: overlap = 0.0           ! Desired length of overlaps
    real(rp) :: lowerOverlap = 0.0      ! Desired length of lower overlaps
    real(rp) :: upperOverlap = 0.0      ! Desired length of lower overlaps
    integer :: overlapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
    integer :: lowerOverlapFamily = PHYQ_Invalid
    integer :: upperOverlapFamily = PHYQ_Invalid
    integer :: noSlaves = 0             ! Number of slave nodes [even]
    integer :: homeModule = l_none      ! Which module to consider [orbital]
    real(rp) :: homeGeodAngle = 0.0     ! Aim for one chunk to start here [orbital]
    logical   :: scanLLSet = .false.    ! True if scan lower limit should be used
    logical   :: scanULSet = .false.    ! True if scan upper limit should be used
    real(rp), dimension(2) :: scanLowerLimit ! Range for bottom of scan
    real(rp), dimension(2) :: scanUpperLimit ! Range for top of scan
    integer   :: criticalModules = l_none ! Which modules must be scanning
    character(len=160), dimension(:), pointer &
      & :: criticalSignals => null() ! Which signals must be on
    real(rp)  :: maxGap = 0.0           ! Length of time/MAFs/orbits allowed for gap
    integer   :: maxGapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
    logical   :: skipL1BCheck = .false. ! Don't check for l1b data probs
  end type ChunkDivideConfig_T

  ! The chunk divide methods are:
  !
  ! Fixed - Ignore the L1B file, just give a fixed set of chunks as described
  !
  ! PE - One chunk centred on the given orbit geodetic angle.
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

  interface dump
    module procedure DUMP_OBSTRUCTIONS
    module procedure DUMP_CRITICALSIGNALS
  end interface
  
  logical, parameter :: CHECKFORSHAREDMAFS = .true.
  logical, parameter :: CHECKFORMAFSINRANGE = .true.
  logical, parameter :: CHECKFORNONNEGOVLPS = .true.

contains ! ===================================== Public Procedures =====

  !----------------------------------------  DestroyChunkDatabase  -----
  subroutine DestroyChunkDatabase ( chunks )
    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use Chunks_m, only: MLSCHUNK_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Warning

    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS
    integer :: STATUS                   ! From deallocate
    integer :: CHUNK                    ! Index

    if ( .not. associated ( chunks ) ) return
    do chunk = 1, size ( chunks )
      call Deallocate_test ( chunks(chunk)%hGridOffsets, &
        & 'chunks(?)%hGridOffsets', ModuleName )
      call Deallocate_test ( chunks(chunk)%hGridTotals, &
        & 'chunks(?)%hGridTotals', ModuleName )
    end do
    deallocate ( chunks, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_DeAllocate // "Chunks" )
  end subroutine DestroyChunkDatabase

  !----------------------------------------  ReduceChunkDatabase  -----
  subroutine ReduceChunkDatabase ( chunks, firstChunk, lastChunk )
    use Chunks_m, only: MLSCHUNK_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate

    type (MLSChunk_T), dimension(:), pointer :: chunks
    integer, intent(in) :: firstChunk, lastChunk
    ! Local variables
    type (MLSChunk_T), dimension(:), pointer :: TEMPDATABASE
    integer :: newSize, status
    ! Executable
    if ( .not. associated ( chunks ) ) return
    newSize = lastChunk - firstChunk + 1
    if ( newSize < 1 ) return
    if ( lastChunk > size(chunks) ) return
    allocate(tempDatabase(newSize), STAT=status)
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "tempDatabase")

    tempDatabase(1:newSize) = chunks(firstChunk:lastChunk)
    call DestroyChunkDatabase ( chunks )
    chunks => tempDatabase

  end subroutine ReduceChunkDatabase

  ! ------------------------------------------------  Chunk Divide -----
  subroutine ChunkDivide ( root, processingRange, l1bInfo, chunks )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Chunks_m, only: DUMP, MLSCHUNK_T
    use Dump_0, only: DUMP
    use EXPR_M, only: EXPR
    use Init_Tables_Module, only: F_OVERLAP, F_LOWEROVERLAP, F_UPPEROVERLAP, &
      & F_MAXLENGTH, F_METHOD, F_NOCHUNKS, &
      & F_HOMEMODULE, F_CRITICALMODULES, F_CRITICALSIGNALS, F_HOMEGEODANGLE, &
      & F_SCANLOWERLIMIT, F_SCANUPPERLIMIT, F_NOSLAVES, F_SKIPL1BCHECK, &
      & FIELD_FIRST, FIELD_LAST, L_EVEN, &
      & L_FIXED, F_MAXGAP, L_ORBITAL, L_PE, S_CHUNKDIVIDE, L_BOTH, L_EITHER
    use Intrinsic, only: L_NONE, FIELD_INDICES, LIT_INDICES, PHYQ_ANGLE, &
      & PHYQ_DIMENSIONLESS, PHYQ_INVALID, PHYQ_LENGTH, PHYQ_MAFS, PHYQ_TIME
    use L1BData, only: DEALLOCATEL1BDATA, L1BDATA_T, NAME_LEN, READL1BDATA, &
      & AssembleL1BQtyName
    use Lexer_core, only: PRINT_SOURCE
    use MLSCommon, only: R8, RP, L1BINFO_T, TAI93_Range_T
    use MLSFiles, only: HDFVERSION_4, HDFVERSION_5, WILDCARDHDFVERSION, &
      & mls_hdf_version
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSSets, only: FINDFIRST
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
      & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, MLSMSG_WARNING
    use MLSNumerics, only: Hunt
    use MLSSignals_m, only: MODULES, ISMODULESPACECRAFT
    use MoreTree, only: GET_BOOLEAN, GET_FIELD_ID, GET_SPEC_ID
    use Output_M, only: BLANKS, OUTPUT
    use String_table, only: GET_STRING, DISPLAY_STRING
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, TOGGLE, SWITCHES
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUBTREE, SUB_ROSA
    use Tree_types, only: N_EQUAL, N_NAMED

    integer, intent(in) :: ROOT    ! Root of the L2CF tree for ChunkDivide
    type( L1BInfo_T ), intent(in) :: L1BINFO
    type( TAI93_Range_T ), intent(in) :: PROCESSINGRANGE
    type( MLSChunk_T ), dimension(:), pointer  :: CHUNKS

    ! Local variables
    type (ChunkDivideConfig_T) :: CONFIG ! Configuration
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    integer :: STATUS                   ! From deallocate
    integer, dimension(2) :: MAFRANGE   ! Processing range in MAFs

    ! For announce_error:
    integer :: ERROR                    ! Error level
    integer :: CHUNK                    ! Loop counter

    integer, parameter :: BadUnits = 1
    integer, parameter :: NotSpecified = BadUnits + 1
    integer, parameter :: Unnecessary = NotSpecified + 1

    ! For timing
    logical :: Timing
    real :: T1

    ! Executable code
    nullify(config%criticalSignals)    ! Just for Sun's compiler

    if ( toggle(gen) ) call trace_begin ( "ChunkDivide", root )

    timing = section_times
    if ( timing ) call time_now ( t1 )

    ! First decode the l2cf instructions
    call ChunkDivideL2CF ( root, config )

    ! For methods other than fixed, we want to survey the L1 data and note the
    ! location of obstructions
    nullify ( obstructions )
    if ( config%method /= l_fixed ) then
      call SurveyL1BData ( processingRange, l1bInfo, config, mafRange,&
      & obstructions )
      if ( index(switches, 'chu') /= 0 ) then
        call output ( 'Requested time range ' )          
        call output ( processingRange%startTime )              
        call output ( ' : ' )    
        call output ( processingRange%endTime, advance='yes' )    
        call output ( 'Corresponding MAF range ' )          
        call output ( mafRange(1) )              
        call output ( ' : ' )
        call output ( mafRange(2), advance='yes' )
      endif
    endif

    ! Now go place the chunks.
    select case ( config%method )
    case ( l_fixed )
      call ChunkDivide_Fixed ( config, chunks )
    case ( l_PE )
      call ChunkDivide_PE ( config, mafRange, l1bInfo, &
        & obstructions, chunks )
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
      if ( index(switches, 'chu') /= 0 ) call Dump_Obstructions ( obstructions )
      deallocate ( obstructions, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'obstructions' )
    end if
    if ( associated(config%criticalSignals) ) then
      if ( index(switches, 'chu') /= 0 ) call Dump_criticalSignals ( config )
      deallocate ( config%criticalSignals, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'config%criticalSignals' )
    end if
    ! Check that no 2 chunks share non-overlapped MAFs
    ! That mafrange(1) <= all_mafs <= mafrange(2)
    ! And that all overlaps are >= 0
    if ( config%method /= l_fixed ) &
      & call CheckChunkSanity ( chunks, mafRange )
    if ( .not. associated(chunks) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &  
        & 'ChunkDivide failed to associate the chunks pointer with a target' )     
    else if ( size(chunks) < 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &  
        & 'ChunkDivide failed to produce any chunks' )     
    endif

    ! Now go through and number the chunks
    do chunk = 1, size ( chunks )
      chunks(chunk)%chunkNumber = chunk
    end do

    if ( toggle(gen) ) call trace_end ( "ChunkDivide" )
    if ( index(switches, 'chu') /= 0 ) call dump ( chunks )

  contains

  ! ========================================== Internal Procedures =====

    !--------------------------------- Add obstruction to database -----
    subroutine AddObstructionToDatabase ( database, item )

      ! Dummy arguments
      type (Obstruction_T), dimension(:), pointer :: DATABASE
      type (Obstruction_T), intent(in) :: ITEM

      ! Local variables
      type (Obstruction_T), dimension(:), pointer :: TEMPDATABASE

      include "addItemToDatabase.f9h"

    end subroutine AddObstructionToDatabase

    !--------------------------------------- Add chunk to database -----
    subroutine AddChunkToDatabase ( database, item )

      ! Dummy arguments
      type (MLSChunk_T), dimension(:), pointer :: DATABASE
      type (MLSChunk_T), intent(in) :: ITEM

      ! Local variables
      type (MLSChunk_T), dimension(:), pointer :: TEMPDATABASE

      include "addItemToDatabase.f9h"

    end subroutine AddChunkToDatabase

    ! ---------------------------------------------- AnnounceError -----
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

    !-------------------------------------------- ChunkDivide_Even -----
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

    !------------------------------------------- ChunkDivide_Fixed -----
    subroutine ChunkDivide_Fixed ( config, chunks )
      type (ChunkDivideConfig_T), intent(in) :: CONFIG
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS

      ! Local variables
      integer :: I                        ! Loop inductor
      integer :: STATUS                   ! From allocate
      integer :: MAXLENGTH                ! nint(config%maxLength)
      integer :: LOWEROVERLAP             ! nint(config%lowerOverlap)
      integer :: UPPEROVERLAP             ! nint(config%upperOverlap)
      integer :: NONONOVERLAP             ! maxLength-2*overlap

      ! Exectuable code
      allocate ( chunks(config%noChunks), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'chunks' )

      maxLength = nint ( config%maxLength )
      lowerOverlap = nint ( config%lowerOverlap )
      upperOverlap = nint ( config%upperOverlap )
      noNonOverlap = maxLength - ( lowerOverlap + upperOverlap )
      do i = 1, config%noChunks
        chunks(i)%firstMAFIndex = max ( (i-1)*maxLength - lowerOverlap, 0 )
        chunks(i)%lastMAFIndex = i*maxLength + upperOverlap - 1
        chunks(i)%noMAFsUpperOverlap = upperOverlap
        if ( i > 1 ) then
          chunks(i)%noMAFsLowerOverlap = lowerOverlap
        else
          chunks(i)%noMAFsLowerOverlap = 0
        end if
      end do

    end subroutine ChunkDivide_Fixed

    !---------------------------------------------- ChunkDivide_PE -----
    subroutine ChunkDivide_PE ( config, mafRange, l1bInfo, &
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
      integer :: MAXLENGTH                ! Max length as integer (MAFs)

      real(r8) :: ANGLEINCREMENT          ! Increment in hunt for homeMAF
      real(r8) :: MAXANGLE                ! Of range in data
      real(r8) :: MAXTIME                 ! Time range in data
      real(r8) :: MINANGLE                ! Of range in data
      real(r8) :: MINTIME                 ! Time range in data
      real(r8) :: TESTANGLE               ! Angle to check for

      real(r8), dimension(:), pointer :: BOUNDARIES ! Used in placing chunks
      real(r8), dimension(:), pointer :: FIELD ! Used in placing chunks
      integer   ::                       l1b_hdf_version
      character(len=NAME_LEN) ::         MAF_start, tp_angle

      ! Executable code

      ! Read in the data we're going to need
      call get_string ( lit_indices(config%homeModule), modNameStr, &
        & strip=.true. )
      if ( LEVEL1_HDFVERSION /= WILDCARDHDFVERSION ) then
        l1b_hdf_version = LEVEL1_HDFVERSION
      else
        l1b_hdf_version = mls_hdf_version(trim(l1bInfo%L1BOAFileName))
        if ( l1b_hdf_version <= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
      endif
      MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
        .false. )
      tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
        & l1b_hdf_version, &
        & .false. )
      call ReadL1BData ( l1bInfo%l1bOAId, trim(tp_angle), &
        & tpGeodAngle, noMAFsRead, flag, hdfVersion=l1b_hdf_version )
      call ReadL1BData ( l1bInfo%l1bOAId, trim(MAF_start), &
        & taiTime, noMAFsRead, flag, hdfVersion=l1b_hdf_version )
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
      if ( tpGeodAngle%dpField(1,1,m1) < 0.0 ) orbit = orbit - 1
      testAngle = config%homeGeodAngle + orbit*360.0
      if ( config%maxLengthFamily == PHYQ_Angle ) then
        angleIncrement = config%maxLength
      else
        angleIncrement = 360.0
      end if

      maxLength = nint ( config%maxLength )

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
        ! Now if this is close enough, and has enough MAFs around it, accept it
        if ( ( abs ( tpGeodAngle%dpField(1,1,home) - &
          & testAngle ) < HomeAccuracy ) &
          & .and. ( home - maxLength/2 >= 1 ) &
          & .and. ( home + maxLength/2 <= noMAFs ) ) exit homeHuntLoop
        ! Otherwise, keep looking
        testAngle = testAngle + angleIncrement
      end do homeHuntLoop
      homeMAF = home + m1 - 1

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OK, now we have a home MAF, get a first cut for the chunks
      ! We work out the chunk ends for each chunk according to how the
      ! maxLength field is specified.
      ! config%maxLengthFamily == PHYQ_MAFs
      noChunks = 1

      ! Allocate the chunk
      allocate ( chunks(noChunks), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'chunks (maxLength/MAFs)' )

      ! Work out its position (in the file)
      chunks(1)%firstMAFIndex = homeMAF - maxLength/2 - 1
      chunks(1)%lastMAFIndex = homeMAF + maxLength/2 - 1

      chunks(1)%noMAFsUpperOverlap = 0
      chunks(1)%noMAFsLowerOverlap = 0

      call DeallocateL1BData ( tpGeodAngle )
      call DeallocateL1BData ( taiTime )

    end subroutine ChunkDivide_PE

    !----------------------------------------- ChunkDivide_Orbital -----
    subroutine ChunkDivide_Orbital ( config, mafRange, l1bInfo, &
      & obstructions, chunks )
      type (ChunkDivideConfig_T), intent(in) :: CONFIG
      integer, dimension(2), intent(in) :: MAFRANGE
      type (L1BInfo_T), intent(in) :: L1BINFO
      type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS

      ! Local parameters
      real(r8), parameter :: HOMEACCURACY = 3.0 ! Try to hit homeGeodAngle within this
      ! (see homeHuntLoop warning below)
      ! Local variables
      type (L1BData_T) :: TAITIME         ! From L1BOA
      type (L1BData_T) :: TPGEODANGLE     ! From L1BOA

      character(len=10) :: MODNAMESTR     ! Home module name as string

      integer :: CHUNK                    ! Loop counter
      integer :: FLAG                     ! From ReadL1B
      integer :: HOME                     ! Index of home MAF in array
      integer :: M1, M2                   ! MafRange + 1
      integer :: NOCHUNKSBELOWHOME        ! Used for placing chunks
      integer :: NOMAFSATORABOVEHOME      ! Fairly self descriptive
      integer :: NOMAFSBELOWHOME          ! Fairly self descriptive
      integer :: NOMAFS                   ! Number of MAFs to consider
      integer :: NOMAFSREAD               ! From ReadL1B
      integer :: ORBIT                    ! Used to locate home
      integer :: STATUS                   ! From allocate etc.
      integer :: NOCHUNKS                 ! Number of chunks
      integer :: MAXLENGTH                ! Max length as integer (MAFs)

      integer, dimension(:), pointer :: NEWFIRSTMAFS ! For thinking about overlaps
      integer, dimension(:), pointer :: NEWLASTMAFS ! For thinking about overlaps

      real(r8) :: ANGLEINCREMENT          ! Increment in hunt for home
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

      integer   ::                       l1b_hdf_version
      character(len=NAME_LEN) ::         MAF_start, tp_angle
      ! Executable code

      ! Read in the data we're going to need
      call get_string ( lit_indices(config%homeModule), modNameStr, strip=.true. )
      if ( LEVEL1_HDFVERSION /= WILDCARDHDFVERSION ) then
        l1b_hdf_version = LEVEL1_HDFVERSION
      else
        l1b_hdf_version = mls_hdf_version(trim(l1bInfo%L1BOAFileName))
        if ( l1b_hdf_version <= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
      endif
      MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
        .false. )
      tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', &
        & l1b_hdf_version, &
        & .false. )
      call ReadL1BData ( l1bInfo%l1bOAId, trim(tp_angle), &
        & tpGeodAngle, noMAFsRead, flag, hdfVersion=l1b_hdf_version )
      call ReadL1BData ( l1bInfo%l1bOAId, trim(MAF_start), &
        & taiTime, noMAFsRead, flag, hdfVersion=l1b_hdf_version )
      noMAFs = mafRange(2) - mafRange(1) + 1
      m1 = mafRange(1) + 1
      m2 = mafRange(2) + 1

      minAngle = minval ( tpGeodAngle%dpField(1,1,m1:m2) )
      maxAngle = maxval ( tpGeodAngle%dpField(1,1,m1:m2) )
      minTime = minval ( taiTime%dpField(1,1,m1:m2) )
      maxTime = maxval ( taiTime%dpField(1,1,m1:m2) )
      if ( index ( switches, 'chu' ) /= 0 ) then
        call output ( 'No MAFs in file: ' )
        call output ( noMAFsRead, advance='yes' )
        call output ( 'MAF time range: ' )
        call output ( minTime )
        call output ( ' : ' )
        call output ( maxTime, advance='yes' )
        call output ( 'Angle range: ' )
        call output ( minAngle )
        call output ( ' : ' )
        call output ( maxAngle, advance='yes' )
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! First try to locate the last MAF before the homeGeodAngle
      orbit = int ( tpGeodAngle%dpField(1,1,m1)/360.0 )
      if ( tpGeodAngle%dpField(1,1,m1) < 0.0 ) orbit = orbit - 1
      testAngle = config%homeGeodAngle + orbit*360.0
      if ( config%maxLengthFamily == PHYQ_Angle ) then
        angleIncrement = config%maxLength
      else
        angleIncrement = 360.0
      end if

      if ( index(switches, 'chu') /= 0 ) then    
        call output ( ' orbit  ', advance='no' )   
        call output ( orbit , advance='no' )      
        call output ( '    testAngle  ', advance='no' )   
        call output ( testAngle , advance='no' )      
        call output ( '    angleIncrement  ', advance='no' )   
        call output ( angleIncrement , advance='yes' )      
      end if
      ! In my opinion (paw) here's what the following loop should do:
      ! Find the 1st MAF within HOMEACCURACY of home_angle
      ! where home_angle has been corrected for the starting orbit number
      ! Afterwards, the preceding MAFs must be divided among one or more
      ! chunks, and the same done with subsequent MAFs
      !
      ! Instead what it actually does is
      ! Find the 1st MAF within HOMEACCURACY of (home_angle + n*angleIncrement)
      ! where home_angle has been corrected for the starting orbit number
      ! In effect the home_angle is set only within an unknown number
      ! of angleIncrements
      ! While there may be few cases in which they don't do about as well
      ! let this be a warning
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
      if ( index(switches, 'chu') /= 0 ) then    
        call output ( 'Test Angle  ' )   
        call output ( testAngle , advance='yes' )      
        call output ( 'Angle(home)  ' )   
        call output ( tpGeodAngle%dpField(1,1,home) , advance='yes' )      
        call output ( 'Home  ' )   
        call output ( home  )      
        call output ( 'Difference  ' )   
        call output (  abs ( tpGeodAngle%dpField(1,1,home) - &
          & testAngle ) , advance='yes' )      
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! OK, now we have a home MAF, get a first cut for the chunks
      ! We work out the chunk ends for each chunk according to how the
      ! maxLength field is specified.
      if ( config%maxLengthFamily == PHYQ_MAFs ) then
        maxLength = nint ( config%maxLength )
        noMAFsBelowHome = home - m1 + 1
        noChunksBelowHome = noMAFsBelowHome / maxLength
        if ( mod ( noMAFsBelowHome, maxLength ) /= 0 ) noChunksBelowHome = noChunksBelowHome + 1
        noMAFsAtOrAboveHome = m2 - home + 1
        if ( config%noChunks == 0 ) then
          ! If user did not request specific number of chunks choose them
          noChunks = noChunksBelowHome + noMAFsAtOrAboveHome / maxLength
          if ( mod ( noMAFsAtOrAboveHome, maxLength ) /= 0 ) &
            & noChunks = noChunks + 1
        else
          ! User requested specific number of chunks
          noChunks = config%noChunks
        end if

        ! Allocate the chunks
        allocate ( chunks(noChunks), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'chunks (maxLength/MAFs)' )

        ! Work out their positions
        do chunk = 1, noChunks
          chunks(chunk)%lastMAFIndex = home + &
            & ( chunk - noChunksBelowHome ) * maxLength - 1
          ! Subtract one to convert from index in array to index in file
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

        noMAFsBelowHome = -999
        noMAFsAtOrAboveHome = -999
        noChunksBelowHome = int ( ( homeV - minV ) / config%maxLength )
        if ( homeV > minV ) noChunksBelowHome = noChunksBelowHome + 1
        if ( config%noChunks == 0 ) then
          ! Choose the number of chunks ourselves
          noChunks = noChunksBelowHome + int ( ( maxV - homeV ) / config%maxLength )
          if ( homeV + config%maxLength * ( noChunks - noChunksBelowHome ) < maxV ) &
            & noChunks = noChunks + 1
        else
          noChunks = config%noChunks
        end if

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

        ! Do some dumping
        if ( index(switches, 'chu') /= 0 ) then
          call output ( ' minV: ' ) 
          call output ( minV  ) 
          call output ( ' maxV: ' ) 
          call output ( maxV , advance='yes' ) 
          call output ( ' homeV: ' ) 
          call output ( homeV , advance='yes' ) 
          call output ( ' noMAFSBelowHome: ' )
          call output ( noMAFsBelowHome, advance='yes' )
          call output ( ' noMAFSAtOrAboveHome: ' )
          call output ( noMAFsAtOrAboveHome, advance='yes' )
          call output ( ' NoChunks: ' ) 
          call output ( NoChunks , advance='yes' ) 
          call output ( ' noChunksBelowHome: ' ) 
          call output ( noChunksBelowHome , advance='yes' ) 
          call dump ( boundaries , 'boundaries' ) 
          call dump ( field, 'field' )
        end if

        call Hunt ( field, boundaries, chunks%lastMAFIndex, start=m1, &
          & allowTopValue=.true., nearest=.true. )
        call Deallocate_test ( boundaries, 'boundaries', ModuleName )
      end if

      ! Now deduce the chunk starts from the ends of their predecessors
      if ( noChunks > 1 ) &
        & chunks(2:noChunks)%firstMAFIndex = &
        & chunks(1:noChunks-1)%lastMAFIndex + 1
      chunks(1)%firstMAFIndex = m1

      ! Now bound the chunks to be within the processing range
      chunks%firstMAFIndex = min ( max ( chunks%firstMAFIndex, m1 ), m2 )
      chunks%lastMAFIndex = min ( max ( chunks%lastMAFIndex, m1 ), m2 )

      ! Now offset these to the index in the file not the array
      ! chunks%firstMAFIndex = chunks%firstMAFIndex + mafRange(1) - 1
      ! chunks%lastMAFIndex = chunks%lastMAFIndex + mafRange(1) - 1 
      chunks%firstMAFIndex = chunks%firstMAFIndex - 1
      chunks%lastMAFIndex = chunks%lastMAFIndex - 1 

      ! If at this point the last two chunks end in the same place, this is
      ! a subtle defect in our chunking algorihtm, lets avoid it
      if ( chunks(noChunks-1)%lastMAFIndex == chunks(noChunks)%lastMAFIndex ) then
        call DeleteChunk ( chunks, noChunks )
        noChunks = noChunks - 1
      end if

      ! Do some dumping
      if ( index(switches, 'chu') /= 0 ) then                      
        call dump ( chunks%lastMAFIndex , 'chunks%lastMAFIndex' )
        call dump ( chunks%firstMAFIndex , 'chunks%firstMAFIndex' )  
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Think about overlaps
      nullify ( newFirstMAFs, newLastMAFs )
      call Allocate_test ( newFirstMAFs, noChunks, 'newFirstMAFs', ModuleName )
      call Allocate_test ( newLastMAFs, noChunks, 'newLastMAFs', ModuleName )
      ! We could split things out to deal with mixed unit, but lets make life easier 
      ! for ourselves.  ChunkDivideL2PC has already insisted that 
      ! lowerOverlapFamily == upperOverlapFamily.
      if ( config%lowerOverlapFamily == PHYQ_MAFs ) then
        newFirstMAFs = max(chunks%firstMAFIndex - nint(config%lowerOverlap), &
          & mafRange(1) )
        newLastMAFs = min(chunks%lastMAFIndex + nint(config%upperOverlap), &
          & mafRange(2) )
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'The bit of code that deals with non-MAF overlaps is probably broken' )
!         ! For angle and time, they are similar enough we'll just do some stuff
!         ! with pointers to allow us to use common code to sort them out
!         ! Note that before we search only over the range in mafRange, here we
!         ! search over the whole file as we can spill over the processing range
!         ! for overlaps.
!         select case ( config%maxLengthFamily )
!         case ( PHYQ_Angle )
!           field => tpGeodAngle%dpField(1,1,:)
!           minV = minAngle
!         case ( PHYQ_Time )
!           field => taiTime%dpField(1,1,:)
!           minV = minTime
!         case ( PHYQ_MAFs)
!         end select
!         call Hunt ( field, field(chunks%firstMAFIndex-1) + &
!           & config%lowerOverlap, newFirstMAFs, allowTopValue=.true. )
!         call Hunt ( field, field(chunks%firstMAFIndex-1) - &
!           & config%upperOverlap, newLastMAFs, allowTopValue=.true. )
!         if ( index(switches, 'chu') /= 0 ) then
!           call dump ( field(chunks%firstMAFIndex-1) + &
!             & config%lowerOverlap , 'fields+lowerOverlap' )
!           call dump ( field(newFirstMAFs) , 'hunted values' )
!           call dump ( field(chunks%firstMAFIndex-1) - &
!             & config%upperOverlap , 'fields-upperOverlap' )
!           call dump ( field(newLastMAFs) , 'hunted values' )
!         end if
!         ! Correct this to be real MAF indices (starting from zero)
!         newFirstMAFs = newFirstMAFs - 1
!         newLastMAFs = newLastMAFs - 1
      end if
      chunks%noMAFsLowerOverlap = chunks%firstMAFIndex - newFirstMAFs
      chunks%noMAFsUpperOverlap = newLastMAFs - chunks%lastMAFIndex
      chunks%firstMAFIndex = newFirstMAFs
      chunks%lastMAFIndex = newLastMAFs
      if ( index(switches, 'chu') /= 0 ) then
        call dump ( newFirstMAFs , 'newFirstMAFs' ) 
        call dump ( newLastMAFs , 'newLastMAFs' ) 
        call dump ( chunks%noMAFsLowerOverlap , 'chunks%noMAFsLowerOverlap' ) 
        call dump ( chunks%noMAFsUpperOverlap , 'chunks%noMAFsUpperOverlap' ) 
      endif
      call Deallocate_test ( newFirstMAFs, 'newFirstMAFs', ModuleName )
      call Deallocate_test ( newLastMAFs, 'newLastMAFs', ModuleName )

      ! Delete any zero length or all overlapped chunks
      call PruneChunks ( chunks )

      if ( index ( switches, 'chu' ) /= 0 ) then
        call output ( 'Before dealing with obstructions', advance='yes' )
        call Dump ( chunks )
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Now think about the obstructions
      call DealWithObstructions ( chunks, obstructions )

      ! Delete any zero length or all overlapped chunks
      call PruneChunks ( chunks )

!       ! Forcibly zero out number of lower (upper) overlaps on 1st (last) chunks
!       noChunks = size ( chunks )
!       chunks(1)%noMAFsLowerOverlap = 0
!       chunks(noChunks)%noMAFsUpperOverlap = 0

      if ( index(switches, 'chu') /= 0 ) then
        call output ( 'After dealing with obstructions', advance='yes' )
        call Dump ( chunks )
      end if

      ! Tidy up
      call DeallocateL1BData ( tpGeodAngle )
      call DeallocateL1BData ( taiTime )

    end subroutine ChunkDivide_Orbital

    !--------------------------------------------- ChunkDivideL2CF -----
    subroutine ChunkDivideL2CF ( sectionRoot, config )
      ! This subroutine identifies, separates, and checks values from the section
      ! of the MLSCF (ChunkDivide) passed to Scan/Divide.
      integer, intent(in) :: SECTIONROOT    ! Root of the ChunkDivide section of the
      ! MLSCF abstract syntax tree
      type (ChunkDivideConfig_T), intent(out) :: CONFIG ! Result of operation

      integer, target, dimension(3) :: NeededForFixed = &
        & (/ f_noChunks, f_maxLength, f_overlap /)
      integer, target, dimension(7) :: NotWantedForFixed = &
        & (/ f_noSlaves, f_homeModule, f_homeGeodAngle, f_scanLowerLimit, &
        &    f_scanUpperLimit, f_criticalModules, f_maxGap /)

      integer, target, dimension(3) :: NeededForPE = &
        & (/ f_maxLength, f_homeModule, f_homeGeodAngle /)
      integer, target, dimension(7) :: NotWantedForPE = &
        & (/ f_noChunks, f_overlap, f_noSlaves, f_scanLowerLimit, &
        &    f_scanUpperLimit, f_criticalModules, f_maxGap /)

      integer, target, dimension(6) :: NeededForOrbital = &
        & (/ f_maxLength, f_overlap, f_homeModule, f_homeGeodAngle, &
        &    f_criticalModules, f_maxGap /)
      integer, target, dimension(1) :: NotWantedForOrbital = &
        & (/ f_noSlaves /)

      integer, target, dimension(6) :: NeededForEven = &
        & (/ f_maxLength, f_overlap, f_maxLength, f_noSlaves, f_maxGap, &
        &    f_criticalModules /)
      integer, target, dimension(3) :: NotWantedForEven = &
        & (/ f_noChunks, f_homeModule, f_homeGeodAngle /)

      ! Local variables
      integer :: FIELDINDEX               ! Tree type
      integer :: fieldValue               ! Node in the tree
      integer :: GSON                     ! A son of son the ChunkDivide section node
      integer :: I                        ! Loop inductor
      integer :: J                        ! Another loop inductor
      integer :: numCriticalSignals       ! How many crit. Sigs
      integer :: ROOT                     ! Root of ChunkDivide command
      integer :: signalsNode              ! Node where crit. Sig. begin
      integer :: SON                      ! A son of the ChunkDivide section node
      integer :: sub_rosa_index
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
          fieldValue = decoration(subtree(2,son)) ! The field's value
        else
          fieldValue = son
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
        case ( f_lowerOverlap )
          config%lowerOverlap = value(1)
          config%lowerOverlapFamily = units(1)
        case ( f_upperOverlap )
          config%upperOverlap = value(1)
          config%upperOverlapFamily = units(1)
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
        case ( f_criticalSignals )
          ! sub_rosa_index = sub_rosa(gson)
          ! call get_string ( sub_rosa_index, config%criticalSignals, &
          !  & strip=.true. )
          signalsNode = son
          numCriticalSignals = nsons(signalsNode) - 1
          call allocate_test(config%criticalSignals, numCriticalSignals, &
            & 'critical signals', ModuleName)
          do j = 2, nsons(signalsNode)    ! Skip name
            call get_string ( sub_rosa ( subtree (j, signalsNode) ), &
              & config%criticalSignals(j-1), strip=.true. )
          end do
        case ( f_maxGap )
          config%maxGap = value(1)
          config%maxGapFamily = units(1)
        case ( f_skipL1BCheck )
          config%skipL1BCheck = get_boolean ( fieldValue )
          if ( config%skipL1BCheck ) &
            & call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'You have elected to skip checking the l1b data for problems' )
        end select
      end do

      ! Now check the sanity of what we've been given, this varies a bit
      ! depending on the method
      select case ( config%method )
      case ( l_fixed )
        needed => NeededForFixed
        notWanted => NotWantedForFixed
      case ( l_PE )
        needed => NeededForPE
        notWanted => NotWantedForPE
      case ( l_orbital )
        needed => NeededForOrbital
        notWanted => NotWantedForOrbital
      case ( l_even )
        needed => NeededForEven
        notWanted => NotWantedForEven
      end select

      ! Some special thought for overlap cases
      if ( any ( got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) ) then
        if ( .not. all ( got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) ) &
          & call AnnounceError ( root, notSpecified, &
          & FindFirst ( .not. got ( (/ f_lowerOverlap, f_upperOverlap /) ) ) )
        if ( got ( f_overlap ) ) &
          & call AnnounceError ( root, unnecessary, f_overlap )
        got ( f_overlap ) = .true.
        ! Insist that upper/lower overlaps be specified in the same units
        if ( config%lowerOverlapFamily /= config%upperOverlapFamily ) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      else
        ! If single overlap specified, copy it into upper/lower
        if ( got ( f_overlap ) ) then
          config%lowerOverlap = config%overlap
          config%lowerOverlapFamily = config%overlapFamily
          config%upperOverlap = config%overlap
          config%upperOverlapFamily = config%overlapFamily
        end if
      end if

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
      select case ( config%method )
      case ( l_PE )
        if ( config%maxLengthFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_maxLength )
      case ( l_orbital )
        if (all(config%maxLengthFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_maxLength )
        if (all(config%lowerOverlapFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_lowerOverlap )
        if (all(config%upperOverlapFamily/=(/PHYQ_MAFs, PHYQ_Angle, PHYQ_Time/))) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      case ( l_fixed, l_even )
        if ( config%maxLengthFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_maxLength )
        if ( config%lowerOverlapFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_lowerOverlap )
        if ( config%upperOverlapFamily /= PHYQ_MAFs ) &
          & call AnnounceError ( root, badUnits, f_upperOverlap )
      end select

      if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Problem with ChunkDivide command' )

      ! That's it, we're all valid now.

    end subroutine ChunkDivideL2CF

    ! --------------------------------- ConvertFlagsToObstructions -----
    subroutine ConvertFlagsToObstructions ( valid, obstructions, mafRange )
      ! This routine takes an array of logicals indicating good/bad data
      ! and converts it into obstruction information.
      logical, dimension(:), intent(in) :: VALID
      type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
      integer, dimension(:), intent(in), optional :: MAFRANGE

      ! Local variables
      integer :: MAF                    ! Loop counter
      integer :: OFFSET                 ! MAF index offset
      logical :: LASTONEVALID           ! Flag
      type (Obstruction_T) :: NEWOBSTRUCTION ! In progrss

      ! Executable code
      lastOneValid = .true.
      offset = 0
      if ( present(mafRange) ) offset = mafRange(1)

      do maf = 1, size(valid)
        if ( valid(maf) .neqv. lastOneValid ) then
          ! A transition either from good to bad or bad to good
          if ( .not. valid(maf) ) then
            ! From good to bad
            newObstruction%range = .true.
            newObstruction%mafs(1) = maf - 1 + offset
          else
            newObstruction%mafs(2) = maf - 2 + offset
            call AddObstructionToDatabase ( obstructions, newObstruction )
          end if
        end if
        lastOneValid = valid(maf)
      end do

      ! Make sure any range at the end gets added
      if ( .not. lastOneValid ) then
        newObstruction%mafs(2) = size(valid) - 1 + offset
        call AddObstructionToDatabase ( obstructions, newObstruction )
      end if
    end subroutine ConvertFlagsToObstructions

    ! --------------------------------------- DealWithObstructions -----
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
          ! (we'll deal with the issue of ranges spilling into overlap regions below)
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
          ! chunks, just modifying existing ones.  Note that this includes
          ! cases where a range starts end ends in a chunk but one or other boundary
          ! is in the overlap regions.

          ! Look for chunks where the range starts in the chunk
          moveEnd: do
            chunk = FindFirst ( &
              & chunks%firstMAFIndex <= firstMAF .and. &
              & chunks%lastMAFIndex >= firstMAF )
            if ( chunk == 0 ) exit moveEnd
            chunks(chunk)%noMAFsUpperOverlap = max ( 0, &
              & chunks(chunk)%noMAFsUpperOverlap - &
              &    ( chunks(chunk)%lastMAFIndex - (firstMAF-1) ) )
            chunks(chunk)%lastMAFIndex = firstMAF - 1
          end do moveEnd

          ! Look for chunks where the range ends in the chunk
          moveStart: do
            chunk = FindFirst ( &
              & chunks%firstMAFIndex <= lastMAF .and. &
              & chunks%lastMAFIndex >= lastMAF )
            if ( chunk == 0 ) exit moveStart
            chunks(chunk)%noMAFsLowerOverlap = max ( 0, &
              & chunks(chunk)%noMAFsLowerOverlap - &
              &    ( (lastMAF+1) - chunks(chunk)%firstMafIndex ) )
            chunks(chunk)%firstMAFIndex = lastMAF + 1
          end do moveStart

        else ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

          ! A wall obstruction
          maf = obstructions(obstruction)%mafs(1)

          ! First look for where the wall occurs in the non overlapped part
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
          wallInLower: do
            chunk = FindFirst ( &
              & ( chunks%firstMAFIndex < maf ) .and. &
              & ( chunks%firstMAFIndex + chunks%noMAFsLowerOverlap > maf ) )
            if ( chunk == 0 ) exit wallInLower
            chunks(chunk)%noMAFsLowerOverlap = &
            & max ( chunks(chunk)%noMAFsLowerOverlap - &
            & ( maf - chunks(chunk)%firstMAFIndex ), 0 )
            chunks(chunk)%firstMAFIndex = maf
          end do wallInLower
          ! Now the upper overlap
          wallInUpper: do
            chunk = FindFirst ( &
              & ( chunks%lastMAFIndex - chunks%noMAFsUpperOverlap < maf ) .and. &
              & ( chunks%lastMAFIndex > maf ) )
            if ( chunk == 0 ) exit wallInUpper
            chunks(chunk)%noMAFsUpperOverlap = &
              & max ( chunks(chunk)%noMAFsUpperOverlap - &
              & ( chunks(chunk)%lastMAFIndex - (maf-1) ), 0 )
            chunks(chunk)%lastMAFIndex = maf - 1
          end do wallInUpper
        end if  ! - - - - - - - - - - - - ! Wall obstructions

      end do                              ! Loop over obstructions

      ! Sort the chunks back into order
      call SortChunks ( chunks )

    end subroutine DealWithObstructions

    ! ------------------------------------------------ DeleteChunk -----
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

    ! ------------------------------------------ DeleteObstruction -----
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

    ! ----------------------------------------- NoteL1BRADChanges -----
    subroutine NoteL1BRADChanges ( obstructions, mafRange, l1bInfo )
      use MLSSignals_m, only:GetModuleFromRadiometer, GetModuleFromSignal, &
        & GetRadiometerFromSignal, GetSignal, GetSignalName, &
        & Signal_T, SIGNALS, MODULES
      use MLSStrings, only: NumStringElements, GetStringElement
      use Parse_signal_m, only: Parse_signal
      ! This routine notes any lack of good data for one of the
      ! signals, and, depending on sensitivity,
      ! augments the database of obstructions
      type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
      type (Obstruction_T) :: NEWOBSTRUCTION ! In progrss
      integer, dimension(2), intent(in) :: MAFRANGE   ! Processing range in MAFs
      type (L1BInfo_T), intent(in) :: L1BINFO

      ! Local variables
      ! This next sets the threshold of what it takes to be an obstruction
      ! At the most sensitive, hair-trigger level, any change in goodness
      ! of any signal would suffice
      ! Otherwise, a dead-zone in any signal must persist for at least a
      ! duration of config%maxGap before we declare it an obstruction
      ! and add it to the database
      logical, parameter :: ANYCHANGEISOBSTRUCTION = .false.
      integer :: critical_index
      integer :: critical_sub_index
      integer, pointer, dimension(:) :: goods_after_gap
      integer, pointer, dimension(:) :: goodness_changes
      logical, dimension(size(signals)) :: good_after_maxgap
      logical, dimension(size(signals)) :: good_signals_now
      logical, dimension(size(signals)) :: good_signals_last
      integer, dimension(size(signals)) :: howlong_nogood
      integer :: i
      integer :: maf              ! mafrange(1) <= maf <= mafrange(2)
      integer :: maf_index   ! 1 <= maf_index <= mafrange(2)-mafrange(1)+1
      integer :: mafset
      integer :: mafset_end
      integer :: mafset_start
      integer, parameter :: MAXMAFSINSET = 250
      integer :: nmafsets
      integer :: num_goods_after_gap
      integer :: num_goodness_changes
      logical, dimension(:), pointer  :: or_valids_buffer
      integer :: SignalIndex
      character(len=40) :: signal_full
      character(len=40) :: signal_str
      integer, pointer :: Signal_Indices(:)         ! Indices in the signals
      logical, dimension(:,:), pointer  :: signals_buffer
      logical :: this_maf_valid
      logical, dimension(:), pointer  :: valids_buffer
      ! Won't check if there are no radiances
      if ( .not. associated(l1bInfo%l1bRadIDs) ) return
      ! Here we will loop over the signals database
      ! Searching for 
      ! (1) mafs where there is no good data for any of the signals
      ! (2) mafs where the good data is switched from one signal to another
      ! (See Construct QuantityTemplates below line 253)
      
      ! (Possibly) time-consuming step:
      ! Read through the l1b radiance files for all signals,
      ! noting which ones are good, which not
      
      nullify ( signals_buffer, goods_after_gap, goodness_changes )
      call allocate_test( &
        & signals_buffer, mafRange(2) - mafRange(1) + 1 , size(signals), &
        & 'signals_buffer', ModuleName)
      call allocate_test( goods_after_gap, mafRange(2) - mafRange(1) + 1,&
        & 'goods_after_gap', ModuleName)
      call allocate_test( goodness_changes, mafRange(2) - mafRange(1) + 1,&
        & 'goodness_changes', ModuleName)
      good_signals_now = .false.   ! Initializing
      do signalIndex=1, size(signals)
        if ( mafRange(2) - mafRange(1) + 1 <= MAXMAFSINSET ) then
          good_signals_now(signalIndex) = &
            & any_good_signaldata ( signalIndex, signals(signalIndex)%sideband, &
            &   l1bInfo, mafRange(1), mafRange(2), &
            &   signals_buffer(:,signalIndex), mafRange )
        else
          nmafsets = (mafRange(2) - mafRange(1))/MAXMAFSINSET + 1
          mafset_end = mafRange(1) - 1   ! A trick--mafset_start is mafRange(1)
          do mafset=1, nmafsets
            mafset_start = mafset_end + 1
            mafset_end = min ( mafRange(2), mafset_end + MAXMAFSINSET )
            good_signals_now(signalIndex) = &
              & any_good_signaldata ( signalIndex, signals(signalIndex)%sideband, &
              &   l1bInfo, mafset_start, mafset_end, &
              &   signals_buffer(:,signalIndex), mafRange )
          enddo
        endif
      enddo
      ! call dump ( signals_buffer, 'signals_buffer' )

      ! Task (1a): Find mafs where there is at least one signal which
      ! changes from either nogood to good or from good to nogood
      ! compared with the last maf
      ! Task (1b): Find mafs where there is at least one signal which
      ! changes from nogood to good after a dead zone 
      ! lasting at least maxGap mafs
      num_goodness_changes = 0
      num_goods_after_gap = 0
      howlong_nogood       = 0
      do maf = mafRange(1), mafRange(2)
        maf_index = maf - mafRange(1) + 1
        do signalIndex=1, size(signals)
          good_signals_now(signalIndex) = signals_buffer(maf_index, signalIndex)
          if ( .not. good_signals_now(signalIndex) ) &
            & howlong_nogood(signalIndex) = howlong_nogood(signalIndex) + 1
          good_after_maxgap(signalIndex) = &
            & good_signals_now(signalIndex) &
            & .and. &
            & ( howlong_nogood(signalIndex) > config%maxGap )
          if ( good_signals_now(signalIndex) ) howlong_nogood(signalIndex) = 0
        enddo
        if ( maf /= mafRange(1) ) then
          if ( Any(good_signals_now .neqv. good_signals_last) ) then
            num_goodness_changes = num_goodness_changes + 1
            goodness_changes(num_goodness_changes) = maf
          endif
        endif
        if ( Any(good_after_maxgap) ) then
          num_goods_after_gap = num_goods_after_gap + 1
          goods_after_gap(num_goods_after_gap) = maf
        endif
        good_signals_last = good_signals_now
      enddo
      
      ! Task (2): Find regions where there is no signal among at least one
      ! of the critical signals
      ! if (config%criticalSignals /= '' ) then
      if ( associated(config%criticalSignals) ) then
        nullify ( valids_buffer )
        call allocate_test( &
          & valids_buffer, mafRange(2) - mafRange(1) + 1 , &
          & 'valids_buffer', ModuleName)
        nullify ( or_valids_buffer )
        call allocate_test( &
          & or_valids_buffer, mafRange(2) - mafRange(1) + 1 , &
          & 'or_valids_buffer', ModuleName)
        valids_buffer = .true.
        if ( index(switches, 'chu') /= 0 ) then
          call output ( 'Checking for critical signals: ')
          do critical_index=1, size(config%criticalSignals)
            call output ( trim(config%criticalSignals(critical_index)), &
              & advance='yes')
          enddo
          call output ( &
            & 'Which signals match the incomplete signal string ' &
            &  // 'weve been given', advance='yes')
        endif
        ! What we're about to do is:
        ! Assume we're given an array of char strings
        ! array = [str_1, str_2, ..]
        ! and each string contains comma-separated signals
        ! str_1 = 'sig_1,sig_2,..'
        ! Then the rule is that for a maf to pass the critical signals test,
        ! for each array element at least one of the signals within that
        ! string must be turned on for that maf
        ! Thus the shorthand formula that we're and-ing the array elements
        ! and or-ing the list elements
        ! You might prefer the formula that a maf is invalid if for any of the
        ! array elements, none of its signals is turned on
        do critical_index=1, size(config%criticalSignals)
          or_valids_buffer = .false.
          do critical_sub_index=1, NumStringElements( &
            & trim(config%criticalSignals(critical_index)), .FALSE.)
            call GetStringElement( &
              & trim(config%criticalSignals(critical_index)), signal_str, &
              & critical_sub_index, .FALSE. )
            ! Which signals match the incomplete signal string we've been given?
            nullify(Signal_Indices)
            call Parse_signal(signal_str, signal_indices)
            if ( .not. associated(Signal_Indices) ) exit
            if ( index(switches, 'chu') /= 0 ) then
              call output ( 'Signal_Indices: ')
              call output ( Signal_Indices, advance='yes')
              do i=1, size(Signal_Indices)
                call GetSignalName(Signal_Indices(i), signal_full)
                call output ( 'Full Signal Name: ')
                call output ( trim(signal_full), advance='yes')
              enddo
            endif
            do maf = mafRange(1), mafRange(2)
              maf_index = maf - mafRange(1) + 1
              this_maf_valid = .false.   ! not valid unless at least one signal
              do i=1, size(Signal_Indices)
                signalIndex = Signal_Indices(i)
                this_maf_valid = &
                  & this_maf_valid .or. signals_buffer(maf_index, signalIndex)
              enddo
              or_valids_buffer(maf_index) = &
                & or_valids_buffer(maf_index) .or. this_maf_valid
            enddo
            call deallocate_test(Signal_Indices, 'Signal_Indices', ModuleName)
          enddo
          valids_buffer = &
            & valids_buffer .and. or_valids_buffer
        enddo
        call ConvertFlagsToObstructions ( valids_buffer, obstructions, mafRange )
        call deallocate_test(valids_buffer, 'valids_buffer', ModuleName)
        call deallocate_test(or_valids_buffer, 'or_valids_buffer', ModuleName)
      endif
      
      ! Depending on sensitivity, add these to Obstructions database
      if ( ANYCHANGEISOBSTRUCTION .and. num_goodness_changes > 0 ) then
        do mafset = 1, num_goodness_changes
          newObstruction%range = .false.
          newObstruction%mafs(1) = goodness_changes(mafset)
          newObstruction%mafs(2) = 0    ! For overzelous Lahey uninitialized checking
          call AddObstructionToDatabase ( obstructions, newObstruction )
        enddo
      elseif ( num_goods_after_gap > 0 ) then
        do mafset = 1, num_goods_after_gap
          newObstruction%range = .false.
          newObstruction%mafs(1) = goods_after_gap(mafset)
          newObstruction%mafs(2) = 0    ! For overzelous Lahey uninitialized checking
          call AddObstructionToDatabase ( obstructions, newObstruction )
        enddo
      endif
      ! OK, we have the mafs where the goodness changes, now what?
      call deallocate_test( signals_buffer, 'signals_buffer', ModuleName )
      call deallocate_test( goods_after_gap, 'goods_after_gap', ModuleName)
      call deallocate_test( goodness_changes, 'goodness_changes', ModuleName)
    end subroutine NoteL1BRADChanges

    ! ----------------------------------------- PruneChunks -----------
    subroutine PruneChunks ( chunks )
      type (MLSChunk_T), dimension(:), pointer :: CHUNKS
      integer :: CHUNK
      ! Now delete chunks that either:
      !  1 - Are nothing but overlap
      !  2 - Have <=0 MAFs
      pruneChunksLoop: do
        chunk = FindFirst ( &
          & ( chunks%noMAFsLowerOverlap + chunks%noMAFsUpperOverlap ) >= &
          & ( chunks%lastMAFIndex - chunks%firstMAFIndex + 1 ) &
          & .or.&
          & ( chunks%firstMAFIndex > chunks%lastMAFIndex ) )
        if ( chunk == 0 ) exit pruneChunksLoop
        call DeleteChunk ( chunks, chunk )
      end do pruneChunksLoop
    end subroutine PruneChunks

    ! ----------------------------------------- PruneObstructions -----
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
              if ( obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
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
                &  obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
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
                exit middleLoop
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

    ! ---------------------------------------------------- SayTime -----
    subroutine SayTime
      real :: T2
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( 'Timing for ChunkDivide = ' )
      call output ( dble(t2-t1), advance='yes' )
      timing = .false.
    end subroutine SayTime

    ! ---------------------------------------------------- CheckChunkSanity -----
    subroutine CheckChunkSanity ( chunks, mafRange )
      type (MLSChunk_T), dimension(:), intent(inout) :: CHUNKS
      integer, dimension(2), intent(in) :: MAFRANGE   ! Processing range in MAFs

      ! Local variables
      integer :: i
      logical :: sharing, outrange, negovlps

      ! Executable statements
      sharing = .false.
      outrange = .false.
      negovlps = .false.
      if ( CHECKFORSHAREDMAFS ) then
        do i = 1, size(chunks) - 1
          sharing = (chunks(i)%lastMAFIndex - chunks(i)%noMAFsUpperOverlap) >= &
            & chunks(i+1)%firstMAFIndex + chunks(i+1)%noMAFsLowerOverlap
          if ( sharing ) exit
        enddo
      endif
      outrange = any(chunks%firstMAFIndex < mafrange(1)) .or. &
        &        any(chunks%lastMAFIndex > mafrange(2))
      negovlps = any(chunks%noMAFsLowerOverlap < 0) .or. &
        &        any(chunks%noMAFsUpperOverlap < 0)
      if ( sharing ) then
        call output ( "Shared non-overlaps MAFs detected ", advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Shared non-overlapped MAFs' )
      elseif ( outrange ) then
        call output ( " MAFs detected outside range: " )
        call output ( mafRange, advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'MAFs outside range' )
      elseif ( negovlps ) then
        call output ( "Negative overlaps detected ", advance='yes' )
        call output ( "(Either tweak ChunkDivide section or get someone " )
        call output ( " to fix the code) ", advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Negative overlaps' )
      endif
    end subroutine CheckChunkSanity

    ! ------------------------------------------------- SortChunks -----
    subroutine SortChunks ( chunks )
      ! Sort the chunks into order of increasing firstMAFIndex
      type (MLSChunk_T), dimension(:), intent(inout) :: CHUNKS

      ! Local variables
      type (MLSChunk_T) :: TEMP
      integer :: I                        ! Loop counters
      integer, dimension(1) :: TOSWAP     ! Index

      ! Executable code
      do i = 1, size(chunks) - 1
        toSwap = minloc ( chunks(i:)%firstMAFIndex + &
          & chunks(i:)%noMAFsLowerOverlap ) + (/ i-1 /)
        if ( toSwap(1) /= i ) then
          temp = chunks(i)
          chunks(i) = chunks(toSwap(1))
          chunks(toSwap(1)) = temp
        end if
      end do
    end subroutine SortChunks

    ! ------------------------------------------- SortObstructions -----
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

    ! ---------------------------------------------- SurveyL1BData -----
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

      integer   ::                       l1b_hdf_version
      character(len=NAME_LEN) ::         MAF_start, tp_alt
      ! Executable code

      if ( LEVEL1_HDFVERSION /= WILDCARDHDFVERSION ) then
        l1b_hdf_version = LEVEL1_HDFVERSION
      else
        l1b_hdf_version = mls_hdf_version(trim(l1bInfo%L1BOAFileName))
        if ( l1b_hdf_version <= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
      endif
      MAF_start = AssembleL1BQtyName ( 'MAFStartTimeTAI', l1b_hdf_version, &
        .false. )
      ! tp_angle = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAngle', l1b_hdf_version, &
      !  .false. )
      ! Read time from the L1BOA file
      call ReadL1BData ( l1bInfo%l1boaId, trim(MAF_start), taiTime, &
        & noMAFsRead, flag, hdfVersion=l1b_hdf_version )

      ! Deduce the first and last MAFs to consider
      call Hunt ( taiTime%dpField(1,1,:), &
        & (/ processingRange%startTime, processingRange%endTime /), &
        & mafRange, allowTopValue=.true., allowBelowValue=.true. )

      ! Check the validity of the MAF range returned
      if ( mafRange(2) == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'L1B data starts after requested processing range' )
      if ( mafRange(1) == taiTime%noMAFs ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'L1B data ends before requested processing range' )
      mafRange = min ( noMAFsRead-1, max ( 0, mafRange - 1 ) )          ! Index from zero
      noMAFS = mafRange(2) - mafRange(1) + 1
      if ( mafRange(2) < mafRange(1) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'mafRange(2) < mafRange(1)' )
      elseif ( mafRange(2) < 1 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'mafRange(2) < 1' )
      endif

      ! Now look through the L1B data, first look for scan problems
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
          tp_alt = AssembleL1BQtyName ( trim(modNameStr)//'.tpGeodAlt', &
            & l1b_hdf_version, .false. )
          if ( .not. modules(mod)%spacecraft .and. &
            & any ( config%criticalModules == &
            &     (/ modules(mod)%name, l_either, l_both /) ) ) then
            ! Read the tangent point altitude
            call ReadL1BData ( l1bInfo%l1boaID, trim(tp_alt), &
              & tpGeodAlt, noMAFsRead, flag, hdfVersion=l1b_hdf_version )
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
        call ConvertFlagsToObstructions ( valid, obstructions )
        call Deallocate_test ( valid, 'valid', ModuleName )
      end if                              ! Consider scan issues

      ! Here we look at radiances and switch changes.
      if ( .not. config%skipL1BCheck) &
        call NoteL1BRADChanges ( obstructions, mafRange, l1bInfo ) 

      ! Sort the obstructions into order; prune them of repeats, overlaps etc.
      call PruneObstructions ( obstructions ) 

      ! Tidy up
      call DeallocateL1BData ( taiTime )

    end subroutine SurveyL1BData
  end subroutine ChunkDivide

  ! -------------------------------------------- Dump Dump_criticalSignals -----
  subroutine Dump_criticalSignals ( config )

    use Output_M, only: OUTPUT

    type (ChunkDivideConfig_T) :: CONFIG ! Configuration
    ! Local variables
    integer :: i                        ! Loop counter
    ! Executable code
    if ( associated ( config%criticalSignals ) ) then
      if ( size(config%criticalSignals) == 0 ) then
        call output ( 'criticalSignals is a zero size array.', advance='yes' )
      else
        call output ( 'Dumping ' )
        call output ( size(config%criticalSignals) )
        call output ( ' criticalSignals:', advance='yes' )
        do i = 1, size(config%criticalSignals)
          call output ( i )
          call output ( ': ' )
          call output ( trim(config%criticalSignals(i)), advance='yes' )
        end do
      end if
    else
      call output ( 'critical Signals is not associated.', advance='yes')
    end if
  end subroutine Dump_criticalSignals

  ! -------------------------------------------- Dump Obstructions -----
  subroutine Dump_Obstructions ( obstructions ) 

    use Output_M, only: OUTPUT

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

! -----------------------------------------------  ANY_GOOD_SIGNALDATA  -----
  function ANY_GOOD_SIGNALDATA ( signal, sideband, l1bInfo, maf, maf2, &
    & good_buffer, mafRange )  &
    & result (answer)
  ! Scalar use:
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  !
  ! Array use
  ! Return one logical value in good_buffer
  ! for each of the_mafs between maf and maf2
  ! Arguments

    use Allocate_Deallocate, only: Deallocate_Test
    use Chunks_m, only: MLSChunk_T
    use L1BData, only: L1BData_T, READL1BDATA, &
      & FindL1BData, AssembleL1BQtyName, PRECISIONSUFFIX, DEALLOCATEL1BDATA
    use MLSCommon, only: L1BInfo_T, RK => R8
    use MLSFiles, only: MLS_HDF_Version
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: GetSignalName

    integer, intent(in)                            :: signal
    integer, intent(in)                            :: sideband
    logical                                        :: answer
    integer, intent(in)                            :: maf
    type (l1bInfo_T), intent(in)                   :: L1bInfo
    integer, optional, intent(in)                  :: maf2
    logical, dimension(:), optional, intent(inout) :: good_buffer
    integer, dimension(2), optional, intent(in)    :: MAFRANGE   ! Processing 
    ! Private                                                   range in MAFs
    integer :: FileID, flag, noMAFs
    character(len=127)  :: namestring
    type (l1bData_T) :: MY_L1BDATA
    integer :: hdfVersion
    integer :: mymaf2
    integer :: the_maf
    integer :: maf_index   ! 1 <= maf_index <= mafrange(2)-mafrange(1)+1

    ! Executable
    answer = .false.
    mymaf2 = maf
    if ( present(maf2) ) mymaf2 = maf2

    ! Set defaults for the good_buffer if present
    if ( present(good_buffer) ) then
      if ( .not. present(mafrange) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'mafRange must be supplied to any_good_signaldata' )
      endif
      good_buffer ( maf - mafRange(1) + 1 : myMAF2 - mafRange(1) + 1 ) = .false.
    end if

    ! Quit if no l1b rad files.
    if ( .not. associated(l1bInfo%l1bRadIDs) ) return

    ! OK, try to find this item in an l1brad file
    hdfVersion = mls_hdf_version ( trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION )
    if ( hdfversion <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
    call GetSignalName ( signal, nameString, &                   
      & sideband=sideband, noChannels=.TRUE. )                     
    nameString = AssembleL1BQtyName ( nameString, hdfVersion, .false. )
    nameString = trim(nameString) // PRECISIONSUFFIX
    fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString, hdfVersion )
    ! If not found, exit appropriately
    if ( fileID <= 0 ) then
      answer = .false.
      return
    end if

    ! OK, we've found it.  Read the data in.
    call ReadL1BData ( fileID , nameString, my_l1bData, noMAFs, flag, &
      & firstMAF=maf, lastMAF=mymaf2, &
      & NeverFail= .true., hdfVersion=hdfVersion )
    ! Quit if the reading failed.
    if ( flag /= 0 ) then
      answer = .false.
      return
    end if

    ! Give detailed or curt response
    if ( present(good_buffer) ) then
      answer = .true.    ! This value should be ignored by the caller
      do the_maf = maf, mymaf2
        maf_index = the_maf - mafRange(1) + 1
        good_buffer(maf_index) = .not. &
          & all (my_l1bData%DpField(:,:, the_maf+1-maf) < 0.0_rk)
      end do
    else
      answer = .not. all (my_l1bData%DpField < 0._rk)
    end if
    call DeallocateL1bData ( my_l1bData )

  end function ANY_GOOD_SIGNALDATA

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ChunkDivide_m

! $Log$
! Revision 2.49  2004/07/31 19:58:28  livesey
! Various bug fixes and clean ups in the light of real data.
!
! Revision 2.47  2004/06/10 00:58:44  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.46  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.45  2003/08/27 20:08:04  livesey
! Removed print statements
!
! Revision 2.44  2003/08/26 19:49:22  livesey
! Removed another print statement.
!
! Revision 2.43  2003/08/26 19:43:42  livesey
! Commented out a print statement
!
! Revision 2.42  2003/08/26 18:05:14  livesey
! Various fixes to the obstruction stuff.  More to come probably.
!
! Revision 2.41  2003/08/25 23:43:08  pwagner
! Added ReduceChunkDatabase
!
! Revision 2.40  2003/08/21 22:51:18  livesey
! Removed yet more print statements!
!
! Revision 2.39  2003/08/21 22:49:22  livesey
! Removed more print statements etc.
!
! Revision 2.38  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.37  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.36  2003/06/09 22:55:23  pwagner
! Improved criticalSignals; some bug fixes
!
! Revision 2.35  2003/06/05 23:37:14  pwagner
! 1st version of criticalSignals in Chunk Divide Orbital
!
! Revision 2.34  2003/05/17 00:06:06  pwagner
! Wont check sidebands for missing radiances
!
! Revision 2.33  2003/05/09 16:43:05  pwagner
! Speedup of L1B radiance check
!
! Revision 2.32  2003/05/07 23:43:05  pwagner
! Optionally may skipL1BCheck
!
! Revision 2.31  2003/04/30 22:06:04  pwagner
! Shouldnt check for good signals if there arent any
!
! Revision 2.30  2003/04/28 23:07:00  pwagner
! Fleshed out notel1brad_changes; not yet tested where needed
!
! Revision 2.29  2003/01/06 20:13:09  livesey
! New split upper/lower overlaps
!
! Revision 2.28  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.27  2002/12/06 01:08:33  pwagner
! Less likely to bomb on single chunks
!
! Revision 2.26  2002/11/13 01:03:11  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.25  2002/11/06 00:21:16  pwagner
! Fixed non-zero starttime prob; some extra checks, printing
!
! Revision 2.24  2002/10/07 23:49:49  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.23  2002/10/07 18:00:10  livesey
! Added call to trace_end, whoops!
!
! Revision 2.22  2002/08/22 01:23:52  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.21  2002/08/04 15:59:45  mjf
! New method "PE" for a single chunk centred on (or near) a given phi.
!
! Revision 2.20  2002/05/24 20:57:34  livesey
! More revisions associated with being able to have a specific number of chunks
!
! Revision 2.19  2002/05/24 16:47:30  livesey
! Allowed user to request specific number of chunks for orbital
!
! Revision 2.18  2001/12/17 23:09:18  livesey
! Now deletes chunks that are nothing but overlap
!
! Revision 2.17  2001/12/16 00:56:43  livesey
! Added dump option
!
! Revision 2.16  2001/11/20 00:10:29  livesey
! Another bug fix
!
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
