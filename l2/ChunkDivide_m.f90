! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ChunkDivide_m

  use EXPR_M, only: EXPR
  use MLSCommon, only: R8, RP, L1BINFO_T, MLSCHUNK_T, TAI93_Range_T
  use Intrinsic, only: L_NONE, FIELD_INDICES
  use Tree, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUBTREE
  use Lexer_core, only: PRINT_SOURCE
  use Tree_types, only: N_EQUAL, N_NAMED
  use Init_Tables_Module, only: F_OVERLAP, F_MAXLENGTH, F_NOCHUNKS, &
    & F_METHOD, F_HOMEMODULE, F_CRITICALMODULES, F_HOMEGEODANGLE, F_SCANLOWERLIMIT, &
    & F_SCANUPPERLIMIT, F_NOSLAVES, FIELD_FIRST, FIELD_LAST, L_EVEN, &
    & L_FIXED, F_MAXGAP, L_ORBITAL, S_CHUNKDIVIDE
  use Units, only: PHYQ_INVALID, PHYQ_LENGTH, PHYQ_MAFS, PHYQ_TIME, PHYQ_ANGLE, &
    & PHYQ_LENGTH, PHYQ_DIMENSIONLESS
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, MLSMSG_WARNING
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
  use Output_M, only: BLANKS, OUTPUT
  use String_table, only: GET_STRING, DISPLAY_STRING
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

    ! Executable code
    
    if ( toggle(gen) ) call trace_begin ("ScanDivide", root )

    timing = section_times
    if ( timing ) call cpu_time ( t1 )

    ! First decode the l2cf instructions
    call ChunkDivideL2CF ( root, config )

    ! For methods other than fixed, we want to survey the L1 data and note the
    ! location of obstructions
    nullify ( obstructions )
    if ( config%method /= l_fixed ) &
      & call SurveyL1BData ( processingRange, l1bInfo, config, obstructions )

    ! Now go place the chunks.
    select case ( config%method )
    case ( l_fixed )
      call ChunkDivide_Fixed ( config, chunks )
    case ( l_orbital )
      call ChunkDivide_Orbital ( config, processingRange, l1bInfo, &
        & obstructions, chunks )
    case ( l_even )
      call ChunkDivide_Even ( config, processingRange, l1bInfo, &
        & obstructions, chunks )
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

  !----------------------------------------- ChunkDivide_Even --------
  subroutine ChunkDivide_Even ( config, processingRange, l1bInfo, &
    & obstructions, chunks )
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local variables
    integer :: I                        ! Loop inductor
    
    ! Exectuable code
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
      chunks(i)%firstMAFIndex = max ( (i-1)*noNonOverlap - overlap, 0 )
      chunks(i)%lastMAFIndex = i*noNonOverlap + overlap - 1
      chunks(i)%noMAFsUpperOverlap = overlap
      if ( i > 1 ) then
        chunks(i)%noMAFsLowerOverlap = overlap
      else
        chunks(i)%noMAFsLowerOverlap = 0
      end if
      chunks(i)%accumulatedMAFs = (i-1)*noNonOverlap
    end do

  end subroutine ChunkDivide_Fixed

  !----------------------------------------- ChunkDivide_Orbital --------
  subroutine ChunkDivide_Orbital ( config, processingRange, l1bInfo, &
    & obstructions, chunks )
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (Obstruction_T), dimension(:), intent(in) :: OBSTRUCTIONS
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local variables
    integer :: I                        ! Loop inductor
    
    ! Exectuable code
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
    logical :: GOT(field_first:field_last) = .false.
    integer :: I                        ! Loop inductor
    integer :: ROOT                     ! Root of ChunkDivide command
    integer :: SON                      ! A son of the ChunkDivide section node
    integer :: GSON                     ! A son of son the ChunkDivide section node
    integer :: UNITS(2)                 ! Units of expression
    real(rp) :: VALUE(2)                ! Value of expression
    integer, dimension(:), pointer :: NEEDED ! Which fields are needed
    integer, dimension(:), pointer :: NOTWANTED ! Which fields are not wanted

    ! Executable code

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
      ! Get value for this field if appropriate
      if ( nsons(son) > 1 ) call expr ( subtree(2,son), units, value )
      select case ( fieldIndex )
      case ( f_method )
        config%method = decoration ( son )
        if ( units(1) /= PHYQ_DimensionLess ) &
          & call AnnounceError ( root, BadUnits, son )
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
        config%homeModule = decoration ( son )
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
        config%criticalModules = decoration ( son )
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
        & call AnnounceError ( root, notSpecified, notWanted(i) )
    end do

    ! Now check the units for various cases
    if ( all ( config%maxGapFamily /= &
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
        call output (' has inappropriate units' )
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

  ! --------------------------------------- Prune Obstructions ----------
  subroutine PruneObstructions ( obstructions )
    ! Dummy arguments
    type(Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

    ! This routine merges overlapping range obstructions and deletes
    ! wall obstructions inside ranges.  Note that it assumes that the
    ! obstructions are sorted.

    ! Local variables
    integer :: I,J                      ! Loop counters
    type (Obstruction_T) :: newObs      ! New Obstruction
    logical :: FOUNDONE                 ! Found at least one

    ! Executable code
    
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
          ! --------------------------- ( Range, range )
          if ( all ( obstructions((/i,j/))%range ) ) then
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
            ! --------------------------- ( Range, wall )
          else if ( obstructions(i)%range .and. .not. obstructions(j)%range ) then
            if ( obstructions(j)%mafs(1) >= obstructions(i)%mafs(1) .and. &
              &  obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) ) then
              ! Delete wall obstruction inside range
              call DeleteObstruction ( obstructions, j )
              foundOne = .true.
              exit middleLoop
            end if
            ! --------------------------- ( Wall, range ) or ( Wall, wall )
          else
            if ( obstructions(i)%mafs(1) == obstructions(j)%mafs(1) ) then
              ! Delete wall obstruction at start of a range or at another wall
              call DeleteObstruction ( obstructions, i )
              foundOne = .true.
              exit Middle Loop
            end if
          end if

        end do innerLoop
      end do middleLoop 
      if ( .not. foundOne ) exit outerLoop
    end do outerLoop

  end subroutine PruneObstructions
  
  ! ----------------------------------------- SortObstructions ----------
  subroutine SortObstructions ( obstructions )
    ! Dummy arguments
    type (Obstruction_T), dimension(:), intent(inout) :: OBSTRUCTIONS

    ! Local variables
    type (Obstruction_T) :: TEMP
    integer :: I                        ! Loop counters
    integer, dimension(1) :: TOSWAP     ! Index

    ! Executable code
    do i = 1, size(obstructions) - 1
      toSwap = maxloc ( obstructions(i:)%mafs(1) ) + (/ i-1 /)
      if ( toSwap(1) /= i ) then
        temp = obstructions(i)
        obstructions(i) = obstructions(toSwap(1))
        obstructions(toSwap(1)) = temp
      end if
    end do
  end subroutine SortObstructions

  ! ------------------------------------------- SayTime -----------------
  subroutine SayTime
    real :: T2
    call cpu_time ( t2 )
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
  subroutine SurveyL1BData ( processingRange, l1bInfo, config, obstructions )
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (ChunkDivideConfig_T), intent(in) :: CONFIG
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS

  end subroutine SurveyL1BData

end module ChunkDivide_m
