module CFM_MLSSetup_m

   use CFM_Tree_Walker, only : Walk_Tree
   use ForwardModelConfig, only: ForwardModelConfig_T, dump
   use SpectroscopyCatalog_m, only: dump
   use QuantityTemplates, only: QuantityTemplate_T
   use VGridsDatabase, only: VGrid_T
   use LEXER_CORE, only: INIT_LEXER
   use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL
   use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE
   use INIT_TABLES_MODULE, only: INIT_TABLES
   use PARSER, only: CONFIGURATION
   use TREE_CHECKER, only: CHECK_TREE
   use H5LIB, ONLY: h5open_f, h5close_f
   ! We don't message here, but don't know why I have to use
   ! MLSMessageConfig, but this is no time to find out
   use MLSMessageModule, only: MLSMessageConfig!, MLSMSG_Allocate, &
!                               MLSMessage, MLSMSG_Error, MLSMSG_Deallocate, &
!                               MLSMSG_Warning
   use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
    & DESTROY_STRING_TABLE
   use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
   use Chunks_m, only: MLSChunk_T
   use L1BData, only: findMaxMaf, NAME_LEN, READL1BDATA, L1BDATA_T, &
                      DEALLOCATEL1BDATA, AssembleL1BQtyName
   use MLSFiles, only: InitializeMLSFile, mls_openFile, &
                       mls_closeFile, GetMLSFileByType
   use MLSCommon, only: MLSFile_T, TAI93_Range_T, r8
   use SDPToolkit, only: PGS_S_SUCCESS
   use Intrinsic, only: l_hdf, l_ascii, lit_indices
   use Hdf, only: DFACC_RDONLY
   use MLSNumerics, only: Hunt
   use ChunkDivide_m, only: ChunkDivideConfig_T, &
                            GetChunkFromTimeRange => CFM_ChunkDivide
   use Init_Tables_Module, only: l_none, phyq_mafs, phyq_time, phyq_angle, &
                                 l_orbital, l_fixed, l_both, l_either, &
                                 l_ghz
   use Allocate_Deallocate, only: allocate_test, deallocate_test
   use MLSSignals_m, only: MODULES
   use String_Table, only: get_string
   use MLSFillValues, only: ISFILLVALUE
   use MLSSets, only: FINDFIRST
   use SDPToolkit, only: pgs_td_utctotai, PGS_S_SUCCESS, PGSTD_E_NO_LEAP_SECS

   use input, only: l1boa

   implicit none

   private
   public :: CFM_MLSSetup, CFM_MLSCleanup
   public :: MLSChunk_T

   type MAFRange_T
      integer, dimension(2) :: L1BCover ! Range found in L1B files
      integer, dimension(2) :: L2Cover  ! Range covered by L2 ProcessingRange
      integer, dimension(2) :: Expanded ! L2Cover plus prior/post overlaps
   end type

   character(len=20), parameter :: moduleName = "CFM_MLSSetup"
   integer, parameter :: CCSDSLen = 27
   type(MLSFile_T), dimension(:), save, pointer :: filedatabase => NULL()
   type(ChunkDivideConfig_T), save :: chunkDivideConfig  ! Using default options

   contains

   subroutine CFM_MLSSetup (startTime, endTime, retVal, filedatabase, &
      fakeChunk, ForwardModelConfigDatabase)

      character(len=CCSDSlen), intent(in) :: startTime, endTime
      integer, intent(out) :: retVal
      type (MLSFile_T), dimension(:), pointer :: filedatabase
      type (MLSChunk_T), intent(out) :: fakeChunk
      type (ForwardModelConfig_T), pointer, optional :: ForwardModelConfigDatabase(:)

      integer :: Root
      integer :: First_Section
      integer :: error
      type (MLSFile_T), target :: l1bfile
      type (TAI93_Range_T) :: processingRange
      type (ForwardModelConfig_T), pointer :: dummy1(:) => NULL()

      !Executables
      retVal = 0
      if (present(ForwardModelConfigDatabase)) nullify (ForwardModelConfigDatabase)

      MLSMessageConfig%useToolkit = .false.
      MLSMessageConfig%logFileUnit = -1

      call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
      call allocate_decl ( ndecls=8000 )
      call allocate_tree ( n_tree=2000000 )
      call init_tables
      call configuration(Root)
      if (Root <= 0) then
         print *, 'A syntax error occurred -- there is no abstract syntax tree'
         retVal = 2
         return
      end if

      call check_tree ( root, error, first_section )
      if (error /= 0) then
         print *, "Error in tree"
         retVal = 3
         return
      end if

      call h5open_f(error)
      if (error /= 0) then
         print *, "Error in initialize hdf5 library"
         retVal = 4
         return
      end if

      if (present(ForwardModelConfigDatabase)) then
         call Walk_Tree ( Root, First_Section, ForwardModelConfigDatabase )
      else
         call Walk_Tree ( Root, First_Section, dummy1 )
         deallocate (dummy1)
      end if

      ! Create filedatabase
      nullify(filedatabase)
      allocate(filedatabase(1), stat=error)
      if (error /= 0) return
         error = InitializeMLSFile(l1bfile, content='l1boa', &
         name=trim(l1boa), &
         shortName='L1BOA', type=l_hdf, access=DFACC_RDONLY)
      if (error == 0) then
         call mls_openFile(L1BFile, error)
         if (error /= PGS_S_SUCCESS ) then
            error = -1
            return
         else
            filedatabase(1) = l1bfile
         end if

         !Populate fakeChunk
         error = pgs_td_utctotai (startTime, processingrange%starttime)
         if (error /= PGS_S_SUCCESS .and. error /= PGSTD_E_NO_LEAP_SECS) then
            print *, "Could not convert UTC Start time to TAI"
            retVal = 1
            return
         end if

         error = pgs_td_utctotai (endtime, processingrange%endtime)
         if (error /= PGS_S_SUCCESS .and. error /= PGSTD_E_NO_LEAP_SECS) then
            print *, "Could not convert UTC End time to TAI"
            retVal = 1
            return
         end if

         ! Initialize ChunkDivideConfig
         chunkDivideConfig%method = l_orbital
         chunkDivideConfig%maxLength = 0  ! chunk can be as big as possible
         chunkDivideConfig%maxLengthFamily = phyq_mafs
         chunkDivideConfig%skipL1BCheck = .true.
         chunkDivideConfig%homeModule = l_ghz
         chunkDivideConfig%criticalModules = l_ghz
         chunkDivideConfig%homeGeodAngle = 0.0_r8
         chunkDivideConfig%maxGap = 2.0_r8
         chunkDivideConfig%maxGapFamily = phyq_mafs
         chunkDivideConfig%maxOrbY = 50000.0_r8 ! Unit is meter
         chunkDivideConfig%scanLowerLimit = (/-20000.0_r8, 10000.0_r8 /) ! unit is meter
         chunkDivideConfig%scanUpperLimit = (/ 40000.0_r8, 200000.0_r8 /) ! unit is meter
         chunkDivideConfig%lowerOverlap = 0.0
         chunkDivideConfig%upperOverlap = 0.0
         chunkDivideConfig%lowerOverlapFamily = phyq_mafs
         chunkDivideConfig%upperOverlapFamily = phyq_mafs
         chunkDivideConfig%noChunks = 1
         fakeChunk = GetChunkFromTimeRange(processingRange, filedatabase, chunkDivideConfig)
!         fakeChunk%lastMafIndex = FindMaxMAF (filedatabase, fakeChunk%firstMafIndex )
!         fakeChunk%lastMafIndex = fakeChunk%lastMafIndex - fakeChunk%firstMafIndex
      end if

   end subroutine CFM_MLSSetup

   subroutine CFM_MLSCleanup

      integer :: error
      integer :: i

      call destroy_char_table
      call destroy_hash_table
      call destroy_string_table
      call destroy_symbol_table
      call deallocate_decl
      call deallocate_tree

      do i = 1, size(filedatabase)
         call mls_closefile(filedatabase(i))
      end do

      if (associated(filedatabase)) deallocate(filedatabase)

      call h5close_f (error)
      if (error /= 0) then
         print *, "Error in finishing up hdf5 library"
         return
      end if

   end subroutine CFM_MLSCleanup

!    type(MLSChunk_T) function GetChunkFromTimeRange &
!       (processingRange, filedatabase, retVal) result(chunk)
!
!       type (TAI93_Range_T), intent(in) :: processingRange
!       type (MLSFile_T), dimension(:), pointer :: filedatabase
!       integer, intent(out) :: retVal
!
!       type (MAFRange_T) :: MAFRange
!       type(Obstruction_T), dimension(:), pointer :: obstructions
!
!       nullify(chunkDivideConfig%criticalSignals)   ! Just for Sun's compiler
!       nullify(obstructions)
!
!       ! not sure how to deal with obstruction, ignore it
!       if (chunkDivideConfig%method /= l_fixed) then
!          call SurveyL1BData ( processingRange, filedatabase, &
!             mafRange, obstructions, retVal )
!          if (retVal /= 0) return
!       end if
!
!       if (chunkdivideconfig%method == l_orbital) then
!          call ChunkDivide_Orbital (mafRange, filedatabase, chunk, retVal)
!       else
!          print *, "chunk divide method is unsupported"
!       end if
!
!    end function
!
!    subroutine SurveyL1BData (processingRange, filedatabase, &
!       mafRange, obstructions, retVal)
!
!       type(TAI93_Range_T), intent(in) :: processingRange
!       type(MLSFile_T), dimension(:), pointer :: filedatabase
!       type(MAFRange_T), intent(out) :: mafRange
!       type(Obstruction_T), dimension(:), pointer :: obstructions
!       integer, intent(out) :: retVal
!
!       type(MLSFile_T), pointer :: l1bFile
!       character(len=NAME_LEN) :: maf_start, tp_alt, tp_orby, tp_angle
!       integer :: flag, noMafs, noMafs2, mod, maf
!       type(L1BData_T) :: taiTime, tpGeodAlt, tpOrbY, tpGeodAngle
!       logical :: THISONEVALID             ! To go into valid
!       logical, dimension(:), pointer :: VALID ! Flag for each MAF
!       logical, dimension(:), pointer :: WASSMOOTHED ! Flag for each MAF
!       logical, dimension(:), pointer :: ANGLEWASSMOOTHED ! Flag for each MAF
!       character(len=10) :: modNameStr   ! Module name
!       real(r8) :: scanMax, scanMin, OrbYMax
!
!       ! Executable
!       L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
!       if ( .not. associated(L1BFile) ) then
!          print *, "Can't make progress in SurveyL1BData without L1BOA files"
!          retVal = 1
!          return
!       end if
!
!       maf_start = AssembleL1BQtyName ('MAFStartTimeTAI', l1bfile%hdfVersion, &
!                                       .false.)
!       ! Read time from the L1BOA file
!       call ReadL1BData (l1bFile, trim(MAF_start), taiTime, noMafs, &
!                         flag, dontPad=.true.)
!
!       call SmoothOutDroppedMafs (taiTime%dpField)
!       ! We shall assume that all l1b files cover this range
!       mafRange%L1BCover(1) = taiTime%firstMaf
!       mafRange%L1BCover(2) = taiTime%noMafs - 1 ! mafs are 0-based indexed
!
!       ! Deduce the first and last MAFs covered the L2 Processing Range
!       call Hunt ( taiTime%dpField(1,1,:), &
!          & (/ processingRange%startTime, processingRange%endTime /), &
!          & mafRange%L2Cover, allowTopValue=.true., allowBelowValue=.true. )
!
!       ! Check the validity of the MAF range returned
!       if (mafRange%L2Cover(2) == 0) then
!          print *, 'L1B data starts after requested processing range'
!          retVal = 1
!          return
!       end if
!
!       if ( mafRange%L2Cover(1) == taiTime%noMAFs ) then
!          print *, 'L1B data ends before requested processing range'
!          retVal = 1
!          return
!       end if
!
!       if ( mafRange%L2Cover(2) < mafRange%L2Cover(1) ) then
!          print *, 'L2 mafRange(2) < L2 mafRange(1)'
!          retVal = 1
!          return
!       elseif ( mafRange%L2Cover(2) < 1 ) then
!          print *, 'L2 mafRange(2) < 1'
!          retVal = 1
!          return
!       endif
!
!       mafRange%Expanded(1) = mafRange%L2Cover(1)
!       mafRange%Expanded(2) = mafRange%L2Cover(2)
!       ! Expand the range covered by L2 with possible overlaps
!       if ( ChunkDivideConfig%allowPriorOverlaps ) then
!          mafRange%Expanded(1) = max(mafRange%L1BCover(1), &
!             mafRange%L2Cover(1) - nint(ChunkDivideConfig%lowerOverlap))
!       end if
!       if ( ChunkDivideConfig%allowPostOverlaps ) then
!          mafRange%Expanded(2) = min( mafRange%L1BCover(2), &
!             mafRange%L2Cover(2) + nint(ChunkDivideConfig%upperOverlap) )
!       endif
!       noMAFS = mafRange%Expanded(2) - mafRange%Expanded(1) + 1
!       ! Now look through the L1B data, first look for scan problems
!       if ( ChunkDivideConfig%criticalModules /= l_none ) then
!          nullify(valid, anglewasSmoothed, wasSmoothed)
!          call Allocate_test ( valid, noMAFs, 'valid', ModuleName )
!          call Allocate_test ( anglewasSmoothed, noMAFs, 'angleWasSmoothed', ModuleName )
!          call Allocate_test ( wasSmoothed, noMAFs, 'wasSmoothed', ModuleName )
!          if ( ChunkDivideConfig%criticalModules == l_both ) then
!             valid = .true.
!          else
!             valid = .false.
!          endif
!
!          do mod=1, size(modules)
!             call get_string (modules(mod)%name, modNameStr, strip=.true.)
!             tp_alt = AssembleL1BQtyName (trim(modNameStr) // '.tpGeodAlt', &
!                l1bfile%hdfVersion, .false.)
!             tp_orby = AssembleL1BQtyName (trim(modNameStr) // '.tpOrbY', &
!                l1bfile%hdfVersion, .false.)
!             tp_angle = AssembleL1BQtyName (trim(modNameStr) // '.tpGeodAngle', &
!                l1bfile%hdfVersion, .false.)
!             if (.not. modules(mod)%spacecraft .and. &
!                & (any(chunkdivideconfig%criticalmodules == (/l_either, l_both/)) &
!                & .or. lit_indices(chunkDivideConfig%criticalModules) == modules(mod)%name)) then
!                ! Read the tangent point altitude
!                call ReadL1BData (l1bfile, trim(tp_alt), tpGeodAlt, noMafs2, flag, &
!                   firstMaf=mafRange%Expanded(1), lastMaf=mafRange%Expanded(2), &
!                   dontPad=.true.)
!                ! Read the out of plane distance
!                call ReadL1BData (l1bfile, trim(tp_orby), tpOrbY, noMafs2, flag, &
!                   firstMaf=mafRange%Expanded(1), lastMaf=mafRange%Expanded(2), &
!                   dontPad=.true.)
!                call ReadL1BData (l1bfile, trim(tp_angle), tpGeodAngle, noMafs2, &
!                   flag, firstMaf=mafRange%Expanded(1), lastMaf=mafRange%Expanded(2), &
!                   dontPad=.true.)
!                call smoothOutDroppedMafs(tpGeodAngle%dpField, angleWasSmoothed, &
!                   monotonize=.true.)
!                call smoothOutDroppedMafs(tpGeodAlt%dpField, wasSmoothed)
!                call smoothOutDroppedMafs(tpOrbY%dpField)
!                wasSmoothed = (wasSmoothed .or. angleWasSmoothed)
!                ! Consider the scan range in each MAF in turn
!                do maf = 1, noMafs
!                   scanMax = maxval(tpGeodAlt%dpField(1,:,maf))
!                   scanMin = minval(tpGeodAlt%dpField(1,:,maf))
!                   orbyMax = maxval(tpOrbY%dpField(1,:,maf))
!                   thisOneValid = ( scanMin >= ChunkDivideConfig%scanLowerLimit(1) &
!                      .and. scanMin <= ChunkDivideConfig%scanLowerLimit(2) ) .and. &
!                      ( scanMax >= ChunkDivideConfig%scanUpperLimit(1) .and. &
!                      scanMax <= ChunkDivideConfig%scanUpperLimit(2) )
!                   if (chunkDivideConfig%maxOrbY > 0.0) then
!                      thisOneValid = thisOneValid .and. orbYMax < chunkDivideConfig%maxOrbY
!                   end if
!                   ! this one is not valid if it's valid only by virtue
!                   ! of having been smoothed
!                   thisOneValid = thisOneValid .and. .not. wasSmoothed(maf)
!                   if ( ChunkDivideConfig%criticalModules == l_both ) then
!                      valid(maf) = valid(maf) .and. thisOneValid
!                   else
!                      valid(maf) = valid(maf) .or. thisOneValid
!                   end if
!                end do   ! Maf loop
!
!                call DeallocateL1BData ( tpgeodalt )
!                call DeallocateL1BData ( tpgeodangle )
!                call DeallocateL1BData ( tporby )
!             end if   ! Consider this module
!          end do    ! Module Loop
!
!          ! Convert this information into obstructions and tidy up.
!          call ConvertFlagsToObstructions (valid, obstructions, retVal, &
!                                          obstructionType='scan')
!          if (retVal /= 0) return
!          call Deallocate_test ( valid, 'valid', ModuleName )
!          call Deallocate_test ( wasSmoothed, 'wasSmoothed', ModuleName )
!          call Deallocate_test ( angleWasSmoothed, 'angleWasSmoothed', ModuleName )
!       end if
!
!       ! Here we look at radiances and switch changes.
!       if (.not. chunkDivideConfig%skipL1BCheck) &
!          call NoteL1BRadChanges(obstructions, mafRange, filedatabase)
!
!       ! Sort the obstructions into order; prune them of repeats, overlaps etc.
!       call PruneObstructions ( obstructions, retVal )
!       if (retval /= 0) return
!
!       call DeallocateL1BData ( taiTime )
!
!    end subroutine
!
!    subroutine PruneObstructions (obstructions, retVal)
!       ! This routine merges overlapping range obstructions and deletes
!       ! wall obstructions inside ranges.  The job is made easier
!       ! by sorting the obstructions into order
!
!       type(Obstruction_T), dimension(:), pointer :: obstructions
!       integer, intent(out) :: retVal
!
!       logical :: foundOne  ! Found at least one
!       integer :: STATUS   ! Flag from allocate
!       integer :: i, j
!       type (Obstruction_T) :: newObs
!
!       ! Executable code
!       retVal = 0
!
!       ! If no obstructions make sure allocate to size zero, not just unassociated pointer
!       if (.not. associated(obstructions)) then
!          allocate(obstructions(0), stat=status)
!          if (status /= 0) then
!             print *, "Out of memory"
!             retVal = 2
!             return
!          end if
!          return
!       end if
!
!       ! Otherwise, do the tidying up
!       outerLoop: do
!          foundOne = .false.
!          i = 0
!          call SortObstructions(obstructions)
!          middleLoop: do
!             i = i + 1
!             if (i >= size(obstructions)) exit
!             j = i
!             innerLoop: do
!                j = j + 1
!                if  (j > size(obstructions)) exit
!                if (all(obstructions((/i,j/))%range)) then
!                   ! --------------------------- ( Range, range )
!                   if ( obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
!                      ! Combine overlapping range obstructions
!                      newObs%range = .true.
!                      newObs%mafs(1) = obstructions(i)%mafs(1)
!                      newObs%mafs(2) = &
!                         & max ( obstructions(i)%mafs(2), obstructions(j)%mafs(2) )
!                      ! Must delete these in order: otherwise
!                      ! if deleted i first where i < j, index would
!                      ! no longer be "j" afterwards
!                      call DeleteObstruction ( obstructions, j, retVal )
!                      if (retVal /= 0) return
!                      call DeleteObstruction ( obstructions, i, retval )
!                      if (retVal /= 0) return
!                      call AddObstructionToDatabase ( obstructions, newObs, retVal )
!                      if (retVal /= 0) return
!                      call SortObstructions ( obstructions )
!                      foundOne = .true.
!                      exit middleLoop
!                   end if
!                else if ( obstructions(i)%range .and. .not. &
!                         obstructions(j)%range ) then
!                   ! --------------------------- ( Range, wall )
!                   if ( obstructions(j)%mafs(1) >= obstructions(i)%mafs(1) .and. &
!                      &  obstructions(j)%mafs(1) <= obstructions(i)%mafs(2) + 1 ) then
!                      ! Delete wall obstruction inside range
!                      call DeleteObstruction ( obstructions, j, retVal )
!                      if (retVal /= 0) return
!                      foundOne = .true.
!                      exit middleLoop
!                   end if
!                else
!                   ! --------------------------- ( Wall, range ) or ( Wall, wall )
!                   ! Becuase the obstructions are in order, we know in the wall, range
!                   ! case that the wall must be at the start of the range, not inside it.
!                   if ( obstructions(i)%mafs(1) == obstructions(j)%mafs(1) ) then
!                      ! Delete wall obstruction at start of a range or at another wall
!                      call DeleteObstruction ( obstructions, i, retVal )
!                      if (retVal /= 0) return
!                      foundOne = .true.
!                      exit middleLoop
!                   end if
!                end if
!
!                ! I'm pretty sure this covers all the possibilities.  It might seem
!                ! not at first glance, but I think the fact that I always re-sort the
!                ! obstructions into order means that the above code does catch everything.
!             end do innerLoop
!          end do middleLoop
!          if (.not. foundOne) exit
!       end do outerLoop
!
!    end subroutine
!
!    subroutine DeleteObstruction ( obstructions, index, retVal )
!       ! Dummy arguments
!       type (Obstruction_T), pointer, dimension(:) :: OBSTRUCTIONS
!       integer, intent(in) :: INDEX
!       integer, intent(out) :: retVal
!
!       ! Local variables
!       type (Obstruction_T), pointer, dimension(:) :: TEMP
!       integer :: STATUS                   ! From allocate
!
!       ! Executable code
!       retVal = 0
!
!       allocate ( temp ( size(obstructions) - 1 ), stat=status )
!       if ( status /= 0 ) then
!          print *, "Out of memory"
!          retVal = 2
!          return
!       end if
!
!       if ( index > 1 ) temp(1:index-1) = obstructions(1:index-1)
!       if ( index < size(obstructions) .and. size(obstructions) > 1 ) &
!          & temp(index:) = obstructions(index+1:)
!
!       deallocate ( obstructions, stat=status )
!       if ( status /= 0 ) then
!          print *, "Invalid memory address"
!          retVal = 3
!          return
!       end if
!
!       obstructions => temp
!
!    end subroutine DeleteObstruction
!
!    subroutine NoteL1BRADChanges (obstructions, mafRange, filedatabase)
!       type(Obstruction_T), dimension(:), pointer :: obstructions
!       type(MAFRange_T) :: mafRange
!       type(MLSFile_T), dimension(:), pointer :: filedatabase
!
!       integer :: noMafs
!
!       ! Won't check if there are no files at all
!       if ( .not. associated(filedatabase) ) return
!       ! Won't check if there are no radiance files
!       if ( .not. any(filedatabase%content == 'l1brad') ) return
!
!       noMafs = mafRange%expanded(2) - mafRange%expanded(1) + 1
!
!       print *, "Not supported NoteL1BRADChanges"
!    end subroutine
!
!    subroutine ChunkDivide_Orbital (mafRange, filedatabase, chunk, retVal)
!       type(MafRange_T), intent(in) :: mafRange
!       type(MLSFile_T), dimension(:), pointer :: filedatabase
!       type(MLSChunk_T), intent(out) :: chunk
!       integer, intent(out) :: retVal
!
!       real(r8), parameter :: homeAccuracy = 3.0 ! Try to hit homeGeodAngle within this
!
!       character(len=10) :: MODNAMESTR     ! Home module name as string
!       type(MLSFile_T), pointer :: l1bFile
!       character(len=name_len) :: maf_start, tp_angle
!       type(L1BData_T) :: taiTime, tpGeodAngle
!       integer :: noMafs, flag, home, noMafsBelowHome, noMafsAtOrAboveHome
!       integer :: m1, m2, mexp1, mexp2, noChunksBelowHome, MAXLENGTH
!       real(r8) :: minAngle, maxAngle, minTime, maxTime, testAngle
!       real(r8) :: angleIncrement
!       integer :: orbit
!       real(r8) :: minV, maxV, homeV
!       real(r8), dimension(:), pointer :: field
!       real(r8) :: boundary
!
!       ! Executable
!       retVal = 0
!
!       ! Read in the data we're going to need
!       call get_string(lit_indices(chunkDivideConfig%homeModule), modNameStr, &
!          strip=.true.)
!       l1bFile => GetMLSFileByType(filedatabase, content='l1boa')
!       if ( .not. associated(L1BFile) ) then
!          print *, "Invalid L1BOA file"
!          retVal = 4
!          return
!       end if
!       maf_start = AssembleL1BQtyName ('MAFStartTimeTAI', l1bFile%hdfVersion, .false.)
!       tp_angle = AssembleL1BQtyName (trim(modNameStr) // '.tpGeodAngle', &
!          l1bFile%hdfVersion, .false.)
!       call ReadL1bData(l1bFile, trim(tp_angle), tpGeodAngle, &
!          noMafs, flag, dontPad=.true.)
!       call smoothOutDroppedMafs(tpGeodAngle%dpField, monotonize=.true.)
!       call ReadL1BData (l1bFile, trim(maf_start), taiTime, &
!          noMafs, flag, dontPad=.true.)
!       call smoothOutDroppedMafs (taiTime%dpField)
!
!       m1 = mafRange%L2Cover(1) + 1
!       m2 = mafRange%L2Cover(2) + 1
!       mexp1 = mafRange%Expanded(1) + 1
!       mexp2 = mafRange%Expanded(2) + 1
!       minAngle = minval (tpGeodAngle%dpField(1,1,m1:m2))
!       maxAngle = maxval (tpGeodAngle%dpfield(1,1,m1:m2))
!       minTime = minval(taiTime%dpField(1,1,m1:m2))
!       maxTime = maxval(taitime%dpfield(1,1,m1:m2))
!
!       ! First try to locate the last MAF before the homeGeodAngle
!       orbit = int(tpGeodAngle%dpField(1,1,m1) / 360.0)
!       if (tpGeodAngle%dpField(1,1,m1) < 0.0) orbit = orbit + 1
!       testAngle = chunkDivideConfig%homeGeodAngle + orbit * 360.0
!       if (chunkDivideConfig%maxLengthFamily == phyq_angle) then
!          angleIncrement = chunkdivideConfig%maxLength
!       else
!          angleIncrement = 360.0
!       end if
!
!       ! In my opinion (paw) here's what the following loop should do:
!       ! Find the 1st MAF within HOMEACCURACY of home_angle
!       ! where home_angle has been corrected for the starting orbit number
!       ! Afterwards, the preceding MAFs must be divided among one or more
!       ! chunks, and the same done with subsequent MAFs
!       !
!       ! Instead what it actually does is
!       ! Find the 1st MAF within HOMEACCURACY of (home_angle + n*angleIncrement)
!       ! where home_angle has been corrected for the starting orbit number
!       ! In effect the home_angle is set only within an unknown number
!       ! of angleIncrements
!       ! While there may be few cases in which they don't do about as well
!       ! let this be a warning
!       do
!          if (testAngle < minAngle) then
!             testAngle = testAngle + angleIncrement
!             cycle
!          end if
!
!          if (testAngle > maxAngle) then
!             print *, 'Unable to establish a home major frame, using the first in your range'
!             home = m1
!             exit
!          end if
!          ! Find MAF which starts before this test angle
!          call Hunt (tpGeodAngle%dpField(1,1,:), testAngle, home, &
!             nearest=.true., allowTopValue=.true.)
!          ! Now if this is close enough, accept it
!          if (abs(tpGeodAngle%dpField(1,1,home) - testAngle) < HomeAccuracy) exit
!          ! Otherwise, keep looking
!          testAngle = testAngle + angleIncrement
!       end do
!
!       ! OK, now we have a home MAF, get a first cut for the chunks
!       ! We work out the chunk ends for each chunk according to how the
!       ! maxLength field is specified.
!       if (chunkDivideConfig%maxLengthFamily == phyq_mafs) then
!          maxLength = nint(ChunkDivideConfig%maxLength)
!          noMafsBelowHome = home - m1
!          noChunksBelowHome = noMafsBelowHome / maxLength
!          if (mod(noMafsBelowHome, maxLength) /= 0) noChunksBelowHome = noChunksBelowHome + 1
!          noMafsAtOrAboveHome = m2 - home + 1
!          ! Subtract one to convert from index in array to index in file
!          chunk%lastMafIndex = home + (1 - noChunksBelowHome) * maxLength - 1
!       else
!          ! For angle and time, they are similar enough we'll just do some stuff
!          ! with pointers to allow us to use common code to sort them out
!          select case (chunkDivideConfig%maxLengthFamily)
!          case (phyq_angle)
!             field => tpGeodAngle%dpField(1,1,:)
!             minV = minAngle
!             maxV = maxAngle
!          case (phyq_time)
!             field => taiTime%dpField(1,1,:)
!             minV = minTime
!             maxV = maxTime
!          end select
!          homeV = field(home)
!
!          noMafsBelowHome = -999
!          noMafsAtOrAboveHome = -999
!          noChunksBelowHome = int((homeV-minV) / ChunkDivideConfig%maxLength)
!          if (homeV > minV) noChunksBelowHome = noChunksBelowHome + 1
!
!          ! When we allow prior overlaps, the first chunk
!          ! sometimes has 1 too many MAFs unless we take extra care
!          if (mexp1 == m1) then
!             boundary = homeV + (1 - noChunksBelowHome) * chunkDivideConfig%maxLength
!          else
!             boundary = homeV + (1 - noChunksBelowHome) * (chunkdivideconfig%maxlength - 1)
!          end if
!          boundary = min(boundary, maxV)
!          boundary = max(boundary, minV)
!
!          call Hunt(field, boundary, chunk%lastMafIndex, start=m1, &
!             allowTopValue=.true., nearest=.true.)
!       end if
!
!       chunk%firstMafIndex = m1
!       chunk%firstMafIndex = min(max(chunk%firstMafIndex, m1), m2)
!       chunk%lastMafIndex = min(max(chunk%lastMafIndex, m1), m2)
!
!       ! Now offset these to the index in the file not the array
!       chunk%firstmafIndex = chunk%firstMafIndex - 1
!       chunk%lastMafIndex = chunk%lastMafIndex - 1
!
!       ! No overlap
!       ! No need to prune chunk
!       ! There is only 1 chunk, how can we deal with obstruction?
!
!       ! Tidy up
!       call DeallocateL1BData(tpGeodAngle)
!       call DeallocateL1BData(taiTime)
!
!    end subroutine
!
!    subroutine smoothOutDroppedMafs (field, wasSmoothed, monotonize)
!       ! detect any fillValues--replace them with nearest neighbor values
!       ! or, optionally, detect and correct any departures from monotone growth
!       real(r8), intent(inout) :: field(:,:,:)
!       logical, dimension(:), optional, intent(out) :: wasSmoothed
!       logical, optional, intent(in) :: monotonize
!
!       integer :: maf, nearest
!       logical :: myMonotonize
!       real(r8) :: lastValue
!
!       myMonotonize = .false.
!       if (present(monotonize)) myMonotonize = monotonize
!       if (present(wasSmoothed)) wasSmoothed = .false.
!       lastValue = field(1,1,1)
!
!       do maf=1, size(field,3)
!          if (myMonotonize) then
!             nearest = max(maf-1,1)
!             if (field(1,1,maf) < lastValue) then
!                if (present(wasSmoothed)) wasSmoothed(maf) = .true.
!                field(:,:,maf) = field(:,:, nearest)
!             else
!                lastValue = field(1,1,maf)
!             end if
!          else if (any(isFillValue(field(:,:,maf)))) then
!             if (present(wasSmoothed)) wasSmoothed(maf) = .true.
!             if (maf == 1) then
!                nearest = findFirst(.not. isFillValue(field(1,1,:)))
!             else
!                nearest = maf - 1
!             end if
!             field(:,:,maf) = field(:,:,nearest)
!          end if
!       end do
!    end subroutine
!
!    subroutine ConvertFlagsToObstructions (valid, obstructions, &
!       retVal, mafRange, obstructionType)
!       ! This routine takes an array of logicals indicating good/bad data
!       ! and converts it into obstruction information.
!       logical, dimension(:), intent(in) :: valid
!       type(Obstruction_T), dimension(:), pointer :: obstructions
!       integer, intent(out) :: retVal
!       integer, dimension(:), intent(in), optional :: mafRange
!       character(len=*), intent(in), optional :: obstructionType
!
!       logical :: LASTONEVALID           ! Flag
!       integer :: MAF                    ! Loop counter
!       type (Obstruction_T) :: NEWOBSTRUCTION ! In progress
!       character(len=64)    :: obstructionTrigger
!       integer :: OFFSET                 ! MAF index offset
!
!       ! Executables
!       retVal = 0
!       lastOneValid = .true.
!       offset = 0
!       if (present(mafRange) ) offset = mafRange(1)
!       obstructionTrigger = 'maf where transition from bad to good made obstruction'
!       if (present(obstructionType)) &
!          obstructionTrigger = 'maf where transition from bad to good ' &
!             // trim(obstructionType) // 'made obstruction'
!
!       do maf = 1, size(valid)
!          if (valid(maf) .neqv. lastOneValid) then
!             ! A transition either from good to bad or bad to good
!             if (.not. valid(maf)) then
!                ! From good to bad
!                newObstruction%range = .true.
!                newObstruction%mafs(1) = maf - 1 + offset
!             else
!                newObstruction%mafs(2) = maf - 2 + offset
!                call AddObstructionToDatabase ( obstructions, newObstruction, retVal )
!                if (retVal /= 0) return
!             end if
!          end if
!
!          lastOneValid = valid(maf)
!       end do
!
!       ! Make sure any range at the end gets added
!       if (.not. lastOneValid) then
!          newObstruction%mafs(2) = size(valid) - 1 + offset
!          call AddObstructionToDatabase ( obstructions, newObstruction, retVal )
!          if (retVal /= 0) return
!       end if
!    end subroutine
!
!    subroutine AddObstructionToDatabase ( database, item, retVal )
!
!       ! Dummy arguments
!       type (Obstruction_T), dimension(:), pointer :: DATABASE
!       type (Obstruction_T), intent(in) :: ITEM
!       integer, intent(out) :: retVal
!
!       ! Local variables
!       type (Obstruction_T), dimension(:), pointer :: TEMPDATABASE
!       integer :: newSize, status
!
!       ! Executables
!       retVal = 0
!       if ( associated(database) ) then ! tree_checker prevents duplicate names
!          newSize=SIZE(database)+1
!       else
!          newSize = 1
!       end if
!
!       allocate(tempDatabase(newSize),STAT=status)
!       if (status /= 0) then
!          print *, "Out of memory in " // ModuleName
!          retVal = 2
!          return
!       end if
!
!       if ( newSize>1 ) tempDatabase(1:newSize-1) = database
!       if (associated(database)) deallocate ( database, stat=status )
!       if (status /= 0) then
!          print *, "Invalid memory address in " // ModuleName
!          retVal = 3
!          return
!       end if
!       database => tempDatabase
!       database(newSize) = item
!
!    end subroutine AddObstructionToDatabase
!
!    subroutine SortObstructions (obstructions)
!       ! Sort the obstructions into order of increasing
!       ! mafs(1) (start/wall MAF index)
!       type(Obstruction_T), dimension(:), intent(inout) :: obstructions
!
!       type(Obstruction_T) :: temp
!       integer :: i
!       integer, dimension(1) :: toswap
!
!       do i=1, size(obstructions) - 1
!          toSwap = minloc(obstructions(i:)%mafs(1)) + (/i-1/)
!          if (toswap(1) /= i) then
!             temp = obstructions(i)
!             obstructions(i) = obstructions(toswap(1))
!             obstructions(toswap(1)) = temp
!          end if
!       end do
!    end subroutine

end module
